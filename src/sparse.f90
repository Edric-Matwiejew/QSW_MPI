!   QSW_MPI -  A package for parallel Quantum Stochastic Walk simulation.


!   Copyright (C) 2019 Edric Matwiejew
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.

!
!   Module: Sparse_Operations
!
!> @brief MPI parallel sparse BLAS operations.
!
module Sparse

    use iso_fortran_env, only: int64, error_unit
    use :: ISO_Precisions
    use :: MPI

    implicit none

    private :: Neighbour_SpMV

    !> @brief Compressed sparse rows (CSR) complex matrix derived type.
    !
    !> @note Supply input Hamiltonians and Lindblad operators in ascending column
    !> order within each row. Distributed matrix-vector products accept unsorted rows.

    type, public  :: CSR

        integer :: rows
        integer :: columns
        character(len=2) :: structure
        integer, dimension(:), pointer :: row_starts => null()
        integer, dimension(:), pointer :: col_indexes => null()
        complex(dp), dimension(:), pointer :: values => null()

        integer, dimension(:), pointer :: local_col_inds => null()
        integer, dimension(:), pointer :: RHS_send_inds => null()
        integer, dimension(:), pointer :: num_send_inds => null()
        integer, dimension(:), pointer :: send_disps => null()
        integer, dimension(:), pointer :: num_rec_inds => null()
        integer, dimension(:), pointer :: rec_disps => null()

        integer :: halo_comm = MPI_COMM_NULL
        integer, allocatable :: halo_send_counts(:), halo_send_disps(:), halo_rec_counts(:), halo_rec_disps(:)
        integer, allocatable :: first_remote(:)

    end type CSR

    !
    !>  @brief Compressed Co-ordinates spase matrix derived type.
    !

    type, public :: COO

        integer :: rows
        integer :: columns
        integer :: nnz = 0
        integer, dimension(:), pointer :: row_indexes => null()
        integer, dimension(:), pointer :: col_indexes => null()
        complex(dp), dimension(:), pointer :: values => null()

    end type COO

    !
    !>  @brief Sparse matrix accumulator, used in COO sum.
    !

    type, private :: SPA
        complex(dp), dimension(:), allocatable :: w
        integer, dimension(:), allocatable :: b
        integer, dimension(:), allocatable :: LS
        integer :: LS_ind = 1
    end type SPA

    real(dp), private, parameter :: spa_tolerance = epsilon(1.0_dp)

    contains

    !
    !   Function: Kronecker_Delta
    !
    !>  @brief Kronecker delta function of integers i and j.
    !
    !>  @details This function is equal to 1 if i = j and 0 otherwise.

    function Kronecker_Delta(i, j)

        integer :: Kronecker_Delta
        integer, intent(in) :: i, j

        Kronecker_Delta = merge(1, 0, i == j)

    end function Kronecker_Delta

    !
    !   Subroutine: Prefix_Sum
    !
    !>  @brief Prefix sum of an integer array.
    !
    !>  @details Perfroms an in-place prefix (or cumulative) sum on an integer
    !>  array.
    !>
    !>> Prefix_Sum([1,2,3]) -> [1,3,6]

    subroutine Prefix_Sum(array)

            integer, dimension(:), intent(inout) :: array

            integer :: i

            do i = 2, size(array, 1)
                array(i) = array(i - 1) + array(i)
            enddo

    end subroutine Prefix_Sum

    !
    ! COO (Compressed Co-ordiantes) subroutines.
    !

    !> @brief Return the first one-based CSR row pointer on this rank.
    !> @details Local row pointers start at one. Global row pointers include
    !> the number of entries on all preceding ranks.
    integer function CSR_Element_Start(local_nnz, rank, comm, local_pointers) result(first)
        integer(int64), intent(in) :: local_nnz
        integer, intent(in) :: rank, comm
        logical, optional, intent(in) :: local_pointers
        logical :: global_pointers
        integer(int64) :: offset
        integer :: ierr

        global_pointers = .true.
        if (present(local_pointers)) global_pointers = .not. local_pointers
        offset = 0_int64
        if (global_pointers) then
            call MPI_Exscan(local_nnz, offset, 1, MPI_INTEGER8, MPI_SUM, comm, ierr)
            if (rank == 0) offset = 0_int64
        endif
        if (offset + local_nnz >= huge(first)) then
            write(error_unit,*) 'CSR row pointer exceeds the 32-bit range on rank', rank
            call MPI_Abort(comm, 1, ierr)
        endif
        first = int(offset) + 1
    end function CSR_Element_Start

    !> @brief Convert a COO matrix to a CSR matrix.

    subroutine COO_to_CSR(A, B, partition_table, MPI_communicator, local_pointers)

        type(COO), intent(inout) :: A
        type(CSR), intent(inout) :: B
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(in) :: MPI_communicator

        logical, optional, intent(in) :: local_pointers
        integer :: lower_bound, upper_bound

        integer :: rank, ierr

        integer :: i

        call mpi_comm_rank(MPI_communicator, rank, ierr)

        lower_bound = partition_table(rank + 1)
        upper_bound = partition_table(rank + 2) - 1

        B%rows = A%rows
        B%columns = A%columns

        allocate(B%row_starts(lower_bound:upper_bound + 1))
        allocate(B%col_indexes(A%nnz))
        allocate(B%values(A%nnz))

        B%row_starts = 0

        do i = 1, A%nnz
            B%values(i) = A%values(i)
            B%col_indexes(i) = A%col_indexes(i)
            B%row_starts(A%row_indexes(i) + 1) = B%row_starts(A%row_indexes(i) + 1) + 1
        enddo

        B%row_starts(lower_bound) = CSR_Element_Start(int(A%nnz, int64), rank, MPI_communicator, local_pointers)
        call prefix_sum(B%row_starts)

    end subroutine COO_to_CSR

    !
    !> @brief Sparse matrix gather, used in COO sum.
    !

    subroutine Gather_SPA(row_spa, A, current_row, comm)

        type(SPA), intent(inout) :: row_spa
        type(COO), intent(inout) :: A
        integer, intent(in) :: current_row, comm
        integer(int64) :: nnz
        integer :: i

        nnz = A%nnz
        do i = 1, row_spa%LS_ind - 1
            if (abs(row_spa%w(row_spa%LS(i))) > spa_tolerance) then
                nnz = nnz + 1
                A%col_indexes(nnz) = row_spa%LS(i)
                A%values(nnz) = row_spa%w(row_spa%LS(i))
            endif
        enddo

        call Check_Local_Nnz(nnz, comm)
        A%row_indexes(A%nnz + 1:nnz) = current_row
        A%nnz = int(nnz)
        row_spa%LS_ind = 1

    end subroutine Gather_SPA

    subroutine Check_Local_Nnz(nnz, comm)
        integer(int64), intent(in) :: nnz
        integer, intent(in) :: comm
        integer :: ierr

        if (nnz >= huge(0)) then
            write(error_unit,*) 'Local sparse matrix exceeds the one-based index range:', nnz
            call MPI_Abort(comm, 1, ierr)
        endif
    end subroutine Check_Local_Nnz

    !
    !> @brief Scatter elements of a sparse matrix, used in COO matrix sum.
    !

    subroutine Scatter_SPA(row_spa, value, col_index, current_row)

        type(SPA), intent(inout) :: row_spa
        complex(dp), intent(in) :: value
        integer, intent(in) :: col_index
        integer, intent(in) :: current_row

        if (row_spa%b(col_index) < current_row) then
            row_spa%w(col_index) = value
            row_spa%b(col_index) = current_row
            row_spa%LS(row_spa%LS_ind) = col_index
            row_spa%LS_ind = row_spa%LS_ind + 1
        else
            row_spa%w(col_index) = row_spa%w(col_index) + value
        endif

    end subroutine Scatter_SPA

    !
    !> @brief Allocate COO arrays prior to a matrix sum.
    !

    subroutine COO_Allocate_Sum(A, B, local_rows, C)

        type(COO), intent(in) :: A
        type(COO), intent(in) :: B
        integer, intent(in) :: local_rows
        type(COO), intent(inout) :: C

        integer(int64) :: estimated_nnz

        estimated_nnz = min(int(A%nnz, int64) + B%nnz, int(A%columns, int64)*local_rows)

        allocate(C%row_indexes(estimated_nnz))
        allocate(C%col_indexes(estimated_nnz))
        allocate(C%values(estimated_nnz))

    end subroutine COO_Allocate_Sum

    !> @brief Deallocate a COO matrix.

    subroutine COO_Deallocate(A)

        type(COO), intent(inout) :: A

        deallocate(A%row_indexes)
        deallocate(A%col_indexes)
        deallocate(A%values)

        A%nnz = 0
        A%rows = 0
        A%columns = 0

    end subroutine COO_Deallocate

    !> @brief Copy COO matrix A to B, deallocate A.

    subroutine COO_Move(A,B)

        type(COO), intent(inout) :: A
        type(COO), intent(inout) :: B

        allocate(B%row_indexes(A%nnz))
        allocate(B%col_indexes(A%nnz))
        allocate(B%values(A%nnz))

        B%rows = A%rows
        B%columns = A%columns
        B%nnz = A%nnz
        B%row_indexes = A%row_indexes(1:A%nnz)
        B%col_indexes = A%col_indexes(1:A%nnz)
        B%values = A%values(1:A%nnz)

        call COO_Deallocate(A)

    end subroutine COO_Move

    !> @brief Sum distributed COO matrices, A + B = C.

    subroutine COO_Sum(A, B, C, partition_table, MPI_communicator)

        type(COO), intent(inout) :: A
        type(COO), intent(inout) :: B
        type(COO), intent(inout) :: C
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(in) :: MPI_communicator

        type(SPA) :: row_spa

        integer :: ka, kb
        integer :: i

        ! MPI
        integer :: rank
        integer :: ierr

        call MPI_comm_rank(MPI_communicator, rank, ierr)

        allocate(row_spa%w(A%columns))
        allocate(row_spa%b(A%columns))
        allocate(row_spa%LS(A%columns))

        row_spa%b = 0
        row_spa%LS = 0

        ka = 1
        kb = 1

        do i = partition_table(rank + 1), partition_table(rank + 2) - 1
            do while (ka <= A%nnz)
                if (A%row_indexes(ka) /= i) exit
                call scatter_spa(row_spa, A%values(ka), A%col_indexes(ka), i)
                ka = ka + 1
            enddo
            do while (kb <= B%nnz)
                if (B%row_indexes(kb) /= i) exit
                call scatter_spa(row_spa, B%values(kb), B%col_indexes(kb), i)
                kb = kb + 1
            enddo
            call Gather_SPA(row_spa, C, i, MPI_communicator)
            row_spa%LS_ind = 1
        enddo

        C%rows = A%rows
        C%columns = A%columns

    end subroutine COO_Sum

    !> @brief Add the Kronecker product of conjg(L) and L to the local rows of B.
    subroutine Add_Jump_Product(L, B, first_row, last_row, C, MPI_communicator)
        type(CSR), intent(in) :: L
        type(COO), intent(in) :: B
        integer, intent(in) :: first_row, last_row, MPI_communicator
        type(COO), intent(out) :: C

        type :: entry_block
            integer, allocatable :: cols(:)
            complex(dp), allocatable :: values(:)
            type(entry_block), pointer :: next => null()
        end type
        type(entry_block), pointer :: head, tail, block
        type(SPA) :: row
        integer, allocatable :: row_ends(:)
        integer(int64) :: local_nnz
        integer :: i, outer, inner, j, k, col, b_pos, capacity, used, offset, n
        complex(dp) :: factor

        C%rows = L%rows**2
        C%columns = L%columns**2
        allocate(row%w(C%columns), row%b(C%columns), row%LS(C%columns))
        allocate(row_ends(first_row - 1:last_row))
        row%b = 0
        row_ends(first_row - 1) = 0
        local_nnz = 0
        b_pos = 1
        capacity = int(min(65536_int64, int(C%columns, int64)*(last_row - first_row + 1)))
        used = capacity
        nullify(head, tail)
        do i = first_row, last_row
            outer = (i - 1) / L%rows + 1
            inner = mod(i - 1, L%rows) + 1
            do j = L%row_starts(outer), L%row_starts(outer + 1) - 1
                col = (L%col_indexes(j) - 1) * L%columns
                factor = conjg(L%values(j))
                do k = L%row_starts(inner), L%row_starts(inner + 1) - 1
                    call Scatter_SPA(row, factor * L%values(k), col + L%col_indexes(k), i)
                enddo
            enddo
            do while (b_pos <= B%nnz)
                if (B%row_indexes(b_pos) /= i) exit
                call Scatter_SPA(row, B%values(b_pos), B%col_indexes(b_pos), i)
                b_pos = b_pos + 1
            enddo
            ! Apply the tolerance to the sum of all contributions to each entry.
            do j = 1, row%LS_ind - 1
                col = row%LS(j)
                if (abs(row%w(col)) > spa_tolerance) then
                    if (used == capacity) then
                        allocate(block)
                        allocate(block%cols(capacity), block%values(capacity))
                        if (associated(tail)) then
                            tail%next => block
                        else
                            head => block
                        endif
                        tail => block
                        used = 0
                    endif
                    used = used + 1
                    tail%cols(used) = col
                    tail%values(used) = row%w(col)
                    local_nnz = local_nnz + 1
                endif
            enddo
            call Check_Local_Nnz(local_nnz, MPI_communicator)
            row_ends(i) = int(local_nnz)
            row%LS_ind = 1
        enddo

        C%nnz = int(local_nnz)
        allocate(C%row_indexes(C%nnz), C%col_indexes(C%nnz), C%values(C%nnz))
        do i = first_row, last_row
            C%row_indexes(row_ends(i-1) + 1:row_ends(i)) = i
        enddo
        offset = 0
        do while (associated(head))
            n = min(capacity, C%nnz - offset)
            C%col_indexes(offset + 1:offset + n) = head%cols(:n)
            C%values(offset + 1:offset + n) = head%values(:n)
            offset = offset + n
            block => head%next
            deallocate(head)
            head => block
        enddo
    end subroutine Add_Jump_Product

    !> @brief Transpose a CSR matrix.

    subroutine CSR_Transpose(A, A_D)

        type(CSR), intent(in) :: A
        type(CSR), intent(out) :: A_D

        integer, allocatable :: next(:)
        integer :: i, j, col, dest, nnz

        A_D%rows = A%columns
        A_D%columns = A%rows
        allocate(A_D%row_starts(A_D%rows + 1))
        nnz = A%row_starts(A%rows + 1) - A%row_starts(1)
        allocate(A_D%col_indexes(nnz), A_D%values(nnz))
        A_D%row_starts = 0
        A_D%row_starts(1) = 1
        do j = A%row_starts(1), A%row_starts(A%rows + 1) - 1
            col = A%col_indexes(j)
            A_D%row_starts(col + 1) = A_D%row_starts(col + 1) + 1
        enddo
        call Prefix_Sum(A_D%row_starts)
        next = A_D%row_starts(:A_D%rows)
        do i = 1, A%rows
            do j = A%row_starts(i), A%row_starts(i + 1) - 1
                col = A%col_indexes(j)
                dest = next(col)
                A_D%col_indexes(dest) = i
                A_D%values(dest) = A%values(j)
                next(col) = dest + 1
            enddo
        enddo

    end subroutine CSR_Transpose

    !
    ! Subroutine: Generate_Partition_Table
    !
    !> @brief Determines a matrix/vector partition scheme.
    !>
    !> @details A given number of rows is divided evenly (if possible) over an
    !> MPI communicator of size *N* nodes. If even division is not possible the
    !> remainder is given to the node of highest rank.
    !>
    !> Array *partition_table* is of size N + 1 and follows a format similar to
    !> that of the *row_start* array in the *CSR* data type:
    !>
    !>> (starting row for rank 0, starting row for rank 1, ..., rows + 1)
    !>
    !> As such the number of rows at a node is given by:
    !>
    !>> partition_table(rank + 2) - partition_table(rank + 1)
    !>
    !> @note Must be called from within an MPI instance.
    !
    !> @todo send row data to MPI communicator from root process.

    subroutine Generate_Partition_Table(rows, partition_table, MPI_communicator)

        integer, intent(in) :: rows !< @param Number of matrix/vector rows.
        integer, intent(out), dimension(:), allocatable :: partition_table !< @param Start/Finish of partitions.
        integer, intent(in) :: MPI_communicator !< @param MPI communicator over which to distribute.

        integer :: i

        ! MPI environment
        integer :: flock ! size of MPI communicator.
        integer :: ierr

        call MPI_comm_size(mpi_communicator, flock, ierr)

        allocate(partition_table(flock + 1))

        do i = 1, flock + 1
            partition_table(i) = int((i - 1)*int(rows, int64)/flock) + 1
        enddo

    end subroutine Generate_Partition_Table

    !> @brief Distribute complex array v row-wise over an MPI communicator.

    subroutine Distribute_Dense_Vector( v, &
                                        partition_table, &
                                        root, &
                                        v_local, &
                                        MPI_communicator)

        complex(dp), dimension(:), intent(in) :: v !< @param Array to distribute.
        integer, dimension(:), intent(in) :: partition_table !< @param MPI communicator partitioning scheme.
        integer, intent(in) :: root !< @param Node from which the array is distributed.
        complex(dp), dimension(:), allocatable, intent(out) :: v_local !< @param Local array partition.
        integer, intent(in) :: mpi_communicator !< @param MPI communicator over which to distribute.

        integer :: lb_send, ub_send, lb_recv, ub_recv

        integer :: i

        ! MPI ENVIRONMENT
        integer :: flock, rank
        integer, dimension(:), allocatable :: requests_out
        integer :: buff_size
        integer :: ierr

        call MPI_comm_size(MPI_communicator, flock, ierr)
        call MPI_comm_rank(MPI_communicator, rank, ierr)

        allocate(requests_out(0:flock - 1))

        if (rank == root) then
            do i = 0, flock - 1
                if (i /= root) then
                    lb_send = partition_table(i + 1)
                    ub_send = partition_table(i + 2) - 1
                    buff_size = ub_send - lb_send + 1

                    call MPI_Isend( v(lb_send:ub_send), buff_size, &
                               MPI_double_complex, i, &
                               rank*10 + i, MPI_communicator, &
                               requests_out(i), ierr)
                endif
            enddo
        endif


        lb_recv = partition_table(rank + 1)
        ub_recv = partition_table(rank + 2) - 1

        allocate(v_local(lb_recv:ub_recv))
        if (rank == root) then
            v_local = v(lb_recv:ub_recv)
        else
            buff_size = ub_recv - lb_recv + 1
            call MPI_recv(  v_local(lb_recv:ub_recv), buff_size, &
                            MPI_double_complex, root, &
                            root*10 + rank, MPI_communicator, &
                            MPI_status_ignore, ierr)
        endif

        call mpi_barrier(MPI_communicator, ierr)

    end subroutine Distribute_Dense_Vector

    !
    !   subroutine: Distribute_CSR_Matrix
    !
    !> @brief Distribute complex CSR matrix over an MPI communicator.
    !>
    !> @details CSR matrix is partitioned following a scheme given by
    !> the *partition_table* array, see @ref generate_partition_table for more.
    !>
    !> @note Must be called from within an MPI instance.

    subroutine Distribute_CSR_Matrix(   A, &
                                        partition_table, &
                                        root, &
                                        A_local, &
                                        MPI_communicator)

        type(CSR), intent(in) :: A  !< @param CSR matrix to distribute.
        integer, dimension(:), intent(in) :: partition_table    !< @param MPI communicator paritioning scheme.
        integer, intent(in) :: root !< @param Node from which A is distributed.
        type(CSR), intent(out) :: A_local !< @param Local CSR parition.
        integer, intent(in) :: MPI_communicator !< @param MPU communicator over which to distribute.

        integer, dimension(:), allocatable :: block_lens, disps
        integer, dimension(:), allocatable :: block_lens_vals, disps_vals

        integer :: lb, ub, lb_vals, ub_vals

        integer :: i

        !MPI ENVIRONMENT
        integer :: rank, flock
        integer :: ierr

        call MPI_comm_size(MPI_communicator, flock, ierr)
        call MPI_comm_rank(MPI_communicator, rank, ierr)


        if (rank == root) then
            A_local%rows = A%rows
            A_local%columns = A%columns
        endif

        call MPI_bcast( A_local%rows, &
                        1, &
                        MPI_integer, &
                        root, &
                        MPI_communicator, &
                        ierr)

        call MPI_bcast( A_local%columns, &
                        1, &
                        MPI_integer, &
                        root, &
                        MPI_communicator, &
                        ierr)

        allocate(block_lens(flock))
        allocate(disps(flock))

        do i = 1, flock
            block_lens(i) = partition_table(i + 1) - partition_table(i) + 1
        enddo

        disps = 0

        do i = 2, flock
            disps(i) = disps(i - 1) + block_lens(i - 1) - 1
        enddo

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2)

        allocate(A_local%row_starts(lb:ub))

        call MPI_scatterv(  A%row_starts, &
                            block_lens, &
                            disps, &
                            MPI_integer, &
                            A_local%row_starts, &
                            block_lens(rank + 1), &
                            MPI_integer, &
                            root, &
                            MPI_communicator, &
                            ierr)

        lb_vals = A_local%row_starts(lb)
        ub_vals = A_local%row_starts(ub) - 1

        allocate(A_local%col_indexes(lb_vals:ub_vals))
        allocate(A_local%values(lb_vals:ub_vals))

        allocate(block_lens_vals(flock))
        allocate(disps_vals(flock))

        if (rank == root) then

            do i = 1, flock

                block_lens_vals(i) = A%row_starts(partition_table(i + 1)) &
                                        - A%row_starts(partition_table(i))

            enddo

        endif

        call MPI_bcast( block_lens_vals, &
                        flock, &
                        MPI_integer, &
                        root, &
                        MPI_communicator, &
                        ierr)

        disps_vals = 0

        do i = 2, flock

            disps_vals(i) = disps_vals(i - 1) + block_lens_vals(i - 1)

        enddo


        call MPI_scatterv(  A%col_indexes, &
                            block_lens_vals, &
                            disps_vals, &
                            MPI_integer, &
                            A_local%col_indexes, &
                            block_lens_vals(rank + 1), &
                            MPI_integer, &
                            root, &
                            MPI_communicator, &
                            ierr)

        call MPI_scatterv(  A%values, &
                            block_lens_vals, &
                            disps_vals, &
                            MPI_double_complex, &
                            A_local%values, &
                            block_lens_vals(rank + 1), &
                            MPI_double_complex, &
                            root, &
                            MPI_communicator, &
                            ierr)

    end subroutine Distribute_CSR_Matrix

    !
    !   subroutine: Gather_Dense_Vector
    !
    !> @brief Gather a 1D double precision complex array distributed over an MPI
    !> communicator.
    !>
    !> @details The array is partitioned following a scheme given by
    !> the *partition_table* array, see @ref generate_partition_table for more.
    !>
    !> @note Must be called from within an MPI instance.
    !>
    !> @todo rewriting Distribute_dense_vector with MPI_scatterv

    subroutine Gather_Dense_Vector( u_local, &
                                    partition_table, &
                                    root, &
                                    u, &
                                    MPI_communicator)

        complex(dp), dimension(:), intent(in) :: u_local !< @param Array to gather.
        integer, dimension(:), intent(in) :: partition_table !< @param MPI communicator partitioning scheme.
        integer, intent(in) :: root !< @param Node to which the array is gathered.
        complex(dp), dimension(:), intent(out) :: u !< @param Gathered array.
        integer, intent(in) :: mpi_communicator !< @param MPI communicator over which to gather.

        integer, dimension(:), allocatable :: block_lens, disps
        integer :: lb, ub
        integer :: i

        !MPI ENVIRONMENT
        integer :: rank, flock
        integer :: ierr

        call MPI_comm_size(MPI_communicator, flock, ierr)
        call MPI_comm_rank(MPI_communicator, rank, ierr)

        allocate(block_lens(flock))
        allocate(disps(flock))

        do i = 1, flock

            block_lens(i) = partition_table(i + 1) - partition_table(i)

        enddo

        disps = 0

        do i = 2, flock

            disps(i) = disps(i - 1) + block_lens(i - 1)

        enddo

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2)

        call MPI_gatherv(   u_local, &
                            block_lens(rank + 1), &
                            MPI_double_complex, &
                            u, &
                            block_lens, &
                            disps, &
                            MPI_double_complex, &
                            root, &
                            MPI_communicator, &
                            ierr)

    end subroutine Gather_Dense_Vector

    !> @brief Distribute a 2 x 2 complex matrix row-wise over a MPI communicator.

    subroutine Scatter_Dense_Matrix(   A, &
                                        partition_table, &
                                        root, &
                                        A_local, &
                                        MPI_communicator)

        complex(dp), dimension(:,:), intent(in) :: A
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(in) :: root
        complex(dp), dimension(:,:), allocatable, intent(out) :: A_local
        integer, intent(in) :: MPI_communicator

        integer :: rec_lb, rec_ub
        integer :: send_lb, send_ub
        integer :: rec_row, rec_col
        integer :: send_row, send_col
        integer :: i

        !MPI ENVIRONMENT
        integer, dimension(:), allocatable :: requests_out
        integer, dimension(MPI_status_size) :: status
        integer :: rank, flock
        integer :: ierr

        call MPI_comm_size(MPI_communicator, flock, ierr)
        call MPI_comm_rank(MPI_communicator, rank, ierr)

        rec_lb = partition_table(rank + 1)
        rec_ub = partition_table(rank + 2) - 1

        rec_row = rec_ub - rec_lb + 1
        rec_col = size(A,2)

        CALL mpi_bcast( rec_col, &
                        1, &
                        MPI_integer, &
                        root, &
                        MPI_communicator, &
                        ierr)

        allocate(requests_out(0:flock-1))


        if (rank == root) then

            send_col = rec_col

            do i = 0, flock - 1

                if (i /= root) then

                send_lb = partition_table(i + 1)
                send_ub = partition_table(i + 2) - 1

                send_row = partition_table(i + 2) - partition_table(i + 1)

                call MPI_send( A(send_lb:send_ub, 1:rec_col), &
                                send_row*send_col, &
                                MPI_double_complex, &
                                i, &
                                i*10, &
                                MPI_communicator, &
                                ierr)

            endif

            enddo

        else

            allocate(A_local(rec_lb:rec_ub, rec_col))

            call MPI_Recv(  A_local(rec_lb:rec_ub, 1:rec_col), &
                            rec_row*rec_col, &
                            MPI_double_complex, &
                            root, &
                            rank*10, &
                            MPI_communicator, &
                            status, &
                            ierr)

        endif

        if (rank == root) then

            allocate(A_local(rec_lb:rec_ub, rec_col))

            A_local = A(rec_lb:rec_ub, :)

        endif

        call mpi_barrier(MPI_communicator, ierr)

    end subroutine Scatter_Dense_Matrix

    !> @brief Gather a row-wise distributed 2 x 2 complex array to MPI rank root.

    subroutine Gather_Dense_Matrix( A_local, &
                                    partition_table, &
                                    root, &
                                    A, &
                                    MPI_communicator)

        complex(dp), dimension(:,:), intent(in) :: A_local
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(in) :: root
        complex(dp), dimension(:,:), intent(out) :: A
        integer, intent(in) :: MPI_communicator

        integer :: rec_lb, rec_ub
        integer :: send_row, send_col
        integer :: rec_row, rec_col
        logical, dimension(:), allocatable :: received
        integer :: i

        !MPI ENVIRONMENT
        integer, dimension(MPI_status_size) :: status
        integer :: rank, flock
        integer :: ierr

        call MPI_comm_size(MPI_communicator, flock, ierr)
        call MPI_comm_rank(MPI_communicator, rank, ierr)

        if (rank /= root) then

            send_row = partition_table(rank + 2) - partition_table(rank + 1)
            send_col = size(A_local, 2)

            call MPI_send(  A_local, &
                            send_row*send_col, &
                            MPI_double_complex, &
                            root, &
                            rank*10, &
                            MPI_communicator, &
                            ierr)

        else

            if (flock > 1) then

                rec_col = size(A_local, 2)

                rec_lb = partition_table(rank + 1)
                rec_ub = partition_table(rank + 2) - 1

                A(rec_lb:rec_ub, :) = A_local

                allocate(received(0:flock - 1))

                received = .false.

                received(rank) = .true.

                do

                    do i = 0, flock - 1

                        if (.not. received(i)) then

                            call MPI_Iprobe(    i, &
                                                i*10, &
                                                MPI_communicator, &
                                                received(i), &
                                                status, &
                                                ierr)

                            if (received (i)) then

                                rec_lb = partition_table(i + 1)
                                rec_ub = partition_table(i + 2) - 1

                                rec_row = rec_ub - rec_lb + 1

                                call MPI_Recv(  A(rec_lb:rec_ub, :), &
                                                rec_row*rec_col, &
                                                MPI_double_complex, &
                                                i, &
                                                i*10, &
                                                MPI_communicator, &
                                                status, &
                                                ierr)
                            endif

                        endif

                    enddo

                    if(all(received)) exit

                enddo

            else

                A = A_local

            endif

        endif

    end subroutine Gather_Dense_Matrix

    !> @brief Merge sort a CSR matrix A by an array of column indexes,
    !> used to form the conjugate transpose of A.


    subroutine Merge_Dagger(column_indexes, &
                            row_indexes, &
                            values, &
                            start, &
                            mid, &
                            finish)

        integer, intent(inout), dimension(:) :: column_indexes
        integer, intent(inout), dimension(:) :: row_indexes
        complex(dp), intent(inout), dimension(:) :: values
        integer, intent(in) :: start
        integer, intent(in) :: mid
        integer, intent(in) :: finish

        integer, dimension(:), allocatable :: col_ind_temp
        integer, dimension(:), allocatable :: row_ind_temp
        complex(dp), dimension(:), allocatable :: val_temp
        integer :: i, j, k

        allocate(col_ind_temp(finish - start + 1))
        allocate(row_ind_temp(finish - start + 1))
        allocate(val_temp(finish - start + 1))

        i = start
        j = mid + 1
        k = 1

        do while (i <= mid .and. j <= finish)

            if (column_indexes(i) <= column_indexes(j)) then
                col_ind_temp(k) = column_indexes(i)
                row_ind_temp(k) = row_indexes(i)
                val_temp(k) = values(i)
                k = k + 1
                i = i + 1
            else
                col_ind_temp(k) = column_indexes(j)
                row_ind_temp(k) = row_indexes(j)
                val_temp(k) = values(j)
                k = k + 1
                j = j+ 1
            endif

        enddo

        do while (i <= mid)
            col_ind_temp(k) = column_indexes(i)
            row_ind_temp(k) = row_indexes(i)
            val_temp(k) = values(i)
            k = k + 1
            i = i + 1
        enddo

        do while (j <= finish)
            col_ind_temp(k) = column_indexes(j)
            row_ind_temp(k) = row_indexes(j)
            val_temp(k) = values(j)
            k = k + 1
            j = j + 1
        enddo

        do i = start, finish
            column_indexes(i) = col_ind_temp(i - start + 1)
            row_indexes(i) = row_ind_temp(i - start + 1)
            values(i) = val_temp(i - start + 1)
        enddo

    end subroutine Merge_Dagger
    
    !> @brief Insertion sort a CSR matrix A by an array of column indexes,
    !> used to form the conjugate transpose of A.

    subroutine Insertion_Sort_Dagger(   column_indexes, &
                                        row_indexes, &
                                        values)

        integer, intent(inout), dimension(:) :: column_indexes
        integer, intent(inout), dimension(:) :: row_indexes
        complex(dp), intent(inout), dimension(:) :: values

        integer :: col_ind_temp
        integer :: row_ind_temp
        complex(dp) :: val_temp

        integer :: i, j

        do i = 2, size(column_indexes)

            col_ind_temp = column_indexes(i)
            row_ind_temp = row_indexes(i)
            val_temp = values(i)

            j = i - 1

            do while (j >= 1)

                if (column_indexes(j) <= col_ind_temp) exit
                    column_indexes(j + 1) = column_indexes(j)
                    row_indexes(j + 1) = row_indexes(j)
                    values(j + 1) = values(j)
                    j = j - 1
            enddo
            column_indexes(j + 1) = col_ind_temp
            row_indexes(j + 1) = row_ind_temp
            values(j + 1) = val_temp

        enddo

    end subroutine Insertion_Sort_Dagger

    !> @brief Merge sort a CSR matrix A by an array of column indexes,
    !> used to form the conjugate transpose of A.

    recursive subroutine Merge_Sort_Dagger( column_indexes, &
                                            row_indexes, &
                                            values, &
                                            start, &
                                            finish)

        integer, intent(inout), dimension(:) :: column_indexes
        integer, intent(inout), dimension(:) :: row_indexes
        complex(dp), intent(inout), dimension(:) :: values
        integer, intent(in) :: start
        integer, intent(in) :: finish

        integer :: mid

        if (start < finish) then
            if (finish - start >= 512) then

                mid = start + (finish - start) / 2

                call Merge_Sort_Dagger( column_indexes, &
                                        row_indexes, &
                                        values, &
                                        start, &
                                        mid)

                call Merge_Sort_Dagger( column_indexes, &
                                        row_indexes, &
                                        values, &
                                        mid + 1, &
                                        finish)

                call Merge_Dagger(  column_indexes, &
                                    row_indexes, &
                                    values, &
                                    start, &
                                    mid, &
                                    finish)

            else
                call insertion_sort_Dagger( column_indexes(start:finish), &
                                            row_indexes(start:finish), &
                                            values(start:finish))
            endif
        endif

    end subroutine Merge_Sort_Dagger

    !> @brief Returns the distributed conjugate transpose of CSR matrix A.

    subroutine CSR_Dagger(A, partition_table, A_T, MPI_communicator, local_pointers)

        type(CSR), intent(in) :: A
        integer, dimension(:), intent(in) :: partition_table
        type(CSR), intent(out) :: A_T
        integer, intent(in) :: MPI_communicator

        logical, optional, intent(in) :: local_pointers
        integer :: lb, ub
        integer :: element_lb_T, element_ub_T

        integer :: nz

        integer, dimension(:), allocatable :: row_indexes, column_indexes
        integer, dimension(:), allocatable :: column_indexes_in
        complex(dp), dimension(:), allocatable :: values

        integer, dimension(:), allocatable :: send_counts, rec_counts
        integer, dimension(:), allocatable :: send_disps, rec_disps

        integer(int64) :: received_nnz

        integer, dimension(:), allocatable :: mapping_disps
        integer, dimension(:), allocatable :: column_indexes_out, row_indexes_out
        complex(dp), dimension(:), allocatable :: values_out

        integer, dimension(:), allocatable :: target_rank

        integer :: i, j

        !MPI_Environment
        integer :: rank
        integer :: flock
        integer :: ierr

        call MPI_comm_size(MPI_communicator, flock, ierr)
        call MPI_comm_rank(MPI_communicator, rank, ierr)

        lb = partition_table(rank + 1)
        ub = partition_table(rank +2) - 1

        nz = size(A%col_indexes)

        A_T%rows = A%rows
        A_T%columns = A%columns

        allocate(column_indexes(A%row_starts(lb):A%row_starts(ub + 1) - 1))
        allocate(row_indexes(A%row_starts(lb):A%row_starts(ub + 1) - 1))
        allocate(values(A%row_starts(lb):A%row_starts(ub + 1) - 1))

        do i = lb, ub
            do j = A%row_starts(i), A%row_starts(i + 1) - 1
                row_indexes(j) = i
            enddo
        enddo

        do i = A%row_starts(lb), A%row_starts(ub + 1) - 1
            column_indexes(i) = A%col_indexes(i)
        enddo

        do i = A%row_starts(lb), A%row_starts(ub + 1) - 1
            values(i) = A%values(i)
        enddo

        allocate(send_counts(flock))

        send_counts = 0

        allocate(target_rank(A%row_starts(lb):A%row_starts(ub + 1) - 1))

        do i = lbound(column_indexes, 1), ubound(column_indexes, 1)

            do j = flock, 1, -1
                if (column_indexes(i) >= partition_table(j)) then
                    send_counts(j) = send_counts(j) + 1
                    target_rank(i) = j
                    exit
                endif
            enddo
        enddo

        allocate(send_disps(flock))

        send_disps(1) = 0

        do i = 2, flock
            send_disps(i) = send_disps(i - 1) + send_counts(i - 1)
        enddo

        allocate(mapping_disps(flock))

        mapping_disps = 0

        allocate(column_indexes_out(A%row_starts(lb):A%row_starts(ub + 1) - 1))
        allocate(row_indexes_out(A%row_starts(lb):A%row_starts(ub + 1) - 1))
        allocate(values_out(A%row_starts(lb):A%row_starts(ub + 1) - 1))

        do i = lb, ub
            do j = A%row_starts(i), A%row_starts(i + 1) - 1

                column_indexes_out(A%row_starts(lb) + send_disps(target_rank(j))  &
                    + mapping_disps(target_rank(j))) = column_indexes(j)

                values_out(A%row_starts(lb) + send_disps(target_rank(j))  &
                    + mapping_disps(target_rank(j))) = conjg(values(j))

                row_indexes_out(A%row_starts(lb) + send_disps(target_rank(j))  &
                    + mapping_disps(target_rank(j))) = row_indexes(j)

                mapping_disps(target_rank(j)) = mapping_disps(target_rank(j)) + 1

            enddo
        enddo

        allocate(rec_counts(flock))

        call MPI_alltoall(  send_counts, &
                            1, &
                            MPI_integer, &
                            rec_counts, &
                            1, &
                            MPI_integer, &
                            MPI_communicator, &
                            ierr)

        received_nnz = sum(int(rec_counts, int64))
        element_lb_T = CSR_Element_Start(received_nnz, rank, MPI_communicator, local_pointers)
        element_ub_T = element_lb_T - 1 + int(received_nnz)

        allocate(column_indexes_in(element_lb_T:element_ub_T))
        allocate(A_T%col_indexes(element_lb_T:element_ub_T))
        allocate(A_T%values(element_lb_T:element_ub_T))

        allocate(rec_disps(flock))

        rec_disps(1) = 0

        do i = 2, flock
            rec_disps(i) = rec_disps(i - 1) + rec_counts(i - 1)
        enddo

        call MPI_alltoallv( column_indexes_out, &
                            send_counts, &
                            send_disps, &
                            MPI_integer, &
                            column_indexes_in, &
                            rec_counts, &
                            rec_disps, &
                            MPI_integer, &
                            MPI_communicator, &
                            ierr)

        call MPI_alltoallv( row_indexes_out, &
                            send_counts, &
                            send_disps, &
                            MPI_integer, &
                            A_T%col_indexes, &
                            rec_counts, &
                            rec_disps, &
                            MPI_integer, &
                            MPI_communicator, &
                            ierr)

        call MPI_alltoallv( values_out, &
                            send_counts, &
                            send_disps, &
                            MPI_double_complex, &
                            A_T%values, &
                            rec_counts, &
                            rec_disps, &
                            MPI_double_complex, &
                            MPI_communicator, &
                            ierr)

        call Merge_Sort_Dagger( column_indexes_in, &
                                A_T%col_indexes, &
                                A_T%values, &
                                1, &
                                size(column_indexes_in))

        allocate(A_T%row_starts(lb:ub+1))

        A_T%row_starts(lb) = element_lb_T
        A_T%row_starts(lb + 1:ub + 1) = 0

        do i = element_lb_T, element_ub_T
            A_T%row_starts(column_indexes_in(i) + 1) = &
                A_T%row_starts(column_indexes_in(i) + 1) + 1
        enddo

        do i = lb + 1, ub + 1
            A_T%row_starts(i) = A_T%row_starts(i) + A_T%row_starts(i - 1)
        enddo

        call MPI_barrier(MPI_communicator, ierr)

    end subroutine CSR_Dagger

    !> @brief Determine which right-hand side entries each MPI rank needs
    !> for sparse matrix multiplication.

    subroutine Reconcile_Communications(A, &
                                        partition_table, &
                                        MPI_communicator)

        type(CSR), intent(inout) :: A
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(in) :: MPI_communicator

        integer :: lb, ub, lb_elements, ub_elements

        integer :: total_rec_inds
        integer, dimension(:), allocatable :: RHS_rec_inds
        integer :: node

        integer :: i, j

        !MPI ENVIRONMENT
        integer :: rank
        integer :: flock
        integer :: ierr

        integer, dimension(:), allocatable :: requires
        integer :: ind

        call mpi_comm_size( mpi_communicator, &
                            flock, &
                            ierr)

        call mpi_comm_rank( mpi_communicator, &
                            rank, &
                            ierr)


        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        allocate(requires(partition_table(flock + 1) - 1))

        requires = 0

        lb_elements = lbound(A%col_indexes, 1)
        ub_elements = ubound(A%col_indexes, 1)

        if (.not. associated(A%num_rec_inds)) then
            allocate(A%num_rec_inds(flock))
        endif

        A%num_rec_inds = 0

        total_rec_inds = 0

        do i = lb_elements, ub_elements
            j = A%col_indexes(i)
            if (j < lb .or. j > ub) requires(j) = 1
        enddo

        do node = 1, flock
            A%num_rec_inds(node) = count(requires(partition_table(node):partition_table(node + 1) - 1) /= 0)
        enddo
        total_rec_inds = sum(A%num_rec_inds)
        if (int(ub, int64) + total_rec_inds > huge(ub)) then
            write(error_unit,*) 'Halo subscripts exceed the 32-bit range on rank', rank
            call MPI_Abort(MPI_communicator, 1, ierr)
        endif

        allocate(RHS_rec_inds(total_rec_inds))
        ind = 0
        do i = 1, size(requires)
            if (requires(i) /= 0) then
                ind = ind + 1
                requires(i) = ind
                RHS_rec_inds(ind) = i
            endif
        enddo

        ! Calculate the offset of the external elements in the 1D receive buffer.

        if (.not. associated(A%rec_disps)) then
            allocate(A%rec_disps(flock))
        endif

        A%rec_disps = 0
        do i = 2, flock
            A%rec_disps(i) = A%rec_disps(i - 1) + A%num_rec_inds(i - 1)
        enddo

        ! Map remote columns to the received entries following the local RHS.

        if (.not. associated(A%local_col_inds)) then
            allocate(A%local_col_inds(lb_elements:ub_elements))
        endif

        do i = lb_elements, ub_elements
            j = A%col_indexes(i)
            A%local_col_inds(i) = j
            if (j < lb .or. j > ub) A%local_col_inds(i) = ub + requires(j)
        enddo

        if (.not. associated(A%num_send_inds)) then
            allocate(A%num_send_inds(flock))
        endif

        call MPI_alltoall(  A%num_rec_inds, &
                            1, &
                            MPI_integer, &
                            A%num_send_inds, &
                            1, &
                            MPI_integer, &
                            MPI_communicator, &
                            ierr)

        if (sum(int(A%num_send_inds, int64)) >= huge(ind)) then
            write(error_unit,*) 'MPI send buffer exceeds the 32-bit range on rank', rank
            call MPI_Abort(MPI_communicator, 1, ierr)
        endif

        if (.not. associated(A%RHS_send_inds)) then
            allocate(A%RHS_send_inds(sum(A%num_send_inds)))
        endif

        if (.not. associated(A%send_disps)) then
            allocate(A%send_disps(flock))
        endif

        ! Calculate the offset of the local send elements in the 1D send buffer.

        A%send_disps = 0
        do i = 2, flock
            A%send_disps(i) = A%send_disps(i - 1) + A%num_send_inds(i - 1)
        enddo

        ! Obtain which RHS indexes to send to each node from the locally
        ! determined receive elements.

        call MPI_alltoallv( RHS_rec_inds, &
                            A%num_rec_inds, &
                            A%rec_disps, &
                            MPI_integer, &
                            A%RHS_send_inds, &
                            A%num_send_inds, &
                            A%send_disps, &
                            MPI_integer, &
                            MPI_communicator, &
                            ierr)

    end subroutine Reconcile_Communications

    !
    !   subroutine: SpMM
    !
    !> @brief MPI parallel CSR sparse-matrix dense-matrix multiplication.
    !>
    !> @details Computes A^(n)*B = C where *A* is a matrix, *B*, *C* are matrices
    !> and *n* is an integer.
    !>
    !> *CSR* matrix *A* and matrix *B* must be paritioned as given by
    !> *partition_table*, see @ref sparse_parts::generate_partition_table for more.
    !>


    subroutine SpMM(A, &
                    n, &
                    B_local, &
                    partition_table, &
                    rank, &
                    C_local, &
                    mpi_communicator)

        type(CSR), intent(inout) :: A !< @param Local CSR array partition.
        integer, intent(in) :: n !< @param Exponent on A.
        complex(dp), dimension(:,:), intent(in) :: B_local !< @param Local input dense matrix partition.
        integer, dimension(:), intent(in) :: partition_table !< @param MPI communicator partition scheme.
        integer, intent(in) :: rank
        complex(dp), &
        dimension(partition_table(rank + 1):partition_table(rank + 2) - 1, size(B_local, 2)), &
            intent(out) :: C_local  !< Local output matrix partition.
        integer, intent(in) :: mpi_communicator !< @param MPI communicator handel.

        complex(dp), dimension(:,:), allocatable :: B_resize
        complex(dp), dimension(:,:), allocatable :: rec_values, send_values

        integer :: num_rec, num_send
        integer :: lb, ub, lb_resize, ub_resize

        integer :: B_col, rhs_type

        integer :: i, j, k, l

        ! MPI environment
        integer :: ierr

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        num_rec = sum(A%num_rec_inds)
        num_send = sum(A%num_send_inds)

        lb_resize = ub + 1
        ub_resize = ub + num_rec

        B_col = size(B_local, 2)

        allocate(B_resize(lb:ub_resize, B_col))

        B_resize(lb:ub,:) = B_local

        allocate(rec_values(B_col, num_rec), send_values(B_col, num_send))
        ! One MPI element contains all right-hand sides for a row.
        call MPI_Type_contiguous(B_col, MPI_double_complex, rhs_type, ierr)
        call MPI_Type_commit(rhs_type, ierr)
        do l = 1, n
            !$omp parallel do
            do i = 1, num_send
                send_values(:, i) = B_resize(A%RHS_send_inds(i), :)
            enddo
            !$omp end parallel do

            call MPI_alltoallv(send_values, A%num_send_inds, A%send_disps, rhs_type, &
                               rec_values, A%num_rec_inds, A%rec_disps, rhs_type, MPI_communicator, ierr)
            do j = 1, B_col
                B_resize(lb_resize:ub_resize, j) = rec_values(j, :)
            enddo

            C_local = 0

            !$omp parallel do
            do i = lb, ub
                do j = A%row_starts(i), A%row_starts(i + 1) - 1
                    do k = 1, B_col
                        C_local(i, k) = A%values(j)*B_resize(A%local_col_inds(j), k) &
                            + C_local(i, k)
                    enddo
                enddo
            enddo
            !$omp end parallel do

            if (l < n) then
                B_resize(lb:ub,:) = C_local
            endif

        enddo

        call MPI_Type_free(rhs_type, ierr)

    end subroutine SpMM

    !> @brief Construct the MPI graph for neighbour exchange and locate the
    !> first remote column in each local row.
    subroutine Prepare_Neighbour_Exchange(A, lb, ub, comm)
        type(CSR), intent(inout) :: A
        integer, intent(in) :: lb, ub, comm
        integer, allocatable :: ranks(:), sources(:), destinations(:)
        integer :: i, j, ierr

        ranks = [(i - 1, i = 1, size(A%num_rec_inds))]
        sources = pack(ranks, A%num_rec_inds > 0)
        destinations = pack(ranks, A%num_send_inds > 0)
        call MPI_Dist_graph_create_adjacent(comm, size(sources), sources, MPI_UNWEIGHTED, &
                size(destinations), destinations, MPI_UNWEIGHTED, MPI_INFO_NULL, .false., A%halo_comm, ierr)
        A%halo_rec_counts = A%num_rec_inds(sources + 1)
        A%halo_rec_disps = A%rec_disps(sources + 1)
        A%halo_send_counts = A%num_send_inds(destinations + 1)
        A%halo_send_disps = A%send_disps(destinations + 1)
        allocate(A%first_remote(lb:ub))
        do i = lb, ub
            do j = A%row_starts(i), A%row_starts(i + 1) - 1
                if (A%local_col_inds(j) > ub) exit
            enddo
            A%first_remote(i) = j
        enddo
    end subroutine Prepare_Neighbour_Exchange

    !> @brief Multiply the leading local entries of each row during the neighbour
    !> exchange, then complete the product with the received vector entries.
    subroutine Neighbour_SpMV(A, lb, ub, u_resize, send_values, rec_values, v_local)
        type(CSR), intent(in) :: A
        integer, intent(in) :: lb, ub
        complex(dp), intent(inout) :: u_resize(lb:)
        complex(dp), contiguous, asynchronous, intent(in) :: send_values(:)
        complex(dp), contiguous, asynchronous, intent(out) :: rec_values(:)
        complex(dp), intent(out) :: v_local(lb:)
        integer :: request, ierr, i, j

        v_local = 0
        call MPI_Ineighbor_alltoallv(send_values, A%halo_send_counts, A%halo_send_disps, MPI_double_complex, &
                rec_values, A%halo_rec_counts, A%halo_rec_disps, MPI_double_complex, A%halo_comm, request, ierr)
        !$omp parallel do
        do i = lb, ub
            do j = A%row_starts(i), A%first_remote(i) - 1
                v_local(i) = A%values(j)*u_resize(A%local_col_inds(j)) + v_local(i)
            enddo
        enddo
        !$omp end parallel do
        call MPI_Wait(request, MPI_STATUS_IGNORE, ierr)
        u_resize(ub + 1:) = rec_values
        !$omp parallel do
        do i = lb, ub
            do j = A%first_remote(i), A%row_starts(i + 1) - 1
                v_local(i) = A%values(j)*u_resize(A%local_col_inds(j)) + v_local(i)
            enddo
        enddo
        !$omp end parallel do
    end subroutine Neighbour_SpMV

    !
    !   subroutine: SpMV_Series
    !
    !> @brief MPI parallel CSR sparse-matrix dense-vector series multiplication.
    !>
    !> @details Each iteration computes v = A*u. The caller may modify the
    !> vector between iterations. If the series ends before *max_it*, call
    !> with both *start_it* and *max_it* set to zero to release the workspace.
    !>
    !> *CSR* matrix *A* and vector *u* must be partitioned as given by
    !> *partition_table*, see @ref sparse_parts::generate_partition_table for more.
    !>

    subroutine SpMV_Series( A, &
                            u_local, &
                            partition_table, &
                            start_it, &
                            current_it, &
                            max_it, &
                            rank, &
                            v_local, &
                            mpi_communicator)

        type(CSR), intent(in) :: A !< @param Local CSR array partition.
        complex(dp), dimension(:), intent(inout) :: u_local !< @param Local input vector partition.
        integer, dimension(:), intent(in) :: partition_table !< @param MPI communicator partition scheme.
        integer, intent(in) :: start_it !< @param starting multiplication index.
        integer, intent(in) :: current_it !< @param Current multiplication index
        integer, intent(in) :: max_it !< @param Final multiplication index.
        integer, intent(in) :: rank
        integer, intent(in) :: mpi_communicator !< @param MPI communicator handel.
        complex(dp), dimension(partition_table(rank + 1):partition_table(rank + 2) - 1), &
            intent(inout) :: v_local !< @param Local output vector partition.

        complex(dp), dimension(:), allocatable, save :: u_resize
        complex(dp), dimension(:), allocatable, save :: send_values
        complex(dp), dimension(:), allocatable, save :: rec_values

        integer :: lb, ub, lb_resize, ub_resize
        integer :: num_send, num_rec

        integer :: i, j

        ! MPI environment
        integer :: ierr

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        num_rec = sum(A%num_rec_inds)
        num_send = sum(A%num_send_inds)

        lb_resize = ub + 1
        ub_resize = ub + num_rec

        if (start_it == 0 .and. max_it == 0) then
            if (allocated(u_resize)) deallocate(u_resize, rec_values, send_values)
            return
        endif
        if (allocated(u_resize) .and. current_it == start_it) then
            if (lbound(u_resize, 1) /= lb .or. ubound(u_resize, 1) /= ub_resize .or. &
                    size(rec_values) /= num_rec .or. size(send_values) /= num_send) then
                deallocate(u_resize, rec_values, send_values)
            endif
        endif
        if (.not. allocated(u_resize)) then
            allocate(u_resize(lb:ub_resize), rec_values(num_rec), send_values(num_send))
        endif

        u_resize(lb:ub) = u_local

        !$omp parallel do
        do i = 1, num_send
            send_values(i) = u_resize(A%RHS_send_inds(i))
        enddo
        !$omp end parallel do

        if (A%halo_comm /= MPI_COMM_NULL) then
            call Neighbour_SpMV(A, lb, ub, u_resize, send_values, rec_values, v_local)
        else
            call MPI_alltoallv(send_values, A%num_send_inds, A%send_disps, MPI_double_complex, &
                               rec_values, A%num_rec_inds, A%rec_disps, MPI_double_complex, MPI_communicator, ierr)
            u_resize(lb_resize:ub_resize) = rec_values
            v_local = 0
            !$omp parallel do
            do i = lb, ub
                do j = A%row_starts(i), A%row_starts(i + 1) - 1
                    v_local(i) = A%values(j)*u_resize(A%local_col_inds(j)) + v_local(i)
                enddo
            enddo
            !$omp end parallel do
        endif

        if (current_it == max_it) then
            deallocate(u_resize)
            deallocate(rec_values)
            deallocate(send_values)
        endif

    end subroutine SpMV_Series

    !> @brief Sort the rows of a distributed CSR matrix.

    subroutine Merge_CSR(   column_indexes, &
                            values, &
                            start, &
                            mid, &
                            finish)

        integer, intent(inout), dimension(:) :: column_indexes
        complex(dp), intent(inout), dimension(:) :: values
        integer, intent(in) :: start
        integer, intent(in) :: mid
        integer, intent(in) :: finish

        integer, dimension(:), allocatable :: col_ind_temp
        complex(dp), dimension(:), allocatable :: val_temp
        integer :: i, j, k

        allocate(col_ind_temp(finish - start + 1))
        allocate(val_temp(finish - start + 1))

        i = start
        j = mid + 1
        k = 1

        do while (i <= mid .and. j <= finish)

            if (column_indexes(i) <= column_indexes(j)) then
                col_ind_temp(k) = column_indexes(i)
                val_temp(k) = values(i)
                k = k + 1
                i = i + 1
            else
                col_ind_temp(k) = column_indexes(j)
                val_temp(k) = values(j)
                k = k + 1
                j = j+ 1
            endif

        enddo

        do while (i <= mid)
            col_ind_temp(k) = column_indexes(i)
            val_temp(k) = values(i)
            k = k + 1
            i = i + 1
        enddo

        do while (j <= finish)
            col_ind_temp(k) = column_indexes(j)
            val_temp(k) = values(j)
            k = k + 1
            j = j + 1
        enddo

        do i = start, finish
            column_indexes(i) = col_ind_temp(i - start + 1)
            values(i) = val_temp(i - start + 1)
        enddo

    end subroutine Merge_CSR

    !> @brief Sort the rows of a distributed CSR matrix.

    subroutine Insertion_Sort_CSR(   column_indexes, &
                                        values)

        integer, intent(inout), dimension(:) :: column_indexes
        complex(dp), intent(inout), dimension(:) :: values

        integer :: col_ind_temp
        complex(dp) :: val_temp

        integer :: i, j

        do i = 2, size(column_indexes)

            col_ind_temp = column_indexes(i)
            val_temp = values(i)

            j = i - 1

            do while (j >= 1)

                if (column_indexes(j) <= col_ind_temp) exit
                    column_indexes(j + 1) = column_indexes(j)
                    values(j + 1) = values(j)
                    j = j - 1
            enddo
            column_indexes(j + 1) = col_ind_temp
            values(j + 1) = val_temp

        enddo

    end subroutine Insertion_Sort_CSR

    !> @brief Sort the rows of a distributed CSR matrix.

    recursive subroutine Merge_Sort_CSR( column_indexes, &
                                            values, &
                                            start, &
                                            finish)

        integer, intent(inout), dimension(:) :: column_indexes
        complex(dp), intent(inout), dimension(:) :: values
        integer, intent(in) :: start
        integer, intent(in) :: finish

        integer :: mid

        if (start < finish) then
            if (finish - start >= 512) then

                mid = start + (finish - start) / 2

                call Merge_Sort_CSR( column_indexes, &
                                        values, &
                                        start, &
                                        mid)

                call Merge_Sort_CSR( column_indexes, &
                                        values, &
                                        mid + 1, &
                                        finish)

                call Merge_CSR(  column_indexes, &
                                    values, &
                                    start, &
                                    mid, &
                                    finish)

            else
                call insertion_sort_CSR( column_indexes(start:finish), &
                                            values(start:finish))
            endif
        endif

    end subroutine Merge_Sort_CSR

    subroutine Sort_CSR(A)

        type(CSR), intent(inout) :: A

        integer :: i

        do i = lbound(A%row_starts,1), ubound(A%row_starts,1) - 1

            call Merge_Sort_CSR(A%col_indexes(A%row_starts(i):A%row_starts(i + 1) - 1), &
                                A%values(A%row_starts(i):A%row_starts(i + 1) - 1), &
                                1, &
                                A%row_starts(i + 1) - A%row_starts(i))

        enddo

    end subroutine Sort_CSR

end module Sparse
