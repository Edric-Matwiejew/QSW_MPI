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
!   Module: Operators
!
!>  @brief Quantum Stochastic Walk operator formation, manipulation and
!>  operator specific sparse operations.
!

module Operators

    use :: ISO_Precisions
    use :: Sparse
    use :: MPI

    implicit none

    contains

    !
    !   Subroutine: Pad_Operator
    !
    !>  @brief Zero-pad a square CSR array.
    !
    !>  @details Returns a zero-padded version of a square CSR array with
    !>  'n_pad' zero-valued bottom rows and left columns.

    subroutine Pad_Operator(    H, &
                                n_pad, &
                                H_pad)

        type(CSR), intent(in) :: H !< @param Square CSR array.
        integer, intent(in) :: n_pad !< @param Padding rows and columns.
        type(CSR), intent(out) :: H_pad !< @param Padded CSR array.

        allocate(   H_pad%row_starts(H%rows + n_pad + 1), &
                    H_pad%col_indexes(size(H%col_indexes)), &
                    H_pad%values(size(H%col_indexes)))

        H_pad%rows = H%rows + n_pad
        H_pad%columns = H_pad%rows

        H_pad%row_starts(1:H%rows + 1) = H%row_starts(1:H%rows + 1)

        H_pad%row_starts(H%rows + 2:H_pad%rows + 1) = H%row_starts(H%rows + 1)

        H_pad%col_indexes(1:size(H%col_indexes)) = &
            H%col_indexes(1:size(H%col_indexes))

        H_pad%values(1:size(H%values)) = H%values(1:size(H%values))

    end subroutine Pad_Operator

    !
    !   Subroutine: Symmetrise_Graph_Weights
    !
    !>  @brief Symmetrise the values of a structurally symmetric CSR array.
    !
    !>  @details Given a CSR array, A, with structurally symmetric non-zero
    !>  values this returns an array 'values' such that when used in place of
    !>  A%values in the CSR datatype,
    !>
    !>> A(i,j) = max(A(i,j), A(j,i)).

    subroutine Symmetrise_Graph_Weights(A, values)

        type(CSR), intent(inout) :: A !< @param Structurally symmetric CSR array.
        complex(dp), dimension(:), intent(inout) :: values !< @param Symmetrised values array

        integer, dimension(:), allocatable :: offset
        integer :: i, j, k

        allocate(offset(A%rows))

        offset = 0
        values = A%values

        do i = 1, A%rows - 1
            do j = A%row_starts(i), A%row_starts(i+1) - 1
                if (A%col_indexes(j) > i) then
                    do k = A%row_starts(A%col_indexes(j)) + offset(A%col_indexes(j)), &
                        A%row_starts(A%col_indexes(j) + 1) - 1
                        if (A%col_indexes(k) > A%col_indexes(j)) then
                                exit
                        elseif (A%col_indexes(k) == i) then
                            offset(A%col_indexes(j)) = k - A%row_starts(A%col_indexes(j)) + 1
                            if (abs(A%values(k)) > abs(A%values(j))) then
                                values(j) = A%values(k)
                                exit
                            else
                                values(k) = A%values(j)
                                exit
                            endif
                        endif
                    enddo
                endif
            enddo
        enddo

    end subroutine Symmetrise_Graph_Weights

    !
    !   Subroutine: Generate_Scattering_Lindblad_Operators
    !
    !>  @brief Construct a QSW Lindblad operator matrix.
    !
    !>  @details Given a CSR transition matrix A this subroutine constructs
    !>  a CSR matrix representing the single-valued Lindblad operator matricies
    !>  as a linear sum.

    subroutine Generate_Scattering_Lindblad_Operators(A, L)

        type(CSR), intent(in) :: A !< @param CSR array.
        type(CSR), intent(out) :: L !< @param CSR array.

        integer, dimension(:), allocatable :: row_offsets

        integer :: i, j

        allocate(L%row_starts(size(A%row_starts)))
        allocate(L%col_indexes(size(A%col_indexes)))
        allocate(L%values(size(A%values)))

        L%rows = A%rows
        L%columns = A%columns
        L%row_starts = A%row_starts

        allocate(row_offsets(L%rows))

        row_offsets = 0

        do i = 1, A%rows
            do j= A%row_starts(i), A%row_starts(i + 1) - 1
            L%values(A%row_starts(A%col_indexes(j)) + row_offsets(A%col_indexes(j))) = &
                sqrt(abs(A%values(j)))
            L%col_indexes(A%row_starts(A%col_indexes(j)) + row_offsets(A%col_indexes(j))) = &
                i
            row_offsets(A%col_indexes(j)) = row_offsets(A%col_indexes(j)) + 1
            enddo
        enddo

        call Sort_CSR(L)

    end subroutine Generate_Scattering_Lindblad_Operators

    !
    !   Subroutine: Generate_Scattering_Superoperator
    !
    !>  @brief: Given a CSR Lindblad operator matrix, lower row partition bound
    !>  and upper row parition bound, a section of a super-operator describing
    !>  a continuous time random walk (CTRW) is constructed.

    subroutine Generate_Scattering_Superoperator(L, lower_bound, upper_bound, B)

        type(CSR), intent(in) :: L !< @param CSR array.
        integer, intent(in) :: upper_bound, lower_bound !< @param lower and upper bounds of the super-operator slice.
        type(CSR), intent(out) :: B !< @param CTRW super-operator row slice.

        integer :: lower_block, upper_block
        complex(dp), dimension(:), allocatable :: value_elements
        complex(dp), dimension(:), allocatable :: value_elements_shifted
        complex(dp), dimension(:), allocatable :: nz_diag, diag_vals
        integer, dimension(:), allocatable :: col_index
        integer, dimension(:), allocatable :: row_start

        integer :: lower_offset, upper_offset

        integer :: N

        integer :: i, j, k, row, indx

        logical :: diag_hit

        lower_block = ceiling(real(lower_bound)/L%rows)
        upper_block = ceiling(real(upper_bound)/L%rows)

        allocate(value_elements(L%rows))
        allocate(nz_diag(L%rows))
        allocate(diag_vals(L%rows**2))

        value_elements = 0
        nz_diag = 0
        diag_vals = 0
        N = L%rows

        do i = 1, L%rows
            do j = L%row_starts(i), L%row_starts(i + 1) - 1
                if (L%col_indexes(j) == i) then
                    row = (i - 1)*N + i
                    nz_diag(i) = L%values(j)
                    diag_vals(row) = abs(L%values(j))**2
                    value_elements(L%col_indexes(j)) = value_elements(L%col_indexes(j)) + &
                        abs(L%values(j))**2
                else
                    value_elements(L%col_indexes(j)) = value_elements(L%col_indexes(j)) + &
                        abs(L%values(j))**2
                endif
            enddo
        enddo

        do i = 1, N
             do j = 1, N
                row = (i - 1)*N + j
                diag_vals(row) = diag_vals(row) - 0.5*(value_elements(i) + value_elements(j))
            enddo
        enddo

        allocate(row_start(L%rows**2 + 1))
        allocate(col_index(L%rows**2))

        row_start(1) = 1
        row_start(2:size(row_start)) = 0

        col_index = 0

        do i = 1, L%rows
            do j = 1, L%rows
                row = (i - 1)*L%rows + j
                if (abs(value_elements(i) + value_elements(j)) > epsilon(0.d0)) then
                    row_start(row + 1) = 1
                    col_index(row) = row
                endif
            enddo
        enddo

        do i = 1, L%rows
            row = (i -1)*L%rows + i
            row_start(row + 1) = row_start(row + 1) + (L%row_starts(i + 1) - L%row_starts(i))
            if (abs(nz_diag(i)) > epsilon(0.d0)) then
                row_start(row + 1) = row_start(row + 1) - 1
            endif
        enddo

        call Prefix_Sum(row_start)

        allocate(B%values(row_start(lower_bound):row_start(upper_bound + 1) - 1))
        allocate(B%col_indexes(row_start(lower_bound):row_start(upper_bound + 1) - 1))

        allocate(value_elements_shifted(0:L%rows - 1))

        value_elements_shifted(0) = value_elements(L%rows)
        value_elements_shifted(1:L%rows - 1) = value_elements(1:L%rows - 1)

        lower_offset = (lower_bound - (lower_block - 1)*L%rows) - 1
        upper_offset = (upper_bound - (upper_block - 1)*L%rows) - L%rows

        diag_hit = .false.
        do i = lower_block, upper_block
            do j = 1 + Kronecker_Delta(i, lower_block)*lower_offset, L%rows &
                + Kronecker_Delta(i, upper_block)*upper_offset

                row = (i - 1)*L%rows + j

                if ((i == j) .and. (L%row_starts(i + 1) - L%row_starts(i) /= 0)) then

                    indx = 0

                    do k = L%row_starts(i), L%row_starts(i + 1) - 1

                        if ((L%col_indexes(k) - 1)*L%rows + L%col_indexes(k) == row) then

                            B%values(row_start(row) + indx) = diag_vals(row)
                            B%col_indexes(row_start(row) + indx) = &
                                (L%col_indexes(k) - 1)*L%rows + L%col_indexes(k)
                            diag_hit = .true.

                        else

                            B%values(row_start(row) + indx) = L%values(k)**2
                            B%col_indexes(row_start(row) + indx) = &
                                (L%col_indexes(k) - 1)*L%rows + L%col_indexes(k)

                        endif

                        indx = indx + 1
                    enddo

                    ! If the Lindblad matrix has zero diagonals.
                    if (.not. diag_hit) then
                       if(abs(diag_vals(row)) > epsilon(0.d0)) then
                        B%values(row_start(row) + indx) = diag_vals(row)
                        B%col_indexes(row_start(row) + indx) = row
                        diag_hit = .false.
                        endif
                     endif
                else

                    B%values(row_start(row)) = diag_vals(row)
                    B%col_indexes(row_start(row)) = row

                endif
            enddo

        enddo

        allocate(B%row_starts(lower_bound:upper_bound + 1))

        B%row_starts(lower_bound:upper_bound + 1) = row_start(lower_bound:upper_bound + 1)
        B%rows = upper_bound - lower_bound + 1
        B%columns = L%columns**2

        call Sort_CSR(B)

    end subroutine Generate_Scattering_Superoperator

    !
    !   Subroutine: Generate_Coherent_Superoperator
    !
    !>  @brief Given a graph Hamiltonian H, lower row bound and upper row bound,
    !>  a row-slice of a super-operator describing a continous time quantum walk
    !>  is constructed.

    subroutine Generate_Coherent_Superoperator( H, &
                                                lower_bound, &
                                                upper_bound, &
                                                B)

        type(CSR), intent(in) :: H !< @param Graph Hamiltonian.
        integer, intent(in) :: lower_bound, upper_bound !< @param Upper and lower row bounds
        type(CSR), intent(out) :: B !< @param CTWQ super-operator row-slice.

        integer :: lower_block, upper_block

        integer, dimension(:), allocatable :: left, right
        complex(dp), dimension(:), allocatable :: diag_val

        integer, dimension(:), allocatable :: row_start

        integer :: lower_offset, upper_offset

        integer :: i, j, row

        lower_block = ceiling(real(lower_bound)/H%rows)
        upper_block = ceiling(real(upper_bound)/H%rows)

        allocate(left(H%rows))
        allocate(right(H%rows))
        allocate(diag_val(H%rows))

        left = 0
        right = 0
        diag_val = 0

        do i = 1, H%rows
            do j = H%row_starts(i), H%row_starts(i + 1) - 1
                if (H%col_indexes(j) < i) then
                    left(i) = left(i) + 1
                elseif (H%col_indexes(j) == i) then
                    diag_val(i) = H%values(j)
                else
                    right(i) = right(i) + 1
                endif
            enddo
        enddo

        allocate(row_start(H%rows**2 + 1))

        row_start(1) = 1
        row_start(2:size(row_start)) = 0

        do i = 1, H%rows
            do j = 1, H%rows
                row = (i - 1)*H%rows + j
                row_start(row + 1) = left(i) + (H%row_starts(j + 1) - H%row_starts(j)) + right(i)
                if ((abs(diag_val(j)) < epsilon(0.d0)) .and. (abs(diag_val(i)) > epsilon(0.d0))) then
                    row_start(row + 1) = row_start(row + 1) + 1
                endif
                if ((abs(diag_val(i) - diag_val(j)) < epsilon(0.d0)) .and. (abs(diag_val(i)) > epsilon(0.d0))) then
                    row_start(row + 1) = row_start(row + 1) - 1
                endif
            enddo
        enddo

        call Prefix_Sum(row_start)

        allocate(B%values(row_start(lower_bound):row_start(upper_bound + 1) - 1))
        allocate(B%col_indexes(row_start(lower_bound):row_start(upper_bound + 1) - 1))

        lower_offset = (lower_bound - (lower_block - 1)*H%rows) - 1
        upper_offset = (upper_bound - (upper_block - 1)*H%rows) - H%rows

        do i = lower_block, upper_block
            do j = 1 + Kronecker_Delta(lower_block, i)*lower_offset, H%rows + Kronecker_Delta(upper_block, i)*upper_offset

                row = (i - 1)*H%rows + j

                B%values(row_start(row):row_start(row) + left(i) - 1) = &
                        & cmplx(0_dp, -1_dp, dp)*(-H%values(H%row_starts(i):H%row_starts(i) + left(i) - 1))
                B%col_indexes(row_start(row):row_start(row) + left(i) - 1) = &
                        & (H%col_indexes(H%row_starts(i):H%row_starts(i) + left(i) - 1) - 1)*H%columns + j

                B%values(row_start(row + 1) - right(i): row_start(row + 1) - 1) = &
                        & cmplx(0_dp, -1_dp, dp)*(-H%values(H%row_starts(i + 1) - right(i):H%row_starts(i + 1) - 1))
                B%col_indexes(row_start(row + 1) - right(i): row_start(row + 1) - 1) = &
                        & (H%col_indexes(H%row_starts(i + 1) - right(i):H%row_starts(i + 1) - 1) - 1)*H%columns + j

                B%values(row_start(row) + left(i): row_start(row) + left(i) + left(j) - 1) = &
                        & cmplx(0_dp, -1_dp, dp)*(H%values(H%row_starts(j):H%row_starts(j) + left(j) - 1))
                B%col_indexes(row_start(row) + left(i): row_start(row) + left(i) + left(j) - 1) = &
                        & H%col_indexes(H%row_starts(j):H%row_starts(j) + left(j) - 1) + (i - 1)*H%columns

                if (abs(diag_val(j) - diag_val(i)) > epsilon(0.d0)) then
                    B%values(row_start(row) + left(i) + left(j)) = cmplx(0_dp, -1_dp, dp)*(diag_val(j) - diag_val(i))
                    B%col_indexes(row_start(row) + left(i) + left(j)) = row
                endif

                B%values(row_start(row + 1) - right(i) - right(j):row_start(row + 1) - right(i) - 1) = &
                        & cmplx(0_dp, -1_dp, dp)*(H%values(H%row_starts(j + 1) - right(j): H%row_starts(j + 1) - 1))
                B%col_indexes(row_start(row + 1) - right(i) - right(j):row_start(row + 1) - right(i) - 1) = &
                        & H%col_indexes(H%row_starts(j + 1) - right(j): H%row_starts(j + 1) - 1) + (i - 1)*H%columns

            enddo
        enddo

        allocate(B%row_starts(lower_bound:upper_bound + 1))

        B%row_starts(lower_bound:upper_bound + 1) = row_start(lower_bound:upper_bound + 1)
        B%rows = upper_bound - lower_bound + 1
        B%columns = H%columns**2

    end subroutine Generate_Coherent_Superoperator

    !
    !   Subroutine: Generate_Graph_Hamiltonian
    !
    !>  @brief Generate a transition matrix from a graph adjacency matrix.
    !
    !>  @details Given the transition rate, gamma, and adjacency matrix, A,
    !>  this constructs a graph transition matrix. If A is non-symmetric the
    !>  transition matrix, B, will also be non-symmetric.

    subroutine Generate_Graph_Hamiltonian(gamma, A, B)

        real(dp), intent(in) :: gamma !< @param System-wide transition rate.
        type(CSR), intent(in) :: A !< @param Graph adjacency matrix.
        type(CSR), intent(out) :: B !< @param Transition matrix.

        complex(dp), dimension(A%rows) :: out_degrees

        integer :: nnz
        integer :: elements

        !SORT
        integer :: i, j, k
        integer :: temp_indx
        complex(dp) :: temp

        allocate(B%row_starts(A%rows + 1))

        B%rows = A%rows
        B%columns = A%columns

        out_degrees = 0

        do i = 1, A%rows
            do j = A%row_starts(i), A%row_starts(i + 1) - 1
                out_degrees(i) = out_degrees(i) + gamma*abs(A%values(j))
            enddo
        enddo

        B%row_starts = A%row_starts

        do i = 1, A%rows
            if (abs(out_degrees(i)) > epsilon(0.d0)) then
                B%row_starts(i + 1 : A%rows + 1) = B%row_starts(i + 1 : A%rows + 1) + 1
            endif
        enddo

        nnz = B%row_starts(B%rows + 1) - 1

        allocate(B%col_indexes(nnz))
        allocate(B%values(nnz))

        do i = 1, B%rows

            elements = A%row_starts(i + 1) - A%row_starts(i)

            if ((elements == 0) .and. (abs(out_degrees(i)) < epsilon(0.d0))) cycle

            B%col_indexes(B%row_starts(i):B%row_starts(i) + elements - 1) = &
                A%col_indexes(A%row_starts(i):A%row_starts(i + 1) - 1)

            B%values(B%row_starts(i):B%row_starts(i) + elements - 1) = &
                -gamma*A%values(A%row_starts(i):A%row_starts(i + 1) - 1)

            if (abs(out_degrees(i)) > epsilon(0.d0)) then
               B%col_indexes(B%row_starts(i + 1) - 1) = i
               B%values(B%row_starts(i + 1) - 1) = out_degrees(i)
            endif

        enddo

        do k = 1, B%rows
            do i = B%row_starts(k) + 1, B%row_starts(k + 1) - 1
                temp = B%values(i)
                temp_indx = B%col_indexes(i)
                j = i - 1
                do while (j >= B%row_starts(k) )
                    if (B%col_indexes(j) <= temp_indx) exit
                        B%values(j + 1) = B%values(j)
                        B%col_indexes(j + 1) = B%col_indexes(j)
                        j = j - 1
                enddo
                B%values(j + 1) = temp
                B%col_indexes(j + 1) = temp_indx
            enddo
        enddo

    end subroutine Generate_Graph_Hamiltonian

    !   Subroutine: Generate_Stochastic_Superoperator
    !
    !>  @bried: Construct a row-slice of a super-operator describing a
    !>  Quantum Stochastic Walk (QSW).
    !
    !>  @details Given a symmetric graph transition matrix, H, a QSW Lindblad
    !>  operator matrix, L, decoherence parameter, omega, lower row bound and
    !>  upper row bound, this constructs a row-slice of a super-operator describing
    !>  a QSW.

    subroutine Generate_Stochastic_Superoperator(H, &
                                    L, &
                                    omega, &
                                    lower_bound, &
                                    upper_bound, &
                                    B)

        type(CSR), intent(in) :: H, L !< @param Transition matrix and Lindblad operator matrix.
        real(dp), intent(in) :: omega !< @param Decoherence parameter.
        integer, intent(in) :: lower_bound, upper_bound !< @param Upper and lower bounds for the row-slice.
        type(CSR), intent(out) :: B !< @ QSW super-operator row-slice.

        integer :: lower_block, upper_block

        integer, dimension(:), allocatable :: H_left, H_right, L_left, L_right
        complex(dp), dimension(:), allocatable :: H_diag_vals, L_diag_vals
        complex(dp), dimension(:), allocatable :: out_degrees
        complex(dp), dimension(:), allocatable :: diag_vals

        integer, dimension(:), allocatable :: row_start

        integer :: lower_offset, upper_offset
        integer :: i, j, row, lower_indx, upper_indx

        integer :: N

        N = H%rows

        lower_block = ceiling(real(lower_bound)/H%rows)
        upper_block = ceiling(real(upper_bound)/H%rows)

        allocate(H_left(N))
        allocate(H_right(N))
        allocate(L_left(N))
        allocate(L_right(N))
        allocate(H_diag_vals(N))
        allocate(L_diag_vals(N))
        allocate(out_degrees(N))
        allocate(diag_vals(N**2))

        H_left = 0
        H_right = 0
        L_left = 0
        L_right = 0
        H_diag_vals = 0
        L_diag_vals = 0
        out_degrees = 0
        diag_vals = 0

        do i = 1, N
            do j = H%row_starts(i), H%row_starts(i + 1) - 1
                if (H%col_indexes(j) < i) then
                    H_left(i) = H_left(i) + 1
                elseif (H%col_indexes(j) == i) then
                    H_diag_vals(i) = -cmplx(0_dp, 1.0_dp - omega, dp)*H%values(j)
                else
                    H_right(i) = H_right(i) + 1
                endif
            enddo
        enddo

        do i = 1, L%rows
            do j = L%row_starts(i), L%row_starts(i + 1) - 1
                if (L%col_indexes(j) < i) then
                    L_left(i) = L_left(i) + 1
                elseif (L%col_indexes(j) > i) then
                    L_right(i) = L_right(i) + 1
                elseif (L%col_indexes(j) == i) then
                    L_diag_vals(i) = L%values(j)
                endif
            enddo
        enddo

        do i = 1, size(L%values, 1)
            out_degrees(l%col_indexes(i)) = out_degrees(l%col_indexes(i)) + &
                omega*l%values(i)*conjg(l%values(i))
        enddo

        do i = 1, N
           row = (i - 1)*N + i
            diag_vals(row) = omega*L_diag_vals(i)*conjg(L_diag_vals(i))
        enddo

        do i = 1, N
             do j = 1, N
                row = (i - 1)*N + j
                diag_vals(row) = diag_vals(row) - 0.5_dp*(out_degrees(i) + out_degrees(j))
            enddo
        enddo

        do i = 1, N
            diag_vals((i - 1)*N + 1:i*N) = &
                diag_vals((i - 1)*N + 1:i*N) + H_diag_vals - H_diag_vals(i)
        enddo

        allocate(row_start(N**2 + 1))

        row_start(1) = 1
        row_start(2:N**2 + 1) = 0

        do i = 1, N
            row = (i - 1)*N + i
            row_start(row + 1) = row_start(row + 1) + L%row_starts(i + 1) - L%row_starts(i)
            if (abs(L_diag_vals(i)) > epsilon(0.d0)) then
                row_start(row + 1) = row_start(row + 1) - 1
            endif
            do j = 1, N
                row = (i - 1)* N + j
                row_start(row + 1) = row_start(row + 1) + H_left(i) &
                    + (H%row_starts(j + 1) - H%row_starts(j)) + H_right(i)
                if (abs(H_diag_vals(j)) > epsilon(0.d0)) then
                    row_start(row + 1) = row_start(row + 1)-1
                endif
                if (abs(diag_vals(row)) > epsilon(0.d0)) then
                    row_start(row + 1) = row_start(row + 1) + 1
                endif
            enddo
        enddo

        call Prefix_Sum(row_start)

        allocate(B%values(row_start(lower_bound):row_start(upper_bound + 1) - 1))
        allocate(B%col_indexes(row_start(lower_bound):row_start(upper_bound + 1) - 1))

        lower_offset = (lower_bound -(lower_block - 1)*N) - 1
        upper_offset = (upper_bound - (upper_block - 1)*N) - N

        do i = lower_block, upper_block
            do j = 1 + Kronecker_Delta(lower_block, i)*lower_offset, N + &
                Kronecker_Delta(upper_block, i)*upper_offset

                row = (i - 1)*H%rows + j

                lower_indx = row_start(row)

                if (abs(diag_vals(row)) > epsilon(0.d0)) then
                    B%values(lower_indx) = diag_vals(row)
                    B%col_indexes(lower_indx) = row
                    lower_indx = lower_indx + 1
                endif

                upper_indx = lower_indx + H_left(i) - 1

                !leftmost H
                B%values(lower_indx:upper_indx) = &
                    & cmplx(0_dp, 1.0_dp - omega, dp)*H%values(H%row_starts(i):H%row_starts(i) + H_left(i) - 1)
                B%col_indexes(lower_indx:upper_indx) = &
                    (H%col_indexes(H%row_starts(i):H%row_starts(i) + H_left(i) - 1) - 1)*H%columns + j

                lower_indx = upper_indx + 1
                upper_indx = lower_indx + H_right(i) - 1

                !rightmost H
                B%values(lower_indx:upper_indx) = &
                        & cmplx(0_dp, 1.0_dp - omega, dp)*H%values(H%row_starts(i + 1) - H_right(i):H%row_starts(i + 1) - 1)
                B%col_indexes(lower_indx:upper_indx) = &
                    & (H%col_indexes(H%row_starts(i + 1) - H_right(i):H%row_starts(i + 1) - 1) - 1)*H%columns + j

                lower_indx = upper_indx + 1
                upper_indx = lower_indx + H_left(j) - 1

                !left H block
                B%values(lower_indx:upper_indx) = &
                    -cmplx(0_dp, 1_dp-omega, dp)*H%values(H%row_starts(j):H%row_starts(j) + H_left(j) - 1)
                B%col_indexes(lower_indx:upper_indx) = &
                    H%col_indexes(H%row_starts(j):H%row_starts(j) + H_left(j) - 1) + (i - 1)*H%columns

                lower_indx = upper_indx + 1
                upper_indx = lower_indx + H_right(j) - 1

                !right H block
                B%values(lower_indx:upper_indx) = &
                        & -cmplx(0_dp, 1.0_dp - omega, dp)*H%values(H%row_starts(j + 1) - H_right(j):H%row_starts(j + 1) - 1)
                B%col_indexes(lower_indx:upper_indx) = &
                        & H%col_indexes(H%row_starts(j + 1) - H_right(j):H%row_starts(j + 1) - 1) + (i - 1)*H%columns

                !leftmost L
                if ((i == j) .and. L_left(i) > 0) then


                    lower_indx = upper_indx + 1
                    upper_indx = lower_indx + L_left(i) - 1

                    B%values(lower_indx:upper_indx) = &
                            & omega*abs(L%values(L%row_starts(i):L%row_starts(i) + L_left(i) - 1))**2
                    B%col_indexes(lower_indx:upper_indx) = &
                            & L%col_indexes(L%row_starts(i):L%row_starts(i) + L_left(i) - 1)*(L%rows + 1) - L%rows

                endif

                !rightmost L
                if ((i == j) .and. L_right(i) > 0) then

                    lower_indx = upper_indx + 1
                    upper_indx = lower_indx + L_right(i) - 1

                    B%values(lower_indx:upper_indx) = &
                            & omega*abs(L%values(L%row_starts(i + 1) - L_right(i):L%row_starts(i + 1) - 1))**2
                    B%col_indexes(lower_indx:upper_indx) = &
                            & L%col_indexes(L%row_starts(i + 1) - L_right(i):L%row_starts(i + 1) - 1)*(L%rows + 1) - L%rows

                endif

            enddo

        enddo

        allocate(B%row_starts(lower_bound:upper_bound + 1))

        B%row_starts = row_start(lower_bound:upper_bound + 1)

        B%rows = upper_bound - lower_bound + 1
        B%columns = H%columns**2

        call Sort_CSR(B)

    end subroutine Generate_Stochastic_Superoperator

    !
    !   Subroutine: Add_Sources_and_Sinks
    !
    !>  @brief Add scattering channels to the Lindblad operator matrix, L, to
    !>  model absorption and emission processes.
    !
    !>  @details Given Lindblad operator matrix, L, sources sites/rates and
    !>  sink sites/rates, this outputs an augmented Lindblad operator of dimension
    !>   N + m, where N is the dimension of L and m is the number of sources and
    !>  sinks, describing absorption and emission from the underlying graph
    !>  structure.

    subroutine Add_Sources_and_Sinks(   L, &
                                        L_aug, &
                                        source_sites, &
                                        source_rates, &
                                        sink_sites, &
                                        sink_rates)

        type(CSR), intent(in) :: L !< @param Lindblad operator matrix.
        type(CSR), intent(out) :: L_aug !< @param Augmented Lindblad operator matrix.
        integer, dimension(:), optional :: source_sites !< @param Graph nodes to which a source is attached.
        real(dp), dimension(:), optional :: source_rates !< @param Transition rates of the attached sources.
        integer, dimension(:), optional :: sink_sites !< @param Graph nodes to which a sink is attached.
        real(dp), dimension(:), optional :: sink_rates !< @param Transition rates of the attached sinks.

        integer :: n_sources
        integer :: n_sinks
        integer :: elements

        integer :: i, indx, j, k
        integer :: temp_indx
        complex(dp) :: temp

        if (present(source_sites)) then
            n_sources = size(source_sites)
        else
            n_sources = 0
        endif

        if (present(sink_sites)) then
            n_sinks = size(sink_sites)
        else
            n_sinks = 0
        endif

        allocate(   L_aug%row_starts(L%rows + n_sources + n_sinks + 1), &
                    L_aug%col_indexes(size(L%col_indexes) + n_sources + n_sinks), &
                    L_aug%values(size(L%col_indexes) + n_sources + n_sinks))

        L_aug%rows = L%rows + n_sinks + n_sources
        L_aug%columns = L_aug%rows

        L_aug%row_starts(1:L%rows + 1) = L%row_starts(1:L%rows + 1)
        L_aug%row_starts(L%rows + 2:L_aug%rows + 1) = L%row_starts(L%rows + 1)

        L_aug%col_indexes = 0

        indx = 1
        do i = 1, n_sources

            L_aug%row_starts(source_sites(i) + 1:L_aug%rows + 1) = &
                L_aug%row_starts(source_sites(i) + 1:L_aug%rows + 1) + 1

            L_aug%col_indexes(L_aug%row_starts(source_sites(i) + 1) - 1) = &
            L%rows + indx

            L_aug%values(L_aug%row_starts(source_sites(i) + 1) - 1) = &
                sqrt(abs(source_rates(i)))
            indx = indx + 1

        enddo

        indx = 1
        do i = 1, n_sinks

            L_aug%row_starts(L%rows + n_sources + indx + 1 : L_aug%rows + 1) = &
                L_aug%row_starts(L%rows + n_sources + indx + 1 : L_aug%rows  + 1) + 1

            L_aug%col_indexes(L_aug%row_starts(L%rows + n_sources + indx + 1) - 1) = &
                sink_sites(indx)

            L_aug%values(L_aug%row_starts(L%rows + n_sources + indx + 1) - 1) = &
                sqrt(abs(sink_rates(indx)))

            indx = indx + 1

        enddo

        do i = 1, L%rows

            elements = L%row_starts(i + 1) - L%row_starts(i)

            L_aug%col_indexes(L_aug%row_starts(i):L_aug%row_starts(i) + elements - 1) = &
                L%col_indexes(L%row_starts(i):L%row_starts(i + 1) - 1)

            L_aug%values(L_aug%row_starts(i):L_aug%row_starts(i) + elements - 1) =&
                L%values(L%row_starts(i):L%row_starts(i + 1) - 1)

        enddo

        do k = 1, L_aug%rows
            do i = L_aug%row_starts(k) + 1, L_aug%row_starts(k + 1) - 1
                temp = L_aug%values(i)
                temp_indx = L_aug%col_indexes(i)
                j = i - 1
                do while (j >= L_aug%row_starts(k) )
                    if (L_aug%col_indexes(j) <= temp_indx) exit
                        L_aug%values(j + 1) = L_aug%values(j)
                        L_aug%col_indexes(j + 1) = L_aug%col_indexes(j)
                        j = j - 1
                enddo
                L_aug%values(j + 1) = temp
                L_aug%col_indexes(j + 1) = temp_indx
            enddo
        enddo

    end subroutine Add_Sources_and_Sinks

    !
    !   Subroutine: Transport_Lindblads
    !
    !>  @brief Creates a Lindblad operator matrix containing only absorption and
    !>  emission channels.
    !
    !>  @detail For the case of a CTQW incorperating absorption and emission
    !>  this creates a Lindblad operator matrix containing only these scattering
    !>  channels.
    !
    ! NOTE: The sink and site variables are mixed up in this subroutine.
    ! as a quick fix I have swappwed them as input arguments.

    subroutine Transport_Lindblads( N, &
                                    sink_sites, &
                                    sink_rates, &
                                    source_sites, &
                                    source_rates, &
                                    L_inout)

        integer, intent(in) :: N !< @param Dimension of the graph adjacency matrix.
        integer, dimension(:), intent(in) :: source_sites !< @param Graph nodes to which a source is attached.
        real(dp), dimension(:), intent(in) :: source_rates !< @param Transition rates of the attached sources.
        integer, dimension(:), intent(in) :: sink_sites !< @param Graph nodes to which a sink is attached.
        real(dp), dimension(:), intent(in) :: sink_rates !< @param Transition rates of the attached sinks.
        type(CSR), intent(out) :: L_inout !< @param Lindblad operator matrix describing only absorption and esmission processes.

        integer :: augs

        integer :: i

        augs = size(source_sites) + size(sink_sites)

        allocate(L_inout%row_starts(N + augs + 1))
        allocate(L_inout%col_indexes(size(source_sites) + size(sink_sites)))
        allocate(L_inout%values(size(source_sites) + size(sink_sites)))

        L_inout%rows = N + augs
        L_inout%columns = N + augs

        ! ARRAY CONSTRUCTOR
        L_inout%col_indexes(1:size(sink_sites)) = [(i + N, i = 1, size(sink_sites))]
        L_inout%values(1:size(sink_sites)) = sqrt(abs(sink_rates))
        L_inout%col_indexes(size(sink_sites) + 1: augs) = source_sites
        L_inout%values(size(sink_sites) + 1: augs) = sqrt(abs(source_rates))

        L_inout%row_starts(2:L_inout%rows + 1) = 0
        L_inout%row_starts(1) = 1

        do i = 1, size(sink_sites)
            L_inout%row_starts(sink_sites(i)+1) =  1
        enddo

        do i = 1, size(source_sites)
            L_inout%row_starts(i + N + size(sink_sites) + 1) = 1
        enddo

        call Prefix_Sum(L_inout%row_starts)

    end subroutine Transport_Lindblads

    !
    !   Subroutine: Generate_Coherent_Transport_Superoperator
    !
    !>  @brief Generate a super-operator slice describing a continous-time
    !>  quantum walk including absorption an emssion processes.
    !
    !>  @detail Given a symmetric graph transition matrix H, a Lindblad operator
    !>  matrix containing only absorption and emission channels, lower row bound,
    !>  and upper row bound, a row slice of a super-operator describing a CTQW
    !>  including absorption and emission processes is constructed.

    subroutine Generate_Coherent_Transport_Superoperator(H, &
                                    L, &
                                    lower_bound, &
                                    upper_bound, &
                                    B)

        type(CSR), intent(in) :: H, L !< @param Symmetric graph transition matrix and Lindblad operator matrix.
        integer, intent(in) :: lower_bound, upper_bound !< @param lower and upper row bounds of the super-operator slice.
        type(CSR), intent(out) :: B !< @param Super-operator row slice.

        integer :: lower_block, upper_block

        integer, dimension(:), allocatable :: H_left, H_right, L_left, L_right
        complex(dp), dimension(:), allocatable :: H_diag_vals, L_diag_vals
        complex(dp), dimension(:), allocatable :: out_degrees
        complex(dp), dimension(:), allocatable :: diag_vals

        integer, dimension(:), allocatable :: row_start

        integer :: lower_offset, upper_offset
        integer :: i, j, row, lower_indx, upper_indx

        integer :: N

        N = H%rows

        lower_block = ceiling(real(lower_bound)/H%rows)
        upper_block = ceiling(real(upper_bound)/H%rows)

        allocate(H_left(N))
        allocate(H_right(N))
        allocate(L_left(N))
        allocate(L_right(N))
        allocate(H_diag_vals(N))
        allocate(L_diag_vals(N))
        allocate(out_degrees(N))
        allocate(diag_vals(N**2))

        H_left = 0
        H_right = 0
        L_left = 0
        L_right = 0
        H_diag_vals = 0
        L_diag_vals = 0
        out_degrees = 0
        diag_vals = 0

        do i = 1, N
            do j = H%row_starts(i), H%row_starts(i + 1) - 1
                if (H%col_indexes(j) < i) then
                    H_left(i) = H_left(i) + 1
                elseif (H%col_indexes(j) == i) then
                    H_diag_vals(i) = -cmplx(0_dp, 1.0_dp, dp)*H%values(j)
                else
                    H_right(i) = H_right(i) + 1
                endif
            enddo
        enddo

        do i = 1, L%rows
            do j = L%row_starts(i), L%row_starts(i + 1) - 1
                if (L%col_indexes(j) < i) then
                    L_left(i) = L_left(i) + 1
                elseif (L%col_indexes(j) > i) then
                    L_right(i) = L_right(i) + 1
                elseif (L%col_indexes(j) == i) then
                    L_diag_vals(i) = L%values(j)
                endif
            enddo
        enddo

        do i = 1, size(L%values, 1)
            out_degrees(L%col_indexes(i)) = out_degrees(L%col_indexes(i)) + &
                L%values(i)*conjg(L%values(i))
        enddo

        do i = 1, N
           row = (i - 1)*N + i
            diag_vals(row) = L_diag_vals(i)*conjg(L_diag_vals(i))
        enddo

        do i = 1, N
             do j = 1, N
                row = (i - 1)*N + j
                diag_vals(row) = diag_vals(row) - 0.5_dp*(out_degrees(i) + out_degrees(j))
            enddo
        enddo

        do i = 1, N
            diag_vals((i - 1)*N + 1:i*N) = &
                diag_vals((i - 1)*N + 1:i*N) + H_diag_vals - H_diag_vals(i)
        enddo

        allocate(row_start(N**2 + 1))

        row_start(1) = 1
        row_start(2:N**2 + 1) = 0

        do i = 1, N
            row = (i - 1)*N + i
            row_start(row + 1) = row_start(row + 1) + L%row_starts(i + 1) - L%row_starts(i)
            if (abs(L_diag_vals(i)) > epsilon(0.d0)) then
                row_start(row + 1) = row_start(row + 1) - 1
            endif
            do j = 1, N
                row = (i - 1)* N + j
                row_start(row + 1) = row_start(row + 1) + H_left(i) &
                    + (H%row_starts(j + 1) - H%row_starts(j)) + H_right(i)
                if (abs(H_diag_vals(j)) > epsilon(0.d0)) then
                    row_start(row + 1) = row_start(row + 1)-1
                endif
                if (abs(diag_vals(row)) > epsilon(0.d0)) then
                    row_start(row + 1) = row_start(row + 1) + 1
                endif
            enddo
        enddo

        call Prefix_Sum(row_start)

        allocate(B%values(row_start(lower_bound):row_start(upper_bound + 1) - 1))
        allocate(B%col_indexes(row_start(lower_bound):row_start(upper_bound + 1) - 1))

        lower_offset = (lower_bound -(lower_block - 1)*N) - 1
        upper_offset = (upper_bound - (upper_block - 1)*N) - N

        do i = lower_block, upper_block
            do j = 1 + Kronecker_Delta(lower_block, i)*lower_offset, N + &
                Kronecker_Delta(upper_block, i)*upper_offset

                row = (i - 1)*H%rows + j

                lower_indx = row_start(row)

                if (abs(diag_vals(row)) > epsilon(0.d0)) then
                    B%values(lower_indx) = diag_vals(row)
                    B%col_indexes(lower_indx) = row
                    lower_indx = lower_indx + 1
                endif

                upper_indx = lower_indx + H_left(i) - 1

                !leftmost H
                B%values(lower_indx:upper_indx) = &
                    & cmplx(0_dp, 1.0_dp, dp)*H%values(H%row_starts(i):H%row_starts(i) + H_left(i) - 1)
                B%col_indexes(lower_indx:upper_indx) = &
                    (H%col_indexes(H%row_starts(i):H%row_starts(i) + H_left(i) - 1) - 1)*H%columns + j

                lower_indx = upper_indx + 1
                upper_indx = lower_indx + H_right(i) - 1

                !rightmost H
                B%values(lower_indx:upper_indx) = &
                        & cmplx(0_dp, 1.0_dp, dp)*H%values(H%row_starts(i + 1) - H_right(i):H%row_starts(i + 1) - 1)
                B%col_indexes(lower_indx:upper_indx) = &
                    & (H%col_indexes(H%row_starts(i + 1) - H_right(i):H%row_starts(i + 1) - 1) - 1)*H%columns + j

                lower_indx = upper_indx + 1
                upper_indx = lower_indx + H_left(j) - 1

                !left H block
                B%values(lower_indx:upper_indx) = -cmplx(0_dp, 1_dp, dp)*H%values(H%row_starts(j):H%row_starts(j) &
                        &  + H_left(j) - 1)
                B%col_indexes(lower_indx:upper_indx) = &
                    H%col_indexes(H%row_starts(j):H%row_starts(j) + H_left(j) - 1) + (i - 1)*H%columns

                lower_indx = upper_indx + 1
                upper_indx = lower_indx + H_right(j) - 1

                !right H block
                B%values(lower_indx:upper_indx) = &
                        & -cmplx(0_dp, 1.0_dp, dp)*H%values(H%row_starts(j + 1) - H_right(j):H%row_starts(j + 1) - 1)
                B%col_indexes(lower_indx:upper_indx) = &
                        & H%col_indexes(H%row_starts(j + 1) - H_right(j):H%row_starts(j + 1) - 1) + (i - 1)*H%columns

                !leftmost L
                if ((i == j) .and. L_left(i) > 0) then


                    lower_indx = upper_indx + 1
                    upper_indx = lower_indx + L_left(i) - 1

                    B%values(lower_indx:upper_indx) = &
                            & abs(L%values(L%row_starts(i):L%row_starts(i) + L_left(i) - 1))**2
                    B%col_indexes(lower_indx:upper_indx) = &
                            & L%col_indexes(L%row_starts(i):L%row_starts(i) + L_left(i) - 1)*(L%rows + 1) - L%rows

                endif

                !rightmost L
                if ((i == j) .and. L_right(i) > 0) then

                    lower_indx = upper_indx + 1
                    upper_indx = lower_indx + L_right(i) - 1

                    B%values(lower_indx:upper_indx) = &
                            & abs(L%values(L%row_starts(i + 1) - L_right(i):L%row_starts(i + 1) - 1))**2
                    B%col_indexes(lower_indx:upper_indx) = &
                            & L%col_indexes(L%row_starts(i + 1) - L_right(i):L%row_starts(i + 1) - 1)*(L%rows + 1) - L%rows

                endif

            enddo

        enddo

        allocate(B%row_starts(lower_bound:upper_bound + 1))

        B%row_starts(lower_bound:upper_bound + 1) = row_start(lower_bound:upper_bound + 1)

        B%rows = upper_bound - lower_bound + 1
        B%columns = H%columns**2

        call Sort_CSR(B)

    end subroutine Generate_Coherent_Transport_Superoperator

    !
    !   Subroutine: Prepare_Super_Operator
    !
    !>  @brief Construct and arbitrary super-operator row slice.
    !
    !>  @details Given decoherence parameter omega, transition matrix H,
    !>  Lindblad operator matrix L, source sites/rates, sink sites/rates and an
    !>  array describing the desired partitioning of the super-operator over a
    !>  MPI communicator, this determines which super-operator construction
    !>  subroutine is appropriate and provides a super-operator row-slice using
    !>  that subroutine.

    subroutine Prepare_Super_Operator(  omega, &
                                        H, &
                                        L, &
                                        source_sites, &
                                        source_rates, &
                                        sink_sites, &
                                        sink_rates, &
                                        partition_table, &
                                        M_local, &
                                        MPI_communicator)

        real(dp), intent(in) :: omega !< @param Decoherence parameter.
        type(CSR), intent(in) :: H !< @param Transition matrix.
        type(CSR), intent(in) :: L !< @param Lindblad operator matrix.
        integer, dimension(:), intent(in) :: source_sites !< @param Graph vertices to which a source is attached.
        real(dp), dimension(:), intent(in) :: source_rates !< @oaram Transition rates of the attached sources.
        integer, dimension(:), intent(in) :: sink_sites !< @param Graph vertices to which a sink is attached.
        real(dp), dimension(:), intent(in) :: sink_rates !< @param Transition rates of the attached sinks.
        integer, dimension(:), allocatable, intent(out) :: partition_table !< @param row-slice paritioning of the super-operator over an MPI communicator.
        type(CSR), intent(out) :: M_local !< @param Super-operator row-slice.
        integer, intent(in) :: MPI_communicator !< @param MPI communicator over which the super-operator will be constructed.

        integer :: L_bound
        integer :: U_bound

        type(CSR) :: H_aug
        type(CSR) :: L_aug

        integer :: N

        real(dp) :: delta = epsilon(omega)

        !MPI ENVIRONMENT
        integer :: rank
        integer :: ierr

        call MPI_comm_rank(MPI_communicator, rank, ierr)

        if ((size(source_sites) > 0) .or. (size(sink_sites) > 0)) then
            N = H%rows + size(source_sites) + size(sink_sites)

            call Generate_Partition_Table(  N**2, &
                                            partition_table, &
                                            MPI_communicator)

            L_bound = partition_table(rank + 1)
            U_bound = partition_table(rank + 2) - 1

            if (omega .le.  (0.0_dp + delta)) then


            call Transport_Lindblads(   H%rows, &
                                        source_sites, &
                                        source_rates, &
                                        sink_sites, &
                                        sink_rates, &
                                        L_aug)

            call Pad_Operator(  H, &
                                size(source_sites) + size(sink_sites), &
                                H_aug)

            call Generate_Coherent_Transport_Superoperator( H_aug, &
                                                    L_aug, &
                                                    L_bound, &
                                                    U_bound, &
                                                    M_local)

            deallocate(H_aug%row_starts, H_aug%col_indexes, H_aug%values)
            deallocate(L_aug%row_starts, L_aug%col_indexes, L_aug%values)

            elseif (omega .ge. (1.0_dp - delta)) then

                call Add_Sources_and_Sinks( L, &
                                            L_aug, &
                                            source_sites = source_sites, &
                                            source_rates = source_rates, &
                                            sink_sites = sink_sites, &
                                            sink_rates = sink_rates)

                call Generate_Scattering_Superoperator( L_aug, &
                                                        L_bound, &
                                                        U_bound, &
                                                        M_local)

            else

                call Add_Sources_and_Sinks( L, &
                                            L_aug, &
                                            source_sites = source_sites, &
                                            source_rates = (1.0_dp/omega)*source_rates, &
                                            sink_sites = sink_sites, &
                                            sink_rates = (1.0_dp/omega)*sink_rates)

                call Pad_Operator(  H, &
                                    size(source_sites) + size(sink_sites), &
                                    H_aug)

                call Generate_Stochastic_Superoperator( H_aug, &
                                                        L_aug, &
                                                        omega, &
                                                        L_bound, &
                                                        U_bound, &
                                                        M_local)

            deallocate(H_aug%row_starts, H_aug%col_indexes, H_aug%values)
            deallocate(L_aug%row_starts, L_aug%col_indexes, L_aug%values)

            endif

        else

            call Generate_Partition_Table(  H%rows**2, &
                                            partition_table, &
                                            MPI_communicator)

            L_bound = partition_table(rank + 1)
            U_bound = partition_table(rank + 2) - 1

            if (omega .le.  (0.0_dp + delta)) then

                call Generate_Coherent_Superoperator(   H, &
                                                        l_bound, &
                                                        u_bound, &
                                                        M_local)

            elseif (omega .ge. (1.0_dp - delta)) then

                call Generate_Scattering_Superoperator( L, &
                                                        L_bound, &
                                                        U_bound, &
                                                        M_local)

            else

                call Generate_Stochastic_Superoperator( H, &
                                                        L, &
                                                        omega, &
                                                        L_bound, &
                                                        U_bound, &
                                                        M_local)
            endif

        endif

    end subroutine Prepare_Super_Operator

!
!    Subroutine: COO_Vectorized_Commutator
!
!>   @breif Returns HSO = (I \otimes H - H^T \otimes).
!
!>   @details Given CSR array H, COO_Vectorized_Commutator returns (I \otimes H - H^T \otimes)
!>   as a COO matrix. This corresponds to the vectorization mapping,
!>
!>>  [H, rho(t)] <-> (I \otimes H - H^T \otimes H)vec(rho(t))


    subroutine COO_Vectorized_Commutator(   H, &
                                            partition_table, &
                                            HSO, &
                                            MPI_communicator)

        type(CSR), intent(in) :: H !< @param CSR array.
        integer, dimension(:), intent(in) :: partition_table !< @param Distributed row partitions.
        type(COO), intent(inout) :: HSO !< @param COO array.
        integer, intent(in) :: MPI_communicator !< @param MPI communicator over which to generate HSO

        type(CSR) :: H_T
        type(COO) :: HSO1, HSO2
        integer :: HSO1_nnz, HSO2_nnz

        integer :: lower_bound, upper_bound, lower_block, upper_block
        integer :: lower_offset, upper_offset
        integer :: local_rows
        integer :: row_index, col_index
        integer :: i, j, k

        ! MPI
        integer :: rank
        integer :: ierr

        call MPI_comm_rank(MPI_communicator, rank, ierr)

        lower_bound = partition_table(rank + 1)
        upper_bound = partition_table(rank + 2) - 1

        lower_block = ceiling(real(lower_bound)/H%rows)
        upper_block = ceiling(real(upper_bound)/H%rows)

        lower_offset = (lower_bound - (lower_block - 1)*H%rows) - 1
        upper_offset = (upper_bound - (upper_block - 1)*H%rows) - H%rows

        local_rows = upper_bound - lower_bound + 1

        HSO1%rows = H%rows**2
        HSO1%columns = H%columns**2

        HSO1_nnz = size(H%col_indexes)*(upper_block - lower_block + 2)
        allocate(HSO1%row_indexes(HSO1_nnz))
        allocate(HSO1%col_indexes(HSO1_nnz))
        allocate(HSO1%values(HSO1_nnz))

        do i = lower_block, upper_block
        do j = 1 + Kronecker_Delta(i, lower_block)*lower_offset, H%rows &
            + Kronecker_Delta(i, upper_block)*upper_offset
                do k = H%row_starts(j), H%row_starts(j + 1) - 1
                    row_index = (i - 1) * H%rows + j
                    col_index = (i - 1) * H%rows + H%col_indexes(k)
                    HSO1%nnz = HSO1%nnz + 1
                    HSO1%row_indexes(HSO1%nnz) = row_index
                    HSO1%col_indexes(HSO1%nnz) = col_index
                    HSO1%values(HSO1%nnz) = H%values(k)
                enddo
            enddo
        enddo

        HSO2%rows = H%rows**2
        HSO2%columns = H%columns**2

        call CSR_Transpose(H, H_T)

        HSO2_nnz = H%rows*maxval(H%row_starts(2:size(H%row_starts)) - &
            H%row_starts(1:size(H%row_starts)-1))*(upper_block - lower_block + 1)

        allocate(HSO2%row_indexes(HSO2_nnz))
        allocate(HSO2%col_indexes(HSO2_nnz))
        allocate(HSO2%values(HSO2_nnz))

        do i = lower_block, upper_block
            do j = 1 + Kronecker_Delta(i, lower_block)*lower_offset, H%rows &
                + Kronecker_Delta(i, upper_block)*upper_offset
                do k = H_T%row_starts(i), H_T%row_starts(i + 1) - 1
                    row_index = (i - 1)*H_T%rows + j
                    col_index = (H_T%col_indexes(k) - 1)*H_T%columns + j
                    HSO2%nnz = HSO2%nnz + 1
                    HSO2%row_indexes(HSO2%nnz) = row_index
                    HSO2%col_indexes(HSO2%nnz) = col_index
                    HSO2%values(HSO2%nnz) = -H_T%values(k)
                enddo
            enddo
        enddo

        call COO_Allocate_Sum(HSO1, HSO2, local_rows, HSO)
        call COO_Sum(HSO1, HSO2, HSO, partition_table, MPI_communicator)

    end subroutine COO_Vectorized_Commutator

!
!>    Subroutine: COO_Vectorized_Dissipator
!
!>    @brief Returns \tau \sum_{k=1}^{K}(L_k^* \otimes L_k - \frac{1}{2}(I \otimes L_k^dagger L_k + L_k^T L_k^* \otimes I)).

!>    @details Given a list of CSR Lindblad operators, this returns the term responsible for environental interactions in
!>    the vectorized Lindblad master equation.

    subroutine COO_Vectorized_Dissipator(   taus, &
                                            Ls, &
                                            partition_table, &
                                            LSO, &
                                            MPI_communicator)

        real(dp), dimension(:), intent(in) :: taus !> @param Lindblad operator coefficients.
        type(CSR), dimension(:), target, intent(in) :: Ls !> @param Array of CSR Lindblad operators.
        integer, dimension(:), intent(in) :: partition_table !> @param !> Distributed row partitions.
        type(COO), intent(inout) :: LSO !> @param COO matrix
        integer, intent(in) :: MPI_communicator !> @param MPI communicator of which to generate LSO.

        integer :: lower_bound, upper_bound
        integer :: lower_block, upper_block
        integer :: lower_offset, upper_offset
        integer :: local_rows

        integer, dimension(:), allocatable :: row_nnzs
        integer :: indx

        integer :: i, j, k, ii, ll

        type(CSR), pointer :: L

        !LSO1
        type(COO) :: LSO1
        integer :: L_row, L_nz
        integer :: current_block
        integer :: col_ind_disp
        integer :: row_nz

        !LSO2
        type(COO) :: LSO2_elements
        type(COO) :: LSO2
        type(CSR) :: L_T
        integer :: lb_indx, ub_indx
        integer :: offset
        integer :: col_index, row_index
        complex(dp) :: value

        !LSO3
        type(COO) :: LSO3_elements
        type(COO) :: LSO3
        integer :: row_indx
        integer :: cnt

        type(COO) :: LSO_temp_1, LSO_temp_2, LSO_temp_3

        ! MPI
        integer :: rank
        integer :: ierr

        call MPI_comm_rank(MPI_communicator, rank, ierr)

        lower_bound = partition_table(rank + 1)
        upper_bound = partition_table(rank + 2) - 1

        lower_block = ceiling(real(lower_bound)/Ls(1)%rows)
        upper_block = ceiling(real(upper_bound)/Ls(1)%rows)

        lower_offset = (lower_bound - (lower_block - 1)*Ls(1)%rows) - 1
        upper_offset = (upper_bound - (upper_block - 1)*Ls(1)%rows) - Ls(1)%rows

        local_rows = upper_bound - lower_bound + 1

        allocate(row_nnzs(lower_bound:upper_bound))


        do ii = 1, size(Ls, 1)

            row_nnzs = 0

            L => Ls(ii)

            LSO1%rows = L%rows**2
            LSO1%columns = L%columns**2
            allocate(LSO1%col_indexes(size(L%col_indexes)**2 + L%rows*LSO1%rows))
            allocate(LSO1%row_indexes(size(L%col_indexes)**2 + L%rows*LSO1%rows))
            allocate(LSO1%values(size(L%col_indexes)**2 + L%rows*LSO1%rows))

            indx = 1

            do i = lower_bound, upper_bound
                L_row = i - (ceiling(i/float(L%rows)) - 1)*L%rows
                L_nz = L%row_starts(L_row + 1) - L%row_starts(L_row)
                current_block = ceiling(real(i)/L%rows)
                col_ind_disp = 0
                row_nz = L%row_starts(current_block + 1) - L%row_starts(current_block)
                do j = 1, L%columns
                   if (L%row_starts(current_block + 1) - L%row_starts(current_block) == 0) exit
                    if (L%col_indexes(L%row_starts(current_block) + col_ind_disp) == j) then
                        row_nnzs(i) = row_nnzs(i) + L%row_starts(L_row + 1) - L%row_starts(L_row)
                        do k = L%row_starts(L_row), L%row_starts(L_row + 1) - 1
                            LSO1%row_indexes(indx) = i
                            LSO1%col_indexes(indx)  = (j - 1)*L%rows + L%col_indexes(k)
                            LSO1%values(indx) = conjg(L%values(L%row_starts(current_block) &
                                                    + col_ind_disp)) * L%values(k)
                            LSO1%nnz = LSO1%nnz + 1
                            indx = indx + 1
                        enddo
                        col_ind_disp = col_ind_disp + 1
                        if (col_ind_disp == row_nz) exit
                    endif
                enddo
            enddo

            call CSR_Transpose(L, L_T)

            allocate(LSO2_elements%row_indexes(L%rows**2))
            allocate(LSO2_elements%col_indexes(L%rows**2))
            allocate(LSO2_elements%values(L%rows**2))

            lb_indx = 0
            ub_indx = 0

            do i = 1, L%rows
                do j = 1, L%rows
                    value = 0
                    do k = L_T%row_starts(i), L_T%row_starts(i + 1) - 1
                        do ll = L_T%row_starts(j), L_T%row_starts(j + 1) - 1
                            if(L_T%col_indexes(ll) == L_T%col_indexes(k)) then
                                value = value + conjg(L_T%values(k))*L_T%values(ll)
                            endif
                            if(L_T%col_indexes(ll) > L_T%col_indexes(k)) exit
                        enddo
                    enddo
                    if (abs(value) > epsilon(real(value))) then
                            LSO2_elements%nnz = LSO2_elements%nnz + 1
                            LSO2_elements%row_indexes(LSO2_elements%nnz) = i
                            LSO2_elements%col_indexes(LSO2_elements%nnz)  = j
                            LSO2_elements%values(LSO2_elements%nnz) = value
                        if(((1 + lower_offset) == i) &
                            .and. (lb_indx == 0)) then
                            lb_indx = LSO2_elements%nnz
                        endif
                    endif

                    if(((L%rows + upper_offset) == i) &
                            .and. (L%rows == j)) then
                        ub_indx = LSO2_elements%nnz
                    endif
                enddo
            enddo

        if (lb_indx == 0) then
            lb_indx = 1
        endif

        allocate(LSO2%row_indexes(LSO2_elements%nnz*(upper_block - lower_block + 1)*2))
        allocate(LSO2%col_indexes(LSO2_elements%nnz*(upper_block - lower_block + 1)*2))
        allocate(LSO2%values(LSO2_elements%nnz*(upper_block - lower_block + 1)*2))

        indx = 1
        do i = lower_block, upper_block
            if (i == lower_block) then
                offset = lb_indx
            else
                offset = 1
            endif
            do j = offset, LSO2_elements%nnz
                row_index = (i-1)*L%rows + LSO2_elements%row_indexes(j)
                col_index = (i-1)*L%columns + LSO2_elements%col_indexes(j)

                LSO2%nnz = LSO2%nnz + 1

                LSO2%row_indexes(LSO2%nnz) = row_index
                LSO2%col_indexes(LSO2%nnz) = col_index
                LSO2%values(LSO2%nnz) = -0.5d0 * LSO2_elements%values(j)
                if ((i == upper_block) .and. (j == ub_indx)) exit

            enddo
        enddo

            call COO_Allocate_Sum(LSO1, LSO2, local_rows, LSO_temp_1)
            call COO_Sum(LSO1, LSO2, LSO_temp_1, partition_table, MPI_communicator)

            call COO_Deallocate(LSO1)
            call COO_Deallocate(LSO2)
            call COO_Deallocate(LSO2_elements)


            allocate(LSO3_elements%row_indexes(L%rows**2))
            allocate(LSO3_elements%col_indexes(L%rows**2))
            allocate(LSO3_elements%values(L%rows**2))

            do i = lower_block, upper_block
                do j = 1, L%rows
                    value = 0
                    do k = L_T%row_starts(i), L_T%row_starts(i + 1) - 1
                        do ll = L_T%row_starts(j), L_T%row_starts(j + 1) - 1
                            if(L_T%col_indexes(ll) == L_T%col_indexes(k)) then
                                value = value + L_T%values(k)*conjg(L_T%values(ll))
                            endif
                            if(L_T%col_indexes(ll) > L_T%col_indexes(k)) exit
                        enddo
                    enddo
                    if (abs(value) > epsilon(real(value))) then
                        LSO3_elements%nnz = LSO3_elements%nnz + 1
                        LSO3_elements%row_indexes(LSO3_elements%nnz) = i
                        LSO3_elements%col_indexes(LSO3_elements%nnz) = j
                        LSO3_elements%values(LSO3_elements%nnz) = value
                    endif
                enddo
            enddo

            allocate(LSO3%row_indexes(LSO3_elements%nnz*L%rows))
            allocate(LSO3%col_indexes(LSO3_elements%nnz*L%rows))
            allocate(LSO3%values(LSO3_elements%nnz*L%rows))

            indx = 1
            do i = lower_block, upper_block
                cnt = 0
                row_indx = indx
                do while(indx <= LSO3_elements%nnz)
                    if (LSO3_elements%row_indexes(indx) /= i) exit
                    row_index = (LSO3_elements%row_indexes(indx) - 1) * L%rows + 1 + lower_offset * &
                        Kronecker_Delta(i, lower_block)
                    col_index = (LSO3_elements%col_indexes(indx) - 1) * L%columns + 1 + lower_offset * &
                        Kronecker_Delta(i, lower_block)
                    value = -0.5d0 * LSO3_elements%values(indx)
                    LSO3%nnz = LSO3%nnz + 1
                    LSO3%row_indexes(LSO3%nnz) = row_index
                    LSO3%col_indexes(LSO3%nnz) = col_index
                    LSO3%values(LSO3%nnz) = value
                    cnt = cnt + 1
                    indx = indx + 1
                enddo
                if (cnt /= 0) then
                    do j = 2 + Kronecker_Delta(i, lower_block)*lower_offset, L%rows &
                        + Kronecker_Delta(i, upper_block)*upper_offset
                    indx = row_indx
                    row_index = (LSO3_elements%row_indexes(indx) - 1) * L%rows + j
                        do k = 1, cnt
                            col_index = (LSO3_elements%col_indexes(indx) - 1) * L%columns + j
                            value = -0.5d0 * LSO3_elements%values(indx)
                            LSO3%nnz = LSO3%nnz + 1
                            LSO3%row_indexes(LSO3%nnz) = row_index
                            LSO3%col_indexes(LSO3%nnz) = col_index
                            LSO3%values(LSO3%nnz) = value
                            indx = indx + 1
                        enddo
                    enddo
                endif
            enddo

            call COO_Deallocate(LSO3_elements)

            call COO_Allocate_Sum(LSO_temp_1, LSO3, local_rows, LSO_temp_2)
            call COO_Sum(LSO_temp_1, LSO3, LSO_temp_2, partition_table, MPI_communicator)

            call COO_Deallocate(LSO_temp_1)
            call COO_Deallocate(LSO3)

            if (size(Ls) > 1) then

                if (ii == 1) then

                    call COO_Move(LSO_temp_2, LSO)
                    LSO%values = taus(ii)*LSO%values

                else

                    call COO_Move(LSO, LSO_temp_3)

                    call COO_Allocate_Sum(LSO_temp_2, LSO_temp_3, local_rows, LSO)
                    call COO_Sum(LSO_temp_2, LSO_temp_3, LSO, partition_table, MPI_communicator)

                    call COO_Deallocate(LSO_temp_2)
                    call COO_Deallocate(LSO_temp_3)

                    LSO%values = taus(ii)*LSO%values

                endif

            else

                call COO_Move(LSO_temp_2, LSO)

                LSO%values = taus(ii)*LSO%values

            endif

        enddo

    end subroutine COO_Vectorized_Dissipator

!>    Subroutine: Generalized_Superoperator
!
!>    @brief Returns the superoperator resulting from vectorization of the Lindblad master equation given
!>    arbitrary Lindblad operators.
!
!>    @details Returns as a distributed CSR array:
!>
!>>   \text{i}(I \otimes H - H^T \otimes I)
!>>   + \tau_k \sum_{k=1}^K(L_k^* \otimes L_k - \frac{1}{2}(I \otimes L_k^\dagger L_k + L_k^T L_k^* \otimes I))

    subroutine Generalized_Superoperator(   taus, &
                                            H, &
                                            Ls, &
                                            partition_table, &
                                            SO, &
                                            MPI_communicator)

        real(dp), dimension(:), intent(in) :: taus !> @param Lindblad operator coefficients.
        type(CSR), intent(in) :: H !> CSR matrix Hamiltonian.
        type(CSR), dimension(:), target, intent(in) :: Ls !> CSR array of Lindblad operators.
        integer, dimension(:), intent(in) :: partition_table !> Distributed row partitions.
        type(CSR), intent(inout) :: SO !> CSR matrix superoperator.
        integer, intent(in) :: MPI_communicator !> MPI communicator over which to generate SO.

        integer :: local_rows

        type(COO) :: HSO
        type(COO) :: LSO

        type(COO) :: SO_COO

        ! MPI
        integer :: rank
        integer :: ierr

        call MPI_comm_rank(MPI_communicator, rank, ierr)

        local_rows = partition_table(rank + 2) - partition_table(rank +1)

        call COO_Vectorized_Commutator( H, &
                                        partition_table, &
                                        HSO, &
                                        MPI_communicator)

        HSO%values = cmplx(0.0d0, -1.0d0,dp) * HSO%values

        call COO_Vectorized_Dissipator( taus, &
                                        Ls, &
                                        partition_table, &
                                        LSO, &
                                        MPI_communicator)
        call COO_Allocate_Sum(HSO, LSO, local_rows, SO_COO)
        call COO_Sum(HSO, LSO, SO_COO, partition_table, MPI_communicator)
        call COO_Deallocate(HSO)
        call COO_Deallocate(LSO)

        call COO_to_CSR(SO_COO, SO, partition_table, MPI_communicator)

    end subroutine Generalized_Superoperator

!>    Subroutine: Demoralized_Superoperator
!
!>    @brief Returns the superoperator resulting from vectorization of the non-moralizing
!>    quantum stochastic walk mater equation.
!
!>    @details Returns, as a CSR array,
!>
!>>    \tilde{\mathcal{L}} = (I \otimes H - H^T \otimes H) + 
!>>    + \omega ((I \otimes H_\text{loc} - H_\text{loc}^T \otimes H_\text{loc})
!>>    + \sum_{k=1}^K(L_k^* \otimes L_k - \frac{1}{2}(I \otimes L_k^\dagger L_k + L_k^T L_k^* \otimes I)))
!>
!>    Where H, H_\text{loc} and L_k are CSR arrays corresponding to a demoralised graph.

    subroutine Demoralized_Superoperator(   omega, &
                                            H, &
                                            H_loc, &
                                            Ls, &
                                            partition_table, &
                                            SO, &
                                            MPI_communicator)

        real(dp), intent(in) :: omega
        type(CSR), intent(in) :: H
        type(CSR), intent(in) :: H_loc
        type(CSR), dimension(:), target, intent(in) :: Ls
        integer, dimension(:), intent(in) :: partition_table
        type(CSR), intent(inout) :: SO
        integer, intent(in) :: MPI_communicator

        type(COO), target :: HSO
        type(COO), target :: HSO_loc
        type(COO), target :: LSO
        type(COO), pointer :: SO_COO_temp => null()
        type(COO), pointer :: SO_COO => null()

        real(dp), dimension(:), allocatable :: taus

        integer :: local_rows

        ! MPI
        integer :: ierr
        integer :: rank

        call MPI_comm_rank(MPI_communicator, rank, ierr)

        local_rows = partition_table(rank + 2) - partition_table(rank + 1)

        allocate(taus(size(Ls, 1)))

        taus = omega

        if ((1.d0-omega) > epsilon(omega)) then

            call COO_Vectorized_Commutator( H, &
                                            partition_table, &
                                            HSO, &
                                            MPI_communicator)

            HSO%values = cmplx(0.0d0, omega - 1.0d0, dp) * HSO%values

        endif

        if (omega > epsilon(omega)) then

            call COO_Vectorized_Dissipator( taus, &
                                            Ls, &
                                            partition_table, &
                                            LSO, &
                                            MPI_communicator)

            if (H_loc%row_starts(H_loc%rows + 1) > 1) then

                call COO_Vectorized_Commutator( H_loc, &
                                                partition_table, &
                                                HSO_loc, &
                                                MPI_communicator)

                HSO_loc%values = -cmplx(0.0d0, omega, dp) * HSO_loc%values

                allocate(SO_COO_temp)
                call COO_Allocate_Sum(HSO_loc, LSO, local_rows, SO_COO_temp)
                call COO_Sum(HSO_loc, LSO, SO_COO_temp, partition_table, MPI_communicator)
                call COO_Deallocate(HSO_loc)
                call COO_Deallocate(LSO)

            else

               SO_COO_temp => LSO

            endif

        endif

        if (((1.d0 - omega) > epsilon(omega)) .and. (omega > epsilon(omega))) then
            allocate(SO_COO)
            call COO_Allocate_Sum(SO_COO_temp, HSO, local_rows, SO_COO)
            call COO_Sum(SO_COO_temp, HSO, SO_COO, partition_table, MPI_communicator)
            call COO_Deallocate(SO_COO_temp)
            call COO_Deallocate(HSO)
        elseif ((1.d0 - omega) > epsilon(omega)) then
            SO_COO => HSO
        elseif (omega > epsilon(omega)) then
            SO_COO => SO_COO_temp
        endif

        call COO_to_CSR(SO_COO, SO, partition_table, MPI_communicator)

        call COO_Deallocate(SO_COO)

    end subroutine Demoralized_Superoperator

    !
    !   Subroutine: Vectorize_Operator
    !
    !>  @brief Vectorise a N x N operator.
    !
    !>  @details Given a dense operator O, the subroutine vectorizes the operator
    !>  such that its vectorized for is distributed over a MPI communicator row-wise
    !>  as described by partition_table.

    subroutine Vectorize_Operator(  O, &
                                    partition_table, &
                                    O_v, &
                                    MPI_communicator)

        complex(dp), dimension(:, :), intent(in) :: O !< @param N x N operator.
        integer, dimension(:), intent(in) :: partition_table !< @param Row-wise partitioning of the vectorized operator.
        complex(dp), dimension(:), allocatable, intent(out) :: O_v !< @param Local partition of the vectorized operator.
        integer, intent(in) :: MPI_communicator !< @param MPI communicator over which the vectorized operator is to be distributed.

        integer :: O_size

        integer :: L_bound, U_bound
        integer :: L_row, U_row
        integer :: L_shift, U_shift

        integer :: i, j
        integer :: indx

        !MPI ENV
        integer :: rank
        integer :: ierr

        call MPI_comm_rank(MPI_communicator, rank, ierr)

        L_bound = partition_table(rank + 1)
        U_bound = partition_table(rank + 2) - 1

        allocate(O_v(L_bound:U_bound))

        O_size = size(O, 1)

        L_row = (L_bound - mod(L_bound, O_size) - &
            Kronecker_Delta(1, mod(L_bound, O_size)+1)*O_size)/O_size + 1

        U_row = (U_bound - mod(U_bound, O_size) - &
            Kronecker_Delta(1, mod(U_bound, O_size)+1)*O_size)/O_size + 1

        L_shift = L_bound - (L_row - 1)*O_size - 1
        U_shift = U_row*O_size - (U_bound)

        do i = L_row, U_row

            do j = 1 + Kronecker_Delta(i, L_row)*L_shift, &
                O_size - Kronecker_Delta(i, U_row)*U_shift

                indx = (i - 1)*size(O, 1) + j

                O_v(indx) = O(i, j)

            enddo

        enddo

    end subroutine Vectorize_Operator

    !
    !   Subroutine: Reshape_Vectorized_Operator_Series
    !
    !>  @breif Reshape a series resulting from a vectorized operator of dimension
    !>  N^2 into a series consisting of N x N matrices.

    subroutine Reshape_Vectorized_Operator_Series( O_v_series, O_series)

        complex(dp), dimension(:,:), intent(in) :: O_v_series !< @param Vectorized operator series.
        complex(dp), dimension(:, :, :), intent(out) :: O_series !< @param reshaped operator series.

        integer :: O_dim
        integer :: steps
        integer :: i, j, k
        integer :: indx

        O_dim = int(sqrt(real(size(O_v_series, 1))))
        steps = size(O_v_series, 2)


        do k = 1, steps
            do i = 1, O_dim
                do j = 1, O_dim
                    indx = (i - 1)*O_dim + j
                    O_series(i, j, k) = O_v_series(indx, k)
                enddo
            enddo
        enddo

    end subroutine Reshape_Vectorized_Operator_Series

    !
    !   Subroutine: Reshape_Vectorized_Operator
    !
    !>  @brief Rehape a vectorized operator of dimension N^2 into a N x N matrix.

    subroutine Reshape_Vectorized_Operator( O_v, O)

        complex(dp), dimension(:), intent(in) :: O_v !< @param Vectorized operator.
        complex(dp), dimension(:, :), intent(out) :: O !< @param Reshaped operator.

        integer :: O_dim
        integer :: i, j
        integer :: indx

        O_dim = int(sqrt(real(size(O_v))))

        do i = 1, O_dim
            do j = 1, O_dim
                indx = (i - 1)*O_dim + j
                O(i, j) = O_v(indx)
            enddo
        enddo

    end subroutine Reshape_Vectorized_Operator

end module Operators
