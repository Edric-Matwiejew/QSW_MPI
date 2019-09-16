subroutine  Super_Operator_Extent( omega, &
                            H_rows, &
                            H_nnz, &
                            H_row_starts, &
                            H_col_indexes, &
                            H_values, &
                            L_nnz, &
                            L_row_starts, &
                            L_col_indexes, &
                            L_values, &
                            n_sources, &
                            source_sites, &
                            source_rates, &
                            n_sinks, &
                            sink_sites, &
                            sink_rates, &
                            flock, &
                            MPI_communicator, &
                            M_nnz_out, &
                            M_rows, &
                            partition_table)

        use :: Operators
        use :: MPI

        implicit none

        real(8), intent(in) :: omega
        integer, intent(in) :: H_rows
        integer, intent(in) :: H_nnz
        integer, dimension(H_rows + 1), intent(in) :: H_row_starts
        integer, dimension(H_nnz), intent(in) :: H_col_indexes
        complex(8), dimension(H_nnz), intent(in) :: H_values
        integer, intent(in) :: L_nnz
        integer, dimension(H_rows + 1), intent(in) :: L_row_starts
        integer, dimension(L_nnz), intent(in) :: L_col_indexes
        complex(8), dimension(L_nnz), intent(in) :: L_values
        integer, intent(in) :: n_sources
        integer, dimension(n_sources), intent(in) :: source_sites
        real(8), dimension(n_sources), target, intent(in) :: source_rates
        integer, intent(in) :: n_sinks
        integer, dimension(n_sinks), intent(in) :: sink_sites
        real(8), dimension(n_sinks), target, intent(in) :: sink_rates
        integer, intent(in) :: flock
        integer, intent(in) :: MPI_communicator
        integer, intent(out) :: M_nnz_out
        integer, intent(out) :: M_rows
        integer, dimension(flock + 1), intent(out) :: partition_table

        type(CSR) :: H, L, M

        integer, dimension(:), allocatable :: source_sites_temp
        real(8), dimension(:), pointer :: source_rates_temp
        integer, dimension(:), allocatable :: sink_sites_temp
        real(8), dimension(:), pointer :: sink_rates_temp
        integer, dimension(:), allocatable :: partition_table_temp

        integer :: ierr
        integer :: i

!f2py  real(kind=8) intent(in) :: omega
!f2py  integer, optional,intent(in),check((len(h_row_starts)-1)>=h_rows),depend(h_row_starts) :: h_rows=(len(h_row_starts)-1)
!f2py  integer, optional,intent(in),check(len(h_col_indexes)>=h_nnz),depend(h_col_indexes) :: h_nnz=len(h_col_indexes)
!f2py  integer dimension(h_rows + 1),intent(in) :: h_row_starts
!f2py  integer dimension(h_nnz),intent(in) :: h_col_indexes
!f2py  complex(kind=8) dimension(h_nnz),intent(in),depend(h_nnz) :: h_values
!f2py  integer, optional,intent(in),check(len(l_col_indexes)>=l_nnz),depend(l_col_indexes) :: l_nnz=len(l_col_indexes)
!f2py  integer dimension(h_rows + 1),intent(in),depend(h_rows) :: l_row_starts
!f2py  integer dimension(l_nnz),intent(in) :: l_col_indexes
!f2py  complex(kind=8) dimension(l_nnz),intent(in),depend(l_nnz) :: l_values
!f2py  integer, optional,intent(in),check(len(source_sites)>=n_sources),depend(source_sites) :: n_sources=len(source_sites)
!f2py  integer dimension(n_sources),intent(in) :: source_sites
!f2py  real(kind=8), target,dimension(n_sources),intent(in),depend(n_sources) :: source_rates
!f2py  integer, optional,intent(in),check(len(sink_sites)>=n_sinks),depend(sink_sites) :: n_sinks=len(sink_sites)
!f2py  integer dimension(n_sinks),intent(in) :: sink_sites
!f2py  real(kind=8), target,dimension(n_sinks),intent(in),depend(n_sinks) :: sink_rates
!f2py  integer intent(in) :: flock
!f2py  integer intent(in) :: mpi_communicator
!f2py  integer intent(out) :: m_nnz_out
!f2py  integer intent(out) :: m_rows
!f2py  integer dimension(flock + 1),intent(out),depend(flock) :: partition_table

        if (source_sites(1) < 0) then
            allocate(source_sites_temp(0))
            allocate(source_rates_temp(0))
        else

            allocate(source_sites_temp(n_sources))

            !$omp parallel do
            do i = 1, n_sources
                source_sites_temp(i) = source_sites(i) + 1
            enddo
            !$omp end parallel do

            source_rates_temp => source_rates

        endif

        if (sink_sites(1) < 0) then
            allocate(sink_sites_temp(0))
            allocate(sink_rates_temp(0))
        else

            allocate(sink_sites_temp(n_sinks))

            !$omp parallel do
            do i = 1, n_sinks
                sink_sites_temp(i) = sink_sites(i) + 1
            enddo
            !$omp end parallel do

            sink_rates_temp => sink_rates

        endif

        allocate(H%row_starts(H_rows + 1))
        allocate(H%col_indexes(H_nnz))
        allocate(H%values(H_nnz))

        H%rows = H_rows
        H%columns = H_rows

        !$omp parallel do
        do i = 1, H_rows + 1
            H%row_starts(i) = H_row_starts(i) + 1
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, H_nnz
            H%col_indexes(i) = H_col_indexes(i) + 1
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, H_nnz
            H%values(i) = H_values(i)
        enddo
        !$omp end parallel do

        allocate(L%row_starts(H_rows + 1))
        allocate(L%col_indexes(L_nnz))
        allocate(L%values(L_nnz))

        L%rows = H_rows
        L%columns = H_rows

        !$omp parallel do
        do i = 1, H_rows + 1
            L%row_starts(i) = L_row_starts(i) + 1
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, L_nnz
            L%col_indexes(i) = L_col_indexes(i) + 1
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, L_nnz
            L%values(i) = L_values(i)
        enddo
        !$omp end parallel do

        call Prepare_Super_Operator(    omega, &
                                        H, &
                                        L, &
                                        source_sites_temp, &
                                        source_rates_temp, &
                                        sink_sites_temp, &
                                        sink_rates_temp, &
                                        partition_table_temp, &
                                        M, &
                                        MPI_communicator)

        M_nnz_out = size(M%col_indexes)
        M_rows = M%columns
        partition_table = partition_table_temp

        deallocate(M%values, M%row_starts, M%col_indexes)
        deallocate(L%values, L%row_starts, L%col_indexes)
        deallocate(H%values, H%row_starts, H%col_indexes)


end subroutine Super_Operator_Extent

subroutine  Super_Operator( omega, &
                            H_rows, &
                            H_nnz, &
                            H_row_starts, &
                            H_col_indexes, &
                            H_values, &
                            L_nnz, &
                            L_row_starts, &
                            L_col_indexes, &
                            L_values, &
                            n_sources, &
                            source_sites, &
                            source_rates, &
                            n_sinks, &
                            sink_sites, &
                            sink_rates, &
                            M_nnz_in, &
                            M_n_row_starts, &
                            flock, &
                            MPI_communicator, &
                            M_row_starts, &
                            M_col_indexes, &
                            M_values, &
                            mu)

        use :: Operators
        use :: Expm
        use :: MPI

        implicit none

        real(8), intent(in) :: omega
        integer, intent(in) :: H_rows
        integer, intent(in) :: H_nnz
        integer, dimension(H_rows + 1), intent(in) :: H_row_starts
        integer, dimension(H_nnz), intent(in) :: H_col_indexes
        complex(8), dimension(H_nnz), intent(in) :: H_values
        integer, intent(in) :: L_nnz
        integer, dimension(H_rows + 1), intent(in) :: L_row_starts
        integer, dimension(L_nnz), intent(in) :: L_col_indexes
        complex(8), dimension(L_nnz), intent(in) :: L_values
        integer, intent(in) :: n_sources
        integer, dimension(n_sources), intent(in) :: source_sites
        real(8), dimension(n_sources), target, intent(in) :: source_rates
        integer, intent(in) :: n_sinks
        integer, dimension(n_sinks), intent(in) :: sink_sites
        real(8), dimension(n_sinks), target, intent(in) :: sink_rates
        integer, intent(in) :: M_n_row_starts
        integer, intent(in) :: M_nnz_in
        integer, intent(in) :: flock
        integer, intent(in) :: MPI_communicator
        integer, dimension(M_n_row_starts), intent(out) :: M_row_starts
        integer, dimension(M_nnz_in), intent(out) :: M_col_indexes
        complex(8), dimension(M_nnz_in), intent(out) :: M_values
        complex(8), intent(out) :: mu

        type(CSR) :: H, L, M

        integer, dimension(:), allocatable :: source_sites_temp
        real(8), dimension(:), pointer :: source_rates_temp
        integer, dimension(:), allocatable :: sink_sites_temp
        real(8), dimension(:), pointer :: sink_rates_temp
        integer, dimension(:), allocatable :: partition_table_temp

        integer :: ierr
        integer :: i

!f2py real(kind=8) intent(in) :: omega
!f2py integer, optional,intent(in),check((len(h_row_starts)-1)>=h_rows),depend(h_row_starts) :: h_rows=(len(h_row_starts)-1)
!f2py integer, optional,intent(in),check(len(h_col_indexes)>=h_nnz),depend(h_col_indexes) :: h_nnz=len(h_col_indexes)
!f2py integer dimension(h_rows + 1),intent(in) :: h_row_starts
!f2py integer dimension(h_nnz),intent(in) :: h_col_indexes
!f2py complex(kind=8) dimension(h_nnz),intent(in),depend(h_nnz) :: h_values
!f2py integer, optional,intent(in),check(len(l_col_indexes)>=l_nnz),depend(l_col_indexes) :: l_nnz=len(l_col_indexes)
!f2py integer dimension(h_rows + 1),intent(in),depend(h_rows) :: l_row_starts
!f2py integer dimension(l_nnz),intent(in) :: l_col_indexes
!f2py complex(kind=8) dimension(l_nnz),intent(in),depend(l_nnz) :: l_values
!f2py integer, optional,intent(in),check(len(source_sites)>=n_sources),depend(source_sites) :: n_sources=len(source_sites)
!f2py integer dimension(n_sources),intent(in) :: source_sites
!f2py real(kind=8), target,dimension(n_sources),intent(in),depend(n_sources) :: source_rates
!f2py integer, optional,intent(in),check(len(sink_sites)>=n_sinks),depend(sink_sites) :: n_sinks=len(sink_sites)
!f2py integer dimension(n_sinks),intent(in) :: sink_sites
!f2py real(kind=8), target,dimension(n_sinks),intent(in),depend(n_sinks) :: sink_rates
!f2py integer intent(in) :: m_nnz_in
!f2py integer intent(in) :: m_n_row_starts
!f2py integer intent(in) :: flock
!f2py integer intent(in) :: mpi_communicator
!f2py integer dimension(m_n_row_starts),intent(out),depend(m_n_row_starts) :: m_row_starts
!f2py integer dimension(m_nnz_in),intent(out),depend(m_nnz_in) :: m_col_indexes
!f2py complex(kind=8) dimension(m_nnz_in),intent(out),depend(m_nnz_in) :: m_values
!f2py complex(kind=8), intent(out) :: mu

        if (source_sites(1) == -1) then
            allocate(source_sites_temp(0))
            allocate(source_rates_temp(0))
        else

            allocate(source_sites_temp(n_sources))

            !$omp parallel do
            do i = 1, n_sources
                source_sites_temp(i) = source_sites(i) + 1
            enddo
            !$omp end parallel do

            source_rates_temp => source_rates

        endif

        if (sink_sites(1) == -1) then
            allocate(sink_sites_temp(0))
            allocate(sink_rates_temp(0))
        else

            allocate(sink_sites_temp(n_sinks))

            !$omp parallel do
            do i = 1, n_sinks
                sink_sites_temp(i) = sink_sites(i) + 1
            enddo
            !$omp end parallel do

            sink_rates_temp => sink_rates

        endif

        allocate(H%row_starts(H_rows + 1))
        allocate(H%col_indexes(H_nnz))
        allocate(H%values(H_nnz))

        H%rows = H_rows
        H%columns = H_rows

        !$omp parallel do
        do i = 1, H_rows + 1
            H%row_starts(i) = H_row_starts(i) + 1
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, H_nnz
            H%col_indexes(i) = H_col_indexes(i) + 1
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, H_nnz
            H%values(i) = H_values(i)
        enddo
        !$omp end parallel do

        allocate(L%row_starts(H_rows + 1))
        allocate(L%col_indexes(L_nnz))
        allocate(L%values(L_nnz))

        L%rows = H_rows
        L%columns = H_rows

        !$omp parallel do
        do i = 1, H_rows + 1
            L%row_starts(i) = L_row_starts(i) + 1
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, L_nnz
            L%col_indexes(i) = L_col_indexes(i) + 1
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, L_nnz
            L%values(i) = L_values(i)
        enddo
        !$omp end parallel do

        call Prepare_Super_Operator(    omega, &
                                        H, &
                                        L, &
                                        source_sites_temp, &
                                        source_rates_temp, &
                                        sink_sites_temp, &
                                        sink_rates_temp, &
                                        partition_table_temp, &
                                        M, &
                                        MPI_communicator)

        call Norm_Reduction(M, mu, MPI_communicator)

        M_row_starts(1:size(M%row_starts)) = M%row_starts
        M_col_indexes(1:size(M%col_indexes)) = M%col_indexes
        M_values(1:size(M%col_indexes)) = M%values

        deallocate(L%values, L%row_starts, L%col_indexes)
        deallocate(H%values, H%row_starts, H%col_indexes)

end subroutine Super_Operator

subroutine rec_a(   M_rows, &
                    M_n_row_starts, &
                    M_n_col_indexes, &
                    M_row_starts, &
                    M_col_indexes, &
                    flock, &
                    partition_table, &
                    MPI_communicator, &
                    M_num_rec_inds, &
                    M_rec_disps, &
                    M_num_send_inds, &
                    M_send_disps)

        use :: Sparse

        implicit none

        integer, intent(in) :: M_rows
        integer, intent(in) :: M_n_row_starts
        integer, intent(in) :: M_n_col_indexes
        integer, dimension(M_n_row_starts), target, intent(in) :: M_row_starts
        integer, dimension(M_n_col_indexes), target, intent(in) :: M_col_indexes
        integer, intent(in) :: flock
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(in) :: MPI_communicator
        integer, dimension(flock), target, intent(out) :: M_num_rec_inds
        integer, dimension(flock), target, intent(out) :: M_rec_disps
        integer, dimension(flock), target, intent(out) :: M_num_send_inds
        integer, dimension(flock), target, intent(out) :: M_send_disps

        type(CSR) :: M

        integer :: rank, ierr

        integer :: lb, ub
        integer :: lb_elements, ub_elements

!f2py integer intent(in) :: m_rows
!f2py integer, optional,intent(in),check(len(m_col_indexes)>=m_n_col_indexes),depend(m_col_indexes) :: m_n_col_indexes=len(m_col_indexes)
!f2py integer, optional,intent(in),check(len(m_row_starts)>=m_n_row_starts),depend(m_row_starts) :: m_n_row_starts=len(m_row_starts)
!f2py integer, target,dimension(m_n_row_starts),intent(in) :: m_row_starts
!f2py integer, target,dimension(m_n_col_indexes),intent(in) :: m_col_indexes
!f2py integer, optional,intent(in),check((len(partition_table)-1)>=flock),depend(partition_table) :: flock=(len(partition_table)-1)
!f2py integer dimension(flock + 1),intent(in) :: partition_table
!f2py integer intent(in) :: mpi_communicator
!f2py integer, target,dimension(flock),intent(out),depend(flock) :: m_num_rec_inds
!f2py integer, target,dimension(flock),intent(out),depend(flock) :: m_rec_disps
!f2py integer, target,dimension(flock),intent(out),depend(flock) :: m_num_send_inds
!f2py integer, target,dimension(flock),intent(out),depend(flock) :: m_send_disps

        call mpi_comm_rank(mpi_communicator, rank, ierr)

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        m%row_starts(lb:ub + 1) => m_row_starts

        lb_elements = m%row_starts(lb)
        ub_elements = m%row_starts(ub + 1) - 1

        M%rows = M_rows
        M%columns = M_rows
        M%col_indexes(lb_elements:ub_elements) => M_col_indexes
        M%num_rec_inds(1:size(M_num_rec_inds)) => M_num_rec_inds
        M%rec_disps(1:size(M_rec_disps)) => M_rec_disps
        M%num_send_inds(1:size(M_num_send_inds)) => M_num_send_inds
        M%send_disps(1:size(M_send_disps)) => M_send_disps

        call Reconcile_Communications_A(M, partition_table, MPI_communicator)

end subroutine rec_a

subroutine rec_b(   M_rows, &
                    M_nnz, &
                    M_n_row_starts, &
                    num_send, &
                    M_row_starts, &
                    M_col_indexes, &
                    M_num_rec_inds, &
                    M_rec_disps, &
                    M_num_send_inds, &
                    M_send_disps, &
                    flock, &
                    partition_table, &
                    MPI_communicator, &
                    M_local_col_inds, &
                    M_RHS_send_inds)

        use :: Sparse

        implicit none

        integer, intent(in) :: M_rows
        integer, intent(in) :: M_nnz
        integer, intent(in) :: M_n_row_starts
        integer, intent(in) :: num_send
        integer, dimension(M_n_row_starts), target, intent(in) :: M_row_starts
        integer, dimension(M_nnz), target, intent(in) :: M_col_indexes
        integer, dimension(flock), target, intent(in) :: M_rec_disps
        integer, dimension(flock), target, intent(in) :: M_num_send_inds
        integer, dimension(flock), target, intent(in) :: M_num_rec_inds
        integer, dimension(flock), target, intent(in) :: M_send_disps
        integer, intent(in) :: flock
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(in) :: MPI_communicator
        integer, dimension(M_nnz), target, intent(out) :: M_local_col_inds
        integer, dimension(num_send), target, intent(out) :: M_RHS_send_inds

        type(CSR) :: M

        integer :: rank, ierr

        integer :: lb, ub
        integer :: lb_elements, ub_elements


!f2py integer intent(in) :: m_rows
!f2py integer, optional,intent(in),check(len(m_col_indexes)>=m_nnz),depend(m_col_indexes) :: m_nnz=len(m_col_indexes)
!f2py integer, optional,intent(in),check(len(m_row_starts)>=m_n_row_starts),depend(m_row_starts) :: m_n_row_starts=len(m_row_starts)
!f2py integer intent(in) :: num_send
!f2py integer, target,dimension(m_n_row_starts),intent(in) :: m_row_starts
!f2py integer, target,dimension(m_nnz),intent(in) :: m_col_indexes
!f2py integer, target,dimension(flock),intent(in) :: m_num_rec_inds
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_rec_disps
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_num_send_inds
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_send_disps
!f2py integer, optional,intent(in),check(len(m_num_rec_inds)>=flock),depend(m_num_rec_inds) :: flock=len(m_num_rec_inds)
!f2py integer dimension(flock + 1),intent(in),depend(flock) :: partition_table
!f2py integer intent(in) :: mpi_communicator
!f2py integer, target,dimension(m_nnz),intent(out),depend(m_nnz) :: m_local_col_inds
!f2py integer, target,dimension(num_send),intent(out),depend(num_send) :: m_rhs_send_inds

        call mpi_comm_rank(mpi_communicator, rank, ierr)

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        M%row_starts(lb:ub + 1) => m_row_starts

        lb_elements = M%row_starts(lb)
        ub_elements = M%row_starts(ub + 1) - 1

        M%rows = M_rows
        M%columns = M_rows
        M%col_indexes(lb_elements:ub_elements) => M_col_indexes
        M%local_col_inds(lb_elements:ub_elements) => M_local_col_inds
        M%num_rec_inds(1:size(M_num_rec_inds)) => M_num_rec_inds
        M%rec_disps(1:size(M_rec_disps)) => M_rec_disps
        M%num_send_inds(1:size(M_num_send_inds)) => M_num_send_inds
        M%send_disps(1:size(M_send_disps)) => M_send_disps
        M%RHS_send_inds(1:size(M_RHS_send_inds)) => M_RHS_send_inds

        call Reconcile_Communications_B(M, partition_table, MPI_communicator)

end subroutine rec_b

subroutine one_norm_series( M_rows, &
                            M_n_col_indexes, &
                            M_n_values, &
                            M_n_local_col_indexes, &
                            M_n_row_starts, &
                            M_sends, &
                            M_row_starts, &
                            M_col_indexes, &
                            M_values, &
                            M_num_rec_inds, &
                            M_rec_disps, &
                            M_num_send_inds, &
                            M_send_disps, &
                            M_local_col_inds, &
                            M_RHS_send_inds, &
                            flock, &
                            partition_table, &
                            MPI_communicator, &
                            one_norm_array, &
                            p)

!f2py integer intent(in) :: m_rows
!f2py integer, optional,intent(in),check(len(m_col_indexes)>=m_n_col_indexes),depend(m_col_indexes) :: m_n_col_indexes=len(m_col_indexes)
!f2py integer, optional,intent(in),check(len(m_values)>=m_n_values),depend(m_values) :: m_n_values=len(m_values)
!f2py integer, optional,intent(in),check(len(m_local_col_inds)>=m_n_local_col_indexes),depend(m_local_col_inds) :: m_n_local_col_indexes=len(m_local_col_inds)
!f2py integer, optional,intent(in),check(len(m_row_starts)>=m_n_row_starts),depend(m_row_starts) :: m_n_row_starts=len(m_row_starts)
!f2py integer, optional,intent(in),check(len(m_rhs_send_inds)>=m_sends),depend(m_rhs_send_inds) :: m_sends=len(m_rhs_send_inds)
!f2py integer, target,dimension(m_n_row_starts),intent(in) :: m_row_starts
!f2py integer, target,dimension(m_n_col_indexes),intent(in) :: m_col_indexes
!f2py complex(kind=8), target,dimension(m_n_values),intent(in) :: m_values
!f2py integer, target,dimension(flock),intent(in) :: m_num_rec_inds
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_rec_disps
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_num_send_inds
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_send_disps
!f2py integer, target,dimension(m_n_local_col_indexes),intent(in) :: m_local_col_inds
!f2py integer, target,dimension(m_sends),intent(in) :: m_rhs_send_inds
!f2py integer, optional,intent(in),check(len(m_num_rec_inds)>=flock),depend(m_num_rec_inds) :: flock=len(m_num_rec_inds)
!f2py integer dimension(flock + 1),intent(in),depend(flock) :: partition_table
!f2py integer intent(in) :: mpi_communicator
!f2py real(kind=8) dimension(9),intent(out) :: one_norm_array
!f2py integer intent(out) :: p

        use :: Sparse
        use :: One_Norms
        use :: MPI
        use :: Expm

        implicit none

        integer, intent(in) :: M_rows
        integer, intent(in) :: M_n_col_indexes
        integer, intent(in) :: M_n_values
        integer, intent(in) :: M_n_local_col_indexes
        integer, intent(in) :: M_n_row_starts
        integer, intent(in) :: M_sends
        integer, dimension(M_n_row_starts), target, intent(in) :: M_row_starts
        integer, dimension(M_n_col_indexes), target, intent(in) :: M_col_indexes
        complex(8), dimension(M_n_values), target, intent(in) :: M_values
        integer, dimension(flock), target, intent(in) :: M_num_rec_inds
        integer, dimension(flock), target, intent(in) :: M_rec_disps
        integer, dimension(flock), target, intent(in) :: M_num_send_inds
        integer, dimension(flock), target, intent(in) :: M_send_disps
        integer, dimension(M_n_local_col_indexes), target, intent(in) :: M_local_col_inds
        integer, dimension(M_sends), target, intent(in) :: M_RHS_send_inds
        integer, intent(in) :: flock
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(in) :: MPI_communicator
        real(8), dimension(9), intent(out) :: one_norm_array
        integer, intent(out) :: p

        type(CSR) :: M
        type(CSR) :: M_T

        integer :: itmax

        integer :: i

        integer :: rank, ierr

        integer :: lb, ub
        integer :: lb_elements, ub_elements

        integer, parameter :: l = 3
        integer, parameter :: pmax = 8

        real(8), dimension(9) :: alphas

        complex(8) :: mu

        call mpi_comm_rank(mpi_communicator, rank, ierr)

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        M%rows = M_rows
        M%columns = M_rows
        M%row_starts(lb:ub + 1) => M_row_starts

        lb_elements = M%row_starts(lb)
        ub_elements = M%row_starts(ub + 1) - 1

        M%col_indexes(lb_elements:ub_elements) => M_col_indexes
        M%values(lb_elements:ub_elements) => M_values
        M%local_col_inds(lb_elements:ub_elements) => M_local_col_inds
        M%num_rec_inds(1:size(M_num_rec_inds)) => M_num_rec_inds
        M%rec_disps(1:size(M_rec_disps)) => M_rec_disps
        M%num_send_inds(1:size(M_num_send_inds)) => M_num_send_inds
        M%send_disps(1:size(M_send_disps)) => M_send_disps
        M%RHS_send_inds(1:size(M_RHS_send_inds)) => M_RHS_send_inds

        call CSR_Dagger(M, partition_table, M_T, MPI_communicator)

        call Reconcile_Communications(  M_T, &
                                        partition_table, &
                                        MPI_communicator)

        itmax = M%columns/l

        one_norm_array = 0

        call One_Norm(  M, &
                        one_norm_array(1), &
                        partition_table, &
                        MPI_communicator)

        one_norm_array(1) = one_norm_array(1)

        p = pmax

        do i = 2, pmax + 1

            call One_Norm_Estimation(   M, &
                                        M_T, &
                                        i, &
                                        l, &
                                        itmax, &
                                        partition_table, &
                                        one_norm_array(i), &
                                        mpi_communicator)

            alphas(i) = one_norm_array(i)**(1_dp/real(i,8))

            if (i >= 3) then
                if((abs((alphas(i - 1) - alphas(i))/alphas(i)) < 0.5)) then
                    p = i - 1
                    exit
                endif
            endif

        enddo

        deallocate(M_T%values, M_T%row_starts, M_T%local_col_inds, M_T%col_indexes, &
            M_T%num_rec_inds, M_T%rec_disps, M_T%num_send_inds, M_T%send_disps, M_T%RHS_send_inds)

end subroutine one_norm_series

subroutine initial_state(   rho0_rows, &
                            M_local_rows, &
                            rho0, &
                            flock, &
                            rank, &
                            partition_table, &
                            MPI_communicator, &
                            rho0_v)

    use :: Operators

    integer, intent(in) :: rho0_rows
    integer, intent(in) :: M_local_rows
    complex(8), dimension(rho0_rows, rho0_rows), intent(in) :: rho0
    integer, intent(in) :: flock
    integer, intent(in) :: rank
    integer, dimension(flock + 1), intent(in) :: partition_table
    integer, intent(in) :: MPI_communicator
    complex(8), dimension(M_local_rows), &
        intent(out) :: rho0_v

    complex(8), dimension(:), allocatable :: rho0_v_temp

    integer :: i

!f2py integer, optional,intent(in),check(shape(rho0,0)==rho0_rows),depend(rho0) :: rho0_rows=shape(rho0,0)
!f2py integer intent(in) :: m_local_rows
!f2py complex(kind=8) dimension(rho0_rows,rho0_rows),intent(in) :: rho0
!f2py integer, optional,intent(in),check((len(partition_table)-1)>=flock),depend(partition_table) :: flock=(len(partition_table)-1)
!f2py integer intent(in) :: rank_bn
!f2py integer dimension(flock + 1),intent(in) :: partition_table
!f2py integer intent(in) :: mpi_communicator
!f2py complex(kind=8) dimension(m_local_rows),intent(out),depend(m_local_rows) :: rho0_v


    call Vectorize_Operator(rho0, &
                            partition_table, &
                            rho0_v_temp, &
                            MPI_communicator)

    !$omp parallel do
    do i = 1, M_local_rows
        rho0_v(i) = rho0_v_temp(i + partition_table(rank + 1) - 1)
    enddo
    !$omp end parallel do

end subroutine initial_state

subroutine step(M_rows, &
                M_n_col_indexes, &
                M_n_values, &
                M_n_local_col_indexes, &
                n_rho0_v, &
                n_rhot_v, &
                M_sends, &
                M_row_starts, &
                M_col_indexes, &
                M_values, &
                M_num_rec_inds, &
                M_rec_disps, &
                M_num_send_inds, &
                M_send_disps, &
                M_local_col_inds, &
                M_RHS_send_inds, &
                mu, &
                t, &
                rho0_v, &
                flock, &
                partition_table, &
                p, &
                one_norm_array, &
                MPI_communicator, &
                rhot_v, &
                target_precision, &
                M_n_row_starts)

        use :: Sparse
        use :: Expm

        implicit none

        integer, intent(in) :: M_rows
        integer, intent(in) :: M_n_row_starts
        integer, intent(in) :: M_n_col_indexes
        integer, intent(in) :: M_n_values
        integer, intent(in) :: M_n_local_col_indexes
        integer, intent(in) :: n_rho0_v
        integer, intent(in) :: n_rhot_v
        integer, intent(in) :: M_sends
        integer, dimension(M_n_row_starts), target, intent(in) :: M_row_starts
        integer, dimension(M_n_col_indexes), target, intent(in) :: M_col_indexes
        complex(8), dimension(M_n_values), target, intent(in) :: M_values
        integer, dimension(flock), target, intent(in) :: M_num_rec_inds
        integer, dimension(flock), target, intent(in) :: M_rec_disps
        integer, dimension(flock), target, intent(in) :: M_num_send_inds
        integer, dimension(flock), target, intent(in) :: M_send_disps
        integer, dimension(M_n_local_col_indexes), target, intent(in) :: M_local_col_inds
        integer, dimension(M_sends), target, intent(in) :: M_RHS_send_inds
        complex(8), intent(inout) :: mu
        real(8), intent(in) :: t
        complex(8), dimension(n_rho0_v), intent(in) :: rho0_v
        integer, intent(in) :: flock
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(inout) :: p
        real(8), dimension(9), intent(inout) :: one_norm_array
        integer, intent(in) :: MPI_communicator
        complex(8), dimension(n_rhot_v), intent(out) :: rhot_v
        character(len=2), intent(in) :: target_precision

        type(CSR) :: M

        integer :: rank, ierr

        integer :: lb, ub
        integer :: lb_elements, ub_elements

!f2py integer intent(in) :: m_rows
!f2py integer, optional,intent(in),check(len(m_col_indexes)>=m_n_col_indexes),depend(m_col_indexes) :: m_n_col_indexes=len(m_col_indexes)
!f2py integer, optional,intent(in),check(len(m_values)>=m_n_values),depend(m_values) :: m_n_values=len(m_values)
!f2py integer, optional,intent(in),check(len(m_local_col_inds)>=m_n_local_col_indexes),depend(m_local_col_inds) :: m_n_local_col_indexes=len(m_local_col_inds)
!f2py integer, optional,intent(in),check(len(rho0_v)>=n_rho0_v),depend(rho0_v) :: n_rho0_v=len(rho0_v)
!f2py integer, optional,intent(in),depend(rho0_v) :: n_rhot_v=len(rho0_v)
!f2py integer, optional,intent(in),check(len(m_rhs_send_inds)>=m_sends),depend(m_rhs_send_inds) :: m_sends=len(m_rhs_send_inds)
!f2py integer, target,dimension(m_n_row_starts),intent(in) :: m_row_starts
!f2py integer, target,dimension(m_n_col_indexes),intent(in) :: m_col_indexes
!f2py complex(kind=8), target,dimension(m_n_values),intent(in),depend(m_n_values) :: m_values
!f2py integer, target,dimension(flock),intent(in) :: m_num_rec_inds
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_rec_disps
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_num_send_inds
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_send_disps
!f2py integer, target,dimension(m_n_local_col_indexes),intent(in) :: m_local_col_inds
!f2py integer, target,dimension(m_sends),intent(in) :: m_rhs_send_inds
!f2py complex(kind=8) intent(in) :: mu
!f2py real(kind=8) intent(in) :: t
!f2py complex(kind=8) dimension(n_rho0_v),intent(in) :: rho0_v
!f2py integer, optional,intent(in),check(len(m_num_rec_inds)>=flock),depend(m_num_rec_inds) :: flock=len(m_num_rec_inds)
!f2py integer dimension(flock + 1),intent(in),depend(flock) :: partition_table
!f2py integer intent(inout) :: p
!f2py real(kind=8) dimension(9),intent(inout) :: one_norm_array
!f2py integer intent(in) :: mpi_communicator
!f2py complex(kind=8) dimension(n_rhot_v),intent(out),depend(n_rhot_v) :: rhot_v
!f2py character*2 intent(in) :: target_precision
!f2py integer, optional,intent(in),check(len(m_row_starts)>=m_n_row_starts),depend(m_row_starts) :: m_n_row_starts=len(m_row_starts)

        call mpi_comm_rank(mpi_communicator, rank, ierr)

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        m%rows = m_rows
        m%columns = m_rows
        m%row_starts(lb:ub + 1) => m_row_starts

        lb_elements = m%row_starts(lb)
        ub_elements = m%row_starts(ub + 1) - 1

        M%col_indexes(lb_elements:ub_elements) => M_col_indexes
        M%values(lb_elements:ub_elements) => M_values
        M%local_col_inds(lb_elements:ub_elements) => M_local_col_inds
        M%num_rec_inds(1:size(M_num_rec_inds)) => M_num_rec_inds
        M%rec_disps(1:size(M_rec_disps)) => M_rec_disps
        M%num_send_inds(1:size(M_num_send_inds)) => M_num_send_inds
        M%send_disps(1:size(M_send_disps)) => M_send_disps
        M%RHS_send_inds(1:size(M_RHS_send_inds)) => M_RHS_send_inds


        call Expm_Multiply( M, &
                            rho0_v, &
                            t, &
                            partition_table, &
                            rhot_v, &
                            MPI_communicator, &
                            one_norm_series = one_norm_array, &
                            p = p, &
                            target_precision = target_precision, &
                            mu = mu)

end subroutine step

subroutine gather_step( M_local_rows, &
                        rhot_v, &
                        flock, &
                        partition_table, &
                        root, &
                        MPI_communicator, &
                        H_aug_rows, &
                        rhot)

    use :: Sparse
    use :: Operators

    implicit none

    integer, intent(in) :: M_local_rows
    complex(8), dimension(M_local_rows), intent(in) :: rhot_v
    integer, intent(in) :: flock
    integer, dimension(flock + 1), intent(in) :: partition_table
    integer, intent(in) :: root
    integer, intent(in) :: MPI_communicator
    integer, intent(in) :: H_aug_rows
    complex(8), dimension(H_aug_rows, H_aug_rows), intent(out) :: rhot

    complex(8), dimension(:), allocatable :: rhot_gathered

!f2py integer, optional,intent(in),check(len(rhot_v)>=m_local_rows),depend(rhot_v) :: m_local_rows=len(rhot_v)
!f2py complex(kind=8) dimension(m_local_rows),intent(in) :: rhot_v
!f2py integer, optional,intent(in),check((len(partition_table)-1)>=flock),depend(partition_table) :: flock=(len(partition_table)-1)
!f2py integer dimension(flock + 1),intent(in) :: partition_table
!f2py integer intent(in) :: root
!f2py integer intent(in) :: mpi_communicator
!f2py integer intent(in) :: h_aug_rows
!f2py complex(kind=8) dimension(h_aug_rows,h_aug_rows),intent(out),depend(h_aug_rows,h_aug_rows) :: rhot

    allocate(rhot_gathered(partition_table(1):partition_table(size(partition_table)) - 1))


        call Gather_Dense_Vector(   rhot_v, &
                                    partition_table, &
                                    root, &
                                    rhot_gathered, &
                                    MPI_communicator)

        call Reshape_Vectorized_Operator(rhot_gathered, rhot)

end subroutine gather_step

subroutine series(  M_rows, &
                    M_n_row_starts,&
                    M_n_col_indexes, &
                    M_n_values, &
                    M_n_local_col_indexes, &
                    n_rho0_v, &
                    n_rhot_v, &
                    M_sends, &
                    M_row_starts, &
                    M_col_indexes, &
                    M_values, &
                    M_num_rec_inds, &
                    M_rec_disps, &
                    M_num_send_inds, &
                    M_send_disps, &
                    M_local_col_inds, &
                    M_RHS_send_inds, &
                    mu, &
                    t_1, &
                    t_2, &
                    rho0_v, &
                    steps, &
                    flock, &
                    partition_table, &
                    p, &
                    one_norm_array, &
                    MPI_communicator, &
                    rhot_v_series, &
                    target_precision)

        use :: Sparse
        use :: Expm

        implicit none

        integer, intent(in) :: M_rows
        integer, intent(in) :: M_n_row_starts
        integer, intent(in) :: M_n_col_indexes
        integer, intent(in) :: M_n_values
        integer, intent(in) :: M_n_local_col_indexes
        integer, intent(in) :: n_rho0_v
        integer, intent(in) :: n_rhot_v
        integer, intent(in) :: M_sends
        integer, dimension(M_n_row_starts), target, intent(in) :: M_row_starts
        integer, dimension(M_n_col_indexes), target, intent(in) :: M_col_indexes
        complex(8), dimension(M_n_values), target, intent(in) :: M_values
        integer, dimension(flock), target, intent(in) :: M_num_rec_inds
        integer, dimension(flock), target, intent(in) :: M_rec_disps
        integer, dimension(flock), target, intent(in) :: M_num_send_inds
        integer, dimension(flock), target, intent(in) :: M_send_disps
        integer, dimension(M_n_local_col_indexes), target, intent(in) :: M_local_col_inds
        integer, dimension(M_sends), target, intent(in) :: M_RHS_send_inds
        complex(8), intent(inout) :: mu
        real(8), intent(in) :: t_1
        real(8), intent(in) :: t_2
        complex(8), dimension(n_rho0_v), intent(in) :: rho0_v
        integer, intent(in) :: steps
        integer, intent(in) :: flock
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(inout) :: p
        real(8), dimension(9), intent(inout) :: one_norm_array
        integer, intent(in) :: MPI_communicator
        complex(8), dimension(n_rhot_v, steps + 1), intent(out) :: rhot_v_series
        character(len = 2), intent(in) :: target_precision

        type(CSR) :: M


        integer :: rank, ierr

        integer :: lb, ub
        integer :: lb_elements, ub_elements

!f2py integer intent(in) :: m_rows
!f2py integer, optional,intent(in),check(len(m_row_starts)>=m_n_row_starts),depend(m_row_starts) :: m_n_row_starts=len(m_row_starts)
!f2py integer, optional,intent(in),check(len(m_col_indexes)>=m_n_col_indexes),depend(m_col_indexes) :: m_n_col_indexes=len(m_col_indexes)
!f2py integer, optional,intent(in),check(len(m_values)>=m_n_values),depend(m_values) :: m_n_values=len(m_values)
!f2py integer, optional,intent(in),check(len(m_local_col_inds)>=m_n_local_col_indexes),depend(m_local_col_inds) :: m_n_local_col_indexes=len(m_local_col_inds)
!f2py integer, optional,intent(in),check(len(rho0_v)>=n_rho0_v),depend(rho0_v) :: n_rho0_v=len(rho0_v)
!f2py integer, optional, intent(in), depend(rhot_v) :: n_rhot_v=len(rho0_v)
!f2py integer, optional,intent(in),check(len(m_rhs_send_inds)>=m_sends),depend(m_rhs_send_inds) :: m_sends=len(m_rhs_send_inds)
!f2py integer, target,dimension(m_n_row_starts),intent(in) :: m_row_starts
!f2py integer, target,dimension(m_n_col_indexes),intent(in) :: m_col_indexes
!f2py complex(kind=8), target,dimension(m_n_values),intent(in) :: m_values
!f2py integer, target,dimension(flock),intent(in) :: m_num_rec_inds
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_rec_disps
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_num_send_inds
!f2py integer, target,dimension(flock),intent(in),depend(flock) :: m_send_disps
!f2py integer, target,dimension(m_n_local_col_indexes),intent(in) :: m_local_col_inds
!f2py integer, target,dimension(m_sends),intent(in) :: m_rhs_send_inds
!f2py real(kind=8) intent(in) :: t_1
!f2py real(kind=8) intent(in) :: t_2
!f2py complex(kind=8) dimension(n_rho0_v),intent(in) :: rho0_v
!f2py integer intent(in) :: steps
!f2py integer, optional,intent(in),check(len(m_num_rec_inds)>=flock),depend(m_num_rec_inds) :: flock=len(m_num_rec_inds)
!f2py integer dimension(flock + 1),intent(in),depend(flock) :: partition_table
!f2py integer intent(inout) :: p
!f2py real(kind=8) dimension(9),intent(inout) :: one_norm_array
!f2py integer intent(in) :: mpi_communicator
!f2py complex(kind=8) dimension(n_rhot_v,steps + 1),intent(out),depend(n_rhot_v,steps) :: rhot_v_series
!f2py character*2 intent(inout) :: target_precision

        call mpi_comm_rank(mpi_communicator, rank, ierr)

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        m%rows = m_rows
        m%columns = m_rows
        m%row_starts(lb:ub + 1) => m_row_starts

        lb_elements = m%row_starts(lb)
        ub_elements = m%row_starts(ub + 1) - 1

        M%col_indexes(lb_elements:ub_elements) => M_col_indexes
        M%values(lb_elements:ub_elements) => M_values
        M%local_col_inds(lb_elements:ub_elements) => M_local_col_inds
        M%num_rec_inds(1:size(M_num_rec_inds)) => M_num_rec_inds
        M%rec_disps(1:size(M_rec_disps)) => M_rec_disps
        M%num_send_inds(1:size(M_num_send_inds)) => M_num_send_inds
        M%send_disps(1:size(M_send_disps)) => M_send_disps
        M%RHS_send_inds(1:size(M_RHS_send_inds)) => M_RHS_send_inds

        call Expm_Multiply_Series(  M, &
                                    rho0_v, &
                                    t_1, &
                                    t_2, &
                                    steps, &
                                    partition_table, &
                                    rhot_v_series, &
                                    MPI_communicator, &
                                    one_norm_series_in = one_norm_array, &
                                    p_ex = p, &
                                    target_precision = target_precision, &
                                    mu = mu)

end subroutine series

subroutine gather_series(   M_local_rows, &
                            steps, &
                            rhot_v_series, &
                            flock, &
                            partition_table, &
                            root, &
                            MPI_communicator, &
                            H_aug_rows, &
                            rhot_series)

    use :: Sparse
    use :: Operators

    implicit none

    integer, intent(in) :: M_local_rows
    integer, intent(in) :: steps
    complex(8), dimension(M_local_rows, steps + 1), intent(in) :: rhot_v_series
    integer, intent(in) :: flock
    integer, dimension(flock + 1), intent(in) :: partition_table
    integer, intent(in) :: root
    integer, intent(in) :: MPI_communicator
    integer, intent(in) :: H_aug_rows
    complex(8), dimension(H_aug_rows, H_aug_rows, steps + 1), intent(out) :: rhot_series

    complex(8), dimension(:,:), allocatable :: rhot_v_series_gathered

!f2py integer, optional,intent(in),check(shape(rhot_v_series,0)==m_local_rows),depend(rhot_v_series) :: m_local_rows=shape(rhot_v_series,0)
!f2py integer, optional,intent(in),check((shape(rhot_v_series,1)-1)==steps),depend(rhot_v_series) :: steps=(shape(rhot_v_series,1)-1)
!f2py complex(kind=8) dimension(m_local_rows,steps + 1),intent(in) :: rhot_v_series
!f2py integer, optional,intent(in),check((len(partition_table)-1)>=flock),depend(partition_table) :: flock=(len(partition_table)-1)
!f2py integer dimension(flock + 1),intent(in) :: partition_table
!f2py integer intent(in) :: root
!f2py integer intent(in) :: mpi_communicator
!f2py integer intent(in) :: h_aug_rows
!f2py complex(kind=8) dimension(h_aug_rows,h_aug_rows,steps + 1),intent(out),depend(h_aug_rows,h_aug_rows,steps) :: rhot_series

    allocate(rhot_v_series_gathered(partition_table(1): &
        partition_table(size(partition_table)) - 1, steps + 1))

    call Gather_Dense_Matrix(   rhot_v_series, &
                                partition_table, &
                                root, &
                                rhot_v_series_gathered, &
                                MPI_communicator)

    call Reshape_Vectorized_Operator_Series( rhot_v_series_gathered, rhot_series)

end subroutine gather_series