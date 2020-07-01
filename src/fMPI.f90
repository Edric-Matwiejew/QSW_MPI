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

subroutine Generalized_Super_Operator_Extent(   taus, &
                                                H_rows, &
                                                H_nnz, &
                                                H_row_starts, &
                                                H_col_indexes, &
                                                H_values, &
                                                L_n, &
                                                L_nnz, &
                                                L_row_total, &
                                                L_nnz_total, &
                                                L_row_starts, &
                                                L_col_indexes, &
                                                L_values, &
                                                flock, &
                                                MPI_communicator, &
                                                M_nnz_out, &
                                                M_rows, &
                                                partition_table)

    use :: iso_precisions
    use :: sparse
    use :: operators
    use :: mpi

    implicit none

    integer, intent(in) :: H_rows
    integer, intent(in) :: H_nnz
    integer, dimension(H_rows + 1), intent(in) :: H_row_starts
    integer, dimension(H_nnz), intent(in) :: H_col_indexes
    complex(dp), dimension(H_nnz), intent(in) :: H_values
    integer, intent(in) :: L_n
    real(dp), dimension(L_n), intent(in) :: taus
    integer, dimension(L_n), intent(in) :: L_nnz
    integer, intent(in) :: L_row_total
    integer, intent(in) :: L_nnz_total
    integer, dimension(L_row_total), intent(in) :: L_row_starts
    integer, dimension(L_nnz_total), intent(in) :: L_col_indexes
    complex(dp), dimension(L_nnz_total), intent(in) :: L_values
    integer, intent(in) :: flock
    integer, intent(in) :: MPI_communicator
    integer, intent(out) :: M_nnz_out
    integer, intent(out) :: M_rows
    integer, dimension(flock + 1), intent(out) :: partition_table

    integer, dimension(:), allocatable :: partition_table_temp
    type(CSR), dimension(L_n) :: Ls
    type(CSR) :: H
    type(CSR) :: M

    integer :: i

    allocate(H%row_starts(H_rows + 1))
    allocate(H%col_indexes(H_nnz))
    allocate(H%values(H_nnz))

    H%rows = H_rows
    H%columns = H_rows
    H%row_starts = H_row_starts + 1
    H%col_indexes = H_col_indexes + 1
    H%values = H_values

    do i = 1, L_n
        Ls(i)%rows = H_rows
        Ls(i)%columns = H_rows
        allocate(Ls(i)%row_starts(H_rows + 1))
        allocate(Ls(i)%col_indexes(L_nnz(i)))
        allocate(Ls(i)%values(L_nnz(i)))
        Ls(i)%row_starts = L_row_starts((i - 1)*(H_rows + 1) + 1: i*(H_rows + 1)) + 1
        Ls(i)%col_indexes = L_col_indexes(sum(L_nnz(1:i - 1)) + 1:sum(L_nnz(1:i))) + 1
        Ls(i)%values = L_values(sum(L_nnz(1:i - 1)) + 1:sum(L_nnz(1:i)))
    enddo

    call generate_partition_table(H_rows**2, partition_table_temp, MPI_communicator)

    call generalized_superoperator( taus, &
                                    H, &
                                    Ls, &
                                    partition_table_temp, &
                                    M, &
                                    MPI_communicator)

    partition_table = partition_table_temp
    M_nnz_out = size(M%col_indexes)
    M_rows = M%columns

end subroutine Generalized_Super_Operator_Extent

subroutine Generalized_Super_Operator(  taus, &
                                        H_rows, &
                                        H_nnz, &
                                        H_row_starts, &
                                        H_col_indexes, &
                                        H_values, &
                                        L_n, &
                                        L_nnz, &
                                        L_row_total, &
                                        L_nnz_total, &
                                        L_row_starts, &
                                        L_col_indexes, &
                                        L_values, &
                                        flock, &
                                        partition_table, &
                                        MPI_communicator, &
                                        M_nnz_in, &
                                        M_n_row_starts, &
                                        M_row_starts, &
                                        M_col_indexes, &
                                        M_values)

    use :: iso_precisions
    use :: sparse
    use :: operators
    use :: mpi

    implicit none

    integer, intent(in) :: H_rows
    integer, intent(in) :: H_nnz
    integer, dimension(H_rows + 1), intent(in) :: H_row_starts
    integer, dimension(H_nnz), intent(in) :: H_col_indexes
    complex(dp), dimension(H_nnz), intent(in) :: H_values
    integer, intent(in) :: L_n
    real(dp), dimension(L_n), intent(in) :: taus
    integer, dimension(L_n), intent(in) :: L_nnz
    integer, intent(in) :: L_row_total
    integer, intent(in) :: L_nnz_total
    integer, dimension(L_row_total), intent(in) :: L_row_starts
    integer, dimension(L_nnz_total), intent(in) :: L_col_indexes
    complex(dp), dimension(L_nnz_total), intent(in) :: L_values
    integer, intent(in) :: flock
    integer, dimension(flock + 1), intent(in) :: partition_table
    integer, intent(in) :: MPI_communicator
    integer, intent(in) :: M_nnz_in
    integer, intent(in) :: M_n_row_starts
    integer, dimension(M_n_row_starts), intent(out) :: M_row_starts
    integer, dimension(M_nnz_in), intent(out) :: M_col_indexes
    complex(dp), dimension(M_nnz_in), intent(out) :: M_values

    type(CSR), dimension(L_n) :: Ls
    type(CSR) :: H
    type(CSR) :: M

    integer :: i

    allocate(H%row_starts(H_rows + 1))
    allocate(H%col_indexes(H_nnz))
    allocate(H%values(H_nnz))

    H%rows = H_rows
    H%columns = H_rows
    H%row_starts = H_row_starts + 1
    H%col_indexes = H_col_indexes + 1
    H%values = H_values

    do i = 1, L_n
        Ls(i)%rows = H_rows
        Ls(i)%columns = H_rows
        allocate(Ls(i)%row_starts(H_rows + 1))
        allocate(Ls(i)%col_indexes(L_nnz(i)))
        allocate(Ls(i)%values(L_nnz(i)))
        Ls(i)%row_starts = L_row_starts((i - 1)*(H_rows + 1) + 1: i*(H_rows + 1)) + 1
        Ls(i)%col_indexes = L_col_indexes(sum(L_nnz(1:i - 1)) + 1:sum(L_nnz(1:i))) + 1
        Ls(i)%values = L_values(sum(L_nnz(1:i - 1)) + 1:sum(L_nnz(1:i)))
    enddo

    call generalized_superoperator( taus, &
                                    H, &
                                    Ls, &
                                    partition_table, &
                                    M, &
                                    MPI_communicator)

    M_row_starts(1:size(M%row_starts)) = M%row_starts
    M_col_indexes(1:size(M%col_indexes)) = M%col_indexes
    M_values(1:size(M%col_indexes)) = M%values

end subroutine Generalized_Super_Operator

subroutine Demoralized_Super_Operator_Extent(   omega, &
                                                H_rows, &
                                                H_nnz, &
                                                H_row_starts, &
                                                H_col_indexes, &
                                                H_values, &
                                                H_loc_nnz, &
                                                H_loc_row_starts, &
                                                H_loc_col_indexes, &
                                                H_loc_values, &
                                                L_n, &
                                                L_nnz, &
                                                L_row_total, &
                                                L_nnz_total, &
                                                L_row_starts, &
                                                L_col_indexes, &
                                                L_values, &
                                                flock, &
                                                MPI_communicator, &
                                                M_nnz_out, &
                                                M_rows, &
                                                partition_table)

    use :: iso_precisions
    use :: sparse
    use :: operators
    use :: mpi

    implicit none

    real(dp), intent(in) :: omega
    integer, intent(in) :: H_rows
    integer, intent(in) :: H_nnz
    integer, dimension(H_rows + 1), intent(in) :: H_row_starts
    integer, dimension(H_nnz), intent(in) :: H_col_indexes
    complex(dp), dimension(H_nnz), intent(in) :: H_values
    integer, intent(in) :: H_loc_nnz
    integer, dimension(H_rows + 1), intent(in) :: H_loc_row_starts
    integer, dimension(H_loc_nnz), intent(in) :: H_loc_col_indexes
    complex(dp), dimension(H_loc_nnz), intent(in) :: H_loc_values
    integer, intent(in) :: L_n
    integer, dimension(L_n), intent(in) :: L_nnz
    integer, intent(in) :: L_row_total
    integer, intent(in) :: L_nnz_total
    integer, dimension(L_row_total), intent(in) :: L_row_starts
    integer, dimension(L_nnz_total), intent(in) :: L_col_indexes
    complex(dp), dimension(L_nnz_total), intent(in) :: L_values
    integer, intent(in) :: flock
    integer, intent(in) :: MPI_communicator
    integer, intent(out) :: M_nnz_out
    integer, intent(out) :: M_rows
    integer, dimension(flock + 1), intent(out) :: partition_table

    integer, dimension(:), allocatable :: partition_table_temp
    type(CSR), dimension(L_n) :: Ls
    type(CSR) :: H
    type(CSR) :: H_loc
    type(CSR) :: M

    integer :: i

    allocate(H%row_starts(H_rows + 1))
    allocate(H%col_indexes(H_nnz))
    allocate(H%values(H_nnz))

    H%rows = H_rows
    H%columns = H_rows
    H%row_starts = H_row_starts + 1
    H%col_indexes = H_col_indexes + 1
    H%values = H_values


    if (H_loc_row_starts(H_rows + 1) == 0) then

        allocate(H_loc%row_starts(H_rows + 1))
        allocate(H_loc%col_indexes(0))
        allocate(H_loc%values(0))

        H_loc%rows = H_rows
        H_loc%columns = H_rows
        H_loc%row_starts = H_loc_row_starts + 1

    else

        allocate(H_loc%row_starts(H_rows + 1))
        allocate(H_loc%col_indexes(H_loc_nnz))
        allocate(H_loc%values(H_loc_nnz))

        H_loc%rows = H_rows
        H_loc%columns = H_rows
        H_loc%row_starts = H_loc_row_starts + 1
        H_loc%col_indexes = H_loc_col_indexes + 1
        H_loc%values = H_loc_values

    endif

    do i = 1, L_n
        Ls(i)%rows = H_rows
        Ls(i)%columns = H_rows
        allocate(Ls(i)%row_starts(H_rows + 1))
        allocate(Ls(i)%col_indexes(L_nnz(i)))
        allocate(Ls(i)%values(L_nnz(i)))
        Ls(i)%row_starts = L_row_starts((i - 1)*(H_rows + 1) + 1: i*(H_rows + 1)) + 1
        Ls(i)%col_indexes = L_col_indexes(sum(L_nnz(1:i - 1)) + 1:sum(L_nnz(1:i))) + 1
        Ls(i)%values = L_values(sum(L_nnz(1:i - 1)) + 1:sum(L_nnz(1:i)))
    enddo

    call generate_partition_table(H_rows**2, partition_table_temp, MPI_communicator)

    call Demoralized_superoperator( omega, &
                                    H, &
                                    H_loc, &
                                    Ls, &
                                    partition_table_temp, &
                                    M, &
                                    MPI_communicator)

    partition_table = partition_table_temp
    M_nnz_out = size(M%col_indexes)
    M_rows = M%columns

end subroutine Demoralized_Super_Operator_Extent

subroutine Demoralized_Super_Operator(  omega, &
                                        H_rows, &
                                        H_nnz, &
                                        H_row_starts, &
                                        H_col_indexes, &
                                        H_values, &
                                        H_loc_nnz, &
                                        H_loc_row_starts, &
                                        H_loc_col_indexes, &
                                        H_loc_values, &
                                        L_n, &
                                        L_nnz, &
                                        L_row_total, &
                                        L_nnz_total, &
                                        L_row_starts, &
                                        L_col_indexes, &
                                        L_values, &
                                        flock, &
                                        partition_table, &
                                        MPI_communicator, &
                                        M_nnz_in, &
                                        M_n_row_starts, &
                                        M_row_starts, &
                                        M_col_indexes, &
                                        M_values)

    use :: iso_precisions
    use :: sparse
    use :: operators
    use :: mpi

    implicit none

    real(dp), intent(in) :: omega
    integer, intent(in) :: H_rows
    integer, intent(in) :: H_nnz
    integer, dimension(H_rows + 1), intent(in) :: H_row_starts
    integer, dimension(H_nnz), intent(in) :: H_col_indexes
    complex(dp), dimension(H_nnz), intent(in) :: H_values
    integer, intent(in) :: H_loc_nnz
    integer, dimension(H_rows + 1), intent(in) :: H_loc_row_starts
    integer, dimension(H_loc_nnz), intent(in) :: H_loc_col_indexes
    complex(dp), dimension(H_loc_nnz), intent(in) :: H_loc_values
    integer, intent(in) :: L_n
    integer, dimension(L_n), intent(in) :: L_nnz
    integer, intent(in) :: L_row_total
    integer, intent(in) :: L_nnz_total
    integer, dimension(L_row_total), intent(in) :: L_row_starts
    integer, dimension(L_nnz_total), intent(in) :: L_col_indexes
    complex(dp), dimension(L_nnz_total), intent(in) :: L_values
    integer, intent(in) :: flock
    integer, dimension(flock + 1), intent(in) :: partition_table
    integer, intent(in) :: MPI_communicator
    integer, intent(in) :: M_nnz_in
    integer, intent(in) :: M_n_row_starts
    integer, dimension(M_n_row_starts), intent(out) :: M_row_starts
    integer, dimension(M_nnz_in), intent(out) :: M_col_indexes
    complex(dp), dimension(M_nnz_in), intent(out) :: M_values

    type(CSR), dimension(L_n) :: Ls
    type(CSR) :: H
    type(CSR) :: H_loc
    type(CSR) :: M

    integer :: i

    allocate(H%row_starts(H_rows + 1))
    allocate(H%col_indexes(H_nnz))
    allocate(H%values(H_nnz))

    H%rows = H_rows
    H%columns = H_rows
    H%row_starts = H_row_starts + 1
    H%col_indexes = H_col_indexes + 1
    H%values = H_values

    if (H_loc_row_starts(H_rows + 1) == 0) then

        allocate(H_loc%row_starts(H_rows + 1))
        allocate(H_loc%col_indexes(0))
        allocate(H_loc%values(0))

        H_loc%rows = H_rows
        H_loc%columns = H_rows
        H_loc%row_starts = H_loc_row_starts + 1

    else

        allocate(H_loc%row_starts(H_rows + 1))
        allocate(H_loc%col_indexes(H_loc_nnz))
        allocate(H_loc%values(H_loc_nnz))

        H_loc%rows = H_rows
        H_loc%columns = H_rows
        H_loc%row_starts = H_loc_row_starts + 1
        H_loc%col_indexes = H_loc_col_indexes + 1
        H_loc%values = H_loc_values

    endif

    do i = 1, L_n
        Ls(i)%rows = H_rows
        Ls(i)%columns = H_rows
        allocate(Ls(i)%row_starts(H_rows + 1))
        allocate(Ls(i)%col_indexes(L_nnz(i)))
        allocate(Ls(i)%values(L_nnz(i)))
        Ls(i)%row_starts = L_row_starts((i - 1)*(H_rows + 1) + 1: i*(H_rows + 1)) + 1
        Ls(i)%col_indexes = L_col_indexes(sum(L_nnz(1:i - 1)) + 1:sum(L_nnz(1:i))) + 1
        Ls(i)%values = L_values(sum(L_nnz(1:i - 1)) + 1:sum(L_nnz(1:i)))
    enddo

    call Demoralized_superoperator( omega, &
                                    H, &
                                    H_loc, &
                                    Ls, &
                                    partition_table, &
                                    M, &
                                    MPI_communicator)

    M_row_starts(1:size(M%row_starts)) = M%row_starts
    M_col_indexes(1:size(M%col_indexes)) = M%col_indexes
    M_values(1:size(M%col_indexes)) = M%values

end subroutine Demoralized_Super_Operator

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

        use :: iso_precisions
        use :: Operators
        use :: MPI

        implicit none

        real(dp), intent(in) :: omega
        integer, intent(in) :: H_rows
        integer, intent(in) :: H_nnz
        integer, dimension(H_rows + 1), intent(in) :: H_row_starts
        integer, dimension(H_nnz), intent(in) :: H_col_indexes
        complex(dp), dimension(H_nnz), intent(in) :: H_values
        integer, intent(in) :: L_nnz
        integer, dimension(H_rows + 1), intent(in) :: L_row_starts
        integer, dimension(L_nnz), intent(in) :: L_col_indexes
        complex(dp), dimension(L_nnz), intent(in) :: L_values
        integer, intent(in) :: n_sources
        integer, dimension(n_sources), intent(in) :: source_sites
        real(dp), dimension(n_sources), target, intent(in) :: source_rates
        integer, intent(in) :: n_sinks
        integer, dimension(n_sinks), intent(in) :: sink_sites
        real(dp), dimension(n_sinks), target, intent(in) :: sink_rates
        integer, intent(in) :: flock
        integer, intent(in) :: MPI_communicator
        integer, intent(out) :: M_nnz_out
        integer, intent(out) :: M_rows
        integer, dimension(flock + 1), intent(out) :: partition_table

        type(CSR) :: H, L, M

        integer, dimension(:), allocatable :: source_sites_temp
        real(dp), dimension(:), pointer :: source_rates_temp
        integer, dimension(:), allocatable :: sink_sites_temp
        real(dp), dimension(:), pointer :: sink_rates_temp
        integer, dimension(:), allocatable :: partition_table_temp

        if (source_sites(1) < 0) then
            allocate(source_sites_temp(0))
            allocate(source_rates_temp(0))
        else

            allocate(source_sites_temp(n_sources))

            source_sites_temp = source_sites + 1
            source_rates_temp => source_rates

        endif

        if (sink_sites(1) < 0) then
            allocate(sink_sites_temp(0))
            allocate(sink_rates_temp(0))
        else

            allocate(sink_sites_temp(n_sinks))

            sink_sites_temp = sink_sites + 1
            sink_rates_temp => sink_rates

        endif

        allocate(H%row_starts(H_rows + 1))
        allocate(H%col_indexes(H_nnz))
        allocate(H%values(H_nnz))

        H%rows = H_rows
        H%columns = H_rows

        H%row_starts = H_row_starts + 1
        H%col_indexes = H_col_indexes + 1
        H%values = H_values

        allocate(L%row_starts(H_rows + 1))
        allocate(L%col_indexes(L_nnz))
        allocate(L%values(L_nnz))

        L%rows = H_rows
        L%columns = H_rows

        L%row_starts = L_row_starts + 1
        L%col_indexes = L_col_indexes + 1
        L%values = L_values

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
                            MPI_communicator, &
                            M_row_starts, &
                            M_col_indexes, &
                            M_values)

        use :: iso_precisions
        use :: Operators
        use :: Expm
        use :: MPI

        implicit none

        real(dp), intent(in) :: omega
        integer, intent(in) :: H_rows
        integer, intent(in) :: H_nnz
        integer, dimension(H_rows + 1), intent(in) :: H_row_starts
        integer, dimension(H_nnz), intent(in) :: H_col_indexes
        complex(dp), dimension(H_nnz), intent(in) :: H_values
        integer, intent(in) :: L_nnz
        integer, dimension(H_rows + 1), intent(in) :: L_row_starts
        integer, dimension(L_nnz), intent(in) :: L_col_indexes
        complex(dp), dimension(L_nnz), intent(in) :: L_values
        integer, intent(in) :: n_sources
        integer, dimension(n_sources), intent(in) :: source_sites
        real(dp), dimension(n_sources), target, intent(in) :: source_rates
        integer, intent(in) :: n_sinks
        integer, dimension(n_sinks), intent(in) :: sink_sites
        real(dp), dimension(n_sinks), target, intent(in) :: sink_rates
        integer, intent(in) :: M_n_row_starts
        integer, intent(in) :: M_nnz_in
        integer, intent(in) :: MPI_communicator
        integer, dimension(M_n_row_starts), intent(out) :: M_row_starts
        integer, dimension(M_nnz_in), intent(out) :: M_col_indexes
        complex(dp), dimension(M_nnz_in), intent(out) :: M_values

        type(CSR) :: H, L, M

        integer, dimension(:), allocatable :: source_sites_temp
        real(dp), dimension(:), pointer :: source_rates_temp
        integer, dimension(:), allocatable :: sink_sites_temp
        real(dp), dimension(:), pointer :: sink_rates_temp
        integer, dimension(:), allocatable :: partition_table_temp

        if (source_sites(1) == -1) then
            allocate(source_sites_temp(0))
            allocate(source_rates_temp(0))
        else

            allocate(source_sites_temp(n_sources))

            source_sites_temp = source_sites + 1
            source_rates_temp => source_rates

        endif

        if (sink_sites(1) == -1) then
            allocate(sink_sites_temp(0))
            allocate(sink_rates_temp(0))
        else

            allocate(sink_sites_temp(n_sinks))

            sink_sites_temp = sink_sites + 1
            sink_rates_temp => sink_rates

        endif

        allocate(H%row_starts(H_rows + 1))
        allocate(H%col_indexes(H_nnz))
        allocate(H%values(H_nnz))

        H%rows = H_rows
        H%columns = H_rows
        H%row_starts = H_row_starts + 1
        H%col_indexes = H_col_indexes + 1
        H%values = H_values

        allocate(L%row_starts(H_rows + 1))
        allocate(L%col_indexes(L_nnz))
        allocate(L%values(L_nnz))

        L%rows = H_rows
        L%columns = H_rows
        L%row_starts = L_row_starts + 1
        L%col_indexes = L_col_indexes + 1
        L%values = L_values

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

        use :: iso_precisions
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

        call Reconcile_Communications(M, partition_table, MPI_communicator)

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

        use :: iso_precisions
        use :: Sparse

        implicit none

        integer, intent(in) :: M_rows
        integer, intent(in) :: M_nnz
        integer, intent(in) :: M_n_row_starts
        integer, intent(in) :: num_send
        integer, dimension(M_n_row_starts), target, intent(in) :: M_row_starts
        integer, intent(in) :: flock
        integer, dimension(M_nnz), target, intent(in) :: M_col_indexes
        integer, dimension(flock), target, intent(in) :: M_rec_disps
        integer, dimension(flock), target, intent(in) :: M_num_send_inds
        integer, dimension(flock), target, intent(in) :: M_num_rec_inds
        integer, dimension(flock), target, intent(in) :: M_send_disps
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(in) :: MPI_communicator
        integer, dimension(M_nnz), target, intent(out) :: M_local_col_inds
        integer, dimension(num_send), target, intent(out) :: M_RHS_send_inds

        type(CSR) :: M

        integer :: rank, ierr

        integer :: lb, ub
        integer :: lb_elements, ub_elements

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

        call Reconcile_Communications(M, partition_table, MPI_communicator)

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

        use :: iso_precisions
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
        complex(dp), dimension(M_n_values), target, intent(in) :: M_values
        integer, intent(in) :: flock
        integer, dimension(flock), target, intent(in) :: M_num_rec_inds
        integer, dimension(flock), target, intent(in) :: M_rec_disps
        integer, dimension(flock), target, intent(in) :: M_num_send_inds
        integer, dimension(flock), target, intent(in) :: M_send_disps
        integer, dimension(M_n_local_col_indexes), target, intent(in) :: M_local_col_inds
        integer, dimension(M_sends), target, intent(in) :: M_RHS_send_inds
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(in) :: MPI_communicator
        real(dp), dimension(9), intent(out) :: one_norm_array
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

        real(dp), dimension(9) :: alphas

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

    use :: iso_precisions
    use :: Operators

    integer, intent(in) :: rho0_rows
    integer, intent(in) :: M_local_rows
    complex(dp), dimension(rho0_rows, rho0_rows), intent(in) :: rho0
    integer, intent(in) :: flock
    integer, intent(in) :: rank
    integer, dimension(flock + 1), intent(in) :: partition_table
    integer, intent(in) :: MPI_communicator
    complex(dp), dimension(M_local_rows), &
        intent(out) :: rho0_v

    complex(dp), dimension(:), allocatable :: rho0_v_temp

    integer :: i

    call Vectorize_Operator(rho0, &
                            partition_table, &
                            rho0_v_temp, &
                            MPI_communicator)

    do i = 1, M_local_rows
        rho0_v(i) = rho0_v_temp(i + partition_table(rank + 1) - 1)
    enddo

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

        use :: iso_precisions
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
        complex(dp), dimension(M_n_values), target, intent(in) :: M_values
        integer, intent(in) :: flock
        integer, dimension(flock), target, intent(in) :: M_num_rec_inds
        integer, dimension(flock), target, intent(in) :: M_rec_disps
        integer, dimension(flock), target, intent(in) :: M_num_send_inds
        integer, dimension(flock), target, intent(in) :: M_send_disps
        integer, dimension(M_n_local_col_indexes), target, intent(in) :: M_local_col_inds
        integer, dimension(M_sends), target, intent(in) :: M_RHS_send_inds
        real(dp), intent(in) :: t
        complex(dp), dimension(n_rho0_v), intent(in) :: rho0_v
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(inout) :: p
        real(dp), dimension(9), intent(inout) :: one_norm_array
        integer, intent(in) :: MPI_communicator
        complex(dp), dimension(n_rhot_v), intent(out) :: rhot_v
        character(len=2), intent(in) :: target_precision

        type(CSR) :: M

        integer :: rank, ierr

        integer :: lb, ub
        integer :: lb_elements, ub_elements

        real(dp) :: start, finish

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


        start = MPI_wtime()
        call Expm_Multiply( M, &
                            rho0_v, &
                            t, &
                            partition_table, &
                            rhot_v, &
                            MPI_communicator, &
                            one_norm_series = one_norm_array, &
                            p = p, &
                            target_precision = target_precision)
        finish = MPI_wtime()

end subroutine step

subroutine gather_step( M_local_rows, &
                        rhot_v, &
                        flock, &
                        partition_table, &
                        root, &
                        MPI_communicator, &
                        H_aug_rows, &
                        rhot)

    use :: iso_precisions
    use :: Sparse
    use :: Operators

    implicit none

    integer, intent(in) :: M_local_rows
    complex(dp), dimension(M_local_rows), intent(in) :: rhot_v
    integer, intent(in) :: flock
    integer, dimension(flock + 1), intent(in) :: partition_table
    integer, intent(in) :: root
    integer, intent(in) :: MPI_communicator
    integer, intent(in) :: H_aug_rows
    complex(dp), dimension(H_aug_rows, H_aug_rows), intent(out) :: rhot

    complex(dp), dimension(:), allocatable :: rhot_gathered

    ! MPI ENVIRONMENT
    integer :: rank
    integer :: ierr

    call mpi_comm_rank(mpi_communicator, rank, ierr)

    allocate(rhot_gathered(partition_table(1):partition_table(size(partition_table)) - 1))


    call Gather_Dense_Vector(   rhot_v, &
                                partition_table, &
                                root, &
                                rhot_gathered, &
                                MPI_communicator)

    if (rank == root) then
        call Reshape_Vectorized_Operator(rhot_gathered, rhot)
    endif

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

        use :: iso_precisions
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
        complex(dp), dimension(M_n_values), target, intent(in) :: M_values
        integer, intent(in) :: flock
        integer, dimension(flock), target, intent(in) :: M_num_rec_inds
        integer, dimension(flock), target, intent(in) :: M_rec_disps
        integer, dimension(flock), target, intent(in) :: M_num_send_inds
        integer, dimension(flock), target, intent(in) :: M_send_disps
        integer, dimension(M_n_local_col_indexes), target, intent(in) :: M_local_col_inds
        integer, dimension(M_sends), target, intent(in) :: M_RHS_send_inds
        real(dp), intent(in) :: t_1
        real(dp), intent(in) :: t_2
        complex(dp), dimension(n_rho0_v), intent(in) :: rho0_v
        integer, intent(in) :: steps
        integer, dimension(flock + 1), intent(in) :: partition_table
        integer, intent(inout) :: p
        real(dp), dimension(9), intent(inout) :: one_norm_array
        integer, intent(in) :: MPI_communicator
        complex(dp), dimension(n_rhot_v, steps + 1), intent(out) :: rhot_v_series
        character(len = 2), intent(in) :: target_precision

        type(CSR) :: M


        integer :: rank, ierr

        integer :: lb, ub
        integer :: lb_elements, ub_elements

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
                                    target_precision = target_precision)

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

    use :: iso_precisions
    use :: Sparse
    use :: Operators

    implicit none

    integer, intent(in) :: M_local_rows
    integer, intent(in) :: steps
    complex(dp), dimension(M_local_rows, steps + 1), intent(in) :: rhot_v_series
    integer, intent(in) :: flock
    integer, dimension(flock + 1), intent(in) :: partition_table
    integer, intent(in) :: root
    integer, intent(in) :: MPI_communicator
    integer, intent(in) :: H_aug_rows
    complex(dp), dimension(H_aug_rows, H_aug_rows, steps + 1), intent(out) :: rhot_series

    complex(dp), dimension(:,:), allocatable :: rhot_v_series_gathered

    ! MPI EXP
    integer :: ierr, rank

    call mpi_comm_rank(mpi_communicator, rank, ierr)

    allocate(rhot_v_series_gathered(partition_table(1): &
        partition_table(size(partition_table)) - 1, steps + 1))

    call Gather_Dense_Matrix(   rhot_v_series, &
                                partition_table, &
                                root, &
                                rhot_v_series_gathered, &
                                MPI_communicator)


    if (rank == root) then
        call Reshape_Vectorized_Operator_Series( rhot_v_series_gathered, rhot_series)
    endif

end subroutine gather_series
