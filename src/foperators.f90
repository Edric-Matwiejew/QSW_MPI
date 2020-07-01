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

subroutine graph(   gamma, &
                    rows, &
                    nnz, &
                    row_starts_in, &
                    col_indexes_in, &
                    values_in, &
                    row_starts_out, &
                    col_indexes_out, &
                    values_out, &
                    nnz_out)

    use :: Sparse
    use :: Operators

    implicit none

    real(8), intent(in) :: gamma
    integer , intent(in) :: rows
    integer, intent(in) :: nnz
    integer, dimension(rows + 1), intent(in) :: row_starts_in
    integer, dimension(nnz), intent(in) :: col_indexes_in
    complex(8), dimension(nnz), intent(in) :: values_in
    integer, dimension(rows + 1), intent(out) :: row_starts_out
    integer, dimension(nnz + rows), intent(out) :: col_indexes_out
    complex(8), dimension(nnz + rows), intent(out) :: values_out
    integer, intent(out) :: nnz_out

    type(CSR) :: A, B

    allocate(A%row_starts(rows + 1))
    allocate(A%col_indexes(nnz))
    allocate(A%values(nnz))

    A%rows = rows
    A%columns = rows

    A%row_starts = row_starts_in + 1
    A%col_indexes = col_indexes_in + 1
    A%values = values_in

    call Generate_Graph_Hamiltonian(gamma, A, B)

    nnz_out = size(B%col_indexes)

    col_indexes_out(1:nnz_out) = B%col_indexes - 1
    values_out(1:nnz_out) = B%values
    row_starts_out = B%row_starts - 1

    deallocate(A%row_starts, A%col_indexes, A%values)

end subroutine graph

subroutine site_lindblads(  rows, &
                            nnz, &
                            row_starts_in, &
                            col_indexes_in, &
                            values_in, &
                            row_starts_out, &
                            col_indexes_out, &
                            values_out, &
                            nnz_out)

    use :: Sparse
    use :: Operators

    implicit none

    integer , intent(in) :: rows
    integer, intent(in) :: nnz
    integer, dimension(rows + 1), intent(in) :: row_starts_in
    integer, dimension(nnz), intent(in) :: col_indexes_in
    complex(8), dimension(nnz), intent(in) :: values_in
    integer, dimension(rows + 1), intent(out) :: row_starts_out
    integer, dimension(nnz + rows), intent(out) :: col_indexes_out
    complex(8), dimension(nnz + rows), intent(out) :: values_out
    integer, intent(out) :: nnz_out

    type(CSR) :: A, B

    allocate(A%row_starts(rows + 1))
    allocate(A%col_indexes(nnz))
    allocate(A%values(nnz))

    A%rows = rows
    A%columns = rows

    A%row_starts = row_starts_in + 1
    A%col_indexes = col_indexes_in + 1
    A%values = values_in

    call Generate_Scattering_Lindblad_Operators(A, B)

    nnz_out = size(B%col_indexes)

    col_indexes_out(1:nnz_out) = B%col_indexes - 1
    values_out(1:nnz_out) = B%values
    row_starts_out = B%row_starts - 1

    deallocate(A%row_starts, A%col_indexes, A%values)

end subroutine site_lindblads

subroutine symmetrise(  rows, &
                        nnz, &
                        row_starts, &
                        col_indexes, &
                        values)

    use :: Sparse
    use :: Operators

    implicit none

    integer , intent(in) :: rows
    integer, intent(in) :: nnz
    integer, dimension(rows + 1), intent(in) :: row_starts
    integer, dimension(nnz), intent(in) :: col_indexes
    complex(8), dimension(nnz), intent(inout), target :: values

    complex(8), dimension(nnz) :: values_temp

    type(CSR) :: A

    allocate(A%row_starts(rows + 1))
    allocate(A%col_indexes(nnz))
    allocate(A%values(nnz))

    A%rows = rows
    A%columns = rows

    A%row_starts = row_starts + 1
    A%col_indexes = col_indexes + 1
    A%values => values

    call Symmetrise_Graph_Weights(A, values_temp)

    values = values_temp

    deallocate(A%row_starts, A%col_indexes)

end subroutine symmetrise
