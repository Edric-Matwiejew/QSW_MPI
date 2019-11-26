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

    integer :: i

    allocate(A%row_starts(rows + 1))
    allocate(A%col_indexes(nnz))
    allocate(A%values(nnz))

    A%rows = rows
    A%columns = rows

    !$omp parallel do
    do i = 1, rows + 1
        A%row_starts(i) = row_starts_in(i) + 1
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, nnz
        A%col_indexes(i) = col_indexes_in(i) + 1
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, nnz
        A%values(i) = values_in(i)
    enddo
    !$omp end parallel do

    call Generate_Graph_Hamiltonian(gamma, A, B)

    nnz_out = size(B%col_indexes)

    !$omp parallel do
    do i = 1, nnz_out
        col_indexes_out(i) = B%col_indexes(i) - 1
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, nnz_out
        values_out(i) = B%values(i)
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, rows + 1
        row_starts_out = B%row_starts - 1
    enddo
    !$omp end parallel do

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

    integer :: i

    allocate(A%row_starts(rows + 1))
    allocate(A%col_indexes(nnz))
    allocate(A%values(nnz))

    A%rows = rows
    A%columns = rows

    !$omp parallel do
    do i = 1, rows + 1
        A%row_starts(i) = row_starts_in(i) + 1
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, nnz
        A%col_indexes(i) = col_indexes_in(i) + 1
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, nnz
        A%values(i) = values_in(i)
    enddo
    !$omp end parallel do

    call Generate_Scattering_Lindblad_Operators(A, B)

    nnz_out = size(B%col_indexes)

    !$omp parallel do
    do i = 1, nnz_out
        col_indexes_out(i) = B%col_indexes(i) - 1
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, nnz_out
        values_out(i) = B%values(i)
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, rows + 1
        row_starts_out = B%row_starts - 1
    enddo
    !$omp end parallel do

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

    type(CSR) :: A, B

    integer :: i

    allocate(A%row_starts(rows + 1))
    allocate(A%col_indexes(nnz))
    allocate(A%values(nnz))

    A%rows = rows
    A%columns = rows

    !$omp parallel do
    do i = 1, rows + 1
        A%row_starts(i) = row_starts(i) + 1
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i = 1, nnz
        A%col_indexes(i) = col_indexes(i) + 1
    enddo
    !$omp end parallel do

    A%values => values

    call Symmetrise_Graph_Weights(A, values_temp)

    !$omp parallel do
    do i = 1, nnz
        values(i) = values_temp(i)
    enddo
    !$omp end parallel do

    deallocate(A%row_starts, A%col_indexes)

end subroutine symmetrise
