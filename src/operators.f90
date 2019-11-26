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
    !   Function: Kronecker_Delta
    !
    !>  @brief Kronecker delta function of integers i and j.
    !
    !>  @details This function is equal to 1 if i = j and 0 otherwise.

    function Kronecker_Delta(i, j)

        integer :: Kronecker_Delta
        integer, intent(in) :: i, j

        Kronecker_Delta = int((real((i+j)-abs(i-j)))/(real((i+j)+abs(i-j))))

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

            do i = 2, size(array)
                array(i) = array(i - 1) + array(i)
            enddo

    end subroutine Prefix_Sum

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

        type(CSR), intent(inout) :: A !< Structurally symmetric CSR array.
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

    subroutine Generate_Scattering_Lindblad_Operators(A, L)

        type(CSR), intent(in) :: A
        type(CSR), intent(out) :: L

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

    subroutine Generate_Scattering_Superoperator(L, lower_bound, upper_bound, B)

        type(CSR), intent(in) :: L
        integer, intent(in) :: upper_bound, lower_bound
        type(CSR), intent(out) :: B

        integer :: lower_block, upper_block
        complex(dp), dimension(:), allocatable :: value_elements
        complex(dp), dimension(:), allocatable :: value_elements_shifted
        complex(dp), dimension(:), allocatable :: nz_diag
        integer, dimension(:), allocatable :: col_index
        integer, dimension(:), allocatable :: row_start

        integer :: lower_offset, upper_offset

        integer :: i, j, k, row, indx

        lower_block = ceiling(real(lower_bound)/L%rows)
        upper_block = ceiling(real(upper_bound)/L%rows)

        allocate(value_elements(L%rows))
        allocate(nz_diag(L%rows))

        value_elements = 0
        nz_diag(i) = 0

        do i = 1, L%rows
            do j = L%row_starts(i), L%row_starts(i + 1) - 1
                value_elements(L%col_indexes(j)) = value_elements(L%col_indexes(j)) &
                    + conjg(L%values(j))*L%values(j)
                if (L%col_indexes(j) == i) then
                    nz_diag(i) = L%values(j)
                endif
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
        value_elements_shifted(1:L%rows - 1) = value_elements(i)

        lower_offset = (lower_bound - (lower_block - 1)*L%rows) - 1
        upper_offset = (upper_bound - (upper_block - 1)*L%rows) - L%rows

        do i = lower_block, upper_block
            do j = 1 + Kronecker_Delta(i, lower_block)*lower_offset, L%rows &
                + Kronecker_Delta(i, upper_block)*upper_offset

                row = (i - 1)*L%rows + j

                if ((i == j) .and. (L%row_starts(i + 1) - L%row_starts(i) /= 0)) then

                    indx = 0

                    do k = L%row_starts(i), L%row_starts(i + 1) - 1

                        if ((L%col_indexes(k) - 1)*L%rows + L%col_indexes(k) == row) then

                            B%values(row_start(row) + indx) =&
                                -0.5_dp*(value_elements_shifted(mod(row, L%rows)) + &
                                & value_elements(ceiling(real(row)/real(L%rows)))) + nz_diag(i)**2
                            B%col_indexes(row_start(row) + indx) = &
                                (L%col_indexes(k) - 1)*L%rows + L%col_indexes(k)

                        else

                            B%values(row_start(row) + indx) = L%values(k)**2
                            B%col_indexes(row_start(row) + indx) = &
                                (L%col_indexes(k) - 1)*L%rows + L%col_indexes(k)

                        endif

                        indx = indx + 1

                    enddo
                else

                    B%values(row_start(row)) = &
                        -0.5_dp*(value_elements_shifted(mod(row, L%rows)) + &
                                & value_elements(ceiling(real(row)/real(L%rows))))
                    B%col_indexes(row_start(row)) = row

                endif
            enddo
        enddo

        allocate(B%row_starts(lower_bound:upper_bound + 1))

        B%row_starts(lower_bound:upper_bound + 1) = row_start(lower_bound:upper_bound + 1)
        B%rows = upper_bound - lower_bound + 1
        B%columns = L%columns**2

    end subroutine Generate_Scattering_Superoperator

    subroutine Generate_Coherent_Superoperator( H, &
                                                lower_bound, &
                                                upper_bound, &
                                                B)

        type(CSR), intent(in) :: H
        integer, intent(in) :: lower_bound, upper_bound
        type(CSR), intent(out) :: B

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

    subroutine Generate_Graph_Hamiltonian(gamma, A, B)

        real(dp), intent(in) :: gamma
        type(CSR), intent(in) :: A
        type(CSR), intent(out) :: B

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

        !$omp parallel do
        do i = 1, A%rows
            out_degrees(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do reduction(+:out_degrees)
        do i = 1, A%rows
            do j = A%row_starts(i), A%row_starts(i + 1) - 1
                out_degrees(i) = out_degrees(i) + gamma*A%values(j)
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, A%rows + 1
            B%row_starts(i) = A%row_starts(i)
        enddo
        !$omp end parallel do

        do i = 1, A%rows
            if (abs(out_degrees(i)) > epsilon(0.d0)) then
                B%row_starts(i + 1 : A%rows + 1) = B%row_starts(i + 1 : A%rows + 1) + 1
            endif
        enddo

        nnz = B%row_starts(B%rows + 1) - 1

        allocate(B%col_indexes(nnz))
        allocate(B%values(nnz))

        !$omp parallel do private(elements)
        do i = 1, B%rows

            elements = A%row_starts(i + 1) - A%row_starts(i)

            if (elements == 0) cycle

            B%col_indexes(B%row_starts(i):B%row_starts(i) + elements - 1) = &
                A%col_indexes(A%row_starts(i):A%row_starts(i + 1) - 1)

            B%values(B%row_starts(i):B%row_starts(i) + elements - 1) = &
                -gamma*A%values(A%row_starts(i):A%row_starts(i + 1) - 1)

            B%col_indexes(B%row_starts(i + 1) - 1) = i

            B%values(B%row_starts(i + 1) - 1) = out_degrees(i)

        enddo
        !$omp end parallel do

        !$omp parallel do private(i, temp, temp_indx, j)
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
        !$omp end parallel do

    end subroutine Generate_Graph_Hamiltonian

    subroutine Generate_Stochastic_Superoperator(H, &
                                    L, &
                                    omega, &
                                    lower_bound, &
                                    upper_bound, &
                                    B)

        type(CSR), intent(in) :: H, L
        real(dp), intent(in) :: omega
        integer, intent(in) :: lower_bound, upper_bound
        type(CSR), intent(out) :: B

        integer :: lower_block, upper_block

        integer, dimension(:), allocatable :: H_left, H_right, L_left, L_right
        complex(dp), dimension(:), allocatable :: H_diag_vals, L_diag_vals
        complex(dp), dimension(:), allocatable :: out_degrees
        complex(dp), dimension(:), allocatable :: diag_vals

        integer, dimension(:), allocatable :: row_start

        integer :: lower_offset, upper_offset
        integer :: i, j, row, lower_indx, upper_indx

        integer :: N

        !sort
        !complex(dp) :: temp
        !integer :: temp_indx
        !integer :: k

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

        !$omp parallel do
        do i = 1, N
            H_left = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            H_right(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            L_left(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            L_right(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            H_diag_vals(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            L_diag_vals(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            out_degrees(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, H%rows**2
            diag_vals(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do private(j)
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
        !$omp end parallel do

        !$omp parallel do private(j)
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
        !$omp end parallel do

        !$omp parallel do reduction(+:out_degrees)
        do i = 1, size(L%values, 1)
            out_degrees(L%col_indexes(i)) = out_degrees(L%col_indexes(i)) + &
                omega*L%values(i)*conjg(L%values(i))
        enddo
        !$omp end parallel do

        !$omp parallel do private(row)
        do i = 1, N
           row = (i - 1)*N + i
            diag_vals(row) = omega*L_diag_vals(i)*conjg(L_diag_vals(i))
        enddo
        !$omp end parallel do

        !$omp parallel do private(j, row)
        do i = 1, N
             do j = 1, N
                row = (i - 1)*N + j
                diag_vals(row) = diag_vals(row) - 0.5_dp*(out_degrees(i) + out_degrees(j))
            enddo
        enddo
        !$omp end parallel do

        do i = 1, N
            diag_vals((i - 1)*N + 1:i*N) = &
                diag_vals((i - 1)*N + 1:i*N) + H_diag_vals - H_diag_vals(i)
        enddo

        allocate(row_start(N**2 + 1))

        row_start(1) = 1

        !$omp parallel do
        do i = 2, N**2 + 1
            row_start(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do private(row, j)
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
        !$omp end parallel do

        call Prefix_Sum(row_start)

        allocate(B%values(row_start(lower_bound):row_start(upper_bound + 1) - 1))
        allocate(B%col_indexes(row_start(lower_bound):row_start(upper_bound + 1) - 1))

        lower_offset = (lower_bound -(lower_block - 1)*N) - 1
        upper_offset = (upper_bound - (upper_block - 1)*N) - N

        !$omp parallel do private(j, row, lower_indx, upper_indx)
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
        !$omp end parallel do

        allocate(B%row_starts(lower_bound:upper_bound + 1))

        B%row_starts = row_start(lower_bound:upper_bound + 1)

        B%rows = upper_bound - lower_bound + 1
        B%columns = H%columns**2

        call Sort_CSR(B)

    end subroutine Generate_Stochastic_Superoperator

    subroutine Add_Sources_and_Sinks(   L, &
                                        L_aug, &
                                        source_sites, &
                                        source_rates, &
                                        sink_sites, &
                                        sink_rates)

        type(CSR), intent(in) :: L
        type(CSR), intent(out) :: L_aug
        integer, dimension(:), optional :: source_sites
        real(dp), dimension(:), optional :: source_rates
        integer, dimension(:), optional :: sink_sites
        real(dp), dimension(:), optional :: sink_rates

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
        L_aug%columns = L_aug%rows + n_sinks + n_sources

        !$omp parallel do
        do i = 1, L%rows + 1
            L_aug%row_starts(i) = L%row_starts(i)
        enddo
        !$omp end parallel do

        L_aug%row_starts(L%rows + 2:L_aug%rows + 1) = L%row_starts(L%rows + 1)

        !$omp parallel do
        do i = 1, size(L_aug%col_indexes)
            L_aug%col_indexes(i) = 0
        enddo
        !$omp end parallel do

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

        !$omp parallel do private(elements)
        do i = 1, L%rows

            elements = L%row_starts(i + 1) - L%row_starts(i)

            L_aug%col_indexes(L_aug%row_starts(i):L_aug%row_starts(i) + elements - 1) = &
                L%col_indexes(L%row_starts(i):L%row_starts(i + 1) - 1)

            L_aug%values(L_aug%row_starts(i):L_aug%row_starts(i) + elements - 1) =&
                L%values(L%row_starts(i):L%row_starts(i + 1) - 1)

        enddo
        !$omp end parallel do

        !$omp parallel do private(i, temp, temp_indx, j)
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
        !$omp end parallel do

    end subroutine Add_Sources_and_Sinks

    ! NOTE: The sink and site variables are mixed up in this subroutine.
    ! as a quick fix I have swappwed them as input arguments.
    subroutine Transport_Lindblads( N, &
                                    sink_sites, &
                                    sink_rates, &
                                    source_sites, &
                                    source_rates, &
                                    L_inout)

        integer, intent(in) :: N
        integer, dimension(:), intent(in) :: source_sites
        real(dp), dimension(:), intent(in) :: source_rates
        integer, dimension(:), intent(in) :: sink_sites
        real(dp), dimension(:), intent(in) :: sink_rates
        type(CSR), intent(out) :: L_inout

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

        !$omp parallel do
        do i = 2, L_inout%rows + 1
            L_inout%row_starts(i) = 0
        enddo
        !$omp end parallel do

        L_inout%row_starts(1) = 1

        do i = 1, size(sink_sites)
            L_inout%row_starts(sink_sites(i)+1) =  1
        enddo

        do i = 1, size(source_sites)
            L_inout%row_starts(i + N + size(sink_sites) + 1) = 1
        enddo

        call Prefix_Sum(L_inout%row_starts)

    end subroutine Transport_Lindblads

    subroutine Generate_Coherent_Transport_Superoperator(H, &
                                    L, &
                                    lower_bound, &
                                    upper_bound, &
                                    B)

        type(CSR), intent(in) :: H, L
        integer, intent(in) :: lower_bound, upper_bound
        type(CSR), intent(out) :: B

        integer :: lower_block, upper_block

        integer, dimension(:), allocatable :: H_left, H_right, L_left, L_right
        complex(dp), dimension(:), allocatable :: H_diag_vals, L_diag_vals
        complex(dp), dimension(:), allocatable :: out_degrees
        complex(dp), dimension(:), allocatable :: diag_vals

        integer, dimension(:), allocatable :: row_start

        integer :: lower_offset, upper_offset
        integer :: i, j, row, lower_indx, upper_indx

        integer :: N

        !sort
        !complex(dp) :: temp
        !integer :: temp_indx
        !integer :: k

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

        !$omp parallel do
        do i = 1, N
            H_left = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
        H_right(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            L_left(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            L_right(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            H_diag_vals(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            L_diag_vals(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            out_degrees(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, H%rows**2
            diag_vals(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do private(j)
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
        !$omp end parallel do

        !$omp parallel do private(j)
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
        !$omp end parallel do

        !$omp parallel do reduction(+:out_degrees)
        do i = 1, size(L%values, 1)
            out_degrees(L%col_indexes(i)) = out_degrees(L%col_indexes(i)) + &
                L%values(i)*conjg(L%values(i))
        enddo
        !$omp end parallel do

        !$omp parallel do private(row)
        do i = 1, N
           row = (i - 1)*N + i
            diag_vals(row) = L_diag_vals(i)*conjg(L_diag_vals(i))
        enddo
        !$omp end parallel do

        !$omp parallel do private(row)
        do i = 1, N
             do j = 1, N
                row = (i - 1)*N + j
                diag_vals(row) = diag_vals(row) - 0.5_dp*(out_degrees(i) + out_degrees(j))
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = 1, N
            diag_vals((i - 1)*N + 1:i*N) = &
                diag_vals((i - 1)*N + 1:i*N) + H_diag_vals - H_diag_vals(i)
        enddo
        !$omp end parallel do

        allocate(row_start(N**2 + 1))

        row_start(1) = 1

        !$omp parallel do
        do i = 2, N**2 + 1
            row_start(i) = 0
        enddo
        !$omp end parallel do

        !$omp parallel do private(row, j)
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
        !$omp end parallel do

        call Prefix_Sum(row_start)

        allocate(B%values(row_start(lower_bound):row_start(upper_bound + 1) - 1))
        allocate(B%col_indexes(row_start(lower_bound):row_start(upper_bound + 1) - 1))

        lower_offset = (lower_bound -(lower_block - 1)*N) - 1
        upper_offset = (upper_bound - (upper_block - 1)*N) - N

        !$omp parallel do private (j, row, lower_indx, upper_indx)
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
        !$omp end parallel do

        allocate(B%row_starts(lower_bound:upper_bound + 1))

        !$omp parallel do
        do i = lower_bound, upper_bound + 1
            B%row_starts(i) = row_start(i)
        enddo
        !$omp end parallel do

        B%rows = upper_bound - lower_bound + 1
        B%columns = H%columns**2

        call Sort_CSR(B)

    end subroutine Generate_Coherent_Transport_Superoperator

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

        real(dp), intent(in) :: omega
        type(CSR), intent(in) :: H
        type(CSR), intent(in) :: L
        integer, dimension(:), intent(in) :: source_sites
        real(dp), dimension(:), intent(in) :: source_rates
        integer, dimension(:), intent(in) :: sink_sites
        real(dp), dimension(:), intent(in) :: sink_rates
        integer, dimension(:), allocatable, intent(out) :: partition_table
        type(CSR), intent(out) :: M_local
        integer, intent(in) :: MPI_communicator

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

    subroutine Vectorize_Operator(  O, &
                                    partition_table, &
                                    O_v, &
                                    MPI_communicator)

        complex(dp), dimension(:, :), intent(in) :: O
        integer, dimension(:), intent(in) :: partition_table
        complex(dp), dimension(:), allocatable, intent(out) :: O_v
        integer, intent(in) :: MPI_communicator

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

        !$omp parallel do private(j, indx)
        do i = L_row, U_row

            do j = 1 + Kronecker_Delta(i, L_row)*L_shift, &
                O_size - Kronecker_Delta(i, U_row)*U_shift

                indx = (i - 1)*size(O, 1) + j

                O_v(indx) = O(i, j)

            enddo

        enddo
        !$omp end parallel do

    end subroutine Vectorize_Operator

    subroutine Reshape_Vectorized_Operator_Series( O_v_series, O_series)

        complex(dp), dimension(:,:), intent(in) :: O_v_series
        complex(dp), dimension(:, :, :), intent(out) :: O_series

        integer :: O_dim
        integer :: steps
        integer :: i, j, k
        integer :: indx

        O_dim = int(sqrt(real(size(O_v_series, 1))))
        steps = size(O_v_series, 2)


        !$omp parallel do private(i, j, indx)
        do k = 1, steps
            do i = 1, O_dim
                do j = 1, O_dim
                    indx = (i - 1)*O_dim + j
                    O_series(i, j, k) = O_v_series(indx, k)
                enddo
            enddo
        enddo
        !$omp end parallel do

    end subroutine Reshape_Vectorized_Operator_Series

    subroutine Reshape_Vectorized_Operator( O_v, O)

        complex(dp), dimension(:), intent(in) :: O_v
        complex(dp), dimension(:, :), intent(out) :: O

        integer :: O_dim
        integer :: i, j
        integer :: indx

        O_dim = int(sqrt(real(size(O_v))))

        !$omp parallel do private(j, indx)
        do i = 1, O_dim
            do j = 1, O_dim
                indx = (i - 1)*O_dim + j
                O(i, j) = O_v(indx)
            enddo
        enddo
        !$omp end parallel do

    end subroutine Reshape_Vectorized_Operator

end module Operators
