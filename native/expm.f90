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

!   Module: Expm
!
!>  @brief Action of the complex matrix exponential on a vector parallalized
!>  using MPI.
!
!>  @deatils This module implements Algorithms 3.2 and 5.2 as described in
!>  "Computing the action of the matrix exponential with an application to
!>  exponential integrators" by Awad H. Al-Mohy and Nicholas J, Higham,
!>  DOI: 10.1137/100788860.

!>  @note There are two omissions. Firstly, optional balancing of the input
!>  matrix is not included. Secondly, the norm of the input matrix is not
!>  minimised via reduction of the Frobenius norm.

module Expm

    use :: ISO_Precisions
    use :: Sparse
    use :: One_Norms
    use :: MPI

    implicit none

    private

    public :: Expm_Multiply, Expm_Multiply_Series

    integer, parameter :: m_max = 55, p_max = 8

    real(qp), parameter :: tol_dp = 2.0_qp**(-53.0_qp)

    real(dp), dimension(100), target :: theta_dp = [2.220446049250313D-16, & 
        2.5809568029717673D-8, 0.000013863478661191213_dp, &
        0.0003397168839976962_dp, 0.002400876357887274_dp, 0.009065656407595102_dp, &
        0.023844555325002736_dp, 0.049912288711153226_dp, 0.08957760203223342_dp, &
        0.14418297616143777_dp, 0.21423580684517107_dp, 0.2996158913811581_dp, &
        0.3997775336316795_dp, 0.5139146936124294_dp, 0.6410835233041199_dp, &
        0.7802874256626575_dp, 0.9305328460786568_dp, 1.0908637192900361_dp, &
        1.2603810606426387_dp, 1.438252596804337_dp, 1.6237159502358214_dp, &
        1.8160778162150857_dp, 2.014710780944616_dp, 2.21904886936509_dp, &
        2.4285825244428265_dp, 2.6428534574594353_dp, 2.861449633934264_dp, &
        3.084000544989162_dp, 3.310172839890271_dp, 3.5396663487436895_dp, &
        3.772210495681751_dp, 4.00756108611804_dp, 4.245497442579696_dp, &
        4.485819859447368_dp, 4.728347345793539_dp, 4.972915626191982_dp, &
        5.219375371084059_dp, 5.467590630524544_dp, 5.717437447572013_dp, &
        5.968802630041848_dp, 6.221582661689891_dp, 6.475682736079984_dp, &
        6.731015898381024_dp, 6.98750228213063_dp, 7.245068429597952_dp, &
        7.503646685788864_dp, 7.763174657377988_dp, 8.02359472893998_dp, &
        8.284853629803916_dp, 8.546902045684934_dp, 8.809694269971322_dp, &
        9.073187890176143_dp, 9.337343505612015_dp, 9.602124472826556_dp, &
        9.8674966757534_dp, 10.133428317897478_dp, 10.399889734191031_dp, &
        10.666853220434106_dp, 10.934292878475777_dp, 11.202184475504577_dp, &
        11.470505316002537_dp, 11.739234125080184_dp, 12.008350942053166_dp, &
        12.277837023246892_dp, 12.547674753126437_dp, 12.817847562946628_dp, &
        13.088339856203294_dp, 13.359136940242902_dp, 13.630224963455024_dp, &
        13.90159085753186_dp, 14.173222284331821_dp, 14.445107586931254_dp, &
        14.717235744490083_dp, 14.989596330594331_dp, 15.262179474771681_dp, &
        15.534975826905704_dp, 15.80797652430087_dp, 16.081173161174043_dp, &
        16.35455776036932_dp,  16.628122747112073_dp, 16.901860924634942_dp, &
        17.175765451524093_dp, 17.449829820647437_dp, 17.72404783953921_dp, &
        17.998413612126303_dp, 18.272921521691813_dp, 18.547566214980513_dp, &
        18.822342587358953_dp, 19.097245768950543_dp, 19.372271111672568_dp, &
        19.647414177108324_dp, 19.922670725152962_dp, 20.19803670337687_dp, &
        20.473508237054766_dp, 20.749081619813076_dp, 21.02475330485181_dp, &
        21.300519896700695_dp, 21.576378143472468_dp, 21.852324929579098_dp, &
        22.128357268879363_dp]

    !real(dp), dimension(55), target :: theta_dp = [2.220446049250313D-16, &
    !    2.5809568029717673D-8, 1.3863478661191213D-5,  3.397168839976962D-4,  &
    !    2.400876357887274D-3, 9.065656407595102D-3,   2.3844555325002736D-2,  &
    !    4.9912288711153226D-2, 8.957760203223342D-2,   1.4418297616143777D-1, &
    !    2.1423580684517107D-1, 2.996158913811581D-1,   3.997775336316795D-1,  &
    !    5.139146936124294D-1, 6.410835233041199D-1,   7.802874256626575D-1,   &
    !    9.305328460786568D-1, 1.0908637192900361_dp,  1.2603810606426387_dp,  &
    !    1.438252596804337_dp, 1.6237159502358214_dp,  1.8160778162150857_dp,  &
    !    2.014710780944616_dp, 2.21904886936509_dp,    2.4285825244428265_dp,  &
    !    2.6428534574594353_dp, 2.861449633934264_dp,   3.084000544989162_dp,  &
    !    3.310172839890271_dp, 3.5396663487436895_dp,  3.772210495681751_dp,   &
    !    4.00756108611804_dp, 4.245497442579696_dp,   4.485819859447368_dp,    &
    !    4.728347345793539_dp, 4.972915626191982_dp,   5.219375371084059_dp,   &
    !    5.467590630524544_dp, 5.717437447572013_dp,   5.968802630041848_dp,   &
    !    6.221582661689891_dp, 6.475682736079984_dp,   6.731015898381024_dp,   &
    !    6.98750228213063_dp, 7.245068429597952_dp,   7.503646685788864_dp,    &
    !    7.763174657377988_dp, 8.02359472893998_dp,    8.284853629803916_dp,   &
    !    8.546902045684934_dp, 8.809694269971322_dp,   9.073187890176143_dp,   &
    !    9.337343505612015_dp, 9.602124472826556_dp,   9.8674966757534_dp]

    real(qp), parameter :: tol_sp = 2.0_dp**(-24.0_dp)

    real(dp), dimension(55), target :: theta_sp = [1.19209D-7, 0.000597886_dp, &
        0.0112339_dp, 0.0511662_dp, 0.130849_dp, 0.249529_dp, 0.401458_dp,     &
        0.580052_dp, 0.779511_dp, 0.995184_dp, 1.22348_dp, 1.46166_dp,         &
        1.70765_dp, 1.95985_dp, 2.21704_dp, 2.47828_dp, 2.74282_dp, 3.01007_dp,&
        3.27956_dp, 3.55093_dp, 3.82386_dp, 4.09811_dp, 4.37347_dp, 4.64978_dp,&
        4.9269_dp, 5.20471_dp, 5.48311_dp, 5.76201_dp, 6.04136_dp, 6.32108_dp, &
        6.60113_dp, 6.88146_dp, 7.16204_dp, 7.44283_dp, 7.7238_dp, 8.00493_dp, &
        8.2862_dp, 8.56759_dp, 8.84908_dp, 9.13065_dp, 9.4123_dp, 9.69402_dp,  &
        9.97579_dp, 10.2576_dp, 10.5394_dp, 10.8213_dp, 11.1032_dp, 11.3852_dp,&
        11.6671_dp, 11.949_dp, 12.231_dp, 12.5129_dp, 12.7949_dp, 13.0769_dp,  &
        13.3588_dp]

    contains

    function Local_Infinity_Norm(B) result(norm)
        complex(dp), intent(in) :: B(:)
        real(dp) :: norm
        integer :: i

        norm = 0
        !$omp parallel do reduction(max:norm)
        do i = 1, size(B)
            norm = max(norm, abs(B(i)))
        enddo
        !$omp end parallel do
    end function Local_Infinity_Norm

    function Infinity_Norm(B, comm) result(norm)
        complex(dp), intent(in) :: B(:)
        integer, intent(in) :: comm
        real(dp) :: local, norm
        integer :: ierr

        local = Local_Infinity_Norm(B)
        call MPI_Allreduce(local, norm, 1, MPI_double, MPI_max, comm, ierr)
    end function Infinity_Norm

    subroutine Taylor_Norms(term, partial_sum, norms, comm)
        complex(dp), intent(in) :: term(:), partial_sum(:)
        real(dp), intent(out) :: norms(2)
        integer, intent(in) :: comm
        real(dp) :: local(2)
        integer :: ierr

        local(1) = Local_Infinity_Norm(term)
        local(2) = Local_Infinity_Norm(partial_sum)
        call MPI_Allreduce(local, norms, 2, MPI_double, MPI_max, comm, ierr)
    end subroutine Taylor_Norms

    subroutine C_m( A, &
                    t, &
                    target_precision, &
                    partition_table, &
                    m_star, &
                    s, &
                    mpi_communicator, &
                    one_norm_series, &
                    p_in)

        type(CSR), intent(inout) :: A
        real(dp), intent(in) :: t
        character(len = 2), intent(in) :: target_precision
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(out) :: m_star, s
        integer, intent(in) :: mpi_communicator
        real(dp), dimension(p_max + 1), optional, intent(inout) :: one_norm_series
        integer, intent(inout), optional :: p_in

        real(dp), dimension(:), pointer :: theta => null()

        real(dp), dimension(2:p_max + 1) :: A_norms
        integer(qp) :: cost, min_cost, scales
        integer :: p, m, k

        if (target_precision == "sp") then
            theta => theta_sp
        elseif (target_precision == "dp") then
            theta => theta_dp
        endif

        if (present(p_in)) then
            p = p_in
        else
            p = 0
        endif

        if (p <= 0 .or. .not. present(one_norm_series)) then
            call One_Norm_Powers(A, partition_table, A_norms, mpi_communicator)
            p = p_max
            if (present(one_norm_series)) then
                one_norm_series(2:) = A_norms
                if (present(p_in)) p_in = p
            endif
        else
            A_norms(2:p+1) = one_norm_series(2:p+1)
        endif

        do k = 2, p + 1
            A_norms(k) = A_norms(k)**(1.0_dp/k)
        enddo
        min_cost = huge(min_cost)
        do k = 2, p
            do m = k*(k - 1) - 1, m_max
                scales = max(1_qp, ceiling(abs(t)*max(A_norms(k), A_norms(k + 1))/theta(m), kind=qp))
                cost = int(m, qp)*scales
                if (cost < min_cost) then
                    min_cost = cost
                    m_star = m
                endif
            enddo
        enddo

        s = int(min_cost/m_star)

    end subroutine C_m

    subroutine Parameters(  A, &
                            t, &
                            target_precision, &
                            partition_table, &
                            m_star, &
                            s, &
                            mpi_communicator, &
                            one_norm_series, &
                            p)

        type(CSR), intent(inout) :: A
        real(dp), intent(in) :: t
        character(len=2), intent(in) :: target_precision
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(out) :: m_star, s
        integer, intent(in) :: mpi_communicator
        real(dp), optional, dimension(p_max + 1), intent(inout) :: one_norm_series
        integer, optional, intent(inout) :: p

        real(dp), dimension(:), pointer :: theta

        real(dp) :: A_norm

        integer(dp) :: m_temp_1, m_temp_2

        integer :: m

        logical :: norm_known

        if (target_precision == "sp") then
            theta => theta_sp
        elseif (target_precision == "dp") then
            theta => theta_dp
        endif

        ! p is zero before norm estimation and -1 when only the 1-norm is known.
        ! A positive p includes estimates through power p + 1.
        norm_known = .false.
        if (present(p) .and. present(one_norm_series)) norm_known = p /= 0
        if (norm_known) then
            A_norm = one_norm_series(1)
        else
            call One_Norm(A, A_norm, partition_table, MPI_communicator)
            if (present(one_norm_series)) then
                one_norm_series(1) = A_norm
                if (present(p)) p = -1
            endif
        endif

        A_norm = abs(t)*A_norm
        if (A_norm <= 2*(theta(m_max)/real(m_max,dp))*(p_max + 3)*p_max) then
            m_star = 1
            m_temp_1 = max(1_dp, ceiling(A_norm/theta(1), kind=dp))
            do m = 2, m_max
                m_temp_2 = m*max(1_dp, ceiling(A_norm/theta(m), kind=dp))
                if (m_temp_2 < m_temp_1) then
                    m_star = m
                    m_temp_1 = m_temp_2
                endif
            enddo
            s = int(m_temp_1/m_star)
        else
            call C_m(A, t, target_precision, partition_table, m_star, s, mpi_communicator, one_norm_series, p)
        endif

    end subroutine  Parameters

    subroutine Expm_Multiply(   A, &
                                B, &
                                t, &
                                partition_table, &
                                C, &
                                mpi_communicator, &
                                one_norm_series, &
                                p, &
                                target_precision)

        type(CSR), intent(inout) :: A
        complex(dp), intent(in), dimension(:) :: B
        real(dp), intent(in) :: t
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(in) :: mpi_communicator
        complex(dp), dimension(:), intent(inout) :: C
        real(dp), dimension(p_max + 1), optional, intent(inout) :: one_norm_series
        integer, optional, intent(inout) :: p
        character(len=2), optional, intent(in) :: target_precision

        complex(dp), dimension(:), allocatable :: B_temp_1, B_temp_2

        integer :: m_star, s

        character(len=2) :: set_target_precision
        real(qp) :: tol
        real(dp) :: norms(2)
        real(qp) :: c_1, c_2

        integer :: i, j, lb, ub

        integer :: rank
        integer :: ierr

        tol = 0

        call mpi_comm_rank(mpi_communicator, rank, ierr)

        if (present(target_precision)) then
            if (target_precision == "sp") then
                set_target_precision = target_precision
                tol = tol_sp
            elseif (target_precision == "dp") then
                set_target_precision = target_precision
                tol = tol_dp
            endif
        else
            set_target_precision = "dp"
            tol = tol_dp
        endif

        lb = lbound(C,1)
        ub = ubound(C,1)

        if (t == 0.0_dp) then
            C = B
            return
        endif
        call Parameters(A, t, set_target_precision, partition_table, m_star, s, mpi_communicator, one_norm_series, p)
        allocate(B_temp_1(lb:ub), B_temp_2(lb:ub))

        B_temp_1(lb:ub) = B(lb:ub)
        C(lb:ub) = B(lb:ub)

        do i = 1, s

            c_1 = Infinity_Norm(B_temp_1, MPI_communicator)

            do j = 1, m_star

                call SpMV_Series(   A, &
                                    B_temp_1, &
                                    partition_table, &
                                    1, &
                                    j, &
                                    m_star, &
                                    rank, &
                                    B_temp_2, &
                                    mpi_communicator)

                B_temp_2(lb:ub) = (t/real(s, dp)/real(j, dp))*B_temp_2(lb:ub)

                C(lb:ub) = C(lb:ub) + B_temp_2(lb:ub)

                call Taylor_Norms(B_temp_2, C, norms, MPI_communicator)
                c_2 = norms(1)
                if (c_1 + c_2 <= tol*norms(2)) then
                    exit
                endif

                B_temp_1(lb:ub) = B_temp_2(lb:ub)

                c_1 = c_2

            enddo

            B_temp_1(lb:ub) = C(lb:ub)

        enddo

        ! Final call to deallocate saved arrays.
        call SpMV_Series(   A, &
                            B_temp_1, &
                            partition_table, &
                            0, &
                            0, &
                            0, &
                            rank, &
                            B_temp_2, &
                            mpi_communicator)

    end subroutine Expm_Multiply

    subroutine Expm_Multiply_Series(    A, &
                                        B, &
                                        t0, &
                                        tq, &
                                        steps, &
                                        partition_table, &
                                        X, &
                                        mpi_communicator, &
                                        one_norm_series_in, &
                                        p_ex, &
                                        target_precision)

        real(dp), intent(in) :: t0, tq
        type(CSR), intent(inout) :: A
        complex(dp), intent(in), dimension(:) :: B
        integer, intent(in) :: steps
        integer, dimension(:), intent(in) :: partition_table
        complex(dp), dimension(:,:), intent(out) :: X
        integer, intent(in) :: mpi_communicator
        real(dp), dimension(p_max + 1), optional, target, intent(inout) :: one_norm_series_in
        integer, optional, target, intent(inout) :: p_ex
        character(len=2), optional, intent(in) :: target_precision

        real(dp), target :: local_norms(p_max + 1)
        integer, target :: local_p
        real(dp), pointer :: one_norm_series(:)
        integer, pointer :: p_in

        integer :: q

        complex(dp), dimension(:), allocatable :: Z, F
        complex(dp), dimension(:,:), allocatable :: K
        real(dp) :: h

        character(len=2) :: set_target_precision
        real(qp) :: tol
        real(dp) :: norms(2)
        real(qp) :: c_1, c_2

        integer :: m_star, s, m_hat
        integer :: d, j, r, d_tilde, p

        ! MPI ENVIRONMENT
        integer :: rank
        integer :: ierr

        integer :: i, kay, indx, lb, ub

        tol = 0
        indx = 0

        call mpi_comm_rank(MPI_communicator, rank, ierr)

        if (present(target_precision)) then
            if (target_precision == "sp") then
                set_target_precision = target_precision
                tol = tol_sp
            elseif (target_precision == "dp") then
                set_target_precision = target_precision
                tol = tol_dp
            endif
        else
            set_target_precision = "dp"
            tol = tol_dp
        endif

        lb = partition_table(rank + 1)
        ub = partition_table(rank + 2) - 1

        allocate(Z(ub - lb + 1))
        allocate(F(ub - lb + 1))

        if (present(one_norm_series_in) .and. present(p_ex)) then
            one_norm_series => one_norm_series_in
            p_in => p_ex
        else
            local_p = 0
            p_in => local_p
            one_norm_series => local_norms
        endif

        call Expm_Multiply( A, &
                            B, &
                            t0, &
                            partition_table, &
                            X(:,1), &
                            mpi_communicator, &
                            one_norm_series = one_norm_series, &
                            p = p_in, &
                            target_precision = set_target_precision)

        if (steps == 0) return

        q = steps
        h = (tq - t0)/real(q, dp)

        call Parameters(A, &
                        tq - t0, &
                        set_target_precision, &
                        partition_table, &
                        m_star, &
                        s, &
                        mpi_communicator, &
                        one_norm_series = one_norm_series, &
                        p = p_in)

        if (q <= s) then

            do kay = 1, q

                call Expm_Multiply( A, &
                                    X(:,kay), &
                                    h, &
                                    partition_table, &
                                    X(:,kay + 1), &
                                    MPI_communicator, &
                                    one_norm_series = one_norm_series, &
                                    p = p_in, &
                                    target_precision = set_target_precision)

            enddo

            return

        endif

        d = q/s
        j = q/d
        r = q - d*j
        d_tilde = d

        ! The time span of each block is at most abs(tq - t0)/s.
        Z = X(:, 1)

        allocate(K(ub - lb + 1, m_star + 1))

        K = 0

        do i = 1, j + 1

            if (i > j) then
               d_tilde = r
            endif

            K(:, 1) = Z
            k(:,2:m_star+1) = 0

            m_hat = 0

            do kay = 1, d_tilde

                F = Z

                c_1 = Infinity_Norm(Z, MPI_communicator)

                do p = 1, m_star

                    if (p > m_hat) then

                        call SpMV_series(   A, &
                                            K(:, p), &
                                            partition_table, &
                                            1, &
                                            p, &
                                            m_star, &
                                            rank, &
                                            K(:, p + 1), &
                                            MPI_communicator)

                        K(:, p + 1) = h*K(:, p+1)/real(p,dp)

                    endif

                    F = F + (real(kay,dp)**real(p, dp))*K(:,p+1)

                    call Taylor_Norms(K(:,p + 1), F, norms, MPI_communicator)
                    c_2 = (real(kay,qp)**real(p, dp))*norms(1)

                    indx = p

                    if (c_1 + c_2 <= tol*norms(2)) then

                        exit

                    endif

                    c_1 = c_2

                enddo

                m_hat = max(m_hat, indx)

                X(:, kay + (i - 1)*d + 1) =  F

            enddo

            if (i <= j) then
                Z = X(:, i*d + 1)
            endif

        enddo

        call SpMV_Series(A, Z, partition_table, 0, 0, 0, rank, F, mpi_communicator)

    end subroutine Expm_Multiply_Series

end module Expm
