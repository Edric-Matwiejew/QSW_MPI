subroutine sub_2(a_size, b_size, a, b, d)

    integer, intent(in) :: a_size
    integer, intent(in) :: b_size
    integer, dimension(a_size), intent(in) :: a
    integer, dimension(b_size), intent(in) :: b
    real, intent(out) :: d

    d = size(a) + size(b)

end subroutine sub_2

subroutine sub_3(a_size, b_size, c_size, a, b, c, d)

    integer, intent(in) :: a_size
    integer, intent(in) :: b_size
    integer, intent(in) :: c_size
    integer, dimension(a_size), intent(in) :: a
    integer, dimension(b_size), intent(in) :: b
    real, dimension(c_size), intent(in) :: c
    real, intent(out) :: d

    d = size(a) + size(b) + size(c)

end subroutine sub_3
