! Utility functions for Fortran MD analysis scripts.
! Written by Arin Khare

module waterbox_utils

    implicit none
    public

contains
    ! Moves values from 2-d output from libgmxfort trj%box(f) into 1-d array of length 3
    subroutine flatten_box(trj_box, f_box)
        implicit none
        real, dimension(1:3, 1:3) :: trj_box
        real, dimension(1:3) :: f_box

        f_box(1) = trj_box(1, 1)
        f_box(2) = trj_box(2, 2)
        f_box(3) = trj_box(3, 3)
    end subroutine flatten_box

    subroutine linreg(x, y, n, m, b)
        implicit none
        integer :: j, n
        real :: sum_x, sum_y, sum_xy, sum_x2, m, b
        real, dimension(n) :: x, y

        sum_x = 0
        sum_y = 0
        sum_xy = 0
        sum_x2 = 0
        do j=1, n
            sum_x = sum_x + x(j)
            sum_y = sum_y + y(j)
            sum_xy = sum_xy + x(j) * y(j)
            sum_x2 = sum_x2 + x(j) ** 2
        enddo

        m = (n*sum_xy - sum_x*sum_y) / (n*sum_x2 - sum_x**2)
        b = (sum_y - m*sum_x) / n

    end subroutine linreg

    ! Squared magnitude of a 3-d vector
    real function get_magnitude2(x)
        implicit none
        real, dimension(1:3) :: x
        integer :: d
        get_magnitude2 = 0
        do d=1, 3
            get_magnitude2 = get_magnitude2 + x(d) ** 2
        enddo
    end function get_magnitude2

    ! Magnitude of a 3-d vector
    real function get_magnitude(x)
        implicit none
        real, dimension(1:3) :: x
        get_magnitude = sqrt(get_magnitude2(x))
    end function get_magnitude

    ! Dot product of two 3-d vectors
    real function get_dot_product(a, b)
        implicit none
        real, dimension(1:3) :: a, b
        integer :: d
        get_dot_product = 0
        do d=1,3
            get_dot_product = get_dot_product + a(d) * b(d)
        enddo
    end function get_dot_product

    ! Adjusts all values in xr (displacement vector) according to the minimum image convention.
    function get_min_image(xr, box)
        implicit none
        real, dimension(1:3) :: box, get_min_image
        real, dimension(1:3) :: xr
        integer :: d
        do d=1, 3
            get_min_image(d) = xr(d) - box(d) * nint(xr(d) / box(d))
        enddo
    end function get_min_image

    ! Distance between two molecules in a box
    real function get_distance(a, b, box)
        implicit none
        real, dimension(1:3) :: a, b, box
        get_distance = get_magnitude(get_min_image(a - b, box))
    end function get_distance

end module waterbox_utils