! Least squares linear regression
! Mandatory args: -f (path to 2-col XVG file), -n (number of datapoints) -nskip
! Written by Arin Khare

program linreg
    implicit none
    integer :: i, n, nskip
    real :: x, y, sum_x, sum_y, sum_xy, sum_x2, m, b
    character(256) :: fpath, arg, n_str, nskip_str

    nskip_str = "0"

    do i=1, command_argument_count()
        call get_command_argument(i, arg)
        if (arg == "-n") then
            call get_command_argument(i + 1, n_str)
        endif
        if (arg == "-f") then
            call get_command_argument(i + 1, fpath)
        endif
        if (arg == "-nskip") then
            call get_command_argument(i + 1, nskip_str)
        endif
    enddo

    read(n_str, *) n
    read(nskip_str, *) nskip
    
    open(10, file=fpath)

    sum_x = 0
    sum_y = 0
    sum_xy = 0
    sum_x2 = 0
    do i=1, nskip + n
        read(10, *) x, y
        if (i <= nskip) then
            cycle
        endif
        sum_x = sum_x + x
        sum_y = sum_y + y
        sum_xy = sum_xy + x * y
        sum_x2 = sum_x2 + x ** 2
    enddo

    m = (n*sum_xy - sum_x*sum_y) / (n*sum_x2 - sum_x**2)
    b = (sum_y - m*sum_x) / n

    write (*, "(A, F0.10)") "m: ", m
    write (*, "(A, F0.10)") "b: ", b

    close(10)

end program linreg