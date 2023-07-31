! Computes the diffusion coefficient of liquid water. Stores the MSD as a time series.
! Written by Arin Khare
! Options: -xtc, -ndx, -num_samples, -sample_size, -o
! Requires unwrapped coordinates if PBC are used.
! This can be done using the command: gmx trjconv -f prd.xtc -pbc nojump -o prd-unwrapped.xtc

program diffusion
    
    use gmxfort_trajectory
    
    implicit none

    real, parameter :: dt = 0.004  ! Nanoseconds

    type(Trajectory) :: trj
    integer :: num_samples, f0_interval, sample, sample_size, f0, f, oxygen, err_status, i, delta
    integer :: count_init, count_final, count_rate, count_max
    real :: dist, time_init, time_final, elapsed_time, slope, y_int
    real, dimension(:), allocatable :: squared_disp, times
    integer, dimension(:), allocatable :: sample_freqs
    character(256) :: arg, traj_file, index_file, num_samples_str, err_iomsg, output_file, sample_size_str

    call system_clock(count_init, count_rate, count_max)
    time_init = count_init * 1.0 / count_rate

    traj_file = "prd-unwrapped.xtc"
    index_file = "prd.ndx"
    num_samples_str = "5"
    sample_size_str = "500"
    output_file = "msd.xvg"

    do i=0, command_argument_count()
        call get_command_argument(i, arg)
        if (arg == "-xtc") then
            call get_command_argument(i + 1, traj_file)
        else if (arg == "-ndx") then
            call get_command_argument(i + 1, index_file)
        else if (arg == "-num_samples") then
            call get_command_argument(i + 1, num_samples_str)
        else if (arg == "-sample_size") then
            call get_command_argument(i + 1, sample_size_str)
        else if (arg == "-o") then
            call get_command_argument(i + 1, output_file)
        endif
    enddo

    call trj%read(traj_file, index_file)

    read(num_samples_str, *) num_samples
    f0_interval = trj%nframes / num_samples
    read(sample_size_str, *) sample_size

    allocate(squared_disp(1:sample_size))
    allocate(times(1:sample_size))
    allocate(sample_freqs(1:sample_size))
    squared_disp = 0.0
    sample_freqs = 0

    write(*, "(I0, A)") num_samples, " samples"
    squared_disp = 0
    sample_freqs = 0
    do sample=1, num_samples
        f0 = (sample - 1) * f0_interval + 1
        write(*, "(A, I0)", advance="no") "\rTaking sample ", sample
        do f=f0, min(trj%nframes, f0 + sample_size)
            delta = (f - f0) + 1
            do oxygen=1, trj%natoms("OW")
                dist = get_squared_dist(trj%x(f, oxygen, "OW"), trj%x(f0, oxygen, "OW"))
                squared_disp(delta) = squared_disp(delta) + dist
            enddo
            sample_freqs(delta) = sample_freqs(delta) + trj%natoms("OW")
        enddo
    enddo
    squared_disp = squared_disp / sample_freqs
    write(*, "(A)") "\ndone."

    write(*, "(A)") "Performing linear regression on MSD timeseries..."
    do delta=1, sample_size
        times(delta) = delta * dt
    enddo
    call linreg(times, squared_disp, sample_size, slope, y_int)
    write(*, "(A)") "done."
    write(*, "(A, F0.7, A)") "Slope: ", slope, "nm/ns"
    write(*, "(A, F0.7, A)") "Y-intercept: ", y_int, "nm"
    write(*, "(A, F0.7, A)") "Diffusivity: ", slope / 6, "µm/ms"

    write(*, "(A)")  "Writing MSD timeseries to "//trim(output_file)//"..."

    open(unit=10, file=output_file, iostat=err_status, iomsg=err_iomsg)
    if (err_status /= 0) then
        write(*, "(A, A)") "Error ", trim(err_iomsg)
        stop
    endif

    do f=1, sample_size
        write (unit=10, fmt="(F10.3, A, F25.17)") (f-1)*dt, "  ", squared_disp(f)
    enddo

    close(unit=10)

    write(*, "(A)") "done."

    call system_clock(count_final, count_rate, count_max)
    time_final = count_final * 1.0 / count_rate
    elapsed_time = time_final - time_init

    write (*, "(A, F0.5, A)") "Wall time: ", elapsed_time, "s"

contains
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

    real function get_squared_dist(r1, r2)
        real, dimension(3) :: r1, r2
        integer :: dim
        get_squared_dist = 0
        do dim=1, 3
            get_squared_dist = get_squared_dist + (r1(dim) - r2(dim)) ** 2
        enddo
    end function get_squared_dist

    ! Moves values from 2-d output from libgmxfort trj%box(f) into 1-d array of length 3
    subroutine get_box(trj_box, f_box)
        implicit none
        real, dimension(1:3, 1:3) :: trj_box
        real, dimension(1:3) :: f_box

        f_box(1) = trj_box(1, 1)
        f_box(2) = trj_box(2, 2)
        f_box(3) = trj_box(3, 3)

    end subroutine get_box

end program diffusion