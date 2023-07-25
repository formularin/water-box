! Computes the diffusion coefficient of liquid water. Stores the MSD as a time series.
! Written by Arin Khare
! Options: -xtc, -ndx, -num_samples, -sample_size, -o

program diffusion
    
    use gmxfort_trajectory
    
    implicit none

    real, parameter :: dt = 0.004  ! Nanoseconds

    type(Trajectory) :: trj
    integer :: num_samples, f0_interval, sample, sample_size, f0, f, oxygen, err_status, i, delta
    integer :: count_init, count_final, count_rate, count_max
    real :: dist, time_init, time_final, elapsed_time
    real, dimension(:), allocatable :: squared_disp
    integer, dimension(:), allocatable :: sample_freqs
    character(256) :: arg, traj_file, index_file, num_samples_str, err_iomsg, output_file, sample_size_str

    call system_clock(count_init, count_rate, count_max)
    time_init = count_init * 1.0 / count_rate

    traj_file = "prd.xtc"
    index_file = "prd.ndx"
    num_samples_str = "100"
    sample_size_str = "10000"
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

    allocate(squared_disp(1:trj%nframes))
    allocate(sample_freqs(1:trj%nframes))
    squared_disp = 0.0
    sample_freqs = 0

    write(*, "(I0, A)") num_samples, " samples"
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
    write(*, "(A)") "\ndone."
    write(*, "(A, A)")  "Writing to ", output_file

    squared_disp = squared_disp / sample_freqs

    open(unit=10, file=output_file, iostat=err_status, iomsg=err_iomsg)
    if (err_status /= 0) then
        write(*, "(A, A)") "Error ", trim(err_iomsg)
        stop
    endif

    do f=1, trj%nframes
        write (unit=10, fmt="(F10.3, A, F25.17)") (f-1)*dt, "  ", squared_disp(f)
    enddo

    close(unit=10)

    write(*, "(A)") "done."

    call system_clock(count_final, count_rate, count_max)
    time_final = count_final * 1.0 / count_rate
    elapsed_time = time_final - time_init

    write (*, "(A, F0.5, A)") "Wall time: ", elapsed_time, "s"

contains
    real function get_squared_dist(r1, r2)
        real, dimension(3) :: r1, r2
        integer :: d
        get_squared_dist = 0
        do d=1, 3
            ! get_squared_dist = get_squared_dist + (dx - box_dims(d) * nint(dx / box_dims(d))) ** 2
            ! Are we supposed to not be using PBC?
            get_squared_dist = get_squared_dist + (r1(d) - r2(d)) ** 2
        enddo
    end function get_squared_dist

end program diffusion