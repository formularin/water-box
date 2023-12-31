! Calculates the diffusion constant of water.
! Written by Arin Khare
! Options: -xtc, -ndx, -num_samples, -sample_size, -calc-tau, -dt, -o
! - If the calc-tau option is specified, the program will also compute the average 
!   time between collisions, assuming the trajectory file contains all frames.
! - Assumes dt is in ns during diffusivity calculation.
! - Requires unwrapped coordinates if PBC are used. This can be done using the command:
!   `gmx trjconv -f spce.xtc -pbc nojump -o spce-unwrapped.xtc`

program diffusion
    
    use gmxfort_trajectory
    use waterbox_utils
    
    implicit none

    real, parameter :: collision_radius = 0.25

    type(Trajectory) :: trj
    logical calc_tau
    integer :: num_samples, f0_interval, sample, sample_size, f0, f, oxygen
    integer :: count_init, count_final, count_rate, count_max, err_status, i, j, delta
    real :: dt, dist, time_init, time_final, elapsed_time, slope, y_int, tau
    real, dimension(:), allocatable :: squared_disp, times
    integer, dimension(:), allocatable :: sample_freqs, collisions
    character(256) :: arg, traj_file, index_file, num_samples_str, err_iomsg, output_file, sample_size_str, dt_str

    call system_clock(count_init, count_rate, count_max)
    time_init = count_init * 1.0 / count_rate

    traj_file = "simulations/water/spce-unwrapped.xtc"
    index_file = "simulations/water/spce.ndx"
    num_samples_str = "10"
    sample_size_str = "500"
    dt_str = "0.004"  ! 4ps
    output_file = "results/msd.xvg"
    calc_tau = .false.

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
        else if (arg == "-dt") then
            call get_command_argument(i + 1, dt_str)
        else if (arg == "-o") then
            call get_command_argument(i + 1, output_file)
        else if (arg == "-calc-tau") then
            calc_tau = .true.
        endif
    enddo

    call trj%read(traj_file, index_file)

    read(num_samples_str, *) num_samples
    f0_interval = trj%nframes / num_samples
    read(sample_size_str, *) sample_size
    read(dt_str, *) dt

    allocate(squared_disp(1:sample_size))
    allocate(times(1:sample_size))
    allocate(sample_freqs(1:sample_size))
    allocate(collisions(1:trj%natoms("OW")))
    squared_disp = 0.0
    sample_freqs = 0

    write(*, "(I0, A)") num_samples, " samples"
    squared_disp = 0
    sample_freqs = 0
    do sample=1, num_samples
        f0 = (sample - 1) * f0_interval + 1
        write(*, "(A, I0)", advance="no") "\rTaking sample ", sample
        do f=f0, min(trj%nframes, f0 + sample_size - 1)
            delta = (f - f0) + 1
            do oxygen=1, trj%natoms("OW")
                dist = get_magnitude2(trj%x(f, oxygen, "OW") - trj%x(f0, oxygen, "OW"))
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
    write(*, "(A, F0.7, A)") "Slope: ", slope, "nm^2/ns"
    write(*, "(A, F0.7, A)") "Y-intercept: ", y_int, "nm^2"
    write(*, "(A, F0.7, A)") "Diffusivity: ", slope / 6, "µm^2/ms"

    write(*, "(A)")  "Writing MSD timeseries to "//trim(output_file)//"..."

    open(unit=10, file=output_file, iostat=err_status, iomsg=err_iomsg)
    if (err_status /= 0) then
        write(*, "(A, A)") "Error ", trim(err_iomsg)
        stop
    endif

    do f=1, sample_size
        write (unit=10, fmt="(F25.17, A, F25.17)") (f-1)*dt, "  ", squared_disp(f)
    enddo

    close(unit=10)

    write(*, "(A)") "done."

    if (calc_tau) then
        write(*, "(A)") "Calculating tau..."
        collisions = 0
        do f=1, trj%nframes
            write(*, "(A, I0)", advance="no") "\rChecking for collisions on frame: ", f
            do i=1, trj%natoms("OW") - 1
                do j=i+1, trj%natoms("OW")
                    dist = get_magnitude2(trj%x(f, i, "OW") - trj%x(f, j, "OW"))
                    if (dist < collision_radius**2) then
                        collisions(i) = collisions(i) + 1
                        collisions(j) = collisions(j) + 1
                    endif
                enddo
            enddo
        enddo
        tau = 0
        do i=1, trj%natoms("OW")
            tau = tau + (dt * trj%nframes) / collisions(i)
        enddo
        tau = tau / trj%natoms("OW")
        write(*, "(A, F0.9)") "\ndone. tau=", tau
    endif

    call system_clock(count_final, count_rate, count_max)
    time_final = count_final * 1.0 / count_rate
    elapsed_time = time_final - time_init

    write (*, "(A, F0.5, A)") "Wall time: ", elapsed_time, "s"

end program diffusion