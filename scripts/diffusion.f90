! Computes the diffusion coefficient of liquid water. Stores the MSD as a time series.
! Written by Arin Khare
! Options: -xtc, -ndx, -num_samples, , -sample_size, -o

program diffusion
    
    use gmxfort_trajectory
    
    implicit none

    real, parameter :: dt = 0.004  ! Nanoseconds

    type(Trajectory) :: trj
    integer :: num_samples, f0_interval, sample, sample_size, f0, f, oxygen, err_status, i, delta
    real :: dist
    real, dimension(:), allocatable :: squared_disp
    integer, dimension(:), allocatable :: sample_freqs
    character(256) :: arg, traj_file, index_file, num_samples_str, err_iomsg, output_file, sample_size_str

    traj_file = "prd.xtc"
    index_file = "prd.ndx"
    num_samples_str = "100"
    sample_size_str = "10000"
    output_file = "msd.xvg"

    do i=0, command_argument_count()
        call get_command_argument(i, arg)
        if (arg .eq. "-xtc") then
            call get_command_argument(i + 1, traj_file)
        else if (arg .eq. "-ndx") then
            call get_command_argument(i + 1, index_file)
        else if (arg .eq. "-num_samples") then
            call get_command_argument(i + 1, num_samples_str)
        else if (arg .eq. "-sample_size") then
            call get_command_argument(i + 1, sample_size_str)
        else if (arg .eq. "-o") then
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

    do sample=1, num_samples
        f0 = (sample - 1) * f0_interval + 1
        do f=f0, min(trj%nframes, f0 + sample_size)
            delta = (f - f0) + 1
            do oxygen=1, trj%natoms("OW")
                dist = get_squared_dist(trj%x(f, oxygen, "OW"), trj%x(f0, oxygen, "OW"), trj%box(f))
                squared_disp(delta) = squared_disp(delta) + dist
            enddo
            ! print *, sample_freqs(delta), trj%natoms("OW")
            sample_freqs(delta) = sample_freqs(delta) + trj%natoms("OW")
        enddo
    enddo

    squared_disp = squared_disp / sample_freqs

    open(unit=10, file=output_file, iostat=err_status, iomsg=err_iomsg)
    if (err_status /= 0) then
        print *, "Error ", trim(err_iomsg)
        stop
    endif

    do f=1, trj%nframes
        write (unit=10, fmt="(F10.3, A, F25.17)") (f-1)*dt, "  ", squared_disp(f)
        ! write (unit=10, fmt="(F10.3, A, I25)") (f-1)*dt, "  ", sample_freqs(f)
    enddo

    close(unit=10)

contains
    real function get_squared_dist(r1, r2, box_dims)
        real, dimension(3) :: r1, r2, box_dims
        real :: dx
        integer :: d
        get_squared_dist = 0
        do d=1, 1
            dx = abs(r1(d) - r2(d))
            get_squared_dist = get_squared_dist + (dx - box_dims(d) * nint(dx / box_dims(d))) ** 2
        enddo
    end function get_squared_dist

end program diffusion