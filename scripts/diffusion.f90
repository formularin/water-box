! Computes the diffusion coefficient of liquid water. Stores the MSD as a time series.
! Written by Arin Khare
! Options: -xtc, -ndx, -num_samples, -sample_size, -o

program diffusion
    
    use gmxfort_trajectory
    
    implicit none

    real, parameter :: dt = 0.004  ! Nanoseconds

    type(Trajectory) :: trj
    integer :: num_samples, f0_interval, sample, sample_size, f0, f, oxygen, err_status, i, delta, d
    integer :: count_init, count_final, count_rate, count_max, side
    real :: dist, time_init, time_final, elapsed_time
    real, dimension(3) :: pos, move, box, pos_0
    real, dimension(:), allocatable :: squared_disp
    real, dimension(:, :), allocatable :: unwrap_terms
    real, dimension(:, :, :), allocatable :: coords
    integer, dimension(:), allocatable :: sample_freqs
    character(256) :: arg, traj_file, index_file, num_samples_str, err_iomsg, output_file, sample_size_str

    call system_clock(count_init, count_rate, count_max)
    time_init = count_init * 1.0 / count_rate

    traj_file = "prd.xtc"
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
    allocate(sample_freqs(1:sample_size))
    allocate(unwrap_terms(1:trj%natoms("OW"), 1:3))
    allocate(coords(1:trj%nframes, 1:trj%natoms("OW"), 1:3))
    squared_disp = 0.0
    sample_freqs = 0

    ! Unwrap coordinates
    do oxygen=1, trj%natoms("OW")
        do d=1, 3
            unwrap_terms(oxygen, d) = 0
        enddo
    enddo

    do f=1, trj%nframes
        do oxygen=1, trj%natoms("OW")
            pos = trj%x(f, oxygen, "OW")
            call get_box(trj%box(f), box)
            if (f == 1) then
                move = (/0, 0, 0/)
            else
                move = pos - trj%x(f - 1, oxygen, "OW")
            endif
            do d=1, 3
                if (move(d) > (box(d) / 2)) then
                    side = 1
                    if (move(d) > 0) then
                        side = -1
                    endif
                    unwrap_terms(oxygen, d) = unwrap_terms(oxygen, d) + box(d) * side
                endif
                coords(f, oxygen, d) = pos(d) + unwrap_terms(oxygen, d)
            enddo
        enddo
    enddo

    write(*, "(I0, A)") num_samples, " samples"
    squared_disp = 0
    sample_freqs = 0
    do sample=1, num_samples
        f0 = (sample - 1) * f0_interval + 1
        write(*, "(A, I0)", advance="no") "\rTaking sample ", sample
        do f=f0, min(trj%nframes, f0 + sample_size)
            delta = (f - f0) + 1
            do oxygen=1, trj%natoms("OW")
                pos_0 = coords(f0, oxygen, 1:3)
                pos = coords(f, oxygen, 1:3)
                dist = get_squared_dist(pos, pos_0)
                squared_disp(delta) = squared_disp(delta) + dist
            enddo
            sample_freqs(delta) = sample_freqs(delta) + trj%natoms("OW")
        enddo
    enddo
    write(*, "(A)") "\ndone."
    write(*, "(A)")  "Writing to "//trim(output_file)

    squared_disp = squared_disp / sample_freqs

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