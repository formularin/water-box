! Creates a probability distribution of the number of waters within a spherical volume at any given time.
! Written by Arin Khare
! Options: -xtc, -ndx, -r, -o
! -r is the radius of the spherical volume (nm).

program count_waters

    use gmxfort_trajectory
    use gmxfort_utils
    use waterbox_utils

    implicit none

    ! File I/O, program timing, command-line args, etc.
    integer :: count_init, count_rate, count_max, count_final, err_status
    character(256) :: arg, traj_file, index_file, output_file, err_iomsg, radius_str
    real :: time_init, time_final, elapsed_time

    ! Program-specific declarations
    type(Trajectory) :: trj
    integer :: i, f, count, max_count
    real :: radius, avg_cell_vol
    ! Index - 1 indicates number of waters
    ! Value indicates probability
    real, dimension(:), allocatable :: distribution
    real, dimension(1:3) :: box, center
    parameter (center = (/1., 1., 1./))

    ! Start timing program
    call system_clock(count_init, count_rate, count_max)
    time_init = count_init * 1.0 / count_rate

    ! Read command-line args
    traj_file = "simulations/water/spce.xtc"
    index_file = "simulations/water/spce.ndx"
    radius_str = "0.3"
    output_file = "results/count_waters.xvg"

    do i=0, command_argument_count()
        call get_command_argument(i, arg)
        if (arg == "-xtc") then
            call get_command_argument(i + 1, traj_file)
        else if (arg == "-ndx") then
            call get_command_argument(i + 1, index_file)
        else if (arg == "-r") then
            call get_command_argument(i + 1, radius_str)
        else if (arg == "-o") then
            call get_command_argument(i + 1, output_file)
        endif
    enddo

    read(radius_str, *) radius

    ! Read XTC file (entire thing loaded into memory)
    call trj%read(traj_file, index_file)

    write(*, "(A)") "Sampling distribution..."
    allocate(distribution(1:trj%natoms("OW")))
    distribution = 0
    max_count = 0
    avg_cell_vol = 0
    do f=1, trj%nframes
        write(*, "(A, I0)", advance="no") "\rOn frame ", f
        call flatten_box(trj%box(f), box)
        avg_cell_vol = avg_cell_vol + box(1) * box(2) * box(3)
        count = 0
        do i=1, trj%natoms("OW")
            if (distance(dble(center), dble(trj%x(f, i, "OW")), dble(trj%box(f))) < radius) then
                count = count + 1
            endif
        enddo
        if (count > max_count) then
            max_count = count
        endif
        distribution(count + 1) = distribution(count + 1) + 1
    enddo
    distribution = distribution / trj%nframes
    avg_cell_vol = avg_cell_vol / trj%nframes
    write(*, "(A)") "\ndone."

    write(*, "(A)") "Writing to "//trim(output_file)//"..."
    open(unit=10, file=output_file, iostat=err_status, iomsg=err_iomsg)
    if (err_status /= 0) then
        write(*, "(A, A)") "Error ", trim(err_iomsg)
        stop
    endif

    write(unit=10, fmt="(A, F0.5)") "# Average cell volume = ", avg_cell_vol
    write(unit=10, fmt="(A, I0)") "# Number of waters = ", trj%natoms("OW")
    do i=1, max_count + 1
        write(unit=10, fmt="(I9, F25.17)") i - 1, distribution(i)
    enddo

    close(unit=10)
    write(*, "(A)") "done."

    ! Print total running time
    call system_clock(count_final, count_rate, count_max)
    time_final = count_final * 1.0 / count_rate
    elapsed_time = time_final - time_init

    write (*, "(A, F0.5, A)") "Wall time: ", elapsed_time, "s"
end program count_waters