! Calculates the fraction of nearest neighbors that are hydrogen bonded.
! Written by Arin Khare
! Options: -xtc, -ndx, -o
! Writes to an XVG file with the following data:
! - The number of water molecules in a comment
! - Two columns (length equal to number of frames in simulation):
!   - One column containing the number of pairs of waters with oxygens within 0.35nm.
!   - One column containing the number of hydrogen bonds.

program hbonds

    use gmxfort_trajectory
    use waterbox_utils

    implicit none

    ! File I/O, prograam timing, command-line args, etc.
    integer :: count_init, count_rate, count_max, count_final, err_status
    integer :: f, i, j
    character(256) :: arg, traj_file, index_file, err_iomsg, output_file
    real :: time_init, time_final, elapsed_time

    ! Program-specifc declarations
    real, parameter :: PI = 3.1415925635897932384626433
    real, parameter :: NEAREST_NEIGHBOR_DIST = 0.35
    real, parameter :: HBOND_ANGLE = PI / 6.

    type(Trajectory) :: trj
    integer, dimension(:), allocatable :: num_nn, num_hbonded
    real, dimension(1:3) :: oo_vec, oh_vec, min_oh_vec, box
    real :: oo_dist, oh_dist, min_oh_dist, dot_prod
    integer :: h_donor
    
    ! Start timing program
    call system_clock(count_init, count_rate, count_max)
    time_init = count_init * 1.0 / count_rate

    ! Read command-line args
    traj_file = "simulations/water/spce.xtc"
    index_file = "simulations/water/spce.ndx"
    output_file = "results/hbonds.xvg"

    do i=1, command_argument_count()
        call get_command_argument(i, arg)
        if (arg == "-xtc") then
            call get_command_argument(i + 1, traj_file)
        else if (arg == "-ndx") then
            call get_command_argument(i + 1, index_file)
        else if (arg == "-o") then
            call get_command_argument(i + 1, output_file)
        endif
    enddo

    ! Read XTC file (entire thing loaded into memory)
    call trj%read(traj_file, index_file)

    ! Get nearest neighbors and hydrogen bonds for all pairs all frames.
    allocate(num_nn(1:trj%nframes))
    allocate(num_hbonded(1:trj%nframes))
    num_nn = 0
    num_hbonded = 0

    do f=1, trj%nframes
        write(*, "(A, I0)", advance="no") "\rOn frame ", f
        call flatten_box(trj%box(f), box)
        do i=1, trj%natoms("OW")
            do j=i + 1, trj%natoms("OW")
                oo_dist = get_distance(trj%x(f, i, "OW"), trj%x(f, j, "OW"), box)
                if (oo_dist < NEAREST_NEIGHBOR_DIST) then
                    num_nn(f) = num_nn(f) + 1

                    ! Assume hydrogen donor is i.
                    oo_vec = get_min_image(trj%x(f, i, "OW") - trj%x(f, j, "OW"), box)

                    ! Get the smallest O-H vector from the 4 possibilities.
                    min_oh_vec = get_min_image(trj%x(f, i, "HW1") - trj%x(f, j, "OW"), box)
                    min_oh_dist = get_magnitude(min_oh_vec)
                    h_donor = i

                    oh_vec = get_min_image(trj%x(f, i, "HW2") - trj%x(f, j, "OW"), box)
                    oh_dist = get_magnitude(oh_vec)
                    if (oh_dist < min_oh_dist) then
                        min_oh_dist = oh_dist
                        min_oh_vec = oh_vec
                    endif

                    oh_vec = get_min_image(trj%x(f, j, "HW1") - trj%x(f, i, "OW"), box)
                    oh_dist = get_magnitude(oh_vec)
                    if (oh_dist < min_oh_dist) then
                        min_oh_dist = oh_dist
                        min_oh_vec = oh_vec
                        h_donor = j
                    endif

                    oh_vec = get_min_image(trj%x(f, j, "HW2") - trj%x(f, i, "OW"), box)
                    oh_dist = get_magnitude(oh_vec)
                    if (oh_dist < min_oh_dist) then
                        min_oh_dist = oh_dist
                        min_oh_vec = oh_vec
                        h_donor = j
                    endif

                    if (h_donor == j) then
                        oo_vec = -oo_vec
                    endif

                    dot_prod = get_dot_product(oh_vec, oo_vec)
                    if (acos(dot_prod / (min_oh_dist * oo_dist)) < HBOND_ANGLE) then
                        num_hbonded(f) = num_hbonded(f) + 1
                    endif
                endif
            enddo
        enddo
    enddo
    write(*, "(A)") "\ndone."
    write(*, "(A)", advance="no") " Writing to file... "

    ! Write output to file.
    open(unit=10, file=output_file, iostat=err_status, iomsg=err_iomsg)
    if (err_status /= 0) then
        write(*, "(A, A)") "Error ", trim(err_iomsg)
        stop
    endif

    write(unit=10, fmt="(A, I5)") "# num_waters=", trj%natoms("OW")
    do f=1, trj%nframes
        write(unit=10, fmt="(I6, I6)") num_nn(f), num_hbonded(f)
    enddo

    close(unit=10)
    write(*, "(A)") "done."

    ! Print total running time
    call system_clock(count_final, count_rate, count_max)
    time_final = count_final * 1.0 / count_rate
    elapsed_time = time_final - time_init

    write(*, "(A, F0.5, A)") "Wall time: ", elapsed_time, "s"
end program hbonds