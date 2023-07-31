! Calculates the O-O radial distribution function of water, including the corners of the box.
! Written by Arin Khare
! Options: -xtc, -ndx, -num_bins, -o
! Note that the max number of bins is 5000
! Outputs an XVG file containing these columns:
! - The right edge of each bin, in nm
! - The frequency of distances for each bin
! - The total volume for each bin (over all frames and oxygens)
! - The number density corresponding to each bin (unnormalized RDF)
! - The normalized RDF value for each bin

program rdf

    use gmxfort_trajectory

    implicit none

    ! File I/O, program timing, command-line args, etc.
    integer :: count_init, count_rate, count_max, count_final, err_status
    character(256) :: arg, traj_file, index_file, output_file
    character(256) :: num_bins_str, err_iomsg
    real*8 :: time_init, time_final, elapsed_time

    ! Program-specific declarations
    real*8, parameter :: dt = 0.004  ! Nanoseconds
    integer, parameter :: nrmax = 5000  ! This must match nrmax of spherenoncube
    real*8, parameter :: pi = 4.d0*atan(1.0)

    type(Trajectory) :: trj
    integer :: num_bins, i, f, j, bin, freq
    real*8 :: dr, max_dist, diag, avg_num_density
    real*8 :: edge, vol, rdf_u, rdf_n
    real*8, dimension(:), allocatable :: volumes, right_bin_edges, rdf_unnorm, rdf_norm
    integer, dimension(:), allocatable :: frequencies
    real*8, dimension(0:nrmax) :: box_volumes
    real*8, dimension(1:3) :: box, disp

    ! Start timing program
    call system_clock(count_init, count_rate, count_max)
    time_init = count_init * 1.0 / count_rate

    ! Read command-line args
    traj_file = "simulations/water/spce.xtc"
    index_file = "simulations/water/spce.ndx"
    num_bins_str = "1000"
    output_file = "results/rdf.xvg"

    do i=0, command_argument_count()
        call get_command_argument(i, arg)
        if (arg == "-xtc") then
            call get_command_argument(i + 1, traj_file)
        else if (arg == "-ndx") then
            call get_command_argument(i + 1, index_file)
        else if (arg == "-num_bins") then
            call get_command_argument(i + 1, num_bins_str)
        else if (arg == "-o") then
            call get_command_argument(i + 1, output_file)
        endif
    enddo

    read(num_bins_str, *) num_bins
    allocate(volumes(1:num_bins))
    allocate(right_bin_edges(1:num_bins))
    allocate(rdf_norm(1:num_bins))
    allocate(rdf_unnorm(1:num_bins))
    allocate(frequencies(1:num_bins))

    ! Read XTC file (entire thing loaded into memory)
    call trj%read(traj_file, index_file)

    write(*, "(A)") "Precomputing volumes..."
    ! Calculate the largest possible distance (diagonal of the biggest cell)
    ! Also save average number density for later use.
    max_dist = 0
    ! Initially holds total volume
    avg_num_density = 0
    do f=1, trj%nframes
        call get_box(trj%box(f), box)
        diag = get_magnitude(box / 2)
        if (diag > max_dist) then
            max_dist = diag
        endif
        avg_num_density = avg_num_density + box(1) * box(2) * box(3)
    enddo
    ! This ensures the binning algorithm works even if there is a distance precisely equal to the maximum.
    max_dist = max_dist + 0.0001
    avg_num_density = trj%natoms("OW") / (avg_num_density / trj%nframes)

    dr = max_dist / num_bins

    ! Precompute volumes for each bin
    volumes = 0
    do f=1, trj%nframes
        call get_box(trj%box(f), box)
        ! box = (/2.3, 2.3, 2.3/)
        call spherenoncube(box, box_volumes, dr)
        do i=1, num_bins
            ! Volume is zero if inner radius is greater than max box dimension
            if (dr * (i - 1) < max(box(1), box(2), box(3))) then
                volumes(i) = volumes(i) + (box_volumes(i) - box_volumes(i - 1)) * trj%natoms("OW")
                ! volumes(i) = volumes(i) + (4. / 3) * pi * ((dr * i) ** 3 - (dr * (i - 1)) ** 3) * trj%natoms("OW")
            endif
        enddo
    enddo
    write(*, "(A)") "done."

    ! Get distances across all frames and create histogram
    write(*, "(A)") "Calculating distances and binning..."
    frequencies = 0
    do f=1, trj%nframes
        write(*, "(A, I0)", advance="no") "\rOn frame ", f
        call get_box(trj%box(f), box)
        ! i and j are oxygens
        do i=1, trj%natoms("OW") - 1
            do j=i + 1, trj%natoms("OW")
                disp = get_min_image(trj%x(f, i, "OW") - trj%x(f, j, "OW"), box)
                bin = int(get_magnitude(disp) / dr) + 1
                frequencies(bin) = frequencies(bin) + 2
            enddo
        enddo
    enddo
    write(*, "(A)") "\ndone."

    ! Write output to file
    write(*, "(A)") "Normalizing and writing to "//trim(output_file)//"..."
    rdf_unnorm = frequencies / volumes
    rdf_norm = rdf_unnorm / avg_num_density

    open(unit=10, file=output_file, iostat=err_status, iomsg=err_iomsg)
    if (err_status /= 0) then
        write(*, "(A, A)") "Error ", trim(err_iomsg)
        stop
    endif

    do i=1, num_bins
        edge = dr * i
        freq = frequencies(i)
        ! vol = volumes(i) / (trj%nframes * trj%natoms("OW"))
        vol = volumes(i)
        rdf_u = rdf_unnorm(i)
        rdf_n = rdf_norm(i)
        write (unit=10, fmt="(F25.17, I9, F35.17, F25.17, F25.17)") edge, freq, vol, rdf_u, rdf_n
    enddo

    close(unit=10)
    write(*, "(A)") "done."

    ! Print total running time
    call system_clock(count_final, count_rate, count_max)
    time_final = count_final * 1.0 / count_rate
    elapsed_time = time_final - time_init

    write (*, "(A, F0.5, A)") "Wall time: ", elapsed_time, "s"

contains

    ! Magnitude of a 3-d vector
    real*8 function get_magnitude(x)
        implicit none
        real*8, dimension(1:3) :: x
        integer :: d
        get_magnitude = 0
        do d=1, 3
            get_magnitude = get_magnitude + x(d) ** 2
        enddo
        get_magnitude = sqrt(get_magnitude)
    end function get_magnitude

    ! Adjusts all values in xr (displacement vector) according to the minimum image convention.
    function get_min_image(xr, box_dims)
        implicit none
        real*8, dimension(1:3) :: box_dims, get_min_image
        real, dimension(1:3) :: xr
        integer :: d
        do d=1, 3
            get_min_image(d) = xr(d) - box_dims(d) * nint(xr(d) / box_dims(d))
        enddo
    end function get_min_image

    ! Moves values from 2-d output from libgmxfort trj%box(f) into 1-d array of length 3
    subroutine get_box(trj_box, f_box)
        implicit none
        real, dimension(1:3, 1:3) :: trj_box
        real*8, dimension(1:3) :: f_box

        f_box(1) = trj_box(1, 1)
        f_box(2) = trj_box(2, 2)
        f_box(3) = trj_box(3, 3)

    end subroutine get_box

end program rdf