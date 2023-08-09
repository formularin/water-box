! Calculates the radial distribution function of water, including the corners of the box.
! Written by Arin Khare
! Options: -xtc, -ndx, -num_bins, -mode, -o
! Mode can be either "OO", "OH," or "HH"
! Note that the max number of bins is 5000
! Outputs an XVG file containing these columns:
! - The right edge of each bin, in nm
! - The frequency of distances for each bin
! - The total volume for each bin (over all frames and oxygens)
! - The number density corresponding to each bin (unnormalized RDF)
! - The normalized RDF value for each bin
! - The standard deviation (error bar) for each bin

program rdf

    use gmxfort_trajectory
    use waterbox_utils

    implicit none

    ! File I/O, program timing, command-line args, etc.
    integer :: count_init, count_rate, count_max, count_final, err_status
    character(256) :: arg, traj_file, index_file, output_file, mode
    character(256) :: num_bins_str, err_iomsg
    real :: time_init, time_final, elapsed_time

    ! Program-specific declarations
    real, parameter :: dt = 0.004  ! Nanoseconds
    integer, parameter :: nrmax = 5000  ! This must match nrmax of spherenoncube
    real, parameter :: pi = 4.d0*atan(1.0)

    type(Trajectory) :: trj
    integer :: num_bins, num_atoms, i, f, j, k, bin, freq
    real :: dr, max_dist, diag, avg_num_density, dist
    real :: edge, vol, rdf_u, rdf_n, err, mean
    real, dimension(:), allocatable :: volumes, right_bin_edges, rdf_unnorm, rdf_norm, stdevs
    integer, dimension(:), allocatable :: frequency_sums
    integer, dimension(:, :), allocatable :: frequencies
    real, dimension(0:nrmax) :: box_volumes
    real, dimension(1:3) :: box
    real, dimension(1:4) :: distances
    integer, dimension(1:4) :: bins

    ! Start timing program
    call system_clock(count_init, count_rate, count_max)
    time_init = count_init * 1.0 / count_rate

    ! Read command-line args
    traj_file = "simulations/water/spce.xtc"
    index_file = "simulations/water/spce.ndx"
    num_bins_str = "1000"
    output_file = "results/rdf.xvg"
    mode = "OO"

    do i=0, command_argument_count()
        call get_command_argument(i, arg)
        if (arg == "-xtc") then
            call get_command_argument(i + 1, traj_file)
        else if (arg == "-ndx") then
            call get_command_argument(i + 1, index_file)
        else if (arg == "-num_bins") then
            call get_command_argument(i + 1, num_bins_str)
        else if (arg == "-mode") then
            call get_command_argument(i + 1, mode)
        else if (arg == "-o") then
            call get_command_argument(i + 1, output_file)
        endif
    enddo

    read(num_bins_str, *) num_bins
    allocate(volumes(1:num_bins))
    allocate(right_bin_edges(1:num_bins))
    allocate(rdf_norm(1:num_bins))
    allocate(rdf_unnorm(1:num_bins))
    allocate(frequency_sums(1:num_bins))
    allocate(stdevs(1:num_bins))

    ! Read XTC file (entire thing loaded into memory)
    call trj%read(traj_file, index_file)

    num_atoms = trj%natoms("OW")
    if (mode == "HH") then
        num_atoms = num_atoms * 2
    endif
    allocate(frequencies(1:num_bins, 1:num_atoms))

    write(*, "(A)") "Precomputing volumes..."
    ! Calculate the largest possible distance (diagonal of the biggest cell)
    ! Also save average number density for later use.
    max_dist = 0
    ! Initially holds total volume
    avg_num_density = 0
    do f=1, trj%nframes
        call flatten_box(trj%box(f), box)
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
        call flatten_box(trj%box(f), box)
        call spherenoncube(box, box_volumes, dr)
        do i=1, num_bins
            ! Volume is zero if inner radius is greater than max box dimension
            if (dr * (i - 1) < max(box(1), box(2), box(3))) then
                volumes(i) = volumes(i) + (box_volumes(i) - box_volumes(i - 1)) * trj%natoms("OW")
            endif
        enddo
    enddo
    write(*, "(A)") "done."

    ! Get distances across all frames and create histogram
    write(*, "(A)") "Calculating distances and binning..."
    frequencies = 0
    do f=1, trj%nframes
        write(*, "(A, I0)", advance="no") "\rOn frame ", f
        call flatten_box(trj%box(f), box)
        ! i and j are oxygens
        do i=1, trj%natoms("OW") - 1
            do j=i + 1, trj%natoms("OW")
                if (mode == "OO") then
                    dist = get_distance(trj%x(f, i, "OW"), trj%x(f, j, "OW"), box)
                    bin = int(dist / dr) + 1
                    frequencies(bin, i) = frequencies(bin, i) + 1
                    frequencies(bin, j) = frequencies(bin, j) + 1
                else if (mode == "OH") then
                    ! i*3 and i*3 - 1 are the Hydrogen atoms for oxygen i
                    ! TODO: figure out how to create groups in the index for each type of hydrogen
                    distances(1) = get_distance(trj%x(f, i, "OW"), trj%x(f, j * 3), box)
                    distances(2) = get_distance(trj%x(f, i, "OW"), trj%x(f, j * 3 - 1), box)
                    distances(3) = get_distance(trj%x(f, j, "OW"), trj%x(f, i * 3), box)
                    distances(4) = get_distance(trj%x(f, j, "OW"), trj%x(f, i * 3 - 1), box)
                    do k=1, 4
                        bin = int(distances(k) / dr) + 1
                        if (k < 3) then
                            frequencies(bin, i) = frequencies(bin, i) + 1
                        else
                            frequencies(bin, j) = frequencies(bin, j) + 1
                        endif
                    enddo
                else if (mode == "HH") then
                    distances(1) = get_distance(trj%x(f, j * 3), trj%x(f, i * 3), box)
                    distances(2) = get_distance(trj%x(f, j * 3), trj%x(f, i * 3 - 1), box)
                    distances(3) = get_distance(trj%x(f, j * 3 - 1), trj%x(f, i * 3), box)
                    distances(4) = get_distance(trj%x(f, j * 3 - 1), trj%x(f, i * 3 - 1), box)
                    do k=1, 4
                        bins(k) = int(distances(k) / dr) + 1
                    enddo
                    ! Update frequencies for each hydrogen atom
                    frequencies(bins(1), j * 2) = frequencies(bins(1), j * 2) + 1
                    frequencies(bins(1), i * 2) = frequencies(bins(1), i * 2) + 1
                    frequencies(bins(2), j * 2) = frequencies(bins(2), j * 2) + 1
                    frequencies(bins(2), i * 2 - 1) = frequencies(bins(2), i * 2 - 1) + 1
                    frequencies(bins(3), j * 2 - 1) = frequencies(bins(1), j * 2 - 1) + 1
                    frequencies(bins(3), i * 2) = frequencies(bins(1), i * 2) + 1
                    frequencies(bins(4), j * 2 - 1) = frequencies(bins(2), j * 2 - 1) + 1
                    frequencies(bins(4), i * 2 - 1) = frequencies(bins(2), i * 2 - 1) + 1
                endif
            enddo
        enddo
    enddo
    write(*, "(A)") "\ndone."

    ! Write output to file
    write(*, "(A)") "Normalizing and writing to "//trim(output_file)//"..."
    do i=1, num_bins
        frequency_sums(i) = sum(frequencies(i, :))
        mean = frequency_sums(i) / num_atoms
        stdevs(i) = sqrt(sum((frequencies(i, :) - mean) ** 2) / (num_atoms - 1))
    enddo
    rdf_unnorm = frequency_sums / volumes
    stdevs = stdevs / volumes
    if (mode == "OO") then
        rdf_norm = rdf_unnorm / avg_num_density
        stdevs = stdevs / avg_num_density
    else
        rdf_norm = rdf_unnorm / (2 * avg_num_density)
        stdevs = stdevs / (2 * avg_num_density)
    endif

    open(unit=10, file=output_file, iostat=err_status, iomsg=err_iomsg)
    if (err_status /= 0) then
        write(*, "(A, A)") "Error ", trim(err_iomsg)
        stop
    endif

    do i=1, num_bins
        edge = dr * i
        freq = frequency_sums(i)
        vol = volumes(i)
        rdf_u = rdf_unnorm(i)
        rdf_n = rdf_norm(i)
        err = stdevs(i)
        write (unit=10, fmt="(F25.17, I9, F35.17, F25.17, F25.17, F25.17)") edge, freq, vol, rdf_u, rdf_n, err
    enddo

    close(unit=10)
    write(*, "(A)") "done."

    ! Print total running time
    call system_clock(count_final, count_rate, count_max)
    time_final = count_final * 1.0 / count_rate
    elapsed_time = time_final - time_init

    write (*, "(A, F0.5, A)") "Wall time: ", elapsed_time, "s"
end program rdf