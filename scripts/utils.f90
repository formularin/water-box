! Utility functions for Fortran MD analysis scripts.
! Written by Arin Khare

module waterbox_utils

    implicit none
    public

contains
    ! Moves values from 2-d output from libgmxfort trj%box(f) into 1-d array of length 3
    subroutine flatten_box(trj_box, f_box)
        implicit none
        real, dimension(1:3, 1:3) :: trj_box
        real, dimension(1:3) :: f_box

        f_box(1) = trj_box(1, 1)
        f_box(2) = trj_box(2, 2)
        f_box(3) = trj_box(3, 3)
    end subroutine flatten_box

end module waterbox_utils