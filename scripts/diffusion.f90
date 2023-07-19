! Computes the mean square displacement and the velocity autocorrelation function of water as time series.
! Written by Arin Khare

program diffusivity
    
    use gmxfort_trajectory
    
    implicit none
    
    type(Trajectory) :: trj
        
    call trj%read("../prd.xtc")

    

end program diffusivity