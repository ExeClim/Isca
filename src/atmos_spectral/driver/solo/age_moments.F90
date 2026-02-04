module age_moments_mod
    implicit none
    public :: get_age_moments
contains

subroutine get_age_moments(nsphum, nsphum_age, previous, grid_tracers, sink, dt_tracers)
    implicit none
        
    integer,                    intent(in)    :: nsphum
    integer,                    intent(in)    :: nsphum_age
    integer,                    intent(in)    :: previous
    real, dimension(:,:,:,:,:), intent(in)    :: grid_tracers
    real, dimension(:,:,:),   intent(in)      :: sink
    real, dimension(:,:,:,:),   intent(inout)   :: dt_tracers

    real :: eps_blowup
    integer :: i

    eps_blowup = 1e-10

    ! Calculates the RHS of age-moment evolution equation
    do i = 2, nsphum_age
        dt_tracers(:,:,:,i)  = dt_tracers(:,:,:,i) +  (i-1) * grid_tracers(:,:,:,previous,i-1) + sink * (grid_tracers(:,:,:,previous,i)/(eps_blowup+grid_tracers(:,:,:,previous,nsphum)))
    end do

end subroutine get_age_moments

end module age_moments_mod

