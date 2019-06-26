module y_originator_mod
    use mpp_domains_mod, only: domain2D
    use spec_mpp_mod,  only: spec_mpp_init, spec_mpp_end, get_grid_domain, get_spec_domain, &
                         grid_domain, spectral_domain, global_spectral_domain, atmosphere_domain
                        !use transforms_mod, only: grid_domain
    implicit none
    type(domain2D), save, public :: var
    public :: grid_domain

    contains

    subroutine rout()
    implicit none
    print *, 'hello'
        
    end subroutine


end module