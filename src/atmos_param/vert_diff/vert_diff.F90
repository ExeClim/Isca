!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vert_diff_mod

!=======================================================================
!
!                         VERTICAL DIFFUSION MODULE
!
!      Routines for computing the tendencies due to vertical diffusion
!
!=======================================================================

use   constants_mod, only:  GRAV, RDGAS, RVGAS, CP_AIR

use         fms_mod, only:  error_mesg, FATAL, uppercase, &
                            write_version_number, stdlog, &
                            mpp_pe, mpp_root_pe 

use   field_manager_mod, only: MODEL_ATMOS, MODEL_LAND, MODEL_ICE
use  tracer_manager_mod, only: query_method, get_number_tracers, &
     get_tracer_index, get_tracer_names, NO_TRACER

use mpp_mod, only: mpp_chksum, stdout

implicit none
private


! public interfaces
!=======================================================================
public :: vert_diff_init,          &
          vert_diff_end,           &
          gcm_vert_diff,               &
          gcm_vert_diff_down,          &
          gcm_vert_diff_up,            &
          vert_diff,                   &
          surf_diff_type

!=======================================================================

! form of interfaces
!=======================================================================


type surf_diff_type

  real, pointer, dimension(:,:) :: dtmass  => NULL(),   &
                                   dflux_t => NULL(),   &
                                   delta_t => NULL(),   &
                                   delta_u => NULL(),   &
                                   delta_v => NULL(), &
                                   sst_miz => NULL()
  real, pointer, dimension(:,:,:) :: dflux_tr => NULL(),& ! tracer flux tendency
                                     delta_tr => NULL()   ! tracer tendency
end type surf_diff_type


real,    allocatable, dimension(:,:,:) :: e_global, f_t_global, f_q_global 

! storage compartment for tracer vert. diffusion options, and for f
! coefficient if necessary
type :: tracer_data_type
   real, pointer :: f(:,:,:) => NULL() ! f coefficient field
   logical :: do_vert_diff
   logical :: do_surf_exch
end type tracer_data_type
! tracer diffusion options and storage for f coefficients
type(tracer_data_type), allocatable :: tracers(:)

      
logical :: do_conserve_energy = .true.
logical :: use_virtual_temp_vert_diff, do_mcm_plev
integer :: sphum, mix_rat

!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: vert_diff.F90,v 19.0 2012/01/06 20:27:29 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201211 $'
logical            :: module_is_initialized = .false.

real :: d608 = 0.

contains

!#######################################################################

subroutine vert_diff_init (Tri_surf, idim, jdim, kdim,    &
                               do_conserve_energy_in,         &
                               use_virtual_temp_vert_diff_in, &
                               do_mcm_plev_in )

 type(surf_diff_type), intent(inout) :: Tri_surf
 integer,              intent(in)    :: idim, jdim, kdim
 logical,              intent(in)    :: do_conserve_energy_in
 logical, optional,    intent(in)    :: use_virtual_temp_vert_diff_in
 logical, optional,    intent(in)    :: do_mcm_plev_in

 integer :: ntprog ! number of prognostic tracers in the atmosphere
 character(len=32)  :: tr_name ! tracer name
 character(len=128) :: scheme  ! tracer diffusion scheme
 integer :: n, logunit

    call write_version_number ( version, tagname )

    !s initialise constants here as rdgas no longer a parameter.
    d608 = (RVGAS-RDGAS)/RDGAS

! get the number of prognostic tracers
    call get_number_tracers( MODEL_ATMOS, num_prog=ntprog)

! get the tracer number for specific humidity
    sphum = get_tracer_index( MODEL_ATMOS, 'sphum')
    mix_rat=get_tracer_index( MODEL_ATMOS, 'mix_rat')
    if(sphum /= NO_TRACER .and. mix_rat /= NO_TRACER) then
      call error_mesg('gcm_vert_diff_init','sphum and mix_rat cannot both'// &
                      'be present in the field_table at the same time', FATAL)
    endif

    logunit = stdlog()
    if (mpp_pe() == mpp_root_pe()) then
      write (logunit,'(a,i12)') 'Tracer number for specific humidity =',sphum
      write (logunit,'(a,i12)') 'Tracer number for mixing ratio      =',mix_rat
    endif

    if(sphum==NO_TRACER) sphum=mix_rat

    if(present(use_virtual_temp_vert_diff_in)) then
      use_virtual_temp_vert_diff = use_virtual_temp_vert_diff_in
    else
      use_virtual_temp_vert_diff = .false.
    endif
    if(present(do_mcm_plev_in)) then
      do_mcm_plev = do_mcm_plev_in
    else
      do_mcm_plev = .false.
    endif

 if (.not. module_is_initialized) then

    if (allocated(  e_global ))    deallocate (  e_global )
    if (allocated(f_t_global ))    deallocate (f_t_global )
    if (allocated(f_q_global ))    deallocate (f_q_global )

    allocate(  e_global (idim, jdim, kdim-1)) ;   e_global = 0.0
    allocate(f_t_global (idim, jdim, kdim-1)) ; f_t_global = 0.0 
    allocate(f_q_global (idim, jdim, kdim-1)) ; f_q_global = 0.0

    module_is_initialized = .true.

 endif

 call alloc_surf_diff_type ( Tri_surf, idim, jdim, ntprog )
 
 do_conserve_energy = do_conserve_energy_in

 ! allocate data storage for tracers
 allocate ( tracers(ntprog) )
 do n = 1,ntprog
    ! skip tracers diffusion if it is turned off in the field table
    tracers(n)%do_vert_diff = .true. 
    if (query_method('diff_vert',MODEL_ATMOS,n,scheme)) then
       tracers(n)%do_vert_diff = (uppercase(scheme) /= 'NONE')
    endif
    ! do not exchange tracer with surface if it is not present in either land or
    ! ice model
    if (n==sphum) then
       tracers(n)%do_vert_diff = .false.
       tracers(n)%do_surf_exch = .false.
    else
       call get_tracer_names ( MODEL_ATMOS, n, tr_name )
       tracers(n)%do_surf_exch = &
            get_tracer_index ( MODEL_LAND, tr_name ) /= NO_TRACER .or.&
            get_tracer_index ( MODEL_ICE,  tr_name ) /= NO_TRACER
    endif
    ! if tracer goes through surface flux, allocate memory to hold f
    ! between downward and upward sweeps
    if(tracers(n)%do_surf_exch)&
         allocate(tracers(n)%f(idim,jdim,kdim-1))
 enddo

 write(logunit,*)'Tracer vertical diffusion properties:'
 do n = 1,ntprog
    call get_tracer_names(MODEL_ATMOS, n, tr_name)
    write(logunit,100)tr_name,tracers(n)%do_vert_diff,tracers(n)%do_surf_exch
 enddo
100 FORMAT('Tracer :',a32,': do_tr_vert_diff=',L1,' : do_tr_surf_exch=',L1)

end subroutine vert_diff_init

!#######################################################################

subroutine alloc_surf_diff_type ( Tri_surf, idim, jdim, ntprog )

type(surf_diff_type), intent(inout) :: Tri_surf
integer,              intent(in)    :: idim, jdim, ntprog

    allocate( Tri_surf%dtmass    (idim, jdim) ) ; Tri_surf%dtmass  = 0.0
    allocate( Tri_surf%dflux_t   (idim, jdim) ) ; Tri_surf%dflux_t = 0.0
    allocate( Tri_surf%delta_t   (idim, jdim) ) ; Tri_surf%delta_t = 0.0
    allocate( Tri_surf%delta_u   (idim, jdim) ) ; Tri_surf%delta_u = 0.0
    allocate( Tri_surf%delta_v   (idim, jdim) ) ; Tri_surf%delta_v = 0.0
    allocate( Tri_surf%sst_miz   (idim, jdim) ) ; Tri_surf%sst_miz = 280.0 !miz
    allocate( Tri_surf%dflux_tr  (idim, jdim, ntprog) ) ; Tri_surf%dflux_tr = 0.0
    allocate( Tri_surf%delta_tr  (idim, jdim, ntprog) ) ; Tri_surf%delta_tr = 0.0

end subroutine alloc_surf_diff_type

!#######################################################################

subroutine dealloc_surf_diff_type ( Tri_surf )

type(surf_diff_type), intent(inout) :: Tri_surf

      deallocate( Tri_surf%dtmass    )
      deallocate( Tri_surf%dflux_t   )
      deallocate( Tri_surf%delta_t   )
      deallocate( Tri_surf%delta_u   )
      deallocate( Tri_surf%delta_v   )
      deallocate( Tri_surf%sst_miz   )!miz
      deallocate( Tri_surf%dflux_tr  )
      deallocate( Tri_surf%delta_tr  )

end subroutine dealloc_surf_diff_type

!#######################################################################

subroutine vert_diff_end

  integer :: n

  if (module_is_initialized) then

    if (allocated(   e_global ))    deallocate (   e_global)
    if (allocated( f_t_global ))    deallocate ( f_t_global)
    if (allocated( f_q_global ))    deallocate ( f_q_global)

    if(allocated(tracers)) then
       do n = 1,size(tracers(:))
          if ( associated(tracers(n)%f) ) deallocate(tracers(n)%f)
       enddo
       deallocate(tracers)
    endif
  endif
  module_is_initialized = .false.


end subroutine vert_diff_end

!#######################################################################

subroutine gcm_vert_diff_down (is, js, delt,                &
                          u, v, t, q, tr,                   &
                          diff_m, diff_t, p_half, p_full,   &
                          z_full, tau_u, tau_v,             &
                          dtau_du, dtau_dv,                 &
                          dt_u, dt_v, dt_t, dt_q, dt_tr,    &
                          dissipative_heat, Tri_surf,       &
                          kbot                              )

integer, intent(in)                        :: is, js
real,    intent(in)                        :: delt
real,    intent(in)   , dimension(:,:,:)   :: u, v, t, q,     &
                                              diff_m, diff_t, &
                                              p_half, p_full, &
                                              z_full
real,    intent(in)   , dimension(:,:,:,:) :: tr
real,    intent(in)   , dimension(:,:)     :: dtau_du, dtau_dv
real,    intent(inout), dimension(:,:)     :: tau_u, tau_v
real,    intent(inout), dimension(:,:,:)   :: dt_u, dt_v, dt_t
real,    intent(in),    dimension(:,:,:)   :: dt_q
real,    intent(inout), dimension(:,:,:,:) :: dt_tr
real,    intent(out)  , dimension(:,:,:)   :: dissipative_heat
type(surf_diff_type), intent(inout)        :: Tri_surf

integer, intent(in)   , dimension(:,:), optional :: kbot

! ---- local vars
real, dimension(size(u,1),size(u,2),size(u,3)) :: &
     tt, mu, nu, e, a, b, c, g, f_tr
real, dimension(size(u,1),size(u,2)) ::           &
     f_t_delt_n1, f_q_delt_n1, f_tr_delt_n1, flux_tr, dflux_dtr, &
     mu_delt_n, nu_n, e_n1, delta_t_n, delta_q_n, delta_tr_n, &
            delta_u_n, delta_v_n
real    :: gcp
integer :: i, j, n, kb, ie, je, ntr, nlev

!-----------------------------------------------------------------------

  if(.not. module_is_initialized) call error_mesg ('gcm_vert_diff_down in vert_diff_mod',  &
      'the initialization routine gcm_vert_diff_init has not been called', &
       FATAL)
    
 ie = is + size(t,1) -1
 je = js + size(t,2) -1
 ntr  = size(tr,4)
 nlev = size(mu,3)
 
 gcp       = GRAV/CP_AIR
 tt  = t + z_full*gcp   ! the vertical gradient of tt determines the
                        ! diffusive flux of temperature

 call compute_mu (p_half, mu)
 call compute_nu (diff_m, p_half, p_full, z_full, t, q, nu) 

!  diffuse u-momentum and v_momentum
 call uv_vert_diff (delt, mu, nu, u, v, dtau_du, dtau_dv, tau_u, tau_v,  &
                    dt_u, dt_v, dt_t, delta_u_n, delta_v_n,         &
                    dissipative_heat, kbot)
                            
!  recompute nu for a different diffusivity
 call compute_nu   (diff_t, p_half, p_full, z_full, t, q, nu)

 ! calculate e, the same for all tracers since their diffusivities are 
 ! the same, and mu_delt_n, nu_n, e_n1
 call compute_e (delt, mu, nu, e, a, b, c, g)
 do j = 1,size(mu,2)
 do i = 1,size(mu,1)
    kb = nlev ; if(present(kbot)) kb=kbot(i,j)
    mu_delt_n(i,j) = mu(i,j,kb  )*delt
         nu_n(i,j) = nu(i,j,kb  )
         e_n1(i,j) = e (i,j,kb-1)
 enddo
 enddo

 do n = 1,ntr
    ! calculate f_tr, f_tr_delt_n1, delta_tr_n for this tracer
    if(.not.tracers(n)%do_vert_diff) cycle ! skip non-diffusive tracers
    call explicit_tend (mu, nu, tr(:,:,:,n), dt_tr(:,:,:,n))
    call compute_f (dt_tr(:,:,:,n), b, c, g, f_tr)
    do j = 1,size(mu,2)
    do i = 1,size(mu,1)
       kb = nlev ; if(present(kbot)) kb=kbot(i,j)
       f_tr_delt_n1(i,j) = f_tr (i,j,kb-1)*delt
       delta_tr_n(i,j)   = dt_tr(i,j,kb,n)*delt
    enddo
    enddo

    ! store information needed by flux_exchange module
    Tri_surf%delta_tr(is:ie,js:je,n) = &
         delta_tr_n(:,:) + mu_delt_n(:,:)*nu_n(:,:)*f_tr_delt_n1(:,:)
    Tri_surf%dflux_tr(is:ie,js:je,n) = -nu_n(:,:)*(1.0 - e_n1(:,:))

    if(tracers(n)%do_surf_exch) then
       ! store f for future use on upward sweep
       tracers(n)%f(is:ie,js:je,:) = f_tr(:,:,:)
    else
       ! upward sweep of tridaigonal solver for tracers that do not exchange 
       ! with surface
       flux_tr  (:,:) = 0.0 ! surface flux of tracer
       dflux_dtr(:,:) = 0.0 ! d(sfc flux)/d(tr atm)
       call diff_surface ( &
            mu_delt_n(:,:), nu_n(:,:), e_n1(:,:), f_tr_delt_n1(:,:), &
            dflux_dtr(:,:), flux_tr(:,:), 1.0, delta_tr_n(:,:) )
       call vert_diff_up ( &
            delt, e(:,:,:), f_tr(:,:,:), delta_tr_n(:,:), dt_tr(:,:,:,n), &
            kbot )
    endif
 enddo

! NOTE: actually e used in the tracer calculations above, and e_global
! calculated in the vert_diff_down_2 below are the same, since they only
! depend on mu and nu.

!  downward sweep of tridiagonal solver for temperature and specific humidity
 call vert_diff_down_2                            & 
         (delt, mu, nu, tt, q, dt_t, dt_q,        &  
         e_global             (is:ie,js:je,:),    &
         f_t_global           (is:ie,js:je,:),    &
         f_q_global           (is:ie,js:je,:),    &
         mu_delt_n, nu_n, e_n1, f_t_delt_n1, f_q_delt_n1, &
         delta_t_n, delta_q_n, kbot)

! store information needed by flux_exchange module

    Tri_surf%delta_t (is:ie,js:je) = delta_t_n + mu_delt_n*nu_n*f_t_delt_n1
    Tri_surf%dflux_t (is:ie,js:je) = -nu_n*(1.0 - e_n1)
    if (sphum/=NO_TRACER) then
       Tri_surf%delta_tr (is:ie,js:je,sphum) = delta_q_n + mu_delt_n*nu_n*f_q_delt_n1
       Tri_surf%dflux_tr (is:ie,js:je,sphum) = -nu_n*(1.0 - e_n1)
    endif
    Tri_surf%dtmass  (is:ie,js:je) = mu_delt_n
    Tri_surf%delta_u (is:ie,js:je) = delta_u_n
    Tri_surf%delta_v (is:ie,js:je) = delta_v_n

!-----------------------------------------------------------------------

end subroutine gcm_vert_diff_down

!#######################################################################

subroutine gcm_vert_diff_up (is, js, delt, Tri_surf, dt_t, dt_q, dt_tr, kbot)

integer, intent(in)                      :: is, js
real,    intent(in)                      :: delt
type(surf_diff_type), intent(in)         :: Tri_surf
real,    intent(out),   dimension(:,:,:) :: dt_t, dt_q
real,    intent(out),   dimension(:,:,:,:) :: dt_tr
integer, intent(in),    dimension(:,:), optional :: kbot

 ! ---- local vars
 integer :: ie, je, n
 real    :: surf_delta_q(size(dt_t,1),size(dt_t,2))

 ie = is + size(dt_t,1) -1
 je = js + size(dt_t,2) -1

! outunit = stdout()
!checksums! write(outunit,'("CHECKSUM::",A32," = ",Z20)')'e_global',mpp_chksum(e_global(is:ie,js:je,:))
!checksums! write(outunit,'("CHECKSUM::",A32," = ",Z20)')'f_t_global',mpp_chksum(f_t_global(is:ie,js:je,:))
!checksums! write(outunit,'("CHECKSUM::",A32," = ",Z20)')'Tri_surf%deta_t',mpp_chksum(Tri_surf%delta_t(is:ie,js:je))

 call vert_diff_up (delt ,                              &
                    e_global          (is:ie,js:je,:) , &
                    f_t_global        (is:ie,js:je,:) , &
                    Tri_surf%delta_t  (is:ie,js:je) ,   &
                    dt_t, kbot )

!checksums! write(outunit,'("CHECKSUM::",A32," = ",Z20)')'dt_t',mpp_chksum(dt_t)

 if(sphum/=NO_TRACER) then
    surf_delta_q = Tri_surf%delta_tr (is:ie,js:je,sphum)
 else
    surf_delta_q = 0.0
 endif

!checksums! write(outunit,'("CHECKSUM::",A32," = ",Z20)')'surf_delta_q',mpp_chksum(surf_delta_q)
!checksums! write(outunit,'("CHECKSUM::",A32," = ",Z20)')'f_q_global',mpp_chksum(f_q_global(is:ie,js:je,:))

 call vert_diff_up (delt ,                              &
                    e_global          (is:ie,js:je,:) , &
                    f_q_global        (is:ie,js:je,:) , &
                    surf_delta_q ,                      &
                    dt_q, kbot )

!checksums! write(outunit,'("CHECKSUM::",A32," = ",Z20)')'dt_q',mpp_chksum(dt_q)

 do n = 1,size(dt_tr,4)
    ! skip tracers if diffusion scheme turned off
    if (tracers(n)%do_vert_diff.and.tracers(n)%do_surf_exch) then
       call vert_diff_up (delt ,                           &
                    e_global           (is:ie,js:je,:) ,   &
                    tracers(n)%f       (is:ie,js:je,:) ,   &
                    Tri_surf%delta_tr  (is:ie,js:je,n) ,   &
                    dt_tr(:,:,:,n), kbot )
    endif
 enddo

end subroutine gcm_vert_diff_up

!#######################################################################

subroutine gcm_vert_diff (delt, u, v, t, q, tr,                    &
                          diff_m, diff_t, p_half, p_full, z_full,  &
                          dtau_du, dtau_dv, dsens_datmos, devap_datmos, &
                          sens, evap, tau_u, tau_v,                &
                          dt_u, dt_v, dt_t, dt_q, dt_tr,           &
                          dissipative_heat, kbot      )

!  one-step diffusion call for gcm in which there is no implicit dependence of 
!    surface fluxes on surface temperature

real,    intent(in)                          :: delt
real,    intent(in)   , dimension(:,:,:)     :: u, v, t, q, p_half, p_full, &
                                                z_full, diff_m, diff_t
real,    intent(in)   , dimension(:,:,:,:)   :: tr
real,    intent(in)   , dimension(:,:)       :: dtau_du, dtau_dv, dsens_datmos, &
                                                devap_datmos
real,    intent(inout), dimension(:,:)       :: tau_u, tau_v, sens, evap
real,    intent(inout), dimension(:,:,:)     :: dt_u, dt_v, dt_t, dt_q
real,    intent(inout), dimension(:,:,:,:)   :: dt_tr
real,    intent(out)  , dimension(:,:,:)     :: dissipative_heat

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(u,1),size(u,2),size(u,3)) :: mu, nu
real, dimension(size(u,1),size(u,2))           :: delta_u_n, delta_v_n


!-----------------------------------------------------------------------

 call compute_mu (p_half, mu)

 call compute_nu (diff_m, p_half, p_full, z_full, t, q, nu) 
 
 call uv_vert_diff (delt, mu, nu, u, v, dtau_du, dtau_dv, tau_u, tau_v, &
                    dt_u, dt_v, dt_t, delta_u_n, delta_v_n,        &
                    dissipative_heat, kbot)
                    
 call compute_nu   (diff_t, p_half, p_full, z_full, t, q, nu)

 call tq_vert_diff (delt, mu, nu, t, q, z_full,  &
                    dsens_datmos, devap_datmos,  &
                    sens, evap, dt_t, dt_q, kbot )

 call tr_vert_diff (delt, mu, nu, tr, dt_tr, kbot )

end subroutine gcm_vert_diff

!#######################################################################

subroutine vert_diff (delt, xi, t, q, diff, p_half, p_full, z_full, &
                      flux, dflux_datmos, factor, dt_xi, kbot)

! one-step diffusion of a single field 

real,    intent(in)                          :: delt
real,    intent(in)   , dimension(:,:,:)     :: xi, t, q, diff, p_half, p_full, z_full
real,    intent(inout), dimension(:,:)       :: flux
real,    intent(in)   , dimension(:,:)       :: dflux_datmos
real,    intent(in)                          :: factor
real,    intent(inout), dimension(:,:,:)     :: dt_xi

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(xi,1),size(xi,2),size(xi,3)  ) :: mu, nu
real, dimension(size(xi,1),size(xi,2),size(xi,3)-1) :: e, f

real, dimension(size(xi,1),size(xi,2))  :: mu_delt_n, nu_n, e_n1,  &
                                           f_delt_n1, delta_xi_n

!-----------------------------------------------------------------------

 call compute_mu    (p_half, mu)

 call compute_nu    (diff, p_half, p_full, z_full, t, q, nu) 

 call vert_diff_down &
     (delt, mu, nu, xi, dt_xi, e, f, mu_delt_n, nu_n, e_n1,  &
      f_delt_n1, delta_xi_n, kbot)

 call diff_surface (mu_delt_n, nu_n, e_n1, f_delt_n1,     &
                    dflux_datmos, flux, factor, delta_xi_n)

 call vert_diff_up (delt, e, f, delta_xi_n, dt_xi, kbot)

end subroutine vert_diff


!#######################################################################

subroutine uv_vert_diff (delt, mu, nu, u, v,  &
                         dtau_du, dtau_dv, tau_u, tau_v, dt_u, dt_v, dt_t, &
                          delta_u_n, delta_v_n, dissipative_heat, kbot )

real,    intent(in)                        :: delt
real,    intent(in)   , dimension(:,:,:)   :: u, v, mu, nu
real,    intent(in)   , dimension(:,:)     :: dtau_du, dtau_dv
real,    intent(inout), dimension(:,:)     :: tau_u, tau_v
real,    intent(inout), dimension(:,:,:)   :: dt_u, dt_v, dt_t
real,    intent(out)  , dimension(:,:,:)   :: dissipative_heat
real,    intent(out)  , dimension(:,:)     :: delta_u_n, delta_v_n

! Note (IH) 
!   delta_u_n = dt_u/delt at lowest model level, and similarly
!   for delta_v_n  -- it is convenient to output them separately

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(u,1),size(u,2)) :: mu_delt_n, nu_n, e_n1,    &
                                        f_u_delt_n1, f_v_delt_n1
                                        
real, dimension(size(u,1),size(u,2),size(u,3)) :: dt_u_temp, dt_v_temp

real, dimension(size(u,1),size(u,2),size(u,3)-1) :: e, f_u, f_v

real    :: half_delt, cp_inv


!-----------------------------------------------------------------------

 half_delt = 0.5*delt
 cp_inv    = 1.0/CP_AIR
 
 if (do_conserve_energy) then
   dt_u_temp = dt_u
   dt_v_temp = dt_v
 endif
 
 call vert_diff_down_2 &
     (delt, mu, nu, u, v, dt_u, dt_v, e, f_u, f_v, &
      mu_delt_n, nu_n, e_n1, f_u_delt_n1, f_v_delt_n1,  &
      delta_u_n, delta_v_n, kbot)        

 call diff_surface (mu_delt_n, nu_n, e_n1, f_u_delt_n1, &
                    dtau_du, tau_u, 1.0, delta_u_n)
 call diff_surface (mu_delt_n, nu_n, e_n1, f_v_delt_n1, &
                    dtau_dv, tau_v, 1.0, delta_v_n)

 call vert_diff_up (delt, e, f_u, delta_u_n, dt_u, kbot)
 call vert_diff_up (delt, e, f_v, delta_v_n, dt_v, kbot)

 if (do_conserve_energy) then
    dt_u_temp = dt_u - dt_u_temp
    dt_v_temp = dt_v - dt_v_temp
    dissipative_heat = - cp_inv*( (u + half_delt*dt_u_temp)*dt_u_temp &
                                 +(v + half_delt*dt_v_temp)*dt_v_temp )
    dt_t = dt_t + dissipative_heat
 else
    dissipative_heat = 0.0
 endif

!-----------------------------------------------------------------------

end subroutine uv_vert_diff

!#######################################################################

subroutine tq_vert_diff (delt, mu, nu, t, q,  z_full, &
                         dsens_datmos, devap_datmos, sens, evap, &
                         dt_t, dt_q, kbot)
                         

real,    intent(in)                        :: delt
real,    intent(in)   , dimension(:,:,:)   :: t, q, z_full, mu, nu
real,    intent(in)   , dimension(:,:)     :: dsens_datmos, devap_datmos
real,    intent(inout), dimension(:,:)     :: sens, evap
real,    intent(inout), dimension(:,:,:)   :: dt_t, dt_q

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(t,1),size(t,2)) :: mu_delt_n, nu_n,          &
                                        e_n1, f_t_delt_n1, f_q_delt_n1, &
                                        delta_t_n, delta_q_n

real, dimension(size(t,1),size(t,2),size(t,3)-1) :: e, f_t, f_q
real, dimension(size(t,1),size(t,2),size(t,3)  ) :: tt

real    :: gcp
!-----------------------------------------------------------------------

 gcp = GRAV/CP_AIR
 tt  = t + z_full*gcp
  
 call vert_diff_down_2 &
     (delt, mu, nu, tt, q, dt_t, dt_q, e, f_t, f_q,    &
      mu_delt_n, nu_n, e_n1, f_t_delt_n1, f_q_delt_n1, &
      delta_t_n, delta_q_n, kbot)


 call diff_surface (mu_delt_n, nu_n, e_n1, f_t_delt_n1,  &
                    dsens_datmos, sens, CP_AIR, delta_t_n)

 call diff_surface (mu_delt_n, nu_n, e_n1, f_q_delt_n1,  &
                    devap_datmos, evap, 1.0, delta_q_n)

 call vert_diff_up (delt, e, f_t, delta_t_n, dt_t, kbot)
 call vert_diff_up (delt, e, f_q, delta_q_n, dt_q, kbot)


!-----------------------------------------------------------------------

end subroutine tq_vert_diff

!#######################################################################

subroutine tr_vert_diff (delt, mu, nu, tr, dt_tr, kbot )

real,    intent(in)                        :: delt
real,    intent(in)   , dimension(:,:,:)   :: mu, nu
real,    intent(in)   , dimension(:,:,:,:) :: tr
real,    intent(inout), dimension(:,:,:,:) :: dt_tr

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(tr,1),size(tr,2)) :: mu_delt_n, nu_n, e_n1

real, dimension(size(tr,1),size(tr,2)) :: f_delt_n1, delta_tr_n
real, dimension(size(tr,1),size(tr,2)) :: dflux_dtr, flux
real, dimension(size(tr,1),size(tr,2),size(tr,3)-1) :: ftr
real, dimension(size(tr,1),size(tr,2),size(tr,3)-1) :: etr
real, dimension(size(tr,1),size(tr,2),size(tr,3)) :: a, b, c, g
integer :: i, j, kb, n, ntr, nlev
character(len=128) :: scheme
!-----------------------------------------------------------------------

 ntr  = size(tr,4) ! number of prognostic tracers

 dflux_dtr = 0.0
 call compute_e (delt, mu, nu, etr, a, b, c, g)
 if (present(kbot)) then
   do j=1,size(tr,2)
   do i=1,size(tr,1)
      kb = kbot(i,j)
      mu_delt_n(i,j) =  mu(i,j,kb  )*delt
           nu_n(i,j) =  nu(i,j,kb  )
           e_n1(i,j) = etr(i,j,kb-1)
   enddo
   enddo
 else
   nlev = size(mu,3)
   mu_delt_n(:,:) =  mu(:,:,nlev  )*delt
        nu_n(:,:) =  nu(:,:,nlev  )
        e_n1(:,:) = etr(:,:,nlev-1)
 endif
  
 do n = 1, ntr
   if ( n == sphum .or. n == mix_rat) cycle
   if (query_method('diff_vert',MODEL_ATMOS,n,scheme)) then
     if(uppercase(trim(scheme)) == 'NONE') cycle
   endif
   call explicit_tend (mu, nu, tr(:,:,:,n), dt_tr(:,:,:,n))
   call compute_f (dt_tr(:,:,:,n), b, c, g, ftr)
   if (present(kbot)) then
     do j=1,size(tr,2)
     do i=1,size(tr,1)
       kb = kbot(i,j)
       f_delt_n1(i,j)  =   ftr(i,j,kb-1)*delt
       delta_tr_n(i,j) = dt_tr(i,j,kb,n)*delt
     enddo
     enddo
   else
      f_delt_n1(:,:) =   ftr(:,:,nlev-1)*delt
     delta_tr_n(:,:) = dt_tr(:,:,nlev  ,n)*delt
   endif
   flux = 0.0
   call diff_surface (mu_delt_n, nu_n, e_n1, f_delt_n1, dflux_dtr, flux, 1.0, delta_tr_n)

! If flux needs to be saved then it should be made a module variable.
! vert_diff_init must allocate it and then call assign_tracer_field
! to set a pointer in tracer_manager_mod. It can be allocated as a
! 3 dimensional array with the 3'd index for tracer number.

   call vert_diff_up (delt, etr, ftr, delta_tr_n, dt_tr(:,:,:,n), kbot)
 end do

!-----------------------------------------------------------------------

end subroutine tr_vert_diff

!#######################################################################

subroutine vert_diff_down &
      (delt, mu, nu, tr, dt_tr, e, f, mu_delt_n, nu_n,  &
       e_n1, f_delt_n1, delta_tr_n, kbot)

!-----------------------------------------------------------------------

real,    intent(in)                         :: delt
real,    intent(in)    , dimension(:,:,:)   :: mu, nu
real,    intent(in)    , dimension(:,:,:)   :: tr
real,    intent(inout) , dimension(:,:,:)   :: dt_tr
real,    intent(out)   , dimension(:,:,:)   :: e
real,    intent(out)   , dimension(:,:,:)   :: f
real,    intent(out)   , dimension(:,:)     :: mu_delt_n, nu_n, e_n1
real,    intent(out)   , dimension(:,:)     :: f_delt_n1, delta_tr_n

integer, intent(in),    dimension(:,:), optional :: kbot

real, dimension(size(tr,1),size(tr,2),size(tr,3)) :: a, b, c, g

integer :: i, j, kb, nlev

!-----------------------------------------------------------------------

 call explicit_tend (mu, nu, tr, dt_tr)

 call compute_e  (delt, mu, nu, e, a, b, c, g)

 call compute_f (dt_tr, b, c, g, f)


 if (present(kbot)) then
    do j=1,size(tr,2)
    do i=1,size(tr,1)
        kb = kbot(i,j)
        mu_delt_n(i,j) =  mu(i,j,kb  )*delt
             nu_n(i,j) =  nu(i,j,kb  )
             e_n1(i,j) =   e(i,j,kb-1)
    enddo
    enddo
    do j=1,size(tr,2)
    do i=1,size(tr,1)
        kb = kbot(i,j)
         f_delt_n1(i,j) =     f(i,j,kb-1)*delt
        delta_tr_n(i,j) = dt_tr(i,j,kb  )*delt
    enddo
    enddo
 else
        nlev = size(mu,3)
        mu_delt_n(:,:) =       mu(:,:,nlev  )*delt
             nu_n(:,:) =       nu(:,:,nlev  )
             e_n1(:,:) =        e(:,:,nlev-1)
        f_delt_n1(:,:) =        f(:,:,nlev-1)*delt
       delta_tr_n(:,:) =    dt_tr(:,:,nlev  )*delt
 endif



!-----------------------------------------------------------------------

end subroutine vert_diff_down

!#######################################################################

subroutine vert_diff_down_2 &
      (delt, mu, nu, xi_1, xi_2, dt_xi_1, dt_xi_2, e, f_1, f_2, &
       mu_delt_n, nu_n, e_n1, f_1_delt_n1, f_2_delt_n1,         &
       delta_1_n, delta_2_n, kbot)

!-----------------------------------------------------------------------

real,    intent(in)                       :: delt
real,    intent(in)    , dimension(:,:,:) :: mu, nu, xi_1, xi_2
real,    intent(in)    , dimension(:,:,:) :: dt_xi_1, dt_xi_2
real,    intent(out)   , dimension(:,:,:) :: e, f_1, f_2
real,    intent(out)   , dimension(:,:)   :: mu_delt_n, nu_n, e_n1,    &
                                             f_1_delt_n1, f_2_delt_n1, &
                                             delta_1_n, delta_2_n

integer, intent(in),    dimension(:,:), optional :: kbot

real, dimension(size(xi_1,1),size(xi_1,2),size(xi_1,3)) :: a, b, c, g, &
                                                      dt_xi_11, dt_xi_22

integer :: i, j, kb, nlev

!-----------------------------------------------------------------------

! local copy of input 
  dt_xi_11 = dt_xi_1
  dt_xi_22 = dt_xi_2

 call explicit_tend (mu, nu, xi_1, dt_xi_11)
 call explicit_tend (mu, nu, xi_2, dt_xi_22)

 call compute_e (delt, mu, nu, e, a, b, c, g)

 call compute_f (dt_xi_11, b, c, g, f_1)
 call compute_f (dt_xi_22, b, c, g, f_2)

 if (present(kbot)) then
    do j=1,size(xi_1,2)
    do i=1,size(xi_1,1)
        kb = kbot(i,j)
        mu_delt_n(i,j)  =      mu(i,j,kb  )*delt
             nu_n(i,j)  =      nu(i,j,kb  )
            e_n1(i,j)  =       e(i,j,kb-1)
     f_1_delt_n1(i,j)  =     f_1(i,j,kb-1)*delt
     f_2_delt_n1(i,j)  =     f_2(i,j,kb-1)*delt
        delta_1_n(i,j)  = dt_xi_11(i,j,kb  )*delt
        delta_2_n(i,j)  = dt_xi_22(i,j,kb  )*delt
    enddo
    enddo
 else
        nlev = size(mu,3)
        mu_delt_n(:,:)  =      mu(:,:,nlev  )*delt
             nu_n(:,:)  =      nu(:,:,nlev  )
            e_n1(:,:)  =       e(:,:,nlev-1)
     f_1_delt_n1(:,:)  =     f_1(:,:,nlev-1)*delt
     f_2_delt_n1(:,:)  =     f_2(:,:,nlev-1)*delt
        delta_1_n(:,:)  = dt_xi_11(:,:,nlev  )*delt
        delta_2_n(:,:)  = dt_xi_22(:,:,nlev  )*delt
 endif



!-----------------------------------------------------------------------

end subroutine vert_diff_down_2

!#######################################################################

subroutine diff_surface (mu_delt, nu, e_n1, f_delt_n1,  &
                         dflux_datmos, flux, factor, delta_xi)

!-----------------------------------------------------------------------

real, intent(in)   , dimension(:,:) :: mu_delt, nu, e_n1, f_delt_n1,  &
                                       dflux_datmos
real, intent(inout), dimension(:,:) :: flux, delta_xi
real, intent(in) :: factor

!-----------------------------------------------------------------------

 real, dimension(size(flux,1),size(flux,2)) :: dflux
 real :: fff

 fff = 1.0/factor

 dflux    = - nu*(1.0 - e_n1)
 delta_xi = delta_xi + mu_delt*nu*f_delt_n1

 delta_xi = (delta_xi + mu_delt*flux*fff)/&
                      (1.0 - mu_delt*(dflux + dflux_datmos*fff))  

 flux     = flux + dflux_datmos*delta_xi


!-----------------------------------------------------------------------

end subroutine diff_surface

!#######################################################################

subroutine vert_diff_up (delt, e, f, delta_xi_n, dt_xi, kbot)

!-----------------------------------------------------------------------

real,    intent(in)                      :: delt
real,    intent(in),    dimension(:,:,:) :: e, f
real,    intent(in) ,   dimension(:,:)   :: delta_xi_n
real,    intent(out),   dimension(:,:,:) :: dt_xi
integer, intent(in),    dimension(:,:), optional :: kbot

integer :: i, j, k, kb, nlev
!-----------------------------------------------------------------------

 if (present(kbot)) then
     do j = 1, size(dt_xi,2)
     do i = 1, size(dt_xi,1)
         kb = kbot(i,j)
         dt_xi(i,j,kb) = delta_xi_n(i,j)/delt
         do k = kb -1, 1, -1
           dt_xi(i,j,k) = e(i,j,k)*dt_xi(i,j,k+1) + f(i,j,k)
         end do
     end do
     end do
 else
    nlev = size(dt_xi,3)
    dt_xi(:,:,nlev) = delta_xi_n/delt
    do k = size(dt_xi,3)-1, 1, -1
      dt_xi(:,:,k) = e(:,:,k)*dt_xi(:,:,k+1) + f(:,:,k)
    end do
 endif

!-----------------------------------------------------------------------

end subroutine vert_diff_up

!#######################################################################

subroutine compute_e (delt, mu, nu, e, a, b, c, g)

!-----------------------------------------------------------------------

real,    intent(in)                       :: delt
real,    intent(in)    , dimension(:,:,:) :: mu, nu
real,    intent(out)   , dimension(:,:,:) :: e, a, b, c, g

integer :: k, nlev

!-----------------------------------------------------------------------

 nlev = size(mu,3)

 a(:,:,1:nlev-1) = - mu(:,:,1:nlev-1)*nu(:,:,2:nlev)*delt
 a(:,:,nlev    ) =   0.0
 c(:,:,2:nlev  ) = - mu(:,:,2:nlev  )*nu(:,:,2:nlev)*delt
 c(:,:,1       ) =   0.0

 b = 1.0 - a - c

 e(:,:,1)   =   - a(:,:,1)/b(:,:,1)
 do  k= 2, nlev - 1
    g(:,:,k) = 1.0/(b(:,:,k) + c(:,:,k)*e(:,:,k-1))
    e(:,:,k) = - a(:,:,k)*g(:,:,k)
 enddo

!-----------------------------------------------------------------------

end subroutine compute_e

!#######################################################################

subroutine compute_f (dt_xi, b, c, g, f)

!-----------------------------------------------------------------------
real,    intent(in)    , dimension(:,:,:) :: dt_xi, b, c, g
real,    intent(out)   , dimension(:,:,:) :: f

integer :: k
!-----------------------------------------------------------------------

 f(:,:,1) =   dt_xi(:,:,1)/b(:,:,1)

 do  k = 2, size(b,3)-1
    f(:,:,k) = (dt_xi(:,:,k) - c(:,:,k)*f(:,:,k-1))*g(:,:,k)
 enddo

!-----------------------------------------------------------------------

end subroutine compute_f

!#######################################################################

subroutine explicit_tend (mu, nu, xi, dt_xi)

!-----------------------------------------------------------------------

real,    intent(in)    , dimension(:,:,:) :: mu, nu, xi
real,    intent(inout) , dimension(:,:,:) :: dt_xi

real, dimension(size(xi,1),size(xi,2),size(xi,3)) :: fluxx

integer :: nlev

!-----------------------------------------------------------------------

 nlev = size(mu,3)

 fluxx(:,:,1)      = 0.0
 fluxx(:,:,2:nlev) = nu(:,:,2:nlev)*(xi(:,:,2:nlev) - xi(:,:,1:nlev-1))

 dt_xi(:,:,1:nlev-1) = dt_xi(:,:,1:nlev-1) +  &
    mu(:,:,1:nlev-1)*(fluxx(:,:,2:nlev) - fluxx(:,:,1:nlev-1))
 dt_xi(:,:,nlev)     = dt_xi(:,:,nlev) - mu(:,:,nlev)*fluxx(:,:,nlev)

!-----------------------------------------------------------------------

end subroutine explicit_tend

!#######################################################################

subroutine compute_mu (p_half, mu)

!-----------------------------------------------------------------------
real,    intent(in)    , dimension(:,:,:) :: p_half
real,    intent(out)   , dimension(:,:,:) :: mu

integer :: nlev
!-----------------------------------------------------------------------

nlev = size(mu,3)

mu(:,:,1:nlev) = GRAV / (p_half(:,:,2:nlev+1) -p_half(:,:,1:nlev))

!-----------------------------------------------------------------------

end subroutine compute_mu


!#######################################################################

subroutine compute_nu (diff, p_half, p_full, z_full, t, q, nu)

!-----------------------------------------------------------------------
real,    intent(in)    , dimension(:,:,:) :: diff, p_half, p_full, &
                                             z_full, t, q
real,    intent(out)   , dimension(:,:,:) :: nu

real, dimension(size(t,1),size(t,2),size(t,3)) :: rho_half, tt
integer :: nlev
!-----------------------------------------------------------------------

nlev = size(nu,3)

if ( use_virtual_temp_vert_diff ) then
  tt = t * (1.0 + d608*q)           ! virtual temperature
else
  tt = t ! Take out virtual temperature effect here to mimic supersource
endif

rho_half(:,:,2:nlev) =  &         ! density at half levels
      2.0*p_half(:,:,2:nlev)/(RDGAS*(tt(:,:,2:nlev)+tt(:,:,1:nlev-1)))

if(do_mcm_plev) then
  nu(:,:,2:nlev) = GRAV*rho_half(:,:,2:nlev)*rho_half(:,:,2:nlev)*diff(:,:,2:nlev)/ &
                    (p_full(:,:,2:nlev)-p_full(:,:,1:nlev-1))
else
  nu(:,:,2:nlev) = rho_half(:,:,2:nlev)*diff(:,:,2:nlev) /  &
                    (z_full(:,:,1:nlev-1)-z_full(:,:,2:nlev))
endif
!-----------------------------------------------------------------------

end subroutine compute_nu

!#######################################################################

end module vert_diff_mod

