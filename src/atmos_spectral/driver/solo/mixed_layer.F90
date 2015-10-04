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

module mixed_layer_mod

!
! Implementation of mixed layer boundary condition
!

use            fms_mod, only: set_domain, write_version_number, &
                              mpp_pe, mpp_root_pe, error_mesg, FATAL, WARNING

use            fms_mod, only: stdlog, check_nml_error, close_file,&
                              open_namelist_file, stdout, file_exist, &
                              read_data, write_data, open_file, &
                              nullify_domain, lowercase

use  field_manager_mod, only: MODEL_ATMOS

use tracer_manager_mod, only: get_tracer_names, get_number_tracers

use      constants_mod, only: HLV, PI, RHO_CP, CP_AIR 

use   diag_manager_mod, only: register_diag_field, send_data

use   time_manager_mod, only: time_type

use     transforms_mod, only: get_deg_lat, grid_domain

use      vert_diff_mod, only: surf_diff_type

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: mixed_layer.F90,v 1.1.2.1 2013/01/24 20:35:37 pjp Exp $'
      
character(len=128) :: tagname= &
'$Name:  $'
character(len=128), parameter :: mod_name='mixed_layer'

!=================================================================================================================================

public :: mixed_layer_init, mixed_layer, mixed_layer_end

!=================================================================================================================================

logical :: evaporation = .true.
real    :: qflux_amp = 0.0
real    :: qflux_width = 16.0  ! width of qflux region in degrees
real    :: depth = 40.0
logical :: load_qflux = .false.


namelist/mixed_layer_nml/ evaporation, qflux_amp, depth, qflux_width, load_qflux

!=================================================================================================================================


logical :: module_is_initialized =.false.
logical :: used

integer :: iter, nhum
integer, dimension(4) :: axes
integer ::                                                                    &
     id_t_surf,            &   ! surface temperature
     id_flux_lhe,          &   ! latent heat flux at surface
     id_flux_oceanq,       &   ! oceanic Q flux 
     id_flux_t                 ! sensible heat flux at surface

real, allocatable, dimension(:,:)   ::                                        &
     ocean_qflux,           &   ! Q-flux 
     rad_lat_2d                 ! latitude in radians 

real, allocatable, dimension(:)   :: deg_lat

real, allocatable, dimension(:,:)   ::                                        &
     gamma_t,               &   ! Used to calculate the implicit
     gamma_q,               &   ! correction to the diffusion in
     fn_t,                  &   ! the lowest layer
     fn_q,                  &   ! 
     en_t,                  &   !
     en_q,                  &   !
     alpha_t,               &   !
     alpha_q,               &   !
     alpha_lw,              &   !
     beta_t,                &   !
     beta_q,                &   !
     beta_lw,               &   !
     t_surf_dependence,     &   !
     corrected_flux,        &   !
     eff_heat_capacity,     &   ! Effective heat capacity
     delta_t_surf               ! Increment in surface temperature

real inv_cp_air

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine mixed_layer_init(is, ie, js, je, num_levels, t_surf, axes, Time)

type(time_type), intent(in)       :: Time
real, intent(out), dimension(:,:) :: t_surf
integer, intent(in), dimension(4) :: axes
integer, intent(in) :: is, ie, js, je, num_levels 

integer :: j
real    :: rad_qwidth
integer:: ierr, io, unit, num_tr, n
character(32) :: tr_name

if(module_is_initialized) return

call write_version_number(version, tagname)

unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=mixed_layer_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'mixed_layer_nml')
enddo
10 call close_file (unit)

if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=mixed_layer_nml)

call get_number_tracers (MODEL_ATMOS, num_prog=num_tr)
do n = 1,num_tr
   call get_tracer_names( MODEL_ATMOS, n, tr_name )
   if(lowercase(tr_name)=='sphum') then
      nhum = n
   endif
enddo

allocate(rad_lat_2d              (is:ie, js:je))
allocate(ocean_qflux             (is:ie, js:je))
allocate(deg_lat                 (js:je))
allocate(gamma_t                 (is:ie, js:je))
allocate(gamma_q                 (is:ie, js:je))
allocate(en_t                    (is:ie, js:je))
allocate(en_q                    (is:ie, js:je))
allocate(fn_t                    (is:ie, js:je))
allocate(fn_q                    (is:ie, js:je))
allocate(alpha_t                 (is:ie, js:je))
allocate(alpha_q                 (is:ie, js:je))
allocate(alpha_lw                (is:ie, js:je))
allocate(beta_t                  (is:ie, js:je))
allocate(beta_q                  (is:ie, js:je))
allocate(beta_lw                 (is:ie, js:je))
allocate(delta_t_surf            (is:ie, js:je))
allocate(eff_heat_capacity       (is:ie, js:je))
allocate(corrected_flux          (is:ie, js:je))
allocate(t_surf_dependence       (is:ie, js:je))
!
!see if restart file exists for the surface temperature
!
if (file_exist('INPUT/mixed_layer.res.nc')) then

   call nullify_domain()
   call read_data(trim('INPUT/mixed_layer.res'), 't_surf',   t_surf, grid_domain)

else if (file_exist('INPUT/swamp.res')) then
         unit = open_file (file='INPUT/swamp.res', &
                           form='native', action='read')
         call read_data (unit, t_surf)
         call close_file (unit)
  call error_mesg('mixed_layer','mixed_layer restart file not found, using swamp restart file', WARNING)
else
  call error_mesg('mixed_layer','mixed_layer restart file not found', WARNING)
endif

id_t_surf = register_diag_field(mod_name, 't_surf',        &
                                axes(1:2), Time, 'surface temperature','K')
id_flux_t = register_diag_field(mod_name, 'flux_t',        &
                                axes(1:2), Time, 'sensible heat flux up at surface','watts/m2')
id_flux_lhe = register_diag_field(mod_name, 'flux_lhe',        &
                                 axes(1:2), Time, 'latent heat flux up at surface','watts/m2')
id_flux_oceanq = register_diag_field(mod_name, 'flux_oceanq',        &
                                 axes(1:2), Time, 'oceanic Q-flux','watts/m2')

! latitude will be needed for oceanic q flux
call get_deg_lat(deg_lat)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*PI/180.
enddo

! calculate ocean Q flux
rad_qwidth = qflux_width*PI/180.
ocean_qflux = qflux_amp*(1-2.*rad_lat_2d**2/rad_qwidth**2) * &
        exp(- ((rad_lat_2d)**2/(rad_qwidth)**2))

! load Q flux 
if (load_qflux) then
  call read_data('INPUT/ocean_qflux.nc', 'ocean_qflux',  ocean_qflux)
endif

inv_cp_air = 1.0 / CP_AIR 

module_is_initialized = .true.

return
end subroutine mixed_layer_init

!=================================================================================================================================

subroutine mixed_layer (                                               &
     Time,                                                             &
     t_surf,                                                           &
     flux_t,                                                           &
     flux_q,                                                           &
     flux_r,                                                           &
     dt,                                                               &
     net_surf_sw_down,                                                 &
     surf_lw_down,                                                     &
     Tri_surf,                                                         &
     dhdt_surf,                                                        &
     dedt_surf,                                                        &
     dedq_surf,                                                        &
     drdt_surf,                                                        &
     dhdt_atm,                                                         &
     dedq_atm)         

! ---- arguments -----------------------------------------------------------
type(time_type), intent(in)       :: Time
real, intent(in),  dimension(:,:) :: &
     net_surf_sw_down, surf_lw_down
real, intent(in), dimension(:,:) :: &
     flux_t,    flux_q,     flux_r
real, intent(inout), dimension(:,:) :: t_surf
real, intent(in), dimension(:,:) :: &
   dhdt_surf, dedt_surf, dedq_surf, &
   drdt_surf, dhdt_atm, dedq_atm  
real, intent(in) :: dt
type(surf_diff_type), intent(inout) :: Tri_surf

if(.not.module_is_initialized) then
  call error_mesg('mixed_layer','mixed_layer module is not initialized',FATAL)
endif

! Need to calculate the implicit changes to the lowest level delta_q and delta_t
! - see the discussion in vert_diff.tech.ps
                                                                                                                                    
! Care is needed to differentiate between the sensible heat flux and the
! diffusive flux of temperature
                                                                                                                                    
gamma_t = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_t + dhdt_atm * inv_cp_air))
gamma_q = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_tr(:,:,nhum) + dedq_atm))
                                                                                                                                 
fn_t = gamma_t * (Tri_surf%delta_t + Tri_surf%dtmass * flux_t * inv_cp_air)
fn_q = gamma_q * (Tri_surf%delta_tr(:,:,nhum) + Tri_surf%dtmass * flux_q)
                                                                                                                                 
en_t = gamma_t * Tri_surf%dtmass * dhdt_surf * inv_cp_air
en_q = gamma_q * Tri_surf%dtmass * dedt_surf
                                                                                                                                    
!
! Note flux_sw doesn't depend on surface or lowest layer values
! Note drdt_atm is not used - should be fixed
!
alpha_t = flux_t * inv_cp_air + dhdt_atm * inv_cp_air * fn_t
alpha_q = flux_q + dedq_atm * fn_q
alpha_lw = flux_r
                                                                                                                                 
beta_t = dhdt_surf * inv_cp_air + dhdt_atm * inv_cp_air * en_t
beta_q = dedt_surf + dedq_atm * en_q
beta_lw = drdt_surf

!
! Implement mixed layer surface boundary condition
!
corrected_flux = - net_surf_sw_down - surf_lw_down + alpha_t * CP_AIR + alpha_lw + ocean_qflux
t_surf_dependence = beta_t * CP_AIR + beta_lw


if (evaporation) then
  corrected_flux = corrected_flux + alpha_q * HLV
  t_surf_dependence = t_surf_dependence + beta_q * HLV
endif

!
! Now update the mixed layer surface temperature using an implicit step
!
eff_heat_capacity = depth * RHO_CP + t_surf_dependence * dt

if (any(eff_heat_capacity .eq. 0.0))  then 
  write(*,*) 'mixed_layer: error', eff_heat_capacity
  call error_mesg('mixed_layer', 'Avoiding division by zero',fatal)
end if

delta_t_surf = - corrected_flux  * dt / eff_heat_capacity

t_surf = t_surf + delta_t_surf
                                                                                                                                    
!
! Finally calculate the increments for the lowest atmospheric layer
!
Tri_surf%delta_t = fn_t + en_t * delta_t_surf
Tri_surf%delta_tr(:,:,nhum) = fn_q + en_q * delta_t_surf


!
! Note:
! When using an implicit step there is not a clearly defined flux for a given timestep
!
if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time)
if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time)
if(id_flux_lhe > 0) used = send_data(id_flux_lhe, HLV * flux_q, Time)
if(id_flux_oceanq > 0)   used = send_data(id_flux_oceanq, ocean_qflux, Time)

end subroutine mixed_layer

!=================================================================================================================================

subroutine mixed_layer_end(t_surf)

real, intent(inout), dimension(:,:) :: t_surf
integer:: unit

if(.not.module_is_initialized) return

! write a restart file for the surface temperature
call nullify_domain()
call write_data(trim('RESTART/mixed_layer.res'), 't_surf',   t_surf, grid_domain)

module_is_initialized = .false.

end subroutine mixed_layer_end

!=================================================================================================================================

end module mixed_layer_mod
