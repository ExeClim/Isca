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
                              nullify_domain, lowercase, field_size

use  field_manager_mod, only: MODEL_ATMOS

use tracer_manager_mod, only: get_tracer_names, get_number_tracers

use      constants_mod, only: HLV, PI, RHO_CP, CP_AIR, KELVIN

use   diag_manager_mod, only: register_diag_field, register_static_field, send_data

use   time_manager_mod, only: time_type     , get_time ! Chung add get_time

use     transforms_mod, only: get_deg_lat, get_deg_lon, grid_domain

use      vert_diff_mod, only: surf_diff_type

use         mpp_domains_mod, only: mpp_get_global_domain !s added to enable qflux reading

! mj know about surface topography
use spectral_dynamics_mod,only: get_surf_geopotential
! mj read SSTs
use interpolator_mod, only: interpolate_type,interpolator_init&
     &,CONSTANT,interpolator
!mj q-flux
use qflux_mod, only: qflux_init,qflux,warmpool


implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: mixed_layer.F90,v 1.1.2.1 2013/01/24 20:35:37 pjp Exp $'

character(len=128) :: tagname= &
'$Name:  $'
character(len=128), parameter :: mod_name='mixed_layer'

!=================================================================================================================================

public :: mixed_layer_init, mixed_layer, mixed_layer_end, albedo_calc

!=================================================================================================================================

logical :: evaporation = .true.
real    :: qflux_amp = 0.0
real    :: qflux_width = 16.0  ! width of qflux region in degrees
logical :: load_qflux = .false.
logical :: time_varying_qflux = .false.
real    :: tconst = 305.0
real    :: delta_T = 40.0
logical :: prescribe_initial_dist = .false.
real    :: albedo_value    = 0.06

!s Surface heat capacity options
real    :: depth           = 40.0,         & !s 2013 implementation
           land_depth      = -1.,          & !mj 2013
           trop_depth      = -1.,          & !mj 2013
           trop_cap_limit  = 15.,          & !mj
           heat_cap_limit  = 60.,          & !mj
           np_cap_factor   =  1.,          & !mj
           albedo_exp      = 2.,           & !mj
           albedo_cntr     = 45.,          & !mj
           albedo_wdth     = 10.,          & !mj
           higher_albedo   = 0.10,         & !mj
           lat_glacier     = 60.,          & !mj
           land_h_capacity_prefactor = 1.0 !s where(land) heat_capcity = land_h_capacity_prefactor * depth * rho_cp

!s Surface albedo options
real    :: land_albedo_prefactor = 1.0 !s where(land) albedo = land_albedo_prefactor * albedo_value

!s Begin mj extra options
integer :: albedo_choice    = 1 ! 1->constant or following 'where(land)', 2->NH or SH step, 3->N-S symmetric step, 4->profile with albedo_exp, 5->tanh with albedo_cntr,albedo_wdth
logical :: do_qflux         = .false. !mj
logical :: do_warmpool      = .false. !mj
logical :: do_read_sst      = .false. !mj
logical :: do_sc_sst        = .false. !mj use specified SSTs
logical :: do_ape_sst       = .false. ! use the AquaPlanet Experiement (APE) sst profile.
logical :: specify_sst_over_ocean_only = .false.
logical :: do_calc_eff_heat_cap = .true. ! assumes specified SST are off the default.

character(len=256) :: sst_file 
character(len=256) :: sst_file_Chung,sst_field_Chung ! Chung 09/2021
character(len=256) :: land_option = 'none'
real,dimension(10) :: slandlon=0,slandlat=0,elandlon=-1,elandlat=-1
!s End mj extra options

character(len=256) :: qflux_file_name  = 'ocean_qflux'
character(len=256) :: qflux_field_name = 'ocean_qflux' !Only used when using a non-time-varying q-flux. Otherwise code assumes field_name = file_name. 

character(len=256) :: ice_file_name  = 'siconc_clim_amip'
real    :: ice_albedo_value = 0.7
real    :: ice_concentration_threshold = 0.5
logical :: update_albedo_from_ice = .false.

logical :: add_latent_heat_flux_anom = .false.
character(len=256) :: flux_lhe_anom_file_name  = 'INPUT/flux_lhe_anom.nc'
character(len=256) :: flux_lhe_anom_field_name = 'flux_lhe_anom'

!!!!! thermodynamic ice !!!!! Chung add 05/2021
logical :: do_prescribe_albedo     = .false. ! prescribe albedo (default ISCA model)
logical :: do_thermodynamic_albedo = .true.  ! use thermodynamic ice by Ian Eisenman, XZ 02/2018
 logical :: do_thickness_lock      = .false. ! fix ice thickness
 character(len=256) :: ice_thickness_file = 'INPUT/h_ice.nc' ! specified ice thickness
 character(len=256) :: ice_thickness_field= 'h_ice' ! (optional) variable name ! Chung add 09/2021
  logical :: do_fraction_lock      = .false. ! fix ice fraction
  character(len=256) :: ice_fraction_file = 'INPUT/a_ice.nc' ! specified ice fraction
  character(len=256) :: ice_fraction_field= 'a_ice' ! (optional) variable name !Chung add 11/2021
 integer  :: num_input_times_h = 72   ! number of times during year in input h_ice file ### ### Chung 09/2021
 ! Sea ice parameters by Ian Eisenman, XZ 02/2018
 logical :: sea_ice_in_mixed_layer = .true.  ! include sea ice model
 logical :: ice_as_albedo_only = .false.     ! no sea ice thermodynamics, but change albedo when ML temp = TFREEZE
 logical :: sfc_melt_from_file = .false.     ! whether to compute ice surface melt based on input surface temperature file
 character(len=256) :: sfc_melt_file = 'INPUT/t_sfc.nc' ! (optional) specified t_surf for ice surface melt
 character(len=256) :: sfc_melt_field= 't_surf' ! (optional) variable name ! Chung add 09/2021
 integer  :: num_input_times = 12 ! 72  ! number of times during year in input t_surf file
 ! parameters for sea ice model, XZ 02/2018
 real :: L_ice = 3e8 ! latent heat of fusion (J/m^3)
 real, parameter  :: k_ice = 2 ! conductivity of ice (W/m/K)
 real, parameter  :: ice_basal_flux_const = 120 ! linear coefficient for heat flux from ocean to ice base (W/m^2/K)
 real :: t_ice_base = 273.15 ! temperature at base of ice, taken to be freshwater freezing point
 real :: t_surf_freeze = 273.15
 real    :: thermodynamic_albedo_ocn = 0.1  ! surface albedo where there is not sea ice, XZ 02/2018
 real    :: thermodynamic_albedo_ice = 0.4  ! surface albedo where there is sea ice, XZ 02/2018
!!!!!

namelist/mixed_layer_nml/ evaporation, depth, qflux_amp, qflux_width, tconst,&
                              delta_T, prescribe_initial_dist,albedo_value,  &
                              land_depth,trop_depth,                         &  !mj
                              trop_cap_limit, heat_cap_limit, np_cap_factor, &  !mj
                              do_qflux,do_warmpool,                          &  !mj
                              albedo_choice,higher_albedo,albedo_exp,        &  !mj
                              albedo_cntr,albedo_wdth,lat_glacier,           &  !mj
                              do_read_sst,do_sc_sst,sst_file,                &  !mj
                              land_option,slandlon,slandlat,                 &  !mj
                              elandlon,elandlat,                             &  !mj
                              land_h_capacity_prefactor,                     &  !s
                              land_albedo_prefactor,                         &  !s
                              load_qflux,qflux_file_name,time_varying_qflux, &
                              update_albedo_from_ice, ice_file_name,         &
                              ice_albedo_value, specify_sst_over_ocean_only, &
                              ice_concentration_threshold,                   &
                              add_latent_heat_flux_anom,flux_lhe_anom_file_name,&
                              flux_lhe_anom_field_name, do_ape_sst, qflux_field_name  ,&
                              !!!!!
                              sst_file_Chung, sst_field_Chung , &
                              !!!!! thermodynamic ice !!!!! Chung
                              do_prescribe_albedo, do_thermodynamic_albedo, &
                              do_thickness_lock,ice_thickness_file,ice_thickness_field, &
                              do_fraction_lock,ice_fraction_file,ice_fraction_field, & !!! Chung 11/2021
                              sea_ice_in_mixed_layer,ice_as_albedo_only, sfc_melt_from_file,sfc_melt_file,sfc_melt_field,&
                              t_surf_freeze, thermodynamic_albedo_ocn, thermodynamic_albedo_ice

!=================================================================================================================================


logical :: module_is_initialized =.false.
logical :: used

integer :: iter, nhum
integer, dimension(4) :: axes
integer ::                                                                    &
     id_t_surf,            &   ! surface temperature
     id_flux_lhe,          &   ! latent heat flux at surface
     id_flux_oceanq,       &   ! oceanic Q flux
     id_flux_t,            &   ! sensible heat flux at surface
     id_heat_cap,          &   ! heat capacity
     id_albedo,            &   ! mj albedo
     id_ice_conc,          &   ! st ice concentration
     id_delta_t_surf,      &
     !!!!! thermodynamic ice !!!!! Chung
     ! Sea ice by Ian Eisenman, XZ 02/2018
     id_h_ice,             &   ! sea ice thickness
     id_a_ice,             &   ! sea ice fractional area
     id_t_ml,              &   ! mixed layer temperature
     id_flux_ice               ! conductive heat flux through ice

real, allocatable, dimension(:,:)   ::                                        &
     ocean_qflux,           &   ! Q-flux
     ice_concentration,     &   ! ice_concentration
     rad_lat_2d,            &   ! latitude in radians
     flux_lhe_anom, flux_q_total, &
     !!!!! thermodynamic ice !!!!! Chung
     ! Sea ice by Ian Eisenman, XZ 02/2018
     t_surf_for_melt            ! if sfc_melt_from_file, current time t_surf from input file (otherwise, =t_surf)

real, allocatable, dimension(:)   :: deg_lat, deg_lon

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
     delta_t_surf,          &   ! Increment in surface temperature
     zsurf,                 &   ! mj know about topography
     land_sea_heat_capacity,&
     sst_new,               &   ! mj input SST
     albedo_initial!,        &

real, allocatable, dimension(:,:)   ::                                        & ! Chung 10/19/2021
     !!!!! thermodynamic ice !!!!! Chung
     ! Sea ice by Ian Eisenman, XZ 02/2018
     dFdt_surf,             &   ! d(corrected_flux)/d(t_surf), for calculation of ice sfc temperature
     delta_t_ml,            &   ! Increment in mixed layer temperature
     delta_h_ice,           &   ! Increment in sea ice thickness
     delta_t_ice,           &   ! Increment in ice (steady-state) surface temperature
     flux_ice                   ! Conductive heat flux through ice

logical, allocatable, dimension(:,:) ::      land_mask

!mj read sst from input file
  type(interpolate_type),save :: sst_interp
  type(interpolate_type),save :: qflux_interp
  type(interpolate_type),save :: ice_interp
  type(interpolate_type),save :: flux_lhe_anom_interp  

!!!!! thermodynamic ice !!!!! Chung add 05/2021
real, allocatable, dimension(:,:,:) ::  input_t_sfc ! Sfc temp read from file for ice sfc melt
real, allocatable, dimension(:,:,:) ::  input_h_ice ! ice thickness read from file for thickness locking
real, allocatable, dimension(:,:,:) ::  input_a_ice ! ice fraction read from file for thickness locking 11/2021
!real, allocatable, dimension(:,:,:) ::  input_tt, input_hh ! Chung 09/2021


real inv_cp_air

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine mixed_layer_init(is, ie, js, je, num_levels, t_surf, bucket_depth, axes, Time, albedo, rad_lonb_2d,rad_latb_2d, land, restart_file_bucket_depth ,&
!!!!! thermodynamic ice !!!!! Chung
h_ice, a_ice,t_ml )

type(time_type), intent(in)       :: Time
real, intent(out), dimension(:,:) :: t_surf, albedo 
real, intent(out), dimension(:,:) :: h_ice, a_ice, t_ml !!!!! thermodynamic Chung
real, intent(out), dimension(:,:,:) :: bucket_depth
integer, intent(in), dimension(4) :: axes
real, intent(in), dimension(:,:) :: rad_lonb_2d, rad_latb_2d
integer, intent(in) :: is, ie, js, je, num_levels

logical, intent(in), dimension(:,:) :: land
logical, intent(in)                 :: restart_file_bucket_depth

integer :: j
real    :: rad_qwidth
integer:: ierr, io, unit, num_tr, n, global_num_lon, global_num_lat
character(32) :: tr_name

! mj shallower ocean in tropics, land-sea contrast
 real :: trop_capacity,land_capacity,lon,lat,loc_cap
 integer :: i,k
integer, dimension(4) :: siz
character(len=12) :: ctmp1='     by     ', ctmp2='     by     '

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
allocate(ice_concentration       (is:ie, js:je))
allocate(flux_lhe_anom           (is:ie, js:je))
allocate(flux_q_total            (is:ie, js:je))
allocate(deg_lat                 (js:je))
allocate(deg_lon                 (is:ie))
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
allocate (albedo_initial         (is:ie, js:je))
allocate(land_sea_heat_capacity  (is:ie, js:je))
allocate(zsurf                   (is:ie, js:je))
allocate(sst_new                 (is:ie, js:je))
allocate(land_mask                 (is:ie, js:je)); land_mask=land
!!!!! thermodynamic ice !!!!! Chung
! Sea ice by Ian Eisenman, XZ 02/2018
allocate(t_surf_for_melt         (is:ie, js:je))
allocate(dFdt_surf               (is:ie, js:je))
allocate(delta_t_ml              (is:ie, js:je))
allocate(delta_h_ice             (is:ie, js:je))
allocate(delta_t_ice             (is:ie, js:je))
allocate(flux_ice                (is:ie, js:je))
allocate(input_t_sfc             (ie-is+1, je-js+1, num_input_times+2))
allocate(input_h_ice             (ie-is+1, je-js+1, num_input_times_h+2))
!allocate(input_tt    (ie-is+1, je-js+1, num_input_times))   ! Chung 09/2021
!allocate(input_hh    (ie-is+1, je-js+1, num_input_times_h)) ! Chung 09/2021
allocate(input_a_ice             (ie-is+1, je-js+1, num_input_times_h+2))

! latitude will be needed for oceanic q flux
!s Moved up slightly so that rad_lat_2d can be used in initial temperature distribution if necessary.

call get_deg_lat(deg_lat)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*PI/180.
enddo
!mj get lon for land_sea_heat_capacity
call get_deg_lon(deg_lon)

!s Adding MiMA options
   if(do_sc_sst) do_read_sst = .true.
   trop_capacity   = trop_depth*RHO_CP
   land_capacity   = land_depth*RHO_CP
   if(trop_capacity .le. 0.) trop_capacity = depth*RHO_CP
   if(land_capacity .le. 0.) land_capacity = depth*RHO_CP
!s End MiMA options

    !mj read fixed SSTs
    if( do_read_sst ) then
!        call interpolator_init( sst_interp, trim(sst_file)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
      !! Chung 09/2021
      call interpolator_init( sst_interp, trim(sst_file_Chung), rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
    endif


!!!!! thermodynamic ice !!!!! Chung
  if ( do_prescribe_albedo .and. do_thermodynamic_albedo)  then
    write(*,*) 'mixed_layer: error, do_prescribe_albedo and do_thermodynamic_albedo cannot be true at the same time'
  elseif  (.not. do_prescribe_albedo .and. .not. do_thermodynamic_albedo) then
    write(*,*) 'mixed_layer: error, do_prescribe_albedo and do_thermodynamic_albedo cannot be false at the same time'
  endif

if (do_thermodynamic_albedo) then
  ! Sea ice by Ian Eisenman, XZ 02/2018
  ! read seasonally-varying input_sfc_melt(lat,lon,time)
  if ( sfc_melt_from_file .and. file_exist(trim(sfc_melt_file)) ) then
    do j = 1, num_input_times
      call read_data( trim(sfc_melt_file),trim(sfc_melt_field) ,input_t_sfc(:,:,j+1), grid_domain, timelevel=j ) !! Chung
!      call read_data( trim(sfc_melt_file),trim(sfc_melt_field),input_tt, grid_domain,timelevel=j ) !! 2021/09/12 !! Chung !!
    end do
!    input_t_sfc(:,:,2:num_input_times+1)=input_tt !! 2021/09/12 !! Chung
    ! wrap edges of seasonal cycle for interpolation
    input_t_sfc(:,:,1)=input_t_sfc(:,:,num_input_times+1)
    input_t_sfc(:,:,num_input_times+2)=input_t_sfc(:,:,2)
  endif

  ! locking ice thickness by Chung 05/2021  ! read seasonally-varying input_h_ice(lat,lon,time)
  if ( do_thickness_lock .and. file_exist(trim(ice_thickness_file)) ) then
    do j = 1, num_input_times_h
      call read_data( trim(ice_thickness_file),trim(ice_thickness_field), input_h_ice(:,:,j+1), grid_domain, timelevel=j )
!      call read_data( trim(ice_thickness_file),trim(ice_thickness_field), input_hh, grid_domain, timelevel=j ) !! 2021/09/25 !! Chung
!      call read_data( trim(ice_thickness_file), 'h_ice', input_hh, grid_domain,timelevel=j ) !! 2021/09/12 !! Chung
    end do
!    input_h_ice(:,:,2:num_input_times_h+1)=input_hh !! 2021/09/12 !! Chung
    ! wrap edges of seasonal cycle for interpolation
    input_h_ice(:,:,1)=input_h_ice(:,:,num_input_times_h+1)
    input_h_ice(:,:,num_input_times_h+2)=input_h_ice(:,:,2)
  endif

  ! locking ice fraction by Chung 11/2021  ! read seasonally-varying input_a_ice(lat,lon,time)
  ! In order to solve the h_ice locking problem (cannot fix 100%)
  if ( do_fraction_lock .and. file_exist(trim(ice_fraction_file)) ) then
    do j = 1, num_input_times_h
      call read_data( trim(ice_fraction_file),trim(ice_fraction_field),input_a_ice(:,:,j+1), grid_domain, timelevel=j )
    end do
    ! wrap edges of seasonal cycle for interpolation
    input_a_ice(:,:,1)=input_a_ice(:,:,num_input_times_h+1)
    input_a_ice(:,:,num_input_times_h+2)=input_a_ice(:,:,2)
  endif

endif !!!!! Chung

!
!see if restart file exists for the surface temperature
!
if (file_exist('INPUT/mixed_layer.res.nc')) then

   call nullify_domain()
   call read_data(trim('INPUT/mixed_layer.res'), 't_surf',   t_surf, grid_domain)

   if (do_thermodynamic_albedo) then
   !!!!! thermodynamic ice !!!!! Chung add 05/2021 ! Sea ice by Ian Eisenman, added by XZ 02/2018
     call read_data(trim('INPUT/mixed_layer.res'), 'h_ice',  h_ice,grid_domain)
     call read_data(trim('INPUT/mixed_layer.res'), 'a_ice',  a_ice,grid_domain)
     call read_data(trim('INPUT/mixed_layer.res'), 't_ml',   t_ml,grid_domain)
   endif
   
   if (restart_file_bucket_depth) then
       call read_data(trim('INPUT/mixed_layer.res'), 'bucket_depth', bucket_depth, grid_domain)
   endif

else if (file_exist('INPUT/swamp.res')) then
         unit = open_file (file='INPUT/swamp.res', &
                           form='native', action='read')
         call read_data (unit, t_surf)
         call close_file (unit)
  call error_mesg('mixed_layer','mixed_layer restart file not found, using swamp restart file', WARNING)

else if( do_read_sst ) then !s Added so that if we are reading sst values then we can restart using them.

!   call interpolator( sst_interp, Time, t_surf, trim(sst_file) )
   call interpolator( sst_interp, Time, t_surf, trim(sst_field_Chung) ) !! Chung 09/2021

elseif (prescribe_initial_dist) then
!  call error_mesg('mixed_layer','mixed_layer restart file not found - initializing from prescribed distribution', WARNING)

    t_surf(:,:) = tconst - delta_T*((3.*sin(rad_lat_2d)**2.)-1.)/3.

else

  call error_mesg('mixed_layer','mixed_layer restart file not found - initializing from lowest model level temp', WARNING)

endif

id_t_surf = register_diag_field(mod_name, 't_surf',        &
                                axes(1:2), Time, 'surface temperature','K')
id_flux_t = register_diag_field(mod_name, 'flux_t',        &
                                axes(1:2), Time, 'sensible heat flux up at surface','watts/m2')
id_flux_lhe = register_diag_field(mod_name, 'flux_lhe',        &
                                 axes(1:2), Time, 'latent heat flux up at surface','watts/m2')
id_flux_oceanq = register_diag_field(mod_name, 'flux_oceanq',        &
                                 axes(1:2), Time, 'oceanic Q-flux','watts/m2')
id_heat_cap = register_static_field(mod_name, 'ml_heat_cap',        &
                                 axes(1:2), 'mixed layer heat capacity','joules/m^2/deg C')
id_delta_t_surf = register_diag_field(mod_name, 'delta_t_surf',        &
                                 axes(1:2), Time, 'change in sst','K')

if (do_prescribe_albedo) then ! == Chung: default code here
  if (update_albedo_from_ice) then
    id_albedo = register_diag_field(mod_name, 'albedo',    &
                                 axes(1:2), Time, 'surface albedo', 'none')
    id_ice_conc = register_diag_field(mod_name, 'ice_conc',    &
                                 axes(1:2), Time, 'ice_concentration', 'none')
  else
    id_albedo = register_static_field(mod_name, 'albedo',    &
                                 axes(1:2), 'surface albedo', 'none')
  endif
endif

if (do_thermodynamic_albedo) then !!!!! thermodynamic ice !!!!! Chung ! sea ice by Ian Eisenman, XZ 02/2018
  id_h_ice    = register_diag_field(mod_name, 'h_ice'   , axes(1:2), Time, 'sea ice thickness','m')
  id_a_ice    = register_diag_field(mod_name, 'a_ice'   , axes(1:2), Time, 'sea ice area','fraction of grid box')
  id_t_ml     = register_diag_field(mod_name, 't_ml'    , axes(1:2), Time, 'mixed layer tempeature','K')
  id_flux_ice = register_diag_field(mod_name, 'flux_ice', axes(1:2), Time, 'conductive heat flux through sea ice','watts/m2')
  id_albedo   = register_diag_field(mod_name, 'albedo'  , axes(1:2), Time, 'surface albedo', 'none')
endif

ocean_qflux = 0.

! load Q flux
if (load_qflux) then

    if (time_varying_qflux) then
       call interpolator_init( qflux_interp, trim(qflux_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
    else

       if(file_exist(trim('INPUT/'//qflux_file_name//'.nc'))) then
         call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat)
         call field_size(trim('INPUT/'//qflux_file_name//'.nc'), trim(qflux_field_name), siz)
         if ( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) then
           call read_data(trim('INPUT/'//qflux_file_name//'.nc'), trim(qflux_field_name), ocean_qflux, grid_domain)
         else
           write(ctmp1(1: 4),'(i4)') siz(1)
           write(ctmp1(9:12),'(i4)') siz(2)
           write(ctmp2(1: 4),'(i4)') global_num_lon
           write(ctmp2(9:12),'(i4)') global_num_lat
           call error_mesg ('get_qflux','Qflux file contains data on a '// &
                  ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
         endif
       else
         call error_mesg('get_qflux','load_qflux="'//trim('True')//'"'// &
                         ' but '//trim(qflux_file_name)//' does not exist', FATAL)
       endif

    endif
endif

!s Adding MiMA options for qfluxes.

if ( do_qflux .or. do_warmpool) then
   call qflux_init
!mj q-flux as in Merlis et al (2013) [Part II]
   if ( do_qflux ) call qflux(rad_latb_2d(1,:),ocean_qflux)
!mj q-flux to create a tropical temperature perturbation
   if ( do_warmpool) call warmpool(rad_lonb_2d(:,1),rad_latb_2d(1,:),ocean_qflux)
endif

!s End MiMA options for qfluxes

if (add_latent_heat_flux_anom) then
       call interpolator_init( flux_lhe_anom_interp, trim(flux_lhe_anom_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
endif



inv_cp_air = 1.0 / CP_AIR

!!!!! Chung !!!!!
if (do_prescribe_albedo) then ! === default ISCA setting  ===
!s Prescribe albedo distribution here so that it will be the same in both two_stream_gray later and rrtmg radiation.

albedo(:,:) = albedo_value

if(trim(land_option) .eq. 'input') then

where(land) albedo = land_albedo_prefactor * albedo

endif

!mj MiMA albedo choices.
select case (albedo_choice)
  case (2) ! higher_albedo northward (lat_glacier>0) or southward (lat_glacier <0 ) of lat_glacier
    do j = 1, size(t_surf,2)
      lat = deg_lat(js+j-1)
      ! mj SH or NH only
      if ( lat_glacier .ge. 0. ) then
         if ( lat > lat_glacier ) then
            albedo(:,j) = higher_albedo
         else
            albedo(:,j) = albedo_value
         endif
      else
         if ( lat < lat_glacier ) then
            albedo(:,j) = higher_albedo
         else
            albedo(:,j) = albedo_value
         endif
      endif
    enddo
  case (3) ! higher_albedo poleward of lat_glacier
    do j = 1, size(t_surf,2)
      lat = deg_lat(js+j-1)
      if ( abs(lat) > lat_glacier ) then
        albedo(:,j) = higher_albedo
      else
        albedo(:,j) = albedo_value
      endif
    enddo
  case (4) ! exponential increase with albedo_exp
     do j = 1, size(t_surf,2)
       lat = abs(deg_lat(js+j-1))
       albedo(:,j) = albedo_value + (higher_albedo-albedo_value)*(lat/90.)**albedo_exp
     enddo
  case (5) ! tanh increase around albedo_cntr with albedo_wdth
     do j = 1, size(t_surf,2)
       lat = abs(deg_lat(js+j-1))
       albedo(:,j) = albedo_value + (higher_albedo-albedo_value)*&
             0.5*(1+tanh((lat-albedo_cntr)/albedo_wdth))
     enddo
!!!  case (101) ! Chung: tanh increase around albedo_cntr with albedo_wdth
!     do j = 1, size(t_surf,2)
!       lat = abs(deg_lat(js+j-1))
!       albedo(:,j) = albedo_value + (higher_albedo-albedo_value)*&
!             0.5*(1+tanh((abs(lat)-albedo_cntr)/albedo_wdth))
!     enddo
!!!!!!
end select

albedo_initial=albedo

if (update_albedo_from_ice) then
    call interpolator_init( ice_interp, trim(ice_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
    call read_ice_conc(Time)
    call albedo_calc(albedo,Time)
else
    if ( id_albedo > 0 ) used = send_data ( id_albedo, albedo )
endif

endif !(do_prescribe_albedo) ! === default ISCA setting  === !!!!! Chung !!!!!

! Note: do_calc_eff_heat_cap is true by default and control when the surface 
! heat capacity is calculated (land and ocean). 
if (do_sc_sst) then
    ! if using specified sst then do not calc the heat capacity
    if (specify_sst_over_ocean_only) then
        ! unless the heat capacity over land is needed
        do_calc_eff_heat_cap = .true.
    else
        do_calc_eff_heat_cap = .false.
    endif
else
    if (do_ape_sst) then
        ! if using specified sst without land, do not calc the heat capacity
        do_calc_eff_heat_cap = .false.
    end if
endif



!s begin surface heat capacity calculation
if (do_calc_eff_heat_cap) then
    land_sea_heat_capacity = depth*RHO_CP
    if(trim(land_option) .ne. 'input') then
         if ( trop_capacity .ne. depth*RHO_CP .or. np_cap_factor .ne. 1. ) then !s Lines above make trop_capacity=depth*RHO_CP if trop_capacity set to be < 0.
            do j=js,je
               lat = deg_lat(j)
               if ( lat .gt. 0. ) then
                  loc_cap = depth*RHO_CP*np_cap_factor
               else
                  loc_cap = depth*RHO_CP
               endif
               if ( abs(lat) .lt. trop_cap_limit ) then
                  land_sea_heat_capacity(:,j) = trop_capacity
               elseif ( abs(lat) .lt. heat_cap_limit ) then
                  land_sea_heat_capacity(:,j) = trop_capacity*(1.-(abs(lat)-trop_cap_limit)/(heat_cap_limit-trop_cap_limit)) + (abs(lat)-trop_cap_limit)/(heat_cap_limit-trop_cap_limit)*loc_cap
               elseif ( lat .gt. heat_cap_limit ) then
                  land_sea_heat_capacity(:,j) = loc_cap
               end if
            enddo
         endif
! mj land heat capacity function of surface topography
         if(trim(land_option) .eq. 'zsurf')then
            call get_surf_geopotential(zsurf)
            where ( zsurf .gt. 10. ) land_sea_heat_capacity = land_capacity
         endif
! mj land heat capacity given through ?landlon, ?landlat
         if(trim(land_option) .eq. 'lonlat')then
            do j=js,je
           lat = deg_lat(j)
               do i=is,ie
                  lon = deg_lon(i)
                  do k=1,size(slandlat)
                     if ( lon .ge. slandlon(k) .and. lon .le. elandlon(k) &
                          &.and. lat .ge. slandlat(k) .and. lat .le. elandlat(k) )then
                        land_sea_heat_capacity(i,j) = land_capacity
                     endif
                  enddo
               enddo
            enddo
         endif
    else  !trim(land_option) .eq. 'input'
        where(land) land_sea_heat_capacity = land_h_capacity_prefactor*land_sea_heat_capacity
    endif !end of if (trim(land_option) .ne. 'input')
endif !end of if(.not.do_sc_sst)

if ( id_heat_cap > 0 ) used = send_data ( id_heat_cap, land_sea_heat_capacity )
!s end surface heat capacity calculation

module_is_initialized = .true.

return
end subroutine mixed_layer_init

!=================================================================================================================================

subroutine mixed_layer (                                               &
     Time,                                                             &
     Time_next,                                                        &
     js, je,                                                           &
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
     dedq_atm,                                                         &
     albedo_out,                                                       &
     albedo, h_ice, a_ice, t_ml) !!!!! thermodynamic ice !!!!! Chung !! is, ie , flux_u

! ---- arguments -----------------------------------------------------------
type(time_type), intent(in)         :: Time, Time_next
integer, intent(in)                 :: js, je
real, intent(in), dimension(:,:)    :: net_surf_sw_down, surf_lw_down
real, intent(in), dimension(:,:)    :: flux_t, flux_q, flux_r
real, intent(inout), dimension(:,:) :: t_surf     
real, intent(inout), dimension(:,:) :: h_ice, a_ice, t_ml !!!!! thermodynamic ice !!!!! Chung
real, intent(in), dimension(:,:)    :: dhdt_surf, dedt_surf, dedq_surf, &
                                       drdt_surf, dhdt_atm, dedq_atm
real, intent(in)                    :: dt
real, intent(out), dimension(:,:)   :: albedo_out   ,albedo !!!!! Chung add
type(surf_diff_type), intent(inout) :: Tri_surf
logical, dimension(size(land_mask,1),size(land_mask,2)) :: land_ice_mask
!!!!! thermodynamic ice !!!!! Chung ! for sfc_melt_from_file, sea ice by Ian Eisenman, XZ 2/10/2018
integer :: i, n, seconds, days
real    :: days_in_year    = 360 ! how many model days in solar year
real    :: day = 0.0

!local variables

integer :: j

if(.not.module_is_initialized) then
  call error_mesg('mixed_layer','mixed_layer module is not initialized',FATAL)
endif

!!!!! Chung !!!!!
if (do_prescribe_albedo) then ! === default ISCA setting  ===

if(update_albedo_from_ice) then
    call read_ice_conc(Time_next)
    land_ice_mask=.false.
    where(land_mask.or.(ice_concentration.gt.ice_concentration_threshold))
        land_ice_mask=.true.
    end where
else
    land_ice_mask=land_mask
endif

call albedo_calc(albedo_out,Time_next)

endif ! === default ISCA setting  ===  !!!!! Chung !!!!!

!s Add latent heat flux anomalies before any of the calculations take place

if (add_latent_heat_flux_anom) then
    call interpolator( flux_lhe_anom_interp, Time, flux_lhe_anom, trim(flux_lhe_anom_file_name) )
    flux_q_total = flux_q + flux_lhe_anom
else
    flux_q_total = flux_q
endif


! Need to calculate the implicit changes to the lowest level delta_q and delta_t
! - see the discussion in vert_diff.tech.ps

! Care is needed to differentiate between the sensible heat flux and the
! diffusive flux of temperature

gamma_t = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_t + dhdt_atm * inv_cp_air))
gamma_q = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_tr(:,:,nhum) + dedq_atm))

fn_t = gamma_t * (Tri_surf%delta_t + Tri_surf%dtmass * flux_t * inv_cp_air)
fn_q = gamma_q * (Tri_surf%delta_tr(:,:,nhum) + Tri_surf%dtmass * flux_q_total)

en_t = gamma_t * Tri_surf%dtmass * dhdt_surf * inv_cp_air
en_q = gamma_q * Tri_surf%dtmass * dedt_surf

!
! Note flux_sw doesn't depend on surface or lowest layer values
! Note drdt_atm is not used - should be fixed
!
alpha_t = flux_t * inv_cp_air + dhdt_atm * inv_cp_air * fn_t
alpha_q = flux_q_total + dedq_atm * fn_q
alpha_lw = flux_r

beta_t = dhdt_surf * inv_cp_air + dhdt_atm * inv_cp_air * en_t
beta_q = dedt_surf + dedq_atm * en_q
beta_lw = drdt_surf

! If time-varying qflux then update value
if(load_qflux.and.time_varying_qflux) then
         call interpolator( qflux_interp, Time, ocean_qflux, trim(qflux_file_name) )

     if(update_albedo_from_ice) then
          where (land_ice_mask) ocean_qflux=0.
     endif

endif

!
! Implement mixed layer surface boundary condition
!
corrected_flux = - net_surf_sw_down - surf_lw_down + alpha_t * CP_AIR + alpha_lw - ocean_qflux
t_surf_dependence = beta_t * CP_AIR + beta_lw
flux_ice = 0 ! Initialize conductive heat flux through sea ice !!!!! Chung

if (evaporation) then
  corrected_flux = corrected_flux + alpha_q * HLV
  t_surf_dependence = t_surf_dependence + beta_q * HLV
endif

!!!!! thermodynamic ice !!!!! Chung

if (do_thermodynamic_albedo .and. .not. do_thickness_lock ) then ! == Chung add

  ! for calculation of ice sfc temperature, sea ice by Ian Eisenman, XZ 02/2018
  dFdt_surf = dhdt_surf + drdt_surf + HLV * (dedt_surf) ! d(corrected_flux)/d(T_surf)

  ! Sea ice by Ian Eisenman, XZ 02/2018
  ! === (3) surface temperature increment is calculated with
  ! explicit forward time step using flux corrected for implicit atm;
  ! next, (4) surface state is updated  ===

  ! ======================================================================
  !
  !          Ocean mixed layer and sea ice model equations
  !
  ! ======================================================================
  !
  ! Mixed layer temperature is evolved is t_ml, and sea ice thickness is
  ! h_ice. Atmosphere cares about surface temparature (t_surf). Where
  ! h_ice=0, t_surf=t_ml, but where h_ice>0, t_surf=t_ice which is the
  ! ice surface temperature (assumed to be in steady-state: "zero layer
  ! model"). So h_ice, t_ml, and t_surf are all saved at each time step
  ! (t_ice does not need to be saved).
  !
  ! This model is equivalent to the Semtner (1976) "zero-layer model",
  ! with several differences (no snow, ice latent heat is same at base
  ! as surface, ice surface temperature calculated by including sensible
  ! and latent heat derivatives in dF/dT_surf, inclusion of frazil
  ! growth).

  if ( sea_ice_in_mixed_layer .and. .not. ice_as_albedo_only ) then ! === use sea ice model ===

    ! calculate values of delta_h_ice, delta_t_ml, and delta_t_ice

    where ( h_ice .le. 0 ) ! = ice-free = [ should be equivalent to use .eq. instead of .le. ]
      delta_t_ml = - ( corrected_flux + ocean_qflux) * dt/(depth*RHO_CP)
      delta_h_ice = 0
    elsewhere ! = ice-covered = ( h>0 )
      delta_t_ml = - ( ice_basal_flux_const * ( t_ml - t_ice_base ) + ocean_qflux) * dt/(depth*RHO_CP)
      delta_h_ice = ( corrected_flux - ice_basal_flux_const * ( t_ml - t_ice_base) ) * dt/L_ice
    endwhere
    where ( t_ml + delta_t_ml .lt. t_surf_freeze ) ! = frazil growth =
      delta_h_ice = delta_h_ice - ( t_ml + delta_t_ml - t_surf_freeze ) *(depth*RHO_CP)/L_ice
      delta_t_ml = t_surf_freeze - t_ml
    endwhere
    where ( ( h_ice .gt. 0 ) .and. ( h_ice + delta_h_ice .le. 0 ) ) ! = complete ablation =
      delta_t_ml = delta_t_ml - ( h_ice + delta_h_ice ) * L_ice/(depth*RHO_CP)
      delta_h_ice = - h_ice
    endwhere
    ! = update surface temperature =
    if ( sfc_melt_from_file ) then ! sfc melt from seasonally varying sfc temp specified from file
      !
      ! linearly interpolate from input file time to model time at each location
      ! interp1( ( (1:num_input_times)-0.5 )*days_in_year/num_input_times,
      ! albedo(x,y,:), day )
      ! = interp1( 1:num_input_times, albedo(x,y,:),
      ! day*num_input_times/days_in_year+0.5 )
      ! find model time
      ! call get_time(Time,seconds,days)
      ! day = days + seconds/86400
      ! ! make sure day is between 0 and days_in_year=360
      ! do while (day .lt. 0)
      !   day = day + days_in_year
      ! end do
      ! do while (day .ge. days_in_year)
      !   day = day - days_in_year
      ! end do
      ! ! find index of nearest input time below
      ! n=floor(day*num_input_times/days_in_year+1.5)
      ! do i = 1, size(input_t_sfc,1)
      !   do j = 1, size(input_t_sfc,2)
      !     t_surf_for_melt(i,j)=input_t_sfc(i,j,n) +
      !     (input_t_sfc(i,j,n+1)-input_t_sfc(i,j,n)) * &
      !        ( (day*num_input_times/days_in_year+1.5)-n )
      !   enddo
      ! enddo
      ! compute surface melt
      where ( h_ice + delta_h_ice .gt. 0 ) ! surface is ice-covered
        ! calculate increment in steady-state ice surface temperature
        flux_ice = k_ice / (h_ice + delta_h_ice) * ( t_ice_base - t_surf )
        delta_t_ice = ( - corrected_flux + flux_ice ) &
                      / ( k_ice / (h_ice + delta_h_ice) + dFdt_surf )
        ! in grid boxes with ice, wherever input t_surf=t_fr, make t_surf=t_fr;
        ! otherwise, let t_surf be whatever it wants (even t_surf>t_fr)
        ! where ( t_surf_for_melt .ge. TFREEZE ) ! surface ablation
        !   delta_t_ice = TFREEZE - t_surf
        ! endwhere
        ! surface is ice-covered, so update t_surf as ice surface temperature
        t_surf = t_surf + delta_t_ice
      elsewhere ! ice-free, so update t_surf as mixed layer temperature
        t_surf = t_ml + delta_t_ml
      endwhere
      !
    else ! no file specifying sfc melt (default)
      ! compute surface melt
      where ( h_ice + delta_h_ice .gt. 0 ) ! surface is ice-covered
        ! calculate increment in steady-state ice surface temperature
        flux_ice = k_ice / (h_ice + delta_h_ice) * ( t_ice_base - t_surf )
        delta_t_ice = ( - corrected_flux + flux_ice ) &
                      / ( k_ice / (h_ice + delta_h_ice) + dFdt_surf )
        where ( t_surf + delta_t_ice .gt. t_surf_freeze ) ! surface ablation
          delta_t_ice = t_surf_freeze - t_surf
        endwhere
        ! surface is ice-covered, so update t_surf as ice surface temperature
        t_surf = t_surf + delta_t_ice
      elsewhere ! ice-free, so update t_surf as mixed layer temperature
        t_surf = t_ml + delta_t_ml
      endwhere
    endif

  else ! === do not use sea ice model: just evolve mixed layer with explicit step===
    ! where ( land_mask .eq. 0 ) ! ocean
    !   delta_t_ml = - ( corrected_flux + ocean_qflux) * dt/(depth*RHO_CP)
    ! elsewhere ! land
    !   delta_t_ml = - ( corrected_flux + ocean_qflux) * dt/(depth_land*RHO_CP)
    ! endwhere
    delta_t_ml = - ( corrected_flux + ocean_qflux) * dt/(depth*RHO_CP)
    delta_h_ice = 0 ! do not evolve sea ice (no ice model)
    ! t_surf and t_ml are equal (they differ only when mixed layer is ice-covered)
    t_surf = t_ml + delta_t_ml
  endif

  ! = update state =
  t_ml = t_ml + delta_t_ml
  h_ice = h_ice + delta_h_ice




!!!!!!!!!!!!! Chung add 2022/10/15
  if (do_fraction_lock) then

    ! linearly interpolate from input file time to model time at each location
    ! find model time
    call get_time(Time,seconds,days)
!    day = days + seconds/86400
    day = real(days) + real(seconds)/86400   !!!!! It is an integer (by Chung 11/02/2021)
    ! make sure day is between 0 and days_in_year=360
    do while (day .lt. 0)
      day = day + days_in_year
    end do
    do while (day .ge. days_in_year)
      day = day - days_in_year
    end do
    ! find index of nearest input time below
    n=floor(day*num_input_times_h/days_in_year+1.5) ! integer:num_input_times_h=72
    do i = 1, size(input_a_ice,1)
      do j = 1, size(input_a_ice,2)
        a_ice(i,j)=input_a_ice(i,j,n)+(input_a_ice(i,j,n+1)-input_a_ice(i,j,n))*((day*num_input_times_h/days_in_year+1.5)-n)
      enddo
    enddo
  else !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  where ( h_ice .gt. 0 ) !!!
    a_ice = 1            !!!
  elsewhere              !!! original code from XZ
    a_ice = 0            !!!
  endwhere               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  endif !!!
!!!!!!!!!!!!! Chung add 2022/10/15




  !if ( sea_ice_in_mixed_layer .and. ice_as_albedo_only ) then ! === model with !only albedo changing when T < Tfr ===
  if ( ice_as_albedo_only ) then ! === model with only albedo changing when T < Tfr ===
    ! a_ice is passed to radiation module for albedo. just set a_ice=1 where   ! T<Tfr.
    where ( t_ml .le. t_surf_freeze)
      a_ice = 1
    elsewhere
      a_ice = 0
    endwhere
  endif

  ! ======================================================================

  ! === (5) save increments in T and Q for lowest atm layer ===

  !Tri_surf%delta_t = fn_t
  !Tri_surf%delta_q = fn_q

  !!!!! Chung add 05/2021 !!!!!
  albedo(:,:) = (1-a_ice(:,:))*thermodynamic_albedo_ocn + a_ice(:,:)*thermodynamic_albedo_ice

endif ! == thermodynamic ice, Chung add 05/2021

!!!!! Chung add 05/2021
if ( do_thickness_lock .and. file_exist(trim(ice_thickness_file)) ) then ! fix ice thickness

  ! linearly interpolate from input file time to model time at each location
  ! find model time
  call get_time(Time,seconds,days)
!  day = days + seconds/86400
  day = real(days) + real(seconds)/86400   !!!!! It is an integer (by Chung 11/02/2021)

  ! make sure day is between 0 and days_in_year=360
  do while (day .lt. 0)
    day = day + days_in_year
  end do
  do while (day .ge. days_in_year)
    day = day - days_in_year
  end do
  ! find index of nearest input time below
  n=floor(day*num_input_times_h/days_in_year+1.5) ! integer: num_input_times_h=72
!  n=floor(day*num_input_times_h/days_in_year+1.4) !!!!! TEST 10/28/2021 Chung
  do i = 1, size(input_h_ice,1)
    do j = 1, size(input_h_ice,2)
      h_ice(i,j)=input_h_ice(i,j,n)+(input_h_ice(i,j,n+1)-input_h_ice(i,j,n))*((day*num_input_times_h/days_in_year+1.5)-n)
!!!!! TEST 10/28/2021 Chung
!      h_ice(i,j)=input_h_ice(i,j,n)+(input_h_ice(i,j,n+1)-input_h_ice(i,j,n))*((day*num_input_times_h/days_in_year+1.4)-n)
    enddo
  enddo

  if ( do_fraction_lock .and. file_exist(trim(ice_fraction_file)) ) then ! fix ice fraction
    n=floor(day*num_input_times_h/days_in_year+1.5)
    do i = 1, size(input_a_ice,1)
      do j = 1, size(input_a_ice,2)
        a_ice(i,j)=input_a_ice(i,j,n)+(input_a_ice(i,j,n+1)-input_a_ice(i,j,n))*((day*num_input_times_h/days_in_year+1.5)-n)
      enddo
    enddo
  else
    where ( h_ice .gt. 0 )
!    where ( h_ice .gt. 0.02 ) !!! Chung 10/2021 changing ice thickness threshold to avoid precision problem
      a_ice = 1
    elsewhere
      a_ice = 0
    endwhere
  endif

  albedo(:,:) = (1-a_ice(:,:))*thermodynamic_albedo_ocn +a_ice(:,:)*thermodynamic_albedo_ice

endif



!if ( .not.(do_thickness_lock) .and. do_fraction_lock .and. file_exist(trim(ice_fraction_file)) ) then ! fix ice fraction !!!!! Chung 2022/09/26
!  n=floor(day*num_input_times_h/days_in_year+1.5)
!  do i = 1, size(input_a_ice,1)
!    do j = 1, size(input_a_ice,2)
!      a_ice(i,j)=input_a_ice(i,j,n)+(input_a_ice(i,j,n+1)-input_a_ice(i,j,n))*((day*num_input_times_h/days_in_year+1.5)-n)
!    enddo
!  enddo
!  albedo(:,:) = (1-a_ice(:,:))*thermodynamic_albedo_ocn+a_ice(:,:)*thermodynamic_albedo_ice
!endif
!!!!! end of Chung's addition




!s Surface heat_capacity calculation based on that in MiMA by mj

if(do_sc_sst) then !mj sst read from input file
     ! read at the new time, as that is what we are stepping to
!     call interpolator( sst_interp, Time_next, sst_new, trim(sst_file) )
     call interpolator( sst_interp, Time_next, sst_new, trim(sst_field_Chung) ) !! Chung 09/2021

     if(specify_sst_over_ocean_only) then
         where (.not.land_ice_mask) delta_t_surf = sst_new - t_surf
         where (.not.land_ice_mask) t_surf = t_surf + delta_t_surf             
     else
         delta_t_surf = sst_new - t_surf
         t_surf = t_surf + delta_t_surf
     endif
end if

if (do_ape_sst) then 
    !
    ! AquaPlanet Experiment protocol (APE) from 
    ! Williams et al 2012 tech report: "The APE Atlas"
    !
    ! use analytic form for setting SST at each timestep.
    ! see appendix equation 1 of Neale and Hoskins 2000
    !     "A standard test for AGCMs including their 
    !      physical parametrizations: I: the proposal"
    ! using a constant longitude.
    do j=js,je   
        if ( (rad_lat_2d(1,j) .gt. -PI/3.) .and. (rad_lat_2d(1,j) .lt. PI/3.) ) then 
            ! between 60N-60S
            sst_new(:,j) = KELVIN+( 27.0*( 1. - (sin( 3./2. * rad_lat_2d(:,j) )**2 ) ))  
            !write(6,*) 'SST profile', rad_lat_2d(1,j)*180/PI, ' j:', j, ' sst:', sst_new(1,j)
        else
            ! from 60N/S to pole
            sst_new(:,j) = KELVIN
            !write(6,*) 'SST is zero', rad_lat_2d(1,j)*180/PI, ' j:', j, ' sst:', sst_new(1,j)
        endif
    enddo
    delta_t_surf = sst_new - t_surf
    t_surf = sst_new
endif


if (do_calc_eff_heat_cap) then
  !s use the land_sea_heat_capacity calculated in mixed_layer_init

    ! Now update the mixed layer surface temperature using an implicit step
    !
    eff_heat_capacity = land_sea_heat_capacity + t_surf_dependence * dt !s need to investigate how this works

    if (any(eff_heat_capacity .eq. 0.0))  then
      write(*,*) 'mixed_layer: error', eff_heat_capacity
      call error_mesg('mixed_layer', 'Avoiding division by zero',fatal)
    end if

    if(do_sc_sst.and.specify_sst_over_ocean_only) then
        where (land_ice_mask) delta_t_surf = - corrected_flux  * dt / eff_heat_capacity
        where (land_ice_mask) t_surf = t_surf + delta_t_surf             
    else
        delta_t_surf = - corrected_flux  * dt / eff_heat_capacity
        t_surf = t_surf + delta_t_surf
    endif

endif !s end of if(do_sc_sst).

!
! Finally calculate the increments for the lowest atmospheric layer
!
Tri_surf%delta_t = fn_t + en_t * delta_t_surf
if (evaporation) Tri_surf%delta_tr(:,:,nhum) = fn_q + en_q * delta_t_surf

!
! Note:
! When using an implicit step there is not a clearly defined flux for a given timestep
! We have taken a time-step, send the values at the next time level.
if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time_next)
if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time_next)
if(id_flux_lhe > 0) used = send_data(id_flux_lhe, HLV * flux_q_total, Time_next)
if(id_flux_oceanq > 0)   used = send_data(id_flux_oceanq, ocean_qflux, Time_next)

if(id_delta_t_surf > 0)   used = send_data(id_delta_t_surf, delta_t_surf, Time_next)

!!!!! thermodynamic ice !!!!! Chung
! Sea ice by Ian Eisenman, XZ 2/10/2018
if(id_h_ice > 0) used = send_data(id_h_ice, h_ice, Time)
if(id_a_ice > 0) used = send_data(id_a_ice, a_ice, Time)
if(id_t_ml > 0) used = send_data(id_t_ml, t_ml, Time)
if(id_flux_ice > 0) used = send_data(id_flux_ice, flux_ice, Time)
if(id_albedo > 0) used = send_data(id_albedo, albedo, Time)

end subroutine mixed_layer

!=================================================================================================================================

subroutine albedo_calc(albedo_inout,Time)

real, intent(out), dimension(:,:) :: albedo_inout
type(time_type), intent(in)       :: Time

albedo_inout=albedo_initial

if(update_albedo_from_ice) then

    where(ice_concentration.gt.ice_concentration_threshold) 
        albedo_inout=ice_albedo_value
    end where

    if ( id_ice_conc > 0 ) used = send_data ( id_ice_conc, ice_concentration, Time )
    if ( id_albedo > 0 ) used = send_data ( id_albedo, albedo_inout, Time )

endif

end subroutine albedo_calc
!=================================================================================================================================

subroutine read_ice_conc(Time)

type(time_type), intent(in)       :: Time


call interpolator( ice_interp, Time, ice_concentration, trim(ice_file_name) )
if ( id_ice_conc > 0 ) used = send_data ( id_ice_conc, ice_concentration, Time )

end subroutine read_ice_conc
!=================================================================================================================================

subroutine mixed_layer_end(t_surf, bucket_depth, restart_file_bucket_depth   , h_ice, a_ice, t_ml ) ! Chung

real, intent(inout), dimension(:,:) :: t_surf  
real, intent(inout), dimension(:,:) :: h_ice, a_ice, t_ml ! Chung ! Added by XZ 02/2018
real, intent(inout), dimension(:,:,:) :: bucket_depth
logical, intent(in)                 :: restart_file_bucket_depth
integer:: unit

if(.not.module_is_initialized) return

! write a restart file for the surface temperature
call nullify_domain()
call write_data(trim('RESTART/mixed_layer.res'), 't_surf',   t_surf, grid_domain)
if (restart_file_bucket_depth) then
    call write_data(trim('RESTART/mixed_layer.res'), 'bucket_depth',   bucket_depth, grid_domain)
endif
!!!!! thermodynamic ice !!!!! Chung
! Sea ice by Ian Eisenman, XZ 2/11/2018
call write_data(trim('RESTART/mixed_layer.res'), 'h_ice', h_ice, grid_domain)
call write_data(trim('RESTART/mixed_layer.res'), 'a_ice', a_ice, grid_domain)
call write_data(trim('RESTART/mixed_layer.res'), 't_ml', t_ml, grid_domain)

module_is_initialized = .false.

end subroutine mixed_layer_end

!=================================================================================================================================

end module mixed_layer_mod
