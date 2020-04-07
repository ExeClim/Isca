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

use   time_manager_mod, only: time_type

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
logical :: do_sc_sst        = .false. !mj
logical :: do_ape_sst  = .false. 
logical :: specify_sst_over_ocean_only = .false.
character(len=256) :: sst_file
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

namelist/mixed_layer_nml/ evaporation, depth, qflux_amp, qflux_width, tconst,&
                              delta_T, prescribe_initial_dist,albedo_value,  &
                              land_depth,trop_depth,                         &  !mj
                              trop_cap_limit, heat_cap_limit, np_cap_factor, &  !mj
			                  do_qflux,do_warmpool,              &  !mj
                              albedo_choice,higher_albedo,albedo_exp,        &  !mj
                              albedo_cntr,albedo_wdth,lat_glacier,           &  !mj
                              do_read_sst,do_sc_sst,do_ape_sst,sst_file,     &  !mj
                              land_option,slandlon,slandlat,                 &  !mj
                              elandlon,elandlat,                             &  !mj
                              land_h_capacity_prefactor,                     &  !s
                              land_albedo_prefactor,                         &  !s
                              load_qflux,qflux_file_name,time_varying_qflux, &
                              update_albedo_from_ice, ice_file_name,         &
                              ice_albedo_value, specify_sst_over_ocean_only, &
                              ice_concentration_threshold,                   &
                              add_latent_heat_flux_anom,flux_lhe_anom_file_name,&
                              flux_lhe_anom_field_name, qflux_field_name

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
     id_delta_t_surf

real, allocatable, dimension(:,:)   ::                                        &
     ocean_qflux,           &   ! Q-flux
     ice_concentration,     &   ! ice_concentration
     rad_lat_2d,            &   ! latitude in radians
     flux_lhe_anom, flux_q_total

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
     albedo_initial

logical, allocatable, dimension(:,:) ::      land_mask

!mj read sst from input file
  type(interpolate_type),save :: sst_interp
  type(interpolate_type),save :: qflux_interp
  type(interpolate_type),save :: ice_interp
  type(interpolate_type),save :: flux_lhe_anom_interp  

real inv_cp_air

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine mixed_layer_init(is, ie, js, je, num_levels, t_surf, bucket_depth, axes, Time, albedo, rad_lonb_2d,rad_latb_2d, land, restart_file_bucket_depth)

type(time_type), intent(in)       :: Time
real, intent(out), dimension(:,:) :: t_surf, albedo
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
!
!see if restart file exists for the surface temperature
!


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
	   call interpolator_init( sst_interp, trim(sst_file)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
	endif




if (file_exist('INPUT/mixed_layer.res.nc')) then

   call nullify_domain()
   call read_data(trim('INPUT/mixed_layer.res'), 't_surf',   t_surf, grid_domain)
   
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

   call interpolator( sst_interp, Time, t_surf, trim(sst_file) )

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
if (update_albedo_from_ice) then
	id_albedo = register_diag_field(mod_name, 'albedo',    &
                                 axes(1:2), Time, 'surface albedo', 'none')
	id_ice_conc = register_diag_field(mod_name, 'ice_conc',    &
                                 axes(1:2), Time, 'ice_concentration', 'none')
else
	id_albedo = register_static_field(mod_name, 'albedo',    &
                                 axes(1:2), 'surface albedo', 'none')
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
end select

albedo_initial=albedo

if (update_albedo_from_ice) then
	call interpolator_init( ice_interp, trim(ice_file_name)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
        call read_ice_conc(Time)
	call albedo_calc(albedo,Time)
else
	if ( id_albedo > 0 ) used = send_data ( id_albedo, albedo )
endif

!s begin surface heat capacity calculation
   if(.not.do_sc_sst.or.(do_sc_sst.and.specify_sst_over_ocean_only)) then
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
     albedo_out)

! ---- arguments -----------------------------------------------------------
type(time_type), intent(in)         :: Time, Time_next
integer, intent(in)                 :: js, je
real, intent(in), dimension(:,:)    :: net_surf_sw_down, surf_lw_down
real, intent(in), dimension(:,:)    :: flux_t, flux_q, flux_r
real, intent(inout), dimension(:,:) :: t_surf
real, intent(in), dimension(:,:)    :: dhdt_surf, dedt_surf, dedq_surf, &
                                       drdt_surf, dhdt_atm, dedq_atm
real, intent(in)                    :: dt
real, intent(out), dimension(:,:)   :: albedo_out
type(surf_diff_type), intent(inout) :: Tri_surf
logical, dimension(size(land_mask,1),size(land_mask,2)) :: land_ice_mask

!local variables

integer :: j

if(.not.module_is_initialized) then
  call error_mesg('mixed_layer','mixed_layer module is not initialized',FATAL)
endif

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

if (evaporation) then
  corrected_flux = corrected_flux + alpha_q * HLV
  t_surf_dependence = t_surf_dependence + beta_q * HLV
endif

!s Surface heat_capacity calculation based on that in MiMA by mj


if(do_sc_sst) then !mj sst read from input file
     ! read at the new time, as that is what we are stepping to
     call interpolator( sst_interp, Time_next, sst_new, trim(sst_file) )

     if(specify_sst_over_ocean_only) then
         where (.not.land_ice_mask) delta_t_surf = sst_new - t_surf
         where (.not.land_ice_mask) t_surf = t_surf + delta_t_surf			 
	 else
	     delta_t_surf = sst_new - t_surf
	     t_surf = t_surf + delta_t_surf
	 endif
end if

if(do_ape_sst) then 
    ! use analytic form for setting SST at each timestep.
    ! see appendix equation 1 of Neale and Hoskins 2000
    !     "A standard test for AGCMs including their 
    !      physical parametrizations: I: the proposal"
    ! using a constant longitude.
    do j=js,je   
        if ( (rad_lat_2d(1,j) .gt. -PI/3.) .and. (rad_lat_2d(1,j) .lt. PI/3.) ) then 
            ! between 60N-60S
            sst_new(:,j) = KELVIN+( 27.0*( 1 - (sin( 3./2. * rad_lat_2d(:,j) )**2 ) ))  
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

if ((.not.do_sc_sst).or.(do_sc_sst.and.specify_sst_over_ocean_only) .or. .not.(do_ape_sst)) then
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

subroutine mixed_layer_end(t_surf, bucket_depth, restart_file_bucket_depth)

real, intent(inout), dimension(:,:) :: t_surf
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

module_is_initialized = .false.

end subroutine mixed_layer_end

!=================================================================================================================================

end module mixed_layer_mod
