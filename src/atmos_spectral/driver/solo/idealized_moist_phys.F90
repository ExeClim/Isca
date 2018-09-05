module idealized_moist_phys_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

use fms_mod, only: write_version_number, file_exist, close_file, stdlog, error_mesg, NOTE, FATAL, read_data, field_size, uppercase, mpp_pe

use           constants_mod, only: grav, rdgas, rvgas, cp_air, PSTD_MKS, dens_h2o !mj cp_air needed for rrtmg !s pstd_mks needed for pref calculation

use        time_manager_mod, only: time_type, get_time, operator( + )

use    vert_turb_driver_mod, only: vert_turb_driver_init, vert_turb_driver, vert_turb_driver_end

use           vert_diff_mod, only: vert_diff_init, gcm_vert_diff_down, gcm_vert_diff_up, vert_diff_end, surf_diff_type

use two_stream_gray_rad_mod, only: two_stream_gray_rad_init, two_stream_gray_rad_down, two_stream_gray_rad_up, two_stream_gray_rad_end

use        cloud_simple_mod, only: cloud_simple_init, cloud_simple_end, cloud_simple

use         mixed_layer_mod, only: mixed_layer_init, mixed_layer, mixed_layer_end, albedo_calc

use         lscale_cond_mod, only: lscale_cond_init, lscale_cond, lscale_cond_end

use qe_moist_convection_mod, only: qe_moist_convection_init, qe_moist_convection, qe_moist_convection_end

use                 ras_mod, only: ras_init, ras_end, ras

use        betts_miller_mod, only: betts_miller, betts_miller_init

use      dry_convection_mod, only: dry_convection_init, dry_convection

use        diag_manager_mod, only: register_diag_field, send_data

use          transforms_mod, only: get_grid_domain

use   spectral_dynamics_mod, only: get_axis_id, get_num_levels, get_surf_geopotential

use        surface_flux_mod, only: surface_flux, gp_surface_flux

use      sat_vapor_pres_mod, only: lookup_es !s Have added this to allow relative humdity to be calculated in a consistent way.

use      damping_driver_mod, only: damping_driver, damping_driver_init, damping_driver_end !s MiMA uses damping

use    press_and_geopot_mod, only: pressure_variables

use         mpp_domains_mod, only: mpp_get_global_domain !s added to enable land reading

use          transforms_mod, only: grid_domain

use tracer_manager_mod, only: get_number_tracers, query_method

use  field_manager_mod, only: MODEL_ATMOS

use rayleigh_bottom_drag_mod, only: rayleigh_bottom_drag_init, compute_rayleigh_bottom_drag

#ifdef RRTM_NO_COMPILE
    ! RRTM_NO_COMPILE not included
#else
    !mj: RRTM radiative scheme
    use rrtmg_lw_init
    use rrtmg_lw_rad
    use rrtmg_sw_init
    use rrtmg_sw_rad
    use rrtm_radiation
    use rrtm_vars
#endif

#ifdef SOC_NO_COMPILE
    ! Socrates not included
#else
use socrates_interface_mod
use soc_constants_mod
#endif

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: idealized_moist_phys.F90,v 1.1.2.1 2013/01/24 18:45:34 pjp Exp $'

character(len=128) :: tagname= &
'$Name:  $'
character(len=10), parameter :: mod_name='atmosphere'

!=================================================================================================================================

public :: idealized_moist_phys_init , idealized_moist_phys , idealized_moist_phys_end

logical :: module_is_initialized =.false.
logical :: turb = .false.
logical :: do_virtual = .false. ! whether virtual temp used in gcm_vert_diff

!s Convection scheme options
character(len=256) :: convection_scheme = 'unset'  !< Use a specific convection scheme.  Valid options
integer, parameter :: UNSET = -1,                & !! are NONE, SIMPLE_BETTS_MILLER, FULL_BETTS_MILLER, DRY
                      NO_CONV = 0,               &
                      SIMPLE_BETTS_CONV = 1,     &
                      FULL_BETTS_MILLER_CONV = 2,&
                      DRY_CONV = 3,              &
                      RAS_CONV = 4
                      
integer :: r_conv_scheme = UNSET  ! the selected convection scheme

logical :: lwet_convection = .false.
logical :: do_bm = .false.
logical :: do_ras = .false.

! Cloud options
logical :: do_cloud_simple = .false.

!s Radiation options
logical :: two_stream_gray = .true.
logical :: do_rrtm_radiation = .false.
logical :: do_socrates_radiation = .false.

!s MiMA uses damping
logical :: do_damping = .false.


logical :: mixed_layer_bc = .false.
logical :: gp_surface = .false. !s Use Schneider & Liu 2009's prescription of lower-boundary heat flux

logical :: do_simple = .false. !s Have added this to enable relative humidity to be calculated correctly below.
real :: roughness_heat = 0.05
real :: roughness_moist = 0.05
real :: roughness_mom = 0.05
real :: land_roughness_prefactor = 1.0

!s options for adding idealised land

character(len=256) :: land_option = 'none'
character(len=256) :: land_file_name  = 'INPUT/land.nc'
character(len=256) :: land_field_name = 'land_mask'

! RG Add bucket
logical :: bucket = .false. 
integer :: future
real :: init_bucket_depth = 1000. ! default large value
real :: init_bucket_depth_land = 20. 
real :: max_bucket_depth_land = 0.15 ! default from Manabe 1969
real :: robert_bucket = 0.04   ! default robert coefficient for bucket depth LJJ
real :: raw_bucket = 0.53       ! default raw coefficient for bucket depth LJJ
! end RG Add bucket

namelist / idealized_moist_phys_nml / turb, lwet_convection, do_bm, do_ras, roughness_heat,  &
                                      do_cloud_simple,                                       &
                                      two_stream_gray, do_rrtm_radiation, do_damping,&
                                      mixed_layer_bc, do_simple,                     &
                                      roughness_moist, roughness_mom, do_virtual,    &
                                      land_option, land_file_name, land_field_name,   & !s options for idealised land
                                      land_roughness_prefactor,               &
                                      gp_surface, convection_scheme,          &
                                      bucket, init_bucket_depth, init_bucket_depth_land, & !RG Add bucket 
                                      max_bucket_depth_land, robert_bucket, raw_bucket, &
                                      do_socrates_radiation


integer, parameter :: num_time_levels = 2 !RG Add bucket - number of time levels added to allow timestepping in this module
real, allocatable, dimension(:,:,:)   :: bucket_depth      ! RG Add bucket
real, allocatable, dimension(:,:    ) :: dt_bucket, filt   ! RG Add bucket

real, allocatable, dimension(:,:)   ::                                        &
     z_surf,               &   ! surface height
     t_surf,               &   ! surface temperature
     q_surf,               &   ! surface moisture
     u_surf,               &   ! surface U wind
     v_surf,               &   ! surface V wind
     rough_mom,            &   ! momentum roughness length for surface_flux
     rough_heat,           &   ! heat roughness length for surface_flux
     rough_moist,          &   ! moisture roughness length for surface_flux
     depth_change_lh,      &   ! tendency in bucket depth due to latent heat transfer     ! RG Add bucket
     depth_change_cond,    &   ! tendency in bucket depth due to condensation rain        ! RG Add bucket
     depth_change_conv,    &   ! tendency in bucket depth due to convection rain          ! RG Add bucket
     gust,                 &   ! gustiness constant
     z_pbl,                &   ! gustiness constant
     flux_t,               &   ! surface sensible heat flux
     flux_q,               &   ! surface moisture flux
     flux_r,               &   ! surface radiation flux
     flux_u,               &   ! surface flux of zonal mom.
     flux_v,               &   ! surface flux of meridional mom.
     drag_m,               &   ! momentum drag coefficient
     drag_t,               &   ! heat drag coefficient
     drag_q,               &   ! moisture drag coefficient
     w_atm,                &   ! wind speed
     ustar,                &   ! friction velocity
     bstar,                &   ! buoyancy scale
     qstar,                &   ! moisture scale
     dhdt_surf,            &   ! d(sensible heat flux)/d(surface temp)
     dedt_surf,            &   ! d(latent heat flux)/d(surface temp)???
     dedq_surf,            &   ! d(latent heat flux)/d(surface moisture)???
     drdt_surf,            &   ! d(upward longwave)/d(surface temp)
     dhdt_atm,             &   ! d(sensible heat flux)/d(atmos.temp)
     dedq_atm,             &   ! d(latent heat flux)/d(atmospheric mixing rat.)
     dtaudv_atm,           &   ! d(stress component)/d(atmos wind)
     dtaudu_atm,           &   ! d(stress component)/d(atmos wind)
     fracland,             &   ! fraction of land in gridbox
     rough,                &   ! roughness for vert_turb_driver
     albedo,               &   !s albedo now defined in mixed_layer_init
     coszen,               &   !s make sure this is ready for assignment in run_rrtmg
     pbltop,               &   !s Used as an input to damping_driver, outputted from vert_turb_driver
     ex_del_m, 		   &   !mp586 for 10m winds and 2m temp
     ex_del_h,		   &   !mp586 for 10m winds and 2m temp
     ex_del_q,		   &   !mp586 for 10m winds and 2m temp
     temp_2m,		   &   !mp586 for 10m winds and 2m temp
     u_10m,		   &   !mp586 for 10m winds and 2m temp
     v_10m,		   &   !mp586 for 10m winds and 2m temp
     q_2m,                 &   ! Add 2m specific humidity
     rh_2m                     ! Add 2m relative humidity

real, allocatable, dimension(:,:,:) ::                                        &
     diff_m,               &   ! momentum diffusion coeff.
     diff_t,               &   ! temperature diffusion coeff.
     tdtlw,                &   ! place holder. appears in calling arguments of vert_turb_driver but not used unless do_edt=.true. -- pjp
     diss_heat,            &   ! heat dissipated by vertical diffusion
     diss_heat_ray,        &   ! heat dissipated by rayleigh bottom drag (used when gp_surface=.True.)
     non_diff_dt_ug,       &   ! zonal wind tendency except from vertical diffusion
     non_diff_dt_vg,       &   ! merid. wind tendency except from vertical diffusion
     non_diff_dt_tg,       &   ! temperature tendency except from vertical diffusion
     non_diff_dt_qg,       &   ! moisture tendency except from vertical diffusion
     conv_dt_tg,           &   ! temperature tendency from convection
     conv_dt_qg,           &   ! moisture tendency from convection
     cond_dt_tg,           &   ! temperature tendency from condensation
     cond_dt_qg                ! moisture tendency from condensation


logical, allocatable, dimension(:,:) ::                                       &
     avail,                &   ! generate surf. flux (all true)
     land,                 &   ! land points (all false)
     coldT,                &   ! should precipitation be snow at this point
     convect                   ! place holder. appears in calling arguments of vert_turb_driver but not used unless do_entrain=.true. -- pjp

real, allocatable, dimension(:,:) ::                                          &
     land_ones                 ! land points (all zeros)

real, allocatable, dimension(:,:) ::                                          &
     klzbs,                &   ! stored level of zero buoyancy values
     cape,                 &   ! convectively available potential energy
     cin,                  &   ! convective inhibition (this and the above are before the adjustment)
     invtau_q_relaxation,  &   ! temperature relaxation time scale
     invtau_t_relaxation,  &   ! humidity relaxation time scale
     rain,                 &   ! Can be resolved or  parameterised
     snow,                 &   !
     precip                    ! cumulus rain  + resolved rain  + resolved snow

real, allocatable, dimension(:,:,:) :: &
     t_ref,          &   ! relaxation temperature for bettsmiller scheme
     q_ref               ! relaxation moisture for bettsmiller scheme

real, allocatable, dimension(:,:) :: &
     net_surf_sw_down,  &   ! net sw flux at surface
     surf_lw_down           ! downward lw flux at surface

integer ::           &
     id_diff_dt_ug,  &   ! zonal wind tendency from vertical diffusion
     id_diff_dt_vg,  &   ! merid. wind tendency from vertical diffusion
     id_diff_dt_tg,  &   ! temperature tendency from vertical diffusion
     id_diff_dt_qg,  &   ! moisture tendency from vertical diffusion
     id_conv_rain,   &   ! rain from convection
     id_cond_rain,   &   ! rain from condensation
     id_precip,      &   ! rain and snow from condensation and convection
     id_conv_dt_tg,  &   ! temperature tendency from convection
     id_conv_dt_qg,  &   ! temperature tendency from convection
     id_cond_dt_tg,  &   ! temperature tendency from condensation
     id_cond_dt_qg,  &   ! temperature tendency from condensation
     id_bucket_depth,      &   ! bucket depth variable for output  - RG Add bucket
     id_bucket_depth_conv, &   ! bucket depth variation induced by convection  - RG Add bucket
     id_bucket_depth_cond, &   ! bucket depth variation induced by condensation  - RG Add bucket
     id_bucket_depth_lh,   &   ! bucket depth variation induced by LH  - RG Add bucket
     id_rh,          & 	 ! Relative humidity
     id_diss_heat_ray,&  ! Heat dissipated by rayleigh bottom drag if gp_surface=.True.
     id_z_tg,        &   ! Relative humidity
     id_cape,        &
     id_cin,	     & 	     
     id_flux_u,      & ! surface flux of zonal mom.
     id_flux_v,      & ! surface flux of meridional mom.
     id_temp_2m,      & !mp586 for 10m winds and 2m temp
     id_u_10m, 	     & !mp586 for 10m winds and 2m temp
     id_v_10m,       & !mp586 for 10m winds and 2m temp
     id_q_2m,        & ! Add 2m specific humidity
     id_rh_2m          ! Add 2m relative humidity

integer, allocatable, dimension(:,:) :: convflag ! indicates which qe convection subroutines are used
real,    allocatable, dimension(:,:) :: rad_lat, rad_lon
real,    allocatable, dimension(:) :: pref, p_half_1d, ln_p_half_1d, p_full_1d,ln_p_full_1d !s pref is a reference pressure profile, which in 2006 MiMA is just the initial full pressure levels, and an extra level with the reference surface pressure. Others are only necessary to calculate pref.
real,    allocatable, dimension(:,:) :: capeflag !s Added for Betts Miller scheme (rather than the simplified Betts Miller scheme).

type(surf_diff_type) :: Tri_surf ! used by gcm_vert_diff
	
!s initialise constants ready to be used in rh_calc	
real :: d622 = 0.
real :: d378 = 0.
	
logical :: used, doing_edt, doing_entrain, do_strat
integer, dimension(4) :: axes
integer :: is, ie, js, je, num_levels, nsphum, dt_integer
real :: dt_real
type(time_type) :: Time_step

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine idealized_moist_phys_init(Time, Time_step_in, nhum, rad_lon_2d, rad_lat_2d, rad_lonb_2d, rad_latb_2d, t_surf_init)
type(time_type), intent(in) :: Time, Time_step_in
integer, intent(in) :: nhum
real, intent(in), dimension(:,:) :: rad_lon_2d, rad_lat_2d, rad_lonb_2d, rad_latb_2d, t_surf_init

integer :: io, nml_unit, stdlog_unit, seconds, days, id, jd, kd
real, dimension (size(rad_lonb_2d,1)-1, size(rad_latb_2d,2)-1) :: sgsmtn !s added for damping_driver

!s added for land reading
integer, dimension(4) :: siz
integer :: global_num_lon, global_num_lat
character(len=12) :: ctmp1='     by     ', ctmp2='     by     '
!s end added for land reading

! Added for RAS
integer :: num_tracers=0,num_ras_tracers=0,n=0
logical :: do_tracers_in_ras = .false.

logical, dimension(:), allocatable :: tracers_in_ras

character(len=80)  :: scheme
! End Added for RAS

if(module_is_initialized) return

call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=idealized_moist_phys_nml, iostat=io)
#else
   if ( file_exist('input.nml') ) then
      nml_unit = open_namelist_file()
      read (nml_unit, idealized_moist_phys_nml, iostat=io)
      call close_file(nml_unit)
   endif
#endif
stdlog_unit = stdlog()
write(stdlog_unit, idealized_moist_phys_nml)

!s initialise variables for rh_calc
d622 = rdgas/rvgas
d378 = 1.-d622

if(do_cloud_simple) then
  call cloud_simple_init()
end if

!s need to make sure that gray radiation and rrtm radiation are not both called.
if(two_stream_gray .and. do_rrtm_radiation) &
   call error_mesg('physics_driver_init','do_grey_radiation and do_rrtm_radiation cannot both be .true.',FATAL)

if(uppercase(trim(convection_scheme)) == 'NONE') then
  r_conv_scheme = NO_CONV
  lwet_convection = .false.
  do_bm           = .false.
  do_ras          = .false.
  call error_mesg('idealized_moist_phys','No convective adjustment scheme used.', NOTE)

else if(uppercase(trim(convection_scheme)) == 'SIMPLE_BETTS_MILLER') then
  r_conv_scheme = SIMPLE_BETTS_CONV
  call error_mesg('idealized_moist_phys','Using Frierson Quasi-Equilibrium convection scheme.', NOTE)
  lwet_convection = .true.
  do_bm           = .false.
  do_ras          = .false.
  

else if(uppercase(trim(convection_scheme)) == 'FULL_BETTS_MILLER') then
  r_conv_scheme = FULL_BETTS_MILLER_CONV
  call error_mesg('idealized_moist_phys','Using Betts-Miller convection scheme.', NOTE)
  do_bm           = .true.
  lwet_convection = .false.
  do_ras          = .false.
  

else if(uppercase(trim(convection_scheme)) == 'RAS') then
  r_conv_scheme = RAS_CONV
  call error_mesg('idealized_moist_phys','Using relaxed Arakawa Schubert convection scheme.', NOTE)
  do_ras          = .true.
  do_bm           = .false.
  lwet_convection = .false.

else if(uppercase(trim(convection_scheme)) == 'DRY') then
  r_conv_scheme = DRY_CONV
  call error_mesg('idealized_moist_phys','Using dry convection scheme.', NOTE)
  lwet_convection = .false.
  do_bm           = .false.
  do_ras          = .false.  

else if(uppercase(trim(convection_scheme)) == 'UNSET') then
  call error_mesg('idealized_moist_phys','determining convection scheme from flags', NOTE)
  if (lwet_convection) then
    r_conv_scheme = SIMPLE_BETTS_CONV
    call error_mesg('idealized_moist_phys','Using Frierson Quasi-Equilibrium convection scheme.', NOTE)
  end if
  if (do_bm) then
    r_conv_scheme = FULL_BETTS_MILLER_CONV
    call error_mesg('idealized_moist_phys','Using Betts-Miller convection scheme.', NOTE)
  end if
  if (do_ras) then
    r_conv_scheme = RAS_CONV
    call error_mesg('idealized_moist_phys','Using  relaxed Arakawa Schubert convection scheme.', NOTE)
  end if    
else
  call error_mesg('idealized_moist_phys','"'//trim(convection_scheme)//'"'//' is not a valid convection scheme.'// &
      ' Choices are NONE, SIMPLE_BETTS, FULL_BETTS_MILLER, RAS, DRY', FATAL)
endif

if(lwet_convection .and. do_bm) &
  call error_mesg('idealized_moist_phys','lwet_convection and do_bm cannot both be .true.',FATAL)
  
if(lwet_convection .and. do_ras) &
  call error_mesg('idealized_moist_phys','lwet_convection and do_ras cannot both be .true.',FATAL)  

if(do_bm .and. do_ras) &
  call error_mesg('idealized_moist_phys','do_bm and do_ras cannot both be .true.',FATAL)  

nsphum = nhum
Time_step = Time_step_in
call get_time(Time_step, seconds, days)
dt_integer   = 86400*days + seconds
dt_real      = float(dt_integer)

call get_grid_domain(is, ie, js, je)
call get_num_levels(num_levels)

allocate(rad_lat     (is:ie, js:je)); rad_lat = rad_lat_2d
allocate(rad_lon     (is:ie, js:je)); rad_lon = rad_lon_2d
allocate (dt_bucket  (is:ie, js:je)); dt_bucket = 0.0         ! RG Add bucket
allocate (filt       (is:ie, js:je)); filt = 0.0              ! RG Add bucket
allocate(bucket_depth (is:ie, js:je, num_time_levels)); bucket_depth = init_bucket_depth        ! RG Add bucket
allocate(depth_change_lh(is:ie, js:je))                       ! RG Add bucket
allocate(depth_change_cond(is:ie, js:je))                     ! RG Add bucket
allocate(depth_change_conv(is:ie, js:je))                     ! RG Add bucket
allocate(z_surf      (is:ie, js:je))
allocate(t_surf      (is:ie, js:je))
allocate(q_surf      (is:ie, js:je)); q_surf = 0.0
allocate(u_surf      (is:ie, js:je)); u_surf = 0.0
allocate(v_surf      (is:ie, js:je)); v_surf = 0.0
allocate(rough_mom   (is:ie, js:je)); rough_mom = roughness_mom
allocate(rough_heat  (is:ie, js:je)); rough_heat = roughness_heat
allocate(rough_moist (is:ie, js:je)); rough_moist = roughness_moist
allocate(gust        (is:ie, js:je)); gust = 1.0
allocate(z_pbl       (is:ie, js:je))
allocate(flux_t      (is:ie, js:je))
allocate(flux_q      (is:ie, js:je))
allocate(flux_r      (is:ie, js:je))
allocate(flux_u      (is:ie, js:je))
allocate(flux_v      (is:ie, js:je))
allocate(drag_m      (is:ie, js:je))
allocate(drag_t      (is:ie, js:je))
allocate(drag_q      (is:ie, js:je))
allocate(w_atm       (is:ie, js:je))
allocate(ustar       (is:ie, js:je))
allocate(bstar       (is:ie, js:je))
allocate(qstar       (is:ie, js:je))
allocate(dhdt_surf   (is:ie, js:je))
allocate(dedt_surf   (is:ie, js:je))
allocate(dedq_surf   (is:ie, js:je))
allocate(drdt_surf   (is:ie, js:je))
allocate(dhdt_atm    (is:ie, js:je))
allocate(dedq_atm    (is:ie, js:je))
allocate(dtaudv_atm  (is:ie, js:je))
allocate(dtaudu_atm  (is:ie, js:je))
allocate(ex_del_m    (is:ie, js:je)) !mp586 added for 10m wind and 2m temp
allocate(ex_del_h    (is:ie, js:je)) !mp586 added for 10m wind and 2m temp
allocate(ex_del_q    (is:ie, js:je)) !mp586 added for 10m wind and 2m temp
allocate(temp_2m     (is:ie, js:je)) !mp586 added for 10m wind and 2m temp
allocate(u_10m       (is:ie, js:je)) !mp586 added for 10m wind and 2m temp
allocate(v_10m       (is:ie, js:je)) !mp586 added for 10m wind and 2m temp
allocate(q_2m        (is:ie, js:je)) ! Add 2m specific humidity
allocate(rh_2m       (is:ie, js:je)) ! Add 2m relative humidity
allocate(land        (is:ie, js:je)); land = .false.
allocate(land_ones   (is:ie, js:je)); land_ones = 0.0
allocate(avail       (is:ie, js:je)); avail = .true.
allocate(fracland    (is:ie, js:je)); fracland = 0.0
allocate(rough       (is:ie, js:je))
allocate(diff_t      (is:ie, js:je, num_levels))
allocate(diff_m      (is:ie, js:je, num_levels))
allocate(diss_heat   (is:ie, js:je, num_levels))
allocate(diss_heat_ray   (is:ie, js:je, num_levels)) !s added for rayleigh_bottom_drag, used when gp_surface=.True.
allocate(tdtlw       (is:ie, js:je, num_levels)); tdtlw = 0.0

allocate(non_diff_dt_ug  (is:ie, js:je, num_levels))
allocate(non_diff_dt_vg  (is:ie, js:je, num_levels))
allocate(non_diff_dt_tg  (is:ie, js:je, num_levels))
allocate(non_diff_dt_qg  (is:ie, js:je, num_levels))

allocate(net_surf_sw_down        (is:ie, js:je))
allocate(surf_lw_down            (is:ie, js:je))
allocate(conv_dt_tg  (is:ie, js:je, num_levels))
allocate(conv_dt_qg  (is:ie, js:je, num_levels))
allocate(cond_dt_tg  (is:ie, js:je, num_levels))
allocate(cond_dt_qg  (is:ie, js:je, num_levels))

allocate(coldT        (is:ie, js:je)); coldT = .false.
allocate(klzbs        (is:ie, js:je))
allocate(cape         (is:ie, js:je))
allocate(cin          (is:ie, js:je))
allocate(invtau_q_relaxation  (is:ie, js:je))
allocate(invtau_t_relaxation  (is:ie, js:je))
allocate(rain         (is:ie, js:je)); rain = 0.0
allocate(snow         (is:ie, js:je)); snow = 0.0
allocate(precip       (is:ie, js:je)); precip = 0.0
allocate(convflag     (is:ie, js:je))
allocate(convect      (is:ie, js:je)); convect = .false.

allocate(t_ref (is:ie, js:je, num_levels)); t_ref = 0.0
allocate(q_ref (is:ie, js:je, num_levels)); q_ref = 0.0

allocate (albedo      (is:ie, js:je)) !s allocate for albedo, to be set in mixed_layer_init.
allocate(coszen       (is:ie, js:je)) !s allocate coszen to be set in run_rrtmg
allocate(pbltop       (is:ie, js:je)) !s allocate coszen to be set in run_rrtmg

allocate(pref(num_levels+1)) !s reference pressure profile, as in spectral_physics.f90 in FMS 2006 and original MiMA.
allocate(p_half_1d(num_levels+1), ln_p_half_1d(num_levels+1))
allocate(p_full_1d(num_levels  ), ln_p_full_1d(num_levels  ))
allocate(capeflag     (is:ie, js:je))

call get_surf_geopotential(z_surf)
z_surf = z_surf/grav

!s initialise the land area
if(trim(land_option) .eq. 'input')then
!s read in land nc file
!s adapted from spectral_init_cond.F90

	   if(file_exist(trim(land_file_name))) then
	     call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat)
	     call field_size(trim(land_file_name), trim(land_field_name), siz)
	     if ( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) then
	       call read_data(trim(land_file_name), trim(land_field_name), land_ones, grid_domain)
	       !s write something to screen to let the user know what's happening.
	     else
	       write(ctmp1(1: 4),'(i4)') siz(1)
	       write(ctmp1(9:12),'(i4)') siz(2)
	       write(ctmp2(1: 4),'(i4)') global_num_lon
	       write(ctmp2(9:12),'(i4)') global_num_lat
	       call error_mesg ('idealized_moist_phys','Land file contains data on a '// &
	              ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
	     endif
	   else
	     call error_mesg('idealized_moist_phys','land_option="'//trim(land_option)//'"'// &
	                     ' but '//trim(land_file_name)//' does not exist', FATAL)
	   endif

	!s convert data in land nc file to land logical array
	where(land_ones > 0.) land = .true.

elseif(trim(land_option) .eq. 'zsurf')then
	!s wherever zsurf is greater than some threshold height then make land = .true.
	where ( z_surf > 10. ) land = .true.
endif


!s Add option to alter surface roughness length over land

if(trim(land_option) .eq. 'input') then

	where(land)
	rough_mom   = land_roughness_prefactor * rough_mom
	rough_heat  = land_roughness_prefactor * rough_heat
	rough_moist = land_roughness_prefactor * rough_moist
	end where

endif

!RG Add bucket - initialise bucket depth
if(bucket) then
where(land)
  bucket_depth(:,:,1)  = init_bucket_depth_land
  bucket_depth(:,:,2)  = init_bucket_depth_land
end where
endif
!RG end Add bucket

!s end option to alter surface roughness length over land


if (gp_surface) then
call rayleigh_bottom_drag_init(get_axis_id(), Time)
axes = get_axis_id()
id_diss_heat_ray = register_diag_field(mod_name, 'diss_heat_ray', &
                   axes(1:3), Time, 'dissipated heat from Rayleigh drag', 'K/s')
endif


!    initialize damping_driver_mod.
      if(do_damping) then
         call pressure_variables(p_half_1d,ln_p_half_1d,pref(1:num_levels),ln_p_full_1d,PSTD_MKS)
	 pref(num_levels+1) = PSTD_MKS
         call damping_driver_init (rad_lonb_2d(:,1),rad_latb_2d(1,:), pref(:), get_axis_id(), Time, & !s note that in the original this is pref(:,1), which is the full model pressure levels and the surface pressure at the bottom. There is pref(:2) in this version with 81060 as surface pressure??
                                sgsmtn)

      endif

if(mixed_layer_bc) then
  ! need an initial condition for the mixed layer temperature
  ! may be overwritten by restart file
  ! choose an unstable initial condition to allow moisture
  ! to quickly enter the atmosphere avoiding problems with the convection scheme
  t_surf = t_surf_init + 1.0

  call mixed_layer_init(is, ie, js, je, num_levels, t_surf, bucket_depth, get_axis_id(), Time, albedo, rad_lonb_2d(:,:), rad_latb_2d(:,:), land, bucket) ! t_surf is intent(inout) !s albedo distribution set here.
  
elseif(gp_surface) then
  albedo=0.0
  call error_mesg('idealized_moist_phys','Because gp_surface=.True., setting albedo=0.0', NOTE)

  call error_mesg('idealized_moist_phys','Note that if grey radiation scheme != Schneider is used, model will seg-fault b/c gp_surface does not define a t_surf, which is required by most grey schemes.', NOTE)

endif

if(turb) then
! need to call vert_diff_init even if using gcm_vert_diff (rather than
! gcm_vert_diff_down) because the variable sphum is not initialized
! otherwise in the vert_diff module
   call vert_diff_init (Tri_surf, ie-is+1, je-js+1, num_levels, .true., do_virtual) !s do_conserve_energy is hard-coded in.
end if

call lscale_cond_init()

axes = get_axis_id()

id_cond_dt_qg = register_diag_field(mod_name, 'dt_qg_condensation',        &
     axes(1:3), Time, 'Moisture tendency from condensation','kg/kg/s')
id_cond_dt_tg = register_diag_field(mod_name, 'dt_tg_condensation',        &
     axes(1:3), Time, 'Temperature tendency from condensation','K/s')
id_cond_rain = register_diag_field(mod_name, 'condensation_rain',          &
     axes(1:2), Time, 'Rain from condensation','kg/m/m/s')
id_precip = register_diag_field(mod_name, 'precipitation',          &
     axes(1:2), Time, 'Precipitation from resolved, parameterised and snow','kg/m/m/s')
id_cape = register_diag_field(mod_name, 'cape',          &
     axes(1:2), Time, 'Convective Available Potential Energy','J/kg')
id_cin = register_diag_field(mod_name, 'cin',          &
     axes(1:2), Time, 'Convective Inhibition','J/kg')
id_flux_u = register_diag_field(mod_name, 'flux_u', &
     axes(1:2), Time, 'Zonal momentum flux', 'Pa')
id_flux_v = register_diag_field(mod_name, 'flux_v', &
     axes(1:2), Time, 'Meridional momentum flux', 'Pa')

if(bucket) then
  id_bucket_depth = register_diag_field(mod_name, 'bucket_depth',            &         ! RG Add bucket
       axes(1:2), Time, 'Depth of surface reservoir', 'm')
  id_bucket_depth_conv = register_diag_field(mod_name, 'bucket_depth_conv',  &         ! RG Add bucket
       axes(1:2), Time, 'Tendency of bucket depth induced by Convection', 'm/s')
  id_bucket_depth_cond = register_diag_field(mod_name, 'bucket_depth_cond',  &         ! RG Add bucket
       axes(1:2), Time, 'Tendency of bucket depth induced by Condensation', 'm/s')
  id_bucket_depth_lh = register_diag_field(mod_name, 'bucket_depth_lh',      &         ! RG Add bucket
       axes(1:2), Time, 'Tendency of bucket depth induced by LH', 'm/s')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! added by mp586 for 10m winds and 2m temperature add mo_profile()!!!!!!!!

id_temp_2m = register_diag_field(mod_name, 'temp_2m',            &         !mp586 add 2m temp
     axes(1:2), Time, 'Air temperature 2m above surface', 'K')
id_u_10m = register_diag_field(mod_name, 'u_10m',                &         !mp586 add 10m wind (u)
     axes(1:2), Time, 'Zonal wind 10m above surface', 'm/s')
id_v_10m = register_diag_field(mod_name, 'v_10m',                &         !mp586 add 10m wind (v)
     axes(1:2), Time, 'Meridional wind 10m above surface', 'm/s')

!!!!!!!!!!!! end of mp586 additions !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

id_q_2m = register_diag_field(mod_name, 'sphum_2m',                  &
     axes(1:2), Time, 'Specific humidity 2m above surface', 'kg/kg')       !Add 2m specific humidity
id_rh_2m = register_diag_field(mod_name, 'rh_2m',                &
     axes(1:2), Time, 'Relative humidity 2m above surface', 'percent')     !Add 2m relative humidity

select case(r_conv_scheme)

case(SIMPLE_BETTS_CONV)
  call qe_moist_convection_init()

case(FULL_BETTS_MILLER_CONV)
  call betts_miller_init()

case(DRY_CONV)
  call dry_convection_init(axes, Time)

case(RAS_CONV)

        !run without startiform cloud scheme

       !---------------------------------------------------------------------
       !    retrieve the number of registered tracers in order to determine 
       !    which tracers are to be convectively transported.
       !---------------------------------------------------------------------

       call get_number_tracers (MODEL_ATMOS, num_tracers= num_tracers)

       allocate (tracers_in_ras(num_tracers))
       ! Instead of finding out which of the tracers need to be advected, manually set this to .false.
       tracers_in_ras = .false.
       do_strat = .false.

       !Commented code not used such that tracers are not advected by RAS. Could implement in future.
       
       ! do n=1, num_tracers
       !   if (query_method ('convection', MODEL_ATMOS, n, scheme)) then
       !    num_ras_tracers = num_ras_tracers + 1
       !    tracers_in_ras(n) = .true.
       !   endif
       ! end do

       ! if (num_ras_tracers > 0) then
       !   do_tracers_in_ras = .true.
       ! else
       !   do_tracers_in_ras = .false.
       ! endif

       !----------------------------------------------------------------------
       !    for each tracer, determine if it is to be transported by convect-
       !    ion, and the convection schemes that are to transport it. set a 
       !    logical flag to .true. for each tracer that is to be transported by
       !    each scheme and increment the count of tracers to be transported
       !    by that scheme.
       !----------------------------------------------------------------------

        call ras_init (do_strat, axes,Time,tracers_in_ras) 

end select

!jp not sure why these diag_fields are fenced when condensation ones above are not...
!if(lwet_convection .or. do_bm) then
   id_conv_dt_qg = register_diag_field(mod_name, 'dt_qg_convection',          &
        axes(1:3), Time, 'Moisture tendency from convection','kg/kg/s')
   id_conv_dt_tg = register_diag_field(mod_name, 'dt_tg_convection',          &
        axes(1:3), Time, 'Temperature tendency from convection','K/s')
   id_conv_rain = register_diag_field(mod_name, 'convection_rain',            &
        axes(1:2), Time, 'Rain from convection','kg/m/m/s')
!endif


if(two_stream_gray) call two_stream_gray_rad_init(is, ie, js, je, num_levels, get_axis_id(), Time, rad_lonb_2d, rad_latb_2d, dt_real)

#ifdef RRTM_NO_COMPILE
    if (do_rrtm_radiation) then
        call error_mesg('idealized_moist_phys','do_rrtm_radiation is .true. but compiler flag -D RRTM_NO_COMPILE used. Stopping.', FATAL)
    endif
#else
    if(do_rrtm_radiation) then
       id=ie-is+1 !s Taking dimensions from equivalend calls in vert_turb_driver_init
       jd=je-js+1
       kd=num_levels
       call rrtmg_lw_ini(cp_air)
       call rrtmg_sw_ini(cp_air)
       call rrtm_radiation_init(axes,Time,id*jd,kd,rad_lonb_2d,rad_latb_2d, Time_step_in)
       id_z_tg = register_diag_field(mod_name, 'interp_t',        &
            axes(1:3), Time, 'temperature interp','T/s')
    endif
#endif

#ifdef SOC_NO_COMPILE
    if (do_socrates_radiation) then
        call error_mesg('idealized_moist_phys','do_socrates_radiation is .true. but compiler flag -D SOC_NO_COMPILE used. Stopping.', FATAL)
    endif
#else
if (do_socrates_radiation) then
    call socrates_init(is, ie, js, je, num_levels, axes, Time, rad_lat, rad_lonb_2d, rad_latb_2d, Time_step_in)
endif
#endif

if(turb) then
   call vert_turb_driver_init (rad_lonb_2d, rad_latb_2d, ie-is+1,je-js+1, &
                 num_levels,get_axis_id(),Time, doing_edt, doing_entrain)

   axes = get_axis_id()
   id_diff_dt_ug = register_diag_field(mod_name, 'dt_ug_diffusion',        &
        axes(1:3), Time, 'zonal wind tendency from diffusion','m/s^2')
   id_diff_dt_vg = register_diag_field(mod_name, 'dt_vg_diffusion',        &
        axes(1:3), Time, 'meridional wind tendency from diffusion','m/s^2')
   id_diff_dt_tg = register_diag_field(mod_name, 'dt_tg_diffusion',        &
        axes(1:3), Time, 'temperature diffusion tendency','K/s')
   id_diff_dt_qg = register_diag_field(mod_name, 'dt_qg_diffusion',        &
        axes(1:3), Time, 'moisture diffusion tendency','kg/kg/s')
endif

   id_rh = register_diag_field ( mod_name, 'rh', &
	axes(1:3), Time, 'relative humidity', 'percent')

end subroutine idealized_moist_phys_init
!=================================================================================================================================
subroutine idealized_moist_phys(Time, p_half, p_full, z_half, z_full, ug, vg, tg, grid_tracers, &
                                previous, current, dt_ug, dt_vg, dt_tg, dt_tracers, mask, kbot)

type(time_type),            intent(in)    :: Time
real, dimension(:,:,:,:),   intent(in)    :: p_half, p_full, z_half, z_full, ug, vg, tg
real, dimension(:,:,:,:,:), intent(in)    :: grid_tracers
integer,                    intent(in)    :: previous, current
real, dimension(:,:,:),     intent(inout) :: dt_ug, dt_vg, dt_tg
real, dimension(:,:,:,:),   intent(inout) :: dt_tracers

real :: delta_t
real, dimension(size(ug,1), size(ug,2), size(ug,3)) :: tg_tmp, qg_tmp, RH,tg_interp, mc, dt_ug_conv, dt_vg_conv

! Simple cloud scheme variabilies to pass to radiation
real, dimension(size(ug,1), size(ug,2), size(ug,3))    :: cfa_rad, reff_rad

real, intent(in) , dimension(:,:,:), optional :: mask
integer, intent(in) , dimension(:,:),   optional :: kbot

real, dimension(1,1,1):: tracer, tracertnd
integer :: nql, nqi, nqa   ! tracer indices for stratiform clouds

if(current == previous) then
   delta_t = dt_real
else
   delta_t = 2*dt_real
endif

if (bucket) then
  dt_bucket = 0.0                ! RG Add bucket
  filt      = 0.0                ! RG Add bucket
endif

rain = 0.0; snow = 0.0; precip = 0.0

select case(r_conv_scheme)

case(SIMPLE_BETTS_CONV)

   call qe_moist_convection ( delta_t,              tg(:,:,:,previous),      &
    grid_tracers(:,:,:,previous,nsphum),        p_full(:,:,:,previous),      &
                          p_half(:,:,:,previous),                coldT,      &
                                 rain,                            snow,      &
                           conv_dt_tg,                      conv_dt_qg,      &
                                q_ref,                        convflag,      &
                                klzbs,                            cape,      &
                                  cin,             invtau_q_relaxation,      &
                  invtau_t_relaxation,                           t_ref)

   tg_tmp = conv_dt_tg + tg(:,:,:,previous)
   qg_tmp = conv_dt_qg + grid_tracers(:,:,:,previous,nsphum)
!  note the delta's are returned rather than the time derivatives

   conv_dt_tg = conv_dt_tg/delta_t
   conv_dt_qg = conv_dt_qg/delta_t
   depth_change_conv = rain/dens_h2o     ! RG Add bucket
   rain       = rain/delta_t
   precip     = rain

   if(id_conv_dt_qg > 0) used = send_data(id_conv_dt_qg, conv_dt_qg, Time)
   if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
   if(id_conv_rain  > 0) used = send_data(id_conv_rain, rain, Time)
   if(id_cape  > 0) used = send_data(id_cape, cape, Time)
   if(id_cin  > 0) used = send_data(id_cin, cin, Time)

case(FULL_BETTS_MILLER_CONV)

   call betts_miller (          delta_t,           tg(:,:,:,previous),       &
    grid_tracers(:,:,:,previous,nsphum),       p_full(:,:,:,previous),       &
                 p_half(:,:,:,previous),                        coldT,       &
                                   rain,                         snow,       &
                             conv_dt_tg,                   conv_dt_qg,       &
                                  q_ref,                     convflag,       &
                                  klzbs,                         cape,       &
                                    cin,                        t_ref,       &
                    invtau_t_relaxation,          invtau_q_relaxation,       &
                               capeflag)

   tg_tmp = conv_dt_tg + tg(:,:,:,previous)
   qg_tmp = conv_dt_qg + grid_tracers(:,:,:,previous,nsphum)
!  note the delta's are returned rather than the time derivatives

   conv_dt_tg = conv_dt_tg/delta_t
   conv_dt_qg = conv_dt_qg/delta_t
   depth_change_conv = rain/dens_h2o     ! RG Add bucket
   rain       = rain/delta_t
   precip     = rain

   if(id_conv_dt_qg > 0) used = send_data(id_conv_dt_qg, conv_dt_qg, Time)
   if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
   if(id_conv_rain  > 0) used = send_data(id_conv_rain, rain, Time)
   if(id_cape  > 0) used = send_data(id_cape, cape, Time)
   if(id_cin  > 0) used = send_data(id_cin, cin, Time)

case(DRY_CONV)
    call dry_convection(Time, tg(:, :, :, previous),                         &
                        p_full(:,:,:,previous), p_half(:,:,:,previous),      &
                        conv_dt_tg, cape, cin)

    tg_tmp = conv_dt_tg*delta_t + tg(:,:,:,previous)
    qg_tmp = grid_tracers(:,:,:,previous,nsphum)

    if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
    if(id_cape  > 0) used = send_data(id_cape, cape, Time)
    if(id_cin  > 0) used = send_data(id_cin, cin, Time)

case(RAS_CONV)

    call ras   (is,   js,     Time,                                                  &  
                tg(:,:,:,previous),   grid_tracers(:,:,:,previous,nsphum),           &
                ug(:,:,:,previous),  vg(:,:,:,previous),    p_full(:,:,:,previous),  &
                p_half(:,:,:,previous), z_half(:,:,:,previous), coldT,  delta_t,     &
                conv_dt_tg,   conv_dt_qg, dt_ug_conv,  dt_vg_conv,                   &
                rain, snow,   do_strat,                                              &                                              
                !OPTIONAL 
                mask,  kbot,                                                         &
                !OPTIONAL OUT
                mc,   tracer(:,:,:), tracer(:,:,:),                          &
               tracer(:,:,:),  tracertnd(:,:,:),                             &
               tracertnd(:,:,:), tracertnd(:,:,:))
                

      !update tendencies - dT and dq are done after cases
      tg_tmp = tg(:,:,:,previous) + conv_dt_tg
      qg_tmp = grid_tracers(:,:,:,previous,nsphum) + conv_dt_qg
      dt_ug = dt_ug + dt_ug_conv
      dt_vg = dt_vg + dt_vg_conv

      precip     = precip + rain + snow

   if(id_conv_dt_qg > 0) used = send_data(id_conv_dt_qg, conv_dt_qg, Time)
   if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
   if(id_conv_rain  > 0) used = send_data(id_conv_rain, precip, Time)


case(NO_CONV)
   tg_tmp = tg(:,:,:,previous)
   qg_tmp = grid_tracers(:,:,:,previous,nsphum)

case default
  call error_mesg('idealized_moist_phys','Invalid convection scheme.', FATAL)

end select


! Add the T and q tendencies due to convection to the timestep
dt_tg = dt_tg + conv_dt_tg
dt_tracers(:,:,:,nsphum) = dt_tracers(:,:,:,nsphum) + conv_dt_qg


! Perform large scale convection
if (r_conv_scheme .ne. DRY_CONV) then
  ! Large scale convection is a function of humidity only.  This is
  ! inconsistent with the dry convection scheme, don't run it!
  rain = 0.0; snow = 0.0
  call lscale_cond (         tg_tmp,                          qg_tmp,        &
             p_full(:,:,:,previous),          p_half(:,:,:,previous),        &
                              coldT,                            rain,        &
                               snow,                      cond_dt_tg,        &
                         cond_dt_qg )

  cond_dt_tg = cond_dt_tg/delta_t
  cond_dt_qg = cond_dt_qg/delta_t
  depth_change_cond = rain/dens_h2o     ! RG Add bucket
  rain       = rain/delta_t
  snow       = snow/delta_t
  precip     = precip + rain + snow

  dt_tg = dt_tg + cond_dt_tg
  dt_tracers(:,:,:,nsphum) = dt_tracers(:,:,:,nsphum) + cond_dt_qg

  if(id_cond_dt_qg > 0) used = send_data(id_cond_dt_qg, cond_dt_qg, Time)
  if(id_cond_dt_tg > 0) used = send_data(id_cond_dt_tg, cond_dt_tg, Time)
  if(id_cond_rain  > 0) used = send_data(id_cond_rain, rain, Time)
  if(id_precip     > 0) used = send_data(id_precip, precip, Time)

endif

! Call the simple cloud scheme in line with SPOOKIE-2 requirements
! Using start of time step variables
! using soecific humidity NOT mixing ratios
if(do_cloud_simple) then

    call cloud_simple(p_half(:,:,:,current), p_full(:,:,:,current),  &
                      Time,                                &
                      tg(:,:,:,previous),                  &
                      grid_tracers(:,:,:,previous,nsphum), &
                      ! outs - 
                      cfa_rad(:,:,:), reff_rad(:,:,:)      & 
                      )

endif


! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

if(two_stream_gray) then
   call two_stream_gray_rad_down(is, js, Time, &
                       rad_lat(:,:),           &
                       rad_lon(:,:),           &
                       p_half(:,:,:,current),  &
                       tg(:,:,:,previous),     &
                       net_surf_sw_down(:,:),  &
                       surf_lw_down(:,:), albedo, &
                       grid_tracers(:,:,:,previous,nsphum))
end if

if(.not.mixed_layer_bc) then

!!$! infinite heat capacity
!    t_surf = surface_temperature_forced(rad_lat)
!!$! no heat capacity:
!!$   t_surf = tg(:,:,num_levels,previous)

!!$! surface temperature has same potential temp. as lowest layer:
!!$  t_surf = surface_temperature(tg(:,:,:,previous), p_full(:,:,:,current), p_half(:,:,:,current))
end if


if(.not.gp_surface) then 
  call surface_flux(                                                        &
                  tg(:,:,num_levels,previous),                              &
 grid_tracers(:,:,num_levels,previous,nsphum),                              &
                  ug(:,:,num_levels,previous),                              &
                  vg(:,:,num_levels,previous),                              &
               p_full(:,:,num_levels,current),                              &
   z_full(:,:,num_levels,current)-z_surf(:,:),                              &
             p_half(:,:,num_levels+1,current),                              &
                                  t_surf(:,:),                              &
                                  t_surf(:,:),                              &
                                  q_surf(:,:),                              & ! is intent(inout)
                                       bucket,                              &     ! RG Add bucket
                    bucket_depth(:,:,current),                              &     ! RG Add bucket
                        max_bucket_depth_land,                              &     ! RG Add bucket
                         depth_change_lh(:,:),                              &     ! RG Add bucket
                       depth_change_conv(:,:),                              &     ! RG Add bucket
                       depth_change_cond(:,:),                              &     ! RG Add bucket
                                  u_surf(:,:),                              &
                                  v_surf(:,:),                              &
                               rough_mom(:,:),                              &
                              rough_heat(:,:),                              &
                             rough_moist(:,:),                              &
                               rough_mom(:,:),                              & ! using rough_mom in place of rough_scale -- pjp
                                    gust(:,:),                              &
                                  flux_t(:,:),                              & ! is intent(out)
                                  flux_q(:,:),                              & ! is intent(out)
                                  flux_r(:,:),                              & ! is intent(out)
                                  flux_u(:,:),                              & ! is intent(out)
                                  flux_v(:,:),                              & ! is intent(out)
                                  drag_m(:,:),                              & ! is intent(out)
                                  drag_t(:,:),                              & ! is intent(out)
                                  drag_q(:,:),                              & ! is intent(out)
                                   w_atm(:,:),                              & ! is intent(out)
                                   ustar(:,:),                              & ! is intent(out)
                                   bstar(:,:),                              & ! is intent(out)
                                   qstar(:,:),                              & ! is intent(out)
                               dhdt_surf(:,:),                              & ! is intent(out)
                               dedt_surf(:,:),                              & ! is intent(out)
                               dedq_surf(:,:),                              & ! is intent(out)
                               drdt_surf(:,:),                              & ! is intent(out)
                                dhdt_atm(:,:),                              & ! is intent(out)
                                dedq_atm(:,:),                              & ! is intent(out)
                              dtaudu_atm(:,:),                              & ! is intent(out)
                              dtaudv_atm(:,:),                              & ! is intent(out)
			        ex_del_m(:,:),				    & ! mp586 for 10m winds and 2m temp
			        ex_del_h(:,:),				    & ! mp586 for 10m winds and 2m temp
			        ex_del_q(:,:),				    & ! mp586 for 10m winds and 2m temp
			         temp_2m(:,:),				    & ! mp586 for 10m winds and 2m temp
			           u_10m(:,:),				    & ! mp586 for 10m winds and 2m temp	
			           v_10m(:,:),				    & ! mp586 for 10m winds and 2m temp
                                    q_2m(:,:),                              & ! Add 2m specific humidity
                                   rh_2m(:,:),                              & ! Add 2m relative humidity
                 	              delta_t,                              &
                                    land(:,:),                              &
                               .not.land(:,:),                              &
                                   avail(:,:)  )

  if(id_flux_u > 0) used = send_data(id_flux_u, flux_u, Time)
  if(id_flux_v > 0) used = send_data(id_flux_v, flux_v, Time)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! added by mp586 for 10m winds and 2m temperature add mo_profile()!!!!!!!!

  if(id_temp_2m > 0) used = send_data(id_temp_2m, temp_2m, Time)    ! mp586 add 2m temp
  if(id_u_10m > 0) used = send_data(id_u_10m, u_10m, Time)          ! mp586 add 10m wind (u)
  if(id_v_10m > 0) used = send_data(id_v_10m, v_10m, Time)          ! mp586 add 10m wind (v)


  !!!!!!!!!!!! end of mp586 additions !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(id_q_2m > 0) used = send_data(id_q_2m, q_2m, Time)           ! Add 2m specific humidity
  if(id_rh_2m > 0) used = send_data(id_rh_2m, rh_2m*1e2, Time)    ! Add 2m relative humidity

endif

! Now complete the radiation calculation by computing the upward and net fluxes.

if(two_stream_gray) then
   call two_stream_gray_rad_up(is, js, Time, &
                     rad_lat(:,:),           &
                     p_half(:,:,:,current),  &
                     t_surf(:,:),            &
                     tg(:,:,:,previous),     &
                     dt_tg(:,:,:), albedo)
end if

#ifdef RRTM_NO_COMPILE
    if (do_rrtm_radiation) then
        call error_mesg('idealized_moist_phys','do_rrtm_radiation is .true. but compiler flag -D RRTM_NO_COMPILE used. Stopping.', FATAL)
    endif
#else
if(do_rrtm_radiation) then
   !need t at half grid
	tg_interp=tg(:,:,:,previous)
   call interp_temp(z_full(:,:,:,current),z_half(:,:,:,current),tg_interp, Time)
   call run_rrtmg(is,js,Time,rad_lat(:,:),rad_lon(:,:),p_full(:,:,:,current),p_half(:,:,:,current),  &
                  albedo,grid_tracers(:,:,:,previous,nsphum),tg_interp,t_surf(:,:),dt_tg(:,:,:),     &
                  coszen,net_surf_sw_down(:,:),surf_lw_down(:,:), cfa_rad(:,:,:), reff_rad(:,:,:),   &     
                  do_cloud_simple )
endif
#endif

#ifdef SOC_NO_COMPILE
    if (do_socrates_radiation) then
        call error_mesg('idealized_moist_phys','do_socrates_radiation is .true. but compiler flag -D SOC_NO_COMPILE used. Stopping.', FATAL)
    endif
#else
if (do_socrates_radiation) then
       ! Socrates interface
       
    call run_socrates(Time, Time+Time_step, rad_lat, rad_lon, tg(:,:,:,previous), grid_tracers(:,:,:,previous,nsphum), t_surf(:,:), p_full(:,:,:,current), &
                      p_half(:,:,:,current),z_full(:,:,:,current),z_half(:,:,:,current), albedo, dt_tg(:,:,:), net_surf_sw_down(:,:), surf_lw_down(:,:), delta_t)

endif
#endif

if(gp_surface) then

	call gp_surface_flux (dt_tg(:,:,:), p_half(:,:,:,current), num_levels)
	
    call compute_rayleigh_bottom_drag( 1,                     ie-is+1, &
                                       1,                     je-js+1, &
                                     Time,                    delta_t, &
                   		     rad_lat(:,:),         dt_ug(:,:,:      ), &
                        dt_vg(:,:,:     ),                             &
                       ug(:,:,:,previous),         vg(:,:,:,previous), &
                     p_half(:,:,:,previous),     p_full(:,:,:,previous), &
                     dt_tg, diss_heat_ray )

	if(id_diss_heat_ray > 0) used = send_data(id_diss_heat_ray, diss_heat_ray, Time)
endif




!----------------------------------------------------------------------
!    Copied from MiMA physics_driver.f90
!    call damping_driver to calculate the various model dampings that
!    are desired.
!----------------------------------------------------------------------
z_pbl(:,:) = pbltop(is:ie,js:je)
if(do_damping) then
     call damping_driver (is, js, rad_lat, Time+Time_step, delta_t,                               &
                             p_full(:,:,:,current), p_half(:,:,:,current),              &
                             z_full(:,:,:,current), z_half(:,:,:,current),              &
                             ug(:,:,:,previous), vg(:,:,:,previous),                    &
                             tg(:,:,:,previous), grid_tracers(:,:,:,previous,nsphum),   &
                             grid_tracers(:,:,:,previous,:),                            &
                             dt_ug(:,:,:), dt_vg(:,:,:), dt_tg(:,:,:),                  &
                             dt_tracers(:,:,:,nsphum), dt_tracers(:,:,:,:),             &
                             z_pbl) !s have taken the names of arrays etc from vert_turb_driver below. Watch ntp from 2006 call to this routine?
endif




if(turb) then

   call vert_turb_driver(            1,                              1, &
                                  Time,                 Time+Time_step, &
                               delta_t, tdtlw(:,:,:),    fracland(:,:), &
                 p_half(:,:,:,current),          p_full(:,:,:,current), &
                 z_half(:,:,:,current),          z_full(:,:,:,current), &
                            ustar(:,:),                     bstar(:,:), &
                            qstar(:,:),                     rough(:,:), &
                          rad_lat(:,:),                   convect(:,:), &
                    ug(:,:,:,current ),             vg(:,:,:,current ), &
                    tg(:,:,:,current ),                                 &
    grid_tracers(:,:,:,current,nsphum),  grid_tracers(:,:,:,current,:), &
                    ug(:,:,:,previous),                                 &
                    vg(:,:,:,previous),             tg(:,:,:,previous), &
   grid_tracers(:,:,:,previous,nsphum), grid_tracers(:,:,:,previous,:), &
                          dt_ug(:,:,:),                   dt_vg(:,:,:), &
                          dt_tg(:,:,:),       dt_tracers(:,:,:,nsphum), &
                   dt_tracers(:,:,:,:),                  diff_t(:,:,:), &
                         diff_m(:,:,:),                      gust(:,:), &
                            z_pbl(:,:) )

      pbltop(is:ie,js:je) = z_pbl(:,:) !s added so that z_pbl can be used subsequently by damping_driver.

!
!! Don't zero these derivatives as the surface flux depends implicitly
!! on the lowest level values
!! However it should be noted that these derivatives do not take into
!! account the change in the Monin-Obukhov coefficients, and so are not
!! very accurate.
!
!!$   dtaudv_atm = 0.0
!!$   dhdt_atm   = 0.0
!!$   dedq_atm   = 0.0

   if(.not.(mixed_layer_bc.or.gp_surface)) then
     call error_mesg('atmosphere','no diffusion implentation for non-mixed layer b.c.',FATAL)
   endif


! We must use gcm_vert_diff_down and _up rather than gcm_vert_diff as the surface flux
! depends implicitly on the surface values

!
! Don't want to do time splitting for the implicit diffusion step in case
! of compensation of the tendencies
!
   non_diff_dt_ug  = dt_ug
   non_diff_dt_vg  = dt_vg
   non_diff_dt_tg  = dt_tg
   non_diff_dt_qg  = dt_tracers(:,:,:,nsphum)

   call gcm_vert_diff_down (1, 1,                                          &
                            delta_t,             ug(:,:,:,previous),       &
                            vg(:,:,:,previous),  tg(:,:,:,previous),       &
                            grid_tracers(:,:,:,previous,nsphum),           &
                            grid_tracers(:,:,:,previous,:), diff_m(:,:,:), &
                            diff_t(:,:,:),          p_half(:,:,:,current), &
                            p_full(:,:,:,current),  z_full(:,:,:,current), &
                            flux_u(:,:),                      flux_v(:,:), &
                            dtaudu_atm(:,:),              dtaudv_atm(:,:), &
                            dt_ug(:,:,:),                    dt_vg(:,:,:), &
                            dt_tg(:,:,:),        dt_tracers(:,:,:,nsphum), &
                            dt_tracers(:,:,:,:),         diss_heat(:,:,:), &
                            Tri_surf)
!
! update surface temperature
!
   if(mixed_layer_bc) then	
   call mixed_layer(                                                       &
                              Time, Time+Time_step,                        &
                              js,                                          & 
                              je,                                          &
                              t_surf(:,:),                                 & ! t_surf is intent(inout)
                              flux_t(:,:),                                 &
                              flux_q(:,:),                                 &
                              flux_r(:,:),                                 &
                                  dt_real,                                 &
                    net_surf_sw_down(:,:),                                 &
                        surf_lw_down(:,:),                                 &
                            Tri_surf,                                      & ! Tri_surf is intent(inout)
                           dhdt_surf(:,:),                                 &
                           dedt_surf(:,:),                                 &
                           dedq_surf(:,:),                                 &
                           drdt_surf(:,:),                                 &
                            dhdt_atm(:,:),                                 &
                            dedq_atm(:,:),                                 &
                              albedo(:,:))
   endif

   call gcm_vert_diff_up (1, 1, delta_t, Tri_surf, dt_tg(:,:,:), dt_tracers(:,:,:,nsphum), dt_tracers(:,:,:,:))

   if(id_diff_dt_ug > 0) used = send_data(id_diff_dt_ug, dt_ug - non_diff_dt_ug, Time)
   if(id_diff_dt_vg > 0) used = send_data(id_diff_dt_vg, dt_vg - non_diff_dt_vg, Time)
   if(id_diff_dt_tg > 0) used = send_data(id_diff_dt_tg, dt_tg - non_diff_dt_tg, Time)
   if(id_diff_dt_qg > 0) used = send_data(id_diff_dt_qg, dt_tracers(:,:,:,nsphum) - non_diff_dt_qg, Time)

endif ! if(turb) then

!s Adding relative humidity calculation so as to allow comparison with Frierson's thesis.
   call rh_calc (p_full(:,:,:,previous),tg_tmp,qg_tmp,RH)
   if(id_rh >0) used = send_data(id_rh, RH*100., Time)


! RG Add bucket
! Timestepping for bucket. 
! NB In tapios github, all physics is still in atmosphere.F90 and this leapfrogging is done there. 
!This part has been included here to avoid editing atmosphere.F90
! Therefore define a future variable locally, but do not feedback any changes to timestepping variables upstream, so as to avoid messing with the model's overall timestepping.
! Bucket diffusion has been cut for this version - could be incorporated later.

if(bucket) then

  if(previous == current) then
    future = num_time_levels + 1 - current
  else
    future = previous
  endif

   ! bucket time tendency
   dt_bucket = depth_change_cond + depth_change_conv - depth_change_lh
   !change in bucket depth in one leapfrog timestep [m]                                 

   ! use the raw filter in leapfrog time stepping

   filt(:,:) = bucket_depth(:,:,previous) - 2.0 * bucket_depth(:,:,current)

   if(previous == current) then
      bucket_depth(:,:,future ) = bucket_depth(:,:,previous) + dt_bucket
      bucket_depth(:,:,current) = bucket_depth(:,:,current ) + robert_bucket &
        *(bucket_depth(:,:,previous) - 2.0*bucket_depth(:,:,current) + bucket_depth(:,:,future)) * raw_bucket
   else
      bucket_depth(:,:,current) = bucket_depth(:,:,current ) + robert_bucket &
        *(bucket_depth(:,:,previous) - 2.0*bucket_depth(:,:,current)) * raw_bucket 
      bucket_depth(:,:,future ) = bucket_depth(:,:,previous) + dt_bucket
      bucket_depth(:,:,current) = bucket_depth(:,:,current) + robert_bucket * bucket_depth(:,:,future) * raw_bucket
   endif

   bucket_depth(:,:,future) = bucket_depth(:,:,future) + robert_bucket * (filt(:,:) + bucket_depth(:,:, future)) &
                           * (raw_bucket - 1.0)  

   where (bucket_depth <= 0.) bucket_depth = 0.

   ! truncate surface reservoir over land points
       where(land .and. (bucket_depth(:,:,future) > max_bucket_depth_land))
            bucket_depth(:,:,future) = max_bucket_depth_land
       end where

   if(id_bucket_depth > 0) used = send_data(id_bucket_depth, bucket_depth(:,:,future), Time)
   if(id_bucket_depth_conv > 0) used = send_data(id_bucket_depth_conv, depth_change_conv(:,:), Time)
   if(id_bucket_depth_cond > 0) used = send_data(id_bucket_depth_cond, depth_change_cond(:,:), Time)
   if(id_bucket_depth_lh > 0) used = send_data(id_bucket_depth_lh, depth_change_lh(:,:), Time)

endif
! end Add bucket section




end subroutine idealized_moist_phys
!=================================================================================================================================
subroutine idealized_moist_phys_end

deallocate (dt_bucket, filt)
if(two_stream_gray)      call two_stream_gray_rad_end
if(lwet_convection)      call qe_moist_convection_end
if(do_ras)               call ras_end

if(turb) then
   call vert_diff_end
   call vert_turb_driver_end
endif
call lscale_cond_end
if(mixed_layer_bc)  call mixed_layer_end(t_surf, bucket_depth, bucket)
if(do_damping) call damping_driver_end

#ifdef SOC_NO_COMPILE
 !No need to end socrates
#else
if(do_socrates_radiation) call run_socrates_end
#endif

end subroutine idealized_moist_phys_end
!=================================================================================================================================

subroutine rh_calc(pfull,T,qv,RH) !s subroutine copied from 2006 FMS MoistModel file moist_processes.f90 (v14 2012/06/22 14:50:00).

        IMPLICIT NONE


        REAL, INTENT (IN),    DIMENSION(:,:,:) :: pfull,T,qv
        REAL, INTENT (OUT),   DIMENSION(:,:,:) :: RH

        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: esat

!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity.
!       This is calculated according to the formula:
!
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!
!       Where epsilon = RDGAS/RVGAS = d622
!
!       and where 1- epsilon = d378
!
!       Note that RH does not have its proper value
!       until all of the following code has been executed.  That
!       is, RH is used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable esat
        CALL LOOKUP_ES(T,esat)						!same as escomp

        !calculate denominator in qsat formula
        if(do_simple) then
          RH(:,:,:) = pfull(:,:,:)
        else
          RH(:,:,:) = pfull(:,:,:)-d378*esat(:,:,:)
        endif

        !limit denominator to esat, and thus qs to epsilon
        !this is done to avoid blow up in the upper stratosphere
        !where pfull ~ esat
        RH(:,:,:) = MAX(RH(:,:,:),esat(:,:,:))

        !calculate RH
        RH(:,:,:)=qv(:,:,:)/(d622*esat(:,:,:)/RH(:,:,:))

        !IF MASK is present set RH to zero
!        IF (present(MASK)) RH(:,:,:)=MASK(:,:,:)*RH(:,:,:)

END SUBROUTINE rh_calc

!=================================================================================================================================

end module idealized_moist_phys_mod
