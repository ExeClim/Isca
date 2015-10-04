module idealized_moist_phys_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

use fms_mod, only: write_version_number, file_exist, close_file, stdlog, error_mesg, FATAL

use           constants_mod, only: grav

use        time_manager_mod, only: time_type, get_time, operator( + )

use    vert_turb_driver_mod, only: vert_turb_driver_init, vert_turb_driver, vert_turb_driver_end

use           vert_diff_mod, only: vert_diff_init, gcm_vert_diff_down, gcm_vert_diff_up, vert_diff_end, surf_diff_type

use two_stream_gray_rad_mod, only: two_stream_gray_rad_init, two_stream_gray_rad_down, two_stream_gray_rad_up, two_stream_gray_rad_end

use         mixed_layer_mod, only: mixed_layer_init, mixed_layer, mixed_layer_end

use         lscale_cond_mod, only: lscale_cond_init, lscale_cond, lscale_cond_end

use qe_moist_convection_mod, only: qe_moist_convection_init, qe_moist_convection, qe_moist_convection_end

use        diag_manager_mod, only: register_diag_field, send_data

use          transforms_mod, only: get_grid_domain

use   spectral_dynamics_mod, only: get_axis_id, get_num_levels, get_surf_geopotential

use        surface_flux_mod, only: surface_flux

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
logical :: lwet_convection = .false.
logical :: two_stream = .true.
logical :: mixed_layer_bc = .false.
real :: roughness_heat = 0.05
real :: roughness_moist = 0.05
real :: roughness_mom = 0.05

namelist / idealized_moist_phys_nml / turb, lwet_convection, roughness_heat, two_stream, mixed_layer_bc, &
                                      roughness_moist, roughness_mom, do_virtual

real, allocatable, dimension(:,:)   ::                                        &
     z_surf,               &   ! surface height
     t_surf,               &   ! surface temperature
     q_surf,               &   ! surface moisture
     u_surf,               &   ! surface U wind
     v_surf,               &   ! surface V wind
     rough_mom,            &   ! momentum roughness length for surface_flux
     rough_heat,           &   ! heat roughness length for surface_flux
     rough_moist,          &   ! moisture roughness length for surface_flux
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
     rough                     ! roughness for vert_turb_driver


real, allocatable, dimension(:,:,:) ::                                        &
     diff_m,               &   ! momentum diffusion coeff.
     diff_t,               &   ! temperature diffusion coeff.
     tdtlw,                &   ! place holder. appears in calling arguments of vert_turb_driver but not used unless do_edt=.true. -- pjp
     diss_heat,            &   ! heat dissipated by vertical diffusion
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
     klzbs,                &   ! stored level of zero buoyancy values
     cape,                 &   ! convectively available potential energy
     cin,                  &   ! convective inhibition (this and the above are before the adjustment)
     invtau_q_relaxation,  &   ! temperature relaxation time scale
     invtau_t_relaxation,  &   ! humidity relaxation time scale
     rain,                 &   !
     snow

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
     id_conv_dt_tg,  &   ! temperature tendency from convection
     id_conv_dt_qg,  &   ! temperature tendency from convection
     id_cond_dt_tg,  &   ! temperature tendency from convection
     id_cond_dt_qg       ! temperature tendency from convection

integer, allocatable, dimension(:,:) :: & convflag ! indicates which qe convection subroutines are used
real,    allocatable, dimension(:,:) :: rad_lat   

type(surf_diff_type) :: Tri_surf ! used by gcm_vert_diff

logical :: used, doing_edt, doing_entrain
integer, dimension(4) :: axes
integer :: is, ie, js, je, num_levels, nsphum, dt_integer
real :: dt_real
type(time_type) :: Time_step

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine idealized_moist_phys_init(Time, Time_step_in, nhum, rad_lat_2d, rad_lonb_2d, rad_latb_2d, t_surf_init)
type(time_type), intent(in) :: Time, Time_step_in
integer, intent(in) :: nhum
real, intent(in), dimension(:,:) :: rad_lat_2d, rad_lonb_2d, rad_latb_2d, t_surf_init

integer :: io, nml_unit, stdlog_unit, seconds, days

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

nsphum = nhum
Time_step = Time_step_in
call get_time(Time_step, seconds, days)
dt_integer   = 86400*days + seconds
dt_real      = float(dt_integer)

call get_grid_domain(is, ie, js, je)
call get_num_levels(num_levels)

allocate(rad_lat     (is:ie, js:je)); rad_lat = rad_lat_2d
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
allocate(land        (is:ie, js:je)); land = .false.
allocate(avail       (is:ie, js:je)); avail = .true.
allocate(fracland    (is:ie, js:je)); fracland = 0.0
allocate(rough       (is:ie, js:je))
allocate(diff_t      (is:ie, js:je, num_levels))
allocate(diff_m      (is:ie, js:je, num_levels))
allocate(diss_heat   (is:ie, js:je, num_levels))
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
allocate(convflag     (is:ie, js:je))
allocate(convect      (is:ie, js:je)); convect = .false.

allocate(t_ref (is:ie, js:je, num_levels)); t_ref = 0.0
allocate(q_ref (is:ie, js:je, num_levels)); q_ref = 0.0

call get_surf_geopotential(z_surf)
z_surf = z_surf/grav

if(mixed_layer_bc) then
  ! need an initial condition for the mixed layer temperature
  ! may be overwritten by restart file
  ! choose an unstable initial condition to allow moisture
  ! to quickly enter the atmosphere avoiding problems with the convection scheme
  t_surf = t_surf_init + 1.0
  call mixed_layer_init(is, ie, js, je, num_levels, t_surf, get_axis_id(), Time) ! t_surf is intent(inout)
endif

if(turb) then
! need to call vert_diff_init even if using gcm_vert_diff (rather than
! gcm_vert_diff_down) because the variable sphum is not initialized
! otherwise in the vert_diff module
   call vert_diff_init (Tri_surf, ie-is+1, je-js+1, num_levels, .true., do_virtual)
end if

call lscale_cond_init()

axes = get_axis_id()

id_cond_dt_qg = register_diag_field(mod_name, 'dt_qg_condensation',        &
     axes(1:3), Time, 'Moisture tendency from condensation','kg/kg/s')
id_cond_dt_tg = register_diag_field(mod_name, 'dt_tg_condensation',        &
     axes(1:3), Time, 'Temperature tendency from condensation','K/s')
id_cond_rain = register_diag_field(mod_name, 'condensation_rain',          &
     axes(1:2), Time, 'Rain from condensation','kg/m/m/s')

if(lwet_convection) then
   call qe_moist_convection_init()
   id_conv_dt_qg = register_diag_field(mod_name, 'dt_qg_convection',          &
        axes(1:3), Time, 'Moisture tendency from convection','kg/kg/s')
   id_conv_dt_tg = register_diag_field(mod_name, 'dt_tg_convection',          &
        axes(1:3), Time, 'Temperature tendency from convection','K/s')
   id_conv_rain = register_diag_field(mod_name, 'convection_rain',            &
        axes(1:2), Time, 'Rain from convection','kg/m/m/s')
endif

if(two_stream) call two_stream_gray_rad_init(is, ie, js, je, num_levels, get_axis_id(), Time)

if(turb) then
   call vert_turb_driver_init (rad_lonb_2d, rad_latb_2d, ie-is+1,je-js+1, &
                 num_levels,get_axis_id(),Time, doing_edt, doing_entrain)

   axes = get_axis_id()
   id_diff_dt_ug = register_diag_field(mod_name, 'dt_ug_diffusion',        &
        axes(1:3), Time, 'zonal wind tendency from diffusion','m/s^2')
   id_diff_dt_vg = register_diag_field(mod_name, 'dt_vg_diffusion',        &
        axes(1:3), Time, 'meridional wind tendency from diffusion','m/s^2')
   id_diff_dt_tg = register_diag_field(mod_name, 'dt_tg_diffusion',        &
        axes(1:3), Time, 'temperature diffusion tendency','T/s')
   id_diff_dt_qg = register_diag_field(mod_name, 'dt_qg_diffusion',        &
        axes(1:3), Time, 'moisture diffusion tendency','T/s')
endif

end subroutine idealized_moist_phys_init
!=================================================================================================================================
subroutine idealized_moist_phys(Time, p_half, p_full, z_half, z_full, ug, vg, tg, grid_tracers, &
                                previous, current, dt_ug, dt_vg, dt_tg, dt_tracers)

type(time_type),            intent(in)    :: Time
real, dimension(:,:,:,:),   intent(in)    :: p_half, p_full, z_half, z_full, ug, vg, tg
real, dimension(:,:,:,:,:), intent(in)    :: grid_tracers
integer,                    intent(in)    :: previous, current
real, dimension(:,:,:),     intent(inout) :: dt_ug, dt_vg, dt_tg
real, dimension(:,:,:,:),   intent(inout) :: dt_tracers

real :: delta_t
real, dimension(size(ug,1), size(ug,2), size(ug,3)) :: tg_tmp, qg_tmp

if(current == previous) then
   delta_t = dt_real
else
   delta_t = 2*dt_real
endif

if (lwet_convection) then
   rain = 0.0; snow = 0.0
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
   rain       = rain/delta_t

   dt_tg = dt_tg + conv_dt_tg
   dt_tracers(:,:,:,nsphum) = dt_tracers(:,:,:,nsphum) + conv_dt_qg

   if(id_conv_dt_qg > 0) used = send_data(id_conv_dt_qg, conv_dt_qg, Time)
   if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
   if(id_conv_rain  > 0) used = send_data(id_conv_rain, rain, Time)

else

   tg_tmp = tg(:,:,:,previous)
   qg_tmp = grid_tracers(:,:,:,previous,nsphum)

endif

rain = 0.0
call lscale_cond (         tg_tmp,                          qg_tmp,        &
           p_full(:,:,:,previous),          p_half(:,:,:,previous),        &
                            coldT,                            rain,        &
                             snow,                      cond_dt_tg,        &
                       cond_dt_qg )
                                                                          
cond_dt_tg = cond_dt_tg/delta_t
cond_dt_qg = cond_dt_qg/delta_t
rain       = rain/delta_t
                                                                             
dt_tg = dt_tg + cond_dt_tg
dt_tracers(:,:,:,nsphum) = dt_tracers(:,:,:,nsphum) + cond_dt_qg
                                                                               
if(id_cond_dt_qg > 0) used = send_data(id_cond_dt_qg, cond_dt_qg, Time)
if(id_cond_dt_tg > 0) used = send_data(id_cond_dt_tg, cond_dt_tg, Time)
if(id_cond_rain  > 0) used = send_data(id_cond_rain, rain, Time)

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

if(two_stream) then
   call two_stream_gray_rad_down(is, js, Time, &
                       rad_lat(:,:),           &
                       p_half(:,:,:,current),  &
                       tg(:,:,:,previous),     &
                       net_surf_sw_down(:,:),  &
                       surf_lw_down(:,:))
end if

if(.not.mixed_layer_bc) then
                                                                               
!!$! infinite heat capacity
!    t_surf = surface_temperature_forced(rad_lat)
!!$! no heat capacity:
!!$   t_surf = tg(:,:,num_levels,previous)
                                                                              
!!$! surface temperature has same potential temp. as lowest layer:
!!$  t_surf = surface_temperature(tg(:,:,:,previous), p_full(:,:,:,current), p_half(:,:,:,current))
end if

call surface_flux(                                                          &
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
                                      delta_t,                              &
                                    land(:,:),                              &
                               .not.land(:,:),                              &
                                   avail(:,:)  )

! Now complete the radiation calculation by computing the upward and net fluxes.

if(two_stream) then
   call two_stream_gray_rad_up(is, js, Time, &
                     rad_lat(:,:),           &
                     p_half(:,:,:,current),  &
                     t_surf(:,:),            &
                     tg(:,:,:,previous),     &
                     dt_tg(:,:,:))
end if

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

   if(.not.mixed_layer_bc) then
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
   call mixed_layer(                                                       &
                              Time,                                        &
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
                            dedq_atm(:,:))


   call gcm_vert_diff_up (1, 1, delta_t, Tri_surf, dt_tg(:,:,:), dt_tracers(:,:,:,nsphum), dt_tracers(:,:,:,:))

   if(id_diff_dt_ug > 0) used = send_data(id_diff_dt_ug, dt_ug - non_diff_dt_ug, Time)
   if(id_diff_dt_vg > 0) used = send_data(id_diff_dt_vg, dt_vg - non_diff_dt_vg, Time)
   if(id_diff_dt_tg > 0) used = send_data(id_diff_dt_tg, dt_tg - non_diff_dt_tg, Time)
   if(id_diff_dt_qg > 0) used = send_data(id_diff_dt_qg, dt_tracers(:,:,:,nsphum) - non_diff_dt_qg, Time)

endif ! if(turb) then

end subroutine idealized_moist_phys
!=================================================================================================================================
subroutine idealized_moist_phys_end

if(two_stream)      call two_stream_gray_rad_end
if(lwet_convection) call qe_moist_convection_end
if(turb) then
   call vert_diff_end
   call vert_turb_driver_end
endif
call lscale_cond_end
if(mixed_layer_bc)  call mixed_layer_end(t_surf)

end subroutine idealized_moist_phys_end
!=================================================================================================================================

end module idealized_moist_phys_mod
