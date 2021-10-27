module shallow_dynamics_mod

!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! This program is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------

use               fms_mod, only: open_namelist_file,   &
                                 open_restart_file,    &
                                 file_exist,           &
                                 check_nml_error,      &
                                 error_mesg,           &
                                 FATAL,                &
                                 write_version_number, &
                                 mpp_pe,               &
                                 mpp_root_pe,          &
                                 read_data,            &
                                 write_data,           &
                                 set_domain,           &
                                 close_file,           &
                                 stdlog

use    time_manager_mod,  only : time_type,      &
                                 get_time,       &
                                 operator(==),   &
                                 operator(-)

use       constants_mod,  only : radius, omega, DEG_TO_RAD

use         transforms_mod, only: transforms_init,         transforms_end,          &
                                  get_grid_boundaries,     horizontal_advection,    &
                                  trans_spherical_to_grid, trans_grid_to_spherical, &  
                                  compute_laplacian,       get_eigen_laplacian,     &
                                  get_sin_lat,             get_cos_lat,             &
                                  get_deg_lon,             get_deg_lat,             &
                                  get_grid_domain,         get_spec_domain,         &
                                  spectral_domain,         grid_domain,             &
                                  vor_div_from_uv_grid,    uv_grid_from_vor_div,    &
                                  area_weighted_global_mean

use   spectral_damping_mod, only: spectral_damping_init,   compute_spectral_damping

use           leapfrog_mod, only: leapfrog

use       fv_advection_mod, only : fv_advection_init, a_grid_horiz_advection

use           stirring_mod, only : stirring, stirring_end

use interpolator_mod, only: interpolate_type,interpolator_init,CONSTANT,interpolator

!======================================================================================
implicit none
private
!======================================================================================

public :: shallow_dynamics_init, shallow_dynamics, shallow_dynamics_end, &
          dynamics_type, grid_type, spectral_type, tendency_type


! version information 
!===================================================================
character(len=128) :: version = '$Id: shallow_dynamics.F90,v 10.0 2003/10/24 22:01:02 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'
!===================================================================

type grid_type
   real, pointer, dimension(:,:,:) :: u=>NULL(), v=>NULL(), vor=>NULL(), div=>NULL(), h=>NULL(), trs=>NULL(), tr=>NULL()
   real, pointer, dimension(:,:)   :: stream=>NULL(), pv=>NULL(), deep_geopot=>NULL()
end type
type spectral_type
   complex, pointer, dimension(:,:,:) :: vor=>NULL(), div=>NULL(), h=>NULL(), trs=>NULL()
end type
type tendency_type
   real, pointer, dimension(:,:) :: u=>NULL(), v=>NULL(), h=>NULL(), trs=>NULL(), tr=>NULL()
end type
type dynamics_type
   type(grid_type)     :: grid
   type(spectral_type) :: spec
   type(tendency_type) :: tend
   integer             :: num_lon, num_lat
   logical             :: grid_tracer, spec_tracer
end type



integer, parameter :: num_time_levels = 2

integer :: is, ie, js, je, ms, me, ns, ne

logical :: module_is_initialized = .false.

real,  allocatable, dimension(:) :: sin_lat, cos_lat, rad_lat, deg_lat, deg_lon, &
                                    coriolis

real,  allocatable, dimension(:,:) :: eigen

integer :: pe, npes


! namelist parameters with default values

integer :: num_lon            = 256 
integer :: num_lat            = 128
integer :: num_fourier        = 85
integer :: num_spherical      = 86
integer :: fourier_inc        = 1
integer :: cutoff_wn          = 30
! (these define a standard T85 model)

logical :: check_fourier_imag = .false.
logical :: south_to_north     = .true.
logical :: triang_trunc   = .true.

real    :: robert_coeff        = 0.04
real    :: robert_coeff_tracer = 0.04
real    :: longitude_origin    = 0.0
real    :: raw_filter_coeff    = 1.0

character(len=64) :: damping_option = 'resolution_dependent'
integer :: damping_order       = 4
real    :: damping_coeff       = 1.e-04
real    :: h_0                 = 3.e04

real    :: u_deep_mag          = 0.
real    :: n_merid_deep_flow   = 3.
real    :: u_upper_mag_init    = 0.

logical :: spec_tracer      = .true.
logical :: grid_tracer      = .true.

!Options for injecting an initial vortex pair
real    :: lon_centre_init_cyc = 0.
real    :: lat_centre_init_cyc = 60.
real    :: lon_centre_init_acyc = 180.
real    :: lat_centre_init_acyc = 60.
real    :: init_vortex_radius_deg = 5.
real    :: init_vortex_vor_f = 0.5
real    :: init_vortex_h_h_0 = 0.1
logical :: add_initial_vortex_pair = .false.
logical :: add_initial_vortex_as_height = .true.

logical :: initial_condition_from_input_file=.false.
character(len=64) :: init_cond_file = 'init_cond_h_vor_div'
character(len=64) :: input_file_div_name = 'div'
character(len=64) :: input_file_height_name = 'height'
character(len=64) :: input_file_vor_name = 'vor'

real, dimension(2) :: valid_range_v = (/-1.e3,1.e3/)

type(interpolate_type),save :: init_cond_interp

namelist /shallow_dynamics_nml/ check_fourier_imag,          &
                          south_to_north, triang_trunc,      &
                          num_lon, num_lat, num_fourier,     &
                          num_spherical, fourier_inc,        &
                          longitude_origin, damping_option,  &
                          damping_order,     damping_coeff,  &
                          robert_coeff, robert_coeff_tracer, &
                          h_0, spec_tracer, grid_tracer,     &
                          valid_range_v, cutoff_wn,          &
                          raw_filter_coeff,                  &
                          u_deep_mag, n_merid_deep_flow,     &
                          u_upper_mag_init,                  &
                          lon_centre_init_cyc,               &
                          lat_centre_init_cyc,               &
                          lon_centre_init_acyc,              &
                          lat_centre_init_acyc,              &                          
                          init_vortex_radius_deg,            &
                          init_vortex_vor_f,                 &
                          init_vortex_h_h_0,                 &
                          add_initial_vortex_pair,           &
                          add_initial_vortex_as_height,      &
                          initial_condition_from_input_file, &
                          init_cond_file,                    &
                          input_file_div_name,               &
                          input_file_height_name,            &
                          input_file_vor_name               


contains

!=======================================================================================

subroutine shallow_dynamics_init (Dyn,  Time, Time_init, dt_real)

type(dynamics_type), intent(inout)  :: Dyn
type(time_type)    , intent(in)     :: Time, Time_init
real               , intent(in)     :: dt_real

integer :: i, j

real,    allocatable, dimension(:)   :: glon_bnd, glat_bnd
real,    allocatable, dimension(:,:) :: rad_lonb_2d, rad_latb_2d
real :: xx, yy, dd, deep_geopot_global_mean, radius_loc_cyc, radius_loc_acyc

integer  :: ierr, io, unit, id_lon, id_lat, id_lonb, id_latb
logical  :: root

! < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > 

call write_version_number (version, tagname)

pe   =  mpp_pe()
root = (pe == mpp_root_pe())

if (file_exist('input.nml')) then
  unit = open_namelist_file ()
  ierr=1
  do while (ierr /= 0)
    read  (unit, nml=shallow_dynamics_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'shallow_dynamics_nml')
  enddo
  10 call close_file (unit)
endif

if (root) write (stdlog(), nml=shallow_dynamics_nml)

call transforms_init(radius, num_lat, num_lon, num_fourier, fourier_inc, num_spherical, &
         south_to_north=south_to_north,   &
         triang_trunc=triang_trunc,       &
         longitude_origin=longitude_origin       )

call get_grid_domain(is,ie,js,je)
call get_spec_domain(ms,me,ns,ne)

Dyn%num_lon = num_lon
Dyn%num_lat = num_lat
Dyn%spec_tracer  = spec_tracer
Dyn%grid_tracer  = grid_tracer

allocate (sin_lat   (js:je))
allocate (cos_lat   (js:je))
allocate (deg_lat   (js:je))
allocate (deg_lon   (is:ie))
allocate (coriolis  (js:je))

call get_deg_lon (deg_lon)
call get_deg_lat (deg_lat)
call get_sin_lat (sin_lat)
call get_cos_lat (cos_lat)

allocate (glon_bnd  (num_lon + 1))
allocate (glat_bnd  (num_lat + 1))
call get_grid_boundaries (glon_bnd, glat_bnd, global=.true.)
allocate (rad_lonb_2d(is:ie+1, js:je+1))
allocate (rad_latb_2d(is:ie+1, js:je+1))

do i = is,ie+1
  rad_lonb_2d(i,:) = glon_bnd(i)
enddo

do j = js,je+1
  rad_latb_2d(:,j) = glat_bnd(j)
enddo

coriolis = 2*omega*sin_lat

call spectral_damping_init(damping_coeff, damping_order, damping_option, cutoff_wn, num_fourier, num_spherical, 1, 0., 0., 0.)

allocate(eigen(ms:me,ns:ne))
call get_eigen_laplacian(eigen)

allocate (Dyn%spec%vor (ms:me, ns:ne, num_time_levels))
allocate (Dyn%spec%div (ms:me, ns:ne, num_time_levels))
allocate (Dyn%spec%h   (ms:me, ns:ne, num_time_levels))

allocate (Dyn%grid%u   (is:ie, js:je, num_time_levels))
allocate (Dyn%grid%v   (is:ie, js:je, num_time_levels))
allocate (Dyn%grid%vor (is:ie, js:je, num_time_levels))
allocate (Dyn%grid%div (is:ie, js:je, num_time_levels))
allocate (Dyn%grid%h   (is:ie, js:je, num_time_levels))

allocate (Dyn%tend%u        (is:ie, js:je))
allocate (Dyn%tend%v        (is:ie, js:je))
allocate (Dyn%tend%h        (is:ie, js:je))
allocate (Dyn%grid%stream   (is:ie, js:je))
allocate (Dyn%grid%pv       (is:ie, js:je))
allocate (Dyn%grid%deep_geopot(is:ie, js:je))


call fv_advection_init(num_lon, num_lat, glat_bnd, 360./float(fourier_inc))
if(Dyn%grid_tracer) then
  allocate(Dyn%Grid%tr (is:ie, js:je, num_time_levels))
  allocate(Dyn%Tend%tr (is:ie, js:je))
endif

if(Dyn%spec_tracer) then
  allocate(Dyn%Grid%trs (is:ie, js:je, num_time_levels))
  allocate(Dyn%Tend%trs (is:ie, js:je))
  allocate(Dyn%Spec%trs (ms:me, ns:ne, num_time_levels))
endif

if( initial_condition_from_input_file ) then
  call interpolator_init( init_cond_interp, trim(init_cond_file)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
endif


do i = is, ie 
    Dyn%grid%deep_geopot(i, js:je) = -2.*omega * u_deep_mag * radius * (1./(1.-n_merid_deep_flow**2.))*(-cos(n_merid_deep_flow*DEG_TO_RAD*deg_lat(js:je))*cos(DEG_TO_RAD*deg_lat(js:je)) - n_merid_deep_flow * (sin(n_merid_deep_flow*DEG_TO_RAD*deg_lat(js:je))*sin(DEG_TO_RAD*deg_lat(js:je))-sin(n_merid_deep_flow*(2.*atan(1.)))))
enddo

deep_geopot_global_mean = area_weighted_global_mean(Dyn%grid%deep_geopot(:,:))
Dyn%grid%deep_geopot(:,:) = Dyn%grid%deep_geopot(:,:)-deep_geopot_global_mean

if(Time == Time_init) then

  if (initial_condition_from_input_file) then 

    call interpolator( init_cond_interp, Dyn%Grid%div(:,:,1), input_file_div_name )
    call interpolator( init_cond_interp, Dyn%Grid%h(:,:,1), input_file_height_name )    
    call interpolator( init_cond_interp, Dyn%Grid%vor(:,:,1), input_file_vor_name )

    Dyn%Grid%h(:,:,1) = Dyn%Grid%h(:,:,1)+h_0 !want to make sure that we keep h_0 consistent, so make sure mean of h input is zero, and add h_0 on afterwards...

  else
      Dyn%Grid%div(:,:,1) = 0.0
      Dyn%Grid%h  (:,:,1) = h_0 - Dyn%grid%deep_geopot(:,:)

    do i = is, ie   
      Dyn%Grid%vor(i,js:je,1) = -((u_upper_mag_init * n_merid_deep_flow)/radius) * sin(DEG_to_RAD * deg_lat(js:je))

      if (add_initial_vortex_pair) then

          do j=js, je

              radius_loc_cyc = ((min((deg_lon(i)-lon_centre_init_cyc)**2., (deg_lon(i)-lon_centre_init_cyc-360.)**2.)+(deg_lat(j)-lat_centre_init_cyc)**2.)**0.5)/init_vortex_radius_deg
              radius_loc_acyc = ((min((deg_lon(i)-lon_centre_init_acyc)**2., (deg_lon(i)-lon_centre_init_acyc-360.)**2.)+(deg_lat(j)-lat_centre_init_acyc)**2.)**0.5)/init_vortex_radius_deg            


              if(radius_loc_cyc.le.1.0 .and. radius_loc_acyc.le.1.0) then
                  call error_mesg('shallow_dynamics','Cannot initialise cyclone and anticyclone in same grid box ', FATAL)
              endif

              if(add_initial_vortex_as_height) then
                  if (radius_loc_cyc.le.2.0) then
                      Dyn%Grid%h(i,j,1) = Dyn%Grid%h(i,j,1) + init_vortex_h_h_0 * -h_0 * exp(-radius_loc_cyc**2.)
                  elseif (radius_loc_acyc.le.2.0) then 
                      Dyn%Grid%h(i,j,1) = Dyn%Grid%h(i,j,1) + init_vortex_h_h_0 * h_0 * exp(-radius_loc_acyc**2.)
                  endif
              else
                  if (radius_loc_cyc.le.1.0) then
                      Dyn%Grid%vor(i,j,1) = init_vortex_vor_f * 2.*omega
                  elseif (radius_loc_acyc.le.1.0) then 
                      Dyn%Grid%vor(i,j,1) = init_vortex_vor_f * -2.*omega
                  endif

              endif

          enddo
          

      endif !add_initial_vortex_pair
    enddo

  endif ! initial_condition_from_input_file
    
  call trans_grid_to_spherical(Dyn%Grid%vor(:,:,1), Dyn%Spec%vor(:,:,1))
  call trans_grid_to_spherical(Dyn%Grid%div(:,:,1), Dyn%Spec%div(:,:,1))
  call trans_grid_to_spherical(Dyn%Grid%h  (:,:,1), Dyn%Spec%h  (:,:,1))

  call uv_grid_from_vor_div   (Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1),  &
                               Dyn%Grid%u  (:,:,1), Dyn%Grid%v  (:,:,1))
  
  if(Dyn%grid_tracer) then
    Dyn%Grid%tr = 0.0
    do j = js, je
      if(deg_lat(j) > 10.0 .and. deg_lat(j) < 20.0) Dyn%Grid%tr(:,j,1) =  1.0
      if(deg_lat(j) > 70.0 )                        Dyn%Grid%tr(:,j,1) = -1.0
    end do
  endif
  
  if(Dyn%spec_tracer) then
    Dyn%Grid%trs = 0.0
    do j = js, je
      if(deg_lat(j) > 10.0 .and. deg_lat(j) < 20.0) Dyn%Grid%trs(:,j,1) =  1.0
      if(deg_lat(j) > 70.0 )                        Dyn%Grid%trs(:,j,1) = -1.0
    end do
    call trans_grid_to_spherical(Dyn%Grid%trs(:,:,1), Dyn%Spec%trs(:,:,1))
  endif
  
else

  call read_restart(Dyn)
  
endif

module_is_initialized = .true.

return
end subroutine shallow_dynamics_init

!========================================================================================

subroutine shallow_dynamics(Time, Time_init, Dyn, previous, current, future, delta_t)

type(time_type)    , intent(in)    :: Time, Time_init
type(dynamics_type), intent(inout) :: Dyn
integer, intent(in   )  :: previous, current, future
real,    intent(in   )  :: delta_t

! < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < >

complex, dimension(ms:me, ns:ne)  :: dt_vors, dt_divs, dt_hs, stream, bs, work
real,    dimension(is:ie, js:je)  :: vorg, bg, h_future, h_dt, dt_vorg
integer :: j

! < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > 

if(.not.module_is_initialized) then
  call error_mesg('shallow_dynamics','dynamics has not been initialized ', FATAL)
endif


do j = js,je
  vorg(:,j) = Dyn%Grid%vor(:,j,current) + coriolis(j) 
end do
Dyn%Tend%u = Dyn%Tend%u + vorg*Dyn%Grid%v(:,:,current)
Dyn%Tend%v = Dyn%Tend%v - vorg*Dyn%Grid%u(:,:,current)

call vor_div_from_uv_grid (Dyn%Tend%u, Dyn%Tend%v, dt_vors, dt_divs)

call horizontal_advection (Dyn%Spec%h(:,:,current),   &
         Dyn%Grid%u(:,:,current), Dyn%Grid%v(:,:,current), Dyn%Tend%h)

Dyn%Tend%h = Dyn%Tend%h - Dyn%Grid%h(:,:,current)*Dyn%Grid%div(:,:,current)

call trans_grid_to_spherical (Dyn%Tend%h, dt_hs)

bg = (Dyn%Grid%h(:,:,current) + Dyn%grid%deep_geopot(:,:) + &
   0.5*(Dyn%Grid%u(:,:,current)**2 + Dyn%Grid%v(:,:,current)**2))

call trans_grid_to_spherical(bg, bs)
dt_divs = dt_divs - compute_laplacian(bs)

call implicit_correction (dt_divs, dt_hs, Dyn%Spec%div, Dyn%Spec%h, &
                            delta_t, previous, current)

call compute_spectral_damping(Dyn%Spec%vor(:,:,previous), dt_vors, delta_t)
call compute_spectral_damping(Dyn%Spec%div(:,:,previous), dt_divs, delta_t)
call compute_spectral_damping(Dyn%Spec%h  (:,:,previous), dt_hs  , delta_t)

call stirring(Time, dt_vors)

call leapfrog(Dyn%Spec%vor , dt_vors  , previous, current, future, delta_t, robert_coeff, raw_filter_coeff)
call leapfrog(Dyn%Spec%div , dt_divs  , previous, current, future, delta_t, robert_coeff, raw_filter_coeff)
call leapfrog(Dyn%Spec%h   , dt_hs    , previous, current, future, delta_t, robert_coeff, raw_filter_coeff)

call trans_spherical_to_grid(Dyn%Spec%vor(:,:,future), Dyn%Grid%vor(:,:,future))
call trans_spherical_to_grid(Dyn%Spec%div(:,:,future), Dyn%Grid%div(:,:,future))
  

call uv_grid_from_vor_div  (Dyn%Spec%vor (:,:,future),  Dyn%Spec%div(:,:,future), &
                            Dyn%Grid%u   (:,:,future),  Dyn%Grid%v  (:,:,future))

call trans_spherical_to_grid (Dyn%Spec%h (:,:,future),  Dyn%Grid%h(:,:,future))

if(minval(Dyn%Grid%v) < valid_range_v(1) .or. maxval(Dyn%Grid%v) > valid_range_v(2)) then
  call error_mesg('shallow_dynamics','meridional wind out of valid range', FATAL)
endif

if(Dyn%spec_tracer) call update_spec_tracer(Dyn%Spec%trs, Dyn%Grid%trs, Dyn%Tend%trs, &
                         Dyn%Grid%u, Dyn%Grid%v, previous, current, future, delta_t)

if(Dyn%grid_tracer) call update_grid_tracer(Dyn%Grid%tr, Dyn%Tend%tr, &
                         Dyn%Grid%u, Dyn%Grid%v, previous, current, future, delta_t)


!  for diagnostics

stream = compute_laplacian(Dyn%Spec%vor(:,:,current), -1) ! for diagnostic purposes
call trans_spherical_to_grid(stream, Dyn%grid%stream)

Dyn%Grid%pv = vorg/(Dyn%Grid%h(:,:,current))

return
end subroutine shallow_dynamics
!================================================================================

subroutine implicit_correction(dt_divs, dt_hs, divs, hs, delta_t, previous, current)

complex, intent(inout), dimension(ms:,ns:)   :: dt_divs, dt_hs
complex, intent(in),    dimension(ms:,ns:,:) :: divs, hs
real   , intent(in)                          :: delta_t
integer, intent(in)                          :: previous, current

real :: xi, mu, mu2

xi    = 0.5 ! centered implicit   (for backwards implicit, set xi = 1.0)

mu  = xi*delta_t
mu2 = mu*mu

dt_hs   = dt_hs + h_0*(divs(:,:,current) - divs(:,:,previous))
dt_divs = dt_divs - eigen*(hs(:,:,current) - hs(:,:,previous))

dt_divs = (dt_divs + mu*eigen*dt_hs)/(1.0 + mu2*eigen*h_0)
dt_hs   =  dt_hs - mu*h_0*dt_divs

return
end subroutine implicit_correction 

!===================================================================================

subroutine update_spec_tracer(tr_spec, tr_grid, dt_tr, ug, vg, &
                              previous, current, future, delta_t)

complex, intent(inout), dimension(ms:me, ns:ne, num_time_levels) :: tr_spec
real   , intent(inout), dimension(is:ie, js:je, num_time_levels) :: tr_grid
real   , intent(inout), dimension(is:ie, js:je                 ) :: dt_tr
real   , intent(in   ), dimension(is:ie, js:je, num_time_levels) :: ug, vg  
real   , intent(in   )  :: delta_t
integer, intent(in   )  :: previous, current, future

complex, dimension(ms:me, ns:ne) :: dt_trs

call horizontal_advection     (tr_spec(:,:,current), ug(:,:,current), vg(:,:,current), dt_tr)
call trans_grid_to_spherical  (dt_tr, dt_trs)
call compute_spectral_damping (tr_spec(:,:,previous), dt_trs, delta_t)
call leapfrog                 (tr_spec, dt_trs, previous, current, future, delta_t, robert_coeff, raw_filter_coeff)
call trans_spherical_to_grid  (tr_spec(:,:,future), tr_grid(:,:,future))

return
end subroutine update_spec_tracer
!==========================================================================

subroutine update_grid_tracer(tr_grid, dt_tr_grid, ug, vg, &
                              previous, current, future, delta_t)

real   , intent(inout), dimension(is:ie, js:je, num_time_levels) :: tr_grid
real   , intent(inout), dimension(is:ie, js:je                 ) :: dt_tr_grid
real   , intent(in   ), dimension(is:ie, js:je, num_time_levels) :: ug, vg

real   , intent(in   )  :: delta_t
integer, intent(in   )  :: previous, current, future

real, dimension(is:ie,js:je) :: tr_current, tr_future

tr_future = tr_grid(:,:,previous) + delta_t*dt_tr_grid
dt_tr_grid = 0.0
call a_grid_horiz_advection (ug(:,:,current), vg(:,:,current), tr_future, delta_t, dt_tr_grid)
tr_future = tr_future + delta_t*dt_tr_grid
tr_current = tr_grid(:,:,current) + &
    robert_coeff_tracer*(tr_grid(:,:,previous) + tr_future - 2.0*tr_grid(:,:,current))
tr_grid(:,:,current) = tr_current
tr_grid(:,:,future)  = tr_future

return
end subroutine update_grid_tracer

!==========================================================================

subroutine read_restart(Dyn)

type(dynamics_type), intent(inout)  :: Dyn

integer :: unit, m, n, nt
real, dimension(ms:me, ns:ne) :: real_part, imag_part

if(file_exist('INPUT/shallow_dynamics.res.nc')) then
  do nt = 1, 2
    call read_data('INPUT/shallow_dynamics.res.nc', 'vors_real', real_part, spectral_domain, timelevel=nt)
    call read_data('INPUT/shallow_dynamics.res.nc', 'vors_imag', imag_part, spectral_domain, timelevel=nt)
    do n=ns,ne
      do m=ms,me
        Dyn%Spec%vor(m,n,nt) = cmplx(real_part(m,n),imag_part(m,n))
      end do
    end do
    call read_data('INPUT/shallow_dynamics.res.nc', 'divs_real', real_part, spectral_domain, timelevel=nt)
    call read_data('INPUT/shallow_dynamics.res.nc', 'divs_imag', imag_part, spectral_domain, timelevel=nt)
    do n=ns,ne
      do m=ms,me
        Dyn%Spec%div(m,n,nt) = cmplx(real_part(m,n),imag_part(m,n))
      end do
    end do
    call read_data('INPUT/shallow_dynamics.res.nc', 'hs_real', real_part, spectral_domain, timelevel=nt)
    call read_data('INPUT/shallow_dynamics.res.nc', 'hs_imag', imag_part, spectral_domain, timelevel=nt)
    do n=ns,ne
      do m=ms,me
        Dyn%Spec%h(m,n,nt) = cmplx(real_part(m,n),imag_part(m,n))
      end do
    end do
    if(Dyn%spec_tracer) then
      call read_data('INPUT/shallow_dynamics.res.nc', 'trs_real', real_part, spectral_domain, timelevel=nt)
      call read_data('INPUT/shallow_dynamics.res.nc', 'trs_imag', imag_part, spectral_domain, timelevel=nt)
      do n=ns,ne
        do m=ms,me
          Dyn%Spec%trs(m,n,nt) = cmplx(real_part(m,n),imag_part(m,n))
        end do
      end do
    endif
    call read_data('INPUT/shallow_dynamics.res.nc', 'u',   Dyn%Grid%u  (:,:,nt), grid_domain, timelevel=nt)
    call read_data('INPUT/shallow_dynamics.res.nc', 'v',   Dyn%Grid%v  (:,:,nt), grid_domain, timelevel=nt)
    call read_data('INPUT/shallow_dynamics.res.nc', 'vor', Dyn%Grid%vor(:,:,nt), grid_domain, timelevel=nt)
    call read_data('INPUT/shallow_dynamics.res.nc', 'div', Dyn%Grid%div(:,:,nt), grid_domain, timelevel=nt)
    call read_data('INPUT/shallow_dynamics.res.nc', 'h',   Dyn%Grid%h  (:,:,nt), grid_domain, timelevel=nt)
    if(Dyn%spec_tracer) then
      call read_data('INPUT/shallow_dynamics.res.nc', 'trs', Dyn%Grid%trs(:,:,nt), grid_domain, timelevel=nt)
    endif
    if(Dyn%grid_tracer) then
      call read_data('INPUT/shallow_dynamics.res.nc', 'tr', Dyn%Grid%tr(:,:,nt), grid_domain, timelevel=nt)
    endif
  end do
else if(file_exist('INPUT/shallow_dynamics.res')) then
  unit = open_restart_file(file='INPUT/shallow_dynamics.res',action='read')

  do nt = 1, 2
    call set_domain(spectral_domain)
    call read_data(unit,Dyn%Spec%vor(:,:, nt))
    call read_data(unit,Dyn%Spec%div(:,:, nt))
    call read_data(unit,Dyn%Spec%h  (:,:, nt))
    if(Dyn%spec_tracer) call read_data(unit,Dyn%Spec%trs(:,:, nt))

    call set_domain(grid_domain)
    call read_data(unit,Dyn%Grid%u   (:,:, nt))
    call read_data(unit,Dyn%Grid%v   (:,:, nt))
    call read_data(unit,Dyn%Grid%vor (:,:, nt))
    call read_data(unit,Dyn%Grid%div (:,:, nt))
    call read_data(unit,Dyn%Grid%h   (:,:, nt))
    if(Dyn%spec_tracer) call read_data(unit,Dyn%Grid%trs(:,:, nt))
    if(Dyn%grid_tracer) call read_data(unit,Dyn%Grid%tr (:,:, nt))
    
  end do
  call close_file(unit)
  
else
    call error_mesg('read_restart', 'restart does not exist', FATAL)
endif
  
return
end subroutine read_restart

!====================================================================

subroutine write_restart(Dyn, previous, current)

type(dynamics_type), intent(in)  :: Dyn
integer, intent(in) :: previous, current

integer :: unit, nt, nn

do nt = 1, 2
  if(nt == 1) nn = previous
  if(nt == 2) nn = current
  call write_data('RESTART/shallow_dynamics.res.nc', 'vors_real',  real(Dyn%Spec%vor(:,:,nn)), spectral_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'vors_imag', aimag(Dyn%Spec%vor(:,:,nn)), spectral_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'divs_real',  real(Dyn%Spec%div(:,:,nn)), spectral_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'divs_imag', aimag(Dyn%Spec%div(:,:,nn)), spectral_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'hs_real',    real(Dyn%Spec%h  (:,:,nn)), spectral_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'hs_imag',   aimag(Dyn%Spec%h  (:,:,nn)), spectral_domain)
  if(Dyn%spec_tracer) then
    call write_data('RESTART/shallow_dynamics.res.nc', 'trs_real',  real(Dyn%Spec%trs(:,:,nn)), spectral_domain)
    call write_data('RESTART/shallow_dynamics.res.nc', 'trs_imag', aimag(Dyn%Spec%trs(:,:,nn)), spectral_domain)
  endif
  call write_data('RESTART/shallow_dynamics.res.nc', 'u',   Dyn%Grid%u  (:,:,nn), grid_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'v',   Dyn%Grid%v  (:,:,nn), grid_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'vor', Dyn%Grid%vor(:,:,nn), grid_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'div', Dyn%Grid%div(:,:,nn), grid_domain)
  call write_data('RESTART/shallow_dynamics.res.nc', 'h',   Dyn%Grid%h  (:,:,nn), grid_domain)
  if(Dyn%spec_tracer) then
    call write_data('RESTART/shallow_dynamics.res.nc', 'trs', Dyn%Grid%trs(:,:,nn), grid_domain)
  endif
  if(Dyn%grid_tracer) then
    call write_data('RESTART/shallow_dynamics.res.nc', 'tr', Dyn%Grid%tr(:,:,nn), grid_domain)
  endif
enddo

!unit = open_restart_file(file='RESTART/shallow_dynamics.res', action='write')

!do n = 1, 2
!  if(n == 1) nn = previous
!  if(n == 2) nn = current
!  
!  call set_domain(spectral_domain)
!  call write_data(unit,Dyn%Spec%vor(:,:, nn))
!  call write_data(unit,Dyn%Spec%div(:,:, nn))
!  call write_data(unit,Dyn%Spec%h  (:,:, nn))
!  if(Dyn%spec_tracer) call write_data(unit,Dyn%Spec%trs(:,:, nn))
!
!  call set_domain(grid_domain)
!  call write_data(unit,Dyn%Grid%u   (:,:, nn))
!  call write_data(unit,Dyn%Grid%v   (:,:, nn))
!  call write_data(unit,Dyn%Grid%vor (:,:, nn))
!  call write_data(unit,Dyn%Grid%div (:,:, nn))
!  call write_data(unit,Dyn%Grid%h   (:,:, nn))
!  if(Dyn%spec_tracer) call write_data(unit,Dyn%Grid%trs(:,:, nn))
!  if(Dyn%grid_tracer) call write_data(unit,Dyn%Grid%tr (:,:, nn))
!  
!end do

!call close_file(unit)

end subroutine write_restart

!====================================================================

subroutine shallow_dynamics_end (Dyn, previous, current)

type(dynamics_type), intent(inout)  :: Dyn
integer, intent(in) :: previous, current

if(.not.module_is_initialized) then
  call error_mesg('shallow_dynamics_end','dynamics has not been initialized ', FATAL)
endif

call write_restart (Dyn, previous, current)

call transforms_end
call stirring_end

module_is_initialized = .false.

return
end subroutine shallow_dynamics_end
!===================================================================================

end module shallow_dynamics_mod
