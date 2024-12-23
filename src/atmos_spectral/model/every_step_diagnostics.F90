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

module every_step_diagnostics_mod

use              fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number

use     time_manager_mod, only: time_type

use    field_manager_mod, only: MODEL_ATMOS

use   tracer_manager_mod, only: get_number_tracers, get_tracer_names

use     diag_manager_mod, only: diag_axis_init, register_diag_field, register_static_field, send_data

use press_and_geopot_mod, only: pressure_variables

use       transforms_mod, only: grid_domain, get_deg_lon, get_deg_lat, get_grid_domain

!===============================================================================================
implicit none
private
!===============================================================================================

public :: every_step_diagnostics_init, every_step_diagnostics, every_step_diagnostics_end

!===============================================================================================

character(len=128), parameter :: version = '$Id: every_step_diagnostics.F90,v 13.0 2006/03/28 21:17:44 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

integer :: id_ps, id_u, id_v, id_t, num_levels, num_tracers
integer, allocatable, dimension(:) :: id_tr

logical :: module_is_initialized = .false.

integer :: two_dt_id_ps, two_dt_id_u, two_dt_id_v, two_dt_id_t
integer, allocatable, dimension(:) :: two_dt_id_tr
integer :: iwt, num_time_steps

real, allocatable, dimension(:,:)     :: two_dt_ps
real, allocatable, dimension(:,:,:)   :: two_dt_u, two_dt_v, two_dt_t
real, allocatable, dimension(:,:,:,:) :: two_dt_tr

integer :: is, ie, js, je

type(time_type) :: Time_save ! every_step_diagnostics_end needs time because it is a required argument of send_data,
                             ! even though the fields are static. When Time is not available to every_step_diagnostics_end
                             ! it uses Time_save instead.
!===============================================================================================

contains

!===============================================================================================

subroutine every_step_diagnostics_init(Time, lon_max, lat_max, num_levels_in, reference_sea_level_press)

type(time_type), intent(in) :: Time
integer, intent(in) :: lon_max, lat_max, num_levels_in
real, intent(in) :: reference_sea_level_press

integer, dimension(2) :: axes_2d
integer, dimension(3) :: axes_3d
integer :: id_lon, id_lat, id_pfull, ntr
real, dimension(2) :: vrange,trange
character(len=128) :: tname, longname, units
character(len=16) :: mod_name = 'dynamics_every'
real, dimension(lon_max) :: lon
real, dimension(lat_max) :: lat
real, dimension(num_levels_in)   :: p_full, ln_p_full
real, dimension(num_levels_in+1) :: p_half, ln_p_half

call write_version_number(version, tagname)

call get_deg_lon(lon)
call get_deg_lat(lat)
id_lon = diag_axis_init('lon_every', lon, 'degrees_E', 'x', 'longitude', set_name=mod_name, Domain2=grid_domain)
id_lat = diag_axis_init('lat_every', lat, 'degrees_N', 'y', 'latitude',  set_name=mod_name, Domain2=grid_domain)
call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, reference_sea_level_press)
p_full = .01*p_full
id_pfull = diag_axis_init('pfull_every', p_full, 'hPa', 'z', 'approx full pressure level', direction=-1, set_name=mod_name)

vrange = (/ -40000., 40000. /)
trange = (/  100., 400. /)
axes_2d = (/ id_lon, id_lat /)
axes_3d = (/ id_lon, id_lat, id_pfull /)

id_ps = register_diag_field(mod_name, 'ps_every', axes_2d, Time, 'surface pressure', 'pascals')
id_u  = register_diag_field(mod_name, 'u_every',  axes_3d, Time, 'zonal wind component', 'm/sec', range=vrange)
id_v  = register_diag_field(mod_name, 'v_every',  axes_3d, Time, 'meridional wind component', 'm/sec', range=vrange)
id_t  = register_diag_field(mod_name, 't_every',  axes_3d, Time, 'temperature', 'deg_k', range=trange)

call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
allocate(id_tr(num_tracers))
do ntr=1,num_tracers
  call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
  id_tr(ntr) = register_diag_field(mod_name, trim(tname)//'_every', axes_3d, Time, longname, units)
enddo

call get_grid_domain(is, ie, js, je)
num_levels = num_levels_in
allocate(two_dt_ps(is:ie, js:je))
allocate(two_dt_u (is:ie, js:je, num_levels))
allocate(two_dt_v (is:ie, js:je, num_levels))
allocate(two_dt_t (is:ie, js:je, num_levels))
allocate(two_dt_tr(is:ie, js:je, num_levels, num_tracers))
allocate(two_dt_id_tr(num_tracers))

two_dt_ps = 0.
two_dt_u  = 0.
two_dt_v  = 0.
two_dt_t  = 0.
two_dt_tr = 0.

two_dt_id_ps = register_static_field(mod_name, '2dt_ps', axes_2d, 'Amplitude of 2*dt wave in surface pressure', 'pascals')
two_dt_id_u  = register_static_field(mod_name, '2dt_u',  axes_3d, 'Amplitude of 2*dt wave in zonal wind', 'm/sec')
two_dt_id_v  = register_static_field(mod_name, '2dt_v',  axes_3d, 'Amplitude of 2*dt wave in meridional wind', 'm/sec')
two_dt_id_t  = register_static_field(mod_name, '2dt_t',  axes_3d, 'Amplitude of 2*dt wave in temperature', 'deg_k')
do ntr=1,num_tracers
  call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
  two_dt_id_tr(ntr) = &
       register_static_field(mod_name, '2dt_'//trim(tname), axes_3d, 'Amplitude of 2*dt wave in '//longname, units)
enddo
iwt = 1
num_time_steps = 0

module_is_initialized = .true.

return
end subroutine every_step_diagnostics_init
!===================================================================================
subroutine every_step_diagnostics(Time, p_surf, u_grid, v_grid, t_grid, tr_grid, time_level)

type(time_type), intent(in) :: Time
real, intent(in), dimension(is:, js:)          :: p_surf
real, intent(in), dimension(is:, js:, :)       :: u_grid, v_grid, t_grid
real, intent(in), dimension(is:, js:, :, :, :) :: tr_grid
integer, intent(in) :: time_level

logical :: used
integer :: i, j, k, ntr
character(len=8) :: err_msg_1, err_msg_2

if(.not.module_is_initialized) then
  call error_mesg('every_step_diagnostics','module is not initialized', FATAL)
endif

Time_save = Time

if(id_ps > 0) used = send_data(id_ps, p_surf, Time)
if(id_u  > 0) used = send_data(id_u,  u_grid, Time)
if(id_v  > 0) used = send_data(id_v,  v_grid, Time)
if(id_t  > 0) used = send_data(id_t,  t_grid, Time)

if(size(tr_grid,5) /= num_tracers) then
  write(err_msg_1,'(i8)') size(tr_grid,5)
  write(err_msg_2,'(i8)') num_tracers
  call error_mesg('every_step_diagnostics','size(tracers)='//err_msg_1//' Should be='//err_msg_2, FATAL)
endif
do ntr=1,num_tracers
  if(id_tr(ntr) > 0) used = send_data(id_tr(ntr), tr_grid(:,:,:,time_level,ntr), Time)
enddo

if(two_dt_id_ps > 0) then
  do j=js,je
    do i=is,ie
      two_dt_ps(i,j) = two_dt_ps(i,j) + iwt*p_surf(i,j)
    enddo
  enddo
endif
if(two_dt_id_u  > 0) then
  do k=1,num_levels
    do j=js,je
      do i=is,ie
        two_dt_u(i,j,k) = two_dt_u(i,j,k) + iwt*u_grid(i,j,k)
      enddo
    enddo
  enddo
endif
if(two_dt_id_v  > 0) then
  do k=1,num_levels
    do j=js,je
      do i=is,ie
        two_dt_v(i,j,k) = two_dt_v(i,j,k) + iwt*v_grid(i,j,k)
      enddo
    enddo
  enddo
endif
if(two_dt_id_t  > 0) then
  do k=1,num_levels
    do j=js,je
      do i=is,ie
        two_dt_t(i,j,k) = two_dt_t(i,j,k) + iwt*t_grid(i,j,k)
      enddo
    enddo
  enddo
endif
do ntr=1,num_tracers
  if(two_dt_id_tr(ntr) > 0) then
    do k=1,num_levels
      do j=js,je
        do i=is,ie
          two_dt_tr(i,j,k,ntr) = two_dt_tr(i,j,k,ntr) + iwt*tr_grid(i,j,k,time_level,ntr)
        enddo
      enddo
    enddo
  endif
enddo

iwt = -1*iwt
num_time_steps = num_time_steps + 1

return
end subroutine every_step_diagnostics
!===================================================================================
subroutine every_step_diagnostics_end(Time_in)
type(time_type), intent(in), optional :: Time_in
logical :: used
integer :: ntr
type(time_type) :: Time

if(.not.module_is_initialized) return

if(present(Time_in)) then
  Time = Time_in
else
  Time = Time_save
endif

two_dt_ps = two_dt_ps / num_time_steps
two_dt_u  = two_dt_u  / num_time_steps
two_dt_v  = two_dt_v  / num_time_steps
two_dt_t  = two_dt_t  / num_time_steps
two_dt_tr = two_dt_tr / num_time_steps

if(two_dt_id_ps > 0) used = send_data(two_dt_id_ps, two_dt_ps, Time)
if(two_dt_id_u  > 0) used = send_data(two_dt_id_u,  two_dt_u,  Time)
if(two_dt_id_v  > 0) used = send_data(two_dt_id_v,  two_dt_v,  Time)
if(two_dt_id_t  > 0) used = send_data(two_dt_id_t,  two_dt_t,  Time)
do ntr=1,num_tracers
  if(two_dt_id_tr(ntr) > 0) used = send_data(two_dt_id_tr(ntr), two_dt_tr(:,:,:,ntr), Time)
enddo

deallocate(id_tr, two_dt_id_tr)
deallocate(two_dt_ps, two_dt_u, two_dt_v, two_dt_t, two_dt_tr)
module_is_initialized = .false.

return
end subroutine every_step_diagnostics_end
!===================================================================================

end module every_step_diagnostics_mod
