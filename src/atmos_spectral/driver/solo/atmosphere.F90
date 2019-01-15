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
 
module atmosphere_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

use                  fms_mod, only: set_domain, write_version_number, field_size, file_exist, stdlog, &
                                    mpp_pe, mpp_root_pe, error_mesg, FATAL, read_data, write_data, nullify_domain

use            constants_mod, only: grav, pi

use         time_manager_mod, only: time_type, get_time, operator( + )

use     press_and_geopot_mod, only: compute_pressures_and_heights

#ifdef COLUMN_MODEL
  use             spec_mpp_mod, only: grid_domain, get_grid_domain, atmosphere_domain
  use               column_mod, only: column_init, column, column_end, &
                                      get_axis_id, column_diagnostics, get_num_levels, &
                                      get_surf_geopotential, get_initial_fields
  use          column_grid_mod, only: get_deg_lon, get_deg_lat, get_grid_boundaries, &
                                      get_lon_max, get_lat_max  
#else
  use    spectral_dynamics_mod, only: spectral_dynamics_init, spectral_dynamics, spectral_dynamics_end, &
                                      get_num_levels, get_axis_id, spectral_diagnostics, get_initial_fields, &
                                      complete_robert_filter, get_surf_geopotential
  use           transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, &
                                      get_deg_lon, get_deg_lat, get_grid_boundaries, grid_domain,    &
                                      spectral_domain, get_grid_domain, get_lon_max, get_lat_max, atmosphere_domain
#endif

use          tracer_type_mod, only: tracer_type

use           hs_forcing_mod, only: hs_forcing_init, hs_forcing, hs_forcing_end

use        field_manager_mod, only: MODEL_ATMOS

use       tracer_manager_mod, only: get_number_tracers

use idealized_moist_phys_mod, only: idealized_moist_phys_init , idealized_moist_phys , idealized_moist_phys_end

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: atmosphere.F90,v 19.0.2.1 2013/01/24 14:54:27 pjp Exp $'
      
character(len=128) :: tagname= &
'$Name:  $'
character(len=10), parameter :: mod_name='atmosphere'

!=================================================================================================================================

public :: atmosphere_init, atmosphere, atmosphere_end, atmosphere_domain

!=================================================================================================================================
logical :: idealized_moist_model = .false.

namelist/atmosphere_nml/ idealized_moist_model

!=================================================================================================================================

integer, parameter :: num_time_levels = 2
integer :: is, ie, js, je, num_levels, num_tracers, nhum
logical :: dry_model, column_model

real, allocatable, dimension(:,:,:,:) :: p_half, p_full
real, allocatable, dimension(:,:,:,:) :: z_half, z_full

type(tracer_type), allocatable, dimension(:) :: tracer_attributes
real, allocatable, dimension(:,:,:,:,:) :: grid_tracers
real, allocatable, dimension(:,:,:    ) :: psg, wg_full
real, allocatable, dimension(:,:,:,:  ) :: ug, vg, tg

real, allocatable, dimension(:,:    ) :: dt_psg
real, allocatable, dimension(:,:,:  ) :: dt_ug, dt_vg, dt_tg
real, allocatable, dimension(:,:,:,:) :: dt_tracers

real, allocatable, dimension(:)   :: deg_lon, deg_lat
real, allocatable, dimension(:,:) :: rad_lon_2d, rad_lat_2d
real, allocatable, dimension(:,:) :: surf_geopotential
real, allocatable, dimension(:,:) :: rad_lonb_2d, rad_latb_2d
real, allocatable, dimension(:)   :: rad_lonb, rad_latb

integer :: previous, current, future
logical :: module_is_initialized =.false.

integer         :: dt_integer
real            :: dt_real
type(time_type) :: Time_step

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine atmosphere_init(Time_init, Time, Time_step_in)

type (time_type), intent(in)  :: Time_init, Time, Time_step_in

integer :: seconds, days, lon_max, lat_max, ntr, nt, i, j, nml_unit, io, stdlog_unit
integer, dimension(4) :: siz
real, dimension(2) :: time_pointers
character(len=64) :: file, tr_name
character(len=256) :: message

if(module_is_initialized) return

call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=atmosphere_nml, iostat=io)
#else  
   if ( file_exist('input.nml') ) then
      nml_unit = open_namelist_file()
      read (nml_unit, atmosphere_nml, iostat=io)
      call close_file(nml_unit)
   endif
#endif
stdlog_unit = stdlog()
write(stdlog_unit, atmosphere_nml)

Time_step = Time_step_in
call get_time(Time_step, seconds, days)
dt_integer   = 86400*days + seconds
dt_real      = float(dt_integer)

call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
allocate (tracer_attributes(num_tracers))
#ifdef COLUMN_MODEL
  call column_init(Time, Time_step, tracer_attributes, dry_model, nhum)
  column_model = .true.
#else
  call spectral_dynamics_init(Time, Time_step, tracer_attributes, dry_model, nhum)
  column_model = .false.
#endif 
call get_grid_domain(is, ie, js, je)
call get_num_levels(num_levels)

allocate (p_half       (is:ie, js:je, num_levels+1, num_time_levels))
allocate (z_half       (is:ie, js:je, num_levels+1, num_time_levels))
allocate (p_full       (is:ie, js:je, num_levels, num_time_levels))
allocate (z_full       (is:ie, js:je, num_levels, num_time_levels))
allocate (wg_full      (is:ie, js:je, num_levels))
allocate (psg          (is:ie, js:je, num_time_levels))
allocate (ug           (is:ie, js:je, num_levels, num_time_levels))
allocate (vg           (is:ie, js:je, num_levels, num_time_levels))
allocate (tg           (is:ie, js:je, num_levels, num_time_levels))
allocate (grid_tracers (is:ie, js:je, num_levels, num_time_levels, num_tracers ))

allocate (dt_psg     (is:ie, js:je))
allocate (dt_ug      (is:ie, js:je, num_levels))
allocate (dt_vg      (is:ie, js:je, num_levels))
allocate (dt_tg      (is:ie, js:je, num_levels))
allocate (dt_tracers (is:ie, js:je, num_levels, num_tracers))

allocate (deg_lon    (is:ie       ))
allocate (rad_lon_2d (is:ie, js:je))
allocate (deg_lat    (       js:je))
allocate (rad_lat_2d (is:ie, js:je))
allocate (rad_lonb_2d(is:ie+1, js:je+1))
allocate (rad_latb_2d(is:ie+1, js:je+1))
allocate (rad_lonb (is:ie+1))
allocate (rad_latb (js:je+1))

p_half = 0.; z_half = 0.; p_full = 0.; z_full = 0.
wg_full = 0.; psg = 0.; ug = 0.; vg = 0.; tg = 0.; grid_tracers = 0.
dt_psg = 0.; dt_ug  = 0.; dt_vg  = 0.; dt_tg  = 0.; dt_tracers = 0.

allocate (surf_geopotential(is:ie, js:je))
call get_surf_geopotential(surf_geopotential)

!--------------------------------------------------------------------------------------------------------------------------------
file = 'INPUT/atmosphere.res.nc'
if(file_exist(trim(file))) then
  call get_lon_max(lon_max)
  call get_lat_max(lat_max)
  call field_size(trim(file), 'ug', siz)
  if(lon_max /= siz(1) .or. lat_max /= siz(2)) then
    write(message,*) 'Resolution of restart data does not match resolution specified on namelist. Restart data: lon_max=', &
                     siz(1),', lat_max=',siz(2),'  Namelist: lon_max=',lon_max,', lat_max=',lat_max
    call error_mesg('atmosphere_init', message, FATAL)
  endif
  call nullify_domain()
! call read_data(trim(file), 'previous', previous)           ! No interface of read_data exists to read integer scalars
! call read_data(trim(file), 'current',  current)            ! No interface of read_data exists to read integer scalars
  call read_data(trim(file), 'time_pointers', time_pointers) ! Getaround for no interface to read integer scalars
  previous = int(time_pointers(1))                           ! Getaround for no interface to read integer scalars
  current  = int(time_pointers(2))                           ! Getaround for no interface to read integer scalars
  do nt=1,num_time_levels
    call read_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'psg', psg(:,:,  nt), grid_domain, timelevel=nt)
    do ntr = 1,num_tracers
      tr_name = trim(tracer_attributes(ntr)%name)
      call read_data(trim(file), trim(tr_name), grid_tracers(:,:,:,nt,ntr), grid_domain, timelevel=nt)
    enddo ! end loop over tracers
  enddo ! end loop over time levels
  call read_data(trim(file), 'wg_full', wg_full, grid_domain)
else
  previous = 1; current = 1
  call get_initial_fields(ug(:,:,:,1), vg(:,:,:,1), tg(:,:,:,1), psg(:,:,1), grid_tracers(:,:,:,1,:))
endif
!--------------------------------------------------------------------------------------------------------------------------------
if(dry_model) then
  call compute_pressures_and_heights(tg(:,:,:,current), psg(:,:,current), surf_geopotential, &
       z_full(:,:,:,current), z_half(:,:,:,current), p_full(:,:,:,current), p_half(:,:,:,current))
  call compute_pressures_and_heights(tg(:,:,:,previous), psg(:,:,previous), surf_geopotential, &
       z_full(:,:,:,previous), z_half(:,:,:,previous), p_full(:,:,:,previous), p_half(:,:,:,previous))
else
  call compute_pressures_and_heights(tg(:,:,:,current), psg(:,:,current), surf_geopotential, &
       z_full(:,:,:,current), z_half(:,:,:,current), p_full(:,:,:,current), p_half(:,:,:,current), &
       grid_tracers(:,:,:,current,nhum))
  call compute_pressures_and_heights(tg(:,:,:,previous), psg(:,:,previous), surf_geopotential, &
       z_full(:,:,:,previous), z_half(:,:,:,previous), p_full(:,:,:,previous), p_half(:,:,:,previous), &
       grid_tracers(:,:,:,previous,nhum))
endif

call get_deg_lon(deg_lon)
do i=is,ie
  rad_lon_2d(i,:) = deg_lon(i)*pi/180.
enddo

call get_deg_lat(deg_lat)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*pi/180.
enddo

call get_grid_boundaries(rad_lonb, rad_latb)

do i = is,ie+1
  rad_lonb_2d(i,:) = rad_lonb(i)
enddo

do j = js,je+1
  rad_latb_2d(:,j) = rad_latb(j)
enddo

if(idealized_moist_model) then
   call idealized_moist_phys_init(Time, Time_step, nhum, rad_lon_2d, rad_lat_2d, rad_lonb_2d, rad_latb_2d, tg(:,:,num_levels,current))
else
   call hs_forcing_init(get_axis_id(), Time, rad_lonb_2d, rad_latb_2d, rad_lat_2d)
endif

module_is_initialized = .true.

return
end subroutine atmosphere_init

!=================================================================================================================================

subroutine atmosphere(Time)
type(time_type), intent(in) :: Time

real    :: delta_t
type(time_type) :: Time_next

if(.not.module_is_initialized) then
  call error_mesg('atmosphere','atmosphere module is not initialized',FATAL)
endif

dt_ug  = 0.0
dt_vg  = 0.0
dt_tg  = 0.0
dt_psg = 0.0
dt_tracers = 0.0

if(current == previous) then
  delta_t = dt_real
else
  delta_t = 2*dt_real
endif

Time_next = Time + Time_step

if(idealized_moist_model) then
   call idealized_moist_phys(Time, p_half, p_full, z_half, z_full, ug, vg, tg, grid_tracers, &
                             previous, current, dt_ug, dt_vg, dt_tg, dt_tracers)
else
   call hs_forcing(1, ie-is+1, 1, je-js+1, delta_t, Time_next, rad_lon_2d, rad_lat_2d, &
                p_half(:,:,:,current ),       p_full(:,:,:,current   ), &
                    ug(:,:,:,previous),           vg(:,:,:,previous  ), &
                    tg(:,:,:,previous), grid_tracers(:,:,:,previous,:), &
                    ug(:,:,:,previous),           vg(:,:,:,previous  ), &
                    tg(:,:,:,previous), grid_tracers(:,:,:,previous,:), &
                 dt_ug(:,:,:         ),        dt_vg(:,:,:           ), &
                 dt_tg(:,:,:         ),   dt_tracers(:,:,:,:), z_full(:,:,:,current))
endif

if(previous == current) then
  future = num_time_levels + 1 - current
else
  future = previous
endif
#ifdef COLUMN_MODEL
call column(Time, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                       tg(:,:,:,future), tracer_attributes, grid_tracers(:,:,:,:,:), future, &
                       dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers, wg_full, &
                       p_full(:,:,:,current), p_half(:,:,:,current), z_full(:,:,:,current))
#else 
  call spectral_dynamics(Time, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                       tg(:,:,:,future), tracer_attributes, grid_tracers(:,:,:,:,:), future, &
                       dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers, wg_full, &
                       p_full(:,:,:,current), p_half(:,:,:,current), z_full(:,:,:,current))
#endif

if(dry_model) then
  call compute_pressures_and_heights(tg(:,:,:,future), psg(:,:,future), surf_geopotential, &
       z_full(:,:,:,future), z_half(:,:,:,future), p_full(:,:,:,future), p_half(:,:,:,future))
else
  call compute_pressures_and_heights(tg(:,:,:,future), psg(:,:,future), surf_geopotential, &
       z_full(:,:,:,future), z_half(:,:,:,future), p_full(:,:,:,future), p_half(:,:,:,future), &
                                     grid_tracers(:,:,:,future,nhum))
endif

#ifdef COLUMN_MODEL
call column_diagnostics(Time_next, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
tg(:,:,:,future), wg_full, grid_tracers(:,:,:,:,:), future)
#else
call spectral_diagnostics(Time_next, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                          tg(:,:,:,future), wg_full, grid_tracers(:,:,:,:,:), future)
#endif

previous = current
current  = future

return
end subroutine atmosphere

!=================================================================================================================================

subroutine atmosphere_end
integer :: ntr, nt
character(len=64) :: file, tr_name

if(.not.module_is_initialized) return

file='RESTART/atmosphere.res'
call nullify_domain()
call write_data(trim(file), 'time_pointers', (/real(previous),real(current)/))
do nt=1,num_time_levels
  call write_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'psg', psg(:,:,  nt), grid_domain)
  do ntr = 1,num_tracers
    tr_name = trim(tracer_attributes(ntr)%name)
    call write_data(trim(file), tr_name, grid_tracers(:,:,:,nt,ntr), grid_domain)
  enddo
enddo
call write_data(trim(file), 'wg_full', wg_full, grid_domain)

deallocate (p_half, z_half, p_full, z_full, wg_full, psg, ug, vg, tg, grid_tracers)
deallocate (dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers)
deallocate (deg_lon, rad_lon_2d, deg_lat, rad_lat_2d)

call set_domain(grid_domain)
if(idealized_moist_model) then
    call idealized_moist_phys_end
else
    call hs_forcing_end
endif
#ifdef COLUMN_MODEL
call column_end(tracer_attributes)
#else
call spectral_dynamics_end(tracer_attributes)
#endif
deallocate(tracer_attributes)

module_is_initialized = .false.

end subroutine atmosphere_end

!=================================================================================================================================

end module atmosphere_mod
