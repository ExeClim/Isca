
module shallow_diagnostics_mod

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

use              fms_mod, only: write_version_number

use       transforms_mod, only: get_grid_boundaries, &
                                get_deg_lon,         &
                                get_deg_lat,         &
                                get_grid_domain,     &
                                get_spec_domain,     & 
                                grid_domain

use     diag_manager_mod, only: diag_axis_init, &
                                register_diag_field, &
                                register_static_field, &
                                send_data

use     time_manager_mod, only: time_type, &
                                get_time


use shallow_physics_mod,  only:  phys_type
use shallow_dynamics_mod, only:  grid_type


implicit none
private

public :: shallow_diagnostics_init, &
          shallow_diagnostics

character(len=84), parameter :: version = '$Id: shallow_diagnostics.F90,v 10.0 2003/10/24 22:01:02 fms Exp $'
character(len=84), parameter :: tagname = '$Name: siena_201207 $'
character(len=8) :: axiset   = 'shallow'
character(len=84) :: mod_name = 'shallow_diagnostics'

logical :: module_is_initialized = .false.

integer :: id_vor, id_stream, id_pv, id_u, id_v, id_div, id_h, id_trs, id_tr

integer :: is, ie, js, je

contains

!-----------------------------------------------------------------------------------------------------------------
subroutine shallow_diagnostics_init(Time, lon_max, lat_max)

type(time_type), intent(in) :: Time
integer, intent(in) :: lon_max, lat_max

real, dimension(lon_max  ) :: lon
real, dimension(lon_max+1) :: lonb
real, dimension(lat_max  ) :: lat
real, dimension(lat_max+1) :: latb

integer, dimension(2) :: axis_2d

integer :: log_unit, id_lonb, id_lon, id_latb, id_lat
integer :: namelist_unit, ierr, io
real    :: rad_to_deg
logical :: used

call write_version_number(version, tagname)

call get_grid_domain(is, ie, js, je)

rad_to_deg = 45./atan(1.)
call get_grid_boundaries(lonb,latb,global=.true.)
call get_deg_lon(lon)
call get_deg_lat(lat)

id_lonb=diag_axis_init('lonb', rad_to_deg*lonb, 'degrees_E', 'x', 'longitude edges', set_name=axiset, Domain2=grid_domain)
id_latb=diag_axis_init('latb', rad_to_deg*latb, 'degrees_N', 'y', 'latitude edges',  set_name=axiset, Domain2=grid_domain)
id_lon =diag_axis_init('lon', lon, 'degrees_E', 'x', 'longitude', set_name=axiset, Domain2=grid_domain, edges=id_lonb)
id_lat =diag_axis_init('lat', lat, 'degrees_N', 'y', 'latitude',  set_name=axiset, Domain2=grid_domain, edges=id_latb)

axis_2d(1) = id_lon
axis_2d(2) = id_lat

id_u      = register_diag_field(mod_name, 'ucomp' , axis_2d, Time, 'u_wind'              , 'm/s'      ) 
id_v      = register_diag_field(mod_name, 'vcomp' , axis_2d, Time, 'v_wind'              , 'm/s'      ) 
id_vor    = register_diag_field(mod_name, 'vor'   , axis_2d, Time, 'relative vorticity'  , '1/s'      )
id_div    = register_diag_field(mod_name, 'div'   , axis_2d, Time, 'divergence'          , '1/s'      )
id_h      = register_diag_field(mod_name, 'h'     , axis_2d, Time, 'geopotential'        , 'm2/s2'    )
id_pv     = register_diag_field(mod_name, 'pv'    , axis_2d, Time, 'potential vorticity' , 's/m2'     )
id_stream = register_diag_field(mod_name, 'stream', axis_2d, Time, 'streamfunction'      , 'm^2/s'    )
id_trs    = register_diag_field(mod_name, 'trs'   , axis_2d, Time, 'spectral tracer'     , 'none'     )
id_tr     = register_diag_field(mod_name, 'tr'    , axis_2d, Time, 'grid tracer'         , 'none'     )

module_is_initialized = .true.

return
end subroutine shallow_diagnostics_init

!--------------------------------------------------------------------------------------------

subroutine shallow_diagnostics(Time, Grid, Phys, time_index)

type(time_type), intent(in) :: Time
type(phys_type), intent(in) :: Phys
type(grid_type), intent(in) :: Grid
integer,         intent(in) :: time_index

logical :: used
   
if(id_u       > 0) used = send_data(id_u      , Grid%u       (:,:, time_index)     , time)
if(id_v       > 0) used = send_data(id_v      , Grid%v       (:,:, time_index)     , time)
if(id_vor     > 0) used = send_data(id_vor    , Grid%vor     (:,:, time_index)     , time)
if(id_div     > 0) used = send_data(id_div    , Grid%div     (:,:, time_index)     , time)
if(id_h       > 0) used = send_data(id_h      , Grid%h       (:,:, time_index)     , time)
if(id_pv      > 0) used = send_data(id_pv     , Grid%pv      (:,:)                 , time)
if(id_stream  > 0) used = send_data(id_stream , Grid%stream  (:,:)                 , time)
if(id_tr      > 0) used = send_data(id_tr     , Grid%tr      (:,:, time_index)     , time)
if(id_trs     > 0) used = send_data(id_trs    , Grid%trs     (:,:, time_index)     , time)

return
end subroutine shallow_diagnostics
!--------------------------------------------------------------------------------------------

end module shallow_diagnostics_mod
