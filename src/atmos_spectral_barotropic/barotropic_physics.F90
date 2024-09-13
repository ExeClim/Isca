module barotropic_physics_mod

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
                                 FATAL, WARNING,       &
                                 write_version_number, &
                                 mpp_pe,               &
                                 mpp_root_pe,          &
                                 fms_init, fms_end,    &
                                 read_data,            &
                                 write_data,           &
                                 set_domain,           &
                                 close_file,           &
                                 stdlog

use         transforms_mod, only: get_sin_lat, get_cos_lat,         &
                                  get_deg_lon, get_deg_lat,         &
                                  get_wts_lat,                      &
                                  get_grid_domain, get_spec_domain, &
                                  grid_domain

use       time_manager_mod, only: time_type

!========================================================================
implicit none
private
!========================================================================

public :: barotropic_physics_init,    &
          barotropic_physics,         &
          barotropic_physics_end,     &
          phys_type

! version information 
!========================================================================
character(len=128) :: version = '$Id: barotropic_physics.F90,v 10.0 2003/10/24 22:00:58 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'
!========================================================================

type phys_type
   real, pointer, dimension(:,:)   :: empty=>NULL()
end type

logical :: module_is_initialized = .false.

integer :: is, ie, js, je


integer :: pe
logical :: root

real, allocatable, dimension(:) :: rad_lat, &
                                   deg_lat, &
                                   sin_lat, &
                                   cos_lat, &
                                   wts_lat

! namelist 
!========================================================================

logical :: empty

namelist /barotropic_physics_nml/ empty
!========================================================================

contains

!========================================================================

subroutine barotropic_physics_init(Phys) 

type(phys_type), intent(inout) :: Phys

integer :: j, unit, ierr, io

call write_version_number (version, tagname)

pe = mpp_pe()
root = (pe == mpp_root_pe())

! read the namelist

if (file_exist('input.nml')) then
  unit = open_namelist_file ()
  ierr=1
  do while (ierr /= 0)
    read  (unit, nml=barotropic_physics_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'barotropic_physics_nml')
  enddo
  10 call close_file (unit)
endif

call get_grid_domain(is,ie,js,je)

allocate ( rad_lat      (js:je) )
allocate ( deg_lat      (js:je) )
allocate ( sin_lat      (js:je) )
allocate ( cos_lat      (js:je) )
allocate ( wts_lat      (js:je) )

call get_wts_lat(wts_lat)
call get_deg_lat(deg_lat)
rad_lat = deg_lat*atan(1.)/45. 
sin_lat = sin(rad_lat)
cos_lat = cos(rad_lat)

module_is_initialized = .true.

return
end subroutine barotropic_physics_init

!=======================================================================

subroutine barotropic_physics(Time, dt_ug, dt_vg, ug, vg,    &
                             delta_t, previous, current, Phys)

real, intent(inout),  dimension(is:ie, js:je)    :: dt_ug, dt_vg
real, intent(in)   ,  dimension(is:ie, js:je, 2) :: ug, vg

real   , intent(in)  :: delta_t
integer, intent(in)  :: previous, current

type(time_type), intent(in)    :: Time
type(phys_type), intent(inout) :: Phys

if(.not.module_is_initialized) call error_mesg('barotropic_physics', &
                       'barotropic_physics is not initialized', FATAL)

! dt_ug = dt_ug +f(ug,vg)
! dt_vg = dt_vg +f(ug,vg)
! Phys%empty = 
  
return
end subroutine barotropic_physics

!======================================================================

subroutine barotropic_physics_end(Phys)

type(phys_type), intent(in) :: Phys

if(.not.module_is_initialized) call error_mesg('barotropic_physics_end', &
                       'barotropic_physics is not initialized', FATAL)

module_is_initialized = .false.
return
end subroutine barotropic_physics_end

!======================================================================

end module barotropic_physics_mod
