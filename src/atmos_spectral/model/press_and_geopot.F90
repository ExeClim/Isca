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

module press_and_geopot_mod

! This module provides utilities that computes half- and full-
! level pressure values, given the surface pressure, and 
! geopotential heights -- assuming a vertical coordinate in which
! the half level values are p_half(k) = pk(k) + bk(k)*surface_pressure,
! where the constants pk, bk define the coordinate levels.
!
! The full-level values are given by the expression recommended by 
! Simmons and Burridge. (See Mon. Weather Review: Vol. 109, No. 4, pp. 758-766)
!  alpha  = 1.0  - p_half(k)*( ln(p_half(k+1)) - ln(p_half(k)) )/(p_half(k+1) - p_half(k))
!  ln(p_full(k)) = ln(p_half(k+1)) - alpha
!
! Geopotentials are computed by assuming isothermal temperatures and in each layer
! integrating the hydrostatic/ideal gas equations exactly

use       fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                         write_version_number

use constants_mod, only: grav, rdgas, rvgas

implicit none

private

public :: press_and_geopot_init, press_and_geopot_end, pressure_variables, half_level_pressures
public :: compute_geopotential, compute_pressures_and_heights, compute_z_bot
 
interface half_level_pressures
   module procedure   half_level_pressures_1d, &
                      half_level_pressures_3d
end interface
 
interface pressure_variables
   module procedure   pressure_variables_1d, &
                      pressure_variables_2d, &
                      pressure_variables_3d
end interface

interface compute_pressures_and_heights
   module procedure   compute_pressures_and_heights_1d, &
                      compute_pressures_and_heights_2d, &
                      compute_pressures_and_heights_3d
end interface

!===============================================================================================

character(len=128), parameter :: version = '$Id: press_and_geopot.F90,v 17.0 2009/07/21 03:00:50 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

!===============================================================================================
 
real, allocatable, dimension(:) :: pk, bk
 
real    :: ln_top_level_factor
integer :: num_levels
logical :: use_virtual_temperature
character(len=64) :: vert_difference_option
 
logical :: module_is_initialized = .false.

contains
 
!------------------------------------------------------------------------------
 
subroutine press_and_geopot_init(pk_in, bk_in, use_virtual_temperature_in, vert_difference_option_in)
 
real,    intent (in), dimension(:)   :: pk_in, bk_in
logical, intent (in)                 :: use_virtual_temperature_in
character(len=*), intent (in)        :: vert_difference_option_in

if(module_is_initialized) return

call write_version_number(version, tagname)

num_levels = size(pk_in,1)-1

allocate(pk (num_levels+1))
allocate(bk (num_levels+1))
pk = pk_in
bk = bk_in

vert_difference_option  = vert_difference_option_in
use_virtual_temperature = use_virtual_temperature_in

ln_top_level_factor = -1.0

module_is_initialized = .true.
 
return
end subroutine press_and_geopot_init
 
!-----------------------------------------------------------------------
 
function half_level_pressures_3d(surface_p) result(p_half)
 
real, intent (in),  dimension(:,:)     :: surface_p
real, dimension(size(surface_p,1),size(surface_p,2),num_levels+1) :: p_half

integer :: k
 
if(.not.module_is_initialized) then
  call error_mesg('half_level_pressures','press_and_geopot_init has not been called', FATAL)
end if

do k=1,num_levels+1
  p_half(:,:,k) = pk(k) + bk(k)*surface_p(:,:)
end do
 
return
end function half_level_pressures_3d
 
!-----------------------------------------------------------------------
 
function half_level_pressures_1d(surface_p) result(p_half)
 
real, intent (in)             :: surface_p
real, dimension(num_levels+1) :: p_half

integer :: k
 
do k=1,num_levels+1
  p_half(k) = pk(k)+bk(k)*surface_p
end do
 
return
end function half_level_pressures_1d
 
!-------------------------------------------------------------------------------------
 
subroutine pressure_variables_3d(p_half, ln_p_half, p_full, ln_p_full, surface_p)
 
real, intent (out), dimension(:,:,:) :: p_half, ln_p_half, p_full, ln_p_full
real, intent (in),  dimension(:,:)   :: surface_p
 
real, dimension(size(p_half,1),size(p_half,2)) :: alpha

integer :: k

if(.not.module_is_initialized) then
  call error_mesg('pressure_variables','press_and_geopot_init has not been called', FATAL)
end if

p_half = half_level_pressures(surface_p)

if(trim(vert_difference_option) == 'simmons_and_burridge') then

  if(pk(1).eq.0.0 .and. bk(1).eq.0.0) then
 
    do k=2,size(p_half,3)
      ln_p_half(:,:,k) = log(p_half(:,:,k))
    end do
 
    do k=2,size(p_half,3)-1
      alpha  = 1.0  - p_half(:,:,k)*(ln_p_half(:,:,k+1) - ln_p_half(:,:,k))/(p_half(:,:,k+1) - p_half(:,:,k))
      ln_p_full(:,:,k) = ln_p_half(:,:,k+1) - alpha
    end do
    ln_p_full(:,:,1) = ln_p_half(:,:,2) + ln_top_level_factor
    ln_p_half(:,:,1) = 0.0
 
  else
 
    do k=1,size(p_half,3)
      ln_p_half(:,:,k) = log(p_half(:,:,k))
    end do
 
    do k=1,size(p_half,3)-1
      alpha  = 1.0  - p_half(:,:,k)*(ln_p_half(:,:,k+1) - ln_p_half(:,:,k))/(p_half(:,:,k+1) - p_half(:,:,k))
      ln_p_full(:,:,k) = ln_p_half(:,:,k+1) - alpha
    end do
 
  end if
  p_full = exp(ln_p_full)

else if(trim(vert_difference_option) == 'mcm') then

  do k = 1,size(p_full,3)
    p_full   (:,:,k) = 0.5*(p_half(:,:,k+1) + p_half(:,:,k))
    ln_p_full(:,:,k) = log(p_full(:,:,k))
  end do
  if(pk(1).eq.0.0 .and. bk(1).eq.0.0) then
    do k=2,size(p_half,3)
      ln_p_half(:,:,k) = log(p_half(:,:,k))
    end do
    ln_p_half(:,:,1) = 0.
  else
    do k=1,size(p_half,3)
      ln_p_half(:,:,k) = log(p_half(:,:,k))
    end do
  end if

else

  call error_mesg('pressure_variables','"'//trim(vert_difference_option)//'"'// &
                  ' is not a valid value for vert_difference_option', FATAL)

end if
 
return
end subroutine pressure_variables_3d
 
!-------------------------------------------------------------------------------------
subroutine compute_z_bot(psg, tg, z_bot, qg)
real, intent(in ), dimension(:,:) :: psg, tg
real, intent(out), dimension(:,:) :: z_bot
real, intent(in ), optional, dimension(:,:) :: qg

real, dimension(size(psg,1), size(psg,2)) ::    p_half_bot,    p_half_nxt
real, dimension(size(psg,1), size(psg,2)) :: ln_p_half_bot, ln_p_half_nxt
real, dimension(size(psg,1), size(psg,2)) :: ln_p_full_bot, alpha, virtual_t

num_levels = size(pk,1) - 1
p_half_bot = pk(num_levels+1) + psg*bk(num_levels+1)
p_half_nxt = pk(num_levels)   + psg*bk(num_levels)
ln_p_half_bot = log(p_half_bot)

if(trim(vert_difference_option) == 'simmons_and_burridge') then
  ln_p_half_nxt = log(p_half_nxt)
  alpha  = 1.0  - p_half_nxt*(ln_p_half_bot - ln_p_half_nxt)/(p_half_bot - p_half_nxt)
  ln_p_full_bot = ln_p_half_bot - alpha
else if(trim(vert_difference_option) == 'mcm') then
  ln_p_full_bot = log(.5*(p_half_bot + p_half_nxt))
endif

if(use_virtual_temperature) then
  if(present(qg)) then
    virtual_t = tg*(1. + (rvgas/rdgas - 1.)*qg)
  else
    call error_mesg('compute_z_bot','qg must be present when use_virtual_temperature=.true.', FATAL)
  endif
else
  virtual_t = tg
endif

z_bot = (rdgas*virtual_t*(ln_p_half_bot - ln_p_full_bot))/grav

return
end subroutine compute_z_bot
!-------------------------------------------------------------------------------------
 
subroutine pressure_variables_1d(p_half,ln_p_half,p_full,ln_p_full,surface_p)
 
real, intent (out), dimension(:) :: p_half, ln_p_half, p_full, ln_p_full
real, intent (in)   :: surface_p
 
real, dimension(1,1)              :: surface_p_2d
real, dimension(1,1,num_levels+1) :: p_half_3d
real, dimension(1,1,num_levels+1) :: ln_p_half_3d
real, dimension(1,1,num_levels)   :: p_full_3d
real, dimension(1,1,num_levels)   :: ln_p_full_3d
 
 
surface_p_2d(1,1) = surface_p
 
call pressure_variables_3d(p_half_3d,ln_p_half_3d,p_full_3d,ln_p_full_3d,surface_p_2d)
 
p_half    =    p_half_3d(1,1,:)
ln_p_half = ln_p_half_3d(1,1,:)
p_full    =    p_full_3d(1,1,:)
ln_p_full = ln_p_full_3d(1,1,:)
 
return
end subroutine pressure_variables_1d
 
!-----------------------------------------------------------------------

subroutine pressure_variables_2d(p_half,ln_p_half,p_full,ln_p_full,surface_p)
 
real, intent (out), dimension(:,:) :: p_half, ln_p_half, p_full, ln_p_full
real, intent (in) , dimension(:)   :: surface_p
 
real, dimension(1,size(surface_p))              :: surface_p_2d
real, dimension(1,size(surface_p),num_levels+1) :: p_half_3d
real, dimension(1,size(surface_p),num_levels+1) :: ln_p_half_3d
real, dimension(1,size(surface_p),num_levels)   :: p_full_3d
real, dimension(1,size(surface_p),num_levels)   :: ln_p_full_3d
 
 
surface_p_2d(1,:) = surface_p
 
call pressure_variables_3d(p_half_3d,ln_p_half_3d,p_full_3d,ln_p_full_3d,surface_p_2d)
 
p_half    =    p_half_3d(1,:,:)
ln_p_half = ln_p_half_3d(1,:,:)
p_full    =    p_full_3d(1,:,:)
ln_p_full = ln_p_full_3d(1,:,:)
 
return
end subroutine pressure_variables_2d

!-----------------------------------------------------------------------

subroutine compute_geopotential(t_grid, ln_p_half, ln_p_full, surf_geopotential, geopot_full, geopot_half, q_grid)

real, intent (in), dimension (:,:,:) :: t_grid, ln_p_half, ln_p_full
real, intent (in), dimension (:,:)   :: surf_geopotential
real, intent(out), dimension (:,:,:) :: geopot_full, geopot_half
real, intent (in), optional, dimension (:,:,:) :: q_grid

real, dimension (size(t_grid,1),size(t_grid,2),size(t_grid,3)  ) :: virtual_t

integer :: ktop, num_levels, k

if(.not.module_is_initialized) then
  call error_mesg('compute_geopotential','press_and_geopot_init has not been called', FATAL)
end if

num_levels = size(t_grid,3)

geopot_half(:,:,num_levels+1) = surf_geopotential

if(pk(1).eq.0.0) then
  ktop = 2
  geopot_half(:,:,1) = 0.0
else
  ktop = 1
endif

if(use_virtual_temperature) then
  if(present(q_grid)) then
    virtual_t = t_grid*(1. + (rvgas/rdgas - 1.)*q_grid)
  else
    call error_mesg('compute_geopotential','q_grid must be present when use_virtual_temperature=.true.', FATAL)
  endif
else
  virtual_t = t_grid
endif

do k=num_levels,ktop,-1
  geopot_half(:,:,k) = geopot_half(:,:,k+1) + rdgas*virtual_t(:,:,k)*(ln_p_half(:,:,k+1) - ln_p_half(:,:,k))
enddo

do k=1,num_levels
  geopot_full(:,:,k) = geopot_half(:,:,k+1) + rdgas*virtual_t(:,:,k)*(ln_p_half(:,:,k+1) - ln_p_full(:,:,k))
enddo

return
end subroutine compute_geopotential

!-----------------------------------------------------------------------

subroutine compute_pressures_and_heights_3d(t_grid, ps_grid, surf_geopotential, z_full, z_half, p_full, p_half, q_grid)

real, intent (in), dimension (:,:,:) :: t_grid
real, intent (in), dimension (:,:  ) :: ps_grid, surf_geopotential
real, intent(in), optional, dimension(:,:,:) :: q_grid

real, intent(out), dimension (size(t_grid,1),size(t_grid,2),size(t_grid,3)  ) :: z_full, p_full
real, intent(out), dimension (size(t_grid,1),size(t_grid,2),size(t_grid,3)+1) :: z_half, p_half

real, dimension(size(t_grid,1),size(t_grid,2),size(t_grid,3)  ) :: ln_p_full
real, dimension(size(t_grid,1),size(t_grid,2),size(t_grid,3)+1) :: ln_p_half

if(.not.module_is_initialized) then
  call error_mesg('compute_pressures_and_heights','press_and_geopot_init has not been called', FATAL)
end if

call pressure_variables (p_half, ln_p_half, p_full, ln_p_full, ps_grid)

call compute_geopotential(t_grid, ln_p_half, ln_p_full, surf_geopotential, z_full, z_half, q_grid)

z_full = z_full/grav
z_half = z_half/grav

return
end subroutine compute_pressures_and_heights_3d

!-----------------------------------------------------------------------
subroutine compute_pressures_and_heights_1d(t_grid, ps_grid, surf_geopotential, z_full, z_half, p_full, p_half, q_grid)

real, intent (in), dimension (:) :: t_grid
real, intent (in)                :: ps_grid, surf_geopotential
real, intent(in), optional, dimension(:) :: q_grid

real, intent(out), dimension (size(t_grid)  ) :: z_full, p_full
real, intent(out), dimension (size(t_grid)+1) :: z_half, p_half
 
real, dimension(1,1) :: ps_grid_2d, surf_geopotential_2d
real, dimension(1,1,size(t_grid))   :: t_grid_3d, z_full_3d, p_full_3d, q_grid_3d
real, dimension(1,1,size(t_grid)+1) :: z_half_3d, p_half_3d
 
t_grid_3d(1,1,:) = t_grid
ps_grid_2d(1,1)  = ps_grid
surf_geopotential_2d(1,1) = surf_geopotential
                      
if(present(q_grid)) then
  q_grid_3d(1,1,:) = q_grid
  call compute_pressures_and_heights_3d(t_grid_3d, ps_grid_2d, surf_geopotential_2d, z_full_3d, z_half_3d, p_full_3d, p_half_3d, &
                                        q_grid_3d)
else
  call compute_pressures_and_heights_3d(t_grid_3d, ps_grid_2d, surf_geopotential_2d, z_full_3d, z_half_3d, p_full_3d, p_half_3d)
endif

z_full = z_full_3d(1,1,:)
z_half = z_half_3d(1,1,:)
p_full = p_full_3d(1,1,:)
p_half = p_half_3d(1,1,:)

end subroutine compute_pressures_and_heights_1d

!-----------------------------------------------------------------------
subroutine compute_pressures_and_heights_2d(t_grid, ps_grid, surf_geopotential, z_full, z_half, p_full, p_half, q_grid)

real, intent (in), dimension (:,:) :: t_grid
real, intent (in), dimension (:)   :: ps_grid, surf_geopotential
real, intent(in), optional, dimension(:,:) :: q_grid

real, intent(out), dimension (size(t_grid,1),size(t_grid,2)  ) :: z_full, p_full
real, intent(out), dimension (size(t_grid,1),size(t_grid,2)+1) :: z_half, p_half

real, dimension(1,size(t_grid,1)) :: ps_grid_2d, surf_geopotential_2d
real, dimension(1,size(t_grid,1),size(t_grid,2))   :: t_grid_3d, z_full_3d, p_full_3d, q_grid_3d
real, dimension(1,size(t_grid,1),size(t_grid,2)+1) :: z_half_3d, p_half_3d
 
t_grid_3d(1,:,:) = t_grid
ps_grid_2d(1,:)  = ps_grid
surf_geopotential_2d(1,:) = surf_geopotential
                      
if(present(q_grid)) then
  q_grid_3d(1,:,:) = q_grid
  call compute_pressures_and_heights_3d(t_grid_3d, ps_grid_2d, surf_geopotential_2d, z_full_3d, z_half_3d, p_full_3d, p_half_3d, &
                                        q_grid_3d)
else
  call compute_pressures_and_heights_3d(t_grid_3d, ps_grid_2d, surf_geopotential_2d, z_full_3d, z_half_3d, p_full_3d, p_half_3d)
endif

z_full = z_full_3d(1,:,:)
z_half = z_half_3d(1,:,:)
p_full = p_full_3d(1,:,:)
p_half = p_half_3d(1,:,:)

end subroutine compute_pressures_and_heights_2d

!-----------------------------------------------------------------------
subroutine press_and_geopot_end

if(.not.module_is_initialized) return

deallocate(pk, bk)
module_is_initialized = .false.

return
end subroutine press_and_geopot_end
!-----------------------------------------------------------------------

end module press_and_geopot_mod
