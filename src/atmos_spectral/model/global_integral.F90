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

module global_integral_mod

use              fms_mod, only: mpp_pe, mpp_root_pe, &
                                write_version_number

use press_and_geopot_mod, only: half_level_pressures

use       transforms_mod, only: area_weighted_global_mean

use        constants_mod, only: grav

use      mpp_domains_mod, only: mpp_global_field

implicit none
private

public :: mass_weighted_global_integral

real :: global_sum_of_wts
logical :: entry_to_logfile_done=.false.
character(len=128), parameter :: version = '$Id: global_integral.F90,v 13.0 2006/03/28 21:17:51 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

contains

!---------------------------------------------------------------------------------------------

function mass_weighted_global_integral(field, surf_press)

!  This function returns the mass weighted vertical integral of field,
!  area averaged over the globe. The units of the result are:
!  (units of field)*(Kg/meters**2)

real :: mass_weighted_global_integral
real, intent(in), dimension(:,:,:) :: field
real, intent(in), dimension(:,:)   :: surf_press
real, dimension(size(field,1), size(field,2), size(field,3)  ) :: dp
real, dimension(size(field,1), size(field,2), size(field,3)+1) :: p_half

real, dimension(size(field,1), size(field,2)) :: vert_integral

integer ::  k, num_levels

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done=.true.
endif


num_levels = size(field,3)
p_half= half_level_pressures(surf_press)
dp = p_half(:,:,2:num_levels+1) - p_half(:,:,1:num_levels)
vert_integral = 0.
do k=1,num_levels
  vert_integral = vert_integral + field(:,:,k)*dp(:,:,k)
enddo
mass_weighted_global_integral = area_weighted_global_mean(vert_integral)/grav

return
end function mass_weighted_global_integral
!---------------------------------------------------------------------------------------------
end module global_integral_mod
