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

module column_initialize_fields_mod

use              fms_mod, only: mpp_pe, mpp_root_pe, write_version_number

use      column_grid_mod, only: area_weighted_global_mean

use        constants_mod, only: rdgas

use         spec_mpp_mod, only: get_grid_domain

implicit none
private

public :: column_initialize_fields

character(len=128), parameter :: version = '$Id: column_initialize_fields.F90,v 0.1 2018/14/11 HH:MM:SS isca Exp $'
character(len=128), parameter :: tagname = '$Name: isca_201811 $'

contains

!-------------------------------------------------------------------------------------------------
subroutine column_initialize_fields(reference_sea_level_press, initial_temperature, surface_wind, &
                                    surf_geopotential, psg, ug, vg, tg)

real,    intent(in) :: reference_sea_level_press
real,    intent(in) :: initial_temperature
real,    intent(in) :: surface_wind

real,    intent(in),  dimension(:,:    ) :: surf_geopotential
real,    intent(out), dimension(:,:    ) :: psg
real,    intent(out), dimension(:,:,:  ) :: ug, vg, tg

real, allocatable, dimension(:,:) :: ln_psg

real :: initial_sea_level_press, global_mean_psg
real :: initial_perturbation   = 1.e-7

integer :: is, ie, js, je, num_levels

call write_version_number(version, tagname)

num_levels = size(ug,3)
call get_grid_domain(is, ie, js, je)
allocate(ln_psg(is:ie, js:je))

initial_sea_level_press = reference_sea_level_press  

ug      = 0.
vg      = 0.
tg      = 0.
psg     = 0.

ug(:,:,num_levels) = surface_wind / sqrt(2.)
vg(:,:,num_levels) = surface_wind / sqrt(2.) 

tg     = initial_temperature
ln_psg = log(initial_sea_level_press) - surf_geopotential/(rdgas*initial_temperature)
psg    = exp(ln_psg)


!  compute and print mean surface pressure
global_mean_psg = area_weighted_global_mean(psg)
if(mpp_pe() == mpp_root_pe()) then
  print '("mean surface pressure=",f9.4," mb")',.01*global_mean_psg
endif

return
end subroutine column_initialize_fields
!================================================================================

end module column_initialize_fields_mod
