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

module spectral_initialize_fields_mod

use              fms_mod, only: mpp_pe, mpp_root_pe, write_version_number

use        constants_mod, only: rdgas

use       transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, vor_div_from_uv_grid, &
                                uv_grid_from_vor_div, get_grid_domain, get_spec_domain, area_weighted_global_mean

implicit none
private

public :: spectral_initialize_fields

character(len=128), parameter :: version = &
'$Id: spectral_initialize_fields.F90,v 17.0 2009/07/21 03:00:42 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: siena_201211 $'

contains

!-------------------------------------------------------------------------------------------------
subroutine spectral_initialize_fields(reference_sea_level_press, triang_trunc, initial_temperature, &
                        surf_geopotential, ln_ps, vors, divs, ts, psg, ug, vg, tg, vorg, divg)

real,    intent(in) :: reference_sea_level_press
logical, intent(in) :: triang_trunc
real,    intent(in) :: initial_temperature

real,    intent(in),  dimension(:,:    ) :: surf_geopotential
complex, intent(out), dimension(:,:    ) :: ln_ps
complex, intent(out), dimension(:,:,:  ) :: vors, divs, ts
real,    intent(out), dimension(:,:    ) :: psg
real,    intent(out), dimension(:,:,:  ) :: ug, vg, tg
real,    intent(out), dimension(:,:,:  ) :: vorg, divg

real, allocatable, dimension(:,:) :: ln_psg

real :: initial_sea_level_press, global_mean_psg
real :: initial_perturbation   = 1.e-7

integer :: ms, me, ns, ne, is, ie, js, je, num_levels

call write_version_number(version, tagname)

num_levels = size(ug,3)
call get_grid_domain(is, ie, js, je)
call get_spec_domain(ms, me, ns, ne)
allocate(ln_psg(is:ie, js:je))

initial_sea_level_press = reference_sea_level_press  

ug      = 0.
vg      = 0.
tg      = 0.
psg     = 0.
vorg    = 0.
divg    = 0.

vors  = (0.,0.)
divs  = (0.,0.)
ts    = (0.,0.)
ln_ps = (0.,0.)

tg     = initial_temperature
ln_psg = log(initial_sea_level_press) - surf_geopotential/(rdgas*initial_temperature)
psg    = exp(ln_psg)

! initial vorticity perturbation used in benchmark code
if(ms <= 1 .and. me >= 1 .and. ns <= 3 .and. ne >= 3) then
  vors(2-ms,4-ns,num_levels  ) = initial_perturbation
  vors(2-ms,4-ns,num_levels-1) = initial_perturbation
  vors(2-ms,4-ns,num_levels-2) = initial_perturbation
endif
if(ms <= 5 .and. me >= 5 .and. ns <= 3 .and. ne >= 3) then
  vors(6-ms,4-ns,num_levels  ) = initial_perturbation
  vors(6-ms,4-ns,num_levels-1) = initial_perturbation
  vors(6-ms,4-ns,num_levels-2) = initial_perturbation
endif
if(ms <= 1 .and. me >= 1 .and. ns <= 2 .and. ne >= 2) then
  vors(2-ms,3-ns,num_levels  ) = initial_perturbation
  vors(2-ms,3-ns,num_levels-1) = initial_perturbation
  vors(2-ms,3-ns,num_levels-2) = initial_perturbation
endif
if(ms <= 5 .and. me >= 5 .and. ns <= 2 .and. ne >= 2) then
  vors(6-ms,3-ns,num_levels  ) = initial_perturbation
  vors(6-ms,3-ns,num_levels-1) = initial_perturbation
  vors(6-ms,3-ns,num_levels-2) = initial_perturbation
endif
call uv_grid_from_vor_div(vors, divs, ug, vg)

!  initial spectral fields (and spectrally-filtered) grid fields

call trans_grid_to_spherical(tg, ts)
call trans_spherical_to_grid(ts, tg)

call trans_grid_to_spherical(ln_psg, ln_ps)
call trans_spherical_to_grid(ln_ps,  ln_psg)
psg = exp(ln_psg)

call vor_div_from_uv_grid(ug, vg, vors, divs, triang=triang_trunc)
call uv_grid_from_vor_div(vors, divs, ug, vg)
call trans_spherical_to_grid(vors, vorg)
call trans_spherical_to_grid(divs, divg)

!  compute and print mean surface pressure
global_mean_psg = area_weighted_global_mean(psg)
if(mpp_pe() == mpp_root_pe()) then
  print '("mean surface pressure=",f9.4," mb")',.01*global_mean_psg
endif

return
end subroutine spectral_initialize_fields
!================================================================================

end module spectral_initialize_fields_mod
