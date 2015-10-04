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

module strat_cloud_mod

use  time_manager_mod, only : time_type
use rad_utilities_mod, only : aerosol_type

implicit none
private
public strat_cloud_init, strat_cloud, strat_cloud_new, strat_cloud_end, strat_cloud_sum, strat_cloud_avg, do_strat_cloud, strat_cloud_restart
Character(len=128) :: Version = '$Id: strat_cloud.F90,v 1.1.2.1 2013/01/24 13:47:35 pjp Exp $'
Character(len=128) :: Tagname = '$Name:  $'

logical, public :: strat_cloud_on = .false.

CONTAINS
!==================================================================================================================
subroutine strat_cloud_init (axes, Time, idim, jdim, kdim, do_legacy_strat_cloud)
integer, intent(in) :: idim, jdim, kdim, axes(4)
type(time_type), intent(in) :: Time
logical, intent(in), optional :: do_legacy_strat_cloud

return
end subroutine strat_cloud_init
!=================================================================================================================
subroutine strat_cloud (Time, is, ie, js, je, dtcloud, pfull, phalf, radturbten2, T, qv, ql, qi ,qa, omega, Mc, diff_t, &
   LAND, ST, SQ, SL, SI, SA, f_snow_berg, rain3d, snow3d, snowclr3d, surfrain, surfsnow, qrat, ahuco, limit_conv_cloud_frac, &
   MASK, qn, Aerosol, SN)
type(time_type), intent (in) :: Time
integer, intent (in) :: is, ie, js, je
real, intent (in) :: dtcloud
real, dimension(:,:,:), intent (in) :: pfull, phalf, T, qv, ql, qi, qa, omega, Mc, diff_t, qrat, ahuco, radturbten2
logical, intent(in) :: limit_conv_cloud_frac
real, dimension(:,:), intent (in) :: LAND
real, dimension(:,:,:), intent (out) :: ST, SQ, SL, SI, SA, rain3d, snow3d, snowclr3d, f_snow_berg
real, dimension(:,:), intent (out) :: surfrain, surfsnow
real, dimension(:,:,:), intent (in), optional :: MASK, qn
type(aerosol_type), intent (in), optional :: Aerosol 
real, dimension(:,:,:), intent (out), optional :: SN

ST = 0.0
SQ = 0.0
SL = 0.0
SI = 0.0
SA = 0.0
rain3d = 0.0
snow3d = 0.0
snowclr3d = 0.0
f_snow_berg = 0.0
surfrain = 0.0
surfsnow = 0.0
if(present(SN)) SN = 0.0
return 
end subroutine strat_cloud
!=================================================================================================================
subroutine strat_cloud_new (Time, is, ie, js, je, dtcloud, pfull, phalf, zhalf, zfull, radturbten2, T_in, qv_in, ql_in, qi_in, &
   qa_in, omega, Mc, diff_t, LAND, ST_out, SQ_out, SL_out, SI_out, SA_out, f_snow_berg, rain3d, snow3d, snowclr3d, surfrain, &
   surfsnow, qrat, ahuco, limit_conv_cloud_frac, Aerosol, MASK3d, qn_in, SN_out, qni_in, SNi_out, lsc_snow, lsc_rain, &
   lsc_snow_size, lsc_rain_size )
type(time_type), intent (in) :: Time
integer, intent (in) :: is,ie,js,je
real, intent (in) :: dtcloud
real, dimension(:,:,:), intent (in) :: pfull, phalf, zhalf, zfull, T_in, qv_in, ql_in
real, dimension(:,:,:), intent (in) :: qi_in, qa_in, omega,Mc, diff_t, qrat, ahuco, radturbten2
logical, intent(in) :: limit_conv_cloud_frac
type(aerosol_type), intent (in) :: Aerosol 
real, dimension(:,:), intent (in) :: LAND
real, dimension(:,:,:), intent (out) :: ST_out, SQ_out, SL_out, SI_out, SA_out, rain3d, snow3d, snowclr3d
real, dimension(:,:,:), intent (out) :: f_snow_berg 
real, dimension(:,:), intent (out) :: surfrain,surfsnow
real, dimension(:,:,:), intent (in), optional :: MASK3d, qn_in, qni_in
real, dimension(:,:,:), intent (out), optional :: SN_out, SNi_out, lsc_snow, lsc_rain, lsc_snow_size, lsc_rain_size

ST_out = 0.0
SQ_out = 0.0
SL_out = 0.0
SI_out = 0.0
SA_out = 0.0
rain3d = 0.0
snow3d = 0.0
snowclr3d = 0.0
f_snow_berg = 0.0
surfrain = 0.0
surfsnow = 0.0
if(present(SN_out)) SN_out = 0.0
if(present(SNi_out)) SNi_out = 0.0
if(present(lsc_snow)) lsc_snow = 0.0
if(present(lsc_rain)) lsc_rain = 0.0
if(present(lsc_snow_size)) lsc_snow_size = 0.0
if(present(lsc_rain_size)) lsc_rain_size = 0.0
return
end subroutine strat_cloud_new
!=================================================================================================================
subroutine strat_cloud_end()

return
end subroutine strat_cloud_end
!=================================================================================================================
subroutine strat_cloud_sum (is, js, ql, qi, cf)
integer, intent(in) :: is, js
real, dimension(:,:,:), intent(in) :: ql, qi, cf

return
end subroutine strat_cloud_sum
!=================================================================================================================
subroutine strat_cloud_avg (is, js, ql, qi, cf, ierr)
integer, intent(in) :: is, js
real, dimension(:,:,:), intent(out) :: ql, qi, cf
integer, intent(out) :: ierr

ierr = 0
return
end subroutine strat_cloud_avg
!=================================================================================================================
function do_strat_cloud ( ) result (answer)
logical :: answer

answer = .false.
end function do_strat_cloud
!=================================================================================================================
subroutine strat_cloud_restart(timestamp)
character(len=*), intent(in), optional :: timestamp

return
end subroutine strat_cloud_restart
!=================================================================================================================
end module strat_cloud_mod
