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

module water_borrowing_mod

use        fms_mod, only: mpp_pe, mpp_root_pe, write_version_number, error_mesg, FATAL

use transforms_mod, only: get_grid_domain

implicit none
private

public :: water_borrowing

character(len=128), parameter :: version = '$Id: water_borrowing.F90,v 10.0 2003/10/24 22:01:01 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'
logical :: entry_to_logfile_done = .false.

contains

subroutine water_borrowing(dt_qg, qg, current, p_half, delta_t)

real,    intent(inout), dimension(:,:,:) :: dt_qg
real,    intent(in),    dimension(:,:,:) :: qg
integer, intent(in)                      :: current
real,    intent(in),    dimension(:,:,:) :: p_half
real,    intent(in)                      :: delta_t

integer :: num_lon, num_lat, num_levels, ibeg, iend, inc, i, j, k, is, ie, js, je
real :: neighboring_water, total_water, ratio
real, dimension(size(qg,1), size(qg,3)) :: dp

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

call get_grid_domain(is, ie, js, je)
num_lon = size(qg,1); num_lat = size(qg,2); num_levels = size(qg,3)

if(is /= 1 .or. ie /= num_lon) then
  call error_mesg('water_borrowing', &
           'hole filling requires that entire latitude circle resides on same processor', FATAL)
endif

do j=1,num_lat
  dp = p_half(:,j,2:num_levels+1) - p_half(:,j,1:num_levels)
  do k=1,num_levels
    if(mod(current,2) == 0) then
      ibeg=1; iend=num_lon; inc=1
    else
      ibeg=num_lon; iend=1; inc=-1
    endif
    do i=ibeg,iend,inc
      if(qg(i,j,k) >= 0.) cycle
      if(i == 1) then
        neighboring_water = qg(num_lon,j,k)*dp(num_lon,k)
      else
        neighboring_water = qg(i-1,j,k)*dp(i-1,k)
      endif
      if(i == num_lon) then
        neighboring_water = neighboring_water + qg(1,j,k)*dp(1,k)
      else
        neighboring_water = neighboring_water + qg(i+1,j,k)*dp(i+1,k)
      endif
      if(k /= 1) then
        neighboring_water = neighboring_water + qg(i,j,k-1)*dp(i,k-1)
      endif
      if(k /= num_levels) then
        neighboring_water = neighboring_water + qg(i,j,k+1)*dp(i,k+1)
      endif
      total_water = neighboring_water + qg(i,j,k)*dp(i,k)
      if(total_water > 0.) then
        ratio = total_water/neighboring_water
        dt_qg(i,j,k) = dt_qg(i,j,k) - qg(i,j,k)/delta_t
        if(i == 1) then
          dt_qg(num_lon,j,k) = dt_qg(num_lon,j,k) + (ratio - 1)*qg(num_lon,j,k)/delta_t
        else
          dt_qg(i-1,j,k) = dt_qg(i-1,j,k) + (ratio - 1)*qg(i-1,j,k)/delta_t
        endif
        if(i == num_lon) then
          dt_qg(1,j,k)   = dt_qg(1,j,k  ) + (ratio - 1)*qg(1,j,k  )/delta_t
        else
          dt_qg(i+1,j,k) = dt_qg(i+1,j,k) + (ratio - 1)*qg(i+1,j,k)/delta_t
        endif
        if(k /= 1) then
          dt_qg(i,j,k-1) = dt_qg(i,j,k-1) + (ratio - 1)*qg(i,j,k-1)/delta_t
        endif
        if(k /= num_levels) then
          dt_qg(i,j,k+1) = dt_qg(i,j,k+1) + (ratio - 1)*qg(i,j,k+1)/delta_t
        endif
      endif
!     if(total_water < 0. .and. neighboring_water > 0.) then

!!!!!!!  Put all neighboring water into water hole and zero out neighboring water (pjp)

!       dt_qg(i,j,k) = dt_qg(i,j,k) + neighboring_water/(dp(i,k)*delta_t)
!       if(i == 1) then
!         dt_qg(num_lon,j,k) = dt_qg(num_lon,j,k) - qg(num_lon,j,k)/delta_t
!       else
!         dt_qg(i-1,j,k) = dt_qg(i-1,j,k) - qg(i-1,j,k)/delta_t
!       endif
!       if(i == num_lon) then
!         dt_qg(1,j,k)   = dt_qg(1,j,k  ) - qg(1,j,k  )/delta_t
!       else
!         dt_qg(i+1,j,k) = dt_qg(i+1,j,k) - qg(i+1,j,k)/delta_t
!       endif
!       if(k /= 1) then
!         dt_qg(i,j,k-1) = dt_qg(i,j,k-1) - qg(i,j,k-1)/delta_t
!       endif
!       if(k /= num_levels) then
!         dt_qg(i,j,k+1) = dt_qg(i,j,k+1) - qg(i,j,k+1)/delta_t
!       endif
!     endif
    enddo
  enddo
enddo

return
end subroutine water_borrowing

end module water_borrowing_mod
