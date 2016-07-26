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

module leapfrog_mod

use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number

!===================================================================================================
implicit none
private
!===================================================================================================

interface leapfrog
  module procedure leapfrog_2d_complex, leapfrog_3d_complex
end interface

interface leapfrog_2level_A
  module procedure leapfrog_2level_A_2d_complex, &
                   leapfrog_2level_A_3d_complex, &
                   leapfrog_2level_A_3d_real
end interface

interface leapfrog_2level_B
  module procedure leapfrog_2level_B_2d_complex, &
                   leapfrog_2level_B_3d_complex, &
                   leapfrog_2level_B_3d_real
end interface

character(len=128), parameter :: version = '$Id leapfrog.f90 $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

public :: leapfrog, leapfrog_2level_A, leapfrog_2level_B

logical :: entry_to_logfile_done = .false.

contains

!================================================================================

subroutine leapfrog_2level_A_3d_complex (a, dt_a, previous, current, future, delta_t, robert_coeff, raw_filter_coeff, prev_curr_part_raw_filter)

complex, intent(inout), dimension(:,:,:,:) :: a
complex, intent(in),    dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_filter_coeff

complex, intent(out), dimension(size(dt_a,1),size(dt_a,2),size(dt_a,3)) :: prev_curr_part_raw_filter


if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

prev_curr_part_raw_filter = a(:,:,:,previous) - 2.0*a(:,:,:,current)

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(prev_curr_part_raw_filter)*raw_filter_coeff
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(prev_curr_part_raw_filter)*raw_filter_coeff
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
endif

return
end subroutine leapfrog_2level_A_3d_complex

!================================================================================

subroutine leapfrog_2level_B_3d_complex (a, part_filt_a, current, future, robert_coeff, raw_filter_coeff)

complex, intent(inout), dimension(:,:,:,:) :: a
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff, raw_filter_coeff

complex, intent(in), dimension(:,:,:) :: part_filt_a

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a(:,:,:,current) = a(:,:,:,current) + robert_coeff * a(:,:,:,future)               * raw_filter_coeff
a(:,:,:,future)  = a(:,:,:,future)  + robert_coeff * (part_filt_a+a(:,:,:,future)) * (raw_filter_coeff-1.0)

return
end subroutine leapfrog_2level_B_3d_complex

!================================================================================

subroutine leapfrog_2level_A_3d_real (a, dt_a, previous, current, future, delta_t, robert_coeff, raw_filter_coeff, prev_curr_part_raw_filter)

real, intent(inout), dimension(:,:,:,:) :: a
real, intent(in),    dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_filter_coeff

real, intent(out), dimension(size(dt_a,1),size(dt_a,2),size(dt_a,3)) :: prev_curr_part_raw_filter


if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

prev_curr_part_raw_filter = a(:,:,:,previous) - 2.0*a(:,:,:,current)

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(prev_curr_part_raw_filter)
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(prev_curr_part_raw_filter)
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
endif

return
end subroutine leapfrog_2level_A_3d_real

!================================================================================

subroutine leapfrog_2level_B_3d_real (a, part_filt_a, current, future, robert_coeff, raw_filter_coeff)

real, intent(inout), dimension(:,:,:,:) :: a
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff, raw_filter_coeff

real, intent(in), dimension(:,:,:) :: part_filt_a


if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a(:,:,:,current) = a(:,:,:,current) + robert_coeff * a(:,:,:,future) * raw_filter_coeff
a(:,:,:,future)  = a(:,:,:,future)  + robert_coeff * (part_filt_a+ a(:,:,:,future)) * (raw_filter_coeff-1.0)

return
end subroutine leapfrog_2level_B_3d_real


!================================================================================

subroutine leapfrog_2level_A_2d_complex (a, dt_a, previous, current, future, delta_t, robert_coeff, raw_filter_coeff, prev_curr_part_raw_filter)

complex, intent(inout), dimension(:,:,:) :: a
complex, intent(in),    dimension(:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_filter_coeff

complex, intent(out), dimension(size(dt_a,1),size(dt_a,2)) :: prev_curr_part_raw_filter

complex,              dimension(size(a,1),size(a,2),1) :: prev_curr_part_raw_filter_tmp


complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d
complex, dimension(size(a,1),size(a,2),1)           :: dt_a_3d

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
dt_a_3d(:,:,1) = dt_a
call leapfrog_2level_A_3d_complex(a_3d, dt_a_3d, previous, current, future, delta_t, robert_coeff, raw_filter_coeff, prev_curr_part_raw_filter_tmp)
a = a_3d(:,:,1,:)
prev_curr_part_raw_filter=prev_curr_part_raw_filter_tmp(:,:,1)

end subroutine leapfrog_2level_A_2d_complex
!================================================================================

subroutine leapfrog_2level_B_2d_complex (a, part_filt_a, current, future, robert_coeff, raw_filter_coeff)

complex, intent(inout), dimension(:,:,:) :: a
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff, raw_filter_coeff
complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d
complex, intent(in), dimension(:,:) :: part_filt_a

complex, dimension(size(a,1),size(a,2),1) :: part_filt_a_3d



if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
part_filt_a_3d(:,:,1) = part_filt_a

call leapfrog_2level_B_3d_complex (a_3d, part_filt_a_3d, current, future, robert_coeff, raw_filter_coeff)
a = a_3d(:,:,1,:)

end subroutine leapfrog_2level_B_2d_complex
!================================================================================

subroutine leapfrog_3d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff, raw_filter_coeff)

complex, intent(inout), dimension(:,:,:,:) :: a
complex, intent(in),    dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_filter_coeff

complex, dimension(size(dt_a,1),size(dt_a,2),size(dt_a,3)) :: prev_curr_part_raw_filter

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

prev_curr_part_raw_filter=a(:,:,:,previous) - 2.0*a(:,:,:,current) !st Defined at the start to get unmodified value of a(:,:,:,current).

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff * (prev_curr_part_raw_filter + a(:,:,:,future ))*raw_filter_coeff
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff * (prev_curr_part_raw_filter                   )*raw_filter_coeff
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t * dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff * a(:,:,:,future)*raw_filter_coeff
endif

a(:,:,:,future ) = a(:,:,:,future ) + robert_coeff * (prev_curr_part_raw_filter + a(:,:,:,future )) * (raw_filter_coeff-1.0) 

!st RAW filter (see e.g. Williams 2011 10.1175/2010MWR3601.1) conserves 3-time-level mean in leap-frog integrations, improving amplitude accuracy of leap-frog scheme from first to third order).

return
end subroutine leapfrog_3d_complex

!================================================================================

subroutine leapfrog_2d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff, raw_filter_coeff)

complex, intent(inout), dimension(:,:,:) :: a
complex, intent(in),    dimension(:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_filter_coeff

complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d
complex, dimension(size(a,1),size(a,2),1)           :: dt_a_3d

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
dt_a_3d(:,:,1) = dt_a
call leapfrog_3d_complex(a_3d, dt_a_3d, previous, current, future, delta_t, robert_coeff, raw_filter_coeff)
a = a_3d(:,:,1,:)

return
end subroutine leapfrog_2d_complex

!================================================================================

end module leapfrog_mod
