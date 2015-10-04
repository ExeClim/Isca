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

module grid_fourier_mod

use       fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                         write_version_number

use constants_mod, only: pi

use fft_mod,       only: fft_init, fft_grid_to_fourier, fft_fourier_to_grid


!    provides one-dimensional grid-to-fourier and fourier-to-grid 
!    transforms on the first dimension of two- or three-dimensional fields.
!
!----------------------------------------------------------------------------

implicit none
private

character(len=128), parameter :: version = '$Id: grid_fourier.F90,v 11.0 2004/09/28 19:30:50 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

public :: grid_fourier_init, grid_fourier_end, trans_grid_to_fourier, trans_fourier_to_grid
public :: get_lon_max, get_longitude_origin, get_deg_lon

interface trans_grid_to_fourier
   module procedure trans_grid_to_fourier_3d, trans_grid_to_fourier_2d
end interface

interface trans_fourier_to_grid
   module procedure trans_fourier_to_grid_3d, trans_fourier_to_grid_2d
end interface

interface verify_fourier_imag
   module procedure verify_fourier_imag_3d, verify_fourier_imag_2d
end interface

integer :: num_lon
real    :: longitude_origin_local
logical :: check_local
logical :: module_is_initialized = .false.

real, allocatable, dimension(:) :: deg_lon

contains

!------------------------------------------------------------------------------------------------------

subroutine grid_fourier_init(num_lon_in, fourier_inc, check, longitude_origin)

integer, intent(in) :: num_lon_in, fourier_inc
logical, intent(in), optional   :: check
real,    intent(in), optional   :: longitude_origin
integer :: i
character(len=8) :: chtmp
real :: total_degrees

if(module_is_initialized) return

call write_version_number(version, tagname)

if(mod(num_lon_in,2) .ne. 0) then
  write(chtmp,'(i8)') num_lon_in
  call error_mesg('grid_fourier_init','num_lon must be even, but num_lon='//chtmp, FATAL)
end if

if(num_lon_in .eq. 0) then
  call error_mesg('grid_fourier_init','num_lon cannot be zero', FATAL)
end if

num_lon = num_lon_in

if(present(check)) then
  check_local = check
else
  check_local = .false.
end if

if(present(longitude_origin)) then
  longitude_origin_local = longitude_origin
else
  longitude_origin_local = 0.0
end if

total_degrees = 360./float(fourier_inc)

if(allocated(deg_lon)) deallocate(deg_lon)
allocate(deg_lon(num_lon))
do i=1,num_lon
  deg_lon(i) = 180*longitude_origin_local/pi &
     + (i-1)*total_degrees/float(num_lon)
  if(deg_lon(i) .ge. total_degrees) then
    deg_lon(i) = deg_lon(i) - total_degrees
  endif
  if(deg_lon(i) .lt. 0.0) then
    deg_lon(i) = deg_lon(i) + total_degrees
  endif
end do

call fft_init(num_lon)

module_is_initialized = .true.

return
end subroutine grid_fourier_init

!------------------------------------------------------------------------------------------------------

function trans_grid_to_fourier_3d(grid) result(fourier)

real,    intent(in),  dimension(:,:,:) :: grid
complex, dimension(0:size(grid,1)/2, size(grid,2), size(grid,3)) :: fourier
real, dimension(size(grid,1)+1,size(grid,2),size(grid,3)) :: grid1
character(len=8) :: chtmp1, chtmp2

if(.not.module_is_initialized) then
  call error_mesg('trans_grid_to_fourier','module grid_fourier not initialized', FATAL)
end if

if(size(grid,1).ne.num_lon) then
  write(chtmp1,'(i8)') size(grid,1)
  write(chtmp2,'(i8)') num_lon
  call error_mesg('trans_grid_to_fourier', &
    'size(grid,1) must equal num_lon, but size(grid,1)='//chtmp1//' num_lon='//chtmp2, FATAL)
end if

grid1(1:size(grid,1),:,:) = grid
grid1(size(grid,1)+1,:,:) = grid(1,:,:)
fourier = fft_grid_to_fourier(grid1)

return
end function trans_grid_to_fourier_3d
!------------------------------------------------------------------------------------------------------

function trans_fourier_to_grid_3d(fourier) result(grid)

complex, intent(in), dimension(0:,:,:) :: fourier
real, dimension(2*(size(fourier,1)-1), size(fourier,2), size(fourier,3)) :: grid
real, dimension(2*size(fourier,1)-1,   size(fourier,2), size(fourier,3)) :: grid1
character(len=8) :: chtmp1, chtmp2

if(.not.module_is_initialized) then
  call error_mesg('trans_fourier_to_grid','module grid_fourier not initialized', FATAL)
end if

if(2*(size(fourier,1)-1) .ne. num_lon) then
  write(chtmp1,'(i8)') size(fourier,1)
  write(chtmp2,'(i8)') num_lon 
  call error_mesg('trans_fourier_to_grid', &
  '2*(size(fourier,1)-1) must equal num_lon, but size(fourier,1)='//chtmp1//' num_lon='//chtmp2, FATAL)
end if

if(check_local)  call verify_fourier_imag(fourier)

grid1 = fft_fourier_to_grid(fourier)
grid = grid1(1:2*(size(fourier,1)-1),:,:)

return
end function trans_fourier_to_grid_3d
!------------------------------------------------------------------------------------------------------

subroutine verify_fourier_imag_3d(fourier)

complex, intent(in), dimension(0:,:,:) :: fourier
integer :: j, k
character(len=56) :: chtmp='aimag(fourier(   ,   ,   ))=                            '

do k=1,size(fourier,3)
  do j=1,size(fourier,2)
    if(aimag(fourier(0,j,k)) .ne. 0.0) then
      chtmp(15:17) = '  0'
      write(chtmp(19:21),'(i3)') j
      write(chtmp(23:25),'(i3)') k
      write(chtmp(29:52),'(1pe24.16)') aimag(fourier(0,j,k))
      call error_mesg('verify_fourier_imag',chtmp//'  It should be zero.', FATAL)
    end if
  end do
end do

if(mod(num_lon,2).eq.0) then
  do k=1,size(fourier,3)
    do j=1,size(fourier,2)
      if(aimag(fourier(num_lon/2,j,k)) .ne. 0.0) then
        write(chtmp(15:17),'(i3)') num_lon/2
        write(chtmp(19:21),'(i3)') j
        write(chtmp(23:25),'(i3)') k
        write(chtmp(29:52),'(1pe24.16)') aimag(fourier(0,j,k))
        call error_mesg('verify_fourier_imag',chtmp//'  It should be zero.', FATAL)
      end if
    end do
  end do
end if

return
end subroutine verify_fourier_imag_3d
!------------------------------------------------------------------------------------------------------

!  The remaining routines are 2d version of the 3d versions described above)

!------------------------------------------------------------------------------------------------------

function trans_fourier_to_grid_2d (fourier) result(grid)

complex, intent(in), dimension(0:,:) :: fourier
complex, dimension(0:size(fourier,1)-1, size(fourier,2), 1) :: fourier_3d

real, dimension(2*(size(fourier,1)-1), size(fourier,2)) :: grid
real, dimension(2*(size(fourier,1)-1), size(fourier,2), 1) :: grid_3d

fourier_3d(:,:,1) = fourier(:,:)
grid_3d = trans_fourier_to_grid_3d(fourier_3d)
grid(:,:) = grid_3d(:,:,1)

return
end function trans_fourier_to_grid_2d
!------------------------------------------------------------------------------------------------------

function trans_grid_to_fourier_2d (grid) result(fourier)

real,    intent(in),  dimension(:,:) :: grid
complex, dimension(0:size(grid,1)/2, size(grid,2)) :: fourier

real, dimension(size(grid,1), size(grid,2), 1) :: grid_3d
complex, dimension(0:size(grid,1)/2, size(grid,2), 1) :: fourier_3d

grid_3d(:,:,1) = grid(:,:)
fourier_3d = trans_grid_to_fourier_3d(grid_3d)
fourier(:,:) = fourier_3d(:,:,1)

return
end function trans_grid_to_fourier_2d

!------------------------------------------------------------------------------------------------------

subroutine verify_fourier_imag_2d(fourier)

complex, intent(in), dimension(0:,:) :: fourier
complex, dimension(0:size(fourier,1)-1, size(fourier,2), 1) :: fourier_3d

fourier_3d(:,:,1) = fourier(:,:)
call verify_fourier_imag_3d(fourier_3d)

return
end subroutine verify_fourier_imag_2d

!------------------------------------------------------------------------------------------------------

subroutine get_lon_max(lon_max_out)

integer, intent (out) :: lon_max_out

if(.not.module_is_initialized) then
  call error_mesg('get_lon_max','module grid_fourier not initialized', FATAL)
end if

lon_max_out = num_lon

return
end subroutine get_lon_max

!------------------------------------------------------------------------------------------------------

subroutine get_longitude_origin(longitude_origin_out)

real, intent (out) :: longitude_origin_out

if(.not.module_is_initialized) then
  call error_mesg('get_longitude_origin','module grid_fourier not initialized', FATAL)
end if

longitude_origin_out = longitude_origin_local

return
end subroutine get_longitude_origin

!------------------------------------------------------------------------------------------------------

subroutine get_deg_lon(deg_lon_out)

real, intent (out), dimension(:) :: deg_lon_out
character(len=8) :: chtmp1, chtmp2

if(.not.module_is_initialized) then
  call error_mesg('get_deg_lon','module grid_fourier not initialized', FATAL)
end if
if(size(deg_lon_out,1).ne.num_lon) then
  write(chtmp1,'(i8)') size(deg_lon_out,1)
  write(chtmp2,'(i8)') num_lon
  call error_mesg('get_deg_lon', &
    'size of deg_lon does not equal num_lon. size(deg_lon)='//chtmp1//' num_lon='//chtmp2, FATAL)
end if

deg_lon_out(:) = deg_lon(:)

return
end subroutine get_deg_lon

!------------------------------------------------------------------------------------------------------

subroutine grid_fourier_end

if(.not.module_is_initialized) return
deallocate(deg_lon)
!call fft_end
module_is_initialized = .false.

return
end subroutine grid_fourier_end

!------------------------------------------------------------------------------------------------------

end module grid_fourier_mod
