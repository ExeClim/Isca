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

module vert_coordinate_mod

!=======================================================================
!
!                          VERT_COORDINATE MODULE
!
!                Sets the coefficients a and b that define 
!                    the generalized vertical coordinate
!
!=======================================================================

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                   write_version_number, stdlog, close_file, &
                   check_nml_error

use constants_mod, only: pi

implicit none
private

public :: compute_vert_coord

!=======================================================================
! subroutine compute_vert_coord (vert_coord_option, scale_heights, surf_res, exponent, p_press, p_sigma, a, b)
!
! output:
!
!    real, dimension(:) :: a, b
!
!          a and b should be dimensioned by the number of interfaces
!                  = 1 + number of levels
!
!          these constants are intended for use in a model in which
!             the pressure at the interfaces will then be given
!
!                     p = a*p_ref + b*p_surf
!
!         where p_ref  is a constant reference pressure 
!           and p_surf is the instantaneous surface pressure
!      
!
!=======================================================================

character(len=128), parameter :: version = &
'$Id: vert_coordinate.F90,v 19.0 2012/01/06 20:35:42 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: siena_201211 $'

character(len=8) :: ch_tmp
logical :: entry_to_logfile_done = .false.

integer, parameter :: max_levels=100
real, dimension(max_levels+1) :: pk, bk

namelist /vert_coordinate_nml/ pk, bk

contains

!=======================================================================

subroutine compute_vert_coord (vert_coord_option, scale_heights, surf_res, exponent, p_press, p_sigma, reference_press, a, b)

character(len=*), intent(in) :: vert_coord_option
real,    intent(in) :: scale_heights, surf_res, exponent, p_press, p_sigma, reference_press
real,    intent(out), dimension(:) :: a, b

real, dimension(size(a,1)) :: a_sigma, b_sigma, a_press, b_press, f
character(len=32) :: chtmp='size(a)=      size(b)=          '
character(len=16) :: chtmp2, chtmp3

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

if(size(a,1).ne.size(b,1)) then
  write(chtmp( 9:12),'(i4)') size(a,1)
  write(chtmp(23:26),'(i4)') size(b,1)
  call error_mesg('compute_vert_coord', 'size(a) does not equal size(b). '//chtmp, FATAL)
endif

if(trim(vert_coord_option) == 'uneven_sigma' .or. trim(vert_coord_option) == 'hybrid') then
  if(scale_heights == 0.) then
    call error_mesg('compute_vert_coord',' zero is an invalid value for scale_heights.', FATAL)
  endif
  if(exponent == 0.) then
    call error_mesg('compute_vert_coord',' zero is an invalid value for exponent.', FATAL)
  endif
  if(surf_res <= 0. .or. surf_res > 1.0 ) then
    write(chtmp2,'(1pe16.8)') surf_res
    call error_mesg('compute_vert_coord', &
     'the namelist parameter surf_res must be < 1.0, but surf_res='//chtmp2, FATAL)
  endif
endif

if(trim(vert_coord_option) == 'hybrid') then
  if(p_sigma < p_press) then 
    write(chtmp2,'(1pe16.8)') p_sigma
    write(chtmp3,'(1pe16.8)') p_press
    call error_mesg('compute_vert_coord', &
     ' p_sigma must be greater than p_press, but p_sigma='//chtmp2//'  p_press='//chtmp3, FATAL)
  endif
endif

!--------------------------------------------------------------------

if(trim(vert_coord_option) == 'input') then
  call read_namelist(a,b)
else if(trim(vert_coord_option) == 'even_sigma') then
  call compute_even_sigma(a, b)
else if(trim(vert_coord_option) == 'uneven_sigma') then
  call compute_uneven_sigma(a, b, scale_heights, surf_res, exponent, .true.)
else if(trim(vert_coord_option) == 'hybrid') then
  call compute_uneven_sigma(a_sigma, b_sigma, scale_heights, surf_res, exponent, .false.)
  call compute_uneven_sigma(b_press, a_press, scale_heights, surf_res, exponent, .false.)
  f = transition(b_sigma, p_sigma, p_press)
  a = a_sigma*f + a_press*(1.0 - f)
  b = b_sigma*f + b_press*(1.0 - f)
  a = reference_press*a
else if(trim(vert_coord_option) == 'mcm') then
  call compute_old_model_sigma(a, b)
else if(trim(vert_coord_option) == 'v197') then
  call compute_v197_sigma(a, b)
else
  call error_mesg('compute_vert_coord','"'//trim(vert_coord_option)//'" is not a valid value for vert_coord_option.', FATAL)
end if

return
end subroutine compute_vert_coord

!-------------------------------------------------------------------------

function transition(p, p_sigma, p_press) result(t)

real, intent(in), dimension(:) :: p
real, intent(in)               :: p_sigma, p_press
real, dimension(size(p,1))       :: t

real    :: x, xx
integer :: k

do k = 1, size(p,1)
  if(p(k) <= p_press) then
    t(k) = 0.0
  else if(p(k) >= p_sigma) then
    t(k) = 1.0
  else
    x  = p(k)    - p_press
    xx = p_sigma - p_press
    t(k) = (sin(0.5*pi*x/xx))**2
  end if
enddo

return
end function transition

!-----------------------------------------------------------------------

subroutine read_namelist (a, b)

real, intent (out), dimension(:) :: a, b
integer :: namelist_unit, ierr, io

pk = 0.
bk = 0.

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=vert_coordinate_nml, iostat=io)
    ierr = check_nml_error(io, 'vert_coordinate_nml')
#else
    namelist_unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(namelist_unit, nml=vert_coordinate_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'vert_coordinate_nml')
    enddo
20  call close_file (namelist_unit)
#endif

if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=vert_coordinate_nml)

if(pk(2)+bk(2) == 0.) then
  call error_mesg('read_namelist in vert_coordinate_mod','No levels specified in namelist vert_coordinate_nml &
                  & or namelist is missing ', FATAL)
endif

if(pk(size(a,1))+bk(size(b,1)) == 0.) then
  call error_mesg('read_namelist in vert_coordinate_mod','Not enough levels specified in namelist vert_coordinate_nml', FATAL)
endif

if(pk(size(a,1)+1)+bk(size(b,1)+1) /= 0.) then
  call error_mesg('read_namelist in vert_coordinate_mod','Too many levels specified in namelist vert_coordinate_nml', FATAL)
endif

a = pk(1:size(a,1))
b = bk(1:size(b,1))

return
end subroutine read_namelist
!-----------------------------------------------------------------------

subroutine compute_even_sigma (a, b)

real, intent (out), dimension(:) :: a, b
integer :: k, num_levels

num_levels = size(a,1)-1

a = 0.0
do k=1,num_levels
  b(k) = float(k - 1)/float(num_levels)
end do
b(num_levels+1) = 1.0

return
end subroutine compute_even_sigma

!-------------------------------------------------------------------------

subroutine compute_uneven_sigma (a, b, scale_heights, surf_res, exponent, zero_top)

real,    intent (out), dimension(:) :: a, b
real,    intent (in) :: scale_heights, surf_res, exponent
logical, intent (in) :: zero_top

integer :: k, num_levels
real :: z, zeta, s2

num_levels = size(a,1)-1

a = 0.0

s2 = 1.0 - surf_res

do k=1,num_levels
  zeta = (1.-(float(k-1)/float(num_levels)))
  z  = surf_res*zeta + s2*(zeta**exponent)
  b(k) = exp(-z*scale_heights)
end do

b(num_levels+1) = 1.0
if(zero_top) b(1) = 0.0

return
end subroutine compute_uneven_sigma

!-------------------------------------------------------------------------
subroutine compute_v197_sigma(a, b)
real, intent (out), dimension(:) :: a, b

integer :: num_levels

num_levels = size(b,1)-1
if(num_levels /= 18) then
  write(ch_tmp,'(i8)') num_levels
  call error_mesg('compute_v197_sigma','num_levels='//ch_tmp//' It must be 18', FATAL)
endif

b = (/0.0, .0089163, .0342936, .0740741, .1262002, .1886145, .2592592, & 
           .3360768, .4170096, .5000000, .5829904, .6639231, .7407407, & 
           .8113854, .8737997, .9259259, .9657064, .9910837, 1.0/)

a = 0.

return
end subroutine compute_v197_sigma
!-------------------------------------------------------------------------
subroutine compute_old_model_sigma(a, b)
real, intent (out), dimension(:) :: a, b

integer :: num_levels

num_levels = size(b,1)-1
if(num_levels /= 14) then
  write(ch_tmp,'(i8)') num_levels
  call error_mesg('compute_old_model_sigma','num_levels='//ch_tmp//' It must be 14', FATAL)
endif
b = (/0.0,.03,.0707,.1311,.2102,.3036,.4062,.5138,.6226,.7284,.8255,.9066,.9640,.9933,1.0/)
a = 0.

return
end subroutine compute_old_model_sigma
!-------------------------------------------------------------------------
end module vert_coordinate_mod
