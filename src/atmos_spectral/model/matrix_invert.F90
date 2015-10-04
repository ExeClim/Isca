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

module matrix_invert_mod

use       fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                         write_version_number

implicit none

public  :: invert
integer, private :: maxmag

character(len=128), parameter :: version = '$Id matrix_invert.f90 $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'
logical :: entry_to_logfile_done = .false.

contains

subroutine invert(matrix, det)

real, intent(inout), dimension(:,:) :: matrix
real, intent(out) :: det

real, dimension(2*size(matrix,1)) :: dd, h
real, dimension(2*size(matrix,1),size(matrix,1)) :: ac, temp
real :: min_det=1.0e-30
character(len=24) :: chtmp
integer :: n, i, j, L, m, k

! *******************************************************************
!
! **** matrix_invert (MATRIX INVERSION AND DETERMINANT)
!
! **** QUESTIONS:  TRIVENI N. UPADHYAY, AUSTIN,TEXAS  X2207,MS 2186
!
! **** PURPOSE:
!        THIS SUBROUTINE INVERTS n by n NONSINGULAR MATRIX AND
!        FINDS IT'S DETERMINANT.
!
! **** ARGUMENTS:
!        matrix : INPUT MATRIX OF DIMENSION n by n TO BE INVERTED
!        n      : DIMENSION OF MATRIX
!        det    : DETERMINANT OF MATRIX
!
! **** PROCEDURE AND LIMITATIONS :
!          THIS SUBROUTINE USES THE METHOD OF ELEMENTARY
!        TRANFORMATIONS TO FIND THE INVERSE OF A MATRIX.
!        THE INPUT MATRIX IS DESTROYED IN COMPUTATION AND THE INVERSE
!        MATRIX TAKES ITS PLACE. FOR NUMERICAL ACCURACY, ELEMENTARY
!        TRANSFORMATIONS ARE PERFORMED IN CONJUNCTION WITH THE
!        'PIVOTAL' ELEMENT METHOD.
!          IF THE INPUT MATRIX IS SINGULAR (DETERMINANT LESS
!        THAN min_det), AN ERROR MESSAGE IS PRINTED OUT AND THE
!        PROGRAM IS TERMINATED.

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

n = size(matrix,1)

!   INITIALIZE

det = 1.0
ac(1:n,:) = matrix(:,:)
do j=1,n
  ac(n+1:2*n,j) = 0.0
  ac(n+j,    j) = 1.0
end do

do k=1,n

! FIND  LARGEST ELEMENT IN THE ROW

  h(k:n) = ac(k,k:n)
  m = n-k+1
  L = max_mag(h(k:n),M)+k

! INTERCHANGE COLUMNS IF THE LARGEST ELEMENT IS NOT THE DIAGONAL ELEMENT.

  if (k-L < 0) then
    do i=k,2*n
      dd(i)   = ac(i,k)
      ac(i,k) = ac(i,L)
      ac(i,L) = dd(i)
    end do
    det = -det
  end if
 
! DIVIDE THE COLUMN BY THE LARGEST ELEMENT
 
  det = det*ac(k,k)
  if (abs(det) < min_det) then
    write(chtmp,'(1pe24.16)') det
    call error_mesg('invert','DETERMINANT OF MATRIX ='//chtmp// &
    & ' THE MAGNITUDE OF THE DETERMINANT IS LESS THAN THE MINIMUM ALLOWED. &
    &  THE INPUT MATRIX APPEARS TO BE SINGULAR.',FATAL)
  endif
  h(k:2*n) = ac(k:2*n,k)/ac(k,k)
  do j=1,n
    temp(k:2*n,j) = h(k:2*n)*ac(k,j)
  end do
  ac(k:2*n,:) = ac(k:2*n,:) - temp(k:2*n,:)
  ac(k:2*n,k) = h(k:2*n)
end do

matrix(1:n,:) = ac(n+1:2*n,:)

return
end subroutine invert

function max_mag(h,m) result(max)

integer, intent(in) :: m
real, intent(in) :: h(m)
integer :: max, i
real :: rmax

max = 0
rmax = abs(h(1))
do i=1,m
  if (abs(h(i)) > rmax ) then
    rmax = abs(h(i))
    max = i-1
  endif
end do
return
end function max_mag

end module matrix_invert_mod
