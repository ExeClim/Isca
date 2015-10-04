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

module gauss_and_legendre_mod

use       fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                         write_version_number

use constants_mod, only: pi

!-----------------------------------------------------------------------
!   computes Gaussian latitudes and associated Legendre polynomials
!
!-----------------------------------------------------------------------

implicit none
private

character(len=128), parameter :: version = '$Id: gauss_and_legendre.F90,v 10.0 2003/10/24 22:01:02 fms Exp $'
character(len=128), parameter :: tagname ='$Name: siena_201211 $'

logical :: entry_to_logfile_done=.false.

public :: compute_legendre, compute_gaussian

contains

!-----------------------------------------------------------------------
subroutine compute_legendre(legendre, num_fourier, fourier_inc,  &
                num_spherical, sin_lat, n_lat)
!-----------------------------------------------------------------------

integer, intent (in) :: num_fourier, fourier_inc, num_spherical, n_lat
real,    intent (in), dimension(n_lat) :: sin_lat

real, intent (out), dimension(0:num_fourier, 0:num_spherical, n_lat)  :: legendre

integer :: j, m, n , fourier_max
real, dimension(0:num_fourier*fourier_inc, 0:num_spherical) :: poly, eps, m2, l2
real, dimension(0:num_fourier*fourier_inc) :: b
real, dimension(n_lat) :: cos_lat

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done=.true.
endif

fourier_max = num_fourier*fourier_inc

do j=1,n_lat
  cos_lat(j)=sqrt(1-sin_lat(j)*sin_lat(j))
end do

do n=0,num_spherical
  do m=0,fourier_max
    m2(m,n)   = m**2
    l2(m,n)   = (m+n)**2
  end do
end do

eps= sqrt((l2 -m2)/(4.0*l2 - 1.0))

do m=1,fourier_max
  b(m) = SQRT(0.5*(2.0*float(m) + 1.0)/float(m))
end do

poly(0,0) = sqrt(0.5)
do j=1,n_lat     
  do m=1,fourier_max
    poly(m,0) = b(m)*cos_lat(j)*poly(m-1,0)
  end do
  do m=0,fourier_max
     poly(m,1) =  sin_lat(j)*poly(m,0)/eps(m,1)
  end do

  do n=2,num_spherical
    do m=0,fourier_max
      poly(m,n) = (sin_lat(j)*poly(m,n-1) - eps(m,n-1)*poly(m,n-2))/eps(m,n)
    end do
  end do

  do n = 0, num_spherical
    do m = 0,num_fourier
      legendre(m,n,j) = poly(m*fourier_inc,n)
    end do
  end do
end do

return
end subroutine compute_legendre

!----------------------------------------------------------------------
subroutine compute_gaussian(sin_hem, wts_hem, n_hem)
!----------------------------------------------------------------------
!
!     reference:
!       press, h. william, et. al., numerical recipes (fortran version),
!       cambridge, england: cambridge university press (1990)
! 
!------------------------------------------------------------------------

integer, intent (in) :: n_hem
real, intent (out), dimension(n_hem) :: sin_hem, wts_hem

real :: converg
integer :: itermax
integer :: i, iter, j, n, nprec
real :: pp, p1, p2, p3, z, z1

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done=.true.
endif

! must use a more relaxed convergence criteria on the
! workstations than that for the cray T90
! fez code is commented out

!if(kind(converg).eq.8) then
!  converg = 1.0E-15
!else if(kind(converg).eq.4) then
!  converg = 1.0E-7
!else
!  call error_mesg('compute_gaussian','dont know what value to use for converg', FATAL)
!end if

! The 2 lines of code below will yeild a different result than the fez code
! when kind(converg)=4. converg is 1.0E-6 instead of 1.0E-7
! This should be investigated further, but it's OK for now because it yeilds
! the same result on the HPCS. -- pjp
nprec = precision(converg)
converg = .1**nprec


itermax = 10

n=2*n_hem
do i=1,n_hem
  z = cos(pi*(i - 0.25)/(n + 0.5))
  do iter=1,itermax
     p1 = 1.0
     p2 = 0.0

     do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
     end do

     pp = n*(z*p1 - p2)/(z*z - 1.0E+00)
     z1 = z
     z  = z1 - p1/pp
     if(ABS(z - z1) .LT. converg) go to 10
  end do
  call error_mesg('compute_gaussian','abscissas failed to converge in itermax iterations', FATAL)

  10  continue

  sin_hem (i)     = z
  wts_hem (i)     = 2.0/((1.0 - z*z)*pp*pp)

end do

return
end subroutine compute_gaussian
!----------------------------------------------------------------------

end module gauss_and_legendre_mod
