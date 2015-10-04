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

module topog_regularization_mod

!  produces regularized topography according to Lindberg and Broccoli,
!  J. of Climate vol 9, no 11 pg. 2641-2659 (1996)

!  Originally coded by Charles Jackson
!  Modified for FMS by Peter Phillipps

use             fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                               write_version_number

use             mpp_mod, only: mpp_chksum

use       constants_mod, only: pi

use      transforms_mod, only: compute_gaussian, compute_legendre, &
                               trans_grid_to_spherical, trans_spherical_to_grid,          &
                               get_sin_lat, get_wts_lat, transforms_are_initialized,      &
                               get_lon_max, get_lat_max, get_num_fourier, get_fourier_inc,&
                               get_num_spherical, get_grid_domain, get_spec_domain,       &
                               grid_domain, spectral_domain, area_weighted_global_mean

use     mpp_domains_mod, only: mpp_global_field

implicit none
private

character(len=128), parameter :: version = &
'$Id: topog_regularization.F90,v 13.0 2006/03/28 21:17:37 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: siena_201211 $'

public :: compute_lambda, regularize

integer, parameter :: itmax=1000
real,    parameter :: tolerance = 1.e-5

integer :: is, ie, js, je, ms, me, ns, ne
integer :: lon_max, lat_max, num_fourier, num_spherical, fourier_inc, nmax

real,    allocatable, dimension(:,:  ) :: smoothed_field_tmp, rough, cost_field
real,    allocatable, dimension(:    ) :: wts_lat_global, sin_lat_global, facm, sin_facm
complex, allocatable, dimension(:,:  ) :: Dnm, anm, bnm, Hnm, DR2, DelAnm, DelBnm

logical :: module_is_initialized=.false.

character(len=8) :: chtmp1, chtmp2

contains

!============================================================================================

subroutine compute_lambda(ocean_topog_smoothing, ocean_mask, unsmoothed_field, lambda, actual_fraction_smoothed)

real,    intent(in) :: ocean_topog_smoothing
logical, intent(in), dimension(:,:) :: ocean_mask
real,    intent(in), dimension(:,:) :: unsmoothed_field
real,    intent(out) :: lambda, actual_fraction_smoothed

real :: lambda_1, lambda_2, fraction_smoothed_1, fraction_smoothed_2
real :: tol_lambda = .001
integer :: it_lambda, itmax_lambda=20

if(.not.module_is_initialized) then
  call topog_regularization_init(ocean_mask)
endif

if(any(shape(unsmoothed_field) /= (/ie-is+1,je-js+1/))) then
  write(chtmp1,'(2i4)') shape(unsmoothed_field)
  write(chtmp2,'(2i4)')  (/ie-is+1,je-js+1/)
  call error_mesg('compute_lambda', &
    'Input argument unsmoothed_field has incorrect dimensions. shape(unsmoothed_field)='//chtmp1//'  Should be '//chtmp2,FATAL)
endif

lambda_1 = 1.e-7
lambda_2 = 2.e-7
call regularize(lambda_1, ocean_mask, unsmoothed_field, smoothed_field_tmp, fraction_smoothed_1)
if(abs(ocean_topog_smoothing-fraction_smoothed_1) < tol_lambda) then
  lambda = lambda_1
  actual_fraction_smoothed = fraction_smoothed_1
  goto 20
endif
call regularize(lambda_2, ocean_mask, unsmoothed_field, smoothed_field_tmp, fraction_smoothed_2)
if(abs(ocean_topog_smoothing-fraction_smoothed_2) < tol_lambda) then
  lambda = lambda_2
  actual_fraction_smoothed = fraction_smoothed_2
  goto 20
endif
if(fraction_smoothed_1 > ocean_topog_smoothing .or. fraction_smoothed_2 > ocean_topog_smoothing) then
  call error_mesg('compute_lambda', &
     'Iterative scheme for computing lambda may not work unless initial values of lambda_1 and lambda_2 are reduced.', FATAL)
endif
lambda_1 = ((fraction_smoothed_2-ocean_topog_smoothing)*lambda_1 + &
            (ocean_topog_smoothing-fraction_smoothed_1)*lambda_2)/(fraction_smoothed_2-fraction_smoothed_1)
if(lambda_1 < 0.) then
  call error_mesg('compute_lambda', &
    'Iterative scheme for finding lambda will not work unless initial values of lambda_1 and lambda_2 are reduced.', FATAL)
endif
call regularize(lambda_1, ocean_mask, unsmoothed_field, smoothed_field_tmp, fraction_smoothed_1)
do it_lambda=1,itmax_lambda
  if(abs(ocean_topog_smoothing-fraction_smoothed_1) < tol_lambda) then
    lambda = lambda_1
    actual_fraction_smoothed = fraction_smoothed_1
    goto 20
  endif
  lambda_2 = ((fraction_smoothed_2-ocean_topog_smoothing)*lambda_1 + &
              (ocean_topog_smoothing-fraction_smoothed_1)*lambda_2)/(fraction_smoothed_2-fraction_smoothed_1)
  if(lambda_2 < 0.) then
    write(chtmp1,'(i8)') it_lambda
    call error_mesg('compute_lambda', &
        'Iterative scheme for finding lambda failed. lambda went negative on iteration number'//chtmp1, FATAL)
  endif
  call regularize(lambda_2, ocean_mask, unsmoothed_field, smoothed_field_tmp, fraction_smoothed_2)
  if(abs(ocean_topog_smoothing-fraction_smoothed_2) < tol_lambda) then
    lambda = lambda_2
    actual_fraction_smoothed = fraction_smoothed_2
    goto 20
  endif
  lambda_1 = ((fraction_smoothed_2-ocean_topog_smoothing)*lambda_1 + &
              (ocean_topog_smoothing-fraction_smoothed_1)*lambda_2)/(fraction_smoothed_2-fraction_smoothed_1)

  call regularize(lambda_1, ocean_mask, unsmoothed_field, smoothed_field_tmp, fraction_smoothed_1)
enddo
call error_mesg('compute_lambda','Cannot converge on a value of lambda. Perhaps more interations are needed.', FATAL)
20 continue

return
end subroutine compute_lambda
!============================================================================================

subroutine regularize(lambda, ocean_mask, unsmoothed_field, smoothed_field, fraction_smoothed)

real,    intent(in) :: lambda
logical, intent(in),  dimension(:,:) :: ocean_mask
real,    intent(in),  dimension(:,:) :: unsmoothed_field
real,    intent(out), dimension(size(ocean_mask,1), size(ocean_mask,2)) :: smoothed_field
real,    intent(out) :: fraction_smoothed

real :: converg, cost, oldcost, lamcost, lamcosti
integer :: m, n, it

if(.not.module_is_initialized) then
  call topog_regularization_init(ocean_mask)
endif

if(any(shape(unsmoothed_field) /= (/ie-is+1,je-js+1/))) then
  write(chtmp1,'(2i4)') shape(unsmoothed_field)
  write(chtmp2,'(2i4)')  (/ie-is+1,je-js+1/)
  call error_mesg('regularize', &
   'Input argument unsmoothed_field has incorrect dimensions. shape(unsmoothed_field)='//chtmp1//'  Should be '//chtmp2,FATAL)
endif

if(is /= 1 .or. ie /= lon_max) then
  call error_mesg('regularize', &
    'subroutine regularize is not yet coded for 2-d decomposition. It was assumed that it will never be needed.',FATAL)
endif

Hnm=cmplx(0.,0.)
do n=ns,nmax
  do m=ms,me
    Hnm(m,n)=1./(1.+lambda*Dnm(m,n)*((n+m)*(n+m+1))**2)
  enddo
enddo

bnm=cmplx(0.,0.)
call trans_grid_to_spherical(unsmoothed_field,bnm)

!    anm (equation 6.3)

anm = cmplx(0.,0.)
do n=ns,nmax
  do m=ms,me
    anm(m,n)=bnm(m,n)/(1.+lambda*((n+m)*(n+m+1))**2)
  enddo
enddo

DelAnm = cmplx(0.,0.)
do n=ns,nmax
  do m=ms,me
    DelAnm(m,n)=(n+m)*(n+m+1)*anm(m,n)
  enddo
enddo

call trans_spherical_to_grid(DelAnm,rough)

converg=1.
cost=0.

DelBnm = cmplx(0.,0.)
do it=1,itmax
  if (abs(converg) < tolerance) goto 60

! rough is zeroed out over land
  where(.not.ocean_mask)
    rough = 0.
  end where

! Do an iteration of smoothing
  call trans_grid_to_spherical(rough,DR2)
  do n=ns,nmax
    do m=ms,me
      DR2(m,n)=(n+m)*(n+m+1)*DR2(m,n)
    enddo
  enddo

  do n=ns,nmax
    do m=max(ms,1),me
      anm(m,n)=(anm(m,n)+Hnm(m,n)*(bnm(m,n)-anm(m,n))-lambda*Hnm(m,n)*DR2(m,n))*sin_facm(m)/facm(m)
    enddo
  enddo
  if(ms == 0) then
    do n=ns,nmax
      anm(0,n)=anm(0,n)+Hnm(0,n)*(bnm(0,n)-anm(0,n))-lambda*Hnm(0,n)*DR2(0,n)
    enddo
  endif

! transform the smoothed field to grid
  call trans_spherical_to_grid(anm,smoothed_field)

! compute cost function (eq. 6.4)
  do n=ns,nmax
    do m=ms,me
      DelAnm(m,n)=(n+m)*(n+m+1)*anm(m,n)
    enddo
  enddo
  call trans_spherical_to_grid(DelAnm,rough)
  where(ocean_mask)
    cost_field = (unsmoothed_field-smoothed_field)**2 + lambda*rough**2
  else where
    cost_field = 0.
  end where
  oldcost=cost
  cost = area_weighted_global_mean(cost_field)

  if (it > 1) then
    converg=(oldcost-cost)/oldcost
  endif
enddo
call error_mesg('regularize','Failure to converge',FATAL)
60 continue

do n=ns,nmax
  do m=ms,me
    DelBnm(m,n)=(n+m)*(n+m+1)*bnm(m,n)
  enddo
enddo
call trans_spherical_to_grid(DelBnm,rough)
where(ocean_mask)
  cost_field = rough**2
else where
  cost_field = 0.
end where
lamcosti = area_weighted_global_mean(cost_field)
call trans_spherical_to_grid(DelAnm,rough)
where(ocean_mask)
  cost_field = rough**2
else where
  cost_field = 0.
end where
lamcost = area_weighted_global_mean(cost_field)
fraction_smoothed = 1. - lamcost/lamcosti

if(mpp_pe() == mpp_root_pe()) then
  print '("Message from subroutine regularize: lambda=",1pe16.8,"  fraction_smoothed=",1pe16.8)',lambda,fraction_smoothed
endif

return
end subroutine regularize
!===================================================================================
subroutine topog_regularization_init(ocean_mask)

logical, intent(in), dimension(:,:) :: ocean_mask
integer :: m, i, j
real,    allocatable, dimension(:,:,:) :: legendre_global, legendre
logical, allocatable, dimension(:,:) :: ocean_mask_global

call write_version_number(version, tagname)

if(.not.transforms_are_initialized()) then
  call error_mesg('topog_regularization_init','Transforms are not initialized',FATAL)
endif

call get_grid_domain(is, ie, js, je)
call get_spec_domain(ms, me, ns, ne)

if(any(shape(ocean_mask) /= (/ie-is+1,je-js+1/))) then
  write(chtmp1,'(2i4)') shape(ocean_mask)
  write(chtmp2,'(2i4)')  (/ie-is+1,je-js+1/)
  call error_mesg('topog_regularization_init', &
     'Input argument ocean_mask has incorrect dimensions. shape(ocean_mask)='//chtmp1//'  Should be '//chtmp2,FATAL)
endif

call get_lon_max(lon_max)
call get_lat_max(lat_max)
call get_num_fourier(num_fourier)
call get_num_spherical(num_spherical)
call get_fourier_inc(fourier_inc)

nmax = min(num_fourier,ne)

allocate(sin_lat_global(lat_max))
allocate(wts_lat_global(lat_max))
call get_sin_lat(sin_lat_global)
call get_wts_lat(wts_lat_global)

allocate(facm    (ms:me))
allocate(sin_facm(ms:me))
do m=ms,me
  facm(m) = pi*m/(2.*num_fourier)
  sin_facm(m) = sin(facm(m))
enddo

allocate(Dnm(ms:me,ns:ne),    Hnm(ms:me,ns:ne),    bnm(ms:me,ns:ne), anm(ms:me,ns:ne))
allocate(DR2(ms:me,ns:ne), DelAnm(ms:me,ns:ne), DelBnm(ms:me,ns:ne))
allocate(smoothed_field_tmp(is:ie,js:je), rough(is:ie,js:je), cost_field(is:ie,js:je))

allocate(legendre_global(0:num_fourier,0:num_spherical,lat_max/2))
allocate(legendre(ms:me, ns:ne, lat_max/2))
call compute_legendre(legendre_global, num_fourier, fourier_inc, num_spherical, -sin_lat_global(1:lat_max/2), lat_max/2)
legendre = legendre_global(ms:me,ns:ne,:)

allocate(ocean_mask_global(lon_max,lat_max))
call mpp_global_field(grid_domain, ocean_mask, ocean_mask_global)
Dnm=cmplx(0.,0.)
do i=1,lon_max
  do j=1,lat_max/2
    if(ocean_mask_global(i,j)) then
      Dnm = Dnm + wts_lat_global(j)*legendre(:,:,j)**2
    endif
  enddo
  do j=lat_max/2+1,lat_max
    if(ocean_mask_global(i,j)) then
      Dnm = Dnm + wts_lat_global(j)*legendre(:,:,lat_max+1-j)**2
    endif
  enddo
enddo
deallocate(legendre_global, legendre, ocean_mask_global)
Dnm = Dnm/lon_max

module_is_initialized = .true.
return

end subroutine topog_regularization_init
!===================================================================================

end module topog_regularization_mod
