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

module spherical_fourier_mod

!Balaji: parallel spectral model using transpose method

use fms_mod,         only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number

use mpp_mod,         only: mpp_error

use mpp_domains_mod, only: domain1D, mpp_get_compute_domains, mpp_get_domain_components, mpp_get_layout

use   constants_mod, only: pi

use spec_mpp_mod,    only: get_grid_domain, grid_domain, get_spec_domain

use spherical_mod,   only: spherical_init, spherical_end,                &
                           compute_lon_deriv_cos, compute_lat_deriv_cos, &
                           compute_laplacian, compute_vor, compute_div,  &
                           get_spherical_wave, get_fourier_wave,         &
                           get_eigen_laplacian, compute_gradient_cos,    &
                           compute_ucos_vcos, compute_vor_div,           &
                           triangular_truncation, rhomboidal_truncation

use gauss_and_legendre_mod, only: compute_legendre, compute_gaussian

!-------------------------------------------------------------------------
!  provides latitudinal transforms from spherical harmonics to 
!       (zonal fourier, latitude (gaussian) grid) space , and 
!       related operations
!       works within the window provided by windows module
!
!       spherical fields are complex, dimension(0:num_fourier, 0:num_spherical)
!         where the zonal wavenumber of (m,n) is M = m*fourier_inc
!             the "meridional" wavenumber is n
!             the "total", 2D, spherical wavenumber L = M+n
!       "fourier fields" are complex, dimension(0:num_fourier, lat_max)
!
!------------------------------------------------------------------------

implicit none
private

character(len=128), parameter :: version = '$Id: spherical_fourier.F90,v 13.0 2006/03/28 21:18:33 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

interface trans_spherical_to_fourier
   module procedure trans_spherical_to_fourier_3d,  &
                    trans_spherical_to_fourier_2d
end interface
interface trans_fourier_to_spherical
   module procedure trans_fourier_to_spherical_3d,  &
                    trans_fourier_to_spherical_2d
end interface

! ---------------------------------------------------

!  public interfaces defined in this module

public :: spherical_fourier_init, spherical_fourier_end
public :: trans_spherical_to_fourier, trans_fourier_to_spherical
public :: get_south_to_north
public :: get_sin_lat, get_cos_lat, get_cosm_lat, get_cosm2_lat
public :: get_deg_lat, get_wts_lat

! ---------------------------------------------------

!  public interfaces carried forward from used modules

public :: spherical_init, spherical_end
public :: compute_lon_deriv_cos, compute_lat_deriv_cos
public :: compute_laplacian, compute_vor, compute_div
public :: get_spherical_wave, get_fourier_wave
public :: get_eigen_laplacian, compute_gradient_cos
public :: compute_ucos_vcos, compute_vor_div
public :: triangular_truncation, rhomboidal_truncation


public :: compute_legendre, compute_gaussian

! ---------------------------------------------------
integer :: fourier_max
integer :: fourier_inc
logical :: south_to_north_local
logical :: make_symmetric_local

integer :: lat_max
integer :: num_fourier
integer :: num_spherical

real, allocatable, dimension(:) :: sin_lat
real, allocatable, dimension(:) :: cos_lat
real, allocatable, dimension(:) :: cosm_lat
real, allocatable, dimension(:) :: cosm2_lat
real, allocatable, dimension(:) :: deg_lat
real, allocatable, dimension(:) :: wts_lat
real, allocatable, dimension(:) :: sin_hem
real, allocatable, dimension(:,:,:) :: legendre
real, allocatable, dimension(:,:,:) :: legendre_wts

type(domain1D),save :: grid_domain_y !used by the trans_s_f anfd f_s routines
integer :: grid_layout(2)
logical, private :: debug=.FALSE.

logical :: module_is_initialized = .false.
integer :: is, ie, js, je, ms, me, ns, ne

contains

!-----------------------------------------------------------------------
subroutine spherical_fourier_init(radius, lat_max_in, num_fourier_in, &
                   fourier_inc_in, num_spherical_in, south_to_north, make_symmetric)
!-----------------------------------------------------------------------

real,    intent(in) :: radius
integer, intent(in) :: lat_max_in
integer, intent(in) :: num_fourier_in
integer, intent(in) :: fourier_inc_in
integer, intent(in) :: num_spherical_in
logical, intent(in), optional :: south_to_north, make_symmetric

call write_version_number(version, tagname)

if(present(south_to_north)) then
  south_to_north_local = south_to_north
else
  south_to_north_local = .true.
end if

if(present(make_symmetric)) then
  make_symmetric_local = make_symmetric
else
  make_symmetric_local = .false.
end if

call get_grid_domain(is, ie, js, je)
call get_spec_domain(ms, me, ns, ne)

lat_max       = lat_max_in

fourier_inc   = fourier_inc_in
num_spherical = num_spherical_in
num_fourier   = num_fourier_in
fourier_max   = num_fourier*fourier_inc

call spherical_init(radius, num_fourier, fourier_inc, num_spherical, make_symmetric=make_symmetric_local)
call define_gaussian
call define_legendre

call mpp_get_domain_components( grid_domain, y=grid_domain_y )
call mpp_get_layout( grid_domain, grid_layout )
module_is_initialized = .true.

return
end subroutine spherical_fourier_init

!-------------------------------------------------------------------------
subroutine trans_spherical_to_fourier_3d(spherical,fourier)
!-------------------------------------------------------------------------

!in parallel, assumes that computational domains are being passed in.

  complex, intent(in),  dimension(ms:,ns:,:) :: spherical
  complex, intent(out), dimension(ms:,:,:,:) :: fourier

  complex :: x_odd, x_even

  integer :: m,n,j,k, num_levels

  integer :: neven, nodd
  integer :: jd, joff, nd, nj
  integer, dimension(grid_layout(2)) :: jstart

  if( size(spherical,1).ne.me-ms+1 )&
       call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_FOURIER: size(spherical,1).NE.me-ms+1.' )
  if( size(spherical,2).NE.ne-ns+1 )&
       call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_FOURIER: size(spherical,2).NE.ne-ns+1.' )
  if( size(fourier,2)*size(fourier,4).NE.lat_max ) &
       call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_FOURIER: size(fourier,2)*size(fourier,4).NE.lat_max.' )

  num_levels = size(spherical,3)
  nd = grid_layout(2)   !also equal to size(fourier,4), could put in an error check
  nj = je-js+1 !this version requires grid domains equal in Y
  if( debug )write( 0,'(a,i2,2i4)' )'pe, nd, nj=', mpp_pe(), nd, nj

  if( mod(ns,2).EQ.0 )then
      neven = ns
      nodd = ns+1
  else
      nodd = ns
      neven = ns+1
  endif

  call mpp_get_compute_domains( grid_domain_y, begin=jstart(:) )
  DOMAIN_LOOP: do jd=1,nd/2+1
!each domain jd has a mirror domain nd+1-jd
!if nd is even we stop at jd=nd/2; if nd is odd we stop at j=nj/2
     if( jd.EQ.nd/2+1 .AND. mod(nd,2).EQ.0 )exit DOMAIN_LOOP
     joff = jstart(jd) - 1
     if( debug )write( 0,'(a,i2,2i4)' )'pe, jd, joff=', mpp_pe(), jd, joff
     do j = 1,nj
        if( jd.EQ.nd/2+1 .AND. j.GT.nj/2 )exit DOMAIN_LOOP
        if(south_to_north_local) then !how long are we going to keep carrying around S2N and N2S?
            do k=1,num_levels
               do m=ms,me
                  x_even = cmplx(0.,0.)
                  x_odd  = cmplx(0.,0.)

                  do n=neven,ne,2
                     x_even = x_even + spherical(m,n,k) * legendre(m,n,j+joff)
                  end do
                  do n=nodd,ne,2
                     x_odd  = x_odd  + spherical(m,n,k) * legendre(m,n,j+joff)
                  end do

                  fourier(m,     j,k,     jd) = x_even - x_odd
                  fourier(m,nj+1-j,k,nd+1-jd) = x_even + x_odd
               end do
            end do
        else
            do k=1,num_levels
               do m=ms,me
                  x_even = cmplx(0.,0.)
                  x_odd  = cmplx(0.,0.)

                  do n=neven,ne,2
                     x_even = x_even + spherical(m,n,k) * legendre(m,n,j+joff)
                  end do
                  do n=nodd,ne,2
                     x_odd  = x_odd  + spherical(m,n,k) * legendre(m,n,j+joff)
                  end do

                  fourier(m,     j,k,     jd) = x_even + x_odd
                  fourier(m,nj+1-j,k,nd+1-jd) = x_even - x_odd
               end do
            end do
        end if
     end do
  end do DOMAIN_LOOP

  return
end subroutine trans_spherical_to_fourier_3d

!-----------------------------------------------------------------------
subroutine trans_fourier_to_spherical_3d(fourier,spherical)
!-----------------------------------------------------------------------

  complex, intent (out), dimension(ms:,ns:,:) :: spherical
  complex, intent (in),  dimension(ms:,:,:,:) :: fourier

  integer :: k,j,m,n
!  complex, dimension(ms:me, size(spherical,3)) :: x_odd, x_even
  complex :: x_odd, x_even

  integer :: num_levels, nd, nj

  integer :: neven, nodd
  integer :: jd, joff
  integer, dimension(grid_layout(2)) :: jstart

  if( size(spherical,1).ne.me-ms+1 )&
       call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_FOURIER: size(spherical,1).NE.me-ms+1.' )
  if( size(spherical,2).NE.ne-ns+1 )&
       call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_FOURIER: size(spherical,2).NE.ne-ns+1.' )
  if( size(fourier,2)*size(fourier,4).NE.lat_max ) &
       call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_FOURIER: size(fourier,2)*size(fourier,4).NE.lat_max.' )

  num_levels = size(spherical,3)
  nd = grid_layout(2)   !also equal to size(fourier,4), could put in an error check
  nj = je-js+1 !this version requires grid domains equal in Y

  if( mod(ns,2).EQ.0 )then
      neven = ns
      nodd = ns+1
  else
      nodd = ns
      neven = ns+1
  endif

  call mpp_get_compute_domains( grid_domain_y, begin=jstart(:) )
  spherical = cmplx(0.,0.)
  DOMAIN_LOOP: do jd=1,nd/2+1
!each domain jd has a mirror domain nd+1-jd
!if nd is even we stop at jd=nd/2; if nd is odd we stop at j=nj/2
     if( jd.EQ.nd/2+1 .AND. mod(nd,2).EQ.0 )exit DOMAIN_LOOP
     joff = jstart(jd) - 1
     do j = 1,nj
        if( jd.EQ.nd/2+1 .AND. j.GT.nj/2 )exit DOMAIN_LOOP
        if(south_to_north_local) then
            do k = 1,num_levels
               do m = ms,me
                  x_even = fourier(m,nj+1-j,k,nd+1-jd) + fourier(m,j,k,jd)
                  x_odd  = fourier(m,nj+1-j,k,nd+1-jd) - fourier(m,j,k,jd)
                  do n=neven,ne,2
                     spherical(m,n,k) = spherical(m,n,k) + x_even*legendre_wts(m,n,j+joff)
                  end do
                  do n=nodd, ne,2                           
                     spherical(m,n,k) = spherical(m,n,k) + x_odd *legendre_wts(m,n,j+joff)
                  end do
               end do
            enddo
        else
            do k = 1,num_levels
               do m = ms,me
                  x_even = fourier(m,j,k,jd) + fourier(m,nj+1-j,k,nd+1-jd)
                  x_odd  = fourier(m,j,k,jd) - fourier(m,nj+1-j,k,nd+1-jd)
                  do n=neven,ne,2
                     spherical(m,n,k) = spherical(m,n,k) + x_even*legendre_wts(m,n,j+joff)
                  end do
                  do n=nodd, ne,2                           
                     spherical(m,n,k) = spherical(m,n,k) + x_odd *legendre_wts(m,n,j+joff)
                  end do
               end do
            enddo
        end if
     end do
  end do DOMAIN_LOOP

  return
end subroutine trans_fourier_to_spherical_3d

!-------------------------------------------------------------------------
subroutine trans_spherical_to_fourier_2d(spherical,fourier) 
!-------------------------------------------------------------------------

complex, intent(in),  dimension(:,:) :: spherical
complex, intent(out), dimension(:,:,:) :: fourier

complex, dimension(size(spherical,1), size(spherical,2), 1) :: spherical_3d
complex, dimension(size(fourier,1),   size(fourier,2),   1, size(fourier,3)) :: fourier_3d

spherical_3d(:,:,1) = spherical(:,:)
call trans_spherical_to_fourier_3d(spherical_3d,fourier_3d)
fourier(:,:,:) = fourier_3d(:,:,1,:)

return
end subroutine trans_spherical_to_fourier_2d

!-----------------------------------------------------------------------
subroutine trans_fourier_to_spherical_2d(fourier,spherical)
!-----------------------------------------------------------------------

complex, intent (in), dimension(:,:,:) :: fourier
complex, intent (out),  dimension(:,:) :: spherical

complex, dimension(size(spherical,1), size(spherical,2), 1) :: spherical_3d
complex, dimension(size(fourier,1),   size(fourier,2),   1, size(fourier,3)) :: fourier_3d

fourier_3d(:,:,1,:) = fourier(:,:,:)
call trans_fourier_to_spherical_3d(fourier_3d,spherical_3d)
spherical(:,:) = spherical_3d(:,:,1)
 
return
end subroutine trans_fourier_to_spherical_2d

!-----------------------------------------------------------------------
subroutine define_legendre
!-----------------------------------------------------------------------

integer j
real, dimension(0:num_fourier,0:num_spherical,lat_max/2) :: legendre_global

allocate(legendre(ms:me,ns:ne,lat_max/2))

allocate(legendre_wts(ms:me,ns:ne,lat_max/2))

call compute_legendre(legendre_global, num_fourier, fourier_inc, num_spherical, sin_hem, lat_max/2)
legendre = legendre_global(ms:me,ns:ne,:)

do j=1,lat_max/2
  legendre_wts(:,:,j) = legendre(:,:,j)*wts_lat(j)
end do

return
end subroutine define_legendre

!----------------------------------------------------------------------
subroutine define_gaussian
!----------------------------------------------------------------------

integer j
real, dimension(lat_max/2) :: wts_hem

allocate (sin_lat(lat_max))
allocate (cos_lat(lat_max))
allocate (cosm_lat(lat_max))
allocate (cosm2_lat(lat_max))
allocate (wts_lat(lat_max))
allocate (deg_lat(lat_max))
allocate (sin_hem(lat_max/2))

call compute_gaussian(sin_hem, wts_hem, lat_max/2)

if(south_to_north_local) then
   sin_lat(1:lat_max/2)   = - sin_hem
else
   sin_lat(1:lat_max/2)   =   sin_hem
end if

do j=1,lat_max/2
  sin_lat(lat_max+1-j) = - sin_lat(j)
  wts_lat(j)           =   wts_hem(j)
  wts_lat(lat_max+1-j) =   wts_hem(j)
end do

cos_lat   = sqrt(1-sin_lat*sin_lat)
cosm_lat  = 1./cos_lat
cosm2_lat = 1./(cos_lat*cos_lat)
deg_lat   = asin(sin_lat)*180.0/pi
  
return
end subroutine define_gaussian

!-----------------------------------------------------------------------
subroutine get_south_to_north(south_to_north_out)
!-----------------------------------------------------------------------

logical, intent (out) :: south_to_north_out

if(.not. module_is_initialized) then
  call error_mesg('get_south_to_north','failed to define package', FATAL)
end if

south_to_north_out = south_to_north_local

return
end subroutine get_south_to_north

!-----------------------------------------------------------------------
subroutine get_sin_lat(sin_lat_out)
!-----------------------------------------------------------------------

real, intent (out), dimension(:) :: sin_lat_out

if(.not. module_is_initialized) then
  call error_mesg('get_sin_lat','failed to define package', FATAL)
end if

if(size(sin_lat_out,1).eq.lat_max) then
    sin_lat_out = sin_lat
else                            !assume grid compute domain
    sin_lat_out = sin_lat(js:je)
end if

return
end subroutine get_sin_lat

!-----------------------------------------------------------------------
subroutine get_cos_lat(cos_lat_out)
!-----------------------------------------------------------------------

real, intent (out), dimension(:) :: cos_lat_out

if(.not. module_is_initialized) then
  call error_mesg('get_cos_lat','failed to define package', FATAL)
end if
  
if(size(cos_lat_out,1).eq.lat_max) then
    cos_lat_out = cos_lat
else                            !assume grid compute domain
    cos_lat_out = cos_lat(js:je)
end if

return
end subroutine get_cos_lat

!-----------------------------------------------------------------------
subroutine get_cosm_lat(cosm_lat_out)
!-----------------------------------------------------------------------

real, intent (out), dimension(:) :: cosm_lat_out

if(.not. module_is_initialized) then
  call error_mesg('get_cosm_lat','failed to define package', FATAL)
end if
  
if(size(cosm_lat_out,1).eq.lat_max) then
    cosm_lat_out = cosm_lat
else                            !assume grid compute domain
    cosm_lat_out = cosm_lat(js:je)
end if

return
end subroutine get_cosm_lat

!-----------------------------------------------------------------------
subroutine get_cosm2_lat(cosm2_lat_out)
!-----------------------------------------------------------------------

real, intent (out), dimension(:) :: cosm2_lat_out

if(.not. module_is_initialized) then 
  call error_mesg('get_cosm2_lat','failed to define package', FATAL)
end if
  
if(size(cosm2_lat_out,1).eq.lat_max) then
    cosm2_lat_out = cosm2_lat
else                            !assume grid compute domain
    cosm2_lat_out = cosm2_lat(js:je)
end if

return
end subroutine get_cosm2_lat

!-----------------------------------------------------------------------
subroutine get_deg_lat(deg_lat_out)
!-----------------------------------------------------------------------

real, intent (out), dimension(:) :: deg_lat_out

if(.not. module_is_initialized) then
  call error_mesg('get_deg_lat','failed to define package', FATAL)
end if
  
if(size(deg_lat_out,1).eq.lat_max) then
    deg_lat_out = deg_lat
else                            !assume grid compute domain
    deg_lat_out = deg_lat(js:je)
end if

return
end subroutine get_deg_lat

!-----------------------------------------------------------------------
subroutine get_wts_lat(wts_lat_out)
!-----------------------------------------------------------------------

real, intent (out), dimension(:) :: wts_lat_out

if(.not. module_is_initialized) then
  call error_mesg('get_wts_lat','failed to define package', FATAL)
end if
  
if(size(wts_lat_out,1).eq.lat_max) then
    wts_lat_out = wts_lat
else                            !assume grid compute domain
    wts_lat_out = wts_lat(js:je)
end if

return
end subroutine get_wts_lat

!----------------------------------------------------------------------
subroutine spherical_fourier_end

if(.not. module_is_initialized) return

deallocate(sin_lat, cos_lat, cosm_lat, cosm2_lat, wts_lat, deg_lat)
deallocate(sin_hem, legendre, legendre_wts)
call spherical_end
module_is_initialized = .false.

return
end subroutine spherical_fourier_end
!----------------------------------------------------------------------
end module spherical_fourier_mod
