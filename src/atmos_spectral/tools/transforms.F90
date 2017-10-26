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

module transforms_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number, stdlog, close_file, check_nml_error

use mpp_mod, only: mpp_chksum, mpp_error, mpp_npes, mpp_sum, mpp_sync, mpp_sync_self, mpp_transmit

use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, mpp_get_domain_components, mpp_get_layout, &
                           mpp_get_pelist, mpp_update_domains, domain1D, XUPDATE, mpp_global_field

use spec_mpp_mod,  only: spec_mpp_init, spec_mpp_end, get_grid_domain, get_spec_domain, &
                         grid_domain, spectral_domain, global_spectral_domain, atmosphere_domain

use constants_mod, only: pi


use spherical_fourier_mod, only: spherical_fourier_init, spherical_fourier_end, &
                                 trans_spherical_to_fourier, trans_fourier_to_spherical, &
                                 get_south_to_north, &
                                 get_sin_lat, get_cos_lat, get_cosm_lat, get_cosm2_lat, &
                                 get_deg_lat, get_wts_lat, &
                                 spherical_init, spherical_end, &
                                 compute_lon_deriv_cos, compute_lat_deriv_cos, &
                                 compute_laplacian, compute_vor, compute_div, &
                                 get_spherical_wave, get_fourier_wave, &
                                 get_eigen_laplacian, compute_gradient_cos, &
                                 compute_ucos_vcos, compute_vor_div, &
                                 triangular_truncation, rhomboidal_truncation, &
                                 compute_legendre, compute_gaussian


use grid_fourier_mod, only: grid_fourier_init, grid_fourier_end, trans_grid_to_fourier, &
                            trans_fourier_to_grid, get_lon_max, get_longitude_origin, get_deg_lon

!---------------------------------------------------------------------------
! This module is initialized by transforms_init
!
! The partial re-initializations required when one changes either the
!   the number of longitudes in the transform grid
!   is provided by reset_num_lon_in_transform
!
! Fields are transformed from the transform grid to spherical harmonics
!   and back by trans_spherical_to_grid and trans_grid_to_spherical
!
! The two utilities, divide_by_cos and divide_by_cos2, provide convenient 
!   ways of dividing grid fields by cos(lat) and cos(lat)**2, 
!   as is often required when transforming to and fro in different contexts,
!   without having to compute or retrieve the gaussian latitudes explicitly.  
!
! trans_filter transforms a grid field to the spectral domain, multiplies
!   by a filter in the spectral domain, and then transforms back
!
!---------------------------------------------------------------------------

implicit none
private

interface trans_spherical_to_grid
   module procedure trans_spherical_to_grid_3d,  &
                    trans_spherical_to_grid_2d
end interface

interface trans_grid_to_spherical
   module procedure trans_grid_to_spherical_3d,  &
                    trans_grid_to_spherical_2d
end interface

interface trans_filter
   module procedure trans_filter_3d,  &
                    trans_filter_2d
end interface

interface divide_by_cos
   module procedure divide_by_cos_3d,  &
                    divide_by_cos_2d
end interface

interface divide_by_cos2
   module procedure divide_by_cos2_3d,  &
                    divide_by_cos2_2d
end interface

interface uv_grid_from_vor_div
   module procedure uv_grid_from_vor_div_2d,  &
                    uv_grid_from_vor_div_3d
end interface

interface vor_div_from_uv_grid
   module procedure vor_div_from_uv_grid_2d,  &
                    vor_div_from_uv_grid_3d
end interface

interface horizontal_advection
   module procedure horizontal_advection_2d,  &
                    horizontal_advection_3d
end interface

character(len=128), parameter :: version=&
'$Id: transforms.F90,v 19.0 2012/01/06 20:36:20 fms Exp $'

character(len=128), parameter :: tagname=&
'$Name: siena_201211 $'

! ---------------------------------------------------

!  public interfaces defined in this module

public :: transforms_init
public :: transforms_end
public :: transforms_are_initialized
public :: reset_num_lon_in_transform
public :: trans_spherical_to_grid
public :: trans_grid_to_spherical
public :: divide_by_cos
public :: divide_by_cos2
public :: trans_filter
public :: get_lat_max
public :: get_triang_trunc
public :: get_num_fourier
public :: get_fourier_inc
public :: get_num_spherical
public :: get_grid_boundaries
public :: vor_div_from_uv_grid
public :: uv_grid_from_vor_div
public :: horizontal_advection
public :: area_weighted_global_mean

! ---------------------------------------------------

!  public interfaces carried forward from used modules

! From spherical_fourier_mod
public :: spherical_fourier_init, spherical_fourier_end
public :: trans_spherical_to_fourier, trans_fourier_to_spherical
public :: get_south_to_north
public :: get_sin_lat, get_cos_lat, get_cosm_lat, get_cosm2_lat
public :: get_deg_lat, get_wts_lat

! From grid_fourier_mod
public :: grid_fourier_init, grid_fourier_end
public :: trans_grid_to_fourier, trans_fourier_to_grid
public :: get_lon_max, get_deg_lon, get_longitude_origin

! From spherical_mod
public :: spherical_init, spherical_end
public :: compute_lon_deriv_cos, compute_lat_deriv_cos
public :: compute_laplacian, compute_vor, compute_div
public :: get_spherical_wave, get_fourier_wave
public :: get_eigen_laplacian, compute_gradient_cos
public :: compute_ucos_vcos, compute_vor_div
public :: triangular_truncation, rhomboidal_truncation

public :: compute_legendre, compute_gaussian

! From spec_mpp_mod
public :: spec_mpp_init, get_grid_domain, get_spec_domain
public :: grid_domain,  spectral_domain,  global_spectral_domain
public :: atmosphere_domain

! ---------------------------------------------------

!  The next six variables are initialized to zero so
!  as to allow a check that the user has supplied them.
integer :: lat_max=0
integer :: num_lon=0
integer :: num_fourier=0
integer :: fourier_inc=0
integer :: num_spherical=0

logical :: south_to_north_local
logical :: triang_trunc_local
logical :: make_symmetric_local
real :: longitude_origin_local

integer :: trunc_fourier
logical :: module_is_initialized = .false.

integer :: npes
integer, dimension(2) :: grid_layout, spectral_layout
integer :: xmaxsize, ymaxsize   !used for dimensioning in transpose routines
logical :: debug=.FALSE.
integer :: ms, me, ns, ne, is, ie, js, je

real, allocatable, dimension(:) :: lat_boundaries_global, lon_boundaries_global
real :: global_sum_of_wts

logical :: check_fourier_imag=.false.
namelist / transforms_nml / check_fourier_imag

contains

!---------------------------------------------------------------------------
subroutine transforms_init(radius,            &
                           lat_max_in,        &
                           num_lon_in,        &
                           num_fourier_in,    &
                           fourier_inc_in,    &
                           num_spherical_in,  &
                           south_to_north,    &
                           triang_trunc,      &
                           longitude_origin,  &
						   make_symmetric)  !
!---------------------------------------------------------------------------

real,    intent(in) :: radius
integer, intent(in) :: lat_max_in
integer, intent(in) :: num_lon_in
integer, intent(in) :: num_fourier_in
integer, intent(in) :: fourier_inc_in
integer, intent(in) :: num_spherical_in

logical, intent(in), optional :: south_to_north, triang_trunc, make_symmetric
real,    intent(in), optional :: longitude_origin
integer :: namelist_unit, ierr, io

real, allocatable :: wts_lat(:)
real :: sum_wts, del_lon
integer :: i, j

if(module_is_initialized) return

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=transforms_nml, iostat=io)
    ierr = check_nml_error(io, 'transforms_nml')
#else
    namelist_unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(namelist_unit, nml=transforms_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'transforms_nml')
    enddo
20  call close_file (namelist_unit)
#endif

call write_version_number(version, tagname)
if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=transforms_nml)

npes = mpp_npes()

!  make local copies of input

lat_max       = lat_max_in
num_lon       = num_lon_in
num_fourier   = num_fourier_in
fourier_inc   = fourier_inc_in
num_spherical = num_spherical_in

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

if(present(triang_trunc)) then
  triang_trunc_local = triang_trunc
else
  triang_trunc_local = .true.
end if

if(present(longitude_origin)) then
  longitude_origin_local = longitude_origin
else
  longitude_origin_local = 0.0
end if

call spec_mpp_init( num_fourier, num_spherical, num_lon, lat_max )
call get_grid_domain(is, ie, js, je)
call get_spec_domain(ms, me, ns, ne)
   
! initialize spherical_fourier (which initializes spherical)
call spherical_fourier_init(radius, lat_max, num_fourier, fourier_inc, num_spherical, &
                            south_to_north=south_to_north_local, make_symmetric=make_symmetric_local)   

trunc_fourier = num_fourier

call grid_fourier_init(num_lon, fourier_inc, &
     check=check_fourier_imag, longitude_origin=longitude_origin_local)

allocate( wts_lat(lat_max) )
call get_wts_lat(wts_lat)
allocate( lon_boundaries_global(num_lon+1) )
allocate( lat_boundaries_global(lat_max+1) )
lat_boundaries_global(1) = -.5*pi
sum_wts = 0.
do j=1,lat_max-1
  sum_wts = sum_wts + wts_lat(j)
  lat_boundaries_global(j+1) = asin(sum_wts-1.)
end do
lat_boundaries_global(lat_max+1) = .5*pi
if (.not. south_to_north_local) then
  lat_boundaries_global(:) = -lat_boundaries_global(:)
end if
del_lon = 2*pi/num_lon
do i=1,num_lon+1
  lon_boundaries_global(i) = longitude_origin_local + (i-1.5)*del_lon
end do

global_sum_of_wts = sum(wts_lat)

call mpp_get_layout( grid_domain, grid_layout )
call mpp_get_layout( spectral_domain, spectral_layout )
call mpp_get_compute_domain( spectral_domain, xmax_size=xmaxsize )
call mpp_get_compute_domain( grid_domain, ymax_size=ymaxsize )
module_is_initialized = .true.

return
end subroutine transforms_init

!--------------------------------------------------------------------------------------------------------------
logical function transforms_are_initialized()
!--------------------------------------------------------------------------------------------------------------

transforms_are_initialized = module_is_initialized
return

end function transforms_are_initialized

!--------------------------------------------------------------------------------------------------------------
subroutine reset_num_lon_in_transform(num_lon_in, trunc_fourier_in, longitude_origin)
!--------------------------------------------------------------------------------------------------------------

integer, intent(in) :: num_lon_in
integer, intent(in) :: trunc_fourier_in
real,    intent(in), optional    :: longitude_origin

num_lon = num_lon_in
trunc_fourier = trunc_fourier_in

if(present(longitude_origin)) then
  longitude_origin_local = longitude_origin
else
  longitude_origin_local = 0.0
end if

if(trunc_fourier > num_fourier) then
  call error_mesg('reset_num_lon_in_transform','trunc_fourier > num_fourier', FATAL)
end if

call grid_fourier_init(num_lon, fourier_inc, &
      check=check_fourier_imag, longitude_origin=longitude_origin_local)
if( npes.GT.1 )call abort()     !need to think about this still

return
end subroutine reset_num_lon_in_transform

!--------------------------------------------------------------------------
 subroutine trans_spherical_to_grid_3d(spherical, grid)
!--------------------------------------------------------------------------

complex, intent (in), dimension (ms:,ns:,:) :: spherical
real, intent(out), dimension (is:,:,:) :: grid
real, dimension(num_lon,size(grid,2),size(grid,3)) :: grid_xglobal

integer(kind=kind(spherical)) :: c1, c2, c3

complex, dimension (0:num_lon/2,          js:je, size(grid,3)) :: fourier_g
complex, dimension (ms:me, je-js+1, size(grid,3), grid_layout(2)) :: fourier_s
logical :: grid_x_is_global, spectral_y_is_global
type(domain1D) :: spectral_domain_y
integer, allocatable :: pelist(:)

if(.not.module_is_initialized) then
  call error_mesg('trans_spherical_to_grid','transforms module is not initialized', FATAL)
end if

if( size(grid,1).ne.ie-is+1 )&
     call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_GRID: size(grid,1).ne.ie-is+1.' )

if( size(grid,2).ne.je-js+1 )&
     call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_GRID: size(grid,2).ne.je-js+1.' )

if( size(spherical,1).ne.me-ms+1 )&
     call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_GRID: size(spherical,1).ne.me-ms+1.' )

if( size(spherical,2).ne.ne-ns+1 )&
     call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_GRID: size(spherical,2).ne.ne-ns+1.' )

if( size(spherical,3).ne.size(grid,3) )&
     call mpp_error( FATAL, 'TRANS_SPHERICAL_TO_GRID: size(spherical,3).ne.size(grid,3.' )

call trans_spherical_to_fourier( spherical, fourier_s )
call mpp_get_compute_domain( spectral_domain, y_is_global=spectral_y_is_global )
if( .NOT.spectral_y_is_global )then
    allocate( pelist(spectral_layout(2)) )
    call mpp_get_domain_components( spectral_domain, y=spectral_domain_y )
    call mpp_get_pelist( spectral_domain_y, pelist )
    call mpp_sum( fourier_s          , size(fourier_s(:,:,:,:)), pelist )
end if

call reverse_transpose_fourier( fourier_s, fourier_g )

fourier_g(trunc_fourier+1:num_lon/2,:,:) = cmplx(0.,0.)
call mpp_get_compute_domain( grid_domain, x_is_global=grid_x_is_global )
if( .NOT.grid_x_is_global )then
    grid_xglobal = trans_fourier_to_grid(fourier_g)
    grid(is:ie,:,:) = grid_xglobal(is:ie,:,:)
else
    grid = trans_fourier_to_grid(fourier_g)
endif

if( debug )then
    c1 = mpp_chksum(spherical)
    c2 = mpp_chksum(fourier_g)
    c3 = mpp_chksum(grid)
    write(  0,'(a,i2,3z18,28i4)' )'S2G: ', mpp_pe(), c1, c2, c3, lbound(spherical), ubound(spherical), &
         lbound(fourier_s), ubound(fourier_s), lbound(fourier_g), ubound(fourier_g), lbound(grid), ubound(grid)
end if

return
end subroutine trans_spherical_to_grid_3d

!--------------------------------------------------------------------------
 subroutine trans_spherical_to_grid_2d(spherical, grid)
!--------------------------------------------------------------------------

complex, intent (in), dimension (:,:)  ::  spherical
real, intent(out), dimension (:,:) :: grid

real, dimension (size(grid,1), size(grid,2), 1) :: grid_3d
complex, dimension(size(spherical,1), size(spherical,2), 1) :: spherical_3d

spherical_3d(:,:,1) = spherical(:,:)
call trans_spherical_to_grid_3d(spherical_3d, grid_3d)
grid(:,:) = grid_3d(:,:,1)

return
end subroutine trans_spherical_to_grid_2d

!-------------------------------------------------------------------------
subroutine trans_grid_to_spherical_3d(grid, spherical, do_truncation)
!-------------------------------------------------------------------------

real, intent(in), dimension (is:,:,:) :: grid
logical, intent(in), optional :: do_truncation
complex, intent(inout), dimension (ms:,ns:,:) :: spherical
real, dimension(num_lon,size(grid,2),size(grid,3)) :: grid_xglobal

logical :: do_truncation_local, grid_x_is_global
complex, dimension (0:num_lon/2,          js:je, size(grid,3)) :: fourier_g
complex, dimension (ms:me, je-js+1, size(grid,3), grid_layout(2)) :: fourier_s

integer(kind=kind(spherical)) :: c1, c2, c3

if(.not.module_is_initialized) then
  call error_mesg('trans_grid_to_spherical','transforms module is not initialized', FATAL)
end if

if( size(grid,1).ne.ie-is+1 )&
     call mpp_error( FATAL, 'TRANSFORMS: size(grid,1).ne.ie-is+1.' )

if( size(grid,2).ne.je-js+1 )&
     call mpp_error( FATAL, 'TRANSFORMS: size(grid,2).ne.je-js+1.' )

if( size(spherical,1).ne.me-ms+1 )&
     call mpp_error( FATAL, 'TRANSFORMS: size(spherical,1).ne.me-ms+1.' )

if( size(spherical,2).ne.ne-ns+1 )&
     call mpp_error( FATAL, 'TRANSFORMS: size(spherical,2).ne.ne-ns+1.' )

if( size(spherical,3).ne.size(grid,3) )&
     call mpp_error( FATAL, 'TRANSFORMS: size(spherical,3).ne.size(grid,3).' )

if(present(do_truncation)) then
  do_truncation_local = do_truncation
else
  do_truncation_local = .true.
end if

call mpp_get_compute_domain( grid_domain, x_is_global=grid_x_is_global )
if( .NOT.grid_x_is_global )then
    grid_xglobal(is:ie,:,:) = grid(is:ie,:,:)
    call mpp_update_domains( grid_xglobal, grid_domain, XUPDATE )
    fourier_g = trans_grid_to_fourier(grid_xglobal)
else
    fourier_g = trans_grid_to_fourier(grid)
endif

if( trunc_fourier.lt.me )fourier_g(trunc_fourier+1:me,:,:) = cmplx(0.,0.)

call transpose_fourier( fourier_g, fourier_s )
call trans_fourier_to_spherical(fourier_s, spherical)

if(do_truncation_local) then
  if(triang_trunc_local) then
    call triangular_truncation(spherical)
  else
    call rhomboidal_truncation(spherical)
  end if
end if

if( debug )then
    c1 = mpp_chksum(spherical)
    c2 = mpp_chksum(fourier_g)
    c3 = mpp_chksum(grid)
    write(  0,'(a,i2,3z18,28i4)' )'G2S: ', mpp_pe(), c1, c2, c3, lbound(spherical), ubound(spherical), &
         lbound(fourier_s), ubound(fourier_s), lbound(fourier_g), ubound(fourier_g), lbound(grid), ubound(grid)
    call mpp_error( FATAL )
end if

return
end subroutine trans_grid_to_spherical_3d

!-------------------------------------------------------------------------
subroutine trans_grid_to_spherical_2d(grid, spherical, do_truncation)
!-------------------------------------------------------------------------

real,    intent(in), dimension (:,:) :: grid
logical, intent(in), optional :: do_truncation
complex, intent(inout), dimension (:,:) :: spherical

real, dimension (size(grid,1), size(grid,2), 1) :: grid_3d
complex, dimension (size(spherical,1), size(spherical,2), 1) :: spherical_3d

grid_3d(:,:,1) = grid(:,:)
spherical_3d(:,:,1) = cmplx(0.0,0.0)
call trans_grid_to_spherical_3d(grid_3d, spherical_3d, do_truncation)
spherical(:,:) = spherical_3d(:,:,1)

return
end subroutine trans_grid_to_spherical_2d

!-------------------------------------------------------------------------
subroutine trans_filter_3d(grid, filter)
!-------------------------------------------------------------------------

real, intent(inout), dimension (:,:,:) :: grid
real, intent(in), optional, dimension (:,:) :: filter

complex, dimension (ms:me, ns:ne, size(grid,3)) :: spherical
integer :: k

if(.not.module_is_initialized) then
  call error_mesg('trans_filter','transforms module is not initialized', FATAL)
end if

call trans_grid_to_spherical(grid, spherical)

if(present(filter)) then
  do k=1, size(grid,3)
    spherical(ms:me,ns:ne,k) = spherical(ms:me,ns:ne,k)*filter
  end do
end if

call trans_spherical_to_grid(spherical, grid)

return

end subroutine trans_filter_3d

!-------------------------------------------------------------------------
subroutine trans_filter_2d(grid, filter)
!-------------------------------------------------------------------------

real, intent(inout), dimension (:,:) :: grid
real, intent(in), optional, dimension (:,:) :: filter

real, dimension(size(grid,1),size(grid,2),1) :: grid_3d

grid_3d(:,:,1) = grid
call trans_filter_3d(grid_3d, filter)
grid = grid_3d(:,:,1)

return
end subroutine trans_filter_2d

!-------------------------------------------------------------------------
subroutine divide_by_cos_3d(grid)
!-------------------------------------------------------------------------

real, intent (inout), dimension(:,:,:) :: grid

real, dimension(size(grid,2)) :: cosm_lat
integer :: i,j,k

if(.not.module_is_initialized) then
  call error_mesg('divide_by_cos','transforms module is not initialized', FATAL)
end if

call get_cosm_lat(cosm_lat)

do k=1,size(grid,3)
  do j=1,size(grid,2)
    do i=1,size(grid,1)
      grid(i,j,k)=grid(i,j,k)*cosm_lat(j)
    end do
  end do
end do

return
end subroutine divide_by_cos_3d

!-------------------------------------------------------------------------
subroutine divide_by_cos2_3d(grid)
!-------------------------------------------------------------------------

real, intent (inout), dimension(:,:,:) :: grid

real, dimension(size(grid,2)) :: cosm2_lat
integer :: i,j,k

if(.not.module_is_initialized) then
  call error_mesg('divide_by_cos2','transforms module is not initialized', FATAL)
end if

call get_cosm2_lat(cosm2_lat)

do k=1,size(grid,3)
  do j=1,size(grid,2)
    do i=1,size(grid,1)
      grid(i,j,k)=grid(i,j,k)*cosm2_lat(j)
    end do
  end do
end do

return
end subroutine divide_by_cos2_3d

!-------------------------------------------------------------------------
subroutine divide_by_cos_2d(grid)
!-------------------------------------------------------------------------

real, intent (inout), dimension(:,:) :: grid

real, dimension(size(grid,1),size(grid,2),1) :: grid_3d

grid_3d(:,:,1) = grid(:,:)
call divide_by_cos_3d(grid_3d)
grid(:,:) = grid_3d(:,:,1)

return
end subroutine divide_by_cos_2d

!-------------------------------------------------------------------------
subroutine divide_by_cos2_2d(grid)
!-------------------------------------------------------------------------

real, intent (inout), dimension(:,:) :: grid

real, dimension(size(grid,1),size(grid,2),1) :: grid_3d

grid_3d(:,:,1) = grid(:,:)
call divide_by_cos2_3d(grid_3d)
grid(:,:) = grid_3d(:,:,1)

return
end subroutine divide_by_cos2_2d

!-------------------------------------------------------------------------
subroutine uv_grid_from_vor_div_2d(vor_spec, div_spec, u_grid, v_grid)
!-------------------------------------------------------------------------

complex, intent(in) , dimension(:,:) :: vor_spec, div_spec
real,    intent(out), dimension(:,:) :: u_grid, v_grid

complex, dimension(size(vor_spec,1), size(vor_spec,2), 1) :: vor_spec_3d, div_spec_3d
real,    dimension(size(u_grid,  1), size(u_grid,  2), 1) ::   u_grid_3d,   v_grid_3d

vor_spec_3d(:,:,1) = vor_spec
div_spec_3d(:,:,1) = div_spec
call uv_grid_from_vor_div_3d(vor_spec_3d, div_spec_3d, u_grid_3d, v_grid_3d)
u_grid = u_grid_3d(:,:,1)
v_grid = v_grid_3d(:,:,1)

return
end subroutine uv_grid_from_vor_div_2d

!-------------------------------------------------------------------------
subroutine uv_grid_from_vor_div_3d(vor_spec, div_spec, u_grid, v_grid)
!-------------------------------------------------------------------------

complex, intent(in) , dimension(:,:,:) :: vor_spec, div_spec
real,    intent(out), dimension(:,:,:) :: u_grid, v_grid

complex , dimension(size(vor_spec,1), size(vor_spec,2), size(vor_spec,3)) :: dx_spec, dy_spec

if(.not.module_is_initialized) then
  call error_mesg('uv_grid_from_vor_div','transforms module is not initialized', FATAL)
end if

call compute_ucos_vcos          (vor_spec, div_spec, dx_spec, dy_spec)
call trans_spherical_to_grid    (dx_spec , u_grid)
call trans_spherical_to_grid    (dy_spec , v_grid)
call divide_by_cos              (u_grid)
call divide_by_cos              (v_grid)

return
end subroutine uv_grid_from_vor_div_3d

!-------------------------------------------------------------------------
subroutine vor_div_from_uv_grid_2d(u_grid, v_grid, vor_spec, div_spec, triang)
!-------------------------------------------------------------------------

real,    intent(in),    dimension(:,:) :: u_grid, v_grid
complex, intent(out  ), dimension(:,:) :: vor_spec, div_spec
logical, intent(in), optional :: triang

complex, dimension(size(vor_spec,1), size(vor_spec,2), 1) :: vor_spec_3d, div_spec_3d
real,    dimension(size(u_grid,  1), size(u_grid,  2), 1) ::   u_grid_3d,   v_grid_3d

u_grid_3d(:,:,1) = u_grid
v_grid_3d(:,:,1) = v_grid
call vor_div_from_uv_grid_3d(u_grid_3d, v_grid_3d, vor_spec_3d, div_spec_3d, triang)
vor_spec = vor_spec_3d(:,:,1)
div_spec = div_spec_3d(:,:,1)

return
end subroutine vor_div_from_uv_grid_2d

!-------------------------------------------------------------------------
subroutine vor_div_from_uv_grid_3d(u_grid, v_grid, vor_spec, div_spec, triang)
!-------------------------------------------------------------------------

real,    intent(in),  dimension(:,:,:) :: u_grid, v_grid
complex, intent(out), dimension(:,:,:) :: vor_spec, div_spec
logical, intent(in), optional :: triang

complex, dimension(size(vor_spec,1), size(vor_spec,2), size(vor_spec,3)) :: dx_spec, dy_spec
real   , dimension(size(u_grid  ,1), size(u_grid  ,2), size(u_grid  ,3)) :: grid_tmp

logical :: do_triang

if(.not.module_is_initialized) then
  call error_mesg('vor_div_from_uv_grid','transforms module is not initialized', FATAL)
end if

if(present(triang)) then
  do_triang = triang
else
  do_triang = .true.
endif

grid_tmp = u_grid
call divide_by_cos(grid_tmp)
call trans_grid_to_spherical(grid_tmp, dx_spec, do_truncation=.false.)

grid_tmp = v_grid
call divide_by_cos(grid_tmp)
call trans_grid_to_spherical(grid_tmp, dy_spec, do_truncation=.false.)

call compute_vor_div(dx_spec, dy_spec, vor_spec, div_spec)

if(do_triang) then
  call triangular_truncation(vor_spec)
  call triangular_truncation(div_spec)
else
  call rhomboidal_truncation(vor_spec)
  call rhomboidal_truncation(div_spec)
endif

return
end subroutine vor_div_from_uv_grid_3d

!-------------------------------------------------------------------------
subroutine horizontal_advection_2d(field_spec, u_grid, v_grid, tendency)
!-------------------------------------------------------------------------

complex, intent(in),    dimension(:,:) :: field_spec
real,    intent(in),    dimension(:,:) :: u_grid, v_grid
real,    intent(inout), dimension(:,:) :: tendency

complex, dimension(size(field_spec,1), size(field_spec,2), 1) :: field_spec_3d
real,    dimension(size(u_grid,    1), size(u_grid,    2), 1) :: u_grid_3d, v_grid_3d
real,    dimension(size(u_grid,    1), size(u_grid,    2), 1) :: tendency_3d

field_spec_3d(:,:,1) = field_spec
u_grid_3d(:,:,1) = u_grid
v_grid_3d(:,:,1) = v_grid
tendency_3d(:,:,1) = tendency
call horizontal_advection_3d(field_spec_3d, u_grid_3d, v_grid_3d, tendency_3d)
tendency = tendency_3d(:,:,1)

return
end subroutine horizontal_advection_2d

!-------------------------------------------------------------------------
subroutine horizontal_advection_3d(field_spec, u_grid, v_grid, tendency)
!-------------------------------------------------------------------------

complex, intent(in),    dimension(:,:,:) :: field_spec
real,    intent(in),    dimension(:,:,:) :: u_grid, v_grid
real,    intent(inout), dimension(:,:,:) :: tendency

real,    dimension(size(u_grid,    1), size(u_grid,    2), size(u_grid,    3)) :: dx_grid, dy_grid 
complex, dimension(size(field_spec,1), size(field_spec,2), size(field_spec,3)) :: dx_spec, dy_spec 

if(.not.module_is_initialized) then
  call error_mesg('horizontal_advection','transforms module is not initialized', FATAL)
end if

call compute_gradient_cos(field_spec, dx_spec, dy_spec)
call trans_spherical_to_grid(dx_spec, dx_grid)
call trans_spherical_to_grid(dy_spec, dy_grid)

call divide_by_cos(dx_grid)
call divide_by_cos(dy_grid)
tendency = tendency - u_grid*dx_grid - v_grid*dy_grid

return
end subroutine horizontal_advection_3d

!-------------------------------------------------------------------------
subroutine get_lat_max(lat_max_out)
!-------------------------------------------------------------------------

integer, intent(out) :: lat_max_out

if(.not.module_is_initialized) then
  call error_mesg('get_lat_max','transforms module is not initialized', FATAL)
end if

lat_max_out = lat_max

return
end subroutine get_lat_max

!-------------------------------------------------------------------------
subroutine get_triang_trunc(triang_trunc_out)
!-------------------------------------------------------------------------

logical, intent(out) :: triang_trunc_out

if(.not.module_is_initialized) then
  call error_mesg('get_triang_trunc','transforms module is not initialized', FATAL)
end if

triang_trunc_out = triang_trunc_local

return
end subroutine get_triang_trunc

!-------------------------------------------------------------------------
subroutine get_num_fourier(num_fourier_out)
!-------------------------------------------------------------------------

integer, intent(out) :: num_fourier_out

if(.not.module_is_initialized) then
  call error_mesg('get_num_fourier','transforms module is not initialized', FATAL)
end if
    
num_fourier_out = num_fourier

return
end subroutine get_num_fourier

!-------------------------------------------------------------------------
subroutine get_fourier_inc(fourier_inc_out)
!-------------------------------------------------------------------------

integer, intent(out) :: fourier_inc_out

if(.not.module_is_initialized) then
  call error_mesg('get_fourier_inc','transforms module is not initialized', FATAL)
end if
    
fourier_inc_out = fourier_inc

return
end subroutine get_fourier_inc

!-------------------------------------------------------------------------
subroutine get_num_spherical(num_spherical_out)
!-------------------------------------------------------------------------

integer, intent(out) :: num_spherical_out

if(.not.module_is_initialized) then
  call error_mesg('get_num_spherical','transforms module is not initialized', FATAL)
end if
    
num_spherical_out = num_spherical

return
end subroutine get_num_spherical

!-------------------------------------------------------------------------
subroutine get_grid_boundaries(lon_boundaries, lat_boundaries,global)
!-------------------------------------------------------------------------

  real, intent(out), dimension(:) :: lon_boundaries, lat_boundaries
  logical,intent(in),optional :: global

  logical :: global_tmp
  character(len=3) :: chtmp1, chtmp2

  if(.not.module_is_initialized) then
     call error_mesg('get_grid_boundaries','transforms module is not initialized', FATAL)
  end if

  if (present(global)) then
     global_tmp = global
  else
     global_tmp = .false.
  endif

  if (.not. global_tmp) then
     if(size(lon_boundaries,1) /= ie-is+2) then
        write(chtmp1,'(i3)') size(lon_boundaries,1)
        write(chtmp2,'(i3)') ie-is+2
        call error_mesg('get_grid_boundaries','size(lon_boundaries) is incorrect. size(lon_boundaries)=' &
                       & //chtmp1//'  Should be'//chtmp2, FATAL)
     endif

     if(size(lat_boundaries,1) /= je-js+2) then
        write(chtmp1,'(i3)') size(lat_boundaries,1)
        write(chtmp2,'(i3)') je-js+2
        call error_mesg('get_grid_boundaries','size(lat_boundaries) is incorrect. size(lat_boundaries)=' &
                       & //chtmp1//'  Should be'//chtmp2, FATAL)
     endif

  else !global call
     if(size(lon_boundaries,1) /= num_lon+1) then
        write(chtmp1,'(i3)') size(lon_boundaries,1)
        write(chtmp2,'(i3)') num_lon+1
        call error_mesg('get_grid_boundaries','size(lon_boundaries) is incorrect. size(lon_boundaries)=' &
                       & //chtmp1//'  Should be'//chtmp2, FATAL)
     endif

     if(size(lat_boundaries,1) /= lat_max+1) then
        write(chtmp1,'(i3)') size(lat_boundaries,1)
        write(chtmp2,'(i3)') lat_max+1
        call error_mesg('get_grid_boundaries','size(lat_boundaries) is incorrect. size(lat_boundaries)=' &
                       & //chtmp1//'  Should be'//chtmp2, FATAL)
     endif
  endif

  if (global_tmp) then
     lon_boundaries = lon_boundaries_global
     lat_boundaries = lat_boundaries_global
  else
     lon_boundaries = lon_boundaries_global(is:ie+1)
     lat_boundaries = lat_boundaries_global(js:je+1)
  endif
  return
end subroutine get_grid_boundaries

!-------------------------------------------------------------------------
subroutine reverse_transpose_fourier( fourier_s, fourier_g )
!-------------------------------------------------------------------------
  complex, intent(in)  :: fourier_s(:,:,:,0:)
  complex, intent(out) :: fourier_g(0:,:,:)
  complex, dimension(xmaxsize*ymaxsize*size(fourier_s,3)) :: get_data
  integer :: i,j,k, jj, jp, jm, pp, pm, nput, nget, jpos
  type(domain1D) :: spectral_domain_x, grid_domain_y
  integer, dimension(0:grid_layout(2)-1) :: pelist, ygridsize, xspecsize, xsbegin, xsend

  if(.not.module_is_initialized) then
    call error_mesg('reverse_transpose_fourier','transforms module is not initialized', FATAL)
  end if

  call mpp_get_domain_components( grid_domain, y=grid_domain_y )
  call mpp_get_domain_components( spectral_domain, x=spectral_domain_x )
  call mpp_get_pelist( grid_domain_y, pelist, jpos )
  call mpp_get_compute_domains( grid_domain_y, size=ygridsize )
  call mpp_get_compute_domains( spectral_domain_x, xsbegin, xsend, xspecsize )
  nput = size(fourier_s,1)*size(fourier_s,2)*size(fourier_s,3)
  fourier_g(ms:me,:,:) = fourier_s(:,:,:,jpos)
  do jj = 1,grid_layout(2)-1
     jp = mod(jpos+jj,grid_layout(2))
     jm = mod(jpos-jj+grid_layout(2),grid_layout(2))
     pp = pelist(jp)
     pm = pelist(jm)
     nget = xspecsize(jm)*ygridsize(jm)*size(fourier_s,3)
     ! Force use of "scalar", integer pointer mpp interface
     call mpp_transmit( put_data=fourier_s(1,1,1,jp), plen=nput, to_pe=pp, &
                        get_data=get_data(1), glen=nget, from_pe=pm )
     nget = 0
     do k = 1,size(fourier_g,3)
        do j = 1,size(fourier_g,2)
           do i = xsbegin(jm),xsend(jm)
              nget = nget + 1
              fourier_g(i,j,k) = get_data(nget)
           end do
        end do
     end do
  end do
  call mpp_sync()
  return
end subroutine reverse_transpose_fourier

!-------------------------------------------------------------------------
subroutine transpose_fourier( fourier_g, fourier_s )
!-------------------------------------------------------------------------
  complex, intent(in)   :: fourier_g(0:,:,:)
  complex, intent(out)  :: fourier_s(:,:,:,0:)
  complex, dimension(xmaxsize*ymaxsize*size(fourier_s,3)) :: put_data
  integer :: i,j,k, ii, ip, im, pp, pm, nput, nget, ipos, jp
  type(domain1D) :: spectral_domain_x, grid_domain_y
  integer, dimension(0:spectral_layout(1)-1) :: pelist, ygridsize, xspecsize, xsbegin, xsend

  if(.not.module_is_initialized) then
    call error_mesg('transpose_fourier','transforms module is not initialized', FATAL)
  end if

  call mpp_get_domain_components( grid_domain, y=grid_domain_y )
  call mpp_get_domain_components( spectral_domain, x=spectral_domain_x )
  call mpp_get_compute_domains( grid_domain_y, size=ygridsize )
  call mpp_get_compute_domains( spectral_domain_x, xsbegin, xsend, xspecsize )
  nget = size(fourier_s,1)*size(fourier_s,2)*size(fourier_s,3)
  call mpp_get_pelist( grid_domain_y, pelist, jp )
  fourier_s(:,:,:,jp) = fourier_g(ms:me,:,:)
  call mpp_get_pelist( spectral_domain_x, pelist, ipos )
  do ii = 1,spectral_layout(1)-1
     ip = mod(ipos+ii,spectral_layout(1))
     im = mod(ipos-ii+spectral_layout(1),spectral_layout(1))
     pp = pelist(ip)
     pm = pelist(im)
     nput = 0
     call mpp_sync_self()
     do k = 1,size(fourier_g,3)
        do j = 1,size(fourier_g,2)
           do i = xsbegin(ip),xsend(ip)
              nput = nput + 1
              put_data(nput) = fourier_g(i,j,k)
           end do
        end do
     end do
     ! Force use of "scalar", integer pointer mpp interface
     call mpp_transmit( put_data=put_data(1), plen=nput, to_pe=pp, &
                        get_data=fourier_s(1,1,1,im), glen=nget, from_pe=pm )
  end do
  call mpp_sync()
  return
end subroutine transpose_fourier

!-------------------------------------------------------------------------
function area_weighted_global_mean(field)
!-------------------------------------------------------------------------
real :: area_weighted_global_mean
real, intent(in), dimension(:,:) :: field
real, dimension(size(field,2))   :: wts_lat
real, dimension(size(field,1), size(field,2)) :: weighted_field_local
real, dimension(num_lon, lat_max) :: weighted_field_global
integer :: j

call get_wts_lat(wts_lat)
do j=1,size(field,2)
  weighted_field_local(:,j) = wts_lat(j)*field(:,j)
enddo

call mpp_global_field(grid_domain, weighted_field_local, weighted_field_global)
area_weighted_global_mean = sum(weighted_field_global)/(global_sum_of_wts*num_lon)

return
end function area_weighted_global_mean

!-------------------------------------------------------------------------
subroutine transforms_end
!-------------------------------------------------------------------------

if(.not.module_is_initialized) return

deallocate(lon_boundaries_global, lat_boundaries_global)
call grid_fourier_end
call spherical_fourier_end
call spec_mpp_end
module_is_initialized = .false.

return
end subroutine transforms_end
!-------------------------------------------------------------------------

end module transforms_mod
