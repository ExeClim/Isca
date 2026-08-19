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

module horiz_interp_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Bruce Wyman </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Performs spatial interpolation between grids.
! </OVERVIEW>

! <DESCRIPTION>
!     This module can interpolate data from any logically rectangular grid
!     to any logically rectangular grid. Four interpolation schems are used here:
!     conservative, bilinear, bicubic and inverse of square distance weighted. 
!     The four interpolation schemes are implemented seperately in 
!     horiz_interp_conserver_mod, horiz_interp_blinear_mod, horiz_interp_bicubic_mod
!     and horiz_interp_spherical_mod. bicubic interpolation requires the source grid
!     is regular lon/lat grid. User can choose the interpolation method in the 
!     public interface horiz_interp_new through optional argument interp_method,
!     with acceptable value "conservative", "bilinear", "bicubic" and "spherical".
!     The default value is "conservative". There is an optional mask field for 
!     missing input data. An optional output mask field may be used in conjunction with
!     the input mask to show where output data exists.
! </DESCRIPTION>

!-----------------------------------------------------------------------
!
!        Performs spatial interpolation between grids.
!
!-----------------------------------------------------------------------

use fms_mod,                    only: write_version_number, fms_error_handler
use mpp_mod,                    only: mpp_error, FATAL, stdout, mpp_min
use constants_mod,              only: pi
use horiz_interp_type_mod,      only: horiz_interp_type, assignment(=)
use horiz_interp_type_mod,      only: CONSERVE, BILINEAR, SPHERICA, BICUBIC
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_init, horiz_interp_conserve
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_new, horiz_interp_conserve_del
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_init, horiz_interp_bilinear
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_new, horiz_interp_bilinear_del
use horiz_interp_bicubic_mod,   only: horiz_interp_bicubic_init, horiz_interp_bicubic
use horiz_interp_bicubic_mod,   only: horiz_interp_bicubic_new, horiz_interp_bicubic_del
use horiz_interp_spherical_mod, only: horiz_interp_spherical_init, horiz_interp_spherical
use horiz_interp_spherical_mod, only: horiz_interp_spherical_new, horiz_interp_spherical_del

 implicit none
 private

!---- interfaces ----

 public   horiz_interp_type, horiz_interp, horiz_interp_new, horiz_interp_del, &
          horiz_interp_init, horiz_interp_end, assignment(=)

! <INTERFACE NAME="horiz_interp_new">
!   <OVERVIEW>
!      Allocates space and initializes a derived-type variable
!      that contains pre-computed interpolation indices and weights.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Allocates space and initializes a derived-type variable
!      that contains pre-computed interpolation indices and weights
!      for improved performance of multiple interpolations between
!      the same grids. This routine does not need to be called if you
!      are doing a single grid-to-grid interpolation.
!   </DESCRIPTION>
!   <IN NAME="lon_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
!      Longitude (in radians) for source data grid. You can pass 1-D lon_in to 
!      represent the geographical longitude of regular lon/lat grid, or just 
!      pass geographical longitude(lon_in is 2-D). The grid location may be 
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="lat_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
!      Latitude (in radians) for source data grid. You can pass 1-D lat_in to 
!      represent the geographical latitude of regular lon/lat grid, or just 
!      pass geographical latitude(lat_in is 2-D). The grid location may be 
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="lon_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Longitude (in radians) for destination data grid. You can pass 1-D lon_out to 
!      represent the geographical longitude of regular lon/lat grid, or just 
!      pass geographical longitude(lon_out is 2-D). The grid location may be 
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="lat_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Latitude (in radians) for destination data grid. You can pass 1-D lat_out to 
!      represent the geographical latitude of regular lon/lat grid, or just 
!      pass geographical latitude(lat_out is 2-D). The grid location may be 
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="verbose" TYPE="integer">
!      Integer flag that controls the amount of printed output.
!      verbose = 0, no output; = 1, min,max,means; = 2, still more
!   </IN>
!   <IN NAME="interp_method" TYPE="character(len=*)" > 
!      interpolation method, = "conservative", using conservation scheme,
!      = "bilinear", using bilinear interpolation, = "spherical",using spherical regrid.
!      = "bicubic", using bicubic interpolation. The default value is "convervative".
!   </IN>
!   <IN NAME = "src_modulo" >
!      Indicate the source data grid is cyclic or not.
!   </IN>
!   <IN NAME = "grid_at_center" >
!      Indicate the data is on the center of grid box or the edge of grid box. 
!      When true, the data is on the center of grid box. default vaule is false.
!      This option is only available when interp_method = "bilinear" or "bicubic".
!   </IN>
!   <OUT NAME="Interp" >
!      A derived-type variable containing indices and weights used for subsequent 
!      interpolations. To reinitialize this variable for a different grid-to-grid 
!      interpolation you must first use the "horiz_interp_del" interface.
!   </OUT>

 interface horiz_interp_new
    module procedure horiz_interp_new_1d     ! Source grid is 1d, destination grid is 1d
    module procedure horiz_interp_new_1d_src ! Source grid is 1d, destination grid is 2d
    module procedure horiz_interp_new_2d     ! Source grid is 2d, destination grid is 2d
    module procedure horiz_interp_new_1d_dst ! Source grid is 2d, destination grid is 1d
 end interface
! </INTERFACE>

! <INTERFACE NAME="horiz_interp">
!
!   <OVERVIEW>
!     Subroutine for performing the horizontal interpolation between two grids.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine for performing the horizontal interpolation between
!     two grids. There are two forms of this interface.
!     Form A requires first calling horiz_interp_new, while Form B
!     requires no initialization.
!   </DESCRIPTION>

!   <IN NAME="Interp" >
!     Derived-type variable containing interpolation indices and weights.
!     Returned by a previous call to horiz_interp_new.
!   </IN>
!   <IN NAME="data_in">
!      Input data on source grid.
!   </IN>
!   <IN NAME="verbose">
!      flag for the amount of print output.
!               verbose = 0, no output; = 1, min,max,means; = 2, still more
!   </IN>
!   <IN NAME="mask_in">
!      Input mask, must be the same size as the input data. The real value of
!      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points 
!      that should not be used or have missing data. It is Not needed for 
!      spherical regrid.
!   </IN>
!   <IN NAME="missing_value" >
!      Use the missing_value to indicate missing data.
!   </IN>
!   <IN NAME="missing_permit">
!      numbers of points allowed to miss for the bilinear interpolation. The value
!      should be between 0 and 3.
!   </IN>
!   <IN NAME="lon_in, lat_in" >
!      longitude and latitude (in radians) of source grid. More explanation can 
!      be found in the documentation of horiz_interp_new.
!   </IN>
!   <IN NAME="lon_out, lat_out" >
!      longitude and latitude (in radians) of destination grid. More explanation can 
!      be found in the documentation of horiz_interp_new.
!   </IN>
!   <OUT NAME="data_out">
!      Output data on destination grid.
!   </OUT>
!   <OUT NAME="mask_out">
!      Output mask that specifies whether data was computed.
!   </OUT>

!   <ERROR MSG="size of input array incorrect" STATUS="FATAL">
!      The input data array does not match the size of the input grid edges
!      specified. If you are using the initialization interface make sure you
!      have the correct grid size.
!   </ERROR>
!   <ERROR MSG="size of output array incorrect" STATUS="FATAL">
!      The output data array does not match the size of the input grid
!      edges specified. If you are using the initialization interface make
!      sure you have the correct grid size.
!   </ERROR>

 interface horiz_interp
    module procedure horiz_interp_base_2d
    module procedure horiz_interp_base_3d
    module procedure horiz_interp_solo_1d
    module procedure horiz_interp_solo_1d_src
    module procedure horiz_interp_solo_2d
    module procedure horiz_interp_solo_1d_dst
    module procedure horiz_interp_solo_old
 end interface
! </INTERFACE>

!-----------------------------------------------------------------------
 character(len=128) :: version = '$Id: horiz_interp.F90,v 19.0.4.1 2012/08/23 18:33:29 Zhi.Liang Exp $'
 character(len=128) :: tagname = '$Name: siena_201211 $'
 logical            :: module_is_initialized = .FALSE.
!-----------------------------------------------------------------------

contains

!#######################################################################
!  <SUBROUTINE NAME="horiz_interp_init">
!  <OVERVIEW>
!     writes version number and tag name to logfile.out
!  </OVERVIEW>
!  <DESCRIPTION>       
!     writes version number and tag name to logfile.out
!  </DESCRIPTION>

  subroutine horiz_interp_init

  if(module_is_initialized) return
  call write_version_number (version, tagname)
  call horiz_interp_conserve_init
  call horiz_interp_bilinear_init
  call horiz_interp_bicubic_init
  call horiz_interp_spherical_init

  module_is_initialized = .true.

  end subroutine horiz_interp_init

!  </SUBROUTINE>

!#######################################################################
!  <SUBROUTINE NAME="horiz_interp_new_1d" INTERFACE="horiz_interp_new">
!  <IN NAME="lon_in" TYPE="real" DIM="(:),(:,:)" UNITS="radians"></IN>
!  <IN NAME="lat_in" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="lon_out" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="lat_out" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="verbose" TYPE="integer, optional"></IN>
!  <IN NAME="interp_method" TYPE="character(len=*),optional"></IN>
!  <IN NAME="src_modulo" TYPE="logical, optional" > </IN>
!  <OUT NAME="Interp" TYPE="type(horiz_interp_type)"></OUT>

!<PUBLICROUTINE INTERFACE="horiz_interp_new">
  subroutine horiz_interp_new_1d (Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                                  interp_method, num_nbrs, max_dist, src_modulo,     &
                                  grid_at_center, mask_in, mask_out)
!</PUBLICROUTINE>

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout)        :: Interp
    real, intent(in),  dimension(:)               :: lon_in , lat_in
    real, intent(in),  dimension(:)               :: lon_out, lat_out
    integer, intent(in),                 optional :: verbose
    character(len=*), intent(in),        optional :: interp_method
    integer, intent(in),                 optional :: num_nbrs
    real,    intent(in),                 optional :: max_dist
    logical, intent(in),                 optional :: src_modulo
    logical, intent(in),                 optional :: grid_at_center
    real, intent(in), dimension(:,:),    optional :: mask_in  ! dummy
    real, intent(out),dimension(:,:),    optional :: mask_out ! dummy
    !-----------------------------------------------------------------------
    real, dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
    real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d
    integer                           :: i, j, nlon_in, nlat_in, nlon_out, nlat_out
    logical                           :: center
    character(len=40)                 :: method
    !-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method

    select case (trim(method))
    case ("conservative")
       Interp%interp_method = CONSERVE
       call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose)
    case ("bilinear")
       Interp%interp_method = BILINEAR
       center = .false.
       if(present(grid_at_center) ) center = grid_at_center
       if(center) then
          nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_out
             lon_dst(i,:) = lon_out(i)
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = lat_out(j)
          enddo

          call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_dst, lat_dst)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
          enddo
          do i = 1, nlon_out
             lon_dst(i,:) = (lon_out(i) + lon_out(i+1)) * 0.5
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = (lat_out(j) + lat_out(j+1)) * 0.5
          enddo
          call horiz_interp_bilinear_new ( Interp, lon_src_1d, lat_src_1d, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst, lat_dst)
       endif
    case ("bicubic")
       Interp%interp_method = BICUBIC
       center = .false.
       if(present(grid_at_center) ) center = grid_at_center
       !No need to expand to 2d, horiz_interp_bicubic_new does 1d-1d
       if(center) then 
          call horiz_interp_bicubic_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
            verbose, src_modulo)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst_1d(nlon_out), lat_dst_1d(nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
          enddo
          do i = 1, nlon_out
             lon_dst_1d(i) = (lon_out(i) + lon_out(i+1)) * 0.5
          enddo
          do j = 1, nlat_out
             lat_dst_1d(j) = (lat_out(j) + lat_out(j+1)) * 0.5
          enddo
          call horiz_interp_bicubic_new ( Interp, lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d)
       endif
    case ("spherical")
       Interp%interp_method = SPHERICA
       nlon_in  = size(lon_in(:));   nlat_in  = size(lat_in(:))
       nlon_out  = size(lon_out(:)); nlat_out = size(lat_out(:))
       allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
       allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
       do i = 1, nlon_in
          lon_src(i,:) = lon_in(i)
       enddo
       do j = 1, nlat_in
          lat_src(:,j) = lat_in(j)
       enddo
       do i = 1, nlon_out
          lon_dst(i,:) = lon_out(i)
       enddo
       do j = 1, nlat_out
          lat_dst(:,j) = lat_out(j)
       enddo
       call horiz_interp_spherical_new ( Interp, lon_src, lat_src, lon_dst, lat_dst, &
            num_nbrs, max_dist, src_modulo)
       deallocate(lon_src, lat_src, lon_dst, lat_dst)
    case default
       call mpp_error(FATAL,'horiz_interp_mod: interp_method should be conservative, bilinear, bicubic, spherical')
    end select

    !-----------------------------------------------------------------------
    Interp%I_am_initialized = .true.

  end subroutine horiz_interp_new_1d
!  </SUBROUTINE>

!#######################################################################

 subroutine horiz_interp_new_1d_src (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                     verbose, interp_method, num_nbrs, max_dist, &
                                     src_modulo, grid_at_center, mask_in, mask_out, is_latlon_out )

   type(horiz_interp_type), intent(inout)        :: Interp
   real, intent(in),  dimension(:)               :: lon_in , lat_in
   real, intent(in),  dimension(:,:)             :: lon_out, lat_out
   integer, intent(in),                 optional :: verbose
   character(len=*), intent(in),        optional :: interp_method
   integer, intent(in),                 optional :: num_nbrs  ! minimum number of neighbors
   real,    intent(in),                 optional :: max_dist
   logical, intent(in),                 optional :: src_modulo
   logical, intent(in),                 optional :: grid_at_center
   real, intent(in), dimension(:,:),    optional :: mask_in
   real, intent(out),dimension(:,:),    optional :: mask_out
   logical, intent(in),                 optional :: is_latlon_out

   real, dimension(:,:), allocatable :: lon_src, lat_src
   real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d
   integer                           :: i, j, nlon_in, nlat_in
   character(len=40)                 :: method
   logical                           :: center
   logical                           :: dst_is_latlon
   !-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'conservative'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      !--- check to see if the destination grid is regular lat-lon grid or not.
      if(PRESENT(is_latlon_out)) then
         dst_is_latlon = is_latlon_out
      else
         dst_is_latlon = is_lat_lon(lon_out, lat_out) 
      end if
      if(dst_is_latlon ) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -.0001) .or. ANY(mask_in > 1.0001)  ) call mpp_error(FATAL, &
                  'horiz_interp_conserve_new_1d_src(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), &
              verbose=verbose )
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if
   case ("conserve_great_circle")
      call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out, clip_method="conserve_great_circle" )
   case ("bilinear")
      Interp%interp_method = BILINEAR
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
         call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
         enddo
         call horiz_interp_bilinear_new ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("bicubic")
      Interp%interp_method = BICUBIC
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
        call horiz_interp_bicubic_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
         enddo
           call horiz_interp_bicubic_new ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("spherical")
      Interp%interp_method = SPHERICA
      nlon_in  = size(lon_in(:));  nlat_in  = size(lat_in(:))
      allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
      do i = 1, nlon_in
         lon_src(i,:) = lon_in(i)
      enddo
      do j = 1, nlat_in
         lat_src(:,j) = lat_in(j)
      enddo
      call horiz_interp_spherical_new ( Interp, lon_src, lat_src, lon_out, lat_out, &
           num_nbrs, max_dist, src_modulo)
      deallocate(lon_src, lat_src)
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical, conserve_great_circle')
   end select

   !-----------------------------------------------------------------------
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_1d_src

!#######################################################################

 subroutine horiz_interp_new_2d (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                 verbose, interp_method, num_nbrs, max_dist, &
                                 src_modulo, mask_in, mask_out, is_latlon_in, is_latlon_out  )
 type(horiz_interp_type), intent(inout)     :: Interp
 real, intent(in),  dimension(:,:)          :: lon_in , lat_in
 real, intent(in),  dimension(:,:)          :: lon_out, lat_out
 integer, intent(in),              optional :: verbose
 character(len=*), intent(in),     optional :: interp_method
 integer, intent(in),              optional :: num_nbrs
 real,    intent(in),              optional :: max_dist
 logical, intent(in),              optional :: src_modulo
 real, intent(in), dimension(:,:), optional :: mask_in
 real, intent(out),dimension(:,:), optional :: mask_out
 logical, intent(in),              optional :: is_latlon_in, is_latlon_out
 logical           :: src_is_latlon, dst_is_latlon
 character(len=40) :: method
!-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      if(PRESENT(is_latlon_in)) then
         src_is_latlon = is_latlon_in
      else
         src_is_latlon = is_lat_lon(lon_in, lat_in)
      end if
      if(PRESENT(is_latlon_out)) then
         dst_is_latlon = is_latlon_out
      else
         dst_is_latlon = is_lat_lon(lon_out, lat_out)
      end if
      if(src_is_latlon .AND. dst_is_latlon) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -.0001) .or. ANY(mask_in > 1.0001)  ) call mpp_error(FATAL, &
              'horiz_interp_conserve_new_2d(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out(:,1), lat_out(1,:), &
              verbose=verbose )
      else if(src_is_latlon) then
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      else if(dst_is_latlon) then
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      else 
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if

   case ("spherical") 
      Interp%interp_method = SPHERICA
      call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                    num_nbrs, max_dist, src_modulo )
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                        verbose, src_modulo )
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select     

!-----------------------------------------------------------------------
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_2d

!#######################################################################
 subroutine horiz_interp_new_1d_dst (Interp, lon_in, lat_in, lon_out, lat_out,   &
      verbose, interp_method, num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in )
   type(horiz_interp_type), intent(inout)     :: Interp
   real, intent(in),  dimension(:,:)          :: lon_in , lat_in
   real, intent(in),  dimension(:)            :: lon_out, lat_out
   integer, intent(in),              optional :: verbose
   character(len=*), intent(in),     optional :: interp_method
   integer, intent(in),              optional :: num_nbrs
   real,    intent(in),              optional :: max_dist
   logical, intent(in),              optional :: src_modulo
   real, intent(in), dimension(:,:), optional :: mask_in
   real, intent(out),dimension(:,:), optional :: mask_out
   logical, intent(in),              optional :: is_latlon_in

   character(len=40) :: method
   !-------------some local variables-----------------------------------------------
   integer                           :: i, j, nlon_out, nlat_out
   real, dimension(:,:), allocatable :: lon_dst, lat_dst
   logical                           :: src_is_latlon
   !-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
   allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
   do i = 1, nlon_out
      lon_dst(i,:) = lon_out(i)
   enddo
   do j = 1, nlat_out
      lat_dst(:,j) = lat_out(j)
   enddo

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      if(PRESENT(is_latlon_in)) then
         src_is_latlon = is_latlon_in
      else
         src_is_latlon = is_lat_lon(lon_in, lat_in)
      end if      

      if(src_is_latlon) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -.0001) .or. ANY(mask_in > 1.0001)  ) call mpp_error(FATAL, &
              'horiz_interp_conserve_new_1d_dst(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, &
              verbose=verbose)
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           verbose, src_modulo )
   case ("spherical")
      Interp%interp_method = SPHERICA
      call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           num_nbrs, max_dist, src_modulo)
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select

   deallocate(lon_dst,lat_dst)

   !-----------------------------------------------------------------------
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_1d_dst

!#######################################################################
! <SUBROUTINE NAME="horiz_interp_base_2d" INTERFACE="horiz_interp">
!   <IN NAME="Interp" TYPE="type(horiz_interp_type)"> </IN>
!   <IN NAME="data_in" TYPE="real" DIM="(:,:),(:,:,:)"> </IN>
!   <IN NAME="lon_in, lat_in" TYPE="real" DIM="(:),(:,:)"> </IN>
!   <IN NAME="lon_out, lat_out" TYPE="real" DIM="(:),(:,:)"> </IN>
!   <IN NAME="missing_value" TYPE="integer, optional" > </IN>
!   <IN NAME="missing_permit" TYPE="integer,optional" > </IN>
!   <IN NAME="verbose" TYPE="integer,optional"> </IN>
!   <IN NAME="mask_in" TYPE="real,optional" DIM="(:,:),(:,:,:)"> </IN>
!   <OUT NAME="data_out" TYPE="real" DIM="(:,:),(:,:,:)"> </OUT>
!   <OUT NAME="mask_out" TYPE="real,optional" DIM="(:,:),(:,:,:)"> </OUT>

!<PUBLICROUTINE INTERFACE="horiz_interp"> 
 subroutine horiz_interp_base_2d ( Interp, data_in, data_out, verbose, &
                                   mask_in, mask_out, missing_value, missing_permit, err_msg )
!</PUBLICROUTINE>
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
      real, intent(in),                   optional :: missing_value
      integer, intent(in),                optional :: missing_permit
   character(len=*), intent(out),         optional :: err_msg
!-----------------------------------------------------------------------
   if(present(err_msg)) err_msg = ''
   if(.not.Interp%I_am_initialized) then
     if(fms_error_handler('horiz_interp','The horiz_interp_type variable is not initialized',err_msg)) return
   endif

   select case(Interp%interp_method)
   case(CONSERVE)
      call horiz_interp_conserve(Interp,data_in, data_out, verbose, mask_in, mask_out)
   case(BILINEAR)
      call horiz_interp_bilinear(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit )
   case(BICUBIC)
      call horiz_interp_bicubic(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit )
   case(SPHERICA)
      call horiz_interp_spherical(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value )
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical')
   end select

   return

 end subroutine horiz_interp_base_2d
! </SUBROUTINE>

!#######################################################################

 subroutine horiz_interp_base_3d ( Interp, data_in, data_out, verbose, mask_in, mask_out, &
      missing_value, missing_permit, err_msg  )
   !-----------------------------------------------------------------------
   !   overload of interface horiz_interp_base_2d
   !   uses 3d arrays for data and mask
   !   this allows for multiple interpolations with one call
   !-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in)           :: Interp
   real, intent(in),  dimension(:,:,:)            :: data_in
   real, intent(out), dimension(:,:,:)            :: data_out
   integer, intent(in),                  optional :: verbose
   real, intent(in),   dimension(:,:,:), optional :: mask_in
   real, intent(out),  dimension(:,:,:), optional :: mask_out
   real, intent(in),                     optional :: missing_value
   integer, intent(in),                  optional :: missing_permit
   character(len=*), intent(out),        optional :: err_msg
   !-----------------------------------------------------------------------
   integer :: n

   if(present(err_msg)) err_msg = ''
   if(.not.Interp%I_am_initialized) then          
     if(fms_error_handler('horiz_interp','The horiz_interp_type variable is not initialized',err_msg)) return
   endif

   do n = 1, size(data_in,3)
      if (present(mask_in))then
         if(present(mask_out)) then
            call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), mask_out(:,:,n), &
                 missing_value, missing_permit )
         else
            call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
      else
         if(present(mask_out)) then
            call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_out=mask_out(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         else
            call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
     endif
   enddo
  
   return
!-----------------------------------------------------------------------
 end subroutine horiz_interp_base_3d

!#######################################################################
!<PUBLICROUTINE INTERFACE="horiz_interp"> 
 subroutine horiz_interp_solo_1d ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                   data_out, verbose, mask_in, mask_out,         &
                                   interp_method, missing_value, missing_permit, &
                                   num_nbrs, max_dist,src_modulo, grid_at_center  )              
!</PUBLICROUTINE>
!-----------------------------------------------------------------------
!   interpolates from a rectangular grid to rectangular grid.
!   interp_method can be the value conservative, bilinear or spherical.
!   horiz_interp_new don't need to be called before calling this routine.

!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:)   :: lon_in , lat_in
      real, intent(in),  dimension(:)   :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real, intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real, intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------
    call horiz_interp_init

    call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo, grid_at_center )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_del ( Interp )
!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d

!#######################################################################

 subroutine horiz_interp_solo_1d_src ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                       data_out, verbose, mask_in, mask_out,         &
                                       interp_method, missing_value, missing_permit, &
                                       num_nbrs, max_dist, src_modulo, grid_at_center )
!-----------------------------------------------------------------------
!
!   interpolates from a uniformly spaced grid to any output grid.
!   interp_method can be the value "onservative","bilinear" or "spherical".
!   horiz_interp_new don't need to be called before calling this routine.
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:)   :: lon_in , lat_in
      real, intent(in),  dimension(:,:) :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real, intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real, intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center

!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: dst_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init
    method = 'conservative'
    if(present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    if(trim(method) == 'conservative') dst_is_latlon = is_lat_lon(lon_out, lat_out)

    if(dst_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               grid_at_center, is_latlon_out = dst_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               grid_at_center, mask_in, mask_out, is_latlon_out = dst_is_latlon)

       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_src


!#######################################################################

 subroutine horiz_interp_solo_2d ( data_in, lon_in, lat_in, lon_out, lat_out, data_out, &
                                   verbose, mask_in, mask_out, interp_method, missing_value,&
                                   missing_permit, num_nbrs, max_dist, src_modulo  )
!-----------------------------------------------------------------------
!
!   interpolates from any grid to any grid. interp_method should be "spherical"
!   horiz_interp_new don't need to be called before calling this routine.
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:,:) :: lon_in , lat_in
      real, intent(in),  dimension(:,:) :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real, intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real, intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: dst_is_latlon, src_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    src_is_latlon = .true.
    if(trim(method) == 'conservative') then
       dst_is_latlon = is_lat_lon(lon_out, lat_out)
       src_is_latlon = is_lat_lon(lon_in, lat_in)
    end if

    if(dst_is_latlon .and. src_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               is_latlon_in=dst_is_latlon, is_latlon_out = dst_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               mask_in, mask_out, &
                               is_latlon_in=dst_is_latlon, is_latlon_out = dst_is_latlon)
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_2d

!#######################################################################

 subroutine horiz_interp_solo_1d_dst ( data_in, lon_in, lat_in, lon_out, lat_out, data_out,    &
                                       verbose, mask_in, mask_out,interp_method,missing_value, &
                                       missing_permit,  num_nbrs, max_dist, src_modulo)
!-----------------------------------------------------------------------
!
!   interpolates from any grid to rectangular longitude/latitude grid. 
!   interp_method should be "spherical".
!   horiz_interp_new don't need to be called before calling this routine.
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in),  dimension(:,:) :: lon_in , lat_in
      real, intent(in),  dimension(:)   :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose  
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real, intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real, intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: src_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method
    src_is_latlon = .true.
    if(trim(method) == 'conservative') src_is_latlon = is_lat_lon(lon_in, lat_in)

    if(src_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               is_latlon_in = src_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               mask_in, mask_out, is_latlon_in = src_is_latlon)

       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_dst

!#######################################################################

 subroutine horiz_interp_solo_old (data_in, wb, sb, dx, dy,  &
                                   lon_out, lat_out, data_out,  &
                                   verbose, mask_in, mask_out)

!-----------------------------------------------------------------------
!       Overloaded version of interface horiz_interp_solo_2
!
! input
!
!   data_in     Global input data stored from west to east (first dimension),
!               south to north (second dimension).  [real, dimension(:,:)]
!
!   wb          Longitude (in radians) that corresponds to western-most
!               boundary of grid box i=1 in array data_in.  [real]
!
!   sb          Latitude (in radians) that corresponds to southern-most
!               boundary of grid box j=1 in array data_in.  [real]
!
!   dx          Grid spacing (in radians) for the longitude axis (first
!               dimension) for the input data.  [real]
!
!   dy          Grid spacing (in radians) for the latitude axis (second
!               dimension) for the input data.  [real]
!
!   lon_out    The longitude edges (in radians) for output data grid boxes.
!               The values are for adjacent grid boxes and must increase in
!               value. If there are MLON grid boxes there must be MLON+1
!               edge values.  [real, dimension(:)]
!
!   lat_out    The latitude edges (in radians) for output data grid boxes.
!               The values are for adjacent grid boxes and may increase or
!               decrease in value. If there are NLAT grid boxes there must
!               be NLAT+1 edge values.  [real, dimension(:)]
!
! OUTPUT
!   data_out    Output data on the output grid defined by grid box
!               edges: blon_out and blat_out.  [real, dimension(:,:)]
!
!-----------------------------------------------------------------------
      real, intent(in),  dimension(:,:) :: data_in
      real, intent(in)                  :: wb, sb, dx, dy
      real, intent(in),  dimension(:)   :: lon_out, lat_out
      real, intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real, intent(in),   dimension(:,:), optional :: mask_in
      real, intent(out),  dimension(:,:), optional :: mask_out
!-----------------------------------------------------------------------
     real, dimension(size(data_in,1)+1)  :: blon_in
     real, dimension(size(data_in,2)+1)  :: blat_in
     integer :: i, j, nlon_in, nlat_in
     real    :: tpi
!-----------------------------------------------------------------------
   call horiz_interp_init
 
   tpi = 2.*pi
   nlon_in = size(data_in,1)
   nlat_in = size(data_in,2)

   do i = 1, nlon_in+1
      blon_in(i) = wb + float(i-1)*dx
   enddo
      if (abs(blon_in(nlon_in+1)-blon_in(1)-tpi) < epsilon(blon_in)) &
              blon_in(nlon_in+1)=blon_in(1)+tpi

   do j = 2, nlat_in
      blat_in(j) = sb + float(j-1)*dy
   enddo
      blat_in(1)         = -0.5*pi
      blat_in(nlat_in+1) =  0.5*pi


   call horiz_interp_solo_1d (data_in, blon_in, blat_in,    &
                              lon_out, lat_out, data_out,   &
                              verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_old

!#######################################################################
! <SUBROUTINE NAME="horiz_interp_del">

!   <OVERVIEW>
!     Deallocates memory used by "horiz_interp_type" variables.
!       Must be called before reinitializing with horiz_interp_new.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Deallocates memory used by "horiz_interp_type" variables.
!     Must be called before reinitializing with horiz_interp_new.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call horiz_interp_del ( Interp )
!   </TEMPLATE>

!   <INOUT NAME="Interp" TYPE="horiz_interp_type">
!     A derived-type variable returned by previous call
!              to horiz_interp_new. The input variable must have
!              allocated arrays. The returned variable will contain
!              deallocated arrays.
!   </INOUT>

! </SUBROUTINE>

 subroutine horiz_interp_del ( Interp )

   type (horiz_interp_type), intent(inout) :: Interp

!-----------------------------------------------------------------------
!  releases space used by horiz_interp_type variables
!  must be called before re-initializing the same variable
!-----------------------------------------------------------------------
   select case(Interp % interp_method) 
   case (CONSERVE)
      call horiz_interp_conserve_del(Interp )
   case (BILINEAR)
      call horiz_interp_bilinear_del(Interp )
   case (BICUBIC)
      call horiz_interp_bicubic_del(Interp )
   case (SPHERICA)
      call horiz_interp_spherical_del(Interp )
   end select

   Interp%I_am_initialized = .false.
!-----------------------------------------------------------------------

 end subroutine horiz_interp_del

 !#####################################################################

! <SUBROUTINE NAME="horiz_interp_end">

!   <OVERVIEW>
!     Dummy routine.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Dummy routine.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call horiz_interp_end
!   </TEMPLATE>

! </SUBROUTINE>

 subroutine horiz_interp_end
 return
 end subroutine horiz_interp_end

 !####################################################################
 function is_lat_lon(lon, lat)
    real, dimension(:,:), intent(in) :: lon, lat
    logical                          :: is_lat_lon
    integer                          :: i, j, nlon, nlat, num

    is_lat_lon = .true.
    nlon = size(lon,1)
    nlat = size(lon,2)
    LOOP_LAT: do j = 1, nlat
       do i = 2, nlon
          if(lat(i,j) .NE. lat(1,j)) then
             is_lat_lon = .false.
             exit LOOP_LAT
          end if
       end do
    end do LOOP_LAT

    if(is_lat_lon) then
       LOOP_LON: do i = 1, nlon
          do j = 2, nlat
             if(lon(i,j) .NE. lon(i,1)) then
                is_lat_lon = .false.
                exit LOOP_LON
             end if
          end do
       end do LOOP_LON
    end if

    num = 0
    if(is_lat_lon) num = 1
    call mpp_min(num)
    if(num == 1) then
       is_lat_lon = .true.
    else
       is_lat_lon = .false.
    end if

    return
 end function is_lat_lon

!#####################################################################

end module horiz_interp_mod

! <INFO>
!   <NOTE>             
!       Has not been checked with grids that do not cover the sphere.
!
!       Has not been checked with the optional mask arguments.
!
!       If a latitude or longitude index cannot be found the tolerance
!       used for making this determination may need to be increased.
!       This can be done by increasing the value of module variable
!       num_iters (default 4).
!   </NOTE>
!   <TESTPROGRAM>  
!     <PRE>
!       program test
!       use horiz_interp_mod
!       implicit none
!       integer, parameter :: nxi=177, nyi=91, nxo=133, nyo=77 ! resolution
!       real :: zi(nxi,nyi), zo(nxo,nyo)                       ! data
!       real :: xi(nxi+1), yi(nyi+1), xo(nxo+1), yo(nyo+1)     ! grid edges
!       real :: pi, tpi, hpi, dx, dy
!     
!       ! constants
!         hpi = acos(0.0)
!          pi = hpi*2.0
!         tpi = hpi*4.0
!     
!       ! grid setup: west to east, south to north
!         dx = tpi/real(nxi); call setaxis (0.,dx,xi);   xi(nxi+1) = xi(1)+tpi
!         dx = tpi/real(nxo); call setaxis (0.,dx,xo);   xo(nxo+1) = xo(1)+tpi
!         dy =  pi/real(nyi); call setaxis (-hpi,dy,yi); yi(nyi+1) = hpi
!         dy =  pi/real(nyo); call setaxis (-hpi,dy,yo); yo(nyo+1) = hpi
!     
!       ! random data on the input grid
!         call random_number (zi)
!     
!       ! interpolate (flipping y-axis)
!         call horiz_interp (zi(:,1:nyi:+1), xi, yi(1:nyi+1:+1), xo, yo(1:nyo+1:+1), zo, verbose=2)
!         call horiz_interp (zi(:,nyi:1:-1), xi, yi(nyi+1:1:-1), xo, yo(1:nyo+1:+1), zo, verbose=2)
!         call horiz_interp (zi(:,nyi:1:-1), xi, yi(nyi+1:1:-1), xo, yo(nyo+1:1:-1), zo, verbose=2)
!         call horiz_interp (zi(:,1:nyi:+1), xi, yi(1:nyi+1:+1), xo, yo(nyo+1:1:-1), zo, verbose=2)
!     
!       contains
!     ! set up a sequence of numbers
!         subroutine setaxis (xo,dx,x)
!         real, intent(in)  :: xo, dx
!         real, intent(out) :: x(:)
!         integer :: i
!           x(1) = xo
!           do i=2,size(x(:))
!             x(i) = x(i-1)+dx
!           enddo
!         end subroutine setaxis
!     
!       end program test
!     </PRE>
!   </TESTPROGRAM>
! </INFO>

#ifdef test_horiz_interp
! T More tests will be added in the future.
program horiz_interp_test

use mpp_mod,          only : mpp_init, mpp_exit, mpp_error, FATAL, stdout, mpp_npes
use mpp_mod,          only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,          only : mpp_pe, mpp_root_pe, NOTE, MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
use mpp_mod,          only : input_nml_file
use mpp_io_mod,       only : mpp_io_init, mpp_io_exit
use mpp_domains_mod,  only : mpp_define_layout, mpp_define_domains, mpp_get_compute_domain
use mpp_domains_mod,  only : mpp_domains_init, domain2d
use fms_mod,          only : file_exist, open_namelist_file, close_file, check_nml_error
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_new, horiz_interp_del
use horiz_interp_mod, only : horiz_interp, horiz_interp_type
use constants_mod,    only : constants_init, PI

implicit none

  integer :: ni_src = 360, nj_src = 180
  integer :: ni_dst = 144, nj_dst = 72

  namelist /test_horiz_interp_nml/ ni_src, nj_src, ni_dst, nj_dst

  real :: lon_src_beg = 0,    lon_src_end = 360
  real :: lat_src_beg = -90,  lat_src_end = 90
  real :: lon_dst_beg = -280, lon_dst_end = 80
  real :: lat_dst_beg = -90,  lat_dst_end = 90
  real :: D2R = PI/180.
  real, parameter :: SMALL = 1.0e-10

  type(domain2d)                    :: domain
  type(horiz_interp_type)           :: Interp
  integer                           :: id1, id2, id3, id4
  integer                           :: isc, iec, jsc, jec, i, j
  integer                           :: nml_unit, io, ierr, layout(2)
  real                              :: dlon_src, dlat_src, dlon_dst, dlat_dst
  real, allocatable, dimension(:)   :: lon1D_src, lat1D_src, lon1D_dst, lat1D_dst
  real, allocatable, dimension(:,:) :: lon2D_src, lat2D_src, lon2D_dst, lat2D_dst
  real, allocatable, dimension(:,:) :: data_src, data1_dst, data2_dst, data3_dst, data4_dst

  call constants_init
  call mpp_init
  call mpp_domains_init
  call mpp_io_init
  call horiz_interp_init

  !--- read namelist
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, test_horiz_interp_nml, iostat=io)
      ierr = check_nml_error(io, 'test_horiz_interp_nml')
#else
  if (file_exist('input.nml')) then
     ierr=1
     nml_unit = open_namelist_file()
     do while (ierr /= 0)
        read(nml_unit, nml=test_horiz_interp_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'test_horiz_interp_nml')
     enddo
10   call close_file(nml_unit)
  endif
#endif

  !--- define domains
  call mpp_define_layout( (/1, ni_dst, 1, nj_dst/), mpp_npes(), layout)
  call mpp_define_domains((/1, ni_dst, 1, nj_dst/), layout, domain)
  call mpp_get_compute_domain(domain,isc,iec,jsc,jec)

  !--- test conservative horiz_interp with a simple test. the source grid is the region
  !    (0:360,-90:90) with grid size ni_src, nj_src ( default 360X180). and the destination 
  !    is the region (-280:80, -90:90) with grid size ni_dstXnj_dst( default 144X72). 
  !    integer checksum and global sum will be printed out for both the 1D and 2D version. 

  allocate(lon2D_src(ni_src+1, nj_src+1), lat2D_src(ni_src+1, nj_src+1) )
  allocate(lon1D_src(ni_src+1), lat1D_src(nj_src+1), data_src(ni_src, nj_src) )

  allocate(lon2D_dst(isc:iec+1, jsc:jec+1), lat2D_dst(isc:iec+1, jsc:jec+1) )
  allocate(lon1D_dst(isc:iec+1), lat1D_dst(jsc:jec+1) )
  allocate(data1_dst(isc:iec, jsc:jec), data2_dst(isc:iec, jsc:jec) )
  allocate(data3_dst(isc:iec, jsc:jec), data4_dst(isc:iec, jsc:jec) )

  ! set up longitude and latitude of source/destination grid.   
  dlon_src = (lon_src_end-lon_src_beg)/ni_src 
  dlat_src = (lat_src_end-lat_src_beg)/nj_src
  dlon_dst = (lon_dst_end-lon_dst_beg)/ni_dst 
  dlat_dst = (lat_dst_end-lat_dst_beg)/nj_dst

  do i = 1, ni_src+1
     lon1D_src(i) = lon_src_beg + (i-1)*dlon_src
  end do

  do j = 1, nj_src+1
     lat1D_src(j) = lat_src_beg + (j-1)*dlat_src
  end do

  do i = isc, iec+1
     lon1D_dst(i) = lon_dst_beg + (i-1)*dlon_dst
  end do

  do j = jsc, jec+1
     lat1D_dst(j) = lat_dst_beg + (j-1)*dlat_dst
  end do

  ! scale grid to radians.
  lon1D_src = lon1D_src * D2R
  lat1D_src = lat1D_src * D2R
  lon1D_dst = lon1D_dst * D2R
  lat1D_dst = lat1D_dst * D2R

  do i = 1, ni_src+1
     lon2D_src(i,:) = lon1D_src(i)
  end do

  do j = 1, nj_src+1
     lat2D_src(:,j) = lat1D_src(j)
  end do

  do i = isc, iec+1
     lon2D_dst(i,:) = lon1D_dst(i)
  end do

  do j = jsc, jec+1
     lat2D_dst(:,j) = lat1D_dst(j)
  end do

  !--- set up the source data
  do j = 1, nj_src
     do i = 1, ni_src
        data_src(i,j) = i + j*0.001
     end do
  end do

  id1 = mpp_clock_id( 'horiz_interp_1dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  id2 = mpp_clock_id( 'horiz_interp_1dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  id3 = mpp_clock_id( 'horiz_interp_2dx1d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  id4 = mpp_clock_id( 'horiz_interp_2dx2d', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )

  ! --- 1dx1d version conservative interpolation
  call mpp_clock_begin(id1)
  call horiz_interp_new(Interp, lon1D_src, lat1D_src, lon1D_dst, lat1D_dst, interp_method = "conservative")
  call horiz_interp(Interp, data_src, data1_dst)
  call horiz_interp_del(Interp)
  call mpp_clock_end(id1)

  ! --- 1dx2d version conservative interpolation
  call mpp_clock_begin(id2)
  call horiz_interp_new(Interp, lon1D_src, lat1D_src, lon2D_dst, lat2D_dst, interp_method = "conservative")
  call horiz_interp(Interp, data_src, data2_dst)
  call horiz_interp_del(Interp)
  call mpp_clock_end(id2)

  ! --- 2dx1d version conservative interpolation
  call mpp_clock_begin(id3)
  call horiz_interp_new(Interp, lon2D_src, lat2D_src, lon1D_dst, lat1D_dst, interp_method = "conservative")
  call horiz_interp(Interp, data_src, data3_dst)
  call horiz_interp_del(Interp)
  call mpp_clock_end(id3)

  ! --- 2dx2d version conservative interpolation
  call mpp_clock_begin(id4)
  call horiz_interp_new(Interp, lon2D_src, lat2D_src, lon2D_dst, lat2D_dst, interp_method = "conservative")
  call horiz_interp(Interp, data_src, data4_dst)
  call horiz_interp_del(Interp)
  call mpp_clock_end(id4)

  !--- compare the data after interpolation between 1-D and 2-D version interpolation
  do j = jsc, jsc
     do i = isc, iec

        if( abs(data1_dst(i,j)-data2_dst(i,j)) > SMALL ) then
           print*, "After interpolation At point (i,j) = (", i, ",", j, "), data1 = ", data1_dst(i,j), &
           ", data2 = ", data2_dst(i,j), ", data1-data2 = ",  data1_dst(i,j) - data2_dst(i,j)
           call mpp_error(FATAL,"horiz_interp_test: data1_dst does not approxiamate data2_dst")
        end if
     end do
  end do

  if(mpp_pe() == mpp_root_pe()) call mpp_error(NOTE,   &
       "The test that verify 1dx2d version horiz_interp can reproduce 1dx1d version of horiz_interp is succesful") 

  do j = jsc, jsc
     do i = isc, iec

        if( abs(data1_dst(i,j)-data3_dst(i,j)) > SMALL ) then
           print*, "After interpolation At point (i,j) = (", i, ",", j, "), data1 = ", data1_dst(i,j), &
           ", data2 = ", data3_dst(i,j), ", data1-data2 = ",  data1_dst(i,j) - data3_dst(i,j)
           call mpp_error(FATAL,"horiz_interp_test: data1_dst does not approxiamate data3_dst")
        end if
     end do
  end do

  if(mpp_pe() == mpp_root_pe()) call mpp_error(NOTE,   &
       "The test that verify 2dx1d version horiz_interp can reproduce 1dx1d version of horiz_interp is succesful") 

  do j = jsc, jsc
     do i = isc, iec

        if( abs(data1_dst(i,j)-data4_dst(i,j)) > SMALL ) then
           print*, "After interpolation At point (i,j) = (", i, ",", j, "), data1 = ", data1_dst(i,j), &
           ", data2 = ", data4_dst(i,j), ", data1-data2 = ",  data1_dst(i,j) - data4_dst(i,j)
           call mpp_error(FATAL,"horiz_interp_test: data1_dst does not approxiamate data4_dst")
        end if
     end do
  end do

  if(mpp_pe() == mpp_root_pe()) call mpp_error(NOTE,   &
       "The test that verify 2dx2d version horiz_interp can reproduce 1dx1d version of horiz_interp is succesful") 

  call mpp_io_exit
  call mpp_exit

end program horiz_interp_test
#endif
