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

module horiz_interp_conserve_mod

  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Bruce Wyman </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

  ! <OVERVIEW>
  !   Performs spatial interpolation between grids using conservative interpolation
  ! </OVERVIEW>

  ! <DESCRIPTION>
  !     This module can conservatively interpolate data from any logically rectangular grid
  !     to any rectangular grid. The interpolation scheme is area-averaging 
  !     conservative scheme. There is an optional mask field for missing input data in both
  !     horiz_interp__conserveinit and horiz_interp_conserve. For efficiency purpose, mask should only be 
  !     kept in horiz_interp_init (will remove the mask in horiz_interp in the future). 
  !     There are 1-D and 2-D version of horiz_interp_conserve_init for 1-D and 2-D grid.
  !     There is a optional argument mask in horiz_interp_conserve_init_2d and no mask should 
  !     to passed into horiz_interp_conserv. optional argument mask will not be passed into
  !     horiz_interp_conserve_init_1d and optional argument mask may be passed into 
  !     horiz_interp_conserve (For the purpose of reproduce Memphis??? results).   
  !     An optional output mask field may be used in conjunction with the input mask to show 
  !     where output data exists.
  ! </DESCRIPTION>

  use mpp_mod,               only: mpp_send, mpp_recv, mpp_pe, mpp_root_pe, mpp_npes
  use mpp_mod,               only: mpp_error, FATAL,  mpp_sync_self 
  use mpp_mod,               only: COMM_TAG_1, COMM_TAG_2
  use fms_mod,               only: write_version_number
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type


  implicit none
  private

  ! public interface


  ! <INTERFACE NAME="horiz_interp_conserve_new">
  !   <OVERVIEW>
  !      Allocates space and initializes a derived-type variable
  !      that contains pre-computed interpolation indices and weights.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !      Allocates space and initializes a derived-type variable
  !      that contains pre-computed interpolation indices and weights
  !      for improved performance of multiple interpolations between
  !      the same grids. 

  !   </DESCRIPTION>
  !   <IN NAME="lon_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
  !      Longitude (in radians) for source data grid. 
  !   </IN>
  !   <IN NAME="lat_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
  !      Latitude (in radians) for source data grid.
  !   </IN>
  !   <IN NAME="lon_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
  !      Longitude (in radians) for destination data grid. 
  !   </IN>
  !   <IN NAME="lat_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
  !      Latitude (in radians) for destination data grid. 
  !   </IN>
  !   <IN NAME="verbose" TYPE="integer, optional" >
  !      flag for the amount of print output.
  !   </IN>
  !   <IN NAME="mask_in" TYPE="real, dimension(:,:),optional">
  !      Input mask.  must be the size (size(lon_in)-1, size(lon. The real value of
  !      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points 
  !      that should not be used or have missing data.
  !   </IN>
  !   <OUT NAME="mask_out" TYPE="real, dimension(:,:),optional">
  !      Output mask that specifies whether data was computed.
  !   </OUT>
  !   <INOUT NAME="Interp" TYPE="type(horiz_interp_type)" >
  !      A derived-type variable containing indices and weights used for subsequent 
  !      interpolations. To reinitialize this variable for a different grid-to-grid 
  !      interpolation you must first use the "horiz_interp_del" interface.
  !   </INOUT>
  interface horiz_interp_conserve_new
     module procedure horiz_interp_conserve_new_1dx1d
     module procedure horiz_interp_conserve_new_1dx2d
     module procedure horiz_interp_conserve_new_2dx1d
     module procedure horiz_interp_conserve_new_2dx2d
  end interface
  ! </INTERFACE>
  public :: horiz_interp_conserve_init 
  public :: horiz_interp_conserve_new, horiz_interp_conserve, horiz_interp_conserve_del

  integer :: pe, root_pe
  !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: horiz_interp_conserve.F90,v 19.0.2.2.2.1 2012/08/23 18:33:29 Zhi.Liang Exp $'
  character(len=128) :: tagname = '$Name: siena_201211 $'
  logical            :: module_is_initialized = .FALSE.

  real, parameter :: range_check_criteria = 0.05

contains

  !#######################################################################
  !  <SUBROUTINE NAME="horiz_interp_conserve_init">
  !  <OVERVIEW>
  !     writes version number and tag name to logfile.out
  !  </OVERVIEW>
  !  <DESCRIPTION>       
  !     writes version number and tag name to logfile.out
  !  </DESCRIPTION>

  subroutine horiz_interp_conserve_init

    if(module_is_initialized) return
    call write_version_number (version, tagname)
    module_is_initialized = .true.

  end subroutine horiz_interp_conserve_init

  !  </SUBROUTINE>

  !#######################################################################
  !<PUBLICROUTINE INTERFACE="horiz_interp_conserve_new">
  subroutine horiz_interp_conserve_new_1dx1d ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, clip_method)
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in),       dimension(:)   :: lon_in , lat_in
    real, intent(in),       dimension(:)   :: lon_out, lat_out
    integer, intent(in),       optional    :: verbose
    character(len=*), intent(in), optional :: clip_method
  !</PUBLICROUTINE>
    !-----------------------------------------------------------------------
    real, dimension(size(lat_out(:))-1,2) :: sph
    real, dimension(size(lon_out(:))-1,2) :: theta
    real, dimension(size(lat_in(:)))      :: slat_in
    real, dimension(size(lon_in(:))-1)    :: dlon_in
    real, dimension(size(lat_in(:))-1)    :: dsph_in
    real, dimension(size(lon_out(:))-1)   :: dlon_out
    real, dimension(size(lat_out(:))-1)   :: dsph_out
    real    :: blon, fac, hpi, tpi, eps
    integer :: num_iters = 4
    integer :: i, j, m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
         iverbose, m2, n2, iter
    logical :: s2n
    character(len=64) :: mesg

    if(PRESENT(clip_method)) call mpp_error(FATAL, &
         'horiz_interp_conserve_new_1dx1d: clip_methos is not implemented, contact developer')
    !-----------------------------------------------------------------------
    iverbose = 0;  if (present(verbose)) iverbose = verbose

    pe      = mpp_pe()
    root_pe = mpp_root_pe()
    !-----------------------------------------------------------------------
    hpi = 0.5*pi
    tpi = 4.*hpi
    Interp%version = 1
    nlon_in = size(lon_in(:))-1;  nlat_in = size(lat_in(:))-1
    nlon_out = size(lon_out(:))-1;  nlat_out = size(lat_out(:))-1

    allocate ( Interp % facj (nlat_out,2), Interp % jlat (nlat_out,2),      &
               Interp % faci (nlon_out,2), Interp % ilon (nlon_out,2),      &
               Interp % area_src (nlon_in, nlat_in),   &
               Interp % area_dst (nlon_out, nlat_out) )

    !-----------------------------------------------------------------------
    !  --- set-up for input grid boxes ---

    do j = 1, nlat_in+1
       slat_in(j) = sin(lat_in(j))
    enddo

    do j = 1, nlat_in
       dsph_in(j) = abs(slat_in(j+1)-slat_in(j))
    enddo

    do i = 1,nlon_in
       dlon_in(i) = abs(lon_in(i+1)-lon_in(i))
    enddo

    !  set south to north flag
    s2n = .true.
    if (lat_in(1) > lat_in(nlat_in+1)) s2n = .false.

    !-----------------------------------------------------------------------
    !  --- set-up for output grid boxes ---

    do n = 1, nlat_out
       dsph_out(n) = abs(sin(lat_out(n+1))-sin(lat_out(n)))
    enddo

    do m = 1,nlon_out
       theta(m,1)  = lon_out(m)
       theta(m,2)  = lon_out(m+1)
       dlon_out(m) = abs(lon_out(m+1)-lon_out(m))
    enddo

    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out
    !***********************************************************************

    !------ set up latitudinal indexing ------
    !------ make sure output grid goes south to north ------

    do n = 1, nlat_out
       if (lat_out(n) < lat_out(n+1)) then
          sph(n,1) = sin(lat_out(n))
          sph(n,2) = sin(lat_out(n+1))
       else
          sph(n,1) = sin(lat_out(n+1))
          sph(n,2) = sin(lat_out(n))
       endif
    enddo

    Interp%jlat = 0
    do n2 = 1, 2         ! looping on grid box edges
       do n = 1, nlat_out   ! looping on output latitudes
          eps = 0.0
          do iter=1,num_iters
             ! find indices from input latitudes
             do j = 1, nlat_in
                if ( (s2n .and. (slat_in(j)-sph(n,n2)) <= eps .and.   &
                     (sph(n,n2)-slat_in(j+1)) <= eps) .or. &
                     (.not.s2n .and. (slat_in(j+1)-sph(n,n2)) <= eps .and.  &
                     (sph(n,n2)-slat_in(j)) <= eps) ) then
                   Interp%jlat(n,n2) = j
                   ! weight with sin(lat) to exactly conserve area-integral
                   fac = (sph(n,n2)-slat_in(j))/(slat_in(j+1)-slat_in(j))
                   if (s2n) then
                      if (n2 == 1) Interp%facj(n,n2) = 1.0 - fac
                      if (n2 == 2) Interp%facj(n,n2) = fac
                   else
                      if (n2 == 1) Interp%facj(n,n2) = fac
                      if (n2 == 2) Interp%facj(n,n2) = 1.0 - fac
                   endif
                   exit
                endif
             enddo
             if ( Interp%jlat(n,n2) /= 0 ) exit
             ! did not find this output grid edge in the input grid
             ! increase tolerance for multiple passes
             eps  = epsilon(sph)*real(10**iter)
          enddo
          ! no match
          if ( Interp%jlat(n,n2) == 0 ) then
             write (mesg,710) n,sph(n,n2)
710          format (': n,sph=',i3,f14.7,40x)
             call mpp_error(FATAL, 'horiz_interp_conserve_mod:no latitude index found'//trim(mesg))
          endif
       enddo
    enddo

    !------ set up longitudinal indexing ------

    Interp%ilon = 0
    do m2 = 1, 2         ! looping on grid box edges
       do m = 1, nlon_out   ! looping on output longitudes
          blon = theta(m,m2)
          if ( blon < lon_in(1)         ) blon = blon + tpi
          if ( blon > lon_in(nlon_in+1) ) blon = blon - tpi
          eps = 0.0
          do iter=1,num_iters
             ! find indices from input longitudes
             do i = 1, nlon_in
                if ( (lon_in(i)-blon) <= eps .and. &
                     (blon-lon_in(i+1)) <= eps ) then
                   Interp%ilon(m,m2) = i
                   fac = (blon-lon_in(i))/(lon_in(i+1)-lon_in(i))
                   if (m2 == 1) Interp%faci(m,m2) = 1.0 - fac
                   if (m2 == 2) Interp%faci(m,m2) = fac
                   exit
                endif
             enddo
             if ( Interp%ilon(m,m2) /= 0 ) exit
             ! did not find this output grid edge in the input grid
             ! increase tolerance for multiple passes
             eps  = epsilon(blon)*real(10**iter)
          enddo
          ! no match
          if ( Interp%ilon(m,m2) == 0 ) then
             print *, 'lon_out,blon,blon_in,eps=',  &
                  theta(m,m2),blon,lon_in(1),lon_in(nlon_in+1),eps
             call mpp_error(FATAL, 'horiz_interp_conserve_mod: no longitude index found')
          endif
       enddo
    enddo

    !  --- area of input grid boxes ---

    do j = 1,nlat_in
       do i = 1,nlon_in
          Interp%area_src(i,j) = dlon_in(i) * dsph_in(j)
       enddo
    enddo

    !  --- area of output grid boxes ---

    do n = 1, nlat_out
       do m = 1, nlon_out
          Interp%area_dst(m,n) = dlon_out(m) * dsph_out(n)
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! this output may be quite lengthy and is not recommended
    ! when using more than one processor
    if (iverbose > 2) then
       write (*,801) (i,Interp%ilon(i,1),Interp%ilon(i,2),  &
            Interp%faci(i,1),Interp%faci(i,2),i=1,nlon_out)
       write (*,802) (j,Interp%jlat(j,1),Interp%jlat(j,2),  &
            Interp%facj(j,1),Interp%facj(j,2),j=1,nlat_out)
801    format (/,2x,'i',4x,'is',5x,'ie',4x,'facis',4x,'facie',  &
            /,(i4,2i7,2f10.5))
802    format (/,2x,'j',4x,'js',5x,'je',4x,'facjs',4x,'facje',  &
            /,(i4,2i7,2f10.5))
    endif
    !-----------------------------------------------------------------------

  end subroutine horiz_interp_conserve_new_1dx1d

  !#######################################################################
  !<PUBLICROUTINE INTERFACE="horiz_interp_conserve_new">
  subroutine horiz_interp_conserve_new_1dx2d ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                               mask_in, mask_out, verbose, clip_method)
    type(horiz_interp_type),        intent(inout) :: Interp
    real, intent(in),              dimension(:)   :: lon_in , lat_in
    real, intent(in),              dimension(:,:) :: lon_out, lat_out
    real, intent(in),    optional, dimension(:,:) :: mask_in
    real, intent(inout), optional, dimension(:,:) :: mask_out
    integer, intent(in), optional                 :: verbose
    character(len=*), intent(in), optional        :: clip_method
  !</PUBLICROUTINE>

    integer :: create_xgrid_1DX2D_order1, get_maxxgrid, maxxgrid
    integer :: create_xgrid_great_circle
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, nxgrid, i, j
    real, dimension(size(lon_in(:))-1, size(lat_in(:))-1) :: mask_src
    integer, allocatable, dimension(:)   :: i_src, j_src, i_dst, j_dst
    real,    allocatable, dimension(:)   :: xgrid_area, clon, clat
    real,    allocatable, dimension(:,:) :: dst_area, lon_src, lat_src
    character(len=64) :: my_clip_method

    if( (size(lon_out,1) .NE. size(lat_out,1)) .OR. (size(lon_out,2) .NE. size(lat_out,2)) )  &
        call mpp_error(FATAL, 'horiz_interp_conserve_mod: size mismatch between lon_out and lat_out')
    nlon_in  = size(lon_in(:)) - 1;  nlat_in  = size(lat_in(:)) - 1
    nlon_out = size(lon_out,1) - 1;  nlat_out = size(lon_out,2) - 1

    mask_src = 1
    if(present(mask_in)) then
       if( (size(mask_in,1) .NE. nlon_in) .OR.  (size(mask_in,2) .NE. nlat_in)) call mpp_error(FATAL, &
         'horiz_interp_conserve_mod: size mismatch between mask_in and lon_in/lat_in')
       mask_src = mask_in
    end if

    maxxgrid = get_maxxgrid()
    allocate(i_src(maxxgrid), j_src(maxxgrid), i_dst(maxxgrid), j_dst(maxxgrid) )
    allocate( xgrid_area(maxxgrid), dst_area(nlon_out, nlat_out) )

    my_clip_method = "conserve_latlon"
    if(present(clip_method)) then
       if(trim(clip_method) .NE. "conserve_latlon" .AND. trim(clip_method) .NE. "conserve_great_circle") &
          call mpp_error(FATAL,'horiz_interp_conserve_mod: clip_method='//trim(clip_method)//' is not a valid option')
       my_clip_method = clip_method
    endif
    if( trim(my_clip_method) == "conserve_latlon") then
       nxgrid = create_xgrid_1DX2D_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, &
                                          mask_src, i_src, j_src, i_dst, j_dst, xgrid_area)
    else
       allocate(lon_src(nlon_in+1,nlat_in+1), lat_src(nlon_in+1,nlat_in+1))
       allocate(clon(maxxgrid), clat(maxxgrid))       
       do j = 1, nlat_in+1
          do i = 1, nlon_in+1
             lon_src(i,j) = lon_in(i)
             lat_src(i,j) = lat_in(j)
          enddo
       enddo
       nxgrid =  create_xgrid_great_circle(nlon_in, nlat_in, nlon_out, nlat_out, lon_src, lat_src, lon_out, lat_out, &
                                          mask_src, i_src, j_src, i_dst, j_dst, xgrid_area, clon, clat, range_check_criteria)
       deallocate(lon_src, lat_src, clon, clat)
    endif
    allocate(Interp%i_src(nxgrid), Interp%j_src(nxgrid) )
    allocate(Interp%i_dst(nxgrid), Interp%j_dst(nxgrid) )
    allocate(Interp%area_frac_dst(nxgrid) )
    Interp%version = 2
    Interp%nxgrid   = nxgrid
    Interp%i_src = i_src(1:nxgrid)+1 ! in C, the starting index is 0
    Interp%j_src = j_src(1:nxgrid)+1
    Interp%i_dst = i_dst(1:nxgrid)+1
    Interp%j_dst = j_dst(1:nxgrid)+1

    ! sum over exchange grid area to get destination grid area
    dst_area = 0
    do i = 1, nxgrid
       dst_area(Interp%i_dst(i), Interp%j_dst(i)) = dst_area(Interp%i_dst(i), Interp%j_dst(i)) + xgrid_area(i)       
    end do    

    do i = 1, nxgrid
       Interp%area_frac_dst(i) = xgrid_area(i)/dst_area(Interp%i_dst(i), Interp%j_dst(i) )
    end do
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out
    if(present(mask_out)) then
       if( (size(mask_out,1) .NE. nlon_out) .OR. (size(mask_out,2) .NE. nlat_out) ) call mpp_error(FATAL, &
         'horiz_interp_conserve_mod: size mismatch between mask_out and lon_out/lat_out')
       mask_out = 0.0
       do i = 1, nxgrid
          mask_out(Interp%i_dst(i),Interp%j_dst(i)) = mask_out(Interp%i_dst(i),Interp%j_dst(i)) + Interp%area_frac_dst(i)
       end do
    end if

    deallocate(i_src, j_src, i_dst, j_dst, xgrid_area, dst_area )

  end subroutine horiz_interp_conserve_new_1dx2d

  !#######################################################################
  !<PUBLICROUTINE INTERFACE="horiz_interp_conserve_new">
  subroutine horiz_interp_conserve_new_2dx1d ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                               mask_in, mask_out, verbose, clip_method)
    type(horiz_interp_type),        intent(inout) :: Interp
    real, intent(in),              dimension(:,:) :: lon_in , lat_in
    real, intent(in),              dimension(:)   :: lon_out, lat_out
    real, intent(in),    optional, dimension(:,:) :: mask_in
    real, intent(inout), optional, dimension(:,:) :: mask_out
    integer, intent(in), optional                 :: verbose
    character(len=*), intent(in),        optional :: clip_method
  !</PUBLICROUTINE>

    integer :: create_xgrid_2DX1D_order1, get_maxxgrid, maxxgrid
    integer :: create_xgrid_great_circle
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, nxgrid, i, j
    real, dimension(size(lon_in,1)-1, size(lon_in,2)-1) :: mask_src
    integer, allocatable, dimension(:)   :: i_src, j_src, i_dst, j_dst
    real,    allocatable, dimension(:)   :: xgrid_area, clon, clat
    real,    allocatable, dimension(:,:) :: dst_area, lon_dst, lat_dst
    character(len=64) :: my_clip_method

    if( (size(lon_in,1) .NE. size(lat_in,1)) .OR. (size(lon_in,2) .NE. size(lat_in,2)) )  &
        call mpp_error(FATAL, 'horiz_interp_conserve_mod: size mismatch between lon_in and lat_in')
    nlon_in  = size(lon_in,1)   - 1;  nlat_in  = size(lon_in,2)   - 1
    nlon_out = size(lon_out(:)) - 1;  nlat_out = size(lat_out(:)) - 1

    mask_src = 1
    if(present(mask_in)) then
       if( (size(mask_in,1) .NE. nlon_in) .OR.  (size(mask_in,2) .NE. nlat_in)) call mpp_error(FATAL, &
         'horiz_interp_conserve_mod: size mismatch between mask_in and lon_in/lat_in')
       mask_src = mask_in
    end if

    maxxgrid = get_maxxgrid()
    allocate(i_src(maxxgrid), j_src(maxxgrid), i_dst(maxxgrid), j_dst(maxxgrid) )
    allocate( xgrid_area(maxxgrid), dst_area(nlon_out, nlat_out) )

    my_clip_method = "conserve_latlon"
    if(present(clip_method)) then
       if(trim(clip_method) .NE. "conserve_latlon" .AND. trim(clip_method) .NE. "conserve_great_circle") &
          call mpp_error(FATAL,'horiz_interp_conserve_mod: clip_method='//trim(clip_method)//' is not a valid option')
       my_clip_method = clip_method
    endif
    if( trim(my_clip_method) == "conserve_latlon") then    
       nxgrid = create_xgrid_2DX1D_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, &
                                       mask_src, i_src, j_src, i_dst, j_dst, xgrid_area)
    else
       allocate(lon_dst(nlon_out+1, nlat_out+1), lat_dst(nlon_out+1, nlat_out+1) )
       allocate(clon(maxxgrid), clat(maxxgrid))   
       do j = 1, nlat_out+1
          do i = 1, nlon_out+1
             lon_dst(i,j) = lon_out(i)
             lat_dst(i,j) = lat_out(j)
          enddo
       enddo
       nxgrid =  create_xgrid_great_circle(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_dst, lat_dst, &
                                          mask_src, i_src, j_src, i_dst, j_dst, xgrid_area, clon, clat, range_check_criteria)
       deallocate(lon_dst, lat_dst, clon, clat)
    endif
    allocate(Interp%i_src(nxgrid), Interp%j_src(nxgrid) )
    allocate(Interp%i_dst(nxgrid), Interp%j_dst(nxgrid) )
    allocate(Interp%area_frac_dst(nxgrid) )
    Interp%version = 2
    Interp%nxgrid   = nxgrid
    Interp%i_src = i_src(1:nxgrid)+1 ! in C, the starting index is 0
    Interp%j_src = j_src(1:nxgrid)+1
    Interp%i_dst = i_dst(1:nxgrid)+1
    Interp%j_dst = j_dst(1:nxgrid)+1

    ! sum over exchange grid area to get destination grid area
    dst_area = 0
    do i = 1, nxgrid
       dst_area(Interp%i_dst(i), Interp%j_dst(i)) = dst_area(Interp%i_dst(i), Interp%j_dst(i)) + xgrid_area(i)       
    end do    

    do i = 1, nxgrid
       Interp%area_frac_dst(i) = xgrid_area(i)/dst_area(Interp%i_dst(i), Interp%j_dst(i) )
    end do
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out
    if(present(mask_out)) then
       if( (size(mask_out,1) .NE. nlon_out) .OR. (size(mask_out,2) .NE. nlat_out) ) call mpp_error(FATAL, &
         'horiz_interp_conserve_mod: size mismatch between mask_out and lon_out/lat_out')
       mask_out = 0.0
       do i = 1, nxgrid
          mask_out(Interp%i_dst(i),Interp%j_dst(i)) = mask_out(Interp%i_dst(i),Interp%j_dst(i)) + Interp%area_frac_dst(i)
       end do
    end if

    deallocate(i_src, j_src, i_dst, j_dst, xgrid_area, dst_area)

  end subroutine horiz_interp_conserve_new_2dx1d

  !#######################################################################
  !<PUBLICROUTINE INTERFACE="horiz_interp_conserve_new">
  subroutine horiz_interp_conserve_new_2dx2d ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                               mask_in, mask_out, verbose, clip_method)
    type(horiz_interp_type),        intent(inout) :: Interp
    real, intent(in),              dimension(:,:) :: lon_in , lat_in
    real, intent(in),              dimension(:,:) :: lon_out, lat_out
    real, intent(in),    optional, dimension(:,:) :: mask_in
    real, intent(inout), optional, dimension(:,:) :: mask_out
    integer, intent(in), optional                 :: verbose
    character(len=*), intent(in),        optional :: clip_method
  !</PUBLICROUTINE>

    integer :: create_xgrid_2DX2D_order1, get_maxxgrid, maxxgrid
    integer :: create_xgrid_great_circle
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, nxgrid, i
    real, dimension(size(lon_in,1)-1, size(lon_in,2)-1) :: mask_src
    integer, allocatable, dimension(:)   :: i_src, j_src, i_dst, j_dst
    real,    allocatable, dimension(:)   :: xgrid_area, clon, clat
    real,    allocatable, dimension(:,:) :: dst_area
    character(len=64) :: my_clip_method

    if( (size(lon_in,1) .NE. size(lat_in,1)) .OR. (size(lon_in,2) .NE. size(lat_in,2)) )  &
        call mpp_error(FATAL, 'horiz_interp_conserve_mod: size mismatch between lon_in and lat_in')
    if( (size(lon_out,1) .NE. size(lat_out,1)) .OR. (size(lon_out,2) .NE. size(lat_out,2)) )  &
        call mpp_error(FATAL, 'horiz_interp_conserve_mod: size mismatch between lon_out and lat_out')
    nlon_in  = size(lon_in,1)  - 1;  nlat_in  = size(lon_in,2)  - 1
    nlon_out = size(lon_out,1) - 1;  nlat_out = size(lon_out,2) - 1

    mask_src = 1
    if(present(mask_in)) then
       if( (size(mask_in,1) .NE. nlon_in) .OR.  (size(mask_in,2) .NE. nlat_in)) call mpp_error(FATAL, &
         'horiz_interp_conserve_mod: size mismatch between mask_in and lon_in/lat_in')
       mask_src = mask_in
    end if

    maxxgrid = get_maxxgrid()
    allocate(i_src(maxxgrid), j_src(maxxgrid), i_dst(maxxgrid), j_dst(maxxgrid) )
    allocate( xgrid_area(maxxgrid), dst_area(nlon_out, nlat_out) )

    my_clip_method = "conserve_latlon"
    if(present(clip_method)) then
       if(trim(clip_method) .NE. "conserve_latlon" .AND. trim(clip_method) .NE. "conserve_great_circle") &
          call mpp_error(FATAL,'horiz_interp_conserve_mod: clip_method='//trim(clip_method)//' is not a valid option')
       my_clip_method = clip_method
    endif
    if( trim(my_clip_method) == "conserve_latlon") then   
       nxgrid = create_xgrid_2DX2D_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, &
                                       mask_src, i_src, j_src, i_dst, j_dst, xgrid_area) 
    else
       allocate(clon(maxxgrid), clat(maxxgrid))   
       nxgrid =  create_xgrid_great_circle(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, &
                                          mask_src, i_src, j_src, i_dst, j_dst, xgrid_area, clon, clat, range_check_criteria)
       deallocate(clon, clat)
    endif

    allocate(Interp%i_src(nxgrid), Interp%j_src(nxgrid) )
    allocate(Interp%i_dst(nxgrid), Interp%j_dst(nxgrid) )
    allocate(Interp%area_frac_dst(nxgrid) )
    Interp%version = 2
    Interp%nxgrid   = nxgrid
    Interp%i_src = i_src(1:nxgrid)+1 ! in C, the starting index is 0
    Interp%j_src = j_src(1:nxgrid)+1
    Interp%i_dst = i_dst(1:nxgrid)+1
    Interp%j_dst = j_dst(1:nxgrid)+1

    ! sum over exchange grid area to get destination grid area
    dst_area = 0
    do i = 1, nxgrid
       dst_area(Interp%i_dst(i), Interp%j_dst(i)) = dst_area(Interp%i_dst(i), Interp%j_dst(i)) + xgrid_area(i)       
    end do    

    do i = 1, nxgrid
       Interp%area_frac_dst(i) = xgrid_area(i)/dst_area(Interp%i_dst(i), Interp%j_dst(i) )
    end do

    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out
    if(present(mask_out)) then
       if( (size(mask_out,1) .NE. nlon_out) .OR. (size(mask_out,2) .NE. nlat_out) ) call mpp_error(FATAL, &
         'horiz_interp_conserve_mod: size mismatch between mask_out and lon_out/lat_out')
       mask_out = 0.0
       do i = 1, nxgrid
          mask_out(Interp%i_dst(i),Interp%j_dst(i)) = mask_out(Interp%i_dst(i),Interp%j_dst(i)) + Interp%area_frac_dst(i)
       end do
    end if

    deallocate(i_src, j_src, i_dst, j_dst, xgrid_area, dst_area )

  end subroutine horiz_interp_conserve_new_2dx2d

  !########################################################################
  ! <SUBROUTINE NAME="horiz_interp_conserve">

  !   <OVERVIEW>
  !      Subroutine for performing the horizontal interpolation between two grids.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Subroutine for performing the horizontal interpolation between two grids. 
  !     horiz_interp_conserve_new must be called before calling this routine.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_conserve ( Interp, data_in, data_out, verbose, mask_in, mask_out)
  !   </TEMPLATE>
  !   
  !   <IN NAME="Interp" TYPE="type(horiz_interp_type)">
  !     Derived-type variable containing interpolation indices and weights.
  !     Returned by a previous call to horiz_interp_new.
  !   </IN>
  !   <IN NAME="data_in" TYPE="real, dimension(:,:)">
  !      Input data on source grid.
  !   </IN>

  !   <IN NAME="verbose" TYPE="integer, optional">
  !      flag for the amount of print output.
  !               verbose = 0, no output; = 1, min,max,means; = 2, still more
  !   </IN>
  !   <IN NAME="mask_in" TYPE="real, dimension(:,:),optional">
  !      Input mask, must be the same size as the input data. The real value of
  !      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points 
  !      that should not be used or have missing data. mask_in will be applied only
  !      when horiz_interp_conserve_new_1d is called. mask_in will be passed into
  !      horiz_interp_conserve_new_2d.
  !   </IN>

  !   <OUT NAME="data_out" TYPE="real, dimension(:,:)">
  !      Output data on destination grid.
  !   </OUT>
  !   <OUT NAME="mask_out" TYPE="real, dimension(:,:),optional">
  !      Output mask that specifies whether data was computed. mask_out will be computed only
  !      when horiz_interp_conserve_new_1d is called. mask_out will be computed in
  !      horiz_interp_conserve_new_2d.
  !   </OUT>

  subroutine horiz_interp_conserve ( Interp, data_in, data_out, verbose, &
       mask_in, mask_out)
    !-----------------------------------------------------------------------
    type (horiz_interp_type), intent(in) :: Interp
    real, intent(in),  dimension(:,:) :: data_in
    real, intent(out), dimension(:,:) :: data_out
    integer, intent(in),                   optional :: verbose
    real, intent(in),   dimension(:,:), optional :: mask_in
    real, intent(out),  dimension(:,:), optional :: mask_out

    !  --- error checking ---
    if (size(data_in,1) /= Interp%nlon_src .or. size(data_in,2) /= Interp%nlat_src) &
         call mpp_error(FATAL, 'horiz_interp_conserve_mod: size of input array incorrect')

    if (size(data_out,1) /= Interp%nlon_dst .or. size(data_out,2) /= Interp%nlat_dst) &
         call mpp_error(FATAL, 'horiz_interp_conserve_mod: size of output array incorrect')

    select case ( Interp%version)
    case (1)
       call horiz_interp_conserve_version1(Interp, data_in, data_out, verbose, mask_in, mask_out)
    case (2)
       if(present(mask_in) .OR. present(mask_out) ) call mpp_error(FATAL,  &
            'horiz_interp_conserve: for version 2, mask_in and mask_out must be passed in horiz_interp_new, not in horiz_interp')
       call horiz_interp_conserve_version2(Interp, data_in, data_out, verbose)     
    end select

  end subroutine horiz_interp_conserve
  ! </SUBROUTINE>

  !##############################################################################
  subroutine horiz_interp_conserve_version1 ( Interp, data_in, data_out, verbose, &
       mask_in, mask_out)
    !-----------------------------------------------------------------------
    type (horiz_interp_type), intent(in) :: Interp
    real, intent(in),  dimension(:,:) :: data_in
    real, intent(out), dimension(:,:) :: data_out
    integer, intent(in),                   optional :: verbose
    real, intent(in),   dimension(:,:), optional :: mask_in
    real, intent(out),  dimension(:,:), optional :: mask_out
    !----------local variables----------------------------------------------------
    integer :: m, n, nlon_in, nlat_in, nlon_out, nlat_out,   &
         miss_in, miss_out, is, ie, js, je,   &
         np, npass, iverbose
    real    :: dsum, wsum, avg_in, min_in, max_in,   &
         avg_out, min_out, max_out, eps, asum,   &
         dwtsum, wtsum, arsum, fis, fie, fjs, fje
    !-----------------------------------------------------------------------
    iverbose = 0;  if (present(verbose)) iverbose = verbose

    eps = epsilon(wtsum)

    nlon_in  = Interp%nlon_src;  nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    if (present(mask_in)) then
       if ( count(mask_in < -.0001 .or. mask_in > 1.0001) > 0 ) &
            call mpp_error(FATAL, 'horiz_interp_conserve_mod: input mask not between 0,1')
    endif

    !-----------------------------------------------------------------------
    !---- loop through output grid boxes ----

    data_out = 0.0
    do n = 1, nlat_out
       ! latitude window
       ! setup ascending latitude indices and weights
       if (Interp%jlat(n,1) <= Interp%jlat(n,2)) then
          js = Interp%jlat(n,1); je = Interp%jlat(n,2)
          fjs = Interp%facj(n,1); fje = Interp%facj(n,2)
       else
          js = Interp%jlat(n,2); je = Interp%jlat(n,1)
          fjs = Interp%facj(n,2); fje = Interp%facj(n,1)
       endif

       do m = 1, nlon_out
          ! longitude window
          is = Interp%ilon(m,1); ie = Interp%ilon(m,2)
          fis = Interp%faci(m,1); fie = Interp%faci(m,2)
          npass = 1
          dwtsum = 0.
          wtsum = 0.
          arsum = 0.

          ! wrap-around on input grid
          ! sum using 2 passes (pass 1: end of input grid)
          if ( ie < is ) then
             ie = nlon_in
             fie = 1.0
             npass = 2
          endif

          do np = 1, npass
             ! pass 2: beginning of input grid
             if ( np == 2 ) then
                is = 1
                fis = 1.0
                ie = Interp%ilon(m,2)
                fie = Interp%faci(m,2)
             endif

             ! summing data*weight and weight for single grid point
             if (present(mask_in)) then
                call data_sum ( data_in(is:ie,js:je), Interp%area_src(is:ie,js:je), &
                     fis, fie, fjs,fje, dwtsum, wtsum, arsum, mask_in(is:ie,js:je)  )
             else if( ASSOCIATED(Interp%mask_in) ) then
                call data_sum ( data_in(is:ie,js:je), Interp%area_src(is:ie,js:je), &
                     fis, fie, fjs,fje, dwtsum, wtsum, arsum, Interp%mask_in(is:ie,js:je)  )
             else
                call data_sum ( data_in(is:ie,js:je), Interp%area_src(is:ie,js:je), &
                     fis, fie, fjs,fje,  dwtsum, wtsum, arsum    )
             endif
          enddo

          if (wtsum > eps) then
             data_out(m,n) = dwtsum/wtsum
             if (present(mask_out)) mask_out(m,n) = wtsum/arsum
          else
             data_out(m,n) = 0.
             if (present(mask_out)) mask_out(m,n) = 0.0
          endif

       enddo
    enddo

    !***********************************************************************
    ! compute statistics: minimum, maximum, and mean
    !-----------------------------------------------------------------------

    if (iverbose > 0) then

       ! compute statistics of input data

       call stats(data_in, Interp%area_src, asum, dsum, wsum, min_in, max_in, miss_in, mask_in)
       ! diagnostic messages
       ! on the root_pe, we can calculate the global mean, minimum and maximum.
       if(pe == root_pe) then
          if (wsum > 0.0) then
             avg_in=dsum/wsum
          else
             print *, 'horiz_interp stats: input area equals zero '
             avg_in=0.0
          endif
          if (iverbose > 1) print '(2f16.11)', 'global sum area_in  = ',  asum, wsum
       endif

       ! compute statistics of output data
       call stats(data_out, Interp%area_dst, asum, dsum, wsum, min_out, max_out, miss_out, mask_out)
       ! diagnostic messages
       if(pe == root_pe) then
          if (wsum > 0.0) then
             avg_out=dsum/wsum
          else
             print *, 'horiz_interp stats: output area equals zero '
             avg_out=0.0
          endif
          if (iverbose > 1) print '(2f16.11)', 'global sum area_out = ',  asum, wsum
       endif
       !---- output statistics ----
       ! the global mean, min and max are calculated on the root pe.
       if(pe == root_pe) then
          write (*,900)
          write (*,901)  min_in ,max_in ,avg_in
          if (present(mask_in))  write (*,903)  miss_in
          write (*,902)  min_out,max_out,avg_out
          if (present(mask_out)) write (*,903)  miss_out
       endif

900    format (/,1x,10('-'),' output from horiz_interp ',10('-'))
901    format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
902    format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
903    format ('          number of missing points = ',i6)

    endif

    !-----------------------------------------------------------------------
  end subroutine horiz_interp_conserve_version1

  !#############################################################################
  subroutine horiz_interp_conserve_version2 ( Interp, data_in, data_out, verbose )
    !-----------------------------------------------------------------------
    type (horiz_interp_type), intent(in) :: Interp
    real,    intent(in),  dimension(:,:) :: data_in
    real,    intent(out), dimension(:,:) :: data_out
    integer, intent(in),        optional :: verbose  
    integer :: i, i_src, j_src, i_dst, j_dst

    data_out = 0.0
    do i = 1, Interp%nxgrid
       i_src = Interp%i_src(i); j_src = Interp%j_src(i)
       i_dst = Interp%i_dst(i); j_dst = Interp%j_dst(i)
       data_out(i_dst, j_dst) = data_out(i_dst, j_dst) + data_in(i_src,j_src)*Interp%area_frac_dst(i)
    end do
    
  end subroutine horiz_interp_conserve_version2

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_conserve_del">

  !   <OVERVIEW>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_new.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_new.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_conserve_del ( Interp )
  !   </TEMPLATE>

  !   <INOUT NAME="Interp" TYPE="horiz_interp_type">
  !     A derived-type variable returned by previous call
  !     to horiz_interp_new. The input variable must have
  !     allocated arrays. The returned variable will contain
  !     deallocated arrays.
  !   </INOUT>

  subroutine horiz_interp_conserve_del ( Interp )

    type (horiz_interp_type), intent(inout) :: Interp

    select case(Interp%version)  
    case (1)
       if(associated(Interp%area_src)) deallocate(Interp%area_src)
       if(associated(Interp%area_dst)) deallocate(Interp%area_dst)
       if(associated(Interp%facj))     deallocate(Interp%facj)
       if(associated(Interp%jlat))     deallocate(Interp%jlat)
       if(associated(Interp%faci))     deallocate(Interp%faci)
       if(associated(Interp%ilon))     deallocate(Interp%ilon)
    case (2)
       if(associated(Interp%i_src)) deallocate(Interp%i_src)
       if(associated(Interp%j_src)) deallocate(Interp%j_src)
       if(associated(Interp%i_dst)) deallocate(Interp%i_dst)
       if(associated(Interp%j_dst)) deallocate(Interp%j_dst)
       if(associated(Interp%area_frac_dst)) deallocate(Interp%area_frac_dst)
    end select

  end subroutine horiz_interp_conserve_del
  ! </SUBROUTINE>

  !#######################################################################
  !---This statistics is for conservative scheme
  subroutine stats ( dat, area, asum, dsum, wsum, low, high, miss, mask )
    real,    intent(in)  :: dat(:,:), area(:,:)
    real,    intent(out) :: asum, dsum, wsum, low, high
    integer, intent(out) :: miss
    real,    intent(in), optional :: mask(:,:)

    integer :: pe, root_pe, npes, p, buffer_int(1)
    real    :: buffer_real(5)

    pe = mpp_pe()
    root_pe = mpp_root_pe()
    npes = mpp_npes()

    ! sum data, data*area; and find min,max on each pe.

    if (present(mask)) then
       asum = sum(area(:,:))
       dsum = sum(area(:,:)*dat(:,:)*mask(:,:))
       wsum = sum(area(:,:)*mask(:,:))
       miss = count(mask(:,:) <= 0.5)
       low  = minval(dat(:,:),mask=mask(:,:) > 0.5)
       high = maxval(dat(:,:),mask=mask(:,:) > 0.5)
    else
       asum = sum(area(:,:))
       dsum = sum(area(:,:)*dat(:,:))
       wsum = sum(area(:,:))
       miss = 0
       low  = minval(dat(:,:))
       high = maxval(dat(:,:))
    endif

    ! other pe send local min, max, avg to the root pe and 
    ! root pe receive these information

    if(pe == root_pe) then
       do p = 1, npes - 1
          ! Force use of "scalar", integer pointer mpp interface
          call mpp_recv(buffer_real(1),glen=5,from_pe=root_pe+p, tag=COMM_TAG_1)
          asum = asum + buffer_real(1)
          dsum = dsum + buffer_real(2)
          wsum = wsum + buffer_real(3)
          low  = min(low, buffer_real(4))
          high = max(high, buffer_real(5))
          call mpp_recv(buffer_int(1),glen=1,from_pe=root_pe+p, tag=COMM_TAG_2)
          miss = miss + buffer_int(1)
       enddo
    else
       buffer_real(1) = asum
       buffer_real(2) = dsum
       buffer_real(3) = wsum
       buffer_real(4) = low
       buffer_real(5) = high
       ! Force use of "scalar", integer pointer mpp interface
       call mpp_send(buffer_real(1),plen=5,to_pe=root_pe, tag=COMM_TAG_1)
       buffer_int(1) = miss
       call mpp_send(buffer_int(1),plen=1,to_pe=root_pe, tag=COMM_TAG_2)
    endif

    call mpp_sync_self()   

  end subroutine stats

  !#######################################################################

  subroutine data_sum( data, area, facis, facie, facjs, facje,  &
       dwtsum, wtsum, arsum, mask )

    !  sums up the data and weights for a single output grid box
    !-----------------------------------------------------------------------
    real, intent(in), dimension(:,:) :: data, area
    real, intent(in)                 :: facis, facie, facjs, facje
    real, intent(inout)              :: dwtsum, wtsum, arsum
    real, intent(in), optional       :: mask(:,:)

    !  fac__ = fractional portion of each boundary grid box included
    !          in the integral
    !  dwtsum = sum(data*area*mask)
    !  wtsum  = sum(area*mask)
    !  arsum  = sum(area)
    !-----------------------------------------------------------------------
    real, dimension(size(area,1),size(area,2)) :: wt
    real    :: asum
    integer :: id, jd
    !-----------------------------------------------------------------------

    id=size(area,1); jd=size(area,2) 

    wt=area
    wt( 1,:)=wt( 1,:)*facis
    wt(id,:)=wt(id,:)*facie
    wt(:, 1)=wt(:, 1)*facjs
    wt(:,jd)=wt(:,jd)*facje

    asum = sum(wt)
    arsum = arsum + asum

    if (present(mask)) then
       wt = wt * mask
       dwtsum = dwtsum + sum(wt*data)
       wtsum =  wtsum + sum(wt)
    else
       dwtsum = dwtsum + sum(wt*data)
       wtsum =  wtsum + asum
    endif
    !-----------------------------------------------------------------------

  end subroutine data_sum


  !#######################################################################

end module horiz_interp_conserve_mod


