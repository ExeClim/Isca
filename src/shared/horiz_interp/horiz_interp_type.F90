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

module horiz_interp_type_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!     define derived data type that contains indices and weights used for subsequent 
!      interpolations.
! </OVERVIEW>

! <DESCRIPTION>
!     define derived data type that contains indices and weights used for subsequent 
!      interpolations.
! </DESCRIPTION>


use mpp_mod, only : mpp_send, mpp_recv, mpp_sync_self, mpp_error, FATAL
use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes
use mpp_mod, only : COMM_TAG_1, COMM_TAG_2

implicit none
private


! parameter to determine interpolation method
 integer, parameter :: CONSERVE = 1
 integer, parameter :: BILINEAR = 2
 integer, parameter :: SPHERICA = 3
 integer, parameter :: BICUBIC  = 4

public :: CONSERVE, BILINEAR, SPHERICA, BICUBIC
public :: horiz_interp_type, stats, assignment(=)

interface assignment(=)
  module procedure horiz_interp_type_eq
end interface

!<PUBLICTYPE >
 type horiz_interp_type
   real,    dimension(:,:), pointer   :: faci =>NULL(), facj =>NULL()   !weights for conservative scheme
   integer, dimension(:,:), pointer   :: ilon =>NULL(), jlat =>NULL()   !indices for conservative scheme
   real,    dimension(:,:), pointer   :: area_src =>NULL()              !area of the source grid
   real,    dimension(:,:), pointer   :: area_dst =>NULL()              !area of the destination grid
   real,    dimension(:,:,:), pointer :: wti =>NULL(),wtj =>NULL()      !weights for bilinear interpolation
                                                                        !wti ist used for derivative "weights" in bicubic 
   integer, dimension(:,:,:), pointer :: i_lon =>NULL(), j_lat =>NULL() !indices for bilinear interpolation 
                                                                        !and spherical regrid
   real,    dimension(:,:,:), pointer :: src_dist =>NULL()              !distance between destination grid and 
                                                                        !neighbor source grid.
   logical, dimension(:,:), pointer   :: found_neighbors =>NULL()       !indicate whether destination grid 
                                                                        !has some source grid around it.
   real                               :: max_src_dist
   integer, dimension(:,:), pointer   :: num_found => NULL()
   integer                            :: nlon_src, nlat_src !size of source grid
   integer                            :: nlon_dst, nlat_dst !size of destination grid
   integer                            :: interp_method      !interpolation method.
                                                            !=1, conservative scheme
                                                            !=2, bilinear interpolation
                                                            !=3, spherical regrid
                                                            !=4, bicubic regrid
   real,    dimension(:,:), pointer   :: rat_x =>NULL(), rat_y =>NULL() !the ratio of coordinates of the dest grid
                                                                        ! (x_dest -x_src_r)/(x_src_l -x_src_r) and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real,    dimension(:), pointer     :: lon_in =>NULL(),  lat_in =>NULL()  !the coordinates of the source grid
   logical                            :: I_am_initialized=.false.
   integer                            :: version                            !indicate conservative interpolation version with value 1 or 2
   !--- The following are for conservative interpolation scheme version 2 ( through xgrid)
   integer                            :: nxgrid                             !number of exchange grid between src and dst grid.
   integer, dimension(:), pointer     :: i_src=>NULL(), j_src=>NULL()       !indices in source grid.
   integer, dimension(:), pointer     :: i_dst=>NULL(), j_dst=>NULL()       !indices in destination grid.
   real,    dimension(:), pointer     :: area_frac_dst=>NULL()              !area fraction in destination grid.
   real,    dimension(:,:), pointer   :: mask_in=>NULL() 
 end type
!</PUBLICTYPE>

contains

!#######################################################################
!---This statistics is for bilinear interpolation and spherical regrid.
 subroutine stats ( dat, low, high, avg, miss, missing_value, mask )
 real,    intent(in)  :: dat(:,:)
 real,    intent(out) :: low, high, avg
 integer, intent(out) :: miss
 real, intent(in), optional :: missing_value
 real,    intent(in), optional :: mask(:,:)

 real :: dsum, npts, buffer_real(3)
 integer :: pe, root_pe, npes, p, buffer_int(2)

   pe = mpp_pe()
   root_pe = mpp_root_pe()
   npes = mpp_npes()

   dsum = 0.0
   miss = 0

   if (present(missing_value)) then
      miss = count(dat(:,:) == missing_value)
      low  = minval(dat(:,:), dat(:,:) /= missing_value)
      high = maxval(dat(:,:), dat(:,:) /= missing_value)
      dsum = sum(dat(:,:), dat(:,:) /= missing_value)
   else if(present(mask)) then
      miss = count(mask(:,:) <= 0.5)
      low  = minval(dat(:,:),mask=mask(:,:) > 0.5)
      high = maxval(dat(:,:),mask=mask(:,:) > 0.5)
      dsum = sum(dat(:,:), mask=mask(:,:) > 0.5)
   else
      miss = 0
      low  = minval(dat(:,:))
      high = maxval(dat(:,:))
      dsum = sum(dat(:,:))
   endif
   avg = 0.0
   
   npts = size(dat(:,:)) - miss
   if(pe == root_pe) then
      do p = 1, npes - 1  ! root_pe receive data from other pe
      ! Force use of "scalar", integer pointer mpp interface
         call mpp_recv(buffer_real(1),glen=3, from_pe=p+root_pe, tag=COMM_TAG_1)
         dsum = dsum + buffer_real(1)
         low  = min(low, buffer_real(2))
         high = max(high, buffer_real(3))
         call mpp_recv(buffer_int(1), glen=2, from_pe=p+root_pe, tag=COMM_TAG_2)
         miss = miss + buffer_int(1)
         npts = npts + buffer_int(2)
      enddo         
      if(npts == 0) then
         print*, 'Warning: no points is valid'
      else
         avg = dsum/real(npts)
      endif
    else   ! other pe send data to the root_pe.
      buffer_real(1) = dsum
      buffer_real(2) = low
      buffer_real(3) = high
      ! Force use of "scalar", integer pointer mpp interface
      call mpp_send(buffer_real(1),plen=3,to_pe=root_pe, tag=COMM_TAG_1)
      buffer_int(1) = miss
      buffer_int(2) = npts
      call mpp_send(buffer_int(1), plen=2, to_pe=root_pe, tag=COMM_TAG_2)
    endif

    call mpp_sync_self()

    return

 end subroutine stats

!#################################################################################################################################
 subroutine horiz_interp_type_eq(horiz_interp_out, horiz_interp_in)
    type(horiz_interp_type), intent(inout) :: horiz_interp_out
    type(horiz_interp_type), intent(in)    :: horiz_interp_in

    if(.not.horiz_interp_in%I_am_initialized) then
      call mpp_error(FATAL,'horiz_interp_type_eq: horiz_interp_type variable on right hand side is unassigned')
    endif

    horiz_interp_out%faci            => horiz_interp_in%faci
    horiz_interp_out%facj            => horiz_interp_in%facj
    horiz_interp_out%ilon            => horiz_interp_in%ilon
    horiz_interp_out%jlat            => horiz_interp_in%jlat
    horiz_interp_out%area_src        => horiz_interp_in%area_src
    horiz_interp_out%area_dst        => horiz_interp_in%area_dst
    horiz_interp_out%wti             => horiz_interp_in%wti
    horiz_interp_out%wtj             => horiz_interp_in%wtj
    horiz_interp_out%i_lon           => horiz_interp_in%i_lon
    horiz_interp_out%j_lat           => horiz_interp_in%j_lat
    horiz_interp_out%src_dist        => horiz_interp_in%src_dist
    horiz_interp_out%found_neighbors => horiz_interp_in%found_neighbors
    horiz_interp_out%max_src_dist    =  horiz_interp_in%max_src_dist
    horiz_interp_out%num_found       => horiz_interp_in%num_found
    horiz_interp_out%nlon_src        =  horiz_interp_in%nlon_src
    horiz_interp_out%nlat_src        =  horiz_interp_in%nlat_src
    horiz_interp_out%nlon_dst        =  horiz_interp_in%nlon_dst
    horiz_interp_out%nlat_dst        =  horiz_interp_in%nlat_dst
    horiz_interp_out%interp_method   =  horiz_interp_in%interp_method
    horiz_interp_out%rat_x           => horiz_interp_in%rat_x
    horiz_interp_out%rat_y           => horiz_interp_in%rat_y
    horiz_interp_out%lon_in          => horiz_interp_in%lon_in
    horiz_interp_out%lat_in          => horiz_interp_in%lat_in
    horiz_interp_out%I_am_initialized = .true.
    horiz_interp_out%i_src           => horiz_interp_in%i_src
    horiz_interp_out%j_src           => horiz_interp_in%j_src
    horiz_interp_out%i_dst           => horiz_interp_in%i_dst
    horiz_interp_out%j_dst           => horiz_interp_in%j_dst
    horiz_interp_out%area_frac_dst   => horiz_interp_in%area_frac_dst
    if(horiz_interp_in%interp_method == CONSERVE) then
       horiz_interp_out%version =  horiz_interp_in%version
       if(horiz_interp_in%version==2) horiz_interp_out%nxgrid = horiz_interp_in%nxgrid
    end if

 end subroutine horiz_interp_type_eq
!#################################################################################################################################

end module horiz_interp_type_mod
