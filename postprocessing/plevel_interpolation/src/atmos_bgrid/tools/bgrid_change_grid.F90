
module bgrid_change_grid_mod

! interpolates data between grids

use bgrid_horiz_mod, only: horiz_grid_type
use         fms_mod, only: error_mesg, FATAL, WARNING, &
                           mpp_pe, mpp_root_pe
implicit none
private

 public :: change_grid

 interface change_grid
    module procedure change_grid_2d, change_grid_2flds_2d, &
                     change_grid_3d, change_grid_2flds_3d
 end interface

! Grid identifiers for the four basic horizontal grids:
!
!     TEMP_GRID = temperature (mass) grid
!     WIND_GRID = velocity (momentum) grid
!     UFLX_GRID = zonal mass flux grid; grid points are located
!                  between temperature points along the x-axis
!     VFLX_GRID = meridional mass flux grid; grid points are located
!                  between temperature points along the y-axis
!
! Grid boxes with the same (i,j) indexing have the following
! proximity to one another:
!
!         VFLX(i,j)   WIND(i,j)
!
!         TEMP(i,j)   UFLX(i,j)

 integer, parameter, public :: TEMP_GRID=101, WIND_GRID=102, &
                               UFLX_GRID=104, VFLX_GRID=105

! identifiers for the type of interpolation weighting
 integer, public, parameter :: AREA = 1, EQUAL = 2

!----------------------------------------------------------------
!------ private data -------

! indexing for interpolation weights array
  integer, parameter :: SOUTH = 1, NORTH = 2

! module variables (used between internal interfaces)
  integer :: ilb, iub, jlb, jub
  real, pointer :: weights(:,:) =>NULL()
  integer :: option

contains

!##############################################################################
! Routine to interpolate one field

subroutine change_grid_3d ( Hgrid, grid_inp, grid_out, data_inp, data_out, &
                            weight, mask_inp )

! Arguments
!    Hgrid    = horizontal grid constants
!    grid_inp = grid identifier for the input grid
!    grid_out = grid identifier for the output grid
!    data_inp = input data
!    data_out = output data
! Optional
!    weight = integer flag that type of interpolation
!               weight=AREA for area weighting (the default) or
!               weight=EQUAL for simple 4-pt (or 2pt) averaging
!    mask_inp = grid box mask (typically for step-coordinate model)

type(horiz_grid_type), intent(in) :: Hgrid
integer, intent(in)  :: grid_inp, grid_out
real,    intent(in)  :: data_inp(Hgrid%ilb:,Hgrid%jlb:,:)
real,    intent(out) :: data_out(Hgrid%ilb:,Hgrid%jlb:,:)
integer, intent(in), optional :: weight
real,    intent(in), optional :: mask_inp(Hgrid%ilb:,Hgrid%jlb:,:)

  call set_options ( Hgrid, grid_inp, grid_out, weight )

  select case (option)
    ! velocity grid to temperature grid
    case (1)
        call vel_to_tmp   ( data_inp, data_out, mask=mask_inp )
    case (2)
        call tmp_to_vel   ( data_inp, data_out, mask=mask_inp )
    case (3)
        call vflx_to_uflx ( data_inp, data_out, mask=mask_inp )
    case (4)
        call uflx_to_vflx ( data_inp, data_out, mask=mask_inp )
    case (5)
        if (present(mask_inp)) call mask_warning
        call tmp_to_uflx  ( data_inp, data_out )
    case (6)
        if (present(mask_inp)) call mask_warning
        call tmp_to_vflx  ( data_inp, data_out )
  end select

  nullify (weights)

end subroutine change_grid_3d

!##############################################################################
! Routine to interpolate two fields (on the same grid) at a time

subroutine change_grid_2flds_3d ( Hgrid, grid_inp, grid_out, data_inp1, data_inp2, &
                                         data_out1, data_out2, weight, mask_inp )

! Arguments
!    Hgrid    = horizontal grid constants
!    grid_inp = grid identifier for the input grids
!    grid_out = grid identifier for the output grids
!    data_inp1 = input data for field 1
!    data_inp2 = input data for field 2
!    data_out1 = output data for field 1
!    data_out2 = output data for field 2
! Optional
!    weight = integer flag that type of interpolation
!               weight=AREA for area weighting (the default) or
!               weight=EQUAL for simple 4-pt (or 2pt) averaging
!    mask_inp = grid box mask (typically for step-coordinate model)

type(horiz_grid_type), intent(in) :: Hgrid
integer, intent(in)  :: grid_inp, grid_out
real,    intent(in)  :: data_inp1(Hgrid%ilb:,Hgrid%jlb:,:), &
                        data_inp2(Hgrid%ilb:,Hgrid%jlb:,:)
real,    intent(out) :: data_out1(Hgrid%ilb:,Hgrid%jlb:,:), &
                        data_out2(Hgrid%ilb:,Hgrid%jlb:,:)
integer, intent(in), optional :: weight
real,    intent(in), optional :: mask_inp(Hgrid%ilb:,Hgrid%jlb:,:)

  call set_options ( Hgrid, grid_inp, grid_out, weight )

  select case (option)
    ! velocity grid to temperature grid
    case (1)
        call vel_to_tmp   ( data_inp1, data_out1, data_inp2, data_out2, mask=mask_inp )
    case (2)
        call tmp_to_vel   ( data_inp1, data_out1, data_inp2, data_out2, mask=mask_inp )
    case (3)
        call vflx_to_uflx ( data_inp1, data_out1, data_inp2, data_out2, mask=mask_inp )
    case (4)
        call uflx_to_vflx ( data_inp1, data_out1, data_inp2, data_out2, mask=mask_inp )
    case (5)
        if (present(mask_inp)) call mask_warning
        call tmp_to_uflx  ( data_inp1, data_out1, data_inp2, data_out2 )
    case (6)
        if (present(mask_inp)) call mask_warning
        call tmp_to_vflx  ( data_inp1, data_out1, data_inp2, data_out2 )
  end select

  nullify (weights)

end subroutine change_grid_2flds_3d

!##############################################################################
! ############## private grid-to-grid interpolation interfaces ##############
!##############################################################################

subroutine vel_to_tmp ( idat1, odat1, idat2, odat2, mask )
real,    intent(in)  :: idat1(ilb:,jlb:,:)
real,    intent(out) :: odat1(ilb:,jlb:,:)
real,    intent(in),  optional :: idat2(ilb:,jlb:,:)
real,    intent(out), optional :: odat2(ilb:,jlb:,:)
real,    intent(in),  optional :: mask (ilb:,jlb:,:)

integer :: i, j, k
real, dimension(ilb:iub) :: savg1, navg1, savg2, navg2, smask, nmask

!--------------------------------------------------------------------
! velocity grid to temperature grid
! south and west boundaries not computed
!----------------------------------------

   do k = 1, size(idat1,3)

      ! compute x-average along southernmost row before starting j-loop
      do i = ilb+1, iub
         savg1(i) = idat1(i-1,jlb,k)+idat1(i,jlb,k)
      enddo
      if (present(idat2)) then
         do i = ilb+1, iub
            savg2(i) = idat2(i-1,jlb,k)+idat2(i,jlb,k)
         enddo
      endif
      if (present(mask)) then
         do i = ilb+1, iub
            smask(i) = (mask(i-1,jlb,k)+mask(i,jlb,k))*0.25
         enddo
      endif

      ! loop over latitudes
      do j = jlb+1, jub

         ! compute x-average for north boxes in 4pt average
         do i = ilb+1, iub
            navg1(i) = idat1(i-1,j,k)+idat1(i,j,k)
         enddo
         if (present(idat2)) then
            do i = ilb+1, iub
               navg2(i) = idat2(i-1,j,k)+idat2(i,j,k)
            enddo
         endif
         if (present(mask)) then
            do i = ilb+1, iub
               nmask(i) = (mask(i-1,j,k)+mask(i,j,k))*0.25
            enddo
         endif

         ! compute 4pt average - then update south average with north average
         do i = ilb+1, iub
            odat1(i,j,k) = savg1(i) * weights(j,SOUTH) + &
                           navg1(i) * weights(j,NORTH)
            savg1(i) = navg1(i)
         enddo
         if (present(idat2)) then
            do i = ilb+1, iub
               odat2(i,j,k) = savg2(i) * weights(j,SOUTH) + &
                              navg2(i) * weights(j,NORTH)
               savg2(i) = navg2(i)
            enddo
         endif
         ! apply weighting for masked grid boxes
         if (present(mask)) then
            do i = ilb+1, iub
               if (smask(i)+nmask(i) > 0.) then
                   odat1(i,j,k) = odat1(i,j,k)/(smask(i)+nmask(i))
                   if (present(idat2)) odat2(i,j,k) = odat2(i,j,k)/(smask(i)+nmask(i))
               else
                   odat1(i,j,k) = 0.
                   if (present(idat2)) odat2(i,j,k) = 0.
               endif
               smask(i) = nmask(i)
            enddo
         endif

      enddo ! end j-loop
   enddo ! end k-loop

 ! bogus data at non-computed grid boxes
   odat1(:,jlb,:) = 1.e20
   odat1(ilb,:,:) = 1.e20
   if (present(odat2)) then
       odat2(:,jlb,:) = 1.e20
       odat2(ilb,:,:) = 1.e20
   endif

end subroutine vel_to_tmp

!##############################################################################

subroutine tmp_to_vel ( idat1, odat1, idat2, odat2, mask )
real,    intent(in)  :: idat1(ilb:,jlb:,:)
real,    intent(out) :: odat1(ilb:,jlb:,:)
real,    intent(in),  optional :: idat2(ilb:,jlb:,:)
real,    intent(out), optional :: odat2(ilb:,jlb:,:)
real,    intent(in),  optional :: mask (ilb:,jlb:,:)

integer :: i, j, k
real, dimension(ilb:iub) :: savg1, navg1, savg2, navg2, smask, nmask

!--------------------------------------------------------------------
! temperature grid to velocity grid
! north and east boundaries not computed

   do k = 1, size(idat1,3)

      ! compute x-average along southernmost row before starting j-loop
      do i = ilb, iub-1
         savg1(i) = idat1(i,jlb,k)+idat1(i+1,jlb,k)
      enddo
      if (present(idat2)) then
         do i = ilb, iub-1
            savg2(i) = idat2(i,jlb,k)+idat2(i+1,jlb,k)
         enddo
      endif
      if (present(mask)) then
         do i = ilb, iub-1
            smask(i) = (mask(i,jlb,k)+mask(i+1,jlb,k))*0.25
         enddo
      endif

      ! loop over latitudes
      do j = jlb, jub-1

         ! compute x-average for north boxes in 4pt average
         do i = ilb, iub-1
            navg1(i) = idat1(i,j+1,k)+idat1(i+1,j+1,k)
         enddo
         if (present(idat2)) then
            do i = ilb, iub-1
               navg2(i) = idat2(i,j+1,k)+idat2(i+1,j+1,k)
            enddo
         endif
         if (present(mask)) then
            do i = ilb, iub-1
               nmask(i) = (mask(i,j+1,k)+mask(i+1,j+1,k))*0.25
            enddo
         endif

         ! compute 4pt average - then update south average with north average
         do i = ilb, iub-1
            odat1(i,j,k) = savg1(i) * weights(j,SOUTH) + &
                           navg1(i) * weights(j,NORTH)
            savg1(i) = navg1(i)
         enddo
         if (present(idat2)) then
            do i = ilb, iub-1
               odat2(i,j,k) = savg2(i) * weights(j,SOUTH) + &
                              navg2(i) * weights(j,NORTH)
               savg2(i) = navg2(i)
            enddo
         endif
         ! apply weighting for masked grid boxes
         if (present(mask)) then
            do i = ilb, iub-1
               if (smask(i)+nmask(i) > 0.) then
                   odat1(i,j,k) = odat1(i,j,k)/(smask(i)+nmask(i))
                   if (present(idat2)) odat2(i,j,k) = odat2(i,j,k)/(smask(i)+nmask(i))
               else
                   odat1(i,j,k) = 0.
                   if (present(idat2)) odat2(i,j,k) = 0.
               endif
               smask(i) = nmask(i)
            enddo
         endif

      enddo ! end j-loop
   enddo ! end k-loop

 ! bogus data at non-computed grid boxes
   odat1(:,jub,:) = 1.e20
   odat1(iub,:,:) = 1.e20
   if (present(odat2)) then
       odat2(:,jub,:) = 1.e20
       odat2(iub,:,:) = 1.e20
   endif

end subroutine tmp_to_vel

!##############################################################################

subroutine vflx_to_uflx ( idat1, odat1, idat2, odat2, mask )
real,    intent(in)  :: idat1(ilb:,jlb:,:)
real,    intent(out) :: odat1(ilb:,jlb:,:)
real,    intent(in),  optional :: idat2(ilb:,jlb:,:)
real,    intent(out), optional :: odat2(ilb:,jlb:,:)
real,    intent(in),  optional :: mask (ilb:,jlb:,:)

integer :: i, j, k
real, dimension(ilb:iub) :: savg1, navg1, savg2, navg2, smask, nmask

!--------------------------------------------------------------------
! velocity grid to temperature grid
! south and east boundaries not computed
!----------------------------------------

   do k = 1, size(idat1,3)

      ! compute x-average along southernmost row before starting j-loop
      do i = ilb, iub-1
         savg1(i) = idat1(i,jlb,k)+idat1(i+1,jlb,k)
      enddo
      if (present(idat2)) then
         do i = ilb, iub-1
            savg2(i) = idat2(i,jlb,k)+idat2(i+1,jlb,k)
         enddo
      endif
      if (present(mask)) then
         do i = ilb, iub-1
            smask(i) = (mask(i,jlb,k)+mask(i+1,jlb,k))*0.25
         enddo
      endif

      ! loop over latitudes
      do j = jlb+1, jub

         ! compute x-average for north boxes in 4pt average
         do i = ilb, iub-1
            navg1(i) = idat1(i,j,k)+idat1(i+1,j,k)
         enddo
         if (present(idat2)) then
            do i = ilb, iub-1
               navg2(i) = idat2(i,j,k)+idat2(i+1,j,k)
            enddo
         endif
         if (present(mask)) then
            do i = ilb, iub-1
               nmask(i) = (mask(i,j,k)+mask(i+1,j,k))*0.25
            enddo
         endif

         ! compute 4pt average - then update south average with north average
         do i = ilb, iub-1
            odat1(i,j,k) = savg1(i) * weights(j,SOUTH) + &
                           navg1(i) * weights(j,NORTH)
            savg1(i) = navg1(i)
         enddo
         if (present(idat2)) then
            do i = ilb, iub-1
               odat2(i,j,k) = savg2(i) * weights(j,SOUTH) + &
                              navg2(i) * weights(j,NORTH)
               savg2(i) = navg2(i)
            enddo
         endif
         ! apply weighting for masked grid boxes
         if (present(mask)) then
            do i = ilb, iub-1
               if (smask(i)+nmask(i) > 0.) then
                   odat1(i,j,k) = odat1(i,j,k)/(smask(i)+nmask(i))
                   if (present(idat2)) odat2(i,j,k) = odat2(i,j,k)/(smask(i)+nmask(i))
               else
                   odat1(i,j,k) = 0.
                   if (present(idat2)) odat2(i,j,k) = 0.
               endif
               smask(i) = nmask(i)
            enddo
         endif

      enddo ! end j-loop
   enddo ! end k-loop

 ! bogus data at non-computed grid boxes
   odat1(:,jlb,:) = 1.e20
   odat1(iub,:,:) = 1.e20
   if (present(odat2)) then
       odat2(:,jlb,:) = 1.e20
       odat2(iub,:,:) = 1.e20
   endif

end subroutine vflx_to_uflx

!##############################################################################

subroutine uflx_to_vflx ( idat1, odat1, idat2, odat2, mask )
real,    intent(in)  :: idat1(ilb:,jlb:,:)
real,    intent(out) :: odat1(ilb:,jlb:,:)
real,    intent(in),  optional :: idat2(ilb:,jlb:,:)
real,    intent(out), optional :: odat2(ilb:,jlb:,:)
real,    intent(in),  optional :: mask (ilb:,jlb:,:)

integer :: i, j, k
real, dimension(ilb:iub) :: savg1, navg1, savg2, navg2, smask, nmask

!--------------------------------------------------------------------
! uflx points to vflx points
! north and west boundaries are not computed

   do k = 1, size(idat1,3)

      ! compute x-average along southernmost row before starting j-loop
      do i = ilb+1, iub
         savg1(i) = idat1(i-1,jlb,k)+idat1(i,jlb,k)
      enddo
      if (present(idat2)) then
         do i = ilb+1, iub
            savg2(i) = idat2(i-1,jlb,k)+idat2(i,jlb,k)
         enddo
      endif
      if (present(mask)) then
         do i = ilb+1, iub
            smask(i) = (mask(i-1,jlb,k)+mask(i,jlb,k))*0.25
         enddo
      endif

      ! loop over latitudes
      do j = jlb, jub-1

         ! compute x-average for north boxes in 4pt average
         do i = ilb+1, iub
            navg1(i) = idat1(i-1,j+1,k)+idat1(i,j+1,k)
         enddo
         if (present(idat2)) then
            do i = ilb+1, iub
               navg2(i) = idat2(i-1,j+1,k)+idat2(i,j+1,k)
            enddo
         endif
         if (present(mask)) then
            do i = ilb+1, iub
               nmask(i) = (mask(i-1,j+1,k)+mask(i,j+1,k))*0.25
            enddo
         endif

         ! compute 4pt average - then update south average with north average
         do i = ilb+1, iub
            odat1(i,j,k) = savg1(i) * weights(j,SOUTH) + &
                           navg1(i) * weights(j,NORTH)
            savg1(i) = navg1(i)
         enddo
         if (present(idat2)) then
            do i = ilb+1, iub
               odat2(i,j,k) = savg2(i) * weights(j,SOUTH) + &
                              navg2(i) * weights(j,NORTH)
               savg2(i) = navg2(i)
            enddo
         endif
         ! apply weighting for masked grid boxes
         if (present(mask)) then
            do i = ilb+1, iub
               if (smask(i)+nmask(i) > 0.) then
                   odat1(i,j,k) = odat1(i,j,k)/(smask(i)+nmask(i))
                   if (present(idat2)) odat2(i,j,k) = odat2(i,j,k)/(smask(i)+nmask(i))
               else
                   odat1(i,j,k) = 0.
                   if (present(idat2)) odat2(i,j,k) = 0.
               endif
               smask(i) = nmask(i)
            enddo
         endif

      enddo ! end j-loop
   enddo ! end k-loop

 ! bogus data at non-computed grid boxes
   odat1(:,jub,:) = 1.e20
   odat1(ilb,:,:) = 1.e20
   if (present(odat2)) then
       odat2(:,jub,:) = 1.e20
       odat2(ilb,:,:) = 1.e20
   endif

end subroutine uflx_to_vflx

!##############################################################################

subroutine tmp_to_uflx ( idat1, odat1, idat2, odat2 )
real,    intent(in)  :: idat1(ilb:,jlb:,:)
real,    intent(out) :: odat1(ilb:,jlb:,:)
real,    intent(in),  optional :: idat2(ilb:,jlb:,:)
real,    intent(out), optional :: odat2(ilb:,jlb:,:)

integer :: i, j, k
real, dimension(ilb:iub) :: avg1, avg2

!--------------------------------------------------------------------
! temperature points to uflx points
! east boundary not computed

   do k = 1, size(odat1,3)
      do j = jlb, jub

         do i = ilb, iub-1
            avg1(i) = (idat1(i,j,k)+idat1(i+1,j,k))*0.5
         enddo
         do i = ilb, iub-1
            odat1(i,j,k) = avg1(i)
         enddo
         if (present(idat2)) then
            do i = ilb, iub-1
               avg2(i) = (idat2(i,j,k)+idat2(i+1,j,k))*0.5
            enddo
            do i = ilb, iub-1
               odat2(i,j,k) = avg2(i)
            enddo
         endif

      enddo ! end j-loop
   enddo ! end k-loop

 ! bogus data at non-computed grid boxes
   odat1(iub,:,:) = 1.e20
   if (present(odat2)) odat2(iub,:,:) = 1.e20

end subroutine tmp_to_uflx

!##############################################################################

subroutine tmp_to_vflx ( idat1, odat1, idat2, odat2 )
real,    intent(in)  :: idat1(ilb:,jlb:,:)
real,    intent(out) :: odat1(ilb:,jlb:,:)
real,    intent(in),  optional :: idat2(ilb:,jlb:,:)
real,    intent(out), optional :: odat2(ilb:,jlb:,:)

integer :: i, j, k
real, dimension(jlb:jub) :: avg1, avg2

!--------------------------------------------------------------------
! temperature points to vflx points
! north boundary not computed

   do k = 1, size(idat1,3)
      do i = ilb, iub

         do j = jlb, jub-1
            avg1(j) = (idat1(i,j  ,k)*weights(j,SOUTH) + &
                       idat1(i,j+1,k)*weights(j,NORTH))*2.0
         enddo
         do j = jlb, jub-1
            odat1(i,j,k) = avg1(j)
         enddo
         if (present(idat2)) then
            do j = jlb, jub-1
               avg2(j) = (idat2(i,j  ,k)*weights(j,SOUTH) + &
                          idat2(i,j+1,k)*weights(j,NORTH))*2.0
            enddo
            do j = jlb, jub-1
               odat2(i,j,k) = avg2(j)
            enddo
         endif

      enddo ! end i-loop
   enddo ! end k-loop

 ! bogus data at non-computed grid boxes
   odat1(:,jub,:) = 1.e20
   if (present(odat2)) odat1(:,jub,:) = 1.e20

end subroutine tmp_to_vflx

!##############################################################################

subroutine set_options ( Hgrid, grid_inp, grid_out, weight )
type(horiz_grid_type), intent(in) :: Hgrid
integer, intent(in)  :: grid_inp, grid_out
integer, intent(in), optional :: weight
logical :: use_area_weight

  option = 0
! valid 4pt interpolations
  if (  grid_inp == WIND_GRID .and. grid_out == TEMP_GRID ) option = 1
  if (  grid_inp == TEMP_GRID .and. grid_out == WIND_GRID ) option = 2
  if (  grid_inp == VFLX_GRID .and. grid_out == UFLX_GRID ) option = 3
  if (  grid_inp == UFLX_GRID .and. grid_out == VFLX_GRID ) option = 4
! valid 2pt interpolations
  if (  grid_inp == TEMP_GRID .and. grid_out == UFLX_GRID ) option = 5
  if (  grid_inp == TEMP_GRID .and. grid_out == VFLX_GRID ) option = 6
! error condition (no need to check later)
  if (option == 0) &
      call error_mesg ('change_grid', 'invalid grid change specified', FATAL)

! weighting flag (area or equal weighting of grid boxes)
  use_area_weight = .true.
  if (present(weight)) then
     if (weight == AREA) then
         use_area_weight = .true.
     else if (weight == EQUAL) then
         use_area_weight = .false.
     else
         call error_mesg ('change_grid', 'invalid value of area_weight', FATAL)
     endif
  endif

! assign pointer for weighting
! will nullify after averaging
  if (use_area_weight) then
     if (option == 1 .or. option == 3) weights => Hgrid%Interp%tmpwts
     if (option == 2 .or. option == 4) weights => Hgrid%Interp%velwts
     if (option == 5)                  weights => Hgrid%Interp%nowts
     if (option == 6)                  weights => Hgrid%Interp%velwts
  else
     weights => Hgrid%Interp%nowts
  endif

! horizontal array limits
  ilb = Hgrid%ilb;  iub = Hgrid%iub
  jlb = Hgrid%jlb;  jub = Hgrid%jub

end subroutine set_options

!##############################################################################
! prints warning message that mask was not used

subroutine mask_warning
   if (mpp_pe() == mpp_root_pe()) &
       call error_mesg ('bgrid_change_grid', 'optional mask argument &
              &not used for this type of grid interpolation', WARNING)
end subroutine mask_warning

!##############################################################################
! overloaded interfaces
!----------------------------------------------------

subroutine change_grid_2d ( Hgrid, grid_inp, grid_out, data_inp, data_out, weight )
type(horiz_grid_type), intent(in) :: Hgrid
integer, intent(in)  :: grid_inp, grid_out
real,    intent(in)  :: data_inp(Hgrid%ilb:,Hgrid%jlb:)
real,    intent(out) :: data_out(Hgrid%ilb:,Hgrid%jlb:)
integer, intent(in), optional :: weight
real :: idat3(size(data_inp,1),size(data_inp,2),1)
real :: odat3(size(data_out,1),size(data_out,2),1)

  idat3(:,:,1) = data_inp
  call change_grid_3d ( Hgrid, grid_inp, grid_out, idat3, odat3, weight )
  data_out = odat3(:,:,1)

end subroutine change_grid_2d

!----------------------------------------------------

subroutine change_grid_2flds_2d ( Hgrid, grid_inp, grid_out, data_inp1, data_inp2, &
                                         data_out1, data_out2, weight )
type(horiz_grid_type), intent(in) :: Hgrid
integer, intent(in)  :: grid_inp, grid_out
real,    intent(in)  :: data_inp1(Hgrid%ilb:,Hgrid%jlb:), &
                        data_inp2(Hgrid%ilb:,Hgrid%jlb:)
real,    intent(out) :: data_out1(Hgrid%ilb:,Hgrid%jlb:), &
                        data_out2(Hgrid%ilb:,Hgrid%jlb:)
integer, intent(in), optional :: weight
real :: idat1_3d(size(data_inp1,1),size(data_inp1,2),1)
real :: idat2_3d(size(data_inp2,1),size(data_inp2,2),1)
real :: odat1_3d(size(data_out1,1),size(data_out1,2),1)
real :: odat2_3d(size(data_out2,1),size(data_out2,2),1)

  idat1_3d(:,:,1) = data_inp1
  idat2_3d(:,:,1) = data_inp2
  call change_grid_2flds_3d ( Hgrid, grid_inp, grid_out, idat1_3d, idat2_3d, &
                               odat1_3d, odat2_3d, weight )
  data_out1 = odat1_3d(:,:,1)
  data_out2 = odat2_3d(:,:,1)

end subroutine change_grid_2flds_2d

!##############################################################################

end module bgrid_change_grid_mod

