
module bgrid_halo_mod

use bgrid_horiz_mod, only: horiz_grid_type, update_np, update_sp, &
                           TGRID, VGRID
use         fms_mod, only: error_mesg, FATAL, mpp_clock_id, &
                           mpp_clock_begin, mpp_clock_end,  &
                           MPP_CLOCK_SYNC, CLOCK_ROUTINE
use mpp_domains_mod, only: mpp_update_domains, SUPDATE, NUPDATE,  &
                                               WUPDATE, EUPDATE

implicit none
private


!------------ public interfaces ------------

public :: update_halo, vel_flux_boundary

interface update_halo
    module procedure  update_halo_2d, update_halo_3d, update_halo_4d
end interface

interface vel_flux_boundary
    module procedure  vel_flux_boundary_2d, vel_flux_boundary_3d
end interface

!-------- public parameters ----------

! possible values for field argument
integer, parameter, public :: TEMP = 21, UWND = 22, VWND = 23, WIND = 24

! possible values for optional halo argument
integer, parameter, public  :: SOUTH = 1, NORTH = 2
integer, parameter, public  ::  WEST = 4,  EAST = 8
integer, parameter, private ::  ALL = SOUTH+NORTH+WEST+EAST ! default

! possible values for optional flags argument
integer, parameter, public :: NOPOLE   = 1
integer, parameter, public :: POLEONLY = 2

!---------------------------------------------------------------
! private timing variables

   integer :: id_update3
   logical :: do_clock_init = .true.

! private module variables (used across module subroutines)

   integer :: domain_flags
   logical :: update_sbnd, update_nbnd, update_wbnd, update_ebnd, &
              no_pole_vel, do_pole_only
   logical :: do_channel

contains

!#######################################################################

 subroutine update_halo_3d (Hgrid, field, data, halos, flags)

!----------------------------------------------------------
!  Halo update for 3-dimensional fields
!----------------------------------------------------------
!  Hgrid = horizontal grid constants
!  field = integer identifier for field
!            possible values: TEMP, UWND, VWND
!  data  = 3-D data
!  halos = identifies which halos should be updated (optional)
!            possible values: NORTH, EAST, WEST, SOUTH
!            default(all): NORTH+EAST+WEST+SOUTH
!  flags = options for handling the pole
!            possible values: NOPOLE, POLEONLY
!----------------------------------------------------------

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer,               intent(in)    :: field
   real,                  intent(inout) :: data(Hgrid%ilb:,Hgrid%jlb:,:)
   integer, optional,     intent(in)    :: halos, flags

   integer :: is, ie, iflags, halosize, xygrid

   if (do_clock_init) call clock_init
   call mpp_clock_begin (id_update3)

!  ----- check dimensions ------

   if (size(data,2) /= Hgrid % jsize)  call error_mesg  &
                    ('update_halo in bgrid_halo_mod',   &
                     'j dimension has wrong size', FATAL)

   if (size(data,1) /= Hgrid % isize)  call error_mesg  &
                    ('update_halo in bgrid_halo_mod',   &
                     'i dimension has wrong size', FATAL)

!  ----- check/set optional flag arguments ----

   call set_domain_flags ( Hgrid, halos, flags )

!  ------ need to determine and check grid -------

   select case (field)
     case (TEMP)
        is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie;  xygrid = TGRID
        ! update non-polar boundaries
        if (.not.do_pole_only) call mpp_update_domains (data, Hgrid%Tmp%Domain, domain_flags)
     case (UWND:VWND)
        is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie;  xygrid = VGRID
        ! update non-polar boundaries
        if (.not.do_pole_only) call mpp_update_domains (data, Hgrid%Tmp%Domain, domain_flags)
     case default
        call error_mesg ('update_halo in bgrid_halo_mod', &
                         'invalid field', FATAL)
   end select

!  ----- update east-west cyclic boundaries (for 1-d decomp only) ----

   halosize = 1
   if (update_wbnd) data(is-halosize:is-1,:,:) = data(ie-halosize+1:ie,:,:)
   if (update_ebnd) data(ie+1:ie+halosize,:,:) = data(is:is+halosize-1,:,:)

!  ------ update south pole ------

   if ( (update_sbnd.or.do_pole_only) .and. update_sp (Hgrid,xygrid) ) then
      call south_boundary_3d (Hgrid, field, data(:,:,:), no_pole_vel)
   endif

!  ------ update north pole ------

   if ( (update_nbnd.or.do_pole_only) .and. update_np (Hgrid,xygrid) ) then
      call north_boundary_3d (Hgrid, field, data(:,:,:), no_pole_vel)
   endif

   call mpp_clock_end (id_update3)

 end subroutine update_halo_3d

!#######################################################################

 subroutine update_halo_4d (Hgrid, field, data, halos, flags)

!----------------------------------------------------------
!  Halo update for 4-dimensional fields
!----------------------------------------------------------
!  Hgrid = horizontal grid constants
!  field = integer identifier for field
!            possible values: TEMP, UWND, VWND, WIND
!            note: if field=WIND, then 4th dim of data must be 2
!  data  = 4-D data
!  halos = identifies which halos should be updated (optional)
!            possible values: NORTH, EAST, WEST, SOUTH
!            default(all): NORTH+EAST+WEST+SOUTH
!  flags = options for handling the pole
!            possible values: NOPOLE, POLEONLY
!----------------------------------------------------------

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer,               intent(in)    :: field
   real,                  intent(inout) :: data(Hgrid%ilb:,Hgrid%jlb:,:,:)
   integer, optional,     intent(in)    :: halos, flags

   integer :: is, ie, iflags, n, halosize, xygrid

   if (do_clock_init) call clock_init
   call mpp_clock_begin (id_update3)

!  ----- check dimensions ------

   if (size(data,2) /= Hgrid % jsize)  call error_mesg  &
                    ('update_halo in bgrid_halo_mod',   &
                     'j dimension has wrong size', FATAL)

   if (size(data,1) /= Hgrid % isize)  call error_mesg  &
                    ('update_halo in bgrid_halo_mod',   &
                     'i dimension has wrong size', FATAL)

   if (field == WIND .and. size(data,4) /= 2) call error_mesg  &
             ('update_halo in bgrid_halo_mod',   &
      '4th dimension must have size 2 for wind components', FATAL)

!  ----- check/set optional flag arguments ----

   call set_domain_flags ( Hgrid, halos, flags )

!  ------ need to determine and check grid -------

   select case (field)
     case (TEMP)
        is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie;  xygrid = TGRID
        ! update non-polar boundaries
        if (.not.do_pole_only) call mpp_update_domains (data, Hgrid%Tmp%Domain, domain_flags)
     case (UWND:WIND)
        is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie;  xygrid = VGRID
        ! update non-polar boundaries
        if (.not.do_pole_only) call mpp_update_domains (data, Hgrid%Tmp%Domain, domain_flags)
     case default
        call error_mesg ('update_halo in bgrid_halo_mod', &
                         'invalid field', FATAL)
   end select

!  ----- update east-west cyclic boundaries (for 1-d decomp only) ----

       halosize = 1
       if (update_wbnd) data(is-halosize:is-1,:,:,:) = data(ie-halosize+1:ie,:,:,:)
       if (update_ebnd) data(ie+1:ie+halosize,:,:,:) = data(is:is+halosize-1,:,:,:)

!  ------ update south pole ------

   do_channel = Hgrid%channel

   if ( (update_sbnd.or.do_pole_only) .and. update_sp (Hgrid,xygrid) ) then
      if (field == WIND) then
          call south_boundary_3d (Hgrid, UWND, data(:,:,:,1), no_pole_vel)
          call south_boundary_3d (Hgrid, VWND, data(:,:,:,2), no_pole_vel)
       else
          do n = 1, size(data,4)
            call south_boundary_3d (Hgrid, field, data(:,:,:,n), no_pole_vel)
          enddo
       endif
   endif

!  ------ update north pole ------

   if ( (update_nbnd.or.do_pole_only) .and. update_np (Hgrid,xygrid) ) then
      if (field == WIND) then
          call north_boundary_3d (Hgrid, UWND, data(:,:,:,1), no_pole_vel)
          call north_boundary_3d (Hgrid, VWND, data(:,:,:,2), no_pole_vel)
      else
          do n = 1, size(data,4)
            call north_boundary_3d (Hgrid, field, data(:,:,:,n), no_pole_vel)
          enddo
      endif
   endif

   call mpp_clock_end (id_update3)

 end subroutine update_halo_4d

!#######################################################################
! updates halos at (and beyond) the north pole

 subroutine north_boundary_3d (Hgrid, field, data, nopole)

   type(horiz_grid_type), intent(in)    :: Hgrid
   integer,               intent(in)    :: field
   real,                  intent(inout) :: data(:,Hgrid%jlb:,:)
   logical,               intent(in)    :: nopole

   integer :: js, je, jeg, halo

      halo = 1  ! assumed halo size = 1

! --- update north pole boundary ---

      select case (field)
         case (TEMP)
!        ---- mass ----
            js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je
            data (:, je+1:je+halo, :) = data (:, je:je-halo+1:-1, :)
         case (UWND)
!        ---- u comp ----
            js = Hgrid % Vel % js;  je = Hgrid % Vel % je;  jeg = Hgrid % Vel % jeg
            if (.not. nopole) then
                if (.not.do_channel) then
                    if ( jeg+1 <= je+halo ) data (:, jeg+1,:) = 0.0
                else
                    if ( jeg+1 <= je+halo ) data (:, jeg+1,:) = data (:, jeg,:)
                endif
            endif
            if ( jeg+2 <= je+halo ) &
            data (:, jeg+2:jeg+halo, :) = data (:, je:je-halo+2:-1, :)
         case (VWND)
!        ---- v comp ----
            js = Hgrid % Vel % js;  je = Hgrid % Vel % je;  jeg = Hgrid % Vel % jeg
            if (.not. nopole) then
                if ( jeg+1 <= je+halo ) data (:, jeg+1, :) = 0.0
            endif
            if ( jeg+2 <= je+halo ) &
            data (:, jeg+2:jeg+halo, :) = - data (:, je:je-halo+2:-1, :)
      end select

 end subroutine north_boundary_3d

!#######################################################################
! updates halos at (and beyond) the south pole

 subroutine south_boundary_3d (Hgrid, field, data, nopole)

   type(horiz_grid_type), intent(in)    :: Hgrid
   integer,               intent(in)    ::  field
   real,                  intent(inout) :: data(:,Hgrid%jlb:,:)
   logical,               intent(in)    :: nopole
      

   integer :: js, je, halo

      halo = 1  ! assumed halo size = 1

! --- update south pole boundary ---

      select case (field)
         case (TEMP)
!        ---- mass ----
            js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je
            data (:, js-1:js-halo:-1, :) = data (:, js:js+halo-1, :)
         case (UWND)
!        ---- u comp ----
            js = Hgrid % Vel % js;  je = Hgrid % Vel % je
            if (.not. nopole) then
                if (.not.do_channel) then
                    data (:, js-1, :) = 0.0
                else
                    data (:, js-1, :) = data (:, js, :)
                endif
            endif
            data (:, js-2:js-halo:-1, :) = data (:, js:js+halo-2, :)
         case (VWND)
!        ---- v comp ----
            js = Hgrid % Vel % js;  je = Hgrid % Vel % je
            if (.not. nopole) data (:, js-1, :) = 0.0
            data (:, js-2:js-halo:-1, :) = - data (:, js:js+halo-2, :)
      end select

 end subroutine south_boundary_3d

!#######################################################################

 subroutine vel_flux_boundary_3d (Hgrid, data)

! zero-out the flux between pole and first velocity row
! do this on both sides of the ppole
! meridional indexing coincides with mass grid
! assumed halo size is one

!  Hgrid = horizontal grid constants
!  data = 3-D data

   type(horiz_grid_type), intent(in)    :: Hgrid
   real,                  intent(inout) :: data(:,Hgrid%jlb:,:)

      if ( update_sp (Hgrid,TGRID) ) then
           data (:, Hgrid%Tmp%js-1, :) = 0.0
           data (:, Hgrid%Tmp%js  , :) = 0.0
      endif

      if ( update_np (Hgrid,TGRID) ) then
           data (:, Hgrid%Tmp%je  , :) = 0.0
           data (:, Hgrid%Tmp%je+1, :) = 0.0
      endif

 end subroutine vel_flux_boundary_3d

!#######################################################################
! overload interfaces
!#######################################################################

 subroutine update_halo_2d (Hgrid, field, data, halos, flags)

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer,               intent(in)    :: field
   real,                  intent(inout) :: data(:,:)
   integer, optional,     intent(in)    :: halos, flags

   real, dimension(size(data,1),size(data,2),1) :: data3

   data3(:,:,1) = data
   call update_halo_3d (Hgrid, field, data3, halos, flags)
   data = data3(:,:,1)

 end subroutine update_halo_2d

!#######################################################################

 subroutine vel_flux_boundary_2d (Hgrid, data)

   type(horiz_grid_type), intent(in)    :: Hgrid
   real,                  intent(inout) :: data(:,Hgrid%jlb:)

   real, dimension(size(data,1),size(data,2),1) :: data3

   data3(:,:,1) = data
   call vel_flux_boundary_3d (Hgrid, data3)
   data = data3(:,:,1)

 end subroutine vel_flux_boundary_2d

!#######################################################################
! private interfaces
!#######################################################################

 subroutine set_domain_flags ( Hgrid, halos, flags )
 type(horiz_grid_type), intent(in) :: Hgrid
 integer, optional,     intent(in) :: halos, flags
 integer :: ihalos, iflags

! sets module variables -- domain_flags,
!                          update_sbnd, update_nbnd,
!                          update_wbnd, update_ebnd,
!                          no_pole_vel, do_pole_only

   ihalos = ALL;  if (present(halos)) ihalos = halos
   iflags = 0;    if (present(flags)) iflags = flags

   domain_flags = 0

   ! south and north boundary
     update_sbnd  = btest(ihalos,0)
     update_nbnd  = btest(ihalos,1)
     if ( update_sbnd ) domain_flags = domain_flags + SUPDATE
     if ( update_nbnd ) domain_flags = domain_flags + NUPDATE
     ! turn off polar update if domain is periodic in y
    !if (Hgrid%Vel%jeg == Hgrid%Tmp%jeg) then
     if (Hgrid%double_periodic) then
        update_sbnd = .false.
        update_nbnd = .false.
     endif
   ! west and east boundary
     update_wbnd  = btest(ihalos,2)
     update_ebnd  = btest(ihalos,3)
     if (Hgrid%decompx) then
        if ( update_wbnd ) domain_flags = domain_flags + WUPDATE
        if ( update_ebnd ) domain_flags = domain_flags + EUPDATE
        update_wbnd = .false. ! turn off ?
        update_ebnd = .false.
     endif

   ! set additional flags related to polar boundary
     no_pole_vel  = btest(iflags,0)    !4)
     do_pole_only = btest(iflags,1)    !5)

   ! might update south or north poles ?
   ! if (do_pole_only) then
   !     update_sbnd = .true.
   !     update_nbnd = .true.
   ! endif

 end subroutine set_domain_flags

!#######################################################################
! initializes mpp performance clock for bgrid halo updates

 subroutine clock_init
   id_update3 = mpp_clock_id ('BGRID: update_halo',   &
                  flags=MPP_CLOCK_SYNC, grain=CLOCK_ROUTINE)
   do_clock_init = .false. 
 end subroutine clock_init

!#######################################################################

end module bgrid_halo_mod

