
module bgrid_horiz_mod

!-----------------------------------------------------------------------
!
!        allocates memory and initializes grid constants
!              for the FMS B-grid dynamical core
!
!-----------------------------------------------------------------------

use mpp_domains_mod, only: mpp_domains_init,           &
                           mpp_define_domains,         &
                           domain1D, domain2D,         &
                           mpp_get_global_domain,      &
                           mpp_get_data_domain,        &
                           mpp_get_compute_domain,     &
                           mpp_get_compute_domains,    &
                           mpp_get_domain_components,  &
                           mpp_get_pelist,             &
                           mpp_define_layout,          &
                           mpp_get_layout,             &
                           CYCLIC_GLOBAL_DOMAIN,       &
                           mpp_global_sum,             &
                           BITWISE_EXACT_SUM

use  constants_mod, only: RADIUS
use        fms_mod, only: error_mesg, FATAL, NOTE,  &
                          mpp_pe, mpp_root_pe, mpp_npes, &
                          write_version_number, stdlog

implicit none
private

!-----------------------------------------------------------------------
!------- public interfaces -------

public  horiz_grid_init,       &
        get_horiz_grid_bound,  &
        get_horiz_grid_size,   &
        update_np, update_sp,  &
        TGRID, VGRID

!-----------------------------------------------------------------------
!             public derived data types
!-----------------------------------------------------------------------

public bgrid_type
public bgrid_interp_type
public horiz_grid_type

!-----------------------------------------------------------------------
!
!    NOTE: all horizontal indexing references global indices
!                    ( 1:nlon, 1:nlat )

type bgrid_type
   integer :: is,  ie,  js,  je         ! compute domain indices
   integer :: isd, ied, jsd, jed        ! data    domain indices
   integer :: isg, ieg, jsg, jeg        ! global  domain indices
   real, pointer, dimension(:)   :: blong=>NULL(), blatg=>NULL()  ! global grid edges
   real, pointer, dimension(:)   ::    dx=>NULL(),   rdx=>NULL(), &
                                     area=>NULL(), rarea=>NULL(), &
                                      tph=>NULL(),   tlm=>NULL()
   real, pointer, dimension(:,:) ::   aph=>NULL(),   alm=>NULL()
   real                          :: dy, rdy
   real                          :: areasum
   type(domain2D) :: Domain, Domain_nohalo
end type

!    is   = starting x-axis index for the compute domain
!    ie   = ending   x-axis index for the compute domain
!    js   = starting y-axis index for the compute domain
!    je   = ending   y-axis index for the compute domain
!    isd  = starting x-axis index for the data    domain
!    ied  = ending   x-axis index for the data    domain
!    jsd  = starting y-axis index for the data    domain
!    jed  = ending   y-axis index for the data    domain
!    isg  = starting x-axis index for the global  domain
!    ieg  = ending   x-axis index for the global  domain
!    jsg  = starting y-axis index for the global  domain
!    jeg  = ending   y-axis index for the global  domain
!
!    dx    = grid spacing for x-axis (in meters)
!    rdx   = reciprocal of dx (1/m)
!    dy    = grid spacing for y-axis (in meters)
!    rdy   = reciprocal of dy (1/m)
!    area  = area of a grid box (in m2)
!    rarea = reciprocal of area (1/m2)
!    areasum = bit-reproducible global sum of area
!
!    tph  = latitude at the center of grid box (in radians)
!    tlm  = longitude at the center of grid box (in radians)
!    aph  = actual latitude at the center of grid box (in radians)
!    alm  = actual longitude at the center of grid box (in radians)
!
!    blong = longitude grid box boundaries along the global x/longitude axis (in radians)
!    blatg = latitude  grid box boundaries along the global y/latitude  axis (in radians)
!
!    Domain        = domain2D variable with halo size = 1
!    Domain_nohalo = domain2D variable with halo size = 0
!                    used for outputing diagnostic fields
!
!-----------------------------------------------------------------------
! interpolation weights: used to move data between model grids

type bgrid_interp_type
   real, pointer, dimension(:,:) :: tmpwts=>NULL(), velwts=>NULL(), nowts=>NULL()
end type bgrid_interp_type

! tmpwts, velwts = interpolation weights for area-weighted 4-point averages
! nowts          = interpolation weights for simple 4-point averages

!-----------------------------------------------------------------------

type horiz_grid_type
   type(bgrid_type)        :: Tmp, Vel
   type(bgrid_interp_type) :: Interp
   integer :: nlon, nlat, isize, jsize
   integer :: ilb, iub, jlb, jub
   logical :: channel, double_periodic, decompx, decompy
   real    :: dlmd, dphd, dlm, dph
   real    :: sb, nb, wb, eb
   real, pointer, dimension(:,:) :: sinphv=>NULL(), tanphv=>NULL()
end type horiz_grid_type

!    Tmp = constants for the temperature/tracer/mass grid
!    Vel = constants for the u/v wind component grid
!
!    sinphv = sine of Vel%aph
!    tanphv = tangent of Vel%aph
!
!    nlon = number of grid points along the global x-axis (no halo points)
!    nlat = number of grid points along the global y-axis (no halo points)
!
!    isize = number of grid points along the x-axis for the current processor
!             (includes halo points)
!    jsize = number of grid points along the y-axis for the current processor
!             (includes halo points)
!
!    ilb  = lower bound x-axis
!    iub  = upper bound x-axis
!    jlb  = lower bound y-axis
!    jub  = upper bound y-axis
!
!    dlm  = grid spacing for x-axis (in radians)
!    dph  = grid spacing for y-axis (in radians)
!    dlmd = grid spacing for x-axis (in degrees of longitude)
!    dphd = grid spacing for y-axis (in degrees of latitude)
!
!    sb = southern boundary of temperature grid (in radians)
!    nb = northern boundary of temperature grid (in radians)
!    wb =  western boundary of temperature grid (in radians)
!    eb =  eastern boundary of temperature grid (in radians)
!
!    channel = channel model, uses grid boundaries: ed,wb,sb,nb.
!              also boundary condition at N/S walls is modified
!              may select f-plane option or double periodic options
!    double_periodic = f-plane channel with periodic boundaries in x and y.
!
!    decompx = x-axis is decomposed across more than one processor
!    decompy = y-axis is decomposed across more than one processor
!-----------------------------------------------------------------------
!-------- public parameters -------------

   integer, parameter :: TGRID = 51, VGRID = 52

!-----------------------------------------------------------------------
!-------- private data ------------

   real,    parameter :: eps=0.0001

!-------- internal parameters --------

!----- untested/unsupported options -----
!         use at you own risk 

   real :: wbd = 0.0  ! western edge in degrees of the first
                      ! longitude (i=1) of temperature grid boxes

   real :: ebd = 360.0 ! eastern edge in degrees of the last
                       ! longitude (i=nlon) of temperature grid boxes

   real :: sbd = -90.0  ! southern edge in degrees of the first
                        ! latitude row (j=1) of temperature grid boxes

   real :: nbd =  90.0  ! northern edge in degrees of the last latitude
                        ! row (j=nlat) of temperature grid boxes

   logical :: do_channel = .false. ! if TRUE, then use sbd,nbd as southern and northern 
                                   ! boundaries of a channel.

   logical :: do_fplane_approx = .false. ! f-plane channel without spherical geometry
                                         ! used in conjunction with do_channel option
   real    :: fphd = 45.                 ! coriolis latitude used for f-plane channel

   logical :: do_double_periodic = .false.  ! domain is cyclic in X and Y
                                            ! this option must be used with do_fplane_approx = true
                                            ! Tmp and Vel grids will have the same size global compute domain
!-----------------------------------------------------------------------

   real :: tph0d = 0.  ! grid/globe transformation (unsupported)
   real :: tlm0d = 0.  ! (set tph0d=tlm0d=0 for no transformation)

!-----------------------------------------------------------------------
 character(len=128) :: version = '$Id: bgrid_horiz.F90,v 11.0 2004/09/28 19:07:28 fms Exp $'
 character(len=128) :: tagname = '$Name: testing $'
 logical :: do_vers = .true.
!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine horiz_grid_init ( Hgrid, nlon, nlat, layout )

!-----------------------------------------------------------------------
!   Hgrid     = horizontal grid constants
!   nlon,nlat = the global horizontal grid resolution
!               number of longitude and latitude grid boxes, respectively
!   layout    = domain decomposition (num X pes by num Y pes),
!                default decomposition along y-axis then x-axis
!-----------------------------------------------------------------------

     type(horiz_grid_type), intent(inout) :: Hgrid
     integer,               intent(in)    :: nlon, nlat
     integer, optional,     intent(in)    :: layout(2)

!-----------------------------------------------------------------------
!------------------- local/private declarations ------------------------

real,    allocatable :: tlmi(:), tphj(:), slat(:), dxj(:), wt(:,:)
integer, allocatable :: xrows(:), yrows(:)

real    :: hpi, dtr
integer :: i, j, ilb, iub, jlb, jub, npes, pe, yflags, nlatv
integer :: is, ie, hs, he, vs, ve
integer :: isd, ied, hsd, hed, vsd, ved
integer :: isg, ieg, hsg, heg, vsg, veg
integer :: domain_layout(2)
logical :: global_x, global_y
type(domain1D) :: Domx, Domy

!------------- parallel interface --------------------------------------
!        ---- domain decomposition -----

      call mpp_domains_init

      npes = mpp_npes()

! write version info to logfile
     if (do_vers) then
       call write_version_number (version, tagname)
       do_vers = .false.
     endif

! error checks and messages

   Hgrid % channel = do_channel
   Hgrid % double_periodic = do_double_periodic

   if (Hgrid%channel) then
       write (stdlog(),'(a)') 'Channel model option has been selected.'
   endif

   if (do_fplane_approx) then
       if (.not.Hgrid%channel) call error_mesg ('horiz_grid_init', &
             'f-plane approximation cannot be used without do_channel=TRUE', FATAL)
       write (stdlog(),'(a)') '... the f-plane approximation will be used'
   endif

   if (do_double_periodic) then
       ! error check
       if (.not.do_fplane_approx) call error_mesg ('horiz_grid_init', &
              'double periodic option cannot be used without do_fplane_approx=TRUE', FATAL)
       yflags = CYCLIC_GLOBAL_DOMAIN
       nlatv  = nlat
       write (stdlog(),'(a)') '... and double periodic boundary conditions will be used'
   else
       yflags = 0
       nlatv  = nlat-1
   endif

! print channel options used
   if (Hgrid%channel) then
       write (stdlog(),'(4x,a4,f10.3,4x,a4,f10.3)') 'WBD=',wbd,'EBD=',ebd
       write (stdlog(),'(4x,a4,f10.3,4x,a4,f10.3)') 'SBD=',sbd,'NBD=',nbd
       if (do_fplane_approx) then
           write (stdlog(),'(4x,a4,f10.3)') 'FPHD=',fphd
       endif
   endif
 
!---- set-up x- & y-axis decomposition -----

 domain_layout = (/ 0, 0 /)
 if (present(layout)) domain_layout = layout
 if (domain_layout(1)+domain_layout(2) == 0) then
     call mpp_define_layout ( (/1,nlon,1,nlat/), npes, domain_layout )
 else
     if (domain_layout(1) == 0) domain_layout(1) = npes/domain_layout(2)
     if (domain_layout(2) == 0) domain_layout(2) = npes/domain_layout(1)
 endif

  if ( domain_layout(1)*domain_layout(2) /= npes ) call error_mesg &
              ('horiz_grid_init', 'number of processors requested not &
                                  &compatible with grid', FATAL )

! flag to indicate axis decomposition
  Hgrid % decompx = domain_layout(1) .gt. 1
  Hgrid % decompy = domain_layout(2) .gt. 1


!    ---- mass/temperature grid domain with halo = 0,1 ----

     call mpp_define_domains ( (/1,nlon,1,nlat/), domain_layout,   &
                               Hgrid % Tmp % Domain,               &
                               xflags = CYCLIC_GLOBAL_DOMAIN,      &
                               yflags = yflags,                    &
                               xhalo = 1, yhalo = 1,               &
          name = 'ATMOSPHERIC (B-GRID) MODEL, temperature grid, halo=1,')
     call mpp_define_domains ( (/1,nlon,1,nlat/), domain_layout,   &
                               Hgrid % Tmp % Domain_nohalo,        &
                               xflags = CYCLIC_GLOBAL_DOMAIN,      &
                               yflags = yflags,                    &
                               xhalo = 0, yhalo = 0                )

!   ---- compute exact decomposition ----
!   ---- compute 2d layout of PEs ----

     allocate ( xrows(domain_layout(1)), yrows(domain_layout(2)) )
     call mpp_get_domain_components ( Hgrid%Tmp%Domain, Domx, Domy )
     call mpp_get_compute_domains   ( Domx, size=xrows )
     call mpp_get_compute_domains   ( Domy, size=yrows )

!    ---- velocity grid may have one less latitude row ----

     if (nlatv == nlat-1) then
         yrows(domain_layout(2)) = yrows(domain_layout(2)) - 1
     endif

!    ---- velocity grid domain with halo = 0,1 ----

     call mpp_define_domains ( (/1,nlon,1,nlatv/), domain_layout,  &
                               Hgrid % Vel % Domain,               &
                               xflags = CYCLIC_GLOBAL_DOMAIN,      &
                               yflags = yflags,                    &
                               xhalo = 1, yhalo = 1,               &
                               xextent = xrows, yextent = yrows    )
     call mpp_define_domains ( (/1,nlon,1,nlatv/), domain_layout,  &
                               Hgrid % Vel % Domain_nohalo,        &
                               xflags = CYCLIC_GLOBAL_DOMAIN,      &
                               yflags = yflags,                    &
                               xhalo = 0, yhalo = 0,               &
                               xextent = xrows, yextent = yrows    )

     deallocate ( xrows, yrows )


!------------- indices for global compute domain -----------------------

     call mpp_get_global_domain ( Hgrid%Tmp%Domain, isg, ieg, hsg, heg )
     call mpp_get_global_domain ( Hgrid%Vel%Domain, isg, ieg, vsg, veg )

     Hgrid % Tmp % isg = isg;   Hgrid % Tmp % ieg = ieg
     Hgrid % Tmp % jsg = hsg;   Hgrid % Tmp % jeg = heg
     Hgrid % Vel % isg = isg;   Hgrid % Vel % ieg = ieg
     Hgrid % Vel % jsg = vsg;   Hgrid % Vel % jeg = veg

!------------- indices for data domain -----------------------

     call mpp_get_data_domain ( Hgrid%Tmp%Domain, isd, ied, hsd, hed )
     call mpp_get_data_domain ( Hgrid%Vel%Domain, isd, ied, vsd, ved )

     Hgrid % Tmp % isd = isd;   Hgrid % Tmp % ied = ied
     Hgrid % Tmp % jsd = hsd;   Hgrid % Tmp % jed = hed
     Hgrid % Vel % isd = isd;   Hgrid % Vel % ied = ied
     Hgrid % Vel % jsd = vsd;   Hgrid % Vel % jed = ved

!------------- indices for computational domain ------------------------

     call mpp_get_compute_domain ( Hgrid%Tmp%Domain, is, ie, hs, he )
     call mpp_get_compute_domain ( Hgrid%Vel%Domain, is, ie, vs, ve )

     Hgrid % Tmp % is = is;   Hgrid % Tmp % ie = ie
     Hgrid % Tmp % js = hs;   Hgrid % Tmp % je = he
     Hgrid % Vel % is = is;   Hgrid % Vel % ie = ie
     Hgrid % Vel % js = vs;   Hgrid % Vel % je = ve

!------------- indices including halo regions --------------------------

     call mpp_get_data_domain ( Hgrid%Tmp%Domain, ilb, iub, jlb, jub )

     Hgrid % ilb = ilb;  Hgrid % iub = iub
     Hgrid % jlb = jlb;  Hgrid % jub = jub

!------------ other values ---------------------------------------------

     Hgrid % nlon = nlon  ! global resolution
     Hgrid % nlat = nlat

     Hgrid % isize = Hgrid % iub - Hgrid % ilb + 1  ! array size for current PE
     Hgrid % jsize = Hgrid % jub - Hgrid % jlb + 1

!-----------------------------------------------------------------------
!    -------- allocate space for arrays ---------

    call alloc_array_space ( ilb, iub, jlb, jub, Hgrid%Tmp )
    call alloc_array_space ( ilb, iub, jlb, jub, Hgrid%Vel )
                 
    allocate ( Hgrid % sinphv(ilb:iub,jlb:jub), &
               Hgrid % tanphv(ilb:iub,jlb:jub)  )

    allocate ( slat (hs-1:he+2), wt (ilb:iub,jlb:jub) ) 

!-----------------------------------------------------------------------
!--------------derived geometrical constants----------------------------

      hpi = acos(0.0)
      dtr = hpi/90.

      if (Hgrid%channel) then
          Hgrid % dphd = (nbd-sbd)/real(nlat)
          Hgrid % dlmd = (ebd-wbd)/real(nlon)
          Hgrid % sb = sbd*dtr
          Hgrid % nb = nbd*dtr
          Hgrid % wb = wbd*dtr
          Hgrid % eb = ebd*dtr
          global_x = abs(ebd-wbd-360.).lt.eps
          global_y = abs(nbd-sbd-180.).lt.eps
      else
          ! global code used for bit-reprocibility
          Hgrid % dphd = 180./real(nlat)
          Hgrid % dlmd = 360./real(nlon)
          Hgrid % sb =  -hpi
          Hgrid % nb =  +hpi
          Hgrid % wb =   0.0
          Hgrid % eb = 360.0*dtr
          global_x = .true.
          global_y = .true.
      endif

      Hgrid % dlm = Hgrid % dlmd*dtr
      Hgrid % dph = Hgrid % dphd*dtr

!     --- dy is the same on both grids ---

      Hgrid % Tmp %  dy = RADIUS * Hgrid % dph
      Hgrid % Vel %  dy = RADIUS * Hgrid % dph
      Hgrid % Tmp % rdy = 1.0 / Hgrid % Tmp % dy
      Hgrid % Vel % rdy = 1.0 / Hgrid % Vel % dy

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---- get grid box boundaries for computational global domain ----

     Hgrid % Tmp % blong(isg-1:ieg+2) = get_grid (nlon+3, Hgrid%wb-Hgrid%dlm, &
                                                                   Hgrid%dlm  )
   ! global cyclic continuity (may not be necessary)
     if (global_x) then
        Hgrid % Tmp % blong(isg-1) = Hgrid % Tmp % blong(ieg)-4.*hpi
        Hgrid % Tmp % blong(ieg+1) = Hgrid % Tmp % blong(isg)+4.*hpi
        Hgrid % Tmp % blong(ieg+2) = Hgrid % Tmp % blong(isg+1)+4.*hpi
     endif

     Hgrid % Tmp % blatg(hsg+1:heg) = get_grid (nlat-1, Hgrid%sb+Hgrid%dph, &
                                                                 Hgrid%dph  )
   ! polar halos
     Hgrid % Tmp % blatg(heg+2) = Hgrid % Tmp % blatg(heg)
     Hgrid % Tmp % blatg(heg+1) = Hgrid % nb
     Hgrid % Tmp % blatg(hsg)   = Hgrid % sb
     Hgrid % Tmp % blatg(hsg-1) = Hgrid % Tmp % blatg(hsg+1)

     do i = is-1, ie+1
        Hgrid % Tmp % tlm(i) = 0.5*(Hgrid%Tmp%blong(i)+Hgrid%Tmp%blong(i+1))
     enddo

     do j = hs-1, he+1
        Hgrid % Tmp % tph(j) = 0.5*(Hgrid%Tmp%blatg(j)+Hgrid%Tmp%blatg(j+1))
     enddo

     if (do_fplane_approx) then
!    --- no spherical geometry ---
         Hgrid % Tmp % dx(hs-1:he+1) = RADIUS * Hgrid % dlm
     else
!    --- sphere ---
         slat(hs-1:he+2) = sin(Hgrid%Tmp%blatg(hs-1:he+2))
         Hgrid % Tmp % dx (hs-1:he+1) = RADIUS * Hgrid % dlm / Hgrid % dph *  &
                                        abs(slat(hs:he+2)-slat(hs-1:he+1))
     endif

!-------------initialize lat/lon at velocity points---------------------

     Hgrid % Vel % blong(isg-1:ieg+2) = get_grid (nlon+3, Hgrid%wb-0.5*Hgrid%dlm, &
                                                                       Hgrid%dlm  )
   ! global cyclic continuity
     if (global_x) then
        Hgrid % Vel % blong(isg-1) = Hgrid % Vel % blong(ieg)-4.*hpi
        Hgrid % Vel % blong(ieg+1) = Hgrid % Vel % blong(isg)+4.*hpi
        Hgrid % Vel % blong(ieg+2) = Hgrid % Vel % blong(isg+1)+4.*hpi
     endif

     Hgrid % Vel % blatg(vsg:veg+1) = get_grid (nlat, Hgrid%sb+0.5*Hgrid%dph, &
                                                                   Hgrid%dph  )
   ! polar halos
     Hgrid % Vel % blatg(vsg-1) = Hgrid % sb
     Hgrid % Vel % blatg(veg+2) = Hgrid % nb

     do i = is-1, ie+1
        Hgrid % Vel % tlm(i) = 0.5*(Hgrid%Vel%blong(i)+Hgrid%Vel%blong(i+1))
     enddo
     do j = vs-1, ve+1
        Hgrid % Vel % tph(j) = 0.5*(Hgrid%Vel%blatg(j)+Hgrid%Vel%blatg(j+1))
     enddo

     if (do_fplane_approx) then
!    --- no spherical geometry ---
         Hgrid % Vel % dx = RADIUS * Hgrid % dlm
     else
!    --- sphere ---
         slat(vs-1:ve+2) = sin(Hgrid%Vel%blatg(vs-1:ve+2))
         Hgrid % Vel % dx (vs-1:ve+1) = RADIUS * Hgrid % dlm / Hgrid % dph *  &
                                        (slat(vs:ve+2)-slat(vs-1:ve+1))
     endif
     Hgrid % Vel % rdx(vs-1:ve+1) = 1.0 / Hgrid % Vel % dx(vs-1:ve+1)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------------- grid box areas ------------------------------------------

      Hgrid % Tmp % area = Hgrid % Tmp % dx * Hgrid % Tmp % dy
      Hgrid % Vel % area = Hgrid % Vel % dx * Hgrid % Vel % dy

! Note: area of the pole on velocity grid is not important because:
!         1) pole is not in compute domain
!         2) current polar boundary condition is u=v=0
!         3) fluxes between pole and sub-pole row are set to zero

!--- reciprocal of area ----

      do j = jlb, jub
         if (Hgrid % Tmp % area(j) > 0.0) Hgrid % Tmp % rarea(j) = 1.0 / Hgrid % Tmp % area(j)
         if (Hgrid % Vel % area(j) > 0.0) Hgrid % Vel % rarea(j) = 1.0 / Hgrid % Vel % area(j)
      enddo

!--- compute bit-reproducible global sum of area ---

      do j = hsd, hed
         wt(:,j) = Hgrid%Tmp%area(j)
      enddo
      Hgrid % Tmp % areasum = mpp_global_sum (Hgrid%Tmp%Domain, wt, flags=BITWISE_EXACT_SUM)

      do j = vsd, ved  ! may not be full array extent
         wt(:,j) = Hgrid%Vel%area(j)
      enddo
      Hgrid % Vel % areasum = mpp_global_sum (Hgrid%Vel%Domain, wt(:,vsd:ved), flags=BITWISE_EXACT_SUM)

      deallocate ( slat, wt )


!--- initialization of weights for grid interpolation ---

     call bgrid_interp_init ( Hgrid )

!-----------------------------------------------------------------------
! unsupported option for transforming the position of poles
! tph0d,tlm0d are the lat,lon position in the transformed grid
! of the point [lat=0,lon=0] ..... I think ????

      if (tph0d > eps .or. tlm0d > eps) then
         ! compute "actual" lat/lon ( aph, aph )
         ! actual lat/lon would be used for the coriolis, radiation, and ...
         ! at temperature points
           call trans_latlon (tph0d,  tlm0d,                        &
                              Hgrid % Tmp % tlm, Hgrid % Tmp % tph, &
                              Hgrid % Tmp % alm, Hgrid % Tmp % aph)
         ! at velocity points
           call trans_latlon (tph0d,  tlm0d,                        &
                              Hgrid % Vel % tlm, Hgrid % Vel % tph, &
                              Hgrid % Vel % alm, Hgrid % Vel % aph)
      else
         ! no transformation
           do j = jlb, jub
           do i = ilb, iub
              Hgrid % Tmp % aph(i,j) = Hgrid % Tmp % tph(j)
              Hgrid % Vel % aph(i,j) = Hgrid % Vel % tph(j)
              Hgrid % Tmp % alm(i,j) = Hgrid % Tmp % tlm(i)
              Hgrid % Vel % alm(i,j) = Hgrid % Vel % tlm(i)
           enddo
           enddo
      endif

!------------- trigometric constants -----------------------------------

      if (do_fplane_approx) then
!     --- channel model ---
          Hgrid % sinphv = sin(fphd*dtr)
          Hgrid % tanphv = tan(fphd*dtr)
      else
!     --- sphere ---
          Hgrid % sinphv = sin(Hgrid % Vel % aph)
          Hgrid % tanphv = tan(Hgrid % Vel % aph)
      endif

!-----------------------------------------------------------------------

end subroutine horiz_grid_init

!##############################################################################
! initializes weights for interpolation between model grids

subroutine bgrid_interp_init ( Hgrid )
type(horiz_grid_type), intent(inout) :: Hgrid

real :: hpi
real :: sph  (1:2*Hgrid%nlat+1), area(-1:2*Hgrid%nlat+2)
real :: areah(Hgrid%Tmp%jsd:Hgrid%Tmp%jed), areav(Hgrid%Vel%jsd:Hgrid%Vel%jed)
integer :: i, j, k
! indexing for interpolation weights array
  integer, parameter :: SOUTH = 1, NORTH = 2

! allocate memory
 allocate ( Hgrid%Interp%tmpwts(Hgrid%jlb:Hgrid%jub,2), &
            Hgrid%Interp%velwts(Hgrid%jlb:Hgrid%jub,2), &
            Hgrid%Interp%nowts (Hgrid%jlb:Hgrid%jub,2)  )

! weights for simple 4-pt averages with equal weighting
  Hgrid%Interp%nowts = 0.25

! equal weighting for non-spherical grid
  if (do_fplane_approx) then
      Hgrid%Interp%tmpwts = 0.25
      Hgrid%Interp%velwts = 0.25
      return
  endif

! create weights that will conserve the quantity being interpolated

! create areas for grid boxes and interpolation
! need sin(lat) every 1/2 delta-lat
  hpi = acos(0.0)
  if (.not.Hgrid%channel) then
      ! original global code (used for bit-reproducibility)
      sph(1) = -1.; sph(2*Hgrid%nlat+1) = 1.
      do j = 2, 2*Hgrid%nlat
        sph(j) = sin(-hpi+real(j-1)*Hgrid%dph*0.5)
      enddo
  else
      sph(1) = sin(Hgrid%sb); sph(2*Hgrid%nlat+1) = sin(Hgrid%nb)
      do j = 2, 2*Hgrid%nlat
        sph(j) = sin(Hgrid%sb+real(j-1)*Hgrid%dph*0.5)
      enddo
  endif
! 1/4 area/radius^2 of grid boxes
  area(1:2*Hgrid%nlat) = 0.5*Hgrid%dlm*(sph(2:2*Hgrid%nlat+1)-sph(1:2*Hgrid%nlat))
! assume halo size of one (fill two 1/2 size grid boxes)
  area(-1:0) = area(2:1:-1)
  area(2*Hgrid%nlat+1:2*Hgrid%nlat+2) = area(2*Hgrid%nlat:2*Hgrid%nlat-1:-1)

! area of temperature grid boxes
  do j = Hgrid%Tmp%jsd, Hgrid%Tmp%jed
    Hgrid%Interp%tmpwts(j,SOUTH) = area(2*j-1)
    Hgrid%Interp%tmpwts(j,NORTH) = area(2*j)
    areah(j) = (Hgrid%Interp%tmpwts(j,SOUTH)+Hgrid%Interp%tmpwts(j,NORTH))*2.0
    ! compute final weights for 4-pt averages
    Hgrid%Interp%tmpwts(j,SOUTH) = Hgrid%Interp%tmpwts(j,SOUTH) / areah(j)
    Hgrid%Interp%tmpwts(j,NORTH) = Hgrid%Interp%tmpwts(j,NORTH) / areah(j)
  enddo

! area of velocity grid boxes
  do j = Hgrid%Vel%jsd, Hgrid%Vel%jed
    Hgrid%Interp%velwts(j,SOUTH) = area(2*j)
    Hgrid%Interp%velwts(j,NORTH) = area(2*j+1)
    areav(j) = (Hgrid%Interp%velwts(j,SOUTH)+Hgrid%Interp%velwts(j,NORTH))*2.0
    ! compute final weights for 4-pt averages
    Hgrid%Interp%velwts(j,SOUTH) = Hgrid%Interp%velwts(j,SOUTH) / areav(j)
    Hgrid%Interp%velwts(j,NORTH) = Hgrid%Interp%velwts(j,NORTH) / areav(j)
  enddo

end subroutine bgrid_interp_init

!#######################################################################

function get_grid (npts, start, space) result (grid)

integer, intent(in) :: npts
real,    intent(in) :: start, space
real                :: grid(npts)
integer :: j

!---- compute equally spaced grid ----

      do j = 1, npts
         grid(j) = start + real(j-1)*space
      enddo

end function get_grid

!#######################################################################

   subroutine trans_latlon (tph0d, tlm0d, tlm, tph, alm, aph)

   real, intent(in)  :: tph0d, tlm0d, tlm(:), tph(:)
   real, intent(out) :: alm(:,:), aph(:,:)

   real :: dtr, pi, tph0, tlm0, stph0, ctph0, ttph0, cc, ee
   real, dimension(size(tph)) :: stph, ctph
   integer :: i, j

!-----------------------------------------------------------------------

    ! scalars
      pi  = 2.*acos(0.0)
      dtr  = pi/180.
      tph0 = tph0d*dtr;  tlm0 = tlm0d*dtr
      stph0 = sin(tph0); ctph0 = cos(tph0)
      ttph0 = stph0/ctph0

!     ------- compute actual lat and lon -------
      stph = sin(tph)
      ctph = cos(tph)
      do j = 1, size(aph,2)
      do i = 1, size(aph,1)
         cc = ctph(j)*cos(tlm(i))
         aph(i,j) = asin(ctph0*stph(j)+stph0*cc)
         ee = cc/(ctph0*cos(aph(i,j))) - tan(aph(i,j))*ttph0
         ee = min(ee,1.)
         if (tlm(i) > pi) then
             alm(i,j) = tlm0-acos(ee)
         else
             alm(i,j) = tlm0+acos(ee)
         endif
      enddo
      enddo

!-----------------------------------------------------------------------

end subroutine trans_latlon

!#######################################################################

subroutine get_horiz_grid_bound ( Hgrid, grid, blon, blat, global )

!       returns the grid box boundaries for either
!           the compute or global grid
!
!   Hgrid = horizontal grid constants
!   grid  = grid identifier, possible values: TGRID, VGRID
!   blon  = longitude edges in radians
!   blat  = latitude edges in radians
!   global = values for compute(F) or global(T) grid?

   type (horiz_grid_type), intent(in)  :: Hgrid
   integer,                intent(in)  ::  grid
   real,                   intent(out) :: blon(:), blat(:)
   logical, optional,      intent(in)  :: global

    select case (grid)
       case (TGRID)
          call horiz_grid_bound ( Hgrid%Tmp, blon, blat, global )
       case (VGRID)
          call horiz_grid_bound ( Hgrid%Vel, blon, blat, global )
       case default
          call error_mesg ('get_horiz_grid_bound', 'invalid grid', FATAL)
    end select

!-----------------------------------------------------------------------

end subroutine get_horiz_grid_bound

!#######################################################################

subroutine get_horiz_grid_size ( Hgrid, grid, nlon, nlat, global )

!     returns the number of longitude and latitude grid boxes
!            for either the compute or global grid
!
!   Hgrid = horizontal grid constants
!   grid  = grid identifier, possible values: TGRID, VGRID
!   nlon  = number longitude grid boxes
!   nlat  = number latitude grid boxes
!   global = values for compute(F) or global(T) grid?

   type (horiz_grid_type), intent(in)  :: Hgrid
   integer,                intent(in)  ::  grid
   integer,                intent(out) :: nlon, nlat
   logical, optional,      intent(in)  :: global

    select case (grid)
       case (TGRID)
          call horiz_grid_size ( Hgrid%Tmp, nlon, nlat, global )
       case (VGRID)
          call horiz_grid_size ( Hgrid%Vel, nlon, nlat, global )
       case default
          call error_mesg ('get_horiz_grid_size', 'invalid grid', FATAL)
    end select

end subroutine get_horiz_grid_size

!#######################################################################

subroutine horiz_grid_bound ( Grid, blon, blat, global )

   type (bgrid_type), intent(in)  :: Grid
   real,              intent(out) :: blon(:), blat(:) 
   logical, optional, intent(in)  :: global

!      private routine that returns the grid box boundaries
!          for either the compute or global grid

    integer :: is, ie, js, je
    logical :: lglobal 

    lglobal = .false.;  if (present(global)) lglobal = global


    if (lglobal) then
      is = Grid % isg; ie = Grid % ieg  ! global grid
      js = Grid % jsg; je = Grid % jeg
    else    
      is = Grid % is ; ie = Grid % ie   ! compute grid
      js = Grid % js ; je = Grid % je
    endif   

     !----- define longitudinal grid box edges -----

      if (size(blon) /= ie-is+2) call error_mesg  &
                         ('get_horiz_grid_bound', &
                          'invalid argument dimension for blon', FATAL)

      blon = Grid % blong (is:ie+1)

     !----- define latitudinal grid box edges -----

      if (size(blat) /= je-js+2) call error_mesg  &
                         ('get_horiz_grid_bound', &
                          'invalid argument dimension for blat', FATAL)

      blat = Grid %  blatg (js:je+1)

!-----------------------------------------------------------------------

end subroutine horiz_grid_bound

!#######################################################################

subroutine horiz_grid_size ( Grid, nlon, nlat, global )

   type (bgrid_type), intent(in)  :: Grid
   integer,                intent(out) :: nlon, nlat
   logical, optional,      intent(in)  :: global

   logical :: lglobal

   lglobal = .false.;  if (present(global)) lglobal = global

  !---- return the size of the requested grid ----

   if (lglobal) then
      nlon = Grid % ieg - Grid % isg + 1   ! global grid
      nlat = Grid % jeg - Grid % jsg + 1
   else
      nlon = Grid % ie - Grid % is + 1     ! compute grid
      nlat = Grid % je - Grid % js + 1
   endif

end subroutine horiz_grid_size

!#######################################################################

function update_np (Hgrid, grid) result (answer)

!   Hgrid = horizontal grid constants
!   grid  = grid identifier, possible values: TGRID, VGRID
!
!   Returns TRUE if the northernmost latitude row
!   is adjacent to the north pole.

  type (horiz_grid_type), intent(in)  :: Hgrid
  integer,                intent(in)  ::  grid
  logical :: answer
  integer :: je, jeg

     answer = .false.
    !if (Hgrid%Tmp%jeg == Hgrid%Vel%jeg) return ! periodic in y
     if (Hgrid%double_periodic) return ! periodic in y

     select case (grid)
       case (TGRID)
          if ( Hgrid%Tmp%je == Hgrid%Tmp%jeg ) answer = .true.
       case(VGRID)
          if ( Hgrid%Vel%je == Hgrid%Vel%jeg ) answer = .true.
       case default
          call error_mesg ('update_np', 'invalid grid', FATAL)
     end select

end function update_np

!#######################################################################

function update_sp (Hgrid, grid) result (answer)

!   Hgrid = horizontal grid constants
!   grid  = grid identifier, possible values: TGRID, VGRID
!
!   Returns TRUE if the southernmost latitude row
!   is adjacent to the south pole.

  type (horiz_grid_type), intent(in)  :: Hgrid
  integer,                intent(in)  ::  grid
  logical :: answer
  integer :: js, jsg

     answer = .false.
    !if (Hgrid%Tmp%jeg == Hgrid%Vel%jeg) return ! periodic in y
     if (Hgrid%double_periodic) return ! periodic in y

     select case (grid)
       case (TGRID)
          if ( Hgrid%Tmp%js == Hgrid%Tmp%jsg ) answer = .true.
       case(VGRID)
          if ( Hgrid%Vel%js == Hgrid%Vel%jsg ) answer = .true.
       case default
          call error_mesg ('update_sp', 'invalid grid', FATAL)
     end select

end function update_sp

!#######################################################################

 subroutine alloc_array_space ( ilb, iub, jlb, jub, Grid )
 integer, intent(in) :: ilb, iub, jlb, jub
 type(bgrid_type), intent(inout) :: Grid
         
      allocate ( Grid % dx   (jlb:jub), &
                 Grid % rdx  (jlb:jub), &
                 Grid % tph  (jlb:jub), &
                 Grid % tlm  (ilb:iub), &
                 Grid % area (jlb:jub), &
                 Grid % rarea(jlb:jub), &
                 Grid % aph  (ilb:iub,jlb:jub), &
                 Grid % alm  (ilb:iub,jlb:jub)  )

      allocate ( Grid % blong (Grid%isg-1:Grid%ieg+2), &
                 Grid % blatg (Grid%jsg-1:Grid%jeg+2)  )

      Grid % dx   = 0.; Grid % rdx   = 0.
      Grid % tph  = 0.; Grid % tlm   = 0.
      Grid % area = 0.; Grid % rarea = 0.

 end subroutine alloc_array_space

!#######################################################################

end module bgrid_horiz_mod

