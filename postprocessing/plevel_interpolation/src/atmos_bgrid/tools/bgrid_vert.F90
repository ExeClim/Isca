
module bgrid_vert_mod

!-----------------------------------------------------------------------
!
!     allocates memory and initializes vertical grid constants
!
!     contains interfaces for computing pressure and height
!
!-----------------------------------------------------------------------

use constants_mod, only: GRAV, RDGAS, RVGAS
use       fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL,       &
                         write_version_number, stdlog,                 &
                         mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                         MPP_CLOCK_SYNC, CLOCK_ROUTINE, CLOCK_LOOP

implicit none
private

!-----------------------------------------------------------------------
!  public derived data type (vert_grid_type)
!  -----------------------------------------
!    nlev  = number of vertical levels (integer)
!    nplev = number of pure pressure levels at the top of the model,
!            will equal zero when not using hybrid coordinate (integer)
!
!    deta  = sigma (aka eta) thicknesses of model layers
!    eta   = sigma (aka eta, bk) values at model layer interfaces (half levels)
!    peta  = reference pressure (pk) values at model layer interfaces
!    dpeta = reference pressure thicknesses of model layers
!
!    **** Note: pfull, phalf, fhalf are based on a SLP=101325. ****
!    pfull = pressure profile at model levels (full levels)
!    phalf = pressure profile at model layer interfaces (half levels)
!    fhalf = profile of geopotental height at half levels
!
!    wta,wtb = weights used to determine values at model levels
!                (based on pressure profile, useful for hybrid coord)
!    
!    psmin  = minimum allowable surface pressure to avoid negative
!             mass in a model layer (important for hybrid coord)
!    hybrid = logical flag (true for hybrid coordinate)
!    pzero  = logical flag (true for pres = 0 at top of model)
!    pref,tref,gamma = reference values for computing fhalf
!
!
!  public interfaces
!  -----------------
!    vert_grid_init        - initializes vert_grid_type data
!
!    compute_pres_depth    - computes pressure depth (mass) of model layers
!    compute_pres_full     - computes pressure at full model levels
!    compute_pres_half     - computes pressure at half model levels
!    compute_height        - computes geopotential height (in meters) at half model levels
!    compute_geop_height   - computes geopotential height (in m2/s2)
!    compute_pres_weights  - computes wta & wtb
!    compute_pressure      - computes pressure at full and half levels
!    compute_height_bottom - computes height at the lowest model level
!
!-----------------------------------------------------------------------
!------- interfaces -------

public  vert_grid_init, compute_geop_height, compute_height,   &
        compute_pres_depth, compute_pres_full, compute_pres_half, &
        compute_pres_weights, compute_pressures, compute_height_bottom

!------- public defined data type -------

public  vert_grid_type

type vert_grid_type
   integer                     :: nlev, nplev
   real, pointer, dimension(:) :: deta =>NULL(),  &
                                  eta =>NULL(),  &
                                  fhalf =>NULL(),  &
                                  dpeta =>NULL(), &
                                  pfull =>NULL(), &
                                  phalf =>NULL(), &
                                  peta =>NULL(), &
                                  wta =>NULL(), &
                                  wtb =>NULL()
   real                        :: pref, tref, gamma, psmin
   logical                     :: hybrid, pzero
end type vert_grid_type

!-----------------------------------------------------------------------

  real, parameter :: d608 = (RVGAS-RDGAS)/RDGAS
  real, parameter :: ginv = 1./GRAV

!------ parameters for eta coordinate reference surface heights --------

  real, parameter :: pref = 101325., tref = 288., gamma = 0.0065

!------ parameters for performance timing of code sections -----

  logical :: do_clock_init = .true.
  integer, dimension(7) :: id
  character(len=16), dimension(7) :: &
            names = (/ 'comp_pres_depth ', &
                       'comp_pres_full  ', &
                       'comp_pressures  ', &
                       'comp_pres_half  ', &
                       'comp_pres_wghts ', &
                       'comp_geop_hght  ', &
                       'comp_height_btm ' /)

!-----------------------------------------------------------------------
 character(len=128) :: version =  '$Id: bgrid_vert.F90,v 11.0 2004/09/28 19:07:49 fms Exp $'
 character(len=128) :: tagname =  '$Name: testing $'
!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine vert_grid_init ( Vgrid, eta, peta, verbose )

!-----------------------------------------------------------------------
!   Vgrid = vertical grid constants
!   eta   = sigma (aka eta,bk) values at model layer interfaces
!             the number of model levels will be "size(eta)-1"
!   peta  = reference pressures (aka pk) at model layer interfaces
!             if specified and non-zero then hybrid coord is used
!             default: peta=0.
!   verbose = controls the amount of printed output
!               possible values are verbose=0,1,2
!               default: verbose=1
!-----------------------------------------------------------------------

  type(vert_grid_type), intent(inout) :: Vgrid
  real,    intent (in)                ::  eta(:)
  real,    intent (in), optional      :: peta(:)
  integer, intent (in), optional      :: verbose

!-----------------------------------------------------------------------
integer :: k, nlev, lverbose
real    :: rgog
real, dimension(size(eta(:)))   :: lphalf, pres
real, dimension(size(eta(:))-1) :: lpfull

!-----------------------------------------------------------------------

   lverbose = 1;  if (present(verbose)) lverbose = verbose

!--------------derived vertical constants-------------------------------

      nlev = size(eta(:)) - 1
      allocate (Vgrid% deta(nlev), Vgrid% eta(nlev+1),  &
                Vgrid%dpeta(nlev), Vgrid%pfull(nlev), Vgrid%peta(nlev+1),  &
                Vgrid%wta  (nlev), Vgrid%wtb  (nlev), &
                Vgrid%phalf(nlev+1), Vgrid% fhalf(nlev+1))

      Vgrid % nlev = nlev

!--------- set-up eta values and hybrid pressure levels ----------
!--------- note: eta(1) and eta(nlev+1) have set values -----
!--------- also note: peta(nlev+1) = 0.0 -----

        Vgrid % eta(1)      = 0.0
        Vgrid % eta(nlev+1) = 1.0
        Vgrid % eta(2:nlev) = eta(2:nlev)

        Vgrid % peta = 0.0
        if (present(peta)) Vgrid % peta(1:nlev) = peta(1:nlev)

      do k = 1, nlev
        Vgrid %  deta(k) =  Vgrid %  eta(k+1) - Vgrid %  eta(k)
        Vgrid % dpeta(k) =  Vgrid % peta(k+1) - Vgrid % peta(k)
        Vgrid % pfull(k) =  0.0
        Vgrid % wta  (k) =  0.0
        Vgrid % wtb  (k) =  0.0
      enddo

!----------- is this a hybrid coordinate ??? -----

      Vgrid % hybrid = .false.

      do k = 1, nlev+1
         if ( Vgrid % peta(k) > 0.0 ) then
              Vgrid % hybrid = .true.
              exit
         endif
      enddo

!----------- find lowest pure pressure level --------------

      Vgrid % nplev = 0

      do k = 1, nlev
         if ( Vgrid % deta(k) > 0.0 ) exit
         Vgrid % nplev = k
      enddo

!     ---- need average pressure in these layers ----

      Vgrid % pzero = .true.

     !if ( Vgrid % nplev >= 1 ) then
         Vgrid % phalf(:) = Vgrid%peta(:) + Vgrid%eta(:)*pref
         if ( Vgrid % phalf(1) <= epsilon(Vgrid % phalf) ) then
            lphalf(1) = 0.0
            lphalf(2:) = log(Vgrid % phalf(2:))
         else
            lphalf(:) = log(Vgrid % phalf(:))
            Vgrid % pzero = .false.
         endif

         do k = 1, nlev
            lpfull(k) = (Vgrid%phalf(k+1)*lphalf(k+1) - Vgrid%phalf(k)*lphalf(k)) &
                        / (Vgrid%phalf(k+1)-Vgrid%phalf(k)) - 1.0
            Vgrid % pfull(k) = exp(lpfull(k))
            Vgrid % wtb  (k) = lphalf(k+1) - lpfull(k)
            Vgrid % wta  (k) = lpfull(k)   - lphalf(k)
         enddo
         if (Vgrid % pzero) Vgrid % wta(1) =Vgrid % wtb(1)

     !endif

!----------- find the minimum allowable surface pressure ------

      Vgrid % psmin = 0.0

      do k = 1, nlev
         if ( Vgrid % deta(k) > 0.0 )  Vgrid % psmin =  &
               max ( Vgrid % psmin, -Vgrid % dpeta(k)/Vgrid % deta(k) )
      enddo


!---------- set-up eta coordinate geopotential heights -------------

      rgog = RDGAS*gamma/GRAV

      do k=1,nlev
        pres(k) = Vgrid%peta(k) + Vgrid%eta(k)*pref
        Vgrid % fhalf(k) = GRAV*tref*(1.0-(pres(k)/pref)**rgog)/gamma
      enddo
        Vgrid % fhalf(nlev+1) = 0.0

        Vgrid % pref  = pref
        Vgrid % tref  = tref
        Vgrid % gamma = gamma

!--- write version info to logfile ---
   call write_version_number (version, tagname)

!--- optional output of coordinate values ----
   if (mpp_pe() == mpp_root_pe()) then
       if (lverbose > 0) then
           write (stdlog(),*) 'Number of vertical levels =', Vgrid%nlev
           write (stdlog(),*) 'Approxiamte model level locations in pascals:'
           write (stdlog(),*) '  Full levels = ', Vgrid%pfull
           write (stdlog(),*) '  Half levels = ', Vgrid%phalf
           if (Vgrid%hybrid) &
           write (stdlog(),*) 'Hybrid coordinate with minimum allowable surface pressure =', Vgrid%psmin
       endif
       if (lverbose > 1) then
           write (stdlog(),*) 'INPUT VALUES:'
           write (stdlog(),*) '  eta =', Vgrid%eta
           write (stdlog(),*) '  peta =', Vgrid%peta
           write (stdlog(),*) 'OTHER VALUES:'
           write (stdlog(),*) '  deta =', Vgrid%deta
           write (stdlog(),*) '  dpeta =', Vgrid%dpeta
           write (stdlog(),*) '  nplev =', Vgrid%nplev
       endif
   endif

! initialize code sections for performance timing 
   if (do_clock_init) then
      do k = 1, size(id(:))
         id(k) = mpp_clock_id ('BGRID: vert ('//trim(names(k))//')', &
                               flags=MPP_CLOCK_SYNC, grain=CLOCK_LOOP)
      enddo
      do_clock_init = .false. 
   endif

!-----------------------------------------------------------------------

end subroutine vert_grid_init

!#######################################################################

subroutine compute_pres_depth (Vgrid, pssl, pdepth)

!--------------------------------------------------------
! Compute the pressure depth (mass) of model layers
!
!   Vgrid = vertical grid constants
!   pssl   = pressure at eta=1
!   pdepth = pressure depth of model layers
!--------------------------------------------------------

  type(vert_grid_type), intent(in)  :: Vgrid
       real           , intent(in)  :: pssl(:,:)
       real           , intent(out) :: pdepth(:,:,:)

       integer :: k, kp, ke

!-----------------------------------------------------------------------
  call mpp_clock_begin (id(1))

  kp = Vgrid % nplev
  ke = Vgrid % nlev

  if (size(pdepth,3) /= ke) call error_mesg (  &
                              'compute_pres_depth in bgrid_vert_mod', &
                              'incorrect dimension 3 for pdepth', FATAL)

! --- check for zero/negative depth layers ---

  if (Vgrid % hybrid) then
      if (minval(pssl) <= Vgrid % psmin) call error_mesg  &
                 ('compute_pres_depth in bgrid_vert_mod',  &
                  'pressure depth <= 0.0', FATAL)
  endif

! --- compute depth ---

  do k = 1, kp
    pdepth(:,:,k) = Vgrid % dpeta(k)
  enddo
    
  do k = kp+1, ke
    pdepth(:,:,k) = Vgrid % dpeta(k) + Vgrid % deta(k) * pssl(:,:)
  enddo

  call mpp_clock_end (id(1))

end subroutine compute_pres_depth

!#######################################################################

subroutine compute_pres_full (Vgrid, pssl, pfull, phalf, dpde)

!--------------------------------------------------------
! Compute the pressure at full model levels (plus options)
!
!   Vgrid = vertical grid constants
!   pssl  = pressure at eta=1
!   pfull = pressure at model levels
!   phalf = pressure at model layer interfaces (half levels)
!   dpde  = pressure depth of model layers
!--------------------------------------------------------

  type(vert_grid_type), intent(in)  :: Vgrid
  real                , intent(in)  :: pssl(:,:)
  real                , intent(out) :: pfull(:,:,:)
  real, optional      , intent(in)  :: phalf(:,:,:), dpde(:,:,:)

   real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)+1) :: ph
   real, dimension(size(pfull,1),size(pfull,2),size(pfull,3))   :: dp
   integer :: k, kp, ke
!-----------------------------------------------------------------------
!      compute the pressure at full model levels
!-----------------------------------------------------------------------
  call mpp_clock_begin (id(2))

  kp = Vgrid % nplev
  ke = Vgrid % nlev

  if (size(pfull,3) /= ke) call error_mesg (  &
                               'compute_pres_full in bgrid_vert_mod', &
                               'incorrect dimension 3 for pfull', FATAL)

!--- set or compute optional arguments ---

  if (present(phalf)) then
      ph = phalf
  else
      call compute_pres_half (Vgrid, pssl, ph)
  endif

  if (present(dpde)) then
      dp = dpde
  else
      call compute_pres_depth (Vgrid, pssl, dp)
  endif

!--- compute p*logp at half levels ---

  if ( Vgrid % pzero ) then
       ph(:,:,1) = 0.0
       ph(:,:,2:ke+1) = ph(:,:,2:ke+1) * log(ph(:,:,2:ke+1))
  else
       ph(:,:,:) = ph(:,:,:) * log(ph(:,:,:))
  endif

!--- compute pressure at full levels ---

  do k = 1, kp
    pfull(:,:,k) = Vgrid % pfull(k)
  enddo

  do k = kp+1, ke
    pfull(:,:,k) = exp( (ph(:,:,k+1)-ph(:,:,k))/dp(:,:,k) - 1.0 )
  enddo
    
  call mpp_clock_end (id(2))

end subroutine compute_pres_full

!#######################################################################

subroutine compute_pressures (Vgrid, pssl, phalf, pfull, dpde, wta, wtb)

!--------------------------------------------------------
! Compute the pressures, layer depths, and weights
!
!   Vgrid   = vertical grid constants
!   pssl    = pressure at eta=1
!   phalf   = pressure at model layer interfaces (half levels)
!   pfull   = pressure at model levels (full levels)
!   dpde    = pressure depth (mass) of model layers
!   wta,wtb = weights for computing data at model levels
!--------------------------------------------------------

  type(vert_grid_type), intent(in)  :: Vgrid
  real                , intent(in)  :: pssl(:,:)
  real                , intent(out) :: phalf(:,:,:), pfull(:,:,:)
  real, optional      , intent(out) :: dpde(:,:,:)
  real, optional      , intent(out) :: wta(:,:,:), wtb(:,:,:)

  real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)+1) :: ph,lph
  real, dimension(size(pfull,1),size(pfull,2),size(pfull,3))   :: dp,lpf
  integer :: k, kp, ke
!-----------------------------------------------------------------------

  kp = Vgrid % nplev
  ke = Vgrid % nlev

  if (size(pfull,3) /= ke) call error_mesg (  &
                               'compute_pressures in bgrid_vert_mod', &
                               'incorrect dimension 3 for pfull', FATAL)

!--- set or compute optional arguments ---

  call compute_pres_half (Vgrid, pssl, ph)
  phalf = ph

  call compute_pres_depth (Vgrid, pssl, dp)
  if (present(dpde))  dpde = dp

! do not include time for previous calls
  call mpp_clock_begin (id(3))

!--- compute p*logp at half levels ---

  if ( Vgrid % pzero ) then
      lph(:,:,1) = 0.0
       ph(:,:,1) = 0.0
      lph(:,:,2:ke+1) = log(ph(:,:,2:ke+1))
       ph(:,:,2:ke+1) = ph(:,:,2:ke+1) * lph(:,:,2:ke+1)
  else
      lph(:,:,:) = log(ph(:,:,:))
       ph(:,:,:) = ph(:,:,:) * lph(:,:,:)
  endif

!--- compute pressure at full levels ---

  do k = 1, kp
    pfull(:,:,k) = Vgrid % pfull(k)
  enddo

! if (present(wta) .or. present(wtb)) then
!   do k = 1, kp
!     lpf(:,:,k) = log(pfull(:,:,k))
!   enddo
! endif

  do k = kp+1, ke
    lpf(:,:,k) = (ph(:,:,k+1)-ph(:,:,k))/dp(:,:,k) - 1.0
    pfull(:,:,k) = exp( lpf(:,:,k) )
  enddo
    
!--- compute weights at full levels ---

   if (present(wtb)) then
      do k = 1, kp
        wtb(:,:,k) = Vgrid % wtb(k)
      enddo
      do k = kp+1, size(wtb,3)
        wtb(:,:,k) = lph(:,:,k+1) - lpf(:,:,k)
      enddo
   endif

   if (present(wta)) then
      do k = 1, kp
        wta(:,:,k) = Vgrid % wta(k)
      enddo
      do k = kp+1, size(wta,3)
        wta(:,:,k) = lpf(:,:,k) - lph(:,:,k)
      enddo
      if (Vgrid % pzero .and. kp == 0) wta(:,:,1) = wtb(:,:,1)
   endif

  call mpp_clock_end (id(3))

end subroutine compute_pressures

!#######################################################################

subroutine compute_pres_half (Vgrid, pssl, phalf)

!-------------------------------------------------------------
! Compute the pressure at the interface between model layers
!
!   Vgrid = vertical grid constants
!   pssl  = pressure at eta=1
!   phalf = pressure at model layer interfaces (half levels)
!-------------------------------------------------------------

  type(vert_grid_type), intent(in)  :: Vgrid
       real           , intent(in)  :: pssl(:,:)
       real           , intent(out) :: phalf(:,:,:)

       integer :: k, kp, ke

!-----------------------------------------------------------------------
  call mpp_clock_begin (id(4))

  kp = Vgrid % nplev + 1
  ke = Vgrid % nlev + 1

  if (size(phalf,3) /= ke) call error_mesg (  &
                               'compute_pres_half in bgrid_vert_mod', &
                               'incorrect dimension 3 for phalf', FATAL)

! pure pressure layers
  do k = 1, kp
    phalf(:,:,k) = Vgrid % peta(k)
  enddo
    
! sigma/pressure layers
  do k = kp+1, ke
    phalf(:,:,k) = Vgrid % peta(k) + Vgrid % eta(k) * pssl(:,:)
  enddo

  call mpp_clock_end (id(4))

end subroutine compute_pres_half

!#######################################################################

subroutine compute_pres_weights ( Vgrid, phalf, pfull, wta, wtb )

!------------------------------------------------------------
! Compute the weights for determining data at model levels
!
!   Vgrid   = vertical grid constants
!   phalf   = pressure at model layer interfaces (half levels)
!   pfull   = pressure at model levels (full levels)
!   wta,wtb = weights for computing data at model levels
!------------------------------------------------------------

 type(vert_grid_type), intent(in)  :: Vgrid
 real, intent(in)  :: phalf(:,:,:), pfull(:,:,:)
 real, intent(out) :: wta(:,:,:), wtb(:,:,:)

 real, dimension(size(phalf,1),size(phalf,2),size(phalf,3)) :: logph
 real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: logpf
 integer :: k, kp, kx, ks

  call mpp_clock_begin (id(5))

   kp = Vgrid % nplev
   kx = size(pfull,3)

   ! start indexing at 2 if ptop=0
   if (Vgrid%pzero) then
       ks = max(2,kp+1)
   else
       ks = kp+1
   endif

   logph(:,:,ks  :kx+1) = log(phalf(:,:,ks  :kx+1))
   logpf(:,:,kp+1:kx  ) = log(pfull(:,:,kp+1:kx  ))

 ! weights for half level below
   do k = 1, kp
     wtb(:,:,k) = Vgrid % wtb(k)
   enddo
   do k = kp+1, kx
     wtb(:,:,k) = logph(:,:,k+1) - logpf(:,:,k)
   enddo

 ! weights for half level above
   if (Vgrid%pzero .and. kp == 0) wta(:,:,1) = wtb(:,:,1)
   do k = 1, kp
     wta(:,:,k) = Vgrid % wta(k)
   enddo
   do k = ks, kx
     wta(:,:,k) = logpf(:,:,k) - logph(:,:,k)
   enddo

  call mpp_clock_end (id(5))

end subroutine compute_pres_weights

!#######################################################################

subroutine compute_geop_height (Vgrid, fssl, vtemp, wta, wtb, &
                                zfull, zhalf, mask)

!-----------------------------------------------------------------
! Compute the geopotential height at full (and half) model levels
!
!   Vgrid   = vertical grid constants
!   fssl    = geopotential height at eta=1
!   vtemp   = layer mean virtual temperature (at full model levels)
!   wta,wtb = weights for computing data at full model levels
!   zfull   = geop height (m2/s2) at model levels (full levels)
!   zhalf   = geop height (m2/s2) at model layer interfaces (half levels)
!   mask    = grid box mask for eta coordinate topography
!------------------------------------------------------------


    type(vert_grid_type), intent(in)   :: Vgrid
    real, intent(in), dimension(:,:)   :: fssl
    real, intent(in), dimension(:,:,:) :: vtemp, wta, wtb
    real, intent(out)                  :: zfull(:,:,:)
    real, intent(out), optional        :: zhalf(:,:,:)
    real, intent(in),  optional        :: mask(:,:,:)

       integer :: k, klev
       real, dimension(size(vtemp,1),size(vtemp,2)) :: zb, zt, rt

!-----------------------------------------------------------------------
  call mpp_clock_begin (id(6))

  klev = Vgrid % nlev

  if (size(zfull,3) /= klev) call error_mesg (  &
                               'compute_geop_height in bgrid_vert_mod', &
                               'incorrect dimension 3 for zfull', FATAL)

  if (present(zhalf)) then
     if (size(zhalf,3) /= klev+1) call error_mesg (  &
                               'compute_geop_height in bgrid_vert_mod', &
                               'incorrect dimension 3 for zhalf', FATAL)
  endif

      zb(:,:) = fssl(:,:)
      if (present(zhalf)) zhalf(:,:,klev+1) = zb(:,:)

!------- vertical integration loop (bottom to top) ----------

      do k = klev, 1, -1

         rt(:,:) = RDGAS * vtemp(:,:,k)
         zt(:,:) = zb(:,:) + rt(:,:) * (wta(:,:,k)+wtb(:,:,k))
         zfull(:,:,k) = zb(:,:) + rt(:,:) * wtb(:,:,k)
         ! eta/step-mountain option
         if (present(mask)) then
            where (mask(:,:,k) < 0.5)
                zt(:,:) = Vgrid % fhalf(k) ! use reference height profile
                zfull(:,:,k) = 0.5*(Vgrid % fhalf(k)+Vgrid % fhalf(k+1))
            endwhere
         endif
         zb(:,:) = zt(:,:)
         if (present(zhalf)) zhalf(:,:,k) = zb(:,:)

      enddo

  call mpp_clock_end (id(6))

end subroutine compute_geop_height

!#######################################################################

subroutine compute_height (Vgrid, fssl, temp, sphum, pfull, phalf, &
                           zfull, zhalf, mask)

!-----------------------------------------------------------------------
! Compute the geopotential height (in meters) at full and half model levels.
!
!   Vgrid   = vertical grid constants
!   fssl    = geopotential height at eta=1
!   temp    = layer mean temperature (at full model levels)
!   sphum   = layer mean specific humidity (at full model levels)
!   pfull   = pressure at model levels (full levels)
!   phalf   = pressure at model layer interfaces (half levels)
!   zfull   = geop height (in meters) at model levels (full levels)
!   zhalf   = geop height (in meters) at model layer interfaces (half levels)
!   mask    = grid box mask for eta coordinate topography
!
! Assumes that specific humidity is used to compute virtual temperature.
!-----------------------------------------------------------------------

   type(vert_grid_type), intent(in)    :: Vgrid
   real, intent(in),  dimension(:,:)   :: fssl
   real, intent(in),  dimension(:,:,:) :: temp, sphum, pfull, phalf
   real, intent(out), dimension(:,:,:) :: zfull, zhalf
   real, intent(in),  optional         :: mask(:,:,:)

   real, dimension(size(temp,1),size(temp,2),size(temp,3)) ::  &
                                                   wta, wtb, vtemp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  call compute_pres_weights ( Vgrid, phalf, pfull, wta, wtb )

  vtemp = temp * (1.0+d608*sphum)  ! WARNING: also computed in bgrid core
                                   !   potential exists for future problems!!!

  if (present(mask)) then
     call compute_geop_height ( Vgrid, fssl, vtemp, wta, wtb, &
                                zfull, zhalf, mask )
  else
     call compute_geop_height ( Vgrid, fssl, vtemp, wta, wtb, &
                                zfull, zhalf )
  endif

  zfull = zfull * ginv
  zhalf = zhalf * ginv

end subroutine compute_height

!#######################################################################

subroutine compute_height_bottom ( Vgrid, pssl, tbot, qbot,  &
                                   zbot, pbot, kbot )

!-----------------------------------------------------------------------
! Compute the height above surface (in meters) of the lowest model level.
!
!   Vgrid   = vertical grid constants
!   pssl    = pressure at eta=1
!   tbot    = mean temperature of lowest model layer
!   qbot    = mean specific humidity of lowest model layer
!
!   zbot    = height above surface (in meters) of the lowest model level
!   pbot    = pressure at lowest model level
!
!   mask    = grid box mask for eta coordinate topography
!

   type(vert_grid_type), intent(in)  :: Vgrid
   real, intent(in),  dimension(:,:) :: pssl, tbot, qbot
   real, intent(out), dimension(:,:) :: zbot, pbot
   integer, intent(in),  optional    :: kbot(:,:)

   real, dimension(size(pssl,1),size(pssl,2)) :: rt, dp, phb, pht, &
                                                 lphb, lpht, lpf
   integer :: i, j, kb

!  ----- pressure at top and bottom interface of bottom level -----

   call mpp_clock_begin (id(7))

 ! compute half level pressure at the "top" and "bottom" of lowest model layer
   if (present(kbot)) then
      do j = 1, size(pssl,2)
      do i = 1, size(pssl,1)
         kb = kbot(i,j)
         pht(i,j) = Vgrid%peta(kb  ) + Vgrid%eta(kb  )*pssl(i,j)
         phb(i,j) = Vgrid%peta(kb+1) + Vgrid%eta(kb+1)*pssl(i,j)
      enddo
      enddo
   else
         kb = Vgrid%nlev
         pht(:,:) = Vgrid%peta(kb+1) + Vgrid%eta(kb+1)*pssl(:,:)
         phb(:,:) = Vgrid%peta(kb  ) + Vgrid%eta(kb  )*pssl(:,:)
   endif

 ! compute log(pressure at lowest level)
   dp = phb - pht
   lphb = log(phb)
   lpht = log(pht)
   lpf = (phb*lphb-pht*lpht)/dp -1
   
 ! compute virtual temperature using specific humidity
   rt = ginv*RDGAS * (tbot * (1.+d608*qbot))
   zbot = rt * (lphb-lpf)
   pbot = exp(lpf)

   call mpp_clock_end (id(7))

end subroutine compute_height_bottom

!#######################################################################

end module bgrid_vert_mod

