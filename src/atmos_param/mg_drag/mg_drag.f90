module mg_drag_mod

!=======================================================================
!         MOUNTAIN GRAVITY WAVE DRAG - PIerrehumbert (1986)            !
!=======================================================================

!-------------------------------------------------------------------
!  Calculates partial tendencies for the zonal and meridional winds
!  due to the effect of mountain gravity wave drag 
!-------------------------------------------------------------------

 use  topography_mod, only: get_topog_stdev

 use         fms_mod, only: mpp_npes, field_size, file_exist, write_version_number, stdlog, &
                            mpp_pe, mpp_root_pe, error_mesg, FATAL, NOTE, read_data, write_data,  &
                            open_namelist_file, close_file, check_nml_error, open_restart_file, mpp_error
 use      fms_io_mod, only: get_restart_io_mode
 use   constants_mod, only: Grav, Kappa, RDgas, cp_air

!-----------------------------------------------------------------------
 implicit none
!-----------------------------------------------------------------------

 private

 character(len=128) :: version = '$Id: mg_drag.f90,v 11.0 2004/09/28 19:19:38 fms Exp $'
 character(len=128) :: tagname = '$Name: lima $'

 real, parameter :: p00 = 1.e5

!---------------------------------------------------------------------
!     Ghprime - array of sub-grid scale mountain height variance
!-----------------------------------------------------------------------

  real, allocatable, dimension(:,:) :: Ghprime
!-----------------------------------------------------------------------
!Contants
!     grav    value of gravity
!     rdgas    universal gas constant for dry air
!     kappa  2/7 (i.e., R/Cp)
!-----------------------------------------------------------------------

 logical :: module_is_initialized = .false.

!---------------------------------------------------------------------
! --- NAMELIST (mg_drag_nml)
!---------------------------------------------------------------------
!     xl_mtn      effective mountain length ( set currently to 100km)
!     acoef       order unity "tunable" parameter
!     gmax    order unity "tunable" parameter 
!             (may be enhanced to increase drag)
!     rho     stand value for density of the air at sea-level (1.13 KG/M**3)
!     low_lev_frac - fraction of atmosphere (from bottom up) considered
!              to be "low-level-layer for base flux calc. and where no
!              wave breaking is allowed.
!     flux_cut_level pressure level (Pa) above which flux divergence is set to zero 
!-----------------------------------------------------------------------

 real :: &
      xl_mtn=1.0e5 &
!     & ,gmax=1.0, acoef=1.0
!  v197 value of gmax = 2.0
      ,gmax=2.0, acoef=1.0, rho=1.13  &
!  v197 value for low-level-layer
      ,low_lev_frac = .23

real  ::  flux_cut_level= 0.0

logical :: do_netcdf_restart = .true.
logical :: do_conserve_energy = .false.
logical :: do_mcm_mg_drag = .false.
character(len=128) :: source_of_sgsmtn = 'input'

    namelist / mg_drag_nml / do_netcdf_restart,  &
                             xl_mtn, gmax, acoef, rho, low_lev_frac, &
                             do_conserve_energy, do_mcm_mg_drag,     &
                             source_of_sgsmtn, flux_cut_level 

 public mg_drag, mg_drag_init, mg_drag_end

 contains

!#############################################################################      

 subroutine mg_drag (is, js, delt, uwnd, vwnd, temp, pfull, phalf, &
                    zfull,zhalf,dtaux,dtauy,dtemp,taubx, tauby, tausf,&
                    kbot)
!===================================================================

! Arguments (intent in)

 integer, intent(in) :: is,js
 real, intent(in)    :: delt
 real, intent(in), dimension (:,:,:) :: &
     &             uwnd, vwnd, temp, pfull, phalf, zfull, zhalf
 integer, intent(in), optional, dimension(:,:)   :: kbot

!
!      INPUT
!      -----
!
!      is,js   - integers containing the starting
!                  i,j indices from the full horizontal grid
!      delt     time step in seconds
!      UWND     Zonal wind (dimensioned IDIM x JDIM x KDIM)
!      VWND     Meridional wind (dimensioned IDIM x JDIM x KDIM)
!      TEMP     Temperature at full model levels
!                   (dimensioned IDIM x JDIM x KDIM)
!      PFULL    Pressure at full model levels
!                   (dimensioned IDIM x JDIM x KDIM)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x KDIM+1)
!      ZHALF    Height at half model levels
!                   (dimensioned IDIM x JDIM x KDIM+1)
!      ZFULL    Height at full model levels
!                   (dimensioned IDIM x JDIM x KDIM+1)
!      KBOT     optional;lowest model level index (integer)
!                   (dimensioned IDIM x JDIM)
!===================================================================
! Arguments (intent out)

 real, intent(out), dimension (:,:) :: taubx, tauby
 real, intent(out), dimension (:,:,:) :: dtaux, dtauy, dtemp, tausf

!      OUTPUT
!      ------

!       TAUBX, TAUBY  base momentum flux componenets - output for diagnostics
!                   (dimensioned IDIM x JDIM)-kg/m/s**2
!                   = -(RHO*U**3/(N*XL))*G(FR) FOR N**2 > 0
!                   =          0               FOR N**2 <=0
!      DTAUX    Tendency of the zonal wind component deceleration 
!                   (dimensioned IDIM x JDIM x KDIM)
!      DTAUY    Tendency of the meridional wind component deceleration 
!                   (dimensioned IDIM x JDIM x KDIM)
!      dtemp    Tendency of temperature due to dissipation of ke
!
!      TAUSF = "CLIPPED" SAT MOMENTUM FLUX ( AT HALF LEVELS below top)
!                  
!===================================================================

!-----------------------------------------------------------------------

!     LLA IS DEFINED AS THE NUMBER OF LEVELS UP FROM THE LOWEST USED
!     TO CALCULATE THE "LOW-LEVEL" AVERAGES.

!     THIS ROUTINE COMPUTES THE DECELERATION OF THE ZONAL WIND AND
!     MERIDIONAL WIND DUE TO MOUNTAIN GRAVITY WAVE DRAG.  THE
!     PARAMETERIZATION WAS DEVELOPED BY R. PIERREHUMBERT AND ADAPTED
!     TO THE SPECTRAL MODEL BY B. STERN.  THE SCHEME IS STRUCTURED TO
!     INCLUDE 4 MAIN (GENERALLY VALID) COMPONENTS
!              1)  CALCULATION OF A BASE MOMENTUM FLUX(TAUB) WHICH IS
!                  A FUNCTION OF LOW-LEVEL->  WINDS, BRUNT-VAISALA FREQ,
!                  AND DENSITY  AS WELL AS THE SUB-GRID SCALE MOUNTAIN
!                  HEIGHT AND EFFECTIVE MOUNTAIN LENGTH.
!              2)  CALCULATION OF A SATURATION MOMENTUM FLUX PROFILE
!                  (TAUS(P) ) - IN GENERAL THIS IS A FUNCTION OF THE
!                  VERTICAL PROFILES OF WINDS, BRUNT VAISALA FREQ AND
!                  DENSITY.
!              3)  DETERMINE THE ACTUAL MOMENTUM FLUX PROFILE.  IT IS
!                  EQUAL TO THE FLUX ENTERING THE LAYER FROM BELOW
!                  BUT CANNOT EXCEED THE SATURATION FLUX IN THAT LAYER.
!              4)  CALCULATE THE DE-CELERATION DUE TO THE DRAG.
!     SATURATION MOMENTUM FLUX PROFILES
!          SCHEME 1:  LINEAR DROP OFF (IN P OR SIGMA) FROM THE BASE
!                     FLUX AT THE BOTTOM OF THE MODEL TO ZERO AT SIGTOP
!                           (^4/86)
!          SCHEME 2:  FUNCTION OF DENSITY, WINDS AND BRUNT VAISALA FREQ.
!                     V SCALE OF WAVE(D)  -
!                       A. FROM WKB THEORY
!                           (^6/87)
!                       B. FROM EXTENSION TO WKB THEORY
!                           (^7/87)
!     THE DECELERATION  IS PROPORTIONAL TO DTAUP/DSIGMA -> MOMENTUM FLUX
!     ABSORPTION WILL TAKE PLACE ONLY IN THOSE REGIONS WHERE TAUP VARIES
!     IN THE VERTICAL - I.E. WAVE BREAKING LAYERS.

!=======================================================================
!  (Intent local)
 real , dimension(size(uwnd,1),size(uwnd,2)) ::  xn, yn, psurf,ptop,taub
 real , dimension(size(uwnd,1),size(uwnd,2),size(uwnd,3)) ::  theta 
 real , dimension(size(uwnd,1),size(uwnd,2),size(uwnd,3)+1) ::  taus
 real vsamp
integer, dimension (size(uwnd,1),size(uwnd,2)) :: ktop, kbtm
integer id, jd, idim, jdim, kdim, kdimm1, kdimp1, ie, je

!              XN,YN  = PROJECTIONS OF "LOW LEVEL" WIND
!                       IN ZONAL & MERIDIONAL DIRECTIONS
!              TAUB = BASE MOMENTUM FLUX
!                   = -(RHO*U**3/(N*XL))*G(FR) FOR N**2 > 0
!                   =          0               FOR N**2 <=0
!              TAUS = SATURATION MOMENTUM FLUX ( AT HALF LEVELS)
!                   = (XN,XY)*(1-(AETA-1)/(SIGTOP-1))*TAUB  - SCHEME 1
!                   = -DENSITY(L)*UMAG(L)*D(L)*GMAX/XL       - SCHEME 2
!                   -> -AETA(L)*PS*UMAG(L)*D(L)*GMAX/XL
!      THETA    POTENTIAL temperature at full model levels
!                   (dimensioned IDIM x JDIM x KDIM)
!      PSURF    Surface pressure 
!                   (dimensioned IDIM x JDIM)
!      PTOP     Pressure at top of low-level layer
!                   (dimensioned IDIM x JDIM)
!      KTOP     Top model level index included in low-level layer
!                   (dimensioned IDIM x JDIM)
!      KBTM     Bottom model level index included in low-level layer
!                   usually the lowest level 
!                   (dimensioned IDIM x JDIM)
!-----------------------------------------------------------------------
!  type loop indicies
 integer i, j, k, kd, kb, kt, kbp1, ktm1 
!-----------------------------------------------------------------------
!  Local variables needed only for code that
!  implements supersource-like gravity wave drag.

integer :: klast, kcrit
real    :: sigtop, small=1.e-10

real,    dimension(size(uwnd,1),size(uwnd,2))              :: ulow, vlow, tlow, thlow
real,    dimension(size(uwnd,1),size(uwnd,2))              :: rlow, zsvar, bvfreq, x
real,    dimension(size(uwnd,1),size(uwnd,2))              :: depth, ave_p
integer, dimension(size(uwnd,1),size(uwnd,2))              :: ntop 
real,    dimension(size(uwnd,1),size(uwnd,2),size(uwnd,3)) :: th, sh_ang, test
real,    dimension(size(uwnd,1),size(uwnd,2),size(uwnd,3)) :: sigma, del_sigma
!real,    dimension(size(uwnd,1),size(uwnd,2),size(uwnd,3)+1) :: sigma_half

!---------------------------------------------------------------------

  idim = size( uwnd, 1 )
  jdim = size( uwnd, 2 )
  kdim = size( uwnd, 3 )
  kdimm1 = kdim - 1
  kdimp1 = kdim + 1

!-----------------------------------------------------------------------

!        CODE VARIABLES     DESCRIPTION

!              XN,YN  = PROJECTIONS OF "LOW LEVEL" WIND
!                       IN ZONAL & MERIDIONAL DIRECTIONS
!              TAUB = BASE MOMENTUM FLUX
!                   = -(RHO*U**3/(N*XL))*G(FR) FOR N**2 > 0
!                   =          0               FOR N**2 <=0
!              TAUS = SATURATION MOMENTUM FLUX ( AT HALF LEVELS)
!                   = (XN,XY)*(1-(AETA-1)/(SIGTOP-1))*TAUB  - SCHEME 1
!                   = -DENSITY(L)*UMAG(L)*D(L)*GMAX/XL       - SCHEME 2
!                   -> -AETA(L)*PS*UMAG(L)*D(L)*GMAX/XL
!              TAUP = MOMENTUM FLUX ( AT HALF LEVELS)
!              TAUP(L) = MIN ( TAUP(L-1),TAUS(L))
!              ULOW = "LOW-LEVEL" WIND MAGNITUDE (M/S)   (= U )
!                    AVERAGE UP TO ^2KM ABOVE SURFACE(LOWEST 1/3 SIGMAS)
!              UMAG = V. PROFILE OF WIND MAGNITUDES-AT HALF LEVS (=U(L))
!              DUDZ = VERTICAL DERIVATIVE OF U(L) WITH RESPECT TO Z
!                     DEFINED AT FULL LEVELS
!              DU2DZ2 = 2ND DERIVATIVE OF U(L) WITH RESPECT TO Z
!                     DEFINED AT HALF LEVELS
!              D = CHARACTERISTI! V. LENGTH SCALE OF WAVES (=D(L))
!                  FOR WKB D(L) = U(L)/N(L)
!                  FOR EXTENDED WKB 1/D**2 = N(L)**2/U(L)**2
!                                            - D2UDZ2(L)/U(L)
!              BNV,BNVK = "LOW-LEVEL",V. PROFILE -  BRUNT VAISALA FREQ(1
!                                                                 (= N,N
!              BNV2,BNVK2 = N**2, N(L)**2
!              HPRIME = Sub-grid scale mountain height 
!                       over local domain (IDIM x JDIM)
!              XL = EFFECTIVE MOUNTAIN LENGTH = (100KM EVERYWHERE)
!              SIGTOP = HIGHEST LEVEL TO WHICH GRAVITY WAVE
!                         MOMENTUM FLUX WILL BE DISTRIBUTED.
!              G = GMAX*FR**2/(FR**2+A**2)
!              	  GMAX = 1.0
!              	  A = 1.0
!=======================================================================

if ( .not.do_mcm_mg_drag ) then

!--- export sub grid scale topography
  ie = is + idim - 1
  je = js + jdim - 1

!-----------------------------------------------------------------------
!     vsamp is a vertical sampling coefficient which serves to amplify
!     the windshear wkb extension term in the calculation of d.
!     it increases this term to adjust for the deficiency of coarse
!     vertical resolution properly resolving the vertical windshear.

!     vsamp = (kdim+63)/kdim
      vsamp = 1.0
!-----------------------------------------------------------------------
!  calculate bottom of low-level layer = lowest level unless kbot is present
    if (present(kbot)) then
       kbtm(:,:) = kbot(:,:)
    else
       kbtm(:,:) = kdim
    endif
!  calculate top of low-level layer, first get surface p from phalf
    if (present(kbot)) then
       do j=1,jdim
       do i=1,idim
         psurf(i,j) = phalf(i,j,kbtm(i,j)+1)
       end do
       end do
    else
       psurf(:,:) = phalf(:,:,kdimp1)
    endif
!     print *,'psurf=', psurf
!  Based on fraction of model atmosphere to be considered "low-level"
!  (input via namelist), find highest model level.

    ptop(:,:) = (1.-low_lev_frac)*psurf(:,:)
    do kd=kdim,1,-1 
         where (pfull(:,:,kd) .ge. ptop(:,:)) 
           ktop(:,:) = kd
         end where
    end do
!  Make sure that low-level layer is at least 2 layer thick
    ktop(:,:) = min(ktop(:,:),(kbtm(:,:)-1) )
!     print *,'ptop=', ptop
!     print *,'ktop=', ktop

!  calculate base flux
    call mgwd_base_flux (is,js,uwnd,vwnd,temp,pfull,phalf,ktop,kbtm,theta, &
         &               xn,yn,taub)

!  split taub in to x and y components
    taubx(:,:) = taub(:,:)*xn(:,:)
    tauby(:,:) = taub(:,:)*yn(:,:)

!  calculate saturation flux profile
    call mgwd_satur_flux (uwnd,vwnd,temp,theta,ktop,kbtm, &
         &                xn,yn,taub,pfull, phalf,zfull,zhalf,vsamp,taus)

!  calculate mountain gravity wave drag tendency contributions
    call mgwd_tend (is,js,xn,yn,taub,phalf,taus,dtaux,dtauy, tausf)

else if ( do_mcm_mg_drag ) then

    if(present(kbot)) then
      call error_mesg ('mg_drag','kbot cannot be present in the calling arguments when using the Manabe Climate Model option',FATAL)
    endif

    do k=1,kdim
      sigma(:,:,k) = pfull(:,:,k)/phalf(:,:,kdimp1)
      del_sigma(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/phalf(:,:,kdimp1)
!     sigma_half(:,:,k) = phalf(:,:,k)/phalf(:,:,kdimp1)
    enddo
!    sigma_half(:,:,kdimp1) = 1.0

    ie = is + idim - 1
    je = js + jdim - 1

    zsvar = (Ghprime(is:ie,js:je))**2

    do k = 1,kdim
      th(:,:,k) = temp(:,:,k) + grav*(zfull(:,:,k)-zhalf(:,:,kdimp1))*kappa/rdgas
    end do

    sigtop = 1.0 - low_lev_frac
    do j = 1,jdim
      do i = 1,idim
        do k = kdim,1,-1
          if (sigma(i,j,k) .lt. sigtop) then
              if ( (sigtop - sigma(i,j,k)) .le. (sigma(i,j,k+1)-sigtop) ) then
                ntop(i,j) = k
              else
                ntop(i,j) = k + 1
              endif
              go to 10
          endif
        end do
        10 continue
      enddo
    enddo

    ulow  = 0.0
    vlow  = 0.0
    tlow  = 0.0
    thlow = 0.0
    depth = 0.0

    do j = 1,jdim
      do i = 1,idim
        do k = ntop(i,j), kdim
           ulow(i,j)  = ulow(i,j)  + del_sigma(i,j,k)* uwnd(i,j,k) 
           vlow(i,j)  = vlow(i,j)  + del_sigma(i,j,k)* vwnd(i,j,k) 
           tlow(i,j)  = tlow(i,j)  + del_sigma(i,j,k)* temp(i,j,k)  
           thlow(i,j) = thlow(i,j) + del_sigma(i,j,k)*th(i,j,k)
           depth(i,j) = depth(i,j) + del_sigma(i,j,k)
        end do
      enddo
    enddo
    ulow  = ulow/depth
    vlow  = vlow/depth
    tlow  = tlow/depth
    thlow = thlow/depth

    do j = 1,jdim
      do i = 1,idim
        ave_p(i,j)  = (phalf(i,j,ntop(i,j)-1) + phalf(i,j,kdim))/2.
        bvfreq(i,j) = (th(i,j,kdim) - th(i,j,ntop(i,j)))/(pfull(i,j,kdim) - pfull(i,j,ntop(i,j)))
        bvfreq(i,j) = -grav*grav*ave_p(i,j)*bvfreq(i,j)/(rdgas*tlow(i,j)*thlow(i,j))  ! thlow should be tlow !IH
      enddo
    enddo

    where(bvfreq > 0.0)
      bvfreq = sqrt(bvfreq)
    elsewhere
      bvfreq = 0.0
    endwhere

!     TK mod: original had rk0*grav*bvfreq*zsvar*ulow .... etc..

    x = grav*bvfreq*zsvar/(rdgas*tlow*xl_mtn)

    rlow = 1.0/sqrt(ulow**2 + vlow**2 + small)
    do k = 1, kdim
      sh_ang(:,:,k) = rlow*(ulow*uwnd(:,:,k) + vlow*vwnd(:,:,k))/      &
                     sqrt(uwnd(:,:,k)**2 + vwnd(:,:,k)**2 + small)
    enddo
    sh_ang = min(sh_ang, 0.99999)
    sh_ang = max(sh_ang,-0.99999)
    sh_ang = acos(sh_ang)

    test = 1.0
    where (sh_ang > 2.*atan(1.0))
      test = 0.0
    endwhere
    do j = 1,jdim
      do i = 1,idim
        do k=ntop(i,j),kdim
          test(i,j,k) = 1.0
        enddo
      enddo
    enddo

    do j = 1,jdim
      do i = 1,idim
        do k = ntop(i,j) - 1, 1, -1
          klast = k
          if (test(i,j,k) /= 1.0 )  go to 20
        end do
        klast = 0
        20 continue
        if ( klast /= 0 )  test(i,j,1:klast) = 0.0
        kcrit = klast + 1
        x(i,j) = x(i,j)/(1.-sigma(i,j,kcrit))  
!       should be  x(i,j) = x(i,j)/(1 - sigma_half(klast)/sigma_half(kdim))
      end do
    end do

    do k = 1,kdim
      dtaux(:,:,k) = - x*ulow*test(:,:,k)
      dtauy(:,:,k) = - x*vlow*test(:,:,k)
    end do

    taub = 0.0
    taubx = 0.0
    tauby = 0.0
    tausf = 0.0

endif

!  calculate temperature tendency due to dissipation of kinetic energy
if (do_conserve_energy) then
  dtemp = -((uwnd+.5*delt*dtaux)*dtaux + (vwnd+.5*delt*dtauy)*dtauy)/cp_air
else
  dtemp = 0.0
endif

return
end subroutine mg_drag
!=======================================================================

!#############################################################################      
 
subroutine mgwd_base_flux (is,js,uwnd,vwnd,temp,pfull,phalf,ktop,kbtm,  &
                          theta,xn,yn,taub)
                                  


!-------------------------------------------------------------------
!  calculates base momentum flux  - taub
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)
 real, intent(in), dimension (:,:,:) :: uwnd, vwnd, temp, pfull, phalf
 integer, intent(in), dimension (:,:) :: ktop, kbtm
 integer, intent(in)   :: is, js
!===================================================================
! Arguments (intent out)
 real, intent(out), dimension (:,:) :: xn, yn, taub
 real , intent(out), dimension (:,:,:) :: theta
!===================================================================
! Arguments (intent inout)
!=======================================================================
!  (Intent local)
real , dimension(size(uwnd,1),size(uwnd,2)) :: sumw, delp, ulow, bnv, &
     &  hprime, fr, g, ubar, vbar, bnv2 
real grav2, xli, a, small
 integer idim, jdim,kdim,ie, je
!-----------------------------------------------------------------------
!  type loop indicies
 integer i, j, k, kb, kt, kbp1, ktm1 
!-----------------------------------------------------------------------
!===================================================================

!-------------------------------------------------------------------
! --- DEFINE CURRENT WINDOW & GET GLOBAL VARIABLES
!-------------------------------------------------------------------

  idim = size( uwnd, 1 )
  jdim = size( uwnd, 2 )
  kdim = size( uwnd, 3 )
  ie = is + idim - 1
  je = js + jdim - 1
  hprime(:,:) = Ghprime(is:ie,js:je)

! define local scalar variables
  xli=1.0/xl_mtn
  grav2=grav*grav
  a = acoef 


!-----------------------------------------------------------------------
!     <><><><><><><><>   base flux code   <><><><><><><><>
!-----------------------------------------------------------------------

!  initialize arrays
        sumw(:,:) = 0.0
        ubar(:,:) = 0.0
        vbar(:,:) = 0.0
        ulow(:,:) = 0.0
        taub(:,:) = 0.0
        xn  (:,:) = 0.0
        yn  (:,:) = 0.0


!     compute low-level averages
!     --------------------------

      do j=1,jdim
        do i=1,idim
          do k=ktop(i,j),kbtm(i,j)
            delp(i,j) = phalf(i,j,k+1)-phalf(i,j,k)
            sumw(i,j) = sumw(i,j) + delp(i,j)
            ubar(i,j) = ubar(i,j) + uwnd(i,j,k)*delp(i,j)
            vbar(i,j) = vbar(i,j) + vwnd(i,j,k)*delp(i,j)
          end do
        end do
      end do
!    print *, 'low-lev aves computed, ubar, vbar =', ubar, vbar

!     calculate projections of low level flow onto wind components (u&v)
!     ------------------------------------------------------------------
        sumw(:,:) = 1./sumw(:,:)
        ubar(:,:) = ubar(:,:) * sumw(:,:)
        vbar(:,:) = vbar(:,:) * sumw(:,:)
        ulow(:,:) =sqrt(ubar(:,:)*ubar(:,:) + vbar(:,:)*vbar(:,:))
        xn(:,:) = ubar(:,:)/(ulow(:,:) + 1.0e-20)
        yn(:,:) = vbar(:,:)/(ulow(:,:) + 1.0e-20)


!     calculate squared brunt vaisala freq
!     ------------------------------------

      theta(:,:,:)=temp(:,:,:)*(pfull(:,:,:)/p00)**(-kappa)
!  v197 uses p* as reference vlues for theta, in above 1000 hPa is used
!      theta(:,:,:)=temp(:,:,:)*(pfull(:,:,:)/ &
!     &             phalf(:,:,kdim+1))**(-kappa)
 
      do j=1,jdim
        do i=1,idim
          kt=ktop(i,j)
          kb=kbtm(i,j)
          bnv2(i,j) = grav2*(pfull(i,j,kt)+pfull(i,j,kb)) &
                 * (theta(i,j,kt)-theta(i,j,kb)) &
              / ( rdgas*(theta(i,j,kt)+theta(i,j,kb)) &
                 * (pfull(i,j,kb)-pfull(i,j,kt)) &
                 *.5*(temp(i,j,kt)+temp(i,j,kb)))
        end do
      end do

!      calculate bnv,fr,g,taub,xn,yn - if n**2>0
!      -----------------------------------------
           small = epsilon(ulow)

           where (bnv2(:,:) .gt. 0.0) 
             bnv(:,:) = sqrt(bnv2(:,:))
             fr (:,:) = bnv(:,:)*hprime(:,:)/(ulow(:,:) + small)
             g  (:,:) = gmax*fr(:,:)*fr(:,:)/(fr(:,:)*fr(:,:)+a*a)
             taub(:,:) = -rho*xli*ulow(:,:)*ulow(:,:)*ulow(:,:) &
     &                 / bnv(:,:)*g(:,:)
           elsewhere
             bnv(:,:) = 0.0
             fr (:,:) = 0.0
             g  (:,:) = 0.0
           endwhere

end subroutine mgwd_base_flux

!#############################################################################      

subroutine mgwd_satur_flux (uwnd,vwnd,temp,theta,ktop,kbtm, &
                           xn,yn,taub,pfull,phalf,zfull,zhalf,vsamp,taus)

!===================================================================
! Arguments (intent in)
 real, intent(in), dimension (:,:,:)  :: &
     &             uwnd, vwnd, temp, theta, pfull, phalf,zfull, zhalf
 real, intent(in), dimension (:,:)    :: xn, yn, taub
 real, intent(in)                     :: vsamp 
 integer, intent(in), dimension (:,:) :: ktop, kbtm
!===================================================================
! Arguments (intent out)
 real, intent(out), dimension (:,:,:) :: taus
!=======================================================================
!  (Intent local)
 real , dimension(size(uwnd,1),size(uwnd,2),size(uwnd,3)) ::  &
     &       dterm, dudz  
 real , dimension(size(uwnd,1),size(uwnd,2),size(uwnd,3)+1) ::  &
     &       umag, bnvk2, d,d2, d2i, d2udz2, extend
 real grav2, xli, small
 integer :: idim, jdim, kdim, kdimm1, kdimp1
!-----------------------------------------------------------------------
!  type loop indicies
 integer i, j, k, kb, kt, kbp1, ktm1 
!-----------------------------------------------------------------------
!  type flux cutoff 
 integer kcut
!=======================================================================


  idim = size( uwnd, 1 )
  jdim = size( uwnd, 2 )
  kdim = size( uwnd, 3 )
  kdimm1 = kdim - 1
  kdimp1 = kdim + 1

! define local scalar variables
  xli=1.0/xl_mtn
  grav2=grav*grav

!-----------------------------------------------------------------------
!     <><><><><><><><>   saturation flux code   <><><><><><><><>
!-----------------------------------------------------------------------

!     scheme 1 - linear profile
!     do 35 l=1,lp1
!     do 35 i=1,idim
!     taus(i,l) = taub(i) * (1-(eta(l)-1)/(sigtop-1) )
!35    continue

!-----------------------------------------------------------------------

!     ********** scheme 2 - wave breaking formulation  **********

!-----------------------------------------------------------------------


!     calculate wind magnitude at 1/2 levels
!     --------------------------------------------

      do k=2,kdim
        umag(:,:,k) =  (0.50*(uwnd(:,:,k-1)+uwnd(:,:,k))*xn(:,:) &
                     + 0.50*(vwnd(:,:,k-1)+vwnd(:,:,k))*yn(:,:))
        umag(:,:,k) = abs( umag(:,:,k) )
      end do


!     set wind magnitude at top of model = to magnitude at top full
!     level.

        umag(:,:,1) = uwnd(:,:,1)*xn(:,:) + vwnd(:,:,1)*yn(:,:)
        umag(:,:,1) = abs( umag(:,:,1) )

!     set wind magnitude at ground = 0.

      do j=1,jdim
        do i=1,idim
          kbp1=kbtm(i,j)+1
          do k=kbp1,kdimp1
            umag(i,j,k) = 0.0
          end do
        end do
      end do

!     set minimum wind magnitude

      small = epsilon (umag)
      where ( umag .lt. small ) umag = 0.0


!      print *, ' umag for sat flux =', umag

!-----------------------------------------------------------------------

!     calculate vertical derivatives of umag, to be used in
!     the extension to the wkb approach for determining d.
!     -- derivative of umag with respect to z is computed at
!        full levels and stored in dudz.
!     dudz(1) is defined using an uncentered difference

         dudz(:,:,1) = (umag(:,:,1)-umag(:,:,2)) &
     &                /(zfull(:,:,1)-zhalf(:,:,2))

      do k=2,kdim
         dudz(:,:,k) = (umag(:,:,k)-umag(:,:,k+1)) &
     &                /(zhalf(:,:,k)-zhalf(:,:,k+1))
      end do

!      print *, ' dudz for sat flux =', dudz



!     assume vertical derivative of umag at the boundaries=0 and
!     compute 2nd derivatives there using uncentered differencing

      do k=2,kdim
         d2udz2(:,:,k) = (dudz(:,:,k)-dudz(:,:,k-1)) &
     &                  /(zfull(:,:,k)-zfull(:,:,k-1))
      end do

!     set d2udz2 = 0 at the top of the atmosphere (original code)
!     set d2udz2 at the top of atm to level 2 value (new code)

!del    d2udz2(:,:,1) = 0.0
        d2udz2(:,:,1) = d2udz2(:,:,2)

      do  j=1,jdim
      do  i=1,idim
        kb=kbtm(i,j)
        kbp1=kb+1
        d2udz2(i,j,kbp1) = dudz(i,j,kb)/(zfull(i,j,kb)-zhalf(i,j,kbp1))
      end do
      end do

!      print *, ' d2udz2 for sat flux =', d2udz2


!-----------------------------------------------------------------------

!     compute wkb extension term for umag > 0
!     ---------------------------------------

         where (umag(:,:,:).gt.0.0) 
            extend(:,:,:) = vsamp*d2udz2(:,:,:)/umag(:,:,:)
         elsewhere
            extend(:,:,:) = 0.0
         endwhere

!      print *, ' wkb exten for sat flux =', extend


!     calculate brunt vaisala frequency at 1/2 levels
!     -----------------------------------------------------

      do k=2,kdim
        bnvk2(:,:,k) =  grav2*(pfull(:,:,k-1)+pfull(:,:,k)) &
     &          * (theta(:,:,k-1)-theta(:,:,k)) &
     &      /     ( rdgas*(theta(:,:,k-1)+theta(:,:,k)) &
     &  * (pfull(:,:,k)-pfull(:,:,k-1))*.5*(temp(:,:,k-1)+temp(:,:,k)) )
      end do


!     keep static stability constant in top & bottom layers of model
!     for taus calculations.

        bnvk2(:,:,1) = bnvk2(:,:,2)


      do j=1,jdim
      do i=1,idim
        kb=kbtm(i,j)
        kbp1=kb+1
        bnvk2(i,j,kbp1)   = bnvk2(i,j,kb)
        bnvk2(i,j,kdimp1) = bnvk2(i,j,kdim)
      end do
      end do

!      print *, ' brunt vaisala for sat flux =', bnvk2


!-----------------------------------------------------------------------

!     calculate d2i (=1/d**2) for umag .gt. 0
!     initialize d2i to a large number, which will result in a very
!     small vertical wavelength (d) where umag = 0.

         where (umag(:,:,:).gt.0.0) 
            d2i(:,:,:) = (bnvk2(:,:,:)/(umag(:,:,:)* &
     &                   umag(:,:,:)) - extend(:,:,:) )
         elsewhere
            d2i(:,:,:) = 1.0e+30
         endwhere

!      print *, ' 1/d**2 for sat flux =', d2i


!     for 1/d**2 approaching 0 calculate d by dividing by a
!     very small but finite number

         where (d2i(:,:,:) .lt. 1.e-30) 
            d(:,:,:) = 1.e+30
         elsewhere
            d2(:,:,:) = 1./d2i(:,:,:)
            d (:,:,:) = sqrt(d2(:,:,:))
         endwhere


!     set d=0 for umag=0.
         where (umag(:,:,:).eq.0.0) 
            d(:,:,:) = 0.0
         endwhere


!-----------------------------------------------------------------------

!      print *, 'd for sat flux =', d

!     calculation of the saturation flux profile for scheme 2
!     -------------------------------------------------------

      do j=1,jdim
        do i=1,idim
          kb=kbtm(i,j)
          kt=ktop(i,j)
          ktm1=kt-1

          do k=2,ktm1
            taus(i,j,k) = -phalf(i,j,k)*umag(i,j,k)*umag(i,j,k) &
     &                  *d(i,j,k)*xli*gmax &
     &                / (0.50*(temp(i,j,k-1)+temp(i,j,k))*rdgas)
          end do

          do k = kt,kdimp1
            taus(i,j,k) = taub(i,j)
          end do

        end do
      end do


!     keep taus profile constant across top model layer (original code)
!     calculate taus profile in top model layer (new code)
 
        taus(:,:,1) = taus(:,:,2)
!del        taus(:,:,1) = -phalf(:,:,1)*umag(:,:,1)*umag(:,:,1) &
!del     &                  *d(:,:,1)*xli*gmax / (temp(:,:,1)*rdgas)


!     do not allow wave breaking for unstable layers
 
       do k = 1,kdimp1
         where ( bnvk2(:,:,k) .lt. 0.0) 
           taus(:,:,k) =  taub(:,:)
         endwhere
       end do


! -------------------------------------------
!        tausat(:,:,1) = 0.             ! use all forcing
!         Instead,  let remaining flux escape above flux_cut_level
      if( flux_cut_level > 0.0 ) then 
         kcut= 1 
         do while( phalf(1,1,kcut) < flux_cut_level )
            kcut= kcut+1
         enddo

        do k= 1, kcut-1
            taus(:,:,k)= taus(:,:,kcut)
        enddo
      endif

end subroutine mgwd_satur_flux

!#############################################################################      

subroutine mgwd_tend (is,js,xn,yn,taub,phalf,taus,dtaux,dtauy,tausf)

!===================================================================
! Arguments (intent in)
 real, intent(in), dimension (:,:,:) :: phalf, taus
 real, intent(in), dimension (:,:) :: xn, yn, taub
 integer, intent(in)   :: is, js
!===================================================================
! Arguments (intent out)
 real, intent(out), dimension (:,:,:) :: dtaux, dtauy, tausf
!=======================================================================
!  (Intent local)
 real , dimension(size(phalf,1),size(phalf,2),size(phalf,3)) ::  dterm 
 real , dimension(size(phalf,1),size(phalf,2),size(phalf,3)+1) ::  taup
 integer kdim, kdimp1
!-----------------------------------------------------------------------
!  type loop indicies
 integer k, kd
!-----------------------------------------------------------------------
!=======================================================================


  kdim = size( dtaux, 3 )
  kdimp1 = kdim + 1


!-----------------------------------------------------------------------
!     <><><><><><><><>   MOMENTUM FLUX CODE   <><><><><><><><>
!-----------------------------------------------------------------------

!     CALCULATE FLUX FROM GROUND UP
!     -----------------------------

        taup (:,:,kdimp1) = taub(:,:)

      do kd=2,kdimp1
        k = kdimp1-kd+1
        tausf(:,:,k)=taup(:,:,k+1)
        taup(:,:,k) = max (taus(:,:,k),taup(:,:,k+1))
      end do

!     ALLOW FLUX TO ESCAPE THE TOP - DO NOT RE-DISTRIBUTE

!     <><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!     <><><><><><><><>   DE-CELERATION CODE   <><><><><><><><>
 
!     CALCULATE DECELERATION TERMS - DTAUX,DTAUY
!     ------------------------------------------

       do k=1,kdim
          dterm(:,:,k) = grav*(taup (:,:,k+1)-taup (:,:,k)) &
     &                     /(phalf(:,:,k+1)-phalf(:,:,k))
 
        dtaux(:,:,k) = xn(:,:)*dterm(:,:,k)
        dtauy(:,:,k) = yn(:,:)*dterm(:,:,k)
       end do

!  print sample output
!            print*, ' mgdrag output for i,j=', is,js
!            print *,'taub = ', taub(is,js)     
!            print *,'taus = ', taus(is,js,:)     
!            print *,'taup = ', taup(is,js,:)     


!     ***********************************************************

end subroutine mgwd_tend

!#######################################################################

  subroutine mg_drag_init( lonb, latb, hprime )

!=======================================================================
! ***** INITIALIZE Mountain Gravity Wave Drag
!=======================================================================

!---------------------------------------------------------------------
! Arguments (Intent in)
!     lonb  = longitude in radians of the grid box edges
!     latb  = latitude  in radians of the grid box edges
!---------------------------------------------------------------------
 real, intent(in), dimension(:) :: lonb, latb
 
!---------------------------------------------------------------------
! Arguments (Intent out - optional)
!     hprime  = array of sub-grid scale mountain heights
!---------------------------------------------------------------------
 real, intent(out), dimension(:,:), optional :: hprime
 
!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------
 integer  ::  ix, iy, unit, io, ierr
 logical  ::  answer
 integer, dimension(4) :: siz
 integer :: global_num_lon, global_num_lat

!=====================================================================

if(module_is_initialized) return

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------
  if( file_exist( 'input.nml' ) ) then
! -------------------------------------
   unit = open_namelist_file()
   ierr = 1
   do while( ierr .ne. 0 )
   read ( unit,  nml = mg_drag_nml, iostat = io, end = 10 ) 
   ierr = check_nml_error(io,'mg_drag_nml')
   end do
10 continue
   call close_file ( unit )
   call get_restart_io_mode(do_netcdf_restart)

! -------------------------------------
  end if

!---------------------------------------------------------------------
! --- Output version
!---------------------------------------------------------------------

  call write_version_number(version, tagname)
  if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=mg_drag_nml)

!---------------------------------------------------------------------
! --- Allocate storage for Ghprime
!---------------------------------------------------------------------

  ix = size(lonb(:)) - 1
  iy = size(latb(:)) - 1

  allocate( Ghprime(ix,iy) ) ; Ghprime = 0.0
  
!-------------------------------------------------------------------
  module_is_initialized = .true.
!---------------------------------------------------------------------
! --- Input hprime
!---------------------------------------------------------------------

  if ( trim(source_of_sgsmtn) == 'computed' ) then
    answer = get_topog_stdev ( lonb, latb, Ghprime )
    if ( .not.answer ) then
      call error_mesg('mg_drag_init','source_of_sgsmtn="'//trim(source_of_sgsmtn)//'"'// &
                      ', but topography data file does not exist', FATAL)
    endif
  else if ( trim(source_of_sgsmtn) == 'input' ) then
    if ( file_exist( 'INPUT/mg_drag.res.nc' ) ) then
       if (mpp_pe() == mpp_root_pe()) call mpp_error ('mg_drag_mod', &
            'Reading NetCDF formatted restart file: INPUT/mg_drag.res.nc', NOTE)
       call read_data ('INPUT/mg_drag.res.nc', 'ghprime', Ghprime)
    else if ( file_exist( 'INPUT/mg_drag.res' ) ) then
       if (mpp_pe() == mpp_root_pe()) call mpp_error ('mg_drag_mod', &
            'Reading native formatted restart file.', NOTE)
      unit = open_restart_file('INPUT/mg_drag.res','read')
      call read_data(unit, Ghprime)
      call close_file(unit)
    else
      call error_mesg ('mg_drag_init','source_of_sgsmtn="'//trim(source_of_sgsmtn)//'"'// &
                       ', but neither ./INPUT/mg_drag.res.nc  or  ./INPUT/mg_drag.res  exists', FATAL)
    endif
  else
    call error_mesg ('mg_drag_init','"'//trim(source_of_sgsmtn)//'"'// &
          ' is not a valid value for source_of_sgsmtn', FATAL)
  endif

! return sub-grid scale topography?
  if (present(hprime)) hprime = Ghprime
 
!=====================================================================
  end subroutine mg_drag_init

!#######################################################################

  subroutine mg_drag_end
  integer :: unit

  if(.not.module_is_initialized) return
  if(do_netcdf_restart) then
     if (mpp_pe() == mpp_root_pe()) call mpp_error ('mg_drag_mod', &
          'Writing NetCDF formatted restart file: RESTART/mg_drag.res.nc', NOTE)
     call write_data('RESTART/mg_drag.res.nc', 'ghprime', ghprime)
  else
     if (mpp_pe() == mpp_root_pe()) call mpp_error ('mg_drag_mod', &
          'Writing native formatted restart file.', NOTE)
     unit = open_restart_file('RESTART/mg_drag.res','write')
     call write_data(unit, Ghprime)
     call close_file(unit)
  endif
  deallocate(ghprime)
  module_is_initialized = .false.

  end subroutine mg_drag_end

!#######################################################################

end module mg_drag_mod
