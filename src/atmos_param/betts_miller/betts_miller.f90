
module betts_miller_mod

!----------------------------------------------------------------------
!use      utilities_mod, only:  file_exist, error_mesg, open_file,  &
!                               check_nml_error, get_my_pe, FATAL,  &
!                               close_file

use            fms_mod, only:  file_exist, error_mesg, open_namelist_file, &
                               check_nml_error, mpp_pe, mpp_root_pe, &
                               FATAL, close_file, write_version_number, stdlog

use sat_vapor_pres_mod, only:  escomp, descomp
use      constants_mod, only:  HLv,HLs,Cp_air,Grav,rdgas,rvgas, CP_VAPOR, & !s cph2ovapr -> CP_VAPOR
                               kappa, es0

implicit none
private
!---------------------------------------------------------------------
!  ---- public interfaces ----

   public  betts_miller, betts_miller_init, lcltabl

!-----------------------------------------------------------------------
!   ---- version number ----

 character(len=128) :: version = '$Id: betts_miller.f90,v 1.1.2.1 2005/05/13 18:26:51 pjp Exp $'
 character(len=128) :: tagname = '$Name:  $'

!-----------------------------------------------------------------------
!   ---- local/private data ----

    real :: d622 = 0.
    real :: d378 = 0.

    logical :: do_init=.true.

!-----------------------------------------------------------------------
!   --- namelist ----

real    :: tau_bm=7200.
real    :: rhbm = .8
logical :: do_simp = .true.
!logical :: do_enadjusttemp = .false.
logical :: do_shallower = .false.
logical :: do_changeqref = .false.
logical :: do_envsat = .false.
logical :: do_taucape = .false.
real    :: capetaubm = 900.
real    :: tau_min = 2400.
real    :: buoyancy_kick = 0.

namelist /betts_miller_nml/  tau_bm, rhbm, do_simp, &
!                            do_enadjusttemp, &
                            do_shallower, do_changeqref, &
                            do_envsat, do_taucape, capetaubm, tau_min, &
			    buoyancy_kick

!-----------------------------------------------------------------------
!           description of namelist variables
!
!  tau_bm    =  betts-miller relaxation timescale (seconds)
!
!  rhbm      = relative humidity that you're relaxing towards
!
!  do_simp = do the simple method where you adjust timescales to make 
!            precip continuous always
!
!  do_enadjusttemp = do the betts-miller (1986) way of conserving energy:
!                    adjust the reference temperature by a fixed amount w/ 
!                    height.  

! ***   if neither of the two above are true, then we instead adjust the 
!       humidity profile so its precipitation is increased to balance the 
!       precip_t change, and hence dries more than the 
! 
!  do_shallower = do the shallow convection scheme where it chooses a smaller
!                 depth such that precipitation is zero
! 
!  do_changeqref = do the shallow convection scheme where if changes the 
!                  profile of both q and T in order make precip zero
! 
!  do_envsat = reference profile is rhbm times saturated wrt environment 
!              (if false, it's rhbm times parcel)
!  
!  do_taucape = scheme where taubm is proportional to CAPE**-1/2
!
!  capetaubm = for the above scheme, the value of CAPE for which 
!              tau = tau_bm
!  
!  tau_min   = minimum relaxation time allowed for the above scheme
!
!  
!
!-----------------------------------------------------------------------

contains

!#######################################################################

   subroutine betts_miller (dt, tin, qin, pfull, phalf, coldT, &
                           rain, snow, tdel, qdel, q_ref, bmflag, &
                           klzbs, cape, cin, t_ref,invtau_bm_t,invtau_bm_q, &
                           capeflag, klcls, mask, conv)

!-----------------------------------------------------------------------
!
!                     Betts-Miller Convection Scheme
!
!-----------------------------------------------------------------------
!
!   input:  dt       time step in seconds
!           tin      temperature at full model levels
!           qin      specific humidity of water vapor at full
!                      model levels
!           pfull    pressure at full model levels
!           phalf    pressure at half (interface) model levels
!           coldT    should precipitation be snow at this point?
!   optional:
!           mask     optional mask (0 or 1.) 
!           conv     logical flag; if true then no betts-miller
!                       adjustment is performed at that grid-point or
!                       model level
!
!  output:  rain     liquid precipitation (kg/m2)
!           snow     frozen precipitation (kg/m2)
!           tdel     temperature tendency at full model levels
!           qdel     specific humidity tendency (of water vapor) at
!                      full model levels
!           bmflag   flag for which routines you're calling
!           klzbs    stored klzb values
!           cape     convectively available potential energy 
!           cin      convective inhibition (this and the above are before the 
!                    adjustment)
!           invtau_bm_t temperature relaxation timescale
!           invtau_bm_q humidity relaxation timescale
!           capeflag a flag that says why cape=0
!
!-----------------------------------------------------------------------
!--------------------- interface arguments -----------------------------

   real   , intent(in) , dimension(:,:,:) :: tin, qin, pfull, phalf
   real   , intent(in)                    :: dt
   logical   , intent(in) , dimension(:,:):: coldT
   real   , intent(out), dimension(:,:)   :: rain,snow, cape, &
       cin, invtau_bm_t, invtau_bm_q, capeflag
   integer, intent(out), dimension(:,:)   :: bmflag, klzbs, klcls
   real   , intent(out), dimension(:,:,:) :: tdel, qdel, q_ref, t_ref
   real   , intent(in) , dimension(:,:,:), optional :: mask
   logical, intent(in) , dimension(:,:,:), optional :: conv
!-----------------------------------------------------------------------
!---------------------- local data -------------------------------------

logical,dimension(size(tin,1),size(tin,2),size(tin,3)) :: do_adjust
logical :: avgbl
   real,dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &
             rin, esat, qsat, desat, dqsat, pmes, pmass
   real,dimension(size(tin,1),size(tin,2))             ::  &
                     hlcp, precip, precip_t
   real,dimension(size(tin,3))                         :: eref, rpc, tpc, &
                                                          tpc1, rpc1

   real                                                ::  & 
       cape1, cin1, tot, deltak, deltaq, qrefint, deltaqfrac, deltaqfrac2, &
       ptopfrac, es, capeflag1, plzb, plcl, small
   integer  i, j, k, ix, jx, kx, klzb, ktop, klcl
!-----------------------------------------------------------------------
!     computation of precipitation by betts-miller scheme
!-----------------------------------------------------------------------

      if (do_init) call error_mesg ('betts_miller',  &
                         'betts_miller_init has not been called.', FATAL)

      ix=size(tin,1)
      jx=size(tin,2)
      kx=size(tin,3)
      avgbl = .false.
      small = 1.e-10

! calculate r
       rin = qin/(1.0 - qin)
       do i=1,ix
          do j=1,jx
             cape1 = 0.
             cin1 = 0.
             tot = 0.
             klzb=0
! the bmflag is written out to show what aspects of the bm scheme is called
! bmflag = 0 is no cape, no convection
! bmflag = 1 is shallow conv, the predicted precip is less than zero
! bmflag = 2 is deep convection
             bmflag(i,j) = 0
             tpc = tin(i,j,:)
             rpc = rin(i,j,:)
! calculate cape, cin, level of zero buoyancy, and parcel properties 
! new code (second order in delta ln p and exact LCL calculation)
             call capecalcnew( kx,  pfull(i,j,:),  phalf(i,j,:),&
                            cp_air, rdgas, rvgas, hlv, kappa, tin(i,j,:), &
                            rin(i,j,:), avgbl, cape1, cin1, tpc, &
                            rpc, klzb, klcl)

! set values for storage
             capeflag(i,j) = capeflag1
             cape(i,j) = cape1
             cin(i,j) = cin1
             klzbs(i,j) = klzb
             klcls(i,j) = klcl
             if(cape1.gt.0.) then 
!             if((tot.gt.0.).and.(cape1.gt.0.)) then 
                bmflag(i,j) = 1
! reference temperature is just that of the parcel all the way up
                t_ref(i,j,:) = tpc
                do k=klzb,kx
! sets reference spec hum to a certain relative hum (change to vapor pressure, 
! multiply by rhbm, then back to spec humid)
                   if(do_envsat) then 
!                      call establ2(es,tin(i,j,k))
                      call escomp(tin(i,j,k),es)
                      es = es*rhbm
!                      rpc(k) = d622*es/(pfull(i,j,k) - d378*es)
                      rpc(k) = rdgas/rvgas*es/pfull(i,j,k)
                      q_ref(i,j,k) = rpc(k)/(1 + rpc(k))
                   else 
                      rpc(k) = rhbm*rpc(k)
!                      eref(k) = rhbm*pfull(i,j,k)*rpc(k)/(d622 + d378*rpc(k))
!                      rpc(k) = d622*eref(k)/(pfull(i,j,k) - d378*eref(k))
                      q_ref(i,j,k) = rpc(k)/(1 + rpc(k))
                   endif 
                end do
! set the tendencies to zero where you don't adjust
! set the reference profiles to be the original profiles (for diagnostic 
! purposes only --  you can think of this as what you're relaxing to in
! areas above the actual convection
                do k=1,max(klzb-1,1)
                   qdel(i,j,k) = 0.0
                   tdel(i,j,k) = 0.0
                   q_ref(i,j,k) = qin(i,j,k)
                   t_ref(i,j,k) = tin(i,j,k)
                end do
! initialize p to zero for the loop
                precip(i,j) = 0.
                precip_t(i,j) = 0.
! makes t_bm prop to (CAPE)**-.5.  Gives a relaxation time of tau_bm when 
! CAPE = sqrt(capetaubm)
                if(do_taucape) then
                   tau_bm = sqrt(capetaubm)*tau_bm/sqrt(cape1)
                   if(tau_bm.lt.tau_min) tau_bm = tau_min
                endif 
                do k=klzb, kx
! relax to reference profiles
                   tdel(i,j,k) = - (tin(i,j,k) - t_ref(i,j,k))/tau_bm*dt
                   qdel(i,j,k) = - (qin(i,j,k) - q_ref(i,j,k))/tau_bm*dt
! Precipitation can be calculated already, based on the change in q on the 
! way up (this isn't altered in the energy conservation scheme).  
                   precip(i,j) = precip(i,j) - qdel(i,j,k)*(phalf(i,j,k+1)- &
                                 phalf(i,j,k))/grav
                   precip_t(i,j)= precip_t(i,j) + cp_air/(hlv+small)*tdel(i,j,k)* &
                                 (phalf(i,j,k+1)-phalf(i,j,k))/grav
                end do
                if ((precip(i,j).gt.0.).and.(precip_t(i,j).gt.0.)) then 
! If precip > 0, then correct energy. 
                   bmflag(i,j) = 2
                   if(precip(i,j).gt.precip_t(i,j)) then
! if the q precip is greater, then lengthen the relaxation timescale on q to
! conserve energy.  qdel is therefore changed.
                      invtau_bm_q(i,j) = precip_t(i,j)/precip(i,j)/tau_bm
                      qdel(i,j,klzb:kx) = tau_bm*invtau_bm_q(i,j)* &
                         qdel(i,j,klzb:kx)
                      precip(i,j) = precip_t(i,j)
                      invtau_bm_t(i,j) = 1./tau_bm
                   else
                      if(do_simp) then
! simple scheme:
! if the t precip is greater, then lengthen the relaxation timescale on t to
! conserve energy.  tdel is therefore changed.
                         invtau_bm_t(i,j) = precip(i,j)/precip_t(i,j)/tau_bm
                         tdel(i,j,klzb:kx) = tau_bm*invtau_bm_t(i,j)* &
                            tdel(i,j,klzb:kx)
                         invtau_bm_q(i,j) = 1./tau_bm
!                      else if(do_enadjusttemp)
                      else
! not simple scheme: shift the reference profile of temperature
! deltak is the energy correction that you make to the temperature reference
! profile
                         deltak = 0.
                         do k=klzb, kx
! Calculate the integrated difference in energy change within each level.
                            deltak = deltak - (tdel(i,j,k) + hlv/cp_air*&
                                     qdel(i,j,k))* &
                                     (phalf(i,j,k+1) - phalf(i,j,k))
                         end do
! Divide by total pressure.
                         deltak = deltak/(phalf(i,j,kx+1) - phalf(i,j,klzb))
! Adjust the reference profile (uniformly with height), and correspondingly 
! the temperature change.
                         t_ref(i,j,klzb:kx) = t_ref(i,j,klzb:kx)+ &
                              deltak*tau_bm/dt
                         tdel(i,j,klzb:kx) = tdel(i,j,klzb:kx) + deltak
                      endif
                   endif
                else if(precip_t(i,j).gt.0.) then
! If precip < 0, then do the shallow conv routine.
! First option: do_shallower = true
! This chooses the depth of convection based on choosing the height that 
! it can make precip zero, i.e., subtract off heights until that precip 
! becomes positive.  
                   if (do_shallower) then
! ktop is the new top of convection.  set this initially to klzb.
                      ktop = klzb
! Work your way down until precip is positive again.
                      do while ((precip(i,j).lt.0.).and.(ktop.le.kx))
                         precip(i,j) = precip(i,j) - qdel(i,j,ktop)* &
                                  (phalf(i,j,ktop) - phalf(i,j,ktop+1))/grav
                         ktop = ktop + 1
                      end do
! since there will be an overshoot (precip is going to be greater than zero 
! once we finish this), the actual new top of convection is somewhere between
! the current ktop, and one level above this.  set ktop to the level above.
                      ktop = ktop - 1
! Adjust the tendencies in the places above back to zero, and the reference 
! profiles back to the original t,q.
                      if (ktop.gt.klzb) then
                         qdel(i,j,klzb:ktop-1) = 0.
                         q_ref(i,j,klzb:ktop-1) = qin(i,j,klzb:ktop-1)
                         tdel(i,j,klzb:ktop-1) = 0.
                         t_ref(i,j,klzb:ktop-1) = tin(i,j,klzb:ktop-1)
                      end if
! Then make the change only a fraction of the new top layer so the precip is 
! identically zero.
! Calculate the fractional penetration of convection through that top layer.  
! This is the amount necessary to make precip identically zero.  
                      if (precip(i,j).gt.0.) then 
                         ptopfrac = precip(i,j)/(qdel(i,j,ktop)* &
                            (phalf(i,j,ktop+1) - phalf(i,j,ktop)))*grav
! Reduce qdel in the top layer by this fraction. 
                         qdel(i,j,ktop) = ptopfrac*qdel(i,j,ktop)
! Set precip to zero
                         precip(i,j) = 0.
! Now change the reference temperature in such a way to make the net 
! heating zero.
!! Reduce tdel in the top layer
                         tdel(i,j,ktop) = ptopfrac*tdel(i,j,ktop)
                         deltak = 0.
                         if (ktop.lt.kx) then
! Integrate temperature tendency up to 1 level below top.
                            do k=ktop,kx
                               deltak = deltak + tdel(i,j,k)* &
                                   (phalf(i,j,k) - phalf(i,j,k+1))
                            end do
! Normalize by the pressure difference.
                            deltak = deltak/(phalf(i,j,kx+1) - phalf(i,j,ktop))
! Subtract this value uniformly from tdel, and make the according change to 
! t_ref.
                            do k=ktop,kx
                               tdel(i,j,k) = tdel(i,j,k) + deltak
                               t_ref(i,j,k) = t_ref(i,j,k) + deltak*tau_bm/dt
                            end do
                         end if
                      else
                         precip(i,j) = 0.
                         qdel(i,j,kx) = 0.
                         q_ref(i,j,kx) = qin(i,j,kx)
                         tdel(i,j,kx) = 0.
                         t_ref(i,j,kx) = tin(i,j,kx)
                         invtau_bm_t(i,j) = 0.
                         invtau_bm_q(i,j) = 0.
                      end if
                   else if(do_changeqref) then
! Change the reference profile of q by a certain fraction so that precip is 
! zero.  This involves calculating the total integrated q_ref dp (this is the
! quantity intqref), as well as the necessary change in q_ref (this is the 
! quantity deltaq).  Then the fractional change in q_ref at each level (the 
! quantity deltaqfrac) is 1-deltaq/intqref.  (have to multiply q_ref by 
! 1-deltaq/intqref at every level)  Then the change in qdel is 
! -deltaq/intqref*q_ref*dt/tau_bm.
! Change the reference profile of T by a uniform amount so that precip is zero.
                      deltak = 0.
                      deltaq = 0.
                      qrefint = 0.
                      do k=klzb,kx
! deltaq = a positive quantity (since int qdel is positive).  It's how 
! much q_ref must be changed by, in an integrated sense.  The requisite 
! change in qdel is this without the factors of tau_bm and dt.
                         deltaq = deltaq - qdel(i,j,k)*tau_bm/dt* &
                                   (phalf(i,j,k) - phalf(i,j,k+1))
! deltak = the amount tdel needs to be changed
                         deltak  = deltak  + tdel(i,j,k)* &
                                   (phalf(i,j,k) - phalf(i,j,k+1))
! qrefint = integrated value of qref
                         qrefint = qrefint - q_ref(i,j,k)* &
                                   (phalf(i,j,k) - phalf(i,j,k+1))
                      end do
! Normalize deltak by total pressure.
                      deltak  = deltak /(phalf(i,j,kx+1) - phalf(i,j,klzb))
! multiplying factor for q_ref is 1 + the ratio
                      deltaqfrac = 1. - deltaq/qrefint
! multiplying factor for qdel adds dt/tau_bm
                      deltaqfrac2 = - deltaq/qrefint*dt/tau_bm
                      precip(i,j) = 0.0
                      do k=klzb,kx
                         qdel(i,j,k) = qdel(i,j,k) + deltaqfrac2*q_ref(i,j,k)
                         q_ref(i,j,k) = deltaqfrac*q_ref(i,j,k)
                         tdel(i,j,k) = tdel(i,j,k) + deltak
                         t_ref(i,j,k) = t_ref(i,j,k) + deltak*tau_bm/dt
                      end do
                   else
                      precip(i,j) = 0.
                      tdel(i,j,:) = 0.
                      qdel(i,j,:) = 0.
                      invtau_bm_t(i,j) = 0.
                      invtau_bm_q(i,j) = 0.
                   end if
! for cases where virtual temp predicts CAPE but precip_t < 0.
                else
                   tdel(i,j,:) = 0.0
                   qdel(i,j,:) = 0.0
                   precip(i,j) = 0.0
                   q_ref(i,j,:) = qin(i,j,:)
                   t_ref(i,j,:) = tin(i,j,:)
                   invtau_bm_t(i,j) = 0.
                   invtau_bm_q(i,j) = 0.
                end if
! if no CAPE, set tendencies to zero.
             else 
                tdel(i,j,:) = 0.0
                qdel(i,j,:) = 0.0
                precip(i,j) = 0.0
                q_ref(i,j,:) = qin(i,j,:)
                t_ref(i,j,:) = tin(i,j,:)
                invtau_bm_t(i,j) = 0.
                invtau_bm_q(i,j) = 0.
             end if
          end do
       end do

       rain = precip
       snow = 0.
   

   end subroutine betts_miller

!#######################################################################

!all new cape calculation.

      subroutine capecalcnew(kx,p,phalf,cp_air,rdgas,rvgas,hlv,kappa,tin,rin,&
                             avgbl,cape,cin,tp,rp,klzb,klcl)

!
!    Input:
!
!    kx          number of levels
!    p           pressure (index 1 refers to TOA, index kx refers to surface)
!    phalf       pressure at half levels
!    cp_air      specific heat of dry air
!    rdgas       gas constant for dry air
!    rvgas       gas constant for water vapor (used in Clausius-Clapeyron, 
!                not for virtual temperature effects, which are not considered)
!    hlv         latent heat of vaporization
!    kappa       the constant kappa
!    tin         temperature of the environment
!    rin         specific humidity of the environment
!    avgbl       if true, the parcel is averaged in theta and r up to its LCL
!
!    Output:
!    cape        Convective available potential energy
!    cin         Convective inhibition (if there's no LFC, then this is set 
!                to zero)
!    tp          Parcel temperature (set to the environmental temperature 
!                where no adjustment)
!    rp          Parcel specific humidity (set to the environmental humidity 
!                where no adjustment, and set to the saturation humidity at 
!                the parcel temperature below the LCL)
!    klzb        Level of zero buoyancy
!    klcl        Lifting condensation level
!
!    Algorithm: 
!    Start with surface parcel. 
!    Calculate the lifting condensation level (uses an analytic formula and a 
!       lookup table).  
!    Average under the LCL if desired, if this is done, then a new LCL must
!       be calculated.  
!    Calculate parcel ascent up to LZB.
!    Calculate CAPE and CIN.  
      implicit none
      integer, intent(in)                    :: kx
      logical, intent(in)                    :: avgbl
      real, intent(in), dimension(:)         :: p, phalf, tin, rin
      real, intent(in)                       :: rdgas, rvgas, hlv, kappa, cp_air
      integer, intent(out)                   :: klzb, klcl
      real, intent(out), dimension(:)        :: tp, rp
      real, intent(out)                      :: cape, cin

      integer            :: k, klfc, klcl2 !  klcl
      logical            :: nocape
      real, dimension(kx)   :: theta
      real                  :: t0, r0, es, rs, theta0, pstar, value, tlcl, &
                               a, b, dtdlnp, d2tdlnp2, thetam, rm, tlcl2, &
                               plcl2, plcl, plzb, small

      pstar = 1.e5
! so we can run dry limit (one expression involves 1/hlv)
      small = 1.e-10

      nocape = .true.
      cape = 0.
      cin = 0.
      plcl = 0.
      plzb = 0.
      klfc = 0
      klcl = 0
      klzb = 0
      tp(1:kx) = tin(1:kx)
      rp(1:kx) = rin(1:kx)

! start with surface parcel (plus a buoyancy kick to temperature)
      t0 = tin(kx) + buoyancy_kick
      r0 = rin(kx)
! calculate the lifting condensation level by the following:
! are you saturated to begin with?  
      call escomp(t0,es)
      rs = rdgas/rvgas*es/p(kx)
      if (r0.ge.rs) then
! if you�re already saturated, set lcl to be the surface value.
         plcl = p(kx)
! the first level where you�re completely saturated.
         klcl = kx
! saturate out to get the parcel temp and humidity at this level
! first order (in delta T) accurate expression for change in temp
         tp(kx) = t0 + (r0 - rs)/(cp_air/(hlv+small) + hlv*rs/rvgas/t0**2.)
         call escomp(tp(kx),es)
         rp(kx) = rdgas/rvgas*es/p(kx)
      else
! if not saturated to begin with, use the analytic expression to calculate the 
! exact pressure and temperature where you�re saturated.  
         theta0 = t0*(pstar/p(kx))**kappa
! the expression that we utilize is 
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! (the division by es00 is necessary because the RHS values are tabulated
! for control moisture content)
! The right hand side of this is only a function of temperature, therefore 
! this is put into a lookup table to solve for temperature.  
         if (r0.gt.0.) then
            value = log(theta0**(-1/kappa)*r0*pstar*rvgas/rdgas/es0)
            call lcltabl(value,tlcl)
            plcl = pstar*(tlcl/theta0)**(1/kappa)
! just in case plcl is very high up
            if (plcl.lt.p(1)) then
               plcl = p(1)
               tlcl = theta0*(plcl/pstar)**kappa
               write (*,*) 'hi lcl'
            end if
            k = kx
         else
! if the parcel sp hum is zero or negative, set lcl to top level
            plcl = p(1)
            tlcl = theta0*(plcl/pstar)**kappa
!            write (*,*) 'zero r0', r0
            do k=1,kx
               tp(k) = theta0*(p(k)/pstar)**kappa
               rp(k) = 0.
! this definition of CIN contains everything below the LCL
               cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            end do
            go to 11
         end if
! calculate the parcel temperature (adiabatic ascent) below the LCL.  
! the mixing ratio stays the same
         do while (p(k).gt.plcl)
            tp(k) = theta0*(p(k)/pstar)**kappa
            call escomp(tp(k),es)
            rp(k) = rdgas/rvgas*es/p(k)
! this definition of CIN contains everything below the LCL
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            k = k-1
         end do
! first level where you're saturated at the level
         klcl = k
	 if (klcl.eq.1) klcl = 2
! do a saturated ascent to get the parcel temp at the LCL.  
! use your 2nd order equation up to the pressure above.  
! moist adaibat derivatives: (use the lcl values for temp, humid, and 
! pressure)
         a = kappa*tlcl + hlv/cp_air*r0
         b = hlv**2.*r0/cp_air/rvgas/tlcl**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
! second order in p (RK2)
! first get temp halfway up 
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)/2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call escomp(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/(p(klcl) + plcl)*2.
         a = kappa*tp(klcl) + hlv/cp_air*rp(klcl)
         b = hlv**2./cp_air/rvgas*rp(klcl)/tp(klcl)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
!         d2tdlnp2 = (kappa + b - 1. - b/tlcl*(hlv/rvgas/tlcl - &
!                   2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*r0/cp_air/ &
!                   (1. + b)
! second order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl) + .5*d2tdlnp2*(log(&
!             p(klcl)/plcl))**2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call escomp(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/p(klcl)
!         write (*,*) 'tp, rp klcl:kx, new', tp(klcl:kx), rp(klcl:kx)
! CAPE/CIN stuff
         if ((tp(klcl).lt.tin(klcl)).and.nocape) then
! if you�re not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(klcl) - &
                 tp(klcl))*log(phalf(klcl+1)/phalf(klcl))
         else
! if you�re buoyant, then add to cape
            cape = cape + rdgas*(tp(klcl) - &
                  tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
! if it�s the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = klcl
            endif
         end if
      end if
! then average the properties over the boundary layer if so desired.  to give 
! a new "parcel".  this may not be saturated at the LCL, so make sure you get 
! to a level where it is before moist adiabatic ascent!
!!!! take out all the below (between the exclamation points) if no avgbl !!!!
      if (avgbl) then
         theta(klcl:kx) = tin(klcl:kx)*(pstar/p(klcl:kx))**kappa
         thetam = 0.
         rm = 0.
         do k=klcl,kx
            thetam = thetam + theta(k)*(phalf(k+1) - phalf(k))
            rm = rm + rin(k)*(phalf(k+1) - phalf(k))
         end do
         thetam = thetam/(phalf(kx+1) - phalf(klcl))
         rm = rm/(phalf(kx+1) - phalf(klcl))
! check if you're saturated at the top level.  if not, then get a new LCL
         tp(klcl) = thetam*(p(klcl)/pstar)**kappa
         call escomp(tp(klcl),es)
         rs = rdgas/rvgas*es/p(klcl)
! if you're not saturated, get a new LCL
         if (rm.lt.rs) then
! reset CIN to zero.  
            cin = 0.
! again, use the analytic expression to calculate the exact pressure and 
! temperature where you�re saturated.  
! the expression that we utilize is 
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00)= log(es/T**(1/kappa))
! (the division by es00 is necessary because the RHS values are tabulated
! for control moisture content)
! The right hand side of this is only a function of temperature, therefore 
! this is put into a lookup table to solve for temperature.  
            value = log(thetam**(-1/kappa)*rm*pstar*rvgas/rdgas/es0)
            call lcltabl(value,tlcl2)
            plcl2 = pstar*(tlcl2/thetam)**(1/kappa)
! just in case plcl is very high up
            if (plcl2.lt.p(1)) then
               plcl2 = p(1)
            end if
            k = kx
! calculate the parcel temperature (adiabatic ascent) below the LCL.  
! the mixing ratio stays the same
            do while (p(k).gt.plcl2) 
               tp(k) = thetam*(p(k)/pstar)**kappa
               call escomp(tp(k),es)
               rp(k) = rdgas/rvgas*es/p(k)
! this definition of CIN contains everything below the LCL
               cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
               k = k-1
            end do
! first level where you�re saturated at the level
            klcl2 = k
	    if (klcl2.eq.1) klcl2 = 2
! do a saturated ascent to get the parcel temp at the LCL.  
! use your 2nd order equation up to the pressure above.  
! moist adaibat derivatives: (use the lcl values for temp, humid, and 
! pressure)
            a = kappa*tlcl2 + hlv/cp_air*rm
            b = hlv**2.*rm/cp_air/rvgas/tlcl2**2.
            dtdlnp = a/(1. + b)
! first order in p
!            tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)
! second order in p (RK2)
! first get temp halfway up 
         tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)/2.
         if ((tp(klcl2).lt.173.16).and.nocape) go to 11
         call escomp(tp(klcl2),es)
         rp(klcl2) = rdgas/rvgas*es/(p(klcl2) + plcl2)*2.
         a = kappa*tp(klcl2) + hlv/cp_air*rp(klcl2)
         b = hlv**2./cp_air/rvgas*rp(klcl2)/tp(klcl2)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)
!            d2tdlnp2 = (kappa + b - 1. - b/tlcl2*(hlv/rvgas/tlcl2 - &
!                          2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*rm/cp_air/ &
!                          (1. + b)
! second order in p
!            tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2) + &
!               .5*d2tdlnp2*(log(p(klcl2)/plcl2))**2.
            call escomp(tp(klcl2),es)
            rp(klcl2) = rdgas/rvgas*es/p(klcl2)
! CAPE/CIN stuff
            if ((tp(klcl2).lt.tin(klcl2)).and.nocape) then
! if you�re not yet buoyant, then add to the CIN and continue
               cin = cin + rdgas*(tin(klcl2) - &
                    tp(klcl2))*log(phalf(klcl2+1)/phalf(klcl2))
            else
! if you�re buoyant, then add to cape
               cape = cape + rdgas*(tp(klcl) - &
                     tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
! if it�s the first time buoyant, then set the level of free convection to k
               if (nocape) then
                  nocape = .false.
                  klfc = klcl2
               endif
            end if
         end if
      end if
!!!! take out all of the above (within the exclamations) if no avgbl !!!!
! then, start at the LCL, and do moist adiabatic ascent by the first order 
! scheme -- 2nd order as well
      do k=klcl-1,1,-1
         a = kappa*tp(k+1) + hlv/cp_air*rp(k+1)
         b = hlv**2./cp_air/rvgas*rp(k+1)/tp(k+1)**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
! second order in p (RK2)
! first get temp halfway up 
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))/2.
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call escomp(tp(k),es)
         rp(k) = rdgas/rvgas*es/(p(k) + p(k+1))*2.
         a = kappa*tp(k) + hlv/cp_air*rp(k)
         b = hlv**2./cp_air/rvgas*rp(k)/tp(k)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
!         d2tdlnp2 = (kappa + b - 1. - b/tp(k+1)*(hlv/rvgas/tp(k+1) - & 
!               2.)*dtdlnp)/(1. + b)*dtdlnp - hlv/cp_air*rp(k+1)/(1. + b)
! second order in p
!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1)) + .5*d2tdlnp2*(log( &
!             p(k)/p(k+1)))**2.
! if you're below the lookup table value, just presume that there's no way 
! you could have cape and call it quits
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call escomp(tp(k),es)
         rp(k) = rdgas/rvgas*es/p(k)
         if ((tp(k).lt.tin(k)).and.nocape) then
! if you�re not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
         elseif((tp(k).lt.tin(k)).and.(.not.nocape)) then
! if you have CAPE, and it�s your first time being negatively buoyant, 
! then set the level of zero buoyancy to k+1, and stop the moist ascent
            klzb = k+1
            go to 11
         else
! if you�re buoyant, then add to cape
            cape = cape + rdgas*(tp(k) - tin(k))*log(phalf(k+1)/phalf(k))
! if it�s the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = k
            endif
         end if
      end do
 11   if(nocape) then 
! this is if you made it through without having a LZB
! set LZB to be the top level.
         plzb = p(1)
         klzb = 0
         klfc = 0
         cin = 0.
         tp(1:kx) = tin(1:kx)
         rp(1:kx) = rin(1:kx)
      end if
!      write (*,*) 'plcl, klcl, tlcl, r0 new', plcl, klcl, tlcl, r0
!      write (*,*) 'tp, rp new', tp, rp
!       write (*,*) 'tp, new', tp
!       write (*,*) 'tin new', tin 
!       write (*,*) 'klcl, klfc, klzb new', klcl, klfc, klzb
      end subroutine capecalcnew

! lookup table for the analytic evaluation of LCL
      subroutine lcltabl(value,tlcl)
!
! Table of values used to compute the temperature of the lifting condensation
! level.  
! 
! the expression that we utilize is 
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! the RHS is tabulated for control moisture content, hence the division 
! by es00 on the LHS
! 
! Gives the values of the temperature for the following range: 
!   starts with -23, is uniformly distributed up to -10.4.  There are a 
! total of 127 values, and the increment is .1.  
!
      implicit none 
      real, intent(in)     :: value
      real, intent(out)    :: tlcl

      integer              :: ival
      real, dimension(127) :: lcltable
      real                 :: v1, v2

      data lcltable/  1.7364512e+02,   1.7427449e+02,   1.7490874e+02, &
      1.7554791e+02,   1.7619208e+02,   1.7684130e+02,   1.7749563e+02, &
      1.7815514e+02,   1.7881989e+02,   1.7948995e+02,   1.8016539e+02, &
      1.8084626e+02,   1.8153265e+02,   1.8222461e+02,   1.8292223e+02, &
      1.8362557e+02,   1.8433471e+02,   1.8504972e+02,   1.8577068e+02, &
      1.8649767e+02,   1.8723077e+02,   1.8797006e+02,   1.8871561e+02, &
      1.8946752e+02,   1.9022587e+02,   1.9099074e+02,   1.9176222e+02, &
      1.9254042e+02,   1.9332540e+02,   1.9411728e+02,   1.9491614e+02, &
      1.9572209e+02,   1.9653521e+02,   1.9735562e+02,   1.9818341e+02, &
      1.9901870e+02,   1.9986158e+02,   2.0071216e+02,   2.0157057e+02, &
      2.0243690e+02,   2.0331128e+02,   2.0419383e+02,   2.0508466e+02, &
      2.0598391e+02,   2.0689168e+02,   2.0780812e+02,   2.0873335e+02, &
      2.0966751e+02,   2.1061074e+02,   2.1156316e+02,   2.1252493e+02, &
      2.1349619e+02,   2.1447709e+02,   2.1546778e+02,   2.1646842e+02, &
      2.1747916e+02,   2.1850016e+02,   2.1953160e+02,   2.2057364e+02, &
      2.2162645e+02,   2.2269022e+02,   2.2376511e+02,   2.2485133e+02, &
      2.2594905e+02,   2.2705847e+02,   2.2817979e+02,   2.2931322e+02, &
      2.3045895e+02,   2.3161721e+02,   2.3278821e+02,   2.3397218e+02, &
      2.3516935e+02,   2.3637994e+02,   2.3760420e+02,   2.3884238e+02, &
      2.4009473e+02,   2.4136150e+02,   2.4264297e+02,   2.4393941e+02, &
      2.4525110e+02,   2.4657831e+02,   2.4792136e+02,   2.4928053e+02, &
      2.5065615e+02,   2.5204853e+02,   2.5345799e+02,   2.5488487e+02, &
      2.5632953e+02,   2.5779231e+02,   2.5927358e+02,   2.6077372e+02, &
      2.6229310e+02,   2.6383214e+02,   2.6539124e+02,   2.6697081e+02, &
      2.6857130e+02,   2.7019315e+02,   2.7183682e+02,   2.7350278e+02, &
      2.7519152e+02,   2.7690354e+02,   2.7863937e+02,   2.8039954e+02, &
      2.8218459e+02,   2.8399511e+02,   2.8583167e+02,   2.8769489e+02, &
      2.8958539e+02,   2.9150383e+02,   2.9345086e+02,   2.9542719e+02, &
      2.9743353e+02,   2.9947061e+02,   3.0153922e+02,   3.0364014e+02, &
      3.0577420e+02,   3.0794224e+02,   3.1014515e+02,   3.1238386e+02, &
      3.1465930e+02,   3.1697246e+02,   3.1932437e+02,   3.2171609e+02, &
      3.2414873e+02,   3.2662343e+02,   3.2914139e+02,   3.3170385e+02 /

      v1 = value
      if (value.lt.-23.0) v1 = -23.0
      if (value.gt.-10.4) v1 = -10.4
      ival = floor(10.*(v1 + 23.0))
      v2 = -230. + ival
      v1 = 10.*v1
      tlcl = (v2 + 1.0 - v1)*lcltable(ival+1) + (v1 - v2)*lcltable(ival+2)


      end subroutine lcltabl

!#######################################################################

   subroutine betts_miller_init ()

!-----------------------------------------------------------------------
!
!        initialization for betts_miller
!
!-----------------------------------------------------------------------

  integer  unit,io,ierr

!----------- read namelist ---------------------------------------------

      if (file_exist('input.nml')) then
         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=betts_miller_nml, iostat=io, end=10)
            ierr = check_nml_error (io,'betts_miller_nml')
         enddo
  10     call close_file (unit)
      endif

!---------- output namelist --------------------------------------------

      call write_version_number ( version, tagname )
      if ( mpp_pe() == mpp_root_pe() ) then
           write (stdlog(),nml=betts_miller_nml)
      endif

   !s initialise here as rdgas no longer a parameter
   d622 = rdgas/rvgas
   d378 = 1.-d622

      do_init=.false.

   end subroutine betts_miller_init

!#######################################################################

end module betts_miller_mod

