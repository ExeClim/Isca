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

module qe_moist_convection_mod
  
  !----------------------------------------------------------------------
  ! This module implements the simple quasi-equilibrium convection scheme
  ! described in Frierson, "The Dynamics of Idealized Convection
  ! Schemes and Their Effect on the Zonally Averaged Tropical
  ! Circulation", J. Atmos. Sci., 64 (2007).
  ! It uses the "shallower" shallow convection scheme described in
  ! that paper, and it incorporates minor modifications mentioned in
  ! O'Gorman and Schneider, "The Hydrological Cycle over a Wide Range
  ! of Climates Simulated with an Idealized GCM", J. Atmos. Sci. 65
  ! (2008). (The modifications primarily concern a consistent use of
  ! virtual temperature effects in the convection scheme and not
  ! approximating the relation between vapor pressure and specific
  ! humidity.) 
  !----------------------------------------------------------------------
  use            fms_mod, only:  file_exist, error_mesg, open_file,  &
                                 check_nml_error, mpp_pe, FATAL,  &
                                 close_file
  use sat_vapor_pres_mod, only:  escomp, descomp
  use      constants_mod, only:  HLv, HLs, Cp_air, Grav, rdgas, rvgas, &
                                 kappa

  implicit none
  private
  !---------------------------------------------------------------------
  !  ---- public interfaces ----

  public  qe_moist_convection, qe_moist_convection_init, qe_moist_convection_end

  !-----------------------------------------------------------------------
  !   ---- version number ----

  character(len=128) :: version = '$Id: qe_moist_convection.F90,v 1.1.2.1 2013/01/24 14:46:13 pjp Exp $'
  character(len=128) :: tag = '$Name:  $'

  !-----------------------------------------------------------------------
  !   ---- local/private data ----

  logical :: do_init=.true.

  !-----------------------------------------------------------------------
  !   --- parameters and defaults (overriden in namelist) ----

  real    :: tau_bm  = 7200.
  real    :: rhbm    = .8
  real    :: Tmin    = 173. ! minimum
  real    :: Tmax    = 335. ! and maximum temperature at LCL
  real    :: val_inc = 0.01
  real    :: val_min = -1.  ! calculated in get_lcl_temp_table_size
  real    :: val_max = 1.   ! calculated in get_lcl_temp_table_size

  real, parameter  :: small = 1.e-10, &  ! to avoid division by 0 in dry limit
                      pref  = 1.e5

  real, allocatable, dimension(:) :: lcl_temp_table

  namelist /qe_moist_convection_nml/  tau_bm, rhbm, Tmin, Tmax, val_inc

  !------------------------------------------------------------------
  !           Description of namelist variables
  !
  !  tau_bm    = Betts-Miller relaxation timescale (seconds)
  !
  !  rhbm      = reference relative humidity toward which convection relaxes
  !
  !  val_inc   = increment in value for the lcl_temp_table
  !
  !  Tmin      = temperature minimum resolvable with lookup table
  !
  !  Tmax      = temperature maximum resolvable with lookup table
  !
  !  val_min   = minimum value passed to get_lcl_temp
  !              calculated in qe_moist_convection_init
  !
  !  val_min   = maximum value passed to get_lcl_temp
  !              calculated in the qe_moist_convection_init
  !-----------------------------------------------------------------------
  
contains


  !######################################################################
  subroutine qe_moist_convection_init ()

    !-----------------------------------------------------------------------
    !
    !        initialization of QE moist convection scheme 
    !
    !-----------------------------------------------------------------------
  
    integer  lcl_temp_table_size, unit, io, ierr

    !----------- read namelist ---------------------------------------------

    if (file_exist('input.nml')) then
       unit = open_file (file='input.nml', action='read')
       ierr = 1
       do while (ierr /= 0)
          read  (unit, nml=qe_moist_convection_nml, iostat=io, end=10)
          ierr = check_nml_error(io, 'qe_moist_convection_nml')
       end do
10     call close_file(unit)
    endif

    !---------- output namelist --------------------------------------------

    unit = open_file (file='logfile.out', action='append')
    if ( mpp_pe() == 0 ) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=qe_moist_convection_nml)
    endif
    call close_file(unit)

    do_init = .false.
    
    ! Calculates the size of the LCL lcl_temp_table with values Tmin Tmax
    call get_lcl_temp_table_size(lcl_temp_table_size)

    ! Generate the lcl_temp_table
    allocate (lcl_temp_table(lcl_temp_table_size))
    call generate_lcl_table(lcl_temp_table)
       
  end subroutine qe_moist_convection_init

  !#######################################################################
  
  subroutine get_lcl_temp_table_size(lcl_temp_table_size)
    
    integer, intent(out)      :: lcl_temp_table_size
    
    call get_val_min_max()
    lcl_temp_table_size = ceiling( (val_max-val_min)/val_inc )
    
  end subroutine get_lcl_temp_table_size
  !#######################################################################
  
  subroutine get_val_min_max()
    
    real      :: esmin, esmax
    
    call escomp(Tmin, esmin)
    call escomp(Tmax, esmax)
    
    val_min = log(esmin/(Tmin**(1.0/kappa))) 
    val_max = log(esmax/(Tmax**(1.0/kappa)))

  end subroutine get_val_min_max

  !##############################################################################

  subroutine generate_lcl_table(lcl_temp_table)
    
    real, intent(out), dimension(:)   :: lcl_temp_table
    real               :: lcl_temp_guess
    integer            :: k
    
    lcl_temp_guess = Tmin
    do k=1, size(lcl_temp_table, 1)
       lcl_temp_table(k) = lcl_temp(val_min + (k-1)*val_inc, lcl_temp_guess)
       lcl_temp_guess     = lcl_temp_table(k)
    end do
    
  end subroutine generate_lcl_table

  !#######################################################################

  subroutine qe_moist_convection (dt, Tin, qin, p_full, p_half, coldT, &
                                  rain, snow, deltaT, deltaq, qref, convflag, &
                                  kLZBs, CAPE, CIN, invtau_q_relaxation, &
                                  invtau_t_relaxation, Tref, kLCLs)

    !-----------------------------------------------------------------------
    !
    !    Minimal quasi-equilibrium (Betts-Miller) convection scheme
    !
    !-----------------------------------------------------------------------
    !
    !   input:  dt       Time step in seconds
    !           Tin      Temperature at full model levels
    !           qin      Specific humidity of water vapor at full
    !                      model levels
    !           p_full   Pressure at full model levels
    !           p_half   Pressure at half (interface) model levels
    !           coldT    Flag indicating whether precipitation should be snow
    !                      (not used)
    !
    !  output:  rain     Liquid precipitation (kg/m2)
    !           snow     Frozen precipitation (kg/m2)
    !           delta_T  Temperature tendency at full model levels
    !           delta_q  Specific humidity tendency (of water vapor) at
    !                      full model levels
    !           convflag Flag indicating what kind of convection occurred:
    !                convflag = 0: no cape, no convection
    !                convflag = 1: shallow conv; predicted precip less than zero
    !                convflag = 2: deep convection
    !           kLZBs    Levels of zero buoyancy
    !           CAPE     Convectively available potential energy 
    !           CIN      Convective inhibition (this and the above are before the 
    !                      adjustment)
    !           invtau_q_relaxation
    !                    Temperature relaxation time scale (1/s)
    !           invtau_t_relaxation
    !                    Humidity relaxation time scale    (1/s)
    !           Tref     Reference temperature profile 
    !
    !--------------------- interface arguments -----------------------------

    real   , intent(in) , dimension(:,:,:) :: tin, qin, p_full, p_half
    real   , intent(in)                    :: dt
    logical, intent(in) , dimension(:,:)   :: coldT
    
    real   , intent(out), dimension(:,:)   :: rain, snow, CAPE, CIN
    real   , intent(out), dimension(:,:)   :: invtau_q_relaxation, invtau_t_relaxation
    integer, intent(out), dimension(:,:)   :: convflag, kLZBs, kLCLs
    real   , intent(out), dimension(:,:,:) :: deltaT, deltaq, qref, Tref

    !-----------------------------------------------------------------------
    !     computation of precipitation by convection scheme
    !----------------------------------------------------------------

    ! Check whether initialization has been completed
    if (do_init) call error_mesg ('qe_moist_convection',  &
         'qe_moist_convection_init has not been called.', FATAL)
    
    ! Call the convection scheme itself
    call SBM_convection_scheme(dt, Tin, qin, p_full, p_half, rain, snow, &
         deltaT, deltaq, kLZBs, CAPE, CIN,invtau_q_relaxation,           &
         invtau_t_relaxation, Tref, qref,                                &
         val_min, val_max, val_inc, lcl_temp_table, convflag, kLCLs)
     
  end subroutine qe_moist_convection

  !#######################################################################
  
  subroutine SBM_convection_scheme(dt, Tin, qin, p_full, p_half, rain,  snow, &
       deltaT, deltaq, kLZBs, CAPE, CIN, invtau_q_relaxation, invtau_t_relaxation,&
       Tref, qref, val_min, val_max, val_inc, lcl_temp_table, convflag, kLCLs)

   
    !-----------------------------------------------------------------------
    !
    !                     SBM Convection Scheme
    !
    !-----------------------------------------------------------------------
    !
    ! Inputs and outputs as for qe_moist_convection, plus arguments
    ! (val_min, val_max, val_inc) pertaining to the lookup of the LCL
    ! temperature from the lookup table lcl_temp_table calculated at
    ! initialization 

    real, intent(in)                       :: dt, val_min, val_max, val_inc
    real, intent(in), dimension(:,:,:)     :: Tin, qin, p_full
    real, intent(in), dimension(:,:,:)     :: p_half
    real, intent(in), dimension(:)         :: lcl_temp_table
    real, intent(out), dimension(:,:)      :: rain, snow, CAPE, CIN
    real, intent(out), dimension(:,:)      :: invtau_q_relaxation, invtau_t_relaxation
    real, intent(out), dimension(:,:,:)    :: deltaT, deltaq, Tref, qref
    integer, intent(out), dimension(:,:)   :: kLZBs, convflag, kLCLs
    
    integer                                :: k_surface, i, j, kLZB, kLCL
    real, dimension(size(Tin, 3))          ::     &
         deltaq_parcel, deltaT_parcel, T_parcel, r_parcel, qref_parcel, Tref_parcel
    real, dimension(size(Tin, 1), size(Tin, 2))            :: Pq
    real, dimension(size(Tin,1), size(Tin,2), size(Tin,3)) :: rin
    real                                   :: cape_parcel, cin_parcel, Pq_parcel, Pt_parcel
    real                                   :: invtau_q_relaxation_parcel, invtau_t_relaxation_parcel

    ! Initialization of parameters and variables
    k_surface = size(Tin, 3)
    deltaq    = 0.
    deltaT    = 0.
    Pq        = 0.
    rin       = qin / (1.0 - qin)

    ! Loop over latitude and longitude
    do i=1, size(Tin, 1)
       do j=1, size(Tin, 2)
          
          ! Definition of variables used for the parcel(i,j)
          deltaq_parcel = deltaq(i, j, :)
          deltaT_parcel = deltaT(i, j, :)
          invtau_q_relaxation = 0.
          invtau_t_relaxation = 0.
          
          convflag(i,j) = 0
          invtau_q_relaxation_parcel = 0.
          invtau_t_relaxation_parcel = 0.
          
          T_parcel = Tin(i,j,:)
          r_parcel = rin(i,j,:)
          
          ! Calculate CAPE (Convective Available Potential Energy) of 
          ! parcel lifted from lowest model level
          call CAPE_calculation(k_surface, p_full(i,j,:), p_half(i,j,:), &
               Tin(i,j,:), rin(i,j,:), kLZB, T_parcel, r_parcel,         &
               cape_parcel, cin_parcel, val_min, val_max, lcl_temp_table,&
               kLCL)
            
          ! Store values
          CAPE(i,j)  = cape_parcel
          CIN(i,j)   = cin_parcel
          kLZBs(i,j) = kLZB
          kLCLs(i,j) = kLCL
          
          ! If CAPE>0, set reference temperature and humidity above and below 
          ! the LZB (Level of Zero Buoyancy) 
          if (cape_parcel .gt. 0) then
             ! moist convection may occur
             convflag(i,j) = 1   ! flag for positive CAPE
             call set_reference_profiles(p_full(i,j,:), qin(i,j,:), Tin(i,j,:), &
                  T_parcel, kLZB, k_surface, r_parcel, deltaq_parcel,   & 
                  deltaT_parcel, qref_parcel, Tref_parcel)
             
             ! Calculate the precipitation rate Pq 
             call Pq_calculation(kLZB, k_surface, qref_parcel, qin(i,j,:),    &
                  p_half(i,j,:), deltaq_parcel, Pq_parcel, dt)
             ! Calculate the humidity change that would be necessary
             ! to balance temperature change by latent heat release 
             call Pt_calculation(kLZB, k_surface, Tref_parcel, Tin(i,j,:),    &
                  p_half(i,j,:), deltaT_parcel, Pt_parcel, dt)
             
             ! If Pq > 0 and Pt > 0, do deep convection
             if ( (Pq_parcel .gt. 0) .and. (Pt_parcel .gt. 0) ) then
                convflag(i,j) = 2 ! deep moist convection

                 call do_deep_convection (kLZB, k_surface,Pt_parcel, dt,p_half(i,j,:),&
                     invtau_q_relaxation_parcel, invtau_t_relaxation_parcel,Pq_parcel,&
                     deltaT_parcel,Tref_parcel,deltaq_parcel)

             else 
                ! Else, if Pq <= 0 and Pt > 0
                if (Pt_parcel .gt. 0) then
                   ! DO shallow convection
                   call do_shallow_convection(kLZB, k_surface,qin(i,j,:),   &
                        qref_parcel, deltaq_parcel, Tin(i,j,:), Tref_parcel, &
                        deltaT_parcel, p_half(i,j,:), Pq_parcel,dt)
                else
                   ! Else, do nothing, and go back to loop over latitude and longitude 
                   Pq_parcel     = 0.
                   call set_profiles_to_full_model_values (1, k_surface, Tin(i,j,:),&
                        qin(i,j,:), Tref_parcel, deltaT_parcel, qref_parcel, &
                        deltaq_parcel)
                end if
             end if
          else
             ! If CAPE < 0, do nothing, and go back to loop over latitude and longitude
             Pq_parcel     = 0.
             call set_profiles_to_full_model_values (1, k_surface, Tin(i,j,:),&
                        qin(i,j,:), Tref_parcel, deltaT_parcel, qref_parcel, &
                        deltaq_parcel)
          end if

          ! Store diagnostics
          deltaT(i,j,:) = deltaT_parcel
          deltaq(i,j,:) = deltaq_parcel
          Pq(i,j)       = Pq_parcel
          qref(i,j,:)   = qref_parcel
          Tref(i,j,:)   = Tref_parcel
          invtau_q_relaxation(i,j)=invtau_q_relaxation_parcel
          invtau_t_relaxation(i,j)=invtau_t_relaxation_parcel


       end do
    end do

    rain = Pq
    snow = 0.
   
  end subroutine sbm_convection_scheme

  !#######################################################################
  
  subroutine CAPE_calculation(k_surface, p_full, p_half, Tin, rin, kLZB, &
       Tp, rp, CAPE, CIN, val_min, val_max, lcl_temp_table, kLCL)
       

    ! Calculates CAPE, CIN, level of zero buoyancy, and parcel properties 
    ! (second order accurate in delta(ln p) and exact LCL calculation)

    integer, intent(in)                  :: k_surface
    real, intent(in), dimension(:)       :: p_full
    real, intent(in), dimension(:)       :: p_half
    real, intent(in), dimension(:)       :: Tin, rin
    real, intent(in)                     :: val_min , val_max
    real, intent(in), dimension(:)       :: lcl_temp_table
    integer, intent(out)                 :: kLZB, kLCL
    real, intent(out), dimension(:)      :: Tp, rp
    real, intent(out)                    :: CAPE, CIN
    
    logical                              :: nocape, saturated, skip
    real                                 :: pLZB, T0, r0, es, rs, pLCL
    integer                              :: kLFC, k
    real, dimension(size(Tin))           :: Tin_virtual
    
    nocape = .true.
    CAPE   = 0.
    CIN    = 0.
    pLZB   = 0.
    kLFC   = 0
    kLZB   = 0
    Tp     = Tin
    rp     = rin
    saturated = .false.
    kLCL = 0
    
    ! Calculation of values to check whether the lowest level is saturated
    ! Calculate the virtual temperature
    do k = 1,k_surface
       Tin_virtual(k) = virtual_temp(Tin(k), rin(k))
    end do
    
    ! Definition of the temperature and the mixing ratio at the surface
    T0 = Tin(k_surface)
    r0 = rin(k_surface)
    
    call escomp(T0,es)
    
    ! Calculates the saturated mixing ratio at the surface
    rs = mixing_ratio(es, p_full(k_surface))
    
    ! Is the lowest level saturated or oversaturated?
    if (r0 .ge. rs) then
       saturated = .true.
    end if
    
    ! Calculation below the lifted condensation level LCL
    call CAPE_below_LCL(saturated, k_surface, p_half, p_full, Tin, Tin_virtual, &
         Tp, T0, r0, rs, rp, rin, CAPE, nocape, CIN, pLZB, kLFC, kLZB,     &
         pLCL, kLCL, skip, val_min, val_max, lcl_temp_table)
    
    ! Calculation above the LCL
    call CAPE_above_LCL(kLCL, kLZB, kLFC, Tp, rp, rin, p_full, nocape, skip,&
         CIN, CAPE, Tin, Tin_virtual, p_half, pLZB)
         
    
  end subroutine CAPE_CALCULATION
  
  !#######################################################################

  subroutine CAPE_below_LCL(saturated, k_surface,  p_half, p_full, Tin, Tin_virtual, &
                       Tp, T0, r0, rs, rp, rin, CAPE, nocape, CIN, pLZB, kLFC, kLZB, &
                       pLCL, kLCL, skip, val_min, val_max, lcl_temp_table)

    logical, intent(in)                :: saturated
    integer, intent(in)                :: k_surface
    real, intent(in), dimension(:)     :: p_half
    real, intent(in), dimension(:)     :: p_full
    real, intent(in), dimension(:)     :: Tin, Tin_virtual
    real, intent(in)                   :: T0
    real, intent(in)                   :: r0, rs
    real, intent(in), dimension(:)     :: rin
    real, intent(in)                   :: val_min, val_max
    real, intent(in), dimension(:)     :: lcl_temp_table
    real, intent(inout), dimension(:)  :: rp, Tp
    real, intent(inout)                :: CAPE, CIN
    logical, intent(inout)             :: nocape
    real, intent(inout)                :: pLZB
    integer, intent(inout)             :: kLFC, kLZB
    logical, intent(out)               :: skip
    real, intent(out)                  :: pLCL
    integer, intent(out)               :: kLCL
    
    real                               :: theta0, value, a, b , dtdlnp, es
    integer                            :: k
    real                               :: TLCL     
    
    skip = .false.
  
    ! If the surface is already (over-)saturated
    if (saturated) then
       pLCL = p_full(k_surface)
       kLCL = k_surface
       
       ! Saturate parcel (wring out excess moisture and change temperature 
       ! correspondinly; the following is the resulting first-order
       ! change in temperature)
       Tp(k_surface) =       & 
            T0 + (r0-rs) / ( (Cp_air/(HLv+small)) + (HLv*rs)/rvgas/T0**2 )
       call escomp(Tp(k_surface), es)
       rp(k_surface) = mixing_ratio(es, p_full(k_surface))
    else
       ! If the lowest level is not saturated, calculate temperature of saturation
       theta0 = Tin(k_surface) * (pref/p_full(k_surface))**kappa
       
       if (r0 .le. 0) then
          ! If the mixing ratio r0 <= 0, LCL is the top of model
          pLCL = p_full(1)
          TLCL = theta0 * (pLCL/pref)**kappa
          skip = .true.
       else
          ! If the mixing ratio r0 > 0, calculate LCL temperature and temperature
          value = log( theta0**(-1/kappa) * pref*r0 / (rdgas/rvgas + r0) )
          
          call get_lcl_temp(lcl_temp_table, value, val_min, val_max, TLCL)
          pLCL = pref * (TLCL/theta0)**(1./kappa)
         
          if (pLCL .lt. p_full(1)) then
             ! If the pLCL is above model domain, use values at top of the model
             pLCL = p_full(1)
             TLCL = theta0 * (pLCL/pref)**kappa
          end if
          
          ! Calculate parcel temperature and CIN below LCL by upward integration
          k   = k_surface
          CIN = 0.
          do while (p_full(k) .gt. pLCL)
             Tp(k) = theta0 * (p_full(k)/pref)**kappa
             call escomp(Tp(k), es)
             
             ! rp is not the actual mixing ratio but will be used
             ! to calculate the reference moisture profile using rhbm
             rp(k) = mixing_ratio(es, p_full(k))
             CIN   = CIN  &
                  + rdgas*( Tin_virtual(k) - virtual_temp(Tp(k), r0) ) &
                    * log(p_half(k+1)/p_half(k))
             k     = k - 1
          end do

          kLCL = k
          
          ! Temperature profile for saturated ascent
          a = kappa * TLCL + (HLv/Cp_air)*r0
          b = (HLv**2) * r0 / (Cp_air*rvgas*TLCL**2)
          dtdlnp = a / (1.0 + b)
          
          ! Second order in p (RK2): first get temperature halfway up
          Tp(kLCL) = TLCL + dtdlnp * log(p_full(kLCL)/pLCL)/2
          if ( (Tp(kLCL) .lt. Tmin) .and. nocape ) then
             skip = .true.
             if (nocape) then
                call set_values_if_nocape (Tin, rin, p_full,&
                     Tp, rp, pLZB, kLZB, kLFC, CIN )
             end if
          else
             call escomp(Tp(kLCL),es)
             rp(kLCL) = mixing_ratio(es, (p_full(kLCL) + pLCL)/2)
             a = kappa * Tp(kLCL) + (HLv/Cp_air) * rp(kLCL)
             b = (HLv**2)*rp(kLCL) / (Cp_air * rvgas * Tp(kLCL)**2)
             dtdlnp = a/(1.0 + b)
             ! Second half of RK2
             Tp(kLCL) = TLCL + dtdlnp * log(p_full(kLCL)/pLCL)
             if ( (Tp(kLCL) .lt. Tmin) .and. nocape) then
                skip = .true.
                if (nocape) then
                   call set_values_if_nocape (Tin, rin, p_full, &
                        Tp, rp, pLZB, kLZB, kLFC, CIN )
                end if
             else
                call escomp(Tp(kLCL), es)
                rp(kLCL) = mixing_ratio(es, p_full(kLCL))
                
                if ((virtual_temp(Tp(kLCL), rp(kLCL)) .lt. Tin_virtual(kLCL)) .and. nocape) then
                   ! If the parcel is not buoyant yet, add to CIN
                   CIN = CIN + rdgas * (Tin_virtual(kLCL) - virtual_temp(Tp(kLCL), rp(kLCL))) &
                        * log(p_half(kLCL+1)/p_half(kLCL))
                else
                   ! If the parcel is buoyant, add to CAPE
                   CAPE = CAPE + rdgas * (virtual_temp(Tp(kLCL), rp(kLCL)) - Tin_virtual(kLCL)) &
                        * log(p_half(kLCL+1)/p_half(kLCL))
                   if (nocape) then
                      ! If it is the first time buoyant
                      nocape = .false.
                      kLFC   = kLCL
                   end if
                end if
             end if
          end if
       end if
    end if
    
  end subroutine CAPE_below_LCL

  !#######################################################################

  subroutine CAPE_above_LCL (kLCL, kLZB, kLFC, Tp, rp, rin, p_full, nocape_, skip, &
                      CIN, CAPE, Tin, Tin_virtual, p_half, pLZB)
    
    
    integer, intent(in)                    :: kLCL
    real, intent(in), dimension(:)         :: rin, p_full, Tin, Tin_virtual
    real, intent(in), dimension(:)         :: p_half
    logical, intent(in)                    :: nocape_, skip
    integer, intent(inout)                 :: kLZB, kLFC
    real, intent(inout), dimension(:)      :: Tp, rp
    real, intent(inout)                    :: CIN, CAPE, pLZB
    
    integer                                :: k
    real                                   :: a, b, dtdlnp, es
    logical                                :: nocape
    
    nocape = nocape_
    
    ! If the mixing ratio r < 0 then skip
    if (skip) then
       if (nocape) then
          call set_values_if_nocape (Tin, rin, p_full, Tp, rp,&
               pLZB, kLZB, kLFC, CIN )
       end if
    else
       ! If the mixing ratio r>0, do moist adiabatic ascent
       ! Loop over k from LCL to top
       do k = kLCL-1, 1, -1
          a      = kappa*Tp(k+1) + (HLv/Cp_air) * rp(k+1)
          b      = (HLv**2) * rp(k+1)/(Cp_air * rvgas * Tp(k+1)**2)
          dtdlnp = a / (1.0 + b)
          
          Tp(k) = Tp (k+1) + dtdlnp * log(p_full(k)/p_full(k+1))/2
          if ( (Tp(k) .lt. Tmin) .and. nocape) then
             if (nocape) then
                call set_values_if_nocape (Tin, rin, p_full,&
                     Tp, rp, pLZB, kLZB, kLFC, CIN )
             end if
             ! Exit the loop over k
             go to 20
          else
             call escomp(Tp(k), es)
             rp(k) = mixing_ratio(es, ( p_full(k) + p_full(k+1) )/2)
             a = kappa * Tp(k) + (HLv/Cp_air)* rp(k)
             b = (HLv**2)*rp(k) / (Cp_air*rvgas*Tp(k)**2)
             dtdlnp = a/( 1.0 + b )
             Tp(k) = Tp(k+1) + dtdlnp * log( p_full(k)/p_full(k+1) )
             
             if ( (Tp(k) .lt. Tmin) .and. nocape ) then
                if (nocape) then
                   call set_values_if_nocape (Tin, rin, p_full,&
                        Tp, rp, pLZB, kLZB, kLFC, CIN )
                end if
                ! Exit the loop over k
                go to 20
             else
                call escomp(Tp(k), es)
                rp(k) = mixing_ratio(es, p_full(k))
                
                if ( (virtual_temp(Tp(k), rp(k)) .lt. Tin_virtual(k) ) .and. nocape) then
                   ! If the parcel is not buoyant and does not yet have CAPE, add to CIN
                   CIN = CIN + rdgas * ( Tin_virtual(k) - virtual_temp(Tp(k), rp(k)) ) &
                        * log(p_half(k+1)/p_half(k))
                else
                   if ( (virtual_temp(Tp(k), rp(k)) .lt. Tin_virtual(k) ) &
                        .and. (.not.nocape) ) then
                      kLZB = k + 1
                      ! Exit the loop over k
                      go to 20
                   else           
                      ! If the parcel is buoyant, add to CAPE
                      CAPE = CAPE + rdgas * (virtual_temp(Tp(k), rp(k)) - Tin_virtual(k)) &
                           * log(p_half(k+1)/p_half(k))
                      ! State that you have CAPE
                      if (nocape) then
                         nocape = .false.
                         kLFC   = k
                      end if
                   end if
                end if
             end if
          end if
       end do
20  end if
    
  end subroutine CAPE_above_LCL

  !#######################################################################
  
  real function mixing_ratio(vapor_pressure, pressure)
    
    ! calculates the mixing ratio from the vapor pressure and pressure
    
    real, intent(in)     :: vapor_pressure, pressure
    
    mixing_ratio = rdgas * vapor_pressure/rvgas/(pressure-vapor_pressure)
    
  end function mixing_ratio
  
  !#######################################################################
  
  real function virtual_temp(temp, r)
    
    ! Calculates the virtual temperature from the temperature and mixing ratio 
    ! consistent with the approximation used in the fms code
    
    real, intent(in)     :: temp         ! temperature
    real, intent(in)     :: r            ! mixing ratio 
    
    real                 :: q            ! specific humidity
    
    q            = r / (1.0 + r)
    virtual_temp = temp * (1.0 + q * (rvgas/rdgas-1.0))
    
  end function virtual_temp
  
  !#######################################################################
  
  subroutine Pq_calculation (kLZB, k_surface, qref, qin, p_half, deltaq, Pq, dt)
    
    integer, intent(in)                 :: kLZB, K_surface
    real, intent(in), dimension(:)      :: qref, qin
    real, intent(in), dimension(:)      :: p_half
    real, intent(in)                    :: dt
    real, intent(inout), dimension(:)   :: deltaq
    real, intent(out)                   :: Pq
    
    integer                             :: k
    
    ! Initialization
    Pq = 0.
    
    ! Calculation of the delta q  and the precipitation
    do k = kLZB, k_surface
       deltaq(k) = - (qin(k) - qref(k)) * dt/tau_bm
       Pq        = Pq + deltaq(k) * (p_half(k)-p_half(k+1))
    end do
    
    Pq = Pq/grav
    
  end subroutine Pq_calculation
  
  !#########################################################################
  
  subroutine Pt_calculation (kLZB, k_surface, Tref, Tin, p_half, deltaT, Pt,dt)
    
    integer, intent(in)                    :: kLZB, k_surface
    real, intent(in), dimension(:)         :: Tref, Tin
    real, intent(in), dimension (:)        :: p_half
    real, intent(in)                       :: dt
    real, intent(inout), dimension(:)      :: deltaT
    real, intent(out)                      :: Pt
    
    integer                                :: k
        
    ! Initialization
    Pt = 0.
    
    ! Calculation of delta T and precipitation
    do k=kLZB, k_surface
       deltaT(k) = -(Tin(k) - Tref(k)) * dt/tau_bm
       Pt        = Pt + (Cp_air/(HLv + small)) * deltaT(k) &
            * (p_half(k+1) - p_half(k))
    end do
    
    Pt = Pt/grav
    
  end subroutine Pt_calculation
  
  !#######################################################################
  
  subroutine set_reference_profiles(p_full, qin, Tin, Tp, kLZB, k_surface, &
       rp, deltaq, deltaT, qref, Tref)
    
    integer, intent(in)                    :: kLZB, k_surface
    real, intent(in), dimension(:)         :: p_full
    real, intent(in), dimension(:)         :: qin
    real, intent(in), dimension(:)         :: Tin, Tp
    real, intent(inout), dimension(:)      :: rp, deltaq, deltaT
    real, intent(out), dimension (:)       :: qref, Tref
    
    integer                                :: k
    real                                   :: eref
    
    ! Initialization
    Tref = Tp
    
    ! Under the LZB
    do k = kLZB, k_surface
       eref       = rhbm * p_full(k) * rp(k) / (rp(k) + (rdgas/rvgas))
       rp(k) = mixing_ratio(eref, p_full(k))
       qref(k)    = rp(k) / (1 + rp(k))
    end do
    
    ! Above the LZB
    k = max(kLZB-1, 1)
    call set_profiles_to_full_model_values (1, k, Tin, qin, Tref, deltaT, qref, &
                        deltaq)
    
  end subroutine set_reference_profiles
  
  !###################################################################################
    
  subroutine do_shallow_convection(kLZB, k_surface, qin, qref, deltaq, Tin, &
       Tref, deltaT, p_half, Pq,dt)
    
    integer, intent(in)                    :: kLZB, k_surface
    real, intent(in), dimension(:)         :: qin
    real, intent(in), dimension(:)         :: Tin
    real, intent(in), dimension(:)         :: p_half
    real, intent(in)                       :: dt
    real, intent(inout), dimension(:)      :: qref, deltaq
    real, intent(inout), dimension(:)      :: Tref, deltaT
    real, intent(inout)                    :: Pq
    
    integer                                :: k_top
    logical                                :: k_zero_precip_found
    
    ! Search for a lower level with Pq > 0
    call level_of_zero_precip(kLZB, k_surface, deltaq, p_half, Pq,  &
         k_zero_precip_found, k_top, deltaT, Tref, qref, Tin, qin)
    
    ! If this level exists   
    if (k_zero_precip_found) then
       ! Change the reference temperature and the LZB
       call change_Tref_LZB_shallowconv(Pq, k_top,k_surface, deltaq, &
            Tref, deltaT, p_half,dt)
       
    else
       ! Else, if this level doesn't exist (k_top = k_surface and P<0),   &
       ! modify reference profiles and do nothing
       if(k_top == kLZB) then
          call set_profiles_to_full_model_values (k_surface, k_surface, Tin, &
               qin, Tref, deltaT, qref, deltaq)

       else
          call set_profiles_to_full_model_values (kLZB, k_top, Tin, qin, Tref, &
               deltaT, qref, deltaq)
       end if
    end if
    
    ! Set the precipitation rate to zero
    Pq = 0.
    
  end subroutine do_shallow_convection
  
  !##########################################################################################
  
  subroutine level_of_zero_precip (kLZB, k_surface, deltaq, p_half, Pq,    &
       k_zero_precip_found, k_zero_precip, deltaT, Tref, qref, Tin, qin)
    
    integer, intent(in)                 :: kLZB, k_surface
    real, intent(in), dimension(:)      :: p_half
    real, intent(in), dimension(:)      :: Tin, qin
    real, intent(inout), dimension(:)   :: deltaq, deltaT, Tref, qref
    real, intent(inout)                 :: Pq
    logical, intent(out)                :: k_zero_precip_found
    integer, intent(out)                :: k_zero_precip
    
    integer                             :: k
    
    ! Current level k
    k    = kLZB
    
    ! Initialization of k_zero_precip_found; by default,
    ! the level of zero precipitation does not exist
    k_zero_precip_found = .false.
    
    ! Calculation of the precipitation up to one level below, until P > 0 
    ! or surface reached 
    do while ( (Pq .lt. 0.) .and. (k .le. k_surface) )
       Pq = Pq - deltaq(k) * (p_half(k) - p_half(k+1))/grav
       k  = k + 1
    end do
    
    ! The level of zero precipitation (if it exists) is 
    ! the level before the while condition is false
    k_zero_precip = k - 1
    
    ! If the level of zero precipitation exists, returns True
    if (Pq .gt. 0.) then
       k_zero_precip_found = .true.
    end if
    
    ! Above k_zero_precip, put original temperature and humidity
    if (k_zero_precip .gt. kLZB) then
       call set_profiles_to_full_model_values (kLZB,k_zero_precip-1 , Tin, qin, &
            Tref, deltaT, qref, deltaq)
    end if
      
  end subroutine level_of_zero_precip
  
  !#########################################################################################


  subroutine change_Tref_LZB_shallowconv(Pq, k_top, k_surface, deltaq, Tref, &
       deltaT, p_half, dt)
    
    real, intent(in)                       :: Pq, dt
    integer, intent(in)                    :: k_top, k_surface
    real, intent(in), dimension(:)         :: p_half
    real, intent(inout), dimension(:)      :: deltaq
    real, intent(inout), dimension(:)      :: Tref, deltaT
    
    integer                                :: k
    real                                   :: c, deltak
    
    ! Below the LZB, put the new Tref and humidity
    
    ! First calculate the coefficient to apply to deltaq(kLZB)
    ! so the precipitation is identically zero
    c = Pq * grav / (deltaq(k_top) * (p_half(k_top+1) - p_half(k_top)))
    
    ! Modify the last fraction of deltaq
    deltaq(k_top) = deltaq(k_top)*c
    
    ! Modify deltaT(kLZB) used in the calculation of delta k
    deltaT(k_top) = deltaT(k_top)*c
    
    ! Calculation of deltak
    deltak = 0.
    do k = k_top,k_surface
       deltak = deltak + deltaT(k) * (p_half(k) - p_half(k+1))
    end do
    deltak = deltak / (p_half(k_surface+1) - p_half(k_top))
    
    ! Modify deltaT and Tref
    if (k_top /= k_surface) then
       deltaT(k_top:k_surface) = deltaT(k_top:k_surface) + deltak
       Tref(k_top:k_surface)   = Tref(k_top:k_surface) + deltak*tau_bm/dt
    end if
    
  end subroutine change_Tref_LZB_shallowconv
  
  !############################################################################################

  subroutine do_deep_convection (kLZB, k_surface,Pt, dt,p_half,invtau_q_relaxation, &
                                 invtau_t_relaxation,Pq,deltaT,Tref,deltaq)

    integer, intent(in)                    :: kLZB, k_surface
    real, intent(in)                       :: Pt, dt
    real, intent(in), dimension(:)         :: p_half
    real, intent(inout)                    :: Pq
    real, intent(inout), dimension(:)      :: deltaT
    real, intent(inout), dimension(:)      :: Tref
    real, intent(inout), dimension(:)      :: deltaq
    real, intent(out)                      :: invtau_q_relaxation, invtau_t_relaxation

    if (Pq.gt.Pt)then
       ! Do deep convection by changing time scales
       call do_change_time_scale_deepconv(kLZB, k_surface, Pt, Pq, deltaq,&
                                  invtau_q_relaxation, invtau_t_relaxation)
    else
       ! Do deep convection by changing the reference temperature profile
       call do_change_Tref_deepconv(kLZB, k_surface, deltaT, deltaq, p_half, Tref,dt)
    end if



  end subroutine do_deep_convection

  !############################################################################################
  
  subroutine do_change_Tref_deepconv(kLZB, k_surface, deltaT, deltaq, p_half, Tref,dt)
    
    integer, intent(in)                    :: kLZB, k_surface
    real, intent(in), dimension(:)         :: deltaq
    real, intent(in), dimension(:)         :: p_half
    real, intent(in)                       :: dt
    real, intent(inout), dimension(:)      :: deltaT
    real, intent(inout), dimension(:)      :: Tref
    
    
    integer                                :: k
    real                                   :: deltak
      
    ! Calculation of deltak: shift of temperature profile necessary 
    ! to conserve enthalpy
    deltak = 0.
    do k=kLZB, k_surface
       deltak = deltak    &
            - (deltaT(k) + (HLv/Cp_air)*deltaq(k)) * (p_half(k+1) - p_half(k))
    end do
    
    ! divide by pressure difference over convective layer
    deltak = deltak / (p_half(k_surface+1) - p_half(kLZB))
    
    ! Modification of the reference temperature profile in convective layer 
    ! (below LZB) by uniform shift to conserve enthalpy
    do k = kLZB, k_surface
       Tref(k)   = Tref(k)   + deltak*tau_bm/dt
       deltaT(k) = deltaT(k) + deltak
    end do
    
  end subroutine do_change_Tref_deepconv

  !##############################################################################

  subroutine do_change_time_scale_deepconv (kLZB, k_surface, Pt, Pq, deltaq,&
                                  invtau_q_relaxation, invtau_t_relaxation)

    integer, intent(in)                     :: kLZB, k_surface
    real, intent(in)                        :: Pt
    real, intent(inout)                     :: Pq
    real, intent(inout), dimension (:)      :: deltaq
    real, intent(out)                       :: invtau_q_relaxation, invtau_t_relaxation

    invtau_q_relaxation = Pt/Pq/tau_bm
    deltaq(kLZB:k_surface) = tau_bm*invtau_q_relaxation*deltaq(kLZB:k_surface)
    Pq = Pt
    invtau_t_relaxation = 1./tau_bm


  end subroutine do_change_time_scale_deepconv

  !##############################################################################

  subroutine set_values_if_nocape (Tin, rin, p_full, Tp, rp, pLZB, kLZB, kLFC, CIN )
    
    real, intent(in), dimension(:)       :: Tin, rin, p_full
    real, intent(inout), dimension(:)      :: Tp, rp
    real, intent(inout)                    :: pLZB
    integer, intent(inout)                 :: kLZB, kLFC
    real, intent(inout)                    :: CIN

    pLZB = p_full(1)
    kLZB = 0
    kLFC = 0
    CIN = 0.
    Tp = Tin
    rp = rin


  end subroutine set_values_if_nocape

  !##############################################################################

  subroutine set_profiles_to_full_model_values (k_1, k_2, Tin, qin, Tref, deltaT, qref, deltaq)

    integer, intent(in)                    :: k_1, k_2
    real, intent(in), dimension(:)         :: Tin, qin                       
    real, intent(inout), dimension(:)      :: Tref, deltaT, qref, deltaq


    ! Set profiles to full model values between levels k_1 and k_2

    Tref   (k_1:k_2) = Tin (k_1:k_2)
    qref   (k_1:k_2) = qin (k_1:k_2)
    deltaT (k_1:k_2) = 0.
    deltaq (k_1:k_2) = 0.


  end subroutine set_profiles_to_full_model_values

  !##############################################################################
  
  subroutine get_lcl_temp(lcl_temp_table, value, val_min, val_max, Tlcl)
    
    ! Returns the approximation of the LCL temperature matching the value
    ! using the lookup table lcl_temp_table
    
    real, intent(in), dimension(:)    :: lcl_temp_table
    real, intent(in)                  :: value, val_min, val_max
    real, intent(out)                 :: Tlcl
    
    real                              :: iv_floor, w_floor, w_ceil
    
    if (value .lt. val_min) then
       write(*,*) 'qe_moist_convection: Value passed to get_lcl_temp: ',value
       call error_mesg ('qe_moist_convection', 'get_lcl_temp: value to low.', FATAL)
    end if
    
    if (value .gt. val_max) then
       write(*,*) 'qe_moist_convection: Value passed to get_lcl_temp: ', value
       call error_mesg ('qe_moist_convection', 'get_lcl_temp: value too high.', FATAL)
    end if
        
    ! Find index of nearest lookup table value below
    iv_floor = floor( (value-val_min)/val_inc ) + 1
    
    ! Linearly interpolate to get the temperature of LCL
    w_floor  = (val_min + (iv_floor-1)*val_inc)
    w_ceil   = (value -w_floor) / val_inc
    Tlcl     = lcl_temp_table(iv_floor+1)*w_ceil - lcl_temp_table(iv_floor)*(w_ceil-1)
      
  end subroutine get_lcl_temp

  !##############################################################################
    
  real function lcl_temp(value, lcl_temp_guess)
    
    real, intent(in)      :: value, lcl_temp_guess
    real, parameter       :: precision = 1.e-7
    integer, parameter    :: max_iter = 100
    real                  :: T, dT
    integer               :: iter 
    
    ! This routine calculates the temperature at the LCL by Newton
    ! iteration. 
    ! 
    ! It solves the nonlinear equation 
    !
    !     value = log[es/T**(1/kappa)]                      (1)
    !
    ! for the temperature. The rationale is as follows:
    !
    ! The potential temperature at the LCL is  the same as at the
    ! parcel origin (dry adiabatic ascent). So theta(parcel) =
    ! theta(LCL). Since mixing ratio r is also conserved in
    ! unsaturated adiabatic ascent, at the LCL we have r(parcel) = r_sat(LCL). So  
    !
    ! T(LCL) = theta(LCL) (p(LCL)/pref)**kappa = theta(parcel) (p(LCL)/pref)**kappa, 
    !
    ! and, at the LCL, 
    !
    ! r(parcel) = r_sat(LCL) =  
    !
    !  rdgas/rvgas * es[T(LCL)]/[p(LCL) - es[T(LCL)]]] 
    !      = rdgas/rvgas *es / [(T(LCL)/theta(parcel))**(1/kappa)*pref - es(LCL)]
    !
    ! Therefore, 
    !
    !            theta**(-1/kappa) * pref * r/(rdgas/rvgas + r) = es/T**(1/kappa).  (2)
    !
    ! On the LHS are parcel properties; on the RHS are properties of the LCL, which are 
    ! a function of temperature only. Given the LHS, we can solve for the temperature at 
    ! the LCL, which is what this function does.
    !
    ! The input value of this function (value) is the logarithm of the LHS of (2), 
    !
    ! value_ = log[theta**(-1/kappa) * pref * r/(rdgas/rvgas + r)]. 
    !
    ! The LCL temperature is the solution of (1), where the
    ! saturation vapor pressure es is a function of temperature. Given
    ! value_, the temperature at which this equation is satisfied is
    ! computed by Newton iteration, up to a given precision.

    ! Initialization
    T    = lcl_temp_guess
    dT   = precision + 1.
    iter = 0
   
    ! Newton iterations to find temperature of LCL
    do while ( (abs(dT) .gt. precision) .and. (iter .lt. max_iter) )
       dT   = lcl_value_difference(T, value) / dlcl_value_difference(T)
       T    = T - dT
       iter = iter + 1
    end do
    
    if (dT .lt. precision) then
       lcl_temp = T
    else
       write(*,*) 'qe_moist_convection: LCL calculation did not converge. Precision not achieved.'
       call error_mesg ('qe_moist_convection', 'lcl_temp', FATAL) 
    end if
    
  end function lcl_temp
  
  !##############################################################################
  
  real function lcl_value_difference (T, value)
      
    real, intent(in)      :: T, value   
    real                  :: es
    
    call escomp(T, es)
    lcl_value_difference = value - log(es * T**(-1/kappa))
    
  end function lcl_value_difference

  !##############################################################################
  
  real function dlcl_value_difference (T)
    
    real, intent(in)       :: T
    
    ! Derivative of the function lcl_value_difference
    dlcl_value_difference = 1/kappa * T**(-1) - HLv/rvgas * T**(-2)
    
  end function dlcl_value_difference
  
  !##############################################################################


  subroutine qe_moist_convection_end
    
    deallocate (lcl_temp_table)
    
  end subroutine qe_moist_convection_end
  
  !##############################################################################
  
end module qe_moist_convection_mod

