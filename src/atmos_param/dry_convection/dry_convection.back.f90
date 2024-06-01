!> Simple Dry Convection Scheme
!!
!! A simple dry convective adjustment parameterisation.
!! Adapted from the dry GCM of Schneider and Walker (2006).
!!
!! Documented and added to the moist-model stack by James Penn.
!!
!! Original code source from Tapio Schneider's FMS repository.
!!
!! ### Namelist parameters
!! @param tau The relaxation timescale to the lapse rate gamma.
!! @param gamma The prescribed lapse rate.
!!              When gamma = 1, parcel lifting temperature is dry potential temperature.
!!
!! @see https://github.com/tapios/fms-idealized
!!
!! @author Tapio Schneider
!! @author James Penn
module dry_convection_mod

  use           fms_mod, only: error_mesg, FATAL, stdlog, mpp_pe,             &
                               mpp_root_pe, open_namelist_file, close_file,   &
                               check_nml_error
  use     constants_mod, only: rdgas, cp_air
  use  time_manager_mod, only: time_type
  use  diag_manager_mod, only: register_diag_field, send_data

!-----------------------------------------------------------------------------

  implicit none
  private

!                             ---  namelist ---
  real :: tau, &            !< relaxation timescale [seconds]
          gamma             !< prescibed lapse rate [non-dim]

  namelist /dry_convection_nml/ tau, gamma

  integer :: i, j, jit, k, num_levels
  real :: cons1 = 0.  !< Potential temperature exponent (R/Cp)

  integer :: id_cape, id_cin, id_lzb, id_lcl, id_tp, id_n_tp, &
       id_dp, id_amb, id_dt

  character(len=14), parameter :: mod_name='dry_convection'
  real :: missing_value = -1.e-10

  public :: dry_convection_init, dry_convection

  contains

!> Initialise the dry convection scheme
!!
!! @param[in] axes The dimensions of the model
!! @param[in] Time Initial time (as a time_type)
    subroutine dry_convection_init(axes, Time)
      integer, intent(in) :: axes(4)
      type(time_type), intent(in) :: Time

      integer :: unit, ierr, io

      unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
         read(unit, nml=dry_convection_nml, iostat=io, end=20)
         ierr = check_nml_error (io, 'dry_convection_nml')
      enddo

20    call close_file (unit)

      if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=dry_convection_nml)

      !s initialise here as rdgas no longer a paramter
      cons1 = rdgas/cp_air  !< Potential temperature exponent (R/Cp)

      id_cape = register_diag_field ( mod_name, 'CAPE', axes(1:2), Time, &
           'CAPE', 'J/kg', missing_value=missing_value)

      id_cin = register_diag_field ( mod_name, 'CIN', axes(1:2), Time, &
           'CIN', 'J/kg', missing_value=missing_value)

      id_dp = register_diag_field ( mod_name, 'dp', axes(1:2), Time, &
           'Pressure interval', 'Pa', missing_value=missing_value)

      id_lzb = register_diag_field ( mod_name, 'LZB', axes(1:2), Time, &
           'Level of zero buoyancy', 'index', missing_value=missing_value)

      id_lcl = register_diag_field ( mod_name, 'LCL', axes(1:2), Time, &
           'Lifting condensation level', 'index', missing_value=missing_value)

      id_tp = register_diag_field ( mod_name, 'parcel_temp', axes(1:3), &
           Time, 'Relaxation temperature', 'K', missing_value=missing_value)

      id_n_tp = register_diag_field ( mod_name, 'nonadj_parcel_temp', &
           axes(1:3), Time, 'Relaxation temperature before adjustment', 'K', &
           missing_value=missing_value)

      id_amb = register_diag_field ( mod_name, 'ambient_temp', axes(1:3), &
           Time, 'Ambient temperature', 'K', missing_value=missing_value)

      id_dt = register_diag_field ( mod_name, 'dt_tg', axes(1:3), &
           Time, 'Temperature tendency', 'K/s', missing_value=missing_value)

    end subroutine dry_convection_init

!> Run a single timestep of the dry convection scheme
!!
!! Taking the current temperature profile of the atmosphere,
!! dry_convection relaxes towards a prescribed lapse rate gamma.
!! This method returns a temperature tendency dt_tg for all points
!! in the model domain.
!!
!! @param[in] Time Current time
!! @param[in] tg   Temperature array
!! @param[in] p_full Pressure values on full levels
!! @param[in] p_half Pressure values on half levels
!! @param[out] dt_tg Calculated temperature tendency
!! @param[out] cape Convective Available Potential Energy
!! @param[out] cin Convective Inhibition
  subroutine dry_convection(Time, tg, p_full, p_half, dt_tg, cape, cin, lzb, lcl)

    type(time_type), intent(in) :: Time

    real, intent(in), dimension(:,:,:) ::                                     &
         tg,               &   ! temperature
         p_full,           &   ! pressure on full model levels
         p_half                ! pressure on half model levels

    real, intent(out), dimension(:,:,:) ::                                    &
         dt_tg                 ! temperature tendency
    real, intent(out), dimension(size(tg,1),size(tg,2)) ::                    &
        cape,              &   !< convectively available potential energy
        cin                    !< convective inhibition
    integer, intent(out), dimension(:,:)   ::                                 &
        lcl,               &   !< lifting condensation level (index)
        lzb                    !< level of zero buoyancy


!                         ---  local variables ---
    real, dimension(size(tg,1),size(tg,2)) ::                                 &
         dp,               &   !< pressure interval from ground to LZB
         ener_int              !< energy integral from ground to LZB

    integer, dimension(size(tg,1),size(tg,2)) ::                              &
         btm                   !< bottom of convecting region

    real, dimension(size(tg,1),size(tg,2), size(tg,3)) ::                     &
         tp,               &   !< parcel lifting temperature
         n_tp,             &   !< parcel lifting temperature before adjustment
         dp_half               !< spacing between half pressure levels

    logical :: used

    num_levels = size(tg,3)

    ! half-level spacings:
    do k=1, num_levels
       dp_half(:,:,k) = p_half(:,:,k+1) - p_half(:,:,k)
    end do

    ! set the lower bound of convecting region (always the lowest level)
    btm = num_levels

    ! calculate convection quantities
    call capecalc(                tg,                          p_full,        &
                              p_half,                         dp_half,        &
                                 btm,                             lzb,        &
                                 lcl,                            cape,        &
                                 cin,                              tp )

    n_tp = tp ! save uncorrected parcel profile

    ! ensure conservation of energy
    ener_int(:,:) = 0.0
    dp(:,:) = 0.0

    do i=1, size(tg,1)
       do j=1, size(tg,2)
         
          do k=1, num_levels
             if(k>=lzb(i,j).and. k<=btm(i,j)) then
                ! in convecting region between ground and LZB
                ener_int(i,j) = ener_int(i,j) +                               &
                     dp_half(i,j,k) * (tg(i,j,k) - tp(i,j,k))
                dp(i,j) = dp(i,j) + dp_half(i,j,k)
             else
                ! ambient in non-convecting regions
                tp(i,j,k) = tg(i,j,k)
             endif
          end do

          ! normalize integral
          ener_int(i,j) = ener_int(i,j) / dp(i,j)

          ! At each level in convecting region, add integral energy
          do k=btm(i,j), lzb(i,j), -1
             tp(i,j,k) = tp(i,j,k) + ener_int(i,j)
          end do
       end do
    end do

    ! temperature tendency
    dt_tg = (tp - tg) / tau


    if (id_cin  > 0) used = send_data ( id_cin,  cin,        Time, 1, 1)
    if (id_cape > 0) used = send_data ( id_cape, cape,       Time, 1, 1)
    if (id_lzb  > 0) used = send_data ( id_lzb,  float(lzb), Time, 1, 1)
    if (id_lcl  > 0) used = send_data ( id_lcl,  float(lcl), Time, 1, 1)
    if (id_tp   > 0) used = send_data ( id_tp,   tp,         Time, 1, 1)
    if (id_n_tp > 0) used = send_data ( id_n_tp, n_tp,       Time, 1, 1)
    if (id_dp   > 0) used = send_data ( id_dp,   dp,         Time, 1, 1)
    if (id_amb  > 0) used = send_data ( id_amb,  tg,         Time, 1, 1)
    if (id_dt   > 0) used = send_data ( id_dt,   dt_tg,      Time, 1, 1)

  end subroutine dry_convection

  !> Calculate convective available potential energy (CAPE)
  !!
  !! Given temperature and pressure profile, calculates the
  !! heights of the and LZB, LVL (in model levels).
  !!
  !! Also calculates CAPE and CIN at each lat-lon point
  !! and the parcel lifting temperature at all lat-lon-p points.
  !!
  !! @param[in] tg      Temperature array
  !! @param[in] p_full  Pressure values on full levels
  !! @param[in] p_half  Pressure values on half levels
  !! @param[in] dp_half Spacing of half pressure levels
  !! @param[in] btm     Level of 'bottom' of convective region
  !!
  !! @param[out] lzb    Level of Zero Buoyancy (top of convection)
  !! @param[out] lcl    Lifting Condensation Level (btm of convection)
  !! @param[out] cape   Convectively Available Potential Energy
  !! @param[out] cin    Convective INhibition
  !! @param[out] tp     Parcel lifting temperature
  subroutine capecalc(tg, p_full, p_half, dp_half, btm, lzb, lcl, cape,       &
       cin, tp)

    integer, intent(in), dimension(:,:) ::                                    &
         btm                   ! level parcel is lifted from

    real, intent(in), dimension(:,:,:) ::                                     &
         tg,               &   ! gridpoint temperature
         p_full,           &   ! gridpoint pressure on full model levels
         p_half,           &   ! gridpoint pressure on half model levels
         dp_half               ! spacing between half pressure levels

    real, intent(out), dimension(:,:) ::                                      &
         cape,             &   ! convectively available potential energy
         cin                   ! convective inhibition

    integer, intent(out), dimension(:,:) ::                                                &
         lzb,              &   ! level of zero buoyancy (top of convection)
         lcl                   ! lifting condensation level (btm of convection)

    real, intent(out), dimension(:,:,:) ::                                    &
         tp                    ! parcel lifting temperature

    real :: zdpkpk             ! pressure spacing

    cape(:,:) = 0.0; cin(:,:) = 0.0
    lzb(:,:) = btm; lcl(:,:) = btm

    tp(:,:,:) = tg(:,:,:)

    do i = 1, size(tg,1)
       do j = 1, size(tg,2)

          ! lift parcel with lapse rate given by gamma
          do k = btm(i,j)-1, 1, -1
             zdpkpk = exp( cons1 * alog(p_full(i,j,k)/p_full(i,j,k+1)))  !< equiv. to (pfull[k]/pfull[k+1])^cons1
             tp(i,j,k) = tp(i,j,k+1) +                                        &
                  gamma * (tp(i,j,k+1)*zdpkpk - tp(i,j,k+1))
          end do

          ! find LCL
          do k = btm(i,j)-1, 1, -1
             if(tp(i,j,k) > tg(i,j,k)) then ! unstable parcel
                if(lzb(i,j) == btm(i,j)) then ! not above a lower cloud
                   ! calculate CAPE
                   cape(i,j) = cape(i,j) +                                    &
                        rdgas*(tp(i,j,k)-tg(i,j,k))*log(p_half(i,j,k+1)/p_half(i,j,k))

                   ! set LCL if parcel is stable at next lower level
                   if(tp(i,j,k+1) < tg(i,j,k+1)) lcl(i,j) = k

                   ! set LZB if parcel is stable at next higher level
                   ! (but limit to the top of the model!)
                   if( (tp(i,j,k-1) < tg(i,j,k-1)) .or. (k == 1)) lzb(i,j) = k
                else ! above cloud level, set parcel temperature to ambient
                   tp(i,j,k) = tg(i,j,k)
                end if
             end if

             if(tp(i,j,k) <= tg(i,j,k)) then ! stable parcel
                if((lzb(i,j) == btm(i,j))) then ! not above a lower cloud
                   if(lcl(i,j) == btm(i,j)) then
                      ! calculate CIN (only if below LCL)
                      cin(i,j) = cin(i,j) -                                   &
                           rdgas*(tp(i,j,k)-tg(i,j,k))*log(p_half(i,j,k+1)/p_half(i,j,k))
                   end if

                else ! set parcel temp to ambient
                   tp(i,j,k) = tg(i,j,k)
                end if
             end if
          end do

          ! if cin > cape, turn off convection by setting
          ! relaxation temperature to ambient temperature
          if(cin(i,j) > cape(i,j)) tp(i,j,:) = tg(i,j,:)

          ! a few checks
          if( (lcl(i,j) /= btm(i,j)) .and. (lzb(i,j) == btm(i,j)))        &
               call error_mesg ('dry_convection','LCL defined, LZB not defined', FATAL)
          if( lcl(i,j) < lzb(i,j) )                                           &
               call error_mesg ('dry_convection','LCL above LZB', FATAL)

          ! if LCL and LZB have not been raised from the ground by the loops above
          ! CAPE and CIN must both be zero.
          if( (lcl(i,j) == btm(i,j)) .and. (lzb(i,j) == btm(i,j))) then
             cape(i,j) = 0.0
             cin(i,j) = 0.0
          end if

       end do
    end do

  end subroutine capecalc

end module dry_convection_mod