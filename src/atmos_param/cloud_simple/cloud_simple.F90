module cloud_simple_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  use            fms_mod, only: stdlog, FATAL, WARNING, NOTE, error_mesg, uppercase
  use   time_manager_mod, only: time_type
  use sat_vapor_pres_mod, only:  compute_qs
  use diag_manager_mod, only: register_diag_field, send_data

  implicit none

  logical ::   do_init = .true.  ! Check if init needs to be run

  real :: zerodegc = 273.15
  integer :: id_cf, id_reff_rad, id_frac_liq, id_qcl_rad, id_rh_in_cf, id_simple_rhcrit, &
             id_tot_cld_amt, id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt
  character(len=14), parameter ::   mod_name_cld = "cloud_simple"

  integer, parameter :: B_SPOOKIE = 1, B_SUNDQVIST = 2, B_SMITH = 3
  integer, private :: cf_diag_formula = B_SPOOKIE
  character(len=32) :: cf_diag_formula_name = 'spookie'
  ! character(len=32) :: diag_cf_method = 'simple_rhcrit' ! 'read_rhcrit', 'linear', 'quadratic'
  logical :: do_simple_rhcrit = .true.
  logical :: do_read_scm_rhcrit = .false.
  logical :: do_linear_diag = .false.
  logical :: do_qcl_with_temp = .false.
  logical :: do_qcl_two_paras = .false.
  logical :: do_cloud_amount_diags = .true.
  logical :: adjust_top = .true.
  logical :: do_add_stratocumulus = .false.
  real, parameter :: FILL_VALUE = -999 ! Fill value for arrays
  real, dimension(100) :: scm_rhcrit = FILL_VALUE   ! Input array for single column critical RH. Max number of levels = 100
  real, dimension(100,2) :: scm_linear_coeffs = FILL_VALUE

  real    ::   simple_cca =  0.0
  real    ::   rhcsfc     = 0.95
  real    ::   rhc700     = 0.7
  real    ::   rhc200     = 0.3
  real    ::   cf_min     = 1e-10

  namelist /cloud_simple_nml/ simple_cca, rhcsfc, rhc700, rhc200, &
                              cf_diag_formula_name, do_simple_rhcrit, &
                              do_linear_diag, scm_linear_coeffs, &
                              do_read_scm_rhcrit, scm_rhcrit, &
                              do_qcl_with_temp, do_qcl_two_paras, &
                              do_cloud_amount_diags, adjust_top, &
                              do_add_stratocumulus

  contains

  !-----------------------------------------------

  subroutine cloud_simple_init (axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    integer :: io ,stdlog_unit

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=cloud_simple_nml, iostat=io)
#else
    if ( file_exist('input.nml') ) then
      nml_unit = open_namelist_file()
      read (nml_unit, cloud_simple_nml, iostat=io)
      call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, cloud_simple_nml)

    ! logical
    if(do_simple_rhcrit) then
      if(do_read_scm_rhcrit) call error_mesg( 'cloud_simple_init', &
        'SETTING DO_READ_SCM_RHCRIT TO FALSE AS DO_SIMPLE_RHCRIT IS .TRUE.', NOTE)
      do_read_scm_rhcrit = .false.
      if(do_linear_diag) call error_mesg( 'cloud_simple_init', &
        'SETTING DO_LINEAR_DIAG TO FALSE AS DO_SIMPLE_RHCRIT IS .TRUE.', NOTE)
      do_linear_diag = .false.
    else if(do_read_scm_rhcrit) then
      if(do_linear_diag) call error_mesg( 'cloud_simple_init', &
        'SETTING DO_LINEAR_DIAG TO FALSE AS DO_READ_SCM_RHCRIT IS .TRUE.', NOTE)
      do_linear_diag = .false.
    else
      if(.NOT. do_linear_diag) call error_mesg( 'cloud_simple_init', &
        'DO_LINEAR_DIAG SHOULD BE .TRUE. AS BOTH DO_SIMPLE_RHCRIT AND DO_READ_SCM_RHCRIT ARE .FALSE.', FATAL)
    endif

    ! Select cloud fraction diag formula
    if(uppercase(trim(cf_diag_formula_name)) == 'SPOOKIE') then
      cf_diag_formula = B_SPOOKIE
      call error_mesg('cloud_simple', 'Using default SPOOKIE cloud fraction diagnostic formula.', NOTE)
    else if(uppercase(trim(cf_diag_formula_name)) == 'SUNDQVIST') then
      cf_diag_formula = B_SUNDQVIST
      call error_mesg('cloud_simple', 'Using Sundqvist (1987) cloud fraction diagnostic formula.', NOTE)
    else if(uppercase(trim(cf_diag_formula_name)) == 'SMITH') then
      cf_diag_formula = B_SMITH
      call error_mesg('cloud_simple', 'Using Smith (1990) cloud fraction diagnostic formula.', NOTE)
    else
      call error_mesg('cloud_simple', '"'//trim(cf_diag_formula_name)//'"'//' is not a valid cloud fraction diagnostic formula.', FATAL)
    endif

    !register diagnostics
    id_cf = &
      register_diag_field ( mod_name_cld, 'cf', axes(1:3), Time, &
      'Cloud fraction for the simple cloud scheme', 'unitless: values 0-1')

    if (do_cloud_amount_diags) then
      id_tot_cld_amt = &
        register_diag_field ( mod_name_cld, 'tot_cld_amt', axes(1:2), Time, &
                            'total cloud amount', 'percent' )

      id_high_cld_amt = &
        register_diag_field ( mod_name_cld, 'high_cld_amt', axes(1:2), Time, &
                            'high cloud amount', 'percent' )

      id_mid_cld_amt = &
        register_diag_field ( mod_name_cld, 'mid_cld_amt', axes(1:2), Time, &
                            'mid cloud amount', 'percent' )

      id_low_cld_amt = &
        register_diag_field ( mod_name_cld, 'low_cld_amt', axes(1:2), Time, &
                            'low cloud amount', 'percent' )
    end if

    id_frac_liq = &
      register_diag_field ( mod_name_cld, 'frac_liq', axes(1:3), Time, &
      'Liquid cloud fraction (liquid, mixed-ice phase, ice)', &
      'unitless: values 0-1')

    id_reff_rad = &
      register_diag_field ( mod_name_cld, 'reff_rad', axes(1:3), Time, &
      'Effective cloud particle radius', &
      'microns')
 
    id_qcl_rad = &
      register_diag_field ( mod_name_cld, 'qcl_rad', axes(1:3), Time, &
      'Specific humidity of cloud liquid', &
      'kg/kg')

    id_rh_in_cf = &
      register_diag_field ( mod_name_cld, 'rh_in_cf', axes(1:3), Time, &
      'RH as a percent', &
      '%')

    if(do_simple_rhcrit .OR. do_read_scm_rhcrit) then
      id_simple_rhcrit = &
        register_diag_field ( mod_name_cld, 'simple_rhcrit', axes(1:3), Time, &
        'RH as a percent', &
        '%')
    end if

    do_init = .false.  !initialisation completed
  end subroutine cloud_simple_init


  subroutine cloud_simple(p_half, p_full, Time,   &
                      temp, q_hum,                &
                      ! outs
                      cf, reff_rad, qcl_rad)
    real       , intent(in), dimension(:,:,:)  :: temp, q_hum, p_full, p_half
    type(time_type) , intent(in)               :: Time
    real       , intent(out), dimension(:,:,:) :: cf, reff_rad, qcl_rad
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: qs, frac_liq, rh_in_cf, simple_rhcrit, theta
    real, dimension(size(temp,1), size(temp,2)) :: tca, high_ca, mid_ca, low_ca
    integer, dimension(size(temp,1), size(temp,2)) :: kdthdp
    integer :: i, j, k, k_surf, k700
    logical :: es_over_liq_and_ice
    real :: strat

    !check initiation has been done - ie read in parameters
    if (do_init) call error_mesg ('cloud_simple',  &
         'cloud_simple_init has not been called.', FATAL)

    if(do_read_scm_rhcrit) then
      if(scm_rhcrit(size(temp,3)).eq.FILL_VALUE) then
        call error_mesg('cloud_simple', 'Input rhcrit must be specified on model pressure levels but not enough levels specified', FATAL)
      endif
      if(scm_rhcrit(size(temp,3)+1).ne.FILL_VALUE) then
        call error_mesg('cloud_simple', 'Input rhcrit must be specified on model pressure levels but too many levels specified', FATAL)
      endif
    end if

    if(do_linear_diag) then
      if(scm_linear_coeffs(size(temp,3),1).eq.FILL_VALUE) then
        call error_mesg('cloud_simple', 'Input linear coefficients must be specified on model pressure levels but not enough levels specified', FATAL)
      endif
      if(scm_linear_coeffs(size(temp,3)+1,1).ne.FILL_VALUE) then
        call error_mesg('cloud_simple', 'Input linear coefficients must be specified on model pressure levels but too many levels specified', FATAL)
      endif
    end if

    ! Get the saturated specific humidity TOTAL (ie ice and vap) ***double check maths!
    call compute_qs(temp, p_full, qs)

    if(do_add_stratocumulus) then
      call calc_theta_dthdp(temp, p_full, p_half, theta, kdthdp, k700)
    end if

    k_surf = size(temp, 3)

    do i=1, size(temp, 1)
      do j=1, size(temp, 2)
        do k=1, size(temp, 3)
          ! caluclate the liquid fraction, effective radius, critical RH for
          ! the simple cloud scheme and cloud fraction.
          ! rh_in_cf is an output diagnostic only for debugging  
          call calc_liq_frac(temp(i,j,k), frac_liq(i,j,k))
          call calc_reff(frac_liq(i,j,k), reff_rad(i,j,k))

          if (do_simple_rhcrit) then
            call calc_rhcrit(p_full(i,j,k), p_full(i,j,k_surf), simple_rhcrit(i,j,k))
          end if
          if(do_read_scm_rhcrit) then
            simple_rhcrit(i,j,k) = scm_rhcrit(k)
          end if
          if(do_simple_rhcrit .OR. do_read_scm_rhcrit) then
            call calc_cf(q_hum(i,j,k), qs(i,j,k), simple_rhcrit(i,j,k), cf(i,j,k), rh_in_cf(i,j,k))
          end if

          if(do_linear_diag) then
            call calc_cf_linear(q_hum(i,j,k), qs(i,j,k), cf(i,j,k), rh_in_cf(i,j,k), scm_linear_coeffs(k,:))
          end if

          if(do_add_stratocumulus) then
            if(k .eq. kdthdp(i,j)) then
              strat = min(1.0, max(0.0, (theta(i,j,k700) - theta(i,j,k_surf)) * 0.057 - 0.5573))
              cf(i,j,k) = max(strat, cf(i,j,k))
            end if
          end if

          if(do_qcl_with_temp) then
            call calc_qcl_rad_with_temp(p_full(i,j,k), temp(i,j,k), cf(i,j,k), qcl_rad(i,j,k))
          else if (do_qcl_two_paras) then
            call calc_qcl_rad_two_paras(p_full(i,j,k), cf(i,j,k), qcl_rad(i,j,k))
          else
            call calc_qcl_rad(p_full(i,j,k), cf(i,j,k), qcl_rad(i,j,k))
          end if
        end do
      end do
    end do

    if (do_cloud_amount_diags) then
      call diag_cloud_amount(cf, p_full, p_half, tca, high_ca, mid_ca, low_ca)
    end if

    !save some diagnotics
    call output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, simple_rhcrit, &
                            tca, high_ca, mid_ca, low_ca, Time)
  end subroutine cloud_simple

  subroutine calc_theta_dthdp(temp, pfull, phalf, theta, kdthdp, k700)
    real,    intent(in),  dimension(:,:,:)  :: temp, pfull, phalf
    real,    intent(out), dimension(:,:,:)  :: theta
    integer, intent(out), dimension(:,:)    :: kdthdp
    integer, intent(out) :: k700
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: p0, dthdp
    real, dimension(size(temp,1), size(temp,2)) :: dthdpmn
    real :: cp, gas_const, premib
    integer :: i, j, k !, kp1, ks

    cp = 1005.0
    gas_const = 287.1
    p0 = 1.0e5
    dthdpmn = 0.0 !-0.125
    kdthdp = 0
    premib = 7.0e4
    dthdp = 0.0

    theta = temp * (p0 / pfull)**(gas_const / cp)

    do k=2, size(temp,3)
      do i=1, size(temp,1)
        do j=1, size(temp,2)
            if (phalf(i,j,k) >= premib) then
              dthdp(i,j,k) = (theta(i,j,k) - theta(i,j,k-1)) / (phalf(i,j,k) - phalf(i,j,k-1)) * 1.0e2
              if (dthdp(i,j,k) < dthdpmn(i,j)) then
                dthdpmn(i,j) = dthdp(i,j,k)
                kdthdp(i,j) = k     ! index of interface of max inversion
              end if
            end if
        end do
      end do
    end do
    k700 = minloc(abs(pfull(0,0,:) - 7.0e4), 1)
    ! write(*,*) 'Model level nearest 700mb is', k700, 'which is', pfull(0,0,k700), 'pascals.'
  end subroutine calc_theta_dthdp

  subroutine calc_liq_frac(temp, frac_liq)
    ! All liquid if above zero and all ice below -40C
    ! linearly interpolate between T=0 and -40C
    real, intent(in)    :: temp
    real, intent(out)   :: frac_liq
   
    if (temp > zerodegc) then
        frac_liq = 1.0
    else if (temp < zerodegc-40.0) then
        frac_liq = 0.0
    else           
        frac_liq = 1.0 - (zerodegc-temp) / 40.0
    end if
  end subroutine calc_liq_frac


  subroutine calc_reff(frac_liq, reff_rad)
    ! the effective cloud radius is bounded between 10 and 20 microns
    real, intent(in) :: frac_liq
    real, intent(out) :: reff_rad

    reff_rad =  10.0 * frac_liq + 20.0 * (1.0 - frac_liq)  !units in microns
  end subroutine calc_reff


  subroutine calc_rhcrit(p_full, p_surf, simple_rhcrit)   !need to check p_full - > p_layer_centres
    !get the RH needed as a threshold for the cloud fraction calc.
    real, intent(in)  :: p_full, p_surf
    real, intent(out) :: simple_rhcrit

    ! Calculate RHcrit as function of pressure
    if ( p_full > 7.0e4 ) then
      simple_rhcrit = rhcsfc - (rhcsfc - rhc700) *        &
                     (p_surf - p_full) / (p_surf - 7.0e4)
    else if ( p_full > 2.0e4 ) then
      simple_rhcrit = rhc700 - (rhc700 - rhc200) *        &
                      (7.0e4 - p_full) / 5.0e4
    else
      simple_rhcrit = rhc200
    endif
  end subroutine calc_rhcrit


  subroutine calc_cf (q_hum, qsat, simple_rhcrit, cf, rh)
    ! Calculate LS (stratiform) cloud fraction 
    ! as a simple linear function of RH
    real, intent(in)  :: q_hum, qsat, simple_rhcrit
    real, intent(out) :: cf, rh
    real :: cca

    rh = q_hum / qsat

    select case(cf_diag_formula)
      case(B_SPOOKIE)
        ! cf = MAX( 0.0, MIN( 1.0, (rh - simple_rhcrit) / (1.0 - simple_rhcrit) ))
        cf = (rh - simple_rhcrit) / (1.0 - simple_rhcrit)
      case(B_SUNDQVIST)
        cf = 1.0 - ((1.0 - MIN(rh,1.0)) / (1.0 - simple_rhcrit))**0.5
      case(B_SMITH)
        cf = 1.0 - (3.0 / sqrt(2.0) * (1.0 - MIN(rh,1.0))/(1.0 - simple_rhcrit))**(2.0/3.0)
      case default
        call error_mesg('cloud_simple', 'invalid cloud fraction diagnostic formula', FATAL)
    end select

    cf = MAX(0.0, MIN(1.0, cf))

    ! include simple convective cloud fraction where present (not currenly used)
    cca = 0.0 ! no convective cloud fraction is calculated
              ! left in for future use
    if (cca > 0.0) then
      cf = MAX(simple_cca, cf)
    end if
  end subroutine calc_cf

  subroutine calc_cf_linear(q_hum, qsat, cf, rh, linear_coeffs)
    ! Calculate LS (stratiform) cloud fraction as a linear function of RH
    real, intent(in)  :: q_hum, qsat
    real, intent(in), dimension(:) :: linear_coeffs
    real, intent(out) :: cf, rh

    rh = q_hum / qsat
    cf = linear_coeffs(1) * rh + linear_coeffs(2)
    cf = MAX(0.0, MIN(1.0, cf))
  end subroutine calc_cf_linear

  subroutine diag_cloud_amount(cf, p_full, p_half, tca, high_ca, mid_ca, low_ca)
    real, intent(in),  dimension(:,:,:) :: cf, p_full, p_half
    real, intent(out), dimension(:,:)   :: tca, high_ca, mid_ca, low_ca
    integer, dimension(size(cf,1),size(cf,2),size(cf,3)) :: ktop, kbot
    real,    dimension(size(cf,1),size(cf,2),size(cf,3)) :: cldamt, cloud
    integer, dimension(size(cf,1),size(cf,2)) :: nclds

    call max_rnd_overlap(cf, p_full, p_half, nclds, ktop, kbot, cldamt)
    call compute_tca_random(nclds, cldamt, tca)
    call expand_cloud(nclds, ktop, kbot, cldamt, cloud)
    call compute_isccp_clds(p_full, cloud, high_ca, mid_ca, low_ca)
  end subroutine diag_cloud_amount

  subroutine calc_qcl_rad(p_full, cf, qcl_rad)
    ! calculate cloud water content
    real , intent(in)   :: p_full, cf
    real , intent(out)  :: qcl_rad
    real :: in_cloud_qcl 

    in_cloud_qcl = 3.0e-4 + (1.0-3.0e-4) * (p_full-2.0e4) / 8.0e4
    in_cloud_qcl = MAX ( 0.0, in_cloud_qcl/1.0e3 ) ! convert to kg/kg
    qcl_rad = cf * in_cloud_qcl
  end subroutine calc_qcl_rad

  subroutine calc_qcl_rad_with_temp(p_full, temp, cf, qcl_rad)
    ! calculate cloud water content
    real , intent(in)   :: p_full, temp, cf
    real , intent(out)  :: qcl_rad
    real :: in_cloud_qcl

    in_cloud_qcl = 0.2 * (temp-220.0) / (280.0-220.0)
    in_cloud_qcl = max(3.0e-4, min(0.2, in_cloud_qcl)) / 1.0e3 ! convert to kg/kg
    qcl_rad = cf * in_cloud_qcl
  end subroutine calc_qcl_rad_with_temp

  subroutine calc_qcl_rad_two_paras(p_full, cf, qcl_rad)
    ! calculate cloud water content
    real , intent(in)   :: p_full, cf
    real , intent(out)  :: qcl_rad
    real :: in_cloud_qcl, qcl_1000, qcl_800, qcl_100

    qcl_1000 = 0.05
    qcl_800 = 0.175
    qcl_100 = 3.0e-4

    if ( p_full >= 8.0e4 ) then
      in_cloud_qcl = qcl_1000 + (qcl_800-qcl_1000) * (1.0e5-p_full) / 2.0e4
    else if ( p_full >= 1.0e4 ) then
      in_cloud_qcl = qcl_100 + (qcl_800-qcl_100) * (p_full-1.0e4) / 7.0e4
    else
      in_cloud_qcl = qcl_100
    end if

    in_cloud_qcl = in_cloud_qcl / 1.0e3 ! convert to kg/kg
    qcl_rad = cf * in_cloud_qcl
  end subroutine calc_qcl_rad_two_paras

  subroutine compute_tca_random(nclds, cldamt, tca)
    ! This subroutine was adapted from AM4 src/atmos_param/clouds/clouds.F90
    integer, intent(in)  :: nclds (:,:)
    real,    intent(in)  :: cldamt(:,:,:)
    real,    intent(out) :: tca   (:,:)
    integer :: i, j, k

    !---- compute total cloud amount assuming that -----
    !       independent clouds overlap randomly
    tca = 1.0
    do i=1,size(cldamt,1)
      do j=1,size(cldamt,2)
        do k = 1,nclds(i,j)
          tca(i,j) = tca(i,j) * (1.0 - cldamt(i,j,k))
        enddo
      enddo
    enddo
    tca = (1.0 - tca) * 1.0e2 ! unit percent
  end subroutine compute_tca_random

  subroutine expand_cloud(nclds, ktop, kbtm, cloud_in, cloud_out)
    integer, intent(in)  :: nclds(:,:), ktop(:,:,:), kbtm(:,:,:)
    real,    intent(in)  :: cloud_in (:,:,:)
    real,    intent(out) :: cloud_out(:,:,:)
    integer :: i, j, n

    cloud_out = 0.0
    do j=1,size(nclds,2)
      do i=1,size(nclds,1)
         do n=1,nclds(i,j)
           cloud_out(i,j,ktop(i,j,n):kbtm(i,j,n)) = cloud_in(i,j,n)
         enddo
      enddo
    enddo
  end subroutine expand_cloud

  subroutine compute_isccp_clds(pfull, cloud, high_ca, mid_ca, low_ca)
    !   define arrays giving the fractional cloudiness for clouds with
    !   tops within the ISCCP definitions of high (10-440 hPa), middle
    !   (440-680 hPa) and low (680-1000 hPa).
    real,  dimension(:,:,:),   intent(in)  :: pfull, cloud
    real,  dimension(:,:),     intent(out) :: high_ca, mid_ca, low_ca
    real,  parameter :: mid_btm = 6.8e4, high_btm = 4.4e4
    ! local array
    integer :: i, j, k

    !---- compute high, middle and low cloud amounts assuming that -----
    !       independent clouds overlap randomly
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    do j=1, size(cloud,2)
      do i=1, size(cloud,1)
        do k = 1, size(cloud,3)
          if (pfull(i,j,k)  <=  high_btm) then
            high_ca(i,j) = high_ca(i,j) * (1. - cloud(i,j,k))
          else if ((pfull(i,j,k) > high_btm) .and. (pfull(i,j,k) <= mid_btm)) then
            mid_ca(i,j) = mid_ca(i,j) * (1. - cloud(i,j,k))
          else if (pfull(i,j,k) > mid_btm ) then
            low_ca(i,j) = low_ca(i,j) * (1. - cloud(i,j,k))
         endif
        enddo
      enddo
    enddo

    high_ca = (1.0 - high_ca) * 1.0e2
    mid_ca = (1.0 - mid_ca) * 1.0e2
    low_ca = (1.0 - low_ca) * 1.0e2
  end subroutine compute_isccp_clds

  subroutine max_rnd_overlap(cf, pfull, phalf, nclds, ktop, kbot, cldamt)
    !max_rnd_overlap returns various cloud specification properties
    !    obtained with the maximum-random overlap assumption.

    ! intent(out) variables:
    !   nclds        Number of (random overlapping) clouds in column
    !   ktop         Level of the top of the cloud
    !   kbot         Level of the bottom of the cloud
    !   cldamt       Cloud amount of condensed cloud [ dimensionless ]
    real,    dimension(:,:,:), intent(in)             :: cf, pfull, phalf
    integer, dimension(:,:),   intent(out)            :: nclds
    integer, dimension(:,:,:), intent(out)            :: ktop, kbot
    real,    dimension(:,:,:), intent(out)            :: cldamt

    ! local variables:
    real, dimension (size(cf,1), size(cf,2), size(cf,3))  :: cldamt_cs
    integer    :: kdim
    integer    :: top_t, bot_t
    integer    :: tmp_top, tmp_bot, nlev
    logical    :: already_in_cloud, cloud_bottom_reached
    real       :: maxcldfrac
    real       :: totcld_bot, max_bot
    real       :: totcld_top, max_top, tmp_val
    integer    :: i, j, k, kc, t
    !--------------------------------------------------------------------
    !   local variables:
    !       kdim              number of model layers
    !       top_t             used temporarily as tag for cloud top index
    !       bot_t             used temporarily as tag for cloud bottom index
    !       tmp_top           used temporarily as tag for cloud top index
    !       tmp_bot           used temporarily as tag for cloud bottom index
    !       nlev              number of levels in the cloud
    !       already_in_cloud  if true, previous layer contained cloud
    !       cloud_bottom_reached
    !                         if true, the cloud-free layer beneath a cloud
    !                         has been reached
    !       maxcldfrac        maximum cloud fraction in any layer of cloud
    !                         [ fraction ]
    !       totcld_bot        total cloud fraction from bottom view
    !       max_bot           largest cloud fraction face from bottom view
    !       totcld_top        total cloud fraction from top view
    !       max_top           largest cloud fraction face from top view
    !       tmp_val           temporary number used in the assigning of top
    !                         [ (kg condensate / m**2) * microns ]
    !       i,j,k,kc,t        do-loop indices
    !----------------------------------------------------------------------

    !---------------------------------------------------------------------
    !    define the number of vertical layers in the model. initialize the
    !    output fields to correspond to the absence of clouds.
    !---------------------------------------------------------------------
    kdim     = size(cf,3)
    nclds    = 0
    ktop     = 1
    kbot     = 1
    cldamt   = 0.0
    !--------------------------------------------------------------------
    !    find the levels with cloud in each column. determine the vertical
    !    extent of each individual cloud, treating cloud in adjacent layers
    !    as components of a multi-layer cloud, and then calculate appropr-
    !    iate values of water paths and effective particle size.
    !--------------------------------------------------------------------
    do j=1,size(cf,2)
      do i=1,size(cf,1)
        ! set a flag indicating that we are searching for the next cloud top.
        already_in_cloud  = .false.
        cloud_bottom_reached = .false.
        ! march down the column.
        do k=1,kdim
          ! find a layer containing cloud in the column.
          if (cf(i,j,k) .gt. cf_min) then
            !--------------------------------------------------------------------
            !    if the previous layer was not cloudy, then a new cloud has been
            !    found. increment the cloud counter, set the flag to indicate the
            !    layer is in a cloud, save its cloud top level, initialize the
            !    values of its ice and liquid contents and fractional area and
            !    effective crystal and drop sizes.
            !--------------------------------------------------------------------
            if (.not. already_in_cloud)  then
              nclds(i,j) = nclds(i,j) + 1
              already_in_cloud = .true.
              cloud_bottom_reached = .false.
              ktop(i,j,nclds(i,j)) = k
              maxcldfrac = 0.0
            endif
            !---------------------------------------------------------------------
            !    add this layer's contributions to the current cloud. total liquid
            !    content, ice content, largest cloud fraction and condensate-
            !    weighted effective droplet and crystal radii are accumulated over
            !    the cloud.
            !---------------------------------------------------------------------
            maxcldfrac = MAX(maxcldfrac, cf(i,j,k))
          endif
          !--------------------------------------------------------------------
          !    when the cloud-free layer below a cloud is reached, or if the
          !    bottom model level is reached, define the cloud bottom level and
          !    set a flag indicating that mean values for the cloud may now be
          !    calculated.
          !--------------------------------------------------------------------
          if (cf(i,j,k) <= cf_min .and. already_in_cloud) then
            cloud_bottom_reached = .true.
            kbot(i,j,nclds(i,j)) = k - 1
          else if (already_in_cloud .and. k == kdim) then
            cloud_bottom_reached = .true.
            kbot(i,j,nclds(i,j)) = kdim
          endif
          !--------------------------------------------------------------------
          ! define the cloud fraction as the largest value of any layer in the cloud.
          !--------------------------------------------------------------------
          if (cloud_bottom_reached) then
            cldamt_cs(i,j,nclds(i,j)) = maxcldfrac
            !----------------------------------------------------------------------
            !    if adjust_top is true, the top and bottom indices of multi-layer
            !    clouds are adjusted to be those that are the most exposed to top
            !    and bottom view.
            !----------------------------------------------------------------------
            if (adjust_top) then
              ! define the cloud thickness.
              nlev = kbot(i,j,nclds(i,j)) - ktop(i,j,nclds(i,j)) + 1
              if (nlev > 1) then
                ! use the current top and bottom as the first guess for the new values.
                tmp_top = ktop(i,j,nclds(i,j))
                tmp_bot = kbot(i,j,nclds(i,j))
                ! initialize local search variables.
                totcld_bot = 0.
                totcld_top = 0.
                max_bot    = 0.
                max_top    = 0.
                !--------------------------------------------------------------------
                !    to find the adjusted cloud top, begin at current top and work
                !    downward. find the layer which is most exposed when viewed from
                !    the top; i.e., the cloud fraction increase is largest for that
                !    layer. the adjusted cloud base is found equivalently, starting
                !    from the actual cloud base and working upwards.
                !--------------------------------------------------------------------
                do t=1,nlev
                  ! find adjusted cloud top.
                  top_t   = ktop(i,j,nclds(i,j)) + t - 1
                  tmp_val = MAX(0., cf(i,j,top_t) - totcld_top)
                  if (tmp_val > max_top) then
                    max_top = tmp_val
                    tmp_top = top_t
                  end if
                  totcld_top = totcld_top + tmp_val
                  ! find adjusted cloud base.
                  bot_t   = kbot(i,j,nclds(i,j)) - t + 1
                  tmp_val = MAX(0., cf(i,j,bot_t) - totcld_bot)
                  if (tmp_val > max_bot) then
                    max_bot = tmp_val
                    tmp_bot = bot_t
                  end if
                  totcld_bot = totcld_bot + tmp_val
                end do
                ! assign tmp_top and tmp_bot as the new ktop and kbot.
                ktop(i,j,nclds(i,j)) = tmp_top
                kbot(i,j,nclds(i,j)) = tmp_bot
              endif  !(nlev > 1)
            endif ! (adjust_top)
            !---------------------------------------------------------------------
            !    reset already_in_cloud and cloud_bottom_reached to indicate that
            !    the current cloud has been exited.
            !---------------------------------------------------------------------
            already_in_cloud     = .false.
            cloud_bottom_reached = .false.
          endif   ! (cloud_bottom_reached)
        end do
      end do
    end do
    !---------------------------------------------------------------------
    !    place cloud properties into physical-space arrays for return to
    !    calling routine. NOTE THAT ALL LEVELS IN A GIVEN CLOUD ARE
    !    ASSIGNED THE SAME PROPERTIES.
    !---------------------------------------------------------------------
    do j=1,size(cf,2)
      do i=1,size(cf,1)
        do kc=1, nclds(i,j)
          !cldamt(i,j,ktop(i,j,kc):kbot(i,j,kc)) = cldamt_cs(i,j,kc)
          cldamt(i,j,kc) = cldamt_cs(i,j,kc)
        end do
      end do
    end do
  end subroutine max_rnd_overlap


  subroutine output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, simple_rhcrit, &
                                tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:,:) :: cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, simple_rhcrit
    real, intent(in), dimension(:,:)   :: tca, high_ca, mid_ca, low_ca
    type(time_type) , intent(in)       :: Time
    real :: used

    if ( id_cf > 0 ) then
      used = send_data ( id_cf, cf, Time)
    endif

    if (do_cloud_amount_diags) then
      if ( id_tot_cld_amt > 0 ) then
        used = send_data ( id_tot_cld_amt, tca, Time)
      endif

      if ( id_high_cld_amt > 0 ) then
        used = send_data ( id_high_cld_amt, high_ca, Time)
      endif

      if ( id_mid_cld_amt > 0 ) then
        used = send_data ( id_mid_cld_amt, mid_ca, Time)
      endif

      if ( id_low_cld_amt > 0 ) then
        used = send_data ( id_low_cld_amt, low_ca, Time)
      endif
    end if

    if ( id_reff_rad > 0 ) then
      used = send_data ( id_reff_rad, reff_rad, Time)
    endif

    if ( id_frac_liq > 0 ) then
      used = send_data ( id_frac_liq, frac_liq, Time)
    endif

    if ( id_qcl_rad > 0 ) then
      used = send_data ( id_qcl_rad, qcl_rad, Time)
    endif

    if ( id_rh_in_cf > 0 ) then
      used = send_data ( id_rh_in_cf, rh_in_cf*100., Time)
    endif

    if ( id_simple_rhcrit > 0 ) then
      used = send_data ( id_simple_rhcrit, simple_rhcrit*100.0, Time)
    endif
  end subroutine output_cloud_diags


  subroutine cloud_simple_end ()
  ! If alloocated are added in init then deallocate them here.
  end subroutine cloud_simple_end

  !-----------------------------------------------

end module cloud_simple_mod
