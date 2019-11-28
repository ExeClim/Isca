module cloud_simple_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  use             fms_mod, only: stdlog, FATAL, WARNING, NOTE, error_mesg, uppercase, &
                                 check_nml_error
  use    time_manager_mod, only: time_type
  use  sat_vapor_pres_mod, only: compute_qs, lookup_es
  use    diag_manager_mod, only: register_diag_field, send_data
  use       constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, HLV, KAPPA, ES0, KELVIN, RADIUS
  use             lcl_mod, only: lcl

  implicit none

  logical ::   do_init = .true.  ! Check if init needs to be run

  real :: zerodegc = 273.15
  integer :: id_cf, id_reff_rad, id_frac_liq, id_qcl_rad, id_rh_in_cf, id_rhcrit
  integer :: id_conv_cf
  ! ----- outputs for EIS, ECTEI and ELF diagnostics ----- !
  integer :: id_theta, id_dthdp, id_lts, id_eis, id_ectei, id_zlcl, &
             id_gamma_850, id_gamma_DL, id_gamma_700, id_z700, &
             id_zinv, id_ELF, id_beta1, id_beta2, id_IS, id_DS, id_alpha, &
             id_low_cld_amt_park, id_marine_strat
  ! ----- outputs for cloud amount diagnostics ----- !
  integer :: id_tot_cld_amt, id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt
  integer :: id_tot_cld_amt_mxr, id_high_cld_amt_mxr, id_mid_cld_amt_mxr, id_low_cld_amt_mxr, &
             id_tot_cld_amt_max, id_high_cld_amt_max, id_mid_cld_amt_max, id_low_cld_amt_max, &
             id_tot_cld_amt_rnd, id_high_cld_amt_rnd, id_mid_cld_amt_rnd, id_low_cld_amt_rnd, &
             id_tot_cld_amt_cam, id_high_cld_amt_cam, id_mid_cld_amt_cam, id_low_cld_amt_cam

  character(len=14), parameter ::   mod_name_cld = "cloud_simple"

  integer, parameter :: B_SPOOKIE=1, B_SUNDQVIST=2, B_LINEAR=3, B_SMITH=4, B_SLINGO=5, B_XR96=6
  integer, private :: cf_diag_formula = B_LINEAR
  character(len=32) :: cf_diag_formula_name = 'linear'
  character(len=32) :: sc_diag_method = 'Park_ELF'
  logical :: do_simple_rhcrit = .false.
  logical :: do_read_scm_rhcrit = .false.
  logical :: do_qcl_with_temp = .false.
  logical :: do_cloud_amount_diags = .true.
  logical :: do_test_overlap = .false.
  logical :: do_cam_cld_cover_diag = .true.
  logical :: do_add_stratocumulus = .false.
  logical :: intermediate_outputs_diags = .false.
  logical :: do_read_ts = .false.
  logical :: do_conv_cld = .false.
  logical :: do_adjust_low_cld = .false.
  logical :: adjust_top = .false.
  real, parameter :: FILL_VALUE = -999.0 ! Fill value for arrays
  real, dimension(100) :: scm_rhcrit = FILL_VALUE   ! Input array for single column critical RH. Max number of levels = 100

  real :: simple_cca = 0.0
  real :: rhcsfc     = 0.95
  real :: rhc700     = 0.7
  real :: rhc200     = 0.3
  real :: cf_min     = 1e-4
  real :: dthdp_min_threshold = -0.05 !K/hPa, which is -0.125 in CESM1.2.1
  real :: pshallow   = 7.5e4 ! copy from am4 diag_cloud.F90
  real :: conv_rain_min = 0.14 ! mm/day, threshold to produce conv cld

  ! Parameters to control the coefficients profile of linear function of RH
  real :: a_surf     = 2.7
  real :: a_top      = 0.35
  real :: b_surf     = -2.2
  real :: b_top      = -0.07
  real :: nx         = 8

  ! For slingo80 scheme
  real :: slingo_rhc_low  = 0.8
  real :: slingo_rhc_mid  = 0.65
  real :: slingo_rhc_high = 0.8

  ! For low clout adjustment
  real :: omega_adj_threshold = -0.1 !Pa/s

  ! For effective cloud droplet radius
  real :: reff_liq = 14.0 ! mm
  real :: reff_ice = 25.0 

  ! For in-cloud liquid water
  real :: qcl_val = 0.2 ! g/kg

  ! Define constants for Earth mass and Newtonian gravational constant
  ! Refer to: https://github.com/Unidata/MetPy/ --> src/metpy/constants.py
  real :: EARTH_MASS = 5.9722e24 ! kg
  ! Refer to: https://physics.nist.gov/cgi-bin/cuu/Value?bg
  real :: GRAV_CONST = 6.674e-11 !m^3 / kg / s^2

  namelist /cloud_simple_nml/ simple_cca, rhcsfc, rhc700, rhc200, &
                              cf_diag_formula_name, &
                              do_read_scm_rhcrit, scm_rhcrit, &
                              do_qcl_with_temp, &
                              do_cloud_amount_diags, adjust_top, do_test_overlap, &
                              do_add_stratocumulus, sc_diag_method, &
                              intermediate_outputs_diags, &
                              dthdp_min_threshold, do_read_ts, &
                              a_surf, a_top, b_surf, b_top, nx, &
                              do_conv_cld, pshallow, cf_min, &
                              slingo_rhc_low, slingo_rhc_mid, slingo_rhc_high, &
                              do_adjust_low_cld, omega_adj_threshold, &
                              conv_rain_min, qcl_val, reff_liq, reff_ice, &
                              do_cam_cld_cover_diag


  contains

  !-----------------------------------------------

  subroutine convective_cloud_init(axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes

    id_conv_cf = register_diag_field (mod_name_cld, 'conv_cf', axes(1:3), Time, &
        'Convective cloud fraction for the simple cloud scheme', 'unitless: values 0-1')
  end subroutine convective_cloud_init

  subroutine marine_strat_cloud_init(axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    character(len=32) :: method_str = ''

    call error_mesg('cloud_simple', 'The stratomulus diagnosis method is '// &
                    uppercase(trim(sc_diag_method)), NOTE)

    method_str = uppercase(trim(sc_diag_method))

    if (method_str(1:3)=='EIS' .or. method_str(1:5)=='ECTEI' .or. method_str(1:4)=='PARK') then
      id_eis = register_diag_field (mod_name_cld, 'eis', axes(1:2), Time, &
                      'estimated inversion strength', 'K')
    end if

    if (method_str(1:5)=='ECTEI' .or. method_str(1:4)=='PARK') then
      id_ectei = register_diag_field (mod_name_cld, 'ectei', axes(1:2), Time, &
                      'estimated cloud top entrainment index', 'K')
    end if

    if (method_str(1:4)=='PARK') then
      id_ELF    = register_diag_field (mod_name_cld, 'ELF', axes(1:2), Time, &
                      'estimated low cloud fraction', '')
    end if

    id_marine_strat = register_diag_field ( mod_name_cld, 'marine_strat', axes(1:3), Time, &
                        'marine low stratus cloud amount', '0-1' )

    id_zlcl = register_diag_field (mod_name_cld, 'zlcl', axes(1:2), Time, &
                  'height of lcl', 'meter')
    id_theta = register_diag_field (mod_name_cld, 'theta', axes(1:3), Time, &
                  'potential temperature', 'K')
    id_lts = register_diag_field (mod_name_cld, 'lts', axes(1:2), Time, &
                  'low-tropospheric stability', 'K')
    if(intermediate_outputs_diags) then
      id_dthdp = register_diag_field (mod_name_cld, 'dthdp', axes(1:3), Time, &
                        'dtheta/dp', 'K/hPa' )
      id_z700 = register_diag_field ( mod_name_cld, 'z700', axes(1:2), Time, &
                        'height of 700mb', 'meter')
      if (method_str(1:3)=='EIS') then
        id_gamma_850 = register_diag_field (mod_name_cld, 'gamma850', axes(1:2), Time, &
                        'moist lapse rate at 850hPa', 'K/m')
      end if
      if (method_str(1:4)=='PARK') then
        id_beta1  = register_diag_field (mod_name_cld, 'beta1', axes(1:2), Time, &
                        'first low-level cloud suppression parameter', '')
        id_beta2  = register_diag_field (mod_name_cld, 'beta2', axes(1:2), Time, &
                        'second low-level cloud suppression parameter', '')
        id_zinv   = register_diag_field (mod_name_cld, 'zinv', axes(1:2), Time, &
                        'height of invesion layer', 'meter')
        id_DS     = register_diag_field (mod_name_cld, 'DS', axes(1:2), Time, &
                        'decoupling strength', 'K')
        id_IS     = register_diag_field (mod_name_cld, 'IS', axes(1:2), Time, &
                        'inversion strength', 'K')
        id_alpha  = register_diag_field (mod_name_cld, 'alpha', axes(1:2), Time, &
                        'decoupling parameter', '')
        id_low_cld_amt_park = register_diag_field ( mod_name_cld, 'low_cld_amt_park', axes(1:2), Time, &
                        'low cloud amount estimated from Park method', 'percent' )
        id_gamma_DL = register_diag_field (mod_name_cld, 'gamma_DL', axes(1:2), Time, &
                        'moist lapse rate at decoupling layer', 'K/m')
        id_gamma_700 = register_diag_field (mod_name_cld, 'gamma700', axes(1:2), Time, &
                        'moist lapse rate at 700hPa', 'K/m')
      end if
    end if
  end subroutine marine_strat_cloud_init

  subroutine cloud_simple_init (axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    integer :: io, ierr, nml_unit, stdlog_unit
    character(len=32) :: method_str = ''

#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=cloud_simple_nml, iostat=io)
    ierr = check_nml_error(io, 'cloud_simple_nml')
#else
    if (file_exist('input.nml')) then
      nml_unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
         read(nml_unit, nml=cloud_simple_nml, iostat=io, end=10)
         ierr = check_nml_error(io, 'cloud_simple_nml')
      enddo
10    call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, cloud_simple_nml)

    ! Select cloud fraction diag formula
    method_str = uppercase(trim(cf_diag_formula_name))
    if(method_str == 'SPOOKIE') then
      cf_diag_formula = B_SPOOKIE
      call error_mesg('cloud_simple', 'Using default SPOOKIE cloud fraction diagnostic formula.', NOTE)
    else if(method_str == 'SUNDQVIST') then
      cf_diag_formula = B_SUNDQVIST
      call error_mesg('cloud_simple', 'Using Sundqvist (1987) cloud fraction diagnostic formula.', NOTE)
    else if(method_str == 'LINEAR') then
      cf_diag_formula = B_LINEAR
      call error_mesg('cloud_simple', 'Using linear cloud fraction diagnostic formula.', NOTE)
    else if(method_str == 'SMITH') then
      cf_diag_formula = B_SMITH
      call error_mesg('cloud_simple', 'Using Smith (1990) cloud fraction diagnostic formula.', NOTE)
    else if(method_str == 'SLINGO') then
      cf_diag_formula = B_SLINGO
      call error_mesg('cloud_simple', 'Using Slingo (1980) cloud fraction diagnostic formula.', NOTE)
    else if(method_str == 'XR96') then
      cf_diag_formula = B_XR96
      call error_mesg('cloud_simple', 'Using Xu and Krueger (1996) cloud fraction diagnostic formula.', NOTE)
    else
      call error_mesg('cloud_simple', '"'//trim(cf_diag_formula_name)//'"'// &
                ' is not a valid cloud fraction diagnostic formula.', FATAL)
    endif

    if (cf_diag_formula .eq. B_SPOOKIE .or. cf_diag_formula .eq. B_SUNDQVIST .or. &
        cf_diag_formula .eq. B_SMITH) then
      do_simple_rhcrit = .true.
    else
      do_simple_rhcrit = .false.
    end if

    !register diagnostics
    id_cf = register_diag_field (mod_name_cld, 'cf', axes(1:3), Time, &
                    'Cloud fraction for the simple cloud scheme', 'unitless: values 0-1')
    id_frac_liq = register_diag_field (mod_name_cld, 'frac_liq', axes(1:3), Time, &
                    'Liquid cloud fraction (liquid, mixed-ice phase, ice)', 'unitless: values 0-1')
    id_reff_rad = register_diag_field (mod_name_cld, 'reff_rad', axes(1:3), Time, &
                    'Effective cloud particle radius', 'microns')
    id_qcl_rad = register_diag_field (mod_name_cld, 'qcl_rad', axes(1:3), Time, &
                    'Specific humidity of cloud liquid', 'kg/kg')
    id_rh_in_cf = register_diag_field (mod_name_cld, 'rh_in_cf', axes(1:3), Time, &
                    'RH as a percent', '%')

    if(do_simple_rhcrit) then
      id_rhcrit = register_diag_field (mod_name_cld, 'rhcrit', axes(1:3), Time, &
                    'RH as a percent', '%')
    end if

    if (do_cloud_amount_diags) then
      id_tot_cld_amt = register_diag_field (mod_name_cld, 'tot_cld_amt', axes(1:2), Time, &
                            'total cloud amount', 'percent')
      id_high_cld_amt = register_diag_field (mod_name_cld, 'high_cld_amt', axes(1:2), Time, &
                            'high cloud amount', 'percent')
      id_mid_cld_amt = register_diag_field (mod_name_cld, 'mid_cld_amt', axes(1:2), Time, &
                            'mid cloud amount', 'percent')
      id_low_cld_amt = register_diag_field (mod_name_cld, 'low_cld_amt', axes(1:2), Time, &
                            'low cloud amount', 'percent')
      if (do_test_overlap) then
        id_tot_cld_amt_mxr = register_diag_field (mod_name_cld, 'tot_cld_amt_mxr', axes(1:2), Time, &
                'total cloud amount', 'percent')
        id_high_cld_amt_mxr = register_diag_field (mod_name_cld, 'high_cld_amt_mxr', axes(1:2), Time, &
                'high cloud amount', 'percent')
        id_mid_cld_amt_mxr = register_diag_field (mod_name_cld, 'mid_cld_amt_mxr', axes(1:2), Time, &
                'mid cloud amount', 'percent')
        id_low_cld_amt_mxr = register_diag_field (mod_name_cld, 'low_cld_amt_mxr', axes(1:2), Time, &
                'low cloud amount', 'percent')
        ! Max overlap
        id_tot_cld_amt_max = register_diag_field (mod_name_cld, 'tot_cld_amt_max', axes(1:2), Time, &
                'total cloud amount', 'percent')
        id_high_cld_amt_max = register_diag_field (mod_name_cld, 'high_cld_amt_max', axes(1:2), Time, &
                'high cloud amount', 'percent')
        id_mid_cld_amt_max = register_diag_field (mod_name_cld, 'mid_cld_amt_max', axes(1:2), Time, &
                'mid cloud amount', 'percent')
        id_low_cld_amt_max = register_diag_field (mod_name_cld, 'low_cld_amt_max', axes(1:2), Time, &
                'low cloud amount', 'percent')
        ! Random overlap
        id_tot_cld_amt_rnd = register_diag_field (mod_name_cld, 'tot_cld_amt_rnd', axes(1:2), Time, &
                'total cloud amount', 'percent')
        id_high_cld_amt_rnd = register_diag_field (mod_name_cld, 'high_cld_amt_rnd', axes(1:2), Time, &
                'high cloud amount', 'percent')
        id_mid_cld_amt_rnd = register_diag_field (mod_name_cld, 'mid_cld_amt_rnd', axes(1:2), Time, &
                'mid cloud amount', 'percent')
        id_low_cld_amt_rnd = register_diag_field (mod_name_cld, 'low_cld_amt_rnd', axes(1:2), Time, &
                'low cloud amount', 'percent')
      end if
      if (do_cam_cld_cover_diag) then
        id_tot_cld_amt_cam = register_diag_field (mod_name_cld, 'tot_cld_amt_cam', axes(1:2), Time, &
                'total cloud amount', 'percent')
        id_high_cld_amt_cam = register_diag_field (mod_name_cld, 'high_cld_amt_cam', axes(1:2), Time, &
                'high cloud amount', 'percent')
        id_mid_cld_amt_cam = register_diag_field (mod_name_cld, 'mid_cld_amt_cam', axes(1:2), Time, &
                'mid cloud amount', 'percent')
        id_low_cld_amt_cam = register_diag_field (mod_name_cld, 'low_cld_amt_cam', axes(1:2), Time, &
                'low cloud amount', 'percent')
      end if
    end if

    if(do_add_stratocumulus) then
      call marine_strat_cloud_init(axes, Time)
    end if

    if (do_conv_cld) then
      call convective_cloud_init(axes, Time)
    end if

    do_init = .false.  !initialisation completed
  end subroutine cloud_simple_init


  subroutine cloud_simple(p_half, p_full, Time, temp, q_hum, z_full, &
                          wg_full, psg, temp_2m, q_2m, rh_2m, precip, klcls, klzbs, &
                          cf, reff_rad, qcl_rad)  ! outs
    real, intent(in),  dimension(:,:,:) :: temp, q_hum, p_full, p_half, z_full, wg_full
    real, intent(in),  dimension(:,:)   :: psg, temp_2m, precip, q_2m, rh_2m
    integer, intent(in), dimension(:,:) :: klcls, klzbs
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: cf, reff_rad, qcl_rad
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: qs, frac_liq, rh_in_cf, &
                                                                 rhcrit, conv_cf, rh_e

    !check initiation has been done - ie read in parameters
    if (do_init) call error_mesg ('cloud_simple',  &
         'cloud_simple_init has not been called.', FATAL)

    ! Get the saturated specific humidity TOTAL (ie ice and vap) ***double check maths!
    call compute_qs(temp, p_full, qs)
    rh_in_cf = q_hum / qs

    call calc_liq_frac(temp, frac_liq)
    call calc_reff(frac_liq, reff_rad)

    if (do_simple_rhcrit) then
      call calc_rhcrit(p_full, rhcrit, scm_rhcrit)
    end if

    conv_cf = 0.0
    if(do_conv_cld) then
      !call calc_convective_cf(p_full, precip, klcls, conv_cf, Time)
      call calc_convective_cf2(p_full, precip, klcls, klzbs, conv_cf, Time)
      !call add_anvil_clouds(p_full, precip, klcls, klzbs, conv_cf, Time)
    end if

    ! rh_e is the effective RH
    rh_e = rh_in_cf * (1.0 - conv_cf)
    call calc_stratiform_cf(p_full, psg, rh_e, q_hum, qs, rhcrit, qcl_rad, cf)

    if(do_adjust_low_cld) then
      call adjust_low_cld(p_full, wg_full, cf)
    end if

    if (do_add_stratocumulus) then
      call add_stratiform_cld(temp, p_full, p_half, z_full, rh_e, q_hum, &
                              temp_2m, q_2m, rh_2m, psg, wg_full, klcls, cf, Time)
    end if

    if(do_conv_cld) then
      call merge_strat_conv_clouds(cf, conv_cf)
    end if

    call calc_qcl_rad(p_full, temp, cf, qcl_rad)

    if (do_cloud_amount_diags) then
      call diag_cloud_amount(cf, p_full, p_half, Time)
      if (do_test_overlap) then
        call diag_cldamt_maxrnd_overlap(cf, p_full, Time)
        call diag_cldamt_max_overlap(cf, p_full, Time)
        call diag_cldamt_random_overlap(cf, p_full, Time)
      end if
      if (do_cam_cld_cover_diag) then
        call cloud_cover_diags_cam(cf, p_full, p_half, Time)
      end if
    end if

    call output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, rhcrit, Time)

  end subroutine cloud_simple

  subroutine adjust_low_cld(p_full, wg_full, cf)
    real, intent(in),    dimension(:,:,:) :: p_full, wg_full
    real, intent(inout), dimension(:,:,:) :: cf
    real :: premib, omega

    premib = 7.0e4 !Pa
    !omega_adj_threshold = -0.1   !Pa/s

    !where(p_full>premib .and. omega_adj_threshold<wg_full .and. wg_full<0.0)
    !  cf = min(1.0, wg_full/omega_adj_threshold) * cf
    !elsewhere (p_full>premib .and. wg_full>=0.0)
    where (p_full>premib .and. wg_full>=0.0)
    !elsewhere (wg_full>=0.0)
      cf = 0.0
      !cf = min(1.0, wg_full/abs(omega_adj_threshold)) * cf
    !elsewhere
    !  cf = cf
    end where
  end subroutine adjust_low_cld


  subroutine estimate_stratiform_cld(method_str, i, j, k, kb, pfull, cf, rh, theta, eis, dthdp, ectei, ELF)
    implicit none
    integer, intent(in) :: i, j, k
    character(len=32), intent(in) :: method_str
    real, intent(in),  dimension(:,:,:) :: rh, theta, pfull, dthdp
    integer, intent(in) :: kb
    real, intent(in),  dimension(:,:)  :: eis, ectei, ELF
    real, intent(out), dimension(:,:,:) :: cf
    real :: strat, rhb_frac
    integer :: k700, k_surf
    
    k_surf = size(pfull, 3)
    k700 = minloc(abs(pfull(i,j,:) - 7.0e4), 1)
  
    if(method_str == 'LTS') then
      strat = min(1.0, max(0.0, (theta(i,j,k700) - theta(i,j,k_surf)) * 0.057 - 0.5573))
      cf(i,j,k) = max(strat, cf(i,j,k))
    else if(method_str == 'SLINGO') then
      strat = min(1.0, max(0.0, -6.67*dthdp(i,j,k) - 0.667))
      rhb_frac = min(1.0, max(0.0, (rh(i,j,kb) - 0.6) / 0.2))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
    else if(method_str == 'SLINGO_NO_RH') then
      strat = min(1.0, max(0.0, -6.67*dthdp(i,j,k) - 0.667))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
    else if(method_str == 'DTHDP') then
      strat = min(1.0, max(0.0, -3.1196*dthdp(i,j,k) - 0.1246))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
    else if(method_str == 'EIS_WOOD') then
      !strat = min(1.0, max(0.0, 0.0221*eis(i,j) + 0.1128))
      strat = min(1.0, max(0.0, 0.06*eis(i,j) + 0.14)) ! Wood and Betherton, 2006
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
    else if(method_str == 'EIS_WOOD_RH') then
      strat = min(1.0, max(0.0, 0.06*eis(i,j)+0.14)) !* (rh(i,j,kb)-0.6)/0.2)) ! Wood and Betherton, 2006
      rhb_frac = min(1.0, max(0.0, (rh(i,j,kb)-0.6) / 0.2))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
    else if(method_str == 'EIS_RH') then
      strat = min(1.0, max(0.0, (0.092*eis(i,j)+ 0.027)*(2.078*rh(i,j,k)-6.45e-19)))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
    else if(method_str == 'ECTEI') then
      ! Kawai, Koshiro and Webb, 2017
      strat = min(1.0, max(0.0, 0.031*ectei(i,j) + 0.39))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
    else if(method_str == 'ECTEI_RH') then
      ! Kawai, Koshiro and Webb, 2017
      strat = min(1.0, max(0.0, 0.031*ectei(i,j) + 0.39))
      rhb_frac = min(1.0, max(0.0, (rh(i,j,kb)-0.6) / 0.2))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
    else if(method_str == 'PARK_ELF') then
      ! Park and Shin, 2019, ACP
      !strat = min(1.0, max(0.0, 0.86*ELF(i,j) + 0.02))
      strat = min(1.0, max(0.0, 1.272*ELF(i,j)-0.366))
      !strat = min(1.0, max(0.0, 1.11*ELF(i,j)-0.107))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
    else
      call error_mesg('cloud_simple', method_str//' is not supported yet!', FATAL)
    end if
  end subroutine estimate_stratiform_cld


  subroutine add_stratiform_cld(temp, p_full, p_half, z_full, rh, q_hum, temp_2m, q_2m, rh_2m, psg, wg_full, klcls, cf, Time)
    implicit none
    real, intent(in),  dimension(:,:,:) :: temp, q_hum, p_full, p_half, z_full, rh, wg_full
    type(time_type),   intent(in)       :: Time
    real, intent(in),  dimension(:,:)   :: temp_2m, q_2m, rh_2m, psg
    integer, intent(in), dimension(:,:) :: klcls
    real, intent(out), dimension(:,:,:) :: cf
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: theta, dthdp, marine_strat
    integer, dimension(size(temp,1), size(temp,2)) :: kdthdp, kinvs
    real,    dimension(size(temp,1), size(temp,2)) :: eis, ectei, ELF, low_ca_park
    real :: strat, rhb_frac, used, omega_pos_threshold
    character(len=32) :: method_str = ''
    integer :: i, j, k, k700, kb, k_surf, kk, nlev

    eis = 0.0
    ectei = 0.0
    ELF = 0.0
    dthdp = 0.0

    call calc_theta_dthdp(temp, temp_2m, p_full, p_half, psg, theta, dthdp, kdthdp)

    method_str = uppercase(trim(sc_diag_method))
    if (method_str(1:3)=='EIS' .or. method_str(1:5)=='ECTEI') then
      call calc_eis(p_full, z_full, temp, temp_2m, psg, klcls, eis, Time)
    end if
    if (method_str(1:5)=='ECTEI') then
      call calc_ectei(p_full, q_hum, q_2m, eis, ectei, Time)
    end if
    if (method_str(1:4)=='PARK') then
      call calc_Park_proxies(p_full, psg, z_full, temp, temp_2m, q_hum, q_2m, rh_2m, klcls, ELF, kinvs, Time)
    end if

    k_surf = size(temp, 3)
    omega_pos_threshold = 0. !1.4*100/3600
    marine_strat = 0.0

    do i=1, size(temp, 1)
      do j=1, size(temp, 2)

        ! =========== Add off-coast marine stratiform clods =========== !
        !if (method_str(1:4)=='PARK') then
        !  kk = kinvs(i,j)
        !  !write(*,*) 'QL kinv,kk',kinvs(i,j), kdthdp(i,j)
        !else
        !  kk = kdthdp(i,j)
        !end if

        kk = kdthdp(i,j)
        kb = min(kk+1, k_surf)
        !do k = kk-1, kb
        do k=kk,kk
        !k = kk
          !if (kk.ne.0 .and. wg_full(i,j,k)>omega_pos_threshold .and. dthdp(i,j,k)<dthdp_min_threshold) then
          if (kk.ne.0 .and. wg_full(i,j,k)>omega_pos_threshold .and. dthdp(i,j,k)<dthdp_min_threshold .and. p_full(i,j,k)>8.5e4) then
          !if (kk.ne.0 .and. dthdp(i,j,k)<dthdp_min_threshold) then
          ! .and. dthdp(i,j,k)<dthdp_min_threshold
          !if (kk.ne.0 .and. wg_full(i,j,k)>omega_pos_threshold .and. p_full(i,j,k)>7.5e4 ) then
            call estimate_stratiform_cld(method_str, i, j, k, kb, p_full, cf, rh, theta, eis, dthdp, ectei, ELF)
            !cf(i,j,k) = min(1.0, max(0.0, cf(i,j,k)*2))
            !cf(i,j,k) = min(1.0, abs(dthdp(i,j,k)/dthdp(i,j,kk))) * cf(i,j,k)
            marine_strat(i,j,k) = min(1.0, max(0.0, cf(i,j,k)))
          end if
        end do

        ! =========== Other stratiform clouds except temperature inversion areas =========== !
        ! This is the first try, the souther ocean is improved a lot,
        ! but subtropical region provide too much low clouds 
        !if (method_str(1:4)=='PARK' .and. p_full(i,j,kk)>7.5e4) then
        
        !if (method_str(1:4)=='PARK') then
        !  kk = kinvs(i,j)
        !end if
        !!if (method_str(1:4)=='PARK' .and. p_full(i,j,kk)>7.5e4 .and. dthdp(i,j,kk)>-0.06) then
        !if (method_str(1:4)=='PARK' .and. p_full(i,j,kk)>7.5e4 &
        !    .and. wg_full(i,j,kk)<0 .and. dthdp(i,j,kk)>-0.06) then
        !  call estimate_stratiform_cld(method_str, i, j, kk, kb, p_full, cf, rh, theta, eis, dthdp, ectei, ELF)
        !  marine_strat(i,j,kk) = min(1.0, max(0.0, cf(i,j,kk)))
        !end if

        !nlev = 0
        !do k=1,k_surf
        !  if(p_full(i,j,k)>7.0e4 .and. k<=klcls(i,j)) then
        !    !call estimate_stratiform_cld(method_str, i, j, k, kb, pfull, cf, rh, theta, eis, dthdp, ectei, ELF)
        !   nlev = nlev + 1
        !  end if
        !end do
        !call estimate_stratiform_cld(method_str, i, j, k, kb, pfull, cf, rh, theta, eis, dthdp, ectei, ELF)
        !nlev = klcls(i,j) - k700 + 1
        !do k=k700,klcls(i,j)
        !  cf(i,j,k) = 1.0 - (1.0 - cf(i,j,k))**(1/nlev) !min(1.0, max(0.0, 1.0-(1.0-cf(i,j,k))**(1/nlev))) ! Random overlap
        !end do
      end do
    end do

    if (id_theta > 0) then
      used = send_data(id_theta, theta, Time)
    end if
    if (id_marine_strat > 0) then
      used = send_data(id_marine_strat, marine_strat, Time)
    end if

    if(intermediate_outputs_diags) then
      if (id_dthdp > 0) then
        used = send_data(id_dthdp, dthdp, Time)
      end if
      if (id_low_cld_amt_park > 0) then
        used = send_data(id_low_cld_amt_park, low_ca_park, Time)
      end if
    end if
  end subroutine add_stratiform_cld

  !subroutine add_stratiform_cld(temp, p_full, p_half, z_full, rh, q_hum, temp_2m, q_2m, rh_2m, psg, wg_full, klcls, cf, Time)
  !  implicit none
  !  real, intent(in),  dimension(:,:,:) :: temp, q_hum, p_full, p_half, z_full, rh, wg_full
  !  type(time_type),   intent(in)       :: Time
  !  real, intent(in),  dimension(:,:)   :: temp_2m, q_2m, rh_2m, psg
  !  integer, intent(in), dimension(:,:) :: klcls
  !  real, intent(out), dimension(:,:,:) :: cf
  !  real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: theta, dthdp, marine_strat
  !  integer, dimension(size(temp,1), size(temp,2)) :: kdthdp
  !  real,    dimension(size(temp,1), size(temp,2)) :: eis, ectei, ELF, low_ca_park
  !  real :: strat, rhb_frac, used, omega_pos_threshold
  !  character(len=32) :: method_str = ''
  !  integer :: i, j, k, k700, kb, k_surf, kk
  !
  !  k_surf = size(temp, 3)
  !  omega_pos_threshold = 1.0*100/3600
  !
  !  call calc_theta_dthdp(temp, temp_2m, p_full, p_half, psg, theta, dthdp, kdthdp)
  !
  !  method_str = uppercase(trim(sc_diag_method))
  !  if (method_str(1:3)=='EIS' .or. method_str(1:5)=='ECTEI') then
  !    call calc_eis(p_full, z_full, temp, temp_2m, psg, klcls, eis, Time)
  !  end if
  !  if (method_str(1:5)=='ECTEI') then
  !    call calc_ectei(p_full, q_hum, q_2m, eis, ectei, Time)
  !  end if
  !  if (method_str(1:4)=='PARK') then
  !    call calc_Park_proxies(p_full, psg, z_full, temp, temp_2m, q_hum, q_2m, rh_2m, klcls, ELF, Time)
  !  end if
  !
  !  marine_strat = 0.0
  !  do i=1, size(temp, 1)
  !    do j=1, size(temp, 2)
  !      kk = kdthdp(i,j)
  !      kb = min(kk+1, k_surf)
  !      !do k=kk-2, kb
  !      k = kk
  !        if (kk.ne.0 .and. wg_full(i,j,k)>omega_pos_threshold .and. dthdp(i,j,k)<dthdp_min_threshold) then
  !          if(method_str == 'LTS') then
  !            strat = min(1.0, max(0.0, (theta(i,j,k700) - theta(i,j,k_surf)) * 0.057 - 0.5573))
  !            cf(i,j,k) = max(strat, cf(i,j,k))
  !          end if
  !          if(method_str == 'SLINGO') then
  !            strat = min(1.0, max(0.0, -6.67*dthdp(i,j,k) - 0.667))
  !            rhb_frac = min(1.0, max(0.0, (rh(i,j,kb) - 0.6) / 0.2))
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
  !          end if
  !          if(method_str == 'SLINGO_NO_RH') then
  !            strat = min(1.0, max(0.0, -6.67*dthdp(i,j,k) - 0.667))
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
  !          end if
  !          if(method_str == 'DTHDP') then
  !            strat = min(1.0, max(0.0, -3.1196*dthdp(i,j,k) - 0.1246))
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
  !          end if
  !          if(method_str == 'EIS_WOOD') then
  !            !strat = min(1.0, max(0.0, 0.0221*eis(i,j) + 0.1128))
  !            strat = min(1.0, max(0.0, 0.06*eis(i,j) + 0.14)) ! Wood and Betherton, 2006
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
  !          end if
  !          if(method_str == 'EIS_WOOD_RH') then
  !            strat = min(1.0, max(0.0, 0.06*eis(i,j)+0.14)) !* (rh(i,j,kb)-0.6)/0.2)) ! Wood and Betherton, 2006
  !            rhb_frac = min(1.0, max(0.0, (rh(i,j,kb)-0.6) / 0.2))
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
  !          end if
  !          if(method_str == 'EIS_RH') then
  !            strat = min(1.0, max(0.0, (0.092*eis(i,j)+ 0.027)*(2.078*rh(i,j,k)-6.45e-19)))
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
  !          end if
  !          if(method_str == 'ECTEI') then
  !            ! Kawai, Koshiro and Webb, 2017
  !            strat = min(1.0, max(0.0, 0.031*ectei(i,j) + 0.39))
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
  !          end if
  !          if(method_str == 'ECTEI_RH') then
  !            ! Kawai, Koshiro and Webb, 2017
  !            strat = min(1.0, max(0.0, 0.031*ectei(i,j) + 0.39))
  !            rhb_frac = min(1.0, max(0.0, (rh(i,j,kb)-0.6) / 0.2))
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
  !          end if
  !          if(method_str == 'PARK_ELF') then
  !            ! Park and Shin, 2019, ACP
  !            !strat = min(1.0, max(0.0, 0.86*ELF(i,j) + 0.02))
  !            strat = min(1.0, max(0.0, 1.272*ELF(i,j)-0.366))
  !            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
  !          end if
  !          cf(i,j,k) = min(1.0, dthdp(i,j,k)/dthdp(i,j,kk)) * cf(i,j,k)
  !          marine_strat(i,j,k) = min(1.0, max(0.0, cf(i,j,k)))
  !        end if
  !      !end do
  !
  !      if(intermediate_outputs_diags .and. method_str=='PARK_ELF') then
  !        !low_ca_park(i,j) = min(1.0, max(0.0, 0.86*ELF(i,j) + 0.02))
  !        low_ca_park(i,j) = min(1.0, max(0.0, 1.272*ELF(i,j)-0.366))
  !      end if
  !
  !    end do
  !  end do
  !
  !  if (id_theta > 0) then
  !    used = send_data(id_theta, theta, Time)
  !  end if
  !  if (id_marine_strat > 0) then
  !    used = send_data(id_marine_strat, marine_strat, Time)
  !  end if
  !
  !  if(intermediate_outputs_diags) then
  !    if (id_dthdp > 0) then
  !      used = send_data(id_dthdp, dthdp, Time)
  !    end if
  !    if (id_low_cld_amt_park > 0) then
  !      used = send_data(id_low_cld_amt_park, low_ca_park, Time)
  !    end if
  !  end if
  !end subroutine add_stratiform_cld

  subroutine calc_theta_dthdp(temp, temp_2m, pfull, phalf, ps, theta, dthdp, kdthdp)
    real,    intent(in),  dimension(:,:,:) :: temp, pfull, phalf
    real,    intent(in),  dimension(:,:)   :: temp_2m, ps
    real,    intent(out), dimension(:,:,:) :: theta, dthdp
    integer, intent(out), dimension(:,:)   :: kdthdp
    real, dimension(size(temp,1), size(temp,2)) :: dthdp_min, theta_0
    real :: premib, pstar
    integer :: i, j, k, kb

    dthdp_min = dthdp_min_threshold  !d_theta / d_p *1e2, lapse rate
    kdthdp = 0
    premib = 8.0e4
    dthdp = 0.0
    pstar = 1.0e5

    kb = size(temp,3)  !bottom level

    do k=1,kb
      theta(:,:,k) =  temp(:,:,k) * (pstar / pfull(:,:,k))**(RDGAS / CP_AIR)
    end do

    do k=1,kb-1
      dthdp(:,:,k) = (theta(:,:,k) - theta(:,:,k+1)) / (phalf(:,:,k) - phalf(:,:,k+1)) * 1.0e2
    end do

    theta_0 = temp_2m * (pstar / ps)**(RDGAS / CP_AIR)
    dthdp(:,:,kb) = (theta(:,:,kb) - theta_0) / (phalf(:,:,kb) - ps) * 1.0e2

    kdthdp = minloc(dthdp, dim=3, mask=(phalf>premib).and.(dthdp<dthdp_min_threshold))

  end subroutine calc_theta_dthdp


  subroutine calc_liq_frac(temp, frac_liq)
    ! All liquid if above zero and all ice below -40C
    ! linearly interpolate between T=0 and -40C
    real, intent(in),  dimension(:,:,:) :: temp
    real, intent(out), dimension(:,:,:) :: frac_liq
    real :: t_max, t_min

    t_max = -5
    t_min = -40

    !where (temp > zerodegc)
    where (temp > zerodegc+t_max)
        frac_liq = 1.0
    elsewhere (temp < zerodegc+t_min)
        frac_liq = 0.0
    elsewhere
        !frac_liq = 1.0 - (zerodegc-temp) / 40.0
        frac_liq = (temp-zerodegc - t_min) / (t_max-t_min)
    end where

  end subroutine calc_liq_frac

  subroutine calc_reff(frac_liq, reff_rad)
    ! The effective cloud radius is bounded between 10 and 20 microns
    real, intent(in),  dimension(:,:,:) :: frac_liq
    real, intent(out), dimension(:,:,:) :: reff_rad

    !reff_rad =  10.0 * frac_liq + 20.0 * (1.0 - frac_liq)  !units in microns
    !reff_rad =  14.0 * frac_liq + 25.0 * (1.0 - frac_liq)  !units in microns
    reff_rad =  reff_liq * frac_liq + reff_ice * (1.0 - frac_liq)  !units in microns

  end subroutine calc_reff

  subroutine calc_rhcrit(p_full, rhcrit, scm_rhcrit)
    !get the RH needed as a threshold for the cloud fraction calc.
    real, intent(in),  dimension(:,:,:) :: p_full
    real, intent(in),  dimension(:)     :: scm_rhcrit
    real, intent(out), dimension(:,:,:) :: rhcrit
    real    :: p_surf
    integer :: i, j, k

    p_surf = 1e5

    if (do_simple_rhcrit) then
      ! Calculate RHcrit as function of pressure
      where (p_full > 7.0e4)
        rhcrit = rhcsfc - (rhcsfc - rhc700) * (p_surf - p_full) / (p_surf - 7.0e4)
      elsewhere (p_full > 2.0e4)
        rhcrit = rhc700 - (rhc700 - rhc200) * (7.0e4 - p_full) / 5.0e4
      elsewhere
        rhcrit = rhc200
      endwhere
    end if

    if (do_read_scm_rhcrit) then
      if(scm_rhcrit(size(p_full,3)) .eq. FILL_VALUE) then
        call error_mesg('cloud_simple', 'Input rhcrit must be specified on model '// &
                        'pressure levels but not enough levels specified', FATAL)
      endif
      if(scm_rhcrit(size(p_full,3)+1) .ne. FILL_VALUE) then
        call error_mesg('cloud_simple', 'Input rhcrit must be specified on model '// &
                        'pressure levels but too many levels specified', FATAL)
      endif

      do i=1, size(p_full,1)
        do j=1, size(p_full,2)
            rhcrit(i,j,:) = scm_rhcrit(1:size(p_full,3))
        end do
      end do
    end if

  end subroutine calc_rhcrit

  subroutine calc_stratiform_cf(pfull, ps, rh, q_hum, qsat, rhcrit, qcl_rad, cf)
    ! Calculate LS (stratiform) cloud fraction as a simple linear function of RH
    real, intent(in),  dimension(:,:,:) :: pfull, rh, q_hum, qsat, rhcrit, qcl_rad
    real, intent(in),  dimension(:,:)   :: ps
    real, intent(out), dimension(:,:,:) :: cf
    real, dimension(size(pfull,1), size(pfull,2), size(pfull,3)) :: rhc
    real :: mid_top, mid_base, p_para, alpha_0, gamma ! For Xu and Krueger (1996)
    integer :: i, j, k

    select case(cf_diag_formula)
      case(B_SPOOKIE)
        cf = (rh - rhcrit) / (1.0 - rhcrit)

      case(B_SUNDQVIST)
        cf = 1.0 - ((1.0 - MIN(rh,1.0)) / (1.0 - rhcrit))**0.5

      case(B_SMITH)
        cf = 1.0 - (3.0 / sqrt(2.0) * (1.0 - MIN(rh,1.0))/(1.0 - rhcrit))**(2.0/3.0)

      case(B_SLINGO)
        mid_base = 8.0e4
        mid_top = 4.0e4
        rhc = 0.0

        where (pfull > mid_base)
          rhc = slingo_rhc_low
        elsewhere (pfull < mid_top)
          rhc = slingo_rhc_high
        elsewhere
          rhc = slingo_rhc_mid
        end where

        where (rh<rhc)
          cf = 0.0
        elsewhere
          cf = ((rh-rhc)/(1-rhc)) ** 2
        end where

      case(B_XR96)
        p_para = 0.25
        alpha_0 = 100.0
        gamma = 0.49

        where (rh.ge.1)
          cf = 1.0
        elsewhere (rh.le.0)
          cf = 0.0
        elsewhere
          cf = rh**p_para * (1.0 - EXP(-alpha_0*qcl_rad / (qsat-q_hum)**gamma))
          ! If only use rh**p_para, then the cld fraction is too large...
        end where

        do i=1, size(cf,1)
          do j=1, size(cf,2)
            do k=1, size(cf,3)
            !write(*,*) 'sum(qcl_rad)=', sum(qcl_rad)
            !if (sum(cf)>0) then
              !write(*,*) 'size=', size(qcl_rad,1), size(qcl_rad,2), size(qcl_rad,3)
              !write(*,*) 'qsat=', qsat(1, 1, 24) * 1e3, 'q_hum=', q_hum(1, 1, 24) * 1e3 , &
              !     'qcl_rad=',qcl_rad(1, 1, 24) * 1e3, 'cf=', cf(1, 1, 24), 'rh=', rh(1, 1, 24)
            !end if
            !write(*,*) 'qcl_rad=', qcl_rad(:, 63, 127) * 1e3
            !write(*,*) 'size qcl_rad=', size(qcl_rad,1), size(qcl_rad,2), size(qcl_rad,3)
            !write(*,*) 'qcl_rad=', qcl_rad * 1e3
            if (isnan(cf(i,j,k))) then
              write(*,*) 'NaN i,j,k=', i, j,k, 'qsat=', qsat(i,j,k) * 1e3, 'q_hum=', q_hum(i,j,k) * 1e3 , &
              'qcl_rad=',qcl_rad(i,j,k) * 1e3, 'cf=', cf(i,j,k), 'rh=', rh(i,j,k)
            end if
            if (qcl_rad(i,j,k)>0) then
              !write(*,*) 'size=', size(qcl_rad,1), size(qcl_rad,2), size(qcl_rad,3)
              write(*,*) 'i,j,k=', i, j,k, 'qsat=', qsat(i,j,k) * 1e3, 'q_hum=', q_hum(i,j,k) * 1e3 , &
                   'qcl_rad=',qcl_rad(i,j,k) * 1e3, 'cf=', cf(i,j,k), 'rh=', rh(i,j,k), &
                   'term2=', (1.0 - EXP(-alpha_0*qcl_rad(i,j,k) / (qsat(i,j,k)-q_hum(i,j,k))**gamma))
            end if
            end do
          end do
        end do

      case(B_LINEAR)
        call calc_cf_linear(pfull, rh, ps, cf)

      case default
        call error_mesg('cloud_simple', 'invalid cloud fraction diagnostic formula', FATAL)
    end select

    cf = MAX(0.0, MIN(1.0, cf))

  end subroutine calc_stratiform_cf

  subroutine calc_cf_linear(p_full, rh, ps, cf)
    ! Calculate LS (stratiform) cloud fraction as a linear function of RH
    real, intent(in),  dimension(:,:,:) :: p_full, rh
    real, intent(in),  dimension(:,:)   :: ps
    real, intent(out), dimension(:,:,:) :: cf
    real, dimension(size(p_full,1), size(p_full,2), size(p_full,3)) :: coeff_a, coeff_b
    integer :: k
    !real :: ps

    !ps = 1.0e5 !Pa
    !coeff_a = a_top + (a_surf-a_top) * exp(1.0 - (ps/p_full)**nx)
    !coeff_b = b_top + (b_surf-b_top) * exp(1.0 - (ps/p_full)**nx)

    do k=1,size(rh,3)
      coeff_a(:,:,k) = a_top + (a_surf-a_top) * exp(1.0 - (ps/p_full(:,:,k))**nx)
      coeff_b(:,:,k) = b_top + (b_surf-b_top) * exp(1.0 - (ps/p_full(:,:,k))**nx)
    end do

    cf = coeff_a * rh + coeff_b
    cf = MAX(0.0, MIN(1.0, cf))
  end subroutine calc_cf_linear

  subroutine calc_qcl_rad(p_full, temp, cf, qcl_rad)
    ! calculate cloud water content
    real, intent(in),  dimension(:,:,:) :: p_full, cf, temp
    real, intent(out), dimension(:,:,:) :: qcl_rad
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: in_cloud_qcl

    if (do_qcl_with_temp) then
      !in_cloud_qcl = 0.2 * (temp-220.0) / (280.0-220.0)
      !in_cloud_qcl = MAX(3.0e-4, MIN(0.2, in_cloud_qcl)) / 1.0e3 ! convert to kg/kg
      in_cloud_qcl =  qcl_val * (temp-220.0) / (280.0-220.0)
      in_cloud_qcl = MAX(3.0e-4, MIN(qcl_val, in_cloud_qcl)) / 1.0e3 ! convert to kg/kg

      !in_cloud_qcl = 0.15 * (temp-220.0) / (280.0-220.0)
      !in_cloud_qcl = MAX(3.0e-4, MIN(0.15, in_cloud_qcl)) / 1.0e3 ! convert to kg/kg
    else
      ! in_cloud_qcl as a function of height
      in_cloud_qcl = 3.0e-4 + (1.0-3.0e-4) * (p_full-2.0e4) / 8.0e4
      in_cloud_qcl = MAX(0.0, in_cloud_qcl/1.0e3) ! convert to kg/kg
    end if

    qcl_rad = cf * in_cloud_qcl
  end subroutine calc_qcl_rad

  subroutine calc_convective_cf(pfull, precip, klcls, conv_cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(in),  dimension(:,:)   :: precip
    integer, intent(in), dimension(:,:) :: klcls
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: conv_cf
    real, dimension(size(pfull,1), size(pfull,2)) :: precip_mm_per_day, plcl
    real    :: a, b, tower_scale_coeff, used
    integer :: k

    call calc_lcls(klcls, pfull=pfull, plcls=plcl)

    precip_mm_per_day = precip * 24.0 * 3600.0  ! change units to mm/day

    a = -0.125 * log(0.14)  !0.246
    b = 0.125
    tower_scale_coeff = 0.25

    conv_cf = 0.0

    do k=1, size(pfull,3)
      where (precip_mm_per_day<0.14)
        conv_cf(:,:,k) = 0.0
      elsewhere (precip_mm_per_day>85.0)
        conv_cf(:,:,k) = 0.8
      elsewhere
        conv_cf(:,:,k) = a + b * log(precip_mm_per_day)
      end where
      ! Convective clouds only exist when above LCL
      where(pfull(:,:,k) > plcl)
        conv_cf(:,:,k) = 0.0
      end where
    end do

    where (pfull < pshallow)
      conv_cf = conv_cf * tower_scale_coeff
    end where
    !!!!!!conv_cf = conv_cf * tower_scale_coeff

    ! Output the diagnostics
    if (id_conv_cf > 0) then
      used = send_data(id_conv_cf, conv_cf, Time)
    endif

  end subroutine calc_convective_cf

  subroutine calc_convective_cf2(pfull, precip, klcls, klzbs, conv_cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(in),  dimension(:,:)   :: precip
    integer, intent(in), dimension(:,:) :: klcls, klzbs
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: conv_cf
    real, dimension(size(pfull,1), size(pfull,2)) :: precip_mm_per_day, plcl2d, conv_cf_tmp
    real    :: convcld_a, convcld_b, used, tower_scale_coeff
    integer :: i, j, k, klcl, ktop, nlayers, k_tower

    call calc_lcls(klcls, pfull=pfull, plcls=plcl2d)

    ! conv_rain_min = 4 mm/day
    precip_mm_per_day = precip * 24.0 * 3600.0  ! change units to mm/day

    !convcld_a = -0.125 * log(0.14)  !0.246, ! 0.001
    convcld_a = -0.125 * log(conv_rain_min)
    convcld_b =  0.125 !0.0418 !
    tower_scale_coeff = 0.25

    conv_cf = 0.0
    !conv_cf_tmp = convcld_a + convcld_b * log(1.0 + precip_mm_per_day)

    where (precip_mm_per_day<conv_rain_min)
        conv_cf_tmp = 0.0
    !elsewhere(precip_mm_per_day>85.0)
    !    conv_cf_tmp = 0.8
    elsewhere
        conv_cf_tmp = convcld_a + convcld_b * log(precip_mm_per_day)
        conv_cf_tmp = min(0.8, max(0.0, conv_cf_tmp))
    end where

    !ktop = 3
    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        klcl = klcls(i,j)
        ktop = klzbs(i,j)
        !if (ktop.ne.0 .and. klcl.ne.ktop) then
        if (ktop.ne.0 .and. pfull(i,j,ktop)<4e4 .and. klcl.ne.ktop) then
          !nlayers = klcl - ktop
          ! Random overlap assumption
          ! conv_cf(i,j,ktop:klcl) = 1.0 - (1.0 - conv_cf_tmp(i,j))**(1.0 / nlayers)
          ! Maxmimum overlap
          !conv_cf(i,j,ktop:klcl) = conv_cf_tmp(i,j)*0.25
          ! Third try:
          conv_cf(i,j,ktop:klcl) = conv_cf_tmp(i,j) !*0.5 !*0.5
          !write(*,*) 'QL', i,j, klcl, ktop, nlayers, conv_cf(i,j,klcl), pfull(i,j,ktop)
          ! k_tower = minloc(abs(pfull(i,j,:) - pshallow), 1)
          ! if (k_tower>ktop .and. k_tower<klcl) then
          !   conv_cf(i,j,ktop:k_tower) = conv_cf_tmp(i,j) * tower_scale_coeff
          ! end if
        end if
      end do
    end do

    !where (pfull < pshallow)
    !  conv_cf = conv_cf * tower_scale_coeff
    !end where

    !conv_cf = min(0.6, max(0.0, conv_cf))

    ! Output the diagnostics
    if (id_conv_cf > 0) then
      used = send_data(id_conv_cf, conv_cf, Time)
    endif

  end subroutine calc_convective_cf2

  subroutine add_anvil_clouds(pfull, precip, klcls, klzbs, conv_cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(in),  dimension(:,:)   :: precip
    integer, intent(in), dimension(:,:) :: klcls, klzbs
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: conv_cf
    real, dimension(size(pfull,1), size(pfull,2)) :: precip_mm_per_day, plcl2d, conv_cf_tmp
    real    :: convcld_a, convcld_b, used, tower_scale_coeff
    integer :: i, j, k, klcl, ktop, nlayers, k_tower

    call calc_lcls(klcls, pfull=pfull, plcls=plcl2d)

    precip_mm_per_day = precip * 24.0 * 3600.0  ! change units to mm/day

    convcld_a = -0.125 * log(0.14)  !0.246, ! 0.001
    convcld_b =  0.125 !0.0418 !
    tower_scale_coeff = 0.25

    conv_cf = 0.0
    conv_cf_tmp = convcld_a + convcld_b * log(1.0 + precip_mm_per_day)

    !ktop = 3
    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        klcl = klcls(i,j)
        ktop = klzbs(i,j)
        if (ktop.ne.0 .and. klcl.ne.ktop) then
          !conv_cf(i,j,ktop:klcl) = conv_cf_tmp(i,j) !*0.5 !*0.5
          k_tower = minloc(abs(pfull(i,j,:) - 4.e4), 1)
          if (k_tower>ktop .and. k_tower<klcl) then
            !conv_cf(i,j,ktop:k_tower) = min(0.6, max(0.0, (conv_cf_tmp(i,j)-0.3)*2))
            conv_cf(i,j,ktop:k_tower) = conv_cf_tmp(i,j) * 0.25 !*0.5
            !write(*,*) 'QL anvil', conv_cf_tmp(i,j)*0.25
          end if
        end if
      end do
    end do

    !conv_cf = min(0.6, max(0.0, conv_cf))

    ! Output the diagnostics
    if (id_conv_cf > 0) then
      used = send_data(id_conv_cf, conv_cf, Time)
    endif

  end subroutine add_anvil_clouds


  subroutine merge_strat_conv_clouds(cf, conv_cf)
    implicit none
    real, intent(in),  dimension(:,:,:) :: conv_cf
    real, intent(out), dimension(:,:,:) :: cf

    cf = max(cf, conv_cf)
    !cf = 1.0 - (1.0 - cf) * (1.0 - conv_cf) ! Random overlap
    ! Refer to CAM5
    !cf = (1-conv_cf) * cf + conv_cf
  end subroutine merge_strat_conv_clouds

  function geopotential_to_height(geopot) result(height)
    ! Calculates the height from geopotential
    ! See https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.geopotential_to_height.html
    ! and https://github.com/Unidata/MetPy --> src/metpy/calc/basic.py

    implicit none
    real, intent(in),  dimension(:,:,:) :: geopot
    real, dimension(size(geopot,1),size(geopot,2),size(geopot,3)):: height, scaled

    scaled = geopot * RADIUS
    height = scaled * RADIUS / (GRAV_CONST * EARTH_MASS - scaled)
  end function geopotential_to_height

  subroutine calc_lcls(klcls, pfull, temp, zfull, ts, ps, rh_surf, plcls, tlcls, zlcls)
    ! Example to call:
    ! call calc_lcls(klcls, pfull=p_full, plcls=plcl2d)
    ! rh_surf in range [0,1]
    implicit none
    integer, intent(in), dimension(:,:) :: klcls
    real, intent(in),  dimension(:,:,:), optional :: temp, pfull, zfull
    real, intent(in),  dimension(:,:),   optional :: rh_surf, ts, ps
    real, intent(out), dimension(:,:),   optional :: plcls, tlcls, zlcls
    integer :: i, j
    !real :: zlcl_l

    do i=1, size(klcls,1)
      do j=1, size(klcls,2)
        if(present(pfull) .and. present(plcls)) then
          plcls(i,j) = pfull(i,j,klcls(i,j))
        end if
        if(present(temp) .and. present(tlcls)) then
          tlcls(i,j) =  temp(i,j,klcls(i,j))
        end if
        if (present(zfull) .and. present(zlcls)) then
          zlcls(i,j) = zfull(i,j,klcls(i,j))
        end if
        if (present(rh_surf) .and. present(ts) .and. present(ps) .and. present(zlcls)) then
          ! Bolton (1980) equation
          ! Refer to: https://journals.ametsoc.org/doi/pdf/10.1175/JAS-D-17-0102.1
          !zlcls(i,j) = max(0.0, CP_AIR/GRAV * (ts(i,j)-55.0 - (1.0/(ts(i,j)-55.0) - log(rh_surf(i,j))/2840.)**(-1)))
          ! The following one is not so correct as it is c_pm not c_p should be used in
          !zlcls(i,j) = CP_AIR/GRAV * (ts(i,j)-55.0 - (1.0/(ts(i,j)-55.0) - log(rh_surf(i,j))/2840.)**(-1))
          !write(*,*) 'QL',  i,j, rh_surf(i,j), zlcls(i,j)

          ! Use the exact LCL formula from D. M. Romps [2017, JAS 74(12)]
          zlcls(i,j) = lcl(ps(i,j), ts(i,j), rh=rh_surf(i,j))
          
          !old_zlcl = CP_AIR/GRAV * (ts(i,j)-55.0 - (1.0/(ts(i,j)-55.0) - log(rh_surf(i,j))/2840.)**(-1))
          !write(*,*) 'QL i,j, zlcl, zlcl_old, diff=',  i, j, zlcls(i,j), old_zlcl, old_zlcl-zlcls(i,j)
          ! if (zlcls(i,j)<0) then
          !   zlcl_l = lcl(ps(i,j), ts(i,j), rh=rh_surf(i,j), return_ldl=.true.)
          !   write(*,*) 'QL i,j, zlcl, zlcl_l, ps, ts, rhs=',  i, j, zlcls(i,j), zlcl_l, ps(i,j), ts(i,j), rh_surf(i,j)
          ! end if
        end if

        if(.not.((present(pfull) .and. present(plcls)) .or. &
                 (present(temp)  .and. present(tlcls)) .or. &
                 (present(zfull) .and. present(zlcls)) .or. &
                 (present(rh_surf) .and. present(ts) .and. present(ps) .and. present(zlcls)))) then
          call error_mesg('calc_lcls in cloud_simple', 'At least one group of '// &
                'pfull(plcls), temp(tlcls) and zfull/rh_surf(zlcls) should exist.', FATAL)
        end if
      end do
    end do

  end subroutine calc_lcls

  subroutine calc_eis(pfull, zfull, temp, ts, ps, klcls, eis, Time)
    ! Estimated inversion stability
    ! Refer to: Wood and Bretherton, 2006, Journal of Climate
    implicit none
    real,    intent(in),  dimension(:,:,:) :: pfull, zfull, temp
    real,    intent(in),  dimension(:,:)   :: ts, ps
    integer, intent(in),  dimension(:,:)   :: klcls
    type(time_type),      intent(in)       :: Time
    real,    intent(out), dimension(:,:)   :: eis
    real, dimension(size(temp,1), size(temp,2)) :: zlcl, z700, Gamma850, LTS
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: zfull_height
    real    :: pstar, T850, es, qs, used
    integer :: k700, i, j

    zfull_height = geopotential_to_height(zfull*GRAV)
    pstar = 1.e5 ! Pa

    do i=1, size(temp,1)
      do j=1, size(temp,2)
        k700 = minloc(abs(pfull(i,j,:) - 7.0e4), 1)
        LTS(i,j) = temp(i,j,k700)*((pstar/pfull(i,j,k700))**(RDGAS/CP_AIR)) - &
                   ts(i,j)*(pstar/ps(i,j))**(RDGAS/CP_AIR)
        T850 = (temp(i,j,k700) + ts(i,j)) / 2.0
        call lookup_es(T850, es)
        qs = 0.622 * es / (8.5e4 - es)
        Gamma850(i,j) = GRAV / CP_AIR * (1.0 - (1.0 + HLV*qs/RDGAS/T850) &
                        / (1.0 + HLV**2 * qs/CP_AIR/RVGAS/T850**2))
        !z700(i,j) = zfull(i,j,k700)
        z700(i,j) = zfull_height(i,j,k700)
      end do
    end do

    call calc_lcls(klcls, zfull=zfull_height, zlcls=zlcl)
    eis = LTS - Gamma850 * (z700 - zlcl)

    ! ----- output diagnositics ------ !
    if(id_eis > 0) then
      used = send_data (id_eis, eis, Time)
    end if
    if(id_lts > 0) then
      used = send_data (id_lts, LTS, Time)
    end if
    if(id_zlcl > 0) then
      used = send_data (id_zlcl, zlcl, Time)
    end if

    if(intermediate_outputs_diags) then
      if(id_z700 > 0) then
        used = send_data (id_z700, z700, Time)
      end if
      if(id_gamma_850 > 0) then
        used = send_data (id_gamma_850, Gamma850, Time)
      end if
    end if
  end subroutine calc_eis

  subroutine calc_ectei(pfull, q_hum, q_surf, eis, ectei, Time)
    ! Estimated Cloud Top Entrainment Index (ECTEI)
    ! Refer to: Eq(3) in Kawai, Koshiro and Webb, 2017, Journal of Climate
    implicit none
    real, intent(in),  dimension(:,:,:) :: pfull, q_hum
    real, intent(in),  dimension(:,:)   :: q_surf, eis
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:)   :: ectei
    real, dimension(size(pfull,1),size(pfull,2)) :: q_700
    integer :: k700, i, j
    real :: k_en, C_qgap, beta, used

    k_en = 0.7
    C_qgap = 0.76
    beta = (1.0 - k_en) * C_qgap

    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        k700 = minloc(abs(pfull(i,j,:) - 7.0e4), 1)
        q_700(i,j) = q_hum(i,j,k700)
      end do
    end do

    ectei = eis - beta * HLV / CP_AIR * (q_surf - q_700)

    if(id_ectei > 0) then
      used = send_data (id_ectei, ectei, Time)
    end if
  end subroutine calc_ectei

  subroutine calc_Park_proxies(pfull, ps, zfull, temp, ts, q_hum, q_surf, rh_surf, klcls, ELF, kinvs, Time)
    ! Refer to: Park and Shin, 2019, Atmospheric Chemistry and Physics
    ! https://www.atmos-chem-phys.net/19/5635/2019/

    implicit none
    real,    intent(in),  dimension(:,:,:) :: pfull, zfull, temp, q_hum
    real,    intent(in),  dimension(:,:)   :: ts, q_surf, ps, rh_surf
    integer, intent(in),  dimension(:,:)   :: klcls
    type(time_type),      intent(in)       :: Time
    real,    intent(out), dimension(:,:)   :: ELF
    integer, intent(out), dimension(:,:)   :: kinvs
    real, dimension(size(temp,1), size(temp,2)) :: plcl, tlcl, zlcl, z700, Gamma_DL, &
                                    Gamma700, LTS, z_ML, zinv, qv_ML, beta2
    ! other paramters
    real, dimension(size(temp,1), size(temp,2)) :: beta1, IS, DS, eis, ectei, alpha, f_para !qs_surf, 
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: zfull_height
    real :: pstar, delta_zs, theta_ML, used
    integer :: k700, i, j

    !call compute_qs(ts, ps, qs_surf)
    !rh_surf = q_surf / qs_surf

    delta_zs = 2750.0 ! meter, constant
    pstar = 1.0e5 ! Pa
    kappa = RDGAS / CP_AIR

    zfull_height = geopotential_to_height(zfull*GRAV)

    !write(*,*) 'QL max of rh_2m', maxval(rh_surf)
    call calc_lcls(klcls, pfull=pfull, temp=temp, ts=ts, ps=ps, rh_surf=rh_surf, &
                   plcls=plcl, tlcls=tlcl, zlcls=zlcl)

    where(zlcl < 0)
      zlcl = 0.0
    end where

    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        k700 = minloc(abs(pfull(i,j,:) - 7.0e4), 1)
        !z700(i,j) = zfull(i,j,k700)
        z700(i,j) = zfull_height(i,j,k700)

        !if (abs(zfull(i,j,k700)-z700(i,j))>5) then
        !  write(*,*) 'QL i,j, z700_g, z700_h=', i, j, zfull(i,j,k700), z700(i,j)
        !end if

        ! Mixed Layer is the LCL
        call calc_moist_lapse_rate(tlcl(i,j), plcl(i,j), Gamma_DL(i,j))
        call calc_moist_lapse_rate(temp(i,j,k700), pfull(i,j,k700), Gamma700(i,j))
        !theta_ML = tlcl(i,j) * (pstar / plcl(i,j))**kappa
        theta_ML = ts(i,j) * (pstar / ps(i,j))**kappa
        LTS(i,j) = temp(i,j,k700) * (pstar / pfull(i,j,k700))**kappa - theta_ML
        qv_ML(i,j) = q_hum(i,j,klcls(i,j))
      end do
    end do

    z_ML = zlcl
    zinv = -LTS/Gamma700 + z700 + delta_zs*(Gamma_DL/Gamma700)

    ! Rest zinv
    where(zinv < z_ML)
      zinv = z_ML
    end where
    where(zinv > z_ML+delta_zs)
      zinv = z_ML + delta_zs
    end where

    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        !kinvs(i,j) = minloc(abs(zinv(i,j)-zfull(i,j,:)), 1)
        kinvs(i,j) = minloc(abs(zinv(i,j)-zfull_height(i,j,:)), 1)
        !write(*,*) 'QL i,j, kinvs=', i,j,kinvs(i,j)
        !write(*,*) 'QL i,j, zinv=', i,j,zinv(i,j)
      end do
    end do

    ! low-level cloud suppression parameters (LCS)
    beta2 = sqrt(zinv*zlcl) / delta_zs
    ! freeze-dry factor (Vavrus and Waliser, 2008)
    f_para = max(0.15, min(1.0, qv_ML/0.003))
    ! Estimated low-cloud fraction (ELF)
    ELF = f_para * (1.0 - beta2)

    ! ----- output diagnostics ----- !
    if(id_ELF>0) then
      used = send_data(id_ELF, ELF, Time)
    end if
    if(id_lts>0) then
      used = send_data(id_lts, LTS, Time)
    end if
    if(id_zlcl>0) then
      used = send_data(id_zlcl, zlcl, Time)
    end if

    if(intermediate_outputs_diags) then
      !============= Other prameters =============!
      beta1 = (zinv + zlcl) / delta_zs
      alpha = (zinv - z_ML) / delta_zs
      IS = (1.0 - alpha) * Gamma_DL * delta_zs
      DS = alpha * Gamma_DL * delta_zs
      eis = LTS + Gamma_DL*z_ML - Gamma700*z700
      call calc_ectei(pfull, q_hum, q_surf, eis, ectei, Time)
      !============= Other prameters =============!

      !Add some diagnostic ouputs
      call output_extra_diags_for_Park_ELF(Time, beta1, beta2, &
              alpha, eis, IS, DS, z700, zinv, Gamma700, Gamma_DL)
    end if
  end subroutine calc_Park_proxies

  subroutine calc_moist_lapse_rate(T, p, Gamma)
    real, intent(in)  :: T, p
    real, intent(out) :: Gamma
    real :: es, qs

    call lookup_es(T, es)
    qs = 0.622 * es / (p - es)
    Gamma = (GRAV/CP_AIR) * (1.0 - (1.0 + HLV*qs/RDGAS/T) / (1.0 + HLV**2 * qs/CP_AIR/RVGAS/T**2))
  end subroutine calc_moist_lapse_rate

  subroutine diag_cloud_amount(cf, p_full, p_half, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full, p_half
    type(time_type),   intent(in)       :: Time
    real,    dimension(size(cf,1),size(cf,2)) :: tca, high_ca, mid_ca, low_ca
    integer, dimension(size(cf,1),size(cf,2),size(cf,3)) :: ktop, kbot
    real,    dimension(size(cf,1),size(cf,2),size(cf,3)) :: cldamt, cloud
    integer, dimension(size(cf,1),size(cf,2)) :: nclds

    call max_rnd_overlap(cf, p_full, p_half, nclds, ktop, kbot, cldamt)
    call compute_tca_random(nclds, cldamt, tca)
    !call expand_cloud(nclds, ktop, kbot, cldamt, cloud)
    call compute_isccp_clds2(p_full, nclds, ktop, cldamt, high_ca, mid_ca, low_ca)

    ! Diagnostics output
    call output_cldamt(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cloud_amount

  subroutine max_rnd_overlap(cf, pfull, phalf, nclds, ktop, kbot, cldamt)
    !max_rnd_overlap returns various cloud specification properties
    !    obtained with the maximum-random overlap assumption.

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

    kdim     = size(cf,3)
    nclds    = 0
    ktop     = 1
    kbot     = 1
    cldamt   = 0.0

    do j=1,size(cf,2)
      do i=1,size(cf,1)
        ! set a flag indicating that we are searching for the next cloud top.
        already_in_cloud  = .false.
        cloud_bottom_reached = .false.
        ! march down the column.
        do k=1,kdim
          ! find a layer containing cloud in the column.
          if (cf(i,j,k) .gt. cf_min) then
            if (.not. already_in_cloud)  then
              nclds(i,j) = nclds(i,j) + 1
              already_in_cloud = .true.
              cloud_bottom_reached = .false.
              ktop(i,j,nclds(i,j)) = k
              maxcldfrac = 0.0
            endif
            maxcldfrac = MAX(maxcldfrac, cf(i,j,k))
          endif

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

  subroutine compute_isccp_clds2(pfull, nclds, ktop, cldamt, high_ca, mid_ca, low_ca)
    real,     dimension(:,:,:), intent(in)  :: pfull, cldamt
    integer,  dimension(:,:),   intent(in)  :: nclds
    integer,  dimension(:,:,:), intent(in)  :: ktop
    real,     dimension(:,:),   intent(out) :: high_ca, mid_ca, low_ca
    real,     parameter :: mid_btm = 6.8e4, high_btm = 4.4e4
    ! local array
    integer :: i, j, k, k_top

    high_ca = 1.0
    mid_ca  = 1.0
    low_ca  = 1.0

    do i=1,size(cldamt,1)
      do j=1,size(cldamt,2)
        do k = 1,nclds(i,j)
          k_top = ktop(i,j,k)
          if (pfull(i,j,k_top)>mid_btm) then
            low_ca(i,j) = low_ca(i,j) * (1.0 - cldamt(i,j,k))
          else if (pfull(i,j,k_top)<high_btm) then
            high_ca(i,j) = high_ca(i,j) * (1.0 - cldamt(i,j,k))
          else
            mid_ca(i,j) = mid_ca(i,j) * (1.0 - cldamt(i,j,k))
          end if
        enddo
      enddo
    enddo

    low_ca  = (1.0 - low_ca)  * 1.0e2 ! unit percent
    mid_ca  = (1.0 - mid_ca)  * 1.0e2
    high_ca = (1.0 - high_ca) * 1.0e2

  end subroutine compute_isccp_clds2


  subroutine cldovrlap(pint, cld, nmxrgn, pmxrgn)
  !subroutine cldovrlap(lchnk   ,ncol    ,pint    ,cld     ,nmxrgn  ,pmxrgn  )
    !  This code is borrowed from CESM.
    !-----------------------------------------------------------------------
    ! Purpose: 
    ! Partitions each column into regions with clouds in neighboring layers.
    ! This information is used to implement maximum overlap in these regions
    ! with random overlap between them.
    ! On output,
    !    nmxrgn contains the number of regions in each column
    !    pmxrgn contains the interface pressures for the lower boundaries of
    !           each region! 
    ! Author: W. Collins
    !-----------------------------------------------------------------------
    !
    ! Input arguments
    !
    !integer, intent(in) :: ncol                ! number of atmospheric columns
    real, intent(in), dimension(:,:,:) :: pint  ! Interface pressure
    real, intent(in), dimension(:,:,:) :: cld   ! Fractional cloud cover
    !
    ! Output arguments
    !
    ! Number of maximally overlapped regions
    integer,  intent(out), dimension(size(cld,1), size(cld,2)) :: nmxrgn 
    real, intent(out), dimension(size(pint,1), size(pint,2), size(pint,3)) :: pmxrgn
    ! Maximum values of pressure for each maximally overlapped region.
    !    0->pmxrgn(i,1) is range of pressure for
    !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
    !    2nd region, etc
    !
    !---------------------------Local variables-----------------------------
    !
    integer :: i, j  ! Lat/Longitude index
    integer :: k     ! Level index
    integer :: n     ! Max-overlap region counter

    real, dimension(size(pint,1), size(pint,2), size(pint,3)) :: pnm  ! Interface pressure

    logical :: cld_found                          ! Flag for detection of cloud
    logical, dimension(size(cld,3)) :: cld_layer  ! Flag for cloud in layer
    !
    !------------------------------------------------------------------------
    !
    integer :: pver, pverp

    pver = size(cld,3)
    pverp = pver + 1
    
    do i = 1,size(cld,1)
      do j = 1,size(cld,2)
        cld_found = .false.
        cld_layer(:) = cld(i,j,:) > 0.0
        pmxrgn(i,j,:) = 0.0
        pnm(i,j,:) = pint(i,j,:) !* 10.0  ! why multiplied by 10?
        n = 1
        do k = 1, pver
          if (cld_layer(k) .and.  .not. cld_found) then
              cld_found = .true.
          else if ( .not. cld_layer(k) .and. cld_found) then
            cld_found = .false.
            if (count(cld_layer(k:pver)) == 0) then
              exit
            endif
            pmxrgn(i,j,n) = pnm(i,j,k)
            n = n + 1
          endif
        end do
        pmxrgn(i,j,n) = pnm(i,j,pverp)
        nmxrgn(i,j) = n
      end do
    end do
  end subroutine cldovrlap

  subroutine cloud_cover_diags_cam(cld, pmid, pint, Time)
    !!! pmid is pfull, and pint is phalf
    real,    intent(in), dimension(:,:,:) :: cld, pmid, pint
    ! Total, low, middle and high random overlap cloud cover
    real, dimension(size(cld,1), size(cld,2)) :: cldtot, cldlow, cldmed, cldhgh 
    
    type(time_type),  intent(in)  :: Time
    real :: used

    !---------------------------Local workspace-----------------------------
    integer :: i,j,k       ! lat/lon,level indices
    integer, dimension(size(cld,1), size(cld,2)) :: irgn    ! Max-overlap region index
    integer :: max_nmxrgn  ! maximum value of nmxrgn over columns
    integer :: ityp        ! Type counter
    real, dimension(size(cld,1), size(cld,2)) :: clrsky      ! Max-random clear sky fraction
    real, dimension(size(cld,1), size(cld,2)) :: clrskymax   ! Maximum overlap clear sky fraction
    
    !------------------------------Parameters-------------------------------
    real, parameter :: plowmax=1.2e5, plowmin=7.0e4  ! Max/min prs for low cloud cover range
    real, parameter :: pmedmax=7.0e4, pmedmin=4.0e4  ! Max/min prs for mid cloud cover range
    real, parameter :: phghmax=4.0e4, phghmin=5.0e3  ! Max/min prs for hgh cloud cover range

    real, dimension(4) :: ptypmin
    real, dimension(4) :: ptypmax
 
    data ptypmin /phghmin, plowmin, pmedmin, phghmin/
    data ptypmax /plowmax, plowmax, pmedmax, phghmax/

    integer, dimension(size(cld,1), size(cld,2)) :: nmxrgn 
    real, dimension(size(pint,1), size(pint,2), size(pint,3)) :: pmxrgn

    ! call the overlap subroutine to obtain the nmxrgn and pmxrgn
    call cldovrlap(pint, cld, nmxrgn, pmxrgn)

    ! Initialize region number
    max_nmxrgn = -1
    do i=1,size(cld,1)
      do j=1,size(cld,2)
        max_nmxrgn = max(max_nmxrgn, nmxrgn(i,j))
      end do
    end do

    do ityp= 1,4
      irgn = 1
      do k =1,max_nmxrgn-1
          do i=1,size(cld,1)
            do j=1,size(cld,2)
              if (pmxrgn(i,j,irgn(i,j)) < ptypmin(ityp) .and. irgn(i,j) < nmxrgn(i,j)) then
                  irgn(i,j) = irgn(i,j) + 1
              end if
            end do
          end do
      end do
      !
      ! Compute cloud amount by estimating clear-sky amounts
      !
      clrsky = 1.0
      clrskymax = 1.0
      do k=1,size(cld,3)
        do i=1,size(cld,1)
          do j=1,size(cld,2)
            if (pmid(i,j,k) >= ptypmin(ityp) .and. pmid(i,j,k) <= ptypmax(ityp)) then
                if (pmxrgn(i,j,irgn(i,j)) < pmid(i,j,k) .and. irgn(i,j) < nmxrgn(i,j)) then
                  irgn(i,j) = irgn(i,j) + 1
                  clrsky(i,j) = clrsky(i,j) * clrskymax(i,j)
                  clrskymax(i,j) = 1.0
                endif
                clrskymax(i,j) = min(clrskymax(i,j), 1.0-cld(i,j,k))
            endif
          end do
        end do
      end do

      if (ityp == 1) cldtot = 1.0 - (clrsky * clrskymax)
      if (ityp == 2) cldlow = 1.0 - (clrsky * clrskymax)
      if (ityp == 3) cldmed = 1.0 - (clrsky * clrskymax)
      if (ityp == 4) cldhgh = 1.0 - (clrsky * clrskymax)
    end do

    ! Write the output diagnostics
    if ( id_tot_cld_amt_cam > 0 ) then
      used = send_data ( id_tot_cld_amt_cam, cldtot*1e2, Time)
    endif
    if ( id_high_cld_amt_cam > 0 ) then
      used = send_data ( id_high_cld_amt_cam, cldlow*1e2, Time)
    endif
    if ( id_mid_cld_amt_cam > 0 ) then
      used = send_data ( id_mid_cld_amt_cam, cldmed*1e2, Time)
    endif
    if ( id_low_cld_amt_cam > 0 ) then
      used = send_data ( id_low_cld_amt_cam, cldhgh*1e2, Time)
    endif

  end subroutine cloud_cover_diags_cam

  subroutine diag_cldamt_maxrnd_overlap(cf, p_full, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full
    type(time_type),   intent(in)       :: Time
    real, dimension(size(cf,1), size(cf,2)) :: tca, high_ca, mid_ca, low_ca
    integer :: i, j, ks, ke !, ks_mid, ke_mid
    logical, dimension(size(cf,3)) :: ind_mid
    real :: mid_btm = 7e4, high_btm = 4e4

    tca = 1.0
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    do i=1,size(cf,1)
      do j=1,size(cf,2)
        ! total cloud amount
        !ks = 1
        !ke = size(cf,3)
        !call max_rnd_overlap_single_lev(cf(i,j,ks:ke), p_full(i,j,ks:ke), tca(i,j))

        ! high cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        call max_rnd_overlap_single_lev(cf(i,j,ks:ke), p_full(i,j,ks:ke), high_ca(i,j))
        !ks_mid = ke + 1

        ! low cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        ke = size(cf,3) !maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        call max_rnd_overlap_single_lev(cf(i,j,ks:ke), p_full(i,j,ks:ke), low_ca(i,j))
        !ke_mid = ks - 1

        ! middle cloud amount
        !ks = ks_mid
        !ke = ke_mid
        ind_mid = high_btm<=p_full(i,j,:) .and. p_full(i,j,:)<=mid_btm
        ks = minloc(p_full(i,j,:), 1, mask=ind_mid)
        ke = maxloc(p_full(i,j,:), 1, mask=ind_mid)
        call max_rnd_overlap_single_lev(cf(i,j,ks:ke), p_full(i,j,ks:ke), mid_ca(i,j))
      enddo
    enddo
    tca = 1.0 - (1.0-high_ca)*(1.0-mid_ca)*(1.0-low_ca)

    ! Diagnostics output
    call output_cldamt_max_random(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cldamt_maxrnd_overlap

  subroutine diag_cldamt_max_overlap(cf, p_full, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full
    type(time_type),   intent(in)       :: Time
    real, dimension(size(cf,1), size(cf,2)) :: tca, high_ca, mid_ca, low_ca
    integer :: i, j, ks, ke
    logical, dimension(size(cf,3)) :: ind_mid
    real :: mid_btm = 7e4, high_btm = 4e4
    tca = 1.0
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    ! total cld amount
    tca = maxval(cf, 3)

    do i=1,size(cf,1)
      do j=1,size(cf,2)
        ! high cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        high_ca(i,j) = maxval(cf(i,j,ks:ke), 1)

        ! low cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        low_ca(i,j) = maxval(cf(i,j,ks:ke), 1)

        ! middle cloud amount
        ind_mid = high_btm<=p_full(i,j,:) .and. p_full(i,j,:)<=mid_btm
        ks = minloc(p_full(i,j,:), 1, mask=ind_mid)
        ke = maxloc(p_full(i,j,:), 1, mask=ind_mid)
        mid_ca(i,j) = maxval(cf(i,j,ks:ke), 1)
      enddo
    enddo

    ! Diagnostics output
    call output_cldamt_max(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cldamt_max_overlap

  subroutine diag_cldamt_random_overlap(cf, p_full, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full
    type(time_type),   intent(in)       :: Time
    real, dimension(size(cf,1), size(cf,2)) :: tca, high_ca, mid_ca, low_ca
    integer :: i, j, ks, ke
    logical, dimension(size(cf,3)) :: ind_mid
    real :: mid_btm = 7e4, high_btm = 4e4

    tca = 1.0
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    do i=1,size(cf,1)
      do j=1,size(cf,2)
        ! total cloud amount
        !ks = 1
        !ke = size(cf,3)
        !call random_overlap_single_lev(cf(i,j,ks:ke), tca(i,j))

        ! high cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        call random_overlap_single_lev(cf(i,j,ks:ke), high_ca(i,j))

        ! low cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        call random_overlap_single_lev(cf(i,j,ks:ke), low_ca(i,j))

        ! middle cloud amount
        ind_mid = high_btm<=p_full(i,j,:) .and. p_full(i,j,:)<=mid_btm
        ks = minloc(p_full(i,j,:), 1, mask=ind_mid)
        ke = maxloc(p_full(i,j,:), 1, mask=ind_mid)
        call random_overlap_single_lev(cf(i,j,ks:ke), mid_ca(i,j))
      enddo
    enddo
    tca = 1.0 - (1.0-high_ca)*(1.0-mid_ca)*(1.0-low_ca)

    ! Diagnostics output
    call output_cldamt_random(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cldamt_random_overlap

  subroutine max_rnd_overlap_single_lev(cf, pfull, cldamt)
    implicit none
    real, dimension(:), intent(in) :: cf, pfull
    real, intent(out) :: cldamt
    ! local variables:
    integer, dimension(size(cf,1)) :: ktop, kbot
    real,    dimension(size(cf,1)) :: cldamt_cs
    integer :: nclds, kdim, k
    logical :: already_in_cloud, cloud_bottom_reached
    real    :: maxcldfrac

    kdim  = size(cf,1)
    nclds = 0
    ktop  = 1
    kbot  = 1
    maxcldfrac = 0.0
  
    ! set a flag indicating that we are searching for the next cloud top.
    already_in_cloud = .false.
    cloud_bottom_reached = .false.
    ! march down the column.
    do k=1,kdim
      ! find a layer containing cloud in the column.
      if (cf(k) .gt. cf_min) then
        if (.not. already_in_cloud) then
          nclds = nclds + 1
          already_in_cloud = .true.
          cloud_bottom_reached = .false.
          ktop(nclds) = k
          maxcldfrac = 0.0
        endif
        maxcldfrac = MAX(maxcldfrac, cf(k))
      endif

      if (cf(k) <= cf_min .and. already_in_cloud) then
        cloud_bottom_reached = .true.
        kbot(nclds) = k - 1
      else if (already_in_cloud .and. k == kdim) then
        cloud_bottom_reached = .true.
        kbot(nclds) = kdim
      endif

      if (cloud_bottom_reached) then
        cldamt_cs(nclds) = maxcldfrac
        already_in_cloud = .false.
        cloud_bottom_reached = .false.
      endif ! (cloud_bottom_reached)
    end do

    ! Random overlap
    cldamt = 1.0
    do k=1, nclds
      cldamt = cldamt * (1-cldamt_cs(k))
    end do
    cldamt = 1.0 - cldamt
  end subroutine max_rnd_overlap_single_lev

  subroutine random_overlap_single_lev(cf, cldamt)
    implicit none
    real, dimension(:), intent(in) :: cf
    real, intent(out) :: cldamt
    integer :: k

    ! Random overlap
    cldamt = 1.0
    do k=1,size(cf,1)
      cldamt = cldamt * (1-cf(k))
    end do
    cldamt = 1.0 - cldamt
  end subroutine random_overlap_single_lev

  subroutine output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, rhcrit, Time)
    real, intent(in), dimension(:,:,:) :: cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, rhcrit
    type(time_type) , intent(in)       :: Time
    real :: used

    if (id_cf > 0) then
      used = send_data (id_cf, cf, Time)
    endif
    if (id_reff_rad > 0) then
      used = send_data (id_reff_rad, reff_rad, Time)
    endif
    if (id_frac_liq > 0) then
      used = send_data (id_frac_liq, frac_liq, Time)
    endif
    if (id_qcl_rad > 0) then
      used = send_data (id_qcl_rad, qcl_rad, Time)
    endif
    if (id_rh_in_cf > 0) then
      used = send_data (id_rh_in_cf, rh_in_cf*100., Time)
    endif
    if (id_rhcrit > 0) then
      used = send_data (id_rhcrit, rhcrit*100.0, Time)
    endif
  end subroutine output_cloud_diags

  subroutine output_extra_diags_for_Park_ELF(Time, beta1, beta2, &
               alpha, eis, IS, DS, z700, zinv, Gamma700, Gamma_DL)

    real, intent(in), dimension(:,:) :: beta1, beta2, alpha, eis, &
                            IS, DS, z700, zinv, Gamma700, Gamma_DL
    type(time_type) , intent(in) :: Time
    real :: used

    if (id_eis>0) then
      ! Notice the eis here is a little different from that in calc_eis
      used = send_data (id_eis, eis, Time)
    endif
    if (id_beta1 > 0) then
      used = send_data (id_beta1, beta1, Time)
    endif
    if (id_beta2 > 0) then
      used = send_data (id_beta2, beta2, Time)
    endif
    if (id_alpha > 0) then
      used = send_data (id_alpha, alpha, Time)
    endif
    if (id_DS > 0) then
      used = send_data(id_DS, DS, Time)
    endif
    if (id_IS > 0) then
      used = send_data (id_IS, IS, Time)
    endif
    if (id_zinv > 0) then
      used = send_data (id_zinv, zinv, Time)
    endif
    if (id_z700 > 0) then
      used = send_data (id_z700, z700, Time)
    endif
    if (id_gamma_700 > 0) then
      used = send_data (id_gamma_700, Gamma700, Time)
    endif
    if (id_gamma_DL > 0) then
      used = send_data (id_gamma_DL, Gamma_DL, Time)
    endif
  end subroutine output_extra_diags_for_Park_ELF

  subroutine output_cldamt(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:) :: tca, high_ca, mid_ca, low_ca
    type(time_type),  intent(in)     :: Time
    real :: used

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
  end subroutine output_cldamt

  subroutine output_cldamt_max_random(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:) :: tca, high_ca, mid_ca, low_ca
    type(time_type),  intent(in)     :: Time
    real :: used

    if ( id_tot_cld_amt > 0 ) then
      used = send_data ( id_tot_cld_amt_mxr, tca*1e2, Time)
    endif
    if ( id_high_cld_amt > 0 ) then
      used = send_data ( id_high_cld_amt_mxr, high_ca*1e2, Time)
    endif
    if ( id_mid_cld_amt > 0 ) then
      used = send_data ( id_mid_cld_amt_mxr, mid_ca*1e2, Time)
    endif
    if ( id_low_cld_amt > 0 ) then
      used = send_data ( id_low_cld_amt_mxr, low_ca*1e2, Time)
    endif
  end subroutine output_cldamt_max_random

  subroutine output_cldamt_max(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:) :: tca, high_ca, mid_ca, low_ca
    type(time_type),  intent(in)     :: Time
    real :: used

    if ( id_tot_cld_amt_max > 0 ) then
      used = send_data ( id_tot_cld_amt_max, tca*1e2, Time)
    endif
    if ( id_high_cld_amt_max > 0 ) then
      used = send_data ( id_high_cld_amt_max, high_ca*1e2, Time)
    endif
    if ( id_mid_cld_amt_max > 0 ) then
      used = send_data ( id_mid_cld_amt_max, mid_ca*1e2, Time)
    endif
    if ( id_low_cld_amt_max > 0 ) then
      used = send_data ( id_low_cld_amt_max, low_ca*1e2, Time)
    endif
  end subroutine output_cldamt_max


  subroutine output_cldamt_random(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:) :: tca, high_ca, mid_ca, low_ca
    type(time_type),  intent(in)     :: Time
    real :: used

    if ( id_tot_cld_amt_rnd > 0 ) then
      used = send_data ( id_tot_cld_amt_rnd, tca*1e2, Time)
    endif
    if ( id_high_cld_amt_rnd > 0 ) then
      used = send_data ( id_high_cld_amt_rnd, high_ca*1e2, Time)
    endif
    if ( id_mid_cld_amt_rnd > 0 ) then
      used = send_data ( id_mid_cld_amt_rnd, mid_ca*1e2, Time)
    endif
    if ( id_low_cld_amt_rnd > 0 ) then
      used = send_data ( id_low_cld_amt_rnd, low_ca*1e2, Time)
    endif
  end subroutine output_cldamt_random

  subroutine cloud_simple_end()
  ! If alloocated are added in init then deallocate them here.
  end subroutine cloud_simple_end

end module cloud_simple_mod
