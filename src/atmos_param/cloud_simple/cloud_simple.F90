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
  use       constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, HLV, KAPPA, ES0, KELVIN
  use    betts_miller_mod, only: lcltabl

  implicit none

  logical ::   do_init = .true.  ! Check if init needs to be run

  real :: zerodegc = 273.15
  integer :: id_cf, id_reff_rad, id_frac_liq, id_qcl_rad, id_rh_in_cf, id_rhcrit
  integer :: id_conv_cf
  ! ----- outputs for EIS, ECTEI and ELF diagnostics ----- !
  integer :: id_theta, id_dthdp, id_lts, id_eis, id_ectei, id_zlcl, &
             id_gamma_850, id_gamma_DL, id_gamma_700, id_z700, &
             id_zinv, id_ELF, id_beta1, id_beta2, id_IS, id_DS, id_alpha, &
             id_low_cld_amt_park
  ! ----- outputs for cloud amount diagnostics ----- !
  integer :: id_tot_cld_amt, id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt

  character(len=14), parameter ::   mod_name_cld = "cloud_simple"

  integer, parameter :: B_SPOOKIE=1, B_SUNDQVIST=2, B_LINEAR=3, B_SMITH=4, B_SLINGO=5, B_XR96=6
  integer, private :: cf_diag_formula = B_LINEAR
  character(len=32) :: cf_diag_formula_name = 'linear'
  character(len=32) :: sc_diag_method = 'Park_ELF'
  logical :: do_simple_rhcrit = .false.
  logical :: do_read_scm_rhcrit = .false.
  logical :: do_qcl_with_temp = .false.
  logical :: do_cloud_amount_diags = .true.
  logical :: adjust_top = .true.
  logical :: do_add_stratocumulus = .false.
  logical :: intermediate_outputs_diags = .false.
  logical :: do_read_ts = .false.
  logical :: do_conv_cld = .false.
  logical :: do_adjust_low_cld = .false.
  real, parameter :: FILL_VALUE = -999.0 ! Fill value for arrays
  real, dimension(100) :: scm_rhcrit = FILL_VALUE   ! Input array for single column critical RH. Max number of levels = 100

  real :: simple_cca = 0.0
  real :: rhcsfc     = 0.95
  real :: rhc700     = 0.7
  real :: rhc200     = 0.3
  real :: cf_min     = 1e-10
  real :: dthdp_min  = -0.125
  real :: pshallow   = 7.5e4 ! copy from am4 diag_cloud.F90

  ! Parameters to control the coefficients profile of linear function of RH
  real :: a_surf     = 2.7
  real :: a_top      = 0.35
  real :: b_surf     = -2.2
  real :: b_top      = -0.07
  real :: nx         = 8

  ! For slingo80 scheme
  real :: slingo_rhc_low = 0.8
  real :: slingo_rhc_mid = 0.65
  real :: slingo_rhc_high = 0.8

  ! For low clout adjustment
  real :: omega_adj_threshold  = -0.1 !Pa/s, -3.6hPa/hour

  namelist /cloud_simple_nml/ simple_cca, rhcsfc, rhc700, rhc200, &
                              cf_diag_formula_name, &
                              do_read_scm_rhcrit, scm_rhcrit, &
                              do_qcl_with_temp, &
                              do_cloud_amount_diags, adjust_top, &
                              do_add_stratocumulus, sc_diag_method, &
                              intermediate_outputs_diags, &
                              dthdp_min, do_read_ts, &
                              a_surf, a_top, b_surf, b_top, nx, &
                              do_conv_cld, pshallow, cf_min, &
                              slingo_rhc_low, slingo_rhc_mid, slingo_rhc_high, &
                              do_adjust_low_cld, omega_adj_threshold

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
    character(len=32) :: tmp_str = ''

    call error_mesg('cloud_simple', 'The stratomulus diagnosis method is '// &
                    uppercase(trim(sc_diag_method)), NOTE)

    tmp_str = uppercase(trim(sc_diag_method))

    if (tmp_str(1:3)=='EIS' .or. tmp_str(1:5)=='ECTEI' .or. tmp_str(1:4)=='PARK') then
      id_eis = register_diag_field (mod_name_cld, 'eis', axes(1:2), Time, &
                      'estimated inversion strength', 'K')
    end if

    if (tmp_str(1:5)=='ECTEI' .or. tmp_str(1:4)=='PARK') then
      id_ectei = register_diag_field (mod_name_cld, 'ectei', axes(1:2), Time, &
                      'estimated cloud top entrainment index', 'K')
    end if

    if (tmp_str(1:4)=='PARK') then
      id_ELF    = register_diag_field (mod_name_cld, 'ELF', axes(1:2), Time, &
                      'estimated low cloud fraction', '')
    end if

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
      if (tmp_str(1:3)=='EIS') then
        id_gamma_850 = register_diag_field (mod_name_cld, 'gamma850', axes(1:2), Time, &
                        'moist lapse rate at 850hPa', 'K/m')
      end if
      if (tmp_str(1:4)=='PARK') then
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

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=cloud_simple_nml, iostat=io)
    ierr = check_nml_error(io, 'cloud_simple_nml')
#else
    if ( file_exist('input.nml') ) then
      nml_unit = open_namelist_file()
      ierr=1; do while (ierr /= 0)
         read (nml_unit, nml=cloud_simple_nml, iostat=io, end=10)
         ierr = check_nml_error(io, 'cloud_simple_nml')
      enddo
10    call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, cloud_simple_nml)

    ! Select cloud fraction diag formula
    if(uppercase(trim(cf_diag_formula_name)) == 'SPOOKIE') then
      cf_diag_formula = B_SPOOKIE
      call error_mesg('cloud_simple', 'Using default SPOOKIE cloud fraction diagnostic formula.', NOTE)
    else if(uppercase(trim(cf_diag_formula_name)) == 'SUNDQVIST') then
      cf_diag_formula = B_SUNDQVIST
      call error_mesg('cloud_simple', 'Using Sundqvist (1987) cloud fraction diagnostic formula.', NOTE)
    else if(uppercase(trim(cf_diag_formula_name)) == 'LINEAR') then
      cf_diag_formula = B_LINEAR
      call error_mesg('cloud_simple', 'Using linear cloud fraction diagnostic formula.', NOTE)
    else if(uppercase(trim(cf_diag_formula_name)) == 'SMITH') then
      cf_diag_formula = B_SMITH
      call error_mesg('cloud_simple', 'Using Smith (1990) cloud fraction diagnostic formula.', NOTE)
    else if(uppercase(trim(cf_diag_formula_name)) == 'SLINGO') then
      cf_diag_formula = B_SLINGO
      call error_mesg('cloud_simple', 'Using Slingo (1980) cloud fraction diagnostic formula.', NOTE)
    else if(uppercase(trim(cf_diag_formula_name)) == 'XR96') then
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
                          wg_full, psg, temp_2m, q_2m, precip, klcls, &
                          cf, reff_rad, qcl_rad)  ! outs
    real, intent(in),  dimension(:,:,:) :: temp, q_hum, p_full, p_half, z_full, wg_full
    real, intent(in),  dimension(:,:)   :: psg, temp_2m, precip, q_2m
    integer, intent(in), dimension(:,:) :: klcls
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
      call calc_convective_cf(p_full, precip, klcls, conv_cf, Time)
      !call calc_convective_cf2(p_full, precip, klcls, conv_cf, Time)
    end if

    ! rh_e is the effective RH
    rh_e = rh_in_cf * (1.0 - conv_cf)

    call calc_stratiform_cf(p_full, psg, rh_e, q_hum, qs, rhcrit, qcl_rad, cf)

    if(do_adjust_low_cld) then
      call adjust_low_cld(p_full, wg_full, cf)
    end if

    if (do_add_stratocumulus) then
      call add_stratiform_cld(temp, p_full, p_half, z_full, rh_e, q_hum, &
                              temp_2m, q_2m, psg, wg_full, klcls, cf, Time)
    end if

    if(do_conv_cld) then
      call merge_strat_conv_clouds(cf, conv_cf)
    end if

    call calc_qcl_rad(p_full, temp, cf, qcl_rad)

    if (do_cloud_amount_diags) then
      call diag_cloud_amount(cf, p_full, p_half, Time)
    end if

    call output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, rhcrit, Time)

  end subroutine cloud_simple

  subroutine adjust_low_cld(p_full, wg_full, cf)
    real, intent(in),    dimension(:,:,:) :: p_full, wg_full
    real, intent(inout), dimension(:,:,:) :: cf
    real :: premib, omega

    premib = 7.0e4 !Pa
    !omega_adj_threshold = -0.1   !Pa/s

    where(p_full>premib .and. omega_adj_threshold<wg_full .and. wg_full<0.0)
      cf = wg_full / omega_adj_threshold * cf
    elsewhere (p_full>premib .and. wg_full>=0.0)
      cf = 0
    elsewhere
      cf = cf
    end where
  end subroutine adjust_low_cld

  subroutine add_stratiform_cld(temp, p_full, p_half, z_full, rh, q_hum, temp_2m, q_2m, psg, wg_full, klcls, cf, Time)
    implicit none
    real, intent(in),  dimension(:,:,:) :: temp, q_hum, p_full, p_half, z_full, rh, wg_full
    type(time_type),   intent(in)       :: Time
    real, intent(in),  dimension(:,:)   :: temp_2m, q_2m, psg
    integer, intent(in), dimension(:,:) :: klcls
    real, intent(out), dimension(:,:,:) :: cf
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: theta, dthdp
    integer, dimension(size(temp,1), size(temp,2)) :: kdthdp
    real,    dimension(size(temp,1), size(temp,2)) :: eis, ectei, ELF, low_ca_park
    real :: strat, rhb_frac, used
    character(len=32) :: tmp_str = ''
    integer :: i, j, k, k700, kb, k_surf, kk

    k_surf = size(temp, 3)

    call calc_theta_dthdp(temp, temp_2m, p_full, p_half, psg, theta, dthdp, kdthdp)

    tmp_str = uppercase(trim(sc_diag_method))
    if (tmp_str(1:3)=='EIS' .or. tmp_str(1:5)=='ECTEI') then
      call calc_eis(p_full, z_full, temp, temp_2m, psg, klcls, eis, Time)
    end if
    if (tmp_str(1:5)=='ECTEI') then
      call calc_ectei(p_full, q_hum, q_2m, eis, ectei, Time)
    end if
    if (tmp_str(1:4)=='PARK') then
      call calc_Park_proxies(p_full, z_full, temp, temp_2m, q_hum, q_2m, klcls, ELF, Time)
    end if

    do i=1, size(temp, 1)
      do j=1, size(temp, 2)
        k = kdthdp(i,j)
        if(k.ne.0 .and. wg_full(i,j,k)>0) then
          kb = min(k+1, k_surf)
          if(uppercase(trim(sc_diag_method)) == 'CAM') then
            strat = min(1.0, max(0.0, (theta(i,j,k700) - theta(i,j,k_surf)) * 0.057 - 0.5573))
            cf(i,j,k) = max(strat, cf(i,j,k))
          end if
          if(uppercase(trim(sc_diag_method)) == 'SLINGO') then
            strat = min(1.0, max(0.0, -6.67*dthdp(i,j,k) - 0.667))
            rhb_frac = min(1.0, max(0.0, (rh(i,j,kb) - 0.6) / 0.2))
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
          end if
          if(uppercase(trim(sc_diag_method)) == 'SLINGO_NO_RH') then
            strat = min(1.0, max(0.0, -6.67*dthdp(i,j,k) - 0.667))
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
          end if
          if(uppercase(trim(sc_diag_method)) == 'DTHDP') then
            strat = min(1.0, max(0.0, -3.1196*dthdp(i,j,k) - 0.1246))
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
          end if
          if(uppercase(trim(sc_diag_method)) == 'EIS_WOOD') then
            !strat = min(1.0, max(0.0, 0.0221*eis(i,j) + 0.1128))
            strat = min(1.0, max(0.0, 0.06*eis(i,j) + 0.14)) ! Wood and Betherton, 2006
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
          end if
          if(uppercase(trim(sc_diag_method)) == 'EIS_WOOD_RH') then
            strat = min(1.0, max(0.0, 0.06*eis(i,j)+0.14)) !* (rh(i,j,kb)-0.6)/0.2)) ! Wood and Betherton, 2006
            rhb_frac = min(1.0, max(0.0, (rh(i,j,kb)-0.6) / 0.2))
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
            !cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
          end if
          if(uppercase(trim(sc_diag_method)) == 'EIS_RH') then
            !if (eis(i,j)>0.0) then
            strat = min(1.0, max(0.0, (0.092*eis(i,j)+ 0.027)*(2.078*rh(i,j,k)-6.45e-19)))
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
          end if
          if(uppercase(trim(sc_diag_method)) == 'ECTEI') then
            ! Kawai, Koshiro and Webb, 2017
            strat = min(1.0, max(0.0, 0.031*ectei(i,j) + 0.39))
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
          end if
          if(uppercase(trim(sc_diag_method)) == 'ECTEI_RH') then
            ! Kawai, Koshiro and Webb, 2017
            strat = min(1.0, max(0.0, 0.031*ectei(i,j) + 0.39))
            rhb_frac = min(1.0, max(0.0, (rh(i,j,kb)-0.6) / 0.2))
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat*rhb_frac))
          end if
          if(uppercase(trim(sc_diag_method)) == 'PARK_ELF') then
            ! Park and Shin, 2019, ACP
            strat = min(1.0, max(0.0, 0.86*ELF(i,j) + 0.02))
            cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
            !cf(i,j,k) = min(1.0, max(0.0, strat))
          end if
          low_ca_park(i,j) = min(1.0, max(0.0, 0.86*ELF(i,j) + 0.02))
        end if
      end do
    end do

    if (id_theta > 0) then
      used = send_data(id_theta, theta, Time)
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

  subroutine calc_theta_dthdp(temp, temp_2m, pfull, phalf, ps, theta, dthdp, kdthdp)
    real,    intent(in),  dimension(:,:,:) :: temp, pfull, phalf
    real,    intent(in),  dimension(:,:)   :: temp_2m, ps
    real,    intent(out), dimension(:,:,:) :: theta, dthdp
    integer, intent(out), dimension(:,:)   :: kdthdp
    real, dimension(size(temp,1), size(temp,2)) :: dthdp_min
    real :: premib
    integer :: i, j, k, kb

    dthdp_min = -0.125  ! d_theta / d_p, lapse rate
    kdthdp = 0
    premib = 7.0e4
    dthdp = 0.0
    kb = size(temp,3)

    do k=1,size(temp,3)
      theta(:,:,k) =  temp(:,:,k) * (ps / pfull(:,:,k))**(RDGAS / CP_AIR)
    end do

    do i=1, size(temp,1)
      do j=1, size(temp,2)
        do k=1, size(temp,3)-1
          dthdp(i,j,k) = (theta(i,j,k) - theta(i,j,k+1)) / (phalf(i,j,k) - phalf(i,j,k+1)) * 1.0e2
          if (phalf(i,j,k) >= premib) then
            if (dthdp(i,j,k) < dthdp_min(i,j)) then
              dthdp_min(i,j) = dthdp(i,j,k)
              kdthdp(i,j) = k     ! index of interface of max inversion
            end if
          end if
        end do

        ! Also check between the bottom layer and the surface
        ! Only perform this check if the criteria were not met above
        if(kdthdp(i,j) .eq. 0) then
          !dthdp(i,j,kb) = (theta(i,j,kb)-temp_2m(i,j)) / (phalf(i,j,kb)-1e5) * 1.0e2
          dthdp(i,j,kb) = (theta(i,j,kb)-temp_2m(i,j)) / (phalf(i,j,kb)-ps(i,j)) * 1.0e2
          if (dthdp(i,j,kb) < dthdp_min(i,j)) then
             dthdp_min(i,j) = dthdp(i,j,kb)
             kdthdp(i,j) = kb     ! index of interface of max inversion
          endif
        endif
      end do
    end do

  end subroutine calc_theta_dthdp


  subroutine calc_liq_frac(temp, frac_liq)
    ! All liquid if above zero and all ice below -40C
    ! linearly interpolate between T=0 and -40C
    real, intent(in),  dimension(:,:,:) :: temp
    real, intent(out), dimension(:,:,:) :: frac_liq

    where (temp > zerodegc)
        frac_liq = 1.0
    elsewhere (temp < zerodegc-40.0)
        frac_liq = 0.0
    elsewhere
        frac_liq = 1.0 - (zerodegc-temp) / 40.0
    end where

  end subroutine calc_liq_frac

  subroutine calc_reff(frac_liq, reff_rad)
    ! The effective cloud radius is bounded between 10 and 20 microns
    real, intent(in),  dimension(:,:,:) :: frac_liq
    real, intent(out), dimension(:,:,:) :: reff_rad

    reff_rad =  10.0 * frac_liq + 20.0 * (1.0 - frac_liq)  !units in microns
    !reff_rad =  10.0 * frac_liq + 30.0 * (1.0 - frac_liq)  !units in microns

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

        cf = min(1.0, cf)

      case(B_XR96)
        p_para = 0.25
        alpha_0 = 100.0
        gamma = 0.49

        where (rh.ge.1)
          cf = 1.0
        elsewhere
          cf = rh**p_para * (1.0 - EXP(-alpha_0*qcl_rad / (qsat-q_hum)**gamma))
        end where

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
      in_cloud_qcl = 0.2 * (temp-220.0) / (280.0-220.0)
      in_cloud_qcl = MAX(3.0e-4, MIN(0.2, in_cloud_qcl)) / 1.0e3 ! convert to kg/kg

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

  subroutine calc_convective_cf2(pfull, precip, klcls, conv_cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(in),  dimension(:,:)   :: precip
    integer, intent(in), dimension(:,:) :: klcls
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: conv_cf
    real, dimension(size(pfull,1), size(pfull,2)) :: precip_mm_per_day, plcl2d, conv_cf_tmp
    real    :: convcld_a, convcld_b, used
    integer :: i, j, k, klcl, ktop, nlayers

    call calc_lcls(klcls, pfull=pfull, plcls=plcl2d)

    precip_mm_per_day = precip * 24.0 * 3600.0  ! change units to mm/day

    convcld_a = -0.125 * log(0.14)  !0.246, ! 0.001
    convcld_b =  0.125 !0.0418 !

    conv_cf = 0.0
    conv_cf_tmp = convcld_a + convcld_b * log(1.0 + precip_mm_per_day)

    ktop = 3
    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        klcl = klcls(i,j)
        nlayers = klcl - ktop + 1
        conv_cf(i,j,ktop:klcl) = 1.0 - (1.0 - conv_cf_tmp(i,j))**(1.0 / nlayers)
      end do
    end do

    ! Output the diagnostics
    if (id_conv_cf > 0) then
      used = send_data(id_conv_cf, conv_cf, Time)
    endif

  end subroutine calc_convective_cf2

  subroutine merge_strat_conv_clouds(cf, conv_cf)
    implicit none
    real, intent(in),  dimension(:,:,:) :: conv_cf
    real, intent(out), dimension(:,:,:) :: cf

    !cf = max(cf, conv_cf)
    cf = 1.0 - (1.0 - cf) * (1.0 - conv_cf) ! Random overlap

  end subroutine merge_strat_conv_clouds

  subroutine calc_lcls(klcls, pfull, temp, zfull, plcls, tlcls, zlcls)
    ! Example to call:
    ! call calc_lcls(klcls, pfull=p_full, plcls=plcl2d)
    implicit none
    integer, intent(in), dimension(:,:) :: klcls
    real, intent(in),  dimension(:,:,:), optional :: temp, pfull, zfull
    real, intent(out), dimension(:,:),   optional :: plcls, tlcls, zlcls
    integer :: i, j

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

        if(.not.((present(pfull) .and. present(plcls)) .or. &
                 (present(temp)  .and. present(tlcls)) .or. &
                 (present(zfull) .and. present(zlcls)))) then
          call error_mesg('calc_lcls in cloud_simple', 'At least one group of '// &
                'pfull(plcls), temp(tlcls) and zfull(zlcls) should exist.', FATAL)
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
    real    :: pstar, T850, es, qs, used
    integer :: k700, i, j

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
        z700(i,j) = zfull(i,j,k700)
      end do
    end do

    call calc_lcls(klcls, zfull=zfull, zlcls=zlcl)
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

  subroutine calc_Park_proxies(pfull, zfull, temp, ts, q_hum, q_surf, klcls, ELF, Time)
    ! Refer to: Park and Shin, 2019, Atmospheric Chemistry and Physics
    ! https://www.atmos-chem-phys.net/19/5635/2019/

    implicit none
    real,    intent(in),  dimension(:,:,:) :: pfull, zfull, temp, q_hum
    real,    intent(in),  dimension(:,:)   :: ts, q_surf
    integer, intent(in),  dimension(:,:)   :: klcls
    type(time_type),      intent(in)       :: Time
    real,    intent(out), dimension(:,:)   :: ELF
    real, dimension(size(temp,1), size(temp,2)) :: plcl, tlcl, zlcl, z700, Gamma_DL, &
                                    Gamma700, LTS, z_ML, zinv, qv_ML, beta2
    ! other paramters
    real, dimension(size(temp,1), size(temp,2)) :: beta1, IS, DS, eis, ectei, alpha, f_para
    real :: pstar, delta_zs, theta_ML, used
    integer :: k700, i, j

    delta_zs = 2750.0 ! meter, constant
    pstar = 1.0e5 ! Pa
    kappa = RDGAS / CP_AIR

    call calc_lcls(klcls, pfull, temp, zfull, plcl, tlcl, zlcl)

    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        k700 = minloc(abs(pfull(i,j,:) - 7.0e4), 1)
        z700(i,j) = zfull(i,j,k700)
        ! Mixed Layer is the LCL
        call calc_moist_lapse_rate(tlcl(i,j), plcl(i,j), Gamma_DL(i,j))
        call calc_moist_lapse_rate(temp(i,j,k700), pfull(i,j,k700), Gamma700(i,j))
        theta_ML = tlcl(i,j) * (pstar / plcl(i,j))**kappa
        LTS(i,j) = temp(i,j,k700) * (pstar / pfull(i,j,k700))**kappa - theta_ML
        qv_ML(i,j) = q_hum(i,j,klcls(i,j))
      end do
    end do

    z_ML = zlcl
    zinv = -LTS/Gamma700 + z700 + delta_zs*(Gamma_DL/Gamma700)

    ! Rest zinv
    where(zinv<z_ML)
      zinv = z_ML
    end where
    where(zinv>z_ML+delta_zs)
      zinv = z_ML + delta_zs
    end where

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
      beta1 = (zinv+zlcl) / delta_zs
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
    call expand_cloud(nclds, ktop, kbot, cldamt, cloud)
    call compute_isccp_clds2(p_full, nclds, ktop, cldamt, high_ca, mid_ca, low_ca)
    !call compute_isccp_clds(p_full, cloud, high_ca, mid_ca, low_ca)

    ! Diagnostics output
    call output_cloud_amount(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cloud_amount

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

  subroutine output_cloud_amount(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:)   :: tca, high_ca, mid_ca, low_ca
    type(time_type) , intent(in)       :: Time
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
  end subroutine output_cloud_amount

  subroutine cloud_simple_end()
  ! If alloocated are added in init then deallocate them here.
  end subroutine cloud_simple_end

end module cloud_simple_mod
