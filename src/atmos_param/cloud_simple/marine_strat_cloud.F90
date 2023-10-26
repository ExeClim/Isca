module marine_strat_cloud_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  use             fms_mod, only: stdlog, FATAL, WARNING, NOTE, error_mesg, &
                                 uppercase, check_nml_error
  use    time_manager_mod, only: time_type
  use  sat_vapor_pres_mod, only: compute_qs, lookup_es
  use    diag_manager_mod, only: register_diag_field, send_data
  use       constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, HLV, KAPPA, RADIUS, TFREEZE
  use             lcl_mod, only: lcl

  implicit none

  character(len=14), parameter :: mod_name = "strat_cloud"

  character(len=32) :: sc_diag_method = 'Park_ELF'
  logical :: intermediate_outputs_diags = .false.
  real :: dthdp_min_threshold = -0.05   ! K/hPa, which is -0.125 in CESM1.2.1

  ! ----- outputs for EIS, ECTEI and ELF diagnostics ----- !
  integer :: id_theta, id_dthdp, id_lts, id_eis, id_ectei, id_zlcl,       &
             id_gamma_850, id_gamma_DL, id_gamma_700, id_z700,            &
             id_zinv, id_ELF, id_beta1, id_beta2, id_IS, id_DS, id_alpha, &
             id_low_cld_amt_park, id_marine_strat

  ! Define constants for Earth mass and Newtonian gravational constant
  ! Refer to: https://github.com/Unidata/MetPy/ --> src/metpy/constants.py
  real :: EARTH_MASS = 5.9722e24  ! kg
  ! Refer to: https://physics.nist.gov/cgi-bin/cuu/Value?bg
  real :: GRAV_CONST = 6.674e-11  ! m^3 / kg / s^2

  ! Linear coefficient for Park_ELF scheme
  real :: park_a = 1.272
  real :: park_b = -0.366

  namelist /marine_strat_cloud_nml/ &
            sc_diag_method, intermediate_outputs_diags, dthdp_min_threshold, &
            park_a, park_b

  contains

  subroutine marine_strat_cloud_init(axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    integer :: io, ierr, nml_unit, stdlog_unit
    character(len=32) :: method_str = ''

#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=marine_strat_cloud_nml, iostat=io)
    ierr = check_nml_error(io, 'marine_strat_cloud_nml')
#else
    if (file_exist('input.nml')) then
      nml_unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
          read(nml_unit, nml=marine_strat_cloud_nml, iostat=io, end=10)
          ierr = check_nml_error(io, 'marine_strat_cloud_nml')
      enddo
10    call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, marine_strat_cloud_nml)

    call error_mesg(mod_name, 'The stratomulus diagnosis method is '// &
                    uppercase(trim(sc_diag_method)), NOTE)

    method_str = uppercase(trim(sc_diag_method))

    if (method_str(1:3)=='EIS' .or. method_str(1:5)=='ECTEI' &
        .or. method_str(1:4)=='PARK') then
      id_eis = register_diag_field (mod_name, 'eis', axes(1:2), Time, &
                  'estimated inversion strength', 'K')
    end if

    if (method_str(1:5)=='ECTEI' .or. method_str(1:4)=='PARK') then
      id_ectei = register_diag_field (mod_name, 'ectei', axes(1:2), Time, &
                      'estimated cloud top entrainment index', 'K')
    end if

    if (method_str(1:4)=='PARK') then
      id_ELF = register_diag_field (mod_name, 'ELF', axes(1:2), Time, &
                  'estimated low cloud fraction', '')
    end if

    id_marine_strat = register_diag_field ( mod_name, 'marine_strat', axes(1:3), Time, &
                        'marine low stratus cloud amount', '0-1' )

    id_zlcl = register_diag_field (mod_name, 'zlcl', axes(1:2), Time, &
                  'height of lcl', 'meter')
    id_theta = register_diag_field (mod_name, 'theta', axes(1:3), Time, &
                  'potential temperature', 'K')
    id_lts = register_diag_field (mod_name, 'lts', axes(1:2), Time, &
                  'low-tropospheric stability', 'K')

    if(intermediate_outputs_diags) then
      id_dthdp = register_diag_field (mod_name, 'dthdp', axes(1:3), Time, &
                        'dtheta/dp', 'K/hPa' )
      id_z700 = register_diag_field ( mod_name, 'z700', axes(1:2), Time, &
                        'height of 700mb', 'meter')

      if (method_str(1:3)=='EIS') then
        id_gamma_850 = register_diag_field (mod_name, 'gamma850', axes(1:2), Time, &
                          'moist lapse rate at 850hPa', 'K/m')
      end if

      if (method_str(1:4)=='PARK') then
        id_beta1  = register_diag_field (mod_name, 'beta1', axes(1:2), Time, &
                        'first low-level cloud suppression parameter', '')
        id_beta2  = register_diag_field (mod_name, 'beta2', axes(1:2), Time, &
                        'second low-level cloud suppression parameter', '')
        id_zinv   = register_diag_field (mod_name, 'zinv', axes(1:2), Time, &
                        'height of invesion layer', 'meter')
        id_DS     = register_diag_field (mod_name, 'DS', axes(1:2), Time, &
                        'decoupling strength', 'K')
        id_IS     = register_diag_field (mod_name, 'IS', axes(1:2), Time, &
                        'inversion strength', 'K')
        id_alpha  = register_diag_field (mod_name, 'alpha', axes(1:2), Time, &
                        'decoupling parameter', '')
        id_low_cld_amt_park = register_diag_field ( mod_name, 'low_cld_amt_park', axes(1:2), Time, &
                        'low cloud amount estimated from Park method', 'percent' )
        id_gamma_DL = register_diag_field (mod_name, 'gamma_DL', axes(1:2), Time, &
                        'moist lapse rate at decoupling layer', 'K/m')
        id_gamma_700 = register_diag_field (mod_name, 'gamma700', axes(1:2), Time, &
                        'moist lapse rate at 700hPa', 'K/m')
      end if
    end if

  end subroutine marine_strat_cloud_init

  subroutine marine_strat_cloud_diag(temp, p_full, p_half, z_full, rh, q_hum, temp_2m, &
                                      q_2m, rh_2m, psg, wg_full, klcls, cf, Time, ocean)
    implicit none
    real, intent(in),  dimension(:,:,:) :: temp, q_hum, p_full, p_half, z_full, rh, wg_full
    type(time_type),   intent(in)       :: Time
    real, intent(in),  dimension(:,:)   :: temp_2m, q_2m, rh_2m, psg
    integer, intent(in), dimension(:,:) :: klcls
    logical, intent(in), dimension(:,:) :: ocean
    real, intent(out), dimension(:,:,:) :: cf
    
    ! local variables
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: theta, dthdp, marine_strat
    integer, dimension(size(temp,1), size(temp,2)) :: kdthdp, kinvs
    real,    dimension(size(temp,1), size(temp,2)) :: eis, ectei, ELF, low_ca_park
    real :: strat, omega_pos_threshold
    logical :: used
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
      call calc_Park_proxies(p_full, psg, z_full, temp, temp_2m, q_hum, &
                             q_2m, rh_2m, klcls, ELF, kinvs, Time)
    end if

    k_surf = size(temp, 3)
    omega_pos_threshold = 0. !1.4*100/3600
    marine_strat = 0.0

    do i=1, size(temp, 1)
      do j=1, size(temp, 2)
        if (ocean(i,j)) then
          ! =========== Add off-coast marine stratiform clouds =========== !
          kk = kdthdp(i,j)

          if (kk .ne. 0) then
            kb = min(kk+1, k_surf)
            do k = kk, kb
              if (wg_full(i,j,k)>omega_pos_threshold .and. &
                  dthdp(i,j,k)<dthdp_min_threshold .and. p_full(i,j,k)>8.0e4) then
                call estimate_stratiform_cld(method_str, i, j, k, kb, p_full, &
                                        cf, rh, theta, eis, dthdp, ectei, ELF)
                marine_strat(i,j,k) = min(1.0, max(0.0, cf(i,j,k)))
              end if
            end do
          endif
        end if
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
  end subroutine marine_strat_cloud_diag

  subroutine estimate_stratiform_cld(method_str, i, j, k, kb, pfull, & 
                                cf, rh, theta, eis, dthdp, ectei, ELF)
    implicit none
    integer, intent(in) :: i, j, k
    integer, intent(in) :: kb
    character(len=32), intent(in) :: method_str
    real, intent(in),  dimension(:,:,:) :: rh, theta, pfull, dthdp
    real, intent(in),  dimension(:,:)   :: eis, ectei, ELF
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

    else if(method_str == 'EIS_WOOD') then
      !strat = min(1.0, max(0.0, 0.0221*eis(i,j) + 0.1128))
      strat = min(1.0, max(0.0, 0.06*eis(i,j) + 0.14)) ! Wood and Betherton, 2006
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))

    else if(method_str == 'ECTEI') then
      ! Kawai, Koshiro and Webb, 2017
      strat = min(1.0, max(0.0, 0.031*ectei(i,j) + 0.39))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
  
    else if(method_str == 'PARK_ELF') then
      ! Park and Shin, 2019, ACP
      ! strat = min(1.0, max(0.0, 1.272*ELF(i,j)-0.366))
      strat = min(1.0, max(0.0, park_a * ELF(i,j) + park_b))
      cf(i,j,k) = min(1.0, max(cf(i,j,k), strat))
   
    else
      call error_mesg('cloud_simple', method_str//' is not supported yet!', FATAL)
   
    end if
  end subroutine estimate_stratiform_cld

  subroutine calc_theta_dthdp(temp, temp_2m, pfull, phalf, ps, theta, dthdp, kdthdp)
    real,    intent(in),  dimension(:,:,:) :: temp, pfull, phalf
    real,    intent(in),  dimension(:,:)   :: temp_2m, ps
    real,    intent(out), dimension(:,:,:) :: theta, dthdp
    integer, intent(out), dimension(:,:)   :: kdthdp
    real, dimension(size(temp,1), size(temp,2)) :: theta_0
    real :: premib, pstar
    integer :: i, j, k, kb

    kdthdp = 0
    premib = 8.0e4
    dthdp = 0.0
    pstar = 1.0e5

    kb = size(temp, 3)  !bottom level
    do k=1,kb
      theta(:,:,k) =  temp(:,:,k) * (pstar / pfull(:,:,k))**(RDGAS / CP_AIR)
    end do

    do k=1,kb-1
      dthdp(:,:,k) = (theta(:,:,k) - theta(:,:,k+1)) / (phalf(:,:,k) - phalf(:,:,k+1)) * 1.0e2
    end do

    theta_0 = temp_2m * (pstar / ps)**(RDGAS / CP_AIR)
    dthdp(:,:,kb) = (theta(:,:,kb) - theta_0) / (phalf(:,:,kb) - ps) * 1.0e2

    kdthdp = minloc(dthdp, dim=3, mask=(pfull>premib).and.(dthdp<dthdp_min_threshold))

  end subroutine calc_theta_dthdp

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
          ! Use the exact LCL formula from D. M. Romps [2017, JAS 74(12)]
          zlcls(i,j) = lcl(ps(i,j), ts(i,j), rh=rh_surf(i,j))
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
    ! Estimated inversion stability (EIS)
    ! Refer to: Wood and Bretherton, 2006, Journal of Climate
    implicit none
    real,    intent(in),  dimension(:,:,:) :: pfull, zfull, temp
    real,    intent(in),  dimension(:,:)   :: ts, ps
    integer, intent(in),  dimension(:,:)   :: klcls
    type(time_type),      intent(in)       :: Time
    real,    intent(out), dimension(:,:)   :: eis
    real, dimension(size(temp,1), size(temp,2)) :: zlcl, z700, Gamma850, LTS
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: zfull_height
    real    :: pstar, T850
    logical :: used
    integer :: k700, i, j

    zfull_height = geopotential_to_height(zfull*GRAV)
    pstar = 1.e5 ! Pa

    do i=1, size(temp,1)
      do j=1, size(temp,2)
        k700 = minloc(abs(pfull(i,j,:) - 7.0e4), 1)
        LTS(i,j) = temp(i,j,k700)*((pstar/pfull(i,j,k700))**(RDGAS/CP_AIR)) - &
                   ts(i,j)*(pstar/ps(i,j))**(RDGAS/CP_AIR)
        T850 = (temp(i,j,k700) + ts(i,j)) / 2.0
        call calc_moist_lapse_rate(T850, 8.5e4, Gamma850(i,j))
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
    real :: k_en, C_qgap, beta
    logical :: used

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

  subroutine calc_Park_proxies(pfull, ps, zfull, temp, ts, q_hum, q_surf, &
                                rh_surf, klcls, ELF, kinvs, Time)
    ! Refer to: Park and Shin, 2019, Atmospheric Chemistry and Physics
    ! Heuristic estimation of low-level cloud fraction over the globe 
    ! based on a decoupling parameterization
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
    real, dimension(size(temp,1), size(temp,2)) :: beta1, IS, DS, eis, ectei, alpha, f_para
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: zfull_height
    real :: pstar, delta_zs, theta_ML
    logical :: used
    integer :: k700, i, j

    delta_zs = 2750.0 ! meter, constant
    pstar = 1.0e5 ! Pa
    kappa = RDGAS / CP_AIR

    zfull_height = geopotential_to_height(zfull*GRAV)

    call calc_lcls(klcls, pfull=pfull, temp=temp, ts=ts, ps=ps, & 
              rh_surf=rh_surf, plcls=plcl, tlcls=tlcl, zlcls=zlcl)

    where(zlcl < 0)
      zlcl = 0.0
    end where

    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        k700 = minloc(abs(pfull(i,j,:) - 7.0e4), 1)
        z700(i,j) = zfull_height(i,j,k700)
        
        ! Mixed Layer is the LCL
        call calc_moist_lapse_rate(tlcl(i,j), plcl(i,j), Gamma_DL(i,j))
        call calc_moist_lapse_rate(temp(i,j,k700), pfull(i,j,k700), Gamma700(i,j))

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
        kinvs(i,j) = minloc(abs(zinv(i,j)-zfull_height(i,j,:)), 1)
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

      !Add some diagnostic ouputs
      call output_extra_diags_for_Park_ELF(Time, beta1, beta2, &
              alpha, eis, IS, DS, z700, zinv, Gamma700, Gamma_DL)
    end if
  end subroutine calc_Park_proxies

  subroutine calc_moist_lapse_rate(T, p, Gamma)
    real, intent(in)  :: T, p
    real, intent(out) :: Gamma
    real :: es, qs
    
    ! Eq(5) in the following paper:
    ! Wood & Bretherton (2006). On the relationship between stratiform low cloud
    ! cover and lower-tropospheric stability. Journal of climate, 19(24), 6425-6432.

    call lookup_es(T, es)
    qs = 0.622 * es / (p - es)
    Gamma = (GRAV/CP_AIR) * (1.0 - (1.0 + HLV*qs/RDGAS/T) / (1.0 + HLV**2 * qs/CP_AIR/RVGAS/T**2))
  end subroutine calc_moist_lapse_rate

  subroutine output_extra_diags_for_Park_ELF(Time, beta1, beta2, &
               alpha, eis, IS, DS, z700, zinv, Gamma700, Gamma_DL)

    real, intent(in), dimension(:,:) :: beta1, beta2, alpha, eis, &
                            IS, DS, z700, zinv, Gamma700, Gamma_DL
    type(time_type) , intent(in) :: Time
    logical :: used

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

  subroutine marine_strat_cloud_end()

  end subroutine marine_strat_cloud_end

end module marine_strat_cloud_mod
