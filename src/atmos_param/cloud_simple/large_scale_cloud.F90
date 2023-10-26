module large_scale_cloud_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  use            fms_mod, only: stdlog, FATAL, WARNING, NOTE, error_mesg, &
                                uppercase, check_nml_error
  use   time_manager_mod, only: time_type
  use sat_vapor_pres_mod, only: compute_qs, lookup_es
  use   diag_manager_mod, only: register_diag_field, send_data
  use      constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, HLV, KAPPA, RADIUS, TFREEZE

  implicit none

  character(len=14), parameter :: mod_name = "ls_cloud"

  integer, parameter :: B_SPOOKIE=1, B_SUNDQVIST=2, B_LINEAR=3, &
                        B_SMITH=4, B_SLINGO=5, B_XR96=6
  integer, private  :: cf_diag_formula = B_LINEAR
  character(len=32) :: cf_diag_formula_name = 'linear'

  logical :: do_simple_rhcrit = .false.
  logical :: do_fitted_rhcrit = .false.
  logical :: do_adjust_cld_by_omega = .false.
  logical :: do_poly_rhcrit = .false.
  logical :: do_freezedry = .false.

  real :: rhcsfc = 0.95
  real :: rhc700 = 0.7
  real :: rhc200 = 0.3

  ! Parameters to control the fitted critical RH profile
  real :: rhc_surf = 0.8
  real :: rhc_top  = 0.4
  real :: n_rhc    = 3.1

  ! Parameters for the freeze-dry problem in polar region
  real :: qv_polar_val    = 0.003  ! kg/kg
  real :: freezedry_power = 2.5

  ! Parameters to control linear coefficient profile
  real :: linear_a_surf = 42
  real :: linear_a_top  = 12
  real :: linear_power  = 8.5

  ! For slingo80 scheme
  real :: slingo_rhc_low  = 0.8
  real :: slingo_rhc_mid  = 0.65
  real :: slingo_rhc_high = 0.8

  ! For cloud adjustment by omega
  real :: omega_adj_threshold = 0.1 !Pa/s
  real :: adj_pres_threshold = 7.0e4 !Pa

  integer :: id_rhcrit

  namelist /large_scale_cloud_nml/ &
            rhcsfc, rhc700, rhc200, &
            do_fitted_rhcrit, do_poly_rhcrit, &
            rhc_surf, rhc_top, n_rhc, &
            cf_diag_formula_name, &
            linear_a_surf, linear_a_top, linear_power, &
            slingo_rhc_low, slingo_rhc_mid, slingo_rhc_high, &
            do_adjust_cld_by_omega, omega_adj_threshold, adj_pres_threshold, &
            do_freezedry, qv_polar_val, freezedry_power

  contains

  subroutine large_scale_cloud_init(axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    integer :: io, ierr, nml_unit, stdlog_unit
    character(len=32) :: method_str = ''

#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=large_scale_cloud_nml, iostat=io)
    ierr = check_nml_error(io, 'large_scale_cloud_nml')
#else
    if (file_exist('input.nml')) then
      nml_unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
        read(nml_unit, nml=large_scale_cloud_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'large_scale_cloud_nml')
      enddo
10    call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, large_scale_cloud_nml)

    ! Select cloud fraction diag formula
    method_str = uppercase(trim(cf_diag_formula_name))

    if(method_str == 'SPOOKIE') then
      cf_diag_formula = B_SPOOKIE
      call error_mesg(mod_name, 'Using default SPOOKIE cloud fraction diagnostic formula.', NOTE)

    else if(method_str == 'SUNDQVIST') then
      cf_diag_formula = B_SUNDQVIST
      call error_mesg(mod_name, 'Using Sundqvist (1989) cloud fraction diagnostic formula.', NOTE)

    else if(method_str == 'LINEAR') then
      cf_diag_formula = B_LINEAR
      call error_mesg(mod_name, 'Using linear cloud fraction diagnostic formula.', NOTE)

    else if(method_str == 'SMITH') then
      cf_diag_formula = B_SMITH
      call error_mesg(mod_name, 'Using Smith (1990) cloud fraction diagnostic formula.', NOTE)

    else if(method_str == 'SLINGO') then
      cf_diag_formula = B_SLINGO
      call error_mesg(mod_name, 'Using Slingo (1980) cloud fraction diagnostic formula.', NOTE)

    else if(method_str == 'XR96') then
      cf_diag_formula = B_XR96
      call error_mesg(mod_name, 'Using Xu and Krueger (1996) cloud fraction diagnostic formula.', NOTE)

    else
      call error_mesg(mod_name, '"'//trim(cf_diag_formula_name)//'"'// &
                ' is not a valid cloud fraction diagnostic formula.', FATAL)
    endif

    if (cf_diag_formula .eq. B_SPOOKIE .or. cf_diag_formula .eq. B_SUNDQVIST &
    .or. cf_diag_formula .eq. B_SMITH) then
      do_simple_rhcrit = .true.
    else
      do_simple_rhcrit = .false.
    end if

    if(do_simple_rhcrit) then
      id_rhcrit = register_diag_field (mod_name, 'rhcrit', axes(1:3), Time, &
                    'critical relative humidity', '%')
    end if

  end subroutine large_scale_cloud_init

  subroutine large_scale_cloud_diag(pfull, ps, rh, q_hum, qsat, qcl_rad, wg_full, cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull, rh, q_hum, qsat, qcl_rad, wg_full
    real, intent(in),  dimension(:,:)   :: ps
    real, intent(out), dimension(:,:,:) :: cf
    type(time_type), intent(in)         :: Time
    real, dimension(size(rh,1), size(rh,2), size(rh,3)) :: rhcrit
    logical :: used

    if (do_simple_rhcrit) then
      call calc_rhcrit(pfull, rhcrit)
    end if

    call calc_large_scale_cld_frac(pfull, ps, rh, q_hum, qsat, rhcrit, qcl_rad, cf)

    if(do_adjust_cld_by_omega) then
      call adjust_cld_by_omega(pfull, wg_full, cf)
    end if

    if(do_freezedry) then
      call freezedry_adjustment(pfull, ps, cf, q_hum)
    end if

    if (id_rhcrit > 0) then
      used = send_data (id_rhcrit, rhcrit*100.0, Time)
    endif

  end subroutine large_scale_cloud_diag

  subroutine calc_rhcrit(p_full, rhcrit)
    !get the RH needed as a threshold for the cloud fraction calc.
    real, intent(in),  dimension(:,:,:) :: p_full
    real, intent(out), dimension(:,:,:) :: rhcrit
    real :: p_surf
    real, dimension(size(p_full,1), size(p_full,2), size(p_full,3)) :: sigma

    real :: rhc1, rhc2, zrhc
    p_surf = 1e5

    if (do_fitted_rhcrit) then
        rhcrit = rhc_top + (rhc_surf - rhc_top) * EXP(1.0 - (p_surf/p_full)**n_rhc)
    else if (do_poly_rhcrit) then
        rhc1 = 0.8
        rhc2 = 1.73
        zrhc = 0.95 !0.85
        sigma = p_full / p_surf
        rhcrit = zrhc - rhc1*sigma * (1.0-sigma) * (1.0 + rhc2*(sigma-0.5))
        !rhcrit = poly_rhc_surf*sigma + poly_rhc_top*(1.0-sigma) + sigma * (1.0-sigma) * (1.0 + rhc2*(sigma-0.5))
    else
      where (p_full > 7.0e4)
        rhcrit = rhcsfc - (rhcsfc - rhc700) * (p_surf - p_full) / (p_surf - 7.0e4)
      elsewhere (p_full > 2.0e4)
        rhcrit = rhc700 - (rhc700 - rhc200) * (7.0e4 - p_full) / 5.0e4
      elsewhere
        rhcrit = rhc200
      endwhere
    end if

  end subroutine calc_rhcrit

  subroutine adjust_cld_by_omega(p_full, wg_full, cf)
    real, intent(in),    dimension(:,:,:) :: p_full, wg_full
    real, intent(inout), dimension(:,:,:) :: cf

    where (p_full>adj_pres_threshold .and. omega_adj_threshold>wg_full .and. wg_full>0.0)
     cf = min(1.0, (omega_adj_threshold-wg_full)/omega_adj_threshold) * cf
    end where

    where (p_full>adj_pres_threshold .and. wg_full>=omega_adj_threshold)
      cf = 0.0
    end where

  end subroutine adjust_cld_by_omega

  subroutine freezedry_adjustment(p_full, psg, cf, q_hum)
    real, intent(in),    dimension(:,:,:) :: p_full, q_hum
    real, intent(in),    dimension(:,:)   :: psg
    real, intent(inout), dimension(:,:,:) :: cf
    integer :: k
    real, dimension(size(cf,1), size(cf,2)) :: qv_k
  
    ! VAVRUS and WALISER, 2008, An Improved Parameterization for Simulating
    !   Arctic Cloud Amount in the CCSM3 Climate Model,
    !   Journal of Climate, DOI: 10.1175/2008JCLI2299.1
    ! 
    ! A difference from the paper is we adjust the clouds not only near the 
    ! surface but through all the levels of atmosphere
  
    do k=1,size(p_full,3)
        qv_k = (p_full(:,:,k) / psg)**freezedry_power * qv_polar_val
        cf(:,:,k) = cf(:,:,k) * MAX(0.15, MIN(1.0, q_hum(:,:,k) / qv_k))
    end do

  end subroutine freezedry_adjustment

  subroutine calc_large_scale_cld_frac(pfull, ps, rh, q_hum, qsat, rhcrit, qcl_rad, cf)
    ! Calculate large scale (stratiform) cloud fraction
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
        ! Refer to: Sundqvist et al., 1989, MWR, Condensation and Cloud Parameterization 
        ! Studies with a Mesoscale Numerical Weather Prediction Model 
        ! https://journals.ametsoc.org/doi/10.1175/1520-0493(1989)117%3C1641:CACPSW%3E2.0.CO%3B2
        !
        ! Eq (3.13) in the paper above

        cf = 1.0 - ((1.0 - MIN(rh,1.0)) / (1.0 - rhcrit))**0.5
      
      case(B_SMITH)
        ! Refer to: Smith (1990). A scheme for predicting layer clouds and their water 
        ! in a general circulation model. QJRMS, 116(492), 435-460.

        cf = 1.0 - (3.0 / sqrt(2.0) * (1.0 - MIN(rh,1.0))/(1.0 - rhcrit))**(2.0/3.0)

      case(B_SLINGO)
        ! Refer to: Slingo (1980). A cloud parametrization scheme derived from
        ! GATE data for use with a numerical model. QJRMS, 106(450), 747-770.

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
        ! Refer to: Xu & Randall (1996). A semiempirical cloudiness parameterization 
        ! for use in climate models. JAS, 53(21), 3084-3102.

        p_para = 0.25
        alpha_0 = 100.0
        gamma = 0.49

        where (rh.ge.1)
          cf = 1.0
        elsewhere (rh.le.0)
          cf = 0.0
        elsewhere
          cf = rh**p_para * (1.0 - EXP(-alpha_0*qcl_rad / (qsat-q_hum)**gamma))
        end where
      
      case(B_LINEAR)
        call calc_cf_linear(pfull, rh, ps, cf)

      case default
        call error_mesg(mod_name, 'invalid cloud fraction diagnostic formula', FATAL)
    end select

    cf = MAX(0.0, MIN(1.0, cf))

  end subroutine calc_large_scale_cld_frac

  subroutine calc_cf_linear(p_full, rh, ps, cf)
    ! Calculate LS (stratiform) cloud fraction as a linear function of RH
    real, intent(in),  dimension(:,:,:) :: p_full, rh
    real, intent(in),  dimension(:,:)   :: ps
    real, intent(out), dimension(:,:,:) :: cf
    real, dimension(size(p_full,1), size(p_full,2), size(p_full,3)) :: coeff_a
    integer :: k 

    do k=1,size(rh,3)
      coeff_a(:,:,k) = linear_a_top + (linear_a_surf-linear_a_top) * &
                        EXP(1.0 - (ps/p_full(:,:,k))**linear_power)
    end do

    cf = coeff_a * (rh - 1.0) + 1.0
    cf = MAX(0.0, MIN(1.0, cf))

  end subroutine calc_cf_linear

  subroutine large_scale_cloud_end()

  end subroutine large_scale_cloud_end

end module large_scale_cloud_mod
