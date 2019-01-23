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
  integer :: id_cf, id_reff_rad, id_frac_liq, id_qcl_rad, id_rh_in_cf, id_simple_rhcrit
  character(len=14), parameter ::   mod_name_cld = "cloud_simple"

  integer, parameter :: B_SPOOKIE = 1, B_SUNDQVIST = 2, B_SMITH = 3
  integer, private :: cf_diag_formula = B_SPOOKIE
  character(len=32) :: cf_diag_formula_name = 'spookie'
  ! character(len=32) :: diag_cf_method = 'simple_rhcrit' ! 'read_rhcrit', 'linear', 'quadratic'
  logical :: do_simple_rhcrit = .true.
  logical :: do_read_scm_rhcrit = .false.
  logical :: do_linear_diag = .false.
  real, parameter :: FILL_VALUE = -999 ! Fill value for arrays
  real, dimension(100) :: scm_rhcrit = FILL_VALUE   ! Input array for single column critical RH. Max number of levels = 100
  real, dimension(100,2) :: scm_linear_coeffs = FILL_VALUE

  real    ::   simple_cca =  0.0
  real    ::   rhcsfc     = 0.95
  real    ::   rhc700     = 0.7
  real    ::   rhc200     = 0.3
  real    ::   rhmsfc     = 0.95
  real    ::   rhm700     = 0.7
  real    ::   rhm200     = 0.3

  namelist /cloud_simple_nml/ simple_cca, rhcsfc, rhc700, rhc200, &
                              rhmsfc, rhm700, rhm200, cf_diag_formula_name, &
                              do_simple_rhcrit, do_linear_diag, scm_linear_coeffs, &
                              do_read_scm_rhcrit, scm_rhcrit

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
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: qs, frac_liq, rh_in_cf, simple_rhcrit
    integer :: i, j, k, k_surf
    logical :: es_over_liq_and_ice

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

          call calc_qcl_rad(p_full(i,j,k), cf(i,j,k), qcl_rad(i,j,k) )
        end do
      end do
    end do

    !save some diagnotics
    call output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, simple_rhcrit, Time)
  end subroutine cloud_simple


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


  subroutine calc_qcl_rad(p_full, cf, qcl_rad)
    ! calculate cloud water content
    real , intent(in)   :: p_full, cf
    real , intent(out)  :: qcl_rad
    real :: in_cloud_qcl 

    in_cloud_qcl = 3.0e-4 + (1.0-3.0e-4) * (p_full-2.0e4) / 8.0e4
    in_cloud_qcl = MAX ( 0.0, in_cloud_qcl/1.0e3 ) ! convert to kg/kg
    qcl_rad = cf * in_cloud_qcl
  end subroutine calc_qcl_rad


  subroutine output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, simple_rhcrit, Time)
    real, intent(in), dimension(:,:,:) :: cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, simple_rhcrit
    type(time_type) , intent(in)       :: Time
    real :: used

    if ( id_cf > 0 ) then
      used = send_data ( id_cf, cf, Time)
    endif

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
