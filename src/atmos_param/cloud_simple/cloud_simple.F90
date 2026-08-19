module cloud_simple_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  use                fms_mod, only: stdlog, FATAL, WARNING, NOTE, error_mesg, &
                                    uppercase, check_nml_error
  use       time_manager_mod, only: time_type
  use     sat_vapor_pres_mod, only: compute_qs, lookup_es
  use       diag_manager_mod, only: register_diag_field, send_data
  use          constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, HLV, KAPPA, RADIUS, TFREEZE
  use                lcl_mod, only: lcl
  use  large_scale_cloud_mod, only: large_scale_cloud_init, large_scale_cloud_diag, &
                                    large_scale_cloud_end
  use marine_strat_cloud_mod, only: marine_strat_cloud_init, marine_strat_cloud_diag, &
                                    marine_strat_cloud_end
  use  cloud_cover_diags_mod, only: cloud_cover_diags_init, cloud_cover_diags, &
                                    cloud_cover_diags_end

  implicit none

  character(len=14), parameter :: mod_name_cld = "cloud_simple"

  logical :: do_init = .true.  ! Check if init needs to be run
  logical :: do_qcl_with_temp = .true.
  logical :: do_cloud_cover_diags = .true.
  logical :: do_add_stratocumulus = .false.

  ! Parameters to control the liquid cloud fraction
  real :: T_max = -5    ! Units in Celcius
  real :: T_min = -40   ! Units in Celcius

  ! For effective cloud droplet radius
  real :: reff_liq = 14.0   ! micron
  real :: reff_ice = 25.0

  ! For in-cloud liquid water
  real :: qcl_val = 0.2     ! g/kg

  ! ----- outputs for baisc cloud diagnostics ----- !
  integer :: id_cf, id_reff_rad, id_frac_liq, id_qcl_rad, id_rh_in_cf

  namelist /cloud_simple_nml/ &
            T_min, T_max, reff_liq, reff_ice, &
            qcl_val, do_qcl_with_temp, &
            do_add_stratocumulus, do_cloud_cover_diags

  contains

  subroutine cloud_simple_init (axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    integer :: io, ierr, nml_unit, stdlog_unit

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
    
    call error_mesg(mod_name_cld, 'Using SimCloud cloud scheme', NOTE)

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

    call large_scale_cloud_init(axes, Time)

    if (do_add_stratocumulus) call marine_strat_cloud_init(axes, Time)

    if (do_cloud_cover_diags) call cloud_cover_diags_init(axes, Time)

    do_init = .false.  !initialisation completed

  end subroutine cloud_simple_init

  ! ====================== Main Cloud Subroutine ====================== !
  subroutine cloud_simple(p_half, p_full, Time, temp, q_hum, z_full, &
                          wg_full, psg, temp_2m, q_2m, rh_2m, klcls, ocean, &
                          cf, reff_rad, qcl_rad)  ! outs
    real, intent(in),  dimension(:,:,:) :: temp, q_hum, p_full, p_half, z_full, wg_full
    real, intent(in),  dimension(:,:)   :: psg, temp_2m, q_2m, rh_2m
    integer, intent(in), dimension(:,:) :: klcls
    logical, intent(in), dimension(:,:) :: ocean
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: cf, reff_rad, qcl_rad
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: qs, frac_liq, rh_in_cf

    !check initiation has been done - ie read in parameters
    if (do_init) call error_mesg ('cloud_simple',  &
         'cloud_simple_init has not been called.', FATAL)

    ! Get the saturated specific humidity TOTAL (ie ice and vap) ***double check maths!
    call compute_qs(temp, p_full, qs)
    rh_in_cf = q_hum / qs

    call calc_liq_frac(temp, frac_liq)
    
    call calc_reff(frac_liq, reff_rad)

    call large_scale_cloud_diag(p_full, psg, rh_in_cf, q_hum, qs, qcl_rad, wg_full, cf, Time)

    if (do_add_stratocumulus) then
      call marine_strat_cloud_diag(temp, p_full, p_half, z_full, rh_in_cf, q_hum, &
              temp_2m, q_2m, rh_2m, psg, wg_full, klcls, cf, Time, ocean)
    end if

    call calc_qcl_rad(p_full, temp, cf, qcl_rad)

    if (do_cloud_cover_diags) call cloud_cover_diags(cf, p_full, p_half, Time)

    call output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, Time)

  end subroutine cloud_simple

  subroutine calc_liq_frac(temp, frac_liq)
    ! All liquid if above T_max and all ice below T_min,
    ! linearly interpolate between T_min and T_max
    real, intent(in),  dimension(:,:,:) :: temp
    real, intent(out), dimension(:,:,:) :: frac_liq

    where (temp > TFREEZE + T_max)
        frac_liq = 1.0
    elsewhere (temp < TFREEZE + T_min)
        frac_liq = 0.0
    elsewhere
        frac_liq = (temp - TFREEZE - T_min) / (T_max - T_min)
    end where

  end subroutine calc_liq_frac

  subroutine calc_reff(frac_liq, reff_rad)
    real, intent(in),  dimension(:,:,:) :: frac_liq
    real, intent(out), dimension(:,:,:) :: reff_rad

    reff_rad =  reff_liq * frac_liq + reff_ice * (1.0 - frac_liq)

  end subroutine calc_reff

  subroutine calc_qcl_rad(p_full, temp, cf, qcl_rad)
    ! calculate cloud water content
    real, intent(in),  dimension(:,:,:) :: p_full, cf, temp
    real, intent(out), dimension(:,:,:) :: qcl_rad
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: in_cloud_qcl

    if (do_qcl_with_temp) then
      in_cloud_qcl =  qcl_val * (temp-220.0) / (280.0-220.0)
      in_cloud_qcl = MAX(3.0e-4, MIN(qcl_val, in_cloud_qcl)) / 1.0e3 ! convert to kg/kg
    else
      ! in_cloud_qcl as a function of height
      in_cloud_qcl = 3.0e-4 + (1.0-3.0e-4) * (p_full-2.0e4) / 8.0e4
      in_cloud_qcl = MAX(0.0, in_cloud_qcl/1.0e3) ! convert to kg/kg
    end if

    qcl_rad = cf * in_cloud_qcl

  end subroutine calc_qcl_rad

  subroutine output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, Time)
    real, intent(in), dimension(:,:,:) :: cf, reff_rad, frac_liq, qcl_rad, rh_in_cf
    type(time_type) , intent(in)       :: Time
    logical :: used

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

  end subroutine output_cloud_diags

  subroutine cloud_simple_end()

    call large_scale_cloud_end()
    if (do_add_stratocumulus)  call marine_strat_cloud_end()
    if (do_cloud_cover_diags)  call cloud_cover_diags_end()

  end subroutine cloud_simple_end

end module cloud_simple_mod
