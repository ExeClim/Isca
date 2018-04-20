MODULE socrates_interface_mod

  ! Socrates calculation interface module
  ! Takes FMS time, spectra, temperature, and pressure
  ! Outputs FMS heating rate, and downwards surface LW and SW
  ! MDH 30/01/18

  !----------

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  ! ExoFMS diagnostics
  USE  diag_manager_mod, ONLY: register_diag_field, send_data

  ! ExoFMS time
  USE time_manager_mod, ONLY: time_type, OPERATOR(+), OPERATOR(-), OPERATOR(/=)

  ! Socrates modules
  USE read_control_mod
  USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
  USE def_spectrum
  use constants_mod, only: grav, rdgas, rvgas, cp_air
  use fms_mod, only: stdlog

  IMPLICIT NONE


  ! Input spectra
  TYPE (StrSpecData) :: spectrum_calc
  TYPE (StrSpecData) :: spectrum_lw, spectrum_lw_hires
  TYPE (StrSpecData) :: spectrum_sw, spectrum_sw_hires

  ! Control options:
  TYPE(StrCtrl) :: control_calc
  TYPE(StrCtrl) :: control_sw, control_sw_hires
  TYPE(StrCtrl) :: control_lw, control_lw_hires

  ! Diagnostic IDs, name, and missing value
  INTEGER :: id_soc_spectral_olr
  INTEGER :: id_soc_olr, id_soc_olr_spectrum_lw, id_soc_surf_spectrum_sw
  INTEGER :: id_soc_heating_sw, id_soc_heating_lw, id_soc_heating_rate
  INTEGER :: id_soc_flux_up_lw, id_soc_flux_down_sw
  INTEGER :: n_soc_bands_lw, n_soc_bands_sw
  INTEGER :: n_soc_bands_lw_hires, n_soc_bands_sw_hires
  CHARACTER(len=10), PARAMETER :: soc_mod_name = 'socrates'
  REAL :: missing_value = -999

  ! Socrates inputs from namelist
  REAL :: stellar_constant = 1370.0
  LOGICAL :: tidally_locked = .TRUE.
  character(len=256) :: lw_spectral_filename='/scratch/sit204/sp_lw_ga7'
  character(len=256) :: lw_hires_spectral_filename='/scratch/sit204/sp_lw_ga7'
  character(len=256) :: sw_spectral_filename='/scratch/sit204/sp_sw_ga7'
  character(len=256) :: sw_hires_spectral_filename='/scratch/sit204/sp_sw_ga7'
    
  NAMELIST/socrates_rad_nml/ stellar_constant, tidally_locked, lw_spectral_filename, lw_hires_spectral_filename, &
                             sw_spectral_filename, sw_hires_spectral_filename



CONTAINS

  SUBROUTINE socrates_init(is, ie, js, je, num_levels, axes, Time, lat)
    !! Initialises Socrates spectra, arrays, and constants

    ! Arguments
    INTEGER, INTENT(in), DIMENSION(4) :: axes
    !! NB axes refers to the handles of the axes defined in fv_diagnostics
    TYPE(time_type), INTENT(in)       :: Time
    INTEGER, INTENT(in)               :: is, ie, js, je, num_levels
    REAL, INTENT(in) , DIMENSION(:,:)   :: lat
    
    integer :: io, stdlog_unit
    !-------------------------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=socrates_rad_nml, iostat=io)
#else
   if ( file_exist('input.nml') ) then
      nml_unit = open_namelist_file()
      read (nml_unit, socrates_rad_nml, iostat=io)
      call close_file(nml_unit)
   endif
#endif
stdlog_unit = stdlog()
write(stdlog_unit, socrates_rad_nml)

    ! Socrates spectral files -- should be set by namelist
    control_lw%spectral_file = lw_spectral_filename
    control_lw_hires%spectral_file = lw_hires_spectral_filename

    control_sw%spectral_file = sw_spectral_filename
    control_sw_hires%spectral_file = sw_hires_spectral_filename

    ! Read in spectral files
    CALL read_spectrum(control_lw%spectral_file,spectrum_lw)
    CALL read_spectrum(control_lw_hires%spectral_file,spectrum_lw_hires)
    CALL read_spectrum(control_sw%spectral_file,spectrum_sw)
    CALL read_spectrum(control_sw_hires%spectral_file,spectrum_sw_hires)
    
    ! Set Socrates configuration
    CALL read_control(control_lw,spectrum_lw)
    CALL read_control(control_lw_hires,spectrum_lw_hires)
    CALL read_control(control_sw,spectrum_sw)
    CALL read_control(control_sw_hires,spectrum_sw_hires)

    ! Specify LW and SW setups
    control_sw%isolir=1
    control_sw_hires%isolir=1
    control_lw%isolir=2
    control_lw_hires%isolir=2


    ! Register diagostic fields
    id_soc_spectral_olr = &
         register_diag_field ( soc_mod_name, 'soc_spectral_olr',(/ axes(1:2), axes(5)/) , Time, &
         'socrates substellar LW OLR spectrum', &
         'watts/m2', missing_value=missing_value          )

    id_soc_olr = &
         register_diag_field ( soc_mod_name, 'soc_olr', axes(1:2), Time, &
         'outgoing longwave radiation', &
         'watts/m2', missing_value=missing_value               )

    id_soc_olr_spectrum_lw = &
         register_diag_field ( soc_mod_name, 'soc_olr_spectrum_lw',(/ axes(1:2), axes(5)/) , Time, &
         'socrates substellar LW OLR spectrum', &
         'watts/m2', missing_value=missing_value               )

    id_soc_surf_spectrum_sw = &
         register_diag_field ( soc_mod_name, 'soc_surf_spectrum_sw',(/ axes(1:2), axes(5)/) , Time, &
         'socrates substellar SW surface spectrum', &
         'watts/m2', missing_value=missing_value               )

    id_soc_heating_lw = &
         register_diag_field ( soc_mod_name, 'soc_heating_lw', axes(1:3), Time, &
         'socrates LW heating rate', &
         'J/s', missing_value=missing_value               )

    id_soc_heating_sw = &
         register_diag_field ( soc_mod_name, 'soc_heating_sw', axes(1:3), Time, &
         'socrates SW heating rate', &
         'J/s', missing_value=missing_value               )

    id_soc_heating_rate = &
         register_diag_field ( soc_mod_name, 'soc_heating_rate', axes(1:3), Time, &
         'socrates total heating rate', &
         'J/s', missing_value=missing_value               )

    id_soc_flux_up_lw = &
         register_diag_field ( soc_mod_name, 'soc_flux_up_lw', axes(1:3), Time, &
         'socrates LW flux up', &
         'watts/m2', missing_value=missing_value               )

    id_soc_flux_down_sw = &
         register_diag_field ( soc_mod_name, 'soc_flux_down_sw', axes(1:3), Time, &
         'socrates SW flux down', &
         'watts/m2', missing_value=missing_value               )

    ! Number of bands
    n_soc_bands_lw = spectrum_lw%dim%nd_band
    n_soc_bands_lw_hires = spectrum_lw_hires%dim%nd_band
    n_soc_bands_sw = spectrum_sw%dim%nd_band
    n_soc_bands_sw_hires = spectrum_sw_hires%dim%nd_band

    ! Print Socrates init data from one processor
    IF (js == 1) THEN
       PRINT*, ' '
       PRINT*, '-----------------------------------'

       PRINT*,'  ____                       _            '
       PRINT*,' / ___|  ___   ___ _ __ __ _| |_ ___  ___ '
       PRINT*,' \___ \ / _ \ / __|  __/ _` | __/ _ \/ __|'
       PRINT*,'  ___) | (_) | (__| | | (_| | ||  __/\__ \'
       PRINT*,' |____/ \___/ \___|_|  \__,_|\__\___||___/'

       PRINT*, ' '

       PRINT*, 'Initialised Socrates v17.03'
       PRINT*, 'Stellar constant = ', stellar_constant
       PRINT*, 'Longwave spectral file = ', TRIM(control_lw%spectral_file), ' WITH ', n_soc_bands_lw, ' bands'
       PRINT*, 'Longwave hires spectral file = ', TRIM(control_lw_hires%spectral_file), ' WITH ', n_soc_bands_lw_hires, ' bands'
       PRINT*, 'Shortwave spectral file = ', TRIM(control_sw%spectral_file), ' WITH ', n_soc_bands_sw, ' bands'
       PRINT*, 'Shortwave hires spectral file = ', TRIM(control_sw_hires%spectral_file), ' WITH ', n_soc_bands_sw_hires, ' bands'
       PRINT*, ' '
       PRINT*, '-----------------------------------'
       PRINT*, ' '
    end if

    return
  end subroutine socrates_init
  ! ==================================================================================


  ! Set up the call to the Socrates radiation scheme
  ! -----------------------------------------------------------------------------
  subroutine socrates_interface(Time_diag, rlat, rlon, soc_mode, hires_mode,      &
       fms_temp, fms_t_surf, fms_p_full, fms_p_half, n_profile, n_layer,        &
       output_heating_rate, fms_net_surf_sw_down, fms_surf_lw_down, fms_stellar_flux )

    use realtype_rd
    use soc_constants_mod
    use read_control_mod
    use socrates_calc_mod
    use compress_spectrum_mod
    use def_spectrum
    use def_dimen,   only: StrDim
    use def_control, only: StrCtrl,  allocate_control,   deallocate_control
    use def_atm,     only: StrAtm,   allocate_atm,       deallocate_atm
    use def_cld,     only: StrCld,   allocate_cld,       deallocate_cld, &
         allocate_cld_prsc,  deallocate_cld_prsc, &
         allocate_cld_mcica, deallocate_cld_mcica
    use def_aer,     only: StrAer,   allocate_aer,       deallocate_aer, &
         allocate_aer_prsc,  deallocate_aer_prsc
    use def_bound,   only: StrBound, allocate_bound,     deallocate_bound
    use def_out,     only: StrOut,                       deallocate_out
    !-----------------------------------------------------------------------
    implicit none

    ! Input time
    type(time_type), intent(in)         :: Time_diag

    INTEGER(i_def), intent(in) :: n_profile, n_layer
    logical, intent(in) :: soc_mode, hires_mode
    INTEGER(i_def) :: nlat

    ! Input arrays
    real(r_def), intent(in) :: fms_temp(:,:,:)
    real(r_def), intent(in) :: fms_p_full(:,:,:)
    real(r_def), intent(in) :: fms_p_half(:,:,:)
    real(r_def), intent(in) :: fms_t_surf(:,:)
    real(r_def), intent(in) :: fms_stellar_flux(:,:)
    real(r_def), intent(in) :: rlon(:,:)
    real(r_def), intent(in) :: rlat(:,:)

    ! Output arrays
    real(r_def), intent(out) :: fms_net_surf_sw_down(:,:)
    real(r_def), intent(out) :: fms_surf_lw_down(:,:)
    real(r_def), intent(out) :: output_heating_rate(:,:,:)
    real(r_def) :: output_heating_rate_lw(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3))
    real(r_def) :: output_heating_rate_sw(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3))
    real(r_def) :: output_soc_flux_up_lw(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3))
    real(r_def) :: output_soc_flux_down_sw(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3))
    real(r_def), allocatable :: output_soc_spectral_olr(:,:,:)

    ! Hi-res output
    INTEGER, PARAMETER :: out_unit=20
    CHARACTER(len=200) :: file_name
    REAL, allocatable :: soc_spectral_olr(:)

    ! Arrays to send to Socrates
    real, dimension(n_layer) :: input_p, input_t, input_mixing_ratio, &
         input_d_mass, input_density, input_layer_heat_capacity, &
         soc_heating_rate, input_o3_mixing_ratio, &
         soc_heating_rate_lw, soc_heating_rate_sw
    real, dimension(0:n_layer) :: input_p_level, input_t_level, soc_flux_direct, &
         soc_flux_down_sw, soc_flux_up_sw, output_flux_net, &
         soc_flux_down_lw, soc_flux_up_lw
    real, dimension(n_profile) :: input_t_surf, input_cos_zenith_angle, input_solar_irrad, &
         input_orog_corr


    ! Socrates options
    logical :: input_l_planet_grey_surface = .true.
    real(r_def) :: input_planet_albedo = 0.06
    real(r_def) :: input_planet_emissivity = 0.5
    integer(i_def) :: input_n_cloud_layer
    integer(i_def) :: input_n_aer_mode
    integer(i_def) :: input_cld_subcol_gen
    integer(i_def) :: input_cld_subcol_req


    ! Dimensions:
    type(StrDim) :: dimen
    type(StrAtm) :: atm_input

    ! Loop variables
    integer(i_def) :: i, j, k, l, n, lon

    !DIAG Diagnostic
    logical :: used

    !----------------------------i


    ! Allocate spectral array sizes
    if (hires_mode == .FALSE.) then
       allocate(soc_spectral_olr(n_soc_bands_lw))
       allocate(output_soc_spectral_olr(size(fms_temp,1),size(fms_temp,2),n_soc_bands_lw))
    else
       allocate(soc_spectral_olr(n_soc_bands_lw_hires))
       allocate(output_soc_spectral_olr(size(fms_temp,1),size(fms_temp,2),n_soc_bands_lw_hires))
    end if

    ! Set array sizes
    input_n_cloud_layer = n_layer
    input_n_aer_mode = n_layer
    input_cld_subcol_gen = n_layer
    input_cld_subcol_req = n_layer

    do lon = 1,size(fms_temp,1)
       do n = 1, size(fms_temp,2)
       
          nlat = n

          !Set input T, p, p_level, and mixing ratio profiles
          input_t = fms_temp(lon,nlat,:)
          input_p = fms_p_full(lon,nlat,:)
          input_p_level = fms_p_half(lon,nlat,:)

          ! TODO: Remove or edit
          input_mixing_ratio = 0.0!
          input_o3_mixing_ratio = 0.0!

          !-------------

          !Default parameters
          input_cos_zenith_angle = 0.7
          input_orog_corr = 0.0
          input_layer_heat_capacity = 29.07

          !Set tide-locked flux - should be set by namelist eventually!
          input_solar_irrad = fms_stellar_flux(lon,nlat)!RESHAPE(fms_stellar_flux, (/n_profile/))
          input_t_surf = fms_t_surf(lon,nlat)!RESHAPE(fms_t_surf, (/n_profile/))


          !--------------
          !Set input t_level by scaling t - NEEDS TO CHANGE!
          DO i = nlat, nlat
             DO k = 0,n_layer
                input_t_level(k) = 0.5*(input_t(k+1)+input_t(k))
             END DO
             input_t_level(n_layer) = input_t(n_layer) + input_t(n_layer) - input_t_level(n_layer-1)
             input_t_level(0) = input_t(1) - (input_t_level(1) - input_t(1))
          END DO

          !Set input dry mass, density, and heat capacity profiles
          DO i=n_layer, 1, -1
             input_d_mass(i) = (input_p_level(i)-input_p_level(i-1))/grav
             input_density(i) = input_p(i)/(rdgas*input_t(i))
             input_layer_heat_capacity(i) = input_d_mass(i)*cp_air
          END DO


          ! Zero heating rate
          soc_heating_rate = 0.0
          soc_heating_rate_lw = 0.0
          soc_heating_rate_sw = 0.0


          ! Test if LW or SW mode
          if (soc_mode == .TRUE.) then
             control_lw%isolir = 2
             CALL read_control(control_lw, spectrum_lw)
             if (hires_mode == .FALSE.) then
                control_calc = control_lw
                spectrum_calc = spectrum_lw
             else
                control_calc = control_lw_hires
                spectrum_calc = spectrum_lw_hires
             end if

          else
             control_sw%isolir = 1
             CALL read_control(control_sw, spectrum_sw)
             if(hires_mode == .FALSE.) then
                control_calc = control_sw
                spectrum_calc = spectrum_sw
             else
                control_calc = control_sw_hires
                spectrum_calc = spectrum_sw_hires
             end if

          end if


          ! Do calculation

          CALL read_control(control_calc, spectrum_calc)

          CALL socrates_calc(Time_diag, control_calc, spectrum_calc,                                          &
               n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
               input_cld_subcol_gen, input_cld_subcol_req,                                  &
               input_p, input_t, input_t_level, input_d_mass, input_density,                &
               input_mixing_ratio, input_o3_mixing_ratio,                                      &
               input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
               input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
               input_layer_heat_capacity,                                                   &
               soc_flux_direct, soc_flux_down_lw, soc_flux_up_lw, soc_heating_rate_lw, soc_spectral_olr)

          ! Set output arrays
          fms_surf_lw_down(lon,nlat) = soc_flux_down_lw(n_layer)
          output_heating_rate(lon,nlat,:) = soc_heating_rate_lw(:)
          output_soc_spectral_olr(lon,nlat,:) = soc_spectral_olr(:)

          if (soc_mode == .TRUE.) then
             fms_surf_lw_down(lon,nlat) = soc_flux_down_lw(n_layer)
             output_soc_flux_up_lw(lon,nlat,:) = soc_flux_up_lw(:)
             output_heating_rate_lw(lon,nlat,:) = soc_heating_rate_lw(:)
          else
             fms_net_surf_sw_down(lon,nlat) = soc_flux_down_lw(n_layer)
             output_soc_flux_down_sw(lon,nlat,:) = soc_flux_down_lw(:)
             output_heating_rate_sw(lon,nlat,:) = soc_heating_rate_lw(:)
          end if

       end do
    end do




    ! Send diagnostics
    if (soc_mode == .TRUE.) then
       used = send_data ( id_soc_heating_lw, output_heating_rate_lw, Time_diag)
       used = send_data ( id_soc_flux_up_lw, output_soc_flux_up_lw, Time_diag)
    else
       used = send_data ( id_soc_heating_sw, output_heating_rate_sw, Time_diag)
       used = send_data ( id_soc_flux_down_sw, output_soc_flux_down_sw, Time_diag)
    endif

    ! Allocate spectral array sizes
    if (hires_mode == .FALSE.) then
       deallocate(soc_spectral_olr)
       deallocate(output_soc_spectral_olr)
    else
       deallocate(soc_spectral_olr)
       deallocate(output_soc_spectral_olr)
    end if

  end subroutine socrates_interface

subroutine run_socrates(Time_diag, rad_lat, rad_lon, temp_in, t_surf_in, p_full_in, p_half_in, &
       temp_tend, net_surf_sw_down, surf_lw_down, delta_t)  

    use soc_constants_mod

    ! Input time
    type(time_type), intent(in)           :: Time_diag
    real, intent(in), dimension(:,:)      :: rad_lat, rad_lon, t_surf_in
    real, intent(in), dimension(:,:,:)   :: temp_in, p_full_in
    real, intent(in), dimension(:,:,:)  :: p_half_in
    real, intent(inout), dimension(:,:,:) :: temp_tend
    real, intent(out), dimension(:,:)   :: net_surf_sw_down, surf_lw_down
    real, intent(in) :: delta_t

    integer(i_def) :: n_profile, n_layer

    real(r_def), dimension(size(temp_in,1), size(temp_in,2)) :: fms_stellar_flux, output_net_surf_sw_down, output_net_surf_lw_down, output_surf_lw_down, t_surf_for_soc, rad_lat_soc, rad_lon_soc
    real(r_def), dimension(size(temp_in,1), size(temp_in,2), size(temp_in,3)) :: tg_tmp_soc, p_full_soc, output_heating_rate_sw, output_heating_rate_lw, output_heating_rate_total
    real(r_def), dimension(size(temp_in,1), size(temp_in,2), size(temp_in,3)+1) :: p_half_soc

    logical :: socrates_hires_mode, soc_lw_mode, used

       !Set tide-locked flux - should be set by namelist!
       fms_stellar_flux = stellar_constant*COS(rad_lat(:,:))*COS(rad_lon(:,:))
       WHERE (fms_stellar_flux < 0.0) fms_stellar_flux = 0.0

       n_profile = INT(1, kind(i_def))
       n_layer   = INT(size(temp_in,3), kind(i_def))
       t_surf_for_soc = REAL(t_surf_in(:,:), kind(r_def))

       ! GCM mode
       socrates_hires_mode = .FALSE.

       ! LW calculation
       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
       soc_lw_mode = .TRUE.
    
       rad_lat_soc = REAL(rad_lat, kind(r_def))
       rad_lon_soc = REAL(rad_lon, kind(r_def))
       tg_tmp_soc =  REAL(temp_in, kind(r_def))
       p_full_soc = REAL(p_full_in, kind(r_def))
       p_half_soc = REAL(p_half_in, kind(r_def))

       CALL socrates_interface(Time_diag, rad_lat_soc, rad_lon_soc, soc_lw_mode, socrates_hires_mode,    &
            tg_tmp_soc, t_surf_for_soc, p_full_soc, p_half_soc, n_profile, n_layer,     &
            output_heating_rate_lw, output_net_surf_sw_down, output_surf_lw_down, fms_stellar_flux )

       tg_tmp_soc = tg_tmp_soc + output_heating_rate_lw*delta_t
       surf_lw_down(:,:) = REAL(output_surf_lw_down(:,:))

       temp_tend(:,:,:) = temp_tend(:,:,:) + real(output_heating_rate_lw)
       
       ! SW calculation
       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
       soc_lw_mode = .FALSE.
       CALL socrates_interface(Time_diag, rad_lat_soc, rad_lon_soc, soc_lw_mode, socrates_hires_mode,    &
            tg_tmp_soc, t_surf_for_soc, p_full_soc, p_half_soc, n_profile, n_layer,     &
            output_heating_rate_sw, output_net_surf_sw_down, output_surf_lw_down, fms_stellar_flux )

       tg_tmp_soc = tg_tmp_soc + output_heating_rate_sw*delta_t
       net_surf_sw_down(:,:) = REAL(output_net_surf_sw_down(:,:))

       temp_tend(:,:,:) = temp_tend(:,:,:) + real(output_heating_rate_sw)
       
       output_heating_rate_total = output_heating_rate_lw + output_heating_rate_sw

       !Sending total heating rates
       used = send_data ( id_soc_heating_rate, output_heating_rate_total, Time_diag)

end subroutine run_socrates  
  
end module socrates_interface_mod
