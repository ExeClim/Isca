MODULE socrates_interface_mod

  ! Socrates calculation interface module
  ! Takes FMS time, spectra, temperature, and pressure
  ! Outputs FMS heating rate, and downwards surface LW and SW
  ! MDH 30/01/18

  !----------

  ! ExoFMS diagnostics
  USE  diag_manager_mod, ONLY: register_diag_field, send_data

  ! ExoFMS time
  USE time_manager_mod, ONLY: time_type, OPERATOR(+), OPERATOR(-), OPERATOR(/=)

  ! Socrates modules
  USE read_control_mod
  USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
  USE def_spectrum


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
  NAMELIST/socrates_nml/ stellar_constant, tidally_locked



CONTAINS

  SUBROUTINE socrates_init(is, ie, js, je, num_levels, axes, Time, lat)
    !! Initialises Socrates spectra, arrays, and constants

    ! Arguments
    INTEGER, INTENT(in), DIMENSION(4) :: axes
    !! NB axes refers to the handles of the axes defined in fv_diagnostics
    TYPE(time_type), INTENT(in)       :: Time
    INTEGER, INTENT(in)               :: is, ie, js, je, num_levels
    REAL, INTENT(in) , DIMENSION(:,:)   :: lat
    !-------------------------------------------------------------------------------------

    ! Socrates spectral files -- should be set by namelist
    control_lw%spectral_file = '/scratch/sit204/Isca/src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'
    control_lw_hires%spectral_file = '/scratch/sit204/Isca/src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'
    !control_lw_hires%spectral_file = '~/Work/spec_file_co2_co'

    control_sw%spectral_file = '/scratch/sit204/Isca/src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'
    control_sw_hires%spectral_file = '/scratch/sit204/Isca/src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'

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

    write(6,*) 'checking arrays passed in'
    write(6,*) n_layer
    write(6,*) n_profile

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

        write(6,*) 'do I make it here?', lon, n, size(fms_temp,1), size(fms_temp,2), size(fms_temp,3)

          !Set input t_level by scaling t - NEEDS TO CHANGE!
          DO i = nlat, nlat
             DO k = 0,n_layer
                input_t_level(k) = 0.5*(input_t(k+1)+input_t(k))
             END DO
             input_t_level(n_layer) = input_t(n_layer) + input_t(n_layer) - input_t_level(n_layer-1)
             input_t_level(0) = input_t(1) - (input_t_level(1) - input_t(1))
          END DO

        write(6,*) 'do I make it past here?'



          !Set input dry mass, density, and heat capacity profiles
          DO i=n_layer, 1, -1
             input_d_mass(i) = (input_p_level(i)-input_p_level(i-1))/23.0
             input_density(i) = input_p(i)/(8.31*input_t(i))!1000.!atm%p(l ,i) / 1000
             !KLUDGE
             input_layer_heat_capacity(i) = input_d_mass(i)*1303.1!17.0*1005.0
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
          control_calc%isolir = 2
          write(6,*) 'Read control'

          CALL read_control(control_calc, spectrum_calc)

          write(6,*) 'made it out of read_control'
        write(6,*) Time_diag, 'Time diag'
        write(6,*) control_calc, 'control calc'
        write(6,*) spectrum_calc, 'spectrum calc'
        write(6,*) ' done 1'
        write(6,*) n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode, input_cld_subcol_gen, input_cld_subcol_req
        write(6,*) ' done 2'        
        write(6,*) size(input_p), size(input_t), size(input_t_level), size(input_d_mass), size(input_density)
                write(6,*) ' done 3'
        write(6,*) input_mixing_ratio, input_o3_mixing_ratio
        write(6,*) ' done 4'        
        write(6,*) size(input_t_surf), size(input_cos_zenith_angle), size(input_solar_irrad), size(input_orog_corr)
        write(6,*) ' done 5'        
        write(6,*) input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity
        write(6,*) ' done 6'        
        write(6,*) input_layer_heat_capacity, soc_flux_direct, soc_flux_down_lw, soc_flux_up_lw, soc_heating_rate_lw, soc_spectral_olr

          write(6,*) 'Socrates calc step'
          CALL socrates_calc(Time_diag, control_calc, spectrum_calc,                                          &
               n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
               input_cld_subcol_gen, input_cld_subcol_req,                                  &
               input_p, input_t, input_t_level, input_d_mass, input_density,                &
               input_mixing_ratio, input_o3_mixing_ratio,                                      &
               input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
               input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
               input_layer_heat_capacity,                                                   &
               soc_flux_direct, soc_flux_down_lw, soc_flux_up_lw, soc_heating_rate_lw, soc_spectral_olr)

          write(6,*) 'made it out of soc calc step'
          write(6,*) 'Output arrays', n_profile, n_layer
          write(6,*) 'Output arrays', size(output_heating_rate,1)          
          write(6,*) 'Output arrays', size(fms_surf_lw_down, 1), size(fms_surf_lw_down, 2)
          write(6,*) 'Output arrays', size(soc_flux_down_lw)

          ! Set output arrays
          fms_surf_lw_down(lon,nlat) = soc_flux_down_lw(n_layer)
          write(6,*) ' done surf lw down'
          output_heating_rate(lon,nlat,:) = soc_heating_rate_lw(:)
          write(6,*) ' done heating rate'          
          output_soc_spectral_olr(lon,nlat,:) = soc_spectral_olr(:)


          write(6,*) 'fluxes etc', soc_mode


          if (soc_mode == .TRUE.) then
             fms_surf_lw_down(lon,nlat) = soc_flux_down_lw(n_layer)
             output_soc_flux_up_lw(lon,nlat,:) = soc_flux_up_lw(:)
             output_heating_rate_lw(lon,nlat,:) = soc_heating_rate_lw(:)
          else
             fms_net_surf_sw_down(lon,nlat) = soc_flux_down_lw(n_layer)

             output_soc_flux_down_sw(lon,nlat,:) = soc_flux_down_lw(:)
             output_heating_rate_sw(lon,nlat,:) = soc_heating_rate_lw(:)
          end if

          write(6,*) 'end inner loop'


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


  end subroutine socrates_interface
end module socrates_interface_mod
