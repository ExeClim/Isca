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
  USE time_manager_mod, ONLY: time_type, OPERATOR(+), OPERATOR(-), OPERATOR(/=), length_of_day, length_of_year, get_time, set_time

  ! Socrates modules
  USE read_control_mod
  USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
  USE def_spectrum
  USE constants_mod, only: grav, rdgas, rvgas, cp_air
  USE fms_mod, only: stdlog, FATAL, WARNING, error_mesg
  USE interpolator_mod, only: interpolate_type  
  USE soc_constants_mod  

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
  INTEGER :: id_soc_flux_up_lw, id_soc_flux_down_sw, id_coszen, id_interp_temp, id_norm_temp
  INTEGER :: n_soc_bands_lw, n_soc_bands_sw
  INTEGER :: n_soc_bands_lw_hires, n_soc_bands_sw_hires
  CHARACTER(len=10), PARAMETER :: soc_mod_name = 'socrates'
  REAL :: missing_value = -999

  type(interpolate_type),save                :: o3_interp, co2_interp            ! use external file for ozone and co2

  REAL :: dt_last !Time of last radiation calculation - used to tell whether it is time to recompute radiation or not
  REAL(r_def), allocatable, dimension(:,:,:) :: tdt_soc_sw_store, tdt_soc_lw_store, thd_sw_flux_down_store, thd_lw_flux_up_store
  REAL(r_def), allocatable, dimension(:,:) :: net_surf_sw_down_store, surf_lw_down_store, coszen_store
  REAL(r_def), allocatable, dimension(:,:,:) :: outputted_soc_spectral_olr

CONTAINS

  SUBROUTINE socrates_init(is, ie, js, je, num_levels, axes, Time, lat, lonb, latb, delta_t_atmos)
    !! Initialises Socrates spectra, arrays, and constants

    USE astronomy_mod, only: astronomy_init
    USE interpolator_mod, only: interpolate_type, interpolator_init, ZERO  
    USE socrates_config_mod

    ! Arguments
    INTEGER, INTENT(in), DIMENSION(4) :: axes
    !! NB axes refers to the handles of the axes defined in fv_diagnostics
    TYPE(time_type), INTENT(in)       :: Time, delta_t_atmos
    INTEGER, INTENT(in)               :: is, ie, js, je, num_levels
    REAL, INTENT(in) , DIMENSION(:,:)   :: lat
    REAL, INTENT(in) , DIMENSION(:,:)   :: lonb, latb
        
    integer :: io, stdlog_unit
    integer :: res, time_step_seconds
    real    :: day_in_s_check

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

    !Initialise astronomy
    call astronomy_init

    !Initialise variables related to radiation timestep

          call get_time(delta_t_atmos,time_step_seconds)

          if (dt_rad .le. 0.) then
              dt_rad = time_step_seconds !Make sure that dt_rad is set if it is not specified in the namelist
          endif

          dt_last = -real(dt_rad) !make sure we are computing radiation at the first time step

            if (dt_rad .gt. time_step_seconds) then
                res=mod(dt_rad, time_step_seconds)

                if (res.ne.0) then
                        call error_mesg( 'socrates_init', &
                         'dt_rad must be an integer multiple of dt_atmos',FATAL)
                endif

                day_in_s_check=length_of_day()
                res=mod(int(day_in_s_check), dt_rad)

                if (res.ne.0) then
                        call error_mesg( 'socrates_init', &
                         'dt_rad does not fit into one day an integer number of times', WARNING)
                endif


            endif

          if(dt_rad_avg .le. 0) dt_rad_avg = dt_rad

    IF (js == 1) THEN

      if (lw_spectral_filename .eq. 'unset') then
       call error_mesg( 'socrates_init', &
                       'lw_spectral_filename is unset, and must point to a valid spectral file',FATAL)
      endif

      if (sw_spectral_filename .eq. 'unset') then
       call error_mesg( 'socrates_init', &
                       'sw_spectral_filename is unset, and must point to a valid spectral file',FATAL)
      endif
    ENDIF

    if (lw_hires_spectral_filename .eq. 'unset') then
       IF (js == 1) THEN
           call error_mesg( 'socrates_init', &
                       'lw_hires_spectral_filename is unset, making equal to lw_spectral_filename',WARNING)
       endif
        lw_hires_spectral_filename = lw_spectral_filename
    endif

    if (sw_hires_spectral_filename .eq. 'unset') then
       IF (js == 1) THEN
           call error_mesg( 'socrates_init', &
                       'sw_hires_spectral_filename is unset, making equal to sw_spectral_filename',WARNING)
       endif
        sw_hires_spectral_filename = sw_spectral_filename
    endif


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

    id_interp_temp = &
         register_diag_field ( soc_mod_name, 'temp_half', (/axes(1),axes(2),axes(4)/), Time, &
         'Interpolated half-level temperatures', &
         'K', missing_value=missing_value               )

    id_norm_temp = &
         register_diag_field ( soc_mod_name, 'temp_norm', (/axes(1),axes(2),axes(3)/), Time, &
         'Normal full-level temperatures', &
         'K', missing_value=missing_value               )


    id_coszen = &
         register_diag_field ( soc_mod_name, 'soc_coszen', axes(1:2), Time, &
         'Cosine(zenith_angle)', &
         'none', missing_value=missing_value               )

      if(do_read_ozone)then
         call interpolator_init (o3_interp, trim(ozone_file_name)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
      endif
      
      if(do_read_co2)then
         call interpolator_init (co2_interp, trim(co2_file_name)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
      endif      

    if (mod((size(lonb,1)-1)*(size(latb,1)-1), chunk_size) .ne. 0) then
    
        call error_mesg( 'socrates_init', &
                     'chunk_size must equally divide number of points per processor, which it currently does not.',FATAL)
    
    endif

    ! Number of bands
    n_soc_bands_lw = spectrum_lw%dim%nd_band
    n_soc_bands_lw_hires = spectrum_lw_hires%dim%nd_band
    n_soc_bands_sw = spectrum_sw%dim%nd_band
    n_soc_bands_sw_hires = spectrum_sw_hires%dim%nd_band

    if (socrates_hires_mode == .True.) then
        allocate(outputted_soc_spectral_olr(size(lonb,1)-1, size(latb,2)-1, n_soc_bands_lw_hires))
    else
        allocate(outputted_soc_spectral_olr(size(lonb,1)-1, size(latb,2)-1, n_soc_bands_lw ))
    endif

    if(store_intermediate_rad) then
        allocate(tdt_soc_sw_store(size(lonb,1)-1, size(latb,2)-1, num_levels))
        allocate(tdt_soc_lw_store(size(lonb,1)-1, size(latb,2)-1, num_levels))

        if (id_soc_flux_up_lw > 0) then
            allocate(thd_lw_flux_up_store(size(lonb,1)-1, size(latb,2)-1, num_levels))
        endif

        if (id_soc_flux_down_sw > 0) then
            allocate(thd_sw_flux_down_store(size(lonb,1)-1, size(latb,2)-1, num_levels))
        endif

        if (id_coszen > 0) then
            allocate(coszen_store(size(lonb,1)-1, size(latb,2)-1))
        endif

        allocate(net_surf_sw_down_store(size(lonb,1)-1, size(latb,2)-1))
        allocate(surf_lw_down_store(size(lonb,1)-1, size(latb,2)-1))
    endif

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
  subroutine socrates_interface(Time_diag, rlat, rlon, soc_mode,  &
       fms_temp, fms_spec_hum, fms_ozone, fms_co2, fms_t_surf, fms_p_full, fms_p_half, fms_z_full, fms_z_half, fms_albedo, n_profile, n_layer,        &
       output_heating_rate, output_soc_flux_down_sw, output_soc_flux_up_lw, fms_net_surf_sw_down, fms_surf_lw_down, fms_coszen, output_soc_spectral_olr, t_half_level_out )

    use realtype_rd
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
    use socrates_config_mod
    !-----------------------------------------------------------------------
    implicit none

    ! Input time
    type(time_type), intent(in)         :: Time_diag

    INTEGER(i_def), intent(in) :: n_profile, n_layer
    logical, intent(in) :: soc_mode
    INTEGER(i_def) :: nlat

    ! Input arrays
    real(r_def), intent(in) :: fms_temp(:,:,:), fms_spec_hum(:,:,:), fms_ozone(:,:,:), fms_co2(:,:,:)
    real(r_def), intent(in) :: fms_p_full(:,:,:)
    real(r_def), intent(in) :: fms_p_half(:,:,:)
    real(r_def), intent(in) :: fms_t_surf(:,:), fms_albedo(:,:)
    real(r_def), intent(in) :: fms_coszen(:,:)
    real(r_def), intent(in) :: rlon(:,:)
    real(r_def), intent(in) :: rlat(:,:)
    real(r_def), intent(in) :: fms_z_full(:,:,:), fms_z_half(:,:,:)


    ! Output arrays
    real(r_def), intent(out) :: fms_net_surf_sw_down(:,:)
    real(r_def), intent(out) :: fms_surf_lw_down(:,:)
    real(r_def), intent(out) :: output_heating_rate(:,:,:)
    real(r_def), intent(out) :: output_soc_flux_up_lw(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3))
    real(r_def), intent(out) :: output_soc_flux_down_sw(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3))
    real(r_def), intent(out), optional :: output_soc_spectral_olr(:,:,:)
    real(r_def)              :: output_heating_rate_lw(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3))
    real(r_def)              :: output_heating_rate_sw(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3))
    real(r_def), intent(out), optional :: t_half_level_out(size(fms_temp,1),size(fms_temp,2),size(fms_temp,3)+1)

    ! Hi-res output
    INTEGER, PARAMETER :: out_unit=20
    CHARACTER(len=200) :: file_name
    REAL :: soc_spectral_olr(n_profile, size(outputted_soc_spectral_olr,3))

    ! Arrays to send to Socrates
    real, dimension(n_profile, n_layer) :: input_p, input_t, input_mixing_ratio, &
         input_d_mass, input_density, input_layer_heat_capacity, &
         soc_heating_rate, input_o3_mixing_ratio, &
         soc_heating_rate_lw, soc_heating_rate_sw, input_co2_mixing_ratio,z_full_reshaped
    real, dimension(n_profile, 0:n_layer) :: input_p_level, input_t_level, soc_flux_direct, &
         soc_flux_down_sw, soc_flux_up_sw, output_flux_net, &
         soc_flux_down_lw, soc_flux_up_lw, z_half_reshaped
    real, dimension(n_profile) :: input_t_surf, input_cos_zenith_angle, input_solar_irrad, &
         input_orog_corr, input_planet_albedo


    ! Socrates options
    integer(i_def) :: input_n_cloud_layer
    integer(i_def) :: input_n_aer_mode
    integer(i_def) :: input_cld_subcol_gen
    integer(i_def) :: input_cld_subcol_req


    ! Dimensions:
    type(StrDim) :: dimen
    type(StrAtm) :: atm_input

    ! Loop variables
    integer(i_def) :: i, j, k, l, n, lon, si, sj, sk

    ! chunking loop variable
    integer(i_def) :: n_chunk_loop, idx_chunk_start, idx_chunk_end, i_chunk, n_profile_chunk

    !DIAG Diagnostic
    logical :: used

    !----------------------------i

    ! Set array sizes
    input_n_cloud_layer = n_layer
    input_n_aer_mode = n_layer
    input_cld_subcol_gen = n_layer
    input_cld_subcol_req = n_layer
    si = size(fms_temp,1)
    sj = size(fms_temp,2)
    sk = size(fms_temp,3)

          !Set input T, p, p_level, and mixing ratio profiles
          input_t = reshape(fms_temp(:,:,:),(/si*sj,sk /))
          input_p = reshape(fms_p_full(:,:,:),(/si*sj,sk /))
          input_p_level = reshape(fms_p_half(:,:,:),(/si*sj,sk+1 /))
          
          if (account_for_effect_of_water == .true.) then
              input_mixing_ratio = reshape(fms_spec_hum(:,:,:) / (1. - fms_spec_hum(:,:,:)),(/si*sj,sk /)) !Mass mixing ratio = q / (1-q)
          else
              input_mixing_ratio = 0.0
          endif
          
          if (account_for_effect_of_ozone == .true.) then
            input_o3_mixing_ratio = reshape(fms_ozone(:,:,:),(/si*sj,sk /))
          else         
            input_o3_mixing_ratio = 0.0
          endif

         input_co2_mixing_ratio = reshape(fms_co2(:,:,:),(/si*sj,sk /))

          !-------------

          !Default parameters
          input_cos_zenith_angle = reshape((fms_coszen(:,:)),(/si*sj /))
          input_orog_corr = 0.0
          input_planet_albedo = reshape(fms_albedo(:,:),(/n_profile /))

          !Set tide-locked flux - should be set by namelist eventually!
          input_solar_irrad = stellar_constant
          input_t_surf = reshape(fms_t_surf(:,:),(/si*sj /))
          z_full_reshaped = reshape(fms_z_full(:,:,:), (/si*sj, sk/))
          z_half_reshaped = reshape(fms_z_half(:,:,:), (/si*sj, sk+1/))

          !--------------
          !Set input t_level by scaling t - NEEDS TO CHANGE!
       if (use_pressure_interp_for_half_levels) then
          DO i = nlat, nlat
             DO k = 0,n_layer
                input_t_level(:,k) = input_t(:,k) + (input_t(:,k+1)-input_t(:,k)) * ((input_p_level(:,k)-input_p(:,k))/(input_p(:,k+1)-input_p(:,k)))
             END DO
!              input_t_level(:,n_layer) = input_t(:,n_layer) + input_t(:,n_layer) - input_t_level(:,n_layer-1)
             input_t_level(:,n_layer) = input_t(:,n_layer) + (input_t(:,n_layer)-input_t(:,n_layer-1)) * ((input_p_level(:,n_layer)-input_p(:,n_layer))/(input_p(:,n_layer)-input_p(:,n_layer-1)))
             
!              input_t_level(:,0) = input_t(:,1) - (input_t_level(:,1) - input_t(:,1))
             input_t_level(:,0) = input_t(:,1) + (input_t(:,2)-input_t(:,1)) * ((input_p_level(:,0)-input_p(:,1))/(input_p(:,2)-input_p(:,1)))
             
          END DO
       else

          call interp_temp(z_full_reshaped,z_half_reshaped,input_t, input_t_level)

       endif

         if (present(t_half_level_out)) then
             t_half_level_out(:,:,:) = reshape(input_t_level,(/si,sj,sk+1 /))
         endif

          !Set input dry mass, density, and heat capacity profiles
          DO i=n_layer, 1, -1
             input_d_mass(:,i) = (input_p_level(:,i)-input_p_level(:,i-1))/grav
             input_density(:,i) = input_p(:,i)/(rdgas*input_t(:,i))
             input_layer_heat_capacity(:,i) = input_d_mass(:,i)*cp_air
          END DO


          ! Zero heating rate
          soc_heating_rate = 0.0
          soc_heating_rate_lw = 0.0
          soc_heating_rate_sw = 0.0


          ! Test if LW or SW mode
          if (soc_mode == .TRUE.) then
             control_lw%isolir = 2
             CALL read_control(control_lw, spectrum_lw)
             if (socrates_hires_mode == .FALSE.) then
                control_calc = control_lw
                spectrum_calc = spectrum_lw
             else
                control_calc = control_lw_hires
                spectrum_calc = spectrum_lw_hires
             end if

          else
             control_sw%isolir = 1
             CALL read_control(control_sw, spectrum_sw)
             if(socrates_hires_mode == .FALSE.) then
                control_calc = control_sw
                spectrum_calc = spectrum_sw
             else
                control_calc = control_sw_hires
                spectrum_calc = spectrum_sw_hires
             end if

          end if


          ! Do calculation
          CALL read_control(control_calc, spectrum_calc)

    n_chunk_loop = (si*sj)/chunk_size
    n_profile_chunk = n_profile / n_chunk_loop

        DO i_chunk=1,n_chunk_loop

            idx_chunk_start = (i_chunk-1)*chunk_size + 1
            idx_chunk_end   = (i_chunk)*chunk_size

          if (soc_mode==.TRUE.) then
 
            CALL socrates_calc(Time_diag, control_calc, spectrum_calc,                      &
               n_profile_chunk, n_layer, input_n_cloud_layer, input_n_aer_mode,             &
               input_cld_subcol_gen, input_cld_subcol_req,                                  &
               input_p(idx_chunk_start:idx_chunk_end,:),                                    & 
               input_t(idx_chunk_start:idx_chunk_end,:),                                    &
               input_t_level(idx_chunk_start:idx_chunk_end,:),                              & 
               input_d_mass(idx_chunk_start:idx_chunk_end,:),                               &
               input_density(idx_chunk_start:idx_chunk_end,:),                              &
               input_mixing_ratio(idx_chunk_start:idx_chunk_end,:),                         &
               input_o3_mixing_ratio(idx_chunk_start:idx_chunk_end,:),                      &
               input_co2_mixing_ratio(idx_chunk_start:idx_chunk_end,:),                     &
               input_t_surf(idx_chunk_start:idx_chunk_end),                                 &
               input_cos_zenith_angle(idx_chunk_start:idx_chunk_end),                       &
               input_solar_irrad(idx_chunk_start:idx_chunk_end),                            &
               input_orog_corr(idx_chunk_start:idx_chunk_end),                              &
               l_planet_grey_surface,                                                       &
               input_planet_albedo(idx_chunk_start:idx_chunk_end),                          &
               input_planet_emissivity,                                                     &
               input_layer_heat_capacity(idx_chunk_start:idx_chunk_end,:),                  &
               soc_flux_direct(idx_chunk_start:idx_chunk_end,:),                            &
               soc_flux_down_lw(idx_chunk_start:idx_chunk_end,:),                           &
               soc_flux_up_lw(idx_chunk_start:idx_chunk_end,:),                             &
               soc_heating_rate_lw(idx_chunk_start:idx_chunk_end,:),                        &
               soc_spectral_olr(idx_chunk_start:idx_chunk_end,:))

          else
            CALL socrates_calc(Time_diag, control_calc, spectrum_calc,                      &
               n_profile_chunk, n_layer, input_n_cloud_layer, input_n_aer_mode,             &
               input_cld_subcol_gen, input_cld_subcol_req,                                  &
               input_p(idx_chunk_start:idx_chunk_end,:),                                    &
               input_t(idx_chunk_start:idx_chunk_end,:),                                    &
               input_t_level(idx_chunk_start:idx_chunk_end,:),                              &
               input_d_mass(idx_chunk_start:idx_chunk_end,:),                               &
               input_density(idx_chunk_start:idx_chunk_end,:),                              &
               input_mixing_ratio(idx_chunk_start:idx_chunk_end,:),                         &
               input_o3_mixing_ratio(idx_chunk_start:idx_chunk_end,:),                      &
               input_co2_mixing_ratio(idx_chunk_start:idx_chunk_end,:),                     &
               input_t_surf(idx_chunk_start:idx_chunk_end),                                 &
               input_cos_zenith_angle(idx_chunk_start:idx_chunk_end),                       &
               input_solar_irrad(idx_chunk_start:idx_chunk_end),                            &
               input_orog_corr(idx_chunk_start:idx_chunk_end),                              &
               l_planet_grey_surface,                                                       &
               input_planet_albedo(idx_chunk_start:idx_chunk_end),                          &
               input_planet_emissivity,                                                     &
               input_layer_heat_capacity(idx_chunk_start:idx_chunk_end,:),                  &
               soc_flux_direct(idx_chunk_start:idx_chunk_end,:),                            &
               soc_flux_down_lw(idx_chunk_start:idx_chunk_end,:),                           &
               soc_flux_up_lw(idx_chunk_start:idx_chunk_end,:),                             &
               soc_heating_rate_lw(idx_chunk_start:idx_chunk_end,:))
          endif

        ENDDO

          ! Set output arrays
          fms_surf_lw_down(:,:) = reshape(soc_flux_down_lw(:,n_layer),(/si,sj/))
          output_heating_rate(:,:,:) = reshape(soc_heating_rate_lw(:,:),(/si,sj,sk /))

          if (soc_mode == .TRUE.) then
              output_soc_spectral_olr(:,:,:) = reshape(soc_spectral_olr(:,:),(/si,sj,int(n_soc_bands_lw,i_def) /))
          endif

          if (soc_mode == .TRUE.) then
             fms_surf_lw_down(:,:) = reshape(soc_flux_down_lw(:,n_layer),(/si,sj /))
             output_soc_flux_up_lw(:,:,:) = reshape(soc_flux_up_lw(:,:),(/si,sj,sk /))
             output_heating_rate_lw(:,:,:) = reshape(soc_heating_rate_lw(:,:),(/si,sj,sk/))
          else
             fms_net_surf_sw_down(:,:) = reshape(soc_flux_down_lw(:,n_layer)-soc_flux_up_lw(:,n_layer),(/si,sj /))
             output_soc_flux_down_sw(:,:,:) = reshape(soc_flux_down_lw(:,:),(/si,sj,sk /))
             output_heating_rate_sw(:,:,:) = reshape(soc_heating_rate_lw(:,:),(/si,sj,sk /))
          end if

  end subroutine socrates_interface

subroutine run_socrates(Time, Time_diag, rad_lat, rad_lon, temp_in, q_in, t_surf_in, p_full_in, p_half_in, z_full_in, z_half_in, albedo_in, &
       temp_tend, net_surf_sw_down, surf_lw_down, delta_t)  

    use astronomy_mod, only: diurnal_solar
    use constants_mod,         only: pi
    use interpolator_mod,only: interpolator
    USE socrates_config_mod    

    ! Input time
    type(time_type), intent(in)           :: Time, Time_diag
    real, intent(in), dimension(:,:)      :: rad_lat, rad_lon, t_surf_in, albedo_in
    real, intent(in), dimension(:,:,:)   :: temp_in, p_full_in, q_in, z_full_in
    real, intent(in), dimension(:,:,:)  :: p_half_in, z_half_in
    real, intent(inout), dimension(:,:,:) :: temp_tend
    real, intent(out), dimension(:,:)   :: net_surf_sw_down, surf_lw_down
    real, intent(in) :: delta_t

    integer(i_def) :: n_profile, n_layer

    real(r_def), dimension(size(temp_in,1), size(temp_in,2)) :: output_net_surf_sw_down, output_net_surf_lw_down, output_surf_lw_down, t_surf_for_soc, rad_lat_soc, rad_lon_soc, albedo_soc
    real(r_def), dimension(size(temp_in,1), size(temp_in,2), size(temp_in,3)) :: tg_tmp_soc, q_soc, ozone_soc, co2_soc, p_full_soc, output_heating_rate_sw, output_heating_rate_lw, output_heating_rate_total, output_soc_flux_down_sw, output_soc_flux_up_lw, z_full_soc
    real(r_def), dimension(size(temp_in,1), size(temp_in,2), size(temp_in,3)+1) :: p_half_soc, t_half_out, z_half_soc

    logical :: soc_lw_mode, used
    integer :: seconds, days, year_in_s
    real :: r_seconds, r_days, r_total_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, dt_rad_radians, day_in_s, r_solday, r_dt_rad_avg
    real, dimension(size(temp_in,1), size(temp_in,2)) :: coszen, fracsun   
    real, dimension(size(temp_in,1), size(temp_in,2), size(temp_in,3)) :: ozone_in, co2_in, thd_sw_flux_down, thd_lw_flux_up
    type(time_type) :: Time_loc

        !check if we really want to recompute radiation
        ! alarm
          call get_time(Time,seconds,days)
                  r_days = real(days)
                  r_seconds = real(seconds)
                  r_total_seconds=r_seconds+(r_days*86400.)
          if(r_total_seconds - dt_last .ge. dt_rad) then
             dt_last = r_total_seconds
          else
             if(store_intermediate_rad) then
                output_heating_rate_sw  = tdt_soc_sw_store
                output_heating_rate_lw  = tdt_soc_lw_store

                if (id_soc_flux_down_sw > 0) then
                    thd_sw_flux_down = thd_sw_flux_down_store
                endif

                if (id_soc_flux_up_lw > 0) then
                    thd_lw_flux_up   = thd_lw_flux_up_store
                endif

                if (id_coszen > 0) then
                    coszen   = coszen_store
                endif

                net_surf_sw_down = real(net_surf_sw_down_store)
                surf_lw_down = real(surf_lw_down_store)
             else
                output_heating_rate_sw = 0.
                output_heating_rate_lw = 0.
                thd_sw_flux_down = 0.
                thd_lw_flux_up   = 0.
                net_surf_sw_down  = 0.
                surf_lw_down  = 0.
             endif
             temp_tend(:,:,:) = temp_tend(:,:,:) + real(output_heating_rate_sw)+real(output_heating_rate_lw)
             output_heating_rate_total = output_heating_rate_sw +output_heating_rate_lw            
            ! Send diagnostics
    if(id_soc_heating_lw > 0) then
      used = send_data ( id_soc_heating_lw, output_heating_rate_lw, Time_diag)
    endif

    if(id_soc_flux_up_lw > 0) then
      used = send_data ( id_soc_flux_up_lw, thd_lw_flux_up, Time_diag)
    endif

    if(id_soc_heating_sw > 0) then
      used = send_data ( id_soc_heating_sw, output_heating_rate_sw, Time_diag)
    endif

    if(id_soc_flux_down_sw > 0) then
      used = send_data ( id_soc_flux_down_sw, thd_sw_flux_down, Time_diag)
    endif

    if(id_soc_heating_rate > 0) then
      used = send_data ( id_soc_heating_rate, output_heating_rate_total, Time_diag)
    endif             

    if(id_coszen > 0) then
      used = send_data ( id_coszen, coszen, Time_diag)
    endif

             return !not time yet
          endif


       !make sure we run perpetual when solday > 0)
          if(solday > 0)then
             Time_loc = set_time(seconds,solday)
          else
             Time_loc = Time
          endif

       !Set tide-locked flux if tidally-locked = .true. Else use diurnal-solar
       !to calculate insolation from orbit!
       if (tidally_locked == .TRUE.) then
           coszen = COS(rad_lat(:,:))*COS(rad_lon(:,:))
           WHERE (coszen < 0.0) coszen = 0.0
       else
       
        ! compute zenith angle
                 call get_time(Time_loc, seconds, days)
                 call get_time(length_of_year(), year_in_s)
                 day_in_s = length_of_day()
         
                 r_seconds=real(seconds)
                 r_days=real(days)
                 r_total_seconds=r_seconds+(r_days*86400.)
         
                 frac_of_day = r_total_seconds / day_in_s

                 if(solday > 0) then
                     r_solday=real(solday)
                     frac_of_year=(r_solday*day_in_s)/year_in_s
                 else
                     frac_of_year = r_total_seconds / year_in_s
                 endif
                 gmt = abs(mod(frac_of_day, 1.0)) * 2.0 * pi
                 time_since_ae = modulo(frac_of_year-equinox_day, 1.0) * 2.0 * pi   
       
           if(do_rad_time_avg) then
	         r_dt_rad_avg=real(dt_rad_avg)
	         dt_rad_radians = (r_dt_rad_avg/day_in_s)*2.0*pi
	         call diurnal_solar(rad_lat, rad_lon, gmt, time_since_ae, coszen, fracsun, rrsun, dt_rad_radians)
           else
	         ! Seasonal Cycle: Use astronomical parameters to calculate insolation
	         call diurnal_solar(rad_lat, rad_lon, gmt, time_since_ae, coszen, fracsun, rrsun)
           end if
           
       endif

       ozone_in = 0.0
       
      !get ozone 
       if(do_read_ozone)then
         call interpolator( o3_interp, Time_diag, p_half_in, ozone_in, trim(ozone_field_name))
       endif

       co2_in = co2_ppmv * 1.e-6

      !get ozone 
       if(do_read_co2)then
         call interpolator( co2_interp, Time_diag, p_half_in, co2_in, trim(co2_field_name))
         co2_in = co2_in * 1.e-6
       endif


       n_profile = INT(size(temp_in,2)*size(temp_in,1), kind(i_def))
       n_layer   = INT(size(temp_in,3), kind(i_def))
       t_surf_for_soc = REAL(t_surf_in(:,:), kind(r_def))
       
       ! LW calculation
       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
       soc_lw_mode = .TRUE.
   
       rad_lat_soc = REAL(rad_lat, kind(r_def))
       rad_lon_soc = REAL(rad_lon, kind(r_def))
       tg_tmp_soc =  REAL(temp_in, kind(r_def))
       q_soc      =  REAL(q_in, kind(r_def))
       ozone_soc  =  REAL(ozone_in, kind(r_def)) 
       co2_soc    =  REAL(co2_in, kind(r_def))      
       p_full_soc = REAL(p_full_in, kind(r_def))
       p_half_soc = REAL(p_half_in, kind(r_def))
       albedo_soc = REAL(albedo_in, kind(r_def))
       z_full_soc = REAL(z_full_in, kind(r_def))
       z_half_soc = REAL(z_half_in, kind(r_def))

       CALL socrates_interface(Time, rad_lat_soc, rad_lon_soc, soc_lw_mode,  &
            tg_tmp_soc, q_soc, ozone_soc, co2_soc, t_surf_for_soc, p_full_soc, p_half_soc, z_full_soc, z_half_soc, albedo_soc, n_profile, n_layer,     &
            output_heating_rate_lw, output_soc_flux_down_sw, output_soc_flux_up_lw,  output_net_surf_sw_down, output_surf_lw_down, coszen, outputted_soc_spectral_olr, t_half_level_out = t_half_out)

       tg_tmp_soc = tg_tmp_soc + output_heating_rate_lw*delta_t
       surf_lw_down(:,:) = REAL(output_surf_lw_down(:,:))
       thd_lw_flux_up = REAL(output_soc_flux_up_lw)

       temp_tend(:,:,:) = temp_tend(:,:,:) + real(output_heating_rate_lw)
       
       ! SW calculation
       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
       soc_lw_mode = .FALSE.
       CALL socrates_interface(Time, rad_lat_soc, rad_lon_soc, soc_lw_mode,  &
            tg_tmp_soc, q_soc, ozone_soc, co2_soc, t_surf_for_soc, p_full_soc, p_half_soc, z_full_soc, z_half_soc, albedo_soc, n_profile, n_layer,     &
            output_heating_rate_sw, output_soc_flux_down_sw, output_soc_flux_up_lw, output_net_surf_sw_down, output_surf_lw_down, coszen )

       tg_tmp_soc = tg_tmp_soc + output_heating_rate_sw*delta_t
       net_surf_sw_down(:,:) = REAL(output_net_surf_sw_down(:,:))
       thd_sw_flux_down = REAL(output_soc_flux_down_sw)


       temp_tend(:,:,:) = temp_tend(:,:,:) + real(output_heating_rate_sw)
       
       output_heating_rate_total = output_heating_rate_lw + output_heating_rate_sw

       if(store_intermediate_rad)then
            tdt_soc_lw_store = output_heating_rate_lw
            tdt_soc_sw_store = output_heating_rate_sw
            if (id_soc_flux_down_sw > 0) then
                thd_sw_flux_down_store = thd_sw_flux_down
            endif

            if (id_soc_flux_up_lw > 0) then
                thd_lw_flux_up_store = thd_lw_flux_up
            endif

            if (id_coszen > 0) then
                coszen_store = coszen
            endif


            net_surf_sw_down_store = real(net_surf_sw_down, kind(r_def))
            surf_lw_down_store = real(surf_lw_down, kind(r_def))
       endif

    !Sending total heating rates
    ! Send diagnostics
    if(id_soc_heating_lw > 0) then
      used = send_data ( id_soc_heating_lw, output_heating_rate_lw, Time_diag)
    endif

    if(id_soc_flux_up_lw > 0) then
      used = send_data ( id_soc_flux_up_lw, thd_lw_flux_up, Time_diag)
    endif

    if(id_soc_heating_sw > 0) then
      used = send_data ( id_soc_heating_sw, output_heating_rate_sw, Time_diag)
    endif

    if(id_soc_flux_down_sw > 0) then
      used = send_data ( id_soc_flux_down_sw, thd_sw_flux_down, Time_diag)
    endif

    if(id_soc_heating_rate > 0) then
      used = send_data ( id_soc_heating_rate, output_heating_rate_total, Time_diag)
    endif

    if(id_coszen > 0) then
      used = send_data ( id_coszen, coszen, Time_diag)
    endif

    if(id_interp_temp > 0) then
      used = send_data ( id_interp_temp, t_half_out, Time_diag)
    endif

    if(id_norm_temp > 0) then
      used = send_data ( id_norm_temp, temp_in, Time_diag)
    endif

end subroutine run_socrates  

subroutine run_socrates_end

    use interpolator_mod, only: interpolator_end
    USE socrates_config_mod    

    if(do_read_ozone) call interpolator_end(o3_interp)
    if(do_read_co2)   call interpolator_end(co2_interp)
    

end subroutine run_socrates_end

!*****************************************************************************************
        subroutine interp_temp(z_full,z_half,temp_in, t_half)
          implicit none

          real(r_def),dimension(:,:),intent(in)  :: z_full,z_half,temp_in
          real(r_def),dimension(size(z_half,1), size(z_half,2)),intent(out) :: t_half

          integer i,k,kend
          real dzk,dzk1,dzk2
          
! note: z_full(kend) = z_half(kend), so there's something fishy
! also, for some reason, z_half(k=1)=0. so we need to deal with k=1 separately
          kend=size(z_full,2)
          do k=2,kend
                do i=1,size(temp_in,1)
                   dzk2 = 1./( z_full(i,k-1)   - z_full(i,k) )
                   dzk  = ( z_half(i,k  )   - z_full(i,k) )*dzk2
                   dzk1 = ( z_full(i,k-1)   - z_half(i,k) )*dzk2 
                   t_half(i,k) = temp_in(i,k)*dzk1 + temp_in(i,k-1)*dzk
                enddo
          enddo
! top of the atmosphere: need to extrapolate. z_half(1)=0, so need to use values
! on full grid
             do i=1,size(temp_in,1)
                !standard linear extrapolation
                !top: use full points, and distance is 1.5 from k=2
                t_half(i,1) = 0.5*(3*temp_in(i,1)-temp_in(i,2))
                !bottom: z=0 => distance is
                !-z_full(kend-1)/(z_full(kend)-z_full(kend-1))
                t_half(i,kend+1) = temp_in(i,kend-1) &
                     + (z_half(i,kend+1) - z_full(i,kend-1))&
                     * (temp_in(i,kend ) - temp_in(i,kend-1))&
                     / (z_full(i,kend  ) - z_full(i,kend-1))
             enddo


        end subroutine interp_temp
!*****************************************************************************************
  
end module socrates_interface_mod
