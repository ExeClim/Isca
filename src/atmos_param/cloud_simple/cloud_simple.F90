module cloud_simple_mod

#ifdef INTERNAL_FILE_NML
    use mpp_mod, only: input_nml_file
#else
    use fms_mod, only: open_namelist_file, close_file
#endif

  use            fms_mod, only: stdlog, FATAL, WARNING, error_mesg
  use   time_manager_mod, only: time_type
  use sat_vapor_pres_mod, only:  compute_qs

  use diag_manager_mod, only: register_diag_field, send_data

  implicit none

  logical ::   do_init = .true.  ! Check if init needs to be run

  real    ::   simple_cca =  0.0

  ! Critical RH (fraction) values - spookie protocol 1 only
  real    ::   rhc_sfc     = 0.95
  real    ::   rhc_base    = 0.7
  real    ::   rhc_top     = 0.3
  ! 
  real    ::   rhmsfc     = 0.95
  real    ::   rhm700     = 0.7
  real    ::   rhm200     = 0.3

  ! Critical RH (fraction) values - spookie protocol 2 only
  real    ::   rh_min_top  = 0.9
  real    ::   rh_min_sfc  = 1.0
  real    ::   rh_min_base = 0.8
  real    ::   rh_max_top  = 1.0
  real    ::   rh_max_sfc  = 1.0
  real    ::   rh_max_base = 1.0

  ! Pressure (Pa.) at cloud bottom and top (very approx)
  real    ::   p_base     = 70000.
  real    ::   p_top      = 20000.

  integer ::  spookie_protocol = 2

  namelist /cloud_simple_nml/  simple_cca, rhc_sfc, rhc_base, rhc_top, &
                               rhmsfc, rhm700, rhm200,                 &
                               rh_min_top, rh_min_sfc, rh_min_base,    &
                               rh_max_top, rh_max_sfc, rh_max_base
      
  real :: zerodegc = 273.15

  integer :: id_cf, id_reff_rad, id_frac_liq, id_qcl_rad, id_rh_in_cf, &
             id_simple_rhcrit, id_rh_min
 
  character(len=14), parameter ::   mod_name_cld = "cloud_simple"

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

   id_simple_rhcrit = &
      register_diag_field ( mod_name_cld, 'simple_rhcrit', axes(1:3), Time, &
      'RH as a percent for spookie protocol 1', &
      '%')

   id_rh_min = &
      register_diag_field ( mod_name_cld, 'rh_min', axes(1:3), Time, &
      'RH as a percent for spookie protocol 2', &
      '%')


    do_init = .false.  !initialisation completed
    
  end subroutine cloud_simple_init

  !-----------------------------------------------

  subroutine cloud_simple(p_half, p_full, Time,   &
                      temp, q_hum,                &
                      ! outs 
                      cf, cca, reff_rad, qcl_rad) 

    real       , intent(in), dimension(:,:,:)  :: temp, q_hum, p_full, p_half
    type(time_type) , intent(in)               :: Time

    real       , intent(out), dimension(:,:,:) :: cf, reff_rad, qcl_rad, cca

    real, dimension(size(temp,1), size(temp, 2), size(temp, 3)) :: qs, frac_liq, rh_in_cf, simple_rhcrit, rh_min, rh_max

    integer :: i, j, k, k_surf

    logical :: es_over_liq_and_ice

    !check initiation has been done - ie read in parameters
    if (do_init) call error_mesg ('cloud_simple',  &
         'cloud_simple_init has not been called.', FATAL)
    
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

          if (spookie_protocol .eq. 1) then
              call calc_rhcrit(p_full(i,j,k), p_full(i,j,k_surf), simple_rhcrit(i,j,k))
          else
              call calc_rh_min_max(p_full(i,j,k), p_full(i,j,k_surf),rh_min(i,j,k), rh_max(i,j,k))
          endif

          call calc_cf(q_hum(i,j,k), qs(i,j,k), cf(i,j,k), cca(i,j,k), rh_in_cf(i,j,k), &
                       simple_rhcrit = simple_rhcrit(i,j,k),                &
                       rh_min = rh_min(i,j,k), rh_max = rh_max(i,j,k) )

          call calc_qcl_rad(p_full(i,j,k), cf(i,j,k), temp(i,j,k), qcl_rad(i,j,k) )
        end do
      end do
    end do

  !save some diagnotics
  call output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, simple_rhcrit, rh_min, Time )  

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
    !  get the RH needed as a threshold for the cloud fraction calc.
    ! This is only requires for spookie_protocol=1
    real, intent(in)  :: p_full, p_surf
    real, intent(out) :: simple_rhcrit

    ! Calculate RHcrit as function of pressure
    if (p_full > p_base  ) then

      simple_rhcrit = rhc_sfc - ( rhc_sfc - rhc_base ) *           &
                     ( p_surf - p_full  ) / ( p_surf - p_base )

    else if ( p_full > p_top ) then 

      simple_rhcrit = rhc_base - ( rhc_base - rhc_top ) *        &
                      (p_base - p_full) / (p_base - p_top)
          
    else
      simple_rhcrit = rhc_top
    endif

  end subroutine calc_rhcrit

  subroutine calc_rh_min_max(p_full, p_surf, rh_min, rh_max)

    real, intent(in)  :: p_full, p_surf
    real, intent(out) :: rh_min, rh_max


    if (p_full > p_base  ) then !surface up to base

      rh_min = rh_min_sfc - ( rh_min_sfc - rh_min_base ) * ( p_surf - p_full  ) / ( p_surf - p_base )
      rh_max = rh_max_sfc - ( rh_max_sfc - rh_max_base ) * ( p_surf - p_full  ) / ( p_surf - p_base )

    else if ( p_full > p_top ) then  ! base up to top

      rh_min = rh_min_base - ( rh_min_base - rh_min_top ) * (p_base - p_full) / (p_base - p_top)
      rh_max = rh_max_base - ( rh_max_base - rh_max_top ) * (p_base - p_full) / (p_base - p_top)
          
    else ! above top

      rh_min = rh_min_top
      rh_max = rh_max_top

    endif


  end subroutine calc_rh_min_max

  subroutine calc_cf(q_hum, qsat, cf, cca, rh, simple_rhcrit, rh_min, rh_max)
    ! Calculate LS (stratiform) cloud fraction 
    ! as a simple linear function of RH

    real, intent(in)            :: q_hum, qsat
    real, intent(in), optional  :: simple_rhcrit, rh_min, rh_max

    real, intent(out) :: cf, rh, cca

    rh = q_hum/qsat

    if (spookie_protocol .eq. 1) then
        cf = (rh - simple_rhcrit ) / ( 1.0 - simple_rhcrit ) 

    else
        cf = (rh - rh_min ) / ( rh_max - rh_min ) 
    end if 

    cf = MAX( 0.0, MIN( 1.0, cf))

    ! include simple convective cloud fraction where present (not currenly used)
    cca = 0.0 ! no convective cloud fraction is calculated
             ! left in for future use

    !cca can not be used in simple clouds as in read_control  
    ! control%i_cloud_representation = ip_cloud_ice_water

    if (cca > 0.0) then
       cf = MAX( simple_cca, cf )
    end if


  end subroutine calc_cf

  subroutine calc_qcl_rad(p_full, cf, temp, qcl_rad)
    ! calculate cloud water content

    real , intent(in)   :: p_full, cf, temp
    real , intent(out)  :: qcl_rad

    real :: in_cloud_qcl 
 
    IF (spookie_protocol .eq. 1) THEN
        ! pressure dependent in_cloud_qcl
        in_cloud_qcl = 3.0e-4 + (1.0-3.0e-4)*(p_full-p_top)/80000.0
        in_cloud_qcl = MAX ( 0.0, in_cloud_qcl/1000.0 ) ! convert to kg/kg
        qcl_rad = cf * in_cloud_qcl
    ELSE
        ! temperatue dependent in_cloud_qcl
        in_cloud_qcl = MIN(0.2, 0.2 * ( temp - 220. ) / ( 280. -220. ))
        in_cloud_qcl = MAX (3.0e-4, in_cloud_qcl/1000.0 ) ! convert to kg/kg
        qcl_rad = cf * in_cloud_qcl
    ENDIF

  end subroutine calc_qcl_rad

  subroutine output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, simple_rhcrit, rh_min, Time)

    real, intent(in), dimension(:,:,:) :: cf, reff_rad, frac_liq, qcl_rad, rh_in_cf
    real, intent(in), dimension(:,:,:), optional :: simple_rhcrit, rh_min

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

    if ( id_rh_min > 0 ) then
      used = send_data ( id_rh_min, rh_min*100.0, Time)
    endif

  end subroutine output_cloud_diags

  subroutine cloud_simple_end ()

  ! If alloocated are added in init then deallocate them here.

  end subroutine cloud_simple_end

  !-----------------------------------------------

end module cloud_simple_mod
