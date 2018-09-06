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
  real    ::   rhcsfc     = 0.95
  real    ::   rhc700     = 0.7
  real    ::   rhc200     = 0.3
  real    ::   rhmsfc     = 0.95
  real    ::   rhm700     = 0.7
  real    ::   rhm200     = 0.3

  namelist /cloud_simple_nml/  simple_cca, rhcsfc, rhc700, rhc200, &
                                           rhmsfc, rhm700, rhm200
  real :: zerodegc = 273.15

  integer :: id_cf_rad, id_reff_rad, id_frac_liq, id_qcl_rad 
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
    id_cf_rad = &
      register_diag_field ( mod_name_cld, 'cf_rad', axes(1:3), Time, &
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

    do_init = .false.  !initialisation completed
    
  end subroutine cloud_simple_init

  !-----------------------------------------------

  subroutine cloud_simple (p_half, p_full, Time,   &
                      temp, q_hum,                 &
                      ! outs 
                      cf_rad, reff_rad, qcl_rad ) 

    real       , intent(in), dimension(:,:,:)  :: temp, q_hum, p_full, p_half
    type(time_type) , intent(in)               :: Time

    real       , intent(out), dimension(:,:,:) :: cf_rad, reff_rad, qcl_rad

    real, dimension(size(temp,1), size(temp, 2), size(temp, 3))    :: qs, frac_liq
    real  ::  simple_rhcrit

    integer :: i, j, k, k_surf

    logical :: es_over_liq_and_ice

          real :: tmp1, tmp2 !remove tmp1 and tmp2 after debugging

    !check initiation has been done - ie read in parameters
    if (do_init) call error_mesg ('cloud_simple',  &
         'cloud_simple_init has not been called.', FATAL)
    
    ! Get the saturated specific humidity TOTAL (ie ice and vap) ***double check maths!
    call compute_qs(temp, p_full, qs, es_over_liq_and_ice=.true.) !qs=qsat in um

    k_surf = size(temp, 3)

    do i=1, size(temp, 1)
      do j=1, size(temp, 2)
        do k=1, size(temp, 3)

          !caluclate the frac_liq 
          call calc_liq_frac(temp(i,j,k), frac_liq(i,j,k))
          call calc_reff(frac_liq(i,j,k), reff_rad(i,j,k))
          call calc_rhcrit(p_full(i,j,k), p_full(i,j,k_surf), simple_rhcrit)
          call calc_cf_rad(q_hum(i,j,k), qs(i,j,k), simple_rhcrit, cf_rad(i,j,k))
          call calc_qcl_rad(p_full(i,j,k), cf_rad(i,j,k), qcl_rad(i,j,k) )
        end do
      end do
    end do

  !save some diagnotics
  call output_cloud_diags(cf_rad, reff_rad, frac_liq, qcl_rad, Time )  

tmp2 = maxval(cf_rad)
tmp1 = maxval(reff_rad)
tmp2 = maxval(cf_rad)


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

    real, intent(in) :: frac_liq
    real, intent(out) :: reff_rad

    reff_rad =  10.0 * frac_liq + 20.0 * (1.0 - frac_liq)  !units in microns

  end subroutine calc_reff

  subroutine calc_rhcrit(p_full, p_surf, simple_rhcrit)   !need to check p_full - > p_layer_centres

    real, intent(in)  :: p_full, p_surf
    real, intent(out) :: simple_rhcrit

    ! Calculate RHcrit as function of pressure
    if (p_full > 70000.0  ) then

      simple_rhcrit = rhcsfc - ( rhcsfc - rhc700 ) *           &
                     ( p_surf - p_full  ) / ( p_surf - 70000.0 )

    else if ( p_full > 20000.0 ) then 

      simple_rhcrit = rhc700 - ( rhc700 - rhc200 ) *        &
                      ( 70000.0 - p_full) / 50000.0
          
    else
      simple_rhcrit = rhc200
    endif

  end subroutine calc_rhcrit

  subroutine calc_cf_rad (q_hum, qsat, simple_rhcrit, cf_rad)
    ! Calculate LS (stratiform) cloud fraction 
    ! as a simple linear function of RH

    real, intent(in)  :: q_hum, qsat, simple_rhcrit
    real, intent(out) :: cf_rad
 
    real :: rh, cca

    rh = q_hum/qsat
    cf_rad = MAX( 0.0, MIN( 1.0, ( rh - simple_rhcrit ) / ( 1.0 - simple_rhcrit ) ))

    ! include simple convective cloud fraction where present

    cca = 0.0 ! no convective cloud fraction is calculated
              ! left in for fture use

    if (cca > 0.0) then
      cf_rad = MAX( simple_cca, cf_rad )
    end if

    

  end subroutine calc_cf_rad

  subroutine calc_qcl_rad(p_full, cf_rad, qcl_rad)
    ! calculate simple water content

    real , intent(in)   :: p_full, cf_rad
    real , intent(out)   :: qcl_rad

    real :: in_cloud_qcl 

    in_cloud_qcl = 3.0e-4 +                                        &
           (1.0-3.0e-4)*(p_full-20000.0)/80000.0

    in_cloud_qcl = MAX ( 0.0, in_cloud_qcl/1000.0 ) ! convert to kg/kg

    qcl_rad = cf_rad * in_cloud_qcl

  end subroutine calc_qcl_rad



  subroutine output_cloud_diags(cf_rad, reff_rad, frac_liq, qcl_rad, Time)

    real, intent(in), dimension(:,:,:) :: cf_rad, reff_rad, frac_liq, qcl_rad

    type(time_type) , intent(in)       :: Time

    real :: used

    if ( id_cf_rad > 0 ) then
      used = send_data ( id_cf_rad, cf_rad, Time)
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


  end subroutine output_cloud_diags

  !-----------------------------------------------


  subroutine cloud_simple_end ()

  ! do we deallocate anything?

  end subroutine cloud_simple_end

  !-----------------------------------------------

end module cloud_simple_mod
