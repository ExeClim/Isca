module cloud_simple_mod

  use            fms_mod, only: stdlog, FATAL, WARNING, error_mesg,  &
                                open_namelist_file, close_file, open_file, &
                                check_nml_error, mpp_pe
  use   time_manager_mod, only: time_type
  use sat_vapor_pres_mod, only: compute_qs
  use      constants_mod, only: KELVIN

  use diag_manager_mod, only: register_diag_field, send_data

  implicit none

  character(len=128) :: version='$Id: cloud_simple.F90,v 1.0 2021/05/11$'
  character(len=128) :: tag='Simple cloud scheme'

  logical ::   do_init = .true. ! update to false after init has been run

  real    ::   cca_lower_limit =  0.0 ! simple convective cloud fraction min
                                      ! the default is zero. Not being used but 
                                      ! could be adapted in future.

  ! There are two testing scenarios developed for SPOOKIE-2 project. Which
  ! is why this code was created. The first was implemented and found not to 
  ! work very well in other models. It has been implemented her for curiousity 
  ! testing but not tested very much yet. The second protocol is what was used 
  ! for the SPOOKIE-2 runs and has been more widely tested. For this reason,
  ! the default spookie_protocol is 2. The 1st protocol could be removed in 
  ! time (once the SPOOKIE-2 project is complete, but is left in the code for 
  ! now in case it is needed.

  integer ::  spookie_protocol = 2 ! default is 2

  ! Critical RH (fraction) values - SPOOKIE-2 protocol version 1
  real    ::   rhc_sfc     = 1.0
  real    ::   rhc_base    = 0.7
  real    ::   rhc_top     = 0.2 ! In the protocol this was 20 % and in 
                                 ! implementation it was 30%. To check in 
                                 ! next round of validation.

  ! Critical RH (fraction) values - SPOOKIE-2 protocol version 2
  ! initial values for RH. Updated in calc_rh_min_max
  real    ::   rh_min_sfc  = 1.0
  real    ::   rh_min_base = 0.8
  real    ::   rh_min_top  = 0.9

  real    ::   rh_max_sfc  = 1.0
  real    ::   rh_max_base = 1.0
  real    ::   rh_max_top  = 1.0

  ! Pressure (Pa) at cloud bottom and top (very approx)
  real    ::   p_base     = 70000.
  real    ::   p_top      = 20000.

  namelist /cloud_simple_nml/  cca_lower_limit, spookie_protocol,      &
                               rhc_sfc, rhc_base, rhc_top,             &
                               rh_min_top, rh_min_sfc, rh_min_base,    &
                               rh_max_top, rh_max_sfc, rh_max_base
      
  integer :: id_cf, id_reff_rad, id_frac_liq, id_qcl_rad, id_rh_in_cf, &
             id_simple_rhcrit, id_rh_min
 
  character(len=14), parameter ::   mod_name_cld = "cloud_simple"

  contains

  !-----------------------------------------------


  subroutine cloud_simple_init (axes, Time)

    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes

    integer :: io, ierr, unit

    unit = open_file (file='input.nml', action='read')
    ierr=1
    do while (ierr /= 0)
        read  (unit, nml=cloud_simple_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'cloud_simple_nml')
    enddo
    10 call close_file (unit)

    unit = open_file (file='logfile.out', action='append')
    if ( mpp_pe() == 0 ) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=cloud_simple_nml)
    endif
    call close_file(unit)

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

   ! rh_in_cf is an output diagnostic only for debugging  

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

    real       , intent(inout), dimension(:,:,:) :: cf, cca, reff_rad, qcl_rad

    real, dimension(size(temp,1), size(temp, 2), size(temp, 3)) :: qs, frac_liq
    real, dimension(size(temp,1), size(temp, 2), size(temp, 3)) :: rh_in_cf
    real, dimension(size(temp,1), size(temp, 2), size(temp, 3)) :: simple_rhcrit
    real, dimension(size(temp,1), size(temp, 2), size(temp, 3)) :: rh_min,rh_max

    integer :: i, j, k, k_surf

    !check initiation has been done
    if (do_init) call error_mesg ('cloud_simple',  &
         'cloud_simple_init has not been called.', FATAL)

    ! Get the saturated specific humidity with respect to water and ice
    ! this is set by the namelist variable sat_vapor_pres_nml
    call compute_qs(temp, p_full, qs) 

    k_surf = size(temp, 3) !set the location of the lowest model level


    ! For future revisions, consider rewriting to remove the loops. 
    do k=1, size(temp, 3)
      do j=1, size(temp, 2)
        do i=1, size(temp, 1)

          ! calculate the liquid fraction
          call calc_liq_frac(temp(i,j,k), frac_liq(i,j,k))

          ! calculate the effective radius
          call calc_reff(frac_liq(i,j,k), reff_rad(i,j,k))

          if (spookie_protocol .eq. 1) then
              ! calculate the critical RH 
              call calc_rhcrit(p_full(i,j,k), p_full(i,j,k_surf), &
                               simple_rhcrit(i,j,k))
          else
              ! calculate the min and max RH 
              call calc_rh_min_max(p_full(i,j,k), p_full(i,j,k_surf), &
                                   rh_min(i,j,k), rh_max(i,j,k))
          endif

          ! calculate the cloud fraction
          call calc_cf(q_hum(i,j,k), qs(i,j,k), cf(i,j,k), cca(i,j,k),         &
                       rh_in_cf(i,j,k), simple_rhcrit = simple_rhcrit(i,j,k),  &
                       rh_min = rh_min(i,j,k), rh_max = rh_max(i,j,k) )

          ! calculate the specific humidity of cloud liquid
          call calc_mixing_ratio(p_full(i,j,k), cf(i,j,k), temp(i,j,k),        &
                                 qcl_rad(i,j,k) )
        end do
      end do
    end do

  !save some diagnotics
  call output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, &
                          simple_rhcrit, rh_min, Time )  

  end subroutine cloud_simple

  subroutine calc_liq_frac(temp, frac_liq)

    real, intent(in)    :: temp
    real, intent(out)   :: frac_liq

    if (temp > KELVIN) then
        ! All liquid if temp above zero
        frac_liq = 1.0
    else if (temp < KELVIN-40.0) then
        ! All ice if temp is below -40C
        frac_liq = 0.0
    else           
        ! linearly interpolate between T=0 and -40C
        frac_liq = 1.0 - (KELVIN - temp) / 40.0
    end if 

  end subroutine calc_liq_frac

  subroutine calc_reff(frac_liq, reff_rad)
    ! the effective cloud radius is bounded between 10 and 20 microns

    real, intent(in) :: frac_liq
    real, intent(out) :: reff_rad

    reff_rad =  10.0 * frac_liq + 20.0 * (1.0 - frac_liq)  
    ! units in microns this will be updated before passing into soc

  end subroutine calc_reff

  subroutine calc_rhcrit(p_full, p_surf, simple_rhcrit)  
    ! Get the RH needed as a threshold for the cloud fraction calc.
    ! This is only requires for spookie_protocol=1
    real, intent(in)  :: p_full, p_surf
    real, intent(out) :: simple_rhcrit

    ! Calculate RHcrit as function of pressure
    if (p_full > p_base) then

      simple_rhcrit = rhc_sfc - (rhc_sfc - rhc_base) *           &
                      (p_surf - p_full) / (p_surf - p_base)

    else if (p_full > p_top) then 

      simple_rhcrit = rhc_base - (rhc_base - rhc_top) *          &
                      (p_base - p_full) / (p_base - p_top)
          
    else
      simple_rhcrit = rhc_top
    endif

  end subroutine calc_rhcrit

  subroutine calc_rh_min_max(p_full, p_surf, rh_min, rh_max)

    real, intent(in)  :: p_full, p_surf
    real, intent(out) :: rh_min, rh_max

    real :: layer

    ! calculate RH min and max as a function of pressure

    if (p_full > p_base) then 
      ! For the layer between the surface and cloud base (default is 700 hpa)

      layer =  (p_surf - p_full) / (p_surf - p_base)

      ! correction step to update initial values
      rh_min = rh_min_sfc - (rh_min_sfc - rh_min_base) * layer
      rh_max = rh_max_sfc - (rh_max_sfc - rh_max_base) * layer

    else if ( p_full > p_top ) then  
      ! For the layer where the cloud is (base up to top)

      layer = (p_base - p_full) / (p_base - p_top)
      rh_min = rh_min_base - ( rh_min_base - rh_min_top ) * layer
      rh_max = rh_max_base - ( rh_max_base - rh_max_top ) * layer
          
    else 
      ! Above the cloud top above top
      rh_min = rh_min_top
      rh_max = rh_max_top
    endif

  end subroutine calc_rh_min_max

  subroutine calc_cf(q_hum, qsat, cf, cca, rh, simple_rhcrit, rh_min, rh_max)
    ! Calculate large scale (stratiform) cloud fraction 
    ! as a simple linear function of RH

    real, intent(in)            :: q_hum, qsat
    real, intent(in), optional  :: simple_rhcrit, rh_min, rh_max

    real, intent(out) :: cf, rh, cca

    ! The environment RH
    rh = q_hum / qsat

    if (spookie_protocol .eq. 1) then
        cf = (rh - simple_rhcrit ) / (1.0 - simple_rhcrit)
    else
        cf = (rh - rh_min) / (rh_max - rh_min) 
    end if 

    cf = MAX(0.0, MIN(1.0, cf))

    ! include simple convective cloud fraction where present
    ! This is currently not being used and array are zeros as
    ! no convective cloud fraction is calculated
    ! left in for future use

    !if (cca > 0.0) then
    !   cf = MAX( cca_lower_limit, cf )
    !end if

  end subroutine calc_cf

  subroutine calc_mixing_ratio(p_full, cf, temp, qcl_rad)

    ! calculate cloud water content

    real , intent(in)   :: p_full, cf, temp
    real , intent(out)  :: qcl_rad ! mixing ratio of cloud liquid

    real :: in_cloud_qcl 
 
    IF (spookie_protocol .eq. 1) THEN
        ! pressure dependent in_cloud_qcl
        ! bounded between:
        !     1 g/kg at 1000hpa
        !  3e-4 g/kg at 200 hpa
        in_cloud_qcl = 3.0e-4 + (1.0 - 3.0e-4) * (p_full - p_top) / 80000.0
        in_cloud_qcl = MAX ( 0.0, in_cloud_qcl) ! in g/kg
    ELSE
        ! temperatue dependent in_cloud_qcl 
        ! bounded between: 
        !     3e-4 g/kg at 220 K 
        !      0.2 g/kg at 280K
        in_cloud_qcl = MIN(0.2, 0.2 * ( temp - 220. ) / ( 280. - 220. ))
        in_cloud_qcl = MAX(3.0e-4, in_cloud_qcl) ! in g/kg
    ENDIF

    qcl_rad = cf * in_cloud_qcl / 1000. ! convert to kg/kg

  end subroutine calc_mixing_ratio

  subroutine output_cloud_diags(cf, reff_rad, frac_liq, qcl_rad, rh_in_cf, &
                                simple_rhcrit, rh_min, Time)

    real, intent(in), dimension(:,:,:) :: cf, reff_rad, frac_liq, qcl_rad, & 
                                          rh_in_cf
    real, intent(in), dimension(:,:,:), optional :: simple_rhcrit, rh_min

    type(time_type) , intent(in)       :: Time

    logical :: used

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

end module cloud_simple_mod
