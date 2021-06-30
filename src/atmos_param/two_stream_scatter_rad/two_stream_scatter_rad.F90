!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Two stream grey radiative transfer code with scattering solved    !!
!! following Pierrehumbert (2010), Principles of Planerary Climate.  !!
!! Chapter 5, Page 353.                                              !!
!!                                                                   !!
!! Forumulates radiative transfer equations as matrix problem and    !!
!! solves via LU decomposition.                                      !!
!!                                                                   !!
!! Code that carries out LU decomposition and solves for a generic   !!
!! system Ax = B was downloaded at the following link:               !!
!! http://jean-pierre.moreau.pagesperso-orange.fr/f_matrices.html    !!
!!                                                                   !!
!! Please direct any questions to Neil Lewis, Univ. Oxford, at       !! 
!! neil.lewis <at> physics.ox.ac.uk                                  !!
!!                                                                   !!
!! Distributed under GNU general public license, available here:     !!
!! http://www.gnu.org/licenses/gpl.txt                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module two_stream_scatter_rad_mod

! ==================================================================================
! ==================================================================================

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file, error_mesg, &
                                    NOTE, FATAL,  uppercase

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks, pstd_mks_earth, seconds_per_sol, orbital_period

   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, length_of_year, length_of_day, &
                                    operator(+), operator(-), operator(/=), get_time

   use astronomy_mod,         only: astronomy_init, diurnal_solar

   use interpolator_mod,      only: interpolate_type, interpolator_init, interpolator, ZERO, interpolator_end

!==================================================================================
implicit none
private
!==================================================================================

! version information

character(len=128) :: version='$Id: two_stream_scatter_rad.F90,v 1.0 2019/01/13$'
character(len=128) :: tag='Two-stream semi-gray atmosphere with scattering'


! public interfaces

public :: two_stream_scatter_rad_init, two_stream_scatter_rad, two_stream_scatter_rad_end
!==================================================================================


! module variables
logical :: initialized     = .false.

! astronomical parameters
real    :: solar_constant  = 1360.0
real    :: del_sol         = 1.4
real    :: del_sw          = 0.0
logical :: do_seasonal     = .false.
integer :: solday          = -10  !s Day of year to run perpetually if do_seasonal=True and solday>0
real    :: equinox_day     = 0.75 !s Fraction of year [0,1] where NH autumn equinox occurs (only really useful if calendar has defined months).
logical :: use_time_average_coszen = .false. !s if .true., then time-averaging is done on coszen so that insolation doesn't depend on timestep
real    :: dt_rad_avg     = -1
real :: diabatic_acce = 1.0 ! artefact from gp scheme in two_stream_grey_rad 

! co2 
real :: carbon_conc = 360.0 ! ppmv

character(len=32) :: sw_optical_depth = 'generic'
integer :: sw_flag 
integer, parameter :: sw_TRANSPARENT = 1
integer, parameter :: sw_GENERIC = 2

! parameters for lw optical depths etc (defaults to BOG... i.e. no scattering)
real :: lw_sca_a = 0.0
real :: lw_sca_b = 0.0
real :: lw_sca_c = 0.0 
real :: lw_abs_a = 0.1627 * 9.80 / 101325.0 ! dtau = a * dsigma   --->   dtau = 1/g * c * dp   where c = ga/p_s 
real :: lw_abs_b = 1997.9 * 9.80 / 101325.0 
real :: lw_abs_c = 0.17   * 9.80 / 101325.0

! parameters for sw optical depths (defaults to BOG... i.e. no scattering or absorption)
real :: sw_sca_a = 0.0
real :: sw_sca_b = 0.0
real :: sw_sca_c = 0.0
real :: sw_abs_a = 0.0
real :: sw_abs_b = 0.0 

real :: gamma = 1.0 ! associated with closure for integrating over all angles 
real :: gammaprime = 1.0 ! " " " (see discussion in Pierrehumbert, 2010, Chapter 5)
real :: g_asym = 0.0 ! assymetry parameter 


real, allocatable, dimension(:,:)   :: insolation, p2,  sw_tau_0 !s albedo now defined in mixed_layer_init
real, allocatable, dimension(:,:)   :: pi_B_surf
real, allocatable, dimension(:,:,:) :: pi_B, tdt_rad, tdt_solar
real, allocatable, dimension(:,:,:) :: lw_up, lw_down, lw_flux, &
                                       sw_up, sw_down, sw_down_direct, sw_flux, &
                                       rad_flux
real, allocatable, dimension(:,:,:) :: lw_tau, sw_tau, lw_del_tau, sw_del_tau!, &
                                       !lw_tau_dummy, sw_tau_dummy
real, allocatable, dimension(:,:,:) :: lw_ss_albedo, sw_ss_albedo
real, allocatable, dimension(:,:,:) :: lw_abs_coeff, lw_scatter_coeff, &
                                       sw_abs_coeff, sw_scatter_coeff
real, allocatable, dimension(:,:,:) :: sw_gammab, sw_gammaone, sw_gammatwo, &
                                       sw_gammaminus, sw_gammaplus 
real, allocatable, dimension(:,:,:) :: lw_gammab, lw_gammaone, lw_gammatwo, &
                                       lw_gammaminus, lw_gammaplus 
real, allocatable, dimension(:,:)   :: olr, net_lw_surf, toa_sw_in, &
                                       toa_sw_out, coszen, fracsun


real, save :: pi, deg_to_rad , rad_to_deg

!extras for reading in co2 concentration
logical                             :: do_read_co2=.false.
type(interpolate_type),save         :: co2_interp           ! use external file for co2
character(len=256)                  :: co2_file='co2'       !  file name of co2 file to read
character(len=256)                  :: co2_variable_name='co2'       !  file name of co2 file to read


namelist/two_stream_scatter_rad_nml/ solar_constant, del_sol, del_sw, &
           do_seasonal, solday, equinox_day,  &
           use_time_average_coszen, dt_rad_avg,&
           diabatic_acce,& !Schneider Liu values 
           lw_abs_a, lw_abs_b, lw_abs_c, &
           lw_sca_a, lw_sca_b, lw_sca_c, &
           sw_abs_a, sw_abs_b, &
           sw_sca_a, sw_sca_b, sw_sca_c, &
           gamma, gammaprime, g_asym, &
           do_read_co2, co2_file, co2_variable_name, carbon_conc, &
           sw_optical_depth

!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_swup_toa, &
           id_net_lw_surf, id_lwdn_sfc, &
           id_lwup_sfc, id_tdt_rad, id_tdt_solar, id_flux_rad, id_flux_lw, &
           id_flux_sw, id_coszen, id_fracsun, &
           id_flux_lw_up, id_flux_lw_down, id_flux_sw_up, id_flux_sw_down, &
           id_flux_sw_down_direct, id_co2

character(len=18), parameter :: mod_name = 'two_stream_scatter'

real :: missing_value = -999.


contains


! ==================================================================================
! ==================================================================================


subroutine two_stream_scatter_rad_init(is, ie, js, je, num_levels, axes, Time, lonb, latb, dt_real)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real ,dimension(:,:),intent(in),optional :: lonb,latb !s Changed to 2d arrays as 2013 interpolator expects this.
real, intent(in)                  :: dt_real !s atmospheric timestep, used for radiation averaging
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=two_stream_scatter_rad_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'two_stream_scatter_rad_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=two_stream_scatter_rad_nml)
endif
call close_file (unit)

pi         = 4. * atan(1.)
deg_to_rad = pi/180.
rad_to_deg = 180./pi

call astronomy_init

if(dt_rad_avg .le. 0) dt_rad_avg = dt_real !s if dt_rad_avg is set to a value in nml then it will be used instead of dt_real

if(do_read_co2)then
   call interpolator_init (co2_interp, trim(co2_file)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
endif

if (uppercase(trim(sw_optical_depth)) == 'TRANSPARENT') then 
  sw_flag = sw_TRANSPARENT
elseif (uppercase(trim(sw_optical_depth)) == 'GENERIC') then
  sw_flag = sw_GENERIC
else
  call error_mesg('two_stream_scatter_rad','"'//trim(sw_optical_depth)//'"'//' is not a valid sw parametrisation.', FATAL)
endif

initialized = .true.

allocate (sw_tau_0         (ie-is+1, je-js+1))
allocate (insolation       (ie-is+1, je-js+1))
allocate (p2               (ie-is+1, je-js+1))

allocate (pi_B_surf        (ie-is+1, je-js+1))

allocate (pi_B             (ie-is+1, je-js+1, num_levels))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_solar        (ie-is+1, je-js+1, num_levels))

allocate (lw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (lw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (lw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down_direct   (ie-is+1, je-js+1, num_levels+1))
allocate (sw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (rad_flux         (ie-is+1, je-js+1, num_levels+1))

allocate (lw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (sw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (lw_del_tau       (ie-is+1, je-js+1, num_levels))
allocate (sw_del_tau       (ie-is+1, je-js+1, num_levels))

allocate (lw_ss_albedo     (ie-is+1, je-js+1, num_levels))
allocate (sw_ss_albedo     (ie-is+1, je-js+1, num_levels))

allocate (lw_abs_coeff     (ie-is+1, je-js+1, num_levels))
allocate (sw_abs_coeff     (ie-is+1, je-js+1, num_levels))
allocate (lw_scatter_coeff (ie-is+1, je-js+1, num_levels))
allocate (sw_scatter_coeff (ie-is+1, je-js+1, num_levels))

allocate (sw_gammab        (ie-is+1, je-js+1, num_levels))
allocate (sw_gammaminus    (ie-is+1, je-js+1, num_levels))
allocate (sw_gammaplus     (ie-is+1, je-js+1, num_levels))
allocate (sw_gammaone      (ie-is+1, je-js+1, num_levels))
allocate (sw_gammatwo      (ie-is+1, je-js+1, num_levels))

allocate (lw_gammab        (ie-is+1, je-js+1, num_levels))
allocate (lw_gammaminus    (ie-is+1, je-js+1, num_levels))
allocate (lw_gammaplus     (ie-is+1, je-js+1, num_levels))
allocate (lw_gammaone      (ie-is+1, je-js+1, num_levels))
allocate (lw_gammatwo      (ie-is+1, je-js+1, num_levels))

allocate (olr              (ie-is+1, je-js+1))
allocate (net_lw_surf      (ie-is+1, je-js+1))
allocate (toa_sw_in        (ie-is+1, je-js+1))
allocate (toa_sw_out       (ie-is+1, je-js+1))
allocate (coszen           (ie-is+1, je-js+1))
allocate (fracsun          (ie-is+1, je-js+1)) !jp from astronomy.f90 : fraction of sun on surface




!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    id_olr = &
    register_diag_field ( mod_name, 'olr', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'W/m2', missing_value=missing_value               )
    id_swdn_sfc = &
    register_diag_field ( mod_name, 'swdn_sfc', axes(1:2), Time, &
               'Absorbed SW at surface', &
               'W/m2', missing_value=missing_value               )
    id_swdn_toa = &
    register_diag_field ( mod_name, 'swdn_toa', axes(1:2), Time, &
               'SW flux down at TOA', &
               'W/m2', missing_value=missing_value               )
    id_swup_toa = &
    register_diag_field ( mod_name, 'swup_toa', axes(1:2), Time, &
                'SW flux up at TOA', &
                'W/m2', missing_value=missing_value               )
    id_lwup_sfc = &
    register_diag_field ( mod_name, 'lwup_sfc', axes(1:2), Time, &
               'LW flux up at surface', &
               'W/m2', missing_value=missing_value               )

    id_lwdn_sfc = &
    register_diag_field ( mod_name, 'lwdn_sfc', axes(1:2), Time, &
               'LW flux down at surface', &
               'W/m2', missing_value=missing_value               )

    id_net_lw_surf = &
    register_diag_field ( mod_name, 'net_lw_surf', axes(1:2), Time, &
               'Net upward LW flux at surface', &
               'W/m2', missing_value=missing_value               )

    id_tdt_rad = &
        register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
               'Temperature tendency due to radiation', &
               'K/s', missing_value=missing_value               )

    id_tdt_solar = &
        register_diag_field ( mod_name, 'tdt_solar', axes(1:3), Time, &
               'Temperature tendency due to solar radiation', &
               'K/s', missing_value=missing_value               )

    id_flux_rad = &
        register_diag_field ( mod_name, 'flux_rad', axes(half), Time, &
               'Total radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw = &
        register_diag_field ( mod_name, 'flux_lw', axes(half), Time, &
               'Net longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw = &
        register_diag_field ( mod_name, 'flux_sw', axes(half), Time, &
               'Net shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )

    id_flux_lw_up = &
        register_diag_field ( mod_name, 'flux_lw_up', axes(half), Time, &
               'Upward longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw_down = &
        register_diag_field ( mod_name, 'flux_lw_down', axes(half), Time, &
               'Downward longwave radiative flux (positive down)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw_up = &
        register_diag_field ( mod_name, 'flux_sw_up', axes(half), Time, &
               'Upward shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw_down = &
        register_diag_field ( mod_name, 'flux_sw_down', axes(half), Time, &
               'Downward (diffuse) shortwave radiative flux (positive down)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw_down_direct = &
        register_diag_field ( mod_name, 'flux_sw_down_direct', axes(half), Time, &
               'Downward direct shortwave radiative flux (positive down)', &
               'W/m^2', missing_value=missing_value               )

    id_coszen  = &
               register_diag_field ( mod_name, 'coszen', axes(1:2), Time, &
                 'cosine of zenith angle', &
                 'none', missing_value=missing_value      )
    id_fracsun  = &
               register_diag_field ( mod_name, 'fracsun', axes(1:2), Time, &
                 'daylight fraction of time interval', &
                 'none', missing_value=missing_value      )

    id_co2  = &
               register_diag_field ( mod_name, 'co2', Time, &
                 'co2 concentration', &
                 'ppmv', missing_value=missing_value      )

return
end subroutine two_stream_scatter_rad_init

! ==================================================================================

subroutine two_stream_scatter_rad (is, js, ie, je, Time_diag, lat, lon, &
                                   p_half, p_full, &
                                   t_surf, t, tdt, net_surf_sw_down, surf_lw_down, albedo, q)


integer, intent(in)                 :: is, js, ie, je
type(time_type), intent(in)         :: Time_diag
real, intent(in), dimension(:,:)    :: lat, lon, albedo
real, intent(out), dimension(:,:)   :: net_surf_sw_down
real, intent(out), dimension(:,:)   :: surf_lw_down
real, intent(in), dimension(:,:,:)  :: t, q,  p_half, p_full
real, intent(in), dimension(:,:)    :: t_surf
real, intent(inout), dimension(:,:,:) :: tdt
integer :: i, j, k, n, dyofyr

integer :: seconds, year_in_s, days
real :: r_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, day_in_s, r_solday, r_total_seconds, r_days, r_dt_rad_avg, dt_rad_radians
logical :: used


real ,dimension(size(q,1),size(q,2),size(q,3)) :: co2f

n = size(t,3)

if(do_read_co2)then
  call interpolator( co2_interp, Time_diag, p_half, co2f, trim(co2_variable_name))
  carbon_conc = maxval(co2f) !Needs maxval just because co2f is a 3d array of a constant, so maxval is just a way to pick out one number
endif

! =================================================================================
! COMPUTE INSOLATION 

! insolation at TOA
if (do_seasonal) then
  ! Seasonal Cycle: Use astronomical parameters to calculate insolation
  call get_time(Time_diag, seconds, days)
  call get_time(length_of_year(), year_in_s)
  r_seconds = real(seconds)
  day_in_s = length_of_day()
  frac_of_day = r_seconds / day_in_s

  if(solday .ge. 0) then
      r_solday=real(solday)
      frac_of_year = (r_solday*day_in_s) / year_in_s
  else
      r_days=real(days)
      r_total_seconds=r_seconds+(r_days*day_in_s)
      frac_of_year = r_total_seconds / year_in_s
  endif

  gmt = abs(mod(frac_of_day, 1.0)) * 2.0 * pi

  time_since_ae = modulo(frac_of_year-equinox_day, 1.0) * 2.0 * pi

  if(use_time_average_coszen) then

     r_dt_rad_avg=real(dt_rad_avg)
     dt_rad_radians = (r_dt_rad_avg/day_in_s)*2.0*pi

     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun, dt_rad_radians)
  else
     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun)
  end if

  insolation = solar_constant * coszen

else
  ! Default: Averaged Earth insolation at all longitudes
  p2          = (1. - 3.*sin(lat)**2)/4.
  insolation  = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * sin(lat))
  coszen      = insolation / solar_constant
end if


! =================================================================================
! COMPUTE SHORTWAVE RADIATION 

! zero fluxes before starting 
sw_down_direct = 0.0
sw_up          = 0.0
sw_down        = 0.0
sw_scatter_coeff = 0.0
sw_abs_coeff     = 0.0
sw_ss_albedo     = 0.0 


select case (sw_flag)
case(sw_TRANSPARENT)
  ! just reflect at surface (direct down is reflected into diffuse up)
  do k = 1, n+1
    sw_down_direct(:,:,k) = insolation(:,:)
  enddo 
  do k = 1, n+1
    sw_up(:,:,k) = albedo(:,:) * sw_down_direct(:,:,n+1)
  enddo 

case(sw_GENERIC)
  ! Absorption and scattering in shortwave applied to both the upward and 
  ! downward beams. Total SW flux is made up of three components: a downward 
  ! direct component, and upward and downward diffuse components. Direct 
  ! radiation impinging on the surface is reflected into the upward diffuse 
  ! beam. 
  ! 
  ! First order differencing is applied to derivatives w.r.t tau. 
  !
  ! Once optical depths have been defined, the problem is solved as a matrix 
  ! equation of the form Ax = B, where x represents the upward and downward 
  ! fluxes (combined into a single vector). A is a banded matrix, and the problem 
  ! is solved by LU decomposition. This is carried out by the routine 
  ! two_stream_scatter_solve. 
  !
  ! Method follows Pierrehumbert (2010): Princples of Planetary Climate. Chapter 
  ! 5, Page 353. 

  
  ! Compute optical properties of atmosphere 
  ! Compute scattering coefficient (currently no contribution from carbon_conc, easily added though...)
  sw_scatter_coeff  = sw_sca_a + sw_sca_b * q 
  ! Compute absorption coefficient 
  sw_abs_coeff = sw_abs_a + sw_abs_b*q 
  ! Compute single scattering albedo 
  where((sw_abs_coeff + sw_scatter_coeff).ne.0.0) 
    sw_ss_albedo = sw_scatter_coeff / (sw_abs_coeff + sw_scatter_coeff)
  endwhere 
  ! Compute SW del_tau 
  sw_del_tau =  1 / grav * (sw_abs_coeff + sw_scatter_coeff) * & 
               (p_half(:,:,1:n) - p_half(:,:,2:n+1) )


  ! Compute SW optical depth 
  ! define optical depth to be zero at surface, following Pierrehumbert(2010)
  sw_tau(:,:,n+1) = 0.0 
  do k = 1, n
    sw_tau(:,:,n+1-k) = sw_tau(:,:,n+2-k) - sw_del_tau(:,:,n+1-k)
  enddo 

  !! gammas (as in Pierrehumbert, 2010)  
  sw_gammaone = gamma *(1.0 - 3./2.*g_asym*sw_ss_albedo) + gammaprime*(1-sw_ss_albedo)
  sw_gammatwo = gamma*(1.0 - 3./2.*g_asym*sw_ss_albedo) - gammaprime*(1-sw_ss_albedo)
  sw_gammab = sw_gammaone - sw_gammatwo 
  do k = 1, n
    sw_gammaplus(:,:,k)  = 0.5*sw_ss_albedo(:,:,k) - gamma*sw_ss_albedo(:,:,k)*3./2.*g_asym*coszen(:,:) 
    sw_gammaminus(:,:,k) = 0.5*sw_ss_albedo(:,:,k) + gamma*sw_ss_albedo(:,:,k)*3./2.*g_asym*coszen(:,:) 
  enddo 

  ! Compute direct downward radiation 
  sw_down_direct(:,:,1) = insolation(:,:)
  do k = 1, n 
    where ((coszen).ne.0.0) 
    sw_down_direct(:,:,k+1) =  sw_down_direct(:,:,k) * &
                             (2*coszen(:,:) + sw_del_tau(:,:,k)) / & 
                             (2*coszen(:,:) - sw_del_tau(:,:,k)) ! << Will this always be non-zero? 
    endwhere 
  enddo 
  
  !!!! Zero blackbody radiation for shortwave !!!!
  pi_B_surf = 0.0
  pi_B = 0.0


  do i = 1, (ie-is+1)
    do j = 1, (je-js+1)
      if (coszen(i,j) .eq. 0.0) then 
        sw_up(i,j,:) = 0.0
        sw_down(i,j,:) = 0.0
      else 
      call two_stream_solver(n, albedo(i,j), coszen(i,j), pi_B_surf(i,j), &
                           pi_B(i,j,:), sw_down_direct(i,j,:), &
                           sw_del_tau(i,j,:)/2, sw_gammaone(i,j,:), &
                           sw_gammatwo(i,j,:), sw_gammab(i,j,:), &
                           sw_gammaplus(i,j,:), sw_gammaminus(i,j,:), &
                           sw_up(i,j,:), sw_down(i,j,:)) ! zero contribution from in-atmosphere emission (see above)
                           ! factor 1/2 in optical depth comes from closure after integrating over all angles (only relevant for diffuse beam)
      endif 
    
    enddo
  enddo 

end select 


! =================================================================================
! LONGWAVE RADIATION
! Supports absorption and scattering in the longwave. 
! 
! First order differencing is applied to derivatives w.r.t tau. 
!
! Once optical depths have been defined, the problem is solved as a matrix 
! equation of the form Ax = B, where x represents the upward and downward 
! fluxes (combined into a single vector). A is a banded matrix, and the problem 
! is solved by LU decomposition. This is carried out by the routine 
! two_stream_scatter_solve. 
!
! Method follows Pierrehumbert (2010): Princples of Planetary Climate. Chapter 
! 5, Page 353. 

! zero fluxes before starting 
lw_up            = 0.0
lw_down          = 0.0 
lw_scatter_coeff = 0.0
lw_abs_coeff     = 0.0
lw_ss_albedo     = 0.0

! Compute optical properties of atmosphere 
! Compute scattering coefficient (currently no contribution from carbon_conc, easily added though...)
lw_scatter_coeff  = lw_sca_a + lw_sca_b * q 
! Compute absorption coefficient 
lw_abs_coeff = lw_abs_a + lw_abs_b * q + lw_abs_c * log(carbon_conc / 360.)
! Compute single scattering albedo 
where ((lw_abs_coeff + lw_scatter_coeff).ne.0.0) 
  lw_ss_albedo = lw_scatter_coeff / (lw_abs_coeff + lw_scatter_coeff)
endwhere 

! Compute lw del_tau 
lw_del_tau = 1 / grav * (lw_abs_coeff + lw_scatter_coeff) * & 
             (p_half(:,:,1:n) - p_half(:,:,2:n+1))

! Compute lw optical depth 
!define optical depth to be zero at surface, following Pierrehumbert (2010)
lw_tau(:,:,n+1) = 0.0 
do k = 1, n
  lw_tau(:,:,n+1-k) = lw_tau(:,:,n+2-k) - lw_del_tau(:,:,n+1-k)
enddo 

!! gammas, as defined in Pierrehumbert (2010)
lw_gammaone = gamma *(1.0 - 3./2.*g_asym*lw_ss_albedo) + gammaprime*(1-lw_ss_albedo)
lw_gammatwo = gamma*(1.0 - 3./2.*g_asym*lw_ss_albedo) - gammaprime*(1-lw_ss_albedo)
lw_gammab = lw_gammaone - lw_gammatwo 
do k = 1, n
  lw_gammaplus(:,:,k)  = 0.5*lw_ss_albedo(:,:,k) - gamma*lw_ss_albedo(:,:,k)*3./2.*g_asym*coszen(:,:) 
  lw_gammaminus(:,:,k) = 0.5*lw_ss_albedo(:,:,k) + gamma*lw_ss_albedo(:,:,k)*3./2.*g_asym*coszen(:,:) 
enddo 

! Compute black body radiation for longwave 
pi_B_surf = stefan * t_surf * t_surf * t_surf * t_surf
pi_B = stefan * t * t * t * t

do i = 1, (ie-is+1)
  do j = 1, (je-js+1)
  call two_stream_solver(n, 0.0, 1.0, pi_B_surf(i,j), &
                         pi_B(i,j,:), sw_down_direct(i,j,:)*0.0, &
                         lw_del_tau(i,j,:)/2, lw_gammaone(i,j,:), &
                         lw_gammatwo(i,j,:), lw_gammab(i,j,:), &
                         lw_gammaplus(i,j,:), lw_gammaminus(i,j,:), &
                         lw_up(i,j,:), lw_down(i,j,:))
  ! zero contribution from direct solar beam 
  ! factor 1/2 comes from integrating over all angles
  enddo
enddo 


! =================================================================================
! Compute net fluxes (positive up) and heating rates 

lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down - sw_down_direct
rad_flux = lw_flux + sw_flux

do k = 1, n
   tdt_rad(:,:,k)   = diabatic_acce * ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do

! =================================================================================
surf_lw_down     = lw_down(:, :, n+1)
toa_sw_in        = sw_down_direct(:, :, 1)
toa_sw_out       = sw_up(:,:,1)
net_surf_sw_down = (sw_down_direct(:,:,n+1)+sw_down(:, :, n+1)) * (1. - albedo)
olr         = lw_up(:,:,1)
net_lw_surf = lw_flux(:, :, n+1)
! =================================================================================


!------- downward lw flux surface -------
if ( id_lwdn_sfc > 0 ) then
  used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
endif
!------- incoming sw flux toa -------
if ( id_swdn_toa > 0 ) then
  used = send_data ( id_swdn_toa, toa_sw_in, Time_diag)
endif
!------- outgoing sw flux toa -------
if ( id_swup_toa > 0 ) then
  used = send_data ( id_swup_toa, toa_sw_out, Time_diag)
endif
!------- downward sw flux surface -------
if ( id_swdn_sfc > 0 ) then
  used = send_data ( id_swdn_sfc, net_surf_sw_down, Time_diag)
endif

!------- cosine of zenith angle ------------
if ( id_coszen > 0 ) then
  used = send_data ( id_coszen, coszen, Time_diag)
endif

!------- carbon dioxide concentration ------------
if ( id_co2 > 0 ) then
  used = send_data ( id_co2, carbon_conc, Time_diag)
endif

!------- outgoing lw flux toa (olr) -------
if ( id_olr > 0 ) then
   used = send_data ( id_olr, olr, Time_diag)
endif
!------- upward lw flux surface -------
if ( id_lwup_sfc > 0 ) then
   used = send_data ( id_lwup_sfc, pi_B_surf, Time_diag)
endif
!------- net upward lw flux surface -------
if ( id_net_lw_surf > 0 ) then
   used = send_data ( id_net_lw_surf, net_lw_surf, Time_diag)
endif
!------- temperature tendency due to radiation ------------
if ( id_tdt_rad > 0 ) then
   used = send_data ( id_tdt_rad, tdt_rad, Time_diag)
endif
!------- temperature tendency due to solar radiation ------------
if ( id_tdt_solar > 0 ) then
   used = send_data ( id_tdt_solar, tdt_solar, Time_diag)
endif
!------- total radiative flux (at half levels) -----------
if ( id_flux_rad > 0 ) then
   used = send_data ( id_flux_rad, rad_flux, Time_diag)
endif
!------- longwave radiative flux (at half levels) --------
if ( id_flux_lw > 0 ) then
   used = send_data ( id_flux_lw, lw_flux, Time_diag)
endif
if ( id_flux_sw > 0 ) then
   used = send_data ( id_flux_sw, sw_flux, Time_diag)
endif
! Upward and downward fluxes 
if ( id_flux_lw_up > 0 ) then
  used = send_data ( id_flux_lw_up, lw_up, Time_diag)
endif
if ( id_flux_lw_down > 0 ) then
  used = send_data ( id_flux_lw_down, lw_down, Time_diag)
endif
if ( id_flux_sw_up > 0 ) then
  used = send_data ( id_flux_sw_up, sw_up, Time_diag)
endif
if ( id_flux_sw_down > 0 ) then
  used = send_data ( id_flux_sw_down, sw_down, Time_diag)
endif
if ( id_flux_sw_down_direct > 0 ) then
  used = send_data ( id_flux_sw_down_direct, sw_down_direct, Time_diag)
endif


return
end subroutine two_stream_scatter_rad

! ==================================================================================


subroutine two_stream_scatter_rad_end

  deallocate (sw_tau_0         )
  deallocate (insolation       )
  deallocate (p2               )
  
  deallocate (pi_B_surf        )
  
  deallocate (pi_B             )
  deallocate (tdt_rad          )
  deallocate (tdt_solar        )
  
  deallocate (lw_up            )
  deallocate (lw_down          )
  deallocate (lw_flux          )
  deallocate (sw_up            )
  deallocate (sw_down          )
  deallocate (sw_down_direct   )
  deallocate (sw_flux          )
  deallocate (rad_flux         )
  
  deallocate (lw_tau           )
  deallocate (sw_tau           )
  deallocate (lw_del_tau       )
  deallocate (sw_del_tau       )
  
  deallocate (lw_ss_albedo     )
  deallocate (sw_ss_albedo     )
  
  deallocate (lw_abs_coeff     )
  deallocate (sw_abs_coeff     )
  deallocate (lw_scatter_coeff )
  deallocate (sw_scatter_coeff )
  
  deallocate (sw_gammab        )
  deallocate (sw_gammaminus    )
  deallocate (sw_gammaplus     )
  deallocate (sw_gammaone      )
  deallocate (sw_gammatwo      )
  
  deallocate (lw_gammab        )
  deallocate (lw_gammaminus    )
  deallocate (lw_gammaplus     )
  deallocate (lw_gammaone      )
  deallocate (lw_gammatwo      )
  
  deallocate (olr              )
  deallocate (net_lw_surf      )
  deallocate (toa_sw_in        )
  deallocate (toa_sw_out       )
  deallocate (coszen           )
  deallocate (fracsun          ) 



if(do_read_co2)call interpolator_end(co2_interp)

end subroutine two_stream_scatter_rad_end

! ==================================================================================

subroutine two_stream_solver(nlevels, alpha_g, coszen, pi_B_surf, pi_B,  I_direct, del_tau, gammaone, & 
  gammatwo, gammab, gammaplus, gammaminus, I_up, I_down)

integer, intent(in) :: nlevels 
real, intent(in)              :: alpha_g, coszen,  pi_B_surf  ! should be zero if longwave 
real, intent(in), dimension(nlevels) :: pi_B,  del_tau, gammaone, gammatwo, &
         gammab, gammaplus, gammaminus 
         !pi_B zero if SW, Lzeroexp zero if longwave  
real, intent(in), dimension(nlevels+1) :: I_direct


real, intent(out), dimension(nlevels+1) :: I_up, I_down  

integer :: N  ! N = 2*(nlevels + 1)
integer, parameter :: MU = 3, ML = 3 ! number of upper and lower diagonals 
real, dimension(2*(nlevels+1), 2*(nlevels+1)) :: A, B
real, dimension(2*(nlevels+1)) :: I_comb, source_term
real, dimension(nlevels) :: I_direct_full 
integer :: i, j, lcounter1, lcounter2, IND, M, k, i1, i2
integer, dimension(2*(nlevels+1)) :: IPVT

N = 2*(nlevels + 1)

I_direct_full = (I_direct(1:nlevels) + I_direct(2:nlevels+1)) / 2 ! move I_direct to full model levels, i.e. I_1 = (I_1/2 + I_3/2) / 2 (approximate, not sure how accurate)


! MAKE MATRIX OF ABSORPTION / SCATTERING COEFFICIENTS - band diagonal matrix 
A = 0.0
! fill main diagonal 
lcounter1 = 1
lcounter2 = 1
do i = 1, N
if (i.eq.2) then 
A(i,i) = 1 ! upper bc
elseif (i.eq.N-1) then 
A(i,i) = 1
else
if (mod(i,2).eq.1) then 
A(i,i) = -1 + gammaone(lcounter1) / 2 * del_tau(lcounter1)
lcounter1 = lcounter1 + 1
else 
A(i,i) = 1 - gammaone(lcounter2) / 2 * del_tau(lcounter2)
lcounter2 = lcounter2 + 1
endif 
endif 
enddo 

! fill first upper and lower diagonal 
lcounter1 = 1
lcounter2 = 1
do i = 1, N-1
if (mod(i,2).eq.1) then
if (i.eq.N-1) then 
A(i,i+1) = - alpha_g ! lower bc
else 
A(i,i+1) = -gammatwo(lcounter1) / 2 * del_tau(lcounter1)
lcounter1 = lcounter1 + 1
endif 
!else 
if (i.ne.1) then ! for i = 1, A(2,1) = 0 for upper bc 
A(i+1, i) = gammatwo(lcounter2) / 2 * del_tau(lcounter2)
lcounter2 = lcounter2 + 1
endif 
endif 
enddo 

! fill second upper and lower diagonal 
lcounter1 = 1
lcounter2 = 1
do i = 1, N-2
if (mod(i,2).eq.1) then 
A(i, i+2) = 1+ gammaone(lcounter1) / 2 * del_tau(lcounter1)
lcounter1 = lcounter1 + 1
else 
A(i+2, i) = -1 - gammaone(lcounter2) / 2 * del_tau(lcounter2)
lcounter2 = lcounter2 + 1
endif 
enddo 

!fill third upper and lower diagonal 
lcounter1 = 1
lcounter2 = 1
do i = 1, N-3
if (mod(i,2).eq.1) then 
A(i, i+3) = - gammatwo(lcounter1) / 2 * del_tau(lcounter1)
lcounter1 = lcounter1 + 1
!else 
A(i+3, i) = gammatwo(lcounter2) / 2 * del_tau(lcounter2)
lcounter2 = lcounter2 + 1
endif 
enddo 




! make banded matrix B from A 
M = ML + MU + 1
do j = 1, N
i1 = max(1,j-MU)
i2 = min(N,j+ML)
do i = i1, i2 
k = i - j + M
B(k,j) = A(i,j)
enddo 
enddo 

! Make source term 
source_term = 0.0
source_term(1) = del_tau(1) * (gammab(1) * pi_B(1) + gammaplus(1) * I_direct_full(1) / coszen)
lcounter1 = 2
lcounter2 = 1
do i = 3, N-2
if (mod(i,2).eq.1) then 
source_term(i) = del_tau(lcounter1) * (gammab(lcounter1) * pi_B(lcounter1) + &
gammaplus(lcounter1) * I_direct_full(lcounter1) / coszen)
lcounter1 = lcounter1 + 1
else
source_term(i) = del_tau(lcounter2) * (-1*gammab(lcounter2) * pi_B(lcounter2) - &
gammaminus(lcounter2) * I_direct_full(lcounter2) / coszen) 
lcounter2 = lcounter2 + 1
endif 
enddo 
source_term(N-1) = alpha_g * I_direct(nlevels+1) + (1 - alpha_g) * pi_B_surf
source_term(N) = del_tau(nlevels) * (-1*gammab(nlevels) * pi_B(nlevels) - gammaminus(nlevels) * I_direct_full(nlevels) / coszen) 


! call LU factorization  
CALL NSBFAC(B,N,N,ML,MU,IPVT,IND)
! solve AX = B for X 
if (IND .ne. 0) then 
write(*,*) 'MATRIX IS ILL-CONDITIONED'
endif 
CALL NSBSLV(B,N,N,ML,MU,IPVT,source_term,I_comb)

lcounter1 = 1
lcounter2 = 1
do i = 1, N
if (mod(i,2).eq.1) then 
I_up(lcounter1) = I_comb(i)
lcounter1 = lcounter1 + 1
else 
I_down(lcounter2) = I_comb(i)
lcounter2 = lcounter2 + 1
endif 
enddo 


end subroutine two_stream_solver

SUBROUTINE NSBFAC(B,LDB,N,ML,MU,IPVT,IND)
!-------------------------------------------------------------------
!     LU factorization of a band matrix (non symmetric) with partial
!     pivoting.
!     INPUTS:
!     B  : banded matrix. The correspondance between full matrix
!     A(i,j) and band matrix B(k,l) is given by following sequence:
!     m=ml+mu+1
!     do j=1,n
!       i1=max(1,j-mu)
!       i2=min(n,j+ml)
!       do i=i1,i2
!         k=i-j+m
!         b(k,j)=a(i,j)
!       enddo
!     enddo
!     LDB : 1st dimension of B in calling program (ldb.ge.2ml+mu+1)
!     N   : size of B                             -----------------
!     ML  : number of lower diagonals
!     MU  : number of upper diagonals
!     OUTPUTS:
!     B   : banded matrix storing the LU elements, the ML first lines
!           of which are used for pivoting.
!     IPVT: integer vector of size N storing the pivoting indices.
!     IND : flag = 0,  B is non singular
!                = k,  B may be singular
!-------------------------------------------------------------------
INTEGER :: N, LDB, ML, MU, IND, OUT
INTEGER :: I, M, J0, J1, JZ, I0, J, JU, K, KP1, L, LM, MM, n_m_one 
REAL :: B(LDB,N),T
INTEGER IPVT(*)
M=ML+MU+1
IND=0

J0=MU+2
J1=MIN(N,M)-1
IF(J1.GE.J0) THEN
DO JZ=J0,J1
I0=M+1-JZ
DO I=I0,ML
B(I,JZ)=0.D0
ENDDO
ENDDO
ENDIF
JZ=J1
JU=0

n_m_one=N-1
IF(n_m_one.GE.1) THEN
DO K=1,n_m_one
KP1=K+1
JZ=JZ+1
IF(JZ.LE.N.AND.ML.GE.1) THEN
DO I=1,ML
B(I,JZ)=0.D0
ENDDO
ENDIF

LM=MIN(ML,N-K)
call IAMAX(B(M,K),LM+1, OUT)
L=OUT+M-1
IPVT(K)=L+K-M
IF(B(L,K).EQ.0.D0) GO TO 10
IF(L.NE.M) THEN
T=B(L,K)
B(L,K)=B(M,K)
B(M,K)=T
ENDIF
T=-1.D0/B(M,K)
CALL SCALE(LM,T,B(M+1,K))

JU=MIN(MAX(JU,MU+IPVT(K)),N)
MM=M
IF(JU.GE.KP1) THEN
DO J=KP1,JU
L=L-1
MM=MM-1
T=B(L,J)
IF(L.NE.MM) THEN
B(L,J)=B(MM,J)
B(MM,J)=T
ENDIF
CALL DAXPY(LM,T,B(M+1,K),B(MM+1,J))
ENDDO
ENDIF
GO TO 20
10     ind=k
20     CONTINUE
ENDDO
ENDIF
IPVT(N)=N
IF(B(M,N).EQ.0.D0) IND=N
END SUBROUTINE NSBFAC

SUBROUTINE DAXPY(N,A,X,Y)
INTEGER N, I
REAL X(*),Y(*),A
DO I=1,N
Y(I)=Y(I)+A*X(I)
ENDDO
END SUBROUTINE DAXPY 

SUBROUTINE IAMAX(A,N,FINAL)
INTEGER N, I, FINAL 
REAL A(*),T
T=0.d0
DO I=1,N
IF(ABS(A(I)).GT.T) THEN
t=abs(A(I))
FINAL=I
ENDIF
enddo
END SUBROUTINE IAMAX 

SUBROUTINE SCALE(N,T,A)
INTEGER N, I
REAL A(*),T
DO I=1,N
A(I)=T*A(I)
ENDDO
END SUBROUTINE SCALE 

SUBROUTINE NSBSLV(A,LDA,N,ML,MU,IPVT,b,x)
!--------------------------------------------------------------------
!     Solve banded linear system Ax = b
!     INPUTS:
!     A   : banded matrix as output of LU factorization by NSBFAC
!           (see storing mode in NSBFAC subroutine).
!     LDA : 1st dimension of A in calling program (lda.ge.2ml+mu+1)
!     N   : order of A                             -----------------
!     ML  : number of lower diagonals
!     MU  : number of upper diagonals
!     IPVT: integer vector of size N storing the pivoting indices
!           as output of NSBFAC.
!     b   : second member vector
!     OUTPUT:
!     x   : solution vector
!---------------------------------------------------------------------

INTEGER :: N, ML, MU, LDA
REAL :: A(LDA,*),B(*),X(*),T
INTEGER IPVT(*)
INTEGER :: I, K, L, M, LM, LA, LB, n_m_one 
DO I=1,N
X(I)=B(I)
ENDDO
M=ML+MU+1
n_m_one=N-1
!     solve L*y = b
IF(ML.NE.0.AND.n_m_one.GE.1) THEN
DO K=1,n_m_one
LM=MIN(ML,N-K)
L=IPVT(K)
T=X(L)
IF(L.NE.K) THEN
X(L)=X(K)
X(K)=T
ENDIF
CALL DAXPY(LM,T,A(M+1,K),X(K+1))
ENDDO
ENDIF
!     solve U*y = x
DO K=N,1,-1
X(K)=X(K)/A(M,K)
LM=MIN(K,M)-1
LA=M-LM
LB=K-LM
T=-X(K)
CALL DAXPY(LM,T,A(LA,K),X(LB))
ENDDO
END SUBROUTINE NSBSLV

end module two_stream_scatter_rad_mod
