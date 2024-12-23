!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module two_stream_gray_rad_mod

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

character(len=128) :: version='$Id: two_stream_gray_rad.F90,v 1.1.2.1 2016/02/23$'
character(len=128) :: tag='Two-stream gray atmosphere'

! Version Details
! [2016/01/28] <James Penn>:  GFDL Astronomy module used to calculate insolation
!   - Diurnal & Seasonal cycle implemented in astronomy_mod.
!   - set `do_seasonal` to .true. in namelist to activate.
!   - default (`do_seasonal = .false.`) remains as permanent global equinox.
!   - tidal locking can be achieved with by setting `2*pi / orbital_period = omega`.
! [2016/02/23] <James Penn>:  Window parameterisation in the longwave stream
!   - reference: Ruth Geen etal, GRL 2016 (supp. information).
!   - Set `wv_exponent` and `solar_exponent` to -1 in namelist to activate.
!   - Moisture and CO2 feedback on SW and LW controlled by `ir_tau...` constants.
!==================================================================================

! public interfaces

public :: two_stream_gray_rad_init, two_stream_gray_rad_down, two_stream_gray_rad_up, two_stream_gray_rad_end
!==================================================================================


! module variables
logical :: initialized     = .false.
real    :: solar_constant  = 1360.0
real    :: del_sol         = 1.4
real    :: del_sw          = 0.0
real    :: ir_tau_eq       = 6.0
real    :: ir_tau_pole     = 1.5
real    :: atm_abs         = 0.0
real    :: odp = 1.0  ! odp = optical depth parameter. Used to simulate GHG concentration, taken from Tapio FMS, added by Nathanael Wong
real    :: sw_diff         = 0.0
real    :: linear_tau      = 0.1
real    :: surf_pres       = 100000.0
real    :: wv_exponent     = 4.0
real    :: solar_exponent  = 4.0
logical :: do_seasonal     = .false.
logical :: tidally_locked  = .true.
real    :: noon_longitude  = 270.0
integer :: solday          = -10  !s Day of year to run perpetually if do_seasonal=True and solday>0
real    :: equinox_day     = 0.75 !s Fraction of year [0,1] where NH autumn equinox occurs (only really useful if calendar has defined months).
logical :: use_time_average_coszen = .false. !s if .true., then time-averaging is done on coszen so that insolation doesn't depend on timestep
real    :: dt_rad_avg     = -1

character(len=32) :: rad_scheme = 'frierson'

integer, parameter :: B_GEEN = 1,  B_FRIERSON = 2, &
                      B_BYRNE = 3, B_SCHNEIDER_LIU=4
integer, private :: sw_scheme = B_FRIERSON
integer, private :: lw_scheme = B_FRIERSON

! constants for RG radiation version
real    :: ir_tau_co2_win  = 0.2150
real    :: ir_tau_wv_win1  = 147.11
real    :: ir_tau_wv_win2  = 1.0814e4
real    :: ir_tau_co2      = 0.1
real    :: ir_tau_wv1      = 23.8
real    :: ir_tau_wv2      = 254.0
real    :: window          = 0.3732
real    :: carbon_conc     = 360.0

! constants for SCHNEIDER_LIU radiation version
real    :: single_albedo      = 0.8
real    :: back_scatter       = 0.398
real    :: lw_tau_0_gp        = 80.0
real    :: sw_tau_0_gp        = 3.0
real    :: lw_tau_exponent_gp = 2.0
real    :: sw_tau_exponent_gp = 1.0
real    :: diabatic_acce = 1.0
real,save :: gp_albedo, Ga_asym, g_asym

! parameters for Byrne and OGorman radiation scheme
real :: bog_a = 0.8678
real :: bog_b = 1997.9
real :: bog_mu = 1.0

!s albedo now defined in mixed_layer_init
real, allocatable, dimension(:,:)   :: insolation, p2, lw_tau_0, sw_tau_0 
real, allocatable, dimension(:,:)   :: b_surf, b_surf_gp
real, allocatable, dimension(:,:,:) :: b, tdt_rad, tdt_solar
real, allocatable, dimension(:,:,:) :: lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux
real, allocatable, dimension(:,:,:) :: lw_tau, sw_tau, lw_dtrans
real, allocatable, dimension(:,:)   :: olr, net_lw_surf, toa_sw_in, coszen, fracsun

! window parameterisation (RG, 2015)
real, allocatable, dimension(:,:,:) :: lw_up_win, lw_down_win, lw_dtrans_win
real, allocatable, dimension(:,:,:) :: b_win, sw_dtrans
real, allocatable, dimension(:,:)   :: sw_wv, del_sol_tau, sw_tau_k, lw_del_tau, lw_del_tau_win

real, save :: pi, deg_to_rad , rad_to_deg

!extras for reading in co2 concentration
logical                             :: do_read_co2=.false.
type(interpolate_type),save         :: co2_interp           ! use external file for co2
character(len=256)                  :: co2_file='co2'       !  file name of co2 file to read
character(len=256)                  :: co2_variable_name='co2'       !  file name of co2 file to read


namelist/two_stream_gray_rad_nml/ solar_constant, del_sol, &
           ir_tau_eq, ir_tau_pole, odp, atm_abs, sw_diff, &
           linear_tau, del_sw, wv_exponent, &
           solar_exponent, do_seasonal, &
           ir_tau_co2_win, ir_tau_wv_win1, ir_tau_wv_win2, &
           ir_tau_co2, ir_tau_wv1, ir_tau_wv2, &
           surf_pres,tidally_locked,noon_longitude, & !Ruizhi add for grey
		       window, carbon_conc, rad_scheme, &
           do_read_co2, co2_file, co2_variable_name, solday, equinox_day, bog_a, bog_b, bog_mu, &
           use_time_average_coszen, dt_rad_avg,&
           diabatic_acce !Schneider Liu values

!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_net_lw_surf, id_lwdn_sfc, id_lwup_sfc, &
           id_tdt_rad, id_tdt_solar, id_flux_rad, id_flux_lw, id_flux_sw, id_coszen, id_fracsun, &
           id_lw_dtrans, id_lw_dtrans_win, id_sw_dtrans, id_co2

character(len=10), parameter :: mod_name = 'two_stream'

real :: missing_value = -999.


contains


! ==================================================================================
! ==================================================================================


subroutine two_stream_gray_rad_init(is, ie, js, je, num_levels, axes, Time, lonb, latb, dt_real)

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
   read  (unit, nml=two_stream_gray_rad_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'two_stream_gray_rad_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=two_stream_gray_rad_nml)
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


if(uppercase(trim(rad_scheme)) == 'GEEN') then
  lw_scheme = B_GEEN
  call error_mesg('two_stream_gray_rad','Using Geen (2015) radiation scheme.', NOTE)
else if(uppercase(trim(rad_scheme)) == 'FRIERSON') then
  lw_scheme = B_FRIERSON
  !call error_mesg('two_stream_gray_rad','Using "'//trim(rad_scheme)//'" radiation scheme.', NOTE)
  call error_mesg('two_stream_gray_rad','Using Frierson (2006) radiation scheme.', NOTE)
else if(uppercase(trim(rad_scheme)) == 'BYRNE') then
  lw_scheme = B_BYRNE
  call error_mesg('two_stream_gray_rad','Using Byrne & OGorman (2013) radiation scheme.', NOTE)
else if(uppercase(trim(rad_scheme)) == 'SCHNEIDER') then
  lw_scheme = B_SCHNEIDER_LIU
  call error_mesg('two_stream_gray_rad','Using Schneider & Liu (2009) radiation scheme for GIANT PLANETS.', NOTE)
else
  call error_mesg('two_stream_gray_rad','"'//trim(rad_scheme)//'"'//' is not a valid radiation scheme.', FATAL)
endif
sw_scheme = lw_scheme


if(lw_scheme == B_SCHNEIDER_LIU) then
	g_asym          = 1 - 2.*back_scatter;                 ! asymmetry factor
	gp_albedo     = ( sqrt(1. - g_asym*single_albedo) - sqrt(1. - single_albedo) )        &
           / ( sqrt(1. - g_asym*single_albedo) + sqrt(1. - single_albedo) );
	Ga_asym         = 2.*sqrt( (1. - single_albedo) * (1. - g_asym*single_albedo) );
endif

if ((lw_scheme == B_BYRNE).or.(lw_scheme == B_GEEN)) then
    if (pstd_mks/=pstd_mks_earth) call error_mesg('two_stream_gray_rad','Pstd_mks and pstd_mks_earth are not the same in the this run, but lw scheme will use pstd_mks_earth because abs coeffs in Byrne and Geen schemes are non-dimensionalized by Earth surface pressure.', NOTE)
endif



initialized = .true.

allocate (b                (ie-is+1, je-js+1, num_levels))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_solar        (ie-is+1, je-js+1, num_levels))

allocate (lw_dtrans        (ie-is+1, je-js+1, num_levels))
allocate (lw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (sw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (lw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (lw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (lw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (rad_flux         (ie-is+1, je-js+1, num_levels+1))

allocate (b_surf           (ie-is+1, je-js+1))
allocate (lw_tau_0         (ie-is+1, je-js+1))
allocate (sw_tau_0         (ie-is+1, je-js+1))
allocate (olr              (ie-is+1, je-js+1))
allocate (net_lw_surf      (ie-is+1, je-js+1))
allocate (toa_sw_in        (ie-is+1, je-js+1))

allocate (insolation       (ie-is+1, je-js+1))
allocate (p2               (ie-is+1, je-js+1))
allocate (coszen           (ie-is+1, je-js+1))
allocate (fracsun          (ie-is+1, je-js+1)) !jp from astronomy.f90 : fraction of sun on surface

!Allocate RG variables only if option selected.
select case(lw_scheme)
case(B_GEEN)
  allocate (b_win            (ie-is+1, je-js+1, num_levels))
  allocate (lw_dtrans_win    (ie-is+1, je-js+1, num_levels))
  allocate (lw_up_win        (ie-is+1, je-js+1, num_levels+1))
  allocate (lw_down_win      (ie-is+1, je-js+1, num_levels+1))
  allocate (lw_del_tau       (ie-is+1, je-js+1))
  allocate (lw_del_tau_win   (ie-is+1, je-js+1))

  allocate (sw_dtrans        (ie-is+1, je-js+1, num_levels))
  allocate (sw_wv            (ie-is+1, je-js+1))
  allocate (del_sol_tau      (ie-is+1, je-js+1))
  allocate (sw_tau_k         (ie-is+1, je-js+1))

case(B_BYRNE)
  allocate (lw_del_tau       (ie-is+1, je-js+1))

case(B_SCHNEIDER_LIU)
	allocate (b_surf_gp      (ie-is+1, je-js+1))
end select


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

  if (lw_scheme.eq.B_GEEN) then
    id_lw_dtrans_win  = &
               register_diag_field ( mod_name, 'lw_dtrans_win', axes(1:3), Time, &
                 'LW window transmission', &
                 'none', missing_value=missing_value      )
  endif
  if (sw_scheme.eq.B_GEEN) then
    id_sw_dtrans  = &
                 register_diag_field ( mod_name, 'sw_dtrans', axes(1:3), Time, &
                   'SW transmission', &
                   'none', missing_value=missing_value      )
  endif
    id_lw_dtrans  = &
               register_diag_field ( mod_name, 'lw_dtrans', axes(1:3), Time, &
                 'LW transmission (non window)', &
                 'none', missing_value=missing_value      )
return
end subroutine two_stream_gray_rad_init

! ==================================================================================

subroutine two_stream_gray_rad_down (is, js, Time_diag, lat, lon, p_half, t,         &
                           net_surf_sw_down, surf_lw_down, albedo, q)

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in), dimension(:,:)    :: lat, lon, albedo
real, intent(out), dimension(:,:)   :: net_surf_sw_down
real, intent(out), dimension(:,:)   :: surf_lw_down
real, intent(in), dimension(:,:,:)  :: t, q,  p_half
integer :: i, j, k, n, dyofyr

integer :: seconds, year_in_s, days
real :: r_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, day_in_s, r_solday, r_total_seconds, r_days, r_dt_rad_avg, dt_rad_radians
logical :: used


real ,dimension(size(q,1),size(q,2),size(q,3)) :: co2f

n = size(t,3)



! albedo(:,:) = albedo_value !s albedo now set in mixed_layer_init.

! =================================================================================
! SHORTWAVE RADIATION

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

else if (sw_scheme==B_SCHNEIDER_LIU) then
  insolation = (solar_constant/pi)*cos(lat)
else if (tidally_locked) then
  ! Tidally locked: insolation depends on latitude and longitude
  do i = 1, size(t, 1)
     do j = 1, size(t, 2)
        insolation(i, j) = solar_constant*cos(lat(i,j))*cos(lon(i,j) - noon_longitude*deg_to_rad)
        if (insolation(i, j) .LT. 0.0) then
           insolation(i, j) = 0.0
        endif
     enddo
  enddo
else
  ! Default: Averaged Earth insolation at all longitudes
  p2          = (1. - 3.*sin(lat)**2)/4.
  insolation  = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * sin(lat))
end if

select case(sw_scheme)
case(B_GEEN)
  ! RG scheme: optical depth a function of wv and co2
  sw_tau_k = 0.
  do k = 1, n
    sw_wv = sw_tau_k + 0.5194
    sw_wv = exp( 0.01887 / (sw_tau_k + 0.009522)                      &
                 + 1.603 / ( sw_wv*sw_wv ) )
    del_sol_tau(:,:) = ( 0.0596 + 0.0029 * log(carbon_conc/360.)       &
                                + sw_wv(:,:) * q(:,:,k) )             &
                     * ( p_half(:,:,k+1) - p_half(:,:,k) ) / p_half(:,:,n+1)
    sw_dtrans(:,:,k) = exp( - del_sol_tau(:,:) )
    sw_tau_k = sw_tau_k + del_sol_tau(:,:)
  end do

  ! compute downward shortwave flux
  sw_down(:,:,1) = insolation(:,:)
  do k = 1, n
     sw_down(:,:,k+1)   = sw_down(:,:,k) * sw_dtrans(:,:,k)
  end do

case(B_FRIERSON, B_BYRNE)
  ! Default: Frierson handling of SW radiation
  ! SW optical thickness
  sw_tau_0    = (1.0 - sw_diff*sin(lat)**2)*atm_abs

  ! compute optical depths for each model level
  do k = 1, n+1
    sw_tau(:,:,k) = sw_tau_0 * (p_half(:,:,k)/pstd_mks)**solar_exponent
  end do

  ! compute downward shortwave flux
  do k = 1, n+1
     sw_down(:,:,k)   = insolation(:,:) * exp(-sw_tau(:,:,k))
  end do

case(B_SCHNEIDER_LIU)
  ! Schneider & Liu 2009 Giant planet scheme
  ! SW optical thickness

  ! compute optical depths for each model level
  do k = 1, n+1
    sw_tau(:,:,k) = sw_tau_0_gp * (p_half(:,:,k)/pstd_mks)**sw_tau_exponent_gp
  end do

  ! compute downward shortwave flux
  do k = 1, n+1
     sw_down(:,:,k)   = insolation(:,:) * (1.0-gp_albedo)* exp(- Ga_asym * sw_tau(:,:,k))
  end do

case default
 call error_mesg('two_stream_gray_rad','invalid radiation scheme',FATAL)

end select

! =================================================================================
! LONGWAVE RADIATION

! longwave source function
b = stefan*t**4

if(do_read_co2)then
  call interpolator( co2_interp, Time_diag, p_half, co2f, trim(co2_variable_name))
  carbon_conc = maxval(co2f) !Needs maxval just because co2f is a 3d array of a constant, so maxval is just a way to pick out one number
endif


select case(lw_scheme)
case(B_GEEN)
  ! split LW in 2 bands: water-vapour window and remaining = non-window
  ! ref: Ruth Geen etal, GRL 2016 (supp. information).
  lw_dtrans_win = 1.
  do k = 1, n
    lw_del_tau    = ( ir_tau_co2 + 0.2023 * log(carbon_conc/360.)                  &
                    + ir_tau_wv1*log(ir_tau_wv2*q(:,:,k) + 1) )                    &
               * ( p_half(:,:,k+1)-p_half(:,:,k) ) / pstd_mks_earth
    lw_dtrans(:,:,k) = exp( - lw_del_tau )
    lw_del_tau_win   = ( ir_tau_co2_win + 0.0954 * log(carbon_conc/360.)           &
                                     + ir_tau_wv_win1*q(:,:,k)                 &
                                     + ir_tau_wv_win2*q(:,:,k)*q(:,:,k) )      &
                  * ( p_half(:,:,k+1)-p_half(:,:,k) ) / pstd_mks_earth
    lw_dtrans_win(:,:,k) = exp( - lw_del_tau_win )
  end do

  ! also compute downward longwave flux for window
  ! Allocate a fraction of the longwave spectrum as window radiation
  b_win = window*b
  b = (1.0 - window)*b
  lw_down_win(:,:,1) = 0.0
  lw_down(:,:,1) = 0.0
  do k = 1,n
    lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
    lw_down_win(:,:,k+1) = lw_down_win(:,:,k)*lw_dtrans_win(:,:,k)              &
                      + b_win(:,:,k)*(1.0 - lw_dtrans_win(:,:,k))
  end do
  lw_down = lw_down + lw_down_win

case(B_BYRNE)
  ! dtau/ds = a*mu + b*q
  ! ref: Byrne, M. P. & O'Gorman, P. A.
  !      Land–ocean warming contrast over a wide range of climates:
  !      Convective quasi-equilibrium theory and idealized simulations.
  !      J. Climate 26, 4000–4106 (2013).

  do k = 1, n
    lw_del_tau    = (bog_a*bog_mu + 0.17 * log(carbon_conc/360.)  + bog_b*q(:,:,k)) * (( p_half(:,:,k+1)-p_half(:,:,k) ) / pstd_mks_earth)
    lw_dtrans(:,:,k) = exp( - lw_del_tau )

  end do

  ! compute downward longwave flux by integrating downward
  lw_down(:,:,1)      = 0.
  do k = 1, n
     lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
  end do

case(B_FRIERSON)
  ! longwave optical thickness function of latitude and pressure
  lw_tau_0 = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*sin(lat)**2
  lw_tau_0 = lw_tau_0 * odp ! scale by optical depth parameter - default 1

  ! compute optical depths for each model level, RUIZHI
  do k = 1, n+1
  !lw_tau(:,:,k) = lw_tau_0 * ( linear_tau * p_half(:,:,k)/pstd_mks     &
  !     + (1.0 - linear_tau) * (p_half(:,:,k)/pstd_mks)**wv_exponent )
  lw_tau(:,:,k) = 0.1 * p_half(:,:,k)/surf_pres &
  + linear_tau*surf_pres/pstd_mks*(p_half(:,:,k)/surf_pres)**wv_exponent
  end do

  ! longwave differential transmissivity
  do k = 1, n
     lw_dtrans(:,:,k) = exp( -(lw_tau(:,:,k+1) - lw_tau(:,:,k)) )
  end do

  ! compute downward longwave flux by integrating downward
  lw_down(:,:,1)      = 0.
  do k = 1, n
     lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
  end do

case(B_SCHNEIDER_LIU)

  ! compute optical depths for each model level
  do k = 1, n+1
  lw_tau(:,:,k) = lw_tau_0_gp * (p_half(:,:,k)/pstd_mks)**lw_tau_exponent_gp
  end do

  ! longwave differential transmissivity
  do k = 1, n
     lw_dtrans(:,:,k) = exp( -(lw_tau(:,:,k+1) - lw_tau(:,:,k)) )
  end do

  ! compute downward longwave flux by integrating downward
  lw_down(:,:,1)      = 0.
  do k = 1, n
     lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
  end do


case default
 call error_mesg('two_stream_gray_rad','invalid radiation scheme',FATAL)

end select

! =================================================================================
surf_lw_down     = lw_down(:, :, n+1)
toa_sw_in        = sw_down(:, :, 1)
net_surf_sw_down = sw_down(:, :, n+1) * (1. - albedo)
! =================================================================================

if(lw_scheme.eq.B_SCHNEIDER_LIU) then
	b_surf_gp=surf_lw_down(:,:)+net_surf_sw_down(:,:)
endif


!------- downward lw flux surface -------
if ( id_lwdn_sfc > 0 ) then
   used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
endif
!------- incoming sw flux toa -------
if ( id_swdn_toa > 0 ) then
   used = send_data ( id_swdn_toa, toa_sw_in, Time_diag)
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

return
end subroutine two_stream_gray_rad_down

! ==================================================================================

subroutine two_stream_gray_rad_up (is, js, Time_diag, lat, p_half, t_surf, t, tdt, albedo)

! Now complete the radiation calculation by computing the upward and net fluxes.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat, albedo
real, intent(in) , dimension(:,:)   :: t_surf
real, intent(in) , dimension(:,:,:) :: t, p_half
real, intent(inout), dimension(:,:,:) :: tdt


integer :: i, j, k, n

logical :: used

n = size(t,3)
b_surf            = stefan*t_surf**4

select case(lw_scheme)
case(B_GEEN)
  ! integrate upward, including window contribution
  lw_up(:,:,n+1)     = b_surf*(1-window)
  lw_up_win(:,:,n+1) = b_surf*window

  do k = n, 1, -1
     lw_up(:,:,k)   = lw_up(:,:,k+1)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
     lw_up_win(:,:,k)   = lw_up_win(:,:,k+1)*lw_dtrans_win(:,:,k) + b_win(:,:,k)*(1.0 - lw_dtrans_win(:,:,k))
  end do
  lw_up = lw_up + lw_up_win

case(B_FRIERSON, B_BYRNE)
  ! compute upward longwave flux by integrating upward
  lw_up(:,:,n+1)    = b_surf
  do k = n, 1, -1
     lw_up(:,:,k)   = lw_up(:,:,k+1)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
  end do

case(B_SCHNEIDER_LIU)
  ! compute upward longwave flux by integrating upward

  lw_up(:,:,n+1)    = b_surf_gp
  do k = n, 1, -1
     lw_up(:,:,k)   = lw_up(:,:,k+1)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
  end do

case default
 call error_mesg('two_stream_gray_rad','invalid radiation scheme',FATAL)

end select

! compute upward shortwave flux (here taken to be constant)
do k = 1, n+1
   sw_up(:,:,k)   = albedo(:,:) * sw_down(:,:,n+1)
end do

! net fluxes (positive up)
lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down
rad_flux = lw_flux + sw_flux

do k = 1, n
   tdt_rad(:,:,k)   = diabatic_acce * ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do

olr         = lw_up(:,:,1)
net_lw_surf = lw_flux(:, :, n+1)

!------- outgoing lw flux toa (olr) -------
if ( id_olr > 0 ) then
   used = send_data ( id_olr, olr, Time_diag)
endif
!------- upward lw flux surface -------
if ( id_lwup_sfc > 0 ) then
   used = send_data ( id_lwup_sfc, b_surf, Time_diag)
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
!------- radiative transmissivities (at half levels) --------
if ( id_lw_dtrans > 0 ) then
   used = send_data ( id_lw_dtrans, lw_dtrans, Time_diag)
endif
if ( id_lw_dtrans_win > 0 ) then
   used = send_data ( id_lw_dtrans_win, lw_dtrans_win, Time_diag)
endif
if ( id_sw_dtrans > 0 ) then
   used = send_data ( id_sw_dtrans, sw_dtrans, Time_diag)
endif

return
end subroutine two_stream_gray_rad_up

! ==================================================================================


subroutine two_stream_gray_rad_end

deallocate (b, tdt_rad, tdt_solar)
deallocate (lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux)
deallocate (b_surf, olr, net_lw_surf, toa_sw_in, lw_tau_0, sw_tau_0)
deallocate (lw_dtrans, lw_tau, sw_tau)
deallocate (insolation, p2) !s albedo

! deallocate RG variables
if (lw_scheme.eq.B_GEEN) then
  deallocate (b_win, lw_dtrans_win, lw_up_win, lw_down_win, lw_del_tau, lw_del_tau_win)
endif
if (sw_scheme.eq.B_GEEN) then
  deallocate (sw_dtrans, sw_wv, del_sol_tau, sw_tau_k)
endif
if (lw_scheme.eq.B_BYRNE) then
  deallocate (lw_del_tau)
endif
if(lw_scheme.eq.B_SCHNEIDER_LIU) then
	deallocate (b_surf_gp)
endif

if(do_read_co2)call interpolator_end(co2_interp)

end subroutine two_stream_gray_rad_end

! ==================================================================================

end module two_stream_gray_rad_mod
