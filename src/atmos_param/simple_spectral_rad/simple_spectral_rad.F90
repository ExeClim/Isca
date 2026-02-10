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

module simple_spectral_rad_mod

! ==================================================================================
! ==================================================================================

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file, error_mesg, &
                                    NOTE, FATAL,  uppercase

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks, pstd_mks_earth, seconds_per_sol, orbital_period

   use    diag_manager_mod,   only: register_diag_field, send_data, diag_axis_init

   use    time_manager_mod,   only: time_type, length_of_year, length_of_day, &
                                    operator(+), operator(-), operator(/=), get_time

   use astronomy_mod,         only: astronomy_init, diurnal_solar

   use interpolator_mod,      only: interpolate_type, interpolator_init, interpolator, ZERO, interpolator_end

!==================================================================================
implicit none
private
!==================================================================================

! version information

character(len=128) :: version='$Id: simple_spectral_rad.F90,v 1.1.1.1 2025/02/28$'
character(len=128) :: tag='Two-stream simple spectral atmosphere'

! Version Details
! [2025/02/28] <Andrew Williams>: initial commit
! I am keeping the simple frierson shortwave scheme, but replacing the
! longwave scheme with a simple spectral scheme where the absorption 
! coefficients are approximated analytically.
!==================================================================================

! public interfaces

public :: simple_spectral_rad_init, simple_spectral_rad_down, simple_spectral_rad_up, simple_spectral_rad_end
!==================================================================================

! Make specific variables public
public :: debug_ssm
  
! module variables
logical :: initialized     = .false.
real    :: solar_constant  = 1360.0
real    :: del_sol         = 1.4
real    :: del_sw          = 0.0
real    :: atm_abs         = 0.0
real    :: sw_diff         = 0.0
real    :: solar_exponent  = 4.0
logical :: do_seasonal     = .false.
integer :: solday          = -10  !s Day of year to run perpetually if do_seasonal=True and solday>0
real    :: equinox_day     = 0.75 !s Fraction of year [0,1] where NH autumn equinox occurs (only really useful if calendar has defined months).
logical :: use_time_average_coszen = .false. !s if .true., then time-averaging is done on coszen so that insolation doesn't depend on timestep
real    :: dt_rad_avg     = -1

! CO2 params
logical :: do_co2          = .false.
real    :: carbon_conc     = 280.0 ! ppmv

! Some physical constants
real    :: planck_h     = 6.626075540e-34    ! Planck's constant
real    :: lightspeed   = 2.99792458e8       ! Speed of light
real    :: boltzmann_k  = 1.38065812e-23     ! Boltzman thermodynamic constant
real    :: stefan_sigma = 5.67051196e-8      ! Stefan-Boltzman constant

! Spectral parameters (todo: expose as namelist options)
real    :: k_rot  = 37 
real    :: l_rot  = 56
real    :: nu_rot = 200 

real    :: k_vr   = 5
real    :: l_vr   = 37
real    :: nu_vr  = 1450

real    :: k_ad   = 5 ! _ad is "additional", that is, in 1700-2500cm-1 range 
real    :: l_ad   = 52
real    :: nu_ad  = 1700

! include H2O continuum 
real    :: k_cnt_1    = 0.0040 ! below 1700cm-1 m^2/kg
real    :: k_cnt_2    = 0.0002 ! above 1700cm-1 m^2/kg
real    :: a          = 0.02  ! K-1, Mlawer et al. 2012, Jeevanjee et al 2018
real    :: Tref       = 260.  ! K
real    :: pvstar_ref = 224.92 ! Pa, = pvstar(Tref) in fitting setup

! CO2 fits 
real    :: k_co2  = 110
real    :: l_co2  = 11.5
real    :: nu_co2 = 667.5

! two-stream diffusivity factor
real    :: D  = 1.5 

! Defining the longwave spectral grid
integer :: nwavenumber     = 41
real    :: wavenumber_min  = 10   ! cm-1
real    :: wavenumber_max  = 2510 ! cm-1
real    :: dwavenumber    ! = (wavenumber_max-wavenumber_min) / nwavenumber    ! cm-1

! Mechanism denial flags
logical :: do_pressure_broadening = .true.

! DEBUGGING FLAGS
logical :: debug_ssm = .false.

real, allocatable, dimension(:,:)     :: insolation, p2, sw_tau_0 !s albedo now defined in mixed_layer_init
real, allocatable, dimension(:,:,:,:) :: lw_tau_nu, lw_tau_h2o_nu, lw_tau_h2o_cont_nu, lw_tau_co2_nu
real, allocatable, dimension(:,:,:,:) :: lw_dtrans_nu
real, allocatable, dimension(:,:,:,:) :: b_nu, lw_up_nu, lw_down_nu
real, allocatable, dimension(:,:,:)   :: b_surf_nu, olr_spec 
real, allocatable, dimension(:,:,:)   :: tdt_rad, tdt_solar
real, allocatable, dimension(:,:,:)   :: lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux
real, allocatable, dimension(:,:,:)   :: sw_tau 
real, allocatable, dimension(:,:)     :: b_surf, olr, net_lw_surf, toa_sw_in, coszen, fracsun

real, allocatable, dimension(:)       :: nu_array, kappa

real, save :: pi, deg_to_rad , rad_to_deg

!extras for reading in co2 concentration
logical                             :: do_read_co2=.false.
type(interpolate_type),save         :: co2_interp           ! use external file for co2
character(len=256)                  :: co2_file='co2'       !  file name of co2 file to read
character(len=256)                  :: co2_variable_name='co2'       !  file name of co2 file to read


namelist/simple_spectral_rad_nml/ solar_constant, del_sol, &
           solar_exponent, do_seasonal, solday, equinox_day, &
           atm_abs, sw_diff, del_sw, &
		   do_co2, carbon_conc, do_read_co2, co2_file, co2_variable_name, & 
           use_time_average_coszen, dt_rad_avg, &
           k_rot, k_vr, k_cnt_1, k_cnt_2, nwavenumber, do_pressure_broadening, debug_ssm

!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_spectral_olr, id_ssm_bins_lw, id_swdn_sfc, id_swdn_toa, id_net_lw_surf, id_lwdn_sfc, id_lwup_sfc, &
           id_tdt_rad, id_tdt_solar, id_flux_rad, id_flux_lw, id_flux_sw, id_coszen, id_fracsun, &
           id_lw_dtrans, id_co2

character(len=15), parameter :: mod_name = 'simple_spectral'

real :: missing_value = -999.


contains


! ==================================================================================
! ==================================================================================


subroutine simple_spectral_rad_init(is, ie, js, je, num_levels, axes, Time, lonb, latb, dt_real)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real ,dimension(:,:),intent(in),optional :: lonb,latb !s Changed to 2d arrays as 2013 interpolator expects this.
real, intent(in)                  :: dt_real !s atmospheric timestep, used for radiation averaging
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/) ! indices for data which is dimensions of lat, lon, phalf
integer :: ierr, io, unit, v
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=simple_spectral_rad_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'simple_spectral_rad_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=simple_spectral_rad_nml)
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

initialized = .true.

allocate (nu_array         (nwavenumber))
allocate (kappa            (nwavenumber))

allocate (b_nu             (ie-is+1, je-js+1, num_levels, nwavenumber))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_solar        (ie-is+1, je-js+1, num_levels))

allocate (lw_tau_h2o_nu        (ie-is+1, je-js+1, num_levels+1, nwavenumber))
allocate (lw_tau_h2o_cont_nu        (ie-is+1, je-js+1, num_levels+1, nwavenumber))

if (do_co2) then
  allocate (lw_tau_co2_nu        (ie-is+1, je-js+1, num_levels+1, nwavenumber))
endif

allocate (lw_dtrans_nu     (ie-is+1, je-js+1, num_levels, nwavenumber))
allocate (lw_tau_nu        (ie-is+1, je-js+1, num_levels+1, nwavenumber))
allocate (sw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (lw_up_nu         (ie-is+1, je-js+1, num_levels+1, nwavenumber))
allocate (lw_down_nu       (ie-is+1, je-js+1, num_levels+1, nwavenumber))
allocate (lw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (lw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (lw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (rad_flux         (ie-is+1, je-js+1, num_levels+1))

allocate (b_surf           (ie-is+1, je-js+1))
allocate (b_surf_nu        (ie-is+1, je-js+1, nwavenumber))
allocate (sw_tau_0         (ie-is+1, je-js+1))
allocate (olr              (ie-is+1, je-js+1))
allocate (olr_spec         (ie-is+1, je-js+1, nwavenumber))
allocate (net_lw_surf      (ie-is+1, je-js+1))
allocate (toa_sw_in        (ie-is+1, je-js+1))

allocate (insolation       (ie-is+1, je-js+1))
allocate (p2               (ie-is+1, je-js+1))
allocate (coszen           (ie-is+1, je-js+1))
allocate (fracsun          (ie-is+1, je-js+1)) !jp from astronomy.f90 : fraction of sun on surface


! Initialize wavenumber increment
dwavenumber     = (wavenumber_max-wavenumber_min) / (nwavenumber - 1)

! Set up wavenumber array for LW calculations
do v = 1, nwavenumber
   nu_array(v) = wavenumber_min + dwavenumber*(v-1)
end do

!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    ! spectral diagostic fields
    id_ssm_bins_lw = diag_axis_init('ssm_bins_lw', nu_array, 'cm^-1', 'n', &
                    'ssm lw spectral bin centers', set_name='ssm_lw_bins')
                    
    id_spectral_olr = &
        register_diag_field ( mod_name, 'olr_spec', (/ axes(1:2), id_ssm_bins_lw/) , Time, &
        'SSM LW OLR spectrum', &
        'watts/m2', missing_value=missing_value )


    ! other diagnostics
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

                 
    if (mpp_pe() == 0) then
        write (*, *) 'Module name used for diagnostics:', mod_name
        write (*, *) 'Diagnostic registration results:'
        write (*, *) 'id_olr = ', id_olr
        write (*, *) 'id_spectral_olr = ', id_spectral_olr
        write (*, *) 'id_swdn_sfc = ', id_swdn_sfc
        write (*, *) 'id_swdn_sfc = ', id_swdn_sfc
        write (*, *) 'id_swdn_toa = ', id_swdn_toa
        write (*, *) 'id_net_lw_surf = ', id_net_lw_surf
        write (*, *) 'id_lwdn_sfc = ', id_lwdn_sfc
        write (*, *) 'id_lwup_sfc = ', id_lwup_sfc
    endif

return
end subroutine simple_spectral_rad_init

! ==================================================================================

subroutine simple_spectral_rad_down (is, js, Time_diag, lat, lon, p_full, p_half, t,         &
                           net_surf_sw_down, surf_lw_down, albedo, q)

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in), dimension(:,:)    :: lat, lon, albedo
real, intent(out), dimension(:,:)   :: net_surf_sw_down
real, intent(out), dimension(:,:)   :: surf_lw_down
real, intent(in), dimension(:,:,:)  :: t, q, p_full   ! num_levels  
real, intent(in), dimension(:,:,:)  :: p_half ! num_levels+1
integer :: i, j, k, n, v, dyofyr

integer :: seconds, year_in_s, days
real :: r_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, day_in_s, r_solday, r_total_seconds, r_days, r_dt_rad_avg, dt_rad_radians
logical :: used

! for debugging
character(len=8) :: dim1, dim2, dim3, dim4 
character(len=20) :: tmp, tmp2, tmp3, tmp4

real ,dimension(size(q,1),size(q,2),size(q,3)) :: co2f

n = size(t,3)
  
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

else
  ! Default: Averaged Earth insolation at all longitudes
  p2          = (1. - 3.*sin(lat)**2)/4.
  insolation  = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * sin(lat))
end if


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

! =================================================================================
! LONGWAVE RADIATION

if(do_read_co2)then
  call interpolator( co2_interp, Time_diag, p_half, co2f, trim(co2_variable_name))
  ! Needs maxval because co2f is a 3d array of a constant, so maxval just picks out one number
  carbon_conc = maxval(co2f) 
endif

! longwave source function
call planckFunction_3D(t, nu_array, b_nu)

! calculate optical depths
if (do_co2) then
  lw_tau_h2o_nu = 0.
  lw_tau_h2o_cont_nu = 0.
  lw_tau_co2_nu = 0.
  call calc_longwave_optical_depth_co2(p_half, t, q, nu_array, lw_tau_co2_nu)
  call calc_longwave_optical_depth_h2o(p_half, p_full, t, q, nu_array, lw_tau_h2o_nu) 
  call calc_longwave_optical_depth_h2o_cont(p_half, p_full, T, q, nu_array, lw_tau_h2o_cont_nu)
  lw_tau_nu = lw_tau_h2o_nu+lw_tau_co2_nu+lw_tau_h2o_cont_nu 
else ! H2O only
  lw_tau_h2o_nu = 0.
  lw_tau_h2o_cont_nu = 0.
  call calc_longwave_optical_depth_h2o(p_half, p_full, t, q, nu_array, lw_tau_h2o_nu)
  call calc_longwave_optical_depth_h2o_cont(p_half, p_full, T, q, nu_array, lw_tau_h2o_cont_nu) 
  lw_tau_nu = lw_tau_h2o_nu + lw_tau_h2o_cont_nu 
endif

! longwave differential transmissivity
! k=1 is model top
do k = 1, n
   lw_dtrans_nu(:,:,k,:) = exp( -(lw_tau_nu(:,:,k+1,:) - lw_tau_nu(:,:,k,:)) ) 
end do

! compute downward longwave flux by integrating downward
lw_down_nu = 0.
do k = 1, n
   lw_down_nu(:,:,k+1,:) = lw_down_nu(:,:,k,:)*lw_dtrans_nu(:,:,k,:) + pi*b_nu(:,:,k,:)*(1. - lw_dtrans_nu(:,:,k,:))
end do

! integrate over wavenumber to get total downwelling flux 
! (is this right?)
lw_down = sum(lw_down_nu, dim=4)*dwavenumber

! =================================================================================
surf_lw_down     = lw_down(:, :, n+1)
toa_sw_in        = sw_down(:, :, 1)
net_surf_sw_down = sw_down(:, :, n+1) * (1. - albedo)
! =================================================================================


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
end subroutine simple_spectral_rad_down

! ==================================================================================

subroutine simple_spectral_rad_up (is, js, Time_diag, lat, p_full, p_half, t_surf, t, tdt, albedo)

! Now complete the radiation calculation by computing the upward and net fluxes.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat, albedo
real, intent(in) , dimension(:,:)   :: t_surf
real, intent(in), dimension(:,:,:)  :: t, p_full ! num_levels  
real, intent(in), dimension(:,:,:)  :: p_half    ! num_levels+1
real, intent(inout), dimension(:,:,:) :: tdt

character(len=8) :: dim1, dim2, dim3, dim4 ! for debugging
character(len=20) :: tmp ! for debugging
integer, dimension(3) :: locn ! for debugging

integer :: i, j, k, v, n

logical :: used

n = size(t,3)

! longwave source function
call planckFunction_2D(t_surf, nu_array, b_surf_nu)

! integrate over wavenumber
b_surf = sum(b_surf_nu, dim=3)*dwavenumber
     
! compute upward longwave flux by integrating upward
lw_up_nu(:,:,n+1,:)    = pi*b_surf_nu

do k = n, 1, -1
    lw_up_nu(:,:,k,:)   = lw_up_nu(:,:,k+1,:)*lw_dtrans_nu(:,:,k,:) + pi*b_nu(:,:,k,:)*(1.0 - lw_dtrans_nu(:,:,k,:))
end do


! integrate over wavenumber to get total downwelling flux
! because delta wavenumber is constant, i can just sum across
! wavenumber, then multiply by dwavenumber (right?)
lw_up = sum(lw_up_nu, dim=4)*dwavenumber

! compute upward shortwave flux (here taken to be constant)
do k = 1, n+1
   sw_up(:,:,k)   = albedo(:,:) * sw_down(:,:,n+1)
end do

! net fluxes (positive up)
lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down
rad_flux = lw_flux + sw_flux

do k = 1, n
   tdt_rad(:,:,k)   = ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )
        
   tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )


   tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do


! ===========================DEBUGGING==========================
if (debug_ssm) then
  locn = minloc(tdt_rad)
  write(dim1, '(I8)') locn(1)
  write(dim2, '(I8)') locn(2)
  write(dim3, '(I8)') locn(3)
  write(tmp, '(E14.6)') minval(tdt_rad)*86400
  call error_mesg('simple_spectral_rad_up', 'minval(tdt_rad*86400)=: '//trim(adjustl(tmp)), NOTE)
  call error_mesg('simple_spectral_rad_up', 'at i,j,k=: '//trim(adjustl(dim1))//'x'//trim(adjustl(dim2))//'x'//trim(adjustl(dim3)), NOTE)
end if


olr         = lw_up(:,:,1)
olr_spec    = lw_up_nu(:,:,1,:)
net_lw_surf = lw_flux(:, :, n+1)
  
!------- outgoing lw flux toa (olr) -------
if ( id_olr > 0 ) then
   used = send_data ( id_olr, olr, Time_diag)
endif
!------- SPECTRAl outgoing lw flux toa (olr_spec) -------
if ( id_spectral_olr > 0 ) then
   used = send_data ( id_spectral_olr, olr_spec, Time_diag)
endif
!------- upward lw flux surface -------
if ( id_lwup_sfc > 0 ) then
   used = send_data ( id_lwup_sfc, pi*b_surf, Time_diag)
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
return
end subroutine simple_spectral_rad_up



! ==================================================================================

subroutine planckFunction_3D(T, nu, B)
  real, dimension(:,:,:), intent(in)    :: T
  real, dimension(:), intent(in)        :: nu
  real, dimension(:,:,:,:), intent(out) :: B
  integer :: k, v
  
  ! W/m2/sr/cm-1
  do v = 1, size(nu) 
    do k = 1, size(T,3)
      B(:,:,k,v) = 100*2*planck_h*((nu(v)*100)**3)*(lightspeed**2) / &
                  (exp( (planck_h * lightspeed * nu(v) * 100) / (boltzmann_k * T(:,:,k))) - 1)

    end do
  end do
end subroutine planckFunction_3D


! ==================================================================================

subroutine planckFunction_2D(Tsfc, nu, Bsfc)
  real, dimension(:,:), intent(in)    :: Tsfc
  real, dimension(:), intent(in)      :: nu
  real, dimension(:,:,:), intent(out) :: Bsfc
  integer :: v
  
  ! W/m2/sr/cm-1
  do v = 1, size(nu)
    Bsfc(:,:,v) = 100*2*planck_h*((nu(v)*100)**3)*(lightspeed**2) / &
                  (exp( (planck_h * lightspeed * nu(v) * 100) / (boltzmann_k * Tsfc)) - 1)

  end do
end subroutine planckFunction_2D


! ==================================================================================
subroutine calc_longwave_optical_depth_h2o_cont(p_half,p_full, T, q, nu, tau_h2o_cont) 
  real, dimension(:,:,:),   intent(in)  :: p_full, T, q
  real, dimension(:,:,:),   intent(in)  :: p_half
  real, dimension(:),       intent(in)  :: nu
  real, dimension(:,:,:,:), intent(out) :: tau_h2o_cont 

  real, dimension(size(T,1),size(T,2),size(T,3))          :: kappa_h2o_cont1
  real, dimension(size(T,1),size(T,2),size(T,3))          :: kappa_h2o_cont2
  real, dimension(size(T,1),size(T,2),size(T,3),size(nu)) :: integrand
  real, dimension(size(T,1),size(T,2),size(T,3))          :: vapor_pressure
  
  integer k,v

  !calc wv pressure
  do k = 1, size(T,3)
    vapor_pressure(:,:,k) = q(:,:,k) * p_full(:,:,k) / ( 0.622 + 0.378 * q(:,:,k) )
  end do
    
  kappa_h2o_cont1 = k_cnt_1*(vapor_pressure/pvstar_ref)*exp(a*(Tref-T))
  kappa_h2o_cont2 = k_cnt_2*(vapor_pressure/pvstar_ref)*exp(a*(Tref-T))
  
  tau_h2o_cont = 0.0
  do v = 1, size(nu)
    if (nu(v) <= 1700) then
      do k = 1, size(T, 3) ! T has num_levels, k=1 is model top
        tau_h2o_cont(:,:,k+1,v) = tau_h2o_cont(:,:,k,v) + kappa_h2o_cont1(:,:,k) * q(:,:,k) * (p_half(:,:,k+1) - p_half(:,:,k)) / grav
      end do
    else if (nu(v) > 1700) then
      do k = 1, size(T, 3) ! T has num_levels, k=1 is model top
        tau_h2o_cont(:,:,k+1,v) = tau_h2o_cont(:,:,k,v) + kappa_h2o_cont2(:,:,k) * q(:,:,k) * (p_half(:,:,k+1) - p_half(:,:,k)) / grav
      end do
    endif
  end do

  tau_h2o_cont = D * tau_h2o_cont

  end subroutine calc_longwave_optical_depth_h2o_cont 
! ==================================================================================

subroutine calc_longwave_optical_depth_h2o(p_half, p_full, T, q, nu, tau_h2o)
  real, dimension(:,:,:),   intent(in)  :: p_full, T, q ! num_levels
  real, dimension(:),       intent(in)  :: nu
  real, dimension(:,:,:),   intent(in)  :: p_half ! num_levels+1
  real, dimension(:,:,:,:), intent(out) :: tau_h2o ! num_levels+1

  ! local variables
  real, dimension(size(nu))                               :: kappa_h2o
  real, dimension(size(T,1),size(T,2),size(T,3),size(nu)) :: kappa_h2o_broadened
  real, dimension(size(T,1),size(T,2),size(T,3),size(nu)) :: integrand
  
  real :: p_ref = 5.0e4 ! p_full in pascal, so be consistent
  integer :: k, v

  ! H2O line absorption
  kappa_h2o = 0.0
  do v = 1, size(nu)
    if (nu(v) <= nu_rot) then
        kappa_h2o(v) = 37
    else if (nu(v) > nu_rot .AND. nu(v) <= 1000.0) then
        kappa_h2o(v) = k_rot * exp(- (nu(v)-nu_rot) / l_rot )
    else if (nu(v) > 1000.0 .AND. nu(v) <= nu_vr) then
        kappa_h2o(v) = k_vr * exp(- (nu_vr-nu(v)) / l_vr )
    else if (nu(v) > nu_vr .AND. nu(v) <= nu_ad) then
        kappa_h2o(v) = 5
    else if (nu(v) > nu_ad .AND. nu(v) <= 2500) then
        kappa_h2o(v) = k_ad * exp(- (nu(v)-nu_ad) / l_ad )
    else if (nu(v) > 2500.) then
        kappa_h2o(v) = 0
    endif
  end do
  
  ! Calculate broadened absorption coefficient
  ! and dtau/dp=integrand, defined at full levels
  do v = 1, size(nu)
    do k = 1, size(q, 3) ! q has num_levels
      if (do_pressure_broadening) then
        kappa_h2o_broadened(:,:,k,v) = kappa_h2o(v) * (p_full(:,:,k)/p_ref)
        integrand(:,:,k,v) = (kappa_h2o_broadened(:,:,k,v)*q(:,:,k)/grav)
      else
        integrand(:,:,k,v) = (kappa_h2o(v)*q(:,:,k)/grav)
      endif
    end do
  end do
  
  ! Calculate optical depth
  ! tau(next_half_level) = tau(last_half_level) + integrand * deltaPressure(between_half_levels) 
  tau_h2o = 0.0
  do v = 1, size(nu)
    do k = 1, size(T, 3) ! t has num_levels, k=1 is model top
      tau_h2o(:,:,k+1,v) = tau_h2o(:,:,k,v) + integrand(:,:,k,v) * (p_half(:,:,k+1) - p_half(:,:,k))
    end do
  end do

  ! multiply by diffusivity factor
  tau_h2o = D * tau_h2o
  
end subroutine calc_longwave_optical_depth_h2o 

! ==================================================================================

subroutine calc_longwave_optical_depth_co2(p_half, T, q, nu, tau_co2)
  real, dimension(:,:,:),   intent(in)  :: p_half ! num_levels+1
  real, dimension(:,:,:),   intent(in)  :: T, q   ! num_levels
  real, dimension(:),       intent(in)  :: nu
  real, dimension(:,:,:,:), intent(out) :: tau_co2

  ! local variables
  real, dimension(size(nu))                                              :: kappa_co2
  real, dimension(size(p_half,1),size(p_half,2),size(p_half,3),size(nu)) :: integrand ! not used???? 
  
  real :: p_ref = 5.0e4 ! p_half in pascal, so be consistent
  real :: q_co2
  integer :: k, v

  ! debugging
  character(len=8) :: dim1, dim2, dim3, dim4 
  character(len=20) :: tmp 

  ! kg/kg co2 from ppmv
  q_co2 = carbon_conc*1.0e-6*44/29

  ! CO2 line absorption
  kappa_co2 = 0.0
  do v = 1, size(nu)
    if (nu(v) > 500. .AND. nu(v) < 850.) then
      kappa_co2(v) = kappa_co2(v) + k_co2 * exp(- abs(nu(v)-nu_co2) / l_co2 ) 
    endif
  end do

  ! Calculate optical depth 
  ! eqn 26 of Jeevanjee and Fueglistaler 2020, assumes pressure broadening,
  ! include diffusivity factor (they forgot)
  tau_co2 = 0.0d0
  do v = 1, size(nu)
    do k = 1, size(p_half, 3)
      tau_co2(:,:,k,v) =  D * kappa_co2(v) * q_co2 * p_half(:,:,k)**2 / (2 * grav * p_ref)
    end do
  end do

end subroutine calc_longwave_optical_depth_co2 

! ==================================================================================

subroutine simple_spectral_rad_end

deallocate (b_nu, tdt_rad, tdt_solar)
deallocate (lw_up, lw_up_nu, lw_down, lw_down_nu, lw_flux, sw_up, sw_down, sw_flux, rad_flux)
deallocate (b_surf, b_surf_nu, olr, olr_spec, net_lw_surf, toa_sw_in, sw_tau_0)
deallocate (lw_dtrans_nu, lw_tau_nu, lw_tau_h2o_nu, sw_tau)
deallocate (insolation, p2) 
  
if (do_co2) then
  deallocate (lw_tau_co2_nu)
endif

if (do_read_co2) call interpolator_end(co2_interp)

end subroutine simple_spectral_rad_end

! ==================================================================================

end module simple_spectral_rad_mod
