module stirring_mod

!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! This program is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------

! Stirring is computed as described in the following paper:

! Vallis, Gerber, Kushner, Cash, 2003: A Mechanism and Simple Dynamical Model of the North Atlantic Oscillation and Annular Modes.
! J. Atmos. Sci., 61, 264-280.

! Stirring is not part of barotropic_physics because barotropic_physics appears to be intended for
! operations that are done completely in grid space. Stirring is computed partly in spectral space.

use    constants_mod, only: pi

use time_manager_mod, only: time_type

use          fms_mod, only: open_namelist_file, check_nml_error, close_file, write_version_number, &
                            stdlog, mpp_pe, mpp_root_pe, file_exist, read_data, write_data, error_mesg, FATAL

use   transforms_mod, only: get_spec_domain, get_grid_domain, trans_spherical_to_grid, trans_grid_to_spherical, &
                            grid_domain, get_lon_max, get_lat_max, get_deg_lon, get_deg_lat, get_grid_boundaries, &
                            get_num_fourier, get_num_spherical, spectral_domain

use diag_manager_mod, only: diag_axis_init, register_static_field, register_diag_field, send_data

implicit none
private

integer :: ms,me,ns,ne,is,ie,js,je
integer :: id_str_amp, id_g_stir_sqr, id_stir
logical :: used
logical, allocatable, dimension(:,:) :: wave_mask   ! wave_mask(m,n) = .true. if spherical wave (m,n) is to be excited
complex, allocatable, dimension(:,:) :: s_stir      ! stirring. Saved from one time step to the next
real,    allocatable, dimension(:,:) :: localize    ! localizes the stirring
real,    allocatable, dimension(:,:) :: g_stir_sqr  ! time mean of g_stir**2 over entire integration
integer, allocatable, dimension(:)   :: seed        ! random number seed
real ::  astir, bstir
integer :: num_steps, num_fourier, num_spherical, nseed

logical :: module_is_initialized = .false.

character(len=128) :: version = '$Id: stirring.F90,v 17.0 2009/07/21 03:00:25 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'

public :: stirring_init, stirring, stirring_end

real :: decay_time=2*86400, amplitude=0.0, lat0=45., widthy=12.
logical :: do_localize=.true.!Default true to allow forcing to be localized in physical space. Set to false to have forcing everywhere.

! Set B to a non-zero value for stirring that has zonal structure.
! The strength of the stirring at latitude=lat0 is: amplitude*(1.0 + B*exp(-.5*((lon-lon0)/widthx)**2))
real :: lon0=180., B=0.0, widthx=45., C=1.0 ! widthx
integer :: n_total_forcing_max = 15 !total wavenumbers LESS THAN this number will be forced
integer :: n_total_forcing_min = 9 !total wavenumbers GREATER THAN this number will be forced
integer :: zonal_forcing_min = 3 !Zonal wavenumbers GREATER THAN this number will be forced, subject to total wavenumber constraints

namelist / stirring_nml / decay_time, amplitude, lat0, lon0, widthy, widthx, B, do_localize, n_total_forcing_max, n_total_forcing_min, zonal_forcing_min

contains

!================================================================================================================================
subroutine stirring_init(dt, Time, id_lon, id_lat, id_lonb, id_latb)
real, intent(in) :: dt
type(time_type), intent(in) :: Time
integer, intent(in) :: id_lon, id_lat, id_lonb, id_latb
real :: xx, kk, rad_to_deg
integer :: i,j,m,n,ierr,io,unit,lon_max,lat_max
real, allocatable, dimension(:) :: ampx, ampy, lon, lat, lonb, latb
real, allocatable, dimension(:,:) :: real_part, imag_part

if(module_is_initialized) return

call write_version_number (version, tagname)

if (file_exist('input.nml')) then
  unit = open_namelist_file ()   
  ierr=1                         
  do while (ierr /= 0)           
    read  (unit, nml=stirring_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'stirring_nml')
  enddo                          
  10 call close_file (unit)      
endif
if(mpp_pe() == mpp_root_pe()) write(stdlog(), nml=stirring_nml)

call get_lon_max(lon_max)
call get_lat_max(lat_max)

allocate(lon (lon_max  )) ; lon  = 0.0
allocate(lat (lat_max  )) ; lat  = 0.0
allocate(lonb(lon_max+1)) ; lonb = 0.0
allocate(latb(lat_max+1)) ; latb = 0.0

call get_deg_lon(lon)          
call get_deg_lat(lat)

module_is_initialized = .true.
if(amplitude == 0.0) return ! stirring does nothing more unless amplitude is non-zero

call get_spec_domain(ms,me,ns,ne)
call get_grid_domain(is,ie,js,je)
call get_num_fourier(num_fourier)
call get_num_spherical(num_spherical)

allocate(wave_mask(ms:me,ns:ne)); wave_mask = .false.
allocate(s_stir(ms:me,ns:ne)); s_stir = cmplx(0.0,0.0)
allocate(ampx(is:ie)); ampx = 0.0
allocate(ampy(js:je)); ampy = 0.0
allocate(localize(is:ie,js:je)); localize = 0.0
allocate(g_stir_sqr(is:ie,js:je)); g_stir_sqr = 0.0

! wave_mask is .true. when  (m+n > 9) .and. (m+n < 15) .and. (m > 3)
do m=(zonal_forcing_min+1),(n_total_forcing_max-1)
  if(m >= ms .and. m <= me) then
    do n=(n_total_forcing_min+1)-m,(n_total_forcing_max-1)-m
      if(n >= ns .and. n <= ne) then
        wave_mask(m,n) = .true.
      endif
    enddo
  endif
enddo

astir = sqrt(1.0 - exp(-2*dt/decay_time))
bstir = exp(-dt/decay_time)

do i=is,ie
  xx = lon(i)-lon0
  ! Make sure xx falls in the range -180. to +180.
  kk = nint(xx/360.)
  xx = xx - 360.*kk
  ampx(i) = (1 + B*exp(-.5*(xx/widthx)**2))
enddo
do j=js,je
  ampy(j) = exp(-.5*((lat(j)-lat0)/widthy)**2)
enddo
if (do_localize) then
    do j=js,je
        do i=is,ie
            localize(i,j) = ampx(i)*ampy(j)
        enddo
    enddo
else
    localize = 1.0
endif

deallocate(ampx, ampy)

num_steps = 0
id_g_stir_sqr = register_static_field('stirring_mod', 'stirring_sqr', (/id_lon,id_lat/), 'stirring sqrared', '1/sec^4')
id_str_amp    = register_static_field('stirring_mod', 'stirring_amp', (/id_lon,id_lat/), 'amplitude of stirring', 'none')
id_stir       = register_diag_field  ('stirring_mod', 'stirring',     (/id_lon,id_lat/), Time, 'stirring', '1/sec^2')
used = send_data(id_str_amp, amplitude*localize)

call random_seed(size=nseed)
allocate(seed(nseed))

if(file_exist('INPUT/stirring.res.nc')) then
  allocate(real_part(ms:me,ns:ne), imag_part(ms:me,ns:ne))
  call read_data('INPUT/stirring.res.nc', 'stir_real', real_part, spectral_domain)
  call read_data('INPUT/stirring.res.nc', 'stir_imag', imag_part, spectral_domain)
  do n=ns,ne
    do m=ms,me
      s_stir(m,n) = cmplx(real_part(m,n),imag_part(m,n))
    end do
  end do
  deallocate(real_part, imag_part)
  call read_data('INPUT/stirring.res.nc', 'ran_nmbr_seed', seed, no_domain=.true.)
  call random_seed(put=seed)
endif

end subroutine stirring_init
!================================================================================================================================
subroutine stirring(Time, dt_vors)
type(time_type), intent(in) :: Time
complex, dimension(ms:me,ns:ne), intent(inout) :: dt_vors
real,    dimension(is:ie,js:je) :: g_stir
complex, dimension(ms:me,ns:ne) :: new_stirring
real,    dimension(0:num_fourier,0:num_spherical,2) :: ran_nmbrs
integer :: i,j,m,n
real :: x,y

if(.not.module_is_initialized) then
  call error_mesg('stirring', 'stirring_init has not been called', FATAL)
end if

if(amplitude == 0.0) return ! stirring does nothing unless amplitude is non-zero

call random_number(ran_nmbrs)

do n=ns,ne
do m=ms,me
  if(wave_mask(m,n)) then
    new_stirring(m,n) = amplitude*astir*cmplx(2*ran_nmbrs(m,n,1)-1, 2*ran_nmbrs(m,n,2)-1)
  else
    new_stirring(m,n) = cmplx(0.0,0.0)
  endif
enddo
enddo
call trans_spherical_to_grid(new_stirring,g_stir)
g_stir = localize*g_stir
call trans_grid_to_spherical(g_stir,new_stirring)
if(ms == 0 .and. ns == 0) then
  new_stirring(0,0)=cmplx(0.0,0.0) ! A non-zero global mean is introduced by the grid space computation, but we don't want it.
endif
s_stir = bstir*s_stir + new_stirring !This is equation A.6 in Vallis et al 2004 - DOI:10.1175/1520-0469(2004)061<0264:AMASDM>2.0.CO;2

dt_vors = dt_vors + s_stir
call trans_spherical_to_grid(s_stir,g_stir)
g_stir_sqr = g_stir_sqr + g_stir*g_stir
num_steps = num_steps + 1
used = send_data(id_stir, g_stir, Time)

end subroutine stirring
!================================================================================================================================
subroutine stirring_end

if(.not.module_is_initialized) return

if(amplitude == 0.0) return ! stirring does nothing unless amplitude is non-zero

g_stir_sqr = g_stir_sqr/num_steps
used = send_data(id_g_stir_sqr, g_stir_sqr)

call write_data('RESTART/stirring.res.nc', 'stir_real',  real(s_stir), spectral_domain)
call write_data('RESTART/stirring.res.nc', 'stir_imag', aimag(s_stir), spectral_domain)
call random_seed(get=seed)
call write_data('RESTART/stirring.res.nc', 'ran_nmbr_seed', seed, no_domain=.true.)

deallocate(wave_mask, s_stir, localize, g_stir_sqr)
module_is_initialized = .false.

end subroutine stirring_end
!================================================================================================================================

end module stirring_mod
