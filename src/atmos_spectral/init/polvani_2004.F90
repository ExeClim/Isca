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

module polvani_2004_mod

! Polvani, L. M., R. K. Scott, and S. J. Thomas, 2004:
! Numerically Converged Solutions of the Global Primitive Equations
! for Testing the Dynamical Core of Atmospheric GCMs.
! Monthly Weather Review, 132, 2539-2552.

use               mpp_mod, only: mpp_error, FATAL

use               fms_mod, only: mpp_pe, mpp_root_pe, stdlog, &
                                 write_version_number, close_file, check_nml_error, open_namelist_file

use         constants_mod, only: PI,    &
                                 GRAV,  & ! GRAV=9.800, not 9.806 as specified by Polvani
                                 RDGAS, & ! RDGAS=287.04, not 287.00 as specified by Polvani
                                 OMEGA, & ! same as specified by Polvani
                                 RADIUS,& ! same as specified by Polvani
                                 KAPPA    ! same as specified by Polvani

use   vert_coordinate_mod, only: compute_vert_coord

use        transforms_mod, only: get_grid_boundaries, get_sin_lat, get_cos_lat, trans_grid_to_spherical, trans_spherical_to_grid,&
                                 get_wts_lat, get_lon_max, get_lat_max, get_deg_lon, get_deg_lat, &
                                 vor_div_from_uv_grid, uv_grid_from_vor_div, get_grid_domain

use  press_and_geopot_mod, only: press_and_geopot_init, pressure_variables

use      diag_manager_mod, only: diag_axis_init, register_static_field, send_data

implicit none
private

character(len=128), parameter :: version = &
'$Id: polvani_2004.F90,v 19.0 2012/01/06 20:35:35 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: siena_201211 $'

public :: polvani_2004

integer, parameter :: num_standard_levels = 8
real, parameter, dimension(num_standard_levels  ) :: z_standard=(/0., 11.e3, 20.e3, 32.e3, 47.e3, 51.e3, 71.e3, 80.e3/)
real, parameter, dimension(num_standard_levels-1) :: lapse_standard=(/-6.5e-3, 0.0, 1.e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3/)
real,            dimension(num_standard_levels  ) :: T_standard

real :: lambda0, phi0, alpha, beta ! These determine the position of the initial temperature perturbation.
real, parameter :: P00=1.e5 ! Used to compute potential temperature. sea_level_press may be the same, but one should not assume that.

real, parameter :: twopi = 2*PI
real :: xx
integer :: kk

real :: H  = 7.340e3      ! meters
real :: z0 = 22.e3        ! meters
real :: delta_z0 = 5.e3   ! meters
real :: z1 = 30.e3        ! meters
real :: u0 = 50.          ! m/sec
real :: perturb_amp = 1.0 ! deg K

namelist / polvani_2004_nml / H, z0, delta_z0, z1, u0, perturb_amp

contains

!=========================================================================================================================

subroutine polvani_2004(sea_level_press, triang_trunc, vert_difference_option, pk, bk, &
                   vors, divs, ts, ln_ps, ug, vg, tg, psg, vorg, divg, surf_geopotential)

real,    intent(in) :: sea_level_press
logical, intent(in) :: triang_trunc 
character(len=*), intent(in) :: vert_difference_option
real,    intent(out), dimension(:) :: pk, bk
complex, intent(out), dimension(:,:,:) :: vors, divs, ts
complex, intent(out), dimension(:,:  ) :: ln_ps
real,    intent(out), dimension(:,:,:) :: ug, vg, tg
real,    intent(out), dimension(:,:  ) :: psg
real,    intent(out), dimension(:,:,:) :: vorg, divg
real,    intent(out), dimension(:,:  ) :: surf_geopotential

integer :: unit, ierr, io, i, j, k, ks, lon_max, lat_max, num_levels
integer :: id_lon, id_lat, id_phalf, id_pfull, id_basic_flow, id_term1_eq10, id_basic_temp, id_pot_temp, id_perturbation
integer :: is, ie, js, je
logical :: used
real, dimension(size(ug,1), size(ug,2)) :: ln_psg
real, dimension(size(ug,3)  ) :: p_full, ln_p_full, T0
real, dimension(size(ug,3)+1) :: p_half, ln_p_half
real :: F, du_dz, dF_dz, ff1, ff2, dzz1_dz, dzz2_dz, dff1_dz, dff2_dz, dTdy, dTdy_prev, ln_slp
real, dimension(size(ug,3)) :: z, zz1, zz2, gmean_term1

! The allocatable arrays below must be allocated for the global domain
real, allocatable, dimension(:)   :: deg_lon, lonb, deg_lat, latb, rad_lon, rad_lat, wts_lat, sin_lat, cos_lat, tan_lat, coriolis
real, allocatable, dimension(:,:) :: basic_flow, term1_eq10, basic_temp, pot_temp
real, allocatable, dimension(:)   :: lon_factor, lat_factor
real, allocatable, dimension(:,:) :: perturbation

!------------------------------------------------------------------------------------------------

unit = open_namelist_file()
ierr=1
do while (ierr /= 0)
  read(unit, nml=polvani_2004_nml, iostat=io, end=20)
  ierr = check_nml_error (io, 'polvani_2004_nml')
enddo
20 call close_file (unit)
call write_version_number(version, tagname)
if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=polvani_2004_nml)

T_standard(1) = 288.15
do ks=2,num_standard_levels
  T_standard(ks) = T_standard(ks-1) + lapse_standard(ks-1)*(z_standard(ks)-z_standard(ks-1))
enddo

call compute_vert_coord ('even_sigma', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pk, bk)

surf_geopotential = 0.0
call press_and_geopot_init(pk, bk, .false., vert_difference_option)
call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, sea_level_press)

call get_lon_max(lon_max)
call get_lat_max(lat_max)
num_levels = size(ug,3)
if(num_levels /= 20) then
  call mpp_error(FATAL,'Number of vertical levels must be 20 when using polvani_2004 initialization')
endif
allocate(rad_lon(lon_max), deg_lon(lon_max), lon_factor(lon_max))
allocate(rad_lat(lat_max), deg_lat(lat_max), lat_factor(lat_max))
allocate(wts_lat(lat_max), sin_lat(lat_max), cos_lat(lat_max), tan_lat(lat_max))
allocate(coriolis(lat_max), latb(lat_max+1), lonb(lon_max+1))
allocate(basic_flow(lat_max,num_levels), term1_eq10(lat_max,num_levels))
allocate(basic_temp(lat_max,num_levels),   pot_temp(lat_max,num_levels))
allocate(perturbation(lon_max,lat_max))
call get_deg_lon(deg_lon)
rad_lon = PI*deg_lon/180.
call get_deg_lat(deg_lat)
rad_lat = PI*deg_lat/180.
call get_wts_lat(wts_lat)
call get_sin_lat(sin_lat)
call get_cos_lat(cos_lat)
tan_lat = sin_lat/cos_lat
call get_grid_boundaries(lonb, latb, global=.true.) ! units = radians
coriolis = 2*OMEGA*sin_lat

ln_slp = log(sea_level_press)
do k=1,num_levels
  z(k) = H*(ln_slp - ln_p_full(k))
  zz1(k) = (z(k)-z0)/delta_z0
  zz2(k) = PI*z(k)/z1
enddo

do k=1,num_levels
  if(z(k) > z_standard(num_standard_levels)) then
    T0(k) = T_standard(num_standard_levels)
    cycle
  endif
  do ks=2,num_standard_levels
    if(z(k) <= z_standard(ks) .and. z(k) >= z_standard(ks-1)) then
      T0(k)=((z_standard(ks)-z(k))*T_standard(ks-1) + (z(k)-z_standard(ks-1))*T_standard(ks))/(z_standard(ks)-z_standard(ks-1))
      exit
    else
      cycle
    endif
  enddo
enddo

do k=1,num_levels
  ff1 = 1 - tanh(zz1(k))**3
  ff2 = sin(zz2(k))
  F = .5*ff1*ff2
  do j=1,lat_max
    if(sin_lat(j) > 0.0) then
      basic_flow(j,k) = u0*F*(sin(PI*sin_lat(j)**2))**3
    else
      basic_flow(j,k) = 0.0
    endif
  enddo
enddo

! The k loop below computes the integral of equation 10 in the paper.
! The inner loop over latitude does the actual integration.
! Prior to that the vertical derivative of F (equation 6) must be computed.
do k=1,num_levels
  ff1 = 1 - tanh(zz1(k))**3
  ff2 = sin(zz2(k))
  dzz1_dz = 1/delta_z0
  dzz2_dz = PI/z1
  dff1_dz = -3 * ( tanh(zz1(k))/cosh(zz1(k)) )**2 * dzz1_dz
  dff2_dz = cos(zz2(k)) * dzz2_dz
  dF_dz = .5*(ff1*dff2_dz + dff1_dz*ff2)

  gmean_term1(k) = 0.0
  do j=1,lat_max
    if(sin_lat(j) > 0.0) then
      du_dz = u0 * (sin(PI*sin_lat(j)**2))**3 * dF_dz
    else
      du_dz = 0.0
    endif
    dTdy = -(H/RDGAS) * (RADIUS*coriolis(j) + 2*basic_flow(j,k)*tan_lat(j)) * du_dz
    if(j == 1) then
!     term1_eq10(j,k) =    .5*wts_lat(j)*dTdy/cos_lat(j)
      term1_eq10(j,k) = (rad_lat(j)-latb(j))*dTdy/cos_lat(j)
    else
!     term1_eq10(j,k) = term1_eq10(j-1,k) + .5*(wts_lat(j-1)*dTdy/cos_lat(j-1) + wts_lat(j)*dTdy/cos_lat(j))
      term1_eq10(j,k) = term1_eq10(j-1,k) + (latb(j)-rad_lat(j-1))*dTdy_prev + (rad_lat(j)-latb(j))*dTdy
    endif
    gmean_term1(k) = gmean_term1(k) + .5*wts_lat(j)*term1_eq10(j,k)
    dTdy_prev = dTdy
  enddo
enddo

! The temperature perturbation is computed below.
lambda0 = 0.0
phi0    = PI/4.
alpha   = 1./3.
beta    = 1./6.
do i=1,lon_max
  xx = rad_lon(i)-lambda0
  ! Make sure xx falls in the range -PI to +PI.
  kk = nint(xx/twopi)
  xx = xx - twopi*kk
  lon_factor(i) = 1./cosh(xx/alpha)**2
enddo
do j=1,lat_max
  lat_factor(j) = 1./cosh((rad_lat(j)-phi0)/beta)**2
enddo
do j=1,lat_max
  do i=1,lon_max
    perturbation(i,j) = perturb_amp*lon_factor(i)*lat_factor(j)
  enddo
enddo

id_phalf = diag_axis_init('phalf',.01*p_half,'hPa','z','approx half pressure level',direction=-1)
id_pfull = diag_axis_init('pfull',.01*p_full,'hPa','z','approx full pressure level',direction=-1,edges=id_phalf)
id_lon   = diag_axis_init('lon',deg_lon,'degrees E','x','longitude')
id_lat   = diag_axis_init('lat',deg_lat,'degrees N','y','latitude')
id_basic_flow   = register_static_field('polvani_2004', 'basic_flow',  (/id_lat,id_pfull/),'initial zonal wind', 'm/sec')
id_term1_eq10   = register_static_field('polvani_2004', 'term1_eq10',  (/id_lat,id_pfull/),'first term of equation 10', 'K')
id_basic_temp   = register_static_field('polvani_2004', 'basic_temp',  (/id_lat,id_pfull/),'initial basic_temp', 'K')
id_pot_temp     = register_static_field('polvani_2004', 'pot_temp',    (/id_lat,id_pfull/),'initial potential basic_temp', 'K')
id_perturbation = register_static_field('polvani_2004', 'perturbation',(/id_lon,id_lat/),  'initial temperature perturbation','K')
if(id_basic_flow > 0) used = send_data(id_basic_flow, basic_flow)
if(id_term1_eq10 > 0) used = send_data(id_term1_eq10, term1_eq10)
do k=1,num_levels
  do j=1,lat_max
    basic_temp(j,k) = term1_eq10(j,k) - gmean_term1(k) + T0(k)
    pot_temp(j,k) = basic_temp(j,k)*(P00/p_full(k))**KAPPA
  enddo
enddo
if(id_basic_temp > 0) used = send_data(id_basic_temp, basic_temp)
if(id_pot_temp   > 0) used = send_data(id_pot_temp,     pot_temp)
if(id_perturbation>0) used = send_data(id_perturbation, perturbation)

call get_grid_domain(is, ie, js, je)

do k=1,num_levels
  do i=1,size(ug,1)
    ug(i,:,k) = basic_flow(js:je,k)
    tg(i,:,k) = basic_temp(js:je,k)
  enddo
enddo
do j=1,size(ug,2)
  do i=1,size(ug,1)
    tg(i,j,:) = tg(i,j,:) + perturbation(is+i-1,js+j-1)
  enddo
enddo
vg = 0.0

! Transform grid fields to spherical waves then transform
! back to ensure that spectral and grid representations match.

call trans_grid_to_spherical(tg, ts)
call trans_spherical_to_grid(ts, tg)
call vor_div_from_uv_grid(ug, vg, vors, divs, triang=triang_trunc)
call uv_grid_from_vor_div(vors, divs, ug, vg)
call trans_spherical_to_grid(vors, vorg)
call trans_spherical_to_grid(divs, divg)
psg = sea_level_press
ln_psg = alog(psg)
call trans_grid_to_spherical(ln_psg, ln_ps)
call trans_spherical_to_grid(ln_ps,  ln_psg)
psg = exp(ln_psg)

deallocate(deg_lon, deg_lat, rad_lon, rad_lat, wts_lat, sin_lat, cos_lat, tan_lat, coriolis)
deallocate(basic_flow, term1_eq10, basic_temp, pot_temp)
deallocate(lon_factor, lat_factor)
deallocate(perturbation)

return
end subroutine polvani_2004
!================================================================================

end module polvani_2004_mod
