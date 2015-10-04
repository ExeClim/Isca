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

module jablonowski_2006_mod

! Christiane Jablonwoski and David Williamson, 2006:
! A baroclinic instability test case for atmospheric model dynamical cores
! Quarterly J. Roy. Met. Soc., Vol. 132, October Part C, No. 621C, 2943-2975

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use               fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, stdlog, &
                                 write_version_number, close_file, check_nml_error

use         constants_mod, only: PI,    &
                                 GRAV,  & ! GRAV=9.800, not 9.806 as specified by Polvani
                                 RDGAS, & ! RDGAS=287.04, not 287.00 as specified by Polvani
                                 OMEGA, & ! same as specified by Polvani
                                 RADIUS,& ! same as specified by Polvani
                                 KAPPA    ! same as specified by Polvani

use   vert_coordinate_mod, only: compute_vert_coord

use        transforms_mod, only: get_grid_boundaries, get_sin_lat, get_cos_lat, trans_grid_to_spherical, trans_spherical_to_grid,&
                                 get_deg_lon, get_deg_lat, &
                                 vor_div_from_uv_grid, uv_grid_from_vor_div, get_grid_domain

use  press_and_geopot_mod, only: press_and_geopot_init, pressure_variables

use      diag_manager_mod, only: diag_axis_init, register_static_field, send_data

implicit none
private

character(len=128), parameter :: version = &
'$Id: jablonowski_2006.F90,v 19.0 2012/01/06 20:35:33 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: siena_201211 $'

public :: jablonowski_2006

real, parameter :: P00=1.e5 ! Used to compute potential temperature. sea_level_press may be the same, but one should not assume that.

real, parameter :: halfpi=.5*PI
real :: xx, sin_latc, cos_latc, nv_surf

real :: n0    = .252 ! sigma level of jet max
real :: U0    =  35. ! max jet speed (m/sec)
real :: nt    = .20  ! sigma level of tropopause
real :: lapse = .005 ! lapse rate below tropopause (deg K/m)
real :: T0    = 288. ! global mean surface temperature (deg K)
real :: Up    =   1. ! amplitude of perturbation (m/sec)
real :: lonc  =  20. ! longitude center of perturbation (deg E)
real :: latc  =  40. ! latitude  center of perturbation (deg N)
real :: deltaT = 4.8e5 ! determines lapse rate above tropopause

namelist / jablonowski_2006_nml / n0, U0, nt, lapse, T0, Up, lonc, latc, deltaT

contains

!=========================================================================================================================

subroutine jablonowski_2006(sea_level_press, triang_trunc, &
                   vert_coord_option, vert_difference_option, scale_heights, surf_res, &
                   p_press, p_sigma, exponent, pk, bk, &
                   vors, divs, ts, ln_ps, ug, vg, tg, psg, vorg, divg, surf_geopotential)

real,    intent(in) :: sea_level_press
logical, intent(in) :: triang_trunc 
character(len=*), intent(in) :: vert_coord_option, vert_difference_option
real,    intent(in) :: scale_heights, surf_res, p_press, p_sigma, exponent
real,    intent(out), dimension(:) :: pk, bk
complex, intent(out), dimension(:,:,:) :: vors, divs, ts
complex, intent(out), dimension(:,:  ) :: ln_ps
real,    intent(out), dimension(:,:,:) :: ug, vg, tg
real,    intent(out), dimension(:,:  ) :: psg
real,    intent(out), dimension(:,:,:) :: vorg, divg
real,    intent(out), dimension(:,:  ) :: surf_geopotential

integer :: unit, ierr, io, i, j, k, num_lon, num_lat, num_levels
integer :: id_lon, id_lat, id_phalf, id_pfull, id_basic_flow, id_basic_temp, id_pot_temp, id_perturbation
logical :: used
real, dimension(size(ug,1), size(ug,2)) :: ln_psg
real, dimension(size(ug,3)  ) :: p_full, ln_p_full, sigma
real, dimension(size(ug,3)+1) :: p_half, ln_p_half
real, dimension(size(ug,3)) :: nv

! The allocatable arrays below must be allocated for the global domain
real, allocatable, dimension(:)   :: deg_lon, deg_lat, rad_lon, rad_lat, sin_lat, cos_lat, coriolis
real, allocatable, dimension(:,:) :: basic_flow, basic_temp, pot_temp
real, allocatable, dimension(:)   :: lon_factor1, lat_factor1, lat_factor2, z_factor1, z_factor2, z_factor3
real, allocatable, dimension(:,:) :: perturbation

!------------------------------------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=jablonowski_2006_nml, iostat=io)
    ierr = check_nml_error(io, 'jablonowski_2006_nml')
#else
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=jablonowski_2006_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'jablonowski_2006_nml')
    enddo
20  call close_file (unit)
#endif
call write_version_number(version, tagname)
if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=jablonowski_2006_nml)

call compute_vert_coord(vert_coord_option, scale_heights, surf_res, exponent, p_press, p_sigma, sea_level_press, pk, bk)

call press_and_geopot_init(pk, bk, .false., vert_difference_option)
call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, sea_level_press)
sigma = p_full/sea_level_press
nv = (sigma - n0)*halfpi

num_lon    = size(ug,1)
num_lat    = size(ug,2)
num_levels = size(ug,3)
allocate(rad_lon(num_lon), deg_lon(num_lon), lon_factor1(num_lon))
allocate(z_factor1(num_levels), z_factor2(num_levels), z_factor3(num_levels))
allocate(rad_lat(num_lat), deg_lat(num_lat), lat_factor1(num_lat), lat_factor2(num_lat))
allocate(sin_lat(num_lat), cos_lat(num_lat), coriolis(num_lat))
allocate(basic_flow(num_lat,num_levels), perturbation(num_lon,num_lat))
allocate(basic_temp(num_lat,num_levels), pot_temp(num_lat,num_levels))

call get_deg_lon(deg_lon)
rad_lon = PI*deg_lon/180.
call get_deg_lat(deg_lat)
rad_lat = PI*deg_lat/180.
call get_sin_lat(sin_lat)
call get_cos_lat(cos_lat)
coriolis = 2*OMEGA*sin_lat

xx = RDGAS*lapse/GRAV
do k=1,num_levels
  z_factor2(k) = U0*cos(nv(k))**1.5
  do j=1,num_lat
    basic_flow(j,k) = z_factor2(k)*sin(2*rad_lat(j))**2
  enddo
  z_factor3(k) = .75*PI*U0*sigma(k)*sin(nv(k))*sqrt(cos(nv(k)))/RDGAS
  if(sigma(k) <= nt) then
    z_factor1(k) = T0*sigma(k)**xx + deltaT*(nt-sigma(k))**5
  else
    z_factor1(k) = T0*sigma(k)**xx
  endif
enddo

do j=1,num_lat
  lat_factor1(j) = 10./63. - 2*sin_lat(j)**6 * (cos_lat(j)**2 + 1./3.)
  lat_factor2(j) = RADIUS*OMEGA*(1.6*cos_lat(j)**3 * (sin_lat(j)**2 + 2./3.) - .25*PI)
enddo

do k=1,num_levels
  do j=1,num_lat
     basic_temp(j,k) = z_factor1(k) + z_factor3(k)*(lat_factor1(j)*2*z_factor2(k) + lat_factor2(j))
  enddo
enddo

nv_surf = (1 - n0)*halfpi
do j=1,num_lat
  surf_geopotential(:,j) = U0*cos(nv_surf)**1.5 * (lat_factor1(j)*U0*cos(nv_surf)**1.5 + lat_factor2(j))
enddo

sin_latc = sin(PI*latc/180)
cos_latc = cos(PI*latc/180)
do i=1,num_lon
  lon_factor1(i) = cos(rad_lon(i)-PI*lonc/180)
enddo
do j=1,num_lat
  do i=1,num_lon
    xx = 10*acos(sin_latc*sin_lat(j) + cos_latc*cos_lat(j)*lon_factor1(i))
    perturbation(i,j) = Up*exp(-xx**2)
  enddo
enddo

id_phalf = diag_axis_init('phalf',.01*p_half,'hPa','z','approx half pressure level',direction=-1)
id_pfull = diag_axis_init('pfull',.01*p_full,'hPa','z','approx full pressure level',direction=-1,edges=id_phalf)
id_lon   = diag_axis_init('lon',deg_lon,'degrees E','x','longitude')
id_lat   = diag_axis_init('lat',deg_lat,'degrees N','y','latitude')
id_basic_flow   = register_static_field('jablonowski_2006', 'basic_flow',  (/id_lat,id_pfull/),'initial zonal wind', 'm/sec')
id_basic_temp   = register_static_field('jablonowski_2006', 'basic_temp',  (/id_lat,id_pfull/),'initial basic_temp', 'K')
id_pot_temp     = register_static_field('jablonowski_2006', 'pot_temp',    (/id_lat,id_pfull/),'initial potential temp', 'K')
id_perturbation = register_static_field('jablonowski_2006', 'perturbation',(/id_lon,id_lat/),  'initial perturbation','m/sec')
if(id_basic_flow > 0) used = send_data(id_basic_flow, basic_flow)
do k=1,num_levels
  do j=1,num_lat
    pot_temp(j,k) = basic_temp(j,k)*(P00/p_full(k))**KAPPA
  enddo
enddo
if(id_basic_temp > 0) used = send_data(id_basic_temp, basic_temp)
if(id_pot_temp   > 0) used = send_data(id_pot_temp,     pot_temp)
if(id_perturbation>0) used = send_data(id_perturbation, perturbation)

do k=1,num_levels
  do j=1,num_lat
    ug(:,j,k) = basic_flow(j,k)
    tg(:,j,k) = basic_temp(j,k)
  enddo
enddo
do j=1,num_lat
  do i=1,num_lon
    ug(i,j,:) = ug(i,j,:) + perturbation(i,j)
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

deallocate(deg_lon, deg_lat, rad_lon, rad_lat, sin_lat, cos_lat, coriolis)
deallocate(basic_flow, basic_temp, pot_temp, perturbation)
deallocate(lon_factor1, lat_factor1, lat_factor2, z_factor1, z_factor2, z_factor3)

return
end subroutine jablonowski_2006
!================================================================================

end module jablonowski_2006_mod
