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

module polvani_2007_mod

! Lorenzo Polvani, and J. Esler:
! Transport and mixing of chemical airmasses in idealized baroclinic life cycles
! J. Geophys. Res., 112, D23102, 10.1029/2007JD008555, 2007.

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use               mpp_mod, only: mpp_error, FATAL

use               fms_mod, only: mpp_pe, mpp_root_pe, stdlog, stdout, lowercase, &
                                 write_version_number, close_file, check_nml_error

use         constants_mod, only: PI,    &
                                 GRAV,  & ! GRAV=9.800, not 9.806 as specified by Polvani
                                 RDGAS, & ! RDGAS=287.04, not 287.00 as specified by Polvani
                                 OMEGA, & ! same as specified by Polvani
                                 RADIUS,& ! same as specified by Polvani
                                 KAPPA    ! same as specified by Polvani

use   vert_coordinate_mod, only: compute_vert_coord

use        transforms_mod, only: get_grid_boundaries, get_sin_lat, get_cos_lat, trans_grid_to_spherical, &
                                 trans_spherical_to_grid, get_wts_lat, get_lon_max, get_lat_max, get_deg_lon, &
                                 get_deg_lat, vor_div_from_uv_grid, uv_grid_from_vor_div, get_grid_domain, grid_domain

use  press_and_geopot_mod, only: press_and_geopot_init, pressure_variables, compute_pressures_and_heights

use      diag_manager_mod, only: diag_axis_init, register_static_field, send_data

use       tracer_type_mod, only: tracer_type

use     field_manager_mod, only: MODEL_ATMOS

use    tracer_manager_mod, only: get_number_tracers, get_tracer_names

implicit none
private

character(len=128), parameter :: version='$Id: polvani_2007.F90,v 19.0 2012/01/06 20:35:37 fms Exp $'
character(len=128), parameter :: tagname='$Name: siena_201211 $'

public :: polvani_2007, polvani_2007_tracer_init, get_polvani_2007_tracers

character(len=12), parameter :: mod_name='polvani_2007'

real,    parameter :: P00=1.e5        ! Used to compute potential temperature.
integer, parameter :: num_iter=10     ! Number of iterations in computation of surface pressure. (Overkill -- only 3 or 4 would do)
integer, parameter :: num_p2007_tr=11 ! Number of tracers associated with polvani_2007 initialization.
                                      ! Additional tracers may be defined elsewhere in the model. There is no conflict. 

integer :: id_LC2_Pstar, id_LC2_Tstar, id_LC2_U1star, id_LC2_U2star, id_LC2_Ustar, id_LC2_lapse00
logical :: used
integer :: i, j, k, lon_max, lat_max, is, ie, js, je, num_levels
integer :: id_LC1_u, id_LC1_pot_t, id_LC1_t, id_perturb, id_Tss, id_Uss, id_LC2_u, id_LC2_pot_t, id_LC2_t
real    :: du_dz, ln_slp, ztmp

real, allocatable, dimension(:,:) :: LC1_u, LC1_t, LC1_pot_t, LC2_u, LC2_t, LC2_pot_t
real, allocatable, dimension(:,:) :: dTdy, Uss, Tss, perturbation
real, allocatable, dimension(:)   :: sin_lat, tan_lat, coriolis, latb, rad_lat, fy1, fy2
real, allocatable, dimension(:)   :: LC2_Tstar, LC2_U1star, LC2_U2star, LC2_Ustar, LC2_Zstar
real, allocatable, dimension(:)   :: deg_lon, lonb, deg_lat, rad_lon, wts_lat, cos_lat
real, allocatable, dimension(:)   :: fz, dfz_dz, lon_factor, lat_factor, z
real, allocatable, dimension(:)   :: LC1_psurf, LC2_psurf
real, allocatable, dimension(:)   :: LC1_p_full, LC1_ln_p_full, LC1_p_half, LC1_ln_p_half
real, allocatable, dimension(:,:) :: LC2_p_full, LC2_ln_p_full, LC2_p_half, LC2_ln_p_half
real, allocatable, dimension(:)   :: trop_strat_boundary

real, allocatable, dimension(:,:,:) :: p2007_tr ! The tracers associated with polvani_2007 initialization.

character(len=15), parameter, dimension(num_p2007_tr) :: p2007_tr_name =   &
(/'polvani_2007_ov','polvani_2007_s1','polvani_2007_s2','polvani_2007_s3', &
  'polvani_2007_t1','polvani_2007_t2','polvani_2007_t3','polvani_2007_un', &
  'polvani_2007_bp','polvani_2007_b0','polvani_2007_bm'/)

character(len=31), parameter, dimension(num_p2007_tr) :: p2007_tr_longname =   &
(/'overworld tracer               ','high middleworld strat tracer  ','center middleworld strat tracer', &
  'low middleworld strat tracer   ','high middleworld trop tracer   ','center middleworld trop tracer ', &
  'low middleworld trop tracer    ','underworld tracer              ', 'high boundary layer tracer     ',&
  'canonical boundary layer tracer','low boundary layer tracer      '/)

! press_diag is the pressure levels used for diagnostic meta data.
! In the LC2 case they are not the actual model pressure levels.
real, allocatable, dimension(:) :: press_diag

character(len=8) :: type_of_init = 'LC1     ' ! 'LC1' or 'LC2'
real :: T_hat = 1.0      ! deg K ------------- Amplitude of temperature perturbation
integer :: m = 6         ! dimensionless ----- Zonal wavenumber of temperature perturbation
real :: theta_hat = 45.  ! degrees latitude -- Latitude center of temperature perturbation
real :: H  = 7.5e3       ! meters
real :: U0 = 45.         ! m/sec
real :: sigma_top = .02  ! dimensionless
real :: zt = 13.e3       ! meters
real :: lapse = -6.5e-3  ! deg K/meter
real :: T0 = 300.        ! deg K
real :: alpha = 10.      ! dimensionless
real :: Us = 45.         ! m/sec
real :: zs = 1.e4        ! meters
real :: theta_s = 35.    ! degrees latitude -- Center of the second term of the zonal wind profile for LC2 ("us" in equation (8) of the paper)
real :: delta_s = 20.    ! degrees latitude -- Width of same

namelist / polvani_2007_nml / type_of_init, T_hat, m, theta_hat, H, U0, &
                              sigma_top, zt, lapse, T0, alpha, Us, zs, theta_s, delta_s

contains

!=================================================================================================================================

subroutine polvani_2007(triang_trunc, vert_difference_option, pk, bk, &
                   vors, divs, ts, ln_ps, ug, vg, tg, psg, vorg, divg, surf_geopotential)

logical, intent(in) :: triang_trunc 
character(len=*), intent(in) :: vert_difference_option
real,    intent(out), dimension(:) :: pk, bk
complex, intent(out), dimension(:,:,:) :: vors, divs, ts
complex, intent(out), dimension(:,:  ) :: ln_ps
real,    intent(out), dimension(:,:,:) :: ug, vg, tg
real,    intent(out), dimension(:,:  ) :: psg
real,    intent(out), dimension(:,:,:) :: vorg, divg
real,    intent(out), dimension(:,:  ) :: surf_geopotential

real, dimension(size(ug,1), size(ug,2)) :: ln_psg
integer :: id_lon, id_lat, id_press, ntmp, id_iter, unit, ierr, io
real :: ln_sigma_top, riter(num_iter)

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=polvani_2007_nml, iostat=io)
    ierr = check_nml_error(io, 'polvani_2007_nml')
#else
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=polvani_2007_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'polvani_2007_nml')
    enddo
20  call close_file (unit)
#endif
call write_version_number(version, tagname)
if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=polvani_2007_nml)

call get_lon_max(lon_max)
call get_lat_max(lat_max)
num_levels = size(ug,3)
ntmp = nint(.46875*lat_max)
if(num_levels /= ntmp) then
  call mpp_error(FATAL,'Number of vertical levels must be ',ntmp, &
                ' when using polvani_2007 initialization with ',lat_max,' latitude rows')
endif
pk = 0.0
bk(1) = sigma_top
ln_sigma_top = alog(bk(1))
do k=2,num_levels
  bk(k) = exp((1.0 - (k-1)/float(num_levels))*ln_sigma_top)
enddo
bk(num_levels+1) = 1.0

surf_geopotential = 0.0

allocate(rad_lon(lon_max), deg_lon(lon_max), lon_factor(lon_max))
allocate(rad_lat(lat_max), deg_lat(lat_max), lat_factor(lat_max), trop_strat_boundary(lat_max))
allocate(wts_lat(lat_max), sin_lat(lat_max), cos_lat(lat_max), tan_lat(lat_max))
allocate(coriolis(lat_max), latb(lat_max+1), lonb(lon_max+1), z(num_levels+1))
allocate(LC1_u(lat_max,num_levels+1), dTdy(lat_max,num_levels+1), Uss(lat_max,num_levels+1))
allocate(LC1_pot_t(lat_max,num_levels+1), LC1_t(lat_max,num_levels+1), Tss(lat_max,num_levels+1))
allocate(LC2_u(lat_max,num_levels+1), LC2_pot_t(lat_max,num_levels+1), LC2_t(lat_max,num_levels+1))
allocate(perturbation(lon_max,lat_max), press_diag(num_levels+1))
allocate(fz(num_levels+1), dfz_dz(num_levels+1), fy1(lat_max), fy2(lat_max))
allocate(LC2_Tstar(lat_max), LC2_U1star(lat_max), LC2_U2star(lat_max), LC2_Ustar(lat_max), LC2_Zstar(lat_max))
allocate(LC1_psurf(lat_max), LC2_psurf(lat_max))
allocate(LC1_p_full(num_levels),   LC1_ln_p_full(num_levels))
allocate(LC1_p_half(num_levels+1), LC1_ln_p_half(num_levels+1))
allocate(LC2_p_full(lat_max,num_levels),   LC2_ln_p_full(lat_max,num_levels))
allocate(LC2_p_half(lat_max,num_levels+1), LC2_ln_p_half(lat_max,num_levels+1))
allocate(p2007_tr(lat_max,num_levels,num_p2007_tr))

call press_and_geopot_init(pk, bk, .false., vert_difference_option)
call pressure_variables(LC1_p_half, LC1_ln_p_half, LC1_p_full, LC1_ln_p_full, P00)
call get_deg_lon(deg_lon)
rad_lon = PI*deg_lon/180.
call get_deg_lat(deg_lat)
rad_lat = PI*deg_lat/180.
trop_strat_boundary = 1500*(55-deg_lat) ! units=meters (only used in the narrow latitude range of about 45 to 52 degrees north)
call get_wts_lat(wts_lat)
call get_sin_lat(sin_lat)
call get_cos_lat(cos_lat)
tan_lat = sin_lat/cos_lat
call get_grid_boundaries(lonb, latb, global=.true.) ! units = radians
coriolis = 2*OMEGA*sin_lat
ln_slp = log(P00)

riter = (/(i,i=1,num_iter)/)
press_diag(1:num_levels) = .01*LC1_p_full
press_diag(num_levels+1) = .01*LC1_p_half(num_levels+1)
id_press = diag_axis_init('p_full',press_diag,'hPa','z','pressure level',direction=-1)
id_lon   = diag_axis_init('lon',deg_lon,'degrees E','x','longitude')
id_lat   = diag_axis_init('lat',deg_lat,'degrees N','y','latitude')
id_iter  = diag_axis_init('iter',riter,'hPa','z','iteration')
id_LC1_u     = register_static_field(mod_name,'LC1_zonal_wind', (/id_lat,id_press/),'initial LC1 zonal wind', 'm/sec')
id_LC1_t     = register_static_field(mod_name,'LC1_temperature',(/id_lat,id_press/),'initial LC1 temperature', 'K')
id_LC1_pot_t = register_static_field(mod_name,'LC1_pot_temp',   (/id_lat,id_press/),'initial LC1 potential temperature','K')
id_Tss       = register_static_field(mod_name,'Tss',            (/id_lat,id_press/),'Tss temperature', 'K')
id_Uss       = register_static_field(mod_name,'Uss',            (/id_lat,id_press/),'Uss zonal wind', 'm/sec')
id_LC2_u     = register_static_field(mod_name,'LC2_zonal_wind', (/id_lat,id_press/),'initial LC2 zonal wind', 'm/sec')
id_LC2_t     = register_static_field(mod_name,'LC2_temperature',(/id_lat,id_press/),'initial LC2 temperature', 'K')
id_LC2_pot_t = register_static_field(mod_name,'LC2_pot_temp',   (/id_lat,id_press/),'initial LC2 potential temperature','K')
id_perturb   = register_static_field(mod_name,'perturbation',   (/id_lon,id_lat/),  'initial temperature perturbation','K')
id_LC2_Pstar = register_static_field(mod_name,'LC2_Pstar',      (/id_lat/),         'initial surface pressure','hPa')
id_LC2_Tstar = register_static_field(mod_name,'LC2_Tstar',      (/id_lat/),         'initial surface temperature','K')
id_LC2_U1star= register_static_field(mod_name,'LC2_U1star',     (/id_lat/),         'initial surface LC1 wind', 'm/sec')
id_LC2_U2star= register_static_field(mod_name,'LC2_U2star',     (/id_lat/),         'initial surface LC2 wind', 'm/sec')
id_LC2_Ustar = register_static_field(mod_name,'LC2_Ustar',      (/id_lat/),         'initial surface wind',     'm/sec')
id_LC2_lapse00=register_static_field(mod_name,'LC2_lapse00',    (/id_lat/),         'initial surface lapse rate','K/m')

call compute_LC1
call compute_LC2
call compute_perturbation

do j=1,lat_max
  do k=1,num_levels+1
    if(k == num_levels+1) then
      LC1_pot_t(j,k) = LC1_t(j,k)
    else
      LC1_pot_t(j,k) = LC1_t(j,k)*(P00/LC1_p_full(k))**KAPPA
    endif
  enddo
enddo
if(id_LC1_u     > 0) used = send_data(id_LC1_u, LC1_u)
if(id_LC1_t     > 0) used = send_data(id_LC1_t, LC1_t)
if(id_LC1_pot_t > 0) used = send_data(id_LC1_pot_t, LC1_pot_t)
if(id_perturb   > 0) used = send_data(id_perturb,  perturbation)

call get_grid_domain(is, ie, js, je)

if(trim(type_of_init) == 'LC1') then
  do k=1,num_levels
    do i=1,size(ug,1)
      ug(i,:,k) = LC1_u(js:je,k)
      tg(i,:,k) = LC1_t(js:je,k)
    enddo
  enddo
  do j=1,size(ug,2)
    psg(:,j) = LC1_psurf(js+j-1)
  enddo
  call compute_p2007_tr(LC1_t,LC1_psurf)
else if(trim(type_of_init) == 'LC2') then
  do k=1,num_levels
    do i=1,size(ug,1)
      ug(i,:,k) = LC2_u(js:je,k)
      tg(i,:,k) = LC2_t(js:je,k)
    enddo
  enddo
  do j=1,size(ug,2)
    psg(:,j) = LC2_psurf(js+j-1)
  enddo
  call compute_p2007_tr(LC2_t,LC2_psurf)
else
  call mpp_error(FATAL,trim(type_of_init)//' is not a valid type_of_init for polvani_2007')
endif

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
ln_psg = alog(psg)
call trans_grid_to_spherical(ln_psg, ln_ps)
call trans_spherical_to_grid(ln_ps,  ln_psg)
psg = exp(ln_psg)

return
end subroutine polvani_2007
!=================================================================================================================================
subroutine compute_LC1

real :: Tr

do k=1,num_levels+1
  if(k == num_levels+1) then
    z(k) = 0.0
  else
    z(k) = H*(ln_slp - LC1_ln_p_full(k))
  endif
  ztmp = z(k)/zt
  fz(k) = ztmp*exp(-0.5*(ztmp**2 - 1))
  dfz_dz(k) = ((1 - ztmp**2)/zt)*exp(-0.5*(ztmp**2 - 1))
enddo

do j=1,lat_max
  if(sin_lat(j) > 0.0) then
    fy1(j) = (sin(PI*sin_lat(j)**2))**3
  else
    fy1(j) = 0.0
  endif
enddo

do k=1,num_levels+1
  do j=1,lat_max
    LC1_u(j,k) = U0*fy1(j)*fz(k)
    du_dz   = U0*fy1(j)*dfz_dz(k)
    dTdy(j,k) = -(H/RDGAS) * (RADIUS*coriolis(j) + 2*LC1_u(j,k)*tan_lat(j)) * du_dz
  enddo
enddo

! The k loop below computes the integral of equation A4 in the paper.
! The inner loop over latitude does the actual integration.
do k=1,num_levels+1
  do j=1,lat_max
    if(j == 1) then
      LC1_t(j,k) = dTdy(j,k)*(rad_lat(j)-latb(j))
    else
      LC1_t(j,k) = LC1_t(j-1,k) + dTdy(j-1,k)*(latb(j)-rad_lat(j-1)) + dTdy(j,k)*(rad_lat(j)-latb(j))
    endif
  enddo
enddo
do k=1,num_levels+1
  if(k == num_levels+1) then
    Tr = T0
  else
    Tr = T0 + lapse/(zt**-alpha + z(k)**-alpha)**(1./alpha)
  endif
  LC1_t(:,k) = Tr + LC1_t(:,k)
enddo

LC1_psurf = P00

end subroutine compute_LC1
!=================================================================================================================================
subroutine compute_perturbation
    
do i=1,lon_max
  lon_factor(i) = cos(m*rad_lon(i))
enddo
ztmp = PI*theta_hat/180
do j=1,lat_max
  lat_factor(j) = 1.0/cosh(m*(rad_lat(j)-ztmp))**2
enddo
do j=1,lat_max 
  do i=1,lon_max
    perturbation(i,j) = T_hat*lon_factor(i)*lat_factor(j)
  enddo
enddo

end subroutine compute_perturbation
!=================================================================================================================================
subroutine compute_LC2

do k=1,num_levels+1
  fz(k) = exp(-z(k)/zs)
  dfz_dz(k) = -fz(k)/zs
enddo

do j=1,lat_max
  fy2(j) = (sin(2*rad_lat(j))**2)*((deg_lat(j)-theta_s)/delta_s)*exp(-((deg_lat(j)-theta_s)/delta_s)**2)
enddo

do k=1,num_levels+1
  do j=1,lat_max
    Uss(j,k) = -Us*fy2(j)*fz(k)
    du_dz    = -Us*fy2(j)*dfz_dz(k)
    dTdy(j,k) = -(H/RDGAS) * (RADIUS*coriolis(j) + 2*Uss(j,k)*tan_lat(j)) * du_dz
  enddo
enddo

! The k loop below computes the integral of equation A4 in the paper.
! The inner loop over latitude does the actual integration.
do k=1,num_levels+1
  do j=1,lat_max
    if(j == 1) then
      Tss(j,k) = dTdy(j,k)*(rad_lat(j)-latb(j))
    else
      Tss(j,k) = Tss(j-1,k) + dTdy(j-1,k)*(latb(j)-rad_lat(j-1)) + dTdy(j,k)*(rad_lat(j)-latb(j))
    endif
  enddo
enddo

if(id_Uss > 0) used = send_data(id_Uss, Uss)
if(id_Tss > 0) used = send_data(id_Tss, Tss)

LC2_u = LC1_u + Uss
LC2_t = LC1_t + Tss

! Surface pressure for LC2 is not constant.
! The paper does not desribe the details of how it is computed, so I concocted my own scheme to do it.
! I have isolated this calculation by putting it in a separate subroutine.
call compute_surf_press

call pressure_variables(LC2_p_half, LC2_ln_p_half, LC2_p_full, LC2_ln_p_full, LC2_psurf)

do j=1,lat_max
  do k=1,num_levels+1
    if(k == num_levels+1) then
      LC2_pot_t(j,k) = LC2_t(j,k)
    else
      LC2_pot_t(j,k) = LC2_t(j,k)*(P00/LC2_p_full(j,k))**KAPPA
    endif
  enddo
enddo
if(id_LC2_u     > 0) used = send_data(id_LC2_u, LC2_u)
if(id_LC2_t     > 0) used = send_data(id_LC2_t, LC2_t)
if(id_LC2_pot_t > 0) used = send_data(id_LC2_pot_t, LC2_pot_t)

return
end subroutine compute_LC2
!=================================================================================================================================
subroutine compute_surf_press

real, allocatable, dimension(:,:) :: Zstar     ! Zstar = H*alog(P00/psurf)
real, allocatable, dimension(:) :: af          ! af = RADIUS*coriolis
real, allocatable, dimension(:) :: dzstar_dy   ! derivative of Zstar with respect to latitude
real, allocatable, dimension(:) :: LC2_lapse00 ! lapse rate at pressure level P00
real, allocatable, dimension(:) :: dlapse00_dy ! derivative of LC2_lapse00 with respect to latitude
real :: rms, c1, c2, e, sqrt_e
integer :: j, itr

allocate(Zstar(lat_max,0:num_iter), af(lat_max), dzstar_dy(lat_max), dlapse00_dy(lat_max), LC2_lapse00(lat_max))
af = RADIUS*coriolis
e = exp(1.0)
sqrt_e = sqrt(e)
c1 = 2*e*(U0/zt)**2
c2 = Us/zs**2

do j=1,lat_max
  if(sin_lat(j) > 0.0) then
    dlapse00_dy(j) = c1*tan_lat(j)*fy1(j)**2 - (af(j) - 2*Us*fy2(j)*tan_lat(j))*c2*fy2(j)
  else
    dlapse00_dy(j) = 0.0
  endif
enddo
dlapse00_dy = -(H/RDGAS)*dlapse00_dy
do j=1,lat_max
  if(j == 1) then
    LC2_lapse00(j) = dlapse00_dy(j)*(rad_lat(j)-latb(j))
  else
    LC2_lapse00(j) = LC2_lapse00(j-1) + dlapse00_dy(j-1)*(latb(j)-rad_lat(j-1)) + dlapse00_dy(j)*(rad_lat(j)-latb(j))
  endif
enddo
LC2_lapse00 = LC2_lapse00 + lapse

! Zstar is computed by an iterative method. The zonal momentum equation for the
! surface winds is integrated from pole to pole each interation. Temperature and
! zonal wind are expressed as a Taylor series expansion about z=0, i.e., Pressure = P00.
Zstar(:,0) = 0.0
do itr=1,num_iter
  do j=1,lat_max
    LC2_Tstar(j)  = LC2_t(j,num_levels+1) + LC2_lapse00(j)*Zstar(j,itr-1)
    LC2_U1star(j) = U0*sqrt_e*fy1(j)*Zstar(j,itr-1)/zt
    LC2_U2star(j) = (Zstar(j,itr-1)/zs - 1)*Us*fy2(j)
    LC2_Ustar(j)  = LC2_U1star(j) + LC2_U2star(j)
    if(sin_lat(j) > 0.0) then
      dzstar_dy(j) = H*LC2_Ustar(j)*(af(j) + LC2_Ustar(j)*tan_lat(j)) / (RDGAS*LC2_Tstar(j))
    else
      dzstar_dy(j) = 0.0
    endif
  enddo
  do j=1,lat_max
    if(j == 1) then
      Zstar(j,itr) = dzstar_dy(j)*(rad_lat(j)-latb(j))
    else
      Zstar(j,itr) = Zstar(j-1,itr) + dzstar_dy(j-1)*(latb(j)-rad_lat(j-1)) + dzstar_dy(j)*(rad_lat(j)-latb(j))
    endif
  enddo
  rms = 0.0
  do j=1,lat_max
    rms = rms + (Zstar(j,itr) - Zstar(j,itr-1))**2
  enddo
  rms = sqrt(rms/lat_max)
  write(stdout(),'("Convergence of surface pressure for polvani_2007 LC2: iteration=",i2,"  rms=",1pe16.8)') itr,rms
enddo

LC2_Zstar = Zstar(:,num_iter)

do j=1,lat_max
  LC2_psurf(j) = P00*exp(-Zstar(j,num_iter)/H)
enddo

if(id_LC2_Pstar  > 0) used = send_data(id_LC2_Pstar,  LC2_psurf)
if(id_LC2_Tstar  > 0) used = send_data(id_LC2_Tstar,  LC2_Tstar )
if(id_LC2_U1star > 0) used = send_data(id_LC2_U1star, LC2_U1star)
if(id_LC2_U2star > 0) used = send_data(id_LC2_U2star, LC2_U2star)
if(id_LC2_Ustar  > 0) used = send_data(id_LC2_Ustar,  LC2_Ustar )
if(id_LC2_lapse00> 0) used = send_data(id_LC2_lapse00,LC2_lapse00)
deallocate(Zstar, af, dzstar_dy, dlapse00_dy, LC2_lapse00)

end subroutine compute_surf_press
!=================================================================================================================================
subroutine polvani_2007_tracer_init(tg, psg, surface_geopotential, tracers)

real, dimension(:,:,:),     intent(in)  :: tg
real, dimension(:,:),       intent(in)  :: psg, surface_geopotential
real, dimension(:,:,:,:,:), intent(out) :: tracers

real,    dimension(is:ie,js:je,num_levels)   :: z_full, p_full, pot_t
real,    dimension(is:ie,js:je,num_levels+1) :: z_half, p_half
integer, dimension(is:ie,js:je,num_levels)   :: world ! overworld=0, high middleworld=1, center middleworld=2, low middleworld=3, underworld=4
logical, dimension(is:ie,js:je,num_levels)   :: strat ! .true. if above tropopause, .false. if below. 

integer :: ntr, num_tracers, ntr_o, ntr_s1, ntr_s2, ntr_s3, ntr_t1, ntr_t2, ntr_t3, ntr_u
character(len=128) :: tname, longname, units

call compute_pressures_and_heights(tg, psg, surface_geopotential, z_full, z_half, p_full, p_half)
do k=1,num_levels
  do j=js,je
    do i=is,je
      pot_t(i,j,k) = tg(i,j,k)*(P00/p_full(i,j,k))**KAPPA
      if(pot_t(i,j,k) >= 380.) then
        world(i,j,k) = 0
      else if(pot_t(i,j,k) >= 350. .and. pot_t(i,j,k) < 380.) then
        world(i,j,k) = 1
      else if(pot_t(i,j,k) > 320. .and. pot_t(i,j,k) < 350.) then
        world(i,j,k) = 2
      else if(pot_t(i,j,k) > 290. .and. pot_t(i,j,k) <= 320.) then
        world(i,j,k) = 3
      else
        world(i,j,k) = 4
      endif
      if(z_full(i,j,k) > trop_strat_boundary(j)) then
        strat(i,j,k) = .true.
      else
        strat(i,j,k) = .false.
      endif
    enddo
  enddo
enddo

ntr_o=0; ntr_s1=0; ntr_s2=0; ntr_s3=0; ntr_t1=0; ntr_t2=0; ntr_t3=0; ntr_u=0
call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
do ntr=1,num_tracers
  call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
  if(trim(lowercase(tname)) == 'polvani_2007_o') then
    ntr_o = ntr
    tracers(:,:,:,:,ntr) = 0.0
  else if(trim(lowercase(tname)) == 'polvani_2007_s1') then
    ntr_s1 = ntr
    tracers(:,:,:,:,ntr) = 0.0
  else if(trim(lowercase(tname)) == 'polvani_2007_s2') then
    ntr_s2 = ntr
    tracers(:,:,:,:,ntr) = 0.0
  else if(trim(lowercase(tname)) == 'polvani_2007_s3') then
    ntr_s3 = ntr
    tracers(:,:,:,:,ntr) = 0.0
  else if(trim(lowercase(tname)) == 'polvani_2007_t1') then
    ntr_t1 = ntr
    tracers(:,:,:,:,ntr) = 0.0
  else if(trim(lowercase(tname)) == 'polvani_2007_t2') then
    ntr_t2 = ntr
    tracers(:,:,:,:,ntr) = 0.0
  else if(trim(lowercase(tname)) == 'polvani_2007_t3') then
    ntr_t3 = ntr
    tracers(:,:,:,:,ntr) = 0.0
  else if(trim(lowercase(tname)) == 'polvani_2007_u') then
    ntr_u = ntr
    tracers(:,:,:,:,ntr) = 0.0
  endif
enddo

do k=1,num_levels
  do j=js,je
    do i=is,je
      if(world(i,j,k) == 0 .and. ntr_o > 0) then
        tracers(i,j,k,:,ntr_o) = 1.0
      else if(world(i,j,k) == 1) then
        if(strat(i,j,k)) then
          if(ntr_s1 > 0) tracers(i,j,k,:,ntr_s1) = 1.0
        else
          if(ntr_t1 > 0) tracers(i,j,k,:,ntr_t1) = 1.0
        endif
      else if(world(i,j,k) == 2) then
        if(strat(i,j,k)) then
          if(ntr_s2 > 0) tracers(i,j,k,:,ntr_s2) = 1.0
        else
          if(ntr_t2 > 0) tracers(i,j,k,:,ntr_t2) = 1.0
        endif
      else if(world(i,j,k) == 3) then
        if(strat(i,j,k)) then
          if(ntr_s3 > 0) tracers(i,j,k,:,ntr_s3) = 1.0
        else
          if(ntr_t3 > 0) tracers(i,j,k,:,ntr_t3) = 1.0
        endif
      else if(world(i,j,k) == 4 .and. ntr_u > 0) then
        tracers(i,j,k,:,ntr_u) = 1.0
      endif
    enddo
  enddo
enddo

end subroutine polvani_2007_tracer_init
!=================================================================================================================================

subroutine compute_p2007_tr(temp, surfp)
real, dimension(lat_max,num_levels), intent(in) :: temp
real, dimension(lat_max),            intent(in) :: surfp

real, dimension(lat_max)              :: surf_geopotential
real, dimension(lat_max,num_levels  ) :: pot_t, z_full, p_full
real, dimension(lat_max,num_levels+1) :: z_half, p_half
integer :: world(lat_max,num_levels) ! overworld=0, high middleworld=1, center middleworld=2, low middleworld=3, underworld=4
logical :: strat(lat_max,num_levels) ! .true. if above tropopause, .false. if below
character(len=15) :: tname
integer :: ntr, id_tr, id_latitude, id_pfull, axes(2)

surf_geopotential = 0.0
call compute_pressures_and_heights(temp, surfp, surf_geopotential, z_full, z_half, p_full, p_half)
do k=1,num_levels
  do j=1,lat_max
    pot_t(j,k) = temp(j,k)*(P00/p_full(j,k))**KAPPA
    if(pot_t(j,k) >= 380.) then
      world(j,k) = 0
    else if(pot_t(j,k) >= 350. .and. pot_t(j,k) < 380.) then
      world(j,k) = 1
    else if(pot_t(j,k) > 320. .and. pot_t(j,k) < 350.) then
      world(j,k) = 2
    else if(pot_t(j,k) > 290. .and. pot_t(j,k) <= 320.) then
      world(j,k) = 3
    else
      world(j,k) = 4
    endif
    if(z_full(j,k) > trop_strat_boundary(j)) then
      strat(j,k) = .true.
    else
      strat(j,k) = .false.
    endif
  enddo
enddo

p2007_tr = 0.0
do ntr=1,num_p2007_tr
  tname = p2007_tr_name(ntr)
  do k=1,num_levels
    do j=1,lat_max

      if(world(j,k) == 0 .and. tname == 'polvani_2007_ov') then
        p2007_tr(j,k,ntr) = 1.0
      else if(world(j,k) == 1) then
        if(strat(j,k)) then
          if(tname == 'polvani_2007_s1') p2007_tr(j,k,ntr) = 1.0
        else
          if(tname == 'polvani_2007_t1') p2007_tr(j,k,ntr) = 1.0
        endif
      else if(world(j,k) == 2) then
        if(strat(j,k)) then
          if(tname == 'polvani_2007_s2') p2007_tr(j,k,ntr) = 1.0
        else
          if(tname == 'polvani_2007_t2') p2007_tr(j,k,ntr) = 1.0
        endif
      else if(world(j,k) == 3) then
        if(strat(j,k)) then
          if(tname == 'polvani_2007_s3') p2007_tr(j,k,ntr) = 1.0
        else
          if(tname == 'polvani_2007_t3') p2007_tr(j,k,ntr) = 1.0
        endif
      else if(world(j,k) == 4 .and. tname == 'polvani_2007_un') then
        p2007_tr(j,k,ntr) = 1.0
      endif

      if(z_full(j,k) < 2500. .and. tname == 'polvani_2007_bp') then
        p2007_tr(j,k,ntr) = 1.0
      endif
      if(z_full(j,k) < 1500. .and. tname == 'polvani_2007_b0') then
        p2007_tr(j,k,ntr) = 1.0
      endif
      if(z_full(j,k) < 500.  .and. tname == 'polvani_2007_bm') then
        p2007_tr(j,k,ntr) = 1.0
      endif

    enddo
  enddo
enddo

id_latitude  = diag_axis_init('latitude' , deg_lat, 'degrees_N', 'y', 'latitude' , set_name=mod_name)
id_pfull     = diag_axis_init('pfull', .01*LC1_p_full, 'hPa'    , 'z', 'pressure' , set_name=mod_name)
axes = (/id_latitude,id_pfull/)
do ntr=1,num_p2007_tr
  id_tr = register_static_field(mod_name, 'initial_'//p2007_tr_name(ntr), axes, p2007_tr_longname(ntr), 'none')
  used = send_data(id_tr, p2007_tr(:,:,ntr))
enddo

end subroutine compute_p2007_tr
!=================================================================================================================================
subroutine get_polvani_2007_tracers(tracers)
real, dimension(:,:,:,:), intent(out) :: tracers
integer :: num_tracers, ntr1, ntr2
character(len=128) :: tname, longname, units

call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
do ntr1=1,num_tracers
  call get_tracer_names(MODEL_ATMOS, ntr1, tname, longname, units)
  do ntr2=1,num_p2007_tr
    if(trim(lowercase(tname)) == p2007_tr_name(ntr2)) then
      do k=1,num_levels
        do i=1,size(tracers,1)
          tracers(i,:,k,ntr1) = p2007_tr(js:je,k,ntr2)
        enddo
      enddo
    endif
  enddo
enddo

end subroutine get_polvani_2007_tracers
!=================================================================================================================================
end module polvani_2007_mod
