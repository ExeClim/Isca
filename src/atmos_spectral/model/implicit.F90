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

module implicit_mod
!||zn, pjk 5/13/99
!following vb's changes to old implicit_mod.
!implicit_init works with global fields, everything else on domain

use              fms_mod, only: error_mesg, FATAL, write_version_number

use        constants_mod, only: rdgas, radius, cp_air

use press_and_geopot_mod, only: pressure_variables

use    matrix_invert_mod, only: invert

use       transforms_mod, only: get_spec_domain

implicit none

private

public :: implicit_init, implicit_correction, implicit_end

interface linear_geopotential
   module procedure   linear_geopotential_1d, linear_geopotential_3d
end interface

interface linear_tp_tendency
   module procedure   linear_tp_tendency_1d, linear_tp_tendency_3d
end interface

logical :: module_is_initialized = .false.

character(len=128), parameter :: version = '$Id: implicit.F90,v 13.0 2006/03/28 21:17:54 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

real,    allocatable, dimension(:)   :: ref_temperature_implicit
real,    allocatable, dimension(:)   :: ref_ln_p_half, ref_ln_p_full, del_ln_p_half, del_ln_p_full
real,    allocatable, dimension(:,:) :: eigen
integer, allocatable, dimension(:,:) :: wavenumber

real :: alpha, ref_surf_p_implicit

real, allocatable, dimension(:,:)   :: div_mat
real, allocatable, dimension(:)     :: h
real, allocatable, dimension(:,:,:) :: wave_matrix
real, allocatable, dimension(:)     :: pk, bk, dpk, dbk

real :: dt = 0.0
real :: xi

integer :: num_levels, num_total_wavenumbers
character(len=64) :: vert_difference_option
integer :: ms, me, ns, ne

contains

!------------------------------------------------------------------------

subroutine implicit_init(pk_in, bk_in, ref_temperature_implicit_in, ref_surf_p_implicit_in, &
                         num_total_wavenumbers_in, eigen_in, wavenumber_in, alpha_in,       &
                         vert_difference_option_in)

real,    intent(in), dimension(:)     :: pk_in, bk_in, ref_temperature_implicit_in
real,    intent(in), dimension(0:,0:) :: eigen_in
integer, intent(in), dimension(0:,0:) :: wavenumber_in

real,    intent(in) :: ref_surf_p_implicit_in, alpha_in
integer, intent(in) :: num_total_wavenumbers_in
character(len=*), intent(in) :: vert_difference_option_in

real, dimension(size(ref_temperature_implicit_in,1)  ) :: p_full_work_1, p_full_work_2, ln_p_full_work_1, ln_p_full_work_2
real, dimension(size(ref_temperature_implicit_in,1)+1) :: p_half_work_1, p_half_work_2, ln_p_half_work_1, ln_p_half_work_2

real :: surface_p_1, eps=1.e-5
integer :: k

if(module_is_initialized) return

call write_version_number(version, tagname)

call get_spec_domain(ms, me, ns, ne)

num_levels = size(ref_temperature_implicit_in,1)
num_total_wavenumbers = num_total_wavenumbers_in

allocate(ref_temperature_implicit(num_levels))
allocate(ref_ln_p_half(num_levels+1))
allocate(ref_ln_p_full(num_levels))
allocate(del_ln_p_half(num_levels+1))
allocate(del_ln_p_full(num_levels))
allocate(div_mat(num_levels,num_levels))
allocate(wave_matrix(num_levels,num_levels,0:num_total_wavenumbers))
allocate(h(num_levels))
allocate(eigen(0:size(eigen_in,1)-1,0:size(eigen_in,2)-1))
allocate(wavenumber(0:size(eigen_in,1)-1,0:size(eigen_in,2)-1))

ref_temperature_implicit = ref_temperature_implicit_in
ref_surf_p_implicit = ref_surf_p_implicit_in
eigen         = eigen_in
wavenumber    = wavenumber_in
alpha         = alpha_in
vert_difference_option = vert_difference_option_in

allocate( pk(size(pk_in(:))),  bk(size(bk_in,1)))
allocate(dpk(num_levels  ), dbk(num_levels  ))

pk = pk_in
bk = bk_in

do k = 1,num_levels
  dpk(k) = pk(k+1) - pk(k)
  dbk(k) = bk(k+1) - bk(k)
end do

call pressure_variables(p_half_work_1, ref_ln_p_half, p_full_work_1, ref_ln_p_full, ref_surf_p_implicit)

!del_ln_p_half = dln_p_half_dps(ref_surf_p_implicit)
!del_ln_p_full = dln_p_full_dps(ref_surf_p_implicit)

!  functions dln_p_half_dps and del_ln_p_full are used nowhere else and
!  they perform very simple calculations, so it only makes the code more
!  confusing to do these calculations in functions.

do k=2,num_levels+1
  del_ln_p_half(k) = bk(k)/(pk(k)+bk(k)*ref_surf_p_implicit)
enddo

if(pk(1).eq.0.0) then
  del_ln_p_half(1) = 1.0/ref_surf_p_implicit
else
  del_ln_p_half(1) = bk(1)/(pk(1)+bk(1)*ref_surf_p_implicit)
endif

surface_p_1 = ref_surf_p_implicit*(1.0 - 0.5*eps)
call pressure_variables(p_half_work_1, ln_p_half_work_1, p_full_work_1, ln_p_full_work_1, surface_p_1)

surface_p_1 = ref_surf_p_implicit*(1.0 + 0.5*eps)
call pressure_variables(p_half_work_2, ln_p_half_work_2, p_full_work_2, ln_p_full_work_2, surface_p_1)

del_ln_p_full = (ln_p_full_work_2 - ln_p_full_work_1)/(eps*ref_surf_p_implicit)

call build_matrix

module_is_initialized = .true.

return
end subroutine implicit_init

!------------------------------------------------------------------------

subroutine build_matrix

real, dimension(num_levels) :: input, zero, dt_t
real, dimension(num_levels+1) :: zero_1
real :: dt_p

real, dimension(num_levels, num_levels) :: gamma, tau
real, dimension(num_levels) :: nu, h1, h2

integer :: k, kk, kkk

zero  = 0.0
zero_1 = 0.0

tau = 0.0
nu = 0.0
do k=1,num_levels
  input    = 0.0
  input(k) = 1.0
  
  call linear_tp_tendency(input, ref_temperature_implicit, ref_surf_p_implicit, &
                          ref_ln_p_half, ref_ln_p_full, dt_p, dt_t)
  nu(k)    = - dt_p
  tau(:,k) = - dt_t

  gamma(:,k) = linear_geopotential(input,zero_1,zero, &
                    ref_temperature_implicit, ref_ln_p_half, ref_ln_p_full)
end do

h1 = pres_grad_funct(ref_temperature_implicit, ref_ln_p_half, ref_ln_p_full, ref_surf_p_implicit)

h2 = linear_geopotential(zero,del_ln_p_half,del_ln_p_full, &
                    ref_temperature_implicit, ref_ln_p_half, ref_ln_p_full)

h = h1 + h2

do k = 1,num_levels
  do kk = 1,num_levels
    div_mat(k,kk)= h(k)*nu(kk)
    do kkk = 1,num_levels
      div_mat(k,kk)=div_mat(k,kk)+gamma(k,kkk)*tau(kkk,kk)   
    enddo
  enddo
enddo

return
end subroutine build_matrix

!-------------------------------------------------------------------------

subroutine build_wave_matrices

real :: factor, det
integer :: k, L

wave_matrix = 0.0
do L = 0,num_total_wavenumbers
  factor = xi*xi*L*(L+1)/radius**2
  do k = 1,num_levels
    wave_matrix(k,k,L) = 1.0
  end do
  wave_matrix(:,:,L) = wave_matrix(:,:,L) + factor*div_mat(:,:)
  call invert(wave_matrix(:,:,L), det)
end do

return
end subroutine build_wave_matrices

!---------------------------------------------------------------------

subroutine implicit_correction(dt_divs, dt_ts, dt_ln_ps, divs, ts, ln_ps, dt_in, previous, current)
!in parallel, wavenumber (0,0) actually means (ms,ns)
!it is necessary to use the (0:,0:) dimensioning here because of the m-n loop below
complex, intent(inout), dimension(0:,0:,:) :: dt_divs, dt_ts
complex, intent(inout), dimension(0:,0:) :: dt_ln_ps
complex, intent(in), dimension(0:,0:,:,:) :: divs, ts
complex, intent(in), dimension(0:,0:,:) :: ln_ps
real, intent(in) :: dt_in
integer,intent(in) :: previous, current

complex, dimension(0:size(divs,1)-1,0:size(divs,2)-1,size(divs,3)) :: dt_ts_temp
complex, dimension(num_levels) ::  work
complex, dimension(0:size(divs,1)-1,0:size(divs,2)-1) :: dt_ps_temp
integer :: m,n,L

if(.not.module_is_initialized) then
  call error_mesg('implicit_correction','failed to initialize implicit', FATAL)
end if

if(dt_in.ne.dt) then 
   dt = dt_in
   xi = dt*alpha
   call build_wave_matrices
end if

call adjust_dt_divs(dt_divs, dt_ts, dt_ln_ps, divs, ts, ln_ps, previous, current)

do n = 0,size(divs,2)-1
  do m = 0,size(divs,1)-1
!||zn, pjk 5/13/99, following vb's dmsm/implicit.f90
    L = wavenumber(m+ms,n+ns)
    if(L.le.num_total_wavenumbers) then
      work = matmul(wave_matrix(:,:,L),dt_divs(m,n,:))
      dt_divs(m,n,:) = work
    end if
  end do
end do

call linear_tp_tendency(dt_divs, ref_temperature_implicit, ref_surf_p_implicit, &
                        ref_ln_p_half, ref_ln_p_full, dt_ps_temp, dt_ts_temp)

dt_ts    = dt_ts    + xi*dt_ts_temp
dt_ln_ps = dt_ln_ps + xi*dt_ps_temp/ref_surf_p_implicit

return
end subroutine implicit_correction

!---------------------------------------------------------------------
subroutine adjust_dt_divs(dt_divs, dt_ts, dt_ln_ps, divs, ts, ln_ps, previous, current)
  !||zn, pjk 5/13/99, following vb's dmsm/implicit.f90. (0:,0:) -> (:,:)
complex, intent(inout), dimension(:,:,:) :: dt_divs
complex, intent(inout), dimension(:,:,:) :: dt_ts
complex, intent(inout), dimension(:,:) :: dt_ln_ps
complex, intent(in), dimension(:,:,:,:) :: divs, ts
complex, intent(in), dimension(:,:,:) :: ln_ps
integer, intent(in) :: previous, current

complex, dimension(size(ts,1),size(ts,2),size(ts,3)) :: divs_temp, ts_temp, geopot, dt_ts_temp, ln_p_full_temp
complex, dimension(size(ts,1),size(ts,2),size(ts,3)+1) :: ln_p_half_temp
complex, dimension(size(ts,1),size(ts,2)) :: ps_temp, dt_ps_temp

integer :: k

divs_temp = divs(:,:,:,previous) - divs(:,:,:,current)

call linear_tp_tendency(divs_temp, ref_temperature_implicit, ref_surf_p_implicit, &
                        ref_ln_p_half, ref_ln_p_full, dt_ps_temp, dt_ts_temp)

dt_ts    = dt_ts    + dt_ts_temp
dt_ln_ps = dt_ln_ps + dt_ps_temp/ref_surf_p_implicit

ts_temp =    ts(:,:,:,previous) -    ts(:,:,:,current) + xi*dt_ts   
ps_temp = ln_ps(:,:,  previous) - ln_ps(:,:,  current) + xi*dt_ln_ps 

ln_p_full_temp = (0.,0.)
ln_p_half_temp = (0.,0.)

geopot = linear_geopotential(ts_temp,ln_p_half_temp,ln_p_full_temp, ref_temperature_implicit, ref_ln_p_half, ref_ln_p_full) 

do k=1,num_levels
  dt_divs(:,:,k)=dt_divs(:,:,k)+eigen(ms:me,ns:ne)*(geopot(:,:,k)+h(k)*ps_temp(:,:)*ref_surf_p_implicit)
enddo

return
end subroutine adjust_dt_divs

!-----------------------------------------------------------------------

function linear_geopotential_3d(del_t, del_ln_p_half, del_ln_p_full, t, ln_p_half, ln_p_full) result (geopot)
!||zn, pjk 5/13/99, following vb's dmsm/implicit.f90. (0:,0:) -> (:,:)

complex, intent (in), dimension (:,:,:) :: del_t
complex, intent (in), dimension (:,:,:) :: del_ln_p_half
complex, intent (in), dimension (:,:,:) :: del_ln_p_full
real,    intent (in), dimension (:)       :: t
real,    intent (in), dimension (:)       :: ln_p_half
real,    intent (in), dimension (:)       :: ln_p_full

complex, dimension (size(del_t,1),size(del_t,2),size(del_t,3))   :: geopot
complex, dimension (size(del_t,1),size(del_t,2),size(del_t,3)+1) :: geopot_half

integer :: k

geopot_half(:,:,num_levels+1) = (0.,0.)

do k=num_levels,2,-1
  geopot_half(:,:,k) = geopot_half(:,:,k+1)                                &
    + rdgas*(del_t(:,:,k)*(ln_p_half(k+1)         - ln_p_half(k))    &
    +t(k)        *(del_ln_p_half(:,:,k+1) - del_ln_p_half(:,:,k)))
enddo

do k=1,num_levels
  geopot(:,:,k) = geopot_half(:,:,k+1)                                     &
    + rdgas*(del_t(:,:,k)*(ln_p_half(k+1)     - ln_p_full(k))        &
    +t(k)        *(del_ln_p_half(:,:,k+1) - del_ln_p_full(:,:,k)))
enddo

return
end function linear_geopotential_3d
!-----------------------------------------------------------------------

function linear_geopotential_1d(del_t,del_ln_p_half,del_ln_p_full, t, ln_p_half, ln_p_full) result (geopot)

real, intent (in), dimension (:) :: del_t
real, intent (in), dimension (:) :: del_ln_p_half
real, intent (in), dimension (:) :: del_ln_p_full
real, intent (in), dimension (:) :: t
real, intent (in), dimension (:) :: ln_p_half
real, intent (in), dimension (:) :: ln_p_full

real, dimension (size(del_t,1)) :: geopot

complex, dimension(0:0, 0:0, size(del_t,1)  ) :: del_t_3d, geopot_3d
complex, dimension(0:0, 0:0, size(del_t,1)+1) :: del_ln_p_half_3d
complex, dimension(0:0, 0:0, size(del_t,1)  ) :: del_ln_p_full_3d

del_t_3d(0,0,:)         = cmplx(del_t,         0.0)
del_ln_p_half_3d(0,0,:) = cmplx(del_ln_p_half, 0.0)
del_ln_p_full_3d(0,0,:) = cmplx(del_ln_p_full, 0.0)

geopot_3d = linear_geopotential(del_t_3d, del_ln_p_half_3d, del_ln_p_full_3d, t, ln_p_half, ln_p_full)

geopot = real(geopot_3d(0,0,:))

return
end function linear_geopotential_1d
!-----------------------------------------------------------------------

function pres_grad_funct(tg, ln_p_half, ln_p_full, p_surf) result(x)

real, intent (in), dimension (:) :: tg, ln_p_half, ln_p_full
real, intent (in)                :: p_surf
real, dimension (size(tg,1)) :: x

real :: dlog_1, dlog_2
integer :: k

if(trim(vert_difference_option) == 'simmons_and_burridge') then
  do k=1,size(ln_p_full,1)
    dlog_1 = ln_p_half(k+1) - ln_p_full(k)
    dlog_2 = ln_p_full(k)   - ln_p_half(k)
    x(k) = rdgas*tg(k) * (bk(k+1)*dlog_1 + bk(k)*dlog_2)/(dpk(k) + dbk(k)*p_surf)
  end do
else if(trim(vert_difference_option) == 'mcm') then
  do k=1,size(ln_p_full,1)
    x(k) = rdgas*tg(k)/p_surf
  end do
endif

return
end function pres_grad_funct
!-----------------------------------------------------------------------

subroutine linear_tp_tendency_3d (div, t_ref, p_surf_ref, ln_p_half_ref, ln_p_full_ref, dt_p_surf, dt_t)

complex, intent (in), dimension (0:,0:,:) :: div
real,    intent (in), dimension (:) :: t_ref, ln_p_half_ref, ln_p_full_ref
real,    intent (in)  :: p_surf_ref

complex, intent (out), dimension (0:,0:)   :: dt_p_surf
complex, intent (out), dimension (0:,0:,:) :: dt_t

real ::  dp, dp_inv, dlog_1, dlog_3, p_full_ref
complex , dimension(0:size(div,1)-1,0:size(div,2)-1) :: dmean, dmean_tot
complex , dimension(0:size(div,1)-1,0:size(div,2)-1,size(div,3)+1) :: vert_vel, temp

real :: kappa
integer :: k

kappa = rdgas/cp_air

dmean_tot = 0.0

if(vert_difference_option == 'simmons_and_burridge') then

  do k=1,num_levels
    dp = dpk(k) + dbk(k)*p_surf_ref
    dp_inv = 1/dp
    dlog_1 = ln_p_half_ref(k+1) - ln_p_full_ref(k)
    dlog_3 = ln_p_half_ref(k+1) - ln_p_half_ref(k)
    dmean = div(:,:,k)*dp
    dt_t(:,:,k) =  - kappa*t_ref(k)*(dmean_tot*dlog_3 + dmean*dlog_1)*dp_inv
    dmean_tot = dmean_tot + dmean
    vert_vel(:,:,k+1) = - dmean_tot
  enddo

else if (vert_difference_option == 'mcm') then

  do k=1,num_levels
    dp = dpk(k) + dbk(k)*p_surf_ref
    dmean = div(:,:,k)*dp
    p_full_ref = 0.5*(pk(k+1) + pk(k)) + 0.5*(bk(k+1) + bk(k))*p_surf_ref
    dt_t(:,:,k) = - (kappa*t_ref(k)/p_full_ref)*(dmean_tot + 0.5*dmean)
    dmean_tot = dmean_tot + dmean
    vert_vel(:,:,k+1) = - dmean_tot
  enddo

endif

dt_p_surf = - dmean_tot

do k=2,num_levels
  vert_vel(:,:,k) = vert_vel(:,:,k) + dmean_tot*bk(k)
enddo

do k=2,num_levels
  temp(:,:,k) = - vert_vel(:,:,k)*(t_ref(k) - t_ref(k-1))
enddo

temp(:,:,1)            = 0.0
temp(:,:,num_levels+1) = 0.0

do k=1,num_levels
  dp = dpk(k) + dbk(k)*p_surf_ref
  dp_inv = 1/dp
  dt_t(:,:,k) = dt_t(:,:,k) + .5*dp_inv*(temp(:,:,k+1)+temp(:,:,k))
enddo

return
end subroutine linear_tp_tendency_3d
!-----------------------------------------------------------------------

subroutine linear_tp_tendency_1d(div, t_ref, p_surf_ref, ln_p_half_ref, ln_p_full_ref, dt_p_surf, dt_t)

real, intent(in), dimension (:) :: div
real, intent(in), dimension (:) :: t_ref, ln_p_half_ref, ln_p_full_ref
real, intent(in)  :: p_surf_ref

real, intent(out) :: dt_p_surf
real, intent(out), dimension (:) :: dt_t

complex, dimension(0:0,0:0,size(div,1)) :: div_3d, dt_t_3d
complex, dimension(0:0,0:0)           :: dt_p_surf_3d

div_3d(0,0,:) = cmplx(div,0.)

call linear_tp_tendency(div_3d, t_ref, p_surf_ref, ln_p_half_ref, ln_p_full_ref, dt_p_surf_3d, dt_t_3d)

dt_t      = real(dt_t_3d(0,0,:))
dt_p_surf = real(dt_p_surf_3d(0,0))

return
end subroutine linear_tp_tendency_1d
!-----------------------------------------------------------------------
subroutine implicit_end

if(.not.module_is_initialized) return

deallocate(ref_temperature_implicit)
deallocate(ref_ln_p_half, ref_ln_p_full, del_ln_p_half, del_ln_p_full)
deallocate(div_mat, wave_matrix, h, eigen, wavenumber, pk,  bk, dpk, dbk)
module_is_initialized = .false.

return
end subroutine implicit_end
!-----------------------------------------------------------------------

end module implicit_mod
