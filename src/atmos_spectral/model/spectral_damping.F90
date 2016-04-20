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

module spectral_damping_mod

use        fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                          write_version_number

use transforms_mod, only: get_eigen_laplacian, get_spec_domain

implicit none

private

interface compute_spectral_damping
   module procedure  compute_spectral_damping_2d, compute_spectral_damping_3d 
end interface

real, allocatable, dimension(:,:) :: damping, damping_vor, damping_div, damping_eddy_sponge
real, allocatable, dimension(:)   :: damping_zmu_sponge, damping_zmv_sponge
logical :: module_is_initialized = .false.
integer :: ms, me, ns, ne, num_levels

! begin LJJ addition via ST
real    :: damping_coeff_vor_local, damping_coeff_div_local, damping_coeff_local
integer :: damping_order_vor_local, damping_order_div_local, damping_order_local
logical :: damping_option_exponential = .false.
! end LJJ addition via ST

character(len=128), parameter :: version = '$Id: spectral_damping.F90,v 17.0 2009/07/21 03:00:53 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

public :: spectral_damping_init, spectral_damping_end, compute_spectral_damping
public :: compute_spectral_damping_vor, compute_spectral_damping_div

contains
!----------------------------------------------------------------------------------------------------------------
subroutine spectral_damping_init (damping_coeff, damping_order, damping_option, cutoff_wn, num_fourier, num_spherical, &
                                  num_levels_in, eddy_sponge_coeff, zmu_sponge_coeff, zmv_sponge_coeff,     &
                                  damping_coeff_vor, damping_order_vor, damping_coeff_div, damping_order_div, &
                                  damping_coeff_r)

real,    intent(in) :: damping_coeff
integer, intent(in) :: cutoff_wn
integer, intent(in) :: damping_order, num_fourier, num_spherical, num_levels_in
real,    intent(in) :: eddy_sponge_coeff, zmu_sponge_coeff, zmv_sponge_coeff
character(len=*), intent(in) :: damping_option
real,    intent(in), optional :: damping_coeff_vor, damping_coeff_div
integer, intent(in), optional :: damping_order_vor, damping_order_div
real,    intent(in), optional :: damping_coeff_r ! linear drag coefficient (units = 1/seconds)

! begin LJJ addition via ST - These vars now defined globally within module
! real    :: damping_coeff_vor_local, damping_coeff_div_local
! integer :: damping_order_vor_local, damping_order_div_local
! end LJJ addition via ST

real, dimension(0:num_fourier, 0:num_spherical) :: eigen, sqrt_eigen

if(module_is_initialized) return

call write_version_number(version, tagname)
num_levels = num_levels_in

call get_spec_domain(ms, me, ns, ne)

allocate(damping_zmu_sponge(0:num_spherical))
allocate(damping_zmv_sponge(0:num_spherical))
allocate(damping_eddy_sponge(0:num_fourier, 0:num_spherical))
allocate(damping    (0:num_fourier, 0:num_spherical))
allocate(damping_vor(0:num_fourier, 0:num_spherical))
allocate(damping_div(0:num_fourier, 0:num_spherical))

call get_eigen_laplacian(eigen)

! begin LJJ addition via ST
sqrt_eigen = sqrt(eigen)
! end LJJ addition via ST

if(present(damping_coeff_vor)) then
  damping_coeff_vor_local = damping_coeff_vor
else
  damping_coeff_vor_local = damping_coeff
endif

if(present(damping_coeff_div)) then
  damping_coeff_div_local = damping_coeff_div
else
  damping_coeff_div_local = damping_coeff
endif

if(present(damping_order_vor)) then
  damping_order_vor_local = damping_order_vor
else
  damping_order_vor_local = damping_order
endif

if(present(damping_order_div)) then
  damping_order_div_local = damping_order_div
else
  damping_order_div_local = damping_order
endif

! begin LJJ addition via ST
damping_coeff_local = damping_coeff
damping_order_local = damping_order
! end LJJ addition via ST

if(trim(damping_option) == 'resolution_dependent') then
  damping    (:,:) = damping_coeff          *((eigen(:,:)/eigen(0,num_spherical-1))**damping_order          )
  damping_vor(:,:) = damping_coeff_vor_local*((eigen(:,:)/eigen(0,num_spherical-1))**damping_order_vor_local)
  damping_div(:,:) = damping_coeff_div_local*((eigen(:,:)/eigen(0,num_spherical-1))**damping_order_div_local)
  ! begin LJJ addition via ST
else if(trim(damping_option) == 'exponential_cutoff') then
   !  implement exponential filter similar to K.S. Smith et al JFM 2002
   !  note typo in equation B6 of paper (minus sign)
   !  note damping coefficients are further modified in the compute routines below

   damping_option_exponential = .true.

   where (eigen(:,:)/eigen(0, cutoff_wn) > 1)
     damping    (:,:) = ((sqrt_eigen(:,:) - sqrt_eigen(0, cutoff_wn) ) & 
        / ( sqrt_eigen(0, num_spherical-1) - sqrt_eigen(0, cutoff_wn))  )**damping_order_local
     damping_vor(:,:) = damping
     damping_div(:,:) = damping
   elsewhere
     damping    (:,:) = 0.0
     damping_vor(:,:) = 0.0
     damping_div(:,:) = 0.0
   end where
! end LJJ addition via ST
else if(trim(damping_option) == 'resolution_independent') then
  damping    (:,:) = damping_coeff          *(eigen(:,:)**damping_order          )
  damping_vor(:,:) = damping_coeff_vor_local*(eigen(:,:)**damping_order_vor_local)
  damping_div(:,:) = damping_coeff_div_local*(eigen(:,:)**damping_order_div_local)
else
  call error_mesg('spectral_damping_init','"'//trim(damping_option)//'" is an invalid value for damping_option', FATAL)
end if

if(present(damping_coeff_r)) then
  damping(:,:) = damping(:,:) + damping_coeff_r
endif

damping_zmu_sponge  =  zmu_sponge_coeff*eigen(0,:)
damping_zmv_sponge  =  zmv_sponge_coeff*eigen(0,:)
damping_eddy_sponge = eddy_sponge_coeff*eigen

module_is_initialized = .true.

return
end subroutine spectral_damping_init

!-----------------------------------------------------------------
!compute_spectral_damping_3d is local
subroutine compute_spectral_damping_3d(spec, dt_spec, current_dt)  

complex, intent(in), dimension(:,:,:) :: spec
real,    intent(in) :: current_dt
complex, intent(inout), dimension (:,:,:) :: dt_spec

real, dimension(size(spec,1),size(spec,2)) :: coeff, damping_eff

integer :: k

if(.not.module_is_initialized) then
  call error_mesg('compute_spectral_damping','initialization of spectral_damping not performed', FATAL)
end if

! begin LJJ addition via ST
if(damping_option_exponential) then
 damping_eff = (exp(log(current_dt*damping_coeff_local+1)*damping(ms:me,ns:ne))-1.0)/current_dt
else
 damping_eff = damping(ms:me, ns:ne)
endif

coeff = 1.0/(1.0 + damping_eff*current_dt)

do k = 1, size(spec,3)
  dt_spec(:,:,k) =   coeff(:,:) * (dt_spec(:,:,k) - damping_eff*spec(:,:,k))
end do
! end LJJ addition via ST

return
end subroutine compute_spectral_damping_3d
!-----------------------------------------------------------------
subroutine compute_spectral_damping_vor(vor, dt_vor, current_dt)

complex, intent(in), dimension(ms:me, ns:ne, num_levels) :: vor
real,    intent(in) :: current_dt
complex, intent(inout), dimension(ms:me, ns:ne, num_levels) :: dt_vor

real, dimension(ms:me, ns:ne) :: coeff, damping_vor_eff

integer :: k, m, n

if(.not.module_is_initialized) then
  call error_mesg('compute_spectral_damping_vor','initialization of spectral_damping not performed', FATAL)
end if

! begin LJJ addition via ST
if(damping_option_exponential) then
 damping_vor_eff = (exp(log(current_dt*damping_coeff_vor_local+1)*damping_vor(ms:me,ns:ne))-1.0)/current_dt
else
 damping_vor_eff = damping_vor(ms:me,ns:ne)
endif

coeff = 1.0/(1.0 + damping_vor_eff*current_dt)

do k = 1, size(vor,3)
  dt_vor(:,:,k) = coeff * (dt_vor(:,:,k) - damping_vor_eff*vor(:,:,k))
end do
! end LJJ addition via ST

do n=ns,ne
  do m=ms,me
    if(m /= 0) then
      dt_vor(m,n,1) = (dt_vor(m,n,1) - damping_eddy_sponge(m,n)*vor(m,n,1))/(1. + damping_eddy_sponge(m,n)*current_dt)
    endif
  enddo
enddo

if(ms == 0) then
  do n=ns,ne
    dt_vor(0,n,1) = (dt_vor(0,n,1) - damping_zmu_sponge(n)*vor(0,n,1))/(1. + damping_zmu_sponge(n)*current_dt)
  enddo
endif

return
end subroutine compute_spectral_damping_vor
!-----------------------------------------------------------------

subroutine compute_spectral_damping_div(div, dt_div, current_dt)
complex, intent(in), dimension(ms:me, ns:ne, num_levels) :: div
real,    intent(in) :: current_dt
complex, intent(inout), dimension (ms:me, ns:ne, num_levels) :: dt_div

real, dimension(ms:me, ns:ne) :: coeff, damping_div_eff

integer :: k, m, n

if(.not.module_is_initialized) then
  call error_mesg('compute_spectral_damping_div','initialization of spectral_damping not performed', FATAL)
end if

! begin LJJ addition via ST
if(damping_option_exponential) then
 damping_div_eff = (exp(log(current_dt*damping_coeff_div_local+1)*damping_div(ms:me,ns:ne))-1.0)/current_dt
else
 damping_div_eff = damping_div(ms:me,ns:ne)
endif

coeff = 1.0/(1.0 + damping_div_eff*current_dt)

do k = 1, size(div,3)
  dt_div(:,:,k) = coeff * (dt_div(:,:,k) - damping_div_eff*div(:,:,k))
end do
! end LJJ addition via ST

do n=ns,ne
  do m=ms,me
    if(m /= 0) then
      dt_div(m,n,1) = (dt_div(m,n,1) - damping_eddy_sponge(m,n)*div(m,n,1))/(1. + damping_eddy_sponge(m,n)*current_dt)
    endif
  enddo
enddo

if(ms == 0) then
  do n=ns,ne
    dt_div(0,n,1) = (dt_div(0,n,1) - damping_zmv_sponge(n)*div(0,n,1))/(1. + damping_zmv_sponge(n)*current_dt)
  enddo
endif

return
end subroutine compute_spectral_damping_div
!-----------------------------------------------------------------
subroutine compute_spectral_damping_2d(spec, dt_spec, current_dt)  

complex, intent(in), dimension(:,:) :: spec
real, intent(in) :: current_dt
complex, intent(inout), dimension(:,:) :: dt_spec

real, dimension(size(spec,1),size(spec,2)) :: coeff, damping_eff

if(.not.module_is_initialized) then
  call error_mesg('compute_spectral_damping','initialization of spectral_damping not performed', FATAL)
end if

! begin LJJ addition via ST
if(damping_option_exponential) then
 damping_eff = (exp(log(current_dt*damping_coeff_local+1)*damping(ms:me,ns:ne))-1.0)/current_dt
else
 damping_eff = damping(ms:me,ns:ne)
endif

coeff = 1.0/(1.0 + damping_eff*current_dt)
dt_spec = coeff*(dt_spec - damping_eff*spec)
! end LJJ addition via ST

return
end subroutine compute_spectral_damping_2d

!-----------------------------------------------------------------
subroutine spectral_damping_end

if(.not.module_is_initialized) return
deallocate(damping, damping_vor, damping_div)
deallocate(damping_eddy_sponge, damping_zmu_sponge, damping_zmv_sponge)
module_is_initialized = .false.

return
end subroutine spectral_damping_end
!-----------------------------------------------------------------

end module spectral_damping_mod
