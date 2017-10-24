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

module spherical_mod

!Balaji: parallel spectral model using transpose method
!initialize spherical operates on global domain
!all other operators attempt to detect whether argument is global (0:M,0:N)
!                                                       or local (ms:me,ns:ne)

use       fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                         write_version_number

use  spec_mpp_mod, only: get_spec_domain

!-------------------------------------------------------------------------
!   provides operations on spectral spherical harmonics fields that do not 
!      require transforms
!
!   spectral fields are complex with horizontal dimensions
!        (0:num_fourier,0:num_spherical)
!       where the zonal wavenumber of (m,n) is M = m*fourier_inc
!             the "meridional" wavenumber is n
!             the "total", 2D, spherical wavenumber L = M+n
!
!-----------------------------------------------------------------------

implicit none
private

character(len=128), parameter :: version = '$Id: spherical.F90,v 13.0 2006/03/28 21:18:30 fms Exp $'
character(len=128), parameter :: tagname ='$Name: siena_201211 $'

interface compute_lon_deriv_cos
   module procedure compute_lon_deriv_cos_2d,  &
                    compute_lon_deriv_cos_3d
end interface

interface compute_lat_deriv_cos
   module procedure compute_lat_deriv_cos_2d,  &
                    compute_lat_deriv_cos_3d
end interface

interface compute_gradient_cos
   module procedure compute_gradient_cos_2d,  &
                    compute_gradient_cos_3d
end interface

interface compute_laplacian
   module procedure compute_laplacian_2d,  &
                    compute_laplacian_3d
end interface

interface compute_ucos_vcos
   module procedure compute_ucos_vcos_2d,  &
                    compute_ucos_vcos_3d
end interface

interface compute_vor_div
   module procedure compute_vor_div_2d,  &
                    compute_vor_div_3d
end interface

interface compute_vor
   module procedure compute_vor_2d,  &
                    compute_vor_3d
end interface

interface compute_div
   module procedure compute_div_2d,  &
                    compute_div_3d
end interface

interface compute_alpha_operator
   module procedure compute_alpha_operator_2d,  &
                    compute_alpha_operator_3d
end interface

interface triangular_truncation
   module procedure triangular_truncation_2d,  &
                    triangular_truncation_3d
end interface

interface rhomboidal_truncation
   module procedure rhomboidal_truncation_2d,  &
                    rhomboidal_truncation_3d
end interface

public :: spherical_init, spherical_end
public :: compute_lon_deriv_cos, compute_lat_deriv_cos
public :: compute_laplacian, compute_vor, compute_div
public :: get_spherical_wave, get_fourier_wave
public :: get_eigen_laplacian, compute_gradient_cos
public :: compute_ucos_vcos, compute_vor_div
public :: triangular_truncation, rhomboidal_truncation

real,    allocatable, dimension(:,:) :: eigen_laplacian, epsilon, r_epsilon
integer, allocatable, dimension(:,:) :: fourier_wave, spherical_wave

integer :: fourier_max, fourier_inc, num_fourier, num_spherical
integer :: ms, me, ns, ne
logical :: module_is_initialized = .false.

logical :: make_symmetric_local

real, allocatable, dimension(:,:) :: coef_uvm
real, allocatable, dimension(:,:) :: coef_uvc
real, allocatable, dimension(:,:) :: coef_uvp
real, allocatable, dimension(:,:) :: coef_alpm
real, allocatable, dimension(:,:) :: coef_alpp
real, allocatable, dimension(:,:) :: coef_dym
real, allocatable, dimension(:,:) :: coef_dx
real, allocatable, dimension(:,:) :: coef_dyp
real, allocatable, dimension(:,:) :: triangle_mask

contains

!--------------------------------------------------------------------------
subroutine spherical_init(radius, num_fourier_in, fourier_inc_in, num_spherical_in, make_symmetric)

real,    intent(in) :: radius
integer, intent(in) :: num_fourier_in, fourier_inc_in, num_spherical_in
logical, intent(in), optional :: make_symmetric

integer :: m,n

if(module_is_initialized) return

call write_version_number(version, tagname)

call get_spec_domain(ms, me, ns, ne)

num_fourier = num_fourier_in
fourier_inc = fourier_inc_in
fourier_max = fourier_inc*num_fourier
num_spherical = num_spherical_in

if(present(make_symmetric)) then
  make_symmetric_local = make_symmetric
else
  make_symmetric_local = .false.
end if

allocate(eigen_laplacian(0:num_fourier,0:num_spherical))
allocate(epsilon        (0:num_fourier,0:num_spherical))
allocate(r_epsilon      (0:num_fourier,0:num_spherical))
allocate(fourier_wave   (0:num_fourier,0:num_spherical))
allocate(spherical_wave (0:num_fourier,0:num_spherical))
allocate(coef_uvm       (0:num_fourier,0:num_spherical))
allocate(coef_uvc       (0:num_fourier,0:num_spherical))
allocate(coef_uvp       (0:num_fourier,0:num_spherical))
allocate(coef_alpm      (0:num_fourier,0:num_spherical))
allocate(coef_alpp      (0:num_fourier,0:num_spherical))
allocate(coef_dym       (0:num_fourier,0:num_spherical))
allocate(coef_dx        (0:num_fourier,0:num_spherical))
allocate(coef_dyp       (0:num_fourier,0:num_spherical))
allocate(triangle_mask  (0:num_fourier,0:num_spherical))

module_is_initialized = .true.

triangle_mask = 1.0
do n=0,num_spherical
  do m=0,num_fourier
    fourier_wave(m,n)   = m*fourier_inc
    spherical_wave(m,n) = fourier_wave(m,n) + n
    if(spherical_wave(m,n).gt.num_spherical-1) triangle_mask(m,n) = 0.0
	if(make_symmetric_local .AND. m.gt.0) triangle_mask(m,n) = 0.0 !GC make zonally symmetric
  end do
end do

epsilon= (spherical_wave**2 - fourier_wave**2)/(4.0*spherical_wave**2 - 1.0)
epsilon = sqrt(epsilon)

eigen_laplacian = spherical_wave*(spherical_wave + 1.0)/(radius*radius)

where (spherical_wave > 0) 
  coef_uvm = -radius*epsilon/spherical_wave
  coef_uvc = -radius*fourier_wave/(spherical_wave*(spherical_wave + 1.0))
else where 
  coef_uvm = 0.0
  coef_uvc = 0.0
end where

coef_uvp(:,0:num_spherical-1) =    &
    -radius*epsilon(:,1:num_spherical)/   &
     (spherical_wave(:,0:num_spherical-1) +1.0)

coef_alpm= (spherical_wave + 1.0)*epsilon/radius
coef_alpp(:,0:num_spherical-1) =    &
     spherical_wave(:,0:num_spherical-1)*epsilon(:,1:num_spherical)/radius

coef_dym = (spherical_wave - 1.0)*epsilon/radius
coef_dx = float(fourier_wave)/radius
coef_dyp(:,0:num_spherical-1) =                                           &
  (spherical_wave(:,0:num_spherical-1) + 2.0)*epsilon(:,1:num_spherical)/radius

return
end subroutine spherical_init

!---------------------------------------------------------------------------
subroutine get_spherical_wave(spherical_wave_out)
!---------------------------------------------------------------------------

integer, intent(out), dimension (:,:) :: spherical_wave_out

if( size(spherical_wave_out,1).EQ.num_fourier+1 .AND. size(spherical_wave_out,2).EQ.num_spherical+1 )then
    spherical_wave_out = spherical_wave
else if( size(spherical_wave_out,1).EQ.me-ms+1 .AND. size(spherical_wave_out,2).EQ.ne-ns+1 )then
    spherical_wave_out = spherical_wave(ms:me,ns:ne)
else
    call error_mesg( 'get_spherical_wave', 'invalid argument size', FATAL )
endif

return 
end subroutine get_spherical_wave

!---------------------------------------------------------------------------
subroutine get_fourier_wave(fourier_wave_out)
!---------------------------------------------------------------------------

integer, intent(out), dimension (:,:) :: fourier_wave_out

if( size(fourier_wave_out,1).EQ.num_fourier+1 .AND. size(fourier_wave_out,2).EQ.num_spherical+1 )then
    fourier_wave_out = fourier_wave
else if( size(fourier_wave_out,1).EQ.me-ms+1 .AND. size(fourier_wave_out,2).EQ.ne-ns+1 )then
    fourier_wave_out = fourier_wave(ms:me,ns:ne)
else
    call error_mesg( 'get_fourier_wave', 'invalid argument size', FATAL )
endif

return
end subroutine get_fourier_wave

!---------------------------------------------------------------------------
subroutine get_eigen_laplacian(eigen_laplacian_out)
!---------------------------------------------------------------------------

real, intent(out), dimension (:,:) :: eigen_laplacian_out

if( size(eigen_laplacian_out,1).EQ.num_fourier+1 .AND. size(eigen_laplacian_out,2).EQ.num_spherical+1 )then
    eigen_laplacian_out = eigen_laplacian
else if( size(eigen_laplacian_out,1).EQ.me-ms+1 .AND. size(eigen_laplacian_out,2).EQ.ne-ns+1 )then
    eigen_laplacian_out = eigen_laplacian(ms:me,ns:ne)
else
    call error_mesg( 'get_eigen_laplacian', 'invalid argument size', FATAL )
endif

return
end subroutine get_eigen_laplacian

!---------------------------------------------------------------------------
function compute_lon_deriv_cos_3d(spherical) result(deriv_lon)
!---------------------------------------------------------------------------

complex, intent(in), dimension (:,:,:) :: spherical
complex, dimension (size(spherical,1), size(spherical,2), size(spherical,3)) :: deriv_lon

integer :: k

if(.not. module_is_initialized ) then
  call error_mesg('compute_lon_deriv','module spherical not initialized', FATAL)
end if

if( size(spherical,1).EQ.num_fourier+1 .AND. size(spherical,2).EQ.num_spherical+1 )then
    do k=1,size(spherical,3)
       deriv_lon(:,:,k) = coef_dx(:,:)*cmplx(-aimag(spherical(:,:,k)),real(spherical(:,:,k)))
    end do
else if( size(spherical,1).EQ.me-ms+1 .AND. size(spherical,2).EQ.ne-ns+1 )then
    do k=1,size(spherical,3)
       deriv_lon(:,:,k) = coef_dx(ms:me,ns:ne)*cmplx(-aimag(spherical(:,:,k)),real(spherical(:,:,k)))
    end do
else
    call error_mesg( 'compute_lon_deriv_cos_3d', 'invalid argument size', FATAL )
endif

return
end function compute_lon_deriv_cos_3d

!--------------------------------------------------------------------------
function compute_lat_deriv_cos_3d(spherical) result(deriv_lat)
!--------------------------------------------------------------------------

complex, intent(in), dimension (:,0:,:) :: spherical
complex, dimension (size(spherical,1), 0:size(spherical,2)-1, size(spherical,3)) :: deriv_lat

integer :: k

if(.not. module_is_initialized ) then
  call error_mesg('compute_lat_deriv','module spherical not initialized', FATAL)
end if

if( size(spherical,2).EQ.num_spherical+1 )then
!could be global domain, or only global in N
    if( size(spherical,1).EQ.num_fourier+1 )then
        deriv_lat(:,0,:) = cmplx(0.,0.)
        do k=1,size(spherical,3)
           deriv_lat(:,1:num_spherical,k) = &  
                - spherical(:,0:num_spherical-1,k)*coef_dym(:,1:num_spherical)
           deriv_lat(:,0:num_spherical-1,k) =  deriv_lat(:,0:num_spherical-1,k)   &
                + spherical(:,1:num_spherical,k)*coef_dyp(:,0:num_spherical-1)
        end do
    else if( size(spherical,1).EQ.me-ms+1 )then
        deriv_lat(:,0,:) = cmplx(0.,0.)
        do k=1,size(spherical,3)
           deriv_lat(:,1:num_spherical,k) = &  
                - spherical(:,0:num_spherical-1,k)*coef_dym(ms:me,1:num_spherical)
           deriv_lat(:,0:num_spherical-1,k) =  deriv_lat(:,0:num_spherical-1,k)   &
                + spherical(:,1:num_spherical,k)*coef_dyp(ms:me,0:num_spherical-1)
        end do
    endif
else if( size(spherical,1).EQ.me-ms+1 .AND. size(spherical,2).EQ.ne-ns+1 )then
!need to write stuff to acquire data at ns-1,ne+1
    call abort()
else
    call error_mesg( 'compute_lat_deriv_cos_3d', 'invalid argument size', FATAL )
endif

return
end function compute_lat_deriv_cos_3d

!-------------------------------------------------------------------
subroutine compute_gradient_cos_3d(spherical, deriv_lon, deriv_lat) 
!-------------------------------------------------------------------

complex, intent(in), dimension (:,:,:) :: spherical
complex, intent(out), dimension (:,:,:) :: deriv_lat
complex, intent(out), dimension (:,:,:) :: deriv_lon

deriv_lon = compute_lon_deriv_cos(spherical)
deriv_lat = compute_lat_deriv_cos(spherical)

return
end subroutine compute_gradient_cos_3d

!----------------------------------------------------------------------
function compute_laplacian_3d(spherical, power) result(laplacian)
!----------------------------------------------------------------------

  complex, intent(in), dimension (:,:,:) :: spherical
  integer, optional :: power

  complex, dimension (size(spherical,1), size(spherical,2), size(spherical,3)) :: laplacian

  integer :: k
  real, dimension(size(spherical,1), size(spherical,2)) :: factor

  if(.not. module_is_initialized ) then
      call error_mesg('compute_laplacian','module spherical not initialized', FATAL)
  end if

  if( size(spherical,1).EQ.num_fourier+1 .AND. size(spherical,2).EQ.num_spherical+1 )then
      if(present(power)) then 
          if(power >= 0) then
              factor = (-eigen_laplacian)**power
          else 
              where (eigen_laplacian .ne. 0.0) 
                  factor = (-eigen_laplacian)**power
              else where
                  factor = 0.0
              end where
          end if
      else
          factor = -eigen_laplacian
      end if
  else if( size(spherical,1).EQ.me-ms+1 .AND. size(spherical,2).EQ.ne-ns+1 )then
      if(present(power)) then 
          if(power >= 0) then
              factor = (-eigen_laplacian(ms:me,ns:ne))**power
          else 
              where (eigen_laplacian(ms:me,ns:ne) .ne. 0.0) 
                  factor = (-eigen_laplacian(ms:me,ns:ne))**power
              else where
                  factor = 0.0
              end where
          end if
      else
          factor = -eigen_laplacian(ms:me,ns:ne)
      end if
  else
      call error_mesg( 'compute_laplacian', 'invalid argument size', FATAL )
  endif

  do k= 1,size(spherical,3)
     laplacian(:,:,k) = spherical(:,:,k)*factor
  end do

  return
end function compute_laplacian_3d

!----------------------------------------------------------------------
subroutine compute_ucos_vcos_3d(vorticity , divergence, u_cos, v_cos)
!----------------------------------------------------------------------

complex, intent(in), dimension (:,0:,:) :: vorticity
complex, intent(in), dimension (:,0:,:) :: divergence
complex, intent(out), dimension (:,0:,:) :: u_cos
complex, intent(out), dimension (:,0:,:) :: v_cos


integer :: k

if(.not. module_is_initialized ) then
  call error_mesg('compute_ucos_vcos','module spherical not initialized', FATAL)
end if

if( size(vorticity,2).EQ.num_spherical+1 )then
!could be global domain, or only global in N
    if( size(vorticity,1).EQ.num_fourier+1 )then
        do k=1,size(vorticity,3)
           u_cos(:,:,k) = coef_uvc(:,:)*                                     &
                cmplx(-aimag(divergence(:,:,k)),real(divergence(:,:,k)))
           v_cos(:,:,k) = coef_uvc(:,:)*                                     &
                cmplx(-aimag(vorticity(:,:,k)),real(vorticity(:,:,k)))

           u_cos(:,1:num_spherical,k) = u_cos(:,1:num_spherical,k) +         &
                coef_uvm(:,1:num_spherical)*vorticity(:,0:num_spherical-1,k)
           v_cos(:,1:num_spherical,k) = v_cos(:,1:num_spherical,k) -         &
                coef_uvm(:,1:num_spherical)*divergence(:,0:num_spherical-1,k)

           u_cos(:,0:num_spherical-1,k) = u_cos(:,0:num_spherical-1,k) -     &
                coef_uvp(:,0:num_spherical-1)*vorticity(:,1:num_spherical,k)
           v_cos(:,0:num_spherical-1,k) = v_cos(:,0:num_spherical-1,k) +     &
                coef_uvp(:,0:num_spherical-1)*divergence(:,1:num_spherical,k)         
        end do
    else if( size(vorticity,1).EQ.me-ms+1 )then
        do k=1,size(vorticity,3)
           u_cos(:,:,k) = coef_uvc(ms:me,:)*                                     &
                cmplx(-aimag(divergence(:,:,k)),real(divergence(:,:,k)))
           v_cos(:,:,k) = coef_uvc(ms:me,:)*                                     &
                cmplx(-aimag(vorticity(:,:,k)),real(vorticity(:,:,k)))

           u_cos(:,1:num_spherical,k) = u_cos(:,1:num_spherical,k) +         &
                coef_uvm(ms:me,1:num_spherical)*vorticity(:,0:num_spherical-1,k)
           v_cos(:,1:num_spherical,k) = v_cos(:,1:num_spherical,k) -         &
                coef_uvm(ms:me,1:num_spherical)*divergence(:,0:num_spherical-1,k)

           u_cos(:,0:num_spherical-1,k) = u_cos(:,0:num_spherical-1,k) -     &
                coef_uvp(ms:me,0:num_spherical-1)*vorticity(:,1:num_spherical,k)
           v_cos(:,0:num_spherical-1,k) = v_cos(:,0:num_spherical-1,k) +     &
                coef_uvp(ms:me,0:num_spherical-1)*divergence(:,1:num_spherical,k)         
        end do
    endif
else if( size(vorticity,1).EQ.me-ms+1 .AND. size(vorticity,2).EQ.ne-ns+1 )then
!need to write stuff to acquire data at ns-1,ne+1
    call abort()
else
    call error_mesg( 'compute_ucos_vcos', 'invalid argument size', FATAL )
endif

return
end subroutine compute_ucos_vcos_3d

!-------------------------------------------------------------------------
subroutine compute_vor_div_3d(u_cos, v_cos, vorticity, divergence)
!-------------------------------------------------------------------------

complex, intent(in), dimension (:,:,:) :: u_cos
complex, intent(in), dimension (:,:,:) :: v_cos
complex, intent(out), dimension (:,:,:) :: vorticity
complex, intent(out), dimension (:,:,:) :: divergence

vorticity  = compute_alpha_operator(v_cos, u_cos, -1)
divergence = compute_alpha_operator(u_cos, v_cos, +1)

return
end subroutine compute_vor_div_3d

!-------------------------------------------------------------------------
function compute_vor_3d(u_cos, v_cos) result(vorticity)
!-------------------------------------------------------------------------

complex, intent(in), dimension (:,:,:) :: u_cos
complex, intent(in), dimension (:,:,:) :: v_cos
complex, dimension (size(u_cos,1), size(u_cos,2), size(u_cos,3)) :: vorticity

vorticity = compute_alpha_operator(v_cos, u_cos, -1)

return
end function compute_vor_3d

!-------------------------------------------------------------------------
function compute_div_3d(u_cos, v_cos) result(divergence)
!-------------------------------------------------------------------------

complex, intent(in), dimension (:,:,:) :: u_cos
complex, intent(in), dimension (:,:,:) :: v_cos
complex, dimension (size(u_cos,1), size(u_cos,2), size(u_cos,3)) :: divergence

divergence = compute_alpha_operator(u_cos, v_cos, +1)

return
end function compute_div_3d

!--------------------------------------------------------------------------------
function compute_alpha_operator_3d(spherical_a, spherical_b, isign) result(alpha)
!--------------------------------------------------------------------------------

complex, intent(in), dimension (:,0:,:) :: spherical_a
complex, intent(in), dimension (:,0:,:) :: spherical_b
integer,intent(in) :: isign

complex, dimension (size(spherical_a,1), 0:size(spherical_a,2)-1, size(spherical_a,3)) :: alpha

integer :: k

if(.not. module_is_initialized ) then
  call error_mesg('compute_vor or div','module spherical not initialized', FATAL)
end if

alpha = cmplx(0.,0.)

if( size(spherical_a,2).EQ.num_spherical+1 )then
!could be global domain, or only global in N
    if( size(spherical_a,1).EQ.num_fourier+1 )then
        do k=1,size(spherical_a,3)
           alpha(:,:,k) = coef_dx(:,:)*    &
                cmplx(-aimag(spherical_a(:,:,k)),real(spherical_a(:,:,k)))
           alpha(:,1:num_spherical,k) = alpha(:,1:num_spherical,k) -  &
                isign*coef_alpm(:,1:num_spherical)  &
                *spherical_b(:,0:num_spherical-1,k)
           alpha(:,0:num_spherical-1,k) = alpha(:,0:num_spherical-1,k) +  &
                isign*coef_alpp(:,0:num_spherical-1)*spherical_b(:,1:num_spherical,k)
        end do
    else if( size(spherical_a,1).EQ.me-ms+1 )then
        do k=1,size(spherical_a,3)
           alpha(:,:,k) = coef_dx(ms:me,:)*    &
                cmplx(-aimag(spherical_a(:,:,k)),real(spherical_a(:,:,k)))
           alpha(:,1:num_spherical,k) = alpha(:,1:num_spherical,k) -  &
                isign*coef_alpm(ms:me,1:num_spherical)  &
                *spherical_b(:,0:num_spherical-1,k)
           alpha(:,0:num_spherical-1,k) = alpha(:,0:num_spherical-1,k) +  &
                isign*coef_alpp(ms:me,0:num_spherical-1)*spherical_b(:,1:num_spherical,k)
        end do
    endif
else if( size(spherical_a,1).EQ.me-ms+1 .AND. size(spherical_a,2).EQ.ne-ns+1 )then
!need to write stuff to acquire data at ns-1,ne+1
    call abort()
else
    call error_mesg( 'compute_alpha_operator_3d', 'invalid argument size', FATAL )
endif

return
end function compute_alpha_operator_3d

!-----------------------------------------------------------------------
subroutine triangular_truncation_3d(spherical, trunc)
!-----------------------------------------------------------------------

complex, intent(inout), dimension (:,:,:) :: spherical
integer, intent(in), optional :: trunc
integer :: k

if( size(spherical,1).EQ.num_fourier+1 .AND. size(spherical,2).EQ.num_spherical+1 )then
    if(present(trunc)) then
        do k=1, size(spherical,3)
           where (spherical_wave > trunc)
               spherical(:,:,k) = cmplx(0.,0.)
           end where
        end do
    else
        do k=1, size(spherical,3)
           spherical(:,:,k) = spherical(:,:,k)*triangle_mask
        end do
    end if
else if( size(spherical,1).EQ.me-ms+1 .AND. size(spherical,2).EQ.ne-ns+1 )then
    if(present(trunc)) then
        do k=1, size(spherical,3)
           where (spherical_wave(ms:me,ns:ne) > trunc)
               spherical(:,:,k) = cmplx(0.,0.)
           end where
        end do
    else
        do k=1, size(spherical,3)
           spherical(:,:,k) = spherical(:,:,k)*triangle_mask(ms:me,ns:ne)
        end do
    end if
else
    call error_mesg( 'triang_trunc', 'invalid argument size', FATAL )
endif

return
end subroutine triangular_truncation_3d

!-----------------------------------------------------------------------------
subroutine rhomboidal_truncation_3d(spherical, trunc_fourier, trunc_spherical)
!-----------------------------------------------------------------------------

complex, intent(inout), dimension (0:,0:,:) :: spherical
integer, intent(in), optional:: trunc_fourier, trunc_spherical

integer :: k,n,m

if( size(spherical,1).EQ.num_fourier+1 .AND. size(spherical,2).EQ.num_spherical+1 )then
    if(present(trunc_fourier).and.present(trunc_spherical)) then
        do k=1,size(spherical,3)
           do n=0,size(spherical,2)-1
              do m=0,size(spherical,1)-1
                 if(fourier_wave(m,n) .GT. trunc_fourier.OR.n .GT. trunc_spherical)  &
                      spherical(m,n,k) = cmplx(0.,0.)
              end do
           end do
        end do
    else
        spherical(:,num_spherical,:) = cmplx(0.,0.)
    end if
else if( size(spherical,1).EQ.me-ms+1 .AND. size(spherical,2).EQ.ne-ns+1 )then
!wavenumbers are numbered from (0,0) which is actually (ms,ns)
    if(present(trunc_fourier).and.present(trunc_spherical)) then
        do k=1,size(spherical,3)
           do n=0,size(spherical,2)-1
              do m=0,size(spherical,1)-1
                 if(fourier_wave(m+ms,n+ns).GT.trunc_fourier .OR. n+ns.GT.trunc_spherical)  &
                      spherical(m,n,k) = cmplx(0.,0.)
              end do
           end do
        end do
    else
!truncate on n=num_spherical, only on domain containing n=num_spherical
        if( ne.EQ.num_spherical )spherical(:,ne-ns,:) = cmplx(0.,0.)
    end if
else
    call error_mesg( 'rhomboid_trunc', 'invalid argument size', FATAL )
endif

return
end subroutine rhomboidal_truncation_3d

!rest are 2d versions of 3d routines

!---------------------------------------------------------------------------
function compute_lon_deriv_cos_2d(spherical) result(deriv_lon)
!---------------------------------------------------------------------------

complex, intent(in),  dimension (:,:) :: spherical
complex, dimension (size(spherical,1), size(spherical,2)) :: deriv_lon

complex, dimension (size(spherical,1), size(spherical,2), 1) :: spherical_3d
complex, dimension (size(spherical,1), size(spherical,2), 1) :: deriv_lon_3d

spherical_3d(:,:,1) = spherical(:,:)
deriv_lon_3d = compute_lon_deriv_cos_3d(spherical_3d)
deriv_lon(:,:) = deriv_lon_3d(:,:,1)


return
end function compute_lon_deriv_cos_2d

!--------------------------------------------------------------------------
function compute_lat_deriv_cos_2d(spherical) result(deriv_lat)
!--------------------------------------------------------------------------

complex, intent(in),  dimension (:,:) :: spherical
complex, dimension (size(spherical,1), size(spherical,2)) :: deriv_lat

complex, dimension (size(spherical,1), size(spherical,2), 1) :: spherical_3d
complex, dimension (size(spherical,1), size(spherical,2), 1) :: deriv_lat_3d

spherical_3d(:,:,1) = spherical(:,:)
deriv_lat_3d = compute_lat_deriv_cos_3d(spherical_3d)
deriv_lat(:,:) = deriv_lat_3d(:,:,1)


return
end function compute_lat_deriv_cos_2d

!-------------------------------------------------------------------
subroutine compute_gradient_cos_2d(spherical, deriv_lon, deriv_lat) 
!-------------------------------------------------------------------

complex, intent(in),  dimension(:,:) :: spherical
complex, intent(out), dimension(size(spherical,1), size(spherical,2)) :: deriv_lat
complex, intent(out), dimension(size(spherical,1), size(spherical,2)) :: deriv_lon

complex, dimension(size(spherical,1), size(spherical,2), 1) :: spherical_3d
complex, dimension(size(spherical,1), size(spherical,2), 1) :: deriv_lat_3d
complex, dimension(size(spherical,1), size(spherical,2), 1) :: deriv_lon_3d

spherical_3d(:,:,1) = spherical(:,:)
call compute_gradient_cos_3d(spherical_3d, deriv_lon_3d, deriv_lat_3d)
deriv_lon(:,:) = deriv_lon_3d(:,:,1)
deriv_lat(:,:) = deriv_lat_3d(:,:,1)

return
end subroutine compute_gradient_cos_2d

!----------------------------------------------------------------------
function compute_laplacian_2d(spherical, power) result(laplacian)
!----------------------------------------------------------------------

complex, intent(in), dimension (:,:) :: spherical
integer, optional :: power

complex, dimension (size(spherical,1), size(spherical,2)) :: laplacian

complex, dimension(size(spherical,1), size(spherical,2), 1) :: spherical_3d
complex, dimension(size(spherical,1), size(spherical,2), 1) :: laplacian_3d

spherical_3d(:,:,1) = spherical(:,:)
laplacian_3d = compute_laplacian_3d(spherical_3d, power)
laplacian(:,:) = laplacian_3d(:,:,1)

return
end function compute_laplacian_2d

!----------------------------------------------------------------------
subroutine compute_ucos_vcos_2d(vorticity , divergence, u_cos, v_cos)
!----------------------------------------------------------------------

complex, intent(in),  dimension (:,:) :: vorticity
complex, intent(in),  dimension (size(vorticity,1), size(vorticity,2)) :: divergence
complex, intent(out), dimension (size(vorticity,1), size(vorticity,2)) :: u_cos
complex, intent(out), dimension (size(vorticity,1), size(vorticity,2)) :: v_cos

complex, dimension (size(vorticity,1), size(vorticity,2), 1) :: vorticity_3d
complex, dimension (size(vorticity,1), size(vorticity,2), 1) :: divergence_3d
complex, dimension (size(vorticity,1), size(vorticity,2), 1) :: u_cos_3d
complex, dimension (size(vorticity,1), size(vorticity,2), 1) :: v_cos_3d

vorticity_3d(:,:,1)  = vorticity(:,:)
divergence_3d(:,:,1) = divergence(:,:)
call compute_ucos_vcos_3d(vorticity_3d, divergence_3d, u_cos_3d, v_cos_3d)
u_cos(:,:) = u_cos_3d(:,:,1)
v_cos(:,:) = v_cos_3d(:,:,1)

return
end subroutine compute_ucos_vcos_2d

!-------------------------------------------------------------------------
subroutine compute_vor_div_2d(u_div_cos, v_div_cos, vorticity, divergence)
!-------------------------------------------------------------------------

complex, intent(in),  dimension (:,:) :: u_div_cos
complex, intent(in),  dimension (size(u_div_cos,1), size(u_div_cos,2)) :: v_div_cos
complex, intent(out), dimension (size(u_div_cos,1), size(u_div_cos,2)) :: vorticity
complex, intent(out), dimension (size(u_div_cos,1), size(u_div_cos,2)) :: divergence

complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: u_div_cos_3d
complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: v_div_cos_3d
complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: vorticity_3d
complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: divergence_3d

u_div_cos_3d(:,:,1) = u_div_cos(:,:)
v_div_cos_3d(:,:,1) = v_div_cos(:,:)
call compute_vor_div_3d(u_div_cos_3d, v_div_cos_3d, vorticity_3d, divergence_3d)
vorticity(:,:)  = vorticity_3d(:,:,1)
divergence(:,:) = divergence_3d(:,:,1)

return
end subroutine compute_vor_div_2d

!-------------------------------------------------------------------------
function compute_vor_2d(u_div_cos, v_div_cos) result(vorticity)
!-------------------------------------------------------------------------

complex, intent(in),  dimension (:,:) :: u_div_cos
complex, intent(in),  dimension (:,:) :: v_div_cos
complex, dimension (size(u_div_cos,1), size(u_div_cos,2)) :: vorticity

complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: u_div_cos_3d
complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: v_div_cos_3d
complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: vorticity_3d

u_div_cos_3d(:,:,1) = u_div_cos(:,:)
v_div_cos_3d(:,:,1) = v_div_cos(:,:)
vorticity_3d = compute_vor_3d(u_div_cos_3d, v_div_cos_3d)
vorticity(:,:) = vorticity_3d(:,:,1)

return
end function compute_vor_2d

!-------------------------------------------------------------------------
function compute_div_2d(u_div_cos, v_div_cos) result(divergence)
!-------------------------------------------------------------------------

complex, intent(in),  dimension (:,:) :: u_div_cos
complex, intent(in),  dimension (:,:) :: v_div_cos
complex, dimension (size(u_div_cos,1), size(u_div_cos,2)) :: divergence

complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: u_div_cos_3d
complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: v_div_cos_3d
complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: divergence_3d

u_div_cos_3d(:,:,1) = u_div_cos(:,:)
v_div_cos_3d(:,:,1) = v_div_cos(:,:)
divergence_3d = compute_div_3d(u_div_cos_3d, v_div_cos_3d)
divergence(:,:) = divergence_3d(:,:,1)

return
end function compute_div_2d

!--------------------------------------------------------------------------------
function compute_alpha_operator_2d(spherical_a, spherical_b, isign) result(alpha)
!--------------------------------------------------------------------------------

complex, intent(in),  dimension (:,:) :: spherical_a
complex, intent(in),  dimension (:,:) :: spherical_b
integer, intent(in) :: isign

complex, dimension (size(spherical_a,1), size(spherical_a,2)) :: alpha

complex, dimension (size(spherical_a,1), size(spherical_a,2), 1) :: spherical_a_3d
complex, dimension (size(spherical_a,1), size(spherical_a,2), 1) :: spherical_b_3d
complex, dimension (size(spherical_a,1), size(spherical_a,2), 1) :: alpha_3d

spherical_a_3d(:,:,1) = spherical_a(:,:)
spherical_b_3d(:,:,1) = spherical_b(:,:)
alpha_3d = compute_alpha_operator_3d(spherical_a_3d, spherical_b_3d, isign)
alpha(:,:) = alpha_3d(:,:,1)

return
end function compute_alpha_operator_2d

!-----------------------------------------------------------------------
subroutine triangular_truncation_2d(spherical, trunc)
!-----------------------------------------------------------------------

complex, intent(inout), dimension (:,:) :: spherical
integer, intent(in), optional :: trunc
complex, dimension (size(spherical,1), size(spherical,2), 1) :: spherical_3d

spherical_3d(:,:,1) = spherical(:,:)
call triangular_truncation_3d(spherical_3d, trunc)
spherical(:,:) = spherical_3d(:,:,1)

return
end subroutine triangular_truncation_2d

!--------------------------------------------------------------------------------
subroutine rhomboidal_truncation_2d(spherical, trunc_fourier, trunc_spherical)
!--------------------------------------------------------------------------------

complex, intent(inout), dimension (:,:) :: spherical
integer, intent(in), optional :: trunc_fourier, trunc_spherical

complex, dimension (size(spherical,1), size(spherical,2), 1) :: spherical_3d

spherical_3d(:,:,1) = spherical(:,:)
call rhomboidal_truncation_3d(spherical_3d, trunc_fourier, trunc_spherical)
spherical(:,:) = spherical_3d(:,:,1)

return
end subroutine rhomboidal_truncation_2d

!-----------------------------------------------------------------------
subroutine spherical_end

if(.not.module_is_initialized) return

deallocate(eigen_laplacian, epsilon, r_epsilon, fourier_wave, spherical_wave)
deallocate(coef_uvm, coef_uvc, coef_uvp, coef_alpm, coef_alpp, coef_dym, coef_dx, coef_dyp, triangle_mask)
module_is_initialized = .false.

return
end subroutine spherical_end
!-----------------------------------------------------------------------

end module spherical_mod
