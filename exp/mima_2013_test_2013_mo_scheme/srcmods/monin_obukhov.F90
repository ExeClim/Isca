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

module monin_obukhov_mod


!=======================================================================
!
!                         MONIN-OBUKHOV MODULE
!
!          Routines for computing surface drag coefficients 
!                 from data at the lowest model level 
!              and for computing the profile of fields 
!           between the lowest model level and the ground
!                  using Monin-Obukhov scaling 
!
!=======================================================================


use constants_mod, only: grav, vonkarm
use mpp_mod,       only: input_nml_file
use fms_mod,       only: error_mesg, FATAL, file_exist,   &
                         check_nml_error, open_namelist_file,      &
                         mpp_pe, mpp_root_pe, close_file, stdlog, &
                         write_version_number

implicit none
private

!=======================================================================
 public monin_obukhov_init
 public monin_obukhov_end
 public mo_drag
 public mo_profile
 public mo_diff
 public stable_mix
!=======================================================================

interface mo_drag
    module procedure  mo_drag_0d, mo_drag_1d, mo_drag_2d
end interface


interface mo_profile
    module procedure  mo_profile_0d,   mo_profile_1d,   mo_profile_2d, &
                      mo_profile_0d_n, mo_profile_1d_n, mo_profile_2d_n
end interface

interface mo_diff
    module procedure  mo_diff_0d_n, mo_diff_0d_1, &
                      mo_diff_1d_n, mo_diff_1d_1, &
                      mo_diff_2d_n, mo_diff_2d_1
end interface

interface stable_mix
    module procedure  stable_mix_0d, stable_mix_1d, &
                      stable_mix_2d, stable_mix_3d
end interface


!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: monin_obukhov.F90,v 19.0 2012/01/06 20:10:48 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201211 $'

!=======================================================================

!  DEFAULT VALUES OF NAMELIST PARAMETERS:

real    :: rich_crit      = 2.0
real    :: drag_min       = 1.e-05          
logical :: neutral        = .false.
integer :: stable_option  = 1
real    :: zeta_trans     = 0.5


namelist /monin_obukhov_nml/ rich_crit, neutral, drag_min, &
                             stable_option, zeta_trans

!=======================================================================

!  MODULE VARIABLES

real, parameter    :: small  = 1.e-04
real               :: b_stab, r_crit, sqrt_drag_min, lambda, rich_trans
logical            :: module_is_initialized = .false.


contains

!=======================================================================

subroutine monin_obukhov_init

integer :: unit, ierr, io, logunit

!------------------- read namelist input -------------------------------

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=monin_obukhov_nml, iostat=io)
      ierr = check_nml_error(io,"monin_obukhov_nml")
#else
      if (file_exist('input.nml')) then
         unit = open_namelist_file ()
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=monin_obukhov_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'monin_obukhov_nml')
         enddo
  10     call close_file (unit)
      endif
#endif

!---------- output namelist to log-------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           logunit = stdlog()
           write (logunit, nml=monin_obukhov_nml)
      endif
      
!----------------------------------------------------------------------

if(rich_crit.le.0.25)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'rich_crit in monin_obukhov_mod must be > 0.25', FATAL)

if(drag_min.le.0.0)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'drag_min in monin_obukhov_mod must be >= 0.0', FATAL)

if(stable_option < 1 .or. stable_option > 2) call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'the only allowable values of stable_option are 1 and 2', FATAL)

if(stable_option == 2 .and. zeta_trans < 0) call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'zeta_trans must be positive', FATAL)

b_stab = 1.0/rich_crit
r_crit = 0.95*rich_crit  ! convergence can get slow if one is 
                         ! close to rich_crit

sqrt_drag_min = 0.0
if(drag_min.ne.0.0) sqrt_drag_min = sqrt(drag_min)

lambda     = 1.0 + (5.0 - b_stab)*zeta_trans   ! used only if stable_option = 2
rich_trans = zeta_trans/(1.0 + 5.0*zeta_trans) ! used only if stable_option = 2

module_is_initialized = .true.

return
end subroutine monin_obukhov_init

!=======================================================================

subroutine monin_obukhov_end

module_is_initialized = .false.

end subroutine monin_obukhov_end

!=======================================================================

subroutine mo_drag_1d &
         (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, &
          u_star, b_star, avail)

real, intent(in)   , dimension(:) :: pt, pt0, z, z0, zt, zq, speed
real, intent(inout), dimension(:) :: drag_m, drag_t, drag_q, u_star, b_star
logical, intent(in), optional, dimension(:) :: avail

logical            :: lavail, avail_dummy(1)
integer            :: n, ier

integer, parameter :: max_iter = 20
real   , parameter :: error=1.e-04, zeta_min=1.e-06, small=1.e-04

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call monin_obukhov_init

n      = size(pt)
lavail = .false.
if(present(avail)) lavail = .true.


if(lavail) then 
   if (count(avail) .eq. 0) return
   call monin_obukhov_drag_1d(grav, vonkarm,               &
        & error, zeta_min, max_iter, small,                         &
        & neutral, stable_option, rich_crit, zeta_trans, drag_min,  &
        & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
        & drag_q, u_star, b_star, lavail, avail, ier)
else
   call monin_obukhov_drag_1d(grav, vonkarm,               &
        & error, zeta_min, max_iter, small,                         &
        & neutral, stable_option, rich_crit, zeta_trans, drag_min,  &
        & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
        & drag_q, u_star, b_star, lavail, avail_dummy, ier)
endif

end subroutine mo_drag_1d


!=======================================================================

subroutine mo_profile_1d(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)

real,    intent(in)                :: zref, zref_t
real,    intent(in) , dimension(:) :: z, z0, zt, zq, u_star, b_star, q_star
real,    intent(out), dimension(:) :: del_m, del_t, del_q
logical, intent(in) , optional, dimension(:) :: avail

logical                            :: dummy_avail(1)
integer                            :: n, ier

! #include "monin_obukhov_interfaces.h"

if(.not. module_is_initialized) call monin_obukhov_init

n = size(z)
if(present(avail)) then

   if (count(avail) .eq. 0) return

   call monin_obukhov_profile_1d(vonkarm, &
        & neutral, stable_option, rich_crit, zeta_trans, &
        & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
        & del_m, del_t, del_q, .true., avail, ier)

else

   call monin_obukhov_profile_1d(vonkarm, &
        & neutral, stable_option, rich_crit, zeta_trans, &
        & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
        & del_m, del_t, del_q, .false., dummy_avail, ier)

endif

end subroutine mo_profile_1d

!=======================================================================

subroutine stable_mix_3d(rich, mix)

real, intent(in) , dimension(:,:,:)  :: rich
real, intent(out), dimension(:,:,:)  :: mix

integer :: n, ier

if(.not. module_is_initialized) call monin_obukhov_init

n = size(rich,1)*size(rich,2)*size(rich,3)
call monin_obukhov_stable_mix(stable_option, rich_crit, zeta_trans, &
     & n, rich, mix, ier)


end subroutine stable_mix_3d

!=======================================================================

subroutine mo_diff_2d_n(z, u_star, b_star, k_m, k_h)

real, intent(in),  dimension(:,:,:) :: z
real, intent(in),  dimension(:,:)   :: u_star, b_star
real, intent(out), dimension(:,:,:) :: k_m, k_h

integer            :: ni, nj, nk, ier
real, parameter    :: ustar_min = 1.e-10

if(.not.module_is_initialized) call monin_obukhov_init

ni = size(z, 1); nj = size(z, 2); nk = size(z, 3)
call monin_obukhov_diff(vonkarm,                           &
          & ustar_min,                                     &
          & neutral, stable_option, rich_crit, zeta_trans, &
          & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

end subroutine mo_diff_2d_n

!=======================================================================
! The following routines are used by the public interfaces above
!=======================================================================

subroutine solve_zeta(rich, z, z0, zt, zq, f_m, f_t, f_q, mask)

real   , intent(in) , dimension(:) :: rich, z, z0, zt, zq
logical, intent(in) , dimension(:) :: mask
real   , intent(out), dimension(:) :: f_m, f_t, f_q


real, parameter    :: error    = 1.e-04
real, parameter    :: zeta_min = 1.e-06
integer, parameter :: max_iter = 20

real    :: max_cor
integer :: iter

real, dimension(size(rich(:))) ::   &
          d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
          ln_z_z0, ln_z_zt, ln_z_zq, zeta,                    &
          phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
          zeta_0, zeta_t, zeta_q, df_m, df_t
          
logical, dimension(size(rich(:))) :: mask_1


z_z0 = z/z0
z_zt = z/zt
z_zq = z/zq
ln_z_z0 = log(z_z0)
ln_z_zt = log(z_zt)
ln_z_zq = log(z_zq)

corr = 0.0
mask_1 = mask

! initial guess

where(mask_1) 
  zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
elsewhere
  zeta = 0.0
end where

where (mask_1 .and. rich >= 0.0)
  zeta = zeta/(1.0 - rich/rich_crit)
end where 

iter_loop: do iter = 1, max_iter

  where (mask_1 .and. abs(zeta).lt.zeta_min) 
    zeta = 0.0
    f_m = ln_z_z0
    f_t = ln_z_zt
    f_q = ln_z_zq
    mask_1 = .false.  ! don't do any more calculations at these pts
  end where
  
  where (mask_1)
    rzeta  = 1.0/zeta
    zeta_0 = zeta/z_z0
    zeta_t = zeta/z_zt
    zeta_q = zeta/z_zq
  elsewhere
    zeta_0 = 0.0
    zeta_t = 0.0
    zeta_q = 0.0
  end where

  call mo_derivative_m(phi_m  , zeta  , mask_1)
  call mo_derivative_m(phi_m_0, zeta_0, mask_1)
  call mo_derivative_t(phi_t  , zeta  , mask_1)
  call mo_derivative_t(phi_t_0, zeta_t, mask_1)
                   
  call mo_integral_m(f_m, zeta, zeta_0, ln_z_z0, mask_1)
  call mo_integral_tq(f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask_1)

  where (mask_1)
    df_m  = (phi_m - phi_m_0)*rzeta
    df_t  = (phi_t - phi_t_0)*rzeta
    rich_1 = zeta*f_t/(f_m*f_m)
    d_rich = rich_1*( rzeta +  df_t/f_t - 2.0 *df_m/f_m) 
    correction = (rich - rich_1)/d_rich  
    corr = min(abs(correction),abs(correction/zeta)) 
      ! the criterion corr < error seems to work ok, but is a bit arbitrary
      !  when zeta is small the tolerance is reduced
  end where
  
  max_cor= maxval(corr)

  if(max_cor > error) then
    mask_1 = mask_1 .and. (corr > error)  
       ! change the mask so computation proceeds only on non-converged points
    where(mask_1) 
      zeta = zeta + correction
    end where
    cycle iter_loop
  else
    return
  end if

end do iter_loop

call error_mesg ('solve_zeta in monin_obukhov_mod',  &
                 'surface drag iteration did not converge', FATAL)

end subroutine solve_zeta

!=======================================================================

subroutine mo_derivative_m(phi_m, zeta, mask)

! the differential similarity function for momentum

real    , intent(out),  dimension(:) :: phi_m
real    , intent(in),   dimension(:) :: zeta
logical , intent(in),   dimension(:) :: mask

logical, dimension(size(zeta(:))) :: stable, unstable
real   , dimension(size(zeta(:))) :: x

stable   = mask .and. zeta >= 0.0
unstable = mask .and. zeta <  0.0

where (unstable) 
  x     = (1 - 16.0*zeta  )**(-0.5)
  phi_m = sqrt(x)  ! phi_m = (1 - 16.0*zeta)**(-0.25)
end where

if(stable_option == 1) then 

  where (stable) 
    phi_m = 1.0 + zeta  *(5.0 + b_stab*zeta)/(1.0 + zeta)
  end where
  
else if(stable_option == 2) then

  where (stable .and. zeta < zeta_trans)
    phi_m = 1 + 5.0*zeta
  end where
  where (stable .and. zeta >= zeta_trans)
    phi_m = lambda + b_stab*zeta
  end where

endif

return
end subroutine mo_derivative_m

!=======================================================================

subroutine mo_derivative_t(phi_t, zeta, mask)

! the differential similarity function for buoyancy and tracers

real    , intent(out),  dimension(:) :: phi_t
real    , intent(in),   dimension(:) :: zeta
logical , intent(in),   dimension(:) :: mask

logical, dimension(size(zeta(:))) :: stable, unstable

stable   = mask .and. zeta >= 0.0
unstable = mask .and. zeta <  0.0

where (unstable) 
  phi_t = (1 - 16.0*zeta)**(-0.5)
end where

if(stable_option == 1) then 

  where (stable) 
    phi_t = 1.0 + zeta*(5.0 + b_stab*zeta)/(1.0 + zeta)
  end where
  
else if(stable_option == 2) then

  where (stable .and. zeta < zeta_trans)
    phi_t = 1 + 5.0*zeta
  end where
  where (stable .and. zeta >= zeta_trans)
    phi_t = lambda + b_stab*zeta
  end where

endif

return
end subroutine mo_derivative_t

!=======================================================================

subroutine mo_integral_tq (psi_t, psi_q, zeta, zeta_t, zeta_q, &
                           ln_z_zt, ln_z_zq, mask)

! the integral similarity function for moisture and tracers

real    , intent(out), dimension(:) :: psi_t, psi_q
real    , intent(in),  dimension(:) :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
logical , intent(in),  dimension(:) :: mask

real, dimension(size(zeta(:))) :: x, x_t, x_q
                               
logical, dimension(size(zeta(:))) :: stable, unstable, &
                                  weakly_stable, strongly_stable

stable   = mask .and. zeta >= 0.0
unstable = mask .and. zeta <  0.0

where(unstable) 

  x     = sqrt(1 - 16.0*zeta)
  x_t   = sqrt(1 - 16.0*zeta_t)
  x_q   = sqrt(1 - 16.0*zeta_q)
  
  psi_t = ln_z_zt - 2.0*log( (1.0 + x)/(1.0 + x_t) )
  psi_q = ln_z_zq - 2.0*log( (1.0 + x)/(1.0 + x_q) )

end where

if( stable_option == 1) then

  where (stable) 
  
    psi_t = ln_z_zt + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_t)) &
       + b_stab*(zeta - zeta_t) 
    psi_q = ln_z_zq + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_q)) &
       + b_stab*(zeta - zeta_q) 
       
  end where
  
else if (stable_option == 2) then

  weakly_stable   = stable .and. zeta <= zeta_trans
  strongly_stable = stable .and. zeta >  zeta_trans

  where (weakly_stable)
    psi_t = ln_z_zt + 5.0*(zeta - zeta_t) 
    psi_q = ln_z_zq + 5.0*(zeta - zeta_q) 
  end where
  
  where(strongly_stable)
    x = (lambda - 1.0)*log(zeta/zeta_trans) + b_stab*(zeta - zeta_trans)
  endwhere
  
  where (strongly_stable .and. zeta_t <= zeta_trans)
    psi_t = ln_z_zt + x + 5.0*(zeta_trans - zeta_t)
  end where
  where (strongly_stable .and. zeta_t > zeta_trans)
    psi_t = lambda*ln_z_zt + b_stab*(zeta  - zeta_t)
  endwhere
  
  where (strongly_stable .and. zeta_q <= zeta_trans)
    psi_q = ln_z_zq + x + 5.0*(zeta_trans - zeta_q)
  end where
  where (strongly_stable .and. zeta_q > zeta_trans)
    psi_q = lambda*ln_z_zq + b_stab*(zeta  - zeta_q)
  endwhere
  
end if

return
end subroutine mo_integral_tq

!=======================================================================

subroutine mo_integral_m (psi_m, zeta, zeta_0, ln_z_z0, mask)

!  the integral similarity function for momentum

real    , intent(out), dimension(:) :: psi_m
real    , intent(in),  dimension(:) :: zeta, zeta_0, ln_z_z0
logical , intent(in),  dimension(:) :: mask

real, dimension(size(zeta(:))) :: x, x_0, x1, x1_0, num, denom, y
                               
logical, dimension(size(zeta(:))) :: stable, unstable, &
                                  weakly_stable, strongly_stable

stable   = mask .and. zeta >= 0.0
unstable = mask .and. zeta <  0.0

where(unstable) 

  x     = sqrt(1 - 16.0*zeta)
  x_0   = sqrt(1 - 16.0*zeta_0)

  x      = sqrt(x)
  x_0    = sqrt(x_0)
  
  x1     = 1.0 + x
  x1_0   = 1.0 + x_0
  
  num    = x1*x1*(1.0 + x*x)
  denom  = x1_0*x1_0*(1.0 + x_0*x_0)
  y      = atan(x) - atan(x_0)
  psi_m  = ln_z_z0 - log(num/denom) + 2*y
  
end where

if( stable_option == 1) then

  where (stable) 
    psi_m = ln_z_z0 + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_0)) &
       + b_stab*(zeta - zeta_0) 
  end where
  
else if (stable_option == 2) then

  weakly_stable   = stable .and. zeta <= zeta_trans
  strongly_stable = stable .and. zeta >  zeta_trans

  where (weakly_stable)
    psi_m = ln_z_z0 + 5.0*(zeta - zeta_0) 
  end where
  
  where(strongly_stable)
    x = (lambda - 1.0)*log(zeta/zeta_trans) + b_stab*(zeta - zeta_trans)
  endwhere
  
  where (strongly_stable .and. zeta_0 <= zeta_trans)
    psi_m = ln_z_z0 + x + 5.0*(zeta_trans - zeta_0)
  end where
  where (strongly_stable .and. zeta_0 > zeta_trans)
    psi_m = lambda*ln_z_z0 + b_stab*(zeta  - zeta_0)
  endwhere
  
end if

return
end subroutine mo_integral_m


!=======================================================================
! The following routines allow the public interfaces to be used
! with different dimensions of the input and output
!
!=======================================================================


subroutine mo_drag_2d &
    (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real, intent(in)   , dimension(:,:) :: z, speed, pt, pt0, z0, zt, zq
real, intent(out)  , dimension(:,:) :: drag_m, drag_t, drag_q
real, intent(inout), dimension(:,:) :: u_star, b_star

integer :: j

do j = 1, size(pt,2)
  call mo_drag_1d (pt(:,j), pt0(:,j), z(:,j), z0(:,j), zt(:,j), zq(:,j), &
                   speed(:,j), drag_m(:,j), drag_t(:,j), drag_q(:,j), &
                   u_star(:,j), b_star(:,j))
end do


return
end subroutine mo_drag_2d

!=======================================================================
subroutine mo_drag_0d &
    (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real, intent(in)    :: z, speed, pt, pt0, z0, zt, zq
real, intent(out)   :: drag_m, drag_t, drag_q, u_star, b_star

real, dimension(1) :: pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
                      drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1

pt_1   (1) = pt
pt0_1  (1) = pt0
z_1    (1) = z
z0_1   (1) = z0
zt_1   (1) = zt
zq_1   (1) = zq
speed_1(1) = speed

call mo_drag_1d (pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
                 drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1)

drag_m = drag_m_1(1)
drag_t = drag_t_1(1)
drag_q = drag_q_1(1)
u_star = u_star_1(1)
b_star = b_star_1(1)

return
end subroutine mo_drag_0d
!=======================================================================

subroutine mo_profile_2d(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real, intent(in)                  :: zref, zref_t
real, intent(in) , dimension(:,:) :: z, z0, zt, zq, u_star, b_star, q_star
real, intent(out), dimension(:,:) :: del_m, del_h, del_q

integer :: j

do j = 1, size(z,2)
  call mo_profile_1d (zref, zref_t, z(:,j), z0(:,j), zt(:,j),         &
                      zq(:,j), u_star(:,j), b_star(:,j), q_star(:,j), &
                      del_m(:,j), del_h (:,j), del_q (:,j))
enddo

return
end subroutine mo_profile_2d

!=======================================================================

subroutine mo_profile_0d(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real, intent(in)  :: zref, zref_t
real, intent(in)  :: z, z0, zt, zq, u_star, b_star, q_star
real, intent(out) :: del_m, del_h, del_q

real, dimension(1) :: z_1, z0_1, zt_1, zq_1, u_star_1, b_star_1, q_star_1, &
                      del_m_1, del_h_1, del_q_1

z_1     (1) = z
z0_1    (1) = z0
zt_1    (1) = zt
zq_1    (1) = zq
u_star_1(1) = u_star
b_star_1(1) = b_star
q_star_1(1) = q_star

call mo_profile_1d (zref, zref_t, z_1, z0_1, zt_1, zq_1, &
                    u_star_1, b_star_1, q_star_1,        &
                    del_m_1, del_h_1, del_q_1)
                    
del_m = del_m_1(1)
del_h = del_h_1(1)
del_q = del_q_1(1)
                    

return
end subroutine mo_profile_0d

!=======================================================================

subroutine mo_profile_1d_n(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)

real,    intent(in),  dimension(:)   :: zref
real,    intent(in) , dimension(:)   :: z, z0, zt, zq, u_star, b_star, q_star
real,    intent(out), dimension(:,:) :: del_m, del_t, del_q
logical, intent(in) , optional, dimension(:) :: avail

integer :: k

do k = 1, size(zref(:))
  if(present(avail)) then
    call mo_profile_1d (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k), avail)
  else 
      call mo_profile_1d (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k))
  endif
enddo

return
end subroutine mo_profile_1d_n

!=======================================================================

subroutine mo_profile_0d_n(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q)

real,    intent(in),  dimension(:) :: zref
real,    intent(in)                :: z, z0, zt, zq, u_star, b_star, q_star
real,    intent(out), dimension(:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
  call mo_profile_0d (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(k), del_t(k), del_q(k))
enddo

return
end subroutine mo_profile_0d_n

!=======================================================================

subroutine mo_profile_2d_n(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q)

real,    intent(in),  dimension(:)     :: zref
real,    intent(in),  dimension(:,:)   :: z, z0, zt, zq, u_star, b_star, q_star
real,    intent(out), dimension(:,:,:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
  call mo_profile_2d (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,:,k), del_t(:,:,k), del_q(:,:,k))
enddo

return
end subroutine mo_profile_2d_n

!=======================================================================

subroutine mo_diff_2d_1(z, u_star, b_star, k_m, k_h)

real, intent(in),  dimension(:,:) :: z, u_star, b_star
real, intent(out), dimension(:,:) :: k_m, k_h

real   , dimension(size(z,1),size(z,2),1) :: z_n, k_m_n, k_h_n

z_n(:,:,1) = z

call mo_diff_2d_n(z_n, u_star, b_star, k_m_n, k_h_n)

k_m = k_m_n(:,:,1)
k_h = k_h_n(:,:,1)

return
end subroutine mo_diff_2d_1


!=======================================================================

subroutine mo_diff_1d_1(z, u_star, b_star, k_m, k_h)

real, intent(in),  dimension(:) :: z, u_star, b_star
real, intent(out), dimension(:) :: k_m, k_h

real, dimension(size(z),1,1) :: z_n, k_m_n, k_h_n
real, dimension(size(z),1)   :: u_star_n, b_star_n

z_n   (:,1,1) = z
u_star_n(:,1) = u_star
b_star_n(:,1) = b_star

call mo_diff_2d_n(z_n, u_star_n, b_star_n, k_m_n, k_h_n)

k_m = k_m_n(:,1,1)
k_h = k_h_n(:,1,1)

return
end subroutine mo_diff_1d_1

!=======================================================================

subroutine mo_diff_1d_n(z, u_star, b_star, k_m, k_h)

real, intent(in),  dimension(:,:) :: z
real, intent(in),  dimension(:)   :: u_star, b_star
real, intent(out), dimension(:,:) :: k_m, k_h

real, dimension(size(z,1),1)            :: u_star2, b_star2
real, dimension(size(z,1),1, size(z,2)) :: z2, k_m2, k_h2

integer :: n

do n = 1, size(z,2)
  z2   (:,1,n) = z(:,n)
enddo
u_star2(:,1) = u_star
b_star2(:,1) = b_star

call mo_diff_2d_n(z2, u_star2, b_star2, k_m2, k_h2)

do n = 1, size(z,2)
  k_m(:,n) = k_m2(:,1,n)
  k_h(:,n) = k_h2(:,1,n)
enddo

return
end subroutine mo_diff_1d_n

!=======================================================================

subroutine mo_diff_0d_1(z, u_star, b_star, k_m, k_h)

real, intent(in)  :: z, u_star, b_star
real, intent(out) :: k_m, k_h

integer            :: ni, nj, nk, ier
real, parameter    :: ustar_min = 1.e-10

if(.not.module_is_initialized) call monin_obukhov_init

ni = 1; nj = 1; nk = 1
call monin_obukhov_diff(vonkarm,                           &
          & ustar_min,                                     &
          & neutral, stable_option, rich_crit, zeta_trans, &
          & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

end subroutine mo_diff_0d_1

!=======================================================================

subroutine mo_diff_0d_n(z, u_star, b_star, k_m, k_h)

real, intent(in),  dimension(:) :: z
real, intent(in)                :: u_star, b_star
real, intent(out), dimension(:) :: k_m, k_h

integer            :: ni, nj, nk, ier
real, parameter    :: ustar_min = 1.e-10

if(.not.module_is_initialized) call monin_obukhov_init

ni = 1; nj = 1; nk = size(z(:))
call monin_obukhov_diff(vonkarm,                           &
          & ustar_min,                                     &
          & neutral, stable_option, rich_crit, zeta_trans, &
          & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

end subroutine mo_diff_0d_n

!=======================================================================

subroutine stable_mix_2d(rich, mix)

real, intent(in) , dimension(:,:)  :: rich
real, intent(out), dimension(:,:)  :: mix

real, dimension(size(rich,1),size(rich,2),1) :: rich_3d, mix_3d

rich_3d(:,:,1) = rich

call stable_mix_3d(rich_3d, mix_3d)

mix = mix_3d(:,:,1)

return
end subroutine stable_mix_2d


!=======================================================================

subroutine stable_mix_1d(rich, mix)

real, intent(in) , dimension(:)  :: rich
real, intent(out), dimension(:)  :: mix

real, dimension(size(rich),1,1) :: rich_3d, mix_3d

rich_3d(:,1,1) = rich

call stable_mix_3d(rich_3d, mix_3d)

mix = mix_3d(:,1,1)

return
end subroutine stable_mix_1d

!=======================================================================

subroutine stable_mix_0d(rich, mix)

real, intent(in) :: rich
real, intent(out) :: mix

real, dimension(1,1,1) :: rich_3d, mix_3d

rich_3d(1,1,1) = rich

call stable_mix_3d(rich_3d, mix_3d)

mix = mix_3d(1,1,1)

return
end subroutine stable_mix_0d
!=======================================================================

end module monin_obukhov_mod

