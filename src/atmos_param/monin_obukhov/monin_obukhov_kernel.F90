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

#include <fms_platform.h>

! -*-F90-*-
! $Id: monin_obukhov_kernel.F90,v 13.0 2006/03/28 21:10:33 fms Exp $

!==============================================================================
! Kernel routine interface
!==============================================================================

module monin_obukhov_inter

! explicit interface to all kernel routines
#include "monin_obukhov_interfaces.h"

end module monin_obukhov_inter

!==============================================================================
! Kernel routines
!==============================================================================

_PURE subroutine monin_obukhov_diff(vonkarm,                &
     & ustar_min,                                     &
     & neutral, stable_option, rich_crit, zeta_trans, &
     & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

  implicit none

  real   , intent(in   )                        :: vonkarm
  real   , intent(in   )                        :: ustar_min ! = 1.e-10
  logical, intent(in   )                        :: neutral
  integer, intent(in   )                        :: stable_option
  real   , intent(in   )                        :: rich_crit, zeta_trans
  integer, intent(in   )                        :: ni, nj, nk
  real   , intent(in   ), dimension(ni, nj, nk) :: z
  real   , intent(in   ), dimension(ni, nj)     :: u_star, b_star
  real   , intent(  out), dimension(ni, nj, nk) :: k_m, k_h
  integer, intent(  out)                        :: ier

  real , dimension(ni, nj) :: phi_m, phi_h, zeta, uss
  integer :: j, k

  logical, dimension(ni) :: mask

  interface
     _PURE subroutine monin_obukhov_derivative_t(stable_option, rich_crit, zeta_trans, &
          & n, phi_t, zeta, mask, ier)

       ! the differential similarity function for buoyancy and tracers
       ! Note: seems to be the same as monin_obukhov_derivative_m?

       integer, intent(in   )                :: stable_option
       real   , intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(  out), dimension(n)  :: phi_t
       real   , intent(in   ), dimension(n)  :: zeta
       logical, intent(in   ), dimension(n)  :: mask  
       integer, intent(  out)                :: ier
     end subroutine monin_obukhov_derivative_t

     _PURE subroutine monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
          & n, phi_m, zeta, mask, ier)

       ! the differential similarity function for momentum

       integer, intent(in   )                :: stable_option
       real   , intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(  out), dimension(n)  :: phi_m
       real   , intent(in   ), dimension(n)  :: zeta
       logical, intent(in   ), dimension(n)  :: mask
       integer, intent(out  )                :: ier
     end subroutine monin_obukhov_derivative_m
  end interface

  ier = 0

  mask = .true.
  uss = max(u_star, ustar_min)

  if(neutral) then
     do k = 1, size(z,3)
        k_m(:,:,k) = vonkarm *uss*z(:,:,k)
        k_h(:,:,k) = k_m(:,:,k)
     end do
  else
     do k = 1, size(z,3)
        zeta = - vonkarm * b_star*z(:,:,k)/(uss*uss)
        do j = 1, size(z,2)
           call monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
                & ni, phi_m(:,j), zeta(:,j), mask, ier)
           call monin_obukhov_derivative_t(stable_option, rich_crit, zeta_trans, &
                & ni, phi_h(:,j), zeta(:,j), mask, ier)
        enddo
        k_m(:,:,k) = vonkarm * uss*z(:,:,k)/phi_m
        k_h(:,:,k) = vonkarm * uss*z(:,:,k)/phi_h
     end do
  endif

end subroutine monin_obukhov_diff
!==============================================================================
_PURE subroutine monin_obukhov_drag_1d(grav, vonkarm,               &
     & error, zeta_min, max_iter, small,                         &
     & neutral, stable_option, rich_crit, zeta_trans, drag_min,  &
     & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
     & drag_q, u_star, b_star, lavail, avail, ier)

  implicit none

  real   , intent(in   )                :: grav     
  real   , intent(in   )                :: vonkarm   
  real   , intent(in   )                :: error    ! = 1.e-04
  real   , intent(in   )                :: zeta_min ! = 1.e-06
  integer, intent(in   )                :: max_iter ! = 20
  real   , intent(in   )                :: small    ! = 1.e-04
  logical, intent(in   )                :: neutral
  integer, intent(in   )                :: stable_option
  real   , intent(in   )                :: rich_crit, zeta_trans, drag_min
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: pt, pt0, z, z0, zt, zq, speed
  real   , intent(inout), dimension(n)  :: drag_m, drag_t, drag_q, u_star, b_star
  logical, intent(in   )                :: lavail ! whether to use provided mask or not
  logical, intent(in   ), dimension(n)  :: avail  ! provided mask 
  integer, intent(out  )                :: ier

  real   , dimension(n) :: rich, fm, ft, fq, zz
  logical, dimension(n) :: mask, mask_1, mask_2
  real   , dimension(n) :: delta_b !!, us, bs, qs
  real                  :: r_crit, sqrt_drag_min
  real                  :: us, bs, qs
  integer               :: i

  interface
     _PURE subroutine monin_obukhov_solve_zeta(error, zeta_min, max_iter, small,  &
          & stable_option, rich_crit, zeta_trans,                           &
          & n, rich, z, z0, zt, zq, f_m, f_t, f_q, mask, ier)

       real   , intent(in   )                :: error    ! = 1.e-04
       real   , intent(in   )                :: zeta_min ! = 1.e-06
       integer, intent(in   )                :: max_iter ! = 20
       real   , intent(in   )                :: small    ! = 1.e-04
       integer, intent(in   )                :: stable_option
       real   , intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(in   ), dimension(n)  :: rich, z, z0, zt, zq
       logical, intent(in   ), dimension(n)  :: mask
       real   , intent(  out), dimension(n)  :: f_m, f_t, f_q
       integer, intent(  out)                :: ier
     end subroutine monin_obukhov_solve_zeta
  end interface

  ier = 0
  r_crit = 0.95*rich_crit  ! convergence can get slow if one is 
                           ! close to rich_crit
  sqrt_drag_min = 0.0
  if(drag_min.ne.0.0) sqrt_drag_min = sqrt(drag_min)

  mask = .true.
  if(lavail) mask = avail

  where(mask) 
     delta_b = grav*(pt0 - pt)/pt0
     rich    = - z*delta_b/(speed*speed + small)
     zz      = max(z,z0,zt,zq)
  elsewhere 
     rich = 0.0
  end where

  if(neutral) then

     do i = 1, n
        if(mask(i)) then
           fm(i)   = log(zz(i)/z0(i))
           ft(i)   = log(zz(i)/zt(i))
           fq(i)   = log(zz(i)/zq(i))
           us   = vonkarm/fm(i)
           bs   = vonkarm/ft(i)
           qs   = vonkarm/fq(i)
           drag_m(i)    = us*us
           drag_t(i)    = us*bs
           drag_q(i)    = us*qs
           u_star(i) = us*speed(i)
           b_star(i) = bs*delta_b(i)
        end if
     enddo

  else

     mask_1 = mask .and. rich <  r_crit
     mask_2 = mask .and. rich >= r_crit

     do i = 1, n
        if(mask_2(i)) then
           drag_m(i)   = drag_min
           drag_t(i)   = drag_min
           drag_q(i)   = drag_min
           us       = sqrt_drag_min
           bs       = sqrt_drag_min
           u_star(i)   = us*speed(i)
           b_star(i)   = bs*delta_b(i)
        end if
     enddo

     call monin_obukhov_solve_zeta (error, zeta_min, max_iter, small, &
          & stable_option, rich_crit, zeta_trans,                     &
          & n, rich, zz, z0, zt, zq, fm, ft, fq, mask_1, ier)

     do i = 1, n
        if(mask_1(i)) then
           us   = max(vonkarm/fm(i), sqrt_drag_min)
           bs   = max(vonkarm/ft(i), sqrt_drag_min)
           qs   = max(vonkarm/fq(i), sqrt_drag_min)
           drag_m(i)   = us*us
           drag_t(i)   = us*bs
           drag_q(i)   = us*qs
           u_star(i)   = us*speed(i)
           b_star(i)   = bs*delta_b(i)
        endif
     enddo

  end if

end subroutine monin_obukhov_drag_1d
!==============================================================================
_PURE subroutine monin_obukhov_solve_zeta(error, zeta_min, max_iter, small,  &
     & stable_option, rich_crit, zeta_trans,                           &
     & n, rich, z, z0, zt, zq, f_m, f_t, f_q, mask, ier)

  implicit none

  real   , intent(in   )                :: error    ! = 1.e-04
  real   , intent(in   )                :: zeta_min ! = 1.e-06
  integer, intent(in   )                :: max_iter ! = 20
  real   , intent(in   )                :: small    ! = 1.e-04
  integer, intent(in   )                :: stable_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: rich, z, z0, zt, zq
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(  out), dimension(n)  :: f_m, f_t, f_q
  integer, intent(  out)                :: ier


  real    :: max_cor
  integer :: iter

  real, dimension(n) ::   &
       d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
       ln_z_z0, ln_z_zt, ln_z_zq, zeta,                    &
       phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
       zeta_0, zeta_t, zeta_q, df_m, df_t

  logical, dimension(n) :: mask_1

  interface
     _PURE subroutine monin_obukhov_derivative_t(stable_option, rich_crit, zeta_trans, &
          & n, phi_t, zeta, mask, ier)

       ! the differential similarity function for buoyancy and tracers
       ! Note: seems to be the same as monin_obukhov_derivative_m?

       integer, intent(in   )                :: stable_option
       real   , intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(  out), dimension(n)  :: phi_t
       real   , intent(in   ), dimension(n)  :: zeta
       logical, intent(in   ), dimension(n)  :: mask  
       integer, intent(  out)                :: ier
     end subroutine monin_obukhov_derivative_t
     _PURE subroutine monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
          & n, phi_m, zeta, mask, ier)

       ! the differential similarity function for momentum

       integer, intent(in   )                :: stable_option
       real   , intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(  out), dimension(n)  :: phi_m
       real   , intent(in   ), dimension(n)  :: zeta
       logical, intent(in   ), dimension(n)  :: mask
       integer, intent(out  )                :: ier
     end subroutine monin_obukhov_derivative_m
     _PURE subroutine monin_obukhov_integral_tq(stable_option, rich_crit, zeta_trans, &
          & n, psi_t, psi_q, zeta, zeta_t, zeta_q, &
          & ln_z_zt, ln_z_zq, mask, ier)

       ! the integral similarity function for moisture and tracers

       integer, intent(in   )                :: stable_option
       real,    intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(  out), dimension(n)  :: psi_t, psi_q
       real   , intent(in)   , dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
       logical, intent(in)   , dimension(n)  :: mask
       integer, intent(  out)                :: ier
     end subroutine monin_obukhov_integral_tq

     _PURE subroutine monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
          & n, psi_m, zeta, zeta_0, ln_z_z0, mask, ier)

       !  the integral similarity function for momentum

       integer, intent(in   )                :: stable_option
       real   , intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(  out), dimension(n)  :: psi_m
       real   , intent(in)   , dimension(n)  :: zeta, zeta_0, ln_z_z0
       logical, intent(in)   , dimension(n)  :: mask
       integer, intent(out)                  :: ier
     end subroutine monin_obukhov_integral_m
     
  end interface


  ier = 0

  z_z0 = z/z0
  z_zt = z/zt
  z_zq = z/zq
  ln_z_z0 = log(z_z0)
  ln_z_zt = log(z_zt)
  ln_z_zq = log(z_zq)

  corr = 0.0
  mask_1 = mask

  ! initial guess

  zeta = 0.0
  where(mask_1) 
     zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
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

     
     zeta_0 = 0.0
     zeta_t = 0.0
     zeta_q = 0.0
     where (mask_1)
        rzeta  = 1.0/zeta
        zeta_0 = zeta/z_z0
        zeta_t = zeta/z_zt
        zeta_q = zeta/z_zq
     end where

     call monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
          & n, phi_m  , zeta  , mask_1, ier)
     call monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
          & n, phi_m_0, zeta_0,  mask_1, ier)
     call monin_obukhov_derivative_t(stable_option, rich_crit, zeta_trans, &
          & n, phi_t  , zeta  , mask_1, ier)
     call monin_obukhov_derivative_t(stable_option, rich_crit, zeta_trans, &
          & n, phi_t_0, zeta_t, mask_1, ier)

     call monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
          & n, f_m, zeta, zeta_0, ln_z_z0, mask_1, ier)
     call monin_obukhov_integral_tq(stable_option, rich_crit, zeta_trans, &
          & n, f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask_1, ier)

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

  ier = 1 ! surface drag iteration did not converge

end subroutine monin_obukhov_solve_zeta
!==============================================================================
_PURE subroutine monin_obukhov_derivative_t(stable_option, rich_crit, zeta_trans, &
     & n, phi_t, zeta, mask, ier)

  ! the differential similarity function for buoyancy and tracers
  ! Note: seems to be the same as monin_obukhov_derivative_m?

  implicit none

  integer, intent(in   )                :: stable_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(  out), dimension(n)  :: phi_t
  real   , intent(in   ), dimension(n)  :: zeta
  logical, intent(in   ), dimension(n)  :: mask  
  integer, intent(  out)                :: ier

  logical, dimension(n) :: stable, unstable
  real                  :: b_stab, lambda

  ier = 0
  b_stab     = 1.0/rich_crit

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

     lambda = 1.0 + (5.0 - b_stab)*zeta_trans

     where (stable .and. zeta < zeta_trans)
        phi_t = 1 + 5.0*zeta
     end where
     where (stable .and. zeta >= zeta_trans)
        phi_t = lambda + b_stab*zeta
     end where

  endif

end subroutine monin_obukhov_derivative_t
!==============================================================================
_PURE subroutine monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
     & n, phi_m, zeta, mask, ier)

  ! the differential similarity function for momentum

  implicit none

  integer, intent(in   )                :: stable_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(  out), dimension(n)  :: phi_m
  real   , intent(in   ), dimension(n)  :: zeta
  logical, intent(in   ), dimension(n)  :: mask
  integer, intent(out  )                :: ier

  logical, dimension(n) :: stable, unstable
  real   , dimension(n) :: x
  real                  :: b_stab, lambda


  ier = 0
  b_stab     = 1.0/rich_crit

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

     lambda = 1.0 + (5.0 - b_stab)*zeta_trans

     where (stable .and. zeta < zeta_trans)
        phi_m = 1 + 5.0*zeta
     end where
     where (stable .and. zeta >= zeta_trans)
        phi_m = lambda + b_stab*zeta
     end where

  endif

end subroutine monin_obukhov_derivative_m
!==============================================================================
_PURE subroutine monin_obukhov_profile_1d(vonkarm, &
     & neutral, stable_option, rich_crit, zeta_trans, &
     & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
     & del_m, del_t, del_q, lavail, avail, ier)

  implicit none

  real   , intent(in   )                :: vonkarm
  logical, intent(in   )                :: neutral
  integer, intent(in   )                :: stable_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real,    intent(in   )                :: zref, zref_t
  real,    intent(in   ), dimension(n)  :: z, z0, zt, zq, u_star, b_star, q_star
  real,    intent(  out), dimension(n)  :: del_m, del_t, del_q
  logical, intent(in   )                :: lavail ! whether to use provided mask or not
  logical, intent(in   ), dimension(n)  :: avail  ! provided mask
  integer, intent(out  )                :: ier

  real, dimension(n) :: zeta, zeta_0, zeta_t, zeta_q, zeta_ref, zeta_ref_t, &
       ln_z_z0, ln_z_zt, ln_z_zq, ln_z_zref, ln_z_zref_t,  &
       f_m_ref, f_m, f_t_ref, f_t, f_q_ref, f_q,           &
       mo_length_inv

  logical, dimension(n) :: mask

  interface
     _PURE subroutine monin_obukhov_integral_tq(stable_option, rich_crit, zeta_trans, &
          & n, psi_t, psi_q, zeta, zeta_t, zeta_q, &
          & ln_z_zt, ln_z_zq, mask, ier)

       ! the integral similarity function for moisture and tracers

       integer, intent(in   )                :: stable_option
       real,    intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(  out), dimension(n)  :: psi_t, psi_q
       real   , intent(in)   , dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
       logical, intent(in)   , dimension(n)  :: mask
       integer, intent(  out)                :: ier
     end subroutine monin_obukhov_integral_tq

     _PURE subroutine monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
          & n, psi_m, zeta, zeta_0, ln_z_z0, mask, ier)

       !  the integral similarity function for momentum

       integer, intent(in   )                :: stable_option
       real   , intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real   , intent(  out), dimension(n)  :: psi_m
       real   , intent(in)   , dimension(n)  :: zeta, zeta_0, ln_z_z0
       logical, intent(in)   , dimension(n)  :: mask
       integer, intent(out)                  :: ier
     end subroutine monin_obukhov_integral_m
  end interface

  ier = 0

  mask = .true.
  if(lavail) mask = avail

  del_m = 0.0  ! zero output arrays
  del_t = 0.0
  del_q = 0.0

  where(mask) 
     ln_z_z0     = log(z/z0)
     ln_z_zt     = log(z/zt)
     ln_z_zq     = log(z/zq)
     ln_z_zref   = log(z/zref)
     ln_z_zref_t = log(z/zref_t)
  endwhere

  if(neutral) then

     where(mask)
        del_m = 1.0 - ln_z_zref  /ln_z_z0
        del_t = 1.0 - ln_z_zref_t/ln_z_zt
        del_q = 1.0 - ln_z_zref_t/ln_z_zq
     endwhere

  else

     mo_length_inv = 0.0
     zeta          = 0.0
     zeta_0        = 0.0
     zeta_t        = 0.0
     zeta_q        = 0.0
     zeta_ref      = 0.0
     zeta_ref_t    = 0.0

     where(mask .and. u_star > 0.0) 
        mo_length_inv = - vonkarm * b_star/(u_star*u_star)
        zeta       = z     *mo_length_inv
        zeta_0     = z0    *mo_length_inv
        zeta_t     = zt    *mo_length_inv
        zeta_q     = zq    *mo_length_inv
        zeta_ref   = zref  *mo_length_inv
        zeta_ref_t = zref_t*mo_length_inv
     endwhere

     call monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
          & n, f_m,     zeta, zeta_0,   ln_z_z0,   mask, ier)
     call monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
          & n, f_m_ref, zeta, zeta_ref, ln_z_zref, mask, ier)

     call monin_obukhov_integral_tq(stable_option, rich_crit, zeta_trans, &
          & n, f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask, ier)
     call monin_obukhov_integral_tq(stable_option, rich_crit, zeta_trans, &
          & n, f_t_ref, f_q_ref, zeta, zeta_ref_t, zeta_ref_t, ln_z_zref_t, ln_z_zref_t,  mask, ier)

     where(mask)
        del_m = 1.0 - f_m_ref/f_m
        del_t = 1.0 - f_t_ref/f_t
        del_q = 1.0 - f_q_ref/f_q
     endwhere

  end if


end subroutine monin_obukhov_profile_1d
!==============================================================================
_PURE subroutine monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
     & n, psi_m, zeta, zeta_0, ln_z_z0, mask, ier)

  !  the integral similarity function for momentum

  implicit none

  integer, intent(in   )                :: stable_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(  out), dimension(n)  :: psi_m
  real   , intent(in)   , dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in)   , dimension(n)  :: mask
  integer, intent(out)                  :: ier

  real                   :: b_stab, lambda

  real, dimension(n) :: x, x_0, x1, x1_0, num, denom, y
  logical, dimension(n) :: stable, unstable, &
       weakly_stable, strongly_stable

  ier = 0

  b_stab     = 1.0/rich_crit

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

     lambda = 1.0 + (5.0 - b_stab)*zeta_trans

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

end subroutine monin_obukhov_integral_m
!==============================================================================
_PURE subroutine monin_obukhov_integral_tq(stable_option, rich_crit, zeta_trans, &
     & n, psi_t, psi_q, zeta, zeta_t, zeta_q, &
     & ln_z_zt, ln_z_zq, mask, ier)

  ! the integral similarity function for moisture and tracers

  implicit none

  integer, intent(in   )                :: stable_option
  real,    intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(  out), dimension(n)  :: psi_t, psi_q
  real   , intent(in)   , dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
  logical, intent(in)   , dimension(n)  :: mask
  integer, intent(  out)                :: ier
  
  real, dimension(n)     :: x, x_t, x_q                              
  logical, dimension(n)  :: stable, unstable, &
                                  weakly_stable, strongly_stable
  real                   :: b_stab, lambda

  ier = 0
  
  b_stab     = 1.0/rich_crit

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

   lambda = 1.0 + (5.0 - b_stab)*zeta_trans

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

end subroutine monin_obukhov_integral_tq
!==============================================================================
_PURE subroutine monin_obukhov_stable_mix(stable_option, rich_crit, zeta_trans, &
     &                              n, rich, mix, ier)

  implicit none

  integer, intent(in   )                 :: stable_option
  real   , intent(in   )                 :: rich_crit, zeta_trans
  integer, intent(in   )                 :: n
  real   , intent(in   ), dimension(n)   :: rich
  real   , intent(  out), dimension(n)   :: mix  
  integer, intent(  out)                 :: ier

  real               :: r, a, b, c, zeta, phi
  real               :: b_stab, rich_trans, lambda
  integer            :: i

  ier = 0

mix = 0.0
b_stab     = 1.0/rich_crit
rich_trans = zeta_trans/(1.0 + 5.0*zeta_trans)

if(stable_option == 1) then

     c = - 1.0
     do i = 1, n
        if(rich(i) > 0.0 .and. rich(i) < rich_crit) then
           r = 1.0/rich(i)
           a = r - b_stab
           b = r - (1.0 + 5.0)
           zeta = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a)
           phi = 1.0 + b_stab*zeta + (5.0 - b_stab)*zeta/(1.0 + zeta)
           mix(i) = 1./(phi*phi)
     endif
  end do
  
else if(stable_option == 2) then
 
  lambda = 1.0 + (5.0 - b_stab)*zeta_trans

  where(rich > 0.0 .and. rich <= rich_trans)
    mix = (1.0 - 5.0*rich)**2
  end where
  where(rich > rich_trans .and. rich < rich_crit)
    mix = ((1.0 - b_stab*rich)/lambda)**2
  end where
  
end if

end subroutine monin_obukhov_stable_mix


#ifdef _TEST_MONIN_OBUKHOV
!==============================================================================
! Unit test
!==============================================================================

program test

  use monin_obukhov_inter

  implicit none
  integer, parameter :: i8 = selected_int_kind(18)
  integer(i8)        :: ier_tot, ier

  real    :: grav, vonkarm, error, zeta_min, small, ustar_min
  real    :: zref, zref_t
  integer :: max_iter

  real    :: rich_crit, zeta_trans, drag_min
  logical :: neutral
  integer :: stable_option

  grav          = 9.80
  vonkarm       = 0.4
  error         = 1.0e-4
  zeta_min      = 1.0e-6
  max_iter      = 20
  small         = 1.0e-4
  neutral       = .false.
  stable_option = 1
  rich_crit     =10.0
  zeta_trans    = 0.5
  drag_min      = 1.0e-5
  ustar_min     = 1.e-10
  
  zref   = 10.
  zref_t = 2.
  

  ier_tot = 0
  call test_drag
  print *,'test_drag                    ier = ', ier
  ier_tot = ier_tot + ier

  call test_stable_mix
  print *,'test_stable_mix              ier = ', ier
  ier_tot = ier_tot + ier

  call test_diff
  print *,'test_diff                    ier = ', ier
  ier_tot = ier_tot + ier

  call test_profile
  print *,'test_profile                 ier = ', ier
  ier_tot = ier_tot + ier

  if(ier_tot/=0) then
     print *, ier_tot, '***ERRORS detected***'
  else
     print *,'No error detected.'
  endif

  CONTAINS

    subroutine test_drag

      integer(i8)        :: w

      integer :: i, ier_l
      integer, parameter :: n = 5
      logical :: avail(n), lavail

      real, dimension(n) :: pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star

      ! potential temperature
      pt     = (/ 268.559120403867, 269.799228886728, 277.443023238556, 295.79192777341, 293.268717243262 /)
      pt0    = (/ 273.42369841804 , 272.551410044203, 278.638168565727, 298.133068766049, 292.898163706587/)
      z      = (/ 29.432779269303, 30.0497139076724, 31.6880000418153, 34.1873479240475, 33.2184943356517/)
      z0     = (/ 5.86144925739178e-05, 0.0001, 0.000641655193293549, 3.23383768877187e-05, 0.07/)
      zt     = (/ 3.69403636275411e-05, 0.0001, 1.01735489109205e-05, 7.63933834969505e-05, 0.00947346982656289/)
      zq     = (/ 5.72575636226887e-05, 0.0001, 5.72575636226887e-05, 5.72575636226887e-05, 5.72575636226887e-05/)
      speed  = (/ 2.9693638452068, 2.43308757772094, 5.69418282305367, 9.5608693754561, 4.35302260074334/)
      lavail = .true.
      avail  = (/.true., .true., .true., .true., .true. /)

      drag_m = 0
      drag_t = 0
      drag_q = 0
      u_star = 0
      b_star = 0

      call monin_obukhov_drag_1d(grav, vonkarm,               &
           & error, zeta_min, max_iter, small,                         &
           & neutral, stable_option, rich_crit, zeta_trans, drag_min,  &
           & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
           & drag_q, u_star, b_star, lavail, avail, ier_l)

      ! check sum results
      w = 0
      w = w + transfer(sum(drag_m), w)
      w = w + transfer(sum(drag_t), w)
      w = w + transfer(sum(drag_q), w)
      w = w + transfer(sum(u_star), w)
      w = w + transfer(sum(b_star), w)

      ! plug in check sum here>>>
#if defined(__INTEL_COMPILER) || defined(_LF95)
#define CHKSUM_DRAG 4466746452959549648
#endif
#if defined(_PGF95)
#define CHKSUM_DRAG 4466746452959549650
#endif


      print *,'chksum test_drag      : ', w, ' ref ', CHKSUM_DRAG
      ier = CHKSUM_DRAG - w

    end subroutine test_drag

    subroutine test_stable_mix

      integer(i8)        :: w

      integer, parameter      :: n = 5
      real   , dimension(n)   :: rich
      real   , dimension(n)   :: mix
      integer                 :: ier_l


      stable_option = 1
      rich_crit     = 10.0
      zeta_trans    =  0.5

      rich = (/1650.92431853365, 1650.9256285137, 77.7636819036559, 1.92806556391324, 0.414767442012442/)


      call monin_obukhov_stable_mix(stable_option, rich_crit, zeta_trans, &
           &                              n, rich, mix, ier_l)

      w = transfer( sum(mix) , w)

      ! plug in check sum here>>>
#if defined(__INTEL_COMPILER) || defined(_LF95)
#define CHKSUM_STABLE_MIX 4590035772608644256
#endif
#if defined(_PGF95)
#define CHKSUM_STABLE_MIX 4590035772608644258
#endif

      print *,'chksum test_stable_mix: ', w, ' ref ', CHKSUM_STABLE_MIX
      ier = CHKSUM_STABLE_MIX - w

    end subroutine test_stable_mix

    !========================================================================

    subroutine test_diff

      integer(i8)        :: w

      integer, parameter             :: ni=1, nj=1, nk=1
      real   , dimension(ni, nj, nk) :: z
      real   , dimension(ni, nj)     :: u_star, b_star
      real   , dimension(ni, nj, nk) :: k_m, k_h
      integer                        :: ier_l

      z      = 19.9982554527751
      u_star = 0.129638955971075
      b_star = 0.000991799765557209

      call monin_obukhov_diff(vonkarm,                        &
           & ustar_min,                                     &
           & neutral, stable_option, rich_crit, zeta_trans, &
           & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier_l)

      w = 0
      w = w + transfer( sum(k_m) , w)
      w = w + transfer( sum(k_h) , w)

      ! plug check sum in here>>>
#if defined(__INTEL_COMPILER) || defined(_LF95) || defined(_PGF95)
#define CHKSUM_DIFF -9222066590093362639
#endif

      print *,'chksum test_diff      : ', w, ' ref ', CHKSUM_DIFF

      ier = CHKSUM_DIFF - w  

    end subroutine test_diff

    !========================================================================

    subroutine test_profile

      integer(i8)        :: w

      integer, parameter :: n = 5
      integer            :: ier_l

      logical :: avail(n)

      real, dimension(n) :: z, z0, zt, zq, u_star, b_star, q_star
      real, dimension(n) :: del_m, del_t, del_q

      z      = (/ 29.432779269303, 30.0497139076724, 31.6880000418153, 34.1873479240475, 33.2184943356517 /)
      z0     = (/ 5.86144925739178e-05, 0.0001, 0.000641655193293549, 3.23383768877187e-05, 0.07/)
      zt     = (/ 3.69403636275411e-05, 0.0001, 1.01735489109205e-05, 7.63933834969505e-05, 0.00947346982656289/)
      zq     = (/ 5.72575636226887e-05, 0.0001, 5.72575636226887e-05, 5.72575636226887e-05, 5.72575636226887e-05/)
      u_star = (/ 0.109462510724615, 0.0932942802513508, 0.223232887323184, 0.290918439028557, 0.260087579361467/)
      b_star = (/ 0.00690834676781433, 0.00428178089592372, 0.00121229800895103, 0.00262353784027441, -0.000570314880866852/)
      q_star = (/ 0.000110861442197537, 9.44983279664197e-05, 4.17643828631936e-05, 0.000133135421415819, 9.36317815993945e-06/) 

      avail = (/ .true., .true.,.true.,.true.,.true. /)

      call monin_obukhov_profile_1d(vonkarm, &
           & neutral, stable_option, rich_crit, zeta_trans, &
           & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
           & del_m, del_t, del_q, .true., avail, ier_l)

      ! check sum results
      w = 0
      w = w + transfer(sum(del_m), w)
      w = w + transfer(sum(del_t), w)
      w = w + transfer(sum(del_q), w)

      ! plug check sum in here>>>
#if defined(__INTEL_COMPILER) || defined(_LF95) 
#define CHKSUM_PROFILE -4596910845317820786
#endif
#if defined(_PGF95)
#define CHKSUM_PROFILE -4596910845317820785
#endif

      print *,'chksum test_profile   : ', w, ' ref ', CHKSUM_PROFILE

      ier = CHKSUM_PROFILE - w

    end subroutine test_profile


end program test

!==============================================================================

#endif 
! _TEST_MONIN_OBUKHOV
