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

! -*-f90-*-
!#include <fms_platform.h>

! $Id: monin_obukhov_interfaces.h,v 13.0 2006/03/28 21:10:30 fms Exp $

  interface
     
     _PURE subroutine monin_obukhov_diff(vonkarm,                &
          & ustar_min,                                     &
          & neutral, stable_option, rich_crit, zeta_trans, &
          & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

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

    end subroutine monin_obukhov_diff

    _PURE subroutine monin_obukhov_drag_1d(grav, vonkarm,               &
          & error, zeta_min, max_iter, small,                         &
          & neutral, stable_option, rich_crit, zeta_trans, drag_min,  &
          & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
          & drag_q, u_star, b_star, lavail, avail, ier)

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
       real   , intent(in   ), dimension(n)  :: pt      ! potential temperature
       real   , intent(in   ), dimension(n)  :: pt0     ! reference potential temperature
       real   , intent(in   ), dimension(n)  :: z       ! height above surface
       real   , intent(in   ), dimension(n)  :: z0      ! roughness height (height at which log wind profile is zero)
       real   , intent(in   ), dimension(n)  :: zt      ! 
       real   , intent(in   ), dimension(n)  :: zq
       real   , intent(in   ), dimension(n)  :: speed
       real   , intent(inout), dimension(n)  :: drag_m
       real   , intent(inout), dimension(n)  :: drag_t
       real   , intent(inout), dimension(n)  :: drag_q
       real   , intent(inout), dimension(n)  :: u_star
       real   , intent(inout), dimension(n)  :: b_star
       logical, intent(in   )                :: lavail ! whether to use provided mask or not
       logical, intent(in   ), dimension(n)  :: avail  ! provided mask 
       integer, intent(out  )                :: ier

     end subroutine monin_obukhov_drag_1d

     _PURE subroutine monin_obukhov_profile_1d(vonkarm, &
          & neutral, stable_option, rich_crit, zeta_trans, &
          & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
          & del_m, del_t, del_q, lavail, avail, ier)

       real   , intent(in   )                :: vonkarm
       logical, intent(in   )                :: neutral
       integer, intent(in   )                :: stable_option
       real   , intent(in   )                :: rich_crit, zeta_trans
       integer, intent(in   )                :: n
       real,    intent(in)                   :: zref, zref_t
       real,    intent(in) , dimension(n)    :: z, z0, zt, zq, u_star, b_star, q_star
       real,    intent(out), dimension(n)    :: del_m, del_t, del_q
       logical, intent(in)                   :: lavail ! whether to use provided mask or not
       logical, intent(in) , dimension(n)    :: avail  ! provided mask
       integer, intent(  out)                :: ier
     end subroutine monin_obukhov_profile_1d

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

     _PURE subroutine monin_obukhov_stable_mix(stable_option, rich_crit, zeta_trans, &
          &                              n, rich, mix, ier)

       integer, intent(in   )                 :: stable_option
       real   , intent(in   )                 :: rich_crit, zeta_trans
       integer, intent(in   )                 :: n
       real   , intent(in   ), dimension(n)   :: rich
       real   , intent(  out), dimension(n)   :: mix  
       integer, intent(  out)                 :: ier      
     end subroutine monin_obukhov_stable_mix

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
