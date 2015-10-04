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

module vert_advection_mod

!-------------------------------------------------------------------------------

use fms_mod, only: error_mesg, FATAL, write_version_number, stdout
use mpp_mod, only: mpp_sum, mpp_max            ,mpp_pe,mpp_sync

implicit none
private

public :: vert_advection, vert_advection_end
! for optional argument: scheme
integer, parameter, public :: SECOND_CENTERED          = 101
integer, parameter, public :: FOURTH_CENTERED          = 102
integer, parameter, public :: FINITE_VOLUME_LINEAR     = 103
integer, parameter, public :: FINITE_VOLUME_PARABOLIC  = 104
integer, parameter, public :: FINITE_VOLUME_PARABOLIC2 = 105
integer, parameter, public :: SECOND_CENTERED_WTS      = 106
integer, parameter, public :: FOURTH_CENTERED_WTS      = 107
integer, parameter, public :: VAN_LEER_LINEAR = FINITE_VOLUME_LINEAR
! for optional argument: form
integer, parameter, public :: FLUX_FORM = 201, ADVECTIVE_FORM = 202
! for optional argument: flags
integer, parameter, public :: WEIGHTED_TENDENCY=1
integer, parameter, public :: OUTFLOW_BOUNDARY=2

character(len=128), parameter :: version = '$Id: vert_advection.F90,v 17.0 2009/07/21 03:00:04 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

logical :: module_is_initialized = .false.

interface vert_advection
   module procedure vert_advection_1d, vert_advection_2d, vert_advection_3d
end interface

! buffers for coefficients used by the parabolic scheme
  real, allocatable :: zwts(:,:,:,:), dzs(:,:,:)
  integer :: nlons = 0, nlats = 0, nlevs = 0

! for cfl diagnostics with finite volume schemes
  real    :: cflmaxc = 0.
  real    :: cflmaxx = 0.
  integer :: cflerr  = 0

contains

!-------------------------------------------------------------------------------

 subroutine vert_advection_3d ( dt, w, dz, r, rdt, mask, scheme, form, flags )

 real, intent(in)                    :: dt
 real, intent(in),  dimension(:,:,:) :: w, dz, r
 real, intent(out), dimension(:,:,:) :: rdt
 real,    intent(in), optional :: mask(:,:,:)
 integer, intent(in), optional :: scheme, form, flags

! INPUT
!   dt  = time step in seconds
!   w   = advecting velocity at the vertical boundaries of the grid boxes
!         does not assume velocities at top and bottom are zero
!         units = [units of dz / second]
!   dz  = depth of model layers in arbitrary units (usually pressure)
!   r   = advected quantity in arbitrary units
!
! OUTPUT
!   rdt = advective tendency for quantity "r" (units depend on optional flags argument)
!           the default units = [units of r / second]
!           if flags=WEIGHTED_TENDENCY then units = [units of r * units of dz / second]
!
! OPTIONAL INPUT
!   mask   = mask for below ground layers,
!            where mask > 0 for layers above ground
!   scheme = differencing scheme, use one of these values:
!               SECOND_CENTERED = second-order centered
!               FOURTH_CENTERED = fourth-order centered
!               SECOND_CENTERED_WTS = second-order centered (assuming unequal level spacing)
!               FOURTH_CENTERED_WTS = fourth-order centered (assuming unequal level spacing)
!               FINITE_VOLUME_LINEAR    = piecewise linear, finite volume (van Leer)
!               VAN_LEER_LINEAR         = same as FINITE_VOLUME_LINEAR
!               FINITE_VOLUME_PARABOLIC = piecewise parabolic, finite volume (PPM)
!               FINITE_VOLUME_PARABOLIC2 = piecewise parabolic, finite volume (PPM)
!                                          using relaxed montonicity constraint (Lin,2003)
!   form   = form of equations, use one of these values:
!               FLUX_FORM      = solves for -d(wr)/dt
!               ADVECTIVE_FORM = solves for -w*d(r)/dt
!   flags  = additional optional flags
!               WEIGHTED_TENDENCY = output tendency (rdt) has units = [units of r * units of dz / second]
!               OUTFLOW_BOUNDARY  = advection is computed at the top and bottom for outflow boundaries
!                                   this option is only valid for finite volume schemes
!
! NOTE
!   size(w,3) == size(dz,3)+1 == size(r,3)+1 == size(rdt,3)+1 == size(mask,3)+1

 real, dimension(size(r,1),size(r,2),size(r,3)) :: slp, r_left, r_right
 real, dimension(size(w,1),size(w,2),size(w,3)) :: flux
 real, dimension(0:3,size(r,1),size(r,2),size(r,3)) :: zwt
 real    :: xx, a, b, rm, r6, rst, wt
 real    :: tt, c1, c2
 real    :: small = 1.e-6
 logical :: test_1
 logical :: do_weighted_dt, do_outflow_bnd
 integer :: i, j, k, ks, ke, kstart, kend, iflags
 integer :: diff_scheme, eqn_form

 real :: cn, rsum, dzsum, dtw
 integer :: kk

 if(.not.module_is_initialized) then
   call write_version_number(version, tagname)
   module_is_initialized = .true.
 endif

 ! set default values for optional arguments
   diff_scheme = VAN_LEER_LINEAR
   eqn_form    = FLUX_FORM
   if (present(scheme)) diff_scheme = scheme
   if (present(form))   eqn_form    = form
 ! set optional flags
   iflags = 0; if (present(flags)) iflags = flags
   do_weighted_dt = btest(iflags,0)
   do_outflow_bnd = btest(iflags,1)
   

 ! note: size(r,3)+1 = size(w,3)
   if (size(w,3) /= size(r,3)+1) &
      call error_handler ('vertical dimension of input arrays inconsistent')

 ! vertical indexing
   ks = 1; ke = size(r,3)

 ! start and end indexing for finite volume fluxes
   kstart = ks+1;  kend = ke
   if (do_outflow_bnd) then
       kstart = ks;  kend = ke+1
   endif

 ! set fluxes at boundaries
 ! most likely w = 0 at these points
 ! but for outflow boundaries set flux to zero

   if (do_outflow_bnd) then
       flux(:,:,ks)   = 0.
       flux(:,:,ke+1) = 0.
   else
       flux(:,:,ks)   = w(:,:,ks)  *r(:,:,ks)
       flux(:,:,ke+1) = w(:,:,ke+1)*r(:,:,ke)
   endif

   select case (diff_scheme)

   !------ 2nd-order centered scheme assuming variable grid spacing ------
      case (SECOND_CENTERED_WTS)
         do k = ks+1, ke
         do j = 1, size(r,2)
         do i = 1, size(r,1)
            wt = dz(i,j,k-1)/(dz(i,j,k-1)+dz(i,j,k))
            rst = r(i,j,k-1) + wt*(r(i,j,k)-r(i,j,k-1))
            flux(i,j,k) = w(i,j,k)*rst
         enddo
         enddo
         enddo

   !------ 2nd-order centered scheme assuming uniform grid spacing ------
      case (SECOND_CENTERED)
         do k = ks+1, ke
         do j = 1, size(r,2)
         do i = 1, size(r,1)
            rst = 0.5*(r(i,j,k)+r(i,j,k-1))
            flux(i,j,k) = w(i,j,k)*rst
         enddo
         enddo
         enddo

   !------ 4th-order centered scheme assuming variable grid spacing ------
      case (FOURTH_CENTERED_WTS)
         call compute_weights ( dz, zwt )
         call slope_z ( r, dz, slp, limit=.false., linear=.false. )
         if (present(mask)) then
          ! second order if adjacent to ground
            do k = ks+2, ke-1
            do j = 1, size(r,2)
            do i = 1, size(r,1)
               if (mask(i,j,k+1) > small) then
                  rst = r(i,j,k-1) + zwt(1,i,j,k)*(r(i,j,k)-r(i,j,k-1)) &
                        - zwt(2,i,j,k)*slp(i,j,k) + zwt(3,i,j,k)*slp(i,j,k-1)
               else
                  rst = r(i,j,k-1) + zwt(0,i,j,k)*(r(i,j,k)-r(i,j,k-1))
               endif
            flux(i,j,k) = w(i,j,k)*rst
            enddo
            enddo
            enddo
         else
          ! no mask: always fourth order
            do k = ks+2, ke-1
            do j = 1, size(r,2)
            do i = 1, size(r,1)
               rst = r(i,j,k-1) + zwt(1,i,j,k)*(r(i,j,k)-r(i,j,k-1)) &
                        - zwt(2,i,j,k)*slp(i,j,k) + zwt(3,i,j,k)*slp(i,j,k-1)
               flux(i,j,k) = w(i,j,k)*rst
            enddo
            enddo
            enddo
         endif
         ! second order at top and bottom
         do j = 1, size(r,2)
         do i = 1, size(r,1)
            wt  = dz(i,j,ks)/(dz(i,j,ks)+dz(i,j,ks+1))
            rst = r(i,j,ks) + wt*(r(i,j,ks+1)-r(i,j,ks))
            flux(i,j,ks+1) = w(i,j,ks+1)*rst
            wt  = dz(i,j,ke-1)/(dz(i,j,ke-1)+dz(i,j,ke))
            rst = r(i,j,ke-1) + wt*(r(i,j,ke)-r(i,j,ke-1))
            flux(i,j,ke) = w(i,j,ke)*rst
         enddo
         enddo

   !------ 4th-order centered scheme assuming uniform grid spacing ------
      case (FOURTH_CENTERED)
         c1 = 7./12.;  c2 = 1./12.
         if (present(mask)) then
          ! second order if adjacent to ground
            do k = ks+2, ke-1
            do j = 1, size(r,2)
            do i = 1, size(r,1)
               if (mask(i,j,k+1) > small) then
                  rst = c1*(r(i,j,k)+r(i,j,k-1)) - c2*(r(i,j,k+1)+r(i,j,k-2))
               else
                  rst = 0.5*(r(i,j,k)+r(i,j,k-1))
               endif
               flux(i,j,k) = w(i,j,k)*rst
            enddo
            enddo
            enddo
         else
          ! no mask: always fourth order
            do k = ks+2, ke-1
            do j = 1, size(r,2)
            do i = 1, size(r,1)
               rst = c1*(r(i,j,k)+r(i,j,k-1)) - c2*(r(i,j,k+1)+r(i,j,k-2))
               flux(i,j,k) = w(i,j,k)*rst
            enddo
            enddo
            enddo
         endif
         ! second order at top and bottom
         do j = 1, size(r,2)
         do i = 1, size(r,1)
            rst = 0.5*(r(i,j,ks+1)+r(i,j,ks  ))
            flux(i,j,ks+1) = w(i,j,ks+1)*rst
            rst = 0.5*(r(i,j,ke  )+r(i,j,ke-1))
            flux(i,j,ke) = w(i,j,ke)*rst
         enddo
         enddo

   !------ finite volume scheme using piecewise linear method ------
      case (FINITE_VOLUME_LINEAR)
       ! slope along the z-axis
         call slope_z ( r, dz, slp )
         do k = kstart, kend
         do j = 1, size(r,2)
         do i = 1, size(r,1)
            if (w(i,j,k) >= 0.) then
               if (k == ks) cycle ! inflow
               cn = dt*w(i,j,k)/dz(i,j,k-1)
               rst = r(i,j,k-1) + 0.5*slp(i,j,k-1)*(1.-cn)
            else
               if (k == ke+1) cycle ! inflow
               cn = -dt*w(i,j,k)/dz(i,j,k)
               rst = r(i,j,k  ) - 0.5*slp(i,j,k  )*(1.-cn)
            endif
            flux(i,j,k) = w(i,j,k)*rst
            if (cn > 1.) cflerr = cflerr+1
            cflmaxc = max(cflmaxc,cn)
            cflmaxx = cflmaxc ! same for linear scheme
         enddo
         enddo
         enddo

   !------ finite volume scheme using piecewise parabolic method (PPM) ------
      case (FINITE_VOLUME_PARABOLIC:FINITE_VOLUME_PARABOLIC2)
         call compute_weights ( dz, zwt )
         call slope_z ( r, dz, slp, linear=.false. )
         do k = ks+2, ke-1
         do j = 1, size(r,2)
         do i = 1, size(r,1)
            r_left(i,j,k) = r(i,j,k-1) + zwt(1,i,j,k)*(r(i,j,k)-r(i,j,k-1)) &
                   - zwt(2,i,j,k)*slp(i,j,k) + zwt(3,i,j,k)*slp(i,j,k-1)        
            r_right(i,j,k-1) = r_left(i,j,k)
            ! coming out of this loop, all we need is r_left and r_right
         enddo
         enddo
         enddo

         ! boundary values  ! masks ???????

         do j = 1, size(r,2)
         do i = 1, size(r,1)
           r_left (i,j,ks+1) = r(i,j,ks+1) - 0.5*slp(i,j,ks+1)
           r_right(i,j,ke-1) = r(i,j,ke-1) + 0.5*slp(i,j,ke-1)

         ! pure upstream advection near boundary
         ! r_left (i,j,ks) = r(i,j,ks)
         ! r_right(i,j,ks) = r(i,j,ks)
         ! r_left (i,j,ke) = r(i,j,ke)
         ! r_right(i,j,ke) = r(i,j,ke)

         ! make linear assumption near boundary
         ! NOTE: slope is zero at ks and ks therefore
         !       this reduces to upstream advection near boundary
           r_left (i,j,ks) = r(i,j,ks) - 0.5*slp(i,j,ks)
           r_right(i,j,ks) = r(i,j,ks) + 0.5*slp(i,j,ks)
           r_left (i,j,ke) = r(i,j,ke) - 0.5*slp(i,j,ke)
           r_right(i,j,ke) = r(i,j,ke) + 0.5*slp(i,j,ke)
         enddo
         enddo

     ! monotonicity constraint

       if (diff_scheme == FINITE_VOLUME_PARABOLIC2) then
         ! limiters from Lin (2003), Equation 6 (relaxed constraint)
           do k = ks, ke
           do j = 1, size(r,2)
           do i = 1, size(r,1)
              r_left (i,j,k) = r(i,j,k) - sign( min(abs(slp(i,j,k)),abs(r_left (i,j,k)-r(i,j,k))), slp(i,j,k) )
              r_right(i,j,k) = r(i,j,k) + sign( min(abs(slp(i,j,k)),abs(r_right(i,j,k)-r(i,j,k))), slp(i,j,k) )
           enddo
           enddo
           enddo
       else
         ! limiters from Colella and Woodward (1984), Equation 1.10
           do k = ks, ke
           do j = 1, size(r,2)
           do i = 1, size(r,1)
              test_1 = (r_right(i,j,k)-r(i,j,k))*(r(i,j,k)-r_left(i,j,k)) <= 0.0
              if (test_1) then
                 r_left(i,j,k)  = r(i,j,k)
                 r_right(i,j,k) = r(i,j,k)
              endif
              if (k == ks .or. k == ke) cycle
              rm = r_right(i,j,k) - r_left(i,j,k)
              a = rm*(r(i,j,k) - 0.5*(r_right(i,j,k) + r_left(i,j,k)))
              b = rm*rm/6.
              if (a >  b) r_left (i,j,k) = 3.0*r(i,j,k) - 2.0*r_right(i,j,k)
              if (a < -b) r_right(i,j,k) = 3.0*r(i,j,k) - 2.0*r_left (i,j,k)
           enddo
           enddo
           enddo
       endif

         ! compute fluxes at interfaces

           tt = 2./3.
           do k = kstart, kend
           do j = 1, size(r,2)
           do i = 1, size(r,1)
              if (w(i,j,k) >= 0.) then
                  if (k == ks) cycle ! inflow
                  cn = dt*w(i,j,k)/dz(i,j,k-1)
                  kk = k-1
                  ! extension for Courant numbers > 1
                  if (cn > 1.) then
                      rsum = 0.
                      dzsum = 0.
                      dtw = dt*w(i,j,k)
                      do while (dzsum+dz(i,j,kk) < dtw)
                         if (kk == 1) then
                             exit
                         endif
                         dzsum = dzsum + dz(i,j,kk)
                          rsum =  rsum +  r(i,j,kk)
                         kk = kk-1
                      enddo
                      xx = (dtw-dzsum)/dz(i,j,kk)
                  else
                      xx = cn
                  endif
                  rm = r_right(i,j,kk) - r_left(i,j,kk)
                  r6 = 6.0*(r(i,j,kk) - 0.5*(r_right(i,j,kk) + r_left(i,j,kk)))
                  if (kk == ks) r6 = 0.
                  rst = r_right(i,j,kk) - 0.5*xx*(rm - (1.0 - tt*xx)*r6)
                  ! extension for Courant numbers > 1
                  if (cn > 1.) rst = (xx*rst + rsum)/cn
              else
                  if (k == ke+1) cycle ! inflow
                  cn = - dt*w(i,j,k)/dz(i,j,k)
                  kk = k
                  ! extension for Courant numbers > 1
                  if (cn > 1.) then
                      rsum = 0.
                      dzsum = 0.
                      dtw = -dt*w(i,j,k)
                      do while (dzsum+dz(i,j,kk) < dtw)
                         if (kk == ks) then
                             exit
                         endif
                         dzsum = dzsum + dz(i,j,kk)
                          rsum =  rsum +  r(i,j,kk)
                         kk = kk+1
                      enddo
                      xx = (dtw-dzsum)/dz(i,j,kk)
                  else
                      xx = cn
                  endif
                  rm = r_right(i,j,kk) - r_left(i,j,kk)
                  r6 = 6.0*(r(i,j,kk) - 0.5*(r_right(i,j,kk) + r_left(i,j,kk)))
                  if (kk == ke) r6 = 0.
                  rst = r_left(i,j,kk) + 0.5*xx*(rm + (1.0 - tt*xx)*r6)
                  ! extension for Courant numbers > 1
                  if (cn > 1.) rst = (xx*rst + rsum)/cn
              endif
              flux(i,j,k) = w(i,j,k)*rst
              if (xx > 1.) cflerr = cflerr+1
              cflmaxx = max(cflmaxx,xx)
              cflmaxc = max(cflmaxc,cn)
           enddo
           enddo
           enddo


      case default
        ! ERROR
          call error_handler ('invalid value for optional argument scheme')
   end select


 ! vertical advective tendency

   select case (eqn_form)
      case (FLUX_FORM)
         if (do_weighted_dt) then
            do k = ks, ke
               rdt (:,:,k) = - (flux(:,:,k+1) - flux (:,:,k)) 
            enddo
         else
            do k = ks, ke
               rdt (:,:,k) = - (flux(:,:,k+1) - flux (:,:,k))  / dz(:,:,k)
            enddo
         endif
      case (ADVECTIVE_FORM)
         if (do_weighted_dt) then
            do k = ks, ke
               rdt (:,:,k) = - (flux(:,:,k+1) - flux (:,:,k) - &
                            r(:,:,k)*(w(:,:,k+1)-w(:,:,k)))
            enddo
         else
            do k = ks, ke
               rdt (:,:,k) = - (flux(:,:,k+1) - flux (:,:,k) - &
                            r(:,:,k)*(w(:,:,k+1)-w(:,:,k))) / dz(:,:,k)
            enddo
         endif
      case default
        ! ERROR
          call error_handler ('invalid value for optional argument form')
   end select


 end subroutine vert_advection_3d

!-------------------------------------------------------------------------------

 subroutine vert_advection_end

  integer :: outunit
  ! deallocate storage
    if (allocated(zwts)) deallocate(zwts)
    if (allocated(dzs))  deallocate(dzs)

  ! cfl diagnostics
    call mpp_max (cflmaxc)
    call mpp_max (cflmaxx)
    call mpp_sum (cflerr) ! integer sum
    if (cflmaxc > 0.) then
        outunit = stdout()
        write (outunit,10) cflmaxc, cflmaxx, cflerr
    endif
 10 format (/,' Vertical advection (atmosphere):',    &
            /,'     maximum CFL =',f10.6,'; ',f10.6,  &
            /,'     number of CFL errors =',i5,/)
        
 end subroutine vert_advection_end

!-------------------------------------------------------------------------------

 subroutine slope_z ( r, dz, slope, limit, linear )
 real, intent(in),  dimension(:,:,:) :: r, dz
 real, intent(out), dimension(:,:,:) :: slope
 logical, intent(in), optional :: limit, linear

!real    :: grad(size(r,1),size(r,2),2:size(r,3))
 real    :: grad(2:size(r,3))
 real    :: rmin, rmax
 integer :: i, j, k, n
 logical :: limiters, dolinear

  limiters = .true.
  if (present(limit))  limiters = limit
  dolinear = .true.
  if (present(linear)) dolinear = linear

  n = size(r,3)

! compute slope (weighted for unequal levels)

  do j = 1, size(r,2)
  do i = 1, size(r,1)

     do k = 2, n
       grad(k) = (r(i,j,k)-r(i,j,k-1))/(dz(i,j,k)+dz(i,j,k-1))
     enddo
     if (dolinear) then
         do k = 2, n-1
           slope(i,j,k) = (grad(k+1)+grad(k))*dz(i,j,k)
         enddo
     else
         do k = 2, n-1
           slope(i,j,k) = (grad(k+1)*(2.*dz(i,j,k-1)+dz(i,j,k)) + &
                           grad(k  )*(2.*dz(i,j,k+1)+dz(i,j,k)))  &
                          *dz(i,j,k)/(dz(i,j,k-1)+dz(i,j,k)+dz(i,j,k+1))
         enddo
     endif
     slope(i,j,1) = 2.*grad(2)*dz(i,j,1)
     slope(i,j,n) = 2.*grad(n)*dz(i,j,n)

   ! apply limiters to slope
     if (limiters) then
        do k = 1, n
          if (k >= 2 .and. k <= n-1) then
            rmin = min(r(i,j,k-1), r(i,j,k), r(i,j,k+1))
            rmax = max(r(i,j,k-1), r(i,j,k), r(i,j,k+1))
            slope(i,j,k) = sign(1.,slope(i,j,k)) *  &
                   min( abs(slope(i,j,k)), 2.*(r(i,j,k)-rmin), 2.*(rmax-r(i,j,k)) )
          else
         !else if (k == 1) then               ! always slope=0
         !  rmin = min(r(i,j,k), r(i,j,k+1))
         !  rmax = max(r(i,j,k), r(i,j,k+1))
         !else if (k == n) then
         !  rmin = min(r(i,j,k-1), r(i,j,k))
         !  rmax = max(r(i,j,k-1), r(i,j,k))
            slope(i,j,k) = 0.
          endif
        enddo
     endif

  enddo
  enddo

 end subroutine slope_z

!-------------------------------------------------------------------------------

 subroutine compute_weights ( dz, zwt )
 real, intent(in),  dimension(:,:,:)    :: dz
 real, intent(out), dimension(0:,:,:,:) :: zwt
 real    :: denom1, denom2, denom3, denom4, num3, num4, x, y
 integer :: i, j, k
 logical :: redo

! check the size of stored coefficients
! need to reallocate if size has changed
   if (nlons /= size(dz,1) .or. nlats /= size(dz,2) .or.  nlevs /= size(dz,3)) then
      if (allocated(zwts)) deallocate(zwts)
      if (allocated(dzs))  deallocate(dzs)
      nlons = size(dz,1)
      nlats = size(dz,2)
      nlevs = size(dz,3)
      allocate (zwts(0:3,nlons,nlats,nlevs))
      allocate (dzs (nlons,nlats,nlevs))
      dzs = -1.
   endif
   
! coefficients/weights for computing values at grid box interfaces
! only recompute coefficients for a column when layer depth has changed

   do j = 1, size(dz,2)
   do i = 1, size(dz,1)

    redo = .false.
    do k=1,size(dz,3)
      if (dz(i,j,k) /= dzs(i,j,k)) then
        redo = .true.
        exit
      endif
    enddo

   if (redo) then
     do k = 3, size(dz,3)-1
       denom1 = 1.0/(dz(i,j,k-1) + dz(i,j,k))
       denom2 = 1.0/(dz(i,j,k-2) + dz(i,j,k-1) + dz(i,j,k) + dz(i,j,k+1))
       denom3 = 1.0/(2*dz(i,j,k-1) +   dz(i,j,k))  
       denom4 = 1.0/(  dz(i,j,k-1) + 2*dz(i,j,k))  
       num3   = dz(i,j,k-2) + dz(i,j,k-1)          
       num4   = dz(i,j,k)   + dz(i,j,k+1)        
       x      = num3*denom3 - num4*denom4        
       y      = 2.0*dz(i,j,k-1)*dz(i,j,k) ! everything up to this point is just
                                          ! needed to compute x1,x1,x3                      
       zwt(0,i,j,k) = dz(i,j,k-1)*denom1                ! = 1/2 in equally spaced case
       zwt(1,i,j,k) = zwt(0,i,j,k) + x*y*denom1*denom2  ! = 1/2 in equally spaced case
       zwt(2,i,j,k) = dz(i,j,k-1)*num3*denom3*denom2    ! = 1/6 ''
       zwt(3,i,j,k) = dz(i,j,k)*num4*denom4*denom2      ! = 1/6 ''
     enddo
     dzs(i,j,:) = dz(i,j,:)
     zwts(0:3,i,j,:) = zwt(0:3,i,j,:)
   else

   ! use previously computed coefficients
     zwt(0:3,i,j,:) = zwts(0:3,i,j,:)
   endif

   enddo
   enddo


 end subroutine compute_weights

!-------------------------------------------------------------------------------
!--------------------------- overloaded versions -------------------------------

 subroutine vert_advection_1d ( dt, w, dz, r, rdt, mask, scheme, form, flags )
 
 real, intent(in)                :: dt
 real, intent(in),  dimension(:) :: w, dz, r
 real, intent(out), dimension(:) :: rdt
 real,    intent(in), optional :: mask(:)
 integer, intent(in), optional :: scheme, form, flags

 real, dimension(1,1,size(r,1)) :: dz3, r3, rdt3, mask3
 real, dimension(1,1,size(w,1)) :: w3

  ! input
    w3 (1,1,:) = w
    dz3(1,1,:) = dz
    r3 (1,1,:) = r

    if (present(mask)) then
       mask3(1,1,:) = mask
       call vert_advection_3d ( dt, w3, dz3, r3, rdt3, mask=mask3, scheme=scheme, form=form, flags=flags )
    else
       call vert_advection_3d ( dt, w3, dz3, r3, rdt3, scheme=scheme, form=form, flags=flags )
    endif

  ! output
    rdt = rdt3(1,1,:)

 end subroutine vert_advection_1d

!-------------------------------------------------------------------------------

 subroutine vert_advection_2d ( dt, w, dz, r, rdt, mask, scheme, form, flags )

 real, intent(in)                  :: dt
 real, intent(in),  dimension(:,:) :: w, dz, r
 real, intent(out), dimension(:,:) :: rdt
 real,    intent(in), optional :: mask(:,:)
 integer, intent(in), optional :: scheme, form, flags

 real, dimension(size(r,1),1,size(r,2)) :: dz3, r3, rdt3, mask3
 real, dimension(size(w,1),1,size(w,2)) :: w3

  ! input
    w3 (:,1,:) = w
    dz3(:,1,:) = dz
    r3 (:,1,:) = r

    if (present(mask)) then
       mask3(:,1,:) = mask
       call vert_advection_3d ( dt, w3, dz3, r3, rdt3, mask=mask3, scheme=scheme, form=form, flags=flags )
    else
       call vert_advection_3d ( dt, w3, dz3, r3, rdt3, scheme=scheme, form=form, flags=flags )
    endif

  ! output
    rdt = rdt3(:,1,:)

 end subroutine vert_advection_2d

!-------------------------------------------------------------------------------

 subroutine error_handler ( message )
 character(len=*), intent(in) :: message

   call error_mesg ('vert_advection', trim(message), FATAL)

!  print *, 'FATAL ERROR in vert_advection'
!  print *, trim(message)
!  stop 111

 end subroutine error_handler

!-------------------------------------------------------------------------------

end module vert_advection_mod

