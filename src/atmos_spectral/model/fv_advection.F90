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

module fv_advection_mod

use         fms_mod, only: mpp_pe, mpp_npes, mpp_root_pe, error_mesg, FATAL, write_version_number

use   constants_mod, only : radius, pi

use mpp_domains_mod, only : mpp_define_domains, mpp_update_domains, &
                            mpp_get_compute_domain, domain2D

implicit none
private

character(len=128), parameter :: version = '$Id: fv_advection.F90,v 13.0 2006/03/28 21:17:47 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

type(domain2D), save, public :: advection_domain

logical :: module_is_initialized = .FALSE.
logical :: monotone      = .true.
integer :: is, ie, js, je, pe, npes, nx, ny, nz

real, allocatable, dimension(:) :: y, yy, c, s, cc, dy, dyy, dyyy, dy_plus, dy_minus
real :: dx

public :: fv_advection_init, fv_advection_end
public :: a_grid_horiz_advection

interface a_grid_horiz_advection
   module procedure a_grid_horiz_advection_3d
   module procedure a_grid_horiz_advection_2d
end interface

!===========================================================================================
contains
!===========================================================================================


subroutine fv_advection_init(nx_in, ny_in, yy_in, degrees_lon, advection_layout)

  integer, intent(in) ::  nx_in, ny_in
  real   , intent(in), dimension(:) :: yy_in
  real   , intent(in) :: degrees_lon
  integer, intent(in), optional :: advection_layout(2)


  integer :: layout(2)

  if (module_is_initialized) return

  call write_version_number(version, tagname)

  pe   = mpp_pe()
  npes = mpp_npes()

  nx = nx_in
  ny = ny_in

  allocate( y(ny), c(ny), s(ny), yy(ny+1), cc(ny+1), dy(-1:ny+2), dyy(ny+1), dyyy(0:ny+1) )
  allocate( dy_plus(0:ny+1), dy_minus(0:ny+1) )

  yy = yy_in

  y  = 0.5*(yy(2:ny+1) + yy(1:ny))

  c  = cos(y)
  s  = sin(y)
  cc = cos(yy)

  dy(1:ny) = yy(2:ny+1) - yy(1:ny)  ! distance between half level points (size of grid boxes)
  dy(-1)    = dy(2)
  dy(0)     = dy(1)
  dy(ny+1) = dy(ny)
  dy(ny+2) = dy(ny-1)

  dyy(2:ny) = y(2:ny) - y(1:ny-1)  ! distance between full points 
  dyy(1)    = 2*(y(1) - yy(1))
  dyy(ny+1) = 2*(yy(ny+1) - y(ny))

  dy_plus (0:ny+1) = dy(0:ny+1)/(dy( 0:ny+1) + dy(1:ny+2))
  dy_minus(0:ny+1) = dy(0:ny+1)/(dy(-1:ny  ) + dy(0:ny+1))

  y    =    y*radius
  yy   =   yy*radius
  dy   =   dy*radius
  dyy  =  dyy*radius

  dx = (degrees_lon/360.)*2.0*pi*radius/float(nx)

!1D decomposition along Y only
  layout = (/1,npes/)    
    
  if( present (advection_layout) ) layout = advection_layout

  call mpp_define_domains( (/1,nx,1,ny/), layout, advection_domain, yhalo=2 )

  module_is_initialized=.TRUE.

  call mpp_get_compute_domain( advection_domain, is, ie, js, je )

return
end subroutine fv_advection_init

!===========================================================================================

subroutine a_grid_horiz_advection_3d(ua, va, q, dt, dq_dt, flux)


real, intent(in)   , dimension(:,js:,:) :: ua, va, q
real, intent(in)                        :: dt
real, intent(inout), dimension(:,js:,:) :: dq_dt
logical, optional, intent(in) :: flux

real, dimension(nx,js-2:je+2,size(q,3)) :: vx, qx
real, dimension(nx,js  :je+1,size(q,3)) :: vc
real, dimension(nx,js  :je  ,size(q,3)) :: uc, div

integer :: i, j, k
integer, dimension(nx) :: ii

logical :: flux_local

if(.not.module_is_initialized) then
  call error_mesg('a_grid_horiz_advection','fv_advection_mod is not initialized', FATAL)
endif

flux_local = .false.
if(present(flux)) flux_local = flux

vx = 0.0
qx = 0.0

do i = 1,nx
  ii(i) = i + nx/2
  if (ii(i) > nx) ii(i) = ii(i) - nx
end do

vx(:, js:je, :) = va(:, js:je, :)
qx(:, js:je, :) = q (:, js:je, :)

call mpp_update_domains(vx, advection_domain)
call mpp_update_domains(qx, advection_domain)

if(js == 1) then
  do i = 1,nx
    vx(i, 0,:) = - vx(ii(i),1,:)
    qx(i, 0,:) =   qx(ii(i),1,:)
    qx(i,-1,:) =   qx(ii(i),2,:)
  end do
endif

if(je == ny) then
  do i = 1,nx
    vx(i,ny+1,:) = - vx(ii(i),ny  ,:)
    qx(i,ny+1,:) =   qx(ii(i),ny  ,:)
    qx(i,ny+2,:) =   qx(ii(i),ny-1,:)
  end do
endif

uc(2:nx,js:je,:)   = 0.5*(ua(1:nx-1, js:je  ,:) + ua(2:nx,js:je,:))
uc(1   ,js:je,:)   = 0.5*(ua(nx    , js:je  ,:) + ua(1   ,js:je,:))

 do k=1,size(vc,3)
   do j=js,je+1
     do i=1,nx
       vc(i,j,k) = 0.5*(vx(i,j-1,k) + vx(i,j,k))
     enddo
   enddo
 enddo

if(.not.flux_local) then 
  do j = js,je
    div(:,j,:) = (vc(:,j+1,:)*cc(j+1) - vc(:,j,:)*cc(j))/(c(j)*dy(j))
  enddo

  do j = js, je
    div(1:nx-1,j,:) = div(1:nx-1,j,:) + (uc(2:nx,j,:) - uc(1:nx-1,j,:))/(c(j)*dx)
    div(nx    ,j,:) = div(nx    ,j,:) + (uc(1   ,j,:) - uc(nx    ,j,:))/(c(j)*dx)
  enddo

  dq_dt = dq_dt + q*div
endif

call advection_sphere_3d(dq_dt, dt, qx, uc, vc, ua, va)

return
end subroutine a_grid_horiz_advection_3d

!===========================================================================================

subroutine a_grid_horiz_advection_2d(ua, va, q, dt, dq_dt, flux)

real, intent(in),    dimension(:,js:) :: ua, va, q
real, intent(in)                      :: dt
real, intent(inout), dimension(:,js:) :: dq_dt
logical, intent(in), optional :: flux

real, dimension(nx,js:je ,1) :: ua_3d, va_3d, q_3d, dq_dt_3d

q_3d     (:,js:je,1)   = q     (:,js:je)
ua_3d    (:,js:je,1)   = ua    (:,js:je)
va_3d    (:,js:je,1)   = va    (:,js:je)
dq_dt_3d (:,js:je,1)   = dq_dt (:,js:je)

if(present(flux)) then
  call a_grid_horiz_advection_3d(ua_3d, va_3d, q_3d, dt, dq_dt_3d, flux = flux)
else 
  call a_grid_horiz_advection_3d(ua_3d, va_3d, q_3d, dt, dq_dt_3d)
endif

dq_dt  =  dq_dt_3d(:,:,1)

return
end subroutine a_grid_horiz_advection_2d

!===========================================================================================

subroutine advection_sphere_3d(dq_dt, dt, q, uc, vc, ua, va)

real, intent(in)   , dimension(:,js-2:,:) :: q
real, intent(in)   , dimension(:,js  :,:) :: vc
real, intent(in)   , dimension(:,js  :,:) :: uc, ua, va
real, intent(in)                          :: dt
real, intent(inout), dimension(:,js  :,:) :: dq_dt

real, dimension(nx, js  :je  , size(q,3)) :: q2
real, dimension(nx, js-2:je+2, size(q,3)) :: q1

integer, dimension(nx) :: ii
integer :: i


call semi_x_3d(q1(:,js:je,:), ua(:,js:je,:), q(:,js :je ,:), 0.5*dt)
q1(:,js:je,:) = q(:,js:je,:) + q1(:,js:je,:)

call semi_y_3d(q2(:,js:je,:), va(:,js:je,:), q(:,js-2:je+2,:), 0.5*dt)
q2(:,js:je,:) = q(:,js:je,:) + q2(:,js:je,:)

call mpp_update_domains(q1, advection_domain)

do i = 1,nx
  ii(i) = i + nx/2
  if (ii(i) > nx) ii(i) = ii(i) - nx
end do

if(js == 1) then
  do i = 1,nx
    q1(i, 0,:) =   q1(ii(i),1,:)
    q1(i,-1,:) =   q1(ii(i),2,:)
  end do
endif

if(je == ny) then
  do i = 1,nx
    q1(i,ny+1,:) =   q1(ii(i),ny  ,:)
    q1(i,ny+2,:) =   q1(ii(i),ny-1,:)
  end do
endif

call vanleer_x_3d     (dq_dt(:,js:je,:), uc(:,js:je,:)  , q2(:,js:je,:)    , dt)

call vanleer_sphere_3d(dq_dt(:,js:je,:), vc(:,js:je+1,:), q1(:,js-2:je+2,:), dt)

return
end subroutine advection_sphere_3d

!===========================================================================================

subroutine vanleer_sphere_3d(dq_dt, vc, q, dt)

real, intent(in),    dimension(:,js-2:,:) :: q
real, intent(in),    dimension(:,js  :,:) :: vc
real, intent(in)                          :: dt
real, intent(inout), dimension(:,js  :,:) :: dq_dt

real, dimension(nx, js  :je+1, size(q,3)) :: flux
real, dimension(nx ,js-1:je+1, size(q,3)) :: s

real, dimension(js-1:je+1) :: dtdy 
real, dimension(js  :je) :: dyc
integer :: j

dtdy(js-1:je+1) = dt/(dy(js-1:je+1)) 
dyc (js  :je) = 1.0/(dy(js:je)*c(js:je))

call slope_sphere(s(:,js-1:je+1,:), q(:,js-2:je+2,:))

do j = js,je+1
  where (vc(:,j,:) >= 0.0) 
    flux(:,j,:) = vc(:,j,:)*cc(j) * &
             (q(:,j-1,:) + 0.5*s(:,j-1,:)*(1.0 - dtdy(j-1)*vc(:,j,:)))
  elsewhere
    flux(:,j,:) = vc(:,j,:)*cc(j) * &
             (q(:,j,:)   - 0.5*s(:,j,:)  *(1.0 + dtdy(j)  *vc(:,j,:)))
  end where
end do

if(js == 1) flux(:,js,: ) = 0.0
if(je == ny) flux(:,je+1,:) = 0.0

do j = js,je
  dq_dt(:,j,:) = dq_dt(:,j,:) - dyc(j)*(flux(:,j+1,:) - flux(:,j,:))
end do

return
end subroutine vanleer_sphere_3d

!===========================================================================================

subroutine vanleer_x_3d(dq_dt, uc, q, dt)

real, intent(in),  dimension(:,js:,:)   :: uc, q
real, intent(in)                        :: dt
real, intent(inout), dimension(:,js:,:) :: dq_dt

real   , dimension(nx+1,js:je,size(q,3)) :: flux
real   , dimension(nx  ,js:je,size(q,3)) :: b, bb, qq, ss, s
integer, dimension(nx  ,js:je,size(q,3)) :: ii

integer :: i, j, k

do j = js,je
  b(:,j,:)  = uc(:,j,:)*dt/(dx*c(j))
end do
bb = b - int(b)

flux = 0.0
do j = js, je   ! try doing one row at a time to avoid unecessary computations
  if(maxval(abs(b(:,j:j,:))) > 1.0 ) &
      call integer_flux_x(flux(1:nx,j:j,:), b(:,j:j,:), q(:,j:j,:))
end do
call slope_x(s, q)
call find_cell_x(ii, b)

do k = 1, size(q,3)
  do j = js,je
    do i = 1, nx
      qq(i,j,k)  = q(ii(i,j,k),j,k)
      ss(i,j,k)  = s(ii(i,j,k),j,k)
    end do
  end do
end do

flux(1:nx,:,:) = flux(1:nx,:,:) + bb*(qq + 0.5*ss*(sign(1.0,bb) - bb))
flux(nx+1,:,:) = flux(1,:,:)

dq_dt = dq_dt - (flux(2:nx+1,:,:) - flux(1:nx,:,:))/dt

do j = js,je
  flux(:,j,:)  = flux(:,j,:)*dt/(dx*c(j))
end do


return
end subroutine vanleer_x_3d

!===========================================================================================

subroutine semi_x_3d(dq, ua, q, dt)

real, intent(in),  dimension(:,js:,:) :: ua, q
real, intent(in)                      :: dt
real, intent(out), dimension(:,js:,:) :: dq

real   , dimension(nx,js:je,size(q,3)) :: b, bb, q_left, q_right
integer, dimension(nx,js:je,size(q,3)) :: ii, i_left, i_right

integer :: i, j, k

do j = js,je
  b(:,j,:)  = ua(:,j,:)*dt/(dx*c(j))
end do

call find_cell_x(ii, b)

i_left  = ii
i_right = i_left + 1
where(i_right.gt.nx) i_right = 1

bb = b - floor(b)

do k = 1, size(q,3)
  do j = js,je
    do i = 1, nx
      q_left (i,j,k) = q(i_left (i,j,k),j,k)
      q_right(i,j,k) = q(i_right(i,j,k),j,k)
    end do
  end do
end do

dq(:,js:je,:) = bb(:,js:je,:)*q_left (:,js:je,:) + (1.0 - bb(:,js:je,:))*q_right(:,js:je,:) &
               - q(:,js:je,:)

return
end subroutine semi_x_3d

!===========================================================================================

subroutine semi_y_3d(dq, va, qx, dt)

real, intent(out), dimension(:,js  :,:) :: dq
real, intent(in),  dimension(:,js  :,:) :: va
real, intent(in),  dimension(:,js-2:,:) :: qx
real, intent(in)                        :: dt

integer :: j

do j = js, je
  where (va(:,j,:) >= 0.0) 
    dq(:,j,:) = va(:,j,:)*dt*(qx(:,j-1,:) - qx(:,j  ,:))/dyy(j)
  elsewhere
    dq(:,j,:) = va(:,j,:)*dt*(qx(:,j  ,:) - qx(:,j+1,:))/dyy(j+1)
  end where
enddo

return
end subroutine semi_y_3d

!===========================================================================================

subroutine find_cell_x(ii,b)
!!dir$ INLINEALWAYS find_cell_x
integer, intent(out),  dimension(:,:,:) :: ii
real   , intent(in) ,  dimension(:,:,:) :: b

integer :: i

do i = 1,nx
  ii(i,:,:) = i-1 
end do

ii = ii - floor(b)

where(ii.gt.nx) ii = ii - nx
where(ii.lt.1 ) ii = ii + nx

return
end subroutine find_cell_x

!===========================================================================================

subroutine slope_x(slope, q)

real, intent(in),  dimension(:,:,:) :: q
real, intent(out), dimension(:,:,:) :: slope

real, dimension(size(q,1),size(q,2),size(q,3)) :: grad, q_min, q_max
integer :: n

n = size(q,1)

grad(2:n,:,:)  = q(2:n,:,:)-q(1:n-1,:,:)
grad(1,:,:)    = q(1,:,:)-q(n,:,:)

slope(1:n-1,:,:) = (grad(2:n,:,:) + grad(1:n-1,:,:))/2
slope(n,:,:) = (grad(1,:,:) + grad(n,:,:))/2

if(monotone) then

  q_min(2:n-1,:,:) = min(q(1:n-2,:,:),q(2:n-1,:,:),q(3:n,:,:))
  q_min(1,:,:)     = min(q(n,:,:)    ,q(1,:,:)    ,q(2,:,:)  )
  q_min(n,:,:)     = min(q(n-1,:,:)  ,q(n,:,:)    ,q(1,:,:)  )

  q_max(2:n-1,:,:) = max(q(1:n-2,:,:),q(2:n-1,:,:),q(3:n,:,:))
  q_max(1,:,:)     = max(q(n,:,:)    ,q(1,:,:)    ,q(2,:,:)  )
  q_max(n,:,:)     = max(q(n-1,:,:)  ,q(n,:,:)    ,q(1,:,:)  )

  slope = sign(1.0,slope) * min( abs(slope), 2.0*(q - q_min), 2.0*(q_max - q))
else
  slope = sign(1.0,slope) * min( abs(slope), 2.0*q )
end if

return
end subroutine slope_x

!===========================================================================================

subroutine integer_flux_x(flux, c, q)

real, intent(out),  dimension(:,:,:) :: flux
real, intent(in) ,  dimension(:,:,:) :: c, q

integer, dimension(size(c,1),size(c,2),size(c,3)) :: ii
integer :: n, i, j, k

n = size(c,1)
ii = int(c)

flux = 0.0

do k = 1, size(c,3)
  do j = 1, size(c,2)

      do i = 1, size(c,1)
        if(ii(i,j,k) >= 1) then
          if(i-ii(i,j,k) >= 1) then
             flux(i,j,k) = sum(q(i-ii(i,j,k):i-1,j,k))
          else
             flux(i,j,k) = sum(q(1:i-1,j,k)) + sum(q(i-ii(i,j,k)+n:n,j,k))
          end if
        else if (ii(i,j,k) <= -1) then
          if(i-1-ii(i,j,k) <= n) then
            flux(i,j,k) = -sum(q(i:i-1-ii(i,j,k),j,k))
          else
            flux(i,j,k) = -sum(q(i:n,j,k)) - sum(q(1:i-1-ii(i,j,k)-n,j,k))
          end if
        end if
      end do

  end do
end do
return
end subroutine integer_flux_x

!===========================================================================================

subroutine slope_sphere(slope, q)

real, intent(in),  dimension(:,js-2:,:) :: q
real, intent(out), dimension(:,js-1:,:) :: slope

real, dimension(nx,js-1:je+1,size(q,3)) :: q_max, q_min
integer :: j

do j = js-1, je+1
  slope(:,j,:) = (q(:,j+1,:) - q(:,j  ,:))*dy_plus (j) &
               + (q(:,j  ,:) - q(:,j-1,:))*dy_minus(j)
end do

if(monotone) then
  q_min = min(q(:,js-2:je,:),q(:,js-1:je+1,:),q(:,js:je+2,:))
  q_max = max(q(:,js-2:je,:),q(:,js-1:je+1,:),q(:,js:je+2,:))

  slope = sign(1.0,slope) * &
        min( abs(slope), 2.0*(q(:,js-1:je+1,:) - q_min), 2.0*(q_max - q(:,js-1:je+1,:)) )
else
  slope = sign(1.0,slope) * min( abs(slope), 2.0*q(:,js-1:je+1,:) )
end if
return
end subroutine slope_sphere

!===========================================================================================

subroutine solid_body(u, v)

!  used for check-out only

real, intent(out), dimension(nx,ny) :: u,v
real :: beta
integer :: i,j

beta = 45.0

do j = 1, ny
  do i = 1, nx
    v(i,j) = - sin(beta*pi/180.0)*sin(2.*pi*float(i)/float(nx))
    u(i,j) =   cos(beta*pi/180.0)*c(j)        &
             + sin(beta*pi/180.0)*cos(2.*pi*float(i)/float(nx))*s(j)
  end do
end do

end subroutine solid_body

!===========================================================================================
subroutine fv_advection_end

if(.not.module_is_initialized) return

deallocate(y, yy, c, s, cc, dy, dyy, dyyy, dy_plus, dy_minus)
module_is_initialized = .false.

return
end subroutine fv_advection_end
!===========================================================================================

end module fv_advection_mod
