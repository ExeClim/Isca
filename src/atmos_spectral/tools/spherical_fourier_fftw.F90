module spherical_fourier_fftw_mod
use fftw3
use fms_mod,         only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number

use mpp_mod,         only: mpp_error

use mpp_domains_mod, only: domain1D, mpp_get_compute_domains, mpp_get_domain_components, mpp_get_layout

use   constants_mod, only: pi

use spec_mpp_mod,    only: get_grid_domain, grid_domain, get_spec_domain

use spherical_mod,   only: spherical_init, spherical_end,                &
                           compute_lon_deriv_cos, compute_lat_deriv_cos, &
                           compute_laplacian, compute_vor, compute_div,  &
                           get_spherical_wave, get_fourier_wave,         &
                           get_eigen_laplacian, compute_gradient_cos,    &
                           compute_ucos_vcos, compute_vor_div,           &
                           triangular_truncation, rhomboidal_truncation

use gauss_and_legendre_mod, only: compute_legendre, compute_gaussian

implicit none
private

character(len=128), parameter :: version = '$Id: spherical_fourier_fftw.F90,v 13.0 2019 fms Exp $'
character(len=128), parameter :: tagname = '$Name: siena_201211 $'

interface trans_spherical_to_fourier_fftw
   module procedure trans_spherical_to_fourier_fftw_3d,  &
                    trans_spherical_to_fourier_fftw_2d
end interface
interface trans_fourier_to_spherical_fftw
   module procedure trans_fourier_to_spherical_fftw_3d,  &
                    trans_fourier_to_spherical_fftw_2d
end interface

interface trans_grid_to_spherical_fftw
    module procedure trans_grid_to_spherical_fftw_3d, &
                     trans_grid_to_spherical_fftw_2d 
end interface

logical :: module_is_initialized = .false.
integer :: is, ie, js, je, ms, me, ns, ne
integer :: fourier_inc, num_fourier, fourier_max, lat_max, &
           num_spherical
logical :: make_symmetric_local

public :: spherical_fourier_fftw_init, trans_grid_to_spherical_fftw

contains

!-----------------------------------------------------------------------
subroutine spherical_fourier_fftw_init(radius, lat_max_in, num_fourier_in, &
                   fourier_inc_in, num_spherical_in, south_to_north, make_symmetric)
!-----------------------------------------------------------------------

real,    intent(in) :: radius
integer, intent(in) :: lat_max_in
integer, intent(in) :: num_fourier_in
integer, intent(in) :: fourier_inc_in
integer, intent(in) :: num_spherical_in
logical, intent(in), optional :: south_to_north, make_symmetric

call write_version_number(version, tagname)


if(present(make_symmetric)) then
  make_symmetric_local = make_symmetric
else
  make_symmetric_local = .false.
end if

call get_grid_domain(is, ie, js, je)
call get_spec_domain(ms, me, ns, ne)

lat_max       = lat_max_in

fourier_inc   = fourier_inc_in
num_spherical = num_spherical_in
num_fourier   = num_fourier_in
fourier_max   = num_fourier*fourier_inc



module_is_initialized = .true.

return
end subroutine spherical_fourier_fftw_init



!!!! work on this one
subroutine trans_grid_to_spherical_fftw_3d(grid, spherical)
  real, intent(inout), dimension (is:,:,:) :: grid
  complex, intent(out), dimension (ms:,ns:,:) :: spherical

  ! call forwardsht_real_fortran()

  ! integer :: num_levels
  ! type(C_PTR) :: p, plan
  ! complex(C_DOUBLE_COMPLEX), pointer :: arr(:,:,:)

  ! num_levels = 25

  ! p = fftw_alloc_complex(int(ie * je * num_levels, C_SIZE_T))
  ! plan = fftw_plan_dft_r2c_3d(ie, je, num_levels, grid, spherical, 0)

  ! call fftw_execute_dft_r2c(plan, grid, spherical)

return
end subroutine trans_grid_to_spherical_fftw_3d

subroutine trans_grid_to_spherical_fftw_2d(grid, spherical)
  real, intent(in), dimension (:,:) :: grid
  complex, intent (out), dimension (:,:)  ::  spherical






  ! integer :: num_levels
  ! type(C_PTR) :: p, plan
  ! complex(C_DOUBLE_COMPLEX), pointer :: arr(:,:,:)

  ! num_levels = 25

  ! p = fftw_alloc_complex(int(ie * je, C_SIZE_T))
  ! plan = fftw_plan_dft_r2c_2d(ie, je, grid, spherical, 0)

  ! call fftw_execute_dft_r2c(plan, grid, spherical)





return
end subroutine trans_grid_to_spherical_fftw_2d












subroutine trans_spherical_to_fourier_fftw_3d(spherical, fourier)
  complex, intent(in),  dimension(ms:,ns:,:) :: spherical
  complex, intent(out), dimension(ms:,:,:,:) :: fourier

end subroutine trans_spherical_to_fourier_fftw_3d

subroutine trans_spherical_to_fourier_fftw_2d(spherical, fourier)
complex, intent(in),  dimension(:,:) :: spherical
complex, intent(out), dimension(:,:,:) :: fourier

end subroutine trans_spherical_to_fourier_fftw_2d



subroutine trans_fourier_to_spherical_fftw_3d(spherical, fourier)
complex, intent(in),  dimension(:,:,:) :: spherical
complex, intent(out), dimension(:,:,:) :: fourier

end subroutine trans_fourier_to_spherical_fftw_3d

subroutine trans_fourier_to_spherical_fftw_2d(spherical, fourier)
complex, intent(in),  dimension(:,:) :: spherical
complex, intent(out), dimension(:,:,:) :: fourier

end subroutine trans_fourier_to_spherical_fftw_2d



end module spherical_fourier_fftw_mod