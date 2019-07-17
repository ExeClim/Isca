module spherical_fourier_fftw_mod
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

interface trans_spherical_to_fourier_fftw
   module procedure trans_spherical_to_fourier_fftw_3d,  &
                    trans_spherical_to_fourier_fftw_2d
end interface
interface trans_fourier_to_spherical_fftw
   module procedure trans_fourier_to_spherical_fftw_3d,  &
                    trans_fourier_to_spherical_fftw_2d
end interface


end module spherical_fourier_fftw_mod