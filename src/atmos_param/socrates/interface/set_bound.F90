! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE set_bound_mod
IMPLICIT NONE
CONTAINS

! Subroutine to set the boundary fields (surface and top-of-atmosphere).
!------------------------------------------------------------------------------
SUBROUTINE set_bound(control, dimen, spectrum, bound, n_profile,               &
  t_rad_surf, cos_zenith_angle, solar_irrad, orog_corr,                        &
  l_planet_grey_surface, planet_albedo, planet_emissivity)

USE rad_pcf
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_bound,    ONLY: StrBound, allocate_bound

use soc_constants_mod, only: i_def, r_def

IMPLICIT NONE

! Control options:
TYPE(StrCtrl),      INTENT(IN)  :: control

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! Boundary properties:
TYPE(StrBound),     INTENT(OUT) :: bound

INTEGER(i_def), INTENT(IN) :: n_profile
!   Number of atmospheric profiles for radiation calculations

REAL(r_def), INTENT(IN) :: t_rad_surf(n_profile)
!   Effective radiative temperature over whole grid-box
REAL(r_def), INTENT(IN) :: cos_zenith_angle(n_profile)
!   Cosine of solar zenith angle
REAL(r_def), INTENT(IN) :: solar_irrad(n_profile)
!   Solar irradiance at top-of-atmosphere (mean over timestep)
REAL(r_def), INTENT(IN) :: orog_corr(n_profile)
!   Orographic correction factor

LOGICAL, INTENT(IN) :: l_planet_grey_surface
!   Set a single grey albedo / emissivity for the surface
REAL(r_def), INTENT(IN) :: planet_albedo(:)
!   Surface albedo used for SW calculations
REAL(r_def), INTENT(IN) :: planet_emissivity
!   Surface emissivity used for LW calculations

! Local variables.
INTEGER :: l,j
!   Loop variables


! Allocate structure for the core radiation code interface
CALL allocate_bound(bound, dimen, spectrum)

! Set the radiative characteristics of the surface.
SELECT CASE (control%isolir)

CASE (ip_solar)
  IF (l_planet_grey_surface) THEN
    do j=1,spectrum%basic%n_band
      bound%rho_alb(1:n_profile, ip_surf_alb_diff, j) &
        = planet_albedo
      bound%rho_alb(1:n_profile, ip_surf_alb_dir,  j) &
        = planet_albedo
    end do
  ELSE
    bound%rho_alb(1:n_profile, ip_surf_alb_diff, 1:spectrum%basic%n_band) = 0.0
    bound%rho_alb(1:n_profile, ip_surf_alb_dir,  1:spectrum%basic%n_band) = 0.0
  END IF

CASE (ip_infra_red)
  ! Surface temperature
  DO l=1, n_profile
    bound%t_ground(l) = t_rad_surf(l)
  END DO

  ! Zero the irrelevant direct albedo.
  bound%rho_alb(1:n_profile, ip_surf_alb_dir, 1:spectrum%basic%n_band) = 0.0

  ! Set the diffuse albedo
  IF (l_planet_grey_surface) THEN
    bound%rho_alb(1:n_profile, ip_surf_alb_diff, 1:spectrum%basic%n_band) &
      = 1.0 - planet_emissivity
  ELSE
    bound%rho_alb(1:n_profile, ip_surf_alb_diff, 1:spectrum%basic%n_band) = 0.0
  END IF

END SELECT


! Set the surface basis functions for a Lambertian surface.
bound%n_brdf_basis_fnc=1
! By defining F_{1,0,0,0} to be 4, rho_alb becomes equal to the diffuse albedo.
bound%f_brdf(1, 0, 0, 0)=4.0


IF (control%l_tile) THEN
  ! Set up the surface tiling variables (none for the minute).
  bound%n_tile=3
  bound%n_point_tile=0
END IF


! Set the incident solar flux and orographic correction factor.
IF (control%isolir == ip_solar) THEN
  DO l=1, n_profile
    IF (cos_zenith_angle(l) > 0.0) THEN
      bound%solar_irrad(l)=solar_irrad(l)
      bound%zen_0(l)=1.0/cos_zenith_angle(l)
    ELSE
      bound%solar_irrad(l)=0.0
      bound%zen_0(l)=1.0
    END IF
    IF (control%l_orog) THEN
      bound%orog_corr(l)=orog_corr(l)
    END IF
  END DO
END IF

END SUBROUTINE set_bound
END MODULE set_bound_mod
