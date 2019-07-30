! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
MODULE set_dimen_mod
IMPLICIT NONE
CONTAINS

! Subroutine to set dimensions for the radiation code.
!------------------------------------------------------------------------------
SUBROUTINE set_dimen(control, dimen, spectrum, n_profile, n_layer,             &
  n_cloud_layer, n_aer_mode, cld_subcol_gen, cld_subcol_req)

USE rad_pcf
USE def_control,  ONLY: StrCtrl
USE def_dimen,    ONLY: StrDim
USE def_spectrum, ONLY: StrSpecData

use soc_constants_mod, only: i_def

IMPLICIT NONE

! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Dimensions:
  TYPE(StrDim),       INTENT(INOUT) :: dimen

! Spectral data:
  TYPE (StrSpecData), INTENT(IN)    :: spectrum

INTEGER(i_def), INTENT(IN) :: n_profile
!   Number of columns to operate on
INTEGER(i_def), INTENT(IN) :: n_layer
!   Number of layers for radiation
INTEGER(i_def), INTENT(IN) :: n_cloud_layer
!   Number of potentially cloudy layers
INTEGER(i_def), INTENT(IN) :: n_aer_mode
!   Number of aerosol modes
INTEGER(i_def), INTENT(IN) :: cld_subcol_gen
!   Number of sub-grid cloud columns generated
INTEGER(i_def), INTENT(IN) :: cld_subcol_req
!   Number of sub-grid cloud columns required


dimen%nd_profile = n_profile
dimen%nd_layer = n_layer


! Cloud
dimen%nd_cloud_type      = 4
dimen%nd_cloud_component = 4
dimen%nd_cloud_representation = 4
SELECT CASE(control%i_cloud)
CASE (ip_cloud_column_max)
  dimen%nd_column = 3 * n_cloud_layer + 2
CASE DEFAULT
  dimen%nd_column = 1
END SELECT
dimen%nd_subcol_gen = cld_subcol_gen
dimen%nd_subcol_req = cld_subcol_req

dimen%id_cloud_top = dimen%nd_layer + 1 - n_cloud_layer
IF (control%l_cloud) THEN
  dimen%nd_layer_clr = dimen%id_cloud_top - 1
ELSE
  dimen%nd_layer_clr = dimen%nd_layer
END IF


! Aerosol
dimen%nd_aerosol_mode = MAX(1,n_aer_mode)

! Arrays for prescribed optical properties (not used here)
dimen%nd_profile_aerosol_prsc   = 1
dimen%nd_profile_cloud_prsc     = 1
dimen%nd_opt_level_aerosol_prsc = 1
dimen%nd_opt_level_cloud_prsc   = 1


! Tiled surface.
IF (control%l_tile) THEN
  dimen%nd_tile_type  = 3
  dimen%nd_point_tile = MAX(1,n_profile)
  dimen%nd_tile       = 3
ELSE
  dimen%nd_tile_type  = 1
  dimen%nd_point_tile = 1
  dimen%nd_tile       = 1
END IF

dimen%nd_viewing_level    = 1
dimen%nd_radiance_profile = 1
dimen%nd_j_profile        = 1
dimen%nd_direction        = 1
dimen%nd_brdf_basis_fnc   = 2
dimen%nd_brdf_trunc       = 1
dimen%nd_flux_profile     = dimen%nd_profile
dimen%nd_channel          = 1
dimen%nd_2sg_profile      = dimen%nd_profile
dimen%nd_source_coeff     = 2
dimen%nd_max_order        = 1
dimen%nd_sph_coeff        = 1


SELECT CASE(control%i_solver)
CASE(ip_solver_mix_app_scat, ip_solver_mix_direct, ip_solver_mix_direct_hogan)
  dimen%nd_overlap_coeff=8
  dimen%nd_region=2
CASE(ip_solver_triple_app_scat, ip_solver_triple, ip_solver_triple_hogan)
  dimen%nd_overlap_coeff=18
  dimen%nd_region=3
CASE DEFAULT
  dimen%nd_overlap_coeff=1
  dimen%nd_region=2
END SELECT


END SUBROUTINE set_dimen
END MODULE set_dimen_mod
