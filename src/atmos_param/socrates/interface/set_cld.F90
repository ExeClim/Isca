! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE set_cld_mod
IMPLICIT NONE
CONTAINS

! Subroutine to set the cloud fields for the core radiation code.
!------------------------------------------------------------------------------
SUBROUTINE set_cld(control, dimen, spectrum, cld, n_profile)

USE rad_pcf
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_cld,      ONLY: StrCld, allocate_cld, allocate_cld_prsc
use soc_constants_mod,   only: i_def, r_def

IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(IN)  :: control

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! Cloud properties:
TYPE(StrCld),       INTENT(OUT) :: cld

INTEGER(i_def), INTENT(IN) :: n_profile
!   Number of atmospheric profiles for radiation calculations


! Allocate structure for the core radiation code interface
CALL allocate_cld(cld, dimen, spectrum)
CALL allocate_cld_prsc(cld, dimen, spectrum)

END SUBROUTINE set_cld
END MODULE set_cld_mod
