! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE socrates_set_aer
IMPLICIT NONE
CONTAINS

! Subroutine to set the aerosol fields for the core radiation code.
!------------------------------------------------------------------------------
SUBROUTINE set_aer(control, dimen, spectrum, aer, n_profile, n_layer, dust)

USE rad_pcf,      ONLY: ip_aersrc_classic_ron
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_aer,      ONLY: StrAer, allocate_aer, allocate_aer_prsc
use soc_constants_mod,   only: i_def, r_def

IMPLICIT NONE

integer :: i_aer, i, l

! Control options:
TYPE(StrCtrl),      INTENT(IN)  :: control

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! Aerosol properties:
TYPE(StrAer),       INTENT(OUT) :: aer

INTEGER(i_def), INTENT(IN) :: n_profile
!   Number of atmospheric profiles for radiation calculations
INTEGER(i_def), INTENT(IN) :: n_layer
!   Number of atmospheric layers for radiation calculations

REAL(r_def), INTENT(IN) :: dust(n_profile, n_layer)
!   Dust mass mixing ratio

! Allocate structure for the core radiation code interface
CALL allocate_aer(aer, dimen, spectrum)
CALL allocate_aer_prsc(aer, dimen, spectrum)

! Only touch mr_source/mix_ratio when aerosol radiative effects are actually
! switched on (control%l_aerosol, set from do_dust_forcing in read_control.F90)
! - for ordinary runs aer%mr_source must be left exactly as
! allocate_aer/allocate_aer_prsc set it, matching the upstream set_aer_mod's
! behaviour. This can NOT be gated on spectrum%aerosol%n_aerosol > 0: the
! standard ga7 spectral files also define aerosol species blocks (inherited
! from the Met Office UM) despite Isca never wiring up their parametrizations,
! so that check is true even for non-Mars, non-dust runs.
IF (control%l_aerosol) THEN
  aer%mr_source = ip_aersrc_classic_ron

  DO i_aer=1, spectrum%aerosol%n_aerosol
    aer%mr_type_index(i_aer)=i_aer
    DO i=1, n_layer
      DO l=1, n_profile
        aer%mix_ratio(l,i,i_aer) = dust(l,i)
      END DO
    END DO
  END DO
END IF

END SUBROUTINE set_aer
END MODULE socrates_set_aer
