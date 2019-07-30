! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE set_control_mod
IMPLICIT NONE
CONTAINS

! Subroutine to set per-timestep options for the core radiation code
!------------------------------------------------------------------------------
SUBROUTINE set_control(control)

USE def_control,  ONLY: StrCtrl, allocate_control

IMPLICIT NONE


! Control options:
TYPE(StrCtrl), INTENT(INOUT) :: control


! Set diagnostic flags (hardwire to false for now)
control%l_flux_up_band            = .FALSE.!.TRUE.
control%l_flux_down_band          = .FALSE.!.TRUE.
control%l_flux_up_clear_band      = .TRUE.!.FALSE.
control%l_flux_down_clear_band    = .TRUE.!.FALSE.

control%l_clear = .TRUE.
control%l_clear = control%l_clear &
             .OR. control%l_flux_up_clear_band &
             .OR. control%l_flux_down_clear_band

END SUBROUTINE set_control
END MODULE set_control_mod
