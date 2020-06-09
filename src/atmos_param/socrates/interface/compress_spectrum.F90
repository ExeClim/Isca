! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE compress_spectrum_mod
IMPLICIT NONE
CONTAINS

! Subroutine to compress the spectral file data to only those required
!------------------------------------------------------------------------------
SUBROUTINE compress_spectrum(con, spec)

USE def_control,         ONLY: StrCtrl
USE def_spectrum,        ONLY: StrSpecData
USE gas_list_pcf
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE


TYPE (StrCtrl),     INTENT(IN)    :: con
TYPE (StrSpecData), INTENT(INOUT) :: spec

! Local variables
INTEGER :: i, j, n_band_absorb
LOGICAL :: l_retain_absorb(spec%gas%n_absorb)
!   Flags for the retention of gases in the spectral file

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COMPRESS_SPECTRUM'


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Search the spectrum to find those gases to be retained.
l_retain_absorb=.FALSE.
DO i=1, spec%gas%n_absorb
  IF (((spec%gas%type_absorb(i) == ip_h2o)     .AND. con%l_h2o    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_co2)     .AND. con%l_co2    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_o3)      .AND. con%l_o3     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_o2)      .AND. con%l_o2     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_n2o)     .AND. con%l_n2o    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_ch4)     .AND. con%l_ch4    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_so2)     .AND. con%l_so2    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cfc11)   .AND. con%l_cfc11  ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cfc12)   .AND. con%l_cfc12  ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cfc113)  .AND. con%l_cfc113 ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cfc114)  .AND. con%l_cfc114 ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_hcfc22)  .AND. con%l_hcfc22 ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_hfc125)  .AND. con%l_hfc125 ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_hfc134a) .AND. con%l_hfc134a) .OR.       &
      ((spec%gas%type_absorb(i) == ip_co)      .AND. con%l_co     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_nh3)     .AND. con%l_nh3    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_tio)     .AND. con%l_tio    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_vo)      .AND. con%l_vo     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_h2)      .AND. con%l_h2     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_he)      .AND. con%l_he     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_na)      .AND. con%l_na     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_k)       .AND. con%l_k      ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_li)      .AND. con%l_li     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_rb)      .AND. con%l_rb     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cs)      .AND. con%l_cs     )) THEN
    l_retain_absorb(i)=.TRUE.
  END IF
END DO

DO i=1, spec%basic%n_band
  n_band_absorb=0
  DO j=1, spec%gas%n_band_absorb(i)
    IF (l_retain_absorb(spec%gas%index_absorb(j, i))) THEN
      n_band_absorb = n_band_absorb + 1
      spec%gas%index_absorb(n_band_absorb, i) = spec%gas%index_absorb(j, i)
    END IF
  END DO
  spec%gas%n_band_absorb(i)=n_band_absorb
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE compress_spectrum
END MODULE compress_spectrum_mod
