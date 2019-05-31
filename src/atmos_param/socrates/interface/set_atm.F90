! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE set_atm_mod
IMPLICIT NONE
CONTAINS

! Subroutine to set the input atmospheric profiles for the core radiation code
!------------------------------------------------------------------------------
SUBROUTINE set_atm(control, dimen, spectrum, atm, n_profile, n_layer,          &
  p_layer, t_layer, t_layer_boundaries, d_mass, density, h2o, o3, co2)

USE rad_pcf
USE def_control,  ONLY: StrCtrl
USE def_dimen,    ONLY: StrDim
USE def_spectrum, ONLY: StrSpecData
USE def_atm,      ONLY: StrAtm, allocate_atm
USE socrates_config_mod, only: co_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, o2_mix_ratio, so2_mix_ratio, cfc11_mix_ratio, cfc12_mix_ratio, cfc113_mix_ratio, hcfc22_mix_ratio, hfc134a_mix_ratio
USE gas_list_pcf, ONLY: ip_h2o, ip_co2, ip_o3, ip_n2o, ip_ch4, ip_o2, ip_so2,  &
  ip_cfc11, ip_cfc12, ip_cfc113, ip_hcfc22, ip_hfc134a, ip_co

use soc_constants_mod, only: i_def, r_def

IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(IN)    :: control

! Dimensions:
TYPE(StrDim),       INTENT(IN)    :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)    :: spectrum

! Atmospheric properties:
TYPE(StrAtm),       INTENT(INOUT) :: atm

INTEGER(i_def), INTENT(IN) :: n_profile
!   Number of atmospheric profiles for radiation calculations
INTEGER(i_def), INTENT(IN) :: n_layer
!   Number of atmospheric layers for radiation calculations

REAL(r_def), INTENT(IN) :: p_layer(n_profile, n_layer)
!   Pressure at layer centres
REAL(r_def), INTENT(IN) :: t_layer(n_profile, n_layer)
!   Temperature at layer centres
REAL(r_def), INTENT(IN) :: t_layer_boundaries(n_profile, 0:n_layer)
!   Temperature at layer boundaries
REAL(r_def), INTENT(IN) :: d_mass(n_profile, n_layer)
!   Mass of layer (kg m-2)
REAL(r_def), INTENT(IN) :: density(n_profile, n_layer)
!   Density of layer (kg m-3)
REAL(r_def), INTENT(IN) :: h2o(n_profile, n_layer)
!   Mass mixing ratio of water vapour
REAL(r_def), INTENT(IN) :: o3(n_profile, n_layer)
!   Mass mixing ratio of ozone
REAL(r_def), INTENT(IN) :: co2(n_profile, n_layer)
!   Mass mixing ratio of co2

! Local variables.
INTEGER :: i, l, i_gas
!   Loop variables


CALL allocate_atm(atm, dimen, spectrum)

! Set up atmosphere grid
atm%n_profile = n_profile
atm%n_layer   = n_layer

! Set the pressures, temperatures, masses (per square metre) and densities
DO i=0, n_layer
  DO l=1, n_profile
    atm%t_level(l, i) = t_layer_boundaries(l, i)
  END DO
END DO
DO i=1, n_layer
  DO l=1, n_profile
    atm%p(l, i)       = p_layer(l, i)
    atm%t(l, i)       = t_layer(l, i)
    atm%mass(l, i)    = d_mass(l, i)
    atm%density(l, i) = density(l, i)
  END DO
END DO

! Set gas mass mixing ratios
DO i_gas=1, spectrum%gas%n_absorb
  SELECT CASE(spectrum%gas%type_absorb(i_gas))

  CASE(ip_co)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = co_mix_ratio
      END DO
    END DO
    
  CASE(ip_co2)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = co2(l,i)
      END DO
    END DO

  CASE(ip_h2o)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = h2o(l,i)
      END DO
    END DO

  CASE(ip_o3)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = o3(l,i)
      END DO
    END DO

  CASE(ip_n2o)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = n2o_mix_ratio
      END DO
    END DO

  CASE(ip_ch4)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = ch4_mix_ratio
      END DO
    END DO

  CASE(ip_o2)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = o2_mix_ratio
      END DO
    END DO
    
  CASE(ip_so2)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = so2_mix_ratio
      END DO
    END DO    

  CASE(ip_cfc11)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = cfc11_mix_ratio
      END DO
    END DO   

  CASE(ip_cfc12)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = cfc12_mix_ratio
      END DO
    END DO   

  CASE(ip_cfc113)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = cfc113_mix_ratio
      END DO
    END DO   

  CASE(ip_hcfc22)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = hcfc22_mix_ratio
      END DO
    END DO   

  CASE(ip_hfc134a)
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = hfc134a_mix_ratio
      END DO
    END DO 
      
  CASE DEFAULT
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = 0.0_r_def
      END DO
    END DO
  END SELECT
END DO

END SUBROUTINE set_atm
END MODULE set_atm_mod
