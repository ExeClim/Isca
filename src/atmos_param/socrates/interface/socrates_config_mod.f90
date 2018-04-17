module socrates_config_mod

use def_control,  only: StrCtrl
use def_spectrum, only: StrSpecData

use soc_constants_mod, only: r_def

implicit none


TYPE (StrCtrl), SAVE :: sw_control, lw_control
TYPE (StrSpecData), SAVE :: sw_spectrum, lw_spectrum

! Set a grey surface albedo
LOGICAL :: l_planet_grey_surface = .TRUE.
REAL(r_def) :: planet_albedo = 0.06
REAL(r_def) :: planet_emissivity = 0.97

! Well mixed gas concentrations (kg / kg)
REAL(r_def) :: co2_mix_ratio = 0.0!2.E-6!!0.01!1e-2!6.002e-5!4
REAL(r_def) :: co_mix_ratio = 1.0!2.E-6!1.0!1.0!1.E-4!0.001!6.002e-2!4
!REAL(r_def) :: n2o_mix_ratio = 0.0!4.945e-07
!REAL(r_def) :: ch4_mix_ratio = 0.0!1.006e-06
!REAL(r_def) :: o2_mix_ratio = 0.2314
!REAL(r_def) :: so2_mix_ratio = 0.0
!REAL(r_def) :: cfc11_mix_ratio = 1.110e-09
!REAL(r_def) :: cfc12_mix_ratio = 2.187e-09
!REAL(r_def) :: cfc113_mix_ratio = 4.826e-10
!REAL(r_def) :: hcfc22_mix_ratio = 6.866e-10
!REAL(r_def) :: hfc134a_mix_ratio = 2.536e-10
end module socrates_config_mod
