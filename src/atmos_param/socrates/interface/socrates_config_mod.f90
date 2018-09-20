module socrates_config_mod

use def_control,  only: StrCtrl
use def_spectrum, only: StrSpecData

use soc_constants_mod, only: r_def

implicit none


TYPE (StrCtrl), SAVE :: sw_control, lw_control
TYPE (StrSpecData), SAVE :: sw_spectrum, lw_spectrum

!  Set a grey surface albedo
LOGICAL :: l_planet_grey_surface = .TRUE.

  ! Socrates inputs from namelist
  
  REAL :: stellar_constant = 1368.22 !Consistent with RRTM default
  LOGICAL :: tidally_locked = .FALSE.
  !Options for annual-mean incoming solar, as prescribed in Frierson's grey scheme.
  LOGICAL :: frierson_solar_rad = .FALSE.
  REAL    :: del_sol         = 1.4
  REAL    :: del_sw          = 0.0
  LOGICAL :: socrates_hires_mode = .FALSE.  !If false then run in 'GCM mode', if true then uses high-res spectral files
  character(len=256) :: lw_spectral_filename='unset'
  character(len=256) :: lw_hires_spectral_filename='unset'
  character(len=256) :: sw_spectral_filename='unset'
  character(len=256) :: sw_hires_spectral_filename='unset'
  logical :: account_for_effect_of_water=.TRUE. !if False then radiation is fed water mixing ratios = 0. If true it's fed mixing ratios based on model specific humidity.
  logical :: account_for_effect_of_ozone=.TRUE. !if False then radiation is fed ozone mixing ratios = 0. If true it's fed mixing ratios based on model ozone field.
  logical :: do_read_ozone = .FALSE. ! read ozone from an external file?
  character(len=256) :: ozone_file_name='ozone' !Name of file containing ozone field - n.b. don't need to include '.nc'
  character(len=256) :: ozone_field_name='ozone' !Name of ozone variable in ozone file
  logical            :: input_o3_file_is_mmr=.true. ! Does the ozone input file contain values as a mass mixing ratio (set to true) or a volume mixing ratio (set to false)?
  logical :: do_read_co2 = .FALSE. ! read ozone from an external file?
  character(len=256) :: co2_file_name='co2' !Name of file containing co2 field - n.b. don't need to include '.nc'
  character(len=256) :: co2_field_name='co2' !Name of co2 variable in co2 file  
  real(r_def) :: input_planet_emissivity = 1.0 !Emissivity of surface. Defined as constant all over surface.
  real :: co2_ppmv = 300. !Default CO2 concentration in PPMV
  logical ::  input_co2_mmr=.false. !Socrates wants input concentrations as mmr not vmr, so need to make sure input data supplied is converted if necessary

  logical :: use_pressure_interp_for_half_levels = .False. !By default (.False.) does linear interpolation in height for half-level temperatures. True does linear interp using pressure. 
      
  ! Incoming radiation options for namelist
  
  integer   :: solday=0  ! if >0, do perpetual run corresponding to day of the year = solday \in [0,days per year]
  logical   :: do_rad_time_avg = .true. ! Average coszen for SW radiation over dt_rad?
  real      :: equinox_day=0.75                ! fraction of the year defining NH autumn equinox \in [0,1]

  ! radiation time stepping and spatial sampling
  integer   :: dt_rad=0                        ! Radiation timestep - every step if dt_rad<dt_atmos
  logical   :: store_intermediate_rad =.true.  ! Keep rad constant over entire dt_rad?
  integer   :: dt_rad_avg = -1                 ! If averaging, over what time? dt_rad_avg=dt_rad if dt_rad_avg<=0

  integer   :: chunk_size = 16 !number of gridpoints to pass to socrates at a time


  ! Well mixed gas concentrations (kg / kg) #Don't know the source of these numbers. Need to check them. e.g. co mix ratio.
  REAL(r_def) :: co_mix_ratio = 0.0
  REAL(r_def) :: n2o_mix_ratio = 4.945e-07
  REAL(r_def) :: ch4_mix_ratio = 1.006e-06
  REAL(r_def) :: o2_mix_ratio = 0.2314
  REAL(r_def) :: so2_mix_ratio = 0.0
  REAL(r_def) :: cfc11_mix_ratio = 1.110e-09
  REAL(r_def) :: cfc12_mix_ratio = 2.187e-09
  REAL(r_def) :: cfc113_mix_ratio = 4.826e-10
  REAL(r_def) :: hcfc22_mix_ratio = 6.866e-10
  REAL(r_def) :: hfc134a_mix_ratio = 2.536e-10


  ! Whether to include radiative effects of particular gases
  
  logical   :: inc_h2o = .TRUE.
!   control%l_h2o            = .TRUE.

  logical   :: inc_co2 = .TRUE.
!  control%l_co2            = .TRUE.
  
  logical   :: inc_co  = .FALSE. 
!   control%l_co             = .FALSE.

  logical   :: inc_o3  = .TRUE.
!  control%l_o3             = .TRUE.

  logical   :: inc_n2o  = .FALSE.
!   control%l_n2o            = .TRUE.

  logical   :: inc_ch4  = .FALSE.
!   control%l_ch4            = .TRUE.

  logical   :: inc_o2   = .FALSE.
!   control%l_o2             = .TRUE.

  logical   :: inc_so2  = .FALSE.
!   control%l_so2            = .FALSE.

  logical   :: inc_cfc11= .FALSE.
!   control%l_cfc11          = .FALSE.

  logical   :: inc_cfc12= .FALSE.
!   control%l_cfc12          = .FALSE.

  logical   :: inc_cfc113=.FALSE.
!   control%l_cfc113         = .FALSE.

  logical   :: inc_hcfc22=.FALSE.
!   control%l_hcfc22         = .FALSE.

  logical   :: inc_hfc134a=.FALSE.
!   control%l_hfc134a        = .FALSE.


  NAMELIST/socrates_rad_nml/ stellar_constant, tidally_locked, lw_spectral_filename, lw_hires_spectral_filename, &
                             sw_spectral_filename, sw_hires_spectral_filename, socrates_hires_mode, &
                             input_planet_emissivity, co2_ppmv, &
                             account_for_effect_of_water, account_for_effect_of_ozone, &
                             do_read_ozone, ozone_file_name, ozone_field_name, input_o3_file_is_mmr, &
                             do_read_co2, co2_file_name, co2_field_name, input_co2_mmr, &                             
                             solday, do_rad_time_avg, equinox_day,  &
                             store_intermediate_rad, dt_rad_avg, dt_rad, &
                             chunk_size, &
                             co_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, &
                             o2_mix_ratio, so2_mix_ratio, cfc11_mix_ratio, &
                             cfc12_mix_ratio, cfc113_mix_ratio, hcfc22_mix_ratio, &
                             hfc134a_mix_ratio, &
                             inc_h2o, inc_co2, inc_co, inc_o3, inc_n2o, inc_ch4, inc_o2, &
                             inc_so2, inc_cfc11, inc_cfc12, inc_cfc113, inc_hcfc22, inc_hfc134a, &
                             use_pressure_interp_for_half_levels,  &
                             frierson_solar_rad, del_sol, del_sw

end module socrates_config_mod
