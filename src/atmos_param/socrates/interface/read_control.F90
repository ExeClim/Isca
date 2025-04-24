! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE read_control_mod
IMPLICIT NONE
CONTAINS

! Subroutine to set input algorithmic options for the core radiation code
!------------------------------------------------------------------------------
SUBROUTINE read_control(control, spectrum, do_clouds)

USE rad_pcf
USE def_control,  ONLY: StrCtrl, allocate_control
USE def_spectrum, ONLY: StrSpecData
USE socrates_config_mod, ONLY: l_planet_grey_surface, inc_h2o, inc_co2, inc_co,      & 
                               inc_o3, inc_n2o, inc_ch4, inc_o2, inc_so2, inc_cfc11, &
                               inc_cfc12, inc_cfc113, inc_hcfc22, inc_hfc134a

IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(INOUT) :: control

! Spectral data:
TYPE (StrSpecData), INTENT(IN)    :: spectrum

LOGICAL, INTENT(IN), OPTIONAL     :: do_clouds


! Local variables.
INTEGER :: i
!   Loop variable


! Eventually these options will be read from a namelist - hardwired for now

! Spectral bands
control%first_band = 1
control%last_band  = spectrum%basic%n_band

! Spectral region-specific options
select case(control%isolir)
case(ip_solar)
  control%i_2stream        = ip_pifm80 ! 16
  control%i_scatter_method = ip_scatter_full ! 1
  control%l_rayleigh       = .TRUE.
  control%l_orog           = .FALSE.
  control%l_solvar         = .FALSE.
  control%l_h2o            = inc_h2o
  control%l_co2            = inc_co2
  control%l_co             = inc_co
  control%l_o3             = inc_o3
  control%l_n2o            = inc_n2o
  control%l_ch4            = inc_ch4
  control%l_o2             = inc_o2
  control%l_so2            = inc_so2
  control%l_cfc11          = inc_cfc11
  control%l_cfc12          = inc_cfc12
  control%l_cfc113         = inc_cfc113
  control%l_hcfc22         = inc_hcfc22
  control%l_hfc134a        = inc_hfc134a
  control%i_st_water       = 5
  control%i_cnv_water      = 5
  control%i_st_ice         = 11
  control%i_cnv_ice        = 11
case(ip_infra_red)
  control%i_2stream        = ip_elsasser
  control%i_scatter_method = ip_scatter_full !ip_scatter_hybrid
  control%l_ir_source_quad = .TRUE.
  control%l_h2o            = inc_h2o
  control%l_co2            = inc_co2
  control%l_co             = inc_co
  control%l_o3             = inc_o3
  control%l_n2o            = inc_n2o
  control%l_ch4            = inc_ch4
  control%l_so2            = inc_so2
  control%l_cfc11          = inc_cfc11
  control%l_cfc12          = inc_cfc12
  control%l_cfc113         = inc_cfc113
  control%l_hcfc22         = inc_hcfc22
  control%l_hfc134a        = inc_hfc134a
  control%i_st_water       = 5
  control%i_cnv_water      = 5
  control%i_st_ice         = 11
  control%i_cnv_ice        = 11
end select

! Angular integration (including algorithmic options):
control%n_channel              = 1
control%i_angular_integration  = ip_two_stream
control%l_rescale              = .FALSE.
control%n_order_forward        = 2
control%l_mixing_ratio         = .TRUE.

! Gaseous absorption
control%l_gas          = .TRUE.
control%l_continuum    = .TRUE.
control%l_cont_gen     = .TRUE.
control%i_gas_overlap  = ip_overlap_k_eqv_scl!ip_overlap_random

! Properties of clouds
if (do_clouds) then
  control%i_cloud_representation = ip_cloud_ice_water
else
  control%i_cloud_representation = ip_cloud_off
end if
control%i_overlap              = ip_max_rand
control%i_inhom                = ip_homogeneous

IF (control%i_cloud_representation == ip_cloud_ice_water) THEN
  control%l_microphysics = .TRUE.
  control%l_cloud        = .TRUE.
  control%l_drop         = .TRUE.
  control%l_ice          = .TRUE.
  IF (control%i_inhom == ip_mcica) THEN
    control%i_cloud = ip_cloud_mcica
    IF ( (control%i_scatter_method == ip_no_scatter_abs) .OR.   &
         (control%i_scatter_method == ip_no_scatter_ext) ) THEN
      control%i_solver       = ip_solver_no_scat
      control%i_solver_clear = ip_solver_no_scat
    ELSE
      control%i_solver       = ip_solver_homogen_direct
      control%i_solver_clear = ip_solver_homogen_direct
    END IF
  ELSE
    IF (control%i_scatter_method == ip_scatter_approx) THEN
      control%i_solver       = ip_solver_mix_app_scat
      control%i_solver_clear = ip_solver_homogen_direct
    ELSE
      control%i_solver       = ip_solver_mix_direct_hogan
      control%i_solver_clear = ip_solver_homogen_direct
    END IF
    IF (control%i_overlap == ip_max_rand) THEN
      control%i_cloud = ip_cloud_mix_max
    ELSE IF (control%i_overlap == ip_exponential_rand) THEN
      control%i_cloud = ip_cloud_part_corr
    ELSE
      control%i_cloud = ip_cloud_mix_random
    END IF
  END IF
ELSE IF (control%i_cloud_representation == ip_cloud_csiw) THEN
  ! Not compatible with control%i_inhom == ip_mcica
  IF (control%i_inhom == ip_mcica) control%i_inhom = ip_homogeneous
  control%l_microphysics = .TRUE.
  control%l_cloud        = .TRUE.
  control%l_drop         = .TRUE.
  control%l_ice          = .TRUE.
  IF (control%i_scatter_method == ip_scatter_approx) THEN
    control%i_solver       = ip_solver_triple_app_scat
    control%i_solver_clear = ip_solver_homogen_direct
  ELSE
    control%i_solver       = ip_solver_triple_hogan
    control%i_solver_clear = ip_solver_homogen_direct
  END IF
  IF (control%i_overlap == ip_max_rand) THEN
    control%i_cloud = ip_cloud_triple
  ELSE IF (control%i_overlap == ip_exponential_rand) THEN
    control%i_cloud = ip_cloud_part_corr_cnv
  ELSE
    control%i_cloud = ip_cloud_mix_random
  END IF
ELSE
  ! No treatment of cloud
  control%l_microphysics = .FALSE.
  control%l_cloud        = .FALSE.
  control%l_drop         = .FALSE.
  control%l_ice          = .FALSE.
  control%i_cloud        = ip_cloud_clear
  IF ( (control%i_scatter_method == ip_no_scatter_abs) .OR.            &
       (control%i_scatter_method == ip_no_scatter_ext) ) THEN
    control%i_solver       = ip_solver_no_scat
    control%i_solver_clear = ip_solver_no_scat
  ELSE
    control%i_solver       = ip_solver_homogen_direct
    control%i_solver_clear = ip_solver_homogen_direct
  END IF
END IF

! Aerosols
control%l_aerosol      = .FALSE.
control%l_aerosol_mode = .FALSE.
control%l_aerosol_ccn  = .FALSE.

! Tiling is not needed if the surface is to be treated as grey.
IF (l_planet_grey_surface) THEN
  control%l_tile=.FALSE.
ELSE
  control%l_tile=.TRUE.
END IF

! Allocate band-by-band control options
CALL allocate_control(control, spectrum)

! Set properties for individual bands.
DO i = 1, spectrum%basic%n_band
  control%map_channel(i)           = 1
  control%weight_band(i)           = 1.0
  control%i_scatter_method_band(i) = control%i_scatter_method
  control%i_gas_overlap_band(i)    = control%i_gas_overlap
  IF (ANY(spectrum%gas%i_scale_fnc(i,:) == ip_scale_ses2)) THEN
    ! If SES2 scaling is used in this band then the overlap must also use SES2:
    control%i_gas_overlap_band(i)  = ip_overlap_mix_ses2
  END IF
END DO

END SUBROUTINE read_control
END MODULE read_control_mod
