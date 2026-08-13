# Mars dust + latent-heating input files

`socrates_mars_dust_test_case.py` requires the following files to be placed in this
folder before it can be run:

* `sp_lw_17_dsa_mars_dust`
* `sp_lw_17_dsa_mars_dust_k`
* `sp_sw_42_dsa_mars_sun_dust`
* `sp_sw_42_dsa_mars_sun_dust_k`
* `t42_mola_mars.nc`
* `cdod_clim_MY24.nc` .. `cdod_clim_MY34.nc` (one file per Mars Year, 24-34)
* `cdod_clim_scenario.nc`
* `cdod_cold.nc`
* `cdod_warm_25.nc`
* `cdod_all_years.nc`
* `cdod_all_years_long.nc`

These are not included in this repository (the dust-climatology files alone are around
100MB). They can be found in Emily Ball's own fork, which this test case is ported
from: <https://github.com/emilyrball/Isca-Mars/tree/bp1_lh_dev/exp/socrates_mars/input>.

## What these files are

* The `sp_*_dsa_mars_dust*` files are Mars-specific Socrates spectral files with dust
  optical properties included, distinct from the plain (non-dust) Mars spectral files
  used by the `socrates_mars` test case.
* `t42_mola_mars.nc` is MOLA-derived Mars topography, also used by `socrates_mars`.
* The `cdod_*.nc` files are column dust optical depth climatologies, used to set the
  reference dust mass mixing ratio via `do_read_cdod`/`cdod_file_name` in
  `socrates_rad_nml`. `cdod_field_name` selects which variable within the chosen file
  to read.

## Physics this exercises

This test case switches on two capabilities from Emily Ball's PhD work (Ball et al.
2021, "The roles of latent heating and dust in the structure and variability of the
northern Martian polar vortex"), both gated behind their own namelist flags and off by
default everywhere else in Isca:

* `do_lscale_cond_lh` (`idealized_moist_phys_nml`) - a CO2-condensation latent-heating
  condensation scheme, an alternative to the standard `do_lscale_cond`.
* `do_dust_forcing` (`socrates_rad_nml`) - builds a Conrath-type Mars dust vertical
  profile from the dust optical depth climatology and feeds it to Socrates as a
  radiatively-active aerosol. `dust_scale` controls the optical-depth-to-mass-mixing-
  ratio conversion.
