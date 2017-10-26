# List of test cases with paper references:

`MiMA`
* CNTL case of the MiMA model described in <https://doi.org/10.1175/JCLI-D-17-0127.1>

`bucket_hydrology`
* As described in Isca paper (Vallis et al., 2017) but without q-fluxes

`frierson`
* Control case of the so-called `Frierson model` described in e.g. <https://doi.org/10.1175/JAS3753.1>

`giant_planet`
* Control case of the Jupiter simulation in <https://doi.org/10.1175/2008JAS2798.1>

`axisymmetric`
* Model setup with zonal-mean dynamics enforced, i.e. a run without eddies. Uses prescribed SSTs and vertical diffusion.

`held_suarez`
* Commonly-used dynamical-core test case described [here][1]

`top_down`
* A general thermal relaxation scheme allowing seasonal variation, dependent on atmospheric and orbital parameters.

`realistic_continents`
* Most complex configuration of Isca v1.0, including realistic continents, realistic topography, a simple ice model, and the option of prescribed AMIP SSTs or prescribed q-fluxes derived from AMIP ssts. Same setup used in Thomson & Vallis, 2017 (submitted), also described in Isca paper.

`variable_co2_concentration`
* Example experiments with a time-varying co2 concentration read-in from an input file. Examples are given with varying co2 in the [Byrne & O'Gorman radiation scheme](https://doi.org/10.1175/JCLI-D-12-00262.1) and RRTM.

# List of Isca features and the location of the appropriate python scripts

Make land-mask:
`Isca/src/extra/python/gfdl/land_generator_fn.py`

Make SST climatology:
`Isca/src/extra/python/scripts/create_amip_sst_timeseries.py`

Make time-varying q-fluxes:
`Isca/src/extra/python/scripts/calculate_qflux/calculate_qflux.py`

Make CO2 concentration input file:
`Isca/src/extra/python/scripts/create_co2_timeseries.py`

Interpolate from model sigma levels to pressure levels:
`Isca/postprocessing/plevel_interpolation/scripts/run_plevel.py`

Explore options for vertical co-ordinates:
`Isca/src/extra/python/scripts/vert_coord_options.py`

Explore options for horizontal resolution:
`Isca/src/extra/python/scripts/resolutions.py`

Change horizontal resolution part-way through a run:
`Isca/src/extra/python/scripts/change_horizontal_resolution_of_restart_file.py`

Check if simulation is spun-up:
`Isca/src/extra/python/scripts/general_spinup_fn.py`

Alter input file so model's linear interpolator preserves monthly-means:
`Isca/src/extra/python/scripts/edit_nc_file_to_preserve_monthly_means.py`

Plot how long a run is taking to calculate each month:
`Isca/src/extra/python/scripts/modified_time_script.py`

[1]: https://doi.org/10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2
