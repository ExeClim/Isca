.. _examples:

Examples
========

Test Cases
----------

Isca provides a number of example test cases, located within ``src/test_cases``, which can be used as a basis for constructing your own model. These examples vary in complexity and a number make use of additional scripts which are also listed below.

- **Axisymmetric:** Model setup with zonal-mean dynamics enforced, i.e. a run without eddies. Uses prescribed SSTs and vertical diffusion.
- **Bucket Hydrology:** As described in Isca paper (Vallis et al., 2017) but without q-fluxes
- **Frierson:** Control case of the so-called `Frierson model` described in e.g. <https://doi.org/10.1175/JAS3753.1>
- **Giant Planet:** Control case of the Jupiter simulation in <https://doi.org/10.1175/2008JAS2798.1>
- **Held Suarez:** Commonly-used dynamical-core test case described [here][1]
- **MiMA:** CNTL case of the MiMA model described in <https://doi.org/10.1175/JCLI-D-17-0127.1>
- **Realistic Continents:** Most complex configuration of Isca v1.0, including realistic continents, realistic topography, a simple ice model, and the option of prescribed AMIP SSTs or prescribed q-fluxes derived from AMIP ssts. Same setup used in Thomson & Vallis, 2017 (submitted), also described in Isca paper.
- **Top Down:** A general thermal relaxation scheme allowing seasonal variation, dependent on atmospheric and orbital parameters.
- **Variable CO2 Concentration:** Example experiments with a time-varying co2 concentration read-in from an input file. Examples are given with varying co2 in the [Byrne & O'Gorman radiation scheme](https://doi.org/10.1175/JCLI-D-12-00262.1) and RRTM.

Auxiliary Scripts
-----------------

A number of these scripts can be used to help construct input files for a model, whilst others may be used for verification or post-processing of model data.

- Create land-mask: ``src/extra/python/gfdl/land_generator_fn.py``
- Make SST climatology: ``src/extra/python/scripts/create_amip_sst_timeseries.py``
- Make time-varying q-fluxes: ``src/extra/python/scripts/calculate_qflux/calculate_qflux.py``
- Make CO2 concentration input file: ``src/extra/python/scripts/create_co2_timeseries.py``
- Interpolate from model sigma levels to pressure levels: ``postprocessing/plevel_interpolation/scripts/run_plevel.py``
- Explore options for vertical co-ordinates: ``src/extra/python/scripts/vert_coord_options.py``
- Display options for horizontal resolution: ``src/extra/python/scripts/resolutions.py``
- Change horizontal resolution part-way through a run: ``src/extra/python/scripts/change_horizontal_resolution_of_restart_file.py``
- Check if simulation is spun-up: ``src/extra/python/scripts/general_spinup_fn.py``
- Alter input file so model's linear interpolator preserves monthly-means: ``src/extra/python/scripts/edit_nc_file_to_preserve_monthly_means.py``
- Plot how long a run is taking to calculate each month: ``src/extra/python/scripts/modified_time_script.py``

[1]: https://doi.org/10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2

.. note::

   These examples provide practical ways to illustrate the usage of the software, and do not necessarily represent the best scientific methods for any particular type of study.

Authors
-------
This documentation was written by Daniel Williams and reviewed by <name>




