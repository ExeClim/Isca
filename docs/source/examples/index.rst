.. _examples:

Examples
========

Test Cases
----------

Isca provides a number of example test cases, located within ``src/test_cases``, which can be used as a basis for constructing your own model. These examples vary in complexity and a number make use of additional scripts which are also listed below.

- **APE Aquaplanet:** Aquaplanet control experiment based off Blackburn et al. 2013. [Blackburn2013]_
- **Axisymmetric:** Model setup with zonal-mean dynamics enforced, i.e. a run without eddies. Uses prescribed SSTs and vertical diffusion.
- **Bucket Hydrology:** As described in the first Isca paper by Vallis et al., 2017  but without q-fluxes. [Vallis2017]_
- **Frierson:** Control case of the so-called `Frierson model` described in Frierson et al. 2006. [Frierson2006a]_
- **Giant Planet:** Control case of the Jupiter simulation defined by Schneider & Liu 2009. [Schneider2009]_
- **Held Suarez:** Commonly-used dynamical-core test case described by Held & Suarez 1994. [Held1994]_
- **MiMA:** CNTL case of the MiMA model described by Jucker & Gerber 2007. [Jucker2017]_
- **Realistic Continents:** Most complex configuration of Isca v1.0, including realistic continents, realistic topography, a simple ice model, and the option of prescribed AMIP SSTs or prescribed q-fluxes derived from AMIP ssts. Same setup used in Thomson & Vallis  2019; also described in the Isca paper. [Thomson_and_Vallis2019]_
- **SOCRATES:** A number of similar test cases all making use of the SOCRATES radiation scheme described in Manners et al. 2015. Example variants include a simple aquaplanet as well as realistic continents (as above) or with a simple cloud scheme. [MannersEtAl2015]_
- **Top Down:** A general thermal relaxation scheme allowing seasonal variation, dependent on atmospheric and orbital parameters.
- **Variable CO2 Concentration:** Example experiments with a time-varying co2 concentration read-in from an input file. Examples are given with varying co2 in the `Byrne & O'Gorman radiation scheme <https://doi.org/10.1175/JCLI-D-12-00262.1>`_ and RRTM. [Byrne2013]_

.. note::

   These examples provide practical ways to illustrate the usage of the software, and do not necessarily represent the best scientific methods for any particular type of study.


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


References
----------

| [Blackburn2013]_
| [Byrne2013]_
| [Frierson2006a]_
| [Held1994]_
| [Jucker2017]_
| [MannersEtAl2015]_
| [Schneider2009]_
| [Thomson_and_Vallis2019]_
| [Vallis2017]_


Authors
-------
This documentation was written by Daniel Williams and reviewed by <name>




