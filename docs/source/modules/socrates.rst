
SOCRATES Radiation Scheme Interface
===================================

Summary
-------
.. This summary is modified from Stephen Thomson's P/R for Socrates: <https://github.com/ExeClim/Isca/pull/61>

SOCRATES (Suite Of Community RAdiative Transfer codes based on Edwards and Slingo) is the radiation scheme used by UK Met Office for Earth and planetary science [MannersEtAl2015]_, which has many significant advantages over RRTM, notably its flexibility in terms of atmospheric composition and the spectral properties of the radiation scheme (e.g. number of bands, etc).

* The code used to integrate Socrates into Isca is contained within the folder ``src/atmos_params/socrates/interface``.
* The Socrates source code itself is **NOT** packed within this Isca repository, and **NEW** users will need to download it from the `Met Office Science Repository <https://code.metoffice.gov.uk/trac/socrates>`_, and put it under ``src/atmos_params/socrates/src/trunk``. Detailed instructions on how to do this are included in the `README.md <https://github.com/ExeClim/Isca/blob/master/exp/test_cases/socrates_test/README.md>`_ for the Socrates test-case.
* The basis of ``socrates_interface`` was coded by Mark Hammond (Univ. of Oxford) and James Manners (Met Office) and modified by Stephen Thomson (Univ. of Exeter) [Thomson_and_Vallis2019]_. Features added include seasonality in the radiation based on Isca's ``astronomy`` package, and the ability to use a radiation timestep != atmospheric timestep.
* Socrates radiation scheme requires ``mass-mixing`` ratios (mmr) for all quantities (e.g. CO2, water vapour etc). This contrasts with RRTM, which wants ``volume-mixing`` ratios (vmr).


Namelist options
---------------- 

The Socrates namelist ``socrates_rad_nml`` can be found at ``src/atmos_param/socrates/interface/socrates_config_mod.f90``.

Radiation options
^^^^^^^^^^^^^^^^^

+----------------------------+---------------+-----------------------------------------------------------------------------------------+
| Name                       | Default value | Description                                                                             |
+============================+===============+=========================================================================================+
|``solday``                  | 0             | If >0, do perpetual run corresponding to day of the year = solday in [0, days per year] |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+
|``do_rad_time_avg``         | True          | Average ``coszen`` for shortwave radiation over ``dt_rad``                              |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+
|``equinox_day``             | 0.75          | fraction of the year defining NH autumn equinox in [0, 1]                               |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+
|``stellar_constant``        | 1368.22       | :math:`Wm^{-2}`, solar constant, consistent with RRTM default                           |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+
|``tidally_locked``          | False         | Tidally locked or not                                                                   |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+
|``frierson_solar_rad``      | False         | Options for annual-mean incoming solar, as prescribed in Frierson's grey scheme         |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+
|``del_sol``                 | 1.4           |                                                                                         |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+
|``del_sw``                  | 0             |                                                                                         |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+
|``input_planet_emissivity`` | 1.0           | Emissivity of surface. Defined as constant all over surface.                            |
+----------------------------+---------------+-----------------------------------------------------------------------------------------+

The following namelist variables set radiation time stepping and spatial sampling:

+----------------------------+---------------+--------------------------------------------------------------------------+
| Name                       | Default value | Description                                                              |
+============================+===============+==========================================================================+
| ``dt_rad``                 | 0             | Radiation timestep - every step if ``dt_rad<dt_atmos``                   |
+----------------------------+---------------+--------------------------------------------------------------------------+
| ``store_intermediate_rad`` | True          | Keep rad constant over entire dt_rad?                                    |
+----------------------------+---------------+--------------------------------------------------------------------------+
| ``dt_rad_avg``             | -1            | If averaging, over what time? ``dt_rad_avg=dt_rad`` if ``dt_rad_avg<=0`` |
+----------------------------+---------------+--------------------------------------------------------------------------+
| ``chunk_size``             | 16            | Number of gridpoints to pass to socrates at a time                       |
+----------------------------+---------------+--------------------------------------------------------------------------+


Spectral files
^^^^^^^^^^^^^^

+--------------------------------+---------------+-----------------------------------------------------+
| Name                           | Default value | Description                                         |
+================================+===============+=====================================================+
| ``socrates_hires_mode``        | False         | If false then run in 'GCM mode', and                |
|                                |               | if true then usehigh-res spectral file              |
+--------------------------------+---------------+-----------------------------------------------------+
| ``lw_hires_spectral_filename`` | 'unset'       |                                                     |
+--------------------------------+---------------+-----------------------------------------------------+
| ``sw_hires_spectral_filename`` | 'unset'       |                                                     |
+--------------------------------+---------------+-----------------------------------------------------+
| ``lw_spectral_filename``       | 'unset'       | Longwave spectral file, which can be found at       |
|                                |               | Socrates source code ``data`` directory (e.g.       |
|                                |               | ``data/spectra/ga7/sp_lw_ga7``) or self-generated   |
+--------------------------------+---------------+-----------------------------------------------------+
| ``sw_spectral_filename``       | 'unset'       | Shortwave spectral file, which can be found at      |
|                                |               | Socrates source code ``data`` directory (e.g.       |
|                                |               | ``data/spectra/ga7/sp_sw_ga7``) or self-generated   |
+--------------------------------+---------------+-----------------------------------------------------+

CO2 and ozone
^^^^^^^^^^^^^

+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| Name                              | Default value | Description                                                                 |
+===================================+===============+=============================================================================+
| ``account_for_effect_of_water``   | True          | - False: radiation is fed water mixing ratios = 0                           |
|                                   |               | - True:  radiation is fed mixing ratios based on model specific humidity.   |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``account_for_effect_of_ozone``   | True          | - False: radiation is fed ozone mixing ratios = 0                           |
|                                   |               | - True:  radiation is fed mixing ratios based on model ozone field          |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``do_read_co2``                   | False         | - Read CO2 from an external file?                                           |
|                                   |               | - If true, needs to specify CO2 file and variable names                     |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``co2_file_name``                 | 'co2'         | Name of file containing CO2 field - n.b. don't need to include '.nc'        |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``co2_field_name``                | 'co2'         | Name of CO2 variable in CO2 file                                            |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``input_co2_mmr``                 | False         | - True if the input file contain values as ``mass mixing ratio``            |
|                                   |               | - False if the input file contain values as ``volume mixing ratio``         |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``co2_ppmv``                      | 300           | Default CO2 concentration in ``ppmv``                                       |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``do_read_ozone``                 | False         | - Read ozone from an external file?                                         |
|                                   |               | - If true, needs to specify ozone file and variable names                   |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``ozone_file_name``               | 'ozone'       | Name of file containing ozone field - n.b. don't need to include '.nc'      |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``ozone_field_name``              | 'ozone'       | Name of ozone variable in ozone file                                        |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``input_o3_file_is_mmr``          | True          | - ``True`` if the input file contain values as ``mass mixing ratio``        |
|                                   |               | - ``False`` if the input file contain values as ``volume mixing ratio``     |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+

Diagnostics
-----------
.. What diagnostics are available for this part of the code.


Relevant modules and subroutines
--------------------------------
.. List the names of relevant modules, subroutines, functions, etc.
.. You can add also code snippets, using Sphinx code formatting


References
----------
[MannersEtAl2015]_
[Thomson_and_Vallis2019]_

Authors
-------
This documentation was written by Qun Liu, peer reviewed by X, and quality controlled by X.
