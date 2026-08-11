
SOCRATES Radiation Scheme Interface
===================================

Summary
-------
.. This summary is modified from Stephen Thomson's PR for Socrates: <https://github.com/ExeClim/Isca/pull/61>

SOCRATES (Suite Of Community RAdiative Transfer codes based on Edwards and Slingo) is the radiation scheme used by UK Met Office for Earth and planetary science [MannersEtAl2015]_, which has many significant advantages over RRTM, notably its flexibility in terms of atmospheric composition and the spectral properties of the radiation scheme (e.g. number of bands, etc).

* The code used to integrate Socrates into Isca is contained within the folder ``src/atmos_params/socrates/interface``.
* The Socrates source code itself is **NOT** packed within this Isca repository, but as of 2026 it is open source (BSD-3-Clause) and hosted at `github.com/MetOffice/socrates <https://github.com/MetOffice/socrates>`_. Isca depends on it via a git submodule pinned to a default, tested version (``src/atmos_param/socrates/src/<version>``, see ``.gitmodules``), and ``SocratesCodeBase``/``SocColumnCodeBase`` will fetch the matching source automatically the first time you compile if it isn't already present -- there is no manual download step for new users any more. See `Getting Socrates source code`_ below for how version selection, fetching and compilation work, and the `README.md <https://github.com/ExeClim/Isca/blob/master/exp/test_cases/socrates_test/README.md>`_ for the Socrates test-case (``exp/test_cases/socrates_test/README.md``) for a worked example.
* The basis of ``socrates_interface`` was coded by Mark Hammond (Univ. of Oxford) and James Manners (Met Office) and modified by Stephen Thomson (Univ. of Exeter) [Thomson_and_Vallis2019]_. Features added include seasonality in the radiation based on Isca's ``astronomy`` package, and the ability to use a ``radiation timestep != atmospheric timestep``.
* Socrates radiation scheme requires ``mass mixing ratios`` for all quantities (e.g. CO2, water vapor etc). This contrasts with RRTM, which wants ``volume mixing ratios``.


Getting Socrates source code
-----------------------------

Socrates' own directory layout and exact file list change slightly from release to release, so rather than hand-maintaining a list of files to compile for each version, Isca derives it automatically from Socrates' own ``make/Mk_src_*``/``Mk_mod_*`` build manifests plus a dependency scan starting from Isca's interface code (see ``isca.socrates_paths`` and ``src/extra/python/scripts/generate_socrates_path_names.py`` for the implementation). A handful of tested versions already have this list committed at ``src/extra/model/socrates/socrates_version_paths/<version>``; for any other version it is generated on the fly the first time you compile.

To compile with Socrates, pick a codebase and, optionally, a version::

    from isca import SocratesCodeBase
    cb = SocratesCodeBase.from_directory(GFDL_BASE, socrates_version='2026.07.1')  # version is optional; defaults to the vendored version
    cb.compile()

The first time this runs for a given version, Isca looks for the matching Socrates source in this order:

1. Already present at ``src/atmos_param/socrates/src/<version>`` (e.g. because you ran ``git submodule update --init`` for the default vendored version).
2. The legacy single-version override ``GFDL_SOC``, if set (points directly at a Socrates checkout, regardless of version -- kept for backwards compatibility with older setups).
3. ``$GFDL_SOC_DIR/<version>``, if ``GFDL_SOC_DIR`` is set and that directory exists.
4. Otherwise, Isca fetches it: a sparse, shallow clone of just Socrates' ``src/`` and ``make/`` directories plus the ``ga7``/``ga3_1`` spectral data Isca's own test cases use (``data/spectra/ga7``, ``data/spectra/ga3_1`` -- together under 15MB, not the ~300MB full repository with all of its example/data files) from ``github.com/MetOffice/socrates`` at the tag matching ``socrates_version``, into ``$GFDL_SOC_DIR/<version>`` (or a default cache under ``$GFDL_WORK`` if ``GFDL_SOC_DIR`` isn't set). If a test case needs a different spectral file, widen the sparse-checkout in that directory yourself (``git -C <dir> sparse-checkout add data/spectra/<whatever>``) or fetch it directly from the GitHub repository.

This means a machine with no pre-existing Socrates checkout and no environment variables set can still go straight from ``cb = SocratesCodeBase.from_directory(...)`` to a working compile, as long as it has network access to GitHub at that point (e.g. on a login node -- if compute nodes are offline, fetch the version you need in advance on a node that does have access, or pre-populate ``GFDL_SOC_DIR``).


Namelist options
---------------- 

The Socrates namelist ``socrates_rad_nml`` can be found at ``src/atmos_param/socrates/interface/socrates_config_mod.f90``.

Radiation options
^^^^^^^^^^^^^^^^^

Here are some options to set incoming radiation:

+----------------------------+----------+-----------------------------------------------------------------------------------------+
| Name                       | Default  | Description                                                                             |
+============================+==========+=========================================================================================+
|``solday``                  | 0        | If ``solday>0``, do perpetual run corresponding to day of the year = solday in          |
|                            |          | [0, days per year]                                                                      |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``do_rad_time_avg``         | True     | Average ``coszen`` for shortwave radiation over ``dt_rad``                              |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``equinox_day``             | 0.75     | Fraction of the year defining NH autumn equinox in ``[0, 1]``                           |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``stellar_constant``        | 1368.22  | solar constant (units: Wm :math:`^{-2}`), consistent with RRTM default                  |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``tidally_locked``          | False    | Tidally locked or not                                                                   |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``frierson_solar_rad``      | False    | Options for annual-mean incoming solar, as prescribed in Frierson's grey scheme         |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``del_sol``                 | 1.4      | Latitudinal variation of shortwave radiation, as prescribed in Frierson's grey scheme   |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``del_sw``                  | 0        | Latitudinal variation of shortwave radiation, as prescribed in Frierson's grey scheme   |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``input_planet_emissivity`` | 1.0      | Emissivity of surface. Defined as constant all over surface.                            |
+----------------------------+----------+-----------------------------------------------------------------------------------------+

The following namelist variables set radiation time stepping and spatial sampling:

+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| Name                       | Default  | Description                                                                                             |
+============================+==========+=========================================================================================================+
| ``dt_rad``                 | 0        | Radiation timestep - every step if ``dt_rad<dt_atmos``                                                  |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``store_intermediate_rad`` | True     | Keep rad constant over entire ``dt_rad``? (recommended)                                                 |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``dt_rad_avg``             | -1       | If averaging, over what time? ``dt_rad_avg=dt_rad`` if ``dt_rad_avg<=0``                                |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``chunk_size``             | 16       | Number of gridpoints to pass to Socrates at a time                                                      |
|                            |          | (see `README.md <https://github.com/ExeClim/Isca/blob/master/exp/test_cases/socrates_test/README.md>`_) |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+

Spectral files
^^^^^^^^^^^^^^

Socrates reads external input spectral files that tell it the number of spectral bands to use, with one file setting the shortwave options, and another file setting the longwave options. Some spectral files have lots of bands, which will make the model run slowly. The default files used in the Met Office's Unified Model-GA7, and also in Isca, can be found within the ``data/spectra`` directory of Socrates source code. For example, it can be here if you put the Socrates source code within the ``trunk`` directory:
::
  src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7 for the longwave
  src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7 for the shortwave

or here if you have set ``GFDL_SOC`` as an environment variable:
::
  $GFDL_SOC/data/spectra/ga7/sp_lw_ga7 for the longwave
  $GFDL_SOC/data/spectra/ga7/sp_sw_ga7 for the shortwave

+--------------------------------+----------+----------------------------------------------------------------------------+
| Name                           | Default  | Description                                                                |
+================================+==========+============================================================================+
| ``socrates_hires_mode``        | False    | If ``False`` then run in 'GCM mode', and                                   |
|                                |          | if ``True`` then use high-res spectral file                                |
+--------------------------------+----------+----------------------------------------------------------------------------+
| ``lw_spectral_filename``       | 'unset'  | Longwave spectral file, which can be found at Socrates source code         |
|                                |          | ``data`` directory (e.g. ``data/spectra/ga7/sp_lw_ga7``) or self-generated |
+--------------------------------+----------+----------------------------------------------------------------------------+
| ``sw_spectral_filename``       | 'unset'  | Shortwave spectral file, which can be found at Socrates source code        |
|                                |          | ``data`` directory (e.g. ``data/spectra/ga7/sp_sw_ga7``) or self-generated |
+--------------------------------+----------+----------------------------------------------------------------------------+
| ``lw_hires_spectral_filename`` | 'unset'  | High-res longwave spectral file if ``socrates_hires_mode`` is ``True``     |
+--------------------------------+----------+----------------------------------------------------------------------------+
| ``sw_hires_spectral_filename`` | 'unset'  | High-res shortwave spectral file if ``socrates_hires_mode`` is ``True``    |
+--------------------------------+----------+----------------------------------------------------------------------------+


CO2, ozone and other gases
^^^^^^^^^^^^^^^^^^^^^^^^^^

To include radiative effects of water vapor, CO2 and ozone, the following switches should be ``True``:

+-------------+---------------+-----------------------------------------------------+
| Name        | Default       | Description                                         |
+=============+===============+=====================================================+
| ``inc_h2o`` | True          | To include radiative effects of water vapor         |
+-------------+---------------+-----------------------------------------------------+
| ``inc_co2`` | True          | To include radiative effects of CO2                 |
+-------------+---------------+-----------------------------------------------------+
| ``inc_o3``  | True          | To include radiative effects of ozone               |
+-------------+---------------+-----------------------------------------------------+

In addition, you need to set their concentrations by specifing them directly or reading from the external files (close attention should be paid to the concentration units):

+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| Name                              | Default       | Description                                                                 |
+===================================+===============+=============================================================================+
| ``account_for_effect_of_water``   | True          | - ``False``: radiation is fed water mixing ratios = 0                       |
|                                   |               | - ``True``: radiation is fed mixing ratios based on model specific humidity |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``account_for_effect_of_ozone``   | True          | - ``False``: radiation is fed ozone mixing ratios = 0                       |
|                                   |               | - ``True``: radiation is fed mixing ratios based on model ozone field       |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``do_read_co2``                   | False         | - Read CO2 from an external file?                                           |
|                                   |               | - If ``True``, needs to specify CO2 file and variable names                 |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``co2_file_name``                 | 'co2'         | Name of file containing CO2 field - n.b. don't need to include ``'.nc'``    |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``co2_field_name``                | 'co2'         | Name of CO2 variable in CO2 file (specified by ``co2_file_name``)           |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``input_co2_mmr``                 | False         | - ``True`` if the input file contain values as ``mass mixing ratio``        |
|                                   |               | - ``False`` if the input file contain values as ``volume mixing ratio``     |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``co2_ppmv``                      | 300           | Default CO2 concentration in ``ppmv``                                       |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``do_read_ozone``                 | False         | - Read ozone from an external file?                                         |
|                                   |               | - If ``True``, needs to specify ozone file and variable names               |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``ozone_file_name``               | 'ozone'       | Name of file containing ozone field - n.b. don't need to include ``'.nc'``  |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``ozone_field_name``              | 'ozone'       | Name of ozone variable in ozone file (specified by ``ozone_file_name``)     |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+
| ``input_o3_file_is_mmr``          | True          | - ``True`` if the input file contain values as ``mass mixing ratio``        |
|                                   |               | - ``False`` if the input file contain values as ``volume mixing ratio``     |
+-----------------------------------+---------------+-----------------------------------------------------------------------------+

To include the radiative effects of other gases, such as CO, CH4, O2, SO2, CFC, etc, first you need to turn on the switches starting with ``inc_`` (default ``False``), then specify the corresponding concentrations through variables ending with ``_mix_ratio`` in the namelist.


Diagnostics
-----------

Diagnostics from Socrates are under module name ``socrates``. The outputs include the temperature tendencies due to LW/SW radiation, LW/SW radiation fluxes at each level, and the fluxes at surface and the top of the atmosphere (TOA).

+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
| Name                     | Description                                         | Units               | Dimension (not including time) |
+==========================+=====================================================+=====================+================================+
|``soc_tdt_lw``            | Socrates temperature tendency due to LW radiation   | Ks :math:`^{-1}`    | (pfull, lat, lon)              |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_tdt_sw``            | Socrates temperature tendency due to SW radiation   | Ks :math:`^{-1}`    | (pfull, lat, lon)              |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_tdt_rad``           | Socrates temperature tendency due to radiation      | Ks :math:`^{-1}`    | (pfull, lat, lon)              |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_flux_lw``           | Socrates net LW flux (positive up)                  | Wm :math:`^{-2}`    | (phalf, lat, lon)              |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_flux_sw``           | Socrates net SW flux (positive up)                  | Wm :math:`^{-2}`    | (phalf, lat, lon)              |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_surf_flux_lw``      | Socrates net LW surface flux (positive up)          | Wm :math:`^{-2}`    | (lat, lon)                     |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_surf_flux_lw_down`` | Socrates LW surface flux down                       | Wm :math:`^{-2}`    | (lat, lon)                     |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_surf_flux_sw``      | Socrates net SW surface flux (positive down)        | Wm :math:`^{-2}`    | (lat, lon)                     |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_surf_flux_sw_down`` | Socrates SW surface flux down                       | Wm :math:`^{-2}`    | (lat, lon)                     |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_olr``               | Socrates TOA LW flux (positive up)                  | Wm :math:`^{-2}`    | (lat, lon)                     |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_toa_sw``            | Socrates net TOA SW flux (positive down)            | Wm :math:`^{-2}`    | (lat, lon)                     |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_toa_sw_down``       | Socrates net TOA SW flux down                       | Wm :math:`^{-2}`    | (lat, lon)                     |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_coszen``            | Socrates cosine (zenith_angle)                      | None                | (lat, lon)                     |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_co2``               | Socrates CO2 concentration (mass mixing ratio)      | kg kg :math:`^{-1}` | (pfull, lat, lon)              |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_ozone``             | Socrates ozone concentration (mass mixing ratio)    | kg kg :math:`^{-1}` | (pfull, lat, lon)              |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+
|``soc_spectral_olr``      | Socrates substellar OLR spectrum                    | Wm :math:`^{-2}`    | (socrates_lw_bins, lat, lon)   |
+--------------------------+-----------------------------------------------------+---------------------+--------------------------------+


Relevant modules and subroutines
--------------------------------

The Socrates radiation scheme is initiatized and called by ``src/atmos_spectral/driver/solo/idealized_moist_phys.F90``.

The major modules/files under ``src/atmos_param/socrates/interface/`` are:

* ``socrates_interface.F90`` and ``socrates_calc.F90``: The Socrates interfaces that initialize/finalize the Socrates, call subroutines to get inputs, set options for radiation and run core radiaiton code, and output the diagnostics.
* ``socrates_config_mod.f90``: module to set the namelist, including the solar radiation options, time-step, and concentrations of CO2, ozone and other well-mixed gases
* ``read_control.F90`` and ``set_control.F90``: The Socrates use the ``StrCtrl`` structure to control the switches for core radiaiton code. For example, if you want to include the effects of CO2, you not only need to provide the value of CO2 concentration, but also need to turn on the switch to tell Socrates to calculate its effect: set ``control%l_co2 = .true.``, where ``control`` is a ``StrCtrl`` structure. Basically, all the logical switches are set in these two files.
* ``set_bound.F90`` and ``set_dimen.F90``: modules to set the boundary fields and dimensions for the radiation code
* ``set_atm.F90``, ``set_aer.F90``, and ``set_cld.F90``: set the input atmospheric profiles, aerosol fields and clouds fields for the core radiation code (currently the aerosol and clouds are not activated)

Other radiation schemes employed in Isca can be found at:

* RRTM: see ``src/atmos_param/rrtm_radiation``
* Two-stream gray radiation: see ``src/atmos_param/two_stream_gray_rad``

References
----------
[MannersEtAl2015]_
[Thomson_and_Vallis2019]_

Authors
-------
This documentation was written by Qun Liu (thanks to Stephen's PR of Socrates), peer reviewed by Stephen Thomson, and quality controlled by Ruth Geen.
