..  DO NOT MODIFY THIS FILE UNLESS YOU ARE A CORE MAINTAINER OF ISCA!

..
    This is a reStructuredText template file for creating
    a new documentation entry for the Isca model.
    
    Please make a copy of this file with the appropriate file name and place it
    to the appropriate location within docs/source/ and start writing.
    Once you are done, remove all the comments from your .rst file.
    
    Here is a guide on reST formatting:
    https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html


RRTM Radiation Scheme Interface
===================================

Summary
-------

RRTM (Rapid Radiative Transfer Model) is a multiband correlated-k radiative transfer scheme, described in [MlawerEtAl1997]_ and [CloughEtAl2005]_.

The correlated-k method, with k being the absorption coefficient, is a means to efficiently calculate radiative transfer over a broad spectral range by collecting wave number intervals with similar spectral properties and by supposing that these spectral properties are correlated from one level to another. A relatively small set of absorption coefficients can then be chosen to be representative of the absorption coefficients for all frequencies, leading to an enormous speed-up over line-by-line calculations and much better accuracy than traditional band methods that more simplistically just group together similar wave numbers. The implementation of this scheme in Isca largely follows that of Jucker and Gerber (2017) in the MiMA model (see [VallisEtAl2018]_).

The fortran code is found in the folder ``src/atmos_param/rrtm_radiation/``
Specifically, the code is driven from ``src/atmos_param/rrtm_radiation/rrtm_radiation.f90``

Namelist options
----------------

The RRTM namelist ``rrtm_radiation_nml`` can be found at ``src/atmos_param/rrtm_radiation/rrtm_radiation.f90``.

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
|``solr_cnst``               | 1368.22  | solar constant (units: Wm :math:`^{-2}`)                                                |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``frierson_solar_rad``      | False    | Options for annual-mean incoming solar, as prescribed in Frierson's grey scheme         |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``del_sol``                 | 0.95     | Latitudinal variation of shortwave radiation, as prescribed in Frierson's grey scheme   |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``del_sw``                  | 0        | Latitudinal variation of shortwave radiation, as prescribed in Frierson's grey scheme   |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``slowdown_rad``            | 1        | Factor to simulate slower seasonal cycle: >1 means faster, <1 slower. NB. This can also |
|                            |          | be achieved by using the ``orbital_period`` option in ``constants_nml``                 |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``use_dyofyr``              | False    | Option to use day of the year to compute the Earth-Sun distance. NB. This is done       |
|                            |          | within RRTM, and assumes 365days/year!  ______________________________________   ______ |
+----------------------------+----------+-----------------------------------------------------------------------------------------+
|``solrad``                  | 1.       | Sets Earth-Sun distance if ``use_dyofyr`` is false.                                     |
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
| ``lonstep``                | 1        | Allows subsampling of fields along longitude for faster radiation calculation                           |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``do_zm_tracers``          | False    | Feed only the zonal mean of tracers to radiation                                                        |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``do_zm_rad``              | False    | Only compute zonal mean radiation                                                                       |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+


The following options allow components of the radiative fluxes to be prescribed from input files:
+----------------------------+-------------+---------------------------------------------------------------------------------------------------------+
| Name                       | Default     | Description                                                                                             |
+============================+=============+=========================================================================================================+
| ``do_read_radiation``      | False       | Read SW and LW radiation in the atmosphere from external file? Surface fluxes are still computed        |
+----------------------------+-------------+---------------------------------------------------------------------------------------------------------+
| ``radiation_file``         | 'radiation' | File name to read radiation from                                                                        |
+----------------------------+-------------+---------------------------------------------------------------------------------------------------------+
| ``rad_missing_value``      | -1.e19.     | Missing value in input files: if <0, replace everything below this value with 0                         |
|                            |             |if >0, replace everything above this value with 0                                                        |
+----------------------------+-------------+---------------------------------------------------------------------------------------------------------+
| ``do_read_sw_flux``        | False       | Read SW surface fluxes from external file?                                                              |
+----------------------------+-------------+---------------------------------------------------------------------------------------------------------+
| ``sw_flux_file``           | 'sw_flux'   | File name to read SW surface fluxes from                                                                |
+----------------------------+-------------+---------------------------------------------------------------------------------------------------------+
| ``do_read_lw_flux``        | False       | Read LW surface fluxes from external file?                                                              |
+----------------------------+-------------+---------------------------------------------------------------------------------------------------------+
| ``lw_flux_file``           | 'lw_flux'   | File name to read LW surface fluxes from                                                                |
+----------------------------+-------------+---------------------------------------------------------------------------------------------------------+


While clouds are not currently incorporated into RRTM within Isca, the following options allow the albedo to be modified based on where precipitation is falling:
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| Name                       | Default  | Description                                                                                             |
+============================+==========+=========================================================================================================+
| ``do_precip_albedo``       | False    | Modify albedo depending on large scale precipitation                                                    |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``precip_albedo_mode``     | 'full'   | Select whether to use precipitation from only large scale condensation ('lscale'),                      |
|                            |          | only convection ('conv') or both ('full')                                                               |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``precip_albedo``          | 0.35     | Set cloud albedo if do_``precip_albedo`` is True                                                        |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``precip_lat``             | 0.0      | Apply precip_albdeo poleward of this latitude                                                           |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+


Some safety boundaries:
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| Name                       | Default  | Description                                                                                             |
+============================+==========+=========================================================================================================+
| ``h2o_lower_limit``        | 2.e-7    | Set lower limit on water vapor to be seen by RRTM, input values below this are replaced by this value   | 
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``temp_lower_limit``       | 100.     | Set lower limit on temperature to be seen by RRTM, input values below this are replaced by this value   |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+			 
| ``temp_upper_limit``       | 370.     | Set upper limit on temperature to be seen by RRTM, input values above this are replaced by this value   |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+


CO2, ozone and other gases
^^^^^^^^^^^^^^^^^^^^^^^^^^

Isca does not include advection or chemistry of ozone and carbon dioxide. Values of these can be prescribed using the following options:
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| Name                       | Default  | Description                                                                                             |
+============================+==========+=========================================================================================================+
| ``do_read_ozone``          | False    | Read ozone from a NetCDF file? NB this is the only way to apply ozone in RRTM                           |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``ozone_file``             | 'ozone'  | File name for ozone input file to read (assumed to be the same as the ozone variable name in NetCDF)    |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``input_o3_file_is_mmr``   | True     | Does the ozone input file contain values as a mass mixing ratio (True) or volume mixing ratio (False)?  |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``co2ppmv``                | 300.     | CO2 concentration in ppmv                                                                               |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``do_read_co2``            | False    | Read CO2 concentraton from a NetCDF file?                                                               |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``co2_file``               | 'co2'    | File name for CO2 input file to read                                                                    |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``co2_variable_name``.     | 'co2'    | Variable name in CO2 input file                                                                         |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+


As a default, RRTM will use the atmospheric specific humidity to calculate radiative fluxes, and it needs to convert this to a volume mixing ratio. Water vapor can also be prescribed from a file or with a constant value: 
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| Name                       | Default  | Description                                                                                             |
+============================+==========+=========================================================================================================+
| ``convert_sphum_to_vmr``   | True     | Convert specific humidity to volume mixing ratio
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``do_read_h2o``            | False    | Read water vapor from an NetCDF file?                                                                   |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``h2o_file``               | 'h2o'    | Filename to read water vapor from (assumed to be the same as the h2o variable name in the NetCDF)       |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``do_fixed_water``         | False    | Feed a fixed value for water vapor to RRTM?                                                             |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``fixed_water``            | 2.e-6    | Fixed value to use if ``do_fixed_water`` is True                                                        |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``fixed_water_pres``       | 10000.   | Apply this fixed value above which pressure level?                                                      |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``fixed_water_lat``        | 90.      | Apply this fixed value equatorward of which latitude?                                                   |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+


Values of secondary gases can also be prescribed:
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| Name                       | Default  | Description                                                                                             |
+============================+==========+=========================================================================================================+
| ``include_secondary_gases``| False    | Use non-zero values for the following secondary gases                                                   |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``ch4_val``                | 0.       | CH4 (Methane)                                                                                           |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``n2o_val``                | 0.       | N2O (Nitrous Oxide)                                                                                     |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``o2_val``                 | 0.       | O2 (Oxygen)                                                                                             |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``cfc11_val``              | 0.       | CFC11 (Trichlorofluoromethane)                                                                          |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``cfc12_val``              | 0.       | CFC12 (Dichlorodifluoromethane)                                                                         |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``cfc22_val``              | 0.       | CFC22 (Chlorodifluoromethane)                                                                           |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+
| ``ccl4_val``               | 0.       | CCl4 (Carbon tetrachloride)                                                                             |
+----------------------------+----------+---------------------------------------------------------------------------------------------------------+


		  
Diagnostics
-----------
.. What diagnostics are available for this part of the code.

Diagnostics from RRTM are under module name ``rrtm_radiation``. The outputs include the temperature tendencies due to LW/SW radiation, LW/SW radiation fluxes at each level, and the fluxes at surface and the top of the atmosphere (TOA).



	 
Relevant modules and subroutines
--------------------------------
.. List the names of relevant modules, subroutines, functions, etc.
.. You can add also code snippets, using Sphinx code formatting



References
----------
..
   Add relevant references. This is done in 2 steps:
   1. Add the reference itself to docs/source/references.rst
   2. Insert the citation key here, e.g. [Vallis2017]_
   
   See the Contributing guide for more info.

[CloughEtAl2005]_
[MlawerEtAl1997]_
[VallisEtAl2018]_

Authors
-------
This documentation was written by Ruth Geen, peer reviewed by Neil Lewis, and quality controlled by Matthew Henry
