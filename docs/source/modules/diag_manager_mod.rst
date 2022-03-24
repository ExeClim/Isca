Diagnostics Manager Module & Saving model output 
==============

Summary
-------

This module handles the writing of diagnostic output to netCDF files. The user can specify which fields should be output and at which temporal resolution (e.g. monthly means, daily means ... ). The source code is located at ``src/shared/diag_manager/diag_manager.F90``. 


Namelist options
----------------

+--------------------------------+----------+-----------------------------------------------------------------------------------------+
| Name                           | Default  | Description                                                                             |
+================================+==========+=========================================================================================+
|``append_pelist_name``          | False    | Decides whether to append the pelist_name to file name                                  |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``mix_snapshot_average_fields`` | False    | Allow both time average and instantaneous fields in the same output file                |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``max_files``                   | 31       | Sets the maximum number of output files allowed                                         |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``max_output_fields``           | 300      | Sets the maximum number of output fields allowed                                        |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``max_input_fields``            | 300      | Sets the maximum number of input fields allowed                                         |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``max_axes``                    | 60       | Sets the maximum number of independent axes                                             |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``do_diag_field_log``           | False    | to log the registration of the data field, put the OPTIONAL parameter                   |
|                                |          | ``do_not_log`` = ``False`` and the namelist variable ``do_diag_field_log`` to ``True``  |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``write_bytes_in_files``        | False    | Write out the number of bytes of data saved to this file                                |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``debug_diag_manager``          | False    | When set to true, bounds are checked in ``send_data`` (send data to output fields)      |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``max_num_axis_sets``           | 25       | Set maximum number of axes for output (e.g. time and space axes)                        |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``use_cmor``                    | False    | Let the ``diag_manager`` know if the missing value (if supplied) should be overridden   |
|                                |          | to be the CMOR standard value of -1.0e20                                                |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``issue_oor_warnings``          | True     | If ``True`` check for values outside the valid range. This range is passed to the       |
|                                |          | ``diag_manager_mod`` via the OPTIONAL variable range in the                             |
|                                |          | ``register_diag_field`` function                                                        |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+
|``oor_warnings_fatal``          | True     | If ``True`` issue a fatal error if any values for the output field are outside the      |
|                                |          | given range                                                                             |
+--------------------------------+----------+-----------------------------------------------------------------------------------------+


Diagnostics
-----------

This part of the code does not have its own diagnostics, but rather handles the saving of all variables. See also ``/src/extra/python/isca/diagtable.py``

Output files 
^^^^^^^^^^^^

In order to save output in Isca, an output file is created first in your experiment runscript (see examples in ``/exp/test_cases/``). Commonly used output timesteps include monthly, daily or x-hourly.

``diag.add_file('atmos_monthly', 30, 'days', time_units='days')``
``diag.add_file('atmos_daily', 1, 'days', time_units='days')``
``diag.add_file('atmos_6_hourly', 6, 'hours', time_units='hours')``

The frequency at which data is saved can be set in the ``main_nml``: ::

	'main_nml': {
		'dt_atmos': 600,    
		'days': 30,    
		'calendar': 'thirty_day'    
	}

For example, set ``'days': 15`` and ``'calendar': 'fifteen_day'``. 


Output fields
^^^^^^^^^^^^^

An output field is created via ``diag.add_field(module, name, time_avg, files)`` in the experiment runscript. 
The default for ``time_avg`` = False, the default for ``files`` = None. 
``time_avg`` is usually set to True for most variables when an output field is defined.

If ``files`` = None, then the diagnostics will be saved to all of the given output files (in our example monthly, daily and 6h). 
An output file can be specified via e.g. ``files=['atmos_6_hourly']`` in 
``diag.add_field('dynamics', 'ucomp', time_avg=True, files=['atmos_6_hourly'])`` if 6h zonal winds shall be saved, but not monthly/daily



Below is a list of commonly saved diagnostics. See the relevant modules for an exhaustive list of available diagnostics. 

+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| Module                   | Name                 | Dimensions              | Description                                 | Units              |
+==========================+======================+=========================+=============================================+====================+
| ``dynamics``             | ``ps``               | (time, lat, lon)        | surface pressure                            | :math:`Pa`         |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``bk``               | (phalf)                 | vertical coordinate sigma values            | :math:`Pa`         |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``pk``               | (phalf)                 | vertical coordinate pressure values         | :math:`Pa`         |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``slp``              | (time, lat, lon)        | sea level pressure                          | :math:`Pa`         |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``height``           | (time, pfull, lat, lon) | geopotential height at full model levels    | :math:`m`          |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``zsurf``            | (lat, lon)              | geopotential height at the surface          | :math:`m`          |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``u_comp``           | (time, pfull, lat, lon) | zonal component of the horizontal winds     | :math:`m/s`        |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``v_comp``           | (time, pfull, lat, lon) | meridional component of the horizontal winds| :math:`m/s`        |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``omega``            | (time, pfull, lat, lon) | vertical velocity                           | :math:`Pa/s`       |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``sphum``            | (time, pfull, lat, lon) | specific humidity                           | :math:`kg/kg`      |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``temp``             | (time, pfull, lat, lon) | temperature                                 | :math:`K`          |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``vor``              | (time, pfull, lat, lon) | vorticity                                   | :math:`1/s`        |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``dynamics``             | ``div``              | (time, pfull, lat, lon) | divergence                                  | :math:`1/s`        |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``atmosphere``           | ``precipitation``    | (time, lat, lon)        | precipitation from resolved, parameterised  | :math:`kg/(m^2 s)` |
|                          |                      |                         | and snow                                    |                    |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``atmosphere``           | ``rh``               | (time, pfull, lat, lon) | relative humidity                           | :math:`\%`         |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``mixed-layer``          | ``t_surf``           | (time, lat, lon)        | surface temperature                         | :math:`K`          |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``mixed-layer``          | ``flux_t``           | (time, lat, lon)        | sensible heat flux at the surface (up)      | :math:`W/m^2`      |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``mixed-layer``          | ``flux_lhe``         | (time, lat, lon)        | latent heat flux at the surface (up)        | :math:`W/m^2`      |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``rrtm_radiation``       | ``flux_sw``          | (time, lat, lon)        | net shortwave flux at the surface (down)    | :math:`W/m^2`      |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+
| ``rrtm_radiation``       | ``flux_lw``          | (time, lat, lon)        | longwave flux at the surface (down only)    | :math:`W/m^2`      |
+--------------------------+----------------------+-------------------------+---------------------------------------------+--------------------+


Relevant modules and subroutines
--------------------------------

The ``diag_manager_mod`` uses several modules and subroutines, including 

* ``diag_axis``
* ``diag_grid``
* ``diag_output``
* ``diag_util``
* ``diag_data``
* ``diag_table``


.. References
.. ----------
.. ..
..    Add relevant references. This is done in 2 steps:
..    1. Add the reference itself to docs/source/references.rst
..    2. Insert the citation key here, e.g. [Vallis2017]_
   
..    See the Contributing guide for more info.

.. None

Authors
-------

This documentation was written by Marianne Pietschnig, peer reviewed by Stephen Thomson and quality controlled by Ross Castle. 
