damping_driver.f90
======================

Summary
-------
The ``damping driver`` module is called by the ``idealized_moist_phys`` module by setting ``do_damping`` to true. It controls the upper level momentum damping in Isca. It controls 4 functions:

1. **Rayleigh friction** which acts at levels ``1 to kbot``. This function is located in ``damping_driver`` itself.
2. A (orographic) **mountain gravity wave drag** module may be called. This is the ``mg_drag`` module.
3. A (non orographic) **convective gravity wave drag** module may be called. This is the ``cg_drag`` module.
4. A **time independent drag** may be called. This function is located in ``damping_driver`` itself.

There is another module references called ``topo_drag`` however this is not available in Isca at present. 

Namelist options
----------------
``trayfric`` - (for Rayleigh friction) damping time in seconds for rayleigh damping momentum in the top model layers, the number of which is specified by ``nlev_rayfric`` (non namelist parameter, automatically determined in code). If ``trayfric`` < 0 then time in days. Default 0.

``do_rayleigh`` - On/Off switch for doing Rayleigh friction. Default False

``sponge_pbottom`` - (for Rayleigh friction) used to calculate ``nlev_rayfric``, it specifies the bottom level where the Rayleigh friction starts. Default 50Pa.

``do_cg_drag`` - On/Off switch for doing mountain gravity wave drag. Default False

``do_topo_drag`` - On/Off switch for doing the topo_drag module which is currently unavailable. Default False. 

``do_mg_drag`` - On/Off switch for doing convective gravity wave drag. Default False.

``do_conserve_energy`` - (for Rayleigh friction) On/Off switch for also calculating the temperature tendency (if True). Default False.

``do_const_drag`` - On/Off switch for doing the constant drag. Default False.

``const_drag_amp`` - (for constant drag) Parameter for adjusting drag (Amplitude?). Default 3.e-04.

``const_drag_off`` - (for constant drag) Parameter for adjusting drag (Offset?). Default 0.

For a typical idealized Earth set up there is no parameterised gravity wave drag and rayleigh damping is only needed to keep the model stable. The namelist would then look like: 
``'do_rayleigh': True,
'trayfric': -0.25,
'sponge_pbottom':  5000.,
'do_conserve_energy': True``

Diagnostics
-----------
There are many diagnostics either in or passed by damping driver. 

**Rayleigh friction:**

+-----------------------+------------------------------------+------------------------+
| Name                  | Description                        | Units                  |
|                       |                                    |                        |
+=======================+====================================+========================+
| ``udt_rdamp``         | u wind tendency                    |:math:`m s^{-2}`        |
+-----------------------+------------------------------------+------------------------+
| ``vdt_rdamp``         | v wind tendency                    |:math:`m s^{-2}`        |
+-----------------------+------------------------------------+------------------------+
| ``tdt_diss_rdamp``    | dissipative heating                |:math:`K s^{-1}`        |
+-----------------------+------------------------------------+------------------------+
| ``diss_heat_rdamp``   | integrated dissipative heating     |:math:`W m^{-2}`        |
+-----------------------+------------------------------------+------------------------+

**Mountain GWD:**

+-----------------------+------------------------------------+------------------------+
| Name                  | Description                        | Units                  |
|                       |                                    |                        |
+=======================+====================================+========================+
| ``udt_gwd``           | u wind tendency                    |:math:`m s^{-2}`        |
+-----------------------+------------------------------------+------------------------+
| ``vdt_gwd``           | v wind tendency                    |:math:`m s^{-2}`        |
+-----------------------+------------------------------------+------------------------+
| ``taubx``             | x base flux                        |:math:`kg m^{-1} s^{-2}`|
+-----------------------+------------------------------------+------------------------+
| ``tauby``             | y base flux                        |:math:`kg m^{-1} s^{-2}`|
+-----------------------+------------------------------------+------------------------+
| ``taus``              | saturation flux                    |:math:`kg m^{-1} s^{-2}`|
+-----------------------+------------------------------------+------------------------+
| ``tdt_diss_gwd``      | dissipative heating                |:math:`K s^{-1}`        |
+-----------------------+------------------------------------+------------------------+
| ``diss_heat_gwd``     | integrated dissipative heating     |:math:`W s^{-2}`        |
+-----------------------+------------------------------------+------------------------+
| ``sgsmtn``            | sub-grid scale topography variance |:math:`m`               |
+-----------------------+------------------------------------+------------------------+

**Convective GWD:**

+-----------------------+------------------------------------+------------------------+
| Name                  | Description                        | Units                  |
|                       |                                    |                        |
+=======================+====================================+========================+
| ``udt_cgwd``          | u wind tendency                    |:math:`m s^{-2}`        |
+-----------------------+------------------------------------+------------------------+

**Constant Drag:**

+-----------------------+------------------------------------+------------------------+
| Name                  | Description                        | Units                  |
|                       |                                    |                        |
+=======================+====================================+========================+
| ``udt_cnstd``         | u wind tendency                    |:math:`m s^{-2}`        |
+-----------------------+------------------------------------+------------------------+

**topo_drag:**

Note: These are not currently available

+-----------------------+------------------------------------+------------------------+
| Name                  | Description                        | Units                  |
|                       |                                    |                        |
+=======================+====================================+========================+
| ``udt_topo``          | u wind tendency                    |:math:`m s^{-2}`        |
+-----------------------+------------------------------------+------------------------+
| ``vdt_topo``          | v wind tendency                    |:math:`m s^{-2}`        |
+-----------------------+------------------------------------+------------------------+

Relevant modules and subroutines
--------------------------------
The code is split into 4 subroutines; ``damping_driver``, ``damping_driver_init``, ``damping_driver_end`` and ``rayleigh``. The ``_init`` and ``_end`` subroutines are for initializing and closing the module. The majority of the ``damping driver`` code is just a switchboard with the exeption of the Rayleigh and constant drag calculations. The calculations for the other drag schemes are given in their own documentation.

**Rayleigh Drag**

Located in the ``rayleigh`` subroutine. This code damps the momentum tomward zero in the upper model levels specified. The zonal/meridional tendency for each grid is calculated my multiplying the zonal/meridional velocity by a factor determined by the pressure and rayleigh parameters. The higher the wind velocity, the more the damping.

The temperature tendency is calculated using the wind velocities, wind tendencies and the heat capacity of air.

**Constant Drag**

Located in the ``damping_driver`` subroutine. This is modeled on Alexander-Dunkerton winter average, it uses a 3rd order polynomial and the constant drag parameters to caclulate a time invariant drag. Earth use only recommended. 


References
----------
Convective Gravity Wave Drag (cg_drag) [Pierrehumbert1986]_

Orographic Gravity Wave Drag (mg_drag) [Alexander1999]_
   
Authors
-------
This documentation was written by Ross Castle, peer reviewed by Stephen Thompson, and quality controlled by Matthew Henry.
