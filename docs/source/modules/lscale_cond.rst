Large Scale Condensation and Precipitation
==========================================

Summary
-------

The module computes the large scale temperature and specific humidity adjustments needed in model layers where the relative humidity exceeds a threshold relative humidity, and returns the mass of rain and snow (or other frozen precipitation) that reaches the ground. The module also outputs the temperature tendency and the specific humidity tendency. Features include the option for the re-evaporation of falling precipitation and energetically consistent adjustment with precipitation type.

See [Frierson2006a]_ (Section 2e) for a detailed description `here <https://doi.org/10.1175/JAS3753.1>`_.

Namelist options
----------------

There are three namelist variables. 

:hc: The relative humidity at which large scale condensation, where :math:`0.0 <= hc <= 1.0`. Default is :math:`hc=1.0`.
:do_evap: The flag for the re-evaporation of moisture in sub-saturated layers below, if ``True`` then re-evaporation is performed. Default is ``False``.
:do_simple: If ``True`` then all precipitation is rain/liquid precipitation, there is no snow/frozen precipitation. Default is ``False``.

Diagnostics
-----------
The diagnostics are not sent out within the call to convection but are sent out from ``idealized_moist_phys.F90`` with the module name ``atmosphere``. The diagnostics relevant to the lscale_cond module are set out below:

+-------------------+--------------------------------------+------------------------------+
| Name              | Description                          | Units                        |
+===================+======================================+==============================+
| cond_dt_qg        | Moisture tendency                    |:math:`kg~kg^{-1}~s^{-1}`     |
+-------------------+--------------------------------------+------------------------------+
| cond_dt_tg        | Temperature tendency                 |:math:`K~s^{-1}`              |
+-------------------+--------------------------------------+------------------------------+
| cond_rain         | Rain from condensation               |:math:`kg~m^{-2}~s^{-1}`      |
+-------------------+--------------------------------------+------------------------------+
| precip            | Rain and Snow from resolved and      |:math:`kg~m^{-2}~s^{-1}`      |
|                   | parameterised condensation/convection|                              |
+-------------------+--------------------------------------+------------------------------+

Relevant Modules and Subroutines
--------------------------------
Relevant modules are: sat_vapor_pres_mod, utilities_mod, constants_mod.

References
----------

The reference list within this docs page are: [Frierson2006a]_

Authors
-------
This documentation was written by Ross Castle and reviewed by Ruth Geen
