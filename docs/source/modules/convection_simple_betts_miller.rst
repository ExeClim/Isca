Convection scheme: Simple Betts-Miller
=================================


Summary
-------
The simplified Betts-Miller (SBM) scheme is a moist adjustment convection scheme developed by [Frierson2007]_ and the motifications described in [OGormanSchneider2008]_. The SBM scheme is Isca's default convection scheme. The key purpose of the routine is output a moisture tendency ``deltaq`` and a temperature tendency ``deltaT``. 


Logic of the routine
-------

1. Perform a particle ascent to calculating the environmental vertical profile (moisture and temperature) and convective available potential energy (CAPE) within a grid cell.

2. Test if the convection scheme should trigger. The trigger is: :math:`CAPE>0`.

3. Calculate the temperature (:math:`T_{ref}`) and moisture reference profiles (:math:`q_{ref}`). Think of these profiles as the desirable parcel properties for the atmosphere to be stable (with respect to convection).

4. Calculate the first guess reference profiles using equations 1-2 from [Frierson2007]_: 

:math:`\partial q = - \frac{q-q_{ref}}{\tau}`,

:math:`\partial T = - \frac{T-T_{ref}}{\tau}`.

These profiles are called *first guess* profiles as they will later be corrected (to conserve enthalpy).

5. Determine if deep convection, shallow convection or no convection occurs in the grid cell (see Section 2 of [Frierson2007]_ for how these are seperated). If the parcel profile is more unstable than the reference profile then it will generate rain (this will occur for both the moisture and temperature profiles to account that condensation can occur due to an increase in moisture or a decrease in temperature).

6. If deep or shallow convection occurs, test if the vertical integral of enthalpy (:math:`H`) is conserved (it usually is not) and adjust the rainfall and tendency equations to conserve enthalpy. Note that :math:`\partial H = -\partial T + \frac{L}{cp}\partial q`.


Namelist options
----------------

There are five optional namelist variables. 

The first two are directly relevant to how the tendency equations and precipitation are calculated

:tau_bm: The relaxation time (:math:`\tau`) with the default value of :math:`7200` :math:`s`.
:rhbm: The critical Relative Humidity (RH) used to calculate the moisture reference profile (:math:`RH_{SBM}`) with the default value of :math:`0.8`.

The final three are used to construct the look up table for calculating the lifting condensation level (LCL).

:Tmin: The minimum look up table value with the default of :math:`173` :math:`K`.
:Tmax: The maximum look up table value with the default of :math:`335` :math:`K`.
:val_inc: The increment used to linearly interpolate within the look up table with the defaul of :math:`0.01`.


Diagnostics
-----------
The diagnostics are not sent out within the call to convection but are sent out from ``idealized_moist_phys.F90`` with the module name ``atmosphere``. The diagnostics relevant to the simple Betts-Miller scheme are in the table below.

+-------------------+----------------------------+------------------------------------+
| Name              | Description                | Units                              |
+===================+============================+====================================+
| dt_qg_convection  | Moisture tendency          |:math:`\frac{kg}{kg}\frac{1}{s}`    |
+-------------------+----------------------------+------------------------------------+
| dt_qg_convection  | Temperature tendency       |:math:`\frac{K}{s}`                 |
+-------------------+----------------------------+------------------------------------+
| convection_rain   | Convective precipitation   |:math:`\frac{kg}{m^2s}`             |
+-------------------+----------------------------+------------------------------------+
| cape              | Convective available       |:math:`\frac{J}{kg}`                |
|                   | potential energy           |                                    |
+-------------------+----------------------------+------------------------------------+
| cin               | Convective inhibition      |:math:`\frac{J}{kg}`                |
+-------------------+----------------------------+------------------------------------+


Relevant modules and subroutines
--------------------------------

The convection scheme is located in: ``src/atmos_param/qe_moist_convection/qe_moist_convection.F90``. The files name reflects this style of a quasi-equilibirum (qe) convection scheme (for more information on convective quasi-equilibirum see [YanoPlant2016]_).

The simple Betts-Miller scheme is initiatised and called by ``/src/atmos_spectral/driver/solo/idealized_moist_phys.F90``.

Relevant routines which are called by the convection scheme are:
```src/shared/sat_vapor_pres/sat_vapor_pres.F90```



References
----------

[Frierson2007]_
[Betts1986]_
[BettsMiller1986]_
[YanoPlant2016]_
[OGormanSchneider2008]_

Authors
----------

This documentation was written by Penelope Maher, peer reviewed by X, and quality controlled by X.
