Convection scheme: Simple Betts-Miller
=================================


Summary
-------
The simplified Betts-Miller (SBM) scheme is a moist adjustment convection scheme developed by [Frierson2007]_ and the modifications described in [OGormanSchneider2008]_. The SBM scheme is Isca's default convection scheme. The key purpose of the routine is to output a moisture tendency ``deltaq``, a temperature tendency ``deltaT``, and a convective rainfall amount ``rain``. The scheme has both deep convection and non-precipitating shallow convection.


Logic of the routine
-------

1. Perform a particle ascent to calculating the environmental vertical profiles of moisture and temperature, and convective available potential energy (CAPE) within a grid cell.

2. Test if the convection scheme should trigger. The trigger is: :math:`CAPE>0`.

3. Calculate the temperature (:math:`T_{ref}`) and moisture reference profiles (:math:`q_{ref}`). Think of these profiles as the desirable parcel properties for the atmosphere to be stable (with respect to convection).

4. Calculate the first guess reference profiles using equations 1-2 from [Frierson2007]_: 

:math:`\partial q = - \frac{q-q_{ref}}{\tau}`,

:math:`\partial T = - \frac{T-T_{ref}}{\tau}`.

These profiles are called *first guess* profiles as they will later be corrected (to conserve enthalpy).

5. Determine if deep convection, shallow convection or no convection occurs in the grid cell (see Section 2 of [Frierson2007]_ for how these are separated). If the parcel profile is unstable, with respect to the reference profile, it will adjust the profiles and generate rain. This adjustment occurs for both the moisture and temperature profiles so that rain is generated when the environment has excess moisture or is cooler compared to the reference profiles).

6. If deep or shallow convection occurs, test if the vertical integral of enthalpy (:math:`H`) is conserved (it usually is not) and adjust the rainfall and tendencies to conserve enthalpy. Note that :math:`\partial H = -\partial T + \frac{L}{cp}\partial q`.

7. The final rainfall amount is given by the enthaly corrected vertical integral of the moisture tendency. The updated tendencies (no longer called *first-guess*) and rainfall are then passed back to ``idealized_moist_phys.F90`` along with other diagnostics.


Namelist options
----------------

There are five optional namelist variables. 

The first two are directly relevant to how the tendency equations and precipitation are calculated

:tau_bm: The relaxation time (:math:`\tau`) with the default value of :math:`7200` :math:`s`.
:rhbm: The critical Relative Humidity (RH) used to calculate the moisture reference profile (:math:`RH_{SBM}`) with the default value of :math:`0.8`.

The final three are used to construct the look up table for calculating the lifting condensation level (LCL).

:Tmin: The minimum look up table value with the default of :math:`173` :math:`K`.
:Tmax: The maximum look up table value with the default of :math:`335` :math:`K`.
:val_inc: The increment used to linearly interpolate within the look up table with the default of :math:`0.01`.


Diagnostics
-----------
The diagnostics are not sent out within the call to convection but are sent out from ``idealized_moist_phys.F90`` with the module name ``atmosphere``. The diagnostics relevant to the simple Betts-Miller scheme are in the table below.

+-------------------+----------------------------+------------------------------------+
| Name              | Description                | Units                              |
+===================+============================+====================================+
| dt_qg_convection  | Moisture tendency          |:math:`kgkg^{-1}s^{-1}`             |
+-------------------+----------------------------+------------------------------------+
| dt_qg_convection  | Temperature tendency       |:math:`Ks^{-1}`                     |
+-------------------+----------------------------+------------------------------------+
| convection_rain   | Convective precipitation   |:math:`kgm^{-2}s^{-1}`              |
+-------------------+----------------------------+------------------------------------+
| cape              | Convective available       |:math:`Jkg^{-1}`                    |
|                   | potential energy           |                                    |
+-------------------+----------------------------+------------------------------------+
| cin               | Convective inhibition      |:math:`Jkg^{-1}`                    |
+-------------------+----------------------------+------------------------------------+


Relevant modules and subroutines
--------------------------------

The convection scheme is located in: ``src/atmos_param/qe_moist_convection/qe_moist_convection.F90``. The files name reflects this style of a quasi-equilibirum (qe) convection scheme (for more information on convective quasi-equilibirum see [YanoPlant2016]_).

The simple Betts-Miller scheme is initialised and called by: ``/src/atmos_spectral/driver/solo/idealized_moist_phys.F90``.

Relevant routines which are called by the convection scheme are:
``src/shared/sat_vapor_pres/sat_vapor_pres.F90``.



References
----------

[Frierson2007]_
[Betts1986]_
[BettsMiller1986]_
[YanoPlant2016]_
[OGormanSchneider2008]_

Authors
----------

This documentation was written by Penelope Maher, peer reviewed by Denis Sergeev , and quality controlled by Ross Castle.
