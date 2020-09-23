Dry dynamical core (Held--Suarez) Forcing  
=======================================================================================

Documentation for ``hs_forcing_mod``. 


Summary
-------
This module is contains the 'physics options' used by Isca when the flag ``idealised_moist_model`` in ``atmosphere_nml`` is set to ``FALSE``. 

The main parametrisations included in this module are summarised in the table below, listed in terms of their subroutine name within the module.

+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``subroutine``                 | Short description                                                                  | Reference            |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``rayleigh_damping``           | | This subroutine contains options to apply a linear relaxation to the zonal and   | [Held1994]_          |
|                                | | meridional winds. The main option is constitutes a linear drag applied in the    |                      |
|                                | | boundary layer (as in Held and Suarez, 1994). It is also possible to relax the   |                      |
|                                | | entire u and v wind fields towards a specified (zonally averaged) profile.       |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``newtonian_damping``          | | This subroutine contains options to apply a linear relaxation to the model       | | [Held1994]_        |
|                                | | temperature field. The main option is a linear relaxation to the Earth-like      | | [Penn2018]_        |
|                                | | temperature profile described in Held and Suarez (1994). There are other options |                      |
|                                | | to read in the relaxation temperature profile from an input file, or to          |                      |
|                                | | configure the relaxation profile so that it produces heating rates appropriate   |                      |
|                                | | for a tidally-locked exoplanet, as in Penn and Vallis (2018).                    |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``top_down_newtonian_damping`` | | An alternative to ``newtonian_damping`` that calculates a radiative-convective   | | [VallisEtAl2018]_  |
|                                | | relaxation profile in terms of astronomical parameters. This option can be used  |                      |
|                                | | in order to obtain forcing by a linear temperature relaxation that evolves with  |                      |
|                                | | time (e.g. due to the inclusion of a seasonal cycle).                            |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``local_heating``              | | This subroutine can be used to add an additional local heating to the relaxation |                      |
|                                | | fields computed by either ``newtonian_damping`` or ``top_down_newtonian_damping``|                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+



Description of forcing 
----------------------


References
----------

| [Held1994]_ 
| [Penn2018]_
| [VallisEtAl2018]_
