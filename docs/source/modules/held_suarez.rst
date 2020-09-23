Dry dynamical core (Held--Suarez) Forcing  
=======================================================================================

Documentation for ``hs_forcing_mod``. 


Summary
-------
This module is contains the 'physics options' used by Isca when the flag ``idealised_moist_model`` in ``atmosphere_nml`` is set to ``FALSE``. 

**In this scenario, Isca is run as a dry model.** That is, no condensible species is included in the model, and consequently the effect of latent heating on the atmospheric dynamics is omitted. 

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

In the next few sections we will describe each of the main options above in detail. A list of namelist parameters for ``hs_forcing_mod``, and their effect and default values, is included at the end of this documentation page. 


Rayleigh damping 
----------------------

The subroutine ``rayleigh_damping`` contains options to apply a linear relaxation to the zonal and meridional winds. 

| **Held and Suarez linear drag**
| The default option for the ``rayleigh_damping`` subroutine (selected with ``relax_to_specified_wind=.FALSE.``) applies a linear drag to the zonal and meridional wind within the model boundary layer. The drag takes the form 

.. math::
   \frac{\partial\mathbf{u}}{\partial t} = -k_{\mathbf{u}}(\sigma)\mathbf{u}

following [Held1994]_. :math:`\mathbf{u}=(u,v)` is the horizontal wind, and :math:`k_{\mathbf{u}}` is a damping coefficient that depends on pressure :math:`p`, via the model sigma coordinate :math:`\sigma=p/p_{\text{s}}`. The damping coefficient is defined as 

.. math:: 
   k_{\mathbf{u}}=k_{f}\max\left(0,\frac{\sigma-\sigma_{b}}{1-\sigma_{\text{b}}}\right) 

where :math:`k_{f}` is the damping rate, and :math:`\sigma_{b}` specifies the boundary layer depth (for :math:`\sigma<\sigma_{b}`, no relaxation is applied). By default, :math:`k_{f}=1\,\text{days}^{-1}` and :math:`\sigma_{b}=0.7`. Both of these options can be changed in the namelist (see table of namelist options at bottom of page).

| **Relaxation to a specified wind**
| If ``relax_to_specified_wind=.TRUE.``, then the ``rayleigh_damping`` subroutine instead relaxes the wind to an equilibrium profile supplied as a netcdf file (with dimensions: pressure, latitude, longitude). In this scenario, the drag takes the form 

.. math::
   \frac{\partial\mathbf{u}}{\partial t} = -k_{f}(\overline{\mathbf{u}}-\overline{\mathbf{u}}_{r})

where :math:`\mathbf{u}_{r}` is the specified relaxation wind profile, and the overline denotes a zonal average. :math:`k_{f}` is the same as above. Note that in this instance, the relaxation is applied at all model levels, not just in the boundary layer. 


Newtonian damping 
----------------------

The subroutine ``newtonian_damping`` contains options to apply a linear relexation to the model temperature field, towards a specified 'equilibrium' temperature field :math:`T_{\text{eq}}` 

.. math::
   \frac{\partial T}{\partial t} = -k_{T}(\theta,\sigma)[T-T_{\text{eq}}(\lambda,\theta,p)]

where :math:`T` is the model temperature field, and :math:`k_{T}` is the thermal damping coefficient, defined as 

.. math::
   k_{T}=k_{a} + (k_{s}-k_{a})\max\left(0,\frac{\sigma-\sigma_{b}}{1-\sigma_{b}}\right)\cos^{4}\theta

following [Held1994]_. :math:`k_{a}` is the damping coefficient in the free atmosphere, and :math:`k_{s}` is the damping coefficient at the surface. :math:`\sigma_{b}` is the boundary layer depth (as above in the Rayleigh damping section). By default, :math:`k_{a}=1/40\,\text{days}^{-1}` and :math:`k_{s}=1/4\,\text{days}^{-1}`. Both coefficients can be changed in the namelist (see namelist options table at bottom of page). 


There are four main options for :math:`T_{\text{eq}}`, which are set with the namelist parameter ``equilibrium_t_option``, input as a string. Each ``equilibrium_t_option`` is summarised in the table below. 

+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``equilibrium_t_option``       | Short description                                                                  | Reference            |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'Held_Suarez'``              | | This option relaxes the model temperature towards the Earth-like                 | [Held1994]_          |
|                                | | radiative-convective equilibrium temperature profile of Held and Suarez (1994).  |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'from_file'``                | | This option relaxes the model temperature towards a temperature profile          |                      |
|                                | | specified via an input file.                                                     |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'EXOPLANET'``                | | ...                                                                              | [Penn2018]_          |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'EXOPLANET2'``               | | ...                                                                              |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+



References
----------

| [Held1994]_ 
| [Penn2018]_
| [VallisEtAl2018]_

Authors
----------
This documentation was written by Neil Lewis, peer reviewed by X, and quality controlled by Y. 