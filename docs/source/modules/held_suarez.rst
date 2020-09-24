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

In the next few sections we will describe each of the main options above in detail. A list of namelist parameters for ``hs_forcing_mod``, and their effect and default values, is included at the end of this documentation page in **Namelist options**. 


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

where :math:`k_{f}` is the damping rate, and :math:`\sigma_{b}` specifies the boundary layer depth (for :math:`\sigma<\sigma_{b}`, no relaxation is applied). By default, :math:`k_{f}=1\,\text{days}^{-1}` and :math:`\sigma_{b}=0.7`. Both of these options can be changed in the namelist.

| **Relaxation to a specified wind**
| If ``relax_to_specified_wind=.TRUE.``, then the ``rayleigh_damping`` subroutine instead relaxes the wind to an equilibrium profile supplied as a netcdf file (with dimensions: pressure, latitude, longitude). In this scenario, the drag takes the form 

.. math::
   \frac{\partial\mathbf{u}}{\partial t} = -k_{f}(\overline{\mathbf{u}}-\overline{\mathbf{u}}_{r})

where :math:`\mathbf{u}_{r}` is the specified relaxation wind profile, and the overline denotes a zonal average. :math:`k_{f}` is the same as above. Note that in this instance, the relaxation is applied at all model levels, not just in the boundary layer. 

The name of the netcdf files are specified with the namelist options ``u_wind_file`` and ``v_wind_file``, and should exclude the ``.nc`` suffix. Note: in the netcdf file, the :math:`u` and :math:`v` variables must be the same name as their respective file names, i.e., ``u_wind_file`` and ``v_wind_file``. The path to the file should be configured in the model's python run script (see the test cases for examples of how to do this). 


Newtonian damping 
----------------------

The subroutine ``newtonian_damping`` contains options to apply a linear relexation to the model temperature field, towards a specified 'equilibrium' temperature field :math:`T_{\text{eq}}` 

.. math::
   \frac{\partial T}{\partial t} = -k_{T}(\theta,\sigma)[T-T_{\text{eq}}(\lambda,\theta,p)]

where :math:`T` is the model temperature field, and :math:`k_{T}` is the thermal damping coefficient, defined as 

.. math::
   k_{T}=k_{a} + (k_{s}-k_{a})\max\left(0,\frac{\sigma-\sigma_{b}}{1-\sigma_{b}}\right)\cos^{4}\theta

following [Held1994]_. :math:`k_{a}` is the damping coefficient in the free atmosphere, and :math:`k_{s}` is the damping coefficient at the surface. :math:`\sigma_{b}` is the boundary layer depth (as above in the Rayleigh damping section). By default, :math:`k_{a}=1/40\,\text{days}^{-1}` and :math:`k_{s}=1/4\,\text{days}^{-1}`. Both coefficients can be changed in the namelist. 


There are four main options for :math:`T_{\text{eq}}`, which are set with the namelist parameter ``equilibrium_t_option``, input as a string. Each ``equilibrium_t_option`` is summarised in the table below. 

+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``equilibrium_t_option``       | Short description                                                                  | Reference            |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'Held_Suarez'``              | | Earth-like radiative-convective equilibrium temperature profile of               | [Held1994]_          |
|                                | |  Held and Suarez (1994).                                                         |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'EXOPLANET'``                | | Equilibrium temperature structure                                                | [Penn2018]_          |
|                                | | suitable for a tidally locked terrestrial planet.                                |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'EXOPLANET2'``               | | Idealised equilibrium profile with an isothermal stratosphere, and a neutrally   |                      |
|                                | | buoyant troposphere.                                                             |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'from_file'``                | | Equilibrium temperature specified by input file.                                 |                      |
+--------------------------------+-------------------------------+----------------------------------------------------+----------------------+

| **Held and Suarez**
| The ``'Held_Suarez'`` equilibrium temperature profile is designed to mimic the temperature sturcutre of the Earth in radiative-conevective equilibrium. 

It takes the following form 

.. math:: 
   T_{\text{eq}} = \max\left\{T_{\text{str}}, \left[T^{\ast} - (\Delta\Theta)_{z}\log\left(\frac{p}{p_{0}}\right)\cos^{2}\theta\right]\left(\frac{p}{p_{0}}\right)^{\kappa}\right\}

where 

.. math:: 
   T_{\text{str}} = T_{\text{strat}} - \epsilon\sin\theta

and 

.. math::
   T^{\ast} = T_0 - (\Delta T)_{y}\sin^{2}\theta-\epsilon\sin\theta

:math:`T_{\text{strat}}` is the isothermal stratosphere temperature, which can be modified to include a latitude dependence via non-zero :math:`\epsilon`. :math:`T^{\ast}` is the surface temperature, set by the surface temperature at the equator :math:`T_{0}`, and an equator-to-pole surface temperature difference :math:`(\Delta T)_{y}` (subject to further modification due to non-zero :math:`\epsilon`). :math:`(\Delta\Theta)_{z}` determines the vertical gradient of potential temperature. :math:`p_{0}` is a reference pressure, and :math:`\kappa=R/c_p`. 

When :math:`\epsilon=0`, the equlibrium temperature profile reduces to that of [Held1994]_. 

By default, :math:`T_{\text{strat}}=200\,\text{K}`, :math:`T_{0}=315\,\text{K}`, :math:`(\Delta T)_{y}=60\,\text{K}`, and :math:`(\Delta \Theta)_{z}=10\,\text{K}`. :math:`p_{0}=1\,\text{bar}`, :math:`\kappa=2/7`, and :math:`\epsilon=0`. Each of these parameters can be set in the namelist. 



| **Tidally locked exoplanet**
| The ``'EXOPLANET'`` option configures the Newtonian damping to produce heating rates characteristic of a tidally locked exoplanet. In this configuration :math:`T_{\text{eq}}` is written 

.. math:: 
   T_{\text{eq}}=\max\left\{T_{\text{str}},\left[T^{\ast}-(\Delta\Theta)_{z}\log\left(\frac{p}{p_{0}}\right)\cos\zeta\right]\left(\frac{p}{p_{0}}\right)^{\kappa}\right\}

where :math:`T_{\text{str}}` is the same as for the ``'Held_Suarez'`` configuration, and 

.. math:: 
   T^{\ast} = T_{0} - (\Delta T)_{y})(1-\cos\zeta)-\epsilon\sin\theta

:math:`\zeta` is the solar zenith angle, calculated by ``astronomy_mod``. To achieve tidal locking, ``calendar='no_calendar'`` should be set in ``main_nml``, and, in ``constants_nml``, ``orbital_period`` and ``omega`` (planetary rotation rate) should be set so that ``orbital_period=1 / omega``. It is also possible to use this configuration to simulate a non-tidally locked planet (i.e., a planet with a diurnal cycle). This is achieved simply by setting ``omega`` and ``orbital_period`` in such a way that the planet is not tidally locked. For example output from experiments run using the ``'EXOPLANET'`` configuration, see [Penn2018]_. 


**Note:** to achieve realistic heating on the nightside, it is **crucial** to set :math:`(\Delta T)_{y}=T_{0}-T_{\text{strat}}` (this is **not default**, as :math:`(\Delta T)_{y}` is also used for the ``'Held_Suarez'`` heating, and the default value is set for use in that configuration). 

Each of the other parameters take the same meaning as for the ``'Held_Suarez'`` heating option (see above). 


| **Neutrally stratified exoplanet**
| When the ``'EXOPLANET2'`` option is selected, the equilibrium temperature profile takes the following form 

.. math::
   T_{\text{eq}} = \max\left[T_{\text{strat}},T_{\text{strat}}\cos\theta\left(\frac{p}{p_{0}}\right)^{\alpha}\right] 

where :math:`T_{\text{strat}}` is an isothermal temperature for the stratosphere (as above), :math:`p_{\text{trop}}` is the tropopause pressure, and the exponent :math:`\alpha` controls the stratification. 

When :math:`\alpha=2/7`, then this equilibrium temperature profile yields a troposphere that is neutrally stratified (if :math:`R/c_p=2/7`), that is, :math:`N=0` (:math:`N` is the buoyancy frequency). 

By default, :math:`T_{\text{strat}}=200\,\text{K}` (as above), :math:`p_{\text{trop}}=0.1\,\text{bar}` and :math:`\alpha=2/7`. These can all be changed in the namelist. 


| **Equilibrium temperature profile from file**
| When the ``'from_file'`` option is specified, then :math:`T_{\text{eq}}` is computed from a specified input file. The file should be a netcdf file, and temperature can be defined as a function of pressure, latitude, and longtiude. 

The zonal mean of :math:`T_{\text{eq}}` provided in the file is taken, and this is used as the relaxation temperature field for the Newtonian damping. 

The name of the netcdf file is specified with the namelist option ``equilibrium_t_file``, and should exclude the ``.nc`` suffix. Note: in the netcdf file, the :math:`T_{\text{eq}}` variable must be the same name as the name of the file, i.e., ``equilibrium_t_file``. The path to the file should be configured in the model's python run script (see the test cases for examples of how to do this). 




Top-down Newtonian damping 
----------------------------
Some text here \dots 



Local heating 
----------------------------
Some text here \dots 



Namelist options
----------------


The namelist options for **hs_forcing_nml** are listed below. 


**Namelist option to choose scheme**

:rad_scheme: String choosing the radiation scheme. Options are ``'FRIERSON'``, ``'BYRNE'``, ``'GEEN'``, ``'SCHNEIDER'``. Default option is ``'FRIERSON'``. 

**Namelist options for Frierson scheme longwave optical depth** 

:ir_tau_eq: Surface longwave optical depth at equator. Default :math:`6.0`. 
:ir_tau_pole: Surface longwave optical depth at pole. Default :math:`1.5`. 
:odp: Frierson optical depth scaling parameter :math:`\kappa`. Default :math:`1.0`. 
:linear_tau: :math:`f_l`. Determines partitioning between linear term and :math:`p^{k}` term in Frierson longwave optical depth. Default :math:`0.1`. 
:wv_exponent: Pressure exponent, :math:`k` in definition of optical depth. Default :math:`4.0`. 


Diagnostics 
----------------

Some diagnostics are listed below\dots 

References
----------

| [Held1994]_ 
| [Penn2018]_
| [VallisEtAl2018]_

Authors
----------
This documentation was written by Neil Lewis, peer reviewed by X, and quality controlled by Y. 