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

A further ``equilibrium_t_option``, ``'top_down'`` is available, which constructs :math:`T_{\text{eq}}` from astronomical solar input and an approximate analytic solution to radiative-convective equations with a specified optical depth, lapse rate, radiative relaxation time, and surface mixed-layer depth. In this scenario, the subroutine ``top_down_newtonian_damping`` is used in place of ``newtonian_damping``. This option will be described in the section **Top down Newtonian damping**. 

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

The subroutine ``top_down_newtonian_damping`` is used when ``equilibrium_t_option`` is set to ``'top_down'``. In this scenario, :math:`T_{\text{eq}}` is constructed using astronomical solar input and an approximate analytic solution to radiative-convective equations with a specified optical depth, lapse rate, radiative relaxation time, and surface mixed-layer depth, following [VallisEtAl2018]_. 

:math:`T_{\text{eq}}` is computed as follows. It is assumed that the atmosphere consists of a troposphere, with a given lapse rate, and a stratosphere that has a small optical depth and is in radiative equilibrium. Then, the a radiative-convective tropopause height is obtained by solving [VallisEtAl2015]_ 

.. math:: 
   H_{\text{T}}=\frac{1}{16\Gamma}\left(CT_{\text{T}}+\sqrt{C^{2}T_{\text{T}}^{2}+32\Gamma\tau_{\text{s}}H_{\text{a}}T_{\text{T}}}\right) 

where :math:`\Gamma` is the lapse rate, :math:`\tau_{\text{s}}`, and :math:`H_{\text{a}}` is the scale height of the main infrared absorber. Each of these are parameters that may be specified by the user in the namelist. :math:`C=\log4` is a constant. 

:math:`T_{\text{T}}` is the temperature at the tropopause, which, given the assumptions decsribed above [VallisEtAl2015]_, is defined in terms of the insolation via 

.. math:: 
   T_{\text{T}} &= 2^{-\frac{1}{4}}T_{\text{e}} \\ 
   T_{\text{e}} &= \left[\frac{(1-\alpha)S}{\sigma}\right]^{\frac{1}{4}}

where :math:`S` is the insolation, which varies in time and space, and depends on astronomical parameters such as the solar constant, orbital radius, obliquity, eccentricity, and rotation rate. :math:`\alpha` is the surface albedo. 

Once the height of the tropopause has been calculated, a surface temperature is calculated using 

.. math:: 
    T_{\text{s}} = T_{\text{T}} + H_{\text{T}}\Gamma 

This surface temperature is then used to force the actual ground temperature via 

.. math:: 
   C_{\text{g}}\frac{\text{d}T_{\text{g}}}{\text{d}t}=\sigma T_{\text{s}}^{4} - \sigma T_{\text{g}}^{4}

where :math:`C_{\text{g}}` is the surface heat capacity. 

Once the ground temperature has been updated, the relaxation temperature profile :math:`T_{\text{eq}}` can be calculated. In the troposphere 

.. math:: T_{\text{trop}} = T_{\text{g}}-\Gamma z 

There then various options available to determine the stratospheric temperature structure. Each of these are selected by setting the namelist parameter ``stratosphere_t_option``. 

**The default option is** obtained with ``stratosphere_t_option='extend_tp'``. In this scenario, :math:`T_{\text{eq}}` is given by  

.. math:: 
   T_{\text{eq}} = \begin{cases} T_{\text{trop}} &\text{if}\  z\leq H_{\text{T}} \\ T_{\text{T}} &\text{if}\  z > H_{\text{T}} \end{cases}

where :math:`T_{\text{T}}` is the tropopause temperature given above. 

Other **experimental options** are as follows:

If ``stratosphere_t_option='c_above_tp'`` then 

.. math:: 
   T_{\text{eq}} = \begin{cases} T_{\text{trop}} &\text{if}\  z\leq H_{\text{T}} \\ T_{\text{str}} &\text{if}\  z > H_{\text{T}} \end{cases}

where :math:`T_{\text{str}} = T_{\text{strat}} - \epsilon\sin\theta` (as in the Newtonian damping configurations above). This option may lead to a temperature discontinuity at the tropopause.

If ``stratosphere_t_option='hs_like'`` then 

.. math:: 
   T_{\text{eq}} = \max\left(T_{\text{trop}}, T_{\text{str}}\right) 

This choice is very similar to that in the Held-Suarez framework, but where the tropospheric radiative-convective equilibrium profile prescribed by Held and Suarez is replaced with :math:`T_{\text{trop}}`. 

Finally, if ``stratosphere_t_option=''`` (or, in fact, is set to anything other than the three preceding options), then it is assumed that the atmosphere is all troposphere, i.e., 

.. math:: 
   T_{\text{eq}} = \max\left(T_{\text{trop}}, 0\right) 



|
|
**For all of the above options**, once :math:`T_{\text{eq}}` has been obtained, the model temperature is relaxed towards the equilibrium temperature in an indentical manner to that in the regular ``newtonian_damping`` subroutine (see above). 




Local heating 
----------------------------
Some text here \dots 



Namelist options
----------------


The namelist options for **hs_forcing_nml** are listed below. 



Diagnostics 
----------------

The diagnostics available from **hs_forcing_mod** are listed below. 

References
----------

| [Held1994]_ 
| [Penn2018]_
| [VallisEtAl2015]_
| [VallisEtAl2018]_

Authors
----------
This documentation was written by Neil Lewis, peer reviewed by X, and quality controlled by Y. 