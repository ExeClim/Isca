
Radiative transfer: Simple gray and multi-band radiation schemes. 
=======================================================================================

Documentation for ``two_stream_gray_rad_mod``. 


Summary
-------
This module includes different simple configurations for solving the two stream radiative transfer equations. There are three semi-gray configurations available (gray in shortwave and longwave), and one intermediate multi-band configuration (gray in shortwave, multi-band in longwave). These configurations are selected by setting the ``rad_scheme`` namelist option. 

+------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``rad_scheme``   | Short description                                                                  | Reference            |
+------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'FRIERSON'``   | | Semi-gray scheme with prescribed longwave and shortwave                          | [Frierson2006a]_     |
|                  | | optical depths.                                                                  |                      |
+------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'BYRNE'``      | | Semi-gray scheme with longwave optical depth dependent on                        | | [Byrne2013]_       |
|                  | | water vapour content and :math:`\text{CO}_{2}` concentration. Shortwave          | | [VallisEtAl2018]_  |
|                  | | optical depth is prescribed.                                                     |                      |
+------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'GEEN'``       | | Multi-band scheme with two longwave bands and one shortwave                      | | [Geen2016]_        |
|                  | | band. One longwave band corresponds to an infrared window                        | | [VallisEtAl2018]_  |
|                  | | region (:math:`8-14\,\mu\,\text{m}`) and the second corresponds to all other     |                      |
|                  | | infrared wavelengths (:math:`>4\,\mu\,\text{m}`). Longwave and shortwave optical |                      |
|                  | | depths depend on water vapour content and :math:`\text{CO}_{2}` concentration.   |                      |
|                  | |                                                                                  |                      |
|                  | | Note: the Geen scheme is currently only stable when run in a fixed SST           |                      |
|                  | | configuration. A version stable with a slab ocean is expected by 30/11/2020.     |                      |
+------------------+-------------------------------+----------------------------------------------------+----------------------+
| ``'SCHNEIDER'``  | | Semi-gray scheme for use in **giant planet** simulations. Longwave               | [Schneider2009]_     |
|                  | | and shortwave optical depths are prescribed. Does not require a                  |                      |
|                  | | surface temperature as input, and allows specification of an                     |                      |
|                  | | interior heat flux.                                                              |                      |
+------------------+-------------------------------+----------------------------------------------------+----------------------+




For each of these schemes diurnally-and-seasonally averaged, and time-dependent options for the incoming solar radiation are available. 

In the next few sections we will describe the main options available to the user for solving the radiative transfer equations, and prescribing the insolation. Namelist options and their defaults, and available output diagnostics are summarised at the end of the documentation. 

Frierson/Byrne schemes 
----------------------

Each of these schemes solve the same equations to compute the longwave and shortwave fluxes. The difference between them lies in the definition of the longwave optical depths. 

| **Longwave**
| Longwave radiative transfer is described by the following two equations:

.. math::
   \frac{\text{d}I_{+}}{\text{d}\tau}=I_{+}-\sigma T^{4} \\
   \frac{\text{d}I_{-}}{\text{d}\tau}=\sigma T^{4} - I_{-}

:math:`I_{+}` and :math:`I_{-}` are the upward and downward longwave fluxes, and :math:`\tau` is the longwave optical depth. At the surface, :math:`I_{+}(p=p_{s}) = \sigma T_{s}^{4}`, and at the top of atmosphere :math:`I_{-}(p=0) = 0`. 

| **Shortwave**
| The downward shortwave flux :math:`I_{-}^{\ast}` is given by 

.. math::
    I_{-}^{\ast} = S\exp(-\tau^{\ast})

where :math:`\tau^{\ast}` is the shortwave optical depth, and :math:`S` is the incoming solar radiation. The implementation of the shortwave band in the Frierson and Byrne schemes assumes that there is no scattering within the atmosphere, and shortwave radiation is only absorbed from the downward stream. All reflection occurs at the surface. The upward flux of shortwave :math:`I_{+}^{\ast}` is not absorbed and escapes through the top of atmosphere, i.e. :math:`I_{+}^{\ast}=\alpha I_{-}^{\ast}(p=p_s)`, where :math:`\alpha` is the surface albedo. Note that the surface albedo is set by ``mixed_layer_mod``.

| **Longwave optical depth** 
| With ``rad_scheme`` set to ``'FRIERSON'``, the longwave optical depth is specified as a function of pressure 

.. math:: 
    \tau = \tau_{0}\left[f_{l}\left(\frac{p}{p_0}\right)+(1-f_{l})\left(\frac{p}{p_0}\right)^{k}\right]

where :math:`p_0=10^{5}\,\text{Pa}` is a reference pressure, :math:`f_{l}` sets how large the linear term is relative to the :math:`p^{k}` term. :math:`k` is a constant exponent (default 4) which determines the pressure dependence of optical depth. The linear term is included to reduce relaxation times in the stratosphere. The surface longwave optical depth :math:`\tau_{0}` is given by 

.. math:: 
    \tau_{0} = \kappa[\tau_{0e} + (\tau_{0p}-\tau_{0e})\sin^{2}\theta]

:math:`\tau_{0e}` is the value of the surface optical depth at the equation, and :math:`\tau_{0p}` is the value at the pole.  :math:`\kappa` is a parameter that can be used to scale the optical depth (default 1). 

When ``rad_scheme`` is set to ``'BYRNE'``, the longwave optical depth is calculated as a function of specific humidity :math:`q`, :math:`\text{CO}_{2}` mixing ratio, and pressure

.. math:: 
    \frac{\text{d}\tau}{\text{d}(p/p_{0})} = a\mu+bq+0.17\log(\text{CO}_{2}/360)

In the equation above, :math:`a` and :math:`b` are constant absorption coefficients, and :math:`\mu` is a scaling parameter intended to represent absorption by well-mixed gases (default value 1). 

| **Shortwave optical depth** 
| Both the Frierson and Byrne schemes define the shortwave optical depth as 

.. math:: 
    \tau^{\ast}= \tau^{\ast}_0\left(\frac{p}{p_0}\right)^{k^{\ast}}

where :math:`k^{\ast}` is a constant exponent that sets the pressure dependence of :math:`\tau^{\ast}`. :math:`\tau_{0}^{\ast}` is the shortwave optical depth at the surface, given by 

.. math:: 
    \tau^{\ast}_{0}=(1-\Delta\tau^{\ast}\sin^{2}\theta)*\tau_{0e}^{\ast}

where :math:`\tau_{0e}^{\ast}` is the surface optical depth at the equator, and :math:`\Delta\tau^{\ast}` determines the variation of shortwave optical depth with latitude.



Schneider giant planet scheme 
-----------------------------

This scheme is suitable to use for giant planet experiments where there is no solid surface at the bottom of the model atmosphere. It is selected by setting ``rad_scheme`` to ``'SCHNEIDER'``.


| **Longwave**
| The Schneider scheme solves the same equations for longwave radiative transfer as the Frierson/Byrne schemes. The only difference is in the lower boundary condition, where energy balance is enforced and the upward thermal radiative flux is set equal to the sum of the downward solar and thermal radiative fluxes, :math:`I_{+}(p=p_{s}) = I_{-}(p=p_{s}) + I_{-}^{\ast}(p=p_{\text{s}})`. An intrinsic heat flux from the planet's deep interior can also be specified, this is done in ``surface_flux_mod``. 

| **Shortwave**
| The downward shortwave flux :math:`I_{-}^{\ast}` is given by 

.. math::
   I_{-}^{\ast} = S(1-\alpha_{\text{b}})\exp(-\Gamma\tau^{\ast})

where 

.. math::
    \Gamma = 2\sqrt{1-\tilde{\omega}}\sqrt{1-\tilde{\omega}\gamma}

and 

.. math:: 
    \alpha_{\text{b}}=\frac{\sqrt{1-\tilde{\omega}\gamma}-\sqrt{1-\tilde{\omega}}}{\sqrt{1-\tilde{\omega}\gamma}+\sqrt{1-\tilde{\omega}}}

is the Bond albedo. Here, :math:`\tilde{\omega}=0.8` is the single-scattering albedo, and :math:`\gamma=1-2f_{b}` is the asymmetry factor, where :math:`f_{b}=0.398` is the fraction of radiation back scattered. 

The only reflection of shortwave radiation occurs at the top of atmosphere where a fraction :math:`(1-\alpha_{\text{b}})` is removed from the incoming solar radiation. 

There is no modelled upward flux of shortwave radiation in the atmosphere :math:`I_{+}^{\ast}=0`. Instead any shortwave radiation that reaches the bottom of the atmosphere is re-added to the upward longwave beam as a flux through the lower boundary. 

| **Longwave optical depth** 
| For the Schneider scheme, longwave optical depth is a simple function of pressure 

.. math:: 
    \tau = \tau_{0,\text{gp}}\left(\frac{p}{p_{0}}\right)^{k_{\text{gp}}}

where :math:`\tau_{0,\text{gp}}=80.0` is the longwave optical depth at pressure :math:`p_{0}` and :math:`k_{\text{gp}}=2.0` is a constant exponent that sets the pressure dependence of :math:`\tau`. 

| **Shortwave optical depth** 
| The shortwave optical depth takes a similar functional form to the longwave optical depth 

.. math:: 
    \tau^{\ast} = \tau^{\ast}_{0,\text{gp}}\left(\frac{p}{p_{0}}\right)^{k^{\ast}_{\text{gp}}}

where :math:`\tau^{\ast}_{0,\text{gp}}=3.0` is the shortwave optical depth at pressure :math:`p_{0}` and :math:`k^{\ast}_{\text{gp}}=1.0` is a constant exponent that sets the pressure dependence of :math:`\tau^{\ast}`. 

| **Note on input parameters for this scheme** 
| At present, each of the input parameters :math:`\tilde{\omega}`, :math:`f_{b}`, :math:`\tau_{0,\text{gp}}`, :math:`k_{\text{gp}}`, :math:`\tau^{\ast}_{0,\text{gp}}`, and  :math:`k^{\ast}_{\text{gp}}` are hardcoded and are not available as namelist options. There default values are given in the description of the scheme above. The code could be easily modified if a user wished to vary these parameters from the namelist. 


Geen scheme
-----------

The Geen scheme provides an intermediate option between gray radiation and more complete descriptions of radiative transfer (e.g., the correlated-:math:`k` schemes SOCRATES and RRTM). It has two infrared bands and one solar band. The shortwave band (:math:`<4\,\mu\,\text{m}`) treats all solar radiation. Two long-wave bands treat absorption: one in the infrared region of the spectral (:math:`8-14\,\mu\,\text{m}`), and the other in all other longwave wavelengths (:math:`<4\,\mu\,\text{m}`, non-window).  

Note: the Geen scheme is currently only stable when run in a fixed SST configuration. A version stable with a slab ocean is expected by 30/11/2020. 

| **Longwave**
| Longwave radiative transfer is described by the following set of equations. In the non-window region:

.. math::
   \frac{\text{d}I_{+}^{\text{nw}}}{\text{d}\tau^{\text{nw}}}=I_{+}^{\text{nw}}-R^{\text{nw}}\sigma T^{4} \\
   \frac{\text{d}I_{-}^{\text{nw}}}{\text{d}\tau^{\text{nw}}}=R^{\text{nw}}\sigma T^{4} - I_{-}^{\text{nw}}

and in the window region: 

.. math::
   \frac{\text{d}I_{+}^{\text{win}}}{\text{d}\tau^{\text{win}}}=I_{+}^{\text{win}}-R^{\text{win}}\sigma T^{4} \\
   \frac{\text{d}I_{-}^{\text{win}}}{\text{d}\tau^{\text{win}}}=R^{\text{win}}\sigma T^{4} - I_{-}^{\text{win}}

The superscripts :math:`^{\text{nw}}` and :math:`^{\text{win}}` refer to the non-window and window regions, respectively. :math:`R^{\text{win}}` and :math:`R^{\text{nw}}=1-R^{\text{win}}` are the fraction irradiances in the non-window and window regions. 


| **Shortwave**
| The Geen scheme solves the same equations for shortwave radiative transfer as the Frierson/Byrne schemes. The only difference is in the specification of the shortwave optical depth (see below). 

As with the Frierson and Byrne schemes, the implementation of the shortwave band assumes that there is no scattering within the atmosphere, and shortwave radiation is only absorbed from the downward stream. All reflection occurs at the surface. The upward flux of shortwave :math:`I_{+}^{\ast}` is not absorbed and escapes through the top of atmosphere. 

| **Longwave optical depth** 
| The longwave optical depth in the Geen scheme is defined seperately for the non-window and window regions. In the window region 

.. math:: 
    \frac{\text{d}\tau^{\text{win}}}{\text{d}(p/p_{0})} = a_{\text{win}}+b_{\text{win}}q+c_{\text{win}}q^{2}+0.0954\log(\text{CO}_{2}/360)

and in the non-window region 

.. math:: 
    \frac{\text{d}\tau^{\text{nw}}}{\text{d}(p/p_{0})} = a_{\text{nw}}+b_{\text{nw}}\log(c_{\text{nw}}q+1)+0.2023\log(\text{CO}_{2}/360)

where the :math:`a`, :math:`b`, and :math:`c`'s are constant absorption coefficients that may be specified in the namelist. The default values for these coefficients were fitted to output from `Santa Barbara DISORT Atmospheric Radiative Tranfer 60` (SBDART) [Ricchiazzi1998]_. 

| **Shortwave optical depth** 
| In the shortwave, the optical depth in the Geen scheme is specified as 

.. math:: 
    \frac{\text{d}\tau^{\ast}}{\text{d}(p/p_{0})} = a_{\text{sw}}+b_{\text{sw}}(\tau^{\ast})q+c_{\text{sw}}\log(\text{CO}_{2}/360)

where 

.. math::
    \log[b_{\text{sw}}(\tau^{\ast})]=\frac{0.01887}{\tau^{\ast}+0.009522}+\frac{1.603}{(\tau^{\ast}+0.5194)^{2}}. 

:math:`a_{\text{sw}}=0.0596` and :math:`c_{\text{sw}}=0.0029`. As it stands, :math:`a_{\text{sw}}`, :math:`b_{\text{sw}}`, and :math:`c_{\text{sw}}` are hardcoded and not modifiable via the namelist. As with the longwave absorption coefficients, the default values for the absorption coefficients were fitted to SBDART output. 

Incoming solar radiation
------------------------

For each of the radiative transfer schemes, diurnally-and-seasonally averaged, and time-dependent options for the incoming solar radiation are avaialable. The user selects whether diurnally-and-seasonally averaged or time-dependent solar forcing is used by setting the namelist option ``do_seasonal``. 

| **Diurnally-and-seasonally averaged insolation** 
| Diurnally-and-seasonally averaged insolation is selected by setting ``do_seasonal = False``. 

For the Frierson, Byrne, and Geen schemes, setting ``do_seasonal = False`` imposes a :math:`P_{2}` (second legendre polynomial) insolation profile, which is designed to approximate the Earth's seasonally averaged insolation distribution. When this option is selected, there is no diurnal cycle (i.e. the forcing is fully time-independent). The incoming solar radiation then takes the following form: 

.. math:: 
    S &= \frac{S_{0}}{4}[1+\Delta_{S}P_{2}(\theta)+\Delta_{\text{sw}}\sin\theta] \\ 
    P_{2} &= (1 - 3\sin^{2}\theta)/4 

where :math:`S_{0}` is the solar constant. :math:`\Delta_{S}` is used to set the amplitude of the :math:`P_{2}` insolation profile between the equator and pole, and :math:`\Delta_{\text{sw}}` (default 0) can be used to further modify this with a :math:`\sin\theta` profile. When :math:`\Delta_{\text{sw}}=0`, the insolation difference between the equator and pole is :math:`\Delta S= 3\Delta_{S}/4\times S_{0}/4`.

For the Schneider giant planet scheme, setting ``do_seasonal = False`` imposes an insolation profile that varies with :math:`\cos\theta`,

.. math:: 
    S = \frac{S_{0}}{\pi}\cos\theta 

which corresponds to the latitudinal distribution of radiation received by a planet with no obliquity (perpetual equinox). As with the simple :math:`P_{2}` insolation for the Frierson, Byrne, and Geen schemes, there is no diurnal cycle. 

| **Time-dependent insolation** 
| All of the radiation schemes contained in this module can be run with time-dependent insolation, selected by setting ``do_seasonal = True``. In this case, ``astronomy_mod`` is used to calculate the zenith angle and planet-star distance at each location and time. This information is then used by the ``two_stream_gray_rad`` module to calculate the top of atmosphere insolation as a function of time. 

``astronomy_mod`` calculates the zenith angle and planet-star distance as a function of the following orbital and planetary parameters: obliquity, eccentricity, semi-major axis, longitude of perihelion (w.r.t NH autumn equinox), orbital period, and rotation rate. These options must be specified to ``astronomy_mod_nml``. An input parameter in ``two_stream_gray_rad``, ``equinox_day`` determines the time of year when NH autumn equinox occurs. 

``two_stream_gray_rad`` then calculates the insolation as 

.. math:: 
    S = S_{0}\cos\zeta\left(\frac{a}{r}\right)^{2}

where :math:`\zeta` is the zenith angle, :math:`a` is the semi-major axis of the orbital ellipse, and :math:`r` is the time-varying planet-star distance. NB: :math:`(a/r)^{2}` is called ``rrsun`` in the code. 

If the user wishes, they may average incoming solar radiation over a period :math:`\Delta t_{\text{avg}}` (units :math:`\text{s}`). This allows the user to have seasonally varying forcing without a diurnal cycle, for example. To achieve diurnally averaged insolation for a planet with the Earth's length of day one would set :math:`\Delta t_{\text{avg}}=86400.0\,\text{s}`. 

There is also an option to run perpetually on one day by setting the namelist variable ``solday``. For example if Northern Hemisphere autumn equinox was set to occur on day :math:`270` of a :math:`360` day year, then one could run a perpetual solstice simulation by setting ``solday=180``. This can be used in conjunction with an appropriate choice for :math:`\text{d}t_{\text{avg}}` to remove the diurnal cycle in such an experiment. 



Namelist options
----------------


The namelist options for **two_stream_gray_rad_nml** are listed below. 


**Namelist option to choose scheme**

:rad_scheme: String choosing the radiation scheme. Options are ``'FRIERSON'``, ``'BYRNE'``, ``'GEEN'``, ``'SCHNEIDER'``. Default option is ``'FRIERSON'``. 

**Namelist options for Frierson scheme longwave optical depth** 

:ir_tau_eq: Surface longwave optical depth at equator. Default :math:`6.0`. 
:ir_tau_pole: Surface longwave optical depth at pole. Default :math:`1.5`. 
:odp: Frierson optical depth scaling parameter :math:`\kappa`. Default :math:`1.0`. 
:linear_tau: :math:`f_l`. Determines partitioning between linear term and :math:`p^{k}` term in Frierson longwave optical depth. Default :math:`0.1`. 
:wv_exponent: Pressure exponent, :math:`k` in definition of optical depth. Default :math:`4.0`. 

**Namelist options for Byrne scheme longwave optical depth** 

:bog_a: Absorption coefficient :math:`a` in Byrne longwave optical depth. Default :math:`0.8678`. 
:bog_b: Absorption coefficient :math:`b` in Byrne longwave optical depth. Default 1997.9. 
:bog_mu: Scaling parameter :math:`\mu` in Byrne longwave optical depth. Default 1.0. 

**Namelist options for Frierson/Byrne scheme shortwave optical depth** 

:atm_abs: Shortwave optical depth at the equator :math:`\tau_{0e}^{\ast}`. Default :math:`0.0`. 
:sw_diff: Amplitude of latitudinal optical depth variation :math:`\Delta\tau^{\ast}`.  Default :math:`0.0`. 
:solar_exponent: Pressure exponent :math:`k^{\ast}`. Default value is :math:`4.0`.

| **Namelist options for Geen multi-band scheme** 
| Note: for the Geen scheme, the input parameters determining the shortwave optical depth are hardcoded and cannot be set from the namelist (see scheme description above for their default values). Thus, only parameters determining the longwave optical depth are listed here. 


:window: Window fraction :math:`R^{\text{win}}` for longwave radiative transfer. Default value :math:`0.3732`. 
:ir_tau_co2_win: Absorption coefficient :math:`a_{\text{win}}` in longwave optical depth. Default value :math:`0.2150`.
:ir_tau_wv_win1: Absorption coefficient :math:`b_{\text{win}}` in longwave optical depth. Default value :math:`147.11`. 
:ir_tau_wv_win2: Absorption coefficient :math:`c_{\text{win}}` in longwave optical depth. Default value :math:`1.0814\times10^{4}`. 
:ir_tau_co2: Absorption coefficient :math:`a_{\text{nw}}` in longwave optical depth. Default value :math:`0.1`. 
:ir_tau_wv1: Absorption coefficient :math:`b_{\text{nw}}` in longwave optical depth. Default value :math:`23.8`. 
:ir_tau_wv2: Absorption coefficient :math:`c_{\text{nw}}` in longwave optical depth. Default value :math:`254.0`. 


| **Namelist options for Schneider giant planet scheme** 
| Note: for the Schneider scheme, the input parameters determining shortwave and longwave optical depth are currently hardcoded and cannot be set from the namelist (see scheme description above for their default values). Thus, these parameters are not listed here. 

:diabatic_acce: Multiplicative scaling factor for temperature tendency due to radiation in Schneider scheme. This can be used to speed up model spin-up. It will speed up the model spin-up if greater than 1. Default value is 1. 

**Namelist options for incoming solar radiation**

:do_seasonal: Sets whether insolation is diurnally-and-seasonally averaged (``FALSE``), or time dependent (``TRUE``). Default ``FALSE``. 
:solar_constant: The solar constant :math:`S_{0}`. Default value: :math:`1360.0\,\text{W}\,\text{m}^{-2}`. 

The following namelist options are used when ``rad_scheme`` is set to ``'FRIERSON'``, ``'BYRNE'`` or ``'GEEN'`` and ``do_seasonal=FALSE``:

:del_sol:  Parameter :math:`\Delta_{S}` determining :math:`P_{2}` insolation amplitude. Default value: :math:`1.4`.  
:del_sw: Parameter :math:`\Delta_{\text{sw}}` defining magnitude of :math:`\sin\theta` modification to :math:`P_{2}` insolation profile. Default value: :math:`0.0`. 

The following namelist options are used when ``do_seasonal=TRUE``: 

:use_time_average_coszen: ``TRUE`` or ``FALSE``. If ``TRUE``, average :math:`\cos\zeta` (:math:`\zeta` is zenith angle) over the period ``dt_rad_avg``. For example, for the Earth's diurnal period, ``use_time_average_coszen=TRUE`` and ``dt_rad_avg=86400.`` would achieve diurnally averaged insolation. 
:dt_rad_avg: Averaging period (units: seconds) for time-dependent insolation :math:`\Delta t_{\text{avg}}`. Default=-1 sets averaging period to model timestep. 
:solday: Day of year to run time-dependent insolation perpetually. If negative, the option to run perpetually on a specific day is not used. Default -10. 
:equinox_day: Fraction of year [0,1] where Northern Hemisphere autumn equinox occurs. Default = 0.75 (e.g. end of September for 360 day year). 

**Namelist options for setting carbon dioxide concentration** 

:do_read_co2: ``TRUE`` or ``FALSE``. If ``TRUE``, reads time-varying :math:`\text{CO}_{2}` concentration from an input file [needs to be 4D (3 spatial dimensions and time), but no spatial variation should be defined (the code only reads in maximum value at a given time)]. Default ``FALSE``. 
:carbon_conc: Prescribed :math:`\text{CO}_{2}` (units: ppmv). Used if ``do_read_co2=FALSE``. Default value :math:`360.0\,\text{ppmv}`. 
:co2_file: Name of :math:`\text{CO}_{2}` file to read (note, should be specified without .nc appendix). Default ``'co2'``. 
:co2_variable_name: Name of :math:`\text{CO}_{2}` variable in :math:`\text{CO}_{2}` file. Default ``'co2'``.  

**Important parameters not set in two_stream_gray_rad_nml**

:pstd_mks: This is used as the reference pressure :math:`P_{0}` in the Frierson/Byrne/Schneider shortwave optical depth, and the Frierson/Schneider longwave optical depth. It is set in ``constants_mod_nml``. Default value is :math:`10^{5}\,\text{Pa}`. Note: this should be changed to :math:`3\times10^{5}\,\text{Pa}` for the giant planet configuration. 



:albedo_value: The surface albedo :math:`\alpha` used by the Frierson/Byrne/Geen schemes at the lower boundary is set in ``mixed_layer_nml``. Default value is :math:`0.06` for a simple homogeneous slab ocean surface. The albedo can vary spatially if land or ice is introduced. For more details, see the documentation for ``mixed_layer_mod``. 

:flux_heat_gp: A prescribed heat flux through the lower boundary can be added for the giant planet case. This is set in ``surface_flux_nml``. Default value is :math:`5.7\,\text{W}\,\text{m}^{-2}`. 

The following astronomical parameters are set in ``astronomy_mod_nml``. They are used if ``do_seasonal=True``. 

:ecc: Orbital eccentricity. Default :math:`0.0`. 
:obliq: Obliquity. Default :math:`23.439` degrees.
:per: Longitude of perihelion (point in orbit when planet is closest to star) with respect to autumnal equinox in Northern Hemisphere. Default :math:`102.932` degrees. 

The following astronomical parameters are set in ``constant_nml``. They will be used to calculate the diurnal period if ``do_seasonal=True``. 

:orbital_period: Orbital period in seconds. Default is :math:`365.25\times86400.0\,\text{s}`. Only used if ``calendar`` is set to ``'no_calendar'`` in ``main_nml``. 
:omega: Planetary rotation rate in :math:`s^{-1}`. Default value is :math:`7.29\times10^{-5}\,\text{s}^{-1}` 



Diagnostics
-----------

These are the diagnostics associated with the ``two_stream_gray_rad`` module. 


+-------------------+-------------------------------------+------------------------------------+
| Name              | Description                         | Units                              |
+===================+=====================================+====================================+
| olr               | Outgoing longwave radiation.        | :math:`\text{W}\,\text{m}^{-2}`    |
+-------------------+-------------------------------------+------------------------------------+
| swdn_sfc          | Absorbed shortwave at surface.      | :math:`\text{W}\,\text{m}^{-2}`    |
+-------------------+-------------------------------------+------------------------------------+
| swdn_toa          | Shortwave flux down at top of       | :math:`\text{W}\,\text{m}^{-2}`    |
|                   | atmosphere.                         |                                    |
+-------------------+-------------------------------------+------------------------------------+
| net_lw_surf       | Net upward longwave flux at         | :math:`\text{W}\,\text{m}^{-2}`    |
|                   | surface.                            |                                    |
+-------------------+-------------------------------------+------------------------------------+
| lwdn_sfc          | Longwave flux down at surface.      | :math:`\text{W}\,\text{m}^{-2}`    |
+-------------------+-------------------------------------+------------------------------------+
| lwup_sfc          | Longwave flux up at surface.        | :math:`\text{W}\,\text{m}^{-2}`    |
+-------------------+-------------------------------------+------------------------------------+
| tdt_rad           | Temperature tendency due to         | :math:`\text{K}\,\text{s}^{-1}`    |
|                   | radiation.                          |                                    |
+-------------------+-------------------------------------+------------------------------------+
| tdt_solar         | Temperature tendency due to         | :math:`\text{K}\,\text{s}^{-1}`    |
|                   | solar radiation.                    |                                    |
+-------------------+-------------------------------------+------------------------------------+
| flux_rad          | Total radiative flux (positive      | :math:`\text{W}\,\text{m}^{-2}`    |
|                   | up).                                |                                    |
+-------------------+-------------------------------------+------------------------------------+
| flux_lw           | Net longwave radiative flux         | :math:`\text{W}\,\text{m}^{-2}`    |
|                   | (positive up).                      |                                    |
+-------------------+-------------------------------------+------------------------------------+
| flux_sw           | Net shortwave radiative flux        | :math:`\text{W}\,\text{m}^{-2}`    |
|                   | (positive up).                      |                                    |
+-------------------+-------------------------------------+------------------------------------+
| coszen            | Cosine of zenith angle.             | none                               |
+-------------------+-------------------------------------+------------------------------------+
| fracsun           | Daylight fraction of time           | none                               |
|                   | interval.                           |                                    |
+-------------------+-------------------------------------+------------------------------------+
| co2               | Carbon dioxide concentration.       | :math:`\text{ppmv}`                |
+-------------------+-------------------------------------+------------------------------------+
| lw_dtrans         | Longwave (non-window)               | none                               |
|                   | transmission.                       |                                    |
+-------------------+-------------------------------------+------------------------------------+
| lw_dtrans_win     | | Longwave window transmission.     | none                               |
|                   | | Note: only for ``'GEEN'`` scheme. |                                    |
+-------------------+-------------------------------------+------------------------------------+
| sw_dtrans         | | Shortwave transmission.           | none                               |
|                   | | Note: only for ``'GEEN'`` scheme. |                                    |
+-------------------+-------------------------------------+------------------------------------+




Relevant modules and subroutines
--------------------------------

Modules relevant to this one include: 

:astronomy_mod: Module that performs astronomical calcuations used for insolation. 
:mixed_layer_mod: Surface albedo is set here. This is also where the surface temperature is updated. 
:surface_flux_mod: An internal heat flux for giant planets can be set here. 
:constants_mod: Planetary rotation rate and orbital period are set here. These are used in the calculations made by ``astronomy_mod``. 


Other radiative transfer schemes are included in the following modules:

:rrtm_radiation: Correlated-:math:`k` scheme tuned for Earth-like applications. 
:socrates_interface_mod: Interface for flexible Met-Office correlated-:math:`k` scheme used for Earth-like and exoplanetary atmospheres. 


References
----------

| [Byrne2013]_ 
| [Frierson2006a]_
| [Geen2016]_
| [Schneider2009]_
| [VallisEtAl2018]_


Authors
----------
This documentation was written by Neil Lewis, peer reviewed by Ruth Geen, and quality controlled by Matthew Henry.

