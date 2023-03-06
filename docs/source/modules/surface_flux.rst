
Surface flux calculation
=======================================================================================

Documentation for ``surface_flux`` and ``gp_surface_flux``, both found at ``Isca/src/coupler/surface_flux.F90``. 


Summary
-------

This module manages the exchange of heat, momentum and moisture at the planet's surface. 
A brief overview of how the exchange of these is calculated and the options available is given first, with a full summary of the namelist parameters below.


Parametrization for planets with a rock/liquid surface
------------------------------------------------------
Fluxes of heat (:math:`S`), momentum (:math:`\mathcal{T}`) and humidity (:math:`E`) are calculated following the standard drag equations set out in [Frierson2006a].

.. math::
    S = \rho_{a} c_{p} C_{S} |\mathbf{v}_{a}|(\Theta_{a} - \Theta_{s})

.. math::
	\mathcal{T} = \rho_{a} C_{\mathcal{T}} |\mathbf{v}_{a}|\mathbf{v}_{a}
	
.. math::
    E = \rho_{a} C_{E} |\mathbf{v}_{a}|(q_{a} - q_{s}^{*}).
	
In the above, :math:`\rho_{a}`, :math:`|\mathbf{v}_{a}`, :math:`\Theta_{a}`, and :math:`q_{a}` are the density, horizontal wind, potential temperature and specific humidity calculated at the lowest model level. :math:`c_{p}` is the heat capacity at constant pressure, :math:`\Theta_{s}` is the surface potential temperature and :math:`q_{s}^{*}` is the saturation specific humidity at the surface temperature. 

:math:`C_{S}`, :math:`C_{\mathcal{T}}` and :math:`C_{E}` are drag coefficients. These are calculated according to a simplified Monin-Obukhov similarity theory [Frierson2006a] in the module ``monin_obukhov_mod``, found in ``Isca/src/atmos_param/monin_obukhov/monin_obukhov.F90``. Roughness lengths for each flux can be specified as namelist parameters in ``idealised_moist_phys``.

Land
````
Land is implemented in Isca in 4 ways, via adjustment of (1) the roughness length, (2) the evaporative flux, (3) the albedo, (4) the mixed layer depth:

(1) Any desired adjustments to roughness length over land are made in ``idealized_moist_phys_mod`` via the namelist parameter ``land_roughness_prefactor`` and are then fed to ``surface_flux_mod``. This parameter prescribes a prefactor to modify the roughness lengths used over land in the Monin-Obukhov calculations. The same modification prefactor is used for all fluxes.

(2) Adjustments to the evaporative flux are made in ``surface_flux_mod``. 2 options are available:

     (a) Prefactors :math:`\alpha` and :math:`\beta` may be prescribed over land to modify the evaporative flux:

     .. math::
         E = \beta \rho_{a} C_{E} |\mathbf{v}_{a}|(q_{a} - \alpha q_{s}^{*}).
	
     :math:`\alpha` is set using ``land_humidity_prefactor`` and :math:`\beta` is set using ``land_evap_prefactor`` in the ``surface_flux_nml``. 


     (b) A simple bucket hydrology may be used following [Manabe1969]. In this case, hydrology is described by prescribing land with a bucket depth :math:`W` which can vary between 0, corresponding to an empty bucket, and a field capacity :math:`W_{FC}`, corresponding to a full bucket. The evaporation equation is then modified over land grid cells according to the bucket depth.

     If :math:`W \lt 0.75 W_{FC}`:
     
     .. math::
         E = \frac{W}{0.75 W_{FC}} \rho_{a} C_{E} |\mathbf{v}_{a}|(q_{a} - q_{s}^{*})
     	
     If :math:`W \ge 0.75 W_{FC}`:
     
     .. math::
         E = \rho_{a} C_{E} |\mathbf{v}_{a}|(q_{a} - q_{s}^{*})
     
     Bucket namelist parameters are input via ``idealized_moist_phys_mod``.

3.& 4. Adjustments to the albedo and mixed layer depth are made via namelist parameters in the module :mixed_layer_mod: found in ``Isca/src/atmos_spectral/driver/solo/mixed_layer.F90``.


Parametrization for gas giants
------------------------------

For modelling gas giants, the namelist option ``flux_heat_gp`` allows the user to prescribe an intrinsic heat flux at 1 bar, appropriate for a giant planet atmosphere where there is no surface [Schneider2009]. In practice, the subroutine ``gp_surface_flux`` will then impose an atmospheric heating :math:`Q_{gp}` at the lowest model level as:

.. math::
    Q_{gp} = \frac{g \gamma}{c_{p} \Delta p}

where :math:`g` is the gravitational acceleration, :math:`\gamma` is the value of ``flux_heat_gp``, :math:`c_{p}` is the heat capacity at constant pressure and :math:`\Delta p` is the thickness of the lowest model pressure level. To accelerate model spin-up, this intrinsic heat flux can be increased via the namelist option ``diabatic_acce``. Alternatively, Isca's column model configuration could be run to equilibrium.
	
To turn on this parametrization, ``gp_surface`` should be set to True in the ``idealized_moist_phys`` namelist. 


Namelist options
----------------

The namelist options for ``surface_flux_nml`` are listed below. These include options to limit evaporation over land, to use more simplified parametrisations, to account for salinity over ocean, and to help with issues resulting from truncation. We note that Isca is developed from the GFDL FMS model (https://github.com/NOAA-GFDL/FMS) and some options remain to enable old model configurations to be supported, as indicated in the table. 

+---------------------------+----------+-----------------------------------------------------------------------------------------+
| Name                      | Default  | Description                                                                             |
+===========================+==========+=========================================================================================+
|``no_neg_q``               | False    | If q_atm_in (specific humidity input) is negative (because of numerical truncation),    |
|                           |          | then override with 0.                                                                   |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``use_virtual_temp``       | True     | If true, use virtual potential temp to calculate the stability of the surface layer.    |
|                           |          | If false, use potential temp.                                                           |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``alt_gustiness``          | False    | If true, use an alternative formulation for gustiness calculation where a minimum bound | 
|                           |          | on the wind speed is used in flux calculations, with the bound equal to ``gust_const``. |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``gust_const``             | 1.0      | Constant for alternative gustiness calculation.                                         |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``gust_min``               | 0.       | Minimum gustiness used when alt_gustiness is set to False.                              |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``old_dtaudv``             | False    | If true, the derivatives of surface wind stress w.r.t. the zonal wind and meridional    |
|                           |          | wind are each approximated by the same value.                                           |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``use_mixing_ratio``       | False    | GFDL LEGACY: If true, then use the saturation surface mixing ratio instead of specific  |
|                           |          | humidity in calculating surface moisture exchange. This provides the capability to run  |
|                           |          | the Manabe Climate form of the surface flux.                                            |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``do_simple``              | False    | If true, then approximate saturation surface specific humidity as:                      |
|                           |          | :math:`q_{sat} = 0.622*e_{sat} / p_{surf}`.                                             |
|                           |          | If false:    :math:`q_{sat} = 0.622*e_{sat} / (p_{surf} - e_{sat})` is used.            |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``ncar_ocean_flux``        | False    | Use NCAR climate model turbulent flux calculation described by [Large2004].             |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``ncar_ocean_flux_orig``   | False    | GFDL LEGACY: Use NCAR climate model turbulent flux calculation described by [Large2004],|
|                           |          | using the original GFDL implementation, which contains a bug in the specification of the|
|                           |          | exchange coefficient for the sensible heat. Not recommended for new experiments.        |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``raoult_sat_vap``         | False    | If true, reduce saturation vapor pressures over ocean to account for seawater salinity. |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``land_humidity_prefactor``| 1.0      | Factor that multiplies the surface specific humidity over land. This is included to make|
|                           |          | land 'dry'. If it is equal to 1, land behaves like ocean. If it is between 0 and 1,     |
|                           |          | this will decrease the evaporative heat flux in areas of land. Note that this can lead  |
|                           |          | to sign changes in the evaporative flux, and we find this becomes unstable over very    |
|                           |          | shallow mixed layer depths.                                                             |
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``land_evap_prefactor``    | 1.0      | Factor that multiplies the evaporative flux over land. This is included to make land    |
|                           |          | 'dry'. If it is equal to 1, land behaves like ocean. If it is between 0 and 1, this will|
|                           |          | decrease the evaporative heat flux in areas of land. This formulation avoids sign       |
|                           |          | changes in the evaporative flux and remains stable over very shallow mixed layer depths.|
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``flux_heat_gp``           | 5.7      | Intrinsic heat flux imposed to describe the lower boundary of a giant planet atmosphere.|
|                           |          | [Schneider2009]. Default is :math:`5.7 Wm^{-2}`, appropriate for Jupiter [Gierasch2000].|            
+---------------------------+----------+-----------------------------------------------------------------------------------------+
|``diabatic_acce``          | 1.0      | Multiplicative scaling factor for the temperature tendency in the giant planet scheme.  |
|                           |          | This can be used to speed up model spin-up. It will speed up the model spin-up if       |
|                           |          | greater than 1. A similar factor can also be applied for giant planet radiation. Note   |
|                           |          | that an alternative way to accelerate spin-up is to run Isca's column model to          |
|                           |          | equilibrium.                                                                            |
+---------------------------+----------+-----------------------------------------------------------------------------------------+



Diagnostics
-----------

No diagnostics are saved within the ``surface_flux`` module itself. Diagnostics related to surface exchange are available via the `atmosphere` diagnostic list in ``idealized_moist_phys.F90``, and the `mixed_layer` diagnostic list in ``mixed_layer.F90``.


Relevant modules and subroutines
--------------------------------

Modules relevant to this one include: 

:idealized_moist_phys: Momentum flux diagnostics are output here. A land roughness prefactor can also be applied in this module, which is used to alter the roughness lengths that are fed to ``surface_flux_mod``.
:mixed_layer_mod: Diagnostics for fluxes of sensible and latent heat are output here. 


References
----------

| [Frierson2006a]_
| [Gierasch2000]_
| [Large2004]_  
| [Manabe1969]_
| [Schneider2009]_
| [VallisEtAl2018]_


Authors
----------
This documentation was written by Ruth Geen, peer reviewed by Marianne Pietschnig, and quality controlled by Ross Castle.

