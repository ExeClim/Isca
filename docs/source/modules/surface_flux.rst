
Surface flux calculation
=======================================================================================

Documentation for ``surface_flux_mod`` and ``gp_surface_flux``, both found at ``Isca/src/coupler/surface_flux.F90``. 


Summary
-------

This module manages the exchange of heat, momentum and moisture at the planet's surface. 
A brief overview of how the exchange of these is calculated and the options available is given first, with a full summary of the Namelist parameters below.

Heat
----
As default, fluxes are calculated as...

Alternatively, the namelist option ``flux_heat_gp`` allows the user to prescribe an intrinsic heat flux at 1 bar, appropriate for a giant planet atmosphere where there is no surface [Schneider2009]. 

Momentum
--------
As default, fluxes are calculated as....

Humidity
--------
As default, fluxes are calculated as...


flux_q    =  rho_drag * land_evap_prefactor * (land_humidity_prefactor*q_surf0 - q_atm)


Namelist options
----------------

The namelist options for **surface_flux_nml** are listed below. These include options to limit evaporation over land, to use more simplified parametrisations, to account for salinity over ocean, and to help with issues resulting from truncation. We note that Isca is developed from the GFDL FMS model (https://github.com/NOAA-GFDL/FMS) and some options remain to enable old model configurations to be supported, as indicated in the table. 

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

| [Gierasch2000]_
| [Large2004]_  
| [Schneider2009]_
| [VallisEtAl2018]_


Authors
----------
This documentation was written by Ruth Geen, peer reviewed by , and quality controlled by .

