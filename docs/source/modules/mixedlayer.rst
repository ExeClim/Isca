Mixed layer module
=====================

Summary
----------------------
This module updates the sea surface temperature (SST) noted as :math:`T_s` below. 

SST boundary condition options
-----------------------

The SST options are:
    - prescribe SST from an input file (``do_sc_sst`` namelist option).
    - prescribe SST to follow an AquaPlanet Experiment protocol (APE) analytic form (``do_ape_sst`` namelist option).
    - calculate SST based on the surface fluxes and mixed layer depth of a **slab ocean**, with the option of including a Q-flux (either analytic or read from a file).

Note that only the final case will generate a closed surface energy budget. Each of these options are discussed below, followed by an outline of the implicit time-stepping process used when SST is not prescribed. 

Input SST file
-----------------------
Set ``do_sc_sst`` to True if you want to specify the SST field. The SST field will be read in from a NetCDF file with a file name specified by the ``sst_file`` variable. 
Using an input SST field is useful, for example,  when you want to add a temperature anomaly. The file is read in during the initialisation (i.e. within the call to ``mixed_layer_init`` from within ``idealized_moist_phys``).
The input file can be time independent (i.e. no diurnal or seasonal cycle or any changes in the SST from one time step to another) or vary with time. More information can be found in the Diagnostics section below.

APE aquaplanet (analytic SST)
-----------------------
The prescribed SST for the APE aquaplanet protocol is given by:

.. math::
    T_s = 27 \left( 1 - \sin^2\left( \frac{3}{2} \lambda \right) \right),

between 60N-60S, equation 1 of Neele and Hoskins 2004 [NealeHoskins2004]_, and 0 deg C poleward of 60N/S.

Slab ocean 
-----------------------
To allow the SST to evolve based on the surface fluxes, the atmosphere can be coupled to a slab ocean, whose depth is specified by the namelist parameter ``depth``. This also allows a closed surface energy budget, useful for e.g. simulations with increased greenhouse gases.

In this case, during the initialisation, if there is no restart file to open, the surface temperature is set to the prescribed initial distribution (``prescribe_initial_dist = True``):

.. math::
    T_s = T_{surf} -\frac{1}{3} dT \left(3\sin(\lambda)^2-1\right),

where the default values we use in the trip test for the Frierson test case (``exp/test_cases/frierson/frierson_test_case.py``) are: :math:`T_{surf} = 285 K` and :math:`dT = 40 K`.

This form of :math:`T_{surf}` is similar to a 2nd legendre polynomial, it is a parabola that maximises at the equator.

Implicit timestepping procedure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The mixed layer module calculates the evolution of surface temperature using an implicit timestep.
Whereas an explicit method uses the current state of the system to calculate the state of the system 
at the next timestep, an implicit method uses the inferred state of the system at the next timestep.

The net flux into the surface is given by::

    net_surf_sw_down + surf_lw_down - flux_r - SH - LH + ocean_qflux

where ``net_surf_sw_down`` is the net shortwave radiation into the surface, ``surf_lw_down`` the downwelling longwave radiation into the surface, and ``flux_r`` the upwelling longwave radiation at the
surface. The sensible heat (``SH``) is the sum of the surface sensible heat flux and the temperature diffusion flux. The latent heat (``LH``) is the sum of the surface latent heat and the 
moisture diffusion flux (only included if ``evaporation`` is true). The ``ocean_flux`` is an optional term and is zero if the Q-flux is not used.

If ``do_calc_eff_heat_cap`` is True, the ``land_sea_heat_capacity`` is used at each timestep to compute the surface temperature using the following equation::

    land_sea_heat_capacity * dTs/dt = - corrected_flux - t_surf_dependence * dTs

The change in surface temperature in this timestep (``dTs/dt``) times the heat capacity of the surface is equal to the net flux into the surface at the given timestep (``- corrected_flux``) plus the change in the net flux 
into the surface caused by the change in surface temperature (``- t_surf_dependence * dTs``).

This is simplified by defining ``eff_heat_capacity`` as ``land_sea_heat_capacity + t_surf_dependence * dt``, and then ``Ts`` is updated using::

    eff_heat_capacity * dTs/dt = - corrected_flux

Optional Q-flux
^^^^^^^^^^^^^^^
The slab ocean model only communicates between grid-boxes in the vertical (i.e. air-sea exchange) but does not represent any horizontal transport (i.e. no north-south or east-west communications between grid cells). 
An idealised horizontal transport can be included using an ocean heat flux (Q-flux). Atmospheric heat transport is more realistic with an ocean heat transport.

Isca is able to calculate an analytic Q-flux, appropriate to an aquaplanet, following Merlis et al 2013 [MerlisEtAl2013] if ``do_qflux`` is True. Additionally, an analytic warmpool can be added  by setting ``do_warmpool``.  The warmpool structure is set in ``Isca/src/atmos_param/qflux/qflux.f90``. Alternatively an arbitrary Q-flux may be read from a file if ``load_qflux`` is True. 

In the case where both a specific SST distribution (e.g. AMIP climatology) and closed surface energy budget are desired, a Q-flux input file can be generated by running a control experiment with the prescribed SST, creating a NetCDF Q-flux file offline and then passing this file to the model via the python interface run script.

1. Run a prescribed experiment (i.e. a control) using either observations, AMIP or similar.

2. Using the prescribed SST field and the surface fluxes from step 1, create the Q-flux file. This is an offline script that is run independently of the model. An example script is shown in: ``src/extra/python/scripts/calculate_qflux/calculate_qflux.py`` but you can create your own script to do this depending on your application.

3. Add the Q-Flux file to the ``inputfiles`` in the python run script (same as you would for ozone, land etc). Then in the ``mixed_layer_nml`` namelist in the python run script set ``load_qflux`` to True, ``qflux_file_name`` to the name of the input file (don't include the .nc extension) and ``qflux_field_name`` is the Q-flux variable name in the file.


See Q-flux options below for namelist options. Note that the Q-flux is only relevant for slab ocean experiments (not fixed or prescribed SST runs). Also note that if the MiMA radiation code is used then the Q-flux is implemented following Merlis et al 2013 [MerlisEtAl2013]_

More information on the method for Q-flux can be found in Russel et al 1985 [RusselEtAl1985]_


Namelist options
----------------

+-------------------+------------------------------------------------------------+---------+
| Option            | Summary                                                    |Default  |
+===================+============================================================+=========+
|``evaporation``    |Switch for surface evaporation.                             |``True`` |
+-------------------+------------------------------------------------------------+---------+
|``depth``          |Mixed layer depth.                                          | ``40.0``|
+-------------------+------------------------------------------------------------+---------+

Q-flux options
^^^^^^^^^^^^^^^^^^^^
If ``do_qflux`` is True, use ``qflux_amp`` and ``qflux_width`` to calculate a time-independent surface Q-flux.

+-------------------+----------------------------------------------------------------+---------+
| Option            | Summary                                                        |Default  |
+===================+================================================================+=========+
|``do_qflux``       | Switch to calculate time-independent Q-flux.                   |``False``|
+-------------------+----------------------------------------------------------------+---------+
|``qflux_amp``      | Amplitude of time-independent Q-flux if ``do_qflux`` is True.  | ``0.0`` |
+-------------------+----------------------------------------------------------------+---------+
|``qflux_width``    | Width of time-independent Q-flux if ``do_qflux`` is True.      | ``16.0``|
+-------------------+----------------------------------------------------------------+---------+

If ``load_qflux`` is True, use input file to load in a time-independent or time-dependent Q-flux.

+----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+
| Option               | Summary                                                                                                                                                                     |Default          |
+======================+=============================================================================================================================================================================+=================+
|``load_qflux``        | Switch to use input file to get Q-flux.                                                                                                                                     | ``False``       |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+
|``qflux_file_name``   | Name of file among input files, from which to get Q-flux.                                                                                                                   | ``ocean_qflux`` |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+
|``qflux_field_name``  | Name of field name in Q-flux file name, from which to get Q-flux. This is only used when ``time_varying_qflux`` is False. Otherwise the code assumes field_name = file_name.| ``ocean_qflux`` |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+
|``time_varying_qflux``| Flag that determines whether input Q-flux file is time dependent.                                                                                                           | ``False``       |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+

Initialize surface temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+-------------------------------+----------------------------------------------------------------------------------+-----------+
| Option                        | Summary                                                                          |Default    |
+===============================+==================================================================================+===========+
|``prescribe_initial_dist``     | Switch to turn on setting the initial surface temperature distribution.          | ``305.0`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``tconst``                     | Initial surface temperature following formula in ``Slab ocean`` section.         | ``305.0`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``delta_T``                    | Initial surface temperature gradient following formula in ``Slab ocean`` section.| ``40.0``  |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``do_read_sst``                | Flag to use fixed SSTs, prescribed from input file (``sst_file``).               | ``False`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``sst_file``                   | Name of file containing fixed SSTs.                                              | ``False`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``do_sc_sst``                  | Flag to use fixed SSTs, prescribed from input file (``sst_file``).               | ``False`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``specify_sst_over_ocean_only``| Flag to specify SSTs only over ocean, only works if ``do_sc_sst`` is True.       | ``False`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``do_ape_sst``                 | Flag to set prescribed SST according to the APE aquaplanet analytic form         | ``False`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``add_latent_heat_flux_anom``  | Flag to add an anomalous latent heat flux                                        | ``False`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+
|``do_warmpool``                | Flag to call warmpool module, which returns ``ocean_qflux``.                     | ``False`` |
+-------------------------------+----------------------------------------------------------------------------------+-----------+

Surface albedo options
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are 5 options for setting the surface albedo, determined by the value of ``albedo_choice``.
    - 1: Surface albedo is a constant (``albedo_value``). 
    - 2: Glacier with higher albedo in one hemisphere only. If ``lat_glacier`` >0, albedo = ``higher_albedo`` North of ``lat_glacier``. If ``lat_glacier`` <0, albedo = ``higher_albedo`` South of ``lat_glacier``.
    - 3: Glacier with higher albedo in both hemispheres. Albedo = ``higher_albedo`` where latitude > ``|lat_glacier|``.
    - 4: Albedo set by ``albedo_value + (higher_albedo - albedo_value) (lat/90) ^ albedo_exp``.
    - 5: Tanh increase around ``albedo_cntr`` with ``albedo_wdth``::
    
        albedo(lat) = albedo_value + (higher_albedo-albedo_value)* 0.5 *(1+tanh((lat-albedo_cntr)/albedo_wdth)).

+-------------------+-----------------------------------------------------------------------------+---------+
| Option            | Summary                                                                     |Default  |
+===================+=============================================================================+=========+
|``albedo_choice``  | Switch to choose surface albedo option described above.                     | ``1``   |
+-------------------+-----------------------------------------------------------------------------+---------+
|``albedo_value``   | Parameter that sets surface albedo depending on albedo choice.              | ``0.06``|
+-------------------+-----------------------------------------------------------------------------+---------+
|``higher_albedo``  | Parameter that sets surface albedo depending on albedo choice.              | ``0.10``|
+-------------------+-----------------------------------------------------------------------------+---------+
|``lat_glacier``    | Parameter that sets the glacier latitude for albedo choices 2 and 3.        | ``60.0``|
+-------------------+-----------------------------------------------------------------------------+---------+
|``albedo_exp``     | Parameter that sets surface albedo latitude dependence for albedo choice 4. | ``2.``  |
+-------------------+-----------------------------------------------------------------------------+---------+
|``albedo_cntr``    | Parameter that sets surface albedo for albedo choice 5.                     | ``45.0``|
+-------------------+-----------------------------------------------------------------------------+---------+
|``albedo_wdth``    | Parameter that sets surface albedo for albedo choice 5.                     | ``10``  |
+-------------------+-----------------------------------------------------------------------------+---------+

Land options
^^^^^^^^^^^^^^^^

There are 4 options for setting up the land, determined by the value of ``land_option``.
    - ``none``: No land.
    - ``input``: Use input file to determine land mask.
    - ``zsurf``: The surface heat capacity is set to ``land_capacity`` where the surface geopotential is greater than 10.
    - ``lonlat``: The surface heat capacity is set to ``land_capacity`` in the longitude / latitude boxes set by [slandlon(k), elandlon(k)] and [slandlat(k), elandlat(k)] for all k's.

+------------------------------+---------------------------------------------------------------------------------------------------------+----------+
| Option                       | Summary                                                                                                 | Default  |
+==============================+=========================================================================================================+==========+
|``land_option``               | Switch to choose land option as described above.                                                        | ``none`` |
+------------------------------+---------------------------------------------------------------------------------------------------------+----------+
|``land_depth``                | Value of land mixed layer depth.                                                                        | ``-1``   |
+------------------------------+---------------------------------------------------------------------------------------------------------+----------+
|``slandlon``                  | Vector determining lower bounds of longitudes for land masses.                                          | ``0``    |
+------------------------------+---------------------------------------------------------------------------------------------------------+----------+
|``slandlat``                  | Vector determining lower bounds of latitudes for land masses.                                           | ``0``    |
+------------------------------+---------------------------------------------------------------------------------------------------------+----------+
|``elandlon``                  | Vector determining higher bounds of longitudes for land masses.                                         | ``-1``   |
+------------------------------+---------------------------------------------------------------------------------------------------------+----------+
|``elandlat``                  | Vector determining higher bounds of latitudes for land masses.                                          | ``-1``   |
+------------------------------+---------------------------------------------------------------------------------------------------------+----------+
|``land_h_capacity_prefactor`` | Factor by which to multiply ocean heat capacity to get land heat capacity if ``input`` option is used.  | ``1.0``  |
+------------------------------+---------------------------------------------------------------------------------------------------------+----------+
|``land_albedo_prefactor``     | Factor by which to multiply ocean albedo to get land albedo if ``input`` option is used.                | ``1.0``  |
+------------------------------+---------------------------------------------------------------------------------------------------------+----------+

Ice options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
+-------------------------------+-------------------------------------------------------------------------------------------------------------+----------------------+
| Option                        | Summary                                                                                                     |Default               |
+===============================+=============================================================================================================+======================+
|``update_albedo_from_ice``     | Flag to set the surface albedo to ``ice_albedo_value`` where there is ice as specified by ``ice_file_name`` | ``False``            |
+-------------------------------+-------------------------------------------------------------------------------------------------------------+----------------------+
|``ice_albedo_value``           | Value for ice albedo when ``update_albedo_from_ice`` is True.                                               | ``0.7``              |
+-------------------------------+-------------------------------------------------------------------------------------------------------------+----------------------+
|``ice_file_name``              | Name of file containing sea ice concentration.                                                              | ``siconc_clim_amip`` |
+-------------------------------+-------------------------------------------------------------------------------------------------------------+----------------------+
|``ice_concentration_threshold``| Value of sea ice concentration above which albedo should be set to ``ice_albedo_value``.                    | ``0.5``              |
+-------------------------------+-------------------------------------------------------------------------------------------------------------+----------------------+
|``ice_file_name``              | Name of file containing sea ice concentration.                                                              | ``siconc_clim_amip`` |
+-------------------------------+-------------------------------------------------------------------------------------------------------------+----------------------+

Diagnostics
-------------------
+---------------------------+-------------------------------------+----------------------------------------------+
| Name                      | Description                         | Units                                        |
+===========================+=====================================+==============================================+
| ``t_surf``                | Surface temperature                 | K                                            |
+---------------------------+-------------------------------------+----------------------------------------------+
| ``delta_t_surf``          | Surface temperature change          | K                                            |
+---------------------------+-------------------------------------+----------------------------------------------+
| ``flux_t``                | Surface sensible heat flux          | :math:`\text{W}\,\text{m}^{-2}`              |
+---------------------------+-------------------------------------+----------------------------------------------+
| ``flux_lhe``              | Surface latent heat flux            | :math:`\text{W}\,\text{m}^{-2}`              |
+---------------------------+-------------------------------------+----------------------------------------------+
| ``flux_oceanq``           | Ocean heat flux                     | :math:`\text{W}\,\text{m}^{-2}`              |
+---------------------------+-------------------------------------+----------------------------------------------+
| ``ice_conc``              | Sea ice concentration               | 0-1                                          |
+---------------------------+-------------------------------------+----------------------------------------------+
| ``albedo``                | Surface albedo                      | 0-1                                          |
+---------------------------+-------------------------------------+----------------------------------------------+
| ``land_sea_heat_capacity``| Heat capacity of land and sea       | :math:`\text{J}\,\text{m}^{-2},\text{K}^{-1}`|
+---------------------------+-------------------------------------+----------------------------------------------+


Relevant modules and subroutines
--------------------------------
.. List the names of relevant modules, subroutines, functions, etc.
.. You can add also code snippets, using Sphinx code formatting

The mixed layer code is located in: ``src/atmos_spectral/driver/solo/mixed_layer.F90``. The name of this file reflects the fact that the code determines the properties of the single layer (either a slab ocean model 
or prescribed SST) below the air-sea interface.

The mixed layer ocean is initialised and called by:  ``src/atmos_spectral/driver/solo/idealized_moist_phys.F90``.

Relevant routines which are called by the mixed layer ocean:
    - The SST input file is read in using the interpolator module found here: ``src/atmos_shared/interpolator/interpolator.F90``.
    - The Q-flux and warmpool components use the Q-flux module: ``src/atmos_param/qflux/qflux.f90``.

References
----------
..
   Add relevant references. This is done in 2 steps:
   1. Add the reference itself to docs/source/references.rst
   2. Insert the citation key here, e.g. [Vallis2017]_
   
   See the Contributing guide for more info.
[Vallis2017]_
[NealeHoskins2004]_
[MerlisEtAl2013]_
[RusselEtAl1985]_

Authors
----------
..
This documentation was written by Matthew Henry and Penelope Maher, peer reviewed by Stephen Thomson, and quality controlled by Ruth Geen.
