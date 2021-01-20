..  DO NOT MODIFY THIS FILE UNLESS YOU ARE A CORE MAINTAINER OF ISCA!

..
    This is a reStructuredText template file for creating
    a new documentation entry for the Isca model.
    
    Please make a copy of this file with the appropriate file name and place it
    to the appropriate location within docs/source/ and start writing.
    Once you are done, remove all the comments from your .rst file.
    
    Here is a guide on reST formatting:
    https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

Moist physics driver: idealized_moist_phys.F90
==============================================

Summary
-------

``idealized_moist_phys.F90`` calls the various modules associated with Isca's moist physics configurations. The specific parameterisations to be used can be selected via namelist input to this module. This is where the radiation, convection, turbulence and land options are set. In addition, timestepping for the bucket hydrology option is managed locally here.

These options allow users to configure a wide range of planets, including notable configurations from the literature [Frierson2006a]_, [Byrne2013]_, [Schneider2009]_, [Jucker2017]_. Users should bear in mind that the full parameter space is vast and not all options may be compatible with one another. 

The fortran file is found in ``Isca/src/atmos_spectral/driver/solo/idealized_moist_phys.F90``


Namelist options
----------------
This module controls a large number of switches for different modules. Further information on what each module is for can be found on the module's documentation page. 

NB. Defaults are generally set to ``False`` to avoid accidental use of modules, but this is unlikely to be a useful configuration! For examples of how parameters are set to achieve different model configurations, see the test cases in ``Isca/exp/test_cases/``. 

Humidity, Condensation and Convection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Four convection schemes, and the option of no convection, are currently available in Isca. The defaults are set so the user should actively select a scheme to use via one of the methods below. Large-scale condensation will be called following the convection scheme in all cases, unless ``DRY_CONV`` is selected.

Method 1 (preferred)
""""""""""""""""""""
The scheme to be used can be selected using the ``convection_scheme`` namelist parameter, which may be set to:

+--------------------------+---------------------------------------------------------------------------+
|Value                     |Effect                                                                     |
+==========================+===========================================================================+
|``SIMPLE_BETTS_CONV``     |Use Frierson Quasi-Equilibrium convection scheme [Frierson2007]_.          |
+--------------------------+---------------------------------------------------------------------------+
|``FULL_BETTS_MILLER_CONV``|Use the Betts-Miller convection scheme [Betts1986]_, [BettsMiller1986]_.   |
+--------------------------+---------------------------------------------------------------------------+
|``RAS_CONV``              |Use the relaxed Arakawa Schubert convection scheme [Moorthi1992]_.         |
+--------------------------+---------------------------------------------------------------------------+
|``DRY_CONV``              |Use the dry convection scheme [Schneider2006]_.                            |
+--------------------------+---------------------------------------------------------------------------+
|``NO_CONV``               |Use no convection scheme.                                                  |
+--------------------------+---------------------------------------------------------------------------+
|``UNSET``                 |Model looks through flags listed under Method 2.                           |
+--------------------------+---------------------------------------------------------------------------+

Using the ``convection_scheme`` option will ensure that flags are set consistently for the selected convection scheme, avoiding user error. The default is ``UNSET``.

Method 2 (legacy)
"""""""""""""""""
Three namelist parameters exist as switches to select the convection options:

+-------------------+----------------------------------------------------------------------------------+---------+
| Option            | Summary                                                                          |Default  |
+===================+==================================================================================+=========+
|``lwet_convection``|If true, this is equivalent to specifying ``SIMPLE_BETTS_CONV`` in Method 1.      |``False``|
+-------------------+----------------------------------------------------------------------------------+---------+
|``do_bm``          |If true, this is equivalent to specifying ``FULL_BETTS_MILLER_CONV`` in Method 1. |``False``|
+-------------------+----------------------------------------------------------------------------------+---------+
|``do_ras``         |If true, this is equivalent to specifying ``RAS_CONV`` in Method 1.               |``False``|
+-------------------+----------------------------------------------------------------------------------+---------+

If multiple flags are set as True, or all are False, an error will be raised. This method exists for compatibility with configurations pre-dating method 1. Method 1 is preferred for new model configurations.

Humidity calculation options
""""""""""""""""""""""""""""
+-------------+------------------------------------------------------------------+---------+
| Option      | Summary                                                          |Default  |
+=============+==================================================================+=========+
|``do_simple``|If true, use a simplification when calculating relative humidity. |``False``|
+-------------+------------------------------------------------------------------+---------+

Radiation
^^^^^^^^^
Two more comprehensive radiation codes are currently included in Isca: RRTM [MlawerEtAl1997]_ and Socrates [MannersEtAl2015]_. In addition, a number of simple radiation parameterisations for the atmospheres of Earth and other planets can be found in ``two_stream_gray_rad.F90``. The radiation scheme is set through the following flags:

+-------------------------+-----------------------------------------------+---------+
| Option                  | Summary                                       |Default  |
+=========================+===============================================+=========+
|``two_stream_gray``      |If true, call the two stream radiation module. |``True`` |
+-------------------------+-----------------------------------------------+---------+
|``do_rrtm_radiation``    |If true, call the RRTM radiation module.       |``False``|
+-------------------------+-----------------------------------------------+---------+
|``do_socrates_radiation``|If true, call the SOCRATES radiation module.   |``False``|
+-------------------------+-----------------------------------------------+---------+

Drag and Turbulence
^^^^^^^^^^^^^^^^^^^

Lower Boundary Heat, Momentum & Humidity Exchange
"""""""""""""""""""""""""""""""""""""""""""""""""
Near the surface, different processes will be appropriate for terrestrial vs. gaseous planets, and for different experimental designs. These are determined by namelist parameters in ``idealized_moist_phys`` as follows.

**Terrestrial Planets**

Isca will evaluate surface heat exchange provided a gaseous option is not specified (see below). Vertical diffusion may be enabled by setting the namelist parameter ``turb`` to ``True``; default is ``False``. This enables calls to ``vert_turb_driver_mod``, ``vert_diff_mod`` and ``mixed_layer_mod``.

Additional namelist parameters further specify the processes called and parameters used:

+----------------------------+-----------------------------------------------------------------------------+---------+
| Option                     | Summary                                                                     |Default  |
+============================+=============================================================================+=========+
|``mixed_layer_bc``          |Dictates whether ``mixed_layer_mod`` is called.                              |``False``|
+----------------------------+-----------------------------------------------------------------------------+---------+
|``do_virtual``              |Selects whether virtual temperature is used in the vertical diffusion module.|``False``|
+----------------------------+-----------------------------------------------------------------------------+---------+
|``roughness_moist``         |Roughness length for use in surface moisture exchange.                       |``0.05`` |
+----------------------------+-----------------------------------------------------------------------------+---------+
|``roughness_mom``           |Roughness length for use in surface momentum exchange.                       |``0.05`` |
+----------------------------+-----------------------------------------------------------------------------+---------+
|``roughness_heat``          |Roughness length for use in surface heat exchange.                           |``0.05`` |
+----------------------------+-----------------------------------------------------------------------------+---------+
|``land_roughness_prefactor``|Multiplier on the above roughness lengths to allow land-ocean contrast.      | ``1.0`` |
+----------------------------+-----------------------------------------------------------------------------+---------+


**Gaseous Planets**

Isca currently includes options to run a Jupiter-type planet. The relevant namelist parameter for the lower boundary physics in this case is:

+----------------------------+------------------------------------------------------------------------------+---------+
| Option                     | Summary                                                                      |Default  |
+============================+==============================================================================+=========+
|``gp_surface``              |Turns on prescription of lower-boundary heat flux following [Schneider2009]_ .|``False``|
+----------------------------+------------------------------------------------------------------------------+---------+


Upper Level Damping
"""""""""""""""""""

At upper levels, damping may be needed to account for subgrid-scale processes that decelerate fast upper-level winds. The module ``damping_driver_mod`` applies Rayleigh friction at upper levels. This is switched on with:

+----------------------------+-----------------------------------------------------------------------------+---------+
| Option                     | Summary                                                                     |Default  |
+============================+=============================================================================+=========+
|``do_damping``              |If true, call ``damping_driver_mod``                                         |``False``|
+----------------------------+-----------------------------------------------------------------------------+---------+

NB: ``damping_driver_mod`` is currently being configured to allow gravity wave drag options.

Land and hydrology
^^^^^^^^^^^^^^^^^^

Land and hydrology processes are predominantly dealt with in ``surface_flux_mod`` and ``mixed_layer_mod``, but land and bucket hydrology options are initialised with the following namelist parameters. We acknowledge that the bucket hydrology is adapted from code from https://github.com/tapios and follows [Manabe1969]_. Land and hydrology options in this module are:

+----------------------------+---------------------------------------------------------------------------------+-------------------+
| Option                     | Summary                                                                         |Default            |
+============================+=================================================================================+===================+
|``land_option``             |Selects how land-mask is defined, a summary of options is given below.           |``none``           |
+----------------------------+---------------------------------------------------------------------------------+-------------------+
|``land_file_name``          |Filename for the input land-mask.                                                |``'INPUT/land.nc'``|
+----------------------------+---------------------------------------------------------------------------------+-------------------+
|``land_field_name``         |Field name in the input land-mask netcdf.                                        |``'land_mask'``    |
+----------------------------+---------------------------------------------------------------------------------+-------------------+
|``bucket``                  |If true, use bucket hydrology.                                                   |``False``          |
+----------------------------+---------------------------------------------------------------------------------+-------------------+
|``init_bucket_depth``       |Value at which to initialise bucket water depth over ocean (large, in :math:`m`).|``1000.``          |
+----------------------------+---------------------------------------------------------------------------------+-------------------+
|``init_bucket_depth_land``  |Value at which to initialise bucket water depth over land.                       |``20.``            |
+----------------------------+---------------------------------------------------------------------------------+-------------------+
|``max_bucket_depth_land``   |Maximum depth of water in bucket over land following initialisation.             |``0.15``           |
+----------------------------+---------------------------------------------------------------------------------+-------------------+
|``robert_bucket``           |Robert coefficient for RAW filter* on bucket leapfrog timestepping.              |``0.04``           |
+----------------------------+---------------------------------------------------------------------------------+-------------------+
|``raw_bucket``              |RAW coefficient for RAW filter* on bucket leapfrog timestepping.                 |``0.53``           |
+----------------------------+---------------------------------------------------------------------------------+-------------------+

`*` Roberts-Asselin-Williams filter, [Williams2011]_

``land_option`` may be set to:

+---------------+------------------------------------------------------------------------------------------+
|Value          | Effect                                                                                   |
+===============+==========================================================================================+
|``'input'``    |Read land mask from input file.                                                           |
+---------------+------------------------------------------------------------------------------------------+
|``'zsurf'``    |Define land where surface geopotential height at model initialisation exceeds a threshold.|
+---------------+------------------------------------------------------------------------------------------+
|``'none'``     | Do not apply a land mask                                                                 |
+---------------+------------------------------------------------------------------------------------------+



									  
Diagnostics
-----------

Diagnostics from this module are output under ``mod_name = 'atmosphere'``. Some diagnostics may only be output when certain namelist options are set, e.g. those associated with the bucket hydrology. Requesting unsaved diagnostics in your diagnostic list will result in those diagnostics not being output, but will not cause a fatal error or affect other diagnostics.


+----------------------+-------------------------------------------------------------+------------------------------------+
| Name                 | Description                                                 | Units                              |
+======================+=============================================================+====================================+
|``dt_ug_diffusion``   | Zonal wind tendency from vertical diffusion                 |  :math:`ms^{-2}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``dt_vg_diffusion``   | Meridional wind tendency from vertical diffusion            |  :math:`ms^{-2}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``dt_tg_diffusion``   | Temperature tendency from vertical diffusion                |  :math:`Ks^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``dt_qg_diffusion``   | Specific humidity tendency from vertical diffusion          |  :math:`kg kg^{-1} s^{-1}`         |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``convection_rain``   | Rain from convection                                        |  :math:`kg m^{-2} s^{-1}`          |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``condensation_rain`` | Rain from large-scale condensation                          |  :math:`kg m^{-2} s^{-1}`          |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``precipitation``     | Precipitation from resolved, parameterised and snow         |  :math:`kg m^{-2} s^{-1}`          |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``dt_tg_convection``  | Temperature tendency from convection                        |  :math:`Ks^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``dt_qg_convection``  | Specific humidity tendency from convection                  |  :math:`kg kg^{-1} s^{-1}`         |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``dt_tg_condensation``| Temperature tendency from condensation                      |  :math:`Ks^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``dt_qg_condensation``| Specific humidity tendency from condensation                |  :math:`kg kg^{-1} s^{-1}`         |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``rh``                | Relative humidity                                           | %                                  |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``cape``              | Convective Available Potential Energy                       |  :math:`J kg^{-1}`                 |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``cin``               | Convective Inhibition                                       |  :math:`J kg^{-1}`                 |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``flux_u``            | Surface zonal wind stress                                   |  :math:`N m^{-2}`                  |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``flux_v``            | Surface meridional wind stress                              |  :math:`N m^{-2}`                  |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``temp_2m``           | Air temperature 2m above surface                            | :math:`K`                          |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``sphum_2m``          | Specific humidity 2m above surface                          |  :math:`kg kg^{-1}`                |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``rh_2m``             | Relative humidity 2m above surface                          | %                                  |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``u_10m``             | Zonal wind 10m above surface                                |  :math:`ms^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``v_10m``             | Meridional wind 10m above surface                           |  :math:`ms^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``diss_heat_ray``     | Heat dissipated by Rayleigh drag in [Schneider2009]_ scheme |  :math:`Ks^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``bucket_depth``      | Depth of surface reservoir                                  |  :math:`m`                         |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``bucket_depth_conv`` | Tendency of bucket depth due to convection                  |  :math:`ms^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``bucket_depth_cond`` | Tendency of bucket depth due to condensation                |  :math:`ms^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+
|``bucket_depth_lh``   | Tendency of bucket depth due to evaporation                 |  :math:`ms^{-1}`                   |
+----------------------+-------------------------------------------------------------+------------------------------------+

	 
Relevant modules and subroutines
--------------------------------

Key physics modules managed from this module include:

* ``vert_turb_driver_mod``
* ``vert_diff_mod`` 
* ``two_stream_gray_rad_mod``
* RRTM: see ``Isca/src/atmos_param/rrtm_radiation/``
* SOCRATES: see ``Isca/src/atmos_param/socrates/``
* ``mixed_layer_mod`` 
* ``lscale_cond_mod``
* ``qe_moist_convection_mod`` 
* ``ras_mod``
* ``betts_miller_mod``
* ``dry_convection_mod``
* ``surface_flux_mod``
* ``damping_driver_mod``
* ``rayleigh_bottom_drag_mod``

References
----------

[Betts1986]_
[BettsMiller1986]_
[Byrne2013]_
[Frierson2006a]_ 
[Frierson2007]_
[Jucker2017]_
[Manabe1969]_
[MannersEtAl2015]_
[MlawerEtAl1997]_
[Moorthi1992]_
[Schneider2006]_
[Schneider2009]_
[Williams2011]_

Authors
-------
This documentation was written by Ruth Geen, peer reviewed by Marianne Pietschnig, and quality controlled by Ross Castle.