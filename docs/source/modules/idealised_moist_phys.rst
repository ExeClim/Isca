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
.. Don't forget to add a concise and informative title.

Summary
-------
.. Add a short abstract on what the relevant part of code does.

``idealized_moist_phys.F90`` calls the various modules associated with Isca's moist physics configurations. The specific parameterisations to be used can be selected via namelist input to this module. This is where the radiation, convection, turbulence and land options are set. In addition, timestepping for the bucket hydrology option is managed locally here.

These options allow users to configure a wide range of planets, including notable configurations from the literature (Frierson et al. 2006, Byrne and O'Gorman, Schneider and Liu, Jucker et al. and others). Users should bear in mind that the full parameter space is vast and not all options may be compatible with one another. 

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
|``SIMPLE_BETTS_CONV``     |Use Frierson Quasi-Equilibrium convection scheme (Frierson et al., 2007).  |
+--------------------------+---------------------------------------------------------------------------+
|``FULL_BETTS_MILLER_CONV``|Use the Betts-Miller convection scheme (Betts and Miller,).                |
+--------------------------+---------------------------------------------------------------------------+
|``RAS_CONV``              |Use the relaxed Arakawa Schubert convection scheme ().                     |
+--------------------------+---------------------------------------------------------------------------+
|``DRY_CONV``              |Use the dry convection scheme (CITATION NEEDED?).                          |
+--------------------------+---------------------------------------------------------------------------+
|``NO_CONV``               |Use no convection scheme.                                                  |
+--------------------------+---------------------------------------------------------------------------+
|``UNSET``                 |Model looks through flags listed under Method 2.                           |
+--------------------------+---------------------------------------------------------------------------+

Using the ``convection_scheme`` option will ensure that flags are set consistently for the selected convection scheme, avoiding user error. The default is ``UNSET``.

Method 2 (legacy)
"""""""""""""""""
Three namelist parameters exist as switches to select the convection options:

+-------------------+------------------------------------------------------------+---------+
| Option            | Summary                                                    |Default  |
+===================+============================================================+=========+
|``lwet_convection``|If true, use Frierson Quasi-Equilibrium convection scheme.  |``False``|
+-------------------+------------------------------------------------------------+---------+
|``do_bm``          |If true, use the Betts-Miller convection scheme.            |``False``|
+-------------------+------------------------------------------------------------+---------+
|``do_ras``         |If true, use the relaxed Arakawa Schubert convection scheme.|``False``|
+-------------------+------------------------------------------------------------+---------+

If multiple flags are set as True, or all are False, an error will be raised. This method exists for compatability with configurations pre-dating method 1. Method 1 is preferred for new model configurations.

Humidity calculation options
""""""""""""""""""""""""""""
+-------------+------------------------------------------------------------------+---------+
| Option      | Summary                                                          |Default  |
+=============+==================================================================+=========+
|``do_simple``|If true, use a simplification when calculating relative humidity. |``False``|
+-------------+------------------------------------------------------------------+---------+

Radiation
^^^^^^^^^
Two more comprehensive radiation codes are currently included in Isca: RRTM () and Socrates (). In addition, a number of simple radiation parameterisations for the atmospheres of Earth and other planets can be found in ``two_stream_gray_rad.F90``. The radiation scheme is set through the following flags:

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

Isca will evaluate surface heat exchange provided a gaseous option is not specified (see below). Vertical diffusion may be enabled via ``turb = True``; default is ``False``. This enables calls to ``vert_turb_driver_mod``, ``vert_diff_mod`` and ``mixed_layer_mod``.

Additonal namelist parameters further specify the processes called and parameters used:

+----------------------------+-----------------------------------------------------------------------------+---------+
| Option                     | Summary                                                                     |Default  |
+============================+=============================================================================+=========+
|``mixed_layer_bc``          |Dictates whether ``mixed_layer_mod`` is called.                              |``False``|
+----------------------------+-----------------------------------------------------------------------------+---------+
|``do_virtual``              |Selects whether virtual temperature is used in the vertical diffusion module.|``False``|
+----------------------------+-----------------------------------------------------------------------------+---------+
|``roughness_moist``         |Roughness length for moisture used in surface heat exchange.                 |``0.05`` |
+----------------------------+-----------------------------------------------------------------------------+---------+
|``roughness_mom``           |Roughness length for momentum used in surface heat exchange.                 |``0.05`` |
+----------------------------+-----------------------------------------------------------------------------+---------+
|``roughness_heat``          |Roughness length for heat used in surface heat exchange.                     |``0.05`` |
+----------------------------+-----------------------------------------------------------------------------+---------+
|``land_roughness_prefactor``|Multiplier on the above roughness lengths to allow land-ocean contrast.      | ``1.0`` |
+----------------------------+-----------------------------------------------------------------------------+---------+


**Gaseous Planets**

Isca currently includes options to run a Jupiter-type planet. The relevant namelist parameter for the lower boundary physics in this case is:

+----------------------------+-----------------------------------------------------------------------------+---------+
| Option                     | Summary                                                                     |Default  |
+============================+=============================================================================+=========+
|``gp_surface``              |Turns on **Schneider & Liu (2009)** prescription of lower-boundary heat flux.|``False``|
+----------------------------+-----------------------------------------------------------------------------+---------+


Upper Level Damping
"""""""""""""""""""

At upper levels, damping may be needed to account for subgrid-scale processes that decelerate fast upper-level winds. ``damping_driver_mod`` allows Rayleigh friction to be applied at upper levels, and is currently being configured to allow gravity wave drag options. This is switched on with:

+----------------------------+-----------------------------------------------------------------------------+---------+
| Option                     | Summary                                                                     |Default  |
+============================+=============================================================================+=========+
|``do_damping``              |If true, call ``damping_driver_mod``                                         |``False``|
+----------------------------+-----------------------------------------------------------------------------+---------+



Land and hydrology
^^^^^^^^^^^^^^^^^^

Land and hydrology processes are predominantly dealt with in ``surface_flux_mod`` and ``mixed_layer_mod``, but land and bucket hydrology options are initialised here. We acknowledge that the bucket hydrology is adapted from code by (TS github), and follows (citation). Land and hydrology options in this module are:

+----------------------------+----------------------------------------------------------------------+-------------------+
| Option                     | Summary                                                              |Default            |
+============================+======================================================================+===================+
|``land_option``             |Selects how land-mask is defined, a summary of options is given below.|``False``          |
+----------------------------+----------------------------------------------------------------------+-------------------+
|``land_file_name``          |Filename for the input land-mask.                                     |``'INPUT/land.nc'``|
+----------------------------+----------------------------------------------------------------------+-------------------+
|``land_field_name``         |Field name in the input land-mask netcdf.                             |``'land_mask'``    |
+----------------------------+----------------------------------------------------------------------+-------------------+
|``bucket``                  |If true, use bucket hydrology.                                        |``False``          |
+----------------------------+----------------------------------------------------------------------+-------------------+
|``init_bucket_depth``       |Value at which to initialise bucket water depth over ocean (large).   |``1000.``          |
+----------------------------+----------------------------------------------------------------------+-------------------+
|``init_bucket_depth_land``  |Value at which to initialise bucket water depth over land.            |``20.``            |
+----------------------------+----------------------------------------------------------------------+-------------------+
|``max_bucket_depth_land``   |Maximum depth of water in bucket over land following intialisation.   |``0.15``           |
+----------------------------+----------------------------------------------------------------------+-------------------+
|``robert_bucket``           |Robert coefficient for RAW filter on bucket leapfrog timestepping.    |``0.04``           |
+----------------------------+----------------------------------------------------------------------+-------------------+
|``raw_bucket``              |RAW coefficient for RAW filter on bucket leapfrog timestepping.       |``0.53``           |
+----------------------------+----------------------------------------------------------------------+-------------------+

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
.. What diagnostics are available for this part of the code.

Diagnostics from this module are output under ``mod_name = 'atmosphere'``. Some diagnostics may only be output when certain namelist options are set, e.g. those associated with the bucket hydrology. Requesting unsaved diagnostics in your diagnostic list will result in those diagnostics not being output, but will not cause a fatal error or affect other diagnostics.


+----------------------+-----------------------------------------------------+------------------------------------+
| Name                 | Description                                         | Units                              |
+======================+=====================================================+====================================+
|``dt_ug_diffusion``   | Zonal wind tendency from vertical diffusion         | ms :math:`^{-2}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``dt_vg_diffusion``   | Meridional wind tendency from vertical diffusion    | ms :math:`^{-2}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``dt_tg_diffusion``   | Temperature tendency from vertical diffusion        | Ks :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``dt_qg_diffusion``   | Specific humidity tendency from vertical diffusion  | kg kg :math:`^{-1}` s :math:`^{-1}`|
+----------------------+-----------------------------------------------------+------------------------------------+
|``convection_rain``   | Rain from convection                                | kg m :math:`^{-2}` s :math:`^{-1}` |
+----------------------+-----------------------------------------------------+------------------------------------+
|``condensation_rain`` | Rain from large-scale condensation                  | kg m :math:`^{-2}` s :math:`^{-1}` |
+----------------------+-----------------------------------------------------+------------------------------------+
|``precipitation``     | Precipitation from resolved, parameterised and snow | kg m :math:`^{-2}` s :math:`^{-1}` |
+----------------------+-----------------------------------------------------+------------------------------------+
|``dt_tg_convection``  | Temperature tendency from convection                | Ks :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``dt_qg_convection``  | Specific humidity tendency from convection          | kg kg :math:`^{-1}` s :math:`^{-1}`|
+----------------------+-----------------------------------------------------+------------------------------------+
|``dt_tg_condensation``| Temperature tendency from convection                | Ks :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``dt_qg_condensation``| Specific humidity tendency from convection          | kg kg :math:`^{-1}` s :math:`^{-1}`|
+----------------------+-----------------------------------------------------+------------------------------------+
|``dt_qg_condensation``| Specific humidity tendency from convection          | kg kg :math:`^{-1}` s :math:`^{-1}`|
+----------------------+-----------------------------------------------------+------------------------------------+
|``rh``                | Relative humidity                                   | %                                  |
+----------------------+-----------------------------------------------------+------------------------------------+
|``cape``              | Convective Avaliable Potential Energy               | J kg :math:`^{-1}`                 |
+----------------------+-----------------------------------------------------+------------------------------------+
|``cin``               | Convective Inhibition                               | J kg :math:`^{-1}`                 |
+----------------------+-----------------------------------------------------+------------------------------------+
|``flux_u``            | Surface zonal wind stress                           | N m :math:`^{-2}`                  |
+----------------------+-----------------------------------------------------+------------------------------------+
|``flux_v``            | Surface meridional wind stress                      | N m :math:`^{-2}`                  |
+----------------------+-----------------------------------------------------+------------------------------------+
|``temp_2m``           | Air temperature 2m above surface                    | K                                  |
+----------------------+-----------------------------------------------------+------------------------------------+
|``sphum_2m``          | Specific humidity 2m above surface                  | kg kg :math:`^{-1}`                |
+----------------------+-----------------------------------------------------+------------------------------------+
|``rh_2m``             | Relative humidity 2m above surface                  | %                                  |
+----------------------+-----------------------------------------------------+------------------------------------+
|``u_10m``             | Zonal wind 10m above surface                        | ms :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``v_10m``             | Meridional wind 10m above surface                   | ms :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``v_10m``             | Meridional wind 10m above surface                   | ms :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``diss_heat_ray``     | Heat dissipated by Rayleigh drag in SL09 scheme     | Ks :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``bucket_depth``      | Depth of surface reservoir                          | m                                  |
+----------------------+-----------------------------------------------------+------------------------------------+
|``bucket_depth_conv`` | Tendency of bucket depth due to convection          | ms :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``bucket_depth_cond`` | Tendency of bucket depth due to condensation        | ms :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+
|``bucket_depth_lh``   | Tendency of bucket depth due to evaporation         | ms :math:`^{-1}`                   |
+----------------------+-----------------------------------------------------+------------------------------------+

	 
Relevant modules and subroutines
--------------------------------
.. List the names of relevant modules, subroutines, functions, etc.
.. You can add also code snippets, using Sphinx code formatting

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
..
   Add relevant references. This is done in 2 steps:
   1. Add the reference itself to docs/source/references.rst
   2. Insert the citation key here, e.g. [Vallis2017]_
   
   See the Contributing guide for more info.

[SchneiderLiu2009]_
[Frierson2007]_