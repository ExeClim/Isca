..  DO NOT MODIFY THIS FILE UNLESS YOU ARE A CORE MAINTAINER OF ISCA!

..
    This is a reStructuredText template file for creating
    a new documentation entry for the Isca model.
    
    Please make a copy of this file with the appropriate file name and place it
    to the appropriate location within docs/source/ and start writing.
    Once you are done, remove all the comments from your .rst file.
    
    Here is a guide on reST formatting:
    https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

idealised_moist_phys.F90
========================
.. Don't forget to add a concise and informative title.

Summary
-------
.. Add a short abstract on what the relevant part of code does.

idealised_moist_phys.F90 calls the various modules associated with Isca's moist physics configurations. The specific parameterisations to be used can be selected via namelist input to this module. This is where the radiation, convection, turbulence and land options are set. In addition, timestepping for the bucket hydrology option is managed locally here.

These options allow users to configure a wide range of planets, including notable configurations from the literature (Frierson et al. 2006, Byrne and O'Gorman, Schneider and Liu, Jucker et al. and others). However, users should bear in mind that the full parameter space is vast and not all options may be compatible with one another. 


Namelist options
----------------

Convection
^^^^^^^^^^
Four convection schemes, and the option of no convection, are currently avaiable in Isca. The defaults are set so the user should actively select a scheme to use via one of the methods below. 

**Method 1 (preferred)**

The scheme to be used can be selected using the ``convection_scheme`` namelist parameter, which may be set to:

#. ``SIMPLE_BETTS_CONV`` - use Frierson Quasi-Equilibrium convection scheme (Frierson et al., 2007).
#. ``FULL_BETTS_MILLER_CONV`` - use the Betts-Miller convection scheme (Betts and Miller,).
#. ``RAS_CONV`` - use the relaxed Arakawa Schubert convection scheme (Arakawa and Schubert, ).
#. ``DRY_CONV`` - use the dry convection scheme (CITATION NEEDED?).
#. ``NO_CONV`` - use no convection scheme.
#. ``UNSET`` - This is the **default**. If ``convection_scheme`` is unset, the model will look through flags for each scheme, see Method 2. 

Using the ``convection_scheme`` option will ensure that flags are set consistently for the selected convection scheme, minimising user error.

**Method 2 (legacy)**

Three namelist parameters exist as switches to select the convection options:

#. ``lwet_convection`` - If true, use Frierson Quasi-Equilibrium convection scheme (Frierson et al., 2007). Default: False
#. ``do_bm`` - If true, use the Betts-Miller convection scheme (Betts and Miller,). Default: False
#. ``do_ras`` - If true, use the relaxed Arakawa Schubert convection scheme (Arakawa and Schubert, ). Default: False

If multiple flags are set as True, or all are False, an error will be raised. This method exists for compatability with configurations pre-dating method 1. Method 1 is preferred for new model configurations.


Radiation
^^^^^^^^^
Two more comprehensive radiation codes are currently included in Isca: RRTM () and Socrates (). In addition, a number of simple radiation parameterisations for the atmospheres of Earth and other planets can be found in ``two_stream_gray_rad.F90``. The radiation scheme is set through the following flags:

#. ``two_stream_gray`` - if true, then the two stream radiation module will be called. Default: True
#. ``do_rrtm_radiation`` - if true, then the RRTM radiation module will be called. Default: False
#. ``do_socrates_radiation`` - if true, then the SOCRATES radiation module will be called. Default: False


Drag and Turbulence
^^^^^^^^^^^^^^^^^^^
turb: selects whether the vertical turbulence modules are called. default: False
do_damping: default: false
mixed_layer_bc: default: false
gp_surface: If true, use the Schneider & Liu (2009) prescription of lower-boundary heat flux. default: false
do_simple: default: false
roughness_moist: default: 0.05
roughness_mom: default: 0.05
roughness_heat: default: 0.05
land_roughness_prefactor: default: 1.0
do_virtual: selects whether virtual temperature is used in gcm_vert_diff. default: False


Land and hydrology
^^^^^^^^^^^^^^^^^^
land_option: default: 'none'
land_file_name: default: 'INPUT/land.nc'
land_field_name: default: 'land_mask'
bucket: default: False
init_bucket_depth: default: 1000.
init_bucket_depth_land: default: 20.
max_bucket_depth_land: default: 0.15
robert_bucket: default: 0.04
raw_bucket: default: 0.53

									  
Diagnostics
-----------
.. What diagnostics are available for this part of the code.

Diagnostics from this module are output under ``mod_name = 'atmosphere'``.


+----------------------+-----------------------------------------------------+------------------------------------+
| Name                 | Description                                         | Units                              |
+----------------------+-----------------------------------------------------+------------------------------------+
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
|``diss_heat_ray``     | Heat dissipated by Rayleigh drag                    | Ks :math:`^{-1}`                   |
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


References
----------
..
   Add relevant references. This is done in 2 steps:
   1. Add the reference itself to docs/source/references.rst
   2. Insert the citation key here, e.g. [Vallis2017]_
   
   See the Contributing guide for more info.

[SchneiderLiu2009]_
[Frierson2007]_