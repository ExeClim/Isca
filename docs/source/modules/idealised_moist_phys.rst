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

+----------+-------------------------------+-------------+
| Name     | Purpose                       | Default     |
+----------+-------------------------------+-------------+
| turb     |selects whether the vertical   | False       |
|          |turbulence modules are called. |             |
+----------+-------------------------------+-------------+


turb: selects whether the vertical turbulence modules are called. default: False

Convection
^^^^^^^^^^
Four convection schemes are currently avaiable in Isca. Only one may be selected. The model can also be run with no convection scheme.

The scheme to be used can be selected using the *convection_scheme* namelist parameter, which may be set to:
1. 'SIMPLE_BETTS_CONV'
2. 'FULL_BETTS_MILLER_CONV'
3. 'DRY_CONV'
4. 'RAS_CONV'
5. 'NO_CONV'

Three namelist parameters exist as switches to provide these four options:
1. lwet_convection: If true, use Frierson Quasi-Equilibrium convection scheme (Frierson et al., 2007). default: False
2. do_bm: If true, use the Betts-Miller convection scheme (Betts and Miller,). default: False
3. do_ras: If true, use the relaxed Arakawa Schubert convection scheme (Arakawa and Schubert, ). default: False
4. If the above are all false, the model will use a dry convection scheme.


Radiation
^^^^^^^^^
two_stream_gray: if true, then the two stream radiation module will be called. default: True
do_rrtm_radiation: if true, then the RRTM radiation module will be called. default: False
do_socrates_radiation: if true, then the SOCRATES radiation module will be called. default: False


Drag and Turbulence
^^^^^^^^^^^^^^^^^^^
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


Relevant modules and subroutines
--------------------------------
.. List the names of relevant modules, subroutines, functions, etc.
.. You can add also code snippets, using Sphinx code formatting

In brief, the module calls, in the following order:


References
----------
..
   Add relevant references. This is done in 2 steps:
   1. Add the reference itself to docs/source/references.rst
   2. Insert the citation key here, e.g. [Vallis2017]_
   
   See the Contributing guide for more info.

[SchneiderLiu2009]_
[Frierson2007]_