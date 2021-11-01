Topography and Land Masks in Isca
==================

This guide covers:

1. How to implement a land mask and/or topography (using the T42 land mask provided).
2. How to create/modify topography to suit your own needs (using the python tools included in Isca).
3. The ``topography.F90`` module, which used to control land masks/topography but has now been superceeded. It still retains some useful functions however.

1 - Implementing Topography in Isca
-----------------------------------

This is relatively simple to do, assuming that the topography/land mask you are using has already been interpolated to the same spectral grid you are running the model at. If running at T42, you can use the one provided (``input/land_masks/era_land_t42.nc``). If not, you can interpolate fairly easily from any high resolution data, or for configurations of Earth continents, please refer to section 2 of this guide.

- The land mask is a binary field the size of the spectral model grid (e.g. 64 x 128 for T42) that contains a ``1`` if that grid box is designated as being land or a ``0`` if it is not. 
- Topography is a separate field, again the size of the spectral model grid, containing the surface height for that grid box in meters.
- Often a 'land mask' will contain both these fields.

**Can I have land without topography?** Yes - For idealized modelling sometimes it is useful to have the land/sea contrasts (different albedo, no ocean exchange and different surface parameters) without the complexity of topography. It is also possible to have topography without land (aquamountains). 

Implementing these options require different namelist values. There are two namelists that need to be told about land and topograhy; ``idealized_moist_phys_nml`` (land mask) and ``spectral_init_cond_nml`` (topography). 

idealized_moist_phys
^^^^^^^^^^^^^^^^^^^^
NOTE: see the :ref:`idealized moist phys <idealized_moist_phys>` documentation for additional guidance on this namelist.

The ``idealized_moist_phys`` module is primarily concerned with the land mask. Its job is to pass the information about the land mask to other modules, which in the first instance are the ``mixed_layer`` and the ``surface_flux`` modules, which deal with ocean and surface exchanges with the atmosphere respectively. Both modules need to know where there is ocean/land so correct physics can be computed in the correct location. In addition to passing this information, ``idealized_moist_phys`` will also use the land mask for some other model features, for example some land based parameters and calculating buckets if using bucket hydrology.

The ``idealized_moist_phys`` namelist options to include the appropriate land mask are:

| ``land_option : 'input'`` (the default is ``'none'``, - this is essentially the on switch)
| ``land_file_name : 'INPUT/era_land_t42.nc'`` (using the example of the provided land file)
| ``land_field_name : 'land_mask'`` (this is the default and doesn't need to be specified if that is correct)

spectral_init_cond
^^^^^^^^^^^^^^^^^^
Now there is a land mask, topography can be included if desired. If not, this step can be ignored.

The namelist options below are to include topography from the same file as the land mask. It can be a different file, as long as the grid size is still correct and the height units is meters. There are multiple other options for including topography, some of which are unused at the moment.

| ``topography_option : 'input'``  (Tell model to get topography from input file)
| ``topog_file_name : 'era_land_t42.nc'`` (Again, we use the provided file as the example)
| ``topog_field_name : 'zsurf'`` (The height field name)
| ``land_field_name : 'land_mask'`` (The land field name)

NOTE: ``zsurf`` and ``land_mask`` are the default values in Isca and these are the names used in the provided land mask so does not need to be included when using it.

This set up is the standard way to use topography in Isca, using the ``spectral_init_cond`` module. The module deals with topography, setting up the shape of the model boundary layer.

Other ``topography_option`` options are available, but not widely used by the Isca team:

- No topography - the default.
- Flat topography - the surface geopotential is 0.
- Interpolated topography - where it will call the ``topography`` module below, however this set up is not currently used.
- Gaussian Topography - Simple Gaussian-shaped mountains are generated from specified parameters.

2 - Creating Custom Topography
------------------------------

We provide python code that allows land files to be generated for a range of scenarios. This code is found in ``src/extra/python/isca/land_generator_fn.py``. Some tweaking of the code will be needed to suit the users requirements.

Land Options
^^^^^^^^^^^^
- 'square' (default) Square block of land with boundaries specified by boundaries keyword, a list of 4 floats in the form [S,N,W,E]
- 'continents_old' Choose continents from the original continent set-up adapted from the Sauliere 2012 paper (Jan 16), including North and South America, Eurasia, and Africa. 
- 'continents' Choose continents from a newer continet set-up allowing addition of India, Australia, and South East Asia.
- If continents keyword is set to 'all' (default), then this will include all possibilities for the given set-up. Alternatively, continents can be set to a list of strings allowing any combination of continents. Instructions for this are available in the code. 

Topography Options:
^^^^^^^^^^^^^^^^^^^
- 'none' (default) Topography set to zero everywhere
- 'sauliere2012' Choose mountains from [Sauliere2012]_ configuration using mountains keyword. Default is 'all', alternatively only 'rockys' or 'tibet' may be specified
- 'gaussian' Use parameters specified in topo_gauss keyword to set up a Gaussian mountain. topo_gauss should be a list in the form: [central_lat, central_lon, radius_degrees, std_dev, height]

3 - Topography Module (topography.F90)
--------------------------------------

Summary
^^^^^^^
The ``topography`` module contains numerous routines for creating land surface topography fields and land-water masks for specified latitude-longitude grids. It does this by interpolating from a high resolution netCDF file, which is designated in the namelist. The module was originally written to work with the 1/6 degree Navy mean topography and water data sets. However, any netCDF file can be used as an input in the namelist, providing that the file contains grid box boundaries, (which should be named in the namelist) and whether degrees or radians is specified in the namelist.

As mentioned above, this module is generally not called anymore, in the normal Isca set up, land masks are actually specified through the ``idealized_moist_phys`` module, and the model topography through the ``spectral_init_cond`` module. The main use of it is to provide the "subgrid topography" when using the orographic gravity wave drag scheme (``mg_drag``).

The fields that can be generated with this module are mean and standard deviation of topography within the specified grid boxes; and land-ocean (or water) mask and fractional area. The interpolation conserves the area weighted average of the input data by using the ``horiz_interp`` module.

Namelist options
^^^^^^^^^^^^^^^^
| ``topog_file`` - The topography file that you wish to use.
| ``water_file`` - The water data file (not commonly used, not provided in the Isca release)

For a typical Earth set up the namelist would simply be:

``'topog_file': navy_topography.nc``

This essentially just points the ``topography`` module to this file when it is asked for by another subroutine, e.g. ``mg_drag``.

Diagnostics
^^^^^^^^^^^
There are no diagnostics available directly through this module. The subgrid topography variance can be obtained through damping driver by asking for ``sgsmtn`` (``mg_drag`` must be turned on).

Relevant subroutines
^^^^^^^^^^^^^^^^^^^^

NOTE: Some subroutines have dimensional variants, e.g. interp_topog has both a 1d and 2d variant. 

**Public Subroutines**

These subroutines are used by the topography module to produce the land-masks etc that are being asked for by the user. They are largely self explanatory.

| ``get_topog_mean`` returns the mean height from a region of the topography file so that that value can be used as the value when interpolating onto a smaller grid.
| ``get_topog_stdev`` returns the standard deviation from a region of the topography file so that that value can be assoiated with the same region when interpolating onto a smaller grid.
| ``get_ocean_frac`` returns the fraction of the land mask that is covered by ocean.
| ``get_ocean_mask`` returns an ocean/land mask
| ``get_water_frac`` returns the fraction of the land mask that is covered by water.
| ``get_water_mask`` returns a water/land mask
| ``gaussian_topog_init`` and ``get_gaussian_topog`` call the gaussian topography module.

**Private Subroutines**

There are other subroutines called by the above. These are listed below:
``open_topog_file``, ``interp_topog``, ``find_indices``, ``input_data``, ``interp_water``, ``determine_ocean_points``, ``read_namelist``.

References
----------
See [Sauliere2012]_ for the topography option in ``land_generator_fn.py``.
   
Authors
-------
This documentation was written by Ross Castle, peer reviewed by Ruth Geen, and quality controlled by Marianne Pietschnig.
