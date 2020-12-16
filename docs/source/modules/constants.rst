
Constants: Atmospheric and Planetary Parameters. 
================================================

Documentation for ``constants.f90``. 

Summary
-------

``constants.F90`` provides a number of controls that set the fundamental parameters of the planet being studied in Isca. By default, Isca models a planet equivalent to earth in terms of physical characteristics such as size, surface gravity and rotational period. In addition, water is set as the default condensate for any moist physics modules that may be used.

Planetary Parameters
""""""""""""""""""""


+-------------------+--------------+----------------------------------+----------------------------------------+
| Name              | Default      | Units                            | Description                            |
+===================+==============+==================================+========================================+
|``radius``         | 6371.e3      | m                                | Radius of planet                       |
+-------------------+--------------+----------------------------------+----------------------------------------+
|``gravity``        | 9.80         | ms :math:`^{-1}`                 | Surface gravitational acceleration     |
+-------------------+--------------+----------------------------------+----------------------------------------+
|``omega``          | 7.2921150e-5 | rad :math:`\cdot`s :math:`^{-1}` | Rotation rate of planet                |
+-------------------+--------------+----------------------------------+----------------------------------------+
|``rotation_period``| -1.0         | s                                | Rotation period of planet              |
+-------------------+--------------+----------------------------------+----------------------------------------+
|``orbital_period`` | 31557600     | s                                | Orbital period of planet               |
+-------------------+--------------+----------------------------------+----------------------------------------+
|``orbital_rate``   |              | rad :math:`\cdot`s :math:`^{-1}` | Orbital rate of planet                 |
+-------------------+--------------+----------------------------------+----------------------------------------+


Dry Atmosphere
""""""""""""""

Changing the dry atmosphere values allows atmospheres with a wide range of compositions to be studied, including both terrestrial planets (e.g. Earth) and gas giants such as Jupiter. In the case of Earth, the reference surface pressure is taken to be the mean sea level pressure.

The dry air gas constant for any homogeneous atmosphere can be calculated from its chemical composition. It is calculated by dividing the universal gas constant by the average molar mass of the atmosphere -- taking Jupiter as an example (approximately 90% hydrogen and 10% helium) we obtain

.. math:: 
    R_{\text{dry}} = \frac{R}{\left(0.9 * 1\right) + \left(0.1 * 4\right)}

which using R=8.314 and converting from g to kg gives 6395 :math:`^{-1}`K :math:`^{-1}`.


**Note:** If the mean surface pressure value is changed here, it is necessary to also set ``reference_sea_level_press`` from the ``spectral_dynamics_nml`` namelist, else the output file will not extend to the pressure specified.

+------------+-----------------------+----------------------------------+----------------------------------------+
| Name       | Default               | Units                            | Description                            |
+============+=======================+==================================+========================================+
|``pstd_mks``| 101325.0              | Pa / Nm :math:`^{-2}`            | Mean (reference) surface pressure (SI) |
+------------+-----------------------+----------------------------------+----------------------------------------+
|``pstd``    | 1.013250E+06          | Ba / dyncm :math:`^{-2}`         | Mean (reference) surface pressure (cgs)|
+------------+-----------------------+----------------------------------+----------------------------------------+
|``rdgas``   | 287.04                | Jkg :math:`^{-1}`K :math:`^{-1}` | Dry air gas constant                   |
+------------+-----------------------+----------------------------------+----------------------------------------+
|``kappa``   | 2/7                   | dimensionless                    | Heat capacity ratio                    |
+------------+-----------------------+----------------------------------+----------------------------------------+
|``cp_air``  | ``rvgas`` / ``kappa`` | Jkg :math:`^{-1}`K :math:`^{-1}` | Dry air heat capcity                   |
+------------+-----------------------+----------------------------------+----------------------------------------+



Moist Atmosphere
""""""""""""""""

With the addition of moist physics, a number of additional namelist parameters can be used to change the primary condensate present in the atmosphere. This allows atmospheres in temperature regimes significantly different from Earth to be modelled, including Titan, which has an active methane cycle.

+--------------+-------------------+----------------------------------+-------------------------------------------+
| Name         | Default           | Units                            | Description                               |
+==============+===================+==================================+===========================================+
|``rvgas``     | 461.50            | Jkg :math:`^{-1}`K :math:`^{-1}` | Vapour gas constant                       |
+--------------+-------------------+----------------------------------+-------------------------------------------+
|``cp_vapor``  | 4 * ``rvgas``     | Jkg :math:`^{-1}`K :math:`^{-1}` | Vapour heat capacity                      |
+--------------+-------------------+----------------------------------+-------------------------------------------+
|``dens_vapor``| 1000              | kgm :math:`^{-3}`                | Density of condensate in the liquid phase |
+--------------+-------------------+----------------------------------+-------------------------------------------+
|``hlv``       | 2.500e6           | Jkg :math:`^{-1}`                | Latent heat of vapourisation of condensate|
+--------------+-------------------+----------------------------------+-------------------------------------------+
|``hlf``       | 3.34e5            | Jkg :math:`^{-1}`                | Latent heat of fusion of condensate       |
+--------------+-------------------+----------------------------------+-------------------------------------------+
|``hls``       | ``hlv`` + ``hlf`` | Jkg :math:`^{-1}`                | Latent heat of sublimation of condensate  |
+--------------+-------------------+----------------------------------+-------------------------------------------+
|``tfreeze``   | 273.16            | K                                | Freezing point of condensate              |
+--------------+-------------------+----------------------------------+-------------------------------------------+
|``tppress``   | 610.78            | Pa / Nm :math:`^{-2}`            | Triple point pressure of condensate       |
+--------------+-------------------+----------------------------------+-------------------------------------------+




Relevant modules and subroutines
--------------------------------

TO BE EDITED!!!!!


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


Authors
-------
This documentation was written by Daniel Williams, peer reviewed by YY, and quality controlled by ZZ.

