
Constants: Atmospheric and Planetary Parameters 
===============================================

Summary
-------

``constants.F90`` provides a number of controls that set the fundamental parameters of the planet being studied - it is used by nearly every other module in Isca. By default, Isca models a planet equivalent to earth in terms of physical characteristics such as size, surface gravity and rotational period. In addition, water is set as the default condensate for any moist physics modules that may be used.

Planetary Parameters
^^^^^^^^^^^^^^^^^^^^

Below are a set of controls that allow the physical size and rotation to be set. This allows us to model a range of planets, including others in the Solar System such as Mars and Jupiter, in addition to exoplanets. Rotation can be set by two different parameters: a rotation rate :math:`\Omega` can be specified; or alternatively a period in seconds can be given that is converted back to a value of :math:`\Omega`.

+-------------------+----------------------------------+---------------------------------------+-------------------------------------------------+
| Name              | Default                          | Units                                 | Description                                     |
+===================+==================================+=======================================+=================================================+
|``radius``         | 6371e3                           | m                                     | Radius of planet                                |
+-------------------+----------------------------------+---------------------------------------+-------------------------------------------------+
|``gravity``        | 9.80                             | ms :math:`^{-2}`                      | Surface gravitational acceleration              |
+-------------------+----------------------------------+---------------------------------------+-------------------------------------------------+
|``omega``          | 7.2921150e-5                     | rad :math:`\cdot` s :math:`^{-1}`     | Rotation rate of planet                         |
+-------------------+----------------------------------+---------------------------------------+-------------------------------------------------+
|``rotation_period``| -1.0 (inactive)                  | s                                     | Rotation period of planet (overrides ``omega`` )|
+-------------------+----------------------------------+---------------------------------------+-------------------------------------------------+
|``orbital_period`` | 31557600                         | s                                     | Orbital period of planet                        |
+-------------------+----------------------------------+---------------------------------------+-------------------------------------------------+
|``orbital_rate``   | :math:`2\pi` / ``orbital_period``|  rad :math:`\cdot` s :math:`^{-1}`    | Orbital rate of planet                          |
+-------------------+----------------------------------+---------------------------------------+-------------------------------------------------+

.. note:: Whilst the rotation and orbital rates are set within this module, other parameters associated with planetary motion such as axial tilt (obliquity) and eccentricity are controlled from the ``astronomy_mod`` module.

Dry Atmosphere
^^^^^^^^^^^^^^

Changing the dry atmosphere values allows atmospheres with a wide range of compositions to be studied, including both terrestrial planets (e.g. Earth) and gas giants such as Jupiter. In the case of Earth, the reference surface pressure is taken to be the mean sea level pressure.

The dry air gas constant for any homogeneous atmosphere can be calculated from its chemical composition. It is calculated by dividing the universal gas constant :math:`R` by the average molar mass of the atmosphere.

+------------+----------------------+-------------------------------------+-------------------------------------------------------+
| Name       | Default              | Units                               | Description                                           |
+============+======================+=====================================+=======================================================+
|``pstd_mks``| 101325.0             | Pa / Nm :math:`^{-2}`               | Mean (reference) surface pressure (SI)                |
+------------+----------------------+-------------------------------------+-------------------------------------------------------+
|``pstd``    | 1.013250e06          | dyn :math:`\cdot` cm :math:`^{-2}`  | Mean (reference) surface pressure (cgs)               |
+------------+----------------------+-------------------------------------+-------------------------------------------------------+
|``rdgas``   | 287.04               | Jkg :math:`^{-1}`K :math:`^{-1}`    | Dry air gas constant                                  |
+------------+----------------------+-------------------------------------+-------------------------------------------------------+
|``kappa``   | 2/7                  | dimensionless                       | Heat capacity ratio ( :math:`\gamma` for an ideal gas)|
+------------+----------------------+-------------------------------------+-------------------------------------------------------+
|``cp_air``  | ``rvgas`` / ``kappa``| Jkg :math:`^{-1}`K :math:`^{-1}`    | Dry air heat capcity                                  |
+------------+----------------------+-------------------------------------+-------------------------------------------------------+

.. note:: If the mean surface pressure value is changed here, it is necessary to also set ``reference_sea_level_press`` from the ``spectral_dynamics_nml`` namelist, else the output file will not extend to the pressure specified.


Moist Atmosphere
^^^^^^^^^^^^^^^^

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

Relevant Modules
----------------
Since this module provides the definition of a number of physical constants, it is used by most other modules that exist within the Isca framework.

Authors
-------
This documentation was written by Daniel Williams, peer reviewed by YY, and quality controlled by ZZ.
