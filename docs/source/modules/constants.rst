
Constants: Atmospheric and Planetary Parameters 
===============================================

Summary
-------

``constants.F90`` provides a number of controls that set the fundamental parameters of the planet being studied - it is used by nearly every other module in Isca. By default, Isca models a planet equivalent to Earth in terms of physical characteristics such as size, surface gravity and rotational period. In addition, water is set as the default condensate for any moist physics modules that may be used.

Planetary Parameters
^^^^^^^^^^^^^^^^^^^^

Below are a set of controls that allow the physical size and rotation of a planet to be set in the namelist. This allows us to model a range of planets, including others in the Solar System such as Mars and Jupiter, in addition to exoplanets. Rotation can be set by two different parameters: a rotation rate ``omega`` can be specified; or alternatively a period in seconds can be given that is converted back to a value of ``omega``. 

+---------------------+----------------------------------+----------------------------------+-----------------------------------------------------------+
| Name                | Default                          | Units                            | Description                                               |
+=====================+==================================+==================================+===========================================================+
|``radius``           | :math:`6371\times 10^3`          | m                                | Radius of planet                                          |
+---------------------+----------------------------------+----------------------------------+-----------------------------------------------------------+
|``gravity``          | :math:`9.80`                     | ms :math:`^{-2}`                 | Surface gravitational acceleration                        |
+---------------------+----------------------------------+----------------------------------+-----------------------------------------------------------+
|``omega``            | :math:`7.2921150\times 10^{-5}`  | rad :math:`\cdot` s :math:`^{-1}`| Rotation rate of planet                                   |
+---------------------+----------------------------------+----------------------------------+-----------------------------------------------------------+
|``orbital_period``   | :math:`31557600`                 | s                                | Orbital period of planet                                  |
+---------------------+----------------------------------+----------------------------------+-----------------------------------------------------------+
|``solar_constant``   | :math:`1368.22`                  | Wm :math:`^{-2}`                 | Stellar irradiance                                        |
+---------------------+----------------------------------+----------------------------------+-----------------------------------------------------------+
|``earthday_multiple``| False                            | n/a                              | Modifies seconds per sol calculation (planetary solar day)|
+---------------------+----------------------------------+----------------------------------+-----------------------------------------------------------+

By default the parameter ``seconds_per_sol`` (where ``earthday_multiple`` is false) is calculated using

.. math:: \texttt{seconds_per_sol} = \left|\frac{2\pi}{\texttt{orbital_rate} - \texttt{omega}}\right| ,

whereas setting ``earthday_multiple`` to true modifies the calculation to

.. math:: \texttt{seconds_per_sol} = 86400 \cdot \frac{\texttt{earth_omega}}{\texttt{omega}} .

.. note:: Whilst the rotation and orbital rates are set within this module, other parameters associated with planetary motion such as axial tilt (obliquity) and eccentricity are controlled from the ``astronomy_mod`` module. 

Atmospheric Parameters
^^^^^^^^^^^^^^^^^^^^^^

Changing the basic atmospheric properties allows atmospheres with a wide range of compositions to be studied, including both terrestrial planets (e.g. Earth/Mars) and gas giants such as Jupiter. In the case of Earth, the reference surface pressure is taken to be the mean sea level pressure.

The dry air gas constant for any homogeneous atmosphere can be calculated from its chemical composition. It is calculated by dividing the universal gas constant :math:`R` by the average molar mass of the atmosphere.

+------------+----------------------------+-------------------------------------+-------------------------------------------------------+
| Name       | Default                    | Units                               | Description                                           |
+============+============================+=====================================+=======================================================+
|``pstd_mks``| :math:`101325.0`           | Pa / Nm :math:`^{-2}`               | Mean (reference) surface pressure (SI)                |
+------------+----------------------------+-------------------------------------+-------------------------------------------------------+
|``pstd``    | :math:`1.013250\times 10^6`| dyn :math:`\cdot` cm :math:`^{-2}`  | Mean (reference) surface pressure (cgs)               |
+------------+----------------------------+-------------------------------------+-------------------------------------------------------+
|``rdgas``   | :math:`287.04`             | Jkg :math:`^{-1}`K :math:`^{-1}`    | Dry air gas constant                                  |
+------------+----------------------------+-------------------------------------+-------------------------------------------------------+
|``kappa``   | :math:`2/7`                | dimensionless                       | Heat capacity ratio ( :math:`\gamma` for an ideal gas)|
+------------+----------------------------+-------------------------------------+-------------------------------------------------------+
|``es0``     | :math:`1.0`                | dimensionless                       | Humidity factor :math:`^\dagger`                      |
+------------+----------------------------+-------------------------------------+-------------------------------------------------------+

:math:`^\dagger`: The humidity factor controls the atmospheric humidity content via the expression for saturation vapor pressure *if* ``do_simple: True`` is set within the ``idealized_moist_phys`` namelist.

.. note:: If the mean surface pressure value is changed here, it is necessary to also set ``reference_sea_level_press`` from the ``spectral_dynamics`` namelist, else the output file will not extend to the pressure specified.

Moist Atmospheres
^^^^^^^^^^^^^^^^^
We hope to be able to introduce controls for moist atmospheres in a future update to Isca. This will allow condensable species other than water to be easily modelled.

Relevant Modules
----------------
Since this module provides the definition of a number of physical constants, it is used by most other modules that exist within the Isca framework.

Authors
-------
This documentation was written by Daniel Williams, peer reviewed by Stephen Thompson, and quality controlled by Ross Castle.
