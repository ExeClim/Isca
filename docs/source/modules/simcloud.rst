Cloud Scheme (SimCloud)
=======================

Summary
-------
Isca provides a simple diagnostic cloud scheme, SimCloud, originally developed and described in [Liu2021]_. The model diagnoses clouds in three different ways: (1) large-scale clouds from the relative humidity; (2) low-level marine stratus clouds from a function of local inversion strength; and (3) a "freeze-dry" adjustment based on a simple function of specific humidity is also available to reduce an excessive cloud bias in polar regions. Cloud properties such as the effective radius of cloud droplet and cloud liquid water content are specified as simple functions of temperature, and parameters of the model may be tuned individually to suit requirements.

The cloud scheme may be found in the directory ``src/atmos_param/cloud_simple``, with tuneable parameters found therein.

A demonstration test case is provided within Isca's test suite of experiments. The scheme is provided on an as-in basis: parameters have been tuned for an Earth-like demonstration so any deviations from such a setup or modifications to the parameterisation are at your own risk. The model should be used with the SOCRATES radiation scheme.

Diagnostics
-----------
There are a number of different diagnostic outputs available with the cloud scheme. These are found across different diagnostic groups which are detailed below.

Basic Diagnostics
^^^^^^^^^^^^^^^^^
These diagnostics are part of the ``cloud_simple`` module and describe the overall state of the cloud.

+--------------------------+-----------------------------------------------------+--------------------+-------------------------------+
| Name                     | Description                                         | Units              | Dimension (not including time)|
+==========================+=====================================================+====================+===============================+
|``cf``                    | Cloud fraction for the simple cloud scheme          | dimensionless (0-1)| (pfull, lat, lon)             |
+--------------------------+-----------------------------------------------------+--------------------+-------------------------------+
|``frac_liq``              | Liquid cloud fraction (liquid, mixed-ice phase, ice)| dimensionless (0-1)| (pfull, lat, lon)             |
+--------------------------+-----------------------------------------------------+--------------------+-------------------------------+
|``reff_rad``              | Effective cloud particle radius                     | :math:`\mu` m      | (pfull, lat, lon)             |
+--------------------------+-----------------------------------------------------+--------------------+-------------------------------+
|``qcl_rad``               | Specific humidity of cloud liquid                   | kg kg :math:`^{-1}`| (phalf, lat, lon)             |
+--------------------------+-----------------------------------------------------+--------------------+-------------------------------+
|``rh_in_cf``              | RH as a percent                                     | :math:`\%`         | (phalf, lat, lon)             |
+--------------------------+-----------------------------------------------------+--------------------+-------------------------------+

Cloud Coverage
^^^^^^^^^^^^^^
The module ``cloud_cover`` provides additional diagnostics to show spatial coverage of clouds in the lat-lon direction; additionally split into different levels of clouds.

+------------------+-----------------------------+------------+-------------------------------+
| Name             | Description                 | Units      | Dimension (not including time)|
+==================+=============================+============+===============================+
|``tot_cld_amt``   | Total cloud amount (%)      | :math:`\%` | (lat, lon)                    |
+------------------+-----------------------------+------------+-------------------------------+
|``high_cld_amt``  | High level cloud amount (%) | :math:`\%` | (lat, lon)                    |
+------------------+-----------------------------+------------+-------------------------------+
|``mid_cld_amt``   | Hid-level cloud amount (%)  | :math:`\%` | (lat, lon)                    |
+------------------+-----------------------------+------------+-------------------------------+
|``low_cld_amt``   | Low level cloud amount (%)  | :math:`\%` | (lat, lon)                    |
+------------------+-----------------------------+------------+-------------------------------+

Marine Stratospheric Clouds
^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are additional diagnostics specifically relating to marine stratospheric clouds available in the module ``strat_cloud``. These are primarily available for test purposes so are not detailed here, however they may be accessed in the routine ``marine_strat_cloud.F90`` within the cloud scheme directory.


Relevant modules and subroutines
--------------------------------
The SOCRATES radiation scheme is setup in ``src/atmos_param/socrates/interface/``.

References
----------
[Liu2021]_

Authors
-------
This documentation was written by Daniel Williams (with thanks to Qun Liu for implementing the model).