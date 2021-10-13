Large Scale Condensation and Precipitation
 ==========================================

 Summary
 -------
     Temperature and specific humidity adjustments are computed in model layers
     where the relative humidity exceeds a threshold relative humidity.


    The features include:

     1) option for the re-evaporation of falling precipitation
     2) energetically consistent adjustment with precipitation type

 Computes the large-scale condensation adjustments for
                  temperature and specific humidity, and returns the
                  mass of rain and snow that reach the ground (in an
                  energetically consistent way)

                    Output

      rain     liquid precipitation (kg/m2)
                  [real, dimension(:,:)]

      snow     frozen precipitation (kg/m2)
                  [real, dimension(:,:)]

      tdel     temperature tendency at full model levels
                  [real, dimension(:,:,nlev)]

      qdel     specific humidity tendency (of water vapor) at
                 full model levels
                  [real, dimension(:,:,nlev)]

 Logic of the routine
 -------

 1. [Frierson2007]_: 

 Namelist options
 ----------------

 There are five optional namelist variables. 

 The first two are directly relevant to how the tendency equations and precipitation are calculated

 :tau_bm: The relaxation time (:math:`\tau`) with the default value of :math:`7200` :math:`s`.
 :rhbm: The critical Relative Humidity (RH) used to calculate the moisture reference profile (:math:`RH_{SBM}`) with the default value of :math:`0.8`.

 The final three are used to construct the look up table for calculating the lifting condensation level (LCL).

 :Tmin: The minimum look up table value with the default of :math:`173` :math:`K`.
 :Tmax: The maximum look up table value with the default of :math:`335` :math:`K`.
 :val_inc: The increment used to linearly interpolate within the look up table with the default of :math:`0.01`.


 &lscale_cond_nml

    hc         relative humidity at which large scale condensation
               occurs, where 0 <= hc <= 1
                  [real, default: hc=1.]

    do_evap    flag for the re-evaporation of moisture in
               sub-saturated layers below, if do_evap=.true. then
               re-evaporation is performed
                  [logical, default: do_evap=.false.]

   do_simple   if true then all precipitation is rain, there is no snow. Added into Isca, not in GFDL

 Diagnostics
 -----------
 The diagnostics are not sent out within the call to convection but are sent out from ``idealized_moist_phys.F90`` with the module name ``atmosphere``. The diagnostics relevant to the simple Betts-Miller scheme are in the table below.
 copied from betts miller
 +-------------------+----------------------------+------------------------------------+
 | Name              | Description                | Units                              |
 +===================+============================+====================================+
 | dt_qg_convection  | Moisture tendency          |:math:`kg~kg^{-1}~s^{-1}`           |
 +-------------------+----------------------------+------------------------------------+
 | dt_tg_convection  | Temperature tendency       |:math:`K~s^{-1}`                    |
 +-------------------+----------------------------+------------------------------------+
 | convection_rain   | Convective precipitation   |:math:`kg~m^{-2}~s^{-1}`            |
 +-------------------+----------------------------+------------------------------------+
 | cape              | Convective available       |:math:`J~kg^{-1}`                   |
 |                   | potential energy           |                                    |
 +-------------------+----------------------------+------------------------------------+
 | cin               | Convective inhibition      |:math:`J~kg^{-1}`                   |
 +-------------------+----------------------------+------------------------------------+


 Relevant modules and subroutines
 --------------------------------
 sat_vapor_pres_mod, utilities_mod, constants_mod

 References
 ----------

 The reference list within this docs page are: [Frierson2007]_, [Betts1986]_, [BettsMiller1986]_, [YanoPlant2016]_, and [OGormanSchneider2008]_. 

     Reference: Manabe, S., (1969). Mon. Wea. Rev. 97, 739-798.

 Authors
 -------
 This documentation was written by Ross Castle
