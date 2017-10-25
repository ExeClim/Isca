import numpy as np
import os
from gfdl.experiment import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']

baseexp = Experiment('held_suarez_test_case', overwrite_data=False)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

baseexp.use_diag_table(diag)

#Turn off the full, slow radiation scheme compilation

baseexp.disable_rrtm()

#Compile model if not already compiled
baseexp.compile()

#Empty the run directory ready to run
baseexp.clear_rundir()

#Define values for the 'core' namelist
baseexp.namelist['main_nml'] = f90nml.Namelist({
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':720,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'thirty_day'
})

baseexp.namelist['atmosphere_nml']['idealized_moist_model'] = False #tells the model to use newtonian relaxation

baseexp.namelist['hs_forcing_nml'] = f90nml.Namelist({

    't_zero': 315.,    # temperature at reference pressure at equator (default 315K)
    't_strat': 200.,   # stratosphere temperature (default 200K)
    'delh': 60.,       # equator-pole temp gradient (default 60K)
    'delv': 10.,       # lapse rate (default 10K)
    'eps': 0.,         # stratospheric latitudinal variation (default 0K)
    'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)

    # negative sign is a flag indicating that the units are days
    'ka':   -40.,      # Constant Newtonian cooling timescale (default 40 days)
    'ks':    -4.,      # Boundary layer dependent cooling timescale (default 4 days)
    'kf':   -1.,       # BL momentum frictional timescale (default 1 days)

    'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
    })

#Lets do a run!
baseexp.runmonth(1, use_restart=False,num_cores=4, light=False)
for i in range(2,121):
    baseexp.runmonth(i, num_cores=4, light=False)
