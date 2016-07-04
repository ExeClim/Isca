# A basic Earth-like setup.
# - Uses nearly all standard parameter values.
# - Grey radiation scheme
# - Diurnal cycle
# - No seasons
# - Dry convection scheme and no humidity tracer

import numpy as np

import gfdl.experiment

import diagtables
import namelists

# create an experiment based on the code version tagged `dry0.1`
exp = gfdl.experiment.Experiment('exp16_almost_locked',
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='dry0.1')

# compiles source code to exp.execdir
exp.disable_rrtm()
exp.compile()

exp.namelist = namelists.moist.copy()

# don't use a calendar, use 360 day "months"
exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*100,
    'calendar': 'no_calendar'
}

omega = 7.2921150e-5
orbital_period = 2*np.pi / (omega*1.2)

exp.namelist['constants_nml'] = {
    'omega': omega,
    'orbital_period': orbital_period
}

exp.namelist['astronomy_nml'] = {
    'ecc': 0.0,
    'obliq': 0.0
}

# add the diag_table setup to the experiment
exp.use_diag_table(diagtables.basic)

# clean up previous runs.
exp.clear_rundir()

# run month 1 from a cold start
exp.run(1, use_restart=False)
for i in range(2, 6):
    # run subsequent months (default is to find the previous month
    # and use that as restart).
    exp.run(i)
