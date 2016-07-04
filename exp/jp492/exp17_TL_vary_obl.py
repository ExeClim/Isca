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
base = gfdl.experiment.Experiment('base_dry0.1',
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='dry0.1')

# compiles source code to exp.execdir
base.disable_rrtm()
base.compile()

base.namelist = namelists.moist.copy()
base.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*100,
    'calendar': 'no_calendar'
}

# add the diag_table setup to the experiment
base.use_diag_table(diagtables.basic)

for obliq in (0.0, 5.0, 10.0, 15.0, 20.0, 25.0):
    exp = base.derive('exp17_obl_%d' % obliq)
    # don't use a calendar, use 360 day "months"

    omega = 7.2921150e-5
    orbital_period = 2*np.pi / omega  # tidally locked Earth

    exp.namelist['constants_nml'] = {
        'omega': omega,
        'orbital_period': orbital_period
    }

    exp.namelist['astronomy_nml'] = {
        'ecc': 0.0,
        'obliq': obliq
    }

    # clean up previous runs.
    exp.clear_rundir()

    # run month 1 from a cold start
    exp.run(1, use_restart=False, num_cores=16)
    for i in range(2, 6):
        # run subsequent months (default is to find the previous month
        # and use that as restart).
        exp.run(i, num_cores=16)
