# A basic Earth-like setup.
# - Uses nearly all standard parameter values.
# - Grey radiation scheme
# - Diurnal cycle
# - No seasons

import numpy as np

import gfdl.experiment

import diagtables
import namelists

# create an experiment based on the code version tagged `exoplan0.3`
base = gfdl.experiment.Experiment('exp13_ref_dry_earth',
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='dry0.1')

# compiles source code to exp.execdir
base.disable_rrtm()
base.compile()
base.diag_table = diagtables.basic.copy()

moist = base.derive('exp15_moist')
moist.namelist = namelists.moist.copy()

dry = base.derive('exp15_dry')
dry.namelist = namelists.dry.copy()


for exp in (dry, moist):
    # don't use a calendar, use 360 day "months"
    exp.namelist['main_nml'] = {
        'dt_atmos': 900,
        'seconds': 86400.0*100,
        'calendar': 'no_calendar'
    }

    omega = 7.2921150e-5
    orbital_period = 360.0*86400.0

    exp.namelist['constants_nml'] = {
        'omega': omega,
        'orbital_period': orbital_period
    }

    exp.namelist['astronomy_nml'] = {
        'ecc': 0.0,
        'obliq': 0.0
    }

    # clean up previous runs.
    exp.clear_rundir()

    # run month 1 from a cold start
    exp.runmonth(1, use_restart=False, num_cores=16)
    for i in range(4):
        year = i + 2
        exp.run(year, num_cores=16)
