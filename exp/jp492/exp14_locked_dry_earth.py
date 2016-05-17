# A basic Earth-like setup.
# - Uses nearly all standard parameter values.
# - Grey radiation scheme
# - Diurnal cycle
# - No seasons
# - Dry convection scheme and no humidity tracer
import sys

import numpy as np
import sh

import gfdl.experiment

import diagtables
import namelists


args = sys.argv[1:]
depths = [int(arg) for arg in args]

# create an experiment based on the code version tagged `dry0.1`
baseexp = gfdl.experiment.Experiment('exp13_ref_dry_earth',
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='dry0.1')

# compiles source code to exp.execdir
baseexp.disable_rrtm()
baseexp.compile()

# add the dry namelist and diag_table to the experiment
baseexp.namelist = namelists.dry.copy()
baseexp.diag_table = diagtables.dry.copy()

# don't use a calendar, use 100 day "months"
baseexp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*10,
    'calendar': 'no_calendar'
}

omega = 7.2921150e-5
orbital_period = 365.25*86400.0  # standard Earth orbital period

baseexp.namelist['constants_nml'] = {
    'omega': omega,
    'orbital_period': orbital_period
}

# no seasons
baseexp.namelist['astronomy_nml'] = {
    'ecc': 0.0,
    'obliq': 0.0
}


spinup = baseexp.derive('exp14_spinup')
spinup.runmonth(1, use_restart=False)
for i in range(2, 6):
    spinup.runmonth(i)

for depth in depths:
    exp = baseexp.derive('exp14_depth_%d' % depth)

    omega = 7.2921150e-5
    orbital_period = 2*np.pi / omega  # tidally lock

    exp.namelist['constants_nml'] = {
        'omega': omega,
        'orbital_period': orbital_period
    }

    exp.namelist['mixed_layer_nml']['depth'] = depth

    # clean up previous runs.
    exp.clear_rundir()

    try:
        # run month 1 from a start based on reference dry earth
        exp.runmonth(1, restart_file=spinup.get_restart_file(3))
        for i in range(2, 26):
            # run subsequent months (default is to find the previous month
            # and use that as restart).
            exp.runmonth(i)
    except sh.ErrorCode_1:
        exp.log.error('Experiment run ended prematurely.')
        continue