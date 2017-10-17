# An example of using the gfdl library with a specific tagged commit of
# the source code

import namelists

import numpy as np

import gfdl.experiment

import sys

P = gfdl.experiment.P

args = sys.argv[1:]
ratios = [int(arg) for arg in args]

#omega = 7.2921150e-5
orbital_period = 365.25*86400.0

#ratios = [360, 180, 90, 45, 30, 15, 4, 2, 1]

# Use one experiment to compile the source.  All other experiments
# then use the same code with different namelist config
baseexp = gfdl.experiment.Experiment('exp2_base',
    repo='git@github.com:ExeClim/GFDLmoistModel.git',
    commit='exoplan0.2')

baseexp.log.info('Running ratios %r' % ratios)
baseexp.compile()

baseexp.namelist = namelists.basic

baseexp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*30,
    'calendar': 'no_calendar'
}

baseexp.namelist['astronomy_nml'] = {
    'ecc': 0.0,
    'obliq': 0.0
}


diag = gfdl.experiment.DiagTable()

#diag.add_file('6hourly', 6*60*60, 'seconds')
diag.add_file('daily', 1, 'days')

diag.add_field('dynamics', 'ps')
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'omega')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')
diag.add_field('dynamics', 'sphum')
diag.add_field('dynamics', 'slp')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')
diag.add_field('two_stream', 'tdt_rad')

diag.add_field('atmosphere', 'rh')

diag.add_field('mixed_layer', 'flux_t')


for ratio in ratios:
    exp = gfdl.experiment.Experiment('exp2_ratio_%d' % ratio)
    exp.clear_rundir()

    omega  = 2*np.pi / orbital_period * ratio

    exp.use_diag_table(diag)
    exp.execdir = baseexp.execdir

    exp.namelist = baseexp.namelist.copy()

    exp.namelist['constants_nml'] = {
        'omega': omega,
        'orbital_period': orbital_period
    }

    exp.runmonth(1, use_restart=False)
    for i in range(2, 161):
        exp.runmonth(i)
