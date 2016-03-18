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
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='exoplan0.3')

baseexp.log.info('Running ratios %r' % ratios)
baseexp.compile()

baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True
baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 25

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
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')
diag.add_field('dynamics', 'sphum')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')
diag.add_field('two_stream', 'lw_dtrans')

diag.add_field('atmosphere', 'rh')


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
    for i in range(2, 41):
        exp.runmonth(i)
