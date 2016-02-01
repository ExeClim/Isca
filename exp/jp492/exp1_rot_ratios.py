import numpy as np

import gfdl.experiment

#omega = 7.2921150e-5
orbital_period = 365.25*86400.0

# Use one experiment to compile the source.  All other experiments
# then use the same code with different namelist config
baseexp = gfdl.experiment.Experiment('rot_base',
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='exoplan0.2')

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


diag = gfdl.experiment.DiagTable()

diag.add_file('6hourly', 6*60*60, 'seconds')
diag.add_file('daily', 1, 'days')

diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')
diag.add_field('two_stream', 'tdt_rad')
diag.add_field('two_stream', 'tdt_solar')

for ratio in [360.0, 180.0, 90.0, 45.0, 15.0, 2.0, 1.0]:
    exp = gfdl.experiment.Experiment('ratio_%d' % ratio)
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
    for i in range(2, 123):
        exp.runmonth(i)
