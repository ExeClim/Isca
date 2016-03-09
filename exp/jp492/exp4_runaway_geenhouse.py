import numpy as np

import gfdl.experiment

#omega = 7.2921150e-5
orbital_period = 365.25*86400.0

#multipliers = [0.5, 0.75, 1.0]
#multipliers = [1.25, 1.5, 2.0, 3.0]
multipliers = [1.0, 1.5, 2.0]

# Use one experiment to compile the source.  All other experiments
# then use the same code with different namelist config
baseexp = gfdl.experiment.Experiment('greenhouse_base',
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='exoplan0.3')

baseexp.compile()


baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False
baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 25
baseexp.namelist['two_stream_gray_rad_nml']['solar_exponent'] = -1
baseexp.namelist['two_stream_gray_rad_nml']['wv_exponent'] = -1

baseexp.namelist['main_nml'] = {
    'dt_atmos': 450,
    'seconds': 86400.0*30,
    'calendar': 'no_calendar'
}

diag = gfdl.experiment.DiagTable()

#diag.add_file('6hourly', 6*60*60, 'seconds')
diag.add_file('daily', 1, 'days')
diag.add_file('hourly', 1, 'hours')

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
diag.add_field('two_stream', 'lw_dtrans_win')
diag.add_field('two_stream', 'sw_dtrans')

diag.add_field('mixed_layer', 't_surf')


for s in multipliers:
    exp = gfdl.experiment.Experiment('greenhouse_eqnox_%03d' % (s*100))
    exp.clear_rundir()

    S0 = 1360.0 * s

    exp.use_diag_table(diag)
    exp.execdir = baseexp.execdir

    exp.namelist = baseexp.namelist.copy()

    exp.namelist['two_stream_gray_rad_nml']['solar_constant'] = S0

    exp.runmonth(1, use_restart=False)
    for i in range(2, 13):
        exp.runmonth(i)
