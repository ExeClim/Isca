# A basic Earth-like setup.
# - Uses nearly all standard parameter values.
# - Grey radiation scheme
# - No seasonal cycle - p2 like insolation profile

import numpy as np

import gfdl.experiment

# Use one experiment to compile the source.  All other experiments
# then use the same code with different namelist config
exp = gfdl.experiment.Experiment('exp5_byrne_rad_scheme',
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='23628d3c1afdf386452fb13dc531f95822eddad2')

exp.compile()

exp.namelist['spectral_dynamics_nml']['num_levels'] = 25

exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False

exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False
exp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'byrne'


exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*30,
    'calendar': 'no_calendar'
}

diag = gfdl.experiment.DiagTable()

#diag.add_file('6hourly', 6, 'hours')
diag.add_file('daily', 1, 'days')

diag.add_field('dynamics', 'ps')
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')
diag.add_field('dynamics', 'sphum')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')
diag.add_field('two_stream', 'tdt_rad')
diag.add_field('two_stream', 'lw_dtrans')

diag.add_field('mixed_layer', 't_surf')
diag.add_field('mixed_layer', 'flux_oceanq')


exp.clear_rundir()

exp.use_diag_table(diag)

exp.runmonth(1, use_restart=False)
for i in range(2, 40):
    exp.runmonth(i)
