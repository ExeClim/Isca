# A basic Earth-like setup.
# - Uses nearly all standard parameter values.
# - Grey radiation scheme
# - No seasonal cycle - p2 like insolation profile

import numpy as np

import gfdl.experiment

# create an experiment based on the code version tagged `exoplan0.3`
# exp = gfdl.experiment.Experiment('dry_earth',
#     repo='git@github.com:jamesp/GFDLmoistModel.git',
#     commit='hot-n-dry',
#     overwrite_data=True)
exp = gfdl.experiment.Experiment('dry_earth', overwrite_data=True)


# compiles source code to exp.execdir
exp.disable_rrtm()
exp.compile()

# setup the namelist:
# - No limits on temperature
# - Initial temperature of 1200K
# - 25 vertical levels (26 half levels)
# - tidally locked at Earth rotation rate
exp.namelist['spectral_dynamics_nml']['num_levels'] = 25
exp.namelist['spectral_dynamics_nml']['valid_range_t'] = [100, 5000]

# Use the Frierson radiation scheme, with day-night sides.
exp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'frierson'
exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True

exp.namelist['astronomy_nml'] = {
    'ecc': 0.0,
    'obliq': 0.0
}

omega = 7.2921150e-5
orbital_period = 2*np.pi / omega
exp.namelist['constants_nml'] = {
    'omega': omega,
    'orbital_period': orbital_period
}

# don't use a calendar, but do use 30 day "months"
exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*5,
    'calendar': 'no_calendar'
}

# Setup a diag_table
# - Record daily and every 6 hours
# - Basic primitive equation prognostic var output
# - Radiation fluxes

diag = gfdl.experiment.DiagTable()

diag.add_file('6hourly', 6, 'hours')
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

# diag.add_field('mixed_layer', 't_surf')
# diag.add_field('mixed_layer', 'flux_oceanq')

# add the diag_table setup to the experiment
exp.use_diag_table(diag)

# clean up previous runs.
exp.clear_rundir()

# run month 1 from a cold start
exp.runmonth(1, use_restart=False)
for i in range(2, 40):
    # run subsequent months (default is to find the previous month
    # and use that as restart).
    exp.runmonth(i)
