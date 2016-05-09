import numpy as np

import gfdl.experiment

import sys

P = gfdl.experiment.P

args = sys.argv[1:]
depths = [int(arg) for arg in args]

MOIST_MODEL = True

# Use one experiment to compile the source.  All other experiments
# then use the same code with different namelist config
baseexp = gfdl.experiment.Experiment('exp12_base',
    repo='/scratch/jp492/GFDLmoistModel',
    commit='b697ef4')

# Create a diag table and set some output files
diag = gfdl.experiment.DiagTable()
#diag.add_file('6hourly', 6*60*60, 'seconds')
diag.add_file('daily', 1, 'days')

baseexp.log.info('Running depths %r' % depths)
baseexp.disable_rrtm()
baseexp.compile()

### CONFIGURE THE EXPERIMENT

# vertical structure
baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 25
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2

# run 30 day 'months', no calendar
baseexp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*30,
    'calendar': 'no_calendar'
}

# No seasons, equinox with diurnal cycle
baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True
baseexp.namelist['astronomy_nml'] = {
    'ecc': 0.0,
    'obliq': 0.0
}


if MOIST_MODEL:
    baseexp.namelist['idealized_moist_phys_nml']['convection_scheme'] = 'betts_miller'
    baseexp.namelist['spectral_dynamics_nml']['initial_sphum'] = 2e-6
    diag.add_field('dynamics', 'sphum')
    diag.add_field('atmosphere', 'rh')

else:
    baseexp.namelist['idealized_moist_phys_nml']['convection_scheme'] = 'dry'
    baseexp.namelist['spectral_dynamics_nml']['initial_sphum'] = 0.0
    # diag.add_field('dry_convection', 'dp')
    # diag.add_field('dry_convection', 'CAPE')
    # diag.add_field('dry_convection', 'CIN')
    # diag.add_field('dry_convection',  'LZB')
    # diag.add_field('dry_convection', 'LCL')
    # diag.add_field('dry_convection', 'dt_tg')
    # diag.add_field('dry_convection', 'parcel_temp')

diag.add_field('dynamics', 'ps')
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'omega')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')

diag.add_field('dynamics', 'slp')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')
diag.add_field('two_stream', 'tdt_rad')

diag.add_field('mixed_layer', 'flux_t')


for depth in depths:
    exp = gfdl.experiment.Experiment('exp12_depth_%d' % depth)
    exp.clear_rundir()


    omega = 7.2921150e-5
    orbital_period = 2*np.pi/omega

    exp.use_diag_table(diag)
    exp.execdir = baseexp.execdir

    exp.namelist = baseexp.namelist.copy()

    exp.namelist['mixed_layer_nml'] = {
        'albedo_value': 0.27,
        'depth': depth,
        #'prescribe_initial_dist': True
        # 'tconst': 285.0,
        # 'delta_T': 40.0,
        'evaporation': MOIST_MODEL,
        'do_qflux': MOIST_MODEL
    }

    exp.namelist['constants_nml'] = {
        'omega': omega,
        'orbital_period': orbital_period
    }

    exp.runmonth(1, use_restart=False)
    for i in range(2, 81):
        exp.runmonth(i)
