import numpy as np

from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress

NCORES = 16

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = DryCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create a diagnostics output file for daily snapshots
diag = DiagTable()
diag.add_file('atmos_daily', 1, 'days', time_units='days')

# Define diag table entries
diag.add_field('dynamics', 'ps', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'bk', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'pk', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'vor', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'div', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'ucomp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'vcomp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'temp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'slp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'omega', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'height', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'height_half', time_avg=True, files=['atmos_daily'])

diag.add_field('hs_forcing', 'teq', time_avg=True, files=['atmos_daily'])
diag.add_field('hs_forcing', 'h_trop', time_avg=True, files=['atmos_daily'])

# define namelist values as python dictionary
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 300,
        'days': 90,
        'calendar': 'no_calendar'
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    'spectral_dynamics_nml': {
        'num_levels': 30,
        'exponent': 2.5,
        'scale_heights': 4,
        'surf_res': 0.1,
        'robert_coeff': 4e-2,
        'do_water_correction': False,
        'vert_coord_option': 'even_sigma',
        'initial_sphum': 0.,
        'valid_range_T': [0, 700]
    },

    # configure the relaxation profile
    'hs_forcing_nml': {
        'equilibrium_t_option' : 'top_down',
        'ml_depth': 10.,
        'spinup_time': 10800,
        'ka': -20.,
        'ks': -5.
    },

    'constants_nml': {
        'orbital_period': 360,
    },

    'astronomy_nml': {
        'obliq' : 15
    },

    'diag_manager_nml': {
        'mix_snapshot_average_fields': False
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    }
})

if __name__=="__main__":

    obls = [15]
    for obl in obls:
        exp = Experiment('top_down_test_obliquity%d' % obl, codebase=cb)
        exp.clear_rundir()

        exp.diag_table = diag
        exp.namelist = namelist.copy()
        exp.namelist['astronomy_nml']['obliq'] = obl

        with exp_progress(exp, description='o%.0f d{day}' % obl):
            exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
        for i in range(2, 21):
            with exp_progress(exp, description='o%.0f d{day}' % s):
                exp.run(i, num_cores=NCORES, overwrite_data=True)