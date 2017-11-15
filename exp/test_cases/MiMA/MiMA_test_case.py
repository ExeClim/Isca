import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

import sys
sys.path.insert(0, '../')

from namelist_basefile import namelist_base

NCORES = 4

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = IscaCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('mima_test_experiment', codebase=cb)

exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

exp.diag_table = diag


#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = namelist_base  # Calls some defaults from test_cases/namelist_basefile.py

exp.update_namelist({
    'idealized_moist_phys_nml': {
        'two_stream_gray': False,
        'do_rrtm_radiation': True,    #Use RRTM radiation, not grey
        'convection_scheme': 'SIMPLE_BETTS_MILLER'     #Use the simple Betts Miller convection scheme          
    },
    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'depth': 100,
        'albedo_value': 0.205,
        'do_qflux': True        
    },

    'qflux_nml': {
        'qflux_amp': 30.0
    },

    'rrtm_radiation_nml': {
        'dt_rad': 7200, #Use long RRTM timestep
    }    
})

#Lets do a run!
exp.run(1, use_restart=False, num_cores=NCORES)
for i in range(2,121):
    exp.run(i, num_cores=NCORES)
