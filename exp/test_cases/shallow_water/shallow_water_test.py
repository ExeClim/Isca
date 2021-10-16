import os

import numpy as np

from isca import ShallowCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 8
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = ShallowCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('shallow_test_experiment', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('shallow_diagnostics', 'ucomp', time_avg=True)
diag.add_field('shallow_diagnostics', 'vcomp', time_avg=True)
diag.add_field('shallow_diagnostics', 'vor', time_avg=True)
diag.add_field('shallow_diagnostics', 'div', time_avg=True)
diag.add_field('shallow_diagnostics', 'h', time_avg=True)
diag.add_field('shallow_diagnostics', 'pv', time_avg=True)
diag.add_field('shallow_diagnostics', 'stream', time_avg=True)
diag.add_field('shallow_diagnostics', 'trs', time_avg=True)
diag.add_field('shallow_diagnostics', 'tr', time_avg=True)


exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos': 1200,
     'calendar': 'no_calendar',
    },

 'atmosphere_nml':{
   'print_interval': 86400,
    },

'fms_io_nml':{
   'threading_write' :'single',
   'fileset_write': 'single'
    },

 'fms_nml':{
   'print_memory_usage':True,
   'domains_stack_size': 200000,
    },

 'shallow_dynamics_nml':{
   'num_lon'             : 256,
   'num_lat'             : 128,
   'num_fourier'         : 85,
   'num_spherical'       : 86,
   'fourier_inc'         : 1,
   'damping_option'      : 'resolution_dependent',
   'damping_order'       : 4,
   'damping_coeff'       : 1.e-04,
   'h_0'                 : 3.e04,
   'grid_tracer'         : True,
   'spec_tracer'         : True,
   'robert_coeff'        : 0.04,
   'robert_coeff_tracer' : 0.04,
    },

 'shallow_physics_nml':{
   'fric_damp_time'  : -50.0,
   'therm_damp_time' : -10.0,
   'del_h'           :  2.e04,
   'h_0'             :  3.e04,
   'h_amp'           :  1.e05,
   'h_lon'           :  90.0,
   'h_lat'           :  25.0,
   'h_width'         :  15.0,
   'itcz_width'      :  4.0,
   'h_itcz'          :  4.e04,
   },
})

#Lets do a run!
if __name__=="__main__":
    cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2,121):
        exp.run(i, num_cores=NCORES)
