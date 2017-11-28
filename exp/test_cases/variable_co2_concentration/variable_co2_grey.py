import os

import numpy as np

import f90nml

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


NCORES = 16
base_dir=os.getcwd()
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
exp = Experiment('variable_co2_grey_test_experiment', codebase=cb)

exp.inputfiles = [os.path.join(base_dir,'input/co2.nc')]

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
diag.add_field('two_stream', 'co2', time_avg=True)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
namelist_name = os.path.join(GFDL_BASE, 'exp/test_cases/namelist_basefile.nml')
nml = f90nml.read(namelist_name)
exp.namelist = nml

exp.update_namelist({
    'main_nml':{
         'days'   : 30,
         'hours'  : 0,
         'minutes': 0,
         'seconds': 0,
         'dt_atmos':720,
         'current_date' : [0001,1,1,0,0,0],
         'calendar' : 'thirty_day'
    },
    
    'idealized_moist_phys_nml': {       
        'two_stream_gray': True,     #Use the grey radiation scheme
        'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use simple Betts miller convection            
    },

    'mixed_layer_nml': {
        'depth': 2.5,                          #Depth of mixed layer used
        'albedo_value': 0.38,                  #Albedo value used      
    },

    'damping_driver_nml': {
        'trayfric': -0.25,              # neg. value: time in *days*
        'sponge_pbottom':  5000., #Bottom of the model's sponge down to 50hPa
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme':  'byrne',        #Select radiation scheme to use
        'equinox_day':  0.75,          #A calendar parameter to get autumn equinox in september, as in the standard earth calendar.
        'do_read_co2':  True, #Read in CO2 timeseries from input file
        'co2_file':  'co2', #Tell model name of co2 input file        
    },

    'spectral_dynamics_nml': {
        'num_levels':25,      #How many model pressure levels to use
        'vert_coord_option':'input',#Use the vertical levels from Frierson 2006
    }
})

#Lets do a run!
exp.run(1, use_restart=False, num_cores=NCORES)
for i in range(2,121):
    exp.run(i, num_cores=NCORES)
