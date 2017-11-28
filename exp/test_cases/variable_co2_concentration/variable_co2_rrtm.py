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
exp = Experiment('variable_co2_rrtm_test_experiment', codebase=cb)

exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),os.path.join(base_dir,'input/co2.nc')]

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
diag.add_field('rrtm_radiation', 'co2', time_avg=True)

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
        'two_stream_gray': False, #Use RRTM, not grey radiation:
        'do_rrtm_radiation':True,
        'convection_scheme': 'FULL_BETTS_MILLER' #Use the full Betts-miller convection scheme
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'depth': 5., #Use shallow mixed-layer depth
        'albedo_value': 0.25, #set albedo value
        'do_qflux' : False, #Do not use prescribed form for q-fluxes            
    },
  
    'damping_driver_nml': {
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
    },

    'rrtm_radiation_nml': {
        'dt_rad': 4320, #Use 4320 as RRTM radiation timestep
        'do_read_co2': True, #Read in CO2 timeseries from input file
        'co2_file': 'co2' #Tell model name of co2 input file        
    },

    'spectral_dynamics_nml': {
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
    }
})    

#Lets do a run!
exp.run(1, use_restart=False, num_cores=NCORES)
for i in range(2,121):
    exp.run(i, num_cores=NCORES)
