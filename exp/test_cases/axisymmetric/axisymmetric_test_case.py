import os

import numpy as np

import f90nml

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


NCORES = 8
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
exp = Experiment('axisymmetric_test_case', codebase=cb)

exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
                      os.path.join(base_dir,'input/sn_1.000_sst.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'bucket_depth', time_avg=True)
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
        'two_stream_gray':  False, #Don't use grey radiation
        'do_rrtm_radiation':  True, #Do use RRTM radiation
        'convection_scheme':  'SIMPLE_BETTS_MILLER', #Use the simple betts-miller convection scheme
    },
    
    'damping_driver_nml': {
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
    },
    
    'rrtm_radiation_nml': {
        'dt_rad': 3600. #Set RRTM radiation timestep to 3600 seconds, meaning it runs every 5 atmospheric timesteps
    },
    
    'mixed_layer_nml': {
        'albedo_value': 0.25,  #Ocean albedo value
        'do_sc_sst': True,  # Use precribed SSTs
        'sst_file': 'sn_1.000_sst', # SSTs to use            
    },

    'spectral_dynamics_nml': {   
# RG The model can be made symmetric by setting the option below to True. However, with no eddies, the equinoctial state, with two (predominantly eddy driven) Hadley cells symmetric about the equator, is not stable. With a mixed layer ocean the model finds other solutions to avoid this state, e.g. keeping the ITCZ off the equator, and the SST becomes very flat in the tropics. Additionally, both with a mixed layer ocean and with fixed SSTs, the cells tend to develop pressure level scale vertical waves. I therefore recommend using prescribed SSTs, and including some vertical diffusion in the free atmosphere to help dissipate vertical waves. For this test case seasonally varying SSTs are provided.     
        'make_symmetric': True, # Make model zonally symmetric
    },  

    'diffusivity_nml': {
        'free_atm_diff':True,#Turn on vertical diffusion in the free atmosphere 
    },
})



#Lets do a run!
exp.run(1, use_restart=False, num_cores=NCORES)
for i in range(2,121):
    exp.run(i, num_cores=NCORES)
    