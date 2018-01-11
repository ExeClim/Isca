import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

import f90nml


NCORES = 4
base_dir = os.path.dirname(os.path.realpath(__file__))
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
exp = Experiment('realistic_continents_qflux_test_experiment', codebase=cb)

exp.inputfiles = [os.path.join(GFDL_BASE,'input/land_masks/era_land_t42.nc'),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
                      os.path.join(base_dir,'input/ami_qflux_ctrl_ice_4320.nc'), os.path.join(base_dir,'input/siconc_clim_amip.nc')]

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
        'two_stream_gray':False, #Don't use grey radiation
        'do_rrtm_radiation':True, #Do use RRTM radiation
        'convection_scheme': 'FULL_BETTS_MILLER', #Use the full betts-miller convection scheme
        'land_option': 'input', #Use land mask from input file
        'land_file_name': 'INPUT/era_land_t42.nc', #Tell model where to find input file
        'land_roughness_prefactor' : 10.0, #How much rougher to make land than ocean
        'roughness_mom' : 2.e-4, #Ocean roughness lengths
        'roughness_heat' : 2.e-4, #Ocean roughness lengths
        'roughness_moist' : 2.e-4, #Ocean roughness lengths        
    },

    'surface_flux_nml': {
        'land_humidity_prefactor': 0.7 #Evaporative resistance over land
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'delta_T' : 0., #Set latitude contrast in initial temperature profile to zero
        'depth' : 20., #Mixed layer depth
        'land_option' : 'input', #Tell mixed layer to get land mask from input file
        'land_h_capacity_prefactor': 0.1, #What factor to multiply mixed-layer depth by over land. 
        'albedo_value' : 0.25, #Ocean albedo value
        'land_albedo_prefactor': 1.3, #What factor to multiply ocean albedo by over land    
        'do_qflux' : False, #Don't use analytic formula for q-fluxes 
        'load_qflux' : True, #Do load q-flux field from an input file
        'time_varying_qflux' : True, #q-flux will be time-varying
        'qflux_file_name' : 'ami_qflux_ctrl_ice_4320', #Name of q-flux input file
        'update_albedo_from_ice' : True, #Use the simple ice model to update surface albedo
        'ice_albedo_value' : 0.7, #What value of albedo to use in regions of ice
        'ice_concentration_threshold' : 0.5, #ice concentration threshold above which to make albedo equal to ice_albedo_value            
    },

    'damping_driver_nml': {
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
    },

    'rrtm_radiation_nml': {
        'dt_rad': 4320, # Use 4320 as RRTM timestep        
    },

    'spectral_dynamics_nml': {
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'ocean_topog_smoothing': 0.8 #Use model's in-built spatial smoothing to smooth topography in order to prevent unwanted aliasing at low horizontal resolution     
    },

    'spectral_init_cond_nml': {
        'topog_file_name': 'era_land_t42.nc', #Name of land input file, which will also contain topography if generated using Isca's `land_file_generator_fn.py' routine.
        'topography_option':'input' #Tell model to get topography from input file
        }
})

#Lets do a run!
if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2,121):
        exp.run(i, num_cores=NCORES)