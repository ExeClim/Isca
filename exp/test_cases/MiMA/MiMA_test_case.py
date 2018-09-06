import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 16

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

cb.compile()#debug=True)  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('mima_test_experiment_with_mcica', codebase=cb)

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
diag.add_field('cloud_simple', 'cf_rad', time_avg=True)
diag.add_field('cloud_simple', 'reff_rad', time_avg=True)
diag.add_field('cloud_simple', 'frac_liq', time_avg=True)
diag.add_field('cloud_simple', 'qcl_rad', time_avg=True)
diag.add_field('rrtm_radiation', 'olr', time_avg=True)

exp.diag_table = diag


#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml': {
        'days'   : 30,
        'hours'  : 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos':600,
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'two_stream_gray': False,
        'do_rrtm_radiation': True,    #Use RRTM radiation, not grey
        'convection_scheme': 'SIMPLE_BETTS_MILLER',     #Use the simple Betts Miller convection scheme
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,
        'do_cloud_simple': True,                
    },
    'cloud_simple_nml': {
        'simple_cca':0.0,
        'rhcsfc': 0.95,
        'rhc700': 0.7,
        'rhc200': 0.3,
        'rhmsfc': 0.95,
        'rhm700': 0.7,
        'rhm200': 0.3,
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'depth': 100,
        'albedo_value': 0.205,
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,
        'do_qflux': True        
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True,
        'construct_table_wrt_liq_and_ice':True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  150.,
        'do_conserve_energy': True,         
    },

    'qflux_nml': {
        'qflux_amp': 30.0
    },

    'rrtm_radiation_nml': {
        'solr_cnst': 1360,  #s set solar constant to 1360, rather than default of 1368.22
        'dt_rad': 7200, #Use long RRTM timestep
        'do_read_ozone':True,
        'ozone_file':'ozone_1990',
        'icld': 1 #cloud overlap method
    },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':40,
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2,
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    }
    
    
})
#Lets do a run!
if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)#, run_idb=True)
    for i in range(2,25):
        exp.run(i, num_cores=NCORES, overwrite_data=True)
