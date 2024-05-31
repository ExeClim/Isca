import os

import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES =16

# Point to code as defined by $GFDL_BASE
cb = GreyCodeBase.from_directory(GFDL_BASE)

base_dir = os.path.dirname(os.path.realpath(__file__))

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('thermo_ice_test_experiment', codebase=cb)


#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'temp_2m', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('mixed_layer', 'albedo', time_avg=True)
diag.add_field('mixed_layer', 'h_thermo_ice', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True)
diag.add_field('two_stream', 'swdn_toa', time_avg=True)



exp.diag_table = diag # register diag table 


#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml': {
        'days'   : 360, # each run lasts one year, and then multiple runs are strung together below (loop on e.g. Line 310)
        'hours'  : 0,   # a different output file is produced for each run (in this case, each year). Data 
        'minutes': 0,   # output at the frequency specified in the diag table 
        'seconds': 0,
        'dt_atmos':900, 
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'two_stream_gray': True, 
        'do_rrtm_radiation': False, 
        'convection_scheme': 'SIMPLE_BETTS_MILLER',                                          
        'do_damping': True,       
        'turb':True,          
        'mixed_layer_bc':True,                
        'do_virtual' :True, 
        'roughness_mom':5.e-3, 
        'roughness_heat':1.e-5, 
        'roughness_moist':1.e-5, 
        'land_roughness_prefactor':1.0,
        'do_simple':False,
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False, 
        'do_diffusivity': True, 
        'do_edt':False, 
        'constant_gust': 1.0,   
        'use_tau': False, 
        'do_entrain':False, 
        'do_stable_bl':False, 
        'do_shallow_conv':False, 
        'do_simple':False 
    },
    
    'diffusivity_nml': {
        'do_entrain':False, 
        'entr_ratio': 0.0, 
        'free_atm_diff':False, 
        'do_simple': False, 
        'parcel_buoy': 0.0, 
        'frac_inner': 0.1,  
        'fixed_depth': False, 
    },

    'surface_flux_nml': {
        'use_virtual_temp': True, 
        'do_simple': False, 
        'old_dtaudv': False,  
        'gust_const':1.0, 
        'land_humidity_prefactor' : 1.0, 
        'land_evap_prefactor': 1.0,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True 
    },


    'mixed_layer_nml': { 
        'depth': 30., 
        'albedo_value': 0.1, 
        'prescribe_initial_dist':True, 
        'tconst' : 305., 
        'delta_T': 60., 
        'evaporation':True, 
        'do_qflux': True,
        'load_qflux':False, 
        'time_varying_qflux' : False, 
        'do_thermo_ice':True, # turn on thermodynamic ice 
        'thermo_ice_albedo':0.55, # ice albedo 
        't_thermo_ice_base':273.15-2., # freezing point 
        't_surf_freeze':273.15-2., # freezing point 
        'do_var_thermo_ice_albedo':False,
        'read_const_correct':False, 
        'read_nudge_out':False, 
        
    },
    
    'qflux_nml':{
        'qflux_amp':30.,
    }, 
    


    'qe_moist_convection_nml': { 
        'rhbm':0.7, 
        'tau_bm':7200., 
        'Tmin':120., 
        'Tmax':360., 
        'val_inc': 0.01, 
        'precision':1.e-6
    },
    
    'lscale_cond_nml': {
        'do_simple':False, 
        'do_evap':False 
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':False, 
                          
    },    
    
    
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  2600.,
        'do_conserve_energy': True,         
    },
    



    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use
        'do_seasonal': True,                
        'use_time_average_coszen':True,
        'dt_rad_avg':86400.,
        'atm_abs': 0.22,                  
        'solar_exponent':2,
        'ir_tau_eq':7.2, 
        'ir_tau_pole':3.6,#1.8, 
        'del_sol':0.98, 
        'solar_constant':1360, 
        'linear_tau':0.2, 
        'odp':1.4, 
        'do_toa_albedo':True,
    },


    'spectral_dynamics_nml': {
        'damping_order': 4, # Yields lap^8 damping 
        'water_correction_limit': 200.e2, 
        'reference_sea_level_press':1.0e5, 
        'num_levels':30, 
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6], 
        'use_virtual_temperature':True,
        'vert_coord_option':'uneven_sigma', 
        'robert_coeff':0.03, 
        # set to T42 resolution 
        'lon_max': 128, 
        'lat_max': 64, 
        'num_fourier': 42, 
        'num_spherical':43, 
        'surf_res': 0.05, 
        'exponent': 3., 
        'scale_heights': 5 
        
    }, 
    
    'spectral_init_cond_nml':{
        'initial_temperature':280.
    }, 
   
    
    
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


    
    
})

# Lets do a run!
if __name__=="__main__":
    
    
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
    
    for i in range(1,11): # run for 10 years 
        exp.run(i, num_cores=NCORES, overwrite_data=False)
        
        

