import os
import numpy as np
from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress

NCORES = 16
NUM_LEVELS = 25

base_dir = os.path.dirname(os.path.realpath(__file__))
cb = SocratesCodeBase.from_directory(GFDL_BASE)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'slp', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf')
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('atmosphere', 'temp_2m', time_avg=True)
diag.add_field('atmosphere', 'q_2m', time_avg=True)
diag.add_field('atmosphere', 'rh_2m', time_avg=True)
diag.add_field('atmosphere', 'u_10m', time_avg=True)
diag.add_field('atmosphere', 'v_10m', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)
diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
diag.add_field('atmosphere', 'precipitation', time_avg=True)

diag.add_field('atmosphere', 'bucket_depth', time_avg=True)
diag.add_field('atmosphere', 'bucket_depth_cond', time_avg=True)
diag.add_field('atmosphere', 'bucket_depth_conv', time_avg=True)
diag.add_field('atmosphere', 'bucket_depth_lh', time_avg=True)

diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'height_half', time_avg=True)

diag.add_field('socrates', 'soc_tdt_lw', time_avg=True)
diag.add_field('socrates', 'soc_tdt_sw', time_avg=True)
diag.add_field('socrates', 'soc_tdt_rad', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_lw_down', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw_down', time_avg=True)

diag.add_field('socrates', 'soc_olr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw_up', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw_down', time_avg=True)

diag.add_field('socrates', 'soc_olr_clr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw_clr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw_up_clr', time_avg=True)
diag.add_field('socrates', 'soc_flux_lw_clr', time_avg=True)
diag.add_field('socrates', 'soc_flux_sw_clr', time_avg=True)

diag.add_field('socrates', 'soc_surf_flux_lw_clr', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw_clr', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_lw_down_clr', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw_down_clr', time_avg=True)

diag.add_field('cloud_simple', 'cf', time_avg=True)
diag.add_field('cloud_simple', 'reff_rad', time_avg=True)
diag.add_field('cloud_simple', 'frac_liq', time_avg=True)
diag.add_field('cloud_simple', 'qcl_rad', time_avg=True)
#diag.add_field('cloud_simple', 'simple_rhcrit', time_avg=True)
diag.add_field('cloud_simple', 'rh_in_cf', time_avg=True)

diag.add_field('cloud_simple', 'tot_cld_amt', time_avg=True)
diag.add_field('cloud_simple', 'high_cld_amt', time_avg=True)
diag.add_field('cloud_simple', 'mid_cld_amt', time_avg=True)
diag.add_field('cloud_simple', 'low_cld_amt', time_avg=True)

diag.add_field('cloud_simple', 'tot_cld_amt_mxr', time_avg=True)
diag.add_field('cloud_simple', 'high_cld_amt_mxr', time_avg=True)
diag.add_field('cloud_simple', 'mid_cld_amt_mxr', time_avg=True)
diag.add_field('cloud_simple', 'low_cld_amt_mxr', time_avg=True)

diag.add_field('cloud_simple', 'tot_cld_amt_max', time_avg=True)
diag.add_field('cloud_simple', 'high_cld_amt_max', time_avg=True)
diag.add_field('cloud_simple', 'mid_cld_amt_max', time_avg=True)
diag.add_field('cloud_simple', 'low_cld_amt_max', time_avg=True)

diag.add_field('cloud_simple', 'tot_cld_amt_rnd', time_avg=True)
diag.add_field('cloud_simple', 'high_cld_amt_rnd', time_avg=True)
diag.add_field('cloud_simple', 'mid_cld_amt_rnd', time_avg=True)
diag.add_field('cloud_simple', 'low_cld_amt_rnd', time_avg=True)

# diag.add_field('cloud_simple', 'low_cld_amt_park', time_avg=True)
# diag.add_field('cloud_simple', 'marine_strat', time_avg=True)
# diag.add_field('cloud_simple', 'eis', time_avg=True)
# diag.add_field('cloud_simple', 'ectei', time_avg=True)
# diag.add_field('cloud_simple', 'lts', time_avg=True)
# diag.add_field('cloud_simple', 'zlcl', time_avg=True)
# diag.add_field('cloud_simple', 'z700', time_avg=True)
# ##diag.add_field('cloud_simple', 'gamma850', time_avg=True)
# diag.add_field('cloud_simple', 'gamma700', time_avg=True)
# diag.add_field('cloud_simple', 'gamma_DL', time_avg=True)
# diag.add_field('cloud_simple', 'theta', time_avg=True)
# diag.add_field('cloud_simple', 'dthdp', time_avg=True)
# #diag.add_field('cloud_simple', 'theta0', time_avg=True)
# #diag.add_field('cloud_simple', 'theta700', time_avg=True)
# diag.add_field('cloud_simple', 'ELF', time_avg=True)
# diag.add_field('cloud_simple', 'beta1', time_avg=True)
# diag.add_field('cloud_simple', 'beta2', time_avg=True)
# diag.add_field('cloud_simple', 'zinv', time_avg=True)
# diag.add_field('cloud_simple', 'alpha', time_avg=True)
# diag.add_field('cloud_simple', 'DS', time_avg=True)
# diag.add_field('cloud_simple', 'IS', time_avg=True)
# #diag.add_field('cloud_simple', 'RH_inv_minus', time_avg=True)
# #diag.add_field('cloud_simple', 'qv_inv_minus', time_avg=True)
# #diag.add_field('cloud_simple', 'qv_inv_plus', time_avg=True)
# #diag.add_field('cloud_simple', 'qs_inv_minus', time_avg=True)
# #diag.add_field('cloud_simple', 'es_inv_minus', time_avg=True)
# #diag.add_field('cloud_simple', 'temp_inv_minus', time_avg=True)
# #diag.add_field('cloud_simple', 'pinv', time_avg=True)
# diag.add_field('cloud_simple', 'conv_cf', time_avg=True)

diag.add_field('mixed_layer', 'albedo', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)    # latent heat flux (up) at surface
diag.add_field('mixed_layer', 'flux_t', time_avg=True)      # sensible heat flux (up) at surface
# additional output options commented out
diag.add_field('socrates', 'soc_flux_lw', time_avg=True)
diag.add_field('socrates', 'soc_flux_sw', time_avg=True)
#diag.add_field('socrates', 'soc_co2', time_avg=True)
#diag.add_field('socrates', 'soc_ozone', time_avg=True)
#diag.add_field('socrates', 'soc_coszen', time_avg=True)
#diag.add_field('socrates', 'soc_spectral_olr', time_avg=True)


inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'), os.path.join(base_dir,'input/sst_clim_amip.nc'),
              os.path.join(GFDL_BASE,'input/land_masks/era_land_t42_filtered.nc'), os.path.join(base_dir, 'input/siconc_clim_amip.nc')]

exp = Experiment('exp_soc_bucket_linear_temp_minus_5_40_radii_14_25_correct', codebase=cb)

exp.clear_rundir()
exp.diag_table = diag
exp.inputfiles = inputfiles

#Define values for the 'core' namelist
exp.namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    'socrates_rad_nml': {
        'stellar_constant':1370.,
        'lw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'),
        'sw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'),
        'do_read_ozone': True,
        'ozone_file_name':'ozone_1990',
        'ozone_field_name':'ozone_1990',
        'dt_rad':3600,
        'store_intermediate_rad':True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
        'tidally_locked':False,
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,
        'two_stream_gray': False,     #Use the grey radiation scheme
        'do_socrates_radiation': True,
        'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use simple Betts miller convection
        'do_cloud_simple': True,
        'do_clear_sky_pass': True,
        'land_option' : 'input',
        'land_file_name' : 'INPUT/era_land_t42_filtered.nc',
        'land_roughness_prefactor' :10.0,
        'roughness_mom': 2.e-04, #Ocean roughness lengths  
        'roughness_heat': 2.e-04, #Ocean roughness lengths  
        'roughness_moist': 2.e-04, #Ocean roughness lengths  
        'bucket':True, #Run with the bucket model
        'init_bucket_depth_land':0.15, 
    },

    'cloud_simple_nml': {
        'cf_diag_formula_name': 'linear',
        'do_read_scm_rhcrit': False,
        #'do_simple_rhcrit': True,
        'do_qcl_with_temp': True,
        'do_adjust_low_cld': False,
        #'omega_adj_threshold': -0.1, 
        'do_cloud_amount_diags': True,
        'do_test_overlap': True,
        'do_add_stratocumulus': False, #True, #,
        #'sc_diag_method': 'Park_ELF',
        #'intermediate_outputs_diags': True,
        #'dthdp_min_threshold': -0.05, #-0.07, #-0.06,
        #'do_read_ts': True,
        'do_conv_cld': False, #False,
        #'conv_rain_min': 1,
        #'pshallow': 8.0e4,
        #'cf_min': 1e-2,
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
        'old_dtaudv': True,
        #'land_humidity_prefactor': 0.6,
        #'land_evap_prefactor': 1,
        'land_humidity_prefactor': 1,
        'land_evap_prefactor': 0.6,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst': 285.,
        'prescribe_initial_dist': True,
        'evaporation': True,
        'depth': 20.0,                          # Depth of mixed layer used
        'land_option': 'input',                 # Tell mixed layer to get land mask from input file
        'land_h_capacity_prefactor': 0.1,       # What factor to multiply mixed-layer depth by over land.
        'albedo_value': 0.11,                    # Ocean albedo value
        'land_albedo_prefactor': 1.3,           # What factor to multiply ocean albedo by over land
        'do_qflux': False,                      # Don't use the prescribed analytical formula for q-fluxes
        'do_read_sst': True,                    # Read in sst values from input file
        'do_sc_sst': True,                      # Do specified ssts (need both to be true)
        'sst_file': 'sst_clim_amip',            # Set name of sst input file
        'specify_sst_over_ocean_only': True,    # Make sure sst only specified in regions of ocean.
        # Copy from realistic_continents namelist
        'update_albedo_from_ice': True,          # Use the simple ice model to update surface albedo
        'ice_albedo_value': 0.7,                # What value of albedo to use in regions of ice
        'ice_concentration_threshold': 0.5,      # ice concentration threshold above which to make albedo equal to ice_albedo_value        
    },
        # Copy from realistic_continents namelist
    'qe_moist_convection_nml': {
        'rhbm': 0.7,
        'Tmin': 160.,
        'Tmax': 350.,
    },

    'lscale_cond_nml': {
        'do_simple': True,
        'do_evap': True,
    },

    'sat_vapor_pres_nml': {
           'do_simple': True,
           'construct_table_wrt_liq_and_ice': True,
           'show_all_bad_values': True,
       },

    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,               # neg. value: time in *days*
        'sponge_pbottom': 150.,         # Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,
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
        'reference_sea_level_press': 1.0e5,
        'num_levels': NUM_LEVELS,      #How many model pressure levels to use
        'valid_range_t': [100.,800.],
        'initial_sphum': [2.e-6],
        'vert_coord_option': 'uneven_sigma',
        'surf_res': 0.1, #0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights': 11.0,
        'exponent': 7.0,
        'robert_coeff': 0.03,
        'ocean_topog_smoothing': 0.8,
    },

    'spectral_init_cond_nml':{
        'topog_file_name': 'era_land_t42_filtered.nc',
        'topography_option': 'input',
    },
})

if __name__=="__main__":
    print(exp.namelist)
    cb.compile(debug=False)
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)#, run_idb=True)
    for i in range(2, 49):
        exp.run(i, num_cores=NCORES, overwrite_data=False)
    '''
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)#, run_idb=True)
    for i in range(2, 25):
        exp.run(i, num_cores=NCORES, overwrite_data=True)
    '''


