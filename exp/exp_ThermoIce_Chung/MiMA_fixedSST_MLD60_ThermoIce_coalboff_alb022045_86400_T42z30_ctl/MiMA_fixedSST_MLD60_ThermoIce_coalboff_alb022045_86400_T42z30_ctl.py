import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32

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
exp = Experiment('MiMA_fixedSST_MLD60_ThermoIce_coalboff_alb022045_86400_T42z30_ctl', codebase=cb)

#exp.inputfiles = [os.path.join(GFDL_BASE,'input/input_Chung/t_surf/t_surf_MiMActl_pen.nc')]
#exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

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
#
diag.add_field('dynamics', 'slp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'wspd', time_avg=True)
diag.add_field('dynamics', 'ucomp_sq', time_avg=True)
diag.add_field('dynamics', 'ucomp_vcomp', time_avg=True)
diag.add_field('dynamics', 'ucomp_omega', time_avg=True)
diag.add_field('dynamics', 'vcomp_sq', time_avg=True)
diag.add_field('dynamics', 'vcomp_omega', time_avg=True)
diag.add_field('dynamics', 'omega_sq', time_avg=True)
diag.add_field('dynamics', 'ucomp_temp', time_avg=True)
diag.add_field('dynamics', 'vcomp_temp', time_avg=True)
diag.add_field('dynamics', 'omega_temp', time_avg=True)
diag.add_field('dynamics', 'temp_sq', time_avg=True)
diag.add_field('dynamics', 'vcomp_vor', time_avg=True)
diag.add_field('rrtm_radiation', 'olr', time_avg=True) #
diag.add_field('rrtm_radiation', 'toa_sw', time_avg=True) #
diag.add_field('mixed_layer', 'albedo', time_avg=True)
diag.add_field('mixed_layer', 'ml_heat_cap', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True) #
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
#diag.add_field('mixed_layer', 'flux_sw', time_avg=True)
#diag.add_field('mixed_layer', 'flux_lw', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)
diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_rad', time_avg=True)
diag.add_field('vert_turb', 'z_pbl', time_avg=True)
diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True)
diag.add_field('rrtm_radiation', 'flux_lw', time_avg=True)
diag.add_field('rrtm_radiation', 'rrtm_albedo', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_sw', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_lw', time_avg=True)
diag.add_field('rrtm_radiation', 'ozone', time_avg=True) #
diag.add_field('rrtm_radiation', 'co2', time_avg=True) #
#

diag.add_field('mixed_layer', 'flux_u', time_avg=True) ##
diag.add_field('mixed_layer', 'h_ice', time_avg=True) ##
diag.add_field('mixed_layer', 'a_ice', time_avg=True) ##
diag.add_field('mixed_layer', 't_ml', time_avg=True) ##
diag.add_field('mixed_layer', 'flux_ice', time_avg=True) ##
#
diag.add_field('atmosphere', 'temp_2m', time_avg=True) #
diag.add_field('atmosphere', 'u_10m', time_avg=True) #
diag.add_field('atmosphere', 'v_10m', time_avg=True) #
diag.add_field('atmosphere', 'q_2m', time_avg=True) #
diag.add_field('atmosphere', 'rh_2m', time_avg=True) #
#
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
        'depth': 60,
#        'albedo_value': 0.1,
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,
        'do_qflux': False,
        'do_prescribe_albedo': False,    ###
        'do_thermodynamic_albedo': True,  ###
        'sea_ice_in_mixed_layer':True,   ### Chung:active thermodynamic ice        
        ### fixed sea ice
        'do_thickness_lock':True,        ### Chung:active thermodynamic ice thickness lock
        'ice_thickness_file':'/home/pochung/Isca/input/input_Chung/Ice/Ice_MiMActl_pen.nc', ###
        'ice_thickness_field':'h_ice_MiMActl_pen', ###
         'do_fraction_lock':True,        ### Chung:active thermodynamic ice thickness lock
         'ice_fraction_file':'/home/pochung/Isca/input/input_Chung/Ice/Ice_MiMActl_pen.nc', ###
         'ice_fraction_field':'a_ice_MiMActl_pen', ###
        ### fixed sea ice
        'thermodynamic_albedo_ocn':0.22,  ### 0.1 --> 0.22 to make it colder
        'thermodynamic_albedo_ice':0.45,   ### 0.4 --> 0.45
        ### fixed SST
        'do_sc_SST': True,
        'sst_file_Chung': '/home/pochung/Isca/input/input_Chung/t_surf/t_surf_mimactl_pen.nc', ###
        'sst_field_Chung': 't_surf_mimactl_pen' ###
#        'sst_file': 't_surf_MiMActl_pen' ###
        ### fixed SST
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
        'do_simple':True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  50.,
        'do_conserve_energy': True,         
    },

    'qflux_nml': {
        'qflux_amp': 0.0 #30.0 !!! Chung:Is this a problem?
    },

    'rrtm_radiation_nml': {
        'solr_cnst': 1360,  #s set solar constant to 1360, rather than default of 1368.22
        'dt_rad': 7200, #Use long RRTM timestep
        'do_read_ozone':True,
#        'ozone_file':'ozone_1990',
        'ozone_file':'/home/pochung/Isca/input/rrtm_input_files/ozone_1990.nc',
        'ozone_field':'ozone_1990',
   #     'do_rad_time_avg':True,
        'dt_rad_avg':86400,  # must be an integer
        ### Chung: coalbedo
        'do_coalb_toa_P2': False ,
   #     'coalb_toa_P2_a0': 0.90,
   #     'del_coalb_toa_P2': 0.0
#        'co2ppmv':1200, # 4xCO2
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
    #####
    'spectral_init_cond_nml':{
        'initial_temperature':288, #
    }, #####

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':30,
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'even_sigma', ###
        'surf_res':0.5,
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    }
#,
#    'vert_coordinate_nml': {
#        'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
#        'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
#       }    
})

exp.set_resolution('T42')
#Lets do a run!
if __name__=="__main__":
    exp.run(731, use_restart=True , num_cores=NCORES ,overwrite_data=True)
    for i in range(732,841):
        exp.run(i, num_cores=NCORES ,overwrite_data=True)
