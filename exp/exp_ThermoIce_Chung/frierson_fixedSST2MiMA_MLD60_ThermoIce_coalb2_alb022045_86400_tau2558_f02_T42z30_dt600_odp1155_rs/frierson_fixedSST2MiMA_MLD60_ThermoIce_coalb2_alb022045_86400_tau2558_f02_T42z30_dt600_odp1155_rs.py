import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32
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
exp = Experiment('frierson_fixedSST2MiMA_MLD60_ThermoIce_coalb2_alb022045_86400_tau2558_f02_T42z30_dt600_odp1155_rs', codebase=cb)

#exp.inputfiles = [os.path.join(GFDL_BASE,'input/input_Chung/t_surf/t_surf_frictl_pen.nc')]

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
diag.add_field('dynamics', 'height', time_avg=True) #
diag.add_field('dynamics', 'height_half', time_avg=True) #
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
diag.add_field('dynamics', 'ucomp_height', time_avg=True) #
diag.add_field('dynamics', 'vcomp_height', time_avg=True) #
diag.add_field('dynamics', 'omega_height', time_avg=True) #
diag.add_field('dynamics', 'sphum_u', time_avg=True) ### 2022/11/02 #
diag.add_field('dynamics', 'sphum_v', time_avg=True) ### 2022/11/02 #
diag.add_field('dynamics', 'sphum_w', time_avg=True) ### 2022/11/02 #
diag.add_field('dynamics', 'EKE', time_avg=True) ### 2022/10/27
diag.add_field('dynamics', 'temp_sq', time_avg=True)
diag.add_field('dynamics', 'vcomp_vor', time_avg=True)

diag.add_field('mixed_layer', 'albedo', time_avg=True)
diag.add_field('mixed_layer', 'ml_heat_cap', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True) #
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_u', time_avg=True) ##
diag.add_field('mixed_layer', 'h_ice', time_avg=True) ##
diag.add_field('mixed_layer', 'a_ice', time_avg=True) ##
diag.add_field('mixed_layer', 't_ml', time_avg=True) ##
diag.add_field('mixed_layer', 'flux_ice', time_avg=True) ##

diag.add_field('two_stream', 'flux_sw', time_avg=True) #
diag.add_field('two_stream', 'flux_lw', time_avg=True) #
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)
diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
diag.add_field('two_stream', 'tdt_rad', time_avg=True)
diag.add_field('vert_turb', 'z_pbl', time_avg=True)
#
diag.add_field('atmosphere', 'temp_2m', time_avg=True) #
diag.add_field('atmosphere', 'u_10m', time_avg=True) #
diag.add_field('atmosphere', 'v_10m', time_avg=True) #
diag.add_field('atmosphere', 'q_2m', time_avg=True) #
diag.add_field('atmosphere', 'rh_2m', time_avg=True) #
#
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True) #####
#
diag.add_field('two_stream', 'coszen', time_avg=True) ### output for offfline kernels
#
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
     'dt_atmos':600,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
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
        'two_stream_gray': True,     #Use grey radiation
        'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use the simple Betts Miller convection scheme from Frierson
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
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,   
        'depth': 60,                          #Depth of mixed layer used
#        'albedo_value': 0.10,                  #Albedo value used            
        'do_prescribe_albedo': False,    ###
        'do_thermodynamic_albedo': True,  ###
        'sea_ice_in_mixed_layer':True,   ### Chung:active thermodynamic ice
        ### fixed sea ice
        'do_thickness_lock':True,        ### Chung:active thermodynamic ice thickness lock
        'ice_thickness_file':'/home/pochung/Isca/input/input_Chung/Ice/Ice_MiMActl_pen.nc', ### ###
        'ice_thickness_field':'h_ice_MiMActl_pen', ###
         'do_fraction_lock':True,        ### Chung:active thermodynamic ice thickness lock
         'ice_fraction_file':'/home/pochung/Isca/input/input_Chung/Ice/Ice_MiMActl_pen.nc', ### ### 2022/09/02
         'ice_fraction_field':'a_ice_MiMActl_pen', ###
        ### fixed sea ice
#        ### fixed input TS for ice
#        'sfc_melt_from_file':True,
#        'sfc_melt_file':'/home/pochung/Isca/input/input_Chung/t_surf/t_surf_frictl_pen.nc', ###
#        'sfc_melt_field':'t_surf_frictl_pen', ###
        ### fixed input TS for ice
        'thermodynamic_albedo_ocn':0.22,   ### 0.1 --> 0.22 to make it colder
        'thermodynamic_albedo_ice':0.45,   ### 0.4 --> 0.45
        ### fixed SST
        'do_sc_SST': True,
        'sst_file_Chung': '/home/pochung/Isca/input/input_Chung/t_surf/t_surf_mimactl_pen.nc', ###
        'sst_field_Chung': 't_surf_mimactl_pen' ###
        ### fixed SST
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },

    'betts_miller_nml': {
       'rhbm': .7   , 
       'do_simp': False, 
       'do_shallower': True, 
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
        'trayfric': -0.25,              # neg. value: time in *days*
        'sponge_pbottom':  5000.,           #Bottom of the model's sponge down to 50hPa (units are Pa)
        'do_conserve_energy': True,             
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use, which in this case is Frierson
        'do_seasonal': True,                #do_seasonal=false uses the p2 insolation profile from Frierson 2006. do_seasonal=True uses the GFDL astronomy module to calculate seasonally-varying insolation.
        'atm_abs': 0.0,                      # default: 0.0        
        'use_time_average_coszen': True, # Chung: If true, average coszen over the period dt_rad_avg
        'dt_rad_avg': 86400.0,            # Chung: units: seconds       --> with seasonal cycle but no diurnal cycle in this case
        'ir_tau_eq':  5.8 ,  # Chung: default (6.0) --> 6.6
        'ir_tau_pole':2.5 ,  # Chung: change the temp gradient (default:1.5) --> 3.3
        ### Chung: coalbedo
        'do_coalb_toa_P2': True   ,     ###
        'coalb_toa_P2_a0': 0.7535    ,
        'del_coalb_toa_P2':-0.0345  ,
        'linear_tau' : 0.2 , ### vertical structure
        ###
        'odp' :1.155 , ### global warming <-- scaling factor to simulate 4xCO2
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
        'num_levels':30,               #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'even_sigma', #Use the vertical levels from Frierson 2006
        'surf_res':0.5,
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    },
#    'vert_coordinate_nml': {
#        'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
#        'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
#       }
})

exp.set_resolution('T42')
#Lets do a run!
if __name__=="__main__":
    exp.run(800, use_restart=True , num_cores=NCORES )
#    exp.run(1, use_restart=True , num_cores=NCORES ,restart_file ='/data/groups/feldl/pochung/Isca_out/frierson_fixedSST2MiMA_MLD60_ThermoIce_coalb2_alb022045_86400_tau2558_f02_T42z30_dt600_odp1155_rs/CTL_res0840.tar.gz')
    for i in range(801,841):
        exp.run(i, num_cores=NCORES )
