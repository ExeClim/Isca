import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE, GFDL_DATA

NCORES = 16
base_dir = os.getcwd()
# a CodeBase can be a directory on the computer,
# useful for iterative development
# cb = IscaCodeBase.from_directory(GFDL_BASE)


cb = IscaCodeBase.from_repo(repo='https://github.com/mp586/Isca.git', commit='9e19192') #looked up commit which was used for original 2xCO2 experiment

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
exp = Experiment('aquaplanet_finalIscaAPqflux_zerointegral', codebase=cb)



#Add any input files that are necessary for a particular experiment.
exp.inputfiles = [os.path.join(GFDL_BASE,'input/aquaplanet/land.nc'),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),os.path.join(GFDL_BASE,'input/aquaplanet/isca_qflux/zero_integral/ocean_qflux.nc')]
#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days') # to output data in two different time resolutions
# diag.add_file('atmos_6_hourly', 6, 'hours', time_units='hours') 

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True) # don't specify files, then it writes to both
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
diag.add_field('dynamics', 'omega', time_avg=True) 
diag.add_field('dynamics', 'height', time_avg=True) # geopotential height 
diag.add_field('atmosphere', 'rh', time_avg=True) 
diag.add_field('dynamics', 'slp', time_avg=True) # sea level pressure
diag.add_field('dynamics', 'zsurf', time_avg=True) # geopotential height at surface
diag.add_field('rrtm_radiation', 'toa_sw',time_avg=True)

diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True) # net SW surface flux
diag.add_field('rrtm_radiation', 'flux_lw', time_avg=True) # net LW surface flux
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True) # latent heat flux (up) at surface
diag.add_field('mixed_layer', 'flux_t', time_avg=True) # sensible heat flux (up) at surface
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True)

# moisture flux 
diag.add_field('dynamics', 'sphum_u', time_avg=True)
diag.add_field('dynamics', 'sphum_v', time_avg=True)
diag.add_field('dynamics', 'sphum_w', time_avg=True)

#MP added on 11 october 2017
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
        'dt_atmos':720,
        'current_date' : [0001,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':2.e-4, #Ocean roughness lengths
        'roughness_heat':2.e-4,
        'roughness_moist':2.e-4,      
        'two_stream_gray':False, #Don't use grey radiation
        'do_rrtm_radiation':True, #Do use RRTM radiation
        'convection_scheme':'SIMPLE_BETTS_MILLER', #Use the simple betts-miller convection scheme
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
        'depth':20., #Mixed layer depth
        'albedo_value': 0.25, #Ocean albedo value
        'do_qflux' : False, #Do not use prescribed qflux formula
        'load_qflux': True,
        'qflux_file_name':'ocean_qflux',
        'time_varying_qflux': True # shouldn't this be false? what is 
        # a non time varying qflux?
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
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,         
    },

    'rrtm_radiation_nml': {
        'do_read_ozone':True,
        'ozone_file':'ozone_1990',
        'solr_cnst' : 1360., #s set solar constant to 1360, rather than default of 1368.22
        'dt_rad': 3600, #Set RRTM radiation timestep to 3600 seconds, meaning it runs every 5 atmospheric timesteps        
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
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    }
})

#Lets do a run!
exp.run(1, use_restart=False, num_cores=NCORES)
for i in range(2,481):
    exp.run(i, num_cores=NCORES)








####### co2 experiment #######
exp = Experiment('aquaplanet_finalIscaAPqflux_zerointegral_2xCO2_spinup_361', codebase=cb)



#Add any input files that are necessary for a particular experiment.
exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),os.path.join(GFDL_BASE,'input/aquaplanet/isca_qflux/zero_integral/ocean_qflux.nc'), os.path.join(GFDL_BASE,'input/co2_doubling.nc')]
#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
# diag.add_file('atmos_6_hourly', 6, 'hours', time_units='hours') 

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
diag.add_field('dynamics', 'omega', time_avg=True) 
diag.add_field('dynamics', 'height', time_avg=True) # geopotential height 
diag.add_field('atmosphere', 'rh', time_avg=True) 
diag.add_field('dynamics', 'slp', time_avg=True) # sea level pressure
diag.add_field('dynamics', 'zsurf', time_avg=True) # geopotential height at surface
diag.add_field('rrtm_radiation', 'toa_sw',time_avg=True)

diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True) # net SW surface flux
diag.add_field('rrtm_radiation', 'flux_lw', time_avg=True) # net LW surface flux
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True) # latent heat flux (up) at surface
diag.add_field('mixed_layer', 'flux_t', time_avg=True) # sensible heat flux (up) at surface
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True)
diag.add_field('rrtm_radiation', 'co2', time_avg=True)
# moisture flux
diag.add_field('dynamics', 'sphum_u', time_avg=True)
diag.add_field('dynamics', 'sphum_v', time_avg=True)
diag.add_field('dynamics', 'sphum_w', time_avg=True)

#MP added on 11 october 2017
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
        'dt_atmos':720,
        'current_date' : [0001,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':2.e-4, #Ocean roughness lengths
        'roughness_heat':2.e-4,
        'roughness_moist':2.e-4,      
        'two_stream_gray':False, #Don't use grey radiation
        'do_rrtm_radiation':True, #Do use RRTM radiation
        'convection_scheme':'SIMPLE_BETTS_MILLER', #Use the simple betts-miller convection scheme
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
        'depth':20., #Mixed layer depth
        'albedo_value': 0.25, #Ocean albedo value
        'do_qflux' : False, #Do not use prescribed qflux formula
        'load_qflux': True,
        'qflux_file_name':'ocean_qflux',
        'time_varying_qflux': True # shouldn't this be false? what is 
        # a non time varying qflux?
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
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,         
    },

    'rrtm_radiation_nml': {
        'do_read_ozone':True,
        'ozone_file':'ozone_1990',
        'do_read_co2': True,
        'co2_file': 'co2_doubling',
        'solr_cnst' : 1360., #s set solar constant to 1360, rather than default of 1368.22
        'dt_rad': 3600, #Set RRTM radiation timestep to 3600 seconds, meaning it runs every 5 atmospheric timesteps        
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
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    }
})

#Lets do a run!
exp.run(1, restart_file='/scratch/mp586/Isca_DATA/aquaplanet_finalIscaAPqflux_zerointegral/restarts/res0361.tar.gz', num_cores=NCORES)
for i in range(2,481):
    exp.run(i, num_cores=NCORES)
