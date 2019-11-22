#!!!This script is for T42 resolution

import numpy as np

from isca import IscaCodeBase, GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress
from ntfy import notify

NCORES = 16 #8 cores for T21, 16 for T42
RESOLUTION = 'T42', 25 #T42 horizontal resolution, 25 levels in pressure

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = GreyCodeBase.from_directory(GFDL_BASE)
# cb = GreyCodeBase.from_repo(repo='git@github.com:sit23/Isca.git', commit='7145be9')

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create a diagnostics output file for daily snapshots
diag = DiagTable()
#for atmos_monthly: 30 'days' in a 'month' multiplied by length of day in seconds
diag.add_file('atmos_monthly', 30*1379678, 'seconds', time_units='days') 
#for atmos_daily: length of 'day' in seconds
#diag.add_file('atmos_daily', 1379678 , 'seconds', time_units='days')
#diag.add_file('atmos_e_daily',1 , 'days', time_units='days')
#diag.add_file('atmos_timestep', 240, 'seconds', time_units='days')

# Define diag table entries
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk', time_avg=True)
diag.add_field('dynamics', 'pk', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)

diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('atmosphere', 'precipitation', time_avg=True)

diag.add_field('atmosphere', 'dt_tg_convection', time_avg=True)

#below checks the momentum tendency being applied by the sponge layer
#should be non-zero over the top 3 model levels
diag.add_field('damping', 'udt_rdamp', time_avg=True)

diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)

diag.add_field('two_stream', 'mars_solar_long', time_avg=True)
diag.add_field('two_stream', 'coszen', time_avg=True)
diag.add_field('two_stream', 'rrsun', time_avg=True)
diag.add_field('two_stream', 'swdn_toa', time_avg=True)
diag.add_field('two_stream', 'time_since_ae', time_avg=True)
diag.add_field('two_stream', 'true_anomaly', time_avg=True)
diag.add_field('two_stream', 'dec', time_avg=True)
diag.add_field('two_stream', 'ang', time_avg=True)

# define namelist values as python dictionary
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 178, #time step in seconds, multiple of length of day
        'days': 0.,
        'seconds': 30.*1379678, #length of time in each .nc file
        'calendar': 'no_calendar'
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':0.,                
        'two_stream_gray': True,     #if 'True' use grey radiation scheme
        'do_lscale_cond': False,
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
        'use_actual_surface_temperatures':False,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'mixed_layer_nml': {
        'tconst' : 120.,
        'prescribe_initial_dist':True,
        'evaporation':False,   
        'albedo_value': 0.3,
    },

    'qe_moist_convection_nml': {
        'rhbm':0.0,
        'tau_bm':3600.,
    },

    'dry_convection_nml': {
        'tau':7200,
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
        'do_simple':True,
        'tcmin':  -223, #Make sure low temperature limit of saturation vapour pressure is low enough that it doesn't cause an error (note that this giant planet has no moisture anyway, so doesn't directly affect calculation.        
        'tcmax': 350.,
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.125,              # neg. value: time in *days*
        'sponge_pbottom':  1.1205,           #Bottom of the model's sponge down to 1.1205 for 3 model layers for Titan (units are Pa)
    },

    'spectral_dynamics_nml': {
        'num_levels': 25,
        'exponent': 2.5,
        'scale_heights': 4,
        'surf_res': 0.1,
        'robert_coeff': 4e-2,
        'do_water_correction': False,
        'vert_coord_option': 'input',
        'initial_sphum': 0.,
        'valid_range_T': [0, 700]
    },

    'vert_coordinate_nml': {
        'bk': [  0.00000000e+00,   1.53008955e-04,   4.63790800e-04,   1.10977640e-03,   2.37044150e-03,   4.74479200e-03,         9.12245300e-03,   1.70677050e-02,   3.12516100e-02,         5.59939500e-02,   9.76165000e-02,   1.63754900e-01,         2.60315150e-01,   3.85974250e-01,   5.28084800e-01,         6.65956600e-01,   7.81088000e-01,   8.65400050e-01,         9.21109250e-01,   9.55343500e-01,   9.75416950e-01,         9.86856800e-01,   9.93269300e-01,   9.96830200e-01,         9.98799150e-01,   1.00000000e+00], 
        'pk': [0.]*26,
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use, which in this case is Frierson
        'do_seasonal': True,
        'atm_abs': 30,
        'sw_diff':0.0,
        'ir_tau_eq':25,        
        'ir_tau_pole':25,
        'linear_tau': 1.0,
        'equinox_day':0.0,
        'use_time_average_coszen':True,
        'solar_constant':15.08, 
        'solar_exponent':0.65    #!!!line additional to Mars script
    },


#     configure the relaxation profile
#     'hs_forcing_nml': {
#         'equilibrium_t_option' : 'top_down',
#         'ml_depth': 10.,
#         'spinup_time': 108000,
#         'ka': -2.,
#         'ks': -2.,
#         'tau_s': 0.2,
#         'calculate_insolation_from_orbit' : True,
#         'do_rayleigh_damping':False,
#         'albedo':0.25,
#         'pure_rad_equil':True,
#         'stratosphere_t_option':'pure_rad_equil',
#         'peri_time': 0.,
#         'h_a': 10.8,
#         'use_olr_from_t_surf':True,
#         'frac_of_year_ae':0.3032894472101727,
#     },


    'spectral_init_cond_nml': { #namelist additional to the Mars script
        'initial_temperature': 200.
    },

    'astronomy_nml': { 
        'ecc':0.054, #eccentricity of Saturn's orbit around the Sun
        'obliq':26.7, #obliquity wrt the Sun
        'use_mean_anom_in_rrsun_calc':True,
        'use_old_r_inv_squared':False
    },

    'constants_nml': {
        'orbital_period': 928523294,  #value calculated from day/year calculator script
        'solar_const':15.08, 
        'radius':2575.0e3,
        'rdgas':290,
        'kappa':0.2727,
        'rotation_period':1377631, #value calculated from day/year calculator script
    },

})


if __name__=="__main__":

    conv_schemes = ['none']

    depths = [2.]

    pers = [93] #longitude of perihelion

    scale = 1.

    for conv in conv_schemes:
        for depth_val in depths:
            for per_value in pers:
                exp = Experiment('grey_titan_T42', codebase=cb) #name of folder in which .nc files are output
                exp.clear_rundir()

                exp.diag_table = diag
                exp.namelist = namelist.copy()
                exp.namelist['constants_nml']['grav']     = scale * 1.354 #surface gravity
                exp.namelist['constants_nml']['pstd']     = scale * 14670000.0 #surface pressure - 100 times bigger than below values
                exp.namelist['constants_nml']['pstd_mks'] = scale * 146700.0 #in Pa
                exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = scale * 146700.0 #in Pa
                exp.namelist['idealized_moist_phys_nml']['convection_scheme'] = conv
                exp.namelist['mixed_layer_nml']['depth'] = depth_val
                exp.namelist['astronomy_nml']['per'] = per_value
                exp.set_resolution(*RESOLUTION)

#            with exp_progress(exp, description='o%.0f d{day}' % scale):
                exp.run(1, use_restart=False, num_cores=NCORES)
                for i in range(2, 111):
#                with exp_progress(exp, description='o%.0f d{day}' % scale):
                    exp.run(i, num_cores=NCORES)
                notify('top down with conv scheme = '+conv+' has completed', 'isca')
