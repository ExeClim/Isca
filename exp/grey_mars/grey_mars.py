import numpy as np

from isca import IscaCodeBase, GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress
from ntfy import notify

NCORES = 16

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
#diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_file('atmos_daily', 1, 'days', time_units='days')

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

diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)

diag.add_field('two_stream', 'mars_solar_long', time_avg=True)
diag.add_field('two_stream', 'coszen', time_avg=True)
diag.add_field('two_stream', 'rrsun', time_avg=True)
diag.add_field('two_stream', 'swdn_toa', time_avg=True)
diag.add_field('two_stream', 'time_since_ae', time_avg=True)
diag.add_field('two_stream', 'true_anomaly', time_avg=True)


# define namelist values as python dictionary
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 150,
        'days': 30.,
        'seconds': 0,
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
        'two_stream_gray': True,     #Use grey radiation
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
        'rh_target': 50.,
        'delta_t_relax': 7200.,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':False,   
        'albedo_value': 0.7,
        'depth':10.,
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
        'do_simple':True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.25,              # neg. value: time in *days*
        'sponge_pbottom':  0.5,           #Bottom of the model's sponge down to 50hPa (units are Pa)
        'do_conserve_energy': True,             
    },

    'spectral_dynamics_nml': {
        'num_levels': 30,
        'exponent': 2.5,
        'scale_heights': 4,
        'surf_res': 0.1,
        'robert_coeff': 4e-2,
        'do_water_correction': False,
        'vert_coord_option': 'even_sigma',
        'initial_sphum': 0.,
        'valid_range_T': [0, 700]
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use, which in this case is Frierson
        'do_seasonal': True,
        'atm_abs': 0.2,
        'sw_diff':0.0,
        'ir_tau_eq':0.2,        
        'ir_tau_pole':0.2,
        'linear_tau': 1.0,
        'equinox_day':0.8032894472101727,
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

    'astronomy_nml': {
        'ecc':0.0935,
        'obliq':25.19,
        'per':109.18420099566217, #to be equivalent to peri_time = 0. in hs_forcing version of mars
        'use_mean_anom_in_rrsun_calc':False,
        'use_old_r_inv_squared':False
    },

    'constants_nml': {
        'orbital_period': 686.980*86400.,
        'solar_const':589.0,
        'radius':3396.0e3,
        'omega':7.088e-5,
        'rdgas':192.0,
        'kappa':0.22727,
    },

})

if __name__=="__main__":

    conv_schemes = ['dry', 'none']

    scales = [ 1.]

    for conv in conv_schemes:
        for scale in scales:
            exp = Experiment('grey_mars_mk16_fixed_year_length_shifted_thin_10m_conv_outputs_shiny_'+conv, codebase=cb)
            exp.clear_rundir()

            exp.diag_table = diag
            exp.namelist = namelist.copy()
            exp.namelist['constants_nml']['grav']     = scale * 3.71
            exp.namelist['constants_nml']['pstd']     = scale * 6100.0
            exp.namelist['constants_nml']['pstd_mks'] = scale * 610.0
            exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = scale * 610.0
            exp.namelist['idealized_moist_phys_nml']['convection_scheme'] = conv

            with exp_progress(exp, description='o%.0f d{day}' % scale):
                exp.run(1, use_restart=False, num_cores=NCORES)
            for i in range(2, 121):
                with exp_progress(exp, description='o%.0f d{day}' % scale):
                    exp.run(i, num_cores=NCORES)
            notify('top down with conv scheme = '+conv+' has completed', 'isca')
