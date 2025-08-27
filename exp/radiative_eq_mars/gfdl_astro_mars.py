import numpy as np

from isca import IscaCodeBase, GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress
from ntfy import notify

NCORES = 16

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = GreyCodeBase.from_directory(GFDL_BASE)
#cb = GreyCodeBase.from_repo(repo='git@github.com:sit23/Isca.git', commit='fa26ea1')

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
diag.add_file('atmos_daily',88440 , 'seconds', time_units='days')

# Define diag table entries
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk', time_avg=True)
diag.add_field('dynamics', 'pk', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'ucomp_vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'slp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
# diag.add_field('dynamics', 'sphum', time_avg=True)
# diag.add_field('atmosphere', 'precipitation', time_avg=True)
# diag.add_field('atmosphere', 'rh', time_avg=True)

diag.add_field('mixed_layer', 't_surf', time_avg=True)
#diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)

diag.add_field('hs_forcing', 'teq', time_avg=True)
diag.add_field('hs_forcing', 'mars_solar_long', time_avg=True)
diag.add_field('hs_forcing', 'true_anom', time_avg=True)
diag.add_field('hs_forcing', 'incoming_sw', time_avg=True)
diag.add_field('hs_forcing', 'dec', time_avg=True)
diag.add_field('hs_forcing', 'rrsun', time_avg=True)
diag.add_field('hs_forcing', 'ang', time_avg=True)

# define namelist values as python dictionary
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 110,
        'days': 0.,
        'seconds': 30.*88440.,
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
        'two_stream_gray': False,     #Use grey radiation
        'do_newtonian_cooling_as_rad': True,     #Use grey radiation        
        'convection_scheme': 'none', #Use the simple Betts Miller convection scheme from Frierson
        'newt_relax_surface': False,
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
        'rh_target': 50.,
        'delta_t_relax': 7200.,
        'use_actual_surface_temperatures':False,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':False,   
        'update_sst_from_fluxes' : False,
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


    # configure the relaxation profile
    'hs_forcing_nml': {
        'equilibrium_t_option' : 'top_down',
        'ml_depth': 2.,
        'spinup_time': 108000,
        'ka': -40.,
        'ks': -2.,
        'tau_s': 0.2,
        'calculate_insolation_from_orbit' : True,
        'do_rayleigh_damping':False,
        'albedo':0.3,
        'pure_rad_equil':True,
        'stratosphere_t_option':'pure_rad_equil',
        'h_a': 10.8,
        'use_olr_from_t_surf':True,
        'equinox_day': 0.,
        'use_gfdl_astronomy':True,
    },

    'astronomy_nml': {
        'ecc':0.0935,
        'obliq':25.19,
        'per':70.85,        
        'use_mean_anom_in_rrsun_calc':True,
        'use_old_r_inv_squared':False        
    },

    'constants_nml': {
        'orbital_period': 59166360, 
        'solar_const':589.0,
        'radius':3396.0e3,
        'rdgas':192.0,
        'kappa':0.22727,
        'rotation_period':88308,
    },

})

if __name__=="__main__":

    conv_schemes = ['none']

    scales = [ 1.]

    for conv in conv_schemes:
        for scale in scales:
            exp = Experiment('radiaitive_eq_mars_40_rad_ka_gfdl_astro_macda_levels_no_conv_scheme', codebase=cb)
            exp.clear_rundir()

            exp.diag_table = diag
            exp.namelist = namelist.copy()
            exp.namelist['constants_nml']['grav']     = scale * 3.71
            exp.namelist['constants_nml']['pstd']     = scale * 6100.0
            exp.namelist['constants_nml']['pstd_mks'] = scale * 610.0
            exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = scale * 610.0
            exp.namelist['idealized_moist_phys_nml']['convection_scheme'] = conv

            notify('top down with conv scheme = '+conv+' has started', 'isca')

            exp.run(1, use_restart=False, num_cores=NCORES)
            for i in range(2, 361):
                exp.run(i, num_cores=NCORES)
            notify('top down with conv scheme = '+conv+' has completed', 'isca')
