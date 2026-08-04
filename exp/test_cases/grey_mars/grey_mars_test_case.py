import os

import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 8

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = GreyCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = GreyCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('grey_mars_test_experiment', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_daily',88440 , 'seconds', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk', time_avg=True)
diag.add_field('dynamics', 'pk', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('two_stream', 'mars_solar_long', time_avg=True)
diag.add_field('two_stream', 'coszen', time_avg=True)
diag.add_field('two_stream', 'rrsun', time_avg=True)
diag.add_field('two_stream', 'swdn_toa', time_avg=True)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist. Physical constants (grav, pstd, radius, rdgas, etc.)
#are Mars' values, taken from exp/grey_mars/grey_mars.py. That script sweeps over convection
#scheme, mixed-layer depth and orbital argument of perihelion; here we fix one illustrative
#combination from that sweep ('dry' convection, 1.0 m mixed-layer depth, per=70.85 degrees).
exp.namelist = namelist = Namelist({
    'main_nml': {
        'dt_atmos': 110,
        'days': 0,
        'seconds': 30*88440,
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
        'do_lscale_cond': False,
        'convection_scheme': 'none',
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
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':False,
        'albedo_value': 0.3,
        'depth': 2.0,
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
        'tcmin':  -223, #Make sure low temperature limit of saturation vapour pressure is low enough that it doesn't cause an error (this planet has no moisture anyway, so doesn't directly affect the calculation).
        'tcmax': 350,
    },

    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.125,              # neg. value: time in *days*
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
        'valid_range_T': [0, 700],
        'reference_sea_level_press': 610.0,
    },

    'vert_coordinate_nml': {
        'bk': [  0.00000000e+00,   1.53008955e-04,   4.63790800e-04,   1.10977640e-03,   2.37044150e-03,   4.74479200e-03,         9.12245300e-03,   1.70677050e-02,   3.12516100e-02,         5.59939500e-02,   9.76165000e-02,   1.63754900e-01,         2.60315150e-01,   3.85974250e-01,   5.28084800e-01,         6.65956600e-01,   7.81088000e-01,   8.65400050e-01,         9.21109250e-01,   9.55343500e-01,   9.75416950e-01,         9.86856800e-01,   9.93269300e-01,   9.96830200e-01,         9.98799150e-01,   1.00000000e+00],
        'pk': [0.]*26,
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use, which in this case is Frierson
        'do_seasonal': True,
        'atm_abs': 0.2,
        'sw_diff':0.0,
        'ir_tau_eq':0.2,
        'ir_tau_pole':0.2,
        'linear_tau': 1.0,
        'equinox_day':0.0,
        'use_time_average_coszen':True,
        'solar_constant':589.0,
    },

    'astronomy_nml': {
        'ecc':0.0935,
        'obliq':25.19,
        'per': 70.85,
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
        'grav': 3.71,
        'pstd': 6100.0,
        'pstd_mks': 610.0,
    },
})

#Lets do a run!
if __name__=="__main__":
    cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2, 65):
        exp.run(i, num_cores=NCORES)
