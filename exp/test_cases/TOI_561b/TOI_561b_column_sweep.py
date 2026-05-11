"""
Single column experiment for TOI-561b: sweep over a range of lw and sw optical depths at Ps=3 bar.
Generates temperature-pressure profiles for each combination.
"""
import os
import numpy as np
from isca import ColumnCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 1

# Compile code
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = ColumnCodeBase.from_directory(GFDL_BASE)
cb.compile()

# Create experiment
exp = Experiment('TOI_561b_column_sweep', codebase=cb)

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_field('column', 'ps', time_avg=True)
diag.add_field('column', 'temp', time_avg=True)
diag.add_field('column', 'sphum', time_avg=True)
diag.add_field('column', 'press', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
exp.diag_table = diag

exp.clear_rundir()

# Set up base namelist for Ps=3 bar
exp.namelist = namelist = Namelist({
    'main_nml': {
        'days'   : 0,#360, # each run lasts one year, and then multiple runs are strung together below (loop on e.g. Line 310)
        'hours'  : 0,   # a different output file is produced for each run (in this case, each year). Data 
        'minutes': 4000,   # output at the frequency specified in the diag table 
        'seconds': 0,
        'dt_atmos':5,#10,#360,#480, # 900s timestep for dynamical core
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },
    'atmosphere_nml': {
        'idealized_moist_model': True
    },
    'column_nml': {
        'lon_max': 1,
        'lat_max': 1,
        'num_levels': 31,
        'initial_sphum': 1e-3,
        'q_decrease_only': True,
    },
    'column_grid_nml': {
        'lat_value': 0.0
    },
    'column_init_cond_nml': {
        'initial_temperature': 2000.,
        'surf_geopotential': 0.0,
        'surface_wind': 5.
    },
    'idealized_moist_phys_nml': {
        'two_stream_gray': True,
        'do_rrtm_radiation': False,
        'convection_scheme': 'DRYADJ',
        'do_damping': True,
        'turb': True,
        'mixed_layer_bc': True,
        'do_virtual': True,
        'roughness_mom': 5.e-3,
        'roughness_heat': 1,
        'roughness_moist': 1.e-5,
        'do_simple': False,
        'do_tj_bl': True,
    },
    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',
        'do_seasonal': False,
        'do_tl': True,
        'do_closein': True,
        'atm_abs': 0.,
        'solar_exponent': 1,
        'ir_tau_eq': 1,
        'ir_tau_pole': 1,
        'sw_tau_eq': 1,
        'sw_tau_pole': 1,
        'del_sol': 1.2,
        'solar_constant': 4715*1361,
        'R_stellar': 0.843*6.957e8,
        'd_stellar': 0.01055*1.496e11,
        'linear_tau': 1.,
        'odp': 1.0
    },
    'mixed_layer_nml': {
        'depth': 5,
        'albedo_value': 0,
        'prescribe_initial_dist': True,
        'tconst': 2000.,
        'delta_T': 0.,
        'evaporation': False,
        'do_qflux': False,
        'load_qflux': False,
    },
    'sat_vapor_pres_nml': {
        'do_simple': True,
        'tcmin_simple': -273,
        'tcmax_simple': 10000,
    },
    'vert_coordinate_nml': {
        'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
        'pk': [0.000000]*26,
    },
    'constants_nml': {
        'omega': 2*np.pi / (0.44656895*86400.),
        'es0': 1.e-6,
        'rdgas': 8.314/18.0153*1e3,
        'kappa': (8.314/18.0153*1e3)/ (2842 - (8.314/18.0153)*1e3),
        'radius': 1.4195*6.371e6,
        'grav': (6.67e-11*2.24*5.972e24)/((1.4195*6.371e6)**2),
        'pstd': 3e5,
        'pstd_mks': 3e6,
    },
})

# Sweep over a range of lw and sw optical depths

if __name__=="__main__":
    lw_tau_list = [1, 3, 5, 10, 15]
    sw_tau_list = [0.5, 1, 2, 3, 5]
    for lw_tau in lw_tau_list:
        for sw_tau in sw_tau_list:
            run_exp = exp.derive(f'TOI_561b/Column_sweep/TOI_561b_column_lw_{lw_tau}_sw_{sw_tau}_Ps_3bar')
            run_exp.namelist['two_stream_gray_rad_nml']['ir_tau_eq'] = lw_tau
            run_exp.namelist['two_stream_gray_rad_nml']['ir_tau_pole'] = lw_tau
            run_exp.namelist['two_stream_gray_rad_nml']['sw_tau_eq'] = sw_tau
            run_exp.namelist['two_stream_gray_rad_nml']['sw_tau_pole'] = sw_tau
            run_exp.namelist['constants_nml']['pstd'] = 3e5
            run_exp.namelist['constants_nml']['pstd_mks'] = 3e6
            run_exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] =3e5
            overwrite=False
            run_exp.run(1, use_restart=False, num_cores=NCORES)
            for i in range(2,11):
                run_exp.run(i,num_cores=NCORES, overwrite_data=overwrite)

