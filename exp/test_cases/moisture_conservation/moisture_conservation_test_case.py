import os

import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 4

# Purpose: an aquaplanet-style world with land everywhere (land_option='all_land',
# no input land mask needed) and bucket hydrology active at every gridpoint, with
# a generous, non-limiting initial bucket depth. Unlike a standard mixed-layer
# aquaplanet - where the ocean is an effectively infinite moisture source/sink -
# every gridpoint here has a finite water reservoir, so the total column water
# (atmosphere + bucket) should be exactly conserved over time in the absence of
# numerical bugs. Used to test whether gcm_vert_diff_down/up's moisture-tendency
# argument aliasing (dt_q vs dt_tracers(:,:,:,nsphum)) causes measurable
# non-conservation.

cb = GreyCodeBase.from_directory(GFDL_BASE)

exp = Experiment('moisture_conservation_test_experiment', codebase=cb)

diag = DiagTable()
diag.add_file('atmos_daily', 1, 'days', time_units='days')

diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'bucket_depth', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)
diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)

exp.diag_table = diag

exp.clear_rundir()

exp.namelist = namelist = Namelist({
    'main_nml': {
        'days': 1,
        'hours': 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos': 600,
        'current_date': [1, 1, 1, 0, 0, 0],
        'calendar': 'thirty_day',
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb': True,
        'mixed_layer_bc': True,
        'do_virtual': False,
        'do_simple': True,
        'roughness_mom': 2.e-4,
        'roughness_heat': 2.e-4,
        'roughness_moist': 2.e-4,
        'two_stream_gray': True,
        'do_rrtm_radiation': False,
        'convection_scheme': 'SIMPLE_BETTS_MILLER',
        'land_option': 'all_land',      # land everywhere, no input file needed
        'bucket': True,
        'init_bucket_depth_land': 1.,   # 1m of water everywhere at t=0 - generous, non-limiting
        'max_bucket_depth_land': 1000., # effectively unbounded, so no overflow/runoff sink during this test
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,
        'do_diffusivity': True,
        'do_simple': True,
        'constant_gust': 0.0,
        'use_tau': False,
    },

    'diffusivity_nml': {
        'do_entrain': False,
        'do_simple': True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True,
    },

    'mixed_layer_nml': {
        'tconst': 285.,
        'prescribe_initial_dist': True,
        'evaporation': True,
        'depth': 5.,                     # shallow - fast surface temperature response, not the focus here
        'land_h_capacity_prefactor': 0.1,
        'albedo_value': 0.25,
        'do_qflux': False,
    },

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
    },

    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,
        'sponge_pbottom': 150.,
        'do_conserve_energy': True,
    },

    'diag_manager_nml': {
        'mix_snapshot_average_fields': False,
    },

    'fms_nml': {
        'domains_stack_size': 600000,
    },

    'fms_io_nml': {
        'threading_write': 'single',
        'fileset_write': 'single',
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,
        'water_correction_limit': 200.e2,
        'do_water_correction': False,   # isolate the diffusion pathway - no separate mass-fixer masking it
        'reference_sea_level_press': 1.0e5,
        'num_levels': 25,
        'valid_range_t': [100., 800.],
        'initial_sphum': [2.e-6],
        'vert_coord_option': 'uneven_sigma',
        'surf_res': 0.2,
        'scale_heights': 11.0,
        'exponent': 7.0,
        'robert_coeff': 0.03,
    },
})

if __name__ == "__main__":
    cb.compile()
    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2, 16):
        exp.run(i, num_cores=NCORES)
