"""
Non-orographic (convective) gravity wave drag test case, using the cg_drag
scheme -- originally implemented by Chaim Garfinkel (based on Alexander &
Dunkerton, JAS 1999), with edits by Stephen Thomson and further development
by Ross Castle (see src/atmos_param/cg_drag/cg_drag.f90 for the parameter
history).

Physics settings (do_rayleigh off, do_cg_drag on, tuned cg_drag_nml) are the
QBO-capable configuration: do_rayleigh is off so the two forcings don't fight
each other. NUM_MONTHS is a short comparison run -- long enough for the
cg_drag forcing to develop and show up in the zonal-mean wind, well short of
the multi-year integration a full QBO cycle would need.
"""
import os

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 16
RESOLUTION = 'T42', 50
NUM_MONTHS = 6

cb = IscaCodeBase.from_directory(GFDL_BASE)

exp = Experiment('cg_drag_qbo_test', codebase=cb)
exp.clear_rundir()

exp.inputfiles = [os.path.join(GFDL_BASE, 'input/rrtm_input_files/ozone_1990.nc')]

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('rrtm_radiation', 'co2', time_avg=True)
diag.add_field('damping', 'udt_cgwd', time_avg=True)  # cg_drag zonal wind tendency
exp.diag_table = diag

exp.namelist = namelist = Namelist({
    'main_nml': {
        'days': 30,
        'hours': 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos': 120,
        'current_date': [1, 1, 1, 0, 0, 0],
        'calendar': 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb': True,
        'mixed_layer_bc': True,
        'do_virtual': False,
        'do_simple': True,
        'roughness_mom': 3.21e-05,
        'roughness_heat': 3.21e-05,
        'roughness_moist': 3.21e-05,
        'two_stream_gray': False,       # Use RRTM, not grey radiation
        'do_rrtm_radiation': True,
        'convection_scheme': 'FULL_BETTS_MILLER'
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,
        'do_diffusivity': True,
        'do_simple': True,
        'constant_gust': 0.0,
        'use_tau': False
    },

    'diffusivity_nml': {
        'do_entrain': False,
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

    'mixed_layer_nml': {
        'tconst': 285.,
        'prescribe_initial_dist': True,
        'evaporation': True,
        'depth': 100.,
        'albedo_value': 0.25,
        'do_qflux': False,
    },

    'betts_miller_nml': {
        'rhbm': .7,
        'do_simp': False,
        'do_shallower': True,
    },

    'lscale_cond_nml': {
        'do_simple': True,
        'do_evap': True
    },

    'sat_vapor_pres_nml': {
        'do_simple': True
    },

    'damping_driver_nml': {
        'do_rayleigh': False,       # off, so cg_drag is the only source of GWD forcing
        'trayfric': -0.5,
        'sponge_pbottom': 50.,
        'do_conserve_energy': True,
        'do_mg_drag': False,
        'do_cg_drag': True,
    },

    'cg_drag_nml': {
        'Bt_0': 0.0043,
        'Bt_nh': 0.0014,
        'Bt_eq': 0.0043,
        'Bt_sh': 0.0014,
        'phi0n': 15.,
        'phi0s': -15.,
        'dphin': 10.,
        'dphis': -10.,
        'flag': 1,
        'Bw': 0.4,
        'Bn': 8.4,
        'cw': 35.0,
        'cwtropics': 35.0,
        'cn': 2.0,
        'kelvin_kludge': 1.0,
        'weighttop': 0.7,
        'weightminus1': 0.28,
        'weightminus2': 0.02,
        'source_level_pressure': 315.e+02,
        'damp_level_pressure': 0.85e+02,
        'cg_drag_freq': 21600
    },

    'rrtm_radiation_nml': {
        'do_read_ozone': True,
        'ozone_file': 'ozone_1990',
        'solr_cnst': 1360.,
        'dt_rad': 3600,
        'do_read_co2': False,
        'co2ppmv': 300,
    },

    'diag_manager_nml': {
        'mix_snapshot_average_fields': False
    },

    'fms_nml': {
        'domains_stack_size': 600000
    },

    'fms_io_nml': {
        'threading_write': 'single',
        'fileset_write': 'single',
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,
        'water_correction_limit': 200.e2,
        'reference_sea_level_press': 1.0e5,
        'num_levels': 50,
        'valid_range_t': [100., 800.],
        'initial_sphum': [2.e-6],
        'vert_coord_option': 'uneven_sigma',
        'surf_res': 0.1,
        'scale_heights': 7.9,
        'exponent': 1.4,
        'robert_coeff': 0.03,
    }
})

exp.set_resolution(*RESOLUTION)

if __name__ == '__main__':
    cb.compile()
    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2, NUM_MONTHS + 1):
        exp.run(i, num_cores=NCORES)
