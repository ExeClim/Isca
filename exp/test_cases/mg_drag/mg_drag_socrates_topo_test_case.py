"""
Orographic (mountain) gravity wave drag test case, using the mg_drag scheme
from Ross Castle's realistic_mgdrag branch, built on the official Socrates
aquaplanet-with-topography test case
(exp/test_cases/socrates_test/socrates_aquaplanet_amip_with_topo.py).

mg_drag_nml settings are Ross's own tuned values, except source_of_sgsmtn.
His own configuration used 'computed', which reads a raw high-resolution
elevation dataset via topography_nml and computes subgrid-topography
variance internally (see src/shared/topography/topography.F90); that raw
elevation dataset isn't included in this repo.

This test case instead uses source_of_sgsmtn='input', reading
input/mg_drag.res.nc (variable 'ghprime') directly -- a real subgrid-
orography-standard-deviation field regridded from ERA5 (see
src/extra/python/scripts/regrid_era5_sdor_to_t42.py for the regridding
method and source dataset).
"""
import os

from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 16
RESOLUTION = 'T42', 40
NUM_MONTHS = 6

base_dir = os.path.dirname(os.path.realpath(__file__))
mg_drag_ghprime_input = os.path.join(base_dir, 'input', 'mg_drag.res.nc')

cb = SocratesCodeBase.from_directory(GFDL_BASE)

exp = Experiment('mg_drag_socrates_topo_test', codebase=cb)
exp.clear_rundir()

exp.inputfiles = [
    os.path.join(GFDL_BASE, 'input/rrtm_input_files/ozone_1990.nc'),
    os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents/input/era-spectral7_T42_64x128.out.nc'),
    os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents/input/sst_clim_amip.nc'),
    os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents/input/siconc_clim_amip.nc'),
    mg_drag_ghprime_input,
]

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('socrates', 'soc_olr', time_avg=True)
# mg_drag diagnostics
diag.add_field('damping', 'udt_gwd', time_avg=True)
diag.add_field('damping', 'vdt_gwd', time_avg=True)
diag.add_field('damping', 'taubx', time_avg=True)
diag.add_field('damping', 'tauby', time_avg=True)
diag.add_field('damping', 'taus', time_avg=True)
diag.add_field('damping', 'tdt_diss_gwd', time_avg=True)
diag.add_field('damping', 'sgsmtn', time_avg=True)
exp.diag_table = diag

exp.namelist = namelist = Namelist({
    'main_nml': {
        'days': 30,
        'hours': 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos': 180,   # finer than the plain topo test case's 600s, needed for mg_drag stability
        'current_date': [1, 1, 1, 0, 0, 0],
        'calendar': 'thirty_day'
    },

    'socrates_rad_nml': {
        'stellar_constant': 1370.,
        'lw_spectral_filename': os.path.join(GFDL_BASE, 'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'),
        'sw_spectral_filename': os.path.join(GFDL_BASE, 'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'),
        'do_read_ozone': True,
        'ozone_file_name': 'ozone_1990',
        'ozone_field_name': 'ozone_1990',
        'dt_rad': 3600,
        'store_intermediate_rad': True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels': False,
        'tidally_locked': False,
        'solday': 90
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
        'two_stream_gray': False,
        'do_socrates_radiation': True,
        'convection_scheme': 'SIMPLE_BETTS_MILLER',
        'do_cloud_simple': False,
        'land_option': 'input',
        'land_file_name': 'INPUT/era-spectral7_T42_64x128.out.nc',
        'land_roughness_prefactor': 10.0,
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
        'old_dtaudv': True,
        'land_humidity_prefactor': 0.7,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'mixed_layer_nml': {
        'tconst': 285.,
        'prescribe_initial_dist': True,
        'evaporation': True,
        'land_option': 'input',
        'land_h_capacity_prefactor': 0.1,
        'albedo_value': 0.25,
        'land_albedo_prefactor': 1.3,
        'do_qflux': False,
        'do_read_sst': True,
        'do_sc_sst': True,
        'sst_file': 'sst_clim_amip',
        'specify_sst_over_ocean_only': True,
    },

    'qe_moist_convection_nml': {
        'rhbm': 0.7,
        'Tmin': 160.,
        'Tmax': 350.
    },

    'lscale_cond_nml': {
        'do_simple': True,
        'do_evap': True
    },

    'sat_vapor_pres_nml': {
        'do_simple': True,
        'construct_table_wrt_liq_and_ice': True
    },

    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,
        'sponge_pbottom': 150.,
        'do_conserve_energy': True,
        'do_cg_drag': False,
        'do_mg_drag': True,
    },

    'mg_drag_nml': {
        'do_netcdf_restart': True,
        'xl_mtn': 1.0e5,
        'gmax': 1.0,
        'acoef': 1.0,
        'rho': 1.13,
        'low_lev_frac': 0.23,
        'do_conserve_energy': True,
        'do_mcm_mg_drag': False,
        'source_of_sgsmtn': 'input',  # reads INPUT/mg_drag.res.nc ('ghprime'), ERA5-derived -- see module docstring
        'flux_cut_level': 0.0
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
        'num_levels': 40,
        'valid_range_t': [100., 800.],
        'initial_sphum': [2.e-6],
        'vert_coord_option': 'uneven_sigma',
        'surf_res': 0.2,
        'scale_heights': 11.0,
        'exponent': 7.0,
        'robert_coeff': 0.03,
        'ocean_topog_smoothing': 0.0
    },

    'spectral_init_cond_nml': {
        'topog_file_name': 'era-spectral7_T42_64x128.out.nc',
        'topography_option': 'input'
    },
})

exp.set_resolution(*RESOLUTION)

if __name__ == '__main__':
    cb.compile(debug=False)
    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2, NUM_MONTHS + 1):
        exp.run(i, num_cores=NCORES)
