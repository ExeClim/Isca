import numpy as np
from gfdl.experiment import Experiment, Namelist, DiagTable

import sys

REPO = 'git@github.com:jamesp/GFDLmoistModel.git'

# create an experiment based on the code version tagged `exoplan0.3`
exp = Experiment('ref_earth_grey', commit='isca1.0', repo=REPO)

exp.disable_rrtm()
exp.compile()

nml = Namelist({})

nml['main_nml'] = {
    'dt_atmos': 600,
    'days': 30,
    'calendar': 'thirty_day',
    'current_date': [2000,1,1,0,0,0]
}

nml['spectral_dynamics_nml'] = {
    #'damping_option'          : 'resolution_dependent',
    'damping_order'           : 4,                      # default: 2
    #'do_mass_correction'      True
    #'do_energy_correction'    True
    #'do_water_correction'     True
    'water_correction_limit'  : 200.e2,                 # default: 0
    #'use_virtual_temperature' False
    #'vert_advect_uv'          : 'second_centered',
    #'vert_advect_t'           : 'second_centered',
    #'longitude_origin'        : 0.0,
    #'robert_coeff'            : .03,                   # default: 0.04
    #'alpha_implicit'          : .5,
    'reference_sea_level_press': 1.0e5,                  # default: 101325
    'valid_range_t'           : [100., 800.],           # default: (100, 500)
    #'initial_state_option'   : 'quiescent'
    'initial_sphum'           : 2.e-6,                  # default: 0
    'vert_coord_option'       : 'uneven_sigma',         # default: 'even_sigma'
    'scale_heights': 6.0,
    'exponent': 7.5,
    'surf_res': 0.5
}

nml['atmosphere_nml'] = {
    'idealized_moist_model': True
}

nml['diag_manager_nml'] = {
    'mix_snapshot_average_fields': False
}

nml['fms_nml'] = {
    'domains_stack_size': 600000                        # default: 0
}

nml['astronomy_nml'] = {

}

nml['fms_io_nml'] = {
    'threading_write': 'single',                         # default: multi
    'fileset_write': 'single',                           # default: multi
}

# from phys.nml
nml['idealized_moist_phys_nml'] = {
    #'two_stream_gray': True,
    #'do_rrtm_radiation': False,
    'convection_scheme': 'SIMPLE_BETTS_MILLER',
    'do_damping': True,
    'turb': True,
    'mixed_layer_bc': True,
    'do_virtual': False,
    'do_simple': True,
    # Roughness Lengths for Monin-Obukhov theory:
    # Baseline 3.21e-5 m
    # Ref:  Heng et al: Mon. Not. R. Astron. Soc [418] (2011)
    #       Frierson et al: J Atmos. Sci [63] (2006)
    # roughness_mom:
    #   Open water: 1e-4m
    #   Urban terrain: 1m
    'roughness_mom': 3.21e-05,             # default: 0.05
    'roughness_heat': 3.21e-05,            # default: 0.05
    'roughness_moist': 3.21e-05            # default: 0.05
}


nml['vert_turb_driver_nml'] = {
   'do_mellor_yamada': False,     # default: True
   'do_diffusivity': True,        # default: False
   'do_simple': True,             # default: False
   #'do_shallow_conv': False,
   #'gust_scheme': 'constant',
   'constant_gust': 0.0,          # default: 1.0
   #'use_tau': False
}

nml['diffusivity_nml'] = {
    'do_entrain': False,          # default: True
    'do_simple': True,            # default: False
    #'frac_inner': 0.1,
    #'rich_crit_pbl':  1.0
}

# nml['monin_obukhov_nml'] = {
#     'neutral': False,
#     'rich_crit': 2.0,
#     'stable_option': 1
# }

nml['surface_flux_nml'] = {
    'use_virtual_temp': False,
    'do_simple': True,
    'old_dtaudv': True
}

# nml['spectral_init_cond_nml'] = {
#     'initial_temperature': 264.0
# }

nml['two_stream_gray_rad_nml'] = {
    'rad_scheme': 'byrne',           # default: Byrne & O'Gorman
    'do_seasonal': False,                # default: False
    'atm_abs': 0.0,                      # default: 0.0
    'bog_a': 1.5,
    'bog_b': 2000
}

nml['mixed_layer_nml'] = {
    'albedo_value': 0.3,
    'depth': 40.0,   # default: 40.0
    'tconst': 305.,    # default: 305.0
    'delta_T': 40.,
    'prescribe_initial_dist': True,
    'evaporation': True,
    'do_qflux': False
}

nml['qe_moist_convection_nml'] = {
    # 'tau_bm': 7200.0,
    # 'rhbm': 0.7,  # default: 0.8
    # 'val_inc': 0.01,
    'Tmin': 160.,
    'Tmax': 350.
}

nml['lscale_cond_nml'] = {
    'do_simple':True,
    'do_evap': True,
    # 'hc': 1.0
}

nml['sat_vapor_pres_nml'] = {
    'do_simple': True
}

nml['damping_driver_nml'] = {
     'do_rayleigh': True,
     'trayfric': -0.5,              # neg. value: time in *days*
     'sponge_pbottom':  1800.,      # 18hPa, approx 4 sigma levels
     'do_conserve_energy': True,
     # 'do_cg_drag': False
}

nml['astro_nml'] = {
    'solr_cnst': 1360.            # default: 1368.22
}


diag = DiagTable()

#diag.add_file('6hourly', 6*60*60, 'seconds', time_units='days')
#diag.add_file('hourly', 1, 'hours', time_units='days')
diag.add_file('daily', 1, 'days', time_units='days')

diag.add_field('dynamics', 'ps')
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')
diag.add_field('dynamics', 'sphum')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')
diag.add_field('two_stream', 'lw_dtrans')

diag.add_field('mixed_layer', 't_surf')


exp.namelist = nml
exp.set_resolution('T85', 25)
exp.diag_table = diag

exp.clear_rundir()
exp.run(1, use_restart=False, num_cores=16)
for i in range(2, 19):
    exp.run(i, num_cores=16)
