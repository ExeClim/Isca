import f90nml

# Adapted from the original core.nml and phys.nml files included in 2006 codebase.

# where the value is the same as the default in the code base, it is commented out
# and therefore not included in the namelist.
# Where the value is different, the code default is shown in an inline comment
basic =  f90nml.Namelist({})

basic['spectral_dynamics_nml'] = {
    #'damping_option'          : 'resolution_dependent',
    'damping_order'           : 4,                      # default: 2
    #'do_mass_correction':      True
    #'do_energy_correction':    True
    #'do_water_correction':     True
    'water_correction_limit'  : 200.e2,                 # default: 0
    #'use_virtual_temperature': False
    #'vert_advect_uv'          : 'second_centered',
    #'vert_advect_t'           : 'second_centered',
    #'longitude_origin'        : 0.0,
    #'robert_coeff'            : .03,                   # default: 0.04
    #'alpha_implicit'          : .5,
    'reference_sea_level_press':1.0e5,                  # default: 101325
    #'lon_max'                 : 128,
    #'lat_max'                 : 64,
    'num_levels'              : 25,                     # default: 18
    #'num_fourier'             : 42,
    #'num_spherical'           : 43,
    #'fourier_inc'             : 1,
    #'triang_trunc'            :True
    'valid_range_t'           : [100., 800.],           # default: (100, 500)
    #'initial_state_option'   : 'quiescent'
    'initial_sphum'           : 2.e-6,                  # default: 0
    'vert_coord_option'       : 'uneven_sigma',         # default: 'even_sigma'
    'surf_res'                : 0.2,                    # default: 0.1
    'scale_heights'           : 11.0,                   # default: 4
    'exponent'                : 7.0,                    # default: 2.5
}

basic['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*30,
    'calendar': 'no_calendar'
}

basic['diag_manager_nml'] = {
    'mix_snapshot_average_fields': False
}

basic['fms_nml'] = {
    'domains_stack_size': 600000                        # default: 0
}

basic['fms_io_nml'] = {
    'threading_write': 'single',                         # default: multi
    'fileset_write': 'single',                           # default: multi
}

# from phys.nml
basic['idealized_moist_phys_nml'] = {
    #'two_stream_gray': True,
    #'do_rrtm_radiation': False,
    'convection_scheme': 'betts_miller',
    'do_damping': True,
    'turb': True,
    #'lwet_convection': False,
    #'do_bm': True,
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


basic['vert_turb_driver_nml'] = {
   'do_mellor_yamada': False,     # default: True
   'do_diffusivity': True,        # default: False
   'do_simple': True,             # default: False
   #'do_shallow_conv': False,
   #'gust_scheme': 'constant',
   'constant_gust': 0.0,          # default: 1.0
   #'use_tau': False
}

basic['diffusivity_nml'] = {
    'do_entrain': False,          # default: True
    'do_simple': True,            # default: False
    #'frac_inner': 0.1,
    #'rich_crit_pbl':  1.0
}

# basic['monin_obukhov_nml'] = {
#     'neutral': False,
#     'rich_crit': 2.0,
#     'stable_option': 1
# }

basic['surface_flux_nml'] = {
    'use_virtual_temp': False,
    'do_simple': True,
    'old_dtaudv': True
}

basic['atmosphere_nml'] = {
    'idealized_moist_model': True
}

# basic['spectral_init_cond_nml'] = {
#     'initial_temperature': 264.0
# }

basic['two_stream_gray_rad_nml'] = {
    'rad_scheme': 'frierson',           # default: frierson
    'do_seasonal': True,                # default: False
    #'linear_tau': 0.1,
    #'solar_constant': 1360.,
    #'solar_exponent': 4.0,
    #'ir_tau_pole': 1.5,
    #'ir_tau_eq': 6.0,
    #'del_sol': 1.4,
    'atm_abs': 0.2                      # default: 0.0
}

basic['mixed_layer_nml'] = {
    'albedo_value': 0.27, # default: 0.06
    'depth': 100.0,   # default: 40.0
    'tconst': 285.,    # default: 305.0
    #'delta_T': 40.,
    'prescribe_initial_dist': True,
    'evaporation': True,
    # 'do_qflux': False
}

basic['qe_moist_convection_nml'] = {
    # 'tau_bm': 7200.0,
    # 'rhbm': 0.7,  # default: 0.8
    # 'val_inc': 0.01,
    'Tmin': 160.,
    'Tmax': 350.
}

basic['lscale_cond_nml'] = {
    'do_simple':True,
    'do_evap': True,
    # 'hc': 1.0
}

basic['sat_vapor_pres_nml'] = {
    'do_simple': True
}

basic['damping_driver_nml'] = {
     'do_rayleigh': True,
     'trayfric': -0.5,              # neg. value: time in *days*
     'sponge_pbottom':  50.,
     'do_conserve_energy': True,
     # 'do_cg_drag': False
}

# basic['rrtm_radiation_nml'] = {
#      'h2o_lower_limit': 2.e-07,
#      'co2ppmv': 300.,
#      'do_read_ozone': True,
#      'ozone_file': 'ozone_1990',
#      'dt_rad':  4500
# }

basic['astro_nml'] = {
    'solr_cnst': 1360.            # default: 1368.22
}

basic['betts_miller_nml'] = {
   # 'tau_bm': 7200.,
   'rhbm': .7   ,   # default: .8
   'do_simp': False,
   'do_shallower': True,
   # 'do_changeqref': False,
   # 'do_envsat': False,
   # 'do_taucape': False,
   # 'capetaubm': 900.,
   # 'tau_min': 2400.
}

# basic['qflux_nml'] = {
#     'qflux_amp': 30.
# }

moist = basic.copy()


dry = basic.copy()

del dry['betts_miller_nml']
dry['idealized_moist_phys_nml']['convection_scheme'] = 'dry'

dry['dry_convection_nml'] = {
    'tau': 14400.0,  # from tapios/fms-idealized
    'gamma': 1.0,
}

dry['lscale_cond_nml']['do_evap'] = False

dry['spectral_dynamics_nml']['initial_sphum'] = 0.0

dry['mixed_layer_nml'].update({
    'evaporation': False,
    'do_qflux': False
})