# define namelist values as python dictionary
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 110,
        'days': 0.,
        'seconds': 30. * 88440.,
        'calendar': 'no_calendar'
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb': True,
        'mixed_layer_bc': True,
        'do_virtual': False,
        'do_simple': True,
        'roughness_mom': 3.21e-05,
        'roughness_heat': 3.21e-05,
        'roughness_moist': 0.,
        'two_stream_gray': True,  # Use grey radiation
        'do_lscale_cond': False,
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,  # default: True
        'do_diffusivity': True,  # default: False
        'do_simple': True,  # default: False
        'constant_gust': 0.0,  # default: 1.0
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
        'use_actual_surface_temperatures': False,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'mixed_layer_nml': {
        'tconst': 285.,
        'prescribe_initial_dist': True,
        'evaporation': False,
        'albedo_value': 0.3,
    },

    'qe_moist_convection_nml': {
        'rhbm': 0.0,
        'tau_bm': 3600.,
    },

    'dry_convection_nml': {
        'tau': 7200,
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
        'do_simple': True,
        'tcmin': -223,
        'tcmax': 350.,
    },

    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.125,  # neg. value: time in *days*
        'sponge_pbottom': 0.5,  # Bottom of the model's sponge down to 50hPa (units are Pa)
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
        'bk': [0.00000000e+00, 1.53008955e-04, 4.63790800e-04, 1.10977640e-03, 2.37044150e-03, 4.74479200e-03,
               9.12245300e-03, 1.70677050e-02, 3.12516100e-02, 5.59939500e-02, 9.76165000e-02, 1.63754900e-01,
               2.60315150e-01, 3.85974250e-01, 5.28084800e-01, 6.65956600e-01, 7.81088000e-01, 8.65400050e-01,
               9.21109250e-01, 9.55343500e-01, 9.75416950e-01, 9.86856800e-01, 9.93269300e-01, 9.96830200e-01,
               9.98799150e-01, 1.00000000e+00],
        'pk': [0.] * 26,
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',  # Select radiation scheme to use, which in this case is Frierson
        'do_seasonal': True,
        'atm_abs': 0.2,
        'sw_diff': 0.0,
        'ir_tau_eq': 0.2,
        'ir_tau_pole': 0.2,
        'linear_tau': 1.0,
        'equinox_day': 0.0,
        'use_time_average_coszen': True,
        'solar_constant': 589.0,
    },
    'astronomy_nml': {
        'ecc': 0.0935,
        'obliq': 25.19,
        'use_mean_anom_in_rrsun_calc': True,
        'use_old_r_inv_squared': False
    },

    'constants_nml': {
        'orbital_period': 59166360,
        'solar_const': 589.0,
        'radius': 3396.0e3,
        'rdgas': 192.0,
        'kappa': 0.22727,
        'rotation_period': 88308,
    },

})
