from isca import DryCodeBase, DiagTable, Namelist, GFDL_BASE


def setup_held_suarez_diag():
    diag = DiagTable()
    diag.add_file('atmos_monthly', 30, 'days', time_units='days')
    # Tell model which diagnostics to write
    diag.add_field('dynamics', 'ps', time_avg=True)
    diag.add_field('dynamics', 'bk')
    diag.add_field('dynamics', 'pk')
    diag.add_field('dynamics', 'ucomp', time_avg=True)
    diag.add_field('dynamics', 'vcomp', time_avg=True)
    diag.add_field('dynamics', 'temp', time_avg=True)
    diag.add_field('dynamics', 'vor', time_avg=True)
    diag.add_field('dynamics', 'div', time_avg=True)
    return diag


def setup_held_suarez_codebase():
    cb = DryCodeBase.from_directory(GFDL_BASE)
    cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase
    return cb


def setup_held_suarez_namelist():
    namelist = Namelist({
        'main_nml': {
            'dt_atmos': 600,
            'days': 30,
            'calendar': 'thirty_day',
            'current_date': [2000, 1, 1, 0, 0, 0]
        },

        'atmosphere_nml': {
            'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
        },

        'spectral_dynamics_nml': {
            'damping_order': 4,  # default: 2
            'water_correction_limit': 200.e2,  # default: 0
            'reference_sea_level_press': 1.0e5,  # default: 101325
            'valid_range_t': [100., 800.],  # default: (100, 500)
            'initial_sphum': 0.0,  # default: 0
            'vert_coord_option': 'uneven_sigma',  # default: 'even_sigma'
            'scale_heights': 6.0,
            'exponent': 7.5,
            'surf_res': 0.5
        },

        # configure the relaxation profile
        'hs_forcing_nml': {
            't_zero': 315.,  # temperature at reference pressure at equator (default 315K)
            't_strat': 200.,  # stratosphere temperature (default 200K)
            'delh': 60.,  # equator-pole temp gradient (default 60K)
            'delv': 10.,  # lapse rate (default 10K)
            'eps': 0.,  # stratospheric latitudinal variation (default 0K)
            'sigma_b': 0.7,  # boundary layer friction height (default p/ps = sigma = 0.7)

            # negative sign is a flag indicating that the units are days
            'ka': -40.,  # Constant Newtonian cooling timescale (default 40 days)
            'ks': -4.,  # Boundary layer dependent cooling timescale (default 4 days)
            'kf': -1.,  # BL momentum frictional timescale (default 1 days)

            'do_conserve_energy': True,  # convert dissipated momentum into heat (default True)
        },

        'diag_manager_nml': {
            'mix_snapshot_average_fields': False
        },

        'fms_nml': {
            'domains_stack_size': 600000  # default: 0
        },

        'fms_io_nml': {
            'threading_write': 'single',  # default: multi
            'fileset_write': 'single',  # default: multi
        }
    })
    return namelist
