from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = DryCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

cb.compile()

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp_name = 'held_suarez_default'
exp = Experiment(exp_name, codebase=cb)


diag = DiagTable()

# create one or more output files
diag.add_file('6hourly', 6*60*60, 'seconds', time_units='days')

# add diagnostic fields to the output files
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'omega   ')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')

exp.diag_table = diag

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 600,
        'days': 30,
        'calendar': 'thirty_day',
        'current_date': [2000,1,1,0,0,0]
    },

    'spectral_dynamics_nml': {
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
        'initial_sphum'           : 0.0,                  # default: 0
        'vert_coord_option'       : 'uneven_sigma',         # default: 'even_sigma'
        'scale_heights': 6.0,
        'exponent': 7.5,
        'surf_res': 0.5
    },

    'hs_forcing_nml': {
        't_zero': 315.,    # temperature at reference pressure at equator (default 315K)
        't_strat': 200.,   # stratosphere temperature (default 200K)
        'delh': 60.,       # equator-pole temp gradient (default 60K)
        'delv': 10.,       # lapse rate (default 10K)
        'eps': 0.,         # stratospheric latitudinal variation (default 0K)
        'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -40.,      # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -4.,      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1.,       # BL momentum frictional timescale (default 1 days)

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    'diag_manager_nml': {
        'mix_snapshot_average_fields': False
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    }
})

exp.namelist = namelist
exp.set_resolution('T42', num_levels=25)

exp.run(1, use_restart=False)
for i in range(2, 13):
  exp.run(i)  # use the restart i-1 by default
