import os
from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

import sys
sys.path.insert(0, os.path.join(GFDL_BASE, 'src/extra/python/scripts'))
from create_pk_heating_input_file import create_polar_and_midlat_heating_file

NCORES = 16
RESOLUTION = 'T42', 60  # T42 horizontal resolution, 60 levels in pressure

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

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp_name = 'polvani_kushner_polar_and_midlat_heating'
exp = Experiment(exp_name, codebase=cb)

# Generate the combined polar + midlatitude wavenumber-2 heating input file
# (all zeros unless the two functions below are given non-default arguments)
# rather than committing a large binary to git - see
# src/extra/python/scripts/create_pk_heating_input_file.py.
input_dir = os.path.join(GFDL_BASE, 'exp/test_cases/polvani_kushner/input')
os.makedirs(input_dir, exist_ok=True)
heating_file_path, heating_var_name = create_polar_and_midlat_heating_file(input_dir)

exp.inputfiles = [heating_file_path]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_daily', 1, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)  # required for plevel interpolation
diag.add_field('dynamics', 'zsurf')                 # static, so can't be time-averaged
diag.add_field('dynamics', 'bk')                    # required for plevel interpolation
diag.add_field('dynamics', 'pk')                    # required for plevel interpolation
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)  # vertical velocity, e.g. for EP flux
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('hs_forcing', 'teq', time_avg=True)
diag.add_field('hs_forcing', 'local_heating', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)

exp.diag_table = diag

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 240,
        'days': 30,
        'calendar': 'thirty_day',
        'current_date': [2000,1,1,0,0,0]
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    'spectral_dynamics_nml': {
        'damping_order'           : 4,                      # default: 2
        'do_water_correction': False,
        'reference_sea_level_press': 1.0e5,                  # default: 101325
        'valid_range_t'           : [50., 800.],           # default: (100, 500)
        'initial_sphum'           : 0.0,                  # default: 0
        'vert_coord_option'       : 'uneven_sigma',         # default: 'even_sigma'
        'scale_heights': 11.0,
        'exponent': 3.0,
        'surf_res': 0.5
    },

    # configure the relaxation profile
    'hs_forcing_nml': {
        't_zero': 315.,      # temperature at reference pressure at equator (default 315K)
        't_strat': 216.65,   # stratosphere temperature - consistent with US standard T at 20km
        'delh': 60.,         # equator-pole temp gradient (default 60K)
        'delv': 10.,         # lapse rate (default 10K)
        'sigma_b': 0.7,      # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -40.,        # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -4.,        # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1.,         # BL momentum frictional timescale (default 1 days)

        # jet-latitude control, following Garfinkel et al. (2013) - off by default here
        'A': 0., # takes values 0, +-5, +-10
        'B': 0., # takes values 0 to 20 in multiples of 4
        'P_opt': 'Option1', # 'Option1' or 'Option2' depending on jet location requirement

        # stratospheric polar vortex, following Polvani & Kushner (2002)
        'equilibrium_t_option': 'Polvani_Kushner',
        'strat_vtx': True,  # set to False for w_vtx=0, i.e. no polar vortex
        'eps': 0.,          # stratospheric latitudinal variation (+-10 in the P-K paper)
        'vtx_gamma': 4.0,   # lapse rate of winter stratospheric cooling (default 4 K/km)
        'z_ozone': 13.,     # height of stratospheric heating source (km)
        'do_conserve_energy': True,  # convert dissipated momentum into heat (default True)
        'sponge_flag': True,         # sponge layer damping in upper levels

        # prescribed polar cap + midlatitude wavenumber-2 heating perturbation
        'local_heating_option': 'from_file',
        'local_heating_file': heating_var_name,
    },

    'diag_manager_nml': {
        'mix_snapshot_average_fields': False
    },

    'fms_nml': {
        'domains_stack_size': 600000  # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',  # default: multi
        'fileset_write': 'single',    # default: multi
    }
})

exp.namelist = namelist
exp.set_resolution(*RESOLUTION)

#Lets do a run!
if __name__ == '__main__':
    exp.run(1, num_cores=NCORES, use_restart=False)
    for i in range(2, 13):
        exp.run(i, num_cores=NCORES)  # use the restart i-1 by default
