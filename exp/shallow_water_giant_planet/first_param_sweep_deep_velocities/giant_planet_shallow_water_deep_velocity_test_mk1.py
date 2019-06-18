import os

import numpy as np

from isca import ShallowCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

import pdb

NCORES = 16
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = ShallowCodeBase.from_directory(GFDL_BASE)

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

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
# diag.add_file('atmos_daily', 1, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('shallow_diagnostics', 'ucomp', time_avg=True)
diag.add_field('shallow_diagnostics', 'vcomp', time_avg=True)
diag.add_field('shallow_diagnostics', 'vor', time_avg=True)
diag.add_field('shallow_diagnostics', 'div', time_avg=True)
diag.add_field('shallow_diagnostics', 'h', time_avg=True)
diag.add_field('shallow_diagnostics', 'pv', time_avg=True)
diag.add_field('shallow_diagnostics', 'pv_corrected', time_avg=True)
diag.add_field('shallow_diagnostics', 'stream', time_avg=True)
diag.add_field('shallow_diagnostics', 'deep_geopot', time_avg=True)
# diag.add_field('shallow_diagnostics', 'trs', time_avg=True)
# diag.add_field('shallow_diagnostics', 'tr', time_avg=True)
diag.add_field('stirring_mod', 'stirring', time_avg=True)
# diag.add_field('stirring_mod', 'stirring_amp', time_avg=True)
diag.add_field('stirring_mod', 'stirring_sqr', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_kin', time_avg=True)
diag.add_field('shallow_diagnostics', 'h_sqd_mean', time_avg=True)

diag.add_field('shallow_diagnostics', 'e_kin_density', time_avg=True)
diag.add_field('shallow_diagnostics', 'eq_geopot', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_kin_real_units', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_pot_real_units', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_tot_real_units', time_avg=True)
diag.add_field('shallow_diagnostics', 'u_rms', time_avg=True)


#Empty the run directory ready to run

#Define values for the 'core' namelist
namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos': 1200,
     'calendar': 'no_calendar',
    },

 'atmosphere_nml':{
   'print_interval': 86400,
    },

'fms_io_nml':{
   'threading_write' :'single',
   'fileset_write': 'single'
    },

 'fms_nml':{
   'print_memory_usage':False,
   'domains_stack_size': 200000,
    },

 'shallow_dynamics_nml':{
   'num_lon'             : 512,
   'num_lat'             : 256,
   'num_fourier'         : 170,
   'num_spherical'       : 171,
   'fourier_inc'         : 1,
   'damping_option'      : 'resolution_dependent',
   'damping_order'       : 4,
   'damping_coeff'       : 1.e-04,
   'h_0'                 : 3.e04,
   'grid_tracer'         : False,
   'spec_tracer'         : False,
   'robert_coeff'        : 0.04,
   'robert_coeff_tracer' : 0.04,
    },

 'shallow_physics_nml':{
   'fric_damp_time'  :  0.0,
   'h_0'             :  3.e04, #Scott and Polvani begin with L_D_polar = 10 radii in size, i.e. depth of 9,910 metres.
   'h_amp'           :  0.,
   'h_itcz'          :  0.,
   },

 'stirring_nml': {
   'B':0.0,
   'do_localize': False,
   },

})

#Lets do a run!
if __name__=="__main__":

    # for ld_value in [10.0, 0.1, 0.05, 0.01]:
    for ld_value in [0.05,]:


        for u_deep_mag_val in [50.]:

            if u_deep_mag_val!=0.:
                u_deep_merid_arr = [27]
                # u_deep_merid_arr = [5,7,9]

            else:
                u_deep_merid_arr = [3]

            for u_deep_merid in u_deep_merid_arr:

                exp = Experiment('giant_planet_fixed_deep_ics_forced_no_rad_damping_ld_'+str(ld_value)+'_udeep_mag_'+str(u_deep_mag_val)+'_u_deep_merid_'+str(int(u_deep_merid))+'_energy_outputs_3', codebase=cb)

                exp.diag_table = diag 
                exp.namelist = namelist 
                exp.clear_rundir()

                rotation_period = ((9.*3600.)+55.*60 + 30.)
                omega = 2.*np.pi/ rotation_period
                radius = 69911e3
                grav = 24.79
                number_ld_in_radius_units = ld_value

                # Model uses geopotential as its height co-ordinate. So depth is h_0/grav.

                equilibrium_geopotential = (2.*number_ld_in_radius_units*omega*radius)**2.
                equilibrium_depth = equilibrium_geopotential/grav

                exp.update_namelist({
                    'shallow_dynamics_nml':{
                        'h_0': equilibrium_geopotential,
                        'u_deep_mag'   : u_deep_mag_val,
                        'n_merid_deep_flow': u_deep_merid,         
                        # 'u_upper_mag_init':  u_deep_mag_val,
                    },
                    'shallow_physics_nml': {
                        'h_0': equilibrium_geopotential, 
                        'therm_damp_time' : rotation_period * 0., #Thermal damping time in Scott and Polvani is 1 rotation period (v_l = 1)
                    },
                    'constants_nml': {
                        'omega': omega,
                        'radius': radius,
                    },
                    'stirring_nml': {
                        'decay_time':10.*rotation_period, #Scott and Polvani use decorrelation time of 10 planetary rotations - i.e. a long time for Jupiter. 
                        'amplitude':3.e-13,
                        'n_total_forcing_max': 85,
                        'n_total_forcing_min': 79,

                    },

                })

                exp.run(1, use_restart=False, num_cores=NCORES)
                for i in range(2,361):
                    exp.run(i, num_cores=NCORES)
