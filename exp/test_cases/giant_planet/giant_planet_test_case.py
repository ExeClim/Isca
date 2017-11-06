import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 4

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = IscaCodeBase.from_directory(GFDL_BASE)

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
exp = Experiment('giant_planet_test_experiment', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

exp.use_diag_table(diag)

#Turn off the full, slow radiation scheme compilation

#Empty the run directory ready to run
exp.clear_rundir()

#s Namelist changes from default values
exp.namelist = namelist = Namelist({
    'main_nml': {
     'days'   : 30,	
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,			
     'dt_atmos':1800,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'no_calendar'
    },

exp.namelist['fms_nml']['domains_stack_size'] = 620000 #Setting size of stack available to model, which needs to be higher than the default when running at high spatial resolution
exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True #Use grey radiation
exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False #Don't use RRTM
exp.namelist['idealized_moist_phys_nml']['convection_scheme'] = 'dry' #Use the dry convection scheme of Schneider & Walker

exp.namelist['idealized_moist_phys_nml']['gp_surface'] = True #Use the giant-planet option for the surface, meaning it's not a mixed layer, and applies energy conservation and a bottom-boundary heat flux
exp.namelist['idealized_moist_phys_nml']['mixed_layer_bc'] = False #Don't use the mixed-layer surface

exp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'Schneider' #Use the Schneider & Liu option for the grey scheme
exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False #Don't use seasonally-varying radiation 

exp.namelist['two_stream_gray_rad_nml']['solar_constant'] = 50.7 #Change solar constant
exp.namelist['two_stream_gray_rad_nml']['diabatic_acce'] = 1.0 #Parameter to artificially accelerate the diabatic processes during spinup. 1.0 performs no such acceleration, >1.0 performs acceleration
exp.namelist['surface_flux_nml']['diabatic_acce'] = 1.0 #Parameter to artificially accelerate the diabatic processes during spinup. 1.0 performs no such acceleration, >1.0 performs acceleration

#Set Jupiter constants
exp.namelist['constants_nml']['radius'] = 69860.0e3
exp.namelist['constants_nml']['grav'] = 26.0
exp.namelist['constants_nml']['omega'] = 1.7587e-4
exp.namelist['constants_nml']['orbital_period'] = 4332.589*86400.
exp.namelist['constants_nml']['PSTD'] = 3.0E+06
exp.namelist['constants_nml']['PSTD_MKS'] = 3.0E+05
exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = 3.0e5
exp.namelist['constants_nml']['rdgas'] = 3605.38

exp.namelist['sat_vapor_pres_nml']['tcmin_simple'] = -223 #Make sure low temperature limit of saturation vapour pressure is low enough that it doesn't cause an error (note that this giant planet has no moisture anyway, so doesn't directly affect calculation.

exp.namelist['spectral_dynamics_nml']['vert_coord_option'] = 'even_sigma' #Use equally-spaced sigma levels
exp.namelist['spectral_dynamics_nml']['surf_res'] = 0.1 #Parameter for setting vertical distribution of sigma levels
exp.namelist['spectral_dynamics_nml']['exponent'] = 2.0 #Parameter for setting vertical distribution of sigma levels
exp.namelist['spectral_dynamics_nml']['scale_heights'] = 5.0 #Number of vertical scale-heights to include
exp.namelist['spectral_dynamics_nml']['valid_range_t'] =[50.,800.] #Valid range of temperatures in Kelvin

exp.namelist['spectral_dynamics_nml']['num_fourier'] = 213 #Number of Fourier modes
exp.namelist['spectral_dynamics_nml']['num_spherical'] = 214 #Number of spherical harmonics in triangular truncation
exp.namelist['spectral_dynamics_nml']['lon_max'] = 1024 #Lon grid points
exp.namelist['spectral_dynamics_nml']['lat_max'] = 320 #Lat grid points
exp.namelist['spectral_dynamics_nml']['num_levels'] = 30 #Number of vertical levels

exp.namelist['spectral_dynamics_nml']['do_water_correction'] = False #Turn off enforced water conservation as model is dry

exp.namelist['spectral_dynamics_nml']['damping_option'] = 'exponential_cutoff' #Use the high-wavenumber filter option
exp.namelist['spectral_dynamics_nml']['damping_order'] = 4
exp.namelist['spectral_dynamics_nml']['damping_coeff'] = 1.3889e-04
exp.namelist['spectral_dynamics_nml']['cutoff_wn'] = 100

#Set initial conditions
exp.namelist['spectral_dynamics_nml']['initial_sphum'] = 0.0 #No initial specific humidity
exp.namelist['spectral_init_cond_nml']['initial_temperature'] = 200. #Lower than normal initial temperature

#Set parameters for near-surface Rayleigh drag
exp.namelist['rayleigh_bottom_drag_nml']= f90nml.Namelist({
	'kf_days':10.0,
	'do_drag_at_surface':True,
	'variable_drag': False
	})

#Set parameters for dry convection scheme
exp.namelist['dry_convection_nml'] = f90nml.Namelist({
    'tau': 21600.,
    'gamma': 1.0, # K/km
})

#Lets do a run!
exp.runmonth(1, use_restart=False,num_cores=32, light=False)
for i in range(2,121):
    exp.runmonth(i, num_cores=32, light=False)
