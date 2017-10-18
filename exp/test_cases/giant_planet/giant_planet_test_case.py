import numpy as np
import os
from gfdl.experiment import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']

baseexp = Experiment('giant_planet_test_experiment', overwrite_data=False)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps')
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation')
diag.add_field('mixed_layer', 't_surf')
diag.add_field('dynamics', 'sphum')
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')

baseexp.use_diag_table(diag)

#Turn off the full, slow radiation scheme compilation

baseexp.disable_rrtm()

#Compile model if not already compiled
baseexp.compile()

#Empty the run directory ready to run
baseexp.clear_rundir()

#s Namelist changes from default values
baseexp.namelist['main_nml'] = f90nml.Namelist({
     'days'   : 30,	
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,			
     'dt_atmos':1800,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'no_calendar'
})

baseexp.namelist['fms_nml']['domains_stack_size'] = 620000
baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False

baseexp.namelist['idealized_moist_phys_nml']['gp_surface'] = True
baseexp.namelist['idealized_moist_phys_nml']['mixed_layer_bc'] = False

baseexp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'Schneider'
baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False

baseexp.namelist['two_stream_gray_rad_nml']['solar_constant'] = 50.7
baseexp.namelist['two_stream_gray_rad_nml']['diabatic_acce'] = 1.0
baseexp.namelist['surface_flux_nml']['diabatic_acce'] = 1.0

baseexp.namelist['constants_nml']['radius'] = 69860.0e3
baseexp.namelist['constants_nml']['grav'] = 26.0
baseexp.namelist['constants_nml']['omega'] = 1.7587e-4
baseexp.namelist['constants_nml']['orbital_period'] = 4332.589*86400.
baseexp.namelist['constants_nml']['PSTD'] = 3.0E+06
baseexp.namelist['constants_nml']['PSTD_MKS'] = 3.0E+05
baseexp.namelist['constants_nml']['rdgas'] = 3605.38

baseexp.namelist['sat_vapor_pres_nml']['tcmin_simple'] = -223

baseexp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = 3.0e5
baseexp.namelist['spectral_dynamics_nml']['vert_coord_option'] = 'even_sigma'
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.1
baseexp.namelist['spectral_dynamics_nml']['exponent'] = 2.0
baseexp.namelist['spectral_dynamics_nml']['scale_heights'] = 5.0
baseexp.namelist['spectral_dynamics_nml']['valid_range_t'] =[50.,800.]

baseexp.namelist['spectral_dynamics_nml']['num_fourier'] = 213
baseexp.namelist['spectral_dynamics_nml']['num_spherical'] = 214
baseexp.namelist['spectral_dynamics_nml']['lon_max'] = 1024
baseexp.namelist['spectral_dynamics_nml']['lat_max'] = 320
baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 30

baseexp.namelist['spectral_dynamics_nml']['do_water_correction'] = False

baseexp.namelist['spectral_dynamics_nml']['damping_option'] = 'exponential_cutoff'
baseexp.namelist['spectral_dynamics_nml']['damping_order'] = 4
baseexp.namelist['spectral_dynamics_nml']['damping_coeff'] = 1.3889e-04
baseexp.namelist['spectral_dynamics_nml']['cutoff_wn'] = 100
baseexp.namelist['spectral_dynamics_nml']['initial_sphum'] = 0.0
baseexp.namelist['spectral_init_cond_nml']['initial_temperature'] = 200.
baseexp.namelist['idealized_moist_phys_nml']['convection_scheme'] = 'dry'

baseexp.namelist['rayleigh_bottom_drag_nml']= f90nml.Namelist({
	'kf_days':10.0,
	'do_drag_at_surface':True,
	'variable_drag': False
	})

baseexp.namelist['dry_convection_nml'] = f90nml.Namelist({
    'tau': 21600.,
    'gamma': 1.0, # K/km
})

#Lets do a run!
baseexp.runmonth(1, use_restart=False,num_cores=32, light=False)
for i in range(2,121):
    baseexp.runmonth(i, num_cores=32, light=False)
