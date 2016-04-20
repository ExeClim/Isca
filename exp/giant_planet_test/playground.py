import numpy as np
import os

from gfdl.experiment import Experiment, DiagTable

baseexp = Experiment('giant_planet_test', overwrite_data=True)

#s Define input files for experiment - by default they are found in exp_dir/input/

#baseexp.inputfiles = []

#s Define srcmods - by default they are found in exp_dir/srcmods/
baseexp.path_names.insert(0, os.path.join(os.getcwd(),'../../src/atmos_param/rayleigh_bottom_drag/rayleigh_bottom_drag.F90'))

diag = DiagTable()

diag.add_file('atmos_daily', 1, 'days', time_units='days')

# Define diag table entries 
diag.add_field('dynamics', 'ps', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'bk', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'pk', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'vor', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'div', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'ucomp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'vcomp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'temp', time_avg=True, files=['atmos_daily'])
diag.add_field('atmosphere', 'rh', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'slp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'omega', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'height', time_avg=True, files=['atmos_daily'])

diag.add_field('two_stream', 'tdt_rad', time_avg=True, files=['atmos_daily'])

diag.add_field('atmosphere', 'convection_rain', time_avg=True, files=['atmos_daily'])
diag.add_field('atmosphere', 'condensation_rain', time_avg=True, files=['atmos_daily'])

diag.add_field('atmosphere', 'diss_heat_ray', time_avg=True, files=['atmos_daily'])

diag.add_field('damping', 'diss_heat_rdamp', time_avg=True, files=['atmos_daily'])



baseexp.disable_rrtm()

baseexp.use_diag_table(diag)

baseexp.compile()

baseexp.clear_rundir()

#s Namelist changes from default values
baseexp.namelist['main_nml'] = {
     'days'   : 30,	
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,			
     'dt_atmos':1800,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'no_calendar'
}


baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False

baseexp.namelist['idealized_moist_phys_nml']['gp_surface'] = True
baseexp.namelist['idealized_moist_phys_nml']['mixed_layer_bc'] = False

baseexp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'Schneider'
baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False

baseexp.namelist['two_stream_gray_rad_nml']['solar_constant'] = 50.7
baseexp.namelist['two_stream_gray_rad_nml']['diabatic_acce'] = 10.0

baseexp.namelist['surface_flux_nml']['diabatic_acce'] = 10.0

baseexp.namelist['betts_miller_nml']['tau_bm'] = 21600.
baseexp.namelist['betts_miller_nml']['rhbm'] = 0.0

baseexp.namelist['constants_nml']['radius'] = 69860.0e3
baseexp.namelist['constants_nml']['grav'] = 26.0
baseexp.namelist['constants_nml']['omega'] = 1.7587e-4
baseexp.namelist['constants_nml']['orbital_period'] = 4332.589*86400.

baseexp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = 3.0e5
baseexp.namelist['spectral_dynamics_nml']['vert_coord_option'] = 'even_sigma'
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.1
baseexp.namelist['spectral_dynamics_nml']['exponent'] = 2.0
baseexp.namelist['spectral_dynamics_nml']['scale_heights'] = 5.0
baseexp.namelist['spectral_dynamics_nml']['valid_range_t'] =[50.,800.]

baseexp.namelist['spectral_dynamics_nml']['num_fourier'] = 85
baseexp.namelist['spectral_dynamics_nml']['num_spherical'] = 86
baseexp.namelist['spectral_dynamics_nml']['lon_max'] = 256
baseexp.namelist['spectral_dynamics_nml']['lat_max'] = 128
baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 30

baseexp.namelist['spectral_dynamics_nml']['damping_option'] = 'exponential_cutoff'
baseexp.namelist['spectral_dynamics_nml']['damping_order'] = 4
baseexp.namelist['spectral_dynamics_nml']['damping_coeff'] = 1.3889e-04
baseexp.namelist['spectral_dynamics_nml']['cutoff_wn'] = 30

baseexp.namelist['spectral_init_cond_nml']['initial_temperature'] = 200.

baseexp.namelist['rayleigh_bottom_drag_nml']['kf_days'] = 10.0
baseexp.namelist['rayleigh_bottom_drag_nml']['sigma_b'] = 0.9
baseexp.namelist['rayleigh_bottom_drag_nml']['rc'] = 0.84
baseexp.namelist['rayleigh_bottom_drag_nml']['H_lambda'] = 1000.0e3



for exp_number in [22]:
    exp = Experiment('giant_planet_test_%d' % exp_number, overwrite_data=False)
    exp.clear_rundir()

    exp.use_diag_table(diag)
    exp.execdir = baseexp.execdir

    exp.inputfiles = baseexp.inputfiles

    exp.namelist = baseexp.namelist.copy()

    exp.runmonth(1, use_restart=False,num_cores=4)
    for i in range(2, 49):
         exp.runmonth(i)
