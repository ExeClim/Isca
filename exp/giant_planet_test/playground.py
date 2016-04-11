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
     'dt_atmos':900,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'no_calendar'
}


baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False

baseexp.namelist['idealized_moist_phys_nml']['gp_surface'] = True
baseexp.namelist['idealized_moist_phys_nml']['mixed_layer_bc'] = False


baseexp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'SCHNEIDER'

baseexp.namelist['constants_nml']['radius'] = 69860.0e3
baseexp.namelist['constants_nml']['grav'] = 26.0
baseexp.namelist['constants_nml']['omega'] = 1.7587e-4
baseexp.namelist['constants_nml']['orbital_period'] = 4332.589*86400.

baseexp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = 3.0e5


for exp_number in [5]:
    exp = Experiment('giant_planet_test_%d' % exp_number)
    exp.clear_rundir()

    exp.use_diag_table(diag)
    exp.execdir = baseexp.execdir

    exp.inputfiles = baseexp.inputfiles

    exp.namelist = baseexp.namelist.copy()

    exp.runmonth(1, use_restart=False)
    for i in range(2, 13):
         exp.runmonth(i)
