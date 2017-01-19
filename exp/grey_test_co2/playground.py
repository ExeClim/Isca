import numpy as np
import os

from gfdl.experiment import Experiment, DiagTable

baseexp = Experiment('amip_derived_qflux_control', overwrite_data=False)

#s Define input files for experiment - by default they are found in exp_dir/input/

baseexp.inputfiles = [os.path.join(os.getcwd(),'input/land.nc'),os.path.join(os.getcwd(),'input/ozone_1990.nc'),
                      os.path.join(os.getcwd(),'input/co2.nc')]

#s Define srcmods - by default they are found in exp_dir/srcmods/

diag = DiagTable()

diag.add_file('atmos_monthly', 1, 'months', time_units='days')
diag.add_file('atmos_6hourly', 6, 'hours', time_units='days')

# Define diag table entries 

#Add variables required by vertical interpolator to ALL files (note that if `files` is not specified in add_field then the variable is added to all files)
diag.add_field('dynamics', 'bk', time_avg=True)
diag.add_field('dynamics', 'pk', time_avg=True)
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)

diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('atmosphere','rh', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)


diag.add_field('dynamics', 'slp', time_avg=True, files=['atmos_monthly'])
diag.add_field('atmosphere', 'precipitation', time_avg=True, files=['atmos_monthly'])


diag.add_field('rrtm_radiation', 'tdt_rad', time_avg=True, files=['atmos_monthly'])
diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True, files=['atmos_monthly'])
diag.add_field('rrtm_radiation', 'flux_lw', time_avg=True, files=['atmos_monthly'])
diag.add_field('rrtm_radiation', 'tdt_sw', time_avg=True, files=['atmos_monthly'])
diag.add_field('rrtm_radiation', 'tdt_lw', time_avg=True, files=['atmos_monthly'])

diag.add_field('vert_turb', 'z_pbl', time_avg=True, files=['atmos_monthly'])

diag.add_field('mixed_layer', 'ml_heat_cap', time_avg=True, files=['atmos_monthly'])
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True, files=['atmos_monthly'])

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
     'dt_atmos':720,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'thirty_day'
}

baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False

baseexp.namelist['idealized_moist_phys_nml']['land_roughness_prefactor'] = 10.0
baseexp.namelist['idealized_moist_phys_nml']['roughness_mom'] = 2.e-4
baseexp.namelist['idealized_moist_phys_nml']['roughness_heat'] = 2.e-4
baseexp.namelist['idealized_moist_phys_nml']['roughness_moist'] = 2.e-4

baseexp.namelist['idealized_moist_phys_nml']['land_option'] = 'input'
baseexp.namelist['idealized_moist_phys_nml']['land_file_name'] = 'INPUT/land.nc'

baseexp.namelist['spectral_init_cond_nml']['topog_file_name'] = 'land.nc'
baseexp.namelist['spectral_init_cond_nml']['topography_option'] = 'input'

baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2

baseexp.namelist['mixed_layer_nml']['depth'] = 20.
baseexp.namelist['mixed_layer_nml']['delta_T'] = 0.

baseexp.namelist['mixed_layer_nml']['do_qflux'] = False

baseexp.namelist['mixed_layer_nml']['land_option'] = 'input'
baseexp.namelist['mixed_layer_nml']['land_h_capacity_prefactor'] = 0.1

baseexp.namelist['two_stream_gray_rad_nml']['do_read_co2'] = True
baseexp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'Byrne'


baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.25
baseexp.namelist['mixed_layer_nml']['land_albedo_prefactor'] = 1.3

baseexp.namelist['surface_flux_nml']['land_humidity_prefactor'] = 0.7

baseexp.namelist['qflux_nml']['qflux_amp'] = 0.0

baseexp.namelist['spectral_dynamics_nml']['ocean_topog_smoothing'] = 0.8

#s End namelist changes from default values

evap_res_name=1

start_month=[2]
run_length= [721]

do_qflux=[False]

for evap_res in do_qflux:
    evap_res_name = evap_res_name+1
    exp = Experiment('grey_test_co2_%d' % evap_res_name, overwrite_data=False)
    exp.clear_rundir()

    exp.use_diag_table(diag)
    exp.execdir = baseexp.execdir

    exp.inputfiles = baseexp.inputfiles
    exp.namelist = baseexp.namelist.copy()

    exp.runmonth(1, use_restart=False,num_cores=8)

