import numpy as np
import os

from gfdl.experiment import Experiment, DiagTable

baseexp = Experiment('co2_test_fresh', overwrite_data=True)

#s Define input files for experiment - by default they are found in exp_dir/input/

baseexp.inputfiles = [os.path.join(os.getcwd(),'input/land.nc'),os.path.join(os.getcwd(),'input/ozone_1990.nc'),os.path.join(os.getcwd(),'input/co2.nc')]

#s Define srcmods - by default they are found in exp_dir/srcmods/
baseexp.path_names.insert(0, os.path.join(os.getcwd(),'srcmods/surface_flux.F90'))

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
diag.add_field('dynamics', 'height_half', time_avg=True, files=['atmos_daily'])

diag.add_field('mixed_layer', 't_surf', time_avg=True, files=['atmos_daily'])

diag.add_field('atmosphere', 'convection_rain', time_avg=True, files=['atmos_daily'])
diag.add_field('atmosphere', 'condensation_rain', time_avg=True, files=['atmos_daily'])

diag.add_field('rrtm_radiation', 'tdt_rad', time_avg=True, files=['atmos_daily'])
diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True, files=['atmos_daily'])
diag.add_field('rrtm_radiation', 'flux_lw', time_avg=True, files=['atmos_daily'])
diag.add_field('rrtm_radiation', 'tdt_sw', time_avg=True, files=['atmos_daily'])
diag.add_field('rrtm_radiation', 'tdt_lw', time_avg=True, files=['atmos_daily'])
diag.add_field('rrtm_radiation', 'co2', time_avg=True, files=['atmos_daily'])
diag.add_field('rrtm_radiation', 'ozone', time_avg=True, files=['atmos_daily'])

diag.add_field('damping', 'udt_rdamp', time_avg=True, files=['atmos_daily'])
diag.add_field('damping', 'vdt_rdamp', time_avg=True, files=['atmos_daily'])
diag.add_field('damping', 'tdt_diss_rdamp', time_avg=True, files=['atmos_daily'])

diag.add_field('vert_turb', 'z_pbl', time_avg=True, files=['atmos_daily'])

diag.add_field('mixed_layer', 'flux_lhe', time_avg=True, files=['atmos_daily'])
diag.add_field('mixed_layer', 'flux_t', time_avg=True, files=['atmos_daily'])
diag.add_field('mixed_layer', 'ml_heat_cap', time_avg=True, files=['atmos_daily'])


baseexp.use_diag_table(diag)

baseexp.compile()

baseexp.clear_rundir()

#s Namelist changes from default values
baseexp.namelist['main_nml'] = {
     'days'   : 1,	
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,			
     'dt_atmos':900,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'thirty_day'
}


# Choose to use JP's new value of omega, as now it doesn't really make any difference. 
#baseexp.namelist['constants_nml'] = {
#    'omega': 7.292e-5
#}



baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = False
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_bm'] = False
baseexp.namelist['idealized_moist_phys_nml']['lwet_convection'] = False

baseexp.namelist['rrtm_radiation_nml']['do_read_co2'] = True
baseexp.namelist['rrtm_radiation_nml']['co2_file'] = 'co2'
baseexp.namelist['rrtm_radiation_nml']['do_read_ozone'] = True


baseexp.namelist['idealized_moist_phys_nml']['land_roughness_prefactor'] = 10.0
baseexp.namelist['idealized_moist_phys_nml']['land_option'] = 'input'
baseexp.namelist['idealized_moist_phys_nml']['land_file_name'] = 'INPUT/land.nc'

baseexp.namelist['mixed_layer_nml']['depth'] = 20.
baseexp.namelist['mixed_layer_nml']['delta_T'] = 0.
baseexp.namelist['mixed_layer_nml']['do_qflux'] = False
baseexp.namelist['mixed_layer_nml']['land_option'] = 'input'
baseexp.namelist['mixed_layer_nml']['land_h_capacity_prefactor'] = 0.1
baseexp.namelist['mixed_layer_nml']['land_albedo_prefactor'] = 1.0

baseexp.namelist['qflux_nml']['qflux_amp'] = 0.0

#baseexp.namelist['fms_nml']['domains_stack_size'] = 800000

#s Using perpetual equinox
#baseexp.namelist['astro_nml']['solday'] = 90.0 
#s End namelist changes from default values


for evap_res in [0.7]:
    evap_res_name = 1
    exp = Experiment('co2_test_fresh_%d' % evap_res_name)
    exp.clear_rundir()

    exp.use_diag_table(diag)
    exp.execdir = baseexp.execdir

    exp.inputfiles = baseexp.inputfiles

    exp.namelist = baseexp.namelist.copy()

    exp.namelist['surface_flux_nml']['land_humidity_prefactor'] = evap_res

    exp.runmonth(1, use_restart=False)
    for i in range(2, 20):
        exp.runmonth(i)
