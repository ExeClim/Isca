import numpy as np

import os

from gfdl import Experiment, DiagTable


baseexp = Experiment('bucket_test', overwrite_data=False)

baseexp.inputfiles = [os.path.join(os.getcwd(),'input/land.nc'), '/scratch/rg419/GFDL_model/GFDLmoistModel/input/ozone_1990.nc',
                      '/scratch/rg419/GFDL_model/GFDLmoistModel/input/ocean_qflux.nc' ]

diag_spinup = DiagTable()
diag_spinup.add_file('atmos_monthly', 30, 'days', time_units='days')
diag_spinup.add_field('dynamics', 'ps', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'bk', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'pk', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'ucomp', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'vcomp', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'temp', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'omega', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'height', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'sphum', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'slp', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('dynamics', 'zsurf', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('atmosphere', 'convection_rain', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('atmosphere', 'condensation_rain', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('atmosphere', 'bucket_depth', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('atmosphere', 'bucket_depth_cond', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('atmosphere', 'bucket_depth_conv', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('atmosphere', 'bucket_depth_lh', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('mixed_layer', 't_surf', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('mixed_layer', 'flux_lhe', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('mixed_layer', 'flux_t', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('rrtm_radiation', 'flux_sw', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('rrtm_radiation', 'flux_lw', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('rrtm_radiation', 'olr', time_avg=True, files=['atmos_monthly'])
diag_spinup.add_field('rrtm_radiation', 'toa_sw', time_avg=True, files=['atmos_monthly'])


baseexp.use_diag_table(diag_spinup)

baseexp.compile()

baseexp.clear_rundir()


baseexp.namelist['main_nml'] = {
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':720,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'thirty_day'
}

baseexp.namelist['mixed_layer_nml']['depth'] = 20.
baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.25 #0.3
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2
baseexp.namelist['spectral_dynamics_nml']['ocean_topog_smoothing'] = 0.8
baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 150.

baseexp.namelist['rrtm_radiation_nml']['dt_rad'] = 3600
baseexp.namelist['idealized_moist_phys_nml']['lwet_convection'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_bm'] = False

if __name__ == '__main__':
    

    #20m ERA land with topography and amip q fluxes
    exp = baseexp.derive('full_qflux_bucket_scale')
    exp.inputfiles = ['/scratch/rg419/GFDL_model/GFDLmoistModel/input/ocean_qflux.nc',
                      '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc', 
                      '/scratch/rg419/GFDL_model/GFDLmoistModel/input/ozone_1990.nc']
    exp.clear_rundir()
    exp.screen_runmonth_prefix = 'full_qflux'
    heat_cap_ratio = 2./20.
    exp.namelist['mixed_layer_nml']['land_h_capacity_prefactor'] = heat_cap_ratio
    exp.namelist['mixed_layer_nml']['land_albedo_prefactor'] = 1.3 #1.2
    exp.namelist['mixed_layer_nml']['load_qflux'] = True
    exp.namelist['mixed_layer_nml']['qflux_file_name'] = 'ocean_qflux'
    exp.namelist['mixed_layer_nml']['time_varying_qflux'] = True
    exp.namelist['idealized_moist_phys_nml']['land_option'] = 'input'
    exp.namelist['idealized_moist_phys_nml']['land_file_name'] = 'INPUT/land.nc'
    
    #Bucket params
    exp.namelist['idealized_moist_phys_nml']['bucket'] = True
    exp.namelist['idealized_moist_phys_nml']['init_bucket_depth_land'] = 1.
    exp.namelist['idealized_moist_phys_nml']['max_bucket_depth_land'] = 2.
    
    exp.namelist['mixed_layer_nml']['land_option'] = 'input'
    exp.namelist['spectral_init_cond_nml']['topography_option'] = 'input'
    exp.namelist['spectral_init_cond_nml']['topog_file_name'] = 'land.nc'
    exp.runmonth(1, use_restart=False, num_cores=16, light=True)
    for i in range(2, 121):  
        exp.runmonth(i, num_cores=16, light=True)
        
    
    