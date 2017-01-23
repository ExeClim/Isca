import numpy as np
import os
from gfdl import Experiment, DiagTable


baseexp = Experiment('rrtm_rot_test', overwrite_data=False)

baseexp.inputfiles = ['/scratch/rg419/GFDL_model/GFDLmoistModel/input/ozone_1990.nc']

diag = DiagTable()

diag.add_file('atmos_monthly', 30, 'days', time_units='days')

diag.add_field('dynamics', 'ps', time_avg=True, files=['atmos_monthly'])
diag.add_field('dynamics', 'bk', time_avg=True, files=['atmos_monthly'])
diag.add_field('dynamics', 'pk', time_avg=True, files=['atmos_monthly'])
diag.add_field('dynamics', 'ucomp', time_avg=True, files=['atmos_monthly'])
diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True, files=['atmos_monthly'])
        
diag.add_file('atmos_6hrly', 6, 'hours')
diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True, files=['atmos_6hrly'])
        
baseexp.use_diag_table(diag)

baseexp.compile()

baseexp.clear_rundir()


baseexp.namelist['main_nml'] = {
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':720,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'no_calendar'
}

baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.25
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2
baseexp.namelist['spectral_dynamics_nml']['ocean_topog_smoothing'] = 0.8

baseexp.namelist['rrtm_radiation_nml']['dt_rad'] = 3600
baseexp.namelist['idealized_moist_phys_nml']['lwet_convection'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_bm'] = False

baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 150.

baseexp.namelist['mixed_layer_nml']['depth'] = 10.

if __name__ == '__main__':
    
    #Vary rotation rates
    rots = [2., 0.5, 1.]
	
    for rot in rots:
        exp = baseexp.derive('rt_%.3f' % rot)
        exp.inputfiles = baseexp.inputfiles
        exp.clear_rundir()
        exp.screen_runmonth_prefix = 'r%.3f' % rot

        omega = rot*7.2921150e-5
        exp.namelist['constants_nml']['omega'] = omega

        exp.runmonth(1, use_restart=False, num_cores=16,light=True)
        for i in range(2, 25):  
            exp.runmonth(i, num_cores=16, light=True)

            
        
    
    
    
    

