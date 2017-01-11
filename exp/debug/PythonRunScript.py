import numpy as np
import pdb
import os
from gfdl.experiment import Experiment, DiagTable

num_cores = 1

run_idb=True

#Step 1: Do a single experiment to compile code then run remaining experiments off this compile.

compiled_exp = Experiment('debug',run_idb=run_idb)

compiled_exp.inputfiles = [os.path.join(os.getcwd(),'/scratch/pm366/GFDL_2013_FMS/GFDLmoistModel/input/ozone_1990.nc')]

#Step 3: Set up the diag table
diag = DiagTable()

#specify time periods
diag.add_file('atmos_daily', 1, 'days', time_units='days')
diag.add_file('atmos_monthly', -1, 'days', time_units='days')


#Both daily and monthly data  - don't specify files and all will be used
diag.add_field('dynamics', 'slp', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)
diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
diag.add_field('atmosphere', 'precipitation', time_avg=True)
#for pressure interpolations
diag.add_field('dynamics', 'bk', time_avg=True)
diag.add_field('dynamics', 'pk', time_avg=True)
diag.add_field('dynamics', 'ps', time_avg=True)


#Output data to test if model is stable after spin-up 
diag.add_field('dynamics', 'sphum', time_avg=True, files=['atmos_monthly'])


compiled_exp.use_diag_table(diag)
compiled_exp.compile()
compiled_exp.clear_rundir()

#namelists are defined here: /scratch/pm366/GFDL_2013_FMS/GFDLmoistModel/src/extra/python/gfdl/templates/


compiled_exp.namelist['main_nml'] = {
         'days'   : 30,	
         'hours'  : 0,
         'minutes': 0,
         'seconds': 0,			
         'dt_atmos':720,
         'current_date' : [0001,1,1,0,0,0],
         'calendar' : 'thirty_day'
    }

#set the vertical resolution. 0 is many lower trop levels and not many strat. Default 0.5.
compiled_exp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2
#set the radiation scheme
compiled_exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = False
compiled_exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = True
#set the convection scheme
compiled_exp.namelist['idealized_moist_phys_nml']['do_bm'] = False
compiled_exp.namelist['idealized_moist_phys_nml']['lwet_convection'] = True #run sbms

compiled_exp.namelist['rrtm_radiation_nml']['dt_rad'] = 4320 #default is 5*900=4500 - divisible by 24 hours could be a problem.
  


#Step 2 now use the compiled experiment to run new ones for each rh and tau


exp_name = ('mima_2013/debug/np{0}/')

exp = Experiment(exp_name, overwrite_data=False)
exp.clear_rundir()

#use existing diag table, namelist and paths
exp.use_diag_table(diag)
exp.execdir = compiled_exp.execdir   #place where compiled executable is
exp.inputfiles = compiled_exp.inputfiles
exp.namelist = compiled_exp.namelist.copy()
exp.runmonth(1, use_restart = False,  overwrite_data=False,num_cores = num_cores, run_idb = run_idb)

