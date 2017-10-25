import numpy as np
import os
from gfdl.experiment import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']

baseexp = Experiment('variable_co2_rrtm', overwrite_data=False)

baseexp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),os.path.join(base_dir,'input/co2.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('rrtm_radiation', 'co2', time_avg=True)

baseexp.use_diag_table(diag)

#Compile model if not already compiled
baseexp.compile()

#Empty the run directory ready to run
baseexp.clear_rundir()

#Define values for the 'core' namelist
baseexp.namelist['main_nml'] = f90nml.Namelist({
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':720,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'thirty_day'
})

#Use RRTM, not grey radiation:
baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = False
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = True

#Use the full Betts-miller convection scheme
baseexp.namelist['idealized_moist_phys_nml']['convection_scheme'] = 'FULL_BETTS_MILLER'


baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 150.  #Setting the lower pressure boundary for the model sponge layer in Pa.
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2 #Parameter that sets the vertical distribution of sigma levels

baseexp.namelist['mixed_layer_nml']['depth'] = 5. #Use shallow mixed-layer depth
baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.25 #set albedo value

baseexp.namelist['mixed_layer_nml']['do_qflux'] = False #Do not use prescribed form for q-fluxes

baseexp.namelist['rrtm_radiation_nml']['solr_cnst'] = 1360. #s set solar constant to 1360, rather than default of 1368.22
baseexp.namelist['rrtm_radiation_nml']['dt_rad'] = 4320 #Use 4320 as RRTM radiation timestep

baseexp.namelist['rrtm_radiation_nml']['do_read_co2'] = True #Read in CO2 timeseries from input file
baseexp.namelist['rrtm_radiation_nml']['co2_file'] = 'co2' #Tell model name of co2 input file

#Lets do a run!
baseexp.runmonth(1, use_restart=False,num_cores=4, light=False)
for i in range(2,121):
    baseexp.runmonth(i, num_cores=4, light=False)
