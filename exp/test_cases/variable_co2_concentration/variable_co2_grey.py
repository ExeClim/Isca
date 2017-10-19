import numpy as np
import os
from gfdl.experiment import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']

baseexp = Experiment('variable_co2_grey', overwrite_data=False)

baseexp.inputfiles = [os.path.join(base_dir,'input/co2.nc')]

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
diag.add_field('two_stream', 'co2', time_avg=True)

baseexp.use_diag_table(diag)

#Turn off the full, slow radiation scheme compilation

baseexp.disable_rrtm()

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

baseexp.namelist['mixed_layer_nml']['depth'] = 2.5                          #Depth of mixed layer used
baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.38                  #Albedo value used

baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 25               #How many model pressure levels to use
baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True     #Use the grey radiation scheme
baseexp.namelist['idealized_moist_phys_nml']['convection_scheme'] ='SIMPLE_BETTS_MILLER' #Use simple Betts miller convection

baseexp.namelist['spectral_dynamics_nml']['vert_coord_option'] = 'input'   #Use the vertical levels from Frierson 2006
baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 5000.           #Bottom of the model's sponge down to 50hPa
baseexp.namelist['damping_driver_nml']['trayfric'] = -0.25                 #Drag timescale for model's sponge

baseexp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'byrne'        #Select radiation scheme to use
baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True          #do_seasonal=false uses the p2 insolation profile from Frierson 2006. do_seasonal=True uses the GFDL astronomy module to calculate seasonally-varying insolation.
baseexp.namelist['two_stream_gray_rad_nml']['equinox_day'] = 0.75          #A calendar parameter to get autumn equinox in september, as in the standard earth calendar.

baseexp.namelist['two_stream_gray_rad_nml']['do_read_co2'] = True #Read in CO2 timeseries from input file
baseexp.namelist['two_stream_gray_rad_nml']['co2_file'] = 'co2' #Tell model name of co2 input file

#Lets do a run!
baseexp.runmonth(1, use_restart=False,num_cores=4, light=False)
for i in range(2,121):
    baseexp.runmonth(i, num_cores=4, light=False)
