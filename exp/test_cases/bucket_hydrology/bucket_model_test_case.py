import numpy as np
import os
from gfdl.experiment import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']

baseexp = Experiment('bucket_test_experiment', overwrite_data=False)

#Add any input files that are necessary for a particular experiment.
baseexp.inputfiles = [os.path.join(base_dir,'input/land.nc'),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'bucket_depth', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

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

#Set physics scheme options
baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = False #Don't use grey radiation
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = True #Do use RRTM radiation
baseexp.namelist['idealized_moist_phys_nml']['convection_scheme'] = 'SIMPLE_BETTS_MILLER' #Use the simple betts-miller convection scheme
baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 150. #Setting the lower pressure boundary for the model sponge layer in Pa.
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2 #Parameter that sets the vertical distribution of sigma levels

baseexp.namelist['idealized_moist_phys_nml']['land_option'] = 'input' #Use land mask from input file
baseexp.namelist['idealized_moist_phys_nml']['land_file_name'] = 'INPUT/land.nc' #Tell model where to find input file

baseexp.namelist['mixed_layer_nml']['depth'] = 20. #Mixed layer depth
baseexp.namelist['mixed_layer_nml']['land_option'] = 'input' #Tell mixed layer to get land mask from input file
baseexp.namelist['mixed_layer_nml']['land_h_capacity_prefactor'] = 0.1 #What factor to multiply mixed-layer depth by over land. 
baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.25 #Ocean albedo value
baseexp.namelist['mixed_layer_nml']['land_albedo_prefactor'] = 1.3 #What factor to multiply ocean albedo by over land
baseexp.namelist['idealized_moist_phys_nml']['roughness_mom'] = 2.e-4 #Ocean roughness lengths
baseexp.namelist['idealized_moist_phys_nml']['roughness_heat'] = 2.e-4 #Ocean roughness lengths
baseexp.namelist['idealized_moist_phys_nml']['roughness_moist'] = 2.e-4 #Ocean roughness lengths

baseexp.namelist['idealized_moist_phys_nml']['bucket'] = True #Run with the bucket model
baseexp.namelist['idealized_moist_phys_nml']['init_bucket_depth_land'] = 1. #Set initial bucket depth over land
baseexp.namelist['idealized_moist_phys_nml']['max_bucket_depth_land'] = 2. #Set max bucket depth over land

baseexp.namelist['mixed_layer_nml']['do_qflux'] = False #Do not use prescribed qflux formula

baseexp.namelist['rrtm_radiation_nml']['solr_cnst'] = 1360. #s set solar constant to 1360, rather than default of 1368.22
baseexp.namelist['rrtm_radiation_nml']['dt_rad'] = 3600 #Set RRTM radiation timestep to 3600 seconds, meaning it runs every 5 atmospheric timesteps

#Lets do a run!
baseexp.runmonth(1, use_restart=False,num_cores=4, light=False)
for i in range(2,121):
    baseexp.runmonth(i, num_cores=4, light=False)
