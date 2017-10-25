import numpy as np
import os
from gfdl import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']

baseexp = Experiment('axisymmetric_test_case', overwrite_data=False)

#Add any input files that are necessary for a particular experiment.
baseexp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
                      os.path.join(base_dir,'input/sn_1.000_sst.nc')]

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

baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.25  #Ocean albedo value

baseexp.namelist['rrtm_radiation_nml']['dt_rad'] = 3600 #Set RRTM radiation timestep to 3600 seconds, meaning it runs every 5 atmospheric timesteps


# RG The model can be made symmetric by setting the option below to True. However, with no eddies, the equinoctial state, with two (predominantly eddy driven) Hadley cells symmetric about the equator, is not stable. With a mixed layer ocean the model finds other solutions to avoid this state, e.g. keeping the ITCZ off the equator, and the SST becomes very flat in the tropics. Additionally, both with a mixed layer ocean and with fixed SSTs, the cells tend to develop pressure level scale vertical waves. I therefore recommend using prescribed SSTs, and including some vertical diffusion in the free atmosphere to help dissipate vertical waves. For this test case seasonally varying SSTs are provided.
baseexp.namelist['spectral_dynamics_nml']['make_symmetric'] = True # Make model zonally symmetric
baseexp.namelist['mixed_layer_nml']['do_sc_sst'] = True  # Use precribed SSTs
baseexp.namelist['mixed_layer_nml']['sst_file'] = 'sn_1.000_sst' # SSTs to use
baseexp.namelist['diffusivity_nml']['free_atm_diff'] = True  #Turn on vertical diffusion in the free atmosphere 

#Lets do a run!
baseexp.runmonth(1, use_restart=False, num_cores=4,light=False)
for i in range(2, 121):  
    baseexp.runmonth(i, num_cores=4, light=False)
    

    
    