""" 
This script configures a column model that uses Isca's columnwise physics routines. 

Useful for testing new convection / radiation parametrizations, as the dynamical core is 
bypassed so the model runs a gazillion times faster (especially if you're only simulating
one column). Can in principle simulate many (in lat and lon) at the same time. 

The wind is prescribed (it needs to be non-zero at the surface to allow for latent and 
sensible surface heat fluxes). Currently the user can set a namelist variable 'surface_wind'
that sets u_surf and v_surf = surface_wind / sqrt(2), so that wind_surf = sqrt(u_surf**2 + 
v_surf**2) = surface_wind. u and v at all other altitudes are set to zero (hardcoded). 

At the moment the model needs to use the vertical turbulent diffusion parameterization in order 
for the mixed layer code to work. This is not very consistent as the u and v wind are prescribed 
and so the u,v tendenency from the diffusion is thrown away. Hence an implicit assumption when 
using the column model is that 'the dynamics' would restore the surface winds to their prescribed 
speed, so that du/dt total is zero. 

The column model is currently initiated as a bit of a hack. The line 

'from isca import ColumnCodeBase'

sets a compiler flag -DCOLUMN_MODEL that tells the model to use the following files: 

atmos_column/column.F90
atmos_column/column_grid.F90
atmos_column/column_init_cond.F90
atmos_column/column_initialize_fields.F90

to initialize the model (including constructing the model grid), do the model timestepping 
(using a leapfrog scheme as before), and  handle input/output. 

Works with either hs_forcing, or the physics packages in idealized_moist_phys. Even when 
multiple columns are simulated, the model can only run on 1 core at the moment (will endeavour 
to fix this as some point). Also, the column model cannot read in topography input files. 

Any questions to Neil Lewis:  
neil.lewis@physics.ox.ac.uk 
"""


import os

import numpy as np

from isca import ColumnCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


# column model only uses 1 core
NCORES = 1

# compile code 
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = ColumnCodeBase.from_directory(GFDL_BASE)
cb.compile() 

# create an Experiment object to handle the configuration of model parameters
exp = Experiment('column_test_exp', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('column', 'ps', time_avg=True)
diag.add_field('column', 'bk')
diag.add_field('column', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('column', 'sphum', time_avg=True)
diag.add_field('column', 'ucomp', time_avg=True)
diag.add_field('column', 'vcomp', time_avg=True)
diag.add_field('column', 'temp', time_avg=True)
diag.add_field('two_stream', 'swdn_toa', time_avg=True)
diag.add_field('atmosphere', 'dt_ug_diffusion', time_avg=True)
diag.add_field('atmosphere', 'dt_vg_diffusion', time_avg=True)
exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':360,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    
    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'column_nml': {
        'lon_max': 1, # number of columns in longitude, default begins at lon=0.0
        'lat_max': 1, # number of columns in latitude, precise 
                      # latitude can be set in column_grid_nml if only 1 lat used. 
        'num_levels': 31,  # number of levels 
        'initial_sphum': 1e-6, 
    },

    'column_grid_nml': { 
        'lat_value': np.rad2deg(np.arcsin(1/np.sqrt(3))) # set latitude to that which causes insolation in frierson p2 radiation to be insolation / 4. 
        #'global_average': True # don't use this option at the moment
    },

    # set initial condition, NOTE: currently there is not an option to read in initial condition from a file. 
    'column_init_cond_nml': {
        'initial_temperature': 264., # initial atmospheric temperature 
        'surf_geopotential': 0.0, # applied to all columns 
        'surface_wind': 5. # as described above 
    },

    'idealized_moist_phys_nml': {
        'do_damping': False, # no damping in column model, surface wind prescribed 
        'turb':True,        # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'mixed_layer_bc':True, # need surface, how is this trying to modify the wind field? ****
        'do_simple': True, # simple RH calculation 
        'roughness_mom': 3.21e-05, # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
        'two_stream_gray': True,     #Use grey radiation
        'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use the simple Betts Miller convection scheme 
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use, which in this case is Frierson
        'do_seasonal': False,                #do_seasonal=false uses the p2 insolation profile from Frierson 2006. do_seasonal=True uses the GFDL astronomy module to calculate seasonally-varying insolation.
        'atm_abs': 0.2,                      # default: 0.0        
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7, # rh criterion for convection 
        'Tmin':160., # min temperature for convection scheme look up tables 
        'Tmax':350.  # max temperature for convection scheme look up tables 
    },
    
    'lscale_cond_nml': {
        'do_simple':True, # only rain 
        'do_evap':False,  # no re-evaporation of falling precipitation 
    },

    'surface_flux_nml': {
        'use_virtual_temp': True, # use virtual temperature for BL stability 
        'do_simple': True,
        'old_dtaudv': True    
    },

    'vert_turb_driver_nml': { # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':False,
        'evaporation':True,   
        'depth': 2.5,                          #Depth of mixed layer used
        'albedo_value': 0.30,                  #Albedo value used             
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True, 
    },

    # define pressure coordinate 
    'vert_coordinate_nml': {
        'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
        'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
       },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },
    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },
    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },
})

#Lets do a run!
if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2,4):
        exp.run(i, num_cores=NCORES)
