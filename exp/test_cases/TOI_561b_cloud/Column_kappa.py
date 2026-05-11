""" 
This script configures a column model that uses Isca's columnwise physics routines. 

Single column configuration of Isca. Please cite McKim et al. (2024, submitted) (https://doi.org/10.22541/essoar.170904795.55675140/v1) if you use the SCM. 

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
n.t.lewis@exeter.ac.uk
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
    'main_nml': {
        'days'   : 0,#360, # each run lasts one year, and then multiple runs are strung together below (loop on e.g. Line 310)
        'hours'  : 0,   # a different output file is produced for each run (in this case, each year). Data 
        'minutes': 4000,   # output at the frequency specified in the diag table 
        'seconds': 0,
        'dt_atmos':5,#10,#360,#480, # 900s timestep for dynamical core
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
        'initial_sphum': 1e-3, 
        'q_decrease_only':True, # constrain q in stratosphere
    },

    'column_grid_nml': { 
        'lat_value': np.rad2deg(np.arcsin(1/np.sqrt(3))) # set latitude to that which causes insolation in frierson p2 radiation to be insolation / 4. 
        #'global_average': True # don't use this option at the moment
    },

    # set initial condition, NOTE: currently there is not an option to read in initial condition from a file (aside from a restart file). 
    'column_init_cond_nml': {
        'initial_temperature': 264., # initial atmospheric temperature 
        'surf_geopotential': 0.0, # applied to all columns 
        'surface_wind': 5. # as described above 
    },

    'idealized_moist_phys_nml': {
        'two_stream_gray': True, # DEFAULT: TRUE
        'do_rrtm_radiation': False, # DEFAULT: FALSE
                                   # Use RRTM radiation, not grey
        'convection_scheme': 'DRYADJ',#DRYADJ',#NONE',#SIMPLE_BETTS_MILLER', # DEAFULT: 'UNSET' 
                                                    # Use the simple Betts Miller convection scheme (Tan et al.)
        'do_damping': True, # DEFAULT: FALSE 
                            # turns on 'damping_driver', which manages the sponge at the top
        'turb':True, # DEFAULT: FALSE
                     # turns on boundary layer diffusion, managed by 'vert_turb_driver' and 'gcm_vert_diff'
        'mixed_layer_bc':True, # DEFAULT: FALSE
                               # turns on mixed layer bc 
        'do_virtual' :True, # DEFAULT: FALSE 
                            # determines whether virtual temperature is used for diffusion module (TRUE in Tan et al. nml)
        'roughness_mom':5.e-3, # DEFAULT: 0.05
        'roughness_heat':1,#1.e-5, # DEFAULT: 0.05
        'roughness_moist':1.e-5, # DEFAULT: 0.05
                                 # Each of these have been set to their value in Tan's nml
        'do_simple':False,
        'do_tj_bl':True,
    }

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use
        'do_seasonal': False,                #do_seasonal=false uses the p2 insolation profile
        'do_tl':True,
        'do_closein':True,
        'atm_abs': 0.,#1.,#22,                      # default: 0.0  
        'solar_exponent':1,
        'ir_tau_eq':1,#4.5, 
        'ir_tau_pole':1,#1.5, 
        'del_sol':1.2, 
        'solar_constant':4715*1361, 
        'R_stellar':0.843*6.957e8,
        'd_stellar':0.01055*1.496e11,
        'linear_tau':1.,#0.2, 
        'odp':1.0
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
        'use_virtual_temp': True, # DEFAULT: TRUE (TRUE in Tan+ nml) 
                                  # use virtual temp to compute stability of surface layer 
        'do_simple': False, # DEFAULT: FALSE (FALSE in Tan+ nml) 
                            # don't simplify computation of surface specific humidity at saturation 
        'old_dtaudv': False, # DEFAULT: FALSE (FALSE in Tan+ nml) 
                             # don't simplify derivative of surface wind stress (would set dwind_stress/du = dwind_stress/dv)  
        'gust_const':1.0, # DEFAULT: 1.0 (1.0 in Tan+ nml) 
                          # Set constant gustiness factor at surface for computation of surface fluxes 
        'tj_coeff':0.001*1e3,#0.0044*1e2, 
        'tj_coeff_m':0.001*1e3,#0.0044*1e2, # 008
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False, # DEFAULT: TRUE 
                                   # Turn off Mellor-Yamada scheme (FALSE in Tan et al. nml)
        'do_diffusivity': False, # DEFAULT: FALSE 
                                # Use non-local k-profile scheme for BL turbulence, as in (TRUE in Tan et al. nml)
        'do_edt':False, # DEFAULT: FALSE
                        # Turn off Bretherton-Grenier scheme (FALSE in Tan et al. nml)
        'constant_gust': 1.0,   # DEFAULT: 1.0
                                # NOTE: not used, constant gustiness is instead set i surface_flux_nml (1.0 in Tan et al. nml)
        'use_tau': False, # DEFAULT: FALSE
                          # use 'future' (t+dt) variables for computation of BL diffusion coefficients (FALSE in Tan et al. nml)
        'do_entrain':False, # DEFAULT: FALSE
                            # turns off BL entrainment parametrisation (FALSE in Tan et al. nml)
        'do_stable_bl':False, # DEFAULT: FALSE 
                              # turns of stable BL parametrisation (FALSE in Tan et al. nml)
        'do_shallow_conv':False, # DEFAULT: FALSE 
                                 # turns of shallow conv in BL (FALSE in Tan et al. nml) 
        'do_simple':False # DEFAULT: FALSE 
                          # turns off virtual temperature in mellor-yamada code (not used), so irrelevant here 
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': { ## DON'T have namelist for this from Tan+ because SLAB ocean is implemented differently there 
        'depth': 5, # DEFAULT: 40
                     # Use 30m mixed layer depth as described in Tan+ paper 
        'albedo_value': 0,#35, # DEFAULT: 0.06
                              # Surface albedo of 0.25 as in Tan+ paper for RRTM
        'prescribe_initial_dist':True, # DEFAULT: FALSE 
                                       # prescribes initial t_s with t_s = tconst - delta_T*(3*sin^2(lat)-1)/3
                                       # where tconst and delta_T are namelist options
        'tconst' : 2000., # DEFAULT: 305 
                         # tconst in expression above 
        'delta_T': 0., # DEFAULT: 40. 
                        # delta_T in expression above 
        'evaporation':False, # DEFAULT: True 
                            # allow surface evaporation and associated latent heat flux 
        'do_qflux': False,  # DEFAULT: False 
                            # don't do prescribed ocean heat transport 
        'load_qflux':False,  
        
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True, # DEFAULT: FALSE -- NEEDS TO BE CHANGED TO FALSE (FALSE in Tan+ nml)
                          # Mistake! This means model is using simplified expression: 
                          # esat = esat_0 * exp[ -(Lv/Rv) * (1/T - 1/T0) ] 
                          # instead of proper lookup table for esat 
                          # Might be inconsistent with using do_simple: FALSE in lscale_cond_nml above?
                          #
                          # Further note: 
                          # Tan+ nml has construct_table_wrt_liq = .true. and construct_table_wrt_liq_and_ice = .true. 
                          # (used when esat computed with lookup table) whereas the Isca defaults are .false. 
                          # however, I don't think the tables that result from this would be used in Isca (instead only the 
                          # default table, computed over water for all T), so probably safe to leave FALSE even if 
                          # do_simple = FALSE. 
        'tcmin_simple': -273,
        'tcmax_simple': 10000,#350.
                          
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
    exp.run(1, use_restart=False, num_cores=NCORES, mpirun_opts='--bind-to socket')
    for i in range(2,11):
        exp.run(i, num_cores=NCORES, mpirun_opts='--bind-to socket')
