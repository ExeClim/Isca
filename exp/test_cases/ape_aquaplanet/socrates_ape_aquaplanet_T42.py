import os
import numpy as np
from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 16  #32 is max for gv3
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = SocratesCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp = Experiment('ape_aqua_soc_low_res', codebase=cb)
exp.clear_rundir()

inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_file('atmos_daily',    1, 'days', time_units='days')

#when not an aquaplanet, you need to interpolate vertically!
diag.add_field('dynamics',    'ps', time_avg=True) #ps
diag.add_field('dynamics',    'bk')
diag.add_field('dynamics',    'pk')
diag.add_field('dynamics',    'zsurf', time_avg=True) 

diag.add_field('mixed_layer', 't_surf',   files=['atmos_monthly'], time_avg=True) #ts
diag.add_field('mixed_layer', 'flux_t',                            time_avg=True) #hfls
diag.add_field('mixed_layer', 'flux_lhe',                          time_avg=True) #hfss - LH is evap if / L_v
diag.add_field('mixed_layer', 'albedo',   files=['atmos_monthly'], time_avg=True) 

diag.add_field('atmosphere', 'precipitation',     time_avg=True) #pr
diag.add_field('atmosphere', 'condensation_rain', time_avg=True) #pr-pc
diag.add_field('atmosphere', 'rh',                files=['atmos_monthly'], time_avg=True) #hur

#surface wind stress
diag.add_field('atmosphere',   'flux_u', files=['atmos_monthly'], time_avg=True) #tauu - zonal component of stress
diag.add_field('atmosphere',   'flux_v', files=['atmos_monthly'], time_avg=True) #tauv

#near surface properties
diag.add_field('atmosphere',   'temp_2m',                        time_avg=True) 
diag.add_field('atmosphere',   'sphum_2m',                       time_avg=True) 
diag.add_field('atmosphere',   'u_10m', files=['atmos_monthly'], time_avg=True) 
diag.add_field('atmosphere',   'v_10m', files=['atmos_monthly'], time_avg=True) 

#common variables
diag.add_field('dynamics',    'sphum',  time_avg=True) #hus
diag.add_field('dynamics',    'ucomp',  time_avg=True) #ua
diag.add_field('dynamics',    'vcomp',  time_avg=True) #va
diag.add_field('dynamics',    'omega',  time_avg=True) #wap
diag.add_field('dynamics',    'temp',   time_avg=True) #ta
diag.add_field('dynamics',    'height', time_avg=True) #zg

#radiative fluxes
diag.add_field('socrates', 'soc_surf_flux_lw',      files=['atmos_monthly'], time_avg=True) #net
diag.add_field('socrates', 'soc_surf_flux_lw_down', files=['atmos_monthly'], time_avg=True) 

diag.add_field('socrates', 'soc_olr',               files=['atmos_monthly'], time_avg=True)

diag.add_field('socrates', 'soc_surf_flux_sw',      files=['atmos_monthly'], time_avg=True) #net
diag.add_field('socrates', 'soc_surf_flux_sw_down', files=['atmos_monthly'], time_avg=True)

diag.add_field('socrates', 'soc_toa_sw',            files=['atmos_monthly'], time_avg=True) #net
diag.add_field('socrates', 'soc_toa_sw_down',       files=['atmos_monthly'], time_avg=True)

#tendencies
diag.add_field('socrates', 'soc_tdt_lw',            files=['atmos_monthly'], time_avg=True)
diag.add_field('socrates', 'soc_tdt_sw',            files=['atmos_monthly'], time_avg=True)
diag.add_field('socrates', 'soc_tdt_rad',           files=['atmos_monthly'], time_avg=True) 

#needed for eddy flux terms
diag.add_field('dynamics', 'ucomp_vcomp', files=['atmos_monthly'], time_avg=True)
diag.add_field('dynamics', 'sphum_v',     files=['atmos_monthly'], time_avg=True)
diag.add_field('dynamics', 'vcomp_temp',  files=['atmos_monthly'], time_avg=True)

exp.diag_table = diag
exp.inputfiles = inputfiles

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':720,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },

    'socrates_rad_nml': {
        'stellar_constant':1365., #solar constant incoming
        'lw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'),
        'sw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'),
        'do_read_ozone': True,
        'ozone_file_name':'ozone_1990',
        'ozone_field_name':'ozone_1990',
        'dt_rad':4320.,
        'solday':90., #turn off seasonal cycle - diurnal by default
        'co2_ppmv':348.0,
        'store_intermediate_rad':True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
        'tidally_locked':False
    }, 

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,            
        'two_stream_gray': False,     #Use the grey radiation scheme
        'do_socrates_radiation': True,
        'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use simple Betts miller convection            
    },


    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,  
        'depth': 10.0,                          #Depth of mixed layer used
        'albedo_value': 0.38,                   #Albedo value used      
        'do_ape_sst': True,
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
           'do_simple':True,
           'construct_table_wrt_liq_and_ice':True
       },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,      
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

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':40,      #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    },

})

#Lets do a run!
if __name__=="__main__":

        cb.compile(debug=False)
        exp.set_resolution('T42')
        exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)#, run_idb=True)

        #for i in range(2,385): #all runs should be 30 years + spin up
        #    exp.run(i, num_cores=NCORES, overwrite_data=True)



