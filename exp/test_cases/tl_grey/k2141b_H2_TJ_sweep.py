import os

import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES =8

# Point to code as defined by $GFDL_BASE
cb = GreyCodeBase.from_directory(GFDL_BASE)

base_dir = os.path.dirname(os.path.realpath(__file__))

#cb.use_single_precision()
cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('tl_grey', codebase=cb)



#Tell model how to write diagnostics
# Begin with yearly, this is modified to 6hrly for the last 10 years of simulation below (Line 314 and Line 338)
diag = DiagTable()
diag.add_file('atmos_10period', 400, 'minutes', time_units='minutes')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('two_stream', 'flux_sw', time_avg=True)
diag.add_field('two_stream', 'flux_lw', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('vert_turb', 'z_pbl', time_avg=True)
diag.add_field('vert_turb', 'diff_t', time_avg=True)

# diag2 = DiagTable()
# diag2.add_file('atmos_monthly', 30, 'days', time_units='days')

# #Tell model which diagnostics to write
# diag2.add_field('dynamics', 'ps', time_avg=True)
# diag2.add_field('dynamics', 'bk')
# diag2.add_field('dynamics', 'pk')
# diag2.add_field('dynamics', 'ucomp', time_avg=True)
# diag2.add_field('dynamics', 'vcomp', time_avg=True)
# diag2.add_field('dynamics', 'temp', time_avg=True)
# diag2.add_field('dynamics', 'sphum', time_avg=True)
# diag2.add_field('dynamics', 'height', time_avg=True)
# diag2.add_field('dynamics', 'ucomp_height', time_avg=True)
# diag2.add_field('dynamics', 'vcomp_height', time_avg=True)
# diag2.add_field('dynamics', 'ucomp_temp', time_avg=True)
# diag2.add_field('dynamics', 'vcomp_temp', time_avg=True)
# diag2.add_field('dynamics', 'sphum_u', time_avg=True)
# diag2.add_field('dynamics', 'sphum_v', time_avg=True)
# diag2.add_field('atmosphere', 'precipitation', time_avg=True)
# diag2.add_field('atmosphere', 'convection_rain', time_avg=True)
# diag2.add_field('two_stream', 'flux_sw', time_avg=True)
# diag2.add_field('two_stream', 'flux_lw', time_avg=True)
# diag2.add_field('mixed_layer', 'flux_lhe', time_avg=True)
# diag2.add_field('mixed_layer', 'flux_t', time_avg=True)
# diag2.add_field('mixed_layer', 't_surf', time_avg=True)
# diag2.add_field('vert_turb', 'z_pbl', time_avg=True)



exp.diag_table = diag # register diag table 


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
    
    'diffusivity_nml': {
        'do_entrain':False, # DEFAULT: TRUE (TRUE in Tan+ nml)
                            # Turn off modification of diffusion coefficient to represent overshooting convection at BL top 
                            # Mistake! - although doesn't matter as 'entr_ratio' in Tan+ nml is 0.0 
                            # which turns off this code anyway 
        'entr_ratio': 0.0, # DEFAULT: 0.0, turns off BL entrainment code, (Tan+ nml is 0.0)
        'free_atm_diff':False, # DEFAULT: FALSE (FALSE in Tan+ nml) 
                               # Turns off diffusion in free atmosphere 
        'do_simple': True, # DEFAULT: FALSE (FALSE in Tan+ nml)
                            # Use virtual temperature for diffusion / when computing stability 
        'parcel_buoy': 0.0, # DEFAULT: 2.0 -- NEEDS TO BE CHANGED to 0.0 (0.0 in Tan+ nml)
                            # Mistake! This is 0.0 in Tan+ nml. 
                            # Making this change will lead to slightly reduced BL depth in unstable conditions. 
        'frac_inner': 0.1,  # DEFAULT 0.1 (0.1 in Tan+ nml)
                            # Determines depth of surface layer in BL 
        'fixed_depth': False, # DEFAULT: FALSE (FALSE in Tan+ nml)
                              # calculate pbl depth using stability criterion 
        'tj_coeff':0.001*1e3,#0.0044*1e2, 
        'tj_coeff_m':0.001*1e3,#0.0044*1e2, # 008
        'tj_bl_pres':85000.,
        'tj_strato_pres':10000./2.
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

    'atmosphere_nml': {
        'idealized_moist_model': True # Use idealized_moist_phys 
    },


    'mixed_layer_nml': { ## DON'T have namelist for this from Tan+ because SLAB ocean is implemented differently there 
        'depth': 5, # DEFAULT: 40
                     # Use 30m mixed layer depth as described in Tan+ paper 
        'albedo_value': 0.,#35, # DEFAULT: 0.06
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
    

#     'qe_moist_convection_nml': { # Namelist for simple betts--miller convection 
#         'rhbm':0.7, # DEFAULT: 0.8 (0.7 in Tan+ nml)
#                     # relative humidity of reference profile used by convection scheme
#         'tau_bm':7200., # DEFAULT: 7200. (7200. in Tan+ nml)
#                         # timescale for convective relaxation 
#         'Tmin':120., # DEFAULT: 173. (120. in Tan+ nml) 
#                      # minimum temperature at LCL 
#         'Tmax':360.,  # DEFAULT: 335 (360. in Tan+ nml) 
#                      # maximum temperature at LCL 
#         'val_inc': 0.01, #DEFAULT: 0.01 (0.01 in Tan+ nml) 
#                      # increment in value for LCL lookup table
#         'precision':1.e-6,#e-6
#     },
    
    'lscale_cond_nml': {
        'do_simple':True, # DEFAULT: FALSE (FALSE in Tan+ nml)
                           # if true then latent heat computation ignores latent heat of sublimation 
        'do_evap':False # DEFAULT: FALSE, (FALSE in Tan+ nml)
                        # turns off re-evaporation of falling precipitation 
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
    
    
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  260.,
        'do_conserve_energy': True,#False,#True,         
    },
    



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
        'solar_constant':4634724., 
        'R_stellar':476700000.,
        'd_stellar':1117512000.,
        'linear_tau':1.,#0.2, 
        'odp':1.0
    },


    'spectral_dynamics_nml': {
        'damping_order': 4, # DEFAULT: 2 
                            # Yields lap^8 damping 
        'water_correction_limit': 200.e2, # DEFAULT: 0. 
                                          # adds upper limit to water correction (which corrects small non-conservation
                                          # of water in the dynamical core)
        'reference_sea_level_press':1.0e5, # DEFAULT: 101325. 
                                           # used to construct hybrid coord and in implicit timestepping 
                                           # note actual mean sea level pressure is set in constants_mod 
        'num_levels':30, # Number of levels corresponding to set below 
        'valid_range_t':[1.,6000.], # just set this to be a wide temperature range 
        'initial_sphum':[0.], # DEFAULT: 0.0 
                                 # start the model with some water in atmosphere 
        'use_virtual_temperature':True,
        'vert_coord_option':'uneven_sigma', # DEFAULT: 'even_sigma' 
                                     # I have a set of hybrid levels I like to use, input below 
        'robert_coeff':0.03, # DEFAULT: 0.04, used in Robert filter for timestepping 
        # set to T42 resolution (default)
        'lon_max': 64,#128,#256, # DEFAULT: 128
        'lat_max': 32,#64,#128, # DEFAULT: 64
        'num_fourier': 21,#42,#85, # DEFAULT: 42
        'num_spherical': 22,#43,#86, # DEFAULT: 43
        
        ### NOTE (Mistake!): the 'input' levels I use are a little different from Tan et al.  
        ### To get levels similar to Tan et al., CHANGE 'input' to 'hybrid' (vert_coord_option) and it will use these:
        'surf_res': 0.05, # DEFAULT: 0.1 
        'exponent': 3., # DEFAULT: 2.5 
        'scale_heights': 5., # DEFAULT: 4.
        # which define p/p0 = exp[-scale_heights*(surf_res*x + (1-surf_res)*x^exponent)] with x evenly spaced on unit interval
        # the levels obtained with these parameters are essentially identical to Tan's... REMEMBER TO CHANGE sponge_pbottom also!
        
    }, 
    
    'spectral_init_cond_nml':{
        'initial_temperature':2000., 
    },
   

   
    
    
    # FMS Framework configuration -- haven't modified these from defaults present in all experiment scripts
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
    
    'constants_nml':{
        'omega': 2*np.pi / ((400./60./24.)*86400.),#7.2921150e-5 / 8., 
        'es0': 1.e-6,
        'rdgas':4000.,#4000.,#296.8,
        'kappa':4000./14000.,#4000./14000.,#296.8/1039,
        'radius':1.51*6.371e6,
        'grav':21.3, 
        'pstd':1.e5, 
        'pstd_mks':1.e6,
    },


    
    
})


if __name__=="__main__":
    

    ps_list = np.array([0.1,1,10]) 
    drag_list = np.array([1e2,1e3]) 
    
    for n_ps, ps in enumerate(ps_list):
        for n_drag, drag in enumerate(drag_list):
            
            if n_drag == 0:
                run_exp = exp.derive('k2-141b_'+str(ps)+'bar_H2_closein_TJ1e2')
            elif n_drag == 1:
                run_exp = exp.derive('k2-141b_'+str(ps)+'bar_H2_closein_TJ1e3')

            run_exp.namelist['constants_nml']['pstd'] = ps*1e5 
            run_exp.namelist['constants_nml']['pstd_mks'] = ps*1e6 
            run_exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = ps*1e5 
            run_exp.namelist['two_stream_gray_rad_nml']['ir_tau_eq'] = ps 
            run_exp.namelist['two_stream_gray_rad_nml']['ir_tau_pole'] = ps 
            run_exp.namelist['surface_flux_nml']['tj_coeff'] = drag*0.001 
            run_exp.namelist['surface_flux_nml']['tj_coeff_m'] = drag*0.001
            run_exp.namelist['diffusivity_nml']['tj_coeff'] = drag*0.001 
            run_exp.namelist['diffusivity_nml']['tj_coeff_m'] = drag*0.001
            run_exp.namelist['diffusivity_nml']['tj_bl_pres'] = 85000.*ps 
            run_exp.namelist['diffusivity_nml']['tj_strato_pres'] = 10000./2.*ps
            run_exp.namelist['damping_driver_nml']['sponge_pbottom'] = 260. * ps
            
            overwrite=False
            run_exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=overwrite)

            for i in range(2,51):
                run_exp.run(i, num_cores=NCORES, overwrite_data=overwrite)
        
        
                  
    
        
        

