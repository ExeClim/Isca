# LHS3844b
# pure co2 atmosphere
# radius: 1.303 R_earth
# surface gravity: 1.325 g_earth
# orbital period: 0.4629 days
# solar constant: 70.3053 S_earth
### Planet configuration
radius = 1.303
surface_g = 1.325
orbital_period = 0.4629
S0 = 70.3053
eccentricity = 0.000
surface_pressure = 1.0e5 # [Pa] reference_sea_level_press
import os
import numpy as np

from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress

NCORES = 16
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = SocratesCodeBase.from_directory(GFDL_BASE)

exp = Experiment('lhs3844b_co2_test2', codebase=cb)
exp.clear_rundir()

# exp.inputfiles = [os.path.join(base_dir,'input/co2.nc')] # maybe the distribution of co2 is not important?

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Write out diagnostics need for vertical interpolation post-processing
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf')

#Tell model which diagnostics to write
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True) #SH
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True) #LH
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)

#temperature tendency - units are K/s
diag.add_field('socrates', 'soc_tdt_lw', time_avg=True) # net flux lw 3d (up - down)
diag.add_field('socrates', 'soc_tdt_sw', time_avg=True)
diag.add_field('socrates', 'soc_tdt_rad', time_avg=True) #sum of the sw and lw heating rates

#net (up) and down surface fluxes
diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_lw_down', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw_down', time_avg=True)
diag.add_field('socrates', 'soc_spectral_olr', time_avg=True)
#net (up) TOA and downard fluxes
diag.add_field('socrates', 'soc_olr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw', time_avg=True) 
diag.add_field('socrates', 'soc_toa_sw_down', time_avg=True)
diag.add_field('socrates', 'soc_flux_lw', time_avg=True)
diag.add_field('socrates', 'soc_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_co2', time_avg=True)

exp.diag_table = diag
sp_test_root = 'src/atmos_param/socrates/src/trunk/data/spectra/sp_test/'
star_gas = 'L/CO2'

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 30,  # check for hot planet
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':60,  # default: 600
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    'socrates_rad_nml': {
        'stellar_constant':S0*1370,
        'lw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/sp_test/sp_lw_b202_Trappist1_CO2_UCL_T41xP22_001_nk20'),
        'sw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/sp_test/sp_sw_b202_LHS3844_CO2_UCL_T41xP22_001_nk20'),
        'dt_rad':300,    # test!!! not fixeds
        'store_intermediate_rad':True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
        'tidally_locked':True,
        'solday':90,
        'inc_o3': False, # only include CO2
        'inc_h2o': False,
        'account_for_effect_of_water': False, # still water in the atm, not just disable the radiation; ask stephen: h2o; check the rain-> small value; email metoffice again
        'account_for_effect_of_ozone': False,
        # 'co2_ppmv': 1e6, # pure CO2
        #'do_read_co2': True,
        #'co2_file_name':'co2',
        #'co2_field_name':'co2' # CO2 change with time
    },
    
    'spectral_init_cond_nml': {
        'initial_temperature': 2000,
    },
    
    'astronomy_nml': {
        'ecc': eccentricity,
        'obliq': 0.0
    },
    
    'constants_nml': {
        'radius': radius*6371.e3, # form?
        'grav': surface_g*9.81,
        'omega': 2.*np.pi/(orbital_period*24.*3600.), # [s^-1]
        'orbital_period': orbital_period*24.*3600., # [s]
        'solar_const': S0*1370., # [W/m^2]
        'rdgas': 188.95, # gas constant for CO2 [J/kg/K]
        'kappa': 2./9. # R/c_p depends on the molecule
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
        'convection_scheme': 'DRY', #Use dry convection           
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
        'old_dtaudv': True,
        'use_actual_surface_temperatures': False, 
    },

    'atmosphere_nml': {
        'idealized_moist_model': True, 
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation': False,                  # Disable surface evaporation
        'depth': 0.5,                          # Depth of mixed layer used
        'albedo_value': 0.38,                  # Albedo value used      
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   # can't be more than 373K ?
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':False,
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True, # disable the calculation in the kernal code
        'do_not_calculate':True, #turn of esat calc altogether (for exoplanets where temperatures will be outside valid range)
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
        'damping_coeff': 2.3148148e-05,  # *10?        
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':surface_pressure,
        'num_levels':40,      #How many model pressure levels to use
        'valid_range_t':[10.,3500.],
        'initial_sphum':[0.], # set to zero
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, # Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 5.0,
        'exponent':2.0,
        'robert_coeff':0.03,
        'do_water_correction': False # disable water correction
    },
    
    'dry_convection_nml': {
        'tau': 120,
        'small': 0.001,
    },

})

#Lets do a run!
if __name__=="__main__":
    cb.compile()
    exp.run(4, num_cores=NCORES, overwrite_data=False)
    overwrite = False
    for i in range(5,13):
        exp.run(i, num_cores=NCORES, overwrite_data=overwrite)
