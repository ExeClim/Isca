import os
import numpy as np
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

dt_atmos = 10
radius = 1.602
surface_g = 18.1/9.81
orbital_period = 3.092926
S0 = 176.2
kappa = 2./7. # pure n2
mole_mass = 28e-3 # kg/mol
rdgas = 8.314/mole_mass
surface_pressure = 0.01e5

NCORES = 16
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = IscaCodeBase.from_directory(GFDL_BASE)
cb.compile()
exp = Experiment('HD219134b_04', codebase=cb)
diag = DiagTable()

diag.add_file('atmos_20days', 20, 'days', time_units='days')
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('two_stream', 'tdt_rad', time_avg=True)
diag.add_field('two_stream', 'olr', time_avg=True)
diag.add_field('two_stream', 'swdn_sfc', time_avg=True)
diag.add_field('two_stream', 'flux_lw', time_avg=True)
diag.add_field('two_stream', 'flux_sw', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)


exp.diag_table = diag
exp.clear_rundir()

exp.namelist = namelist = Namelist({
    'main_nml':{
        'days'   : 20,
        'hours'  : 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos':dt_atmos,
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },
    'spectral_dynamics_nml':{
        'damping_option': 'resolution_dependent',
        'damping_order': 4,
        'damping_coeff': 9.259259e-04,
        'do_mass_correction': True,
        'do_energy_correction': True,
        'do_water_correction': True,
        'use_virtual_temperature': True,
        'vert_advect_uv': 'second_centered',
        'vert_advect_t': 'second_centered',
        'longitude_origin': 0.,
        'robert_coeff': .04,
        'alpha_implicit': .5,
        'lon_max': 128,
        'lat_max': 64,
        'num_levels': 30,
        'num_fourier': 42,
        'num_spherical': 43,
        'fourier_inc': 1,
        'triang_trunc': True,
        #'topography_option': 'flat',
        'valid_range_t': [0.0,3500.],
        'vert_coord_option': 'uneven_sigma',
        'surf_res': 0.1,
        'scale_heights': 5.0,
        'exponent': 2.0,
        'reference_sea_level_press':surface_pressure,
    },
    'fms_nml':{
        'domains_stack_size': 600000,
    },
    'fms_io_nml':{
        'threading_write': 'single',
        'fileset_write': 'single', 
    },
    'vert_turb_driver_nml':{
        'do_mellor_yamada': False,
        'do_shallow_conv': False,
        'gust_scheme': 'constant',
        'constant_gust ': 1.0,
        'use_tau': True,
        'do_molecular_diffusion': False,
    },
    'diag_manager_nml':{
        'mix_snapshot_average_fields': False,
        'max_input_fields': 400,
        'max_output_fields': 500,
    },
    'diffusivity_nml':{
        'pbl_mcm': False,
        'free_atm_diff': False,
        'entr_ratio': 0.0,
        'parcel_buoy': 0.0,
        'do_virtual_non_mcm': True,
        'fixed_depth': False,
        'frac_inner': 0.1,
    },
    'monin_obukhov_nml':{
        'neutral': False,
        'rich_crit': 2.0,
        'stable_option': 1,
    },
    'surface_flux_nml':{
        'use_virtual_temp': True,
    },
    #ruizhi add the namelist below
    # to set grey radiation, using frierson rad_scheme
    'two_stream_gray_rad_nml':{
        'rad_scheme': 'frierson',
        'solar_constant': S0*1370.,
        # set tauSW
        'atm_abs': 0.0, # to set tauSW always=0
        # set tauLW
        # lw_tau_0 = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*sin(lat)**2
        # set eq = pole ???
        # set ir_tau_eq to ???
        'ir_tau_eq': 1.0, # default: 6.0
        'ir_tau_pole': 1.0,  # default: 1.5 
        #*    For stratospheric converence, add a linear_tau = 0.1.
        #*    tau(p,ps) = 0.1*(p/ps) + tau0*(ps/1 bar) *(p/ps)**2 ,
        #*    where tau0 = {0.2,1,5}
        #(old) lw_tau(:,:,k) = lw_tau_0 * (linear_tau * p_half(:,:,k)/pstd_mks     &
        #  + (1.0 - linear_tau) * (p_half(:,:,k)/pstd_mks)**wv_exponent )
        #(new) lw_tau(:,:,k) = 0.1 * p_half(:,:,k)/surf_pres &
        #  + linear_tau*surf_pres/pstd_mks*(p_half(:,:,k)/surf_pres)**wv_exponent 
        'surf_pres': surface_pressure, # ruizhi add it
        'linear_tau': 1, # actually tau_0 in Daniel's eq, [0.2,1,5]
        'wv_exponent': 2.0, # 4.0 for Earth, Feng suggests 3.0 for dry planets
        'tidally_locked': True,
        'noon_longitude': 270.0,
    },
    # to set dry convection
    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,            
        'two_stream_gray': True,     #Use the grey radiation scheme
        'do_rrtm_radiation':  False, #Don't use RRTM
        'convection_scheme': 'DRY', #Use dry convection           
    },
    'dry_convection_nml': {
        'tau': dt_atmos,
        'small': 0.001,
    },
    # to set planet constants
    'constants_nml': {
        'radius': radius*6371.e3, # form?
        'grav': surface_g*9.81,
        'omega': 2.*np.pi/(orbital_period*24.*3600.), # [s^-1]
        'orbital_period': orbital_period*24.*3600., # [s]
        'solar_const': S0*1370., # [W/m^2]
        'rdgas': rdgas, # gas constant for N2 [J/kg/K]
        'kappa': kappa # R/c_p depends on the molecule
    },
    # to have mixed layer output ?
    'mixed_layer_nml': {
        'tconst' : 310.,                       # can it be higher?
        'prescribe_initial_dist':True,
        #'evaporation': False,                  # Disable surface evaporation
        'depth': 0.5,                          # Depth of mixed layer used
        'albedo_value': 0.0,                   # set to zero for future tests
    },
    'atmosphere_nml': {
        'idealized_moist_model': True
    },
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':False,
    },
    'sat_vapor_pres_nml': {
        'do_simple':True, # disable the calculation in the kernal code
        'do_not_calculate':True, #turn of esat calc altogether (for exoplanets where temperatures will be outside valid range)
    },
})

if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2,21):
        exp.run(i, num_cores=NCORES)