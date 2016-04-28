import numpy as np

from gfdl.experiment import Experiment, DiagTable

exp = Experiment('playground', overwrite_data=True)


diag = DiagTable()

diag.add_file('6hourly', 6*60*60, 'seconds', time_units='days')
diag.add_file('hourly', 1, 'hours', time_units='days')

diag.add_field('dynamics', 'ps')
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')
diag.add_field('dynamics', 'sphum')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')
diag.add_field('two_stream', 'lw_dtrans')
# diag.add_field('two_stream', 'lw_dtrans_win')
# diag.add_field('two_stream', 'sw_dtrans')

diag.add_field('mixed_layer', 't_surf')

diag.add_field('dry_convection', 'dp')
diag.add_field('dry_convection', 'CAPE')
diag.add_field('dry_convection', 'CIN')
diag.add_field('dry_convection',  'LZB')
diag.add_field('dry_convection', 'LCL')
diag.add_field('dry_convection', 'dt_tg')
diag.add_field('dry_convection', 'parcel_temp')



exp.use_diag_table(diag)

#exp.clear_workdir()
exp.disable_rrtm()
exp.compile()

exp.clear_rundir()

exp.namelist['main_nml'] = {
    'dt_atmos': 900,
#    'seconds': 5*86400,
    'days': 5,
    'calendar': 'no_calendar',
#    'current_date': [0001,1,1,0,0,0]
}

exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True
exp.namelist['spectral_dynamics_nml']['num_levels'] = 25
exp.namelist['idealized_moist_phys_nml']['convection_scheme'] = 'betts_miller'
exp.namelist['idealized_moist_phys_nml']['lwet_convection'] = False
exp.namelist['idealized_moist_phys_nml']['do_bm'] = False

exp.namelist['dry_convection_nml'] = {
    'tau': 86400.0*10,
    'gamma': 1.0, # K/km
}

exp.namelist['lscale_cond_nml'] = {
    'do_simple': True,
    'do_evap': False,
    'hc': 1.0,
}

exp.namelist['mixed_layer_nml'] = {
    'albedo_value': 0.27,
    'depth': 10.0,
    #'prescribe_initial_dist': True
    # 'tconst': 285.0,
    # 'delta_T': 40.0,
    'evaporation': False,
    'do_qflux': False
}

exp.runmonth(1, use_restart=False)

#for i, scheme in enumerate(('frierson', 'geen', 'byrne')):
    # exp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = scheme
    # exp.runmonth(i, use_restart=False)