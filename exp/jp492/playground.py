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
diag.add_field('two_stream', 'lw_dtrans_win')
diag.add_field('two_stream', 'sw_dtrans')

diag.add_field('mixed_layer', 't_surf')



exp.use_diag_table(diag)

exp.compile()

exp.clear_rundir()

exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*10,
    'calendar': 'no_calendar'
}

exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
# turn on Ruth's radiation scheme:
exp.namelist['two_stream_gray_rad_nml']['solar_exponent'] = -1
exp.namelist['two_stream_gray_rad_nml']['wv_exponent'] = -1
exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False
exp.namelist['spectral_dynamics_nml']['num_levels'] = 25

exp.runmonth(1, use_restart=False)