import numpy as np

import mima

exp = mima.Experiment('playground', overwrite_data=True)

diag = mima.DiagTable()

diag.add_file('6hourly', 6*60*60, 'seconds', time_units='days')
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')

exp.use_diag_table(diag)

exp.compile()

exp.clear_rundir()

exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*30,     # approximately 1 earth month of integration
    'calendar': 'no_calendar'
}

exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True
exp.namelist['spectral_dynamics_nml']['num_levels'] = 25

exp.runmonth(1, use_restart=False)