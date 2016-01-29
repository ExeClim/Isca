import numpy as np

import mima

exp = mima.Experiment('playground', overwrite_data=True)


exp.diag_table_file = exp.diag_table_file+'.no_calendar'
#exp.clear_workdir()
exp.compile()

exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*30,     # approximately 1 earth month of integration
    'calendar': 'no_calendar'
}

exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True
exp.namelist['spectral_dynamics_nml']['num_levels'] = 25

exp.clear_rundir()


#omega = 7.2921150e-5
orbital_period = 360*86400.0

for ratio in [360.0, 180.0, 90.0, 45.0, 15.0, 1.0]:
    exp1 = mima.Experiment('ratio_%d' % ratio)
    omega  = 2*np.pi / orbital_period * ratio

    exp1.diag_table_file = exp.diag_table_file
    exp1.execdir = exp.execdir

    exp1.namelist = exp.namelist.copy()

    exp1.namelist['constants_nml'] = {
        'omega': omega,
        'orbital_period': orbital_period
    }

    exp1.runmonth(1, use_restart=False)
    for i in range(2, 60):
        exp1.runmonth(i)