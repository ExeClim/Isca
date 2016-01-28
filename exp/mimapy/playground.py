import numpy as np

import mima

exp = mima.Experiment('playground', overwrite_data=True)


exp.diag_table_file = exp.diag_table_file+'.no_calendar'
#exp.clear_workdir()
exp.compile()


exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*2,     # approximately 1 earth month of integration
    'calendar': 'no_calendar'
}

omega = 7.2921150e-5
exp.namelist['constants_nml'] = {
    'grav': 10.0,
    'omega': omega,
}


exp.clear_rundir()

exp.runmonth(1, use_restart=False)

exp.namelist['constants_nml'].update({
    'omega': 7.2921150e-5*2.0
})

exp.runmonth(2, use_restart=False)

exp.namelist['constants_nml'] = {
    'grav': 10.0,
    'omega': omega,
    'orbital_period': 2*np.pi / omega
}

exp.runmonth(3, use_restart=False)