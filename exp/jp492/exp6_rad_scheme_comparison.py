import numpy as np
import sh

from gfdl.experiment import Experiment, DiagTable, P

exp = Experiment('exp6_rad_scheme_comparison')

diag = DiagTable()

diag.add_file('daily', 1, 'days', time_units='days')
#diag.add_file('hourly', 1, 'hours', time_units='days')

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

diag.add_field('atmosphere', 'rh')

exp.use_diag_table(diag)

exp.compile()
exp.clear_rundir()

exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*100,
    'calendar': 'no_calendar'
}

exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False

exp.namelist['spectral_dynamics_nml']['num_levels'] = 25

base_datadir = exp.datadir

for scheme in ('frierson', 'geen', 'byrne'):
    exp.datadir = P(base_datadir, scheme)
    exp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = scheme
    exp.runmonth(0, use_restart=False)
    # copy the data to the base data directory and rename to the scheme
    sh.cp(P(base_datadir, scheme, 'run0', 'daily.nc'), P(base_datadir, '%s.nc' % scheme))
