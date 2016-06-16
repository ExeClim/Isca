import gfdl
import namelists, diagtables

import numpy as np
import sh

from exp22_atmos_mass_x19 import baseexp, mass_ratios

NUM_CORES = 16

#baseexp = gfdl.experiment.Experiment('playground')


#baseexp.disable_rrtm()
#baseexp.compile()

#baseexp.namelist = namelists.moist.copy()

#baseexp.namelist['spectral_dynamics_nml']['valid_range_t'] =  [100., 800.]

# baseexp.namelist['main_nml'] = {
#     'dt_atmos': 900,
#     'seconds': 86400.0*100,
#     'calendar': 'no_calendar'
# }

# baseexp.namelist['astronomy_nml'] = {
#     'ecc': 0.0,
#     'obliq': 0.0
# }

# baseexp.namelist['mixed_layer_nml']['depth'] = 25.0
# baseexp.namelist['mixed_layer_nml']['do_qflux'] = False

# baseexp.use_diag_table(diagtables.moist)



if __name__ == '__main__':
    baseexp.log.info('Running atmospheric masses %r' % (mass_ratios))

    experiments = []
    for m in mass_ratios:
        exp = baseexp.derive('exp23_m%.1f' % (m))
        exp.screen_runmonth_prefix = 'm%.1f' % (m)

        exp.diag_table.add_field('two_stream', 'lw_dtrans')
        exp.diag_table.add_field('two_stream', 'sw_dtrans')
        exp.diag_table.add_field('two_stream', 'tdt_solar')
        exp.diag_table.add_field('two_stream', 'flux_rad')

        # standard Earth rotation and orbital rate
        omega = 7.29e-5
        orbital_period = 365.25*86400.0

        exp.namelist['constants_nml'] = {
            'omega': omega,
            'orbital_period': orbital_period
        }

        exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = 1.0e5*m

        # no shortwave absorption in the atmosphere
        exp.namelist['two_stream_gray_rad_nml']['atm_abs'] = 0.0

        experiments.append(exp)

    for exp in experiments:
        exp.clear_rundir()
        try:
            exp.run(1, use_restart=False, num_cores=NUM_CORES)
        except:
            exp.log.error('Experiment run ended prematurely.')

    for i in range(2, 11):
        for exp in experiments:
            try:
                exp.run(i, num_cores=NUM_CORES)
            except:
                exp.log.error('Experiment run ended prematurely.')