# An earth-like setup
# Obliquity

import sys
import numpy as np

import gfdl.experiment

import diagtables
import namelists


P = gfdl.experiment.P

args = sys.argv[1:]
ratios = [float(arg) for arg in args]

omega = 7.2921150e-5
earth_orbital_period = 365.25*86400.0
NUM_CORES = 16

#ratios = [360, 180, 90, 45, 30, 15, 4, 2, 1]

# Use one experiment to compile the source.  All other experiments
# then use the same code with different namelist config
baseexp = gfdl.experiment.Experiment('base_dry0.1',
    repo='git@github.com:jamesp/GFDLmoistModel.git',
    commit='dry0.1')

baseexp.log.info('Running ratios %r' % ratios)
baseexp.disable_rrtm()
baseexp.compile()

baseexp.namelist = namelists.moist.copy()

baseexp.namelist['spectral_dynamics_nml']['valid_range_t'] =  [100., 800.]

baseexp.namelist['main_nml'] = {
    'dt_atmos': 450,
    'seconds': 86400.0*100,
    'calendar': 'no_calendar'
}

baseexp.namelist['astronomy_nml'] = {
    'ecc': 0.0,
    'obliq': 23.4
}

baseexp.namelist['mixed_layer_nml']['depth'] = 25.0
baseexp.namelist['mixed_layer_nml']['do_qflux'] = False

baseexp.use_diag_table(diagtables.moist)

if __name__ == '__main__':
    for ratio in ratios:
        exp = baseexp.derive('exp20_ratio_%.3f' % ratio)
        exp.clear_rundir()
        exp.screen_runmonth_prefix = 'r%.3f' % ratio

        #omega  = 2*np.pi / orbital_period
        p_omega = omega * ratio

        exp.namelist['constants_nml'] = {
            'omega': p_omega,
            'orbital_period': earth_orbital_period
        }

        exp.runmonth(1, use_restart=False, num_cores=NUM_CORES)
        for i in range(2, 21):  # 81, 161
            exp.runmonth(i, num_cores=NUM_CORES)
