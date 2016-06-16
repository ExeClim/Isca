import numpy as np
import sh

from exp19_tidally_locked_moist import baseexp, omega

NUM_CORES = 16

rot_ratios = [0.01, 0.1, 0.5, 1.0, 2.0]
mass_ratios = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]

if __name__ == '__main__':
    baseexp.log.info('Running rot_ratios %r and atmospheric masses %r' % (rot_ratios, mass_ratios))

    experiments = []
    for r in rot_ratios:
        for m in mass_ratios:
            exp = baseexp.derive('exp22_r%.1f_m%.1f' % (r, m))
            exp.screen_runmonth_prefix = 'r%.1f:m%.1f' % (r, m)

            #omega  = 2*np.pi / orbital_period
            p_omega = omega * r
            orbital_period = 2*np.pi / p_omega

            exp.namelist['constants_nml'] = {
                'omega': p_omega,
                'orbital_period': orbital_period
            }

            exp.namelist['spectral_dynamics_nml'] = {
                'reference_sea_level_press': 1.0e5*m,
            }

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

