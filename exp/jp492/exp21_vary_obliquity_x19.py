import numpy as np

from exp19_tidally_locked_moist import baseexp, omega

NUM_CORES = 16

ratios = [0.1, 0.5, 1.0, 2.0]
obls = [0.0, 5.0, 10.0, 20.0, 25.0]

if __name__ == '__main__':
    baseexp.log.info('Running ratios %r and obliquities %r' % (ratios, obls))

    experiments = []
    for r in ratios:
        for o in obls:
            exp = baseexp.derive('exp21_r%.1f_o%.1f' % (r, o))
            exp.screen_runmonth_prefix = 'r%.1f:o%.1f' % (r, o)

            #omega  = 2*np.pi / orbital_period
            p_omega = omega * r
            orbital_period = 2*np.pi / p_omega

            exp.namelist['constants_nml'] = {
                'omega': p_omega,
                'orbital_period': orbital_period
            }

            exp.namelist['astronomy_nml'] = {
                'ecc': 0.0,
                'obliq': o
            }

            experiments.append(exp)

for exp in experiments:
    exp.clear_rundir()
    exp.run(1, use_restart=False, num_cores=NUM_CORES)

for i in range(2, 11):
    for exp in experiments:
        exp.run(i, num_cores=NUM_CORES)

