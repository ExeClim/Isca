# Moist Earth
# - 360 day orbital period
# - Rotation rate: 361 rotations / year

import numpy as np
import exputil

exp = exputil.new_experiment('exp101_std_earth',
        tag='dry0.1', namelist='basic', diagtable='basic')

secs_per_day = 86400.0
days_per_year = 360.0
revs_per_year = days_per_year + 1
omega = 2*np.pi/secs_per_day*(revs_per_year/days_per_year)

exp.update_namelist({
    'constants_nml': {
        'omega': omega,
        'orbital_period': secs_per_day*days_per_year,
    },
    'two_stream_gray_rad_nml': {
        'atm_abs': 0.0,     # no SW absorption in the atmosphere
    },
    'mixed_layer_nml': {
        'albedo_value': 0.31,     #  needs to be higher to account for no atm_abs
        'depth': 40.0,
        'prescribe_initial_dist': False,  # use constant temp from init cond
    },
    'spectral_init_cond_nml' :{
        'initial_temperature': 285.0
    },
    'astronomy_nml': {
        'ecc': 0.0,
        'obliq': 23.4,
        'per': 0.0        # begin at autumn equinox
    },

    'main_nml': {
        'dt_atmos': 900,
        'seconds': 30*secs_per_day,
        'calendar': 'no_calendar'
    }
})

if __name__ == '__main__':
    exp.clear_rundir()
    exp.run(1, use_restart=False, num_cores=16)
    for i in range(1, 3*5):
        exp.run(i+1, num_cores=16)