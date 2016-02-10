# This script can be run to check consistency after making changes to the source.
# It performs the following steps:
#
# 1. Compiles the current state of src
# 2. Runs with a 360 day calendar (grey rad, seasonal)
# 3. Runs with a 360 day calendar (grey rad, nonseasonal)
# 4. Runs with no calendar
#
# Other tests to be implemented

import os
import sys

import numpy as np

from gfdl.experiment import Experiment, DiagTable


diag = DiagTable()

diag.add_file('6hourly', 6*60*60, 'seconds')

diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')


# 1. Compile the code
testbase = Experiment('test', overwrite_data=True)

# # get the directory of this file
# srcdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'src')
# testbase.srcdir = srcdir

testbase.compile()



# 2. Run with calendar
caltest = Experiment('caltest', overwrite_data=True)
caltest.clear_workdir()

caltest.execdir = testbase.execdir

caltest.namelist['main_nml'].update({
    'dt_atmos': 900,
    'days': 2
})

caltest.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
caltest.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
caltest.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True
caltest.namelist['spectral_dynamics_nml']['num_levels'] = 25

caltest.use_diag_table(diag)

caltest.runmonth(1, use_restart=False)
# TODO: assert here that netcdf files created are valid


# 3. Run with calendar, no seasonal cycle
caltest.clear_rundir()

caltest.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False

caltest.runmonth(2, use_restart=False)

#TODO: assert here that month 2 is different to month 1

caltest.rm_workdir()



# 4. Runs with NO_CALENDAR
nocaltest = Experiment('nocaltest', overwrite_data=True)
nocaltest.clear_workdir()

nocaltest.execdir = testbase.execdir

nocaltest.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*2,
    'calendar': 'no_calendar'
}



nocaltest.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
nocaltest.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
nocaltest.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True
nocaltest.namelist['spectral_dynamics_nml']['num_levels'] = 25

nocaltest.use_diag_table(diag)

nocaltest.runmonth(1, use_restart=False)


nocaltest.rm_workdir()
