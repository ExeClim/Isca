import numpy as np

from isca import IscaCodeBase, GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress
from util import get_core_list
import sh
import os
import namelist

# from ntfy import notify

NCORES = 16
CORE_LIST = get_core_list()

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = GreyCodeBase.from_directory(GFDL_BASE)
cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create a diagnostics output file for daily snapshots
diag = DiagTable()
diag.add_file('atmos_daily', 88440, 'seconds', time_units='days')

# Define diag table entries
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk', time_avg=True)
diag.add_field('dynamics', 'pk', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)

diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('atmosphere', 'precipitation', time_avg=True)

diag.add_field('atmosphere', 'dt_tg_convection', time_avg=True)

diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)

diag.add_field('two_stream', 'mars_solar_long', time_avg=True)
diag.add_field('two_stream', 'coszen', time_avg=True)
diag.add_field('two_stream', 'rrsun', time_avg=True)
diag.add_field('two_stream', 'swdn_toa', time_avg=True)
diag.add_field('two_stream', 'time_since_ae', time_avg=True)
diag.add_field('two_stream', 'true_anomaly', time_avg=True)
diag.add_field('two_stream', 'dec', time_avg=True)
diag.add_field('two_stream', 'ang', time_avg=True)



if __name__ == "__main__":

    conv_schemes = ['none']

    depths = [2.]

    pers = [70.85]

    scale = 1.

    for conv in conv_schemes:
        for depth_val in depths:
            for per_value in pers:
                for core in CORE_LIST:
                    exp = Experiment('grey_mars_mk36_per_value' + str((per_value)) + '_' + conv + '_mld_' + str(
                        depth_val) + '_cores__' + core,
                                     codebase=cb)
                    exp.clear_rundir()
                    exp.diag_table = diag
                    exp.namelist = namelist.copy()
                    exp.namelist['constants_nml']['grav'] = scale * 3.71
                    exp.namelist['constants_nml']['pstd'] = scale * 6100.0
                    exp.namelist['constants_nml']['pstd_mks'] = scale * 610.0
                    exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = scale * 610.0
                    exp.namelist['idealized_moist_phys_nml']['convection_scheme'] = conv
                    exp.namelist['mixed_layer_nml']['depth'] = depth_val
                    exp.namelist['astronomy_nml']['per'] = per_value

                    # with exp_progress(exp, description='o%.0f d{day}' % scale):
                    exp.run(1, use_restart=False, num_cores=NCORES)
                    for i in range(2, 241):
                        #                with exp_progress(exp, description='o%.0f d{day}' % scale):
                        exp.run(i, num_cores=NCORES)
    #                notify('top down with conv scheme = '+conv+' has completed', 'isca')
