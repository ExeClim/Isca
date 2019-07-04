from util_scale import *
from util_output import *
import time
from constants import GFDL_TRACE
from isca import GFDL_WORK
import sh
import os


def main():
    max_cores, min_cores, codebase_name, res_list = parse_arguments()
    core_list = get_core_list(max_cores, min_cores)
    res_list = get_resolution_list(res_list)
    # cb, diag, namelist = setup_experiment(codebase_name)

    for i, ncores in enumerate(core_list):
        for j, resolution in enumerate(res_list):
            print(ncores, resolution)
            exp_name = f'trace_{codebase_name}_{ncores}_{resolution[0]}'
            # run_experiment(ncores, cb, diag, namelist, resolution, exp_name)


if __name__ == '__main__':
    main()
