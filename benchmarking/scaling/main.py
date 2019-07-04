from util_scale import *
from util_output import *
import time
from constants import GFDL_BENCH
import sh

if __name__ == "__main__":
    max_cores, min_cores, codebase_name, res_list = parse_arguments()
    core_list = get_core_list(max_cores, min_cores)
    res_list = get_resolution_list(res_list)
    cb, diag, namelist = setup_experiment(codebase_name)

    for i, ncores in enumerate(core_list):
        for j, resolution in enumerate(res_list):
            print(ncores, resolution)
            exp_name = f'{codebase_name}_{ncores}_{resolution[0]}'
            print(exp_name)
            start = time.time()
            run_experiment(ncores, cb, diag, namelist, resolution, exp_name)
            end = time.time()
            time_delta = end - start
            data = [ncores, resolution[0], 'Total', time_delta]
            write_to_csvfile(f'{GFDL_BENCH}/{exp_name}', data)
