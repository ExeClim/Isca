import os

import numpy as np
from field_table_write import write_ft
from run_exp import create_exp_obj
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
import sys
n_moments = 2
horizontal_resolution = "T85"
NCORES = 64
sp = int(sys.argv[1])
spinup_exp_name = f"RT85_HYB_50_TR1_{sp}" #"real_T85_sst_spinup"
res_file = "/home/philbou/scratch/isca_data/RT85_sst_sponge_200/restarts/res0074.tar.gz"
delta_sst = 0
namelist = "namelist_basefile_T85"

    
exp = create_exp_obj(spinup_exp_name,delta_sst, n_moments ,namelist,
                horizontal_resolution, vertical_resolution=50,dt_atm=720,
                trayfric=-1,sponge=sp)

N_restart = 0
n_end = 361
exp.run(N_restart + 1,
        num_cores=NCORES, overwrite_data=True,
        use_restart=False, restart_file= res_file)

for i in range(N_restart + 2, n_end):
            exp.run(i, num_cores=NCORES, overwrite_data=True)