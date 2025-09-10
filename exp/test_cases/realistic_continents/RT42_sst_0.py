import os

import numpy as np
from field_table_write import write_ft
from run_exp import create_exp_obj
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

n_moments = 2
horizontal_resolution = "T42"
NCORES = 32

exp_name = "RT42_sst_0_bucket" #"real_T42_sst_spinup"


N_restart = 0
res_file = f"/home/philbou/scratch/isca_data/RT42_sst_0_bucket/restarts/res0{N_restart}.tar.gz"

delta_sst = 0
n_end = 361
namelist = "namelist_basefile_T42"
exp = create_exp_obj(exp_name,delta_sst, n_moments ,namelist,
                horizontal_resolution, vertical_resolution=45)


exp.run(N_restart + 1,
        num_cores=NCORES, overwrite_data=True,
        use_restart=False, restart_file= res_file)

for i in range(N_restart + 2, n_end):
            exp.run(i, num_cores=NCORES, overwrite_data=True)