import os

import numpy as np
from field_table_write import write_ft
from run_exp import start_sst_exp
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

n_moments = 2
horizontal_resolution = "T85"


exp_name = "RTEST85_sst_2" #"real_T85_sst_spinup"

res_file = "/home/philbou/scratch/isca_data/RTEST85_sst_spinup/restarts/res0002.tar.gz"

delta_sst = 2
namelist = "namelist_basefile_T85"
start_sst_exp(exp_name,delta_sst, n_moments,2 ,namelist,
                horizontal_resolution, vertical_resolution=45,
                use_restart = True,restart_file = res_file,n_start_months=1)
    
    