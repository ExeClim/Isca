import os

import numpy as np
from field_table_write import write_ft
from run_exp import start_sst_exp
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

n_moments = 2
horizontal_resolution = "T85"


spinup_exp_name = "RTEST85_sst_spinup" #"real_T85_sst_spinup"
res_file = " "
delta_sst = 0
namelist = "namelist_basefile_T85"
start_sst_exp(spinup_exp_name,delta_sst, n_moments,1 ,namelist,
                horizontal_resolution, vertical_resolution=45,
                use_restart = False,restart_file = res_file,n_start_months=1)
    
    