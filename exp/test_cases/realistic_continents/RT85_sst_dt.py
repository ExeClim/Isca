import os

import numpy as np
from field_table_write import write_ft
from run_exp import start_sst_exp
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

n_moments = 2
horizontal_resolution = "T85"


spinup_exp_name = "RT85_sst_dt" #"real_T85_sst_spinup"
res_file = "/home/philbou/scratch/isca_data/RT85_sst_dt/restarts/res0235.tar.gz"
delta_sst = 0
namelist = "namelist_basefile_T85"
start_sst_exp(spinup_exp_name,delta_sst, n_moments,240 ,namelist,
                horizontal_resolution, vertical_resolution=45,
                use_restart = True,restart_file = res_file,n_start_months=236,
                dt_atm=360,sponge = 150)
    
    