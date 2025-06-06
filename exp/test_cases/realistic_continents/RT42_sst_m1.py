import os

import numpy as np
from field_table_write import write_ft
from run_exp import start_sst_exp
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

n_moments = 2
horizontal_resolution = "T42"


exp_name = "RT42_sst_m1" #"real_T42_sst_spinup"

res_file = "/home/philbou/scratch/isca_data/RT42_sst_m2/restarts/res0162.tar.gz"

delta_sst = -1
namelist = "namelist_basefile_T42"
start_sst_exp(exp_name,delta_sst, n_moments,360 ,namelist,
                horizontal_resolution, vertical_resolution=40,
                use_restart = False,restart_file = res_file,n_start_months=1)
    
    