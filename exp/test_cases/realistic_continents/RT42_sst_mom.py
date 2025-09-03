import os

import numpy as np
from field_table_write import write_ft
from run_exp import start_sst_exp
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

n_moments = 5
horizontal_resolution = "T42"


exp_name = "RT42_sst_0_5mom" #"real_T42_sst_spinup"

res_file = "/home/philbou/scratch/isca_data/RT42_sst_0_5mom/restar1ts/res0089.tar.gz"

delta_sst = 0
namelist = "namelist_basefile_T42"
start_sst_exp(exp_name,delta_sst, n_moments,120 ,namelist,
                horizontal_resolution, vertical_resolution=40,
                use_restart = True,restart_file = res_file,n_start_months=90)
    
    