import os

import numpy as np
from field_table_write import write_ft
from run_exp import create_exp_obj
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
import sys
import datetime
n_moments = 2
horizontal_resolution = "T170"
NCORES = 64

spinup_exp_name = "RT170_HYB_ADPT_" #"real_T85_sst_spinup"

delta_sst = 0
namelist = "namelist_basefile_T170"


dt = 360
original_dt = 360
min_dt = 60
n_start_month = 1
n_end_month = 361
restart = False
NCORES = 64                   # number of cores
current_n = n_start_month


exp = create_exp_obj(spinup_exp_name,delta_sst, n_moments ,namelist,
                horizontal_resolution, vertical_resolution=50,
                dt_atm=dt)


while current_n < n_end_month:
    try:
        for i in range(current_n , n_end_month):
            if current_n ==0: restart = False
            else: restart = True
            exp.run(i, num_cores=NCORES, overwrite_data=True, use_restart = restart)
            print(f"{datetime.datetime.now()} : month {current_n} completed at dt = {dt}",file=sys.stdout, flush=True)
            current_n+=1
    except:
        if dt/2 < min_dt:
            raise RuntimeError(f"dt too small ({dt/2}), cannot continue month {i}")   
        else:   
            print(f"{datetime.datetime.now()} : dt -> dt/2")
            dt = dt/2
            success = False
            exp.update_namelist({'main_nml' : {'dt_atmos' : dt}})
            print(f"{datetime.datetime.now()} : namelist updated to dt = {dt}",file=sys.stdout, flush=True)
            if current_n == n_start_month: # First month did not work
                pass 
            else: # any other month
                while success == False:
                    try:
                        # run for 2 months
                        for j in range(current_n,current_n + 2):
                            
                            exp.run(j, num_cores=NCORES, overwrite_data=True)
                            print(f"{datetime.datetime.now()} : month {i} completed at dt = {dt}",file=sys.stdout, flush=True)
                            current_n+=1
                    except:
                        if dt/2 < min_dt:
                            raise RuntimeError(f"dt too small ({dt/2}), cannot continue month {i}")  
                        print(f"{datetime.datetime.now()} : dt -> dt/2")
                        dt = dt/2
                dt = original_dt
                exp.update_namelist({'main_nml' : {'dt_atmos' : dt}})
                print(f"{datetime.datetime.now()} : namelist updated to dt = {dt}",file=sys.stdout, flush=True)

