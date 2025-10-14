import os

import numpy as np
from field_table_write import write_ft
from run_exp import create_exp_obj
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
import datetime
import sys 

def update_nml(exp,dt):
    exp.update_namelist({'main_nml' : {'dt_atmos' : dt}})
    print(f"{datetime.datetime.now()} : namelist updated to dt = {dt}",file=sys.stdout, flush=True)

n_moments = 2
horizontal_resolution = "T42"
NCORES = 32
delta_sst = int(sys.argv[1])
n_start_month = int(sys.argv[2])

if delta_sst < 0:
        delta_sst_name = f"m{abs(delta_sst)}"
else:
        delta_sst_name = str(delta_sst)



exp_name = f"R{horizontal_resolution}_sst_{delta_sst_name}_bucket" #"real_T42_sst_spinup"

dt = 720
original_dt = 720
min_dt = 30

n_end_month = 361
if n_start_month != 1:  
    restart = True
else: restart= False
current_n = n_start_month + 1
n_end = 361


spinup_exp_name = exp_name
res_file = f"/home/philbou/scratch/isca_data/{spinup_exp_name}/restarts/res{n_start_month-1:04d}.tar.gz"

namelist = f"namelist_basefile_{horizontal_resolution}"
exp = create_exp_obj(spinup_exp_name,delta_sst, n_moments ,namelist,
                horizontal_resolution, vertical_resolution=45,
                dt_atm=dt)
"""
exp.update_namelist({
    'mixed_layer_nml': {  
        'albedo_choice' : 1, 
        'albedo_cntr_lat' : 25,
        'albedo_wdth_lat' : 10,
        'albedo_cntr_lon' : 5,
        'albedo_wdth_lon' :25,
        'higher_albedo' : 0.4,
        'albedo_value' : 0.25}})
"""
exp.run(n_start_month, num_cores=NCORES, overwrite_data=True,use_restart=restart, restart_file=res_file)
print(f"{datetime.datetime.now()} : month {n_start_month} completed at dt = {dt}",file=sys.stdout, flush=True)

while current_n < n_end_month:
    try:
        for i in range(current_n , n_end_month):
            if current_n ==0: restart = False
            else: restart = True
            exp.run(i, num_cores=NCORES, overwrite_data=True)
            print(f"{datetime.datetime.now()} : month {current_n} completed at dt = {dt}",file=sys.stdout, flush=True)
            current_n+=1
    except:
        if dt/2 < min_dt:
            raise RuntimeError(f"dt too small ({dt/2}), cannot continue month {i}")   
        else:   
            dt = dt/2
            print(f"{datetime.datetime.now()} : dt -> dt/2",file=sys.stdout, flush=True)
            update_nml(exp,dt)
            
            # Run for 3 months at the smaller dt
            success = False
            # if first month did not work
            if current_n == n_start_month: 
                raise RuntimeError(f"Crashed at first month, PICK A SMALLER DT")  # Let it crash, you fool
            # any other month
            else: 
                success = False
                while success == False:
                    try:
                        # run for 2 months
                        for j in range(current_n,current_n + 2):
                            exp.run(j, num_cores=NCORES, overwrite_data=True)
                            print(f"{datetime.datetime.now()} : month {current_n} completed at dt = {dt}",file=sys.stdout, flush=True)
                            current_n+=1
                        success = True
                    except:
                        if dt/2 < min_dt:
                            raise RuntimeError(f"dt too small ({dt/2}), cannot continue month {i}")  
                        print(f"{datetime.datetime.now()} : dt -> dt/2")
                        dt = dt/2
                        update_nml(exp,dt)
                print(f"{datetime.datetime.now()} : dt/2 -> OG dt")
                dt = original_dt
                update_nml(exp,dt)


