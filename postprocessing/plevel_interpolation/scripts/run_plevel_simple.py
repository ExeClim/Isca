from netCDF4 import Dataset  
from plevel_fn import plevel_call, daily_average, join_files, two_daily_average, monthly_average
import sys
import os
import time
import pdb
import subprocess

start_time=time.time()
base_dir='/home/sthomson/datadir_isca/'
exp_name_list = ['project_7_ice_albedo_False_t42']
avg_or_daily_list=['monthly']
start_file=1
end_file=7
nfiles=(end_file-start_file)+1

level_set='standard' #Default is the standard levels used previously. ssw_diagnostics are the ones blanca requested for MiMa validation
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 

try:
    out_dir
except:
    out_dir = base_dir

plevs={}
var_names={}

if level_set=='standard':

    plevs['monthly']=' -p "1000 10000 25000 50000 85000 92500"'
    plevs['6hourly']=' -p "1000 10000 25000 50000 85000 92500"'
    plevs['daily']  =' -p "1000 10000 25000 50000 85000 92500"'
    
    var_names['monthly']='-a slp height'
    var_names['6hourly']='-a slp height'
    var_names['daily']='-a slp height'
    file_suffix='_interp'

for exp_name in exp_name_list:
    for n in range(nfiles):
        for avg_or_daily in avg_or_daily_list:
            print(n+start_file)

            nc_file_in = base_dir+'/'+exp_name+'/run%04d'%(n+start_file)+'/atmos_'+avg_or_daily+'.nc'
            nc_file_out = out_dir+'/'+exp_name+'/run%04d'%(n+start_file)+'/atmos_'+avg_or_daily+file_suffix+'.nc'

            if not os.path.isfile(nc_file_out):
                plevel_call(nc_file_in,nc_file_out, var_names = var_names[avg_or_daily], p_levels = plevs[avg_or_daily], mask_below_surface_option=mask_below_surface_set)

print('execution time', time.time()-start_time)

