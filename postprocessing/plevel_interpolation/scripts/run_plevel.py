from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess



start_time=time.time()
base_dir='/scratch/sit204/Data_2013/'
exp_name_list=['giant_drag_exp_chai_values_1']
avg_or_daily_list=['monthly']
start_file=1100
end_file=1200
nfiles=(end_file-start_file)+1

plevs={}

#plevs['monthly']=' -p "3 16 51 138 324 676 1000 1266 2162 3407 5014 6957 9185 10000 11627 14210 16864 19534 20000 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 70000 73915 78095 82510 85000 87175 92104 97312"'

plevs['monthly']=' -p "3680 14720 24830 34880 44910 54920 64940 74940 84950 94960 100000 104960 114960 124970 134970 144970 154970 164970 174980 184980 194980 204980 214980 224980 234980 244980 254980 264980 274980 284990 294990"'

#plevs['6hourly']=' -p "25000 85000 92500"'
#plevs['6hourly']=' -p "1000 10000"'

var_names={}

var_names['monthly']='-a'
#var_names['6hourly']='ucomp ps height vor t_surf vcomp omega'
var_names['6hourly']='ucomp vcomp temp'


for exp_name in exp_name_list:
	for n in range(nfiles):
		for avg_or_daily in avg_or_daily_list:
			print n+start_file
			nc_file_in = base_dir+'/'+exp_name+'/run'+str(n+start_file)+'/atmos_'+avg_or_daily+'.nc'
			nc_file_out = base_dir+'/'+exp_name+'/run'+str(n+start_file)+'/atmos_'+avg_or_daily+'_interp_new_model_lev.nc'

			if not os.path.isfile(nc_file_out):
				plevel_call(nc_file_in,nc_file_out, var_names = var_names[avg_or_daily], p_levels = plevs[avg_or_daily])

print 'execution time', time.time()-start_time


