# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num
import pdb

nfiles=50

sst_all=np.zeros((nfiles,12,180,360))

for file_tick in range(nfiles):

	filename='amipbc_sst_360x180_19'+str(file_tick+50)


	resolution_file = Dataset('/scratch/sit204/sst_amip_files/1950_1999/'+filename+'.nc', 'r')

	lons = resolution_file.variables['longitude'][:]
	lats = resolution_file.variables['latitude'][:]

	sst_in = resolution_file.variables['tosbcs'][:]

	sst_all[file_tick,:,:,:]=sst_in

sst_in=np.mean(sst_all,axis=0)


lonbs = resolution_file.variables['bounds_longitude'][:]
latbs = resolution_file.variables['bounds_latitude'][:]

nlon=lons.shape[0]
nlat=lats.shape[0]

nlonb=lonbs.shape[0]
nlatb=latbs.shape[0]

lonbs_adjusted=np.zeros(nlonb+1)
latbs_adjusted=np.zeros(nlatb+1)

lonbs_adjusted[0:nlonb]=lonbs[:,0]
lonbs_adjusted[nlonb]=lonbs[-1,1]

latbs_adjusted[0:nlatb]=latbs[:,0]
latbs_adjusted[nlatb]=latbs[-1,1]

day_number = resolution_file.variables['time'][:]

time_arr = day_number_to_date(day_number, calendar_type = 'gregorian', units_in = 'days since 1870-1-1')

time_arr_adj=np.arange(15,360,30)

p_full=None
p_half=None

npfull=None
nphalf=None

#Find grid and time numbers

ntime=time_arr.day.shape[0]

#Output it to a netcdf file. 
file_name='sst_clim_amip_test.nc'
variable_name='sst_clim_amip'

number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['npfull']=npfull
number_dict['nphalf']=nphalf
number_dict['ntime']=ntime

time_units='days since 0000-01-01 00:00:00.0'

cts.output_to_file(sst_in,lats,lons,latbs_adjusted,lonbs_adjusted,p_full,p_half,time_arr_adj,time_units,file_name,variable_name,number_dict)


