# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num
import pdb
import create_timeseries as cts

resolution_file = Dataset('/scratch/sit204/ozone_cmip5_files/Ozone_CMIP5_ACC_SPARC_1990-1999_T3M_O3.nc', 'r')

lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]

ozone_in = resolution_file.variables['O3'][:]

zm_ozone_in = np.mean(ozone_in,axis=3)

lons = [0.]
lonb_temp = [0., 360.]

latb_temp=np.zeros(lats.shape[0]+1)

for tick in np.arange(1,lats.shape[0]):
    latb_temp[tick]=(lats[tick-1]+lats[tick])/2.

latb_temp[0]=-90.
latb_temp[-1]=90.

latbs=latb_temp
lonbs=lonb_temp

lats[0]=-90.+0.01
lats[-1]=90.-0.01


nlon=1#lons.shape[0]
nlat=lats.shape[0]

nlonb=len(lonbs)
nlatb=latbs.shape[0]

p_full = resolution_file.variables['plev'][:]

time_arr = resolution_file.variables['time'][:]

date_arr=np.mod((time_arr-time_arr[0]),12)
date_arr_new=np.zeros(12)

ozone_new=np.zeros((12,p_full.shape[0],nlat))

mol_to_kg_factor=(3.*16.)/(28.97)

#time average ozone by month
for month in np.arange(0,12):
    time_idx=date_arr==month
    ozone_new[month,:,:]=np.mean(zm_ozone_in[time_idx,:,:],axis=0)*(mol_to_kg_factor)
    date_arr_new[month]=np.mean(date_arr[time_idx])

p_half_temp=np.zeros(p_full.shape[0]+1)

p_half_temp[0]=1200.

for tick in np.arange(1,p_full.shape[0]):

    p_half_temp[tick]=(p_full[tick-1]+p_full[tick])/2.


p_half=p_half_temp

time_arr=(date_arr_new+0.5)*30.


#DO FLIPPING

p_full=p_full[::-1]
p_half=p_half[::-1]

ozone_new=ozone_new[:,::-1,:]


#Find grid and time numbers

nlon=len(lons)
nlat=len(lats)

nlonb=len(lonbs)
nlatb=len(latbs)

npfull=len(p_full)
nphalf=len(p_half)

ntime=len(time_arr)


#Output it to a netcdf file. 
file_name='ozone_1990_cmip5_test.nc'
variable_name='ozone_1990_cmip5'

number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['npfull']=npfull
number_dict['nphalf']=nphalf
number_dict['ntime']=ntime

time_units='days since 0000-01-01 00:00:00.0'

cts.output_to_file(ozone_new,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,variable_name,number_dict)

