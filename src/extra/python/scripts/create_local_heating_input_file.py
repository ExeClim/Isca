# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import date2num
import xarray as xar
import pdb
import create_timeseries as cts
from finite_difference import diffz

grav = 9.81
c_p = 287.04 / (2./7.)

#Read in radiation-scheme output file(s) Will need to use open_mfdataset if more than one.
resolution_file = xar.open_dataset('~/Desktop/full_continents_global_monthly_control_test.nc', decode_times=False)

lons = resolution_file['lon']
lats = resolution_file['lat']


#Set up output file variables. 

latb_temp=np.zeros(lats.shape[0]+1)

for tick in np.arange(1,lats.shape[0]):
    latb_temp[tick]=(lats[tick-1]+lats[tick])/2.

lonb_temp=np.zeros(lons.shape[0]+1)

for tick in np.arange(1,lons.shape[0]):
    lonb_temp[tick]=(lons[tick-1]+lons[tick])/2.

latb_temp[0]=-90.
latb_temp[-1]=90.

lonb_temp[0]=0.
lonb_temp[-1]=360.

latbs=latb_temp
lonbs=lonb_temp

nlon=lons.shape[0]
nlat=lats.shape[0]

nlonb=len(lonbs)
nlatb=latbs.shape[0]

p_full = resolution_file['pfull'][::-1]

p_half = resolution_file['phalf'][::-1]

time_arr = resolution_file['time'].values

nphalf =p_half.shape[0] 

heating_rate = np.zeros((nphalf, nlat, nlon))


# DEFINE A HEATING RATE HERE.


heating_rate_new = np.zeros((12, heating_rate.shape[1], nlat, nlon))

#Arrange 12 months to be 30 days apart
time_arr_new = np.arange(0,12,1)*30.
time_arr = time_arr_new


#Repeat 1 month 12 times just because I only have 1 input file. This would need replacing with different values each month.
for t in range(len(time_arr)):
    heating_rate_new[t,...] = heating_rate

heating_rate = heating_rate_new[:,:,...]

date_arr=np.mod((time_arr-time_arr[0]),12)
date_arr_new=np.zeros(12)

#Find grid and time numbers

npfull=p_full.shape[0]
nphalf=p_half.shape[0]

ntime=len(time_arr)

#Output it to a netcdf file. 
file_name='full_continents_heating_input_file.nc'
variable_name='heating_rate'

number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['npfull']=npfull
number_dict['nphalf']=nphalf
number_dict['ntime']=ntime

#Line that would need changing depending on year-on-year repeating heating or timeseries heating. 
time_units='days since 0000-01-01 00:00:00.0'

cts.output_to_file(heating_rate[:,::-1,...],lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,variable_name,number_dict, dims=total_flux.dims)

