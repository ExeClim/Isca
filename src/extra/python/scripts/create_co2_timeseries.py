# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num
import sys
import create_timeseries as cts

#create grid
manual_grid_option=False

lons,lats,lonbs,latbs,nlon,nlat,nlonb,nlatb=cts.create_grid(manual_grid_option)

p_full,p_half,npfull,nphalf=cts.create_pressures()

#create times
is_climatology=False
num_years=10	

time_arr,day_number,ntime,time_units=cts.create_time_arr(num_years,is_climatology)

#create time series based on times
co2 = np.zeros((ntime, npfull, nlat, nlon))

for tick in np.arange(0,len(day_number)):
    co2[tick,...] = 300.*(day_number[tick]/360.+1.) #Some scenario in dimensionless units. 1.e-6 is to convert from ppmv. 

#Output it to a netcdf file. 
file_name='co2_test.nc'
variable_name='co2'

number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['npfull']=npfull
number_dict['nphalf']=nphalf
number_dict['ntime']=ntime

cts.output_to_file(co2,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,variable_name,number_dict)



