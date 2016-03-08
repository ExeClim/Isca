# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num

#create grid
lons = [0.]
lonbs = [0., 360.]

lats = [0.]
latbs = [-90., 90.]

p_full = [1000.]

#create times
#day_number = np.arange(0,10)

day_number = [16.0,45.0,75.0,105.0,136.0,166.0,197.0,228.0,258.0,289.0,319.0,350.0]

time_arr = day_number_to_date(day_number)

#Find grid and time numbers

nlon=len(lons)
nlat=len(lats)

nlonb=len(lonbs)
nlatb=len(latbs)

npfull=len(p_full)
ntime=len(time_arr)


#create time series based on times

co2 = np.zeros((ntime, npfull, nlat, nlon))

for tick in np.arange(0,len(day_number)):
    co2[tick,...] = 300.*(day_number[tick]/360.+1.)*1.e-6 #Some scenario in dimensionless units. 1.e-6 is to convert from ppmv. 


#Output it to a netcdf file. 

co2_file = Dataset('co2.nc', 'w', format='NETCDF3_CLASSIC')

lat = co2_file.createDimension('lat', nlat)
lon = co2_file.createDimension('lon', nlon)

latb = co2_file.createDimension('latb', nlatb)
lonb = co2_file.createDimension('lonb', nlonb)

pfull = co2_file.createDimension('pfull', npfull)
time = co2_file.createDimension('time', 0) #s Key point for a year-independent climatology file is to have the length of the time axis 0, or 'unlimited'. This seems necessary to get the code to run properly. 

latitudes = co2_file.createVariable('lat','d',('lat',))
longitudes = co2_file.createVariable('lon','d',('lon',))

latitudebs = co2_file.createVariable('latb','d',('latb',))
longitudebs = co2_file.createVariable('lonb','d',('lonb',))

pfulls = co2_file.createVariable('pfull','d',('pfull',))
times = co2_file.createVariable('time','d',('time',))

latitudes.units = 'degrees_N'.encode('utf-8')
latitudes.cartesian_axis = 'Y'
latitudes.edges = 'latb'
latitudes.long_name = 'latitude'

longitudes.units = 'degrees_E'.encode('utf-8')
longitudes.cartesian_axis = 'X'
longitudes.edges = 'lonb'
longitudes.long_name = 'longitude'

latitudebs.units = 'degrees_N'.encode('utf-8')
latitudebs.cartesian_axis = 'Y'
latitudebs.long_name = 'latitude edges'

longitudebs.units = 'degrees_E'.encode('utf-8')
longitudebs.cartesian_axis = 'X'
longitudebs.long_name = 'longitude edges'

pfulls.units = 'hPa'
pfulls.cartesian_axis = 'Z'
pfulls.positive = 'down'
pfulls.long_name = 'full pressure level'

times.units = 'days since 0000-01-01 00:00:00.0'
times.calendar = 'THIRTY_DAY_MONTHS'
times.calendar_type = 'THIRTY_DAY_MONTHS'
times.cartesian_axis = 'T'
times.climatology = '1979-01-01 00:00:00, 1998-01-01 00:00:00'

co2_array_netcdf = co2_file.createVariable('co2','f4',('time','pfull','lat','lon',))

latitudes[:] = lats
longitudes[:] = lons

latitudebs[:] = latbs
longitudebs[:] = lonbs

pfulls[:]     = p_full
times[:]     = date2num(time_arr,units='days since 0001-01-01 00:00:00.0',calendar='360_day')

co2_array_netcdf[:] = co2

co2_file.close()

