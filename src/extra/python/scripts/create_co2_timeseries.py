# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num

#create grid
t_res = 42

resolution_file = Dataset('/scratch/sit204/FMS2013/GFDLmoistModel/land_file_generator/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')

lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]

lonbs = resolution_file.variables['lonb'][:]
latbs = resolution_file.variables['latb'][:]

#lons = [0.]
#lonbs = [0., 360.]

#lats = [-10., 0., 10.]
#latbs = [-90., -5., 5., 90.]

p_full = [950., 300.]

#create times
day_number = np.arange(0,10)

time_arr = day_number_to_date(day_number)

#Find grid and time numbers

#nlon=len(lons)
#nlat=len(lats)

#nlonb=len(lonbs)
#nlatb=len(latbs)

nlon=lons.shape[0]
nlat=lats.shape[0]

nlonb=lonbs.shape[0]
nlatb=latbs.shape[0]

npfull=len(p_full)
ntime=len(time_arr)


#create time series based on times

co2 = np.zeros((ntime, npfull, nlat, nlon))

for tick in np.arange(0,len(day_number)):
	co2[tick,...] = 300.*1.e-6*day_number[tick]


#Output it to a netcdf file. 

co2_file = Dataset('co2.nc', 'w', format='NETCDF3_CLASSIC')

lat = co2_file.createDimension('lat', nlat)
lon = co2_file.createDimension('lon', nlon)

latb = co2_file.createDimension('latb', nlatb)
lonb = co2_file.createDimension('lonb', nlonb)

pfull = co2_file.createDimension('pfull', npfull)
time = co2_file.createDimension('time', ntime)

latitudes = co2_file.createVariable('lat','d',('lat',))
longitudes = co2_file.createVariable('lon','d',('lon',))

latitudebs = co2_file.createVariable('latb','d',('latb',))
longitudebs = co2_file.createVariable('lonb','d',('lonb',))

pfulls = co2_file.createVariable('pfull','f4',('pfull',))
times = co2_file.createVariable('time','f4',('time',))

latitudes.units = 'degree'.encode('utf-8')
latitudes.cartesian_axis = 'Y'

longitudes.units = 'degree'.encode('utf-8')
longitudes.cartesian_axis = 'X'

latitudebs.units = 'degree'.encode('utf-8')
latitudebs.cartesian_axis = 'Y'

longitudebs.units = 'degree'.encode('utf-8')
longitudebs.cartesian_axis = 'X'


pfulls.units = 'hPa'
pfulls.cartesian_axis = 'Z'

times.units = 'days since 0000-01-01 00:00:00.0'
times.calendar = 'THIRTY_DAY_MONTHS'
times.calendar_type = 'THIRTY_DAY_MONTHS'
times.cartesian_axis = 'T'

co2_array_netcdf = co2_file.createVariable('co2','f4',('time','pfull','lat','lon',))

latitudes[:] = lats
longitudes[:] = lons

latitudebs[:] = latbs
longitudebs[:] = lonbs

pfulls[:]     = p_full
times[:]     = date2num(time_arr,units='days since 0001-01-01 00:00:00.0',calendar='360_day')

co2_array_netcdf[:] = co2

co2_file.close()







