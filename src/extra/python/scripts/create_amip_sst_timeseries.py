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



#DO FLIPPING




#Find grid and time numbers

ntime=time_arr.day.shape[0]

#Output it to a netcdf file. 

sst_file = Dataset('sst_clim_amip.nc', 'w', format='NETCDF3_CLASSIC')

lat = sst_file.createDimension('lat', nlat)
lon = sst_file.createDimension('lon', nlon)

latb = sst_file.createDimension('latb', nlatb+1)
lonb = sst_file.createDimension('lonb', nlonb+1)

time = sst_file.createDimension('time', 0) #s Key point for a year-independent climatology file is to have the length of the time axis 0, or 'unlimited'. This seems necessary to get the code to run properly. 

latitudes = sst_file.createVariable('lat','d',('lat',))
longitudes = sst_file.createVariable('lon','d',('lon',))

latitudebs = sst_file.createVariable('latb','d',('latb',))
longitudebs = sst_file.createVariable('lonb','d',('lonb',))

times = sst_file.createVariable('time','d',('time',))

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


times.units = 'days since 0000-01-01 00:00:00.0'
times.calendar = 'THIRTY_DAY_MONTHS'
times.calendar_type = 'THIRTY_DAY_MONTHS'
times.cartesian_axis = 'T'

sst_array_netcdf = sst_file.createVariable('sst_clim_amip','f4',('time', 'lat','lon',))

latitudes[:] = lats
longitudes[:] = lons

latitudebs[:] = latbs_adjusted
longitudebs[:] = lonbs_adjusted

times[:]     = time_arr_adj #date2num(time_arr,units='days since 0001-01-01 00:00:00.0',calendar='360_day')

sst_array_netcdf[:] = sst_in[:,:,:]

sst_file.close()

