# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num

#create grid
#lons = [0.]
#lonbs = [0., 360.]

#lats = [0.]
#latbs = [-90., 90.]

t_res=42

resolution_file = Dataset('/scratch/sit204/FMS2013/GFDLmoistModel/land_file_generator/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')

lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]

lonbs = resolution_file.variables['lonb'][:]
latbs = resolution_file.variables['latb'][:]

nlon=lons.shape[0]
nlat=lats.shape[0]

nlonb=lonbs.shape[0]
nlatb=latbs.shape[0]


#p_full = [1000.]

p_full = [0.017223, 0.078654, 0.180749, 0.375337, 0.71314, 1.253787, 2.060447, 3.193908, 4.707791, 6.646065, 9.04322, 11.926729, 15.321008, 19.251999, 23.751706, 28.862245, 34.639281, 41.154864, 48.49985, 56.786085, 66.14858, 76.747823, 88.772421, 102.442154, 118.011585, 135.77429, 156.067832, 179.279569, 205.853424, 236.297753, 271.194477, 311.209677, 357.105858, 409.756139, 470.160652, 539.465462, 618.984358, 710.223896, 814.912138, 935.031579, 1100.]

p_half = [0., 0.046817, 0.115544, 0.254986, 0.510238, 0.937468, 1.599324, 2.558913, 3.874132, 5.593891, 7.756986, 10.393581, 13.528609, 17.186225, 21.394476, 26.18965, 31.62002, 37.748946, 44.657459, 52.446517, 61.23914, 71.182614, 82.450941, 95.247663, 109.809174, 126.408616, 145.360464, 167.025894, 191.81904, 220.214288, 252.754733, 290.062007, 332.84765, 381.926291, 438.230874, 502.830259, 576.949493, 661.993148, 759.57212, 871.53435, 1000., 1200.]


#create times
#day_number = np.arange(0,10)

day_number = [3.0,45.0,75.0,105.0,136.0,166.0,197.0,228.0,258.0,289.0,319.0,350.0]

time_arr = day_number_to_date(day_number)

#Find grid and time numbers

nlon=len(lons)
nlat=len(lats)

nlonb=len(lonbs)
nlatb=len(latbs)

npfull=len(p_full)
nphalf=len(p_half)

ntime=len(time_arr)


#create time series based on times

co2 = np.zeros((ntime, npfull, nlat, nlon))

for tick in np.arange(0,len(day_number)):
    co2[tick,...] = 300.*(day_number[tick]/360.+1.) #Some scenario in dimensionless units. 1.e-6 is to convert from ppmv. 
#    co2[tick,...] = 300. #Some scenario in dimensionless units. 1.e-6 is to convert from ppmv. 

#Output it to a netcdf file. 

co2_file = Dataset('co2.nc', 'w', format='NETCDF3_CLASSIC')

lat = co2_file.createDimension('lat', nlat)
lon = co2_file.createDimension('lon', nlon)

latb = co2_file.createDimension('latb', nlatb)
lonb = co2_file.createDimension('lonb', nlonb)

pfull = co2_file.createDimension('pfull', npfull)
phalf = co2_file.createDimension('phalf', nphalf)
time = co2_file.createDimension('time', 0) #s Key point for a year-independent climatology file is to have the length of the time axis 0, or 'unlimited'. This seems necessary to get the code to run properly. 

latitudes = co2_file.createVariable('lat','d',('lat',))
longitudes = co2_file.createVariable('lon','d',('lon',))

latitudebs = co2_file.createVariable('latb','d',('latb',))
longitudebs = co2_file.createVariable('lonb','d',('lonb',))

pfulls = co2_file.createVariable('pfull','d',('pfull',))
phalfs = co2_file.createVariable('phalf','d',('phalf',))

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

phalfs.units = 'hPa'
phalfs.cartesian_axis = 'Z'
phalfs.positive = 'down'
phalfs.long_name = 'half pressure level'


times.units = 'days since 0000-01-01 00:00:00.0'
times.calendar = 'THIRTY_DAY_MONTHS'
times.calendar_type = 'THIRTY_DAY_MONTHS'
times.cartesian_axis = 'T'
times.climatology = '1979-01-01 00:00:00, 1998-01-01 00:00:00'

co2_array_netcdf = co2_file.createVariable('co2','f4',('time','pfull', 'lat','lon',))

latitudes[:] = lats
longitudes[:] = lons

latitudebs[:] = latbs
longitudebs[:] = lonbs

pfulls[:]     = p_full
phalfs[:]     = p_half

times[:]     = date2num(time_arr,units='days since 0001-01-01 00:00:00.0',calendar='360_day')

co2_array_netcdf[:] = co2

co2_file.close()

