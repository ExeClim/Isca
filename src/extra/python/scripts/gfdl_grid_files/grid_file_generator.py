import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# specify resolution
t_res = 42

#read in grid from approriate file
resolution_file = Dataset('t'+str(t_res)+'_atmos_daily.nc', 'r', format='NETCDF3_CLASSIC')

lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]

lonsb = resolution_file.variables['lonb'][:]
latsb = resolution_file.variables['latb'][:]

nlon=lons.shape[0]
nlat=lats.shape[0]

nlonb=lonsb.shape[0]
nlatb=latsb.shape[0]


output_file = Dataset('t'+str(t_res)+'.nc', 'w', format='NETCDF3_CLASSIC')

lat = output_file.createDimension('lat', nlat)
lon = output_file.createDimension('lon', nlon)
latb = output_file.createDimension('latb', nlatb)
lonb = output_file.createDimension('lonb', nlonb)

latitudes = output_file.createVariable('lat','f4',('lat',))
longitudes = output_file.createVariable('lon','f4',('lon',))
latitudesb = output_file.createVariable('latb','f4',('latb',))
longitudesb = output_file.createVariable('lonb','f4',('lonb',))

latitudes[:] = lats
longitudes[:] = lons
latitudesb[:] = latsb
longitudesb[:] = lonsb


output_file.close()
