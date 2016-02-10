import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# specify resolution
t_res = 42

#read in grid from approriate file
resolution_file = Dataset('gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')

lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]

lonb = resolution_file.variables['lonb'][:]
latb = resolution_file.variables['latb'][:]

nlon=lons.shape[0]
nlat=lats.shape[0]

land_array = np.zeros((nlat,nlon))


#make 2d arrays of latitude and longitude
lon_array, lat_array = np.meshgrid(lons,lats)
lonb_array, latb_array = np.meshgrid(lonb,latb)



#s for land use this kind of indexing approach. Then it will be super easy.
idx = (20. < lat_array) & (lat_array < 60) & (20. < lon_array) & (lon_array < 60)

land_array[idx] = 1.0



land_file = Dataset('land.nc', 'w', format='NETCDF3_CLASSIC')
lat = land_file.createDimension('lat', nlat)
lon = land_file.createDimension('lon', nlon)

latitudes = land_file.createVariable('latitude','f4',('lat',))
longitudes = land_file.createVariable('longitude','f4',('lon',))
land_array_netcdf = land_file.createVariable('land_mask','f4',('lat','lon',))

latitudes[:] = lats
longitudes[:] = lons
land_array_netcdf[:] = land_array

land_file.close()


lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lon_array, lat_array)

plt.figure()
m.contour(xi,yi,land_array)
cs = m.contourf(xi,yi,land_array, cmap=plt.get_cmap('RdBu_r'))
plt.xticks(np.linspace(0,360,13))
plt.yticks(np.linspace(-90,90,7))
cb = plt.colorbar(cs, shrink=0.5, extend='both')

plt.show()

