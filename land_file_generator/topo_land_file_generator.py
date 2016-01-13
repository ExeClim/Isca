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

topo_array = np.zeros((nlat,nlon))
land_array = np.zeros((nlat,nlon))

#make 2d arrays of latitude and longitude
lon_array, lat_array = np.meshgrid(lons,lats)
lonb_array, latb_array = np.meshgrid(lonb,latb)



#s for land use this kind of indexing approach. Then it will be super easy.

central_lat = 40. 
central_lon = 40.
radius_degrees = 20.
std_dev = 10.
height = 3500.


rsqd_array = np.sqrt((lon_array - central_lon)**2.+(lat_array - central_lat)**2.)

idx = (rsqd_array < radius_degrees) 

topo_array[idx] = height* np.exp(-(rsqd_array[idx]**2.)/(2.*std_dev**2.))
land_array[idx] = 1.0


topo_file = Dataset('land.nc', 'w', format='NETCDF3_CLASSIC')
lat = topo_file.createDimension('lat', nlat)
lon = topo_file.createDimension('lon', nlon)

latitudes = topo_file.createVariable('latitude','f4',('lat',))
longitudes = topo_file.createVariable('longitude','f4',('lon',))
topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
land_array_netcdf = topo_file.createVariable('land_mask','f4',('lat','lon',))

latitudes[:] = lats
longitudes[:] = lons
topo_array_netcdf[:] = topo_array
land_array_netcdf[:] = land_array

topo_file.close()

lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lon_array, lat_array)


plt.figure()
m.contour(xi,yi,topo_array)
cs = m.contourf(xi,yi,topo_array, cmap=plt.get_cmap('RdBu_r'))
plt.xticks(np.linspace(0,360,13))
plt.yticks(np.linspace(-90,90,7))
cb = plt.colorbar(cs, shrink=0.5, extend='both')

plt.show()

