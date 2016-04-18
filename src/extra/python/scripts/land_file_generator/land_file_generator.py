import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# specify resolution
t_res = 42

#read in grid from approriate file
resolution_file = Dataset('../gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')

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
#idx = (20. < lat_array) & (lat_array < 60) & (30. < lon_array) & (lon_array < 130)
#land_array[idx] = 1.0

#make continents

idx_na = (103.-43./40.*(lon_array-180) < lat_array) & ((lon_array-180)*43./50. -51.8 < lat_array) &( lat_array < 60.)
idx_sa = (737.-7.2*(lon_array-180) < lat_array) & ((lon_array-180)*10./7. + -212.1 < lat_array) &( lat_array < -22./45*(lon_array-180) +65.9)
idx_ea = (17. <= lat_array) & (lat_array < 60.) & (-5. < lon_array) & ( 43./40.*lon_array -101.25 < lat_array)
idx_ea_neg = (17. <= lat_array) & (lat_array < 60.) & (355. < lon_array) 
idx_af = (lat_array < 17.) & (-52./27.*lon_array + 7.37 < lat_array) & (52./38.*lon_array -65.1 < lat_array) 
idx_af_neg = (lat_array < 17.) & (-52./27.*(lon_array-360) + 7.37 < lat_array)
idx_world = idx_na + idx_sa + idx_ea + idx_ea_neg + idx_af + idx_af_neg
land_array[idx_world] = 1.0



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


