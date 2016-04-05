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

#for land use this kind of indexing approach. Then it will be super easy.

#s parameters from Sauliere 2012 footnote 1. Do Rocky's first.
h_0 = 2670.
central_lon = 247.5
central_lat = 40.
L_1 = 7.5
L_2 = 20.
gamma_1 = 42.
gamma_2 = 42.

delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2


h_arr_rockys = h_0 * np.exp(-(delta_1**2. + delta_2**2.))

idx = (h_arr_rockys / h_0 > 0.05) #s make sure exponentials are cut at some point - use the value from p70 of Brayshaw's thesis. 



topo_array[idx] = h_arr_rockys[idx]

#s parameters from Sauliere 2012 footnote 1. Do `Tibet' second
h_0 = 5700.
central_lon = 82.5
central_lat = 28
L_1 = 12.5
L_2 = 12.5
gamma_1 = -49.5
gamma_2 = -18.

delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2


h_arr_tibet_no_amp = np.exp(-(delta_1**2.))*(1./delta_2)*np.exp(-0.5*(np.log(delta_2))**2.)

maxval = np.nanmax(h_arr_tibet_no_amp) #For some reason my maximum value of h_arr_tibet_no_amp > 1. Renormalise so h_0 sets amplitude. 

h_arr_tibet = (h_arr_tibet_no_amp/maxval)*h_0

idx = (h_arr_tibet / h_0 > 0.05) 


topo_array[idx] = h_arr_tibet[idx]

#Add back land from Sauliere too. 

idx_na = (103.-43./40.*(lon_array-180) < lat_array) & ((lon_array-180)*43./50. -51.8 < lat_array) &( lat_array < 60.)
idx_sa = (737.-7.2*(lon_array-180) < lat_array) & ((lon_array-180)*10./7. + -212.1 < lat_array) &( lat_array < -22./45*(lon_array-180) +65.9)
idx_ea = (17. <= lat_array) & (lat_array < 60.) & (-5. < lon_array) & ( 43./40.*lon_array -101.25 < lat_array)
#idx_ea = (17. <= lat_array) & (lat_array < 60.) & (-5. < lon_array) & (137. > lon_array)
idx_ea_neg = (17. <= lat_array) & (lat_array < 60.) & (355. < lon_array) 
idx_af = (lat_array < 17.) & (-52./27.*lon_array + 7.37 < lat_array) & (52./38.*lon_array -65.1 < lat_array) 
idx_af_neg = (lat_array < 17.) & (-52./27.*(lon_array-360) + 7.37 < lat_array)
idx_world = idx_na + idx_sa + idx_ea + idx_ea_neg + idx_af + idx_af_neg
land_array[idx_world] = 1.0


idx = (land_array == 0.) & (topo_array != 0.)

topo_array[idx] = 0. 




topo_file = Dataset('land_world_mountains.nc', 'w', format='NETCDF3_CLASSIC')
lat = topo_file.createDimension('lat', nlat)
lon = topo_file.createDimension('lon', nlon)

latitudes = topo_file.createVariable('lat','f4',('lat',))
longitudes = topo_file.createVariable('lon','f4',('lon',))
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
m.contour(xi,yi,land_array)
cs = m.contourf(xi,yi,topo_array, cmap=plt.get_cmap('RdBu_r'))
plt.xticks(np.linspace(0,360,13))
plt.yticks(np.linspace(-90,90,7))
cb = plt.colorbar(cs, shrink=0.5, extend='both')

plt.show()

