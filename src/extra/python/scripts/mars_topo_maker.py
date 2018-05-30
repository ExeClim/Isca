import xarray as xar
import gauss_grid as gg
import numpy as np
import mpl_toolkits.basemap as basemap
import matplotlib.pyplot as plt
from netCDF4 import Dataset

num_lat_out = 64
num_lon_out = 128


mola_original = xar.open_dataset('/scratch/sit204/planet_topo_files/mars/mola32.nc')

longitudes_out  = np.arange(0., 360., (360./num_lon_out))
latitudes_out  = gg.gaussian_latitudes(int(num_lat_out/2.))[0]  

lon_array, lat_array = np.meshgrid(longitudes_out, latitudes_out)

output_array = basemap.interp(mola_original.alt.values[::-1,:], mola_original.longitude.values, mola_original.latitude.values[::-1], lon_array, lat_array, order=1)


nlat = latitudes_out.shape[0]
nlon = longitudes_out.shape[0]

#Write land and topography arrays to file
topo_filename = './t42_mola_mars.nc'
topo_file = Dataset(topo_filename, 'w', format='NETCDF3_CLASSIC')
lat = topo_file.createDimension('lat', nlat)
lon = topo_file.createDimension('lon', nlon)
latitudes = topo_file.createVariable('lat','f4',('lat',))
longitudes = topo_file.createVariable('lon','f4',('lon',))
topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
latitudes[:] = latitudes_out
longitudes[:] = longitudes_out
topo_array_netcdf[:] = output_array
topo_file.close()

