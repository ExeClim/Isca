import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#from isca import GFDL_BASE # MP added
import os # MP added

# specify resolution
t_res = '1_deg'

#read in grid from approriate file
# resolution_file = Dataset('t'+str(t_res)+'_atmos_daily.nc', 'r', format='NETCDF3_CLASSIC') # MP commented

GFDL_BASE = os.environ['GFDL_BASE']
resolution_file = Dataset(os.path.join(GFDL_BASE,'input/sst_clim_amip.nc'), 'r', format='NETCDF3_CLASSIC') # MP added

# resolution_file = Dataset('/scratch/mp586/Isca/input/sst_clim_amip.nc', 'r', format='NETCDF3_CLASSIC') # MP added

lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]

lonsb = resolution_file.variables['lonb'][:]
latsb = resolution_file.variables['latb'][:]

nlon=lons.shape[0]
nlat=lats.shape[0]

nlonb=lonsb.shape[0]
nlatb=latsb.shape[0]

# MP commented
# output_file = Dataset('t'+str(t_res)+'.nc', 'w', format='NETCDF3_CLASSIC')

# MP added
output_file = Dataset(os.path.join(GFDL_BASE,'src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc'), 'w', format='NETCDF3_CLASSIC')


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
