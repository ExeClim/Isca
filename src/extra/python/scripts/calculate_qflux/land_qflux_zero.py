import numpy as np
import xarray
from xarray import ufuncs as xruf
import time
from scipy import stats
from mpl_toolkits.basemap import shiftgrid
import matplotlib.pyplot as plt
import os
from netCDF4 import Dataset


GFDL_BASE = os.environ['GFDL_BASE']


landfile=Dataset(os.path.join(GFDL_BASE,'input/two_continents/land.nc'),mode='r')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land.nc'),mode='r')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/sqland_plus_antarctica/land.nc'),mode='r')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/aquaplanet/land.nc'),mode='r')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/square_South_America/land.nc'))
# landfile=Dataset(os.path.join(GFDL_BASE,'input/all_continents/land.nc'))

landmask=landfile.variables['land_mask'][:]

#input file
dsin = Dataset(os.path.join(GFDL_BASE,'input/aquaplanet/ocean_qflux.nc'))

#output file
dsout = Dataset(os.path.join(GFDL_BASE,'input/two_continents/ocean_qflux.nc'), "w", format="NETCDF3_CLASSIC")

#Copy dimensions
for dname, the_dim in dsin.dimensions.iteritems():
    print dname, len(the_dim)
    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)


# Copy variables
for v_name, varin in dsin.variables.iteritems():
    
    if v_name == 'ocean_qflux':
        outVar = dsout.createVariable('ocean_qflux', varin.datatype, varin.dimensions)
        print varin.datatype
    
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})

        qflux_in = varin[:]
        qflux_out = qflux_in[:]

        for i in range(0,12):
            qflux_i = qflux_in[i,:,:]
            qflux_i[landmask == 1.] = 0. # 2.0
            qflux_out[i,:,:] = qflux_i
            print(i)
            print(np.mean(qflux_i))
        
        outVar[:] = np.asarray(qflux_out)

    else:
        
        outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
        print varin.datatype
    
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        outVar[:] = varin[:]

# close the output file
dsout.close()







