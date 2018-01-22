# adapted from https://gist.github.com/guziy/8543562

from netCDF4 import Dataset
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from scipy import interpolate
import os 
GFDL_BASE = os.environ['GFDL_BASE']


# landfile=Dataset(os.path.join(GFDL_BASE,'input/two_continents/land.nc'),mode='r')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land.nc'),mode='r')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/sqland_plus_antarctica/land.nc'),mode='r')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/aquaplanet/land.nc'),mode='r')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/square_South_America/land.nc'))
landfile=Dataset(os.path.join(GFDL_BASE,'input/all_continents/land.nc'))


landmask=landfile.variables['land_mask'][:]

#input file
dsin = Dataset('/scratch/mp586/Isca_DATA/flat_continents_simple2mbucket_qflux/restarts/extract_res241/mixed_layer.res.nc')

#output file
dsout = Dataset('/scratch/mp586/Isca_DATA/flat_continents_simple2mbucket_qflux/restarts/tocompile_241/mixed_layer.res.nc', "w", format="NETCDF3_CLASSIC") # This is to differentiate from _sqland, because that was NH SH symmetric in each month, not only in the annual mean! --> bug!

#Copy dimensions
for dname, the_dim in dsin.dimensions.iteritems():
    print dname, len(the_dim)
    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)


# Copy variables
for v_name, varin in dsin.variables.iteritems():
    
    if v_name == 'bucket_depth':
        outVar = dsout.createVariable('bucket_depth', varin.datatype, varin.dimensions)
        print varin.datatype
    
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})

        bdepth_0 = varin[0,0,:,:]
        bdepth_1 = varin[0,1,:,:]

        bdepth_0[landmask == 1.] = 2.00000 # 2.0
        bdepth_1[landmask == 1.] = 2.00000 # 2.0

        
        outVar[0,0,:,:] = np.asarray(bdepth_0)
        outVar[0,1,:,:] = np.asarray(bdepth_1)

    else:
        
        outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
        print varin.datatype
    
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        outVar[:] = varin[:]

# close the output file
dsout.close()








