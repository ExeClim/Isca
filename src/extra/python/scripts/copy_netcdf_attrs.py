'''Function to copy netcdf attributes to a new files'''

import numpy as np
import xarray as xr
import os
from netcdftime import utime
from netCDF4 import Dataset, date2num
import pdb


def copy_netcdf_attrs(dsin, dsout, copy_vars = True):
    
    #Copy dimensions
    for dname, the_dim in list(dsin.dimensions.items()):
        print((dname, len(the_dim)))
        if dname!='Time':
            dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
        else:
            dsout.createDimension(dname, None)            

    if copy_vars:
        # Copy variables
        for v_name, varin in list(dsin.variables.items()):
            outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
            print((varin.datatype))
    
            # Copy variable attributes
            for k in varin.ncattrs():
                if k!='_FillValue':
                    outVar.setncatts({k: varin.getncattr(k)})
#                 else:
#                     outVar.setncatts({k: 1.0E20})                
    
            outVar[:] = varin[:]
    
    return dsout #returns without closing
    
if __name__=="__main__":

    dsin= Dataset('./spectral_dynamics_old_mod_213.res.nc', 'a', format='NETCDF3_CLASSIC')
    dsout= Dataset('./spectral_dynamics_old_rmod_213.res.nc', 'w', format='NETCDF3_CLASSIC')    

    copy_netcdf_attrs(dsin, dsout, copy_vars = True)

    dsout.close()