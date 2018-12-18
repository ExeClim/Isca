import numpy as np
import xarray
from xarray import ufuncs as xruf
import time
from scipy import stats
from mpl_toolkits.basemap import shiftgrid
import matplotlib.pyplot as plt
import os
from netCDF4 import Dataset
import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
# import plotting_routines
from plotting_routines_kav7 import * # isca and gfdl have 0:04 and 0:03 

GFDL_BASE = os.environ['GFDL_BASE']
GFDL_DATA = os.environ['GFDL_DATA']

model = input('Enter model ')
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
exp_name= input('Enter data directory name as string ')
runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number ')
testdir = model_data + '/' + exp_name

landmask_name = input('Which landmask? ')
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landmask_name+'/land.nc'),mode='r')
landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]

# for specified lats

landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

# get sst climatology from chosen experiment

[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')

i = 1
if model=='isca':
    runnr="{0:04}".format(i)
elif model=='gfdl':
	runnr="{0:03}".format(i)
filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
nc = Dataset(filename,mode='r')




#input file
dsin = Dataset(os.path.join(GFDL_BASE,'input/sst_clim_amip.nc'))

#output file
# dsout = Dataset(os.path.join(GFDL_BASE,'input/'+landmask_name+'/prescribed_ssts_control.nc'), "w", format="NETCDF3_CLASSIC")
# for control, the input data directory was Isca/two_continents_newbucket_finalIscaAPqflux_landqfluxzero_zerointegral_with6hrly run 121-481
# for perturbed, the input data directory was Isca/two_continents_newbucket_finalIscaAPqflux_landqfluxzero_zerointegral_with6hrly_2xCO2_spinup_361 run 120-480

dsout = Dataset(os.path.join(GFDL_BASE,'input/'+landmask_name+'/ssts_from_twoCs_0qflux.nc'), "w", format="NETCDF3_CLASSIC")

#Copy dimensions
for dname, the_dim in dsin.dimensions.iteritems():
    if (dname == 'lat') or (dname == 'lon') or (dname == 'latb') or (dname == 'lonb'):
        print dname
        dsout.createDimension(dname, len(nc.variables[dname][:]))
    else:
        print dname, len(the_dim)
        dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)


# Copy variables
for v_name, varin in dsin.variables.iteritems():
    
    if v_name == 'sst_clim_amip':
#        outVar = dsout.createVariable('prescribed_ssts_control', varin.datatype, varin.dimensions)
        outVar = dsout.createVariable('ssts_from_twoCs_0qflux', varin.datatype, varin.dimensions)
        print v_name, varin.datatype
    
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})

        tsurf_in = varin[:]
        tsurf_out = tsurf_in[:]

        tsurf_out = tsurf_month_avg[:]
        
        outVar[:] = np.asarray(tsurf_out)

    elif v_name in ['lat','lon','latb','lonb']:

        outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
        print v_name, varin.datatype
    
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        
        outVar[:] = np.asarray(nc.variables[v_name][:])        

    else:
        outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
        print v_name, varin.datatype
    
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        outVar[:] = varin[:]

# close the output file
dsout.close()
