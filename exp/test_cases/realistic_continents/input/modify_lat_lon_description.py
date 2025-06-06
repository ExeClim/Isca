import xarray as xar
path = "/home/philbou/Isca/exp/test_cases/realistic_continents/input/era_land_t42.nc"

ds = xar.open_mfdataset(path, decode_times=False,chunks=None)  
ds['lat'].attrs['units'] = 'degrees_N'
ds['lat'].attrs['cartesian_axis'] = 'Y'
ds['lat'].attrs['edges'] = 'latb'
ds['lat'].attrs['long_name'] = 'latitude'

ds['lon'].attrs['units'] = 'degrees_E'
ds['lon'].attrs['cartesian_axis'] = 'X'
ds['lon'].attrs['edges'] = 'lonb'
ds['lon'].attrs['long_name'] = 'longitude'

ds.to_netcdf("/home/philbou/Isca/exp/test_cases/realistic_continents/input/mod_era_land_t42.nc")