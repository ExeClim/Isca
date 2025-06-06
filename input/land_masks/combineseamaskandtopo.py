import xarray as xar
filename = "ERA5_land_sea_mask.nc"
ds_lsm = xar.open_dataset(filename)
filename = "ERA5_geopotential.nc"
ds_z = xar.open_dataset(filename)
g = 9.80665
zsurf = ds_z.z/g
lsm = ds_lsm.lsm

combined_set = xar.Dataset({"land_mask": lsm, "zsurf" :zsurf.isel(valid_time = 0)})
combined_set.to_netcdf("ERA5_land.nc")