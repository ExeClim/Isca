import xarray as xar
import mpl_toolkits.basemap as basemap
import gauss_grid as gg
import numpy as np
from netCDF4 import Dataset
import pdb

def linear_interpolate_for_regrid(lon_list_in_grid, lat_list_in_grid, lon_list_out_grid, lat_list_out_grid, input_array):

    lat_length, lon_length = input_array.shape
    
    lon_array, lat_array = np.meshgrid(lon_list_out_grid,lat_list_out_grid)

    output_array = np.zeros((lat_list_out_grid.shape[0], lon_list_out_grid.shape[0]))

    output_array = basemap.interp(input_array, lon_list_in_grid, lat_list_in_grid, lon_array, lat_array, order=1)
            
    return output_array

def reorder_lats(latitudes, data_in):
    dlat = np.gradient(latitudes)
    if np.max(dlat) < 0:
        latitudes_out = latitudes[::-1]
        data_out = data_in[::-1,:]
    else:
        latitudes_out = latitudes
        data_out = data_in

    return latitudes_out, data_out

def linear_interpolation_regridding(land_mask_in, topog_in, lat_in_lsm, lat_in_z, lon_out, lat_out, lsm_in, z_in):

    regrided_lsm = linear_interpolate_for_regrid(land_mask_in.longitude.values, lat_in_lsm, lon_out, lat_out, lsm_in)

    regrided_z   = linear_interpolate_for_regrid(topog_in.longitude.values,     lat_in_z,   lon_out, lat_out, z_in)

    regrided_lsm_out = np.round(regrided_lsm)
    
    return regrided_lsm_out, regrided_z

def nearest_neighbour_regrid(land_mask_in, topog_in, lat_in_lsm, lat_in_z, lon_out, lat_out, lsm_in, z_in):

    lons_grid, lats_grid = np.meshgrid(land_mask_in.longitude.values, lat_in_lsm)

    lons_grid_out, lats_grid_out = np.meshgrid(lon_out, lat_out)

    lsm_out = np.zeros_like(lons_grid_out)
    z_out   = np.zeros_like(lons_grid_out)

    for j in range(len(lat_out)):
        print(j)
        for i in range(len(lon_out)):

            distance = (lat_out[j]-lats_grid)**2. + (lon_out[i]-lons_grid)**2.    
            ste = np.where(distance == np.min(distance))
#             if len(ste[0])>1:
#                 ste = ste[0]
#                 pdb.set_trace()
            pdb.set_trace()
            lsm_out[j,i] = lsm_in[ste]
            z_out[j,i] = z_in[ste]            
            
    pdb.set_trace()

def find_equiv(x,y):
    #for two coordinate axes x and y find the index in y representing the closest value to x
    x1=np.zeros(len(x))
    for i in range(0,len(x)):
        x1[i] = np.argmin(abs(y-x[i])) 
    x1=x1.astype(int)
    return x1

def nearest_neighbour_regrid_ruth(land_mask_in, topog_in, lat_in_lsm, lat_in_z, lon_out, lat_out, lsm_in, z_in):

    lons_equiv = find_equiv(lon_out,land_mask_in.longitude.values)
    lats_equiv = find_equiv(lat_out,lat_in_lsm)

    lons_grid_out, lats_grid_out = np.meshgrid(lon_out, lat_out)

    lsm_out = np.zeros_like(lons_grid_out)
    z_out   = np.zeros_like(lons_grid_out)

    for j in range(len(lats_equiv)):
        for i in range(len(lons_equiv)):
            lsm_out[j,i] = lsm_in[lats_equiv[j],lons_equiv[i]]
            z_out[j,i]   = z_in[lats_equiv[j],lons_equiv[i]]
    
    idx = (lsm_out == 0.) & (z_out != 0.)
    z_out[idx] = 0.     
    
    return lsm_out, z_out

def main():
    target_grid=170
    do_linear_interp_regrid = False
    do_nearest_neighbour_regrid = False
    do_nearest_neighbour_regrid_ruth = True

    n_lons = {21:64, 42:128, 85:256, 170:512, 213: 1024}
    n_lats = {21:32, 42:64,  85:128, 170:256, 213: 320}

    land_mask_in = xar.open_dataset('./era_interim_lsm_0_125_deg_grid.nc')
    topog_in     = xar.open_dataset('./era_interim_surf_geopotential_0_125_deg_grid.nc')

    if target_grid=='input':
        input_grid_file = xar.open_dataset('/scratch/sit204/collab/laura_wilcox/files_from_laura_13_12_17/HadGEM_clim/psl_clim.nc')
        
        lon_out = input_grid_file.longitude.values
        lat_out = input_grid_file.latitude.values

        latbs_out = np.zeros(lat_out.shape[0]+1)
        delta_lat = 180./lat_out.shape[0]
        delta_lon = 360./lon_out.shape[0]

        for i in range(len(lat_out)):
            latbs_out[i] = lat_out[i]-0.5*delta_lat
            
    else:
        lat_out, latbs_out = gg.gaussian_latitudes(int(n_lats[target_grid]*0.5))

        delta_lon = 360./n_lons[target_grid]

        lon_out = np.arange(0.,360., delta_lon)
        
    lon_b_out = np.zeros(lon_out.shape[0]+1)

    for i in range(len(lon_out)):
        lon_b_out[i] = lon_out[i]-0.5*delta_lon
    
    lon_b_out[-1] = lon_b_out[-2]+delta_lon

    lat_in_lsm, lsm_in = reorder_lats(land_mask_in.latitude.values, land_mask_in['lsm'].values.squeeze())
    lat_in_z, z_in = reorder_lats(topog_in.latitude.values, topog_in['z'].values.squeeze())

    z_in = z_in /9.8 #rescale surface geopotential to get into a height

    if do_linear_interp_regrid:
        regrided_lsm_out, regrided_z = linear_interpolation_regridding(land_mask_in, topog_in, lat_in_lsm, lat_in_z, lon_out, lat_out, lsm_in, z_in)
        file_name_snippet = '_lininterp_'
    elif do_nearest_neighbour_regrid:
        regrided_lsm_out, regrided_z = nearest_neighbour_regrid(land_mask_in, topog_in, lat_in_lsm, lat_in_z, lon_out, lat_out, lsm_in, z_in)
        file_name_snippet = '_nns_'
    elif do_nearest_neighbour_regrid_ruth:
        regrided_lsm_out, regrided_z = nearest_neighbour_regrid_ruth(land_mask_in, topog_in, lat_in_lsm, lat_in_z, lon_out, lat_out, lsm_in, z_in)        
        file_name_snippet = '_nnr_'     

    if target_grid=='input':
        topo_filename = './era_land'+file_name_snippet+'_input.nc'
    else:
        topo_filename = './era_land'+file_name_snippet+'t'+str(target_grid)+'.nc'
            
    topo_file = Dataset(topo_filename, 'w', format='NETCDF3_CLASSIC')
    lat = topo_file.createDimension('lat', lat_out.shape[0])
    lon = topo_file.createDimension('lon', lon_out.shape[0])
    latitudes = topo_file.createVariable('lat','f4',('lat',))
    longitudes = topo_file.createVariable('lon','f4',('lon',))
    topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
    land_array_netcdf = topo_file.createVariable('land_mask','f4',('lat','lon',))
    latitudes[:] = lat_out
    longitudes[:] = lon_out
    topo_array_netcdf[:] = regrided_z
    land_array_netcdf[:] = regrided_lsm_out
    topo_file.close()
    
if __name__=="__main__":
    main()
