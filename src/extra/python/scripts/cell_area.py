import numpy as np
from netCDF4 import Dataset

def cell_area_all(t_res,base_dir, radius=6376.0e3):
    """read in grid from approriate file, and return 2D array of grid cell areas in metres**2."""
    resolution_file = Dataset(base_dir+'src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')

    lons = resolution_file.variables['lon'][:]
    lats = resolution_file.variables['lat'][:]

    lonb = resolution_file.variables['lonb'][:]
    latb = resolution_file.variables['latb'][:]

    area_array,xsize_array,ysize_array = cell_area_calculate(lons, lats, lonb, latb, radius)


    return area_array,xsize_array,ysize_array

def cell_area(t_res,base_dir):
    """wrapper for cell_area_all, such that cell_area only returns area array, and not xsize_array and y_size_array too."""
    area_array,xsize_array,ysize_array = cell_area_all(t_res,base_dir)
    return area_array

def cell_area_calculate(lons, lats, lonb, latb, radius):

    nlon=lons.shape[0]
    nlat=lats.shape[0]

    area_array = np.zeros((nlat,nlon))
    area_array_2 = np.zeros((nlat,nlon))
    xsize_array = np.zeros((nlat,nlon))
    ysize_array = np.zeros((nlat,nlon))

    for i in np.arange(len(lons)):
        for j in np.arange(len(lats)):
            xsize_array[j,i] = radius*np.absolute(np.radians(lonb[i+1]-lonb[i])*np.cos(np.radians(lats[j])))
            ysize_array[j,i] = radius*np.absolute(np.radians(latb[j+1]-latb[j]))
            area_array[j,i] = xsize_array[j,i]*ysize_array[j,i]
            area_array_2[j,i] = (radius**2.)*np.absolute(np.radians(lonb[i+1]-lonb[i]))*np.absolute(np.sin(np.radians(latb[j+1]))-np.sin(np.radians(latb[j])))

    return area_array_2,xsize_array,ysize_array

def cell_area_from_xar(dataset, lat_name='lat', lon_name = 'lon', latb_name='latb', lonb_name='lonb', radius=6376.0e3):

    lats = dataset[lat_name].values
    lons = dataset[lon_name].values

    try:
        latb = dataset[latb_name].values
        lonb = dataset[lonb_name].values
    except KeyError:
        delta_lat=(lats[1]-lats[0])
        if np.all((lats[1:10]-lats[0:9]) == delta_lat):
            latb = np.zeros((len(lats)+1))

            for latb_idx in range(len(lats)):
                latb[latb_idx] = lats[latb_idx]-delta_lat / 2.
            latb[-1] = lats[-1] + delta_lat / 2.

        delta_lon=(lons[1]-lons[0])
        if np.all((lons[1:10]-lons[0:9]) == delta_lon):
            lonb = np.zeros((len(lons)+1))

            for lonb_idx in range(len(lons)):
                lonb[lonb_idx] = lons[lonb_idx]-delta_lon / 2.
            lonb[-1] = lons[-1] + delta_lon / 2.  
            
            

    area_array,xsize_array,ysize_array = cell_area_calculate(lons, lats, lonb, latb, radius)

    return area_array,xsize_array,ysize_array





if __name__ == "__main__":

    # specify resolution
    t_res = 42
    # specify base dir
    base_dir= '/scratch/sit204/FMS2013/GFDLmoistModel/'
    #return area_array
    area_array=cell_area(t_res,base_dir)

