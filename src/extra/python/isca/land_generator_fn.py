# Function to allow land to be generated for a range of scenarios

# Land Options:
# 'square' (default) Square block of land with boundaries specified by boundaries keyword, a list of 4 floats in the form [S,N,W,E]
# 'continents_old' Choose continents from the original continent set-up adapted from the Sauliere 2012 paper (Jan 16), including North and South America, Eurasia, and Africa. 
# 'continents' Choose continents from a newer continet set-up allowing addition of India, Australia, and South East Asia.
#    If continents keyword is set to 'all' (default), then this will include all possibilities for the given set-up
#    Alternatively, continents can be set to a list of strings containing options: 
#    NA - North America
#    SA - South America
#    EA - Eurasia
#    AF - Africa
#    OZ - Australia
#    IN - India
#    SEA - South East Asia

# Topography Options:
#'none' (default) Topography set to zero everywhere
#'sauliere2012' Choose mountains from Sauliere 2012 configuration using mountains keyword. Default is 'all', alternatively only 'rockys' or 'tibet' may be specified
#'gaussian' Use parameters specified in topo_gauss keyword to set up a Gaussian mountain. topo_gauss should be a list in the form: [central_lat,central_lon,radius_degrees,std_dev,height]

# Topography boundary options:
# If waterworld keyword is set to False (default), then topography can only be non-zero on continents - important as topography has a Gaussian structure and tends exponentially to zero.
# If waterworld keyword is set to True, aquamountains are possible - extra work needed here to deal with exponential issues!

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os

def write_land(exp,land_mode='square',boundaries=[20.,60.,20.,60.],continents=['all'],topo_mode='none',mountains=['all'],topo_gauss=[40.,40.,20.,10.,3500.],waterworld=False):

# Common features of set-ups
    # specify resolution
    t_res = 42
    #read in grid from approriate file
    GFDL_BASE = os.environ['GFDL_BASE']
    resolution_file = Dataset(GFDL_BASE + 'src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')
    lons = resolution_file.variables['lon'][:]
    lats = resolution_file.variables['lat'][:]
    lonb = resolution_file.variables['lonb'][:]
    latb = resolution_file.variables['latb'][:]
    nlon=lons.shape[0]
    nlat=lats.shape[0]
    topo_array = np.zeros((nlat,nlon))
    land_array = np.zeros((nlat,nlon))

    #make 2d arrays of latitude and longitude
    lon_array, lat_array = np.meshgrid(lons,lats)
    lonb_array, latb_array = np.meshgrid(lonb,latb)
    
    #create dictionary for continents
    cont_dic = {'NA':0, 'SA':1, 'EA':2, 'AF':3, 'OZ':4, 'IN':5, 'SEA':6}

# Firstly determine the land set-up to be used    
    # 1) Set-up in which a square of land is included
    if land_mode=='square':
        idx = (boundaries[0] <= lat_array) & (lat_array < boundaries[1]) & (boundaries[2] < lon_array) & (boundaries[3] > lon_array)
        land_array[idx] = 1.0
    
    # 2) Set-up in which some or all of 'original' continents are included
    elif land_mode=='continents_old':  #Older configuration of continents: Addition of India and SE Asia required some restructuring. This may be removed once obsolete.
        idx_c = np.zeros((4,nlat,nlon), dtype=bool)
        idx_c[0,:,:] = (103.-43./40.*(lon_array-180) < lat_array) & ((lon_array-180)*43./50. -51.8 < lat_array) &( lat_array < 60.)   #North America
        idx_c[1,:,:] = (737.-7.2*(lon_array-180) < lat_array) & ((lon_array-180)*10./7. + -212.1 < lat_array) &( lat_array < -22./45*(lon_array-180) +65.9)   #South America
        eurasia_pos = (17. <= lat_array) & (lat_array < 60.) & (-5. < lon_array) & ( 43./40.*lon_array -101.25 < lat_array)
        eurasia_neg =  (17. <= lat_array) & (lat_array < 60.) & (355. < lon_array)
        idx_c[2,:,:] = eurasia_pos + eurasia_neg    #Eurasia
        africa_pos = (lat_array < 17.) & (-52./27.*lon_array + 7.37 < lat_array) & (52./38.*lon_array -65.1 < lat_array) 
        africa_neg = (lat_array < 17.) & (-52./27.*(lon_array-360) + 7.37 < lat_array)
        idx_c[3,:,:] = africa_pos + africa_neg   #Africa
        
        if 'all' in continents:
            idx = idx_c[0,:,:] + idx_c[1,:,:] + idx_c[2,:,:] + idx_c[3,:,:] 
            land_array[idx] = 1.            
        else:
            idx = np.zeros((nlat,nlon), dtype=bool)
            for cont in continents:
                idx = idx + idx_c[cont_dic[cont],:,:]
                land_array[idx] = 1.0
                
    # 2) Set-up in which some or all of 'new' continents are included
    elif land_mode=='continents':  
        idx_c = np.zeros((7,nlat,nlon), dtype=bool)
        idx_c[0,:,:] = (103.-43./40.*(lon_array-180) < lat_array) & ((lon_array-180)*43./50. -51.8 < lat_array) &( lat_array < 60.)   #North America
        idx_c[1,:,:] = (737.-7.2*(lon_array-180) < lat_array) & ((lon_array-180)*10./7. + -212.1 < lat_array) &( lat_array < -22./45*(lon_array-180) +65.9)   #South America
        eurasia_pos = (23. <= lat_array) & (lat_array < 60.) & (-8. < lon_array) & ( 43./40.*lon_array -101.25 < lat_array) 
        eurasia_neg =  (23. <= lat_array) & (lat_array < 60.) & (352. < lon_array)
        idx_c[2,:,:] = eurasia_pos + eurasia_neg    #Eurasia
        africa_pos = (lat_array < 23.) & (-52./27.*lon_array + 7.59 < lat_array) & (52./38.*lon_array -65.1 < lat_array)
        africa_neg = (lat_array < 23.) & (-52./27.*(lon_array-360) + 7.59 < lat_array)
        idx_c[3,:,:] = africa_pos + africa_neg   #Africa
        idx_c[4,:,:] = (lat_array > - 35.) & (lat_array < -17.) & (lon_array > 115.) & (lon_array < 150.) #Australia
        idx_c[5,:,:] = (lat_array < 23.) & (-15./8.*lon_array + 152 < lat_array) & (15./13.*lon_array - 81 < lat_array) #India
        idx_c[6,:,:] = (lat_array < 23.) & ( 43./40.*lon_array -101.25 < lat_array) & (-14./13.*lon_array +120 < lat_array)      #South East Asia
        if 'all' in continents:
            idx = idx_c[0,:,:] + idx_c[1,:,:] + idx_c[2,:,:] + idx_c[3,:,:] + idx_c[4,:,:] + idx_c[5,:,:] + idx_c[6,:,:] 
            land_array[idx] = 1.
        else:
            idx = np.zeros((nlat,nlon), dtype=bool)
            for cont in continents:
                idx = idx + idx_c[cont_dic[cont],:,:]
                land_array[idx] = 1.
                
    elif land_mode=='none':  
        land_array = np.zeros((nlat,nlon))
        
    
# Next produce a topography array
    if topo_mode == 'none':
        topo_array = np.zeros((nlat,nlon))
        
    elif topo_mode == 'sauliere2012':
        # Rockys from Sauliere 2012
        h_0 = 2670.
        central_lon = 247.5
        central_lat = 40.
        L_1 = 7.5
        L_2 = 20.
        gamma_1 = 42.
        gamma_2 = 42.
        delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
        delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2
        h_arr_rockys = h_0 * np.exp(-(delta_1**2. + delta_2**2.))
        idx_rockys = (h_arr_rockys / h_0 > 0.05) #s make sure exponentials are cut at some point - use the value from p70 of Brayshaw's thesis. 

        #Tibet from Sauliere 2012
        h_0 = 5700.
        central_lon = 82.5
        central_lat = 28
        L_1 = 12.5
        L_2 = 12.5
        gamma_1 = -49.5
        gamma_2 = -18.
        delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
        delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2
        h_arr_tibet_no_amp = np.exp(-(delta_1**2.))*(1./delta_2)*np.exp(-0.5*(np.log(delta_2))**2.)
        maxval = np.nanmax(h_arr_tibet_no_amp) #For some reason my maximum value of h_arr_tibet_no_amp > 1. Renormalise so h_0 sets amplitude. 
        h_arr_tibet = (h_arr_tibet_no_amp/maxval)*h_0
        idx_tibet = (h_arr_tibet / h_0 > 0.05)

        if 'all' in mountains:
            topo_array[idx_rockys] = h_arr_rockys[idx_rockys]
            topo_array[idx_tibet] =  h_arr_tibet[idx_tibet]
        elif 'rockys' in mountains:
            topo_array[idx_rockys] = h_arr_rockys[idx_rockys]
        elif 'tibet' in mountains:
            topo_array[idx_tibet] =  h_arr_tibet[idx_tibet]
        else:
            print('No valid mountain options detected for Sauliere 2012 topography')

            
    elif topo_mode == 'gaussian':
        #Options to define simple Gaussian Mountain
        central_lat = topo_gauss[0]
        central_lon = topo_gauss[1]
        radius_degrees = topo_gauss[2]
        std_dev = topo_gauss[3]
        height = topo_gauss[4]
        rsqd_array = np.sqrt((lon_array - central_lon)**2.+(lat_array - central_lat)**2.)
        #generalise to ellipse - needs checking but may be useful later (RG)
        #ax_rot = 1. #gradient of new x axis
        #ax_rat = 2. #axis ratio a**2/b**2
        #rsqd_array = np.sqrt((lon_array - central_lon + ax_rot*(lat_array - central_lat))**2.+ ax_rat*(lat_array - central_lat - ax_rot*(lon_array - central_lon))**2.)*np.cos(np.arctan(ax_rot))
        #divide by factor of cos(atan(m)) to account for change in coords
        idx = (rsqd_array < radius_degrees) 
        topo_array[idx] = height* np.exp(-(rsqd_array[idx]**2.)/(2.*std_dev**2.))
        
    else:
        print('Invalid topography option given')
        

    if waterworld != True:      #Leave flexibility to allow aquamountains!
        idx = (land_array == 0.) & (topo_array != 0.)
        topo_array[idx] = 0. 


    #Write land and topography arrays to file
    topo_filename = GFDL_BASE + 'exp/' + exp + '/input/land.nc'
    topo_file = Dataset(topo_filename, 'w', format='NETCDF3_CLASSIC')
    lat = topo_file.createDimension('lat', nlat)
    lon = topo_file.createDimension('lon', nlon)
    latitudes = topo_file.createVariable('lat','f4',('lat',))
    longitudes = topo_file.createVariable('lon','f4',('lon',))
    topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
    land_array_netcdf = topo_file.createVariable('land_mask','f4',('lat','lon',))
    latitudes[:] = lats
    longitudes[:] = lons
    topo_array_netcdf[:] = topo_array
    land_array_netcdf[:] = land_array
    topo_file.close()
    print('Output written to: ' + topo_filename)


    #Show configuration on screen to allow checking
    lon_0 = lons.mean()
    lat_0 = lats.mean()
    m = Basemap(lat_0=lat_0,lon_0=lon_0)
    xi, yi = m(lon_array, lat_array)
    plt.figure()
    if land_mode != 'none':
        m.contour(xi,yi,land_array)
    if topo_mode != 'none':
        cs = m.contourf(xi,yi,topo_array, cmap=plt.get_cmap('RdBu_r'))
        cb = plt.colorbar(cs, shrink=0.5, extend='both')
    plt.xticks(np.linspace(0,360,13))
    plt.yticks(np.linspace(-90,90,7))
    plt.show()



if __name__ == "__main__":

    write_land('test',land_mode='continents')

