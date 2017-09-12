import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# specify resolution
t_res = 42

lat_0 = 16.

simple_or_complex_continents = 'complex'

if simple_or_complex_continents == 'simple':
    exp_name='/scratch/sit204/FMS2013/GFDLmoistModel/exp/simple_continents_post_princeton_qflux_control/'
    base_qflux_file_name='simp_p_prin_qflux'
    base_sea_ice_file_name = None
    file_prefix='simp_pp_lon_'
elif simple_or_complex_continents == 'complex':
    exp_name='/scratch/sit204/FMS2013/GFDLmoistModel/exp/annual_mean_ice_princeton_qflux_control/'
    base_qflux_file_name='ami_qflux_ctrl_ice_4320'
    base_sea_ice_file_name = '/scratch/sit204/Data_2013/annual_mean_ice_princeton_qflux_control_1/run001/atmos_monthly_test.nc'
    file_prefix='ami_lon_'
elif simple_or_complex_continents == 'aquaplanet':
    exp_name='/scratch/sit204/FMS2013/GFDLmoistModel/exp/aquaplanet_qflux_control/'
    base_qflux_file_name='aquaplanet_qflux_zm'
    base_sea_ice_file_name = None
    file_prefix='aq_pl_lon_'
else:
    raise NotImplementedError('option not known')

seasonal_sst_anom='always' #Either 'always', 'djf', 'mam', 'jja', or 'son'

do_shielded_anomaly=False

month_dict={'jan':0, 'feb':1, 'mar':2, 'apr':3, 'may':4, 'jun':5, 'jul':6, 'aug':7, 'sep':8, 'oct':9, 'nov':10, 'dec':11}

#read in grid from approriate file
resolution_file = Dataset(exp_name+'input/'+base_qflux_file_name+'.nc', 'r', format='NETCDF3_CLASSIC')

lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]

lonbs = resolution_file.variables['lonb'][:]
latbs = resolution_file.variables['latb'][:]

seasonal_qflux = resolution_file.variables[base_qflux_file_name][:]

time_in = resolution_file.variables['time'][:]

nlon=lons.shape[0]
nlat=lats.shape[0]

nlonb=lonbs.shape[0]
nlatb=latbs.shape[0]

ntime_in=time_in.shape[0]


area_array = np.zeros((nlat,nlon))

for i in np.arange(len(lons)):
    for j in np.arange(len(lats)):
        area_array[j,i] = np.absolute(np.radians(lonbs[i+1]-lonbs[i])*np.cos(np.radians(lats[j]))*np.radians(latbs[j+1]-latbs[j]))




#warmpool_lat_centre = 0.
#warmpool_lon_centre = 240.

warmpool_loc_list=[
# [30.,165.], [30.,180.], [30.,195.],
# [50.,165.], [50.,180.], [50.,195.],
# [30.,300.], [30.,315.], [30.,330.],
# [50.,315.], [50.,330.], [50.,345.],
# [10.,165.], [10.,180.], [10.,195.], [10.,210.],
# [10.,315.], [10.,330.]
# [0.,165.], [0.,180.], [0.,195.], [0.,210.],
# [0.,330.], [0.,345.],
#[0.,225.], [0.,240.], [0.,255.], [0.,270.]
# [40., 165.], [40., 180.], [40., 195.]

[0., 330.], [10., 330.], [30., 330.], [40.,330.], [50., 330.]


]

#warmpool_loc_list=[[0.,210.]]

#warmpool_loc_list=[ [30.,330.]]
#warmpool_loc_list=[ [45.,345.]]

for file_values in warmpool_loc_list:

    warmpool_array = np.zeros((nlat,nlon))
    warmpool_shield_array = np.zeros((nlat,nlon))

    warmpool_lat_centre=file_values[0]
    warmpool_lon_centre=file_values[1]

    print((warmpool_lat_centre,warmpool_lon_centre))

    warmpool_width = 7.5
    warmpool_width_lon = 7.5
    warmpool_amp = 200.


    for j in np.arange(len(lats)):
         lat = 0.5*(latbs[j+1] + latbs[j])
         lat = (lat-warmpool_lat_centre)/warmpool_width
    #     if( np.absolute(lat) <= 1.0 ):
         for i in np.arange(len(lons)):
              lon = 0.5*(lonbs[i+1] + lonbs[i])
              lon = (lon-warmpool_lon_centre)/warmpool_width_lon
              if( lat**2.+lon**2. <= 1.0 ):
                  warmpool_array[j,i] = (1.-lat**2.-lon**2.)*warmpool_amp
              if( do_shielded_anomaly and ((lat**2.+lon**2. > 1.0) and ( lat**2.+lon**2. <= 2.0 )) ):
                  warmpool_shield_array[j,i] = (1.-(np.sqrt(lat**2.+lon**2.)-1.5)**2)

    warmpool_integral = np.sum(area_array*warmpool_array)
    warmpool_shield_integral = np.sum(area_array*warmpool_shield_array)

    land_array  = np.zeros((nlat,nlon))
    lon_array, lat_array = np.meshgrid(lons,lats)

    try:
        land_file = Dataset(exp_name+'/input/land.nc', 'r', format='NETCDF3_CLASSIC')
        land_array = land_file.variables['land_mask'][:]
    except:
        print('No land file')


    if base_sea_ice_file_name is not None:
        sea_ice_file = Dataset(base_sea_ice_file_name, 'r', format='NETCDF3_CLASSIC')

        ice_array = sea_ice_file.variables['albedo'][:]

        ice_array_one_month=ice_array[0,:,:].squeeze()
        ice_array_one_month=np.round(ice_array_one_month, decimals=2)

        ice_idx = (ice_array_one_month == 0.7)        

        ice_mask=np.zeros_like(ice_array_one_month)

        ice_mask[ice_idx]=1.0

    else:
        ice_mask = np.zeros_like(land_array)



    ocean_array = np.ones((nlat,nlon))-land_array

    if not do_shielded_anomaly:

        cooling_idx = (ocean_array == 1.0) & (warmpool_array == 0.) & (ice_mask == 0.0)

        cooling_array = np.zeros((nlat,nlon))
        cooling_array[cooling_idx] = 1.

        warmpool_array[cooling_idx]=-1.*warmpool_integral/np.sum(area_array[cooling_idx])

    else:
        shield_scaling = warmpool_integral/(-1.*warmpool_shield_integral)
        warmpool_array=warmpool_array+(warmpool_shield_array*shield_scaling)



    warmpool_integral_adj = np.sum(area_array*warmpool_array)

    print(('warmpool_integral', warmpool_integral))
    print(('warmpool_integral_adj', warmpool_integral_adj))



    warmpool_array_total=np.copy(seasonal_qflux)
#    warmpool_array_total=np.zeros_like(seasonal_qflux)

    for t in np.arange(ntime_in):
        if seasonal_sst_anom=='always':
            warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition=''
        elif seasonal_sst_anom=='djf':
            if t==0 or t==1 or t==11:
                warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition='_djf'
        elif seasonal_sst_anom=='mam':
            if t==2 or t==3 or t==4:
                warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition='_mam'
        elif seasonal_sst_anom=='jja':
            if t==5 or t==6 or t==7:
                warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition='_jja'
        elif seasonal_sst_anom=='son':
            if t==8 or t==9 or t==10:
                warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition='_son'
        elif seasonal_sst_anom=='djf_mid':
            if t==0:
                warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition='_djf_mid'
        elif seasonal_sst_anom=='mam_mid':
            if t==3:
                warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition='_mam_mid'
        elif seasonal_sst_anom=='jja_mid':
            if t==6:
                warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition='_jja_mid'
        elif seasonal_sst_anom=='son_mid':
            if t==9:
                warmpool_array_total[t,:,:]=warmpool_array_total[t,:,:]+warmpool_array
            file_name_seasonal_addition='_son_mid'


    plt.figure()

    lon_0 = lons.mean()
    lat_0 = lats.mean()

    lon_array, lat_array = np.meshgrid(lons,lats)
    lonb_array, latb_array = np.meshgrid(lonbs,latbs)

    m = Basemap(lat_0=lat_0,lon_0=lon_0)
    xi, yi = m(lon_array, lat_array)

    cs = m.contourf(xi,yi, warmpool_array_total[0,:,:], cmap=plt.get_cmap('RdBu_r'))
    plt.xticks(np.linspace(0,360,13))
    plt.yticks(np.linspace(-90,90,7))
    cb = plt.colorbar(cs, shrink=0.5, extend='both')

    description_string=file_prefix+str(int(warmpool_lon_centre))+'_lat_'+str(int(warmpool_lat_centre))+'_a_'+str(int(warmpool_amp))
    
    var_name_out=description_string+file_name_seasonal_addition

    warmpool_file = Dataset(var_name_out+'.nc', 'w', format='NETCDF3_CLASSIC')
    lat = warmpool_file.createDimension('lat', nlat)
    lon = warmpool_file.createDimension('lon', nlon)

    latb = warmpool_file.createDimension('latb', nlatb)
    lonb = warmpool_file.createDimension('lonb', nlonb)

    time = warmpool_file.createDimension('time', 0)

    latitudes = warmpool_file.createVariable('lat','f4',('lat',))
    longitudes = warmpool_file.createVariable('lon','f4',('lon',))

    latitudebs = warmpool_file.createVariable('latb','f4',('latb',))
    longitudebs = warmpool_file.createVariable('lonb','f4',('lonb',))

    times = warmpool_file.createVariable('time','d',('time',))

    times.units = 'days since 0000-01-01 00:00:00.0'
    times.calendar = 'THIRTY_DAY_MONTHS'
    times.calendar_type = 'THIRTY_DAY_MONTHS'
    times.cartesian_axis = 'T'

    latitudes.units = 'degrees_N'.encode('utf-8')
    latitudes.cartesian_axis = 'Y'
    latitudes.edges = 'latb'
    latitudes.long_name = 'latitude'

    longitudes.units = 'degrees_E'.encode('utf-8')
    longitudes.cartesian_axis = 'X'
    longitudes.edges = 'lonb'
    longitudes.long_name = 'longitude'

    latitudebs.units = 'degrees_N'.encode('utf-8')
    latitudebs.cartesian_axis = 'Y'
    latitudebs.long_name = 'latitude edges'

    longitudebs.units = 'degrees_E'.encode('utf-8')
    longitudebs.cartesian_axis = 'X'
    longitudebs.long_name = 'longitude edges'


    warmpool_array_netcdf = warmpool_file.createVariable(var_name_out,'f4',('time','lat','lon',))

    ice_mask_netcdf = warmpool_file.createVariable('ice_mask','f4',('lat','lon',))


    latitudes[:] = lats
    longitudes[:] = lons

    latitudebs[:] = latbs
    longitudebs[:] = lonbs

    times[:] = time_in
    warmpool_array_netcdf[:] = warmpool_array_total
    ice_mask_netcdf[:] = ice_mask

    warmpool_file.close()



