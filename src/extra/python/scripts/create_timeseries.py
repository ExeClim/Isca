# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num
import sys
import pdb
import os

__author__='Stephen Thomson'

def create_grid(manual_grid_option):

    if(manual_grid_option):

        lons =  [0.]
        lonbs = [0., 360.]

        lats =  [0.]
        latbs = [-90., 90.]

        nlon=len(lons)
        nlat=len(lats)

        nlonb=len(lonbs)
        nlatb=len(latbs)

    else:
        t_res=42

        try:
            GFDL_BASE        = os.environ['GFDL_BASE']
        except Exception as e:
            print('Environment variable GFDL_BASE must be set')
            exit(0)

        resolution_file = Dataset(GFDL_BASE+'/src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')

        lons = resolution_file.variables['lon'][:]
        lats = resolution_file.variables['lat'][:]

        lonbs = resolution_file.variables['lonb'][:]
        latbs = resolution_file.variables['latb'][:]

        nlon=lons.shape[0]
        nlat=lats.shape[0]

        nlonb=lonbs.shape[0]
        nlatb=latbs.shape[0]

    return lons,lats,lonbs,latbs,nlon,nlat,nlonb,nlatb


def create_pressures():

    p_full = [300.,900.]
    p_half=[0.,600.,1200.]

#    p_full = [0.017223, 0.078654, 0.180749, 0.375337, 0.71314, 1.253787, 2.060447, 3.193908, 4.707791, 6.646065, 9.04322, 11.926729, 15.321008, 19.251999, 23.751706, 28.862245, 34.639281, 41.154864, 48.49985, 56.786085, 66.14858, 76.747823, 88.772421, 102.442154, 118.011585, 135.77429, 156.067832, 179.279569, 205.853424, 236.297753, 271.194477, 311.209677, 357.105858, 409.756139, 470.160652, 539.465462, 618.984358, 710.223896, 814.912138, 935.031579, 1100.]

#    p_half = [0., 0.046817, 0.115544, 0.254986, 0.510238, 0.937468, 1.599324, 2.558913, 3.874132, 5.593891, 7.756986, 10.393581, 13.528609, 17.186225, 21.394476, 26.18965, 31.62002, 37.748946, 44.657459, 52.446517, 61.23914, 71.182614, 82.450941, 95.247663, 109.809174, 126.408616, 145.360464, 167.025894, 191.81904, 220.214288, 252.754733, 290.062007, 332.84765, 381.926291, 438.230874, 502.830259, 576.949493, 661.993148, 759.57212, 871.53435, 1000., 1200.]
    
    if(np.min(p_half)!=0.):
        print('Must have minimum p_half level = 0., as otherwise model data will be missing near the top levels.')
        sys.exit(0)

    if(np.max(p_half)<=1000.):
        print('Must have maximum p_half level > 1000., as otherwise vertical interpolation will fail near the surface when p_surf>1000.')
        sys.exit(0)

    npfull=len(p_full)
    nphalf=len(p_half)

    return p_full,p_half,npfull,nphalf


def create_time_arr(num_years,is_climatology,time_spacing):

    if(is_climatology):
        if(num_years!=1.):
            print('note that for climatology file only one year is required, so setting num_years=1.')
        num_days=360.
        num_years=1.
#        time_spacing=num_days//10
        day_number = np.linspace(0,num_days,time_spacing+1)[1:]-(num_days/(2.*time_spacing))
        


        time_units='days since 0000-01-01 00:00:00.0'
        print('when creating a climatology file, the year of the time units must be zero. This is how the model knows it is a climatology.')
    else:
        num_days=num_years*360.
#        time_spacing=num_years
        day_number = np.linspace(0,num_days,time_spacing+1)
        time_units='days since 0001-01-01 00:00:00.0'

    half_spacing = (day_number[1]-day_number[0])/2.
    lower_time_bounds = day_number - half_spacing
    upper_time_bounds = day_number + half_spacing

    time_bounds=np.zeros((len(lower_time_bounds),2))
    time_bounds[:,0]=lower_time_bounds
    time_bounds[:,1]=upper_time_bounds

    time_arr = day_number_to_date(day_number)
    ntime=len(time_arr)

    return time_arr,day_number,ntime,time_units,time_bounds


def output_multiple_variables_to_file(data_dict,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,number_dict, time_bounds=None):
    """Default is to accept multiple data variables to put in the file.
    Input format in data_dict = {variable_name:variable_array}.
    Currently only accepts all 2D fields or all 3D fields. 
    Could be updated in future to accept a mix."""

    output_file = Dataset(file_name, 'w', format='NETCDF3_CLASSIC')

    if p_full is None:
        is_thd=False
    else:
        is_thd=True

    if time_arr is None:
        add_time = False
    else:
        add_time = True

    lat = output_file.createDimension('lat', number_dict['nlat'])
    lon = output_file.createDimension('lon', number_dict['nlon'])

    use_latbs=False
    use_lonbs=False
    
    if latbs is not None:
        latb = output_file.createDimension('latb', number_dict['nlatb'])
        latitudebs = output_file.createVariable('latb','d',('latb',))
        use_latbs=True

    if lonbs is not None:    
        lonb = output_file.createDimension('lonb', number_dict['nlonb'])
        longitudebs = output_file.createVariable('lonb','d',('lonb',))
        use_lonbs=True

    if is_thd:
        pfull = output_file.createDimension('pfull', number_dict['npfull'])
        phalf = output_file.createDimension('phalf', number_dict['nphalf'])

    if add_time:
        time = output_file.createDimension('time', 0) #s Key point is to have the length of the time axis 0, or 'unlimited'. This seems necessary to get the code to run properly. 

    latitudes = output_file.createVariable('lat','d',('lat',))
    longitudes = output_file.createVariable('lon','d',('lon',))
    
    if is_thd:
        pfulls = output_file.createVariable('pfull','d',('pfull',))
        phalfs = output_file.createVariable('phalf','d',('phalf',))

    if add_time:
        times = output_file.createVariable('time','d',('time',))

    latitudes.units = 'degrees_N'.encode('utf-8')
    latitudes.cartesian_axis = 'Y'
    latitudes.long_name = 'latitude'

    longitudes.units = 'degrees_E'.encode('utf-8')
    longitudes.cartesian_axis = 'X'
    longitudes.long_name = 'longitude'

    if use_latbs:
        latitudes.edges = 'latb'
        latitudebs.units = 'degrees_N'.encode('utf-8')
        latitudebs.cartesian_axis = 'Y'
        latitudebs.long_name = 'latitude edges'

    if use_lonbs:
        longitudes.edges = 'lonb'
        longitudebs.units = 'degrees_E'.encode('utf-8')
        longitudebs.cartesian_axis = 'X'
        longitudebs.long_name = 'longitude edges'

    if is_thd:
        pfulls.units = 'hPa'
        pfulls.cartesian_axis = 'Z'
        pfulls.positive = 'down'
        pfulls.long_name = 'full pressure level'

        phalfs.units = 'hPa'
        phalfs.cartesian_axis = 'Z'
        phalfs.positive = 'down'
        phalfs.long_name = 'half pressure level'

    if add_time:
        times.units = time_units
        times.calendar = 'THIRTY_DAY_MONTHS'
        times.calendar_type = 'THIRTY_DAY_MONTHS'
        times.cartesian_axis = 'T'

    if time_bounds is not None:

        vertex_dimension = output_file.createDimension('nv', 2) #s 
        vertices = output_file.createVariable('nv','d',('nv',))
        vertices[:] = [1.,2.]
        time_bounds_file = output_file.createVariable('time_bounds','d',('time','nv'))

        time_bounds_file.long_name = 'time axis boundaries'
        time_bounds_file.units     = time_units
        time_bounds_file[:] = time_bounds

        times.bounds = 'time_bounds'

    latitudes[:] = lats
    longitudes[:] = lons

    if use_latbs:
        latitudebs[:] = latbs
    if use_lonbs:
        longitudebs[:] = lonbs

    if is_thd:
        pfulls[:]     = p_full
        phalfs[:]     = p_half

    if add_time:
        if type(time_arr[0])!=np.float64 and type(time_arr[0])!=np.int64 :
            times[:]     = date2num(time_arr,units='days since 0001-01-01 00:00:00.0',calendar='360_day')
        else:
            times[:]     = time_arr

    for variable_name in data_dict.keys():
        if is_thd:
            if add_time:
                three_d_dims = ('time','pfull', 'lat','lon',)
            else:
                three_d_dims = ('pfull', 'lat','lon',)            

            output_array_netcdf = output_file.createVariable(variable_name,'f4',three_d_dims)
        else:
            if add_time:
                two_d_dims = ('time','lat','lon',)
            else:
                two_d_dims = ('lat','lon',)              
            output_array_netcdf = output_file.createVariable(variable_name,'f4',two_d_dims)

        output_array_netcdf[:] = data_dict[variable_name]

    output_file.close()

def output_to_file(data,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,variable_name,number_dict, time_bounds=None):

    """Special interface for script wanting to output 1 variable only."""

    data_dict_to_send = {variable_name:data}
    output_multiple_variables_to_file(data_dict_to_send,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,number_dict,time_bounds=time_bounds)
