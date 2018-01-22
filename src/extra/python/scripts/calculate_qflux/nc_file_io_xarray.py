import xarray as xar
from netCDF4 import Dataset
import numpy as np
import calendar_calc as cal
import os
import sys
import pdb
import create_timeseries as cts
from scipy import stats

__author__='Stephen Thomson'

def read_data( base_dir, exp_name, start_file, end_file, avg_or_daily, use_interpolated_pressure_level_data, model='fms13', file_name=None):

    if model=='fms13':

        possible_format_strs = [[base_dir+'/'+exp_name+'/run%03d' % m for m in range(start_file, end_file+1)],
                                [base_dir+'/'+exp_name+'/run%04d' % m for m in range(start_file, end_file+1)],
                                [base_dir+'/'+exp_name+'/run%d'   % m for m in range(start_file, end_file+1)]]

        if(use_interpolated_pressure_level_data):
            if avg_or_daily == 'monthly':
#                 extra='_interp.nc'
                extra='_interp_new_height.nc'

            else:
                extra='_interp_new_height_temp.nc'
        else:
            extra='.nc'

        thd_string = '/atmos_'+avg_or_daily+extra

        for format_str_files in possible_format_strs:
            files_temp = format_str_files

            thd_files = [s + thd_string for s in files_temp]

            thd_files_exist=[os.path.isfile(s) for s in thd_files]

            if thd_files_exist[0]:
                break

            if not thd_files_exist[0] and possible_format_strs.index(format_str_files)==(len(possible_format_strs)-1):
                raise EOFError('EXITING BECAUSE NO APPROPRIATE FORMAT STR', [thd_files[elem] for elem in [0] if not thd_files_exist[elem]])


        print(thd_files[0])

        if not(all(thd_files_exist)):
            raise EOFError('EXITING BECAUSE OF MISSING FILES', [thd_files[elem] for elem in range(len(thd_files_exist)) if not thd_files_exist[elem]])

        thd_file_size = [os.path.getsize(s) for s in thd_files]

        mode_file_size = stats.mode(thd_file_size).mode[0]

        thd_files_too_small = [thd_file_size[s] < 0.75*mode_file_size for s in range(len(thd_files))]

        if np.any(thd_files_too_small):
            raise EOFError('EXITING BECAUSE OF FILE TOO SMALL', [thd_files[elem] for elem in range(len(thd_files_exist)) if thd_files_too_small[elem]])


        size_list = init(thd_files[0])

        try:
            da_3d = xar.open_mfdataset(    thd_files,
                    decode_times=False,  # no calendar so tell netcdf lib
                    # choose how data will be broken down into manageable
                    # chunks.
        #                        concat_dim='time',
                    chunks={'time': size_list['ntime'],
                            'lon': size_list['nlons']//4,
                            'lat': size_list['nlats']//2})
        except:
            da_3d = xar.open_mfdataset(    thd_files,
                    decode_times=False,  # no calendar so tell netcdf lib
                    # choose how data will be broken down into manageable
                    # chunks.
        #                        concat_dim='time',
                    chunks={'xofyear': size_list['ntime'],
                            'lon': size_list['nlons']//4,
                            'lat': size_list['nlats']//2})        

        names_dict = {'xofyear':'time'}

        for name in names_dict.keys():
            try:                        
                da_3d.rename({name:names_dict[name]}, inplace=True)
            except ValueError:
                pass       

        time_arr = da_3d.time
        date_arr = cal.day_number_to_date(time_arr)

    da_3d.coords['dayofyear_ax'] = (('dayofyear_ax'),np.unique(date_arr.dayofyear))
    da_3d.coords['months_ax'] = (('months_ax'),np.unique(date_arr.month))
    da_3d.coords['seasons_ax'] = (('seasons_ax'),np.arange(4))
    da_3d.coords['years_ax'] = (('years_ax'),date_arr.year)
    da_3d.coords['all_time_ax'] = (('all_time_ax'),np.arange(1))

    da_3d.coords['dayofyear'] = (('time'),date_arr.dayofyear)
    da_3d.coords['months'] = (('time'),date_arr.month)
    da_3d.coords['years'] = (('time'),date_arr.year)

    seasons_arr = cal.month_to_season(date_arr.month, avg_or_daily)

    da_3d.coords['seasons'] = (('time'), seasons_arr)

    two_months_arr = cal.month_to_two_months(date_arr.month, avg_or_daily)
    
    da_3d.coords['two_months'] = (('time'),two_months_arr)
    da_3d.coords['two_months_ax'] = (('two_months_ax'),np.unique(two_months_arr))



    da_3d.coords['all_time'] = (('time'),time_arr/time_arr)

    da_3d.coords['seq_months'] = (('time'),date_arr.month+12.*((date_arr.year-np.min(date_arr.year))))

    da_3d.coords['seq_seasons'] = (('time'),cal.recurring_to_sequential(seasons_arr))
    
    da_3d.coords['seq_all_time'] = (('time'),list(range(len(time_arr))))    

    da_3d.coords['seq_days'] = (('time'),cal.recurring_to_sequential(date_arr.dayofyear))

    da_3d.coords['seq_seasons_ax'] = (('seq_seasons_ax'),np.mod(np.min(da_3d.seq_seasons.values),4)+np.arange(len(np.unique(da_3d.seq_seasons.values))))

    da_3d.attrs['exp_name']=exp_name
    da_3d.attrs['start_file']=start_file
    da_3d.attrs['end_file']=end_file
    da_3d.attrs['data_type']=avg_or_daily
    try:
        da_3d['precipitation']
    except KeyError:
        try:
           print('aggregating rain')
           da_3d['convection_rain']
           da_3d['condensation_rain']
        except KeyError:
           print('no precip output present')
        else:        
            da_3d['precip']=(('time','lat','lon'),da_3d['convection_rain']+da_3d['condensation_rain'])
            print('done aggregating rain')




    thd_data = da_3d

    return thd_data, time_arr, size_list

def init( nc_file_init):
    "Uses the first nc file to read longitudes, lats etc."
    fh_init = Dataset(nc_file_init, mode='r')
    
    # Initialise variables that will be used throughout
    
    try:
        lons = fh_init.variables['lon'][:]
        lats = fh_init.variables['lat'][:]
    except KeyError:
        try:
            lons = fh_init.variables['g0_lon_3'][:]
            lats = fh_init.variables['g0_lat_2'][:]
        except KeyError:
            lons = fh_init.variables['longitude'][:]
            lats = fh_init.variables['latitude'][:]

    try:
        time_init = fh_init.variables['time'][:] #This is specific and needs deleting once ntime has been found.
    except KeyError:
        try:        
            time_init = fh_init.variables['initial_time0_hours'][:] #This is specific and needs deleting once ntime has been found.
        except KeyError:
            time_init = fh_init.variables['xofyear'][:] #This is specific and needs deleting once ntime has been found.        
    try:
        pfull = fh_init.variables['pfull'][:]
    except KeyError:
        try:
            pfull = fh_init.variables['lv_IBSL1'][:]
        except KeyError:
            pfull=None
    
    #Initialise nlons, nlats, etc
    nlons=len(lons)
    nlats=len(lats)
    ntime=len(time_init)
    try:
        nlevs=len(pfull)
    except TypeError:
        nlevs=0
    
    # Delete any unwanted variables
    del time_init
    fh_init.close()
    
    size_list={}
    size_list['nlons']=nlons
    size_list['nlats']=nlats
    size_list['nlevs']=nlevs
    size_list['ntime']=ntime
    
    return    size_list

def read_land( base_dir,exp_name,land_present, use_interpolated_pressure_level_data, size_list,land_file='input/land.nc', lats_in = None):
    "Function for reading in land mask."
    
    #Create thd (3D) and twd (2D) data stores.
    land_array=np.zeros((size_list['nlats'],size_list['nlons']))
    topo_array=np.zeros((size_list['nlats'],size_list['nlons']))
    
    if land_file is not None:
        if(land_present or use_interpolated_pressure_level_data):
            nc_file = base_dir+'exp/'+exp_name+land_file
            print(nc_file)
            try:
                fh = Dataset(nc_file, mode='r')
            except (OSError, RuntimeError):
                nc_file = base_dir+exp_name+land_file
                fh = Dataset(nc_file, mode='r')
            if land_present:
                try:
                    land_array=fh.variables['land_mask'][:]
                except KeyError:
                    land_array=fh.variables['lsm'][:]
                    land_array=land_array.squeeze()   
                    try:         
                        lat_land = fh.variables['latitude'][:]
                    except:
                        lat_land = fh.variables['lat'][:]
                        
                
                    if lats_in[0]!=lat_land[0]:
                        land_array=land_array[::-1,:]
                    
            if use_interpolated_pressure_level_data:
                topo_array=fh.variables['zsurf'][:]
            fh.close()

    return land_array, topo_array

def output_nc_file(dataset, field_name, model_params, output_dict):

    #create grid
    try:
        lons = dataset.lon.values
        lats = dataset.lat.values
        latbs = dataset.latb.values
        lonbs = dataset.lonb.values
        nlon = lons.shape[0]
        nlat = lats.shape[0]
        nlatb = latbs.shape[0]
        nlonb = lonbs.shape[0]
    except:
        lons,lats,lonbs,latbs,nlon,nlat,nlonb,nlatb=cts.create_grid(output_dict['manual_grid_option'])

    if output_dict['is_thd']:
        p_full,p_half,npfull,nphalf=cts.create_pressures()
    else:
        p_full=None
        p_half=None
        npfull=None
        nphalf=None

    #create times
    try:
        output_dict['num_years']
    except KeyError:
        is_climatology=True
    else:
        if output_dict['num_years']==1:
            is_climatology=True
        else:
            is_climatology=False
        num_years=output_dict['num_years']

    time_arr,day_number,ntime,time_units,time_bounds=cts.create_time_arr(num_years,is_climatology,output_dict['time_spacing_days'])

    #Output it to a netcdf file. 
    file_name=output_dict['file_name']
    variable_name=output_dict['var_name']

    number_dict={}
    number_dict['nlat']=nlat
    number_dict['nlon']=nlon
    number_dict['nlatb']=nlatb
    number_dict['nlonb']=nlonb
    number_dict['npfull']=npfull
    number_dict['nphalf']=nphalf
    number_dict['ntime']=ntime

    data_out=dataset[field_name].load().data

    cts.output_to_file(data_out,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,variable_name,number_dict,time_bounds)


def output_data(dataset, model_params, variable_name_prefix, category_name):

    data_type_folder_name = dataset.data_type

    translate_table = str.maketrans(dict.fromkeys('!@#$/'))

    if dataset.dataset_id == 'diff':
        time_folder_name=str(dataset.attrs['start_file_1'])+'_'+str(dataset.attrs['end_file_1'])+'_minus_'+str(dataset.attrs['start_file_2'])+'_'+str(dataset.attrs['end_file_2'])
        directory='/scratch/sit204/derived_data/diffs/'+dataset.exp_name.translate(translate_table)+'/'+time_folder_name+'/'+data_type_folder_name+'/' 
    else:
        time_folder_name=str(dataset.start_file)+'_'+str(dataset.end_file)
        directory='/scratch/sit204/derived_data/exps/'+dataset.exp_name.translate(translate_table)+'/'+time_folder_name+'/'+data_type_folder_name+'/' 
                
    if not os.path.exists(directory):
        os.makedirs(directory)

    output_file_name = variable_name_prefix+'_'+category_name+'.nc'

    output_dataset = xar.Dataset(coords = dataset.coords)
    
    for key in list(dataset.keys()):
        if variable_name_prefix in key:
            output_dataset[key] = (dataset[key].dims, dataset[key].values)


    output_dataset.to_netcdf(path=directory+output_file_name)
