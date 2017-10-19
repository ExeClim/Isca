import xarray as xar
from netCDF4 import Dataset
import numpy as np
import calendar_calc as cal
import os
import sys
import pdb
import create_timeseries as cts
from scipy import stats

def read_data( base_dir, exp_name, start_file, end_file, avg_or_daily, topo_present, model='fms13', file_name=None):

    if model=='fms13':

        files_temp=[base_dir+'/'+exp_name+'/run%03d' % m for m in range(start_file, end_file+1)]
        if(topo_present):
            if avg_or_daily == 'monthly':
#                 extra='_interp.nc'
                extra='_interp_new_height.nc'            

            else:
                extra='_interp_new_height_temp.nc'            
        else:
            extra='.nc'

        thd_string = '/atmos_'+avg_or_daily+extra

        thd_files = [s + thd_string for s in files_temp]

        thd_files_exist=[os.path.isfile(s) for s in thd_files]
        
        thd_file_size = [os.path.getsize(s) for s in thd_files]

        mode_file_size = stats.mode(thd_file_size).mode[0]
        
        thd_files_too_small = [thd_file_size[s] < 0.75*mode_file_size for s in range(len(thd_files))]

        print(thd_files[0])

        if not(all(thd_files_exist)):
            raise EOFError('EXITING BECAUSE OF MISSING FILES', [thd_files[elem] for elem in range(len(thd_files_exist)) if not thd_files_exist[elem]])

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





    elif model=='jra55_blanca':

        thd_files=base_dir+'/'+exp_name+'/'+file_name+'.nc'
        print(thd_files)

        size_list = init(thd_files)

        try:
            da_3d = xar.open_mfdataset(    thd_files,
                        decode_times=False,  # no calendar so tell netcdf lib
                        # choose how data will be broken down into manageable
                        # chunks.
            #                        concat_dim='time',
                        chunks={'initial_time0_hours': size_list['ntime']//64,
                                'lon': size_list['nlons']//36,
                                'lat': size_list['nlats']})
        except ValueError:
            da_3d = xar.open_mfdataset(    thd_files,
                    decode_times=False,  # no calendar so tell netcdf lib
                    # choose how data will be broken down into manageable
                    # chunks.
        #                        concat_dim='time',
                    chunks={'initial_time0_hours': size_list['ntime']//64,
                            'g0_lon_3': size_list['nlons']//36,
                            'g0_lat_2': size_list['nlats']})
                            
        names_dict = {'initial_time0_hours':'time', 'PRMSL_GDS0_MSL':'slp', 'UGRD_GDS0_ISBL':'ucomp', 'lv_ISBL1':'pfull', 'HGT_GDS0_ISBL':'height', 'g0_lon_3':'lon', 'g0_lat_2':'lat'}

        for name in names_dict.keys():
            try:                        
                da_3d.rename({name:names_dict[name]}, inplace=True)
            except ValueError:
                pass            

        time_arr = da_3d.time
        date_arr = cal.day_number_to_date(time_arr, calendar_type='standard', units_in='hours since 1800-01-01 00:00:00')

    elif model=='jra55':

        thd_files=base_dir+'/'+exp_name+'/'+file_name+'.nc'
        print(thd_files)

        size_list = init(thd_files)        
        if size_list['ntime']<64:
            chunks_use = {'time': size_list['ntime'],
                        'lon': size_list['nlons']//36,
                        'lat': size_list['nlats']}
        else:
            chunks_use = {'time': size_list['ntime']//64,
                        'lon': size_list['nlons']//36,
                        'lat': size_list['nlats']}
        
        da_3d = xar.open_mfdataset(    thd_files,
                decode_times=False,  # no calendar so tell netcdf lib
                # choose how data will be broken down into manageable
                # chunks.
    #                        concat_dim='time',
                chunks=chunks_use)

        names_dict = {'var2':'slp', 'var33':'ucomp', 'var7':'height', 'var34':'vcomp', 'var11':'temp', 'var39':'omega', 'lev':'pfull'}

        for name in names_dict.keys():
            try:                        
                da_3d.rename({name:names_dict[name]}, inplace=True)
                if name=='lev':
                    da_3d.coords['pfull'] = (('pfull'), da_3d.pfull.values/100.)
            except ValueError:
                pass            

        time_arr = da_3d.time
        date_arr = cal.day_number_to_date(time_arr, calendar_type='standard', units_in='hours since 1958-01-01 00:00:00')

    elif model=='hadisst':

        thd_files=base_dir+'/'+exp_name+'/'+file_name+'.nc'
        print(thd_files)

        size_list = init(thd_files)

        try:
            da_3d = xar.open_mfdataset(    thd_files,
                        decode_times=False,  # no calendar so tell netcdf lib
                        # choose how data will be broken down into manageable
                        # chunks.
            #                        concat_dim='time',
                        chunks={'time': size_list['ntime']//8,
                                'longitude': size_list['nlons']//36,
                                'latitude': size_list['nlats']})
            da_3d.rename({'sst':'t_surf', 'latitude':'lat', 'longitude':'lon'}, inplace=True)
                                
        except:
            da_3d = xar.open_mfdataset(    thd_files,
                        decode_times=False,  # no calendar so tell netcdf lib
                        # choose how data will be broken down into manageable
                        # chunks.
            #                        concat_dim='time',
                        chunks={'time': size_list['ntime']//8,
                                'lon': size_list['nlons']//36,
                                'lat': size_list['nlats']})        
            da_3d.rename({'sst':'t_surf'}, inplace=True)

        time_arr = da_3d.time

        if 'jra' in thd_files and '55' in thd_files:

            date_arr = cal.day_number_to_date(time_arr, calendar_type=da_3d.time.attrs['calendar'], units_in=da_3d.time.attrs['units'])

            names_dict = {'var2':'slp', 'var33':'ucomp', 'var7':'height', 'var34':'vcomp', 'var11':'temp', 'var39':'omega', 'lev':'pfull'}
        
            for name in names_dict.keys():
                try:                        
                    da_3d.rename({name:names_dict[name]}, inplace=True)
                    if name=='lev':
                        da_3d.coords['pfull'] = (('pfull'), da_3d.pfull.values/100.)
                except ValueError:
                    pass                                
                   
        else:
            date_arr = cal.day_number_to_date(time_arr, calendar_type='standard', units_in='days since 1870-01-01 00:00:00')        
            nao_pc1  = xar.open_mfdataset( base_dir+ '/jra_55/xarray_derived_fields/nao_jra_55_pc1_monthly.nc')
            snao_pc1 = xar.open_mfdataset( base_dir+ '/jra_55/xarray_derived_fields/snao_jra_55_pc1_monthly.nc')

            da_3d['nao_pc1']=(('time'),nao_pc1.nao_pc1)
            da_3d['snao_pc1']=(('time'),snao_pc1.snao_pc1)
        
#         ucomp_250_file  = xar.open_mfdataset( base_dir+ '/jra_55/xarray_derived_fields/ucomp_250_lon_280_40_eof_1958_2016_monthly.nc')
# 
#         idx_up_to_2012 = np.where(ucomp_250_file['years'].values < 2013)
# 
#         da_3d['ucomp_250_lon_280_40_winter_pc1']=(('time'),ucomp_250_file['ucomp_250_lon_280_40_winter_pc1'][idx_up_to_2012])
#         da_3d['ucomp_250_lon_280_40_summer_pc1']=(('time'),ucomp_250_file['ucomp_250_lon_280_40_summer_pc1'][idx_up_to_2012])        
        
#         da_3d['lon'] = (('lon'), np.mod(da_3d.lon.values+360.,360.))


    elif model=='era_interim':

        thd_files=base_dir+'/'+exp_name+'/'+file_name+'.nc'
        print(thd_files)

        size_list = init(thd_files)

        da_3d = xar.open_mfdataset(    thd_files,
                decode_times=False,  # no calendar so tell netcdf lib
                # choose how data will be broken down into manageable
                # chunks.
    #                        concat_dim='time',
                chunks={'time': size_list['ntime']//4,
                        'longitude': size_list['nlons']//36,
                        'latitude': size_list['nlats']})
        da_3d.rename({ 'latitude':'lat', 'longitude':'lon'}, inplace=True)

        try:
            da_3d.rename({'level':'pfull'}, inplace=True)
        except KeyError:                    
            da_3d.coords['pfull'] = (('pfull'), [100.])

        try:        
            temp_data = da_3d['t'].values
            temp_data = temp_data[:,np.newaxis,:,:]
            da_3d['temp'] = (('time', 'pfull', 'lat', 'lon'), temp_data)
        
            vcomp_data = da_3d['v'].values
            vcomp_data = vcomp_data[:,np.newaxis,:,:]
            da_3d['vcomp'] = (('time', 'pfull', 'lat', 'lon'), vcomp_data)
        except KeyError:
            try:
                da_3d.rename({'u':'ucomp', 'v':'vcomp'}, inplace=True)
            except ValueError:
                da_3d.rename({'u':'ucomp'}, inplace=True)
                
        time_arr = da_3d.time
        date_arr = cal.day_number_to_date(time_arr, calendar_type='standard', units_in='hours since 1900-01-01 00:00:00')

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

def read_land( base_dir,exp_name,land_present,topo_present,size_list,land_file='input/land.nc', lats_in = None):
    "Function for reading in land mask."
    
    #Create thd (3D) and twd (2D) data stores.
    land_array=np.zeros((size_list['nlats'],size_list['nlons']))
    topo_array=np.zeros((size_list['nlats'],size_list['nlons']))
    
    if land_file is not None:
        if(land_present or topo_present):
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
                    
            if topo_present:
                topo_array=fh.variables['zsurf'][:]
            fh.close()

    return land_array, topo_array

def diff_time_adjust(dataset,dataset_2):

    ntime_1=len(dataset.time.to_index())
    ntime_2=len(dataset_2.time.to_index())

    if ntime_1!=ntime_2:
        print('FATAL - total number of days in dataset_1 and dataset_2 are different. In dataset 1 and 2 there are ',ntime_1,ntime_2,' times, respectively.')
        sys.exit()

    time_arr_1=dataset.time.to_index()
    time_arr_2=dataset_2.time.to_index()

    doy_arr_1=dataset.dayofyear.to_index()
    doy_arr_2=dataset_2.dayofyear.to_index()

    if(all(time_arr_1!=time_arr_2)):

        if doy_arr_1[0]==doy_arr_2[0]:
#            del dataset['years']
#            del dataset['seq_months']
#            del dataset['seq_seasons']
            del dataset_2['years']
            del dataset_2['seq_months']
            del dataset_2['seq_seasons']
            del dataset_2['years_ax']
            del dataset_2['seq_seasons_ax']

            print('adjusting times in dataset_2 to those in dataset_1 to make diff possible.')

            dataset_2.coords['time']=(('time'),time_arr_1)

            dataset_2.coords['years']=(('time'),dataset.years.to_index())
            dataset_2.coords['seq_months']=(('time'),dataset.seq_months.to_index())
            dataset_2.coords['seq_seasons']=(('time'),dataset.seq_seasons.to_index())
            dataset_2.coords['years_ax'] = (('years_ax'),dataset.years_ax.to_index())
            dataset_2.coords['seq_seasons_ax'] = (('seq_seasons_ax'),dataset.seq_seasons_ax.to_index())



        else:
            print('FATAL - dates in DIFF file not lined up. First dayofyear in dataset ',doy_arr_1[0],' First dayofyear in dataset_2 ', doy_arr_2[0])
            sys.exit()

    dataset_diff=dataset-dataset_2

    return dataset_diff

def output_nc_file(dataset, field_name, model_params, output_dict):

    #create grid

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

def load_in_existing_data(dataset, model_params, variable_name_prefix, category_name, data_type_input = None):

    if data_type_input is None:
        data_type_folder_name = dataset.data_type
    else:
        data_type_folder_name = data_type_input

    translate_table = str.maketrans(dict.fromkeys('!@#$/'))

    if dataset.dataset_id == 'diff':
        time_folder_name=str(dataset.attrs['start_file_1'])+'_'+str(dataset.attrs['end_file_1'])+'_minus_'+str(dataset.attrs['start_file_2'])+'_'+str(dataset.attrs['end_file_2'])
        directory='/scratch/sit204/derived_data/diffs/'+dataset.exp_name.translate(translate_table)+'/'+time_folder_name+'/'+data_type_folder_name+'/' 
    else:
        time_folder_name=str(dataset.start_file)+'_'+str(dataset.end_file)
        directory='/scratch/sit204/derived_data/exps/'+dataset.exp_name.translate(translate_table)+'/'+time_folder_name+'/'+data_type_folder_name+'/' 

    input_file_name = variable_name_prefix+'_'+category_name+'.nc'

    try:
        input_dataset = xar.open_dataset(directory + input_file_name)
    except ValueError:
        input_dataset = xar.open_dataset(directory + input_file_name, decode_times=False)
        
    for key in list(input_dataset.keys()):
        if variable_name_prefix in key:
            dataset[key] = (input_dataset[key].dims, input_dataset[key].values)

    input_dataset.close()

def read_in_datasets(exp_details):

    #***********Read data*************************
    dataset, time_arr, size_list = read_data( exp_details.base_dir, exp_details.exp_name, exp_details.start_file, exp_details.end_file, exp_details.avg_or_daily, exp_details.topo_present, exp_details.model, exp_details.file_name)
    
    land_array, topo_array = read_land(exp_details.input_dir, exp_details.base_exp_name, exp_details.land_present, exp_details.topo_present, size_list, exp_details.land_file, dataset.lat.values)
    dataset['land'] = (('lat','lon'),land_array)
    dataset['topo'] = (('lat','lon'),topo_array)

    dataset.attrs['dataset_id'] = 'first'

    if exp_details.avg_or_daily != 'monthly':
        try:
            dataset_monthly, time_arr_monthly, size_list_monthly = read_data( exp_details.base_dir, exp_details.exp_name, exp_details.start_file, exp_details.end_file, 'monthly', exp_details.topo_present, exp_details.model, exp_details.file_name)
            dataset_monthly['land'] = (('lat','lon'),land_array)
            dataset_monthly['topo'] = (('lat','lon'),topo_array)
        except EOFError:
            dataset_monthly=None
    else:
        dataset_monthly=None

    if exp_details.exp_name_2 is None:
        print('running analysis on one experiment only')
        one_exp=True
        dataset_2 = None
        dataset_monthly_2 = None
        dataset_diff = None
        dataset_monthly_diff = None
    else:
        dataset_2, time_arr2, size_list = read_data( exp_details.base_dir, exp_details.exp_name_2, exp_details.start_file_2, exp_details.end_file_2, exp_details.avg_or_daily, exp_details.topo_present, exp_details.model, exp_details.file_name)
        dataset_2['land'] = (('lat','lon'),land_array)
        dataset_2['topo'] = (('lat','lon'),topo_array)
        one_exp=False        

        dataset_2.attrs['dataset_id'] = 'second'

        if exp_details.avg_or_daily != 'monthly':
            dataset_monthly_2, time_arr_monthly_2, size_list_monthly_2 = read_data( exp_details.base_dir, exp_details.exp_name_2, exp_details.start_file_2, exp_details.end_file_2,'monthly', exp_details.topo_present, exp_details.model, exp_details.file_name)
            dataset_monthly_2['land'] = (('lat','lon'),land_array)
            dataset_monthly_2['topo'] = (('lat','lon'),topo_array)
        else:
            dataset_monthly_2=None

        if((time_arr2 == time_arr).to_index().all()):
            dataset_diff=diff_time_adjust(dataset, dataset_2)
            dataset_diff['land'] =(('lat','lon'),land_array)
            dataset_diff['topo'] =(('lat','lon'),topo_array)
            dataset_diff.attrs['exp_name']='diff_'+exp_details.exp_name+'_minus_'+exp_details.exp_name_2

            dataset_diff.attrs['dataset_id'] = 'diff'
            
            if dataset.data_type == dataset_2.data_type:
                dataset_diff.attrs['data_type']  = dataset.data_type
            else:
                raise AttributeError('How can data with different time intervals be safely compared? Not yet implemented.')
        
            dataset_diff.attrs['start_file_1']= dataset.start_file
            dataset_diff.attrs['end_file_1']  = dataset.end_file
            dataset_diff.attrs['start_file_2']= dataset_2.start_file
            dataset_diff.attrs['end_file_2']  = dataset_2.end_file

            if exp_details.avg_or_daily != 'monthly':
                dataset_monthly_diff=diff_time_adjust(dataset_monthly, dataset_monthly_2)
                dataset_monthly_diff['land'] =(('lat','lon'),land_array)
                dataset_monthly_diff.attrs['exp_name']='diff_'+exp_details.exp_name+'_minus_'+exp_details.exp_name_2
            else:
                dataset_monthly_diff=None
    #************End Read data********************

    return dataset, dataset_2, dataset_diff, dataset_monthly, dataset_monthly_2, dataset_monthly_diff, one_exp

def close_datasets(dataset, dataset_2, dataset_diff, dataset_monthly, dataset_monthly_2, dataset_monthly_diff):

    dataset_list = [dataset, dataset_2, dataset_diff, dataset_monthly, dataset_monthly_2, dataset_monthly_diff]

    for dataset_obj in dataset_list:
        if dataset_obj is not None:
            dataset_obj.close()

class exp_object(object):
    def __init__(self, input_dir, base_dir, land_file, base_exp_name, exp_name, exp_name_2, land_present, topo_present, avg_or_daily, start_file, end_file, start_file_2, end_file_2, model, file_name):
        self.input_dir = input_dir
        self.base_dir = base_dir
        self.land_file = land_file
        self.base_exp_name = base_exp_name
        self.exp_name = exp_name
        self.exp_name_2 = exp_name_2
        self.land_present = land_present
        self.topo_present = topo_present
        self.avg_or_daily = avg_or_daily
        self.start_file = start_file
        self.end_file = end_file
        self.start_file_2 = start_file_2
        self.end_file_2 = end_file_2
        self.model = model
        self.file_name = file_name
        
def create_exp_object(input_dir, base_dir, land_file, base_exp_name, exp_name, exp_name_2, land_present, topo_present, avg_or_daily, start_file, end_file, start_file_2, end_file_2, model='fms13', file_name = None):

    final_exp_object = exp_object(input_dir, base_dir, land_file, base_exp_name, exp_name, exp_name_2, land_present, topo_present, avg_or_daily, start_file, end_file, start_file_2, end_file_2, model, file_name)
    
    return final_exp_object

