# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num
import pdb
import create_timeseries as cts
import xarray as xar
from mpl_toolkits.basemap import shiftgrid
import matplotlib.pyplot as plt

def add_sst_anomaly(sst_in, anomaly_type=None):
 
    
    if anomaly_type=='el_nino':
        el_nino_dataset_location = '/scratch/sit204/HadISST/hadisst_el_nino_3p4_analysis_with_rolling_mean.nc'
        dataset_el_nino          = xar.open_dataset(el_nino_dataset_location)
    
        start_month_idx = 472
        end_month_idx = 483
        sst_anom_orig_order = dataset_el_nino['t_surf_anom'][start_month_idx:(end_month_idx+1),...]
        
        print(('Selecting months between ' + str(dataset_el_nino.time.values[start_month_idx]) + ' to ' + str(dataset_el_nino.time.values[end_month_idx])))


        lat_range_keep = [-30., 15.]
        lon_range_keep = [-180., -70.]
        taper_length = 15.

        masked_sst_anom_orig_order = apply_lat_lon_mask(sst_anom_orig_order, lat_range_keep, lon_range_keep, taper_length)


        masked_sst_anom_orig_order = masked_sst_anom_orig_order.fillna(0.)

        lon_shift = np.around(dataset_el_nino.lon.mean().values/180.)*180.

        if lon_shift==0.:
        
            lon_shift_practical = lon_shift+0.001 
            #Added a little extra so that it doesn't start the grid at -0.5, but rather starts it at +0.5, as in the bcs data.
        
            sst_anom_orig_order_shift_lon,lons = shiftgrid(lon_shift_practical,masked_sst_anom_orig_order,dataset_el_nino.lon.values)
        else: 
            sst_anom_orig_order_shift_lon = masked_sst_anom_orig_order
        
        sst_anom_all_shifted = sst_anom_orig_order_shift_lon[:,::-1,:]

        print(' reordering months to line up with input file')
                
        start_month = dataset_el_nino.months.values[start_month_idx]
        end_month = dataset_el_nino.months.values[end_month_idx]
        
        sst_anom = np.zeros_like(sst_anom_orig_order)

        for month_tick in range(12):
        
            orig_order_idx = np.mod(month_tick - (start_month-1), 12)       
            sst_anom[month_tick, ...] = sst_anom_all_shifted[orig_order_idx]
        
        

        sst_with_anomaly = sst_in + sst_anom

    return sst_with_anomaly, lons

def apply_lat_lon_mask( unmasked_input, lat_range, lon_range_in, taper_length, power = 5):

    width         = np.abs(lon_range_in[1]-lon_range_in[0])
    central_point = (lon_range_in[1]+lon_range_in[0])/2.
    
    lon_range = [-width/2., width/2.]

    zeros_mask = xar.DataArray(np.zeros_like(unmasked_input.values), coords=unmasked_input.coords, dims=unmasked_input.dims)
    untapered_mask = xar.DataArray(np.zeros_like(unmasked_input.values), coords=unmasked_input.coords, dims=unmasked_input.dims)

    lat_array = unmasked_input.lat.values
    lon_array = unmasked_input.lon.values


    for lat_tick in range(len(lat_array)):
        
        lat_value = lat_array[lat_tick]
        
        for lon_tick in range(len(lon_array)):
            
            lon_value = lon_array[lon_tick]
            
            if (lat_value < lat_range[1] and lat_value > lat_range[0]) and (lon_value < lon_range[1] and lon_value > lon_range[0]):
                zeros_mask[:,lat_tick, lon_tick] = 1.
                untapered_mask[:,lat_tick, lon_tick] = 1.
                #All of those points inside the un-tapered box have the mask set to 1.
                
            elif (lat_value < lat_range[1]+taper_length and lat_value > lat_range[0]-taper_length) and (lon_value < lon_range[1]+taper_length and lon_value > lon_range[0]-taper_length) and not (lat_value < lat_range[1] and lat_value > lat_range[0]):    
             
                min_distance = np.min([np.abs(lat_value-lat_range[1]), np.abs(lat_value - lat_range[0])])
                zeros_mask[:,lat_tick, lon_tick] = (1. -((min_distance)/(taper_length)))**power
                #This are the points that are within the larger tapered box, but in a latitude band outside of the untapered box.

                if not (lon_value < lon_range[1] and lon_value > lon_range[0]) and not (lat_value < lat_range[1] and lat_value > lat_range[0]) :
                    min_distance = np.min([np.abs(lon_value-lon_range[1]), np.abs(lon_value - lon_range[0])])
                    zeros_mask[:,lat_tick, lon_tick] = zeros_mask[:,lat_tick, lon_tick] * (1. -((min_distance)/(taper_length)))**power
                    #These are the corners - so a more longitudinally restricted version of the above

            elif (lon_value < lon_range[1]+taper_length and lon_value > lon_range[0]-taper_length) and (lat_value < lat_range[1] and lat_value > lat_range[0]):    
             
                min_distance = np.min([np.abs(lon_value-lon_range[1]), np.abs(lon_value - lon_range[0])])
                zeros_mask[:,lat_tick, lon_tick] = (1. -((min_distance)/(taper_length)))**power
                #The final set of points is within the larger box, but within the latitude range of the original box.


        #Added a little extra so that it doesn't start the grid at -0.5, but rather starts it at +0.5, as in the bcs data.
     
    zeros_mask_shifted,lons = shiftgrid(-180.-central_point+0.5, zeros_mask, zeros_mask.lon.values)                    
    untapered_mask_shifted ,lons = shiftgrid(-180.-central_point+0.5, untapered_mask, untapered_mask.lon.values)                    
 
    final_mask = xar.DataArray(zeros_mask_shifted, coords=unmasked_input.coords, dims=unmasked_input.dims)
    final_untapered_mask = xar.DataArray(untapered_mask_shifted, coords=unmasked_input.coords, dims=unmasked_input.dims)

    masked_sst = unmasked_input * final_mask


    return masked_sst



def main():

    base_directory='/scratch/sit204/sst_amip_files/'


    amip_data_version='amip_data_version_1_1_0' #s 'amip_data_version_1_1_0' or 'amip_data_version_1_0_0'

    output_name_list  ={'tosbcs':'sst', 'siconc':'siconc'}
    #Note that we are using the bcs (boundary conditions) input4MIPs files, as instructed.
    # The theory is that by using the bcs files (which are mid-month values) the time-average
    # of the interpolated bcs files should be equal to the time-average data provided in 'tos'
    # files, not the 'tosbcs'. See http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amip2bcs.php
    # and http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amipbc_dwnld_files/360x180/v1.0.0/nc/readme_nc

    add_anomaly = False
#     anomaly_type='el_nino'
    months_to_include='all'

    for variable_name in list(output_name_list.keys()):

        if amip_data_version=='amip_data_version_1_0_0':
            nfiles=50
            folder_name='/1950_1999/'
            filename_prefix='amipbc_sst_360x180_19'
            sst_all=np.zeros((nfiles,12,180,360))
            do_annual_mean=True
        elif amip_data_version=='amip_data_version_1_1_0':
            nfiles=1
            folder_name=''
            filename_prefix=variable_name+'_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-0_gs1x1_187001-201512'
            do_annual_mean=False
        elif amip_data_version=='hadgem_t_surf':
            nfiles=1
            folder_name=''
            filename_prefix='ts_clim'
            do_annual_mean=False        

        for file_tick in range(nfiles):

            if nfiles!=1:
                filename=filename_prefix+str(file_tick+50)
            else:
                filename=filename_prefix

            resolution_file = Dataset(base_directory+amip_data_version+'/'+folder_name+'/'+filename+'.nc', 'r')

            try:
                lons = resolution_file.variables['longitude'][:]
                lats = resolution_file.variables['latitude'][:]
            except KeyError:
                lons = resolution_file.variables['lon'][:]
                lats = resolution_file.variables['lat'][:]

            sst_in = resolution_file.variables[variable_name][:]

            try:
                sst_all[file_tick,:,:,:]=sst_in
            except NameError:
                sst_all=sst_in
            except IndexError:
                sst_all=sst_in

        try:
            no_latb_lonb = False
            lonbs = resolution_file.variables['bounds_longitude'][:]
            latbs = resolution_file.variables['bounds_latitude'][:]
        except KeyError:
            try:
                lonbs = resolution_file.variables['lon_bnds'][:]
                latbs = resolution_file.variables['lat_bnds'][:]
            except:
                no_latb_lonb=True

        nlon=lons.shape[0]
        nlat=lats.shape[0]

        if not no_latb_lonb:
            nlonb=lonbs.shape[0]
            nlatb=latbs.shape[0]

            lonbs_adjusted=np.zeros(nlonb+1)
            latbs_adjusted=np.zeros(nlatb+1)

            lonbs_adjusted[0:nlonb]=lonbs[:,0]
            lonbs_adjusted[nlonb]=lonbs[-1,1]

            latbs_adjusted[0:nlatb]=latbs[:,0]
            latbs_adjusted[nlatb]=latbs[-1,1]
        else:
            latbs_adjusted = None
            lonbs_adjusted = None

        try:
            day_number = resolution_file.variables['time'][:]
        except:
            day_number = np.ones(12)   

        time_arr = day_number_to_date(day_number, calendar_type = 'gregorian', units_in = 'days since 1870-1-1')
        time_arr_adj=np.arange(15,360,30)

        annual_mean_name=''

        if len(sst_all.shape)==4:
            sst_in=np.mean(sst_all,axis=0)
        else:
            sst_in=np.zeros((12,nlat,nlon))
            
            if months_to_include=='all':
                for month_tick in np.arange(1,13,1):
                    month_idx = np.where(time_arr.month==month_tick)[0]
                    sst_in[month_tick-1,:,:]=np.mean(sst_all[month_idx,:,:],axis=0)
            
            elif months_to_include=='DJF':
                djf_idx = np.where(np.logical_or(np.logical_or(time_arr.month==1, time_arr.month==2), time_arr.month==12))
                djf_mean = np.mean(sst_all[djf_idx[0],...], axis=0)
                for month_tick in np.arange(1,13,1):
                    sst_in[month_tick-1,...]=djf_mean
                annual_mean_name='_djf'
            
            elif months_to_include=='only_month_available':
                for month_tick in np.arange(1,13,1):
                    month_idx = np.where(time_arr.month==month_tick)[0]
                    sst_in[month_tick-1,:,:]=sst_all

        if do_annual_mean:
            sst_in_am=np.mean(sst_in,axis=0)
            sst_in=np.zeros((12,nlat,nlon))
            for month_tick in np.arange(1,13,1):
                sst_in[month_tick-1,:,:]=sst_in_am
            annual_mean_name='_am'

        if add_anomaly and variable_name=='tosbcs':
            sst_in, shifted_lons = add_sst_anomaly(sst_in, anomaly_type)
            anom_name = '_'+anomaly_type
        else:
            anom_name = ''

        p_full=None
        p_half=None

        npfull=None
        nphalf=None

        #Find grid and time numbers

        ntime=time_arr.day.shape[0]
        if not no_latb_lonb:        
            nlonb=lonbs_adjusted.shape[0]
            nlatb=latbs_adjusted.shape[0]

        #Output it to a netcdf file. 
        variable_name=output_name_list[variable_name]+annual_mean_name+'_clim_'+amip_data_version[0:5]+anom_name        
        file_name=variable_name+'_'+amip_data_version+'.nc'


        number_dict={}
        number_dict['nlat']=nlat
        number_dict['nlon']=nlon
        number_dict['npfull']=npfull
        number_dict['nphalf']=nphalf
        number_dict['ntime']=ntime

        if not no_latb_lonb:
            number_dict['nlatb']=nlatb
            number_dict['nlonb']=nlonb

        time_units='days since 0000-01-01 00:00:00.0'

        cts.output_to_file(sst_in,lats,lons,latbs_adjusted,lonbs_adjusted,p_full,p_half,time_arr_adj,time_units,file_name,variable_name,number_dict)

if __name__=="__main__":
    main()



