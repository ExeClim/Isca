"""A method for calculating seasonally-varying qfluxes, as described in Russell et al 1985 DOI:10.1016/0377-0265(85)90022-3"""
    
import numpy as np
import xarray
from xarray import ufuncs as xruf
import time
from scipy import stats
from mpl_toolkits.basemap import shiftgrid
import area_average as aav
import nc_file_io_xarray as io
import matplotlib.pyplot as plt
import os
import pdb

__author__='Stephen Thomson'


def qflux_calc(dataset, model_params, output_file_name, ice_file_name=None, groupby_name='months'):

    if groupby_name=='months':
        time_varying_ice = ice_mask_calculation(dataset, dataset.land, ice_file_name)
        upper_ocean_heat_content(dataset, model_params, time_varying_ice)
        net_surf_energy_flux(dataset, model_params)
        deep_ocean_heat_content(dataset, model_params)
        ocean_transport(dataset, model_params)

        output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':output_file_name+'.nc', 'var_name':output_file_name} #Have specified that var name is the same as file name as this is what the fortran assumes.
        
    elif groupby_name=='dayofyear':
        time_varying_ice = ice_mask_calculation(dataset, dataset.land, ice_file_name)
        upper_ocean_heat_content(dataset, model_params, time_varying_ice, dayofyear_or_months='dayofyear')
        net_surf_energy_flux(dataset, model_params)
        deep_ocean_heat_content(dataset, model_params)
        ocean_transport(dataset, model_params)

        output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':output_file_name+'.nc', 'var_name':output_file_name}    
        
    elif groupby_name=='all_time':
        time_varying_ice = ice_mask_calculation(dataset, dataset.land, ice_file_name, dayofyear_or_months=groupby_name)
        upper_ocean_heat_content(dataset, model_params, time_varying_ice, dayofyear_or_months=groupby_name)
        net_surf_energy_flux(dataset, model_params, dayofyear_or_months=groupby_name)
        deep_ocean_heat_content(dataset, model_params, dayofyear_or_months=groupby_name)
        ocean_transport(dataset, model_params, dayofyear_or_months=groupby_name)
        regrid_in_time(dataset, groupby_name)

        output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':output_file_name+'.nc', 'var_name':output_file_name}            
        
    io.output_nc_file(dataset,'masked_ocean_transport', model_params, output_dict)

def time_gradient(data_in, delta_t):

    data_out=np.gradient(data_in, delta_t)

    return data_out
    
def ice_mask_calculation(dataset, land_array, ice_file_name, dayofyear_or_months='months'):

    try:
        ice_climatology=dataset['ice_conc'].groupby(dayofyear_or_months).mean('time').load()
        ice_array=ice_climatology.values
        ice_idx=ice_array !=0.
        ice_array[ice_idx]=1.0
        time_varying_ice=True
        print('have gotten ice concentration from climatology')
    except KeyError:

        try:
            ice_data_temp         = xarray.open_dataset(ice_file_name, decode_times=False)
            albedo_constant_value = ice_data_temp['albedo']
            albedo_array          = albedo_constant_value.values.squeeze()
            ice_idx               = np.round(albedo_array, decimals=2) == 0.7
            ice_array=np.zeros_like(albedo_array)
            ice_array[ice_idx]    = 1.0
            time_varying_ice      = False
            print('have gotten ice concentration from one month albedo')
        except TypeError:
            ice_array=np.zeros_like(land_array)
            time_varying_ice = False
            print('no ice climatology')

    land_ice_mask=np.zeros_like(ice_array)

    if time_varying_ice:
        dataset['ice_mask']=((dayofyear_or_months+'_ax','lat','lon'),ice_array)

        for month_tick in range(12):
            land_ice_mask_temp=land_array.values+dataset['ice_mask'].values[month_tick,:,:]
            two_idx=land_ice_mask_temp>1.
            land_ice_mask_temp[two_idx]=1.0
            land_ice_mask[month_tick,:,:]=land_ice_mask_temp

        dataset['land_ice_mask']=((dayofyear_or_months+'_ax','lat','lon'),land_ice_mask)

    else:
        dataset['ice_mask']=(('lat','lon'),ice_array)
        land_ice_mask_temp=land_array.values+dataset['ice_mask'].values
        two_idx=land_ice_mask_temp>1.
        land_ice_mask_temp[two_idx]=1.0
        land_ice_mask[:,:]=land_ice_mask_temp
        dataset['land_ice_mask']=(('lat','lon'),land_ice_mask)

    return time_varying_ice

def upper_ocean_heat_content(dataset, model_params, time_varying_ice, dayofyear_or_months='months'):
    """Calculating upper-ocean heat content assuming a constant mixed layer depth, unlike Russel 1985, who have a seasonally-varying mixed layer depth.

    Note that dayofyear_or_months was designed only so that the time derrivatives of surface temperature could be calculated
    in a more time-resolved way with dayofyear data, then time-averaged onto months, with months used throughout the rest of
    of the code. However, if all_time option is used, this negates the need for a rate-of-change calulation anyway.
    
     """
    print('doing upper ocean heat content and rate of change calculation')
    sst_data=dataset['t_surf'].groupby(dayofyear_or_months).mean('time').load()

    if dayofyear_or_months=='months':
        dataset['sst_clim']=(('months_ax','lat','lon'),sst_data)
        sst_clim=dataset['sst_clim']
    else:
        dataset['sst_clim_'+dayofyear_or_months]=((dayofyear_or_months+'_ax','lat','lon'),sst_data)
        sst_clim=dataset['sst_clim_'+dayofyear_or_months]
        dataset['sst_clim']=((dayofyear_or_months+'_ax','lat','lon'),sst_data)

#    weighted_sst_data=model_params['ocean_rho']*model_params['ocean_cp']*model_params['ml_depth']*sst_data*(1.0-dataset['land'])
    weighted_sst_data=model_params['ocean_rho']*model_params['ocean_cp']*model_params['ml_depth']*sst_clim*(1.0-dataset['land'])

    shape1=np.shape(weighted_sst_data)

    ny=shape1[1]
    nx=shape1[2]
    d_weighted_sst_data_dt=np.zeros_like(weighted_sst_data)

    if dayofyear_or_months=='dayofyear':
        delta_t=model_params['day_length']
    elif dayofyear_or_months=='months':
        delta_t=model_params['day_length']*30.

    if dayofyear_or_months!='all_time':
        for x in range(nx):
            for y in range(ny):
                d_weighted_sst_data_dt[:,y,x]=time_gradient(weighted_sst_data[:,y,x], delta_t)

    for time_tick in range(shape1[0]):
        if time_varying_ice:
            d_weighted_sst_data_dt[time_tick,:,:]=d_weighted_sst_data_dt[time_tick,:,:]*(1.0-dataset['land_ice_mask'].values[time_tick,:,:])
        else:
            d_weighted_sst_data_dt[time_tick,:,:]=d_weighted_sst_data_dt[time_tick,:,:]*(1.0-dataset['land_ice_mask'].values)


    if dayofyear_or_months=='dayofyear':
        
        dataset['d_weighted_sst_data_dt_days']=(('dayofyear_ax','lat','lon'), d_weighted_sst_data_dt)
        months_on_dayofyear_ax = dataset.months.groupby('dayofyear').mean('time').values
        dataset.coords['months_on_dayofyear_ax']=(('dayofyear_ax'),months_on_dayofyear_ax)

        monthly_values = dataset['d_weighted_sst_data_dt_days'].groupby('months_on_dayofyear_ax').mean('dayofyear_ax')
        dataset['d_weighted_sst_data_dt_months_from_days']=(('months_ax','lat','lon'), monthly_values)    

    elif dayofyear_or_months=='months':
        dataset['d_weighted_sst_data_dt']=(('months_ax','lat','lon'), d_weighted_sst_data_dt)    

    elif dayofyear_or_months=='all_time':
        dataset['d_weighted_sst_data_dt']=(('all_time_ax','lat','lon'), d_weighted_sst_data_dt)    

def net_surf_energy_flux(dataset, model_params, dayofyear_or_months='months'):
    """Calculates the net surface energy flux to be used in q-flux calcuation, but also calcualtes a scaling factor such that the annual average of the area-averaged surface flux is zero."""

    print('doing net surf energy flux')

    flux_sw_data=dataset['flux_sw'].groupby(dayofyear_or_months).mean('time')
    dataset['flux_sw_clim']=((dayofyear_or_months+'_ax','lat','lon'),flux_sw_data)

    flux_lw_data=dataset['flux_lw'].groupby(dayofyear_or_months).mean('time')
    dataset['flux_lw_clim']=((dayofyear_or_months+'_ax','lat','lon'),flux_lw_data)

    flux_t_data=dataset['flux_t'].groupby(dayofyear_or_months).mean('time')
    dataset['flux_t_clim']=((dayofyear_or_months+'_ax','lat','lon'),flux_t_data)

    flux_lhe_data=dataset['flux_lhe'].groupby(dayofyear_or_months).mean('time')
    dataset['flux_lhe_clim']=((dayofyear_or_months+'_ax','lat','lon'),flux_lhe_data)


    aav.area_average(dataset, 'flux_sw_clim', model_params, land_ocean_all='ocean_non_ice', axis_in=dayofyear_or_months+'_ax')
    aav.area_average(dataset, 'flux_lw_clim', model_params, land_ocean_all='ocean_non_ice', axis_in=dayofyear_or_months+'_ax')
    aav.area_average(dataset, 'sigma_sb_sst_clim', model_params, land_ocean_all='ocean_non_ice', axis_in=dayofyear_or_months+'_ax')
    aav.area_average(dataset, 'flux_t_clim', model_params, land_ocean_all='ocean_non_ice', axis_in=dayofyear_or_months+'_ax')
    aav.area_average(dataset, 'flux_lhe_clim', model_params, land_ocean_all='ocean_non_ice', axis_in=dayofyear_or_months+'_ax')

    scaling_factor_old=(((dataset['sigma_sb_sst_clim_area_av_ocean_non_ice']+dataset['flux_t_clim_area_av_ocean_non_ice']+dataset['flux_lhe_clim_area_av_ocean_non_ice']-dataset['flux_lw_clim_area_av_ocean_non_ice'])/dataset['flux_sw_clim_area_av_ocean_non_ice'])).mean(dayofyear_or_months+'_ax')

    scaling_factor=(((dataset['sigma_sb_sst_clim_area_av_ocean_non_ice'].mean(dayofyear_or_months+'_ax')+dataset['flux_t_clim_area_av_ocean_non_ice'].mean(dayofyear_or_months+'_ax')+dataset['flux_lhe_clim_area_av_ocean_non_ice'].mean(dayofyear_or_months+'_ax')-dataset['flux_lw_clim_area_av_ocean_non_ice'].mean(dayofyear_or_months+'_ax'))/dataset['flux_sw_clim_area_av_ocean_non_ice'].mean(dayofyear_or_months+'_ax')))

#    scaling_factor=1.0

    print('using scale factor for SW of '+ str(scaling_factor))

    net_surf_energy_fl=(scaling_factor*dataset['flux_sw_clim']+dataset['flux_lw_clim']-(model_params['sigma_sb']*dataset['sst_clim']**4.0)-dataset['flux_t_clim']-dataset['flux_lhe_clim'])*(1.0-dataset['land_ice_mask'])

    dataset['net_surf_energy_fl']=((dayofyear_or_months+'_ax','lat','lon'), net_surf_energy_fl)
    aav.area_average(dataset, 'net_surf_energy_fl', model_params, land_ocean_all='ocean_non_ice', axis_in=dayofyear_or_months+'_ax')
    
def deep_ocean_heat_content(dataset, model_params, dayofyear_or_months='months'):

    print('doing deep ocean heat content')
    aav.area_average(dataset, 'd_weighted_sst_data_dt', model_params, land_ocean_all='ocean_non_ice', axis_in=dayofyear_or_months+'_ax')
    d_deep_ocean_dt=(dataset['net_surf_energy_fl_area_av_ocean_non_ice']-dataset['d_weighted_sst_data_dt_area_av_ocean_non_ice'])

    dataset['d_deep_ocean_dt']=(dataset['net_surf_energy_fl_area_av_ocean_non_ice'].dims, d_deep_ocean_dt)    

def ocean_transport(dataset, model_params, dayofyear_or_months='months'):

    ocean_transport=dataset['d_weighted_sst_data_dt']+dataset['d_deep_ocean_dt']-dataset['net_surf_energy_fl']
    masked_ocean_transport=ocean_transport*(1.0-dataset['land_ice_mask'])
    masked_land_transport=ocean_transport*(dataset['land_ice_mask'])


    dataset['ocean_transport']=((dayofyear_or_months+'_ax','lat','lon'), ocean_transport)    
    dataset['masked_ocean_transport']=((dayofyear_or_months+'_ax','lat','lon'), masked_ocean_transport)
    dataset['masked_land_transport']=((dayofyear_or_months+'_ax','lat','lon'), masked_land_transport)    

def regrid_in_time(dataset, groupby_name):
    """ Simple routine to repeat arrays in all_time-averaging situations, so that monthly data still output."""

    ocean_heat_flux_shape = dataset['masked_ocean_transport'].shape

    if ocean_heat_flux_shape[0]!=12 and groupby_name=='all_time':
    
        dataset_to_repeat = dataset['masked_ocean_transport'].load()
        
        dataset_to_output = np.zeros((12, ocean_heat_flux_shape[1], ocean_heat_flux_shape[2]))
        
        for i in range(12):
            dataset_to_output[i,...] = dataset_to_repeat
        
        dataset['masked_ocean_transport'] = (('months_ax','lat','lon'), dataset_to_output)

def check_surface_flux_dims(dataset):
    ''' This surface flux checker is designed to decide if we're using grey rad or not. If we're using grey rad then the definition
    of flux_sw and flux_lw are different to RRTM. The script was written to use RRTM output, so it changes variable names etc to be 
    equivalent to RRTM definitions.
    '''

    flux_dims = dataset['flux_sw'].dims

    if 'phalf' in flux_dims:
        dataset.rename({'flux_sw':'flux_sw'+'_3d'}, inplace=True)
        max_pressure = dataset.phalf.max()
        flux_at_bottom_phalf_level = dataset['flux_sw_3d'].sel(phalf=max_pressure)
        new_dims = ('time','lat','lon')
        dataset['flux_sw'] = (new_dims, flux_at_bottom_phalf_level)

    flux_dims_lw = dataset['flux_lw'].dims

    if 'phalf' in flux_dims_lw:
        dataset.rename({'flux_lw':'flux_lw'+'_3d'}, inplace=True)
        try:
            # Script assumes flux_lw is the surface lw down (i.e. not a net flux). This is the case with RRTM, but with
            # grey radiation 'flux_lw' is the net lw flux in 3D. So we take the lwdn_sfc output from grey rad and rename it
            # flux_lw. 
            dataset['lwdn_sfc']
            dataset.rename({'lwdn_sfc':'flux_lw'}, inplace=True)
        except:
            #If lwdn_sfc is not available, then we re-calculate it from flux_lw by adding back sigma*t_surf**4, then call it flux_lw
            print('lwdn_sfc not present when using grey radiation, so re-calculating it from flux_lw.')
            max_pressure = dataset.phalf.max()
            lwdn_sfc = dataset.flux_lw_3d.sel(phalf=max_pressure) + sigma_sb*dataset.t_surf**4.
            new_dims = ('time','lat','lon')
            dataset['flux_lw'] = (new_dims, lwdn_sfc)

if __name__ == "__main__":

    import nc_file_io_xarray as io
    import set_and_get_params as sagp

    try:
        GFDL_BASE        = os.environ['GFDL_BASE']
        GFDL_WORK        = os.environ['GFDL_WORK']
        GFDL_DATA        = os.environ['GFDL_DATA']        
    except Exception as e:
        print('Environment variables GFDL_BASE, GFDL_WORK, GFDL_DATA must be set')
        exit(0)


    input_dir=GFDL_BASE
    base_dir=GFDL_DATA
    land_file='input/land.nc'
    base_exp_name='annual_mean_ice_post_princeton_fixed_sst/' #Folder containing the python script and input files that ran the experiment
    exp_name='annual_mean_ice_post_princeton_fixed_sst_1' #Folder within the data directory where the files can be found
#    ice_file_name=base_dir+'annual_mean_ice_albedo_change_test_mk2_4320_dt_rad_4/'+'run360/'+'atmos_monthly.nc'
    ice_file_name = '/scratch/sit204/data_isca/realistic_continents_fixed_sst_test_experiment_albedo/run0001/atmos_daily.nc'
    output_file_name='ami_test_interp' #Proposed name of your output qflux file. Will also be qflux field name in q-flux netcdf file as the fortran assumes file-name = field name. No need to add '.nc' or any file paths in this variable as otherwise they will end up in the field name too. Output file will be stored in the same directory as this script.

    start_file=240
    end_file=360
    land_present=True
    use_interpolated_pressure_level_data = False #Conditions the script on whether to expect data on sigma levels (if False) or pressure levels (if True). Script should be insensitive to this choice if both sets of files exist. 

    #Set time increments of input files (e.g. `monthly` for `atmos_monthly` files.
    avg_or_daily='monthly'

    #Set the time frequency of output data. Valid options are 'months', 'all_time' or 'dayofyear'.
    time_divisions_of_qflux_to_be_calculated='months'

    model_params = sagp.model_params_set(input_dir, delta_t=720., ml_depth=20., res=42)

    dataset, time_arr, size_list = io.read_data( base_dir,exp_name,start_file,end_file,avg_or_daily,use_interpolated_pressure_level_data)

    land_array, topo_array = io.read_land(input_dir,base_exp_name,land_present,use_interpolated_pressure_level_data,size_list,land_file)
    dataset['land'] = (('lat','lon'),land_array)

    check_surface_flux_dims(dataset)
    
    qflux_calc(dataset, model_params, output_file_name, ice_file_name, groupby_name=time_divisions_of_qflux_to_be_calculated)


