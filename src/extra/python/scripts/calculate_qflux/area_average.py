import set_and_get_params as sagp
import numpy as np
import pdb

__author__='Stephen Thomson'

def area_average(dataset, variable_name, model_params, land_ocean_all='all', level=None, axis_in='time', lat_range = None):

    print('performing area average on ',variable_name, 'of type ', land_ocean_all)

    if (variable_name[0:9]=='hc_scaled'):
        variable_name_use=variable_name[10:]
    elif(variable_name[0:8]=='sigma_sb'):
        variable_name_use=variable_name[9:]
    else:
        variable_name_use=variable_name

    if(level!=None):
        data_input=dataset[variable_name_use].sel(pfull=level, method='nearest')
    else:
        data_input=dataset[variable_name_use]


    if (variable_name[0:9]=='hc_scaled'):
        data_to_average=data_input*dataset['ml_heat_cap']/model_params['delta_t']
    elif(variable_name[0:8]=='sigma_sb'):
        data_to_average=model_params['sigma_sb']*data_input**4.
    else:
        data_to_average=data_input

    try:
        grid_area=dataset['grid_cell_area']
    except KeyError:
        sagp.get_grid_sizes(dataset,model_params)
        grid_area=dataset['grid_cell_area']


    if(land_ocean_all == 'land'):
        scaled_grid_area=grid_area*(dataset['land'])

        multiplied=scaled_grid_area*data_to_average
        average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))

    elif(land_ocean_all == 'ocean'):
        scaled_grid_area=grid_area*(1.-dataset['land'])

        multiplied=scaled_grid_area*data_to_average
        average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))

    elif(land_ocean_all == 'ocean_non_ice'):
        scaled_grid_area=grid_area*(1.-dataset['land_ice_mask'])

        multiplied=scaled_grid_area*data_to_average
        average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))

    elif(land_ocean_all == 'all'):
        multiplied=grid_area*data_to_average
        average=multiplied.sum(('lat','lon'))/grid_area.sum(('lat','lon'))

    elif(land_ocean_all == 'qflux_area'):
        scaled_grid_area=grid_area*(dataset['qflux_area'])

        multiplied=scaled_grid_area*data_to_average
        average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))

    elif(land_ocean_all[3:] == 'eur'):
        scaled_grid_area=grid_area*(dataset[land_ocean_all])

        multiplied=scaled_grid_area*data_to_average
        average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))

    elif(land_ocean_all == 'lat_range'):
        scaled_grid_area = grid_area.where((dataset.lat > lat_range[0]) & (dataset.lat < lat_range[1]))
        
        multiplied = scaled_grid_area * data_to_average
        average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))      
    else:
        print('invalid area-average option: ',land_ocean_all)
        return

    new_var_name=variable_name+'_area_av_'+land_ocean_all
    dataset[new_var_name]=((axis_in), average)
    
def european_area_av(dataset, model_params, eur_area_av_input):

    variables_list=eur_area_av_input['variables_list']
    try:
            levels_list  = eur_area_av_input['levels_list']
    except KeyError:
        levels_list  = None

    lats=dataset.lat
    lons=dataset.lon

    lon_array, lat_array = np.meshgrid(lons,lats)

    idx_nw_eur =     (45. <= lat_array) & (lat_array < 60.) & (-5. < lon_array) & (lon_array < 27.5)
    idx_nw_eur_neg = (45. <= lat_array) & (lat_array < 60.) & (np.mod(-5.,360.) < lon_array)

    idx_sw_eur =     (30. <= lat_array) & (lat_array < 45.) & (-5. < lon_array) & (lon_array < 27.5)
    idx_sw_eur_neg = (30. <= lat_array) & (lat_array < 45.) & (np.mod(-5.,360.) < lon_array)

    idx_ne_eur = (45. <= lat_array) & (lat_array < 60.) & (27.5 < lon_array) & (lon_array < 60.)
    idx_se_eur = (30. <= lat_array) & (lat_array < 45.) & (27.5 < lon_array) & (lon_array < 60.)

    idx_all_eur =     (35. <= lat_array) & (lat_array < 60.) & (-10. < lon_array) & (lon_array < 40.)
    idx_all_eur_neg = (35. <= lat_array) & (lat_array < 60.) & (np.mod(-10.,360.) < lon_array)


    land_nw_eur=np.zeros_like(dataset.land)
    land_sw_eur=np.zeros_like(dataset.land)

    land_ne_eur=np.zeros_like(dataset.land)
    land_se_eur=np.zeros_like(dataset.land)
    
    land_all_eur=np.zeros_like(dataset.land)
    
    land_nw_eur[idx_nw_eur]=1.0
    land_nw_eur[idx_nw_eur_neg]=1.0

    land_sw_eur[idx_sw_eur]=1.0
    land_sw_eur[idx_sw_eur_neg]=1.0

    land_ne_eur[idx_ne_eur]=1.0
    land_se_eur[idx_se_eur]=1.0

    land_all_eur[idx_all_eur]=1.0
    land_all_eur[idx_all_eur_neg]=1.0

    dataset['nw_eur']=(('lat','lon'), land_nw_eur)
    dataset['sw_eur']=(('lat','lon'), land_sw_eur)
    dataset['ne_eur']=(('lat','lon'), land_ne_eur)
    dataset['se_eur']=(('lat','lon'), land_se_eur)

    dataset['al_eur']=(('lat','lon'), land_all_eur)
    

    for i in range(np.shape(variables_list)[0]):
        var_name=variables_list[i]
        if levels_list!=None:
            level_in=levels_list[i]
        else:
            level_in=None

        area_average(dataset, var_name, model_params, land_ocean_all='nw_eur',level=level_in)
        area_average(dataset, var_name, model_params, land_ocean_all='sw_eur',level=level_in)
        area_average(dataset, var_name, model_params, land_ocean_all='ne_eur',level=level_in)
        area_average(dataset, var_name, model_params, land_ocean_all='se_eur',level=level_in)
        area_average(dataset, var_name, model_params, land_ocean_all='al_eur',level=level_in)

def qflux_area_av(dataset, model_params, qflux_area_av_input):

    qflux_area=np.zeros_like(dataset.land)

    variables_list     = qflux_area_av_input['variables_list']

    warmpool_lat_centre= qflux_area_av_input['lat_centre']
    warmpool_lon_centre= qflux_area_av_input['lon_centre']

    warmpool_width     = qflux_area_av_input['width']
    warmpool_width_lon = qflux_area_av_input['width_lon']

    lats=dataset.lat
    lons=dataset.lon

    latbs=dataset.latb
    lonbs=dataset.lonb


    for j in np.arange(len(lats)):
         lat = 0.5*(latbs[j+1] + latbs[j])
         lat = (lat-warmpool_lat_centre)/warmpool_width
         for i in np.arange(len(lons)):
              lon = 0.5*(lonbs[i+1] + lonbs[i])
              lon = (lon-warmpool_lon_centre)/warmpool_width_lon
              if( lat**2.+lon**2. <= 1.0 ):
                  qflux_area[j,i]=1.0

    dataset['qflux_area']=(('lat','lon'), qflux_area)

    for i in range(np.shape(variables_list)[0]):
        var_name=variables_list[i]
        original_axes_of_data = dataset[var_name].dims
        list_of_dims = list(original_axes_of_data)
        try:
            list_of_dims.remove('lat')        
            list_of_dims.remove('lon')            
        except:
            print('original data doesnt have lat lon as dimension - problem!')
        tuple_of_dims = tuple(list_of_dims)
        
        area_average(dataset, var_name, model_params, land_ocean_all='qflux_area', axis_in = tuple_of_dims)


