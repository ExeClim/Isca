import numpy as np
import pdb
import sys
sys.path.append('../')
import cell_area as carea

__author__='Stephen Thomson'

def model_params_set(input_dir,res=42, delta_t=900., radius=6376.0e3, day_length=86400., ocean_rho=1.035e3, ocean_cp=3989.24495292815, ml_depth=10., sigma_sb=5.6734e-8, pref=1000.,kappa=2./7., g=9.8, rct=6376.0e3, omega=7.292e-5, rdgas=287.04, lat_name='lat', lon_name='lon', latb_name='latb', lonb_name='lonb'):

    model_params={}

    model_params['input_dir']=input_dir
    model_params['res']=res
    model_params['delta_t']=delta_t
    model_params['planet_radius']=radius
    model_params['sigma_sb']=sigma_sb
    model_params['day_length']=day_length
    model_params['ocean_rho']=ocean_rho
    model_params['ocean_cp']=ocean_cp
    model_params['ml_depth']=ml_depth
    model_params['pref']=pref
    model_params['kappa']=kappa
    model_params['g']=g #Graviational acceleration
    model_params['rct']=rct # Rgct is the R/g constant
    model_params['rgct']=model_params['rct']/model_params['g']
    model_params['omega'] = omega #Planet rotation rate omega
    model_params['rdgas']=rdgas
    model_params['cp_air']=rdgas/kappa
    model_params['lat'] = lat_name
    model_params['lon'] = lon_name
    model_params['latb'] = latb_name
    model_params['lonb'] = lonb_name

    return model_params

def get_grid_sizes(dataset,model_params):

    area,x,y = carea.cell_area_from_xar(dataset, lat_name=model_params['lat'], lon_name = model_params['lon'], latb_name=model_params['latb'], lonb_name=model_params['lonb'], radius = model_params['planet_radius'])

#     area,x,y=carea.cell_area_all(model_params['res'],model_params['input_dir'])

#     area=(model_params['planet_radius']**2.)*area
#     x=(model_params['planet_radius'])*x
#     y=(model_params['planet_radius'])*y

    dataset['grid_cell_area']=(('lat','lon'),area)
    dataset['grid_cell_area_unscaled']=(('lat','lon'),area/(model_params['planet_radius']**2.))
    dataset['grid_cell_size_lon']=(('lat','lon'),x)
    dataset['grid_cell_size_lat']=(('lat','lon'),y)


