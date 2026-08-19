# -*- coding: utf-8 -*-s
from typing import NoReturn
import numpy as np
import pdb
import create_timeseries as cts
import xarray as xar
import gauss_grid as gg
import matplotlib.pyplot as plt
import windspharm as wsp
import pdb

def convert_to_vor_div(u_in, v_in, lat_arr, planet_radius):
    """convert spherical polar velocities to vor and div"""

    uwnd, uwnd_info = wsp.tools.prep_data(u_in, 'yx')
    vwnd, vwnd_info = wsp.tools.prep_data(v_in, 'yx')

    # It is also required that the latitude dimension is north-to-south. Again the
    # bundled tools make this easy.
    lat_1d_ordered, uwnd, vwnd = wsp.tools.order_latdim(lat_arr[:,0], uwnd, vwnd)

    # Create a VectorWind instance to handle the computation of streamfunction and
    # velocity potential.
    w = wsp.standard.VectorWind(uwnd, vwnd, rsphere=planet_radius, gridtype='gaussian')

    # Compute the streamfunction and velocity potential. Also use the bundled
    # tools to re-shape the outputs to the 4D shape of the wind components as they
    # were read off files.
    vor = w.vorticity()
    div = w.divergence()
    # sf, vp = w.sfvp()
    vor = wsp.tools.recover_data(vor, uwnd_info)
    div = wsp.tools.recover_data(div, uwnd_info)

    return vor[::-1,:], div[::-1,:] #need to reverse latitude reordering

def set_u_v_height_field(lon_in, lat_in, lonb_in, latb_in, epsilon, alpha, beta, m, r_0, planet_radius, northern_hemisphere=True):
    """Configure an initial condition for u, v and h given some
    balance condition. Use parameters and gradient-wind balance for Saturn
    from 10.1016/j.icarus.2017.06.006"""

    deformation_scale = 3200e3 #p62 of Rostami et al 2017
    f_0 = 3.2e-4
    timescale = (f_0)**-1
    velocity_scale = deformation_scale/timescale

    lat_rad_2d = np.deg2rad(lat_in)
    lon_rad_2d = np.deg2rad(lon_in)    

    if northern_hemisphere:
        r_array = (planet_radius * (np.pi/2. - lat_rad_2d))/deformation_scale #non-dim
    else:
        r_array = (planet_radius * (np.pi/2. + lat_rad_2d))/deformation_scale #non-dim        

    v = np.zeros_like(lat_in)
    u = epsilon * ((r_array - r_0)**alpha)* np.exp(-m*((r_array-r_0)**beta))

    v_si_units = v * velocity_scale
    u_si_units = u * velocity_scale

    if northern_hemisphere:
        grad_geopot = ((u_si_units**2)/(r_array* deformation_scale)) + (f_0*np.sin(lat_rad_2d)*u_si_units)
    else:
        #I've changed the sign of the coriolis term here. Clearly this isn't really happening, but in this funny radial coordinate system the sign of u would change for the opposite hemisphere, thus necessitating the sign change.
        grad_geopot = ((u_si_units**2)/(r_array* deformation_scale)) - (f_0*np.sin(lat_rad_2d)*u_si_units)        

    geopotential = np.zeros_like(grad_geopot)

    if northern_hemisphere:
        for lat_idx in range(1, len(lat_rad_2d[:,0])):
            geopotential[lat_idx,:] =  geopotential[lat_idx-1,:] + 0.5*(grad_geopot[lat_idx-1,:]+grad_geopot[lat_idx,:])*(r_array[lat_idx]-r_array[lat_idx-1])
    else:
        r_array_opposite = r_array[::-1,:]
        grad_geopot_opposite = grad_geopot[::-1,:]
        for lat_idx in range(1, len(lat_rad_2d[:,0])):
            geopotential[lat_idx,:] =  geopotential[lat_idx-1,:] + 0.5*(grad_geopot_opposite[lat_idx-1,:]+grad_geopot_opposite[lat_idx,:])*(r_array_opposite[lat_idx]-r_array_opposite[lat_idx-1])    

        geopotential = geopotential[::-1,:]    

    #we want to pass a geopotential field that has an area-mean of zero. This is because we want to preserve the mean geopotential that the model sets as its h_0 parameter.

    delta_lat_arr = np.diff(latb_in, axis=0)[:,0:-1]

    area_array = np.cos(np.deg2rad(lat_in))*np.deg2rad(delta_lat_arr)

    area_av_geopot = np.sum(geopotential*area_array)/np.sum(area_array)

    geopotential_av_removed = geopotential-area_av_geopot

    area_av_final = np.sum(geopotential_av_removed*area_array)/np.sum(area_array)

    print(f'old mean = {area_av_geopot}, final area_av geopot = {area_av_final}')

    geopotential_si_units = geopotential_av_removed * deformation_scale

    h_0 = (deformation_scale*f_0)**2.

    return u_si_units, v_si_units, geopotential_si_units, h_0, grad_geopot

nlat=128
nlon=256

latitudes, latitude_bounds_2  = gg.gaussian_latitudes(int(nlat/2))
latitude_bounds = [latitude_bound[0] for latitude_bound in latitude_bounds_2] + [latitude_bounds_2[-1][1]]

longitudes = np.linspace(0., 360., nlon, endpoint=False)
delta_lon = longitudes[1]-longitudes[0]
longitude_bounds = [lon_val-(0.5*delta_lon) for lon_val in longitudes] + [np.max(longitudes)+(0.5*delta_lon)]
time_arr_adj=None

lon_array_2d, lat_array_2d = np.meshgrid(longitudes, latitudes)
lonb_array_2d, latb_array_2d = np.meshgrid(longitude_bounds, latitude_bounds)

#Note that in the following we're making the initial condition symmetric about the equator. This is because if you only set the initial conditions in the northern hemisphere then you end up needing a very large set of latitudinal functions to get that level of asymmetry, and the code gets very upset when translating that to a finite spectral representation. Making it symmetric gets rid of this problem, at least to some extent.

epsilon = 0.15*2.
alpha = 0.42
beta = 1.3
r_0 = 0.
m_param = 1. 
planet_radius = 55000e3

u_array_vortex, v_array_vortex, height_array_vortex, h_0, grad_geopot_vortex = set_u_v_height_field(lon_array_2d, lat_array_2d,lonb_array_2d, latb_array_2d, epsilon, alpha, beta, m_param, r_0, planet_radius)

u_array_vortex_sp, v_array_vortex_sp, height_array_vortex_sp, h_0_sp, grad_geopot_vortex_sp = set_u_v_height_field(lon_array_2d, lat_array_2d,lonb_array_2d, latb_array_2d, epsilon, alpha, beta, m_param, r_0, planet_radius, northern_hemisphere=False)

epsilon = 0.08
alpha = 0.
beta = 2.
r_0 = 3.37
m_param = 3. 
planet_radius = 55000e3

u_array_jet, v_array_jet, height_array_jet, h_0, grad_geopot_jet = set_u_v_height_field(lon_array_2d, lat_array_2d,lonb_array_2d, latb_array_2d, epsilon, alpha, beta, m_param, r_0, planet_radius)

u_array_jet_sp, v_array_jet_sp, height_array_jet_sp, h_0_sp, grad_geopot_jet_sp = set_u_v_height_field(lon_array_2d, lat_array_2d,lonb_array_2d, latb_array_2d, epsilon, alpha, beta, m_param, r_0, planet_radius, northern_hemisphere=False)


u_array_total = u_array_vortex+u_array_vortex_sp + u_array_jet+u_array_jet_sp
v_array_total = v_array_vortex+v_array_vortex_sp + v_array_jet+v_array_jet_sp
height_array_total = height_array_vortex+height_array_vortex_sp + height_array_jet+height_array_jet_sp
grad_geopot_total = grad_geopot_vortex + grad_geopot_vortex_sp + grad_geopot_jet + grad_geopot_jet_sp

vor_array, div_array = convert_to_vor_div(u_array_total, v_array_total, height_array_total, planet_radius)


p_full=None
p_half=None

npfull=None
nphalf=None

#Output it to a netcdf file. 

file_name='rostami_t85_jet_and_vortex_mk7_gg.nc'


number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlat+1
number_dict['nlonb']=nlon+1 
number_dict['npfull']=npfull
number_dict['nphalf']=nphalf
number_dict['ntime']=None

data_dict = {
    'vor':    vor_array,
    'height': height_array_total,
    'div':    div_array,
    'ucomp': u_array_total, 
    'vcomp': v_array_total,
    'grad_geopot': grad_geopot_total
}


time_units=None

cts.output_multiple_variables_to_file(data_dict,latitudes,longitudes,latitude_bounds,longitude_bounds,p_full,p_half,time_arr_adj,time_units,file_name,number_dict)

print(f'Must set h_0 parameter in code to be {h_0}')



