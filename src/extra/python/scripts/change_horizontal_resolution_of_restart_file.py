"""Script for changing the horizontal resolution of an FMS restart file"""
import xarray as xar
import numpy as np
import gauss_grid as gg
import scipy.interpolate as scinterp
import pdb
import mpl_toolkits.basemap as basemap
import matplotlib.pyplot as plt
import sh
from netCDF4 import Dataset
import copy_netcdf_attrs as cna
import tempfile
import shutil
import os

def linear_interpolate_for_regrid(lon_list_in_grid, lat_list_in_grid, lon_list_out_grid, lat_list_out_grid, input_array):

    time_length, z_length, lat_length, lon_length = input_array.shape
    
    lon_array, lat_array = np.meshgrid(lon_list_out_grid,lat_list_out_grid)


    output_array = np.zeros((time_length, z_length, lat_list_out_grid.shape[0], lon_list_out_grid.shape[0]))

    for tim in range(time_length):
        for z in range(z_length):
            input_array_2d = np.squeeze(input_array[tim,z,...])            
            output_array[tim,z,...] = basemap.interp(input_array_2d, lon_list_in_grid, lat_list_in_grid, lon_array, lat_array, order=1)
            
    return output_array

def populate_new_spherical_harmonic_field(x_in, y_in, x_out, y_out, input_array):

    time_length, z_length, y_length, x_length = input_array.shape

    output_array = np.zeros((time_length, z_length, y_out.shape[0], x_out.shape[0]))

    output_array[:,:,0:y_length, 0:x_length] = input_array
    
    return output_array


def process_input_file(file_name, atmosphere_or_spectral_dynamics, num_fourier_out, num_x_out, num_y_out):

    num_spherical_out = num_fourier_out + 1

    dataset = xar.open_dataset(file_name+'.res.nc')

    Time_in = dataset.Time.values #Just 1 and 2
    Time_out = Time_in
        
    dataset_out = dataset.copy(deep=True)

    if atmosphere_or_spectral_dynamics=='atmosphere':
        x_axis_1_in = dataset.xaxis_1.values #Just a single value - 1.0
        x_axis_2_in = dataset.xaxis_2.values #Bizarrely, this is the number of vertical levels+1

        y_axis_1_in = dataset.yaxis_1.values #Just a single value - 1.0
        y_axis_2_in = dataset.yaxis_2.values #Number of spherical harmonics + 1 (87)

        z_axis_1_in = dataset.zaxis_1.values #Just a single value - 1.0
        z_axis_2_in = dataset.zaxis_2.values #Number of full levels (30)

        longitudes_in  = np.arange(0., 360., (360./x_axis_2_in.shape[0]))
        latitudes_in  = gg.gaussian_latitudes(int(y_axis_2_in.shape[0]/2))[0]        
        
        axes_out = {'xaxis_1':x_axis_1_in, 'xaxis_2':np.arange(1.,num_x_out+1), 'yaxis_1':y_axis_1_in, 'yaxis_2':np.arange(1.,num_y_out+1), 'zaxis_1':dataset.zaxis_1.values, 'zaxis_2':dataset.zaxis_2.values, 'Time':Time_in}        

        y_axis_name = 'yaxis_2'
        x_axis_name = 'xaxis_2'
        
    elif atmosphere_or_spectral_dynamics=='spectral_dynamics':    
        x_axis_1_in = dataset.xaxis_1.values #Just a single value - 1.0
        x_axis_2_in = dataset.xaxis_2.values #Bizarrely, this is the number of vertical levels+1
        x_axis_3_in = dataset.xaxis_3.values #Number of spherical harmonics (e.g. 86)
        x_axis_4_in = dataset.xaxis_4.values #Number of gridpoints in physical space

        y_axis_1_in = dataset.yaxis_1.values #Just a single value - 1.0
        y_axis_2_in = dataset.yaxis_2.values #Number of spherical harmonics + 1 (87)
        y_axis_3_in = dataset.yaxis_3.values #Number of gridpoints in physical space

        z_axis_1_in = dataset.zaxis_1.values #Just a single value - 1.0
        z_axis_2_in = dataset.zaxis_2.values #Number of full levels (30)

        longitudes_in  = np.arange(0., 360., (360./x_axis_4_in.shape[0]))
        latitudes_in  = gg.gaussian_latitudes(int(y_axis_3_in.shape[0]/2))[0]

        axes_out = {'xaxis_1':x_axis_1_in, 'xaxis_2':x_axis_2_in, 'xaxis_3':np.arange(1.,num_spherical_out+1), 'xaxis_4':np.arange(1.,num_x_out+1), 'yaxis_1':y_axis_1_in, 'yaxis_2':np.arange(1.,num_spherical_out+2), 'yaxis_3':np.arange(1.,num_y_out+1), 'zaxis_1':dataset.zaxis_1.values, 'zaxis_2':dataset.zaxis_2.values, 'Time':Time_in}

        y_axis_name = 'yaxis_3'
        x_axis_name = 'xaxis_4'

    longitudes_out = np.arange(0., 360., (360./num_x_out))
    latitudes_out = gg.gaussian_latitudes(int(num_y_out/2))[0]    

    for var in list(dataset_out.data_vars.keys()):
        dataset_out.__delitem__(var)

    for coord in list(dataset_out.coords.keys()):
        dataset_out[coord] = axes_out[coord]
        dataset_out[coord].attrs = dataset[coord].attrs

    for var in list(dataset.data_vars.keys()):

        var_dims = dataset[var].dims

        if var_dims[2:4] == (y_axis_name, x_axis_name):
            print((var, 'physical grid'))
            new_var = linear_interpolate_for_regrid(longitudes_in, latitudes_in, longitudes_out, latitudes_out, dataset[var].load().values)
            dataset_out[var] = (dataset[var].dims, new_var)
    
        elif var_dims[2:4] == ('yaxis_2', 'xaxis_3'):    
            print((var, 'spectral grid'))
            new_var = populate_new_spherical_harmonic_field(x_axis_2_in, y_axis_2_in, axes_out['xaxis_3'], axes_out['yaxis_2'], dataset[var].values)
            dataset_out[var] = ((dataset[var].dims, new_var))
        else:
            print((var, 'neither'))
            dataset_out[var] = ((dataset[var].dims, dataset[var].values))
            
        dataset_out[var].attrs = dataset[var].attrs

    out_file_name = file_name+'_mod_'+str(num_fourier_out)+'_onescript.res.nc'
    dataset_out.to_netcdf(path='./temp_'+out_file_name, format='NETCDF3_CLASSIC', engine='scipy')
    
    remove_fill_value_attribute('./temp_'+out_file_name, out_file_name)
    
    os.remove('./temp_'+out_file_name)
    
    return out_file_name

def join_into_cpio(atmosphere_file_name='./atmosphere.res.nc', spectral_dynamics_file_name='./spectral_dynamics.res.nc', atmos_model_file_name='./atmos_model.res', restart_file_out_name='./res_mod'):

    temp_folder_name = tempfile.mkdtemp()    
    
    shutil.move(atmosphere_file_name, temp_folder_name+'/atmosphere.res.nc')
    shutil.move(spectral_dynamics_file_name, temp_folder_name+'/spectral_dynamics.res.nc')    
    shutil.copyfile(atmos_model_file_name, temp_folder_name+'/atmos_model.res')

    state_files_out = ['atmosphere.res.nc', 'spectral_dynamics.res.nc', 'atmos_model.res']

    cwd = os.getcwd()
    
    os.chdir(temp_folder_name) #s have to move to temporary folder as cpio cannot cope with absolute file references, as otherwise when you delete the temporary folder, cpio will go looking for the temporary folder when it's extracted.
    
    sh.cpio('-ov', _in='\n'.join(state_files_out), _out=restart_file_out_name)

    shutil.move(restart_file_out_name, cwd)    

    os.chdir(cwd)

    shutil.move(temp_folder_name+'/atmosphere.res.nc', atmosphere_file_name)
    shutil.move(temp_folder_name+'/spectral_dynamics.res.nc', spectral_dynamics_file_name)    

    shutil.rmtree(temp_folder_name)
    
def remove_fill_value_attribute(in_file_name, out_file_name):

    dsin = Dataset(in_file_name,  'a', format='NETCDF3_CLASSIC')
    dsout= Dataset(out_file_name, 'w', format='NETCDF3_CLASSIC')    
    cna.copy_netcdf_attrs(dsin, dsout, copy_vars = True)
    dsout.close()
    
if __name__=="__main__":

    #Specify the number of fourier modes and lon and lat dimensions for the output
    num_fourier_out = 85
    num_x_out = 256
    num_y_out = 128

    #Specify the name of the input files that you want to regrid
    atmosphere_file_name = 'atmosphere_old'
    spectral_dynamics_file_name = 'spectral_dynamics_old'
    atmos_model_file_name = 'atmos_model.res'
    
    #Specify the name of the output cpio archive    
    restart_file_out_name = 'res_85_onescript'
    
    #Regridding atmosphere file
    atmosphere_out_file_name = process_input_file(atmosphere_file_name,        'atmosphere',        num_fourier_out, num_x_out, num_y_out)
    #regridding spectral dynamics file
    spectral_out_file_name   = process_input_file(spectral_dynamics_file_name, 'spectral_dynamics', num_fourier_out, num_x_out, num_y_out)
    
    #merging into a single archive
    join_into_cpio(atmosphere_out_file_name, spectral_out_file_name, atmos_model_file_name, restart_file_out_name=restart_file_out_name)
    
