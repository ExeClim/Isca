import xarray as xar
import numpy as np
import pdb
import matplotlib.pyplot as plt
import create_timeseries as cts


def adjust_data(input_data):
	""" 
	adjust_data solves the matrix problem Ax=b, where A is a matrix with recurring entries of 0.125, 0.75, 0.125,
	x is 12 unknowns, and b is a column vector of 12 monthly average data points. The purpose of solving for x is
	that the values for x will produce mid-month values that, when linearly interpolated by a model's interpolator
	will reproduce monthly means equal to b. This method is described in detail here:
	http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amip2bcs.php.
	"""

	multiply_matrix=np.asarray([
	[0.75, 0.125, 0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.125],
	[0.125,0.75,  0.125,0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.   ],
	[0.,   0.125, 0.75, 0.125,0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.   ],
	[0.,   0.,    0.125,0.75, 0.125,0.,   0.,   0.,   0.,   0.,   0.,   0.   ],
	[0.,   0.,    0.,   0.125,0.75, 0.125,0.,   0.,   0.,   0.,   0.,   0.   ],
	[0.,   0.,    0.,   0.,   0.125,0.75, 0.125,0.,   0.,   0.,   0.,   0.   ],
	[0.,   0.,    0.,   0.,   0.,   0.125,0.75, 0.125,0.,   0.,   0.,   0.   ],
	[0.,   0.,    0.,   0.,   0.,   0.,   0.125,0.75, 0.125,0.,   0.,   0.   ],
	[0.,   0.,    0.,   0.,   0.,   0.,   0.,   0.125,0.75, 0.125,0.,   0.   ],
	[0.,   0.,    0.,   0.,   0.,   0.,   0.,   0.,   0.125,0.75, 0.125,0.   ],
	[0.,   0.,    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.125,0.75, 0.125],
	[0.125,0.,    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.125,0.75 ]
	])

	real_monthly_means = input_data

	solution = np.linalg.solve(multiply_matrix, real_monthly_means)

	return solution

def perform_adj(data_array):
	"""
	Function for looping over spatial points in data array.
	"""

	output_array=np.zeros_like(data_array)

	nt, ny, nx = np.shape(data_array)

	for j in range(ny):
		print(j)

		for i in range(nx):
			output_array[:,j,i] = adjust_data(data_array[:,j,i])

	return output_array

def output_to_file(dataset, output_array, output_filename, variable_name):
	"""
	Simple function for outputting adjusted data as a netcdf file.
	"""
	nt, ny, nx = np.shape(output_array)

	p_full=None
	p_half=None

	npfull=None
	nphalf=None

	#Find grid and time numbers

	lats=dataset.lat.values
	lons=dataset.lon.values
	latbs=dataset.latb.values
	lonbs=dataset.lonb.values
	time_arr=dataset.time.values

	ntime=nt
	nlonb=dataset.lonb.values.shape[0]
	nlatb=dataset.latb.values.shape[0]


	#Output it to a netcdf file. 


	number_dict={}
	number_dict['nlat']=ny
	number_dict['nlon']=nx
	number_dict['nlatb']=nlatb
	number_dict['nlonb']=nlonb
	number_dict['npfull']=npfull
	number_dict['nphalf']=nphalf
	number_dict['ntime']=ntime

	time_units='days since 0000-01-01 00:00:00.0'

	cts.output_to_file(output_array,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,output_file_name,variable_name,number_dict)

def run_steps(input_file_name, input_variable_name, output_file_name, output_variable_name):
	"""
	Function for performing the necessary steps to execute the program.
	"""
	
	dataset=xar.open_dataset(input_file_name, decode_times=False)

	data_array=dataset[input_variable_name].values
	
	adjusted_data_array = perform_adj(data_array)

	output_to_file(dataset, adjusted_data_array, output_file_name, output_variable_name)


if __name__=="__main__":

	input_file_name='/scratch/sit204/FMS2013/GFDLmoistModel/src/extra/python/scripts/nc_files_11_01_17/ami_qflux_ctrl_ice_4320.nc'
	input_variable_name='ami_qflux_ctrl_ice_4320'

	output_file_name='ami_qflux_ctrl_ice_4320m.nc'
	output_variable_name='ami_qflux_ctrl_ice_4320m'

	run_steps(input_file_name, input_variable_name, output_file_name, output_variable_name)


