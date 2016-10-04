# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num
import pdb
import create_timeseries as cts

base_directory='/scratch/sit204/sst_amip_files/'

amip_data_version='amip_data_version_1_1_0' #s 'amip_data_version_1_1_0' or 'amip_data_version_1_0_0'

output_name_list  ={'tosbcs':'sst','siconc':'siconc'}

for variable_name in output_name_list.iterkeys():

	if amip_data_version=='amip_data_version_1_0_0':
		nfiles=50
		folder_name='/1950_1999/'
		filename_prefix='amipbc_sst_360x180_19'
		sst_all=np.zeros((nfiles,12,180,360))
	elif amip_data_version=='amip_data_version_1_1_0':
		nfiles=1
		folder_name=''
		filename_prefix=variable_name+'_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-0_gs1x1_187001-201512'
		do_annual_mean=True

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
		lonbs = resolution_file.variables['bounds_longitude'][:]
		latbs = resolution_file.variables['bounds_latitude'][:]
	except KeyError:
		lonbs = resolution_file.variables['lon_bnds'][:]
		latbs = resolution_file.variables['lat_bnds'][:]

	nlon=lons.shape[0]
	nlat=lats.shape[0]

	nlonb=lonbs.shape[0]
	nlatb=latbs.shape[0]

	lonbs_adjusted=np.zeros(nlonb+1)
	latbs_adjusted=np.zeros(nlatb+1)

	lonbs_adjusted[0:nlonb]=lonbs[:,0]
	lonbs_adjusted[nlonb]=lonbs[-1,1]

	latbs_adjusted[0:nlatb]=latbs[:,0]
	latbs_adjusted[nlatb]=latbs[-1,1]

	day_number = resolution_file.variables['time'][:]

	time_arr = day_number_to_date(day_number, calendar_type = 'gregorian', units_in = 'days since 1870-1-1')
	time_arr_adj=np.arange(15,360,30)

	if len(sst_all.shape)==4:
		sst_in=np.mean(sst_all,axis=0)
	else:
		sst_in=np.zeros((12,180,360))
		for month_tick in np.arange(1,13,1):
			month_idx = np.where(time_arr.month==month_tick)[0]
			sst_in[month_tick-1,:,:]=np.mean(sst_all[month_idx,:,:],axis=0)

	annual_mean_name=''

	if do_annual_mean:
		sst_in_am=np.mean(sst_in,axis=0)
		sst_in=np.zeros((12,180,360))
		for month_tick in np.arange(1,13,1):
			sst_in[month_tick-1,:,:]=sst_in_am
		annual_mean_name='_am'

	p_full=None
	p_half=None

	npfull=None
	nphalf=None

	#Find grid and time numbers

	ntime=time_arr.day.shape[0]
	nlonb=lonbs_adjusted.shape[0]
	nlatb=latbs_adjusted.shape[0]


	#Output it to a netcdf file. 
	file_name=output_name_list[variable_name]+annual_mean_name+'_clim_amip_'+amip_data_version+'.nc'
	variable_name=output_name_list[variable_name]+'_clim_amip'

	number_dict={}
	number_dict['nlat']=nlat
	number_dict['nlon']=nlon
	number_dict['nlatb']=nlatb
	number_dict['nlonb']=nlonb
	number_dict['npfull']=npfull
	number_dict['nphalf']=nphalf
	number_dict['ntime']=ntime

	time_units='days since 0000-01-01 00:00:00.0'

	cts.output_to_file(sst_in,lats,lons,latbs_adjusted,lonbs_adjusted,p_full,p_half,time_arr_adj,time_units,file_name,variable_name,number_dict)


