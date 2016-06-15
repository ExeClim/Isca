#Calculate vertical integrals of area mean, annual mean specific humidity for whole atmosphere and levels above 100hPa (strat wv)
#Could set a threshold to determine spin-up end, e.g. require changes of less than 1% seems plausible. Need more data to confirm this is appropriate however

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from cell_area import cell_area
import pdb
import os

def q_spinup(run_fol, months, plt_dir):

    #personalise
    #model directory
    model_dir = '/scratch/sit204/FMS2013/GFDLmoistModel/'
    #data directory
    data_dir = '/scratch/sit204/Data_2013/'
    #file name
    file_name='atmos_monthly.nc'
    #time-resolution of plotting
    group='months'

    years=int(np.ceil(months/12.))

    #get cell areas and pressure thicknesses
    rundata = xr.open_dataset( data_dir + run_fol + '/run1/'+file_name,
     	        decode_times=False)  # no calendar so tell netcdf lib

    area = cell_area(42, model_dir)
    area_xr = xr.DataArray(area, [('lat', rundata.lat ), ('lon', rundata.lon)])
    dp = xr.DataArray( np.diff(rundata.phalf), [('pfull',rundata.pfull) ])

    files_temp=[data_dir+'/'+run_fol+'/run%d/' % m for m in range(1, months+1)]

    names = [s + file_name for s in files_temp]

        #read data into xarray 
    print 'opening dataset'
    rundata = xr.open_mfdataset( names,
     	        decode_times=False,  # no calendar so tell netcdf lib
	        # choose how data will be broken down into manageable chunks.
	        chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    time_arr = rundata.time

    rundata.coords['months'] = time_arr // 30 + 1
    rundata.coords['years'] = (time_arr // 360) +1

    q_yr = (rundata.sphum).groupby(group).mean(('time'))

    #take area mean of q
    q_av = q_yr*area_xr
    q_avs = q_av.sum(('lat','lon'))/area_xr.sum(('lat','lon'))



    #integrate over pressure levels above 100hPa and over whole atmosphere
#    q_strat = (q_avs[:,0:24]*dp[0:24]*100).sum(('pfull'))/9.8
    min_id=np.min(np.where(q_avs.pfull.to_index() < 100.))
    max_id=np.max(np.where(q_avs.pfull.to_index() < 100.))+1

    q_strat = (q_avs[:,min_id:max_id]*dp[min_id:max_id]*100).sum(('pfull'))/9.8
    q_vint = (q_avs*dp).sum(('pfull'))/9.8

    print 'plotting'

    #plot timeseries of these
    plt.figure(1)
    plt.plot(q_strat,label=run_fol)
    plt.xlabel(group)
    plt.ylabel('Stratospheric area mean specific humidity, kg/m^2')
    plt.legend()
    plt.savefig(plt_dir + '/qstrat_su_num_years_'+str(years)+'_gp_by_'+group+'.png')

    plt.figure(2)
    plt.plot(q_vint,label=run_fol)
    plt.xlabel(group)
    plt.ylabel('Vertically integrated area mean specific humidity, kg/m^2')
    plt.legend()
    plt.savefig(plt_dir + '/qvint_su_num_years_'+str(years)+'_gp_by_'+group+'.png')

    return q_strat, q_vint



if __name__ == "__main__":

        #set run name
        run_fol = 'warmpool_exp_mk1_5'
        plt_dir = '/scratch/sit204/plots/exps/'+run_fol
	if not os.path.exists(plt_dir):
		os.makedirs(plt_dir)
        #number of years to read
	months=40*12
        #return integral of area mean q over stratosphere and whole atmosphere
	q_strat, q_vint = q_spinup(run_fol, months, plt_dir)


        #set run name
        run_fol_2 = 'warmpool_exp_mk2_1'
        plt_dir_2 = '/scratch/sit204/plots/exps/'+run_fol_2
	if not os.path.exists(plt_dir_2):
		os.makedirs(plt_dir_2)

	q_strat_2, q_vint_2 = q_spinup(run_fol_2, months, plt_dir_2)

	plt.figure()

	plt.plot(q_strat,label='old_vert_res')
	plt.plot(q_strat_2,label='new_vert_res')

	plt.legend()
	plt.show()


