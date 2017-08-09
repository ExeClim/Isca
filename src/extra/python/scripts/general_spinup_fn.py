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

def q_spinup(run_fol, var_to_integrate, start_month, end_month, plt_dir):

    #personalise
    #model directory
    model_dir = '/scratch/sit204/FMS2013/GFDLmoistModel/'
    #data directory
    data_dir = '/scratch/sit204/Data_2013/'
    #file name
    file_name='atmos_monthly.nc'
    #time-resolution of plotting
    group='months'
    scaling=1.
    t_resolution=42
    nlon=128
    nlat=64
    gravity=9.8

    years=int(np.ceil((end_month-start_month)/12.))

    #get cell areas and pressure thicknesses

    try:
        files_temp=[data_dir+'/'+run_fol+'/run%03d/' % m for m in range(start_month, end_month+1)]

        names = [s + file_name for s in files_temp]
        print names[0]
        rundata = xr.open_dataset(names[0],
                     decode_times=False)  # no calendar so tell netcdf lib
    except RuntimeError:
        files_temp=[data_dir+'/'+run_fol+'/run%d/' % m for m in range(start_month, end_month+1)]

        names = [s + file_name for s in files_temp]

        rundata = xr.open_dataset(names[0],
                     decode_times=False)  # no calendar so tell netcdf lib

    area = cell_area(t_resolution, model_dir)
    area_xr = xr.DataArray(area, [('lat', rundata.lat ), ('lon', rundata.lon)])
    dp = xr.DataArray( np.diff(rundata.phalf), [('pfull',rundata.pfull) ])

        #read data into xarray 
    print 'opening dataset'
    rundata = xr.open_mfdataset( names,
                 decode_times=False,  # no calendar so tell netcdf lib
            # choose how data will be broken down into manageable chunks.
            chunks={'time': 30, 'lon': nlon//4, 'lat': nlat//2})

    time_arr = rundata.time

    rundata.coords['months'] = time_arr // 30 + 1
    rundata.coords['years'] = ((time_arr // 360) +1)*scaling - 6.

    q_yr = (rundata[var_to_integrate]).groupby(group).mean(('time'))

    #take area mean of q
    q_av = q_yr*area_xr
    q_avs = q_av.sum(('lat','lon'))/area_xr.sum(('lat','lon'))

    try:
        q_avs.pfull
    except AttributeError:
        print 'data is 2d'
        q_vint=q_avs
        q_vint.load()
        q_strat=q_vint
    else:
        #integrate over pressure levels above 100hPa and over whole atmosphere
    #    q_strat = (q_avs[:,0:24]*dp[0:24]*100).sum(('pfull'))/gravity
        min_id=np.min(np.where(q_avs.pfull.to_index() < 100.))
        max_id=np.max(np.where(q_avs.pfull.to_index() < 100.))+1

        q_strat = (q_avs[:,min_id:max_id]*dp[min_id:max_id]*100).sum(('pfull'))/gravity
        q_vint = (q_avs*dp).sum(('pfull'))/gravity

        q_strat.load()
        q_vint.load()

    time_arr=q_vint.months.values

    rundata.close()

    return q_strat, q_vint, time_arr



if __name__ == "__main__":

    start_month_offset=[0, 0]

#    exp_list=['without_dry_convection_bug/giant_drag_exp_chai_values_1', 'giant_drag_exp_chai_values_1',
#          'without_dry_convection_bug/giant_drag_exp_chai_values_1_bar_damping_1', 'giant_drag_exp_chai_values_1_bar_damping_1']

#    label_arr=['without bug chai 1','with bug chai 1',
#           'without bug chai 1 bar 1','with bug chai 1 bar 1']

    exp_list=['no_ice_flux_lhe_exps_fixed_sst_1', 'no_ice_flux_q_exps_qflux_control_1']
    
    label_arr=['fixed DJF sst no ice', 'q-flux control DJF no ice']
           
    exp_name=exp_list

    #number of years to read
    start_month_arr=[1, 1]
    end_month_arr=[360, 163]

    len_list=[len(start_month_offset), len(exp_list), len(label_arr), len(start_month_arr), len(end_month_arr)]

    if not all(x==len_list[0] for x in len_list):
        raise IndexError, "Input arrays to routine are not all the same length"


    variable_to_integrate='temp'

    plt.figure()

    for exp_number in exp_list:
        #set run name
        run_fol = str(exp_number)
        plt_dir = '/scratch/sit204/plots/exps/'+run_fol
        if not os.path.exists(plt_dir):
            os.makedirs(plt_dir)
        idx=exp_list.index(exp_number)
        print 'running ', exp_number
        #return integral of area mean q over stratosphere and whole atmosphere
        q_strat, q_vint, time = q_spinup(run_fol, variable_to_integrate, start_month_arr[idx]+start_month_offset[idx], end_month_arr[idx]+start_month_offset[idx], plt_dir)
        plt.plot(time,q_vint,label=label_arr[idx])
#         plt.plot(time,q_strat,label='strat '+label_arr[idx])
        

    plt.xlabel('time (months)')
    plt.ylabel('Global average '+variable_to_integrate)
    plt.legend(loc='upper left')
    plt.show()
