import os
import numpy as np
import matplotlib.pyplot as plt
import pdb
import sys
from datetime import datetime


def calculate_month_run_time(exp_dir_list, plot_against_wall_time=True, file_to_use_for_timing = 'logfile.000000.out'):
    """A script that takes a list of experiment names as input, and plots the time taken to run each month in that experiment vs the wall time. """

    try:
        GFDL_DATA        = os.environ['GFDL_DATA']
    except Exception as e:
        print('Environment variables GFDL_DATA must be set')
        sys.exit(0)


    for exp_dir in exp_dir_list:

        exp_dir_full = GFDL_DATA+'/'+exp_dir+'/'

        #Finds all the months for particular experiment
        months_to_check=os.listdir(exp_dir_full)
        try:
            months_to_check.remove('restarts')
        except:
            pass 
        months_to_check.sort()
        
        try:
            months_to_check.remove('.DS_Store')
            months_to_check.remove('._.DS_Store')
        except ValueError:
            pass

        delta_t_arr=np.zeros(len(months_to_check)-1)
        end_t_arr=[]

        for month in np.arange(len(months_to_check)-1)+1:
            #Calculates the time between the current month's folder being modified, and the previous month's folder being modified.
            delta_t = os.path.getctime(exp_dir_full+months_to_check[month]+'/'+file_to_use_for_timing)-os.path.getctime(exp_dir_full+months_to_check[month-1]+'/'+file_to_use_for_timing)
                       
            #Converts this time to minutes from seconds:
            delta_t_arr[month-1]=delta_t/60.
            #Saves time of modification as python datetime object:
            end_t_arr.append(datetime.fromtimestamp(os.path.getmtime(exp_dir_full+months_to_check[month])))

        month_num_arr = [int(months_to_check[num].replace('run', '')) for num in range(len(months_to_check))]

        months_idx_to_remove = [num for num in np.where(np.abs(delta_t_arr) > 10.*np.mean(delta_t_arr))[0]]

        print(('removing anomalously long delta_t for months ', [month_num_arr[month] for month in months_idx_to_remove], [delta_t_arr[month] for month in months_idx_to_remove]))
        delta_t_arr[np.where(np.abs(delta_t_arr) > 10.*np.mean(delta_t_arr))] = np.nan

        #Plots results for particular experiment
        if plot_against_wall_time:
            plt.plot(end_t_arr,delta_t_arr, label=exp_dir)   
            plt.xlabel('Wall time (GMT)')                 
        else:
            plt.plot(month_num_arr[:-1], delta_t_arr, label=exp_dir)                
            plt.xlabel('Month number')
        
    plt.legend()
    plt.ylabel('Wall time elapsed per month (minutes)')

if __name__=="__main__":

    exp_dir_list = ['bog_fixed_sst_low_bog_a_ocean_topog_85']

    calculate_month_run_time(exp_dir_list, plot_against_wall_time=False, file_to_use_for_timing='git_hash_used.txt')

    plt.show()

