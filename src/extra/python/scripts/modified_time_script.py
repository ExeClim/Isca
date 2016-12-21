import os
import numpy as np
import matplotlib.pyplot as plt
import pdb
from datetime import datetime


def calculate_month_run_time(exp_dir_list):
	"""A script that takes a list of experiment names as input, and plots the time taken to run each month in that experiment vs the wall time. """

	try:
	    GFDL_DATA        = os.environ['GFDL_DATA']
	except Exception, e:
	    print('Environment variables GFDL_DATA must be set')
	    exit(0)


	for exp_dir in exp_dir_list:

		exp_dir_full = GFDL_DATA+'/'+exp_dir+'/'

		#Finds all the months for particular experiment
		months_to_check=os.listdir(exp_dir_full) 

		delta_t_arr=np.zeros(len(months_to_check)-1)
		end_t_arr=[]

		for month in np.arange(len(months_to_check)-1)+1:
			#Calculates the time between the current month's folder being modified, and the previous month's folder being modified.
			delta_t = os.path.getctime(exp_dir_full+months_to_check[month])-os.path.getctime(exp_dir_full+months_to_check[month-1])
			#Converts this time to minutes from seconds:
			delta_t_arr[month-1]=delta_t/60.
			#Saves time of modification as python datetime object:
			end_t_arr.append(datetime.fromtimestamp(os.path.getctime(exp_dir_full+months_to_check[month])))

		#Plots results for particular experiment
		plt.plot(end_t_arr,delta_t_arr, label=exp_dir)

	plt.legend()
	plt.xlabel('Wall time (GMT)')
	plt.ylabel('Wall time elapsed per month (minutes)')

if __name__=="__main__":

	exp_dir_list = ['simple_continents_post_princeton_qflux_anoms_1']

	calculate_month_run_time(exp_dir_list)

	plt.show()

