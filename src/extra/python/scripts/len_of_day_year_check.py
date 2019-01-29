import numpy as np
import matplotlib.pyplot as plt
import texttable as tt
import pdb

""" The idea of this script is to try to find orbital periods and rotation rates that avoid leap years in planetary calendars. Like the way that we have 360-day calendars for Earth."""

planet_to_check = 'titan'

if planet_to_check=='mars':
    n_day_range = np.arange(660, 680, 1)
    len_day_range = np.arange(88000, 89000, 1)
    target_number_sols = 668.5991
    target_length_solar_day =  (86400.)+(39.*60.)+35.
elif planet_to_check =='earth':
    n_day_range = np.arange(350, 370, 1)
    len_day_range = np.arange(86000, 88000, 1)    
    target_number_sols = 365.25
    target_length_solar_day =  86400    
elif planet_to_check =='titan':
    n_day_range = np.arange(655, 675, 1)
    len_day_range = np.arange(1377550, 1377750, 1)       

    target_length_sidereal_day = 15.945*86400.
    target_length_orbital_period = 10759.22*86400.

    target_length_solar_day = target_length_sidereal_day / (1.-target_length_sidereal_day/target_length_orbital_period)
    target_number_sols = target_length_orbital_period / target_length_solar_day


n_days = np.shape(n_day_range)[0]
n_len_day = np.shape(len_day_range)[0]

n_day_arr = np.zeros((n_days, n_len_day))
n_len_day_arr = np.zeros((n_days, n_len_day))

for al in range(n_len_day):
    for ts in range(n_days):
        n_day_arr[ts, al] = n_day_range[ts]
        n_len_day_arr[ts,al] = len_day_range[al]

integer_values = n_day_arr * n_len_day_arr / (n_day_arr -1) #equation here is essentially number of sols * length of sol = number of sidereal days * length of sidereal day. Here integer_values = length_of_sol, n_day_arr is number of sidereal days, n_len_day_arr is length of sidereal day, and n_day_arr - 1 = number of sols.

res_arr = np.mod(integer_values,1) #we are checking here that the length of the solar day is an integer number of seconds.
where_integer = np.where(res_arr==0.0)

cs = plt.contourf(n_len_day_arr, n_day_arr, integer_values)
fig = plt.gcf()

cb = fig.colorbar(cs)

table_to_print = tt.Texttable()

headings = ['i,   j', 'nsol', 'nsid', 'lsol', 'lsid', 'orbital_period']
table_to_print.header(headings)

names = ["{0:3d}, {1:3d}".format(where_integer[0][i], where_integer[1][i]) for i in range(len(where_integer[0]))]

nsol = ["{0:6.2f}".format(-1+n_day_arr[where_integer[0][i], where_integer[1][i]]) for i in range(len(where_integer[0]))] #n_day_arr - 1 = number of sols.
nsid = ["{0:6.2f}".format(n_day_arr[where_integer[0][i], where_integer[1][i]]) for i in range(len(where_integer[0]))] #n_day_arr is number of sidereal days
lsid = ["{0:6.2f}".format(n_len_day_arr[where_integer[0][i], where_integer[1][i]]) for i in range(len(where_integer[0]))] #n_len_day_arr is length of sidereal day
lsol = ["{0:6.2f}".format(integer_values[where_integer[0][i], where_integer[1][i]]) for i in range(len(where_integer[0]))] #integer_values = length_of_sol
orbital_period_1  = ["{0:6.2f}".format(float(nsol[i])*float(lsol[i])) for i in range(len(nsol))] #Checking orbital period = number_of_sols * length_of_sol
orbital_period_2  = ["{0:6.2f}".format(float(nsid[i])*float(lsid[i])) for i in range(len(nsid))] #checking orbital period = number of sidereal days * length of sidereal day

if orbital_period_1 != orbital_period_2:
    raise ValueError('orbital period from sidereal day not equal to orbital period from solar day')

for row in zip(names, nsol, nsid, lsol, lsid, orbital_period_1):
    table_to_print.add_row(row)

s = table_to_print.draw()
print('Table for '+planet_to_check, ' with target number of sols = {0:6.3f} and length of solar day = {1:9.2f}'.format(target_number_sols, target_length_solar_day))
print(s)
#For Isca, set the following:

#diag.add_file('atmos_daily',length_of_solar_day , 'seconds', time_units='days')

    # 'main_nml': {
    #     'dt_atmos': whatever_you_like,
    #     'days': 0.,
    #     'seconds': number_of_solar_days_to_put_in_one_file*length_of_solar_day,
    #     'calendar': 'no_calendar'
    # },

    # 'constants_nml': {
    #     'orbital_period': orbital_period_from_this_script,
    #     'rotation_period':length_of_SIDEREAL_day,
    # },

# By doing this, the length of a year will be an integer number of solar and sidereal days.