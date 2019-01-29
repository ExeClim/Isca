import numpy as np
import matplotlib.pyplot as plt
import texttable as tt
import pdb

planet_to_check = 'earth'

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
    target_number_sols = 674.7707745374726
    target_length_solar_day =  1377648.0


n_days = np.shape(n_day_range)[0]
n_len_day = np.shape(len_day_range)[0]

n_day_arr = np.zeros((n_days, n_len_day))
n_len_day_arr = np.zeros((n_days, n_len_day))

for al in range(n_len_day):
    for ts in range(n_days):
        n_day_arr[ts, al] = n_day_range[ts]
        n_len_day_arr[ts,al] = len_day_range[al]

integer_values = n_day_arr * n_len_day_arr / (n_day_arr -1)
res_arr = np.mod(integer_values,1)
where_integer = np.where(res_arr==0.0)

cs = plt.contourf(n_len_day_arr, n_day_arr, integer_values)
fig = plt.gcf()

cb = fig.colorbar(cs)

table_to_print = tt.Texttable()

headings = ['i,   j', 'nsol', 'nsid', 'lsol', 'lsid', 'orbital_period']
table_to_print.header(headings)

names = ["{0:3d}, {1:3d}".format(where_integer[0][i], where_integer[1][i]) for i in range(len(where_integer[0]))]

nsol = ["{0:6.2f}".format(-1+n_day_arr[where_integer[0][i], where_integer[1][i]]) for i in range(len(where_integer[0]))]
nsid = ["{0:6.2f}".format(n_day_arr[where_integer[0][i], where_integer[1][i]]) for i in range(len(where_integer[0]))]
lsid = ["{0:6.2f}".format(n_len_day_arr[where_integer[0][i], where_integer[1][i]]) for i in range(len(where_integer[0]))]
lsol = ["{0:6.2f}".format(integer_values[where_integer[0][i], where_integer[1][i]]) for i in range(len(where_integer[0]))]
orbital_period_1  = ["{0:6.2f}".format(float(nsol[i])*float(lsol[i])) for i in range(len(nsol))]
orbital_period_2  = ["{0:6.2f}".format(float(nsid[i])*float(lsid[i])) for i in range(len(nsid))]

if orbital_period_1 != orbital_period_2:
    raise ValueError('orbital period from sidereal day not equal to orbital period from solar day')

for row in zip(names, nsol, nsid, lsol, lsid, orbital_period_1):
    table_to_print.add_row(row)

s = table_to_print.draw()
print('Table for '+planet_to_check, ' with target number of sols = {0:6.3f} and length of solar day = {1:9.2f}'.format(target_number_sols, target_length_solar_day))
print(s)