import numpy as np
import matplotlib.pyplot as plt
import pdb

n_day_range = np.arange(660, 680, 1)
len_day_range = np.arange(88000, 89000, 1)

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
#cn = plt.contour(albedo_range, tau_s_range, t_surf_arr,levels=[210., 260.], colors=['k'])
#plt.clabel(cn, inline=1, fontsize=10.)
fig = plt.gcf()

cb = fig.colorbar(cs)


