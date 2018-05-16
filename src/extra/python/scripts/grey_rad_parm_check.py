import numpy as np
import matplotlib.pyplot as plt
import pdb

sigma=5.67e-8
solar_constant = 589.0


tau_s_range = np.arange(0.01, 2.0, 0.01)
albedo_range = np.arange(0.1, 0.8, 0.01)

n_tau_s = np.shape(tau_s_range)[0]
n_albedo = np.shape(albedo_range)[0]

tau_s_arr = np.zeros((n_tau_s, n_albedo))
albedo_arr = np.zeros((n_tau_s, n_albedo))

for al in range(n_albedo):
    for ts in range(n_tau_s):
        tau_s_arr[ts, al] = tau_s_range[ts]
        albedo_arr[ts,al] = albedo_range[al]

def calc_t_surf(albedo_in, tau_s_in):
    olr = (1.-albedo_in)* solar_constant

    t_surf = (olr*(tau_s_in+1.)/(2.*sigma))**0.25

    return t_surf

t_surf_arr = calc_t_surf(albedo_arr, tau_s_arr)

cs = plt.contourf(albedo_range, tau_s_range, t_surf_arr)
cn = plt.contour(albedo_range, tau_s_range, t_surf_arr,levels=[210., 260.], colors=['k'])
plt.clabel(cn, inline=1, fontsize=10.)
fig = plt.gcf()

cb = fig.colorbar(cs)

t_surf_value = calc_t_surf(0.25, 0.2)
t_surf_shiny = calc_t_surf(0.7, 0.2)

