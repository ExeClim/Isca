#Function to load interpolated data on pressure levels for multiple months
#Copied from RG on 29/03/16

import subprocess

def plevel_call(nc_file_in,nc_file_out, var_names = '-a', p_levels='default'):
#    for mn in mn_range:
        interper = './plevel.sh'
        nc_file = ' -i '+nc_file_in
        out_file = ' -o '+nc_file_out
        if p_levels == 'model':
            plev = ' -p "2 9 18 38 71 125 206 319 471 665 904 1193 1532 1925 2375 2886 3464 4115 4850 5679 6615 7675 8877 10244 11801 13577 15607 17928 20585 23630 27119 31121 35711 40976 47016 53946 61898 71022 81491 93503" '
            command = interper + nc_file + out_file + plev + var_names
	elif p_levels=='default':
            command = interper + nc_file + out_file + ' ' + var_names
        else:
	    plev=p_levels
	    command = interper + nc_file + out_file + plev +' '+ var_names
        print command
        subprocess.call([command], shell=True)

