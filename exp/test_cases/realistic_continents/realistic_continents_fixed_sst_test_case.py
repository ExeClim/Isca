import os

import numpy as np
from field_table_write import write_ft
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 64

RESOLUTION = "T85", 40
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = IscaCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

#cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase
n_moments =  4

field_table_name = "field_table_age_" + str(n_moments)
write_ft(field_table_name,n_moments)

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('real_T85_mimasm',ext_field_table=field_table_name ,codebase=cb)

#Add any input files that are necessary for a particular experiment.
exp.inputfiles = [os.path.join(base_dir,'input/ssss_ERA5_land_T85.nc'),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
                      os.path.join(base_dir,'input/sst_clim_amip.nc'), os.path.join(base_dir,'input/siconc_clim_amip.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 6, 'hours', time_units='days')


diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True)
diag.add_field('mixed_layer', 'corr_flux', time_avg=True)

diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'wspd', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

diag.add_field('vert_turb', 'z_pbl', time_avg=True)

for ind in range(n_moments):
    name = f"sphum_age_{ind+1}"
    diag.add_field('dynamics', name, time_avg=True)

diag.add_field('atmosphere', 'cape', time_avg=True)
diag.add_field('atmosphere', 'dt_qg_convection', time_avg=True)
diag.add_field('atmosphere', 'dt_qg_condensation', time_avg=True)
diag.add_field('atmosphere', 'dt_sink', time_avg=True)
diag.add_field('atmosphere', 'dt_qg_diffusion', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)

diag.add_field('rrtm_radiation', 'olr', time_avg=True)
diag.add_field('rrtm_radiation', 'toa_sw', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_rad', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_sw_rad', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_lw_rad', time_avg=True)
diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True)
diag.add_field('rrtm_radiation', 'flux_lw', time_avg=True)

exp.diag_table = diag


#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
namelist_name = os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents/namelist_basefile.nml')
nml = f90nml.read(namelist_name)
exp.namelist = nml


def exp_spacing(start, end, num_points, exponent=2.0):
    x = np.linspace(0, 1, num_points)
    return start + (end - start) * x**exponent if num_points > 1 else np.array([start])

def generate_atmospheric_profile(N=41, p_surf=1000):
    bk = np.zeros(N)
    ak = np.zeros(N)

    # Define the layer sizes dynamically
    n_bottom = N // 4
    n_middle = N // 2
    n_top = N - (n_bottom + n_middle)  # Ensure total sum is N
    # Bottom Layer (1000-800)
    ak[:n_bottom] = np.zeros(n_bottom)
    bk[:n_bottom] = np.linspace(1,0.8,n_bottom)#exp_spacing(1, 0.8, n_bottom, 1.7)

    # Middle Layer
    bk[n_bottom:n_bottom + n_middle] = np.linspace(0.76, 0, n_middle)

    # Split middle layer ak allocation to ensure even spread
    mid_half = n_middle // 2 -1
    ak[n_bottom:n_bottom + mid_half] = np.linspace(0,150,mid_half)#exp_spacing(0, 150, mid_half, 1)
    print(mid_half)
    ak[n_bottom + mid_half:n_bottom + n_middle] = exp_spacing(170, 220, n_middle - mid_half, 0.6)

    # Top Layer (200-5)
    ak[-n_top:] = exp_spacing(200, 5, n_top, 0.6)
    bk[-n_top:] = np.zeros(n_top)

    # Compute pressure levels
    pk = ak + bk * p_surf

    # Reverse arrays to match expected order
    bk = bk[::-1]
    ak = ak[::-1]
    pk = pk[::-1]
    
    return ak, bk, pk

# Example usage
n = 41
psurf = 1000
ak, bk, pk = generate_atmospheric_profile(N=n, p_surf=psurf)



exp.update_namelist({
    'mixed_layer_nml': {  
        'do_qflux' : False, #Don't use the prescribed analytical formula for q-fluxes
        'do_read_sst' : True, #Read in sst values from input file
        'do_sc_sst' : True, #Do specified ssts (need both to be true)
        'sst_file' : 'sst_clim_amip', #Set name of sst input file
        'specify_sst_over_ocean_only' : True, #Make sure sst only specified in regions of ocean.
    }})
"""    ,
    'vert_coordinate_nml': {
        'bk': bk,
        'pk': 100*ak,
       }
})"""

exp.set_resolution(*RESOLUTION)
#Lets do a run!
if __name__=="__main__":
    """    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2,241):
        exp.run(i, num_cores=NCORES)"""
        
    # Start from restart file and run for 2 weeks (50 years = 624 weeks) 
    #res_file = "/home/philbou/projects/def-rfajber/philbou/restarts/realistic_continents_T85/res_1200_reallistic_continents_T85.tar.gz"
    #res_file = "/home/philbou/projects/def-rfajber/philbou/restarts/realistic_continents_T85/restart_T85.tar.gz"
    res_file = "/home/philbou/scratch/isca_data/realistic_continents_T85_year/restarts/res0030.tar.gz"
    res_file = "/home/philbou/scratch/isca_data/realistic_conts_T85_sp42/restarts/res0120.tar.gz"
    exp.run(1,restart_file = res_file,use_restart=False, num_cores=NCORES,overwrite_data = True)
    for i in range(2,121):
        exp.run(i, num_cores=NCORES,overwrite_data = True)
        