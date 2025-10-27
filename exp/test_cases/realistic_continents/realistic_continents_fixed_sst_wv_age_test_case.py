import os

import numpy as np

import f90nml

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


# functions to modify field table to add new passive tracers
def get_field(n):
    base = {
            "atmos_mod": "sphum_age",
            "longname": "sphum times",
            "units": "sec (kg/kg)",
            "numerical_representation": "grid",
            "hole_filling": "off",
            "advect_vert": "finite_volume_parabolic",
            "robert_filter": "on",
            "profile_type": ["fixed", "surface_value=0.0"]
        }
    base["atmos_mod"] += f"_{n+1}"
    base["longname"] += f" {n+1}-th moment "
    return base

def write_ft(name,n_moments):
    add_dir = "/src/extra/model/isca/"
    field_table_dir = GFDL_BASE + add_dir

    sphum ={
            "atmos_mod": "sphum",
            "longname": "specific humidity",
            "units": "kg/kg",
            "numerical_representation": "grid",
            "hole_filling": "off",
            "advect_vert": "finite_volume_parabolic",
            "robert_filter": "on",
            "profile_type": ["fixed", "surface_value=0.0"]
            }
        
    
    # Define the file name without an extension
    output_file = field_table_dir + name

    # Open the file in write mode
    with open(output_file, 'w') as file:
        # Write sphum tracer
        file.write(f'"TRACER",')
        for key, value in sphum.items():
            if type(value) == str:
                file.write(f'"{key}", "{value}"\n')
            elif type(value) == list:
                file.write(f'"{key}", "{value[0]}", "{value[1]}"\n')
        file.write('/ \n')

         # write WV age distribution moments tracers
        for ind in range(n_moments):
            field = get_field(ind)
            file.write(f'"TRACER",')
            for key, value in field.items():
                if type(value) == str:
                    file.write(f'"{key}", "{value}"\n')
                elif type(value) == list:
                    file.write(f'"{key}", "{value[0]}", "{value[1]}"\n')
            file.write('/ \n')




NCORES = 4
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

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics


n_moments = 2
field_table_name = "field_table"
write_ft(field_table_name,n_moments)
        
          
exp = Experiment('realistic_continents_fixed_sst_wv_age_test_experiment', codebase=cb)

#Add any input files that are necessary for a particular experiment.
exp.inputfiles = [os.path.join(GFDL_BASE,'input/land_masks/era_land_t42.nc'),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
                      os.path.join(base_dir,'input/sst_clim_amip.nc'), os.path.join(base_dir,'input/siconc_clim_amip.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')


#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf') #need at least ps, pk, bk and zsurf to do vertical interpolation onto plevels from sigma
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

# WV age moments diagnostics
for ind in range(n_moments):
    name = f"sphum_age_{ind+1}"
    diag.add_field('dynamics', name, time_avg=True)

exp.diag_table = diag


#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
namelist_name = os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents/namelist_basefile.nml')
nml = f90nml.read(namelist_name)
exp.namelist = nml

exp.update_namelist({
    'mixed_layer_nml': {  
        'do_qflux' : False, #Don't use the prescribed analytical formula for q-fluxes
        'do_read_sst' : True, #Read in sst values from input file
        'do_sc_sst' : True, #Do specified ssts (need both to be true)
        'sst_file' : 'sst_clim_amip', #Set name of sst input file
        'specify_sst_over_ocean_only' : True, #Make sure sst only specified in regions of ocean.
    }
})

#Lets do a run!
if __name__=="__main__":
    cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2,121):
        exp.run(i, num_cores=NCORES)
