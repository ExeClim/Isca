import os

import numpy as np
from field_table_write import write_ft
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


A_hyb = (1/100) * np.array([
    0.000000, 20.006149, 43.297810, 75.346230, 115.082146, 161.897491,
    215.896912, 350.138184, 539.651489,
    828.398987, 1026.366943, 1271.644531, 1575.537842, 1952.054443,
    2418.549805, 2996.526611, 3712.626221, 4599.856934, 5699.114746,
    6998.388184, 8507.411133, 10181.707031, 11883.089844, 13442.915039,
    14736.354492, 15689.206055, 16266.609375, 16465.003906, 16297.620117,
    15791.597656, 14985.269531, 13925.519531, 12665.294922, 11261.230469,
    9771.406250, 8253.210938, 6761.339844, 5345.917969, 4050.718750,
    2911.570312, 1954.804688, 1195.890625, 638.148438, 271.625000,
    72.062500, 0.000000])

# B coefficients (unitless)
B_hyb = np.array([
    0.000000, 0.000000,0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000100,
    0.000673, 0.003163, 0.009292, 0.020319, 0.036975, 0.059488,
    0.087895, 0.122004, 0.161442, 0.205703, 0.254189, 0.306235,
    0.361145, 0.418202, 0.476688, 0.535887, 0.595084, 0.653565,
    0.710594, 0.765405, 0.817167, 0.864956, 0.907716, 0.946485162,1.00000000e+00
])

A_hyb_v2 = np.array([0.00000000e+00, 7.36780167e-02, 2.65109360e-01, 7.97618508e-01,
       2.05185151e+00, 4.60448742e+00, 9.17678547e+00, 1.65052547e+01,
       2.18396397e+01, 2.71740246e+01, 3.43231888e+01, 4.14723549e+01,
       5.04006462e+01, 5.93289337e+01, 6.98336868e+01, 8.03384399e+01,
       9.20988922e+01, 1.03859337e+02, 1.29142944e+02, 1.47363545e+02,
       1.56892061e+02, 1.62666094e+02, 1.64650039e+02, 1.62976201e+02,
       1.57915977e+02, 1.49852695e+02, 1.39255195e+02, 1.26652949e+02,
       1.12612305e+02, 9.77140625e+01, 8.25321094e+01, 6.76133984e+01,
       5.34591797e+01, 4.05071875e+01, 2.91157031e+01, 1.95480469e+01,
       1.19589063e+01, 6.38148438e+00, 2.71625000e+00, 7.20625000e-01,
       0.00000000e+00])

B_hyb_v2 = np.array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e-04,
       6.73000000e-04, 3.16300000e-03, 9.29200000e-03, 2.03190000e-02,
       3.69750000e-02, 5.94880000e-02, 8.78950000e-02, 1.22004000e-01,
       1.61442000e-01, 2.05703000e-01, 2.54189000e-01, 3.06235000e-01,
       3.61145000e-01, 4.18202000e-01, 4.76688000e-01, 5.35887000e-01,
       5.95084000e-01, 6.53565000e-01, 7.10594000e-01, 7.65405000e-01,
       8.17167000e-01, 8.64956000e-01, 9.07716000e-01, 9.46485162e-01,
       1.00000000e+00])

B_default = np.array([0.00000000e+00, 7.36780203e-05, 2.65109353e-04, 7.97618530e-04,
                    2.05185148e-03, 4.60448721e-03, 9.17678513e-03, 1.65052544e-02,
                    2.18396392e-02, 2.71740239e-02, 3.43231894e-02, 4.14723530e-02,
                    5.04006445e-02, 5.93289323e-02, 6.98336884e-02, 8.03384408e-02,
                    9.20988917e-02, 1.03859335e-01, 1.29142940e-01, 1.55455455e-01,
                    1.82168067e-01, 2.08807260e-01, 2.35068962e-01, 2.60807067e-01,
                    2.86007226e-01, 3.10755134e-01, 3.35205764e-01, 3.59556884e-01,
                    3.84028047e-01, 4.08845514e-01, 4.34231699e-01, 4.60398972e-01,
                    4.87546176e-01, 5.15857577e-01, 5.45503139e-01, 5.76640010e-01,
                    6.09414339e-01, 6.43963873e-01, 6.80420518e-01, 7.18912899e-01,
                    7.59568930e-01, 8.02518070e-01, 8.47893596e-01, 8.95834148e-01,
                    9.46485162e-01, 1.00000000e+00])
A_default =  np.array([0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
                    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                    0.000000])


def start_sst_exp(exp_name,delta_sst, n_moments,n_simulation_months,namelist,
                  horizontal_resolution, vertical_resolution=40,
                  use_restart = False,restart_file ="None",n_start_months = 1,dt_atm = 720,sponge = 150,trayfric = -0.5):
    
    if horizontal_resolution == "T42":
        NCORES = 32
        name_era_land = 'input/era_land_t42.nc'
    elif horizontal_resolution == "T85":
        NCORES = 64
        name_era_land = 'input/era_land_t85.nc'
        
    RESOLUTION = horizontal_resolution, vertical_resolution
    base_dir = os.path.dirname(os.path.realpath(__file__))
    # a CodeBase can be a directory on the computer,
    # useful for iterative development
    cb = IscaCodeBase.from_directory(GFDL_BASE)
    #cb.compile(debug = False)
    field_table_name = "field_table_age_" + str(n_moments)
    write_ft(field_table_name,n_moments)
    
    exp = Experiment(exp_name,ext_field_table=field_table_name ,codebase=cb)

    sst_file_name = f'input/sst_clim_amip({delta_sst}).nc'
    exp.inputfiles = [os.path.join(base_dir,name_era_land),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
                      os.path.join(base_dir,sst_file_name), os.path.join(base_dir,'input/siconc_clim_amip.nc')]

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

    #diag.add_field('vert_turb', 'z_pbl', time_avg=True)

    for ind in range(n_moments):
        name = f"sphum_age_{ind+1}"
        diag.add_field('dynamics', name, time_avg=True)

    diag.add_field('atmosphere', 'cape', time_avg=True)
    diag.add_field('atmosphere', 'dt_qg_convection', time_avg=True)
    diag.add_field('atmosphere', 'dt_qg_condensation', time_avg=True)
    diag.add_field('atmosphere', 'dt_sink', time_avg=True)
    diag.add_field('atmosphere', 'dt_tracer', time_avg=True)
    diag.add_field('atmosphere', 'dt_q', time_avg=True)
    diag.add_field('atmosphere', 'dt_tracer_diff', time_avg=True)
    diag.add_field('atmosphere', 'dt_qg_diffusion', time_avg=True)
    diag.add_field('atmosphere', 'rh', time_avg=True)
    diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
    diag.add_field('atmosphere', 'convection_rain', time_avg=True)
    
    diag.add_field('rrtm_radiation', 'olr', time_avg=True)
    diag.add_field('rrtm_radiation', 'toa_sw', time_avg=True)
    diag.add_field('rrtm_radiation', 'tdt_rad', time_avg=True)
    diag.add_field('rrtm_radiation', 'tdt_sw', time_avg=True)
    diag.add_field('rrtm_radiation', 'tdt_lw', time_avg=True)
    diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True)
    diag.add_field('rrtm_radiation', 'flux_lw', time_avg=True)

    exp.diag_table = diag


    #Empty the run directory ready to run
    exp.clear_rundir()

    #Define values for the 'core' namelist
    namelist_name = os.path.join(GFDL_BASE, f'exp/test_cases/realistic_continents/{namelist}.nml')
    nml = f90nml.read(namelist_name)
    exp.namelist = nml


    exp.update_namelist({
    'mixed_layer_nml': {  
        'do_qflux' : False, #Don't use the prescribed analytical formula for q-fluxes
        'do_read_sst' : True, #Read in sst values from input file
        'do_sc_sst' : True, #Do specified ssts (need both to be true)
        'sst_file' : f'sst_clim_amip({delta_sst})', #Set name of sst input file
        'specify_sst_over_ocean_only' : True, #Make sure sst only specified in regions of ocean.
    },
    'main_nml' : {
        'dt_atmos' : dt_atm
        },
    'damping_driver_nml': { "sponge_pbottom" : sponge , "trayfric" : trayfric}
  ,
    'vert_coordinate_nml': {
            'bk' : B_default,
            'pk' : A_default,
            }
    })

    exp.set_resolution(*RESOLUTION)
    
    n_max = n_simulation_months + 1
    # Run first month with or without the restart_file
    exp.run(n_start_months,restart_file = restart_file,use_restart=use_restart, num_cores=NCORES,overwrite_data = True)
    for i in range(n_start_months + 1,n_max):
        exp.run(i, num_cores=NCORES,overwrite_data = True)
        
    return "Experiment ran succesfully"