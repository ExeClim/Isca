import mima

exp = mima.Experiment('grey_nonseasonal2')

exp.clear_workdir()
exp.compile()

exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False


exp.namelist.update({
    'spectral_dynamics_nml': {
        'num_levels': 40
    },
})

exp.runmonth(1, num_cores=8)