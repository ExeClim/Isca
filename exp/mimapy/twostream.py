import mima

exp = mima.Experiment('grey_nonseasonal')

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

exp.runmonth(1, use_restart=False)
for i in range(2, 80):
    exp.runmonth(i, num_cores=8)