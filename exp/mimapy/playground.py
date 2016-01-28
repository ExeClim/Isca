import mima

exp = mima.Experiment('playground', overwrite_data=True)

#exp.clear_workdir()
exp.compile()


exp.namelist['constants_nml'] = {
    'grav': 10.0
}


exp.clear_rundir()

exp.runmonth(1, use_restart=False)