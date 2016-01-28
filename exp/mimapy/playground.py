import mima

exp = mima.Experiment('playground', overwrite_data=True)

#exp.clear_workdir()
exp.compile()


#exp.namelist['constants_nml'] = {'omega': 6.5e-5}


exp.clear_rundir()

exp.runmonth(1, use_restart=False)