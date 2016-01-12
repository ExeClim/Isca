import mima
import os

nonseasonal = mima.Experiment('grey_nonseasonal')

nonseasonal.clear_workdir()
nonseasonal.compile()

nonseasonal.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
nonseasonal.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
nonseasonal.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False

# nonseasonal.namelist.update({
#     'spectral_dynamics_nml': {
#         'num_levels': 40
#     },
# })

nonseasonal.runmonth(1, use_restart=False)
for i in range(2, 82):
    nonseasonal.runmonth(i)

#nonseasonal.runmonth(81, overwrite_data=True)


# Create the same experiment, but with the seasonal cycle switched on
seasonal = mima.Experiment('grey_seasonal')
seasonal.clear_workdir()
seasonal.execdir = nonseasonal.execdir  # use the same executable as above

seasonal.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
seasonal.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
seasonal.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True

seasonal.runmonth(1, restart_file=nonseasonal.get_restart_file(60))
for i in range(2, 40):
	seasonal.runmonth(i)




