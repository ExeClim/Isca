import numpy as np

import mima

orbit_rate = (2*np.pi) / (360*86400)

for i, omega in enumerate(orbit_rate*np.array([1, 2, 10, 50, 360, 720])):
	exp = mima.Experiment('greyrot_omega_%d' % (i+1))

	# use the templated constants.F90 and change omega
	exp.use_template_file('constants.F90', {'omega': omega})

	exp.compile()
	exp.clear_workdir()

	exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
	exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
	exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = False
	exp.namelist['spectral_dynamics_nml']['num_levels'] = 20


	#restart_file = mima.Experiment('grey_nonseasonal').get_restart_file(60))

	exp.runmonth(1, use_restart=False)
	for i in range(100):
		exp.runmonth(i)
