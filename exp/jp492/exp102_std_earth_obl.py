from exp101_standard_earth import exp as baseexp

for obl in [0.0, 5.0, 10.0, 20.0]:
    exp = baseexp.derive('exp102_obl_%d' % obl)
    exp.update_namelist({
        'astronomy_nml': {
            'obliq': obl,
        }
    })

    exp.clear_rundir()
    exp.run(1, use_restart=False, num_cores=16)
    for i in range(1, 2*5):
        exp.run(i+1, num_cores=16)