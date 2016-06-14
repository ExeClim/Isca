import gfdl.experiment


basic = gfdl.experiment.DiagTable()

basic.add_file('hourly', 1, 'hours')
basic.add_field('two_stream', 'coszen', files=['hourly'])

basic.add_file('daily', 1, 'days')

basic.add_field('dynamics', 'ps', files=['daily'])
basic.add_field('dynamics', 'bk', files=['daily'])
basic.add_field('dynamics', 'pk', files=['daily'])
basic.add_field('dynamics', 'ucomp', files=['daily'])
basic.add_field('dynamics', 'vcomp', files=['daily'])
basic.add_field('dynamics', 'omega', files=['daily'])
basic.add_field('dynamics', 'temp', files=['daily'])
basic.add_field('dynamics', 'vor', files=['daily'])
basic.add_field('dynamics', 'div', files=['daily'])
basic.add_field('dynamics', 'sphum', files=['daily'])

basic.add_field('two_stream', 'olr', files=['daily'])
basic.add_field('two_stream', 'flux_sw', files=['daily'])
basic.add_field('two_stream', 'flux_lw', files=['daily'])
basic.add_field('two_stream', 'tdt_rad', files=['daily'])

basic.add_field('mixed_layer', 't_surf', files=['daily'])

basic.add_file('hourly', 1, 'hours')
basic.add_field('two_stream', 'coszen', files=['hourly'])

moist = basic.copy()

dry = gfdl.experiment.DiagTable()

dry.add_file('daily', 1, 'days')
#dry.add_file('6hourly', 6, 'hours')

dry.add_field('dynamics', 'ps')
dry.add_field('dynamics', 'bk')
dry.add_field('dynamics', 'pk')
dry.add_field('dynamics', 'ucomp')
dry.add_field('dynamics', 'vcomp')
dry.add_field('dynamics', 'temp')
dry.add_field('dynamics', 'vor')
dry.add_field('dynamics', 'div')
dry.add_field('dynamics', 'sphum')

dry.add_field('two_stream', 'olr')
dry.add_field('two_stream', 'flux_sw')
dry.add_field('two_stream', 'flux_lw')
dry.add_field('two_stream', 'tdt_rad')

dry.add_field('mixed_layer', 't_surf')

dry.add_field('dry_convection', 'CAPE')
dry.add_field('dry_convection', 'CIN')
dry.add_field('dry_convection',  'LZB')
dry.add_field('dry_convection', 'LCL')


dry.add_file('hourly', 1, 'hours')

dry.add_field('two_stream', 'coszen', files=['hourly'])