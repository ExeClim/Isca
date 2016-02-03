The phys namelist (phys.nml) is sometimes really hard to work out what the flags refer to. Here is a list of the flags and what they mean.

| group name   |  name list  |  What is it?  | comment  |
|---|---|---|---|
| idealized_moist_phys_nml   | do_rrtm_radiation   | radiation flag    | turn mima radiation on  |   
| idealized_moist_phys_nml   | two_stream_gray     | radiation flag    | do grey radiation     |   
| idealized_moist_phys_nml   | lwet_convection     | convection flag   |  turn on Simple Betts-Miller (SBMS |   
| idealized_moist_phys_nml   | do_bm               | convection flag   |  turn on Betts-Miller BMS |   
| qe_moist_convection_nml    | tau_bm              | convection flag for SBMS| default is 7200.0 |
| qe_moist_convection_nml    | rhbm                | convection flag for SBMS| default is 0.7    |
| | | |
| | | |
