The phys namelist (phys.nml) is sometimes really hard to work out what the flags refer to. Here is a list of the flags and what they mean.

| group name   |  name list  |  What is it?  | comment  |
|---|---|---|---|
| idealized_moist_phys_nml   | do_rrtm_radiation   | radiation flag    | turn mima radiation on. Default false |   
| idealized_moist_phys_nml   | two_stream_gray     | radiation flag    | do grey radiation. Default true |   
| idealized_moist_phys_nml   | lwet_convection     | convection flag   |  turn on Simple Betts-Miller (SBMS |   
| idealized_moist_phys_nml   | do_bm               | convection flag   |  turn on Betts-Miller BMS |   
| idealized_moist_phys_nml   | turb                | turb flag         |  turn on vert_turb_driver. Default false |   
| idealized_moist_phys_nml   | do_damping          | damping flag      |  turn on damping_driver. Default false |   
| idealized_moist_phys_nml   | mixed_layer_bc      | ml flag           |  turn on mixed_layer. Default false |   
| idealized_moist_phys_nml   | do_virtual          |                   |  use virtual temp in gcm_vert_diff. Default false |
| idealized_moist_phys_nml   | do_simple           |                   |  Default false. RH calculation choice? |
| idealized_moist_phys_nml   | roughness_***       |                   |  Set roughness lengths, m. Default 0.05 |
| idealized_moist_phys_nml   | land_option         | land flag         | Choose land input. Default 'none. 'input' read from file. 'zsurf' base on topography height. |
| idealized_moist_phys_nml   | land_roughness_prefactor | land flag         |  Scale factor for roughness lengths over land for 'input' option. |
| idealized_moist_phys_nml   | land_file_name      | land flag         |  Land file name (also set in runscript?) |
| idealized_moist_phys_nml   | land_field_name     | land flag         |  Land field name  |
| idealized_moist_phys_nml   | mixed_layer_bc      | ml flag           |  turn on mixed_layer. Default false |   
| qe_moist_convection_nml    | tau_bm              | convection flag for SBMS| default is 7200.0 |
| qe_moist_convection_nml    | rhbm                | convection flag for SBMS| default is 0.7    |
| spectral_dynamics_nml      | surf_res            | trop vs strat distribution of sigma levels | Larger surf_res gives more levels in strat vs trop when using uneven-sigma option. 0.2 new recommended value.
| | | |
