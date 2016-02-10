|Module name   	| Field name  	|  Field long name 	|  Dimension (not including time) 	| Description (if needed)  	|
|---	          |---	          |---	              |---	                              |---	                      |
|dynamics |	u_comp	|u wind	| (pfull, lat, lon)	| zonal component of the horizontal winds |
|dynamics |	v_comp	|v wind	| (pfull, lat, lon)	| meridional component of the horizontal winds|
|dynamics |	ps	    |surface pressure	| (lat, lon)	| slab |
|dynamics |	slp	    |sea level pressure	| (lat,lon)	|
|dynamics |	bk	    |vertical coordinate sigma values	| (phalf)	|
|dynamics |	pk	    |vertical coordinate pressure values	| (phalf)	|
|dynamics |	vor	    |Vorticity	| (pfull, lat, lon)	|
|dynamics |	div	    |Divergence	| (pfull, lat, lon)	|
|dynamics |	temp	  |atmospheric tempertaure	| (pfull, lat, lon)	|
|dynamics |	omega	  |vertical velocity	| (pfull, lat, lon)	dp/dt|
|mixed_layer | t_surf |	surface temperature	| (lat, lon) | slab |
|mixed_layer | albedo | surface albedo | (lat, lon) | static |
|mixed_layer | ml_heat_cap | mixed layer heat capacity | (lat, lon) | 
|atmosphere |	rh	  | relative humidity |	(pfull, lat, lon)	|
|atmosphere |	convection_rain	| Rain from convection scheme	| (lat, lon)	|
|atmosphere |	condensation_rain |	Large scale resolved rain	| (lat, lon)	|
|atmosphere |	precipitation |	Sum of resolved rain, resolved snow and parameterised rain 	| (lat, lon)	|
|rrtm_radiation |  z_thalf 	|   	|   	|   	|
|rrtm_radiation |	tdt_rad	| Temperature tendency due to radiation | 	(pfull, lat, lon)	|
|rrtm_radiation |	flux_sw	| Net short wave surface flux	| (lat, lon)	|
|rrtm_radiation |	flux_lw	| Long wave surface flux	| (lat, lon)	|
|rrtm_radiation |	rrtm_albedo	 |Interactive albedo	| (lat, lon)	|
|rrtm_radiation | 	tdt_sw	 | Temperature tendency due to SW radiation	| (pfull, lat, lon)	|
|rrtm_radiation |	tdt_lw	| Temperature tendency due to LW radiation	| (pfull, lat, lon)	|
|damping |	udt_rdamp |	u wind tendency for Rayleigh damping	| (pfull, lat, lon)	|
|damping |	vdt_rdamp |	v wind tendency for Rayleigh damping	| (pfull, lat, lon)	|
|damping |	tdt_diss_rdamp |	Dissipative heating from Rayleigh damping	| (pfull, lat, lon)	|
|vert_turb|	z_pbl	| depth of planetary boundary layer	| (lat, lon)| 	
