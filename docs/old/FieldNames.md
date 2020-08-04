|Module name   	| Field name  	|  Field long name 	|  Dimension (not including time) 	| Description (if needed)  	| Units |
|---	          |---	          |---	              |---	                              |---	                      |--- |
|dynamics |	ucomp	|zonal wind u	| (pfull, lat, lon)	| zonal component of the horizontal winds | m/sec |
|dynamics |	vcomp	|meridional wind v	| (pfull, lat, lon)	| meridional component of the horizontal winds| m/sec |
|dynamics |	ps	    |surface pressure	| (lat, lon)	| slab | pascals |
|dynamics |	slp	    |sea level pressure	| (lat,lon)	| | pascals |
|dynamics |	bk	    |vertical coordinate sigma values	| (phalf)	| if pk = 0, then bk = sigma | (dimensionless) |
|dynamics |	pk	    |vertical coordinate pressure values	| (phalf)	| if bk = 0, then pk = pressure (i.e. phalf, BUT NB bk in Pa whereas phalf (& pfull) in hPa!) | pascals |
|dynamics |	pres_full   |pressure at full model levels	| (pfull, lat, lon)	| | pascals |
|dynamics |	pres_half   |pressure at half model levels	| (phalf, lat, lon)	| | pascals |
|dynamics |	height      |geopotential height at full model levels	| (pfull, lat, lon)	| | m |
|dynamics |	height_half |geopotential height at half model levels	| (phalf, lat, lon)	| | m |
|dynamics |	vor	    |Vorticity	| (pfull, lat, lon)	| | sec**-1 |
|dynamics |	div	    |Divergence	| (pfull, lat, lon)	| | sec**-1 |
|dynamics |	temp	  |atmospheric tempertaure	| (pfull, lat, lon)	| | deg_k |
|dynamics |	sphum	  | specific humidity	| (pfull, lat, lon)	| | |
|dynamics |	omega	  |vertical velocity	| (pfull, lat, lon)	| dp/dt | Pa/sec |
|dynamics |	wspd	| horizontal wind speed | (pfull, lat, lon)	| sqrt(u^2 + v^2) | m/sec |
|dynamics |	ucomp_sq	|zonal wind squared	| (pfull, lat, lon)	|  | (m/sec)**2 |
|dynamics |	ucomp_vcomp	|zonal wind * meridional wind	| (pfull, lat, lon)	|  | (m/sec)**2 |
|dynamics |	ucomp_omega	|zonal wind * vertical wind	| (pfull, lat, lon)	|  | m*Pa/sec**2 |
|dynamics |	vcomp_sq	|meridional wind squared	| (pfull, lat, lon)	|  | (m/sec)**2 |
|dynamics |	vcomp_omega	|meridional wind * vertical wind	| (pfull, lat, lon)	|  | m*Pa/sec**2 |
|dynamics |	omega_sq	|omega squared	| (pfull, lat, lon)	| | (Pa/sec)**2 |
|dynamics |	ucomp_temp	|zonal wind * temperature	| (pfull, lat, lon)	|  | m*K/sec |
|dynamics |	vcomp_temp	|meridional wind * temperature	| (pfull, lat, lon)	|  | m*K/sec |
|dynamics |	omega_temp	|dp/dt * temperature	| (pfull, lat, lon)	|  | Pa*K/sec |
|dynamics |	temp_sq	  |temp squared	| (pfull, lat, lon)	| | deg_k**2 |
|dynamics |	vcomp_vor	|meridional wind * vorticity	| (pfull, lat, lon)	|  | m/sec**2 |
|mixed_layer | t_surf |	surface temperature	| (lat, lon) | slab | |
|mixed_layer | albedo | surface albedo | (lat, lon) | static | |
|mixed_layer | ml_heat_cap | mixed layer heat capacity | (lat, lon) |  | |
|mixed_layer |	flux_lhe	| latent heat flux up at surface	| (lat,lon) | | |
|mixed_layer |	flux_sw	| Net shortwave radiative flux (positive up)	| (phalf,lat,lon) | | |
|mixed_layer |	flux_lw	| Net longwave radiative flux (positive up)	 | (phalf,lat,lon) | | |
|atmosphere |	rh	  | relative humidity |	(pfull, lat, lon)	| | |
|atmosphere |	convection_rain	| Rain from convection scheme	| (lat, lon)	| | |
|atmosphere |	condensation_rain |	Large scale resolved rain	| (lat, lon)	| | |
|atmosphere |	precipitation |	Sum of resolved rain, resolved snow and parameterised rain 	| (lat, lon)	| | |
|rrtm_radiation |  z_thalf 	|   	|   	|   	| |
|rrtm_radiation |	flux_sw	| Net short wave surface flux	| (lat, lon)	| | |
|rrtm_radiation |	flux_lw	| Long wave surface flux	| (lat, lon)	| | |
|rrtm_radiation |	rrtm_albedo	 |Interactive albedo	| (lat, lon)	| | |
|rrtm_radiation | 	tdt_sw	 | Temperature tendency due to SW radiation	| (pfull, lat, lon)	| | |
|rrtm_radiation |	tdt_lw	| Temperature tendency due to LW radiation	| (pfull, lat, lon)	| | |
|rrtm_radiation, two_stream |	tdt_rad	| Temperature tendency due to radiation | 	(pfull, lat, lon)	| | |
|damping, hs_forcing |	udt_rdamp |	u wind tendency for Rayleigh damping	| (pfull, lat, lon)	| | m/s2 |
|damping, hs_forcing |	vdt_rdamp |	v wind tendency for Rayleigh damping	| (pfull, lat, lon)	| | m/s2 |
|damping, hs_forcing |	tdt_diss_rdamp |	Dissipative heating from Rayleigh damping	| (pfull, lat, lon)	| | deg/sec |
|hs_forcing |  	tdt_ndamp	| Heating due to newtonian damping | (pfull,lat,lon) | | deg/sec |
|hs_forcing |  	local_heating	| Local heating  | (pfull,lat,lon) | | deg/sec |
|hs_forcing |  	tdt	        | Total heating: newtonian damping + local heating | (pfull,lat,lon) | Also includes dissipative heating if applicable | deg/sec |
|hs_forcing |   	teq	| Equilibration temperature | (pfull,lat,lon) | | deg_K |
|vert_turb|	z_pbl	| depth of planetary boundary layer	| (lat, lon)| 	 | |


The following variables need to be saved to do pressure interpolation: bk, pk and ps. These should be saved for each output frequency (i.e daily and monthly if both are used).
