

      module rrtm_vars
!
!    Modeling an idealized Moist Atmosphere (MiMA)
!    Copyright (C) 2015  Martin Jucker
!    https://github.com/mjucker/MiMA
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or 
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!
!   RRTM_VARS:
!   Contains all variables needed to 
!   run the RRTM code, version for GCMs (hence the 'G'),
!   other than astronomy, i.e. all variables needed
!   for radiation that are not within astro.f90
! 
!   external modules
        use parkind, only         : im => kind_im, rb => kind_rb
        use interpolator_mod, only: interpolate_type
!
!  rrtm_radiation variables
!
        implicit none

        logical                                    :: rrtm_init=.false.    ! has radiation been initialized?
        type(interpolate_type),save                :: o3_interp            ! use external file for ozone
        type(interpolate_type),save                :: h2o_interp           ! use external file for water vapor
        type(interpolate_type),save                :: co2_interp           ! use external file for co2
        type(interpolate_type),save                :: rad_interp           ! use external file for radiation
        type(interpolate_type),save                :: fsw_interp           ! use external file for SW fluxes
        type(interpolate_type),save                :: flw_interp           ! use external file for SLW fluxes
        integer(kind=im)                           :: ncols_rrt,nlay_rrt   ! RRTM field sizes
                                                                           ! ncols_rrt = (size(lon)/lonstep*
                                                                           !             size(lat)
                                                                           ! nlay_rrt  = size(pfull)
        ! gas volume mixing ratios, dimensions (ncols_rrt,nlay=nlevels)
        ! vmr = mass mixing ratio [mmr,kg/kg] scaled with molecular weight [g/mol]
        ! any modification of the units must be accounted for in rrtm_xw_rad.nomica.f90, where x=(s,l)
        real(kind=rb),allocatable,dimension(:,:)   :: h2o                  ! specific humidity [kg/kg]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: o3                   ! ozone [mmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: co2                  ! CO2 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: zeros                ! place holder for any species set
                                                                           !  to zero
        real(kind=rb),allocatable,dimension(:,:)   :: ones                 ! place holder for secondary species
        ! the following species are only set if use_secondary_gases=.true.
        real(kind=rb),allocatable,dimension(:,:)   :: ch4                  ! CH4 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: n2o                  ! N2O [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: o2                   ! O2  [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc11                ! CFC11 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc12                ! CFC12 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: cfc22                ! CFC22 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: ccl4                 ! CCL4 [vmr]
                                                                           ! dimension (ncols_rrt x nlay_rrt)
        real(kind=rb),allocatable,dimension(:,:)   :: emis                 ! surface LW emissivity per band
                                                                           ! dimension (ncols_rrt x nbndlw)
                                                                           ! =1 for black body
        ! clouds stuff
        !  cloud & aerosol optical depths, cloud and aerosol specific parameters. Set to zero
        real(kind=rb),allocatable,dimension(:,:,:) :: taucld,tauaer, sw_zro, zro_sw
        ! heating rates and fluxes, zenith angle when in-between radiation time steps
        real(kind=rb),allocatable,dimension(:,:)   :: sw_flux,lw_flux,zencos, olr, toa_sw! surface and TOA fluxes, cos(zenith angle) 
                                                                            ! dimension (lon x lat)
        real(kind=rb),allocatable,dimension(:,:,:) :: tdt_rad               ! heating rate [K/s]
                                                                            ! dimension (lon x lat x pfull)
        real(kind=rb),allocatable,dimension(:,:,:) :: tdt_sw_rad,tdt_lw_rad ! SW, LW radiation heating rates,
                                                                            ! diagnostics only [K/s]
                                                                            ! dimension (lon x lat x pfull)
        real(kind=rb),allocatable,dimension(:,:,:) :: t_half                ! temperature at half levels [K]
                                                                            ! dimension (lon x lat x phalf)
                                                                            

        real(kind=rb),allocatable,dimension(:,:)   :: rrtm_precip           ! total time of precipitation
                                                                            ! between radiation steps to
                                                                            ! determine precip_albedo
                                                                            ! dimension (lon x lat)
        integer(kind=im)                           :: num_precip            ! number of times precipitation
                                                                            ! has been summed in rrtm_precip
        real(kind=rb)                              :: dt_last               ! time of last radiation calculation
                                                                            ! used for alarm
!---------------------------------------------------------------------------------------------------------------
! some constants
        real(kind=rb)      :: daypersec=1./86400.,deg2rad   !RG: daypersec=1./86400. left in when conversion to non-specific day length made as this only converts heatrates from RRTM from K/day to K/sec
! no clouds in the radiative scheme
        integer(kind=im) :: icld=0,idrv=0, &
             inflglw=0,iceflglw=0,liqflglw=0, &
             iaer=0
!---------------------------------------------------------------------------------------------------------------
!                                namelist values
!---------------------------------------------------------------------------------------------------------------
! input files: file names are always given without '.nc', which is always assumed
!  the field to be read within the file needs to have the same name as the file
! OBS! reading of radiation is still experimental at this point!
        logical            :: do_read_radiation=.false.       ! read SW and LW radiation in the atmosphere from
                                                              !  external file? Surface fluxes are still computed
        character(len=256) :: radiation_file='radiation'      !  file name to read radiation
        real(kind=rb)      :: rad_missing_value=-1.e19        !   missing value in input files:
                                                              !    if <0, replace everything below this value with 0
                                                              !    if >0, replace everything above this value with 0
        logical            :: do_read_sw_flux=.false.         ! read SW surface fluxes from external file?
        character(len=256) :: sw_flux_file='sw_flux'          !  file name to read fluxes
        logical            :: do_read_lw_flux=.false.         ! read LW surface fluxes from external file?
        character(len=256) :: lw_flux_file='lw_flux'          !  file name to read fluxes
        logical            :: do_read_ozone=.false.           ! read ozone from an external file?
                                                              !  this is the only way to get ozone into the model
        character(len=256) :: ozone_file='ozone'              !  file name of ozone file to read
        logical            :: input_o3_file_is_mmr=.true.     ! Does the ozone input file contain values as a mass mixing ratio (set to true) or a volume mixing ratio (set to false)?
        logical            :: do_read_h2o=.false.             ! read water vapor from an external file?
        character(len=256) :: h2o_file='h2o'                  !  file name of h2o file to read
        logical            :: do_read_co2=.false.             ! read co2 concentration from an external file?
        character(len=256) :: co2_file='co2'                  !  file name of co2 file to read
        character(len=256) :: co2_variable_name='co2'         !  field name of co2 file to read

! secondary gases (CH4,N2O,O2,CFC11,CFC12,CFC22,CCL4)
        logical            :: include_secondary_gases=.false. ! non-zero values for above listed secondary gases?
        real(kind=rb)      :: ch4_val  = 0.                   !  if .true., value for CH4 vmr
        real(kind=rb)      :: n2o_val  = 0.                   !                       N2O vmr
        real(kind=rb)      :: o2_val   = 0.                   !                       O2 vmr
        real(kind=rb)      :: cfc11_val= 0.                   !                       CFC11
        real(kind=rb)      :: cfc12_val= 0.                   !                       CFC12
        real(kind=rb)      :: cfc22_val= 0.                   !                       CFC22
        real(kind=rb)      :: ccl4_val = 0.                   !                       CCL4
! some safety boundaries
        real(kind=rb)      :: h2o_lower_limit = 2.e-7         ! never use smaller than this in radiative scheme
        real(kind=rb)      :: temp_lower_limit = 100.         ! never go below this in radiative scheme
        real(kind=rb)      :: temp_upper_limit = 370.         ! never go above this in radiative scheme
! primary gases: CO2 and H2O
        real(kind=rb)      :: co2ppmv=300.                    ! CO2 ppmv concentration
        logical            :: do_fixed_water = .false.        ! feed fixed value for water vapor to RRTM?
        real(kind=rb)      :: fixed_water = 2.e-06            ! if so, what value? [kg/kg]
        real(kind=rb)      :: fixed_water_pres = 100.e02      ! if so, above which pressure level? [hPa]
        real(kind=rb)      :: fixed_water_lat  = 90.          ! if so, equatorward of which latitude? [deg]
        logical            :: do_zm_tracers=.false.           ! Feed only the zonal mean of tracers to radiation
        logical            :: convert_sphum_to_vmr=.true.     ! Model is fed sphum, but RRTM wants vmr. Set to true to make this conversion. May want false if using do_read_h2o and input file is a vmr.             
! radiation time stepping and spatial sampling
        integer(kind=im)   :: dt_rad=0                        ! Radiation time step - every step if dt_rad<dt_atmos
        logical            :: store_intermediate_rad =.true.  ! Keep rad constant over entire dt_rad?
                                                              ! Else only heat radiatively at every dt_rad
        logical            :: do_rad_time_avg =.true.         ! Average coszen for SW radiation over dt_rad?
        integer(kind=im)   :: dt_rad_avg = -1                 ! If averaging, over what time? dt_rad_avg=dt_rad if dt_rad_avg<=0
        integer(kind=im)   :: lonstep=1                       ! Subsample fields along longitude
                                                              !  for faster radiation calculation

!!!!!! mp586 added for annual mean insolation!!!!!

		logical            :: frierson_solar_rad =.false.
		real(kind=rb)	   :: del_sol = 0.95 ! frierson 2006 default = 1.4, but 0.95 gets the curve closer to the annual mean insolation 
		real(kind=rb)	   :: del_sw = 0.0 !frierson 2006 default 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! some fancy radiation tweaks  
        real(kind=rb)      :: slowdown_rad = 1.0              ! factor do simulate slower seasonal cycle: >1 means faster, <1 slower
        logical            :: do_zm_rad=.false.               ! Only compute zonal mean radiation
        logical            :: do_precip_albedo=.false.        ! Modify albedo depending on large scale
                                                              !  precipitation (crude cloud parameterization)
        real(kind=rb)      :: precip_albedo=0.35              ! If so, what's the cloud albedo?
        real(kind=rb)      :: precip_lat = 0.0                ! If so, poleward of which latitude should it be applied?
        character(len=14)  :: precip_albedo_mode = 'full'     ! If so, use 
                                                              !  full precipitation ('full')
                                                              !  only large scale condensation ('lscale')
                                                              !  only convection ('conv')

        logical            :: use_dyofyr = .false.            ! use day of the year to compute Earth-Sun distance?
                                                              !  this is done within RRTM, and assumes 365days/year!
        real(kind=rb)      :: solrad=1.0                      ! distance Earth-Sun [AU] if use_dyofyr=.false.
        integer(kind=im)   :: solday=0                        ! if >0, do perpetual run corresponding to
                                                              !  day of the year = solday \in [0,days per year]
        real(kind=rb)      :: equinox_day=0.75                ! fraction of the year defining NH autumn equinox \in [0,1]
        real(kind=rb)      :: solr_cnst= 1368.22              ! solar constant [W/m2]
!-------------------------------------------------s--------------------------------------------------------------
!
!-------------------- diagnostics fields -------------------------------

        integer :: id_tdt_rad, id_tdt_sw, id_tdt_lw, id_coszen, id_flux_sw, id_flux_lw, id_olr, id_toa_sw, id_albedo,id_ozone, id_co2, id_fracday, id_half_level_temp, id_full_level_temp
        character(len=14), parameter :: mod_name_rad = 'rrtm_radiation' !s changed parameter name from mod_name to mod_name_rad as compiler objected, presumably because mod_name also defined in idealized_moist_physics.F90 after use rrtm_vars is included. 
        real :: missing_value = -999.

!---------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------

        namelist/rrtm_radiation_nml/ include_secondary_gases, do_read_ozone, ozone_file, input_o3_file_is_mmr, &
             &do_read_h2o, h2o_file, convert_sphum_to_vmr, ch4_val, n2o_val, o2_val, cfc11_val, cfc12_val, cfc22_val, ccl4_val, &
             &do_read_radiation, radiation_file, rad_missing_value, &
             &do_read_sw_flux, sw_flux_file, do_read_lw_flux, lw_flux_file,&
             &h2o_lower_limit,temp_lower_limit,temp_upper_limit,co2ppmv, &
             &do_fixed_water,fixed_water,fixed_water_pres,fixed_water_lat, &
             &slowdown_rad, &
             &store_intermediate_rad, do_rad_time_avg, frierson_solar_rad, del_sol, del_sw, dt_rad, dt_rad_avg, & !mp586 added frierson_solar_rad, del_sol, del_sw for annual mean insolation
             &lonstep, do_zm_tracers, do_zm_rad, &
             &do_precip_albedo, precip_albedo_mode, precip_albedo, precip_lat,&
             &do_read_co2, co2_file, co2_variable_name, use_dyofyr, solrad, &
             &solday, equinox_day,solr_cnst

      end module rrtm_vars
!*****************************************************************************************
!*****************************************************************************************
      module rrtm_radiation
        use parkind, only : im => kind_im, rb => kind_rb
        use constants_mod,         only: pi, wtmozone, wtmh2o, gas_constant, rdgas
        implicit none
    
      contains

!*****************************************************************************************
        subroutine rrtm_radiation_init(axes,Time,ncols,nlay,lonb,latb, Time_step)
!
! Initialize diagnostics, allocate variables, set constants
!
! Modules
          use rrtm_vars
          use astronomy_mod,         only: astronomy_init
          use parrrtm, only:          nbndlw
          use parrrsw, only:          nbndsw
          use diag_manager_mod, only: register_diag_field, send_data
          use interpolator_mod, only: interpolate_type, interpolator_init, &
                                      &CONSTANT, ZERO,INTERP_WEIGHTED_P
          use fms_mod, only:          open_namelist_file, check_nml_error,  &
                                      &mpp_pe, mpp_root_pe, close_file, &
                                      &write_version_number, stdlog, &
                                      &error_mesg, NOTE, WARNING, FATAL
          use time_manager_mod, only: time_type, length_of_day, get_time
	  use transforms_mod,   only: get_grid_domain
! Local variables
          implicit none
          
          integer, intent(in), dimension(4) :: axes
          type(time_type), intent(in)       :: Time, Time_step
          integer(kind=im),intent(in)       :: ncols,nlay
          real(kind=rb),dimension(:,:),intent(in),optional :: lonb,latb !s Changed to 2d arrays as 2013 interpolator expects this. 

          integer :: i,k,seconds, time_step_seconds,res

          integer :: ierr, io, unit, is, ie, js, je

          real :: day_in_s_check


! read namelist and copy to logfile
          unit = open_namelist_file ( )
          ierr=1
          do while (ierr /= 0)
             read  (unit, nml=rrtm_radiation_nml, iostat=io, end=10)
             ierr = check_nml_error (io, 'rrtm_radiation_nml')
          enddo
10        call close_file (unit)
          
          !call write_version_number ( version, tagname )
          if ( mpp_pe() == mpp_root_pe() ) then
             write (stdlog(), nml=rrtm_radiation_nml)
          endif
          call close_file (unit)
!----
!------------ initialize diagnostic fields ---------------

          id_tdt_rad = &
               register_diag_field ( mod_name_rad, 'tdt_rad', axes(1:3), Time, &
                 'Temperature tendency due to radiation', &
                 'K/s', missing_value=missing_value               )
          id_tdt_sw = &
               register_diag_field ( mod_name_rad, 'tdt_sw', axes(1:3), Time, &
                 'Temperature tendency due to SW radiation', &
                 'K/s', missing_value=missing_value               )
          id_tdt_lw = &
               register_diag_field ( mod_name_rad, 'tdt_lw', axes(1:3), Time, &
                 'Temperature tendency due to LW radiation', &
                 'K/s', missing_value=missing_value               )
          id_coszen  = &
               register_diag_field ( mod_name_rad, 'coszen', axes(1:2), Time, &
                 'cosine of zenith angle', &
                 'none', missing_value=missing_value               )
          id_flux_sw = &
               register_diag_field ( mod_name_rad, 'flux_sw', axes(1:2), Time, &
                 'Net SW surface flux', &
                 'W/m2', missing_value=missing_value               )
          id_flux_lw = &
               register_diag_field ( mod_name_rad, 'flux_lw', axes(1:2), Time, &
                 'LW surface flux', &
                 'W/m2', missing_value=missing_value               )
	      id_olr = &
	           register_diag_field ( mod_name_rad, 'olr', axes(1:2), Time, &
	             'Outgoing LW radiation', &
	             'W/m2', missing_value=missing_value               )
	      id_toa_sw = &
	           register_diag_field ( mod_name_rad, 'toa_sw', axes(1:2), Time, &
	             'Net TOA SW flux', &
	             'W/m2', missing_value=missing_value               )
          id_albedo  = &
               register_diag_field ( mod_name_rad, 'rrtm_albedo', axes(1:2), Time, &
                 'Interactive albedo', &
                 'none', missing_value=missing_value               )
          id_ozone   = &
               register_diag_field ( mod_name_rad, 'ozone', axes(1:3), Time, &
                 'Ozone', &
                 'mmr', missing_value=missing_value               )
          id_co2   = &
               register_diag_field ( mod_name_rad, 'co2', axes(1:3), Time, &
                 'Co2', &
                 'mmr', missing_value=missing_value               )
          id_fracday  = &
               register_diag_field ( mod_name_rad, 'fracday', axes(1:2), Time, &
                 'fracday', &
                 'none', missing_value=missing_value               )
          id_half_level_temp = &
               register_diag_field ( mod_name_rad, 't_half_rrtm',(/axes(1),axes(2),axes(4)/) , Time, &
                 'Half level temperatures used by RRTM', &
                 'K', missing_value=missing_value               )                 
          id_full_level_temp = &
               register_diag_field ( mod_name_rad, 't_full_rrtm',axes(1:3) , Time, &
                 'Full level temperatures used by RRTM', &
                 'K', missing_value=missing_value               )                     
! 
!------------ make sure namelist choices are consistent -------
! this does not work at the moment, as dt_atmos from coupler_mod induces a circular dependency at compilation
!          if(dt_rad .le. dt_atmos .and. store_intermediate_rad)then
!             call error_mesg ( 'rrtm_gases_init', &
!                  ' dt_rad <= dt_atmos, for conserving memory, I am setting store_intermediate_rad=.false.', &
!                  WARNING)
!             store_intermediate_rad = .false.
!          endif
!          if(dt_rad .gt. dt_atmos .and. .not.store_intermediate_rad)then
!             call error_mesg( 'rrtm_gases_init', &
!                  ' dt_rad > dt_atmos, but store_intermediate_rad=.false. might cause time steps with zero radiative forcing!', &
!                  WARNING)
!          endif

          if(do_read_radiation .and. do_read_sw_flux .and. do_read_lw_flux) then
             if(do_read_ozone) call error_mesg( 'rrtm_gases_init', &
                  'SETTING DO_READ_OZONE TO FALSE AS DO_READ_RADIATION AND DO_READ_?W_FLUX ARE .TRUE.', NOTE)
             do_read_ozone = .false.
             if(do_read_h2o) call error_mesg( 'rrtm_gases_init', &
                  'SETTING DO_READ_H2O TO FALSE AS DO_READ_RADIATION AND DO_READ_?W_FLUX ARE .TRUE.', NOTE)
             do_read_h2o   = .false.
             if(do_precip_albedo) call error_mesg( 'rrtm_gases_init', &
                  'SETTING DO_PRECIP_ALBEDO TO FALSE AS DO_READ_RADIATION AND DO_READ_?W_FLUX ARE .TRUE.', NOTE)
             do_precip_albedo = .false.
             if(do_read_co2) call error_mesg( 'rrtm_gases_init', &
                  'SETTING DO_READ_CO2 TO FALSE AS DO_READ_RADIATION AND DO_READ_?W_FLUX ARE .TRUE.', NOTE)
             do_read_co2   = .false.
          endif

!------------ set some constants and parameters -------

          deg2rad = acos(0.)/90.

          ncols_rrt = ncols/lonstep
          nlay_rrt  = nlay

          call get_time(Time_step,time_step_seconds)

          if (dt_rad .le. 0.) then
              dt_rad = time_step_seconds !Make sure that dt_rad is set if it is not specified in the namelist
          endif

          dt_last = -real(dt_rad) !make sure we are computing radiation at the first time step



	    if (dt_rad .gt. time_step_seconds) then
	        res=mod(dt_rad, time_step_seconds)

		if (res.ne.0) then
			call error_mesg( 'rrtm_gases_init', &
        	         'dt_rad must be an integer multiple of dt_atmos', FATAL)
		endif

		day_in_s_check=length_of_day()
	        res=mod(int(day_in_s_check), dt_rad)

		if (res.ne.0) then
			call error_mesg( 'rrtm_gases_init', &
        	         'dt_rad does not fit into one day an integer number of times', WARNING)
		endif


	    endif

          if(dt_rad_avg .le. 0) dt_rad_avg = dt_rad

!------------ allocate arrays to be used later  -------

          call get_grid_domain(is, ie, js, je)
          allocate(t_half(size(lonb,1)-1,size(latb,2)-1,nlay+1)) !s changed all size(latb) to size(latb,2) as latb now 2d in 2013, where it was 1d in MiMA.

          if(.not. do_read_radiation .or. .not. do_read_sw_flux .and. .not. do_read_lw_flux)then
             allocate(h2o(ncols_rrt,nlay_rrt),o3(ncols_rrt,nlay_rrt), &
                  co2(ncols_rrt,nlay_rrt))
             allocate(ones(ncols_rrt,nlay_rrt), &
                  zeros(ncols_rrt,nlay_rrt))
             allocate(emis(ncols_rrt,nbndlw))
             allocate(taucld(nbndlw,ncols_rrt,nlay_rrt), &
                  tauaer(ncols_rrt,nlay_rrt,nbndlw))
             allocate(sw_zro(nbndsw,ncols_rrt,nlay_rrt), &
                  zro_sw(ncols_rrt,nlay_rrt,nbndsw))
             if(id_coszen > 0)allocate(zencos (size(lonb,1)-1,size(latb,2)-1))
             
             ! gases
             h2o   = 0. !this will be set by the water vapor tracer
             o3    = 0. !this will be set by an input file if do_read_ozone=.true.
             co2   = co2ppmv*1.e-6 ! convert ppmv
             zeros = 0. ! gases and clouds
             ones  = 1. ! gases and clouds

             emis  = 1. !black body: 1.0
             
             ! absorption
             taucld = 0.
             tauaer = 0.
             ! clouds
             sw_zro = 0.
             zro_sw = 0.
          endif !run RRTM?

          if(do_read_radiation)then
             call interpolator_init (rad_interp, trim(radiation_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
          endif

          if(do_read_sw_flux)then
             call interpolator_init (fsw_interp, trim(sw_flux_file)//'.nc'  , lonb, latb, data_out_of_bounds=(/CONSTANT/))
          endif

          if(do_read_lw_flux)then
             call interpolator_init (flw_interp, trim(lw_flux_file)//'.nc'  , lonb, latb, data_out_of_bounds=(/CONSTANT/))
          endif

          if(do_read_ozone)then
             call interpolator_init (o3_interp, trim(ozone_file)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
          endif

          if(do_read_h2o)then
             call interpolator_init (h2o_interp, trim(h2o_file)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
          endif

          if(do_read_co2)then
             call interpolator_init (co2_interp, trim(co2_file)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
          endif

          if(store_intermediate_rad .or. id_flux_sw > 0) &
               allocate(sw_flux(size(lonb,1)-1,size(latb,2)-1))
          if(store_intermediate_rad .or. id_flux_lw > 0) &
               allocate(lw_flux(size(lonb,1)-1,size(latb,2)-1))
	      if(id_olr > 0) &
	           allocate(olr(size(lonb,1)-1,size(latb,2)-1))
	      if(id_toa_sw > 0) &
	           allocate(toa_sw(size(lonb,1)-1,size(latb,2)-1))
          if(do_precip_albedo)allocate(rrtm_precip(size(lonb,1)-1,size(latb,2)-1))
          if(store_intermediate_rad .or. id_tdt_rad > 0)&
               allocate(tdt_rad(size(lonb,1)-1,size(latb,2)-1,nlay))
          if(id_tdt_sw > 0)allocate(tdt_sw_rad(size(lonb,1)-1,size(latb,2)-1,nlay)) 
          if(id_tdt_lw > 0)allocate(tdt_lw_rad(size(lonb,1)-1,size(latb,2)-1,nlay)) 

          if(do_precip_albedo)then
             rrtm_precip = 0.
             num_precip  = 0
          endif

	  call astronomy_init

          if(solday .gt. 0)then
             call error_mesg( mod_name_rad, &
                  ' running perpetual simulation', NOTE)
          endif


          rrtm_init=.true.

        end subroutine rrtm_radiation_init
!*****************************************************************************************
        subroutine interp_temp(z_full,z_half,t, Time)
          use rrtm_vars
          use diag_manager_mod, only: send_data
          use time_manager_mod,only:  time_type
          implicit none

          real(kind=rb),dimension(:,:,:),intent(in)  :: z_full,z_half,t
          type(time_type)               ,intent(in)          :: Time

          integer i,j,k,kend
          real dzk,dzk1,dzk2
          logical :: used

          
! note: z_full(kend) = z_half(kend), so there's something fishy
! also, for some reason, z_half(k=1)=0. so we need to deal with k=1 separately
          kend=size(z_full,3)
          do k=2,kend
             do j=1,size(t,2)
                do i=1,size(t,1)
                   dzk2 = 1./( z_full(i,j,k-1)   - z_full(i,j,k) )
                   dzk  = ( z_half(i,j,k  )   - z_full(i,j,k) )*dzk2
                   dzk1 = ( z_full(i,j,k-1)   - z_half(i,j,k) )*dzk2 
                   t_half(i,j,k) = t(i,j,k)*dzk1 + t(i,j,k-1)*dzk
                enddo
             enddo
          enddo
! top of the atmosphere: need to extrapolate. z_half(1)=0, so need to use values on full grid
          do j=1,size(t,2)
             do i=1,size(t,1)
                !standard linear extrapolation
                !top: use full points, and distance is 1.5 from k=2
                t_half(i,j,1) = 0.5*(3*t(i,j,1)-t(i,j,2))
                !bottom: z=0 => distance is -z_full(kend-1)/(z_full(kend)-z_full(kend-1))
                t_half(i,j,kend+1) = t(i,j,kend-1) &
                     + (z_half(i,j,kend+1) - z_full(i,j,kend-1))&
                     * (t     (i,j,kend  ) - t     (i,j,kend-1))&
                     / (z_full(i,j,kend  ) - z_full(i,j,kend-1))
             enddo
          enddo    


        end subroutine interp_temp
!*****************************************************************************************
!*****************************************************************************************
        subroutine run_rrtmg(is,js,Time,lat,lon,p_full,p_half,albedo,q,t,t_surf_rad,tdt,coszen,flux_sw,flux_lw)
!
! Driver for RRTMG radiation scheme.
! Prepares all inputs, calls SW and LW radiation schemes, 
!  transforms outputs back into FMS form
!
! Modules
          use fms_mod, only:         error_mesg, FATAL
          use mpp_mod, only:         mpp_pe,mpp_root_pe
          use rrtmg_lw_rad, only:    rrtmg_lw
          use rrtmg_sw_rad, only:    rrtmg_sw
          use astronomy_mod,         only: diurnal_solar
          use rrtm_vars
          use time_manager_mod,only: time_type,get_time,set_time, length_of_year, length_of_day
          use interpolator_mod,only: interpolator

          use diag_manager_mod, only: register_diag_field, send_data
          use time_manager_mod,only:  time_type
          use transforms_mod,only:    area_weighted_global_mean
!---------------------------------------------------------------------------------------------------------------
! In/Out variables
          implicit none

          integer, intent(in)                               :: is, js          ! index range for each CPU

          type(time_type),intent(in)                        :: Time            ! global time in calendar

          real(kind=rb),dimension(:,:,:),intent(in)         :: p_full,p_half   ! pressure, full and half levels
                                                                               ! dimension (lat x lon x p*)
          real(kind=rb),dimension(:,:,:),intent(in)         :: q               ! water vapor mixing ratio [g/g]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:,:),intent(in)         :: t               ! temperature [K]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:),intent(in)           :: lat,lon         ! latitude, longitude
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(in)           :: albedo          ! surface albedo
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(in)           :: t_surf_rad      ! surface temperature [K]
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:,:),intent(inout)      :: tdt             ! heating rate [K/s]
                                                                               ! dimension (lat x lon x pfull)
          real(kind=rb),dimension(:,:),intent(out)          :: coszen          ! cosine of zenith angle
                                                                               ! dimension (lat x lon)
          real(kind=rb),dimension(:,:),intent(out),optional :: flux_sw,flux_lw ! surface fluxes [W/m2]
                                                                               ! dimension (lat x lon)
                                                                               ! need to have both or none!
!---------------------------------------------------------------------------------------------------------------
! Local variables
          integer k,j,i,ij,j1,i1,ij1,kend,dyofyr,seconds,days
          integer si,sj,sk,locmin(3)
          real(kind=rb),dimension(size(q,1),size(q,2),size(q,3)) :: o3f
          real(kind=rb),dimension(size(q,1),size(q,2),size(q,3)) :: co2f,co2f_temp
          real(kind=rb),dimension(ncols_rrt,nlay_rrt) :: pfull,tfull,fracday&
               , hr,hrc, swhr, swhrc
          real(kind=rb),dimension(size(tdt,1),size(tdt,2),size(tdt,3)) :: tdt_rrtm
          real(kind=rb),dimension(ncols_rrt,nlay_rrt+1) :: uflx, dflx, uflxc, dflxc&
               ,swuflx, swdflx, swuflxc, swdflxc
          real(kind=rb),dimension(size(q,1)/lonstep,size(q,2),size(q,3)  ) :: swijk,lwijk
          real(kind=rb),dimension(size(q,1)/lonstep,size(q,2)) :: swflxijk,lwflxijk
          real(kind=rb),dimension(ncols_rrt,nlay_rrt+1):: phalf,thalf
          real(kind=rb),dimension(ncols_rrt)   :: tsrf,cosz_rr,albedo_rr
          real(kind=rb) :: dlon,dlat,dj,di 
          type(time_type) :: Time_loc
          real(kind=rb),dimension(size(q,1),size(q,2)) :: albedo_loc
          real(kind=rb),dimension(size(q,1),size(q,2),size(q,3)) :: q_tmp, h2o_vmr
          real(kind=rb),dimension(size(q,1),size(q,2)) :: fracsun
          real(kind=rb),dimension(size(q,1),size(q,2)) :: p2 !mp586 addition for annual mean insolation

	  integer :: year_in_s
          real :: r_seconds, r_days, r_total_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, dt_rad_radians, day_in_s, r_solday, r_dt_rad_avg



! debug
          integer :: indx2(2),indx(3),ii,ji,ki
          logical :: used
!---------------------------------------------------------------------------------------------------------------

          if(.not. rrtm_init)&
               call error_mesg('run_rrtm','module not initialized', FATAL)


!check if we really want to recompute radiation (alarm,input file(s))
! alarm
          call get_time(Time,seconds,days)
		  r_days = real(days)
		  r_seconds = real(seconds)
		  r_total_seconds=r_seconds+(r_days*86400.)
          if(r_total_seconds - dt_last .ge. dt_rad) then
             dt_last = r_total_seconds
          else
             if(store_intermediate_rad)then
                tdt_rrtm = tdt_rad
                flux_sw = sw_flux
                flux_lw = lw_flux
             else
                tdt_rrtm = 0.
                flux_sw  = 0.
                flux_lw  = 0.
             endif
             tdt = tdt + tdt_rrtm
             call write_diag_rrtm(Time,is,js)
             return !not time yet
          endif
!make sure we run perpetual when solday > 0)
          if(solday > 0)then
             Time_loc = set_time(seconds,solday)
          elseif(slowdown_rad .ne. 1.0)then
			 seconds = days*86400 + seconds
             Time_loc = set_time(int(seconds*slowdown_rad))
          else
             Time_loc = Time
          endif

!!!!! mp586 addition for annual mean insolation !!!!!
!!!! following https://github.com/sit23/Isca/blob/master/src/atmos_param/socrates/interface/socrates_interface.F90#L888 !!!!

       	if (frierson_solar_rad) then
            p2     = (1. - 3.*sin(lat(:,:))**2)/4.
            coszen = 0.25 * (1.0 + del_sol * p2 + del_sw * sin(lat(:,:)))
            rrsun  = 1 ! needs to be set, set to 1 so that stellar_radiation is unchanged in socrates_interface
       else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! compute zenith angle
!  this is also an output, so need to compute even if we read radiation from file
	     call get_time(Time_loc, seconds, days)
	     call get_time(length_of_year(), year_in_s)
	     day_in_s = length_of_day()
		 
	     r_seconds=real(seconds)
		 r_days=real(days)
		 r_total_seconds=r_seconds+(r_days*86400.)
		 
	     frac_of_day = r_total_seconds / day_in_s

         if(solday > 0) then
             r_solday=real(solday)
             frac_of_year=(r_solday*day_in_s)/year_in_s
	     else
	         frac_of_year = r_total_seconds / year_in_s
         endif
	     gmt = abs(mod(frac_of_day, 1.0)) * 2.0 * pi
	     time_since_ae = modulo(frac_of_year-equinox_day, 1.0) * 2.0 * pi

          if(do_rad_time_avg) then
	     r_dt_rad_avg=real(dt_rad_avg)
	     dt_rad_radians = (r_dt_rad_avg/day_in_s)*2.0*pi
	     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun,dt_rad_radians)
          else
	     ! Seasonal Cycle: Use astronomical parameters to calculate insolation
	     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun)
          end if

   		end if !mp586 addition for annual mean insolation

! input files: only deal with case where we don't need to call radiation at all
          if(do_read_radiation .and. do_read_sw_flux .and. do_read_lw_flux) then
             call interpolator( rad_interp, Time_loc, p_half, tdt_rrtm, trim(radiation_file))
             call interpolator( fsw_interp, Time_loc, flux_sw, trim(sw_flux_file))
             call interpolator( flw_interp, Time_loc, flux_lw, trim(lw_flux_file))
             ! there might be missing values due to surface topography, which would
             !  put in weird values. This is still work in progress, and cannot be 
             !  used safely!
             !if( rad_missing_value .lt. 0. )then
             !   where( tdt_rrtm .lt. rad_missing_value ) tdt_rrtm = 0.
             !   where( flux_sw  .lt. rad_missing_value ) flux_sw  = 0.
             !   where( flux_lw  .lt. rad_missing_value ) flux_lw  = 0.
             ! else
             !   where( tdt_rrtm .gt. rad_missing_value ) tdt_rrtm = 0.
             !   where( flux_sw  .gt. rad_missing_value ) flux_sw  = 0.
             !   where( flux_lw  .gt. rad_missing_value ) flux_lw  = 0.
             !endif
             tdt = tdt + tdt_rrtm
             tdt_rad = tdt_rrtm
             sw_flux = flux_sw
             lw_flux = flux_lw
             call write_diag_rrtm(Time_loc,is,js)
             return !we're done here
          endif
!---------------------------------------------------------------------------------------------
! we know now that we want to run radiation
          
          si=size(tdt,1)
          sj=size(tdt,2)
          sk=size(tdt,3)

          if(.not. use_dyofyr) dyofyr=0 !use solrad instead of day of year

          !get ozone 
          if(do_read_ozone)then
             call interpolator( o3_interp, Time_loc, p_half, o3f, trim(ozone_file))
             if (input_o3_file_is_mmr) then
                 o3f = o3f * (1000. * gas_constant / rdgas ) / wtmozone !RRTM expects all abundances to be volume mixing ratio. So if input file is mass mixing ratio, it must be converted to volume mixing ratio using the molar masses of dry air and ozone. 
                 ! Molar mass of dry air calculated from gas_constant / rdgas, and converted into g/mol from kg/mol by multiplying by 1000. This conversion is necessary because wtmozone is in g/mol.
             endif 
          endif

          !get co2
          if(do_read_co2)then
             call interpolator( co2_interp, Time, p_half, co2f, trim(co2_variable_name))
	     co2f_temp = co2f*1.e-6
             co2f = co2f_temp
          else
             co2f=co2ppmv*1.e-6
          endif

          !interactive albedo: zonal mean of precipitation
          if(do_precip_albedo .and. num_precip>0)then
             where ( abs(lat) < precip_lat*3.14159265/180. ) rrtm_precip = 0.
             do i=1,size(albedo,1)
                albedo_loc(i,:) = albedo(i,:) + (precip_albedo - albedo(i,:))&
                     &*sum(rrtm_precip,1)/size(rrtm_precip,1)/num_precip
             enddo
             rrtm_precip = 0.
             num_precip  = 0
          else
             albedo_loc = albedo
          endif
!---------------------------------------------------------------------------------------------------------------
          !Compute zonal means if that's what we want to feed to RRTM
          if(do_zm_tracers)then
             do i=1,size(q,1)
                q_tmp(i,:,:) = sum(q,1)/size(q,1)
             enddo
          else
             q_tmp = q
          endif
!---------------------------------------------------------------------------------------------------------------  
          !! water vapor stuff
          ! read water vapor
          if(do_read_h2o)then
             call interpolator( h2o_interp, Time_loc, p_half, q_tmp, trim(h2o_file))
          endif

          ! fixed water vapor
          if(do_fixed_water)then
             do j=1,size(lat,2)
                do i=1,size(lat,1)
                   if( abs(lat(i,j)) <= fixed_water_lat )then
                      where( p_full(i,j,:) <= fixed_water_pres ) q_tmp(i,j,:) = fixed_water
                   endif
                enddo
             enddo
          endif


          if(convert_sphum_to_vmr) then
              h2o_vmr = (q_tmp/(1.-q_tmp))*(1000. * gas_constant / rdgas)/wtmh2o
                 ! Convert sphum to vmr using q/1-q to get mmr, then divide by molar mass ratio.
                 ! Molar mass of dry air calculated from gas_constant / rdgas, and converted into g/mol from kg/mol by multiplying by 1000. This conversion is necessary because wtmh2o is in g/mol.              
          else
              h2o_vmr = q_tmp
          endif


!---------------------------------------------------------------------------------------------------------------
          !RRTM's first pressure level is at the surface - need to inverse order
          !also, RRTM's pressures are in hPa
          !reshape arrays
          pfull = reshape(p_full(1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))*0.01
          phalf = reshape(p_half(1:si:lonstep,:,sk+1:1:-1),(/ si*sj/lonstep,sk+1 /))*0.01
          !for RRTM, we need the top level to be greater than 0
          if(minval(phalf(:,sk+1)) .le. 0.) &
               &phalf(:,sk+1) = pfull(:,sk)*0.5
          tfull = reshape(t     (1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))
          thalf = reshape(t_half(1:si:lonstep,:,sk+1:1:-1),(/ si*sj/lonstep,sk+1 /))
          h2o   = reshape(h2o_vmr (1:si:lonstep,:,sk  :1:-1),(/ si*sj/lonstep,sk   /))
          if(do_read_ozone)o3 = reshape(o3f(1:si:lonstep,:,sk :1:-1),(/ si*sj/lonstep,sk  /))
          if(do_read_co2)co2 = reshape(co2f(1:si:lonstep,:,sk :1:-1),(/ si*sj/lonstep,sk  /))

         
          cosz_rr   = reshape(coszen    (1:si:lonstep,:),(/ si*sj/lonstep /))
          albedo_rr = reshape(albedo_loc(1:si:lonstep,:),(/ si*sj/lonstep /))
          tsrf      = reshape(t_surf_rad(1:si:lonstep,:),(/ si*sj/lonstep /))
          
!---------------------------------------------------------------------------------------------------------------
! now actually run RRTM
!
          swhr = 0.
          swdflx = 0.
          swuflx = 0.
          ! make sure we don't go beyond 'dangerous' values for radiation. Same as Merlis spectral_am2rad
          h2o   = max(h2o  , h2o_lower_limit)
          tfull = max(tfull, temp_lower_limit)
          tfull = min(tfull, temp_upper_limit)
          thalf = max(thalf, temp_lower_limit)
          thalf = min(thalf, temp_upper_limit)

          
          !h2o=h2o_lower_limit
          ! SW seems to have a problem with too small coszen values. 
          ! anything lower than 0.01 (about 15min) is set to zero
!          where(cosz_rr < 1.e-2)cosz_rr=0.
          
          if(include_secondary_gases)then
             call rrtmg_sw &
                  (ncols_rrt, nlay_rrt , icld     , iaer         , &
                  pfull     , phalf    , tfull    , thalf        , tsrf         , &
                  h2o       , o3       , co2      , ch4_val*ones , n2o_val*ones , o2_val*ones , &
                  albedo_rr , albedo_rr, albedo_rr, albedo_rr, &
                  cosz_rr   , solrad   , dyofyr   , solr_cnst, &
                  inflglw   , iceflglw , liqflglw , &
                  ! cloud parameters
                  zeros     , taucld   , sw_zro   , sw_zro   , sw_zro , &
                  zeros     , zeros    , 10*ones  , 10*ones  , &
                  tauaer    , zro_sw   , zro_sw   , zro_sw    , &
                  ! output
                  swuflx    , swdflx   , swhr     , swuflxc  , swdflxc, swhrc)
          else
             call rrtmg_sw &
                  (ncols_rrt, nlay_rrt , icld     , iaer     , &
                  pfull     , phalf    , tfull    , thalf    , tsrf , &
                  h2o       , o3       , co2      , zeros    , zeros, zeros, &
                  albedo_rr , albedo_rr, albedo_rr, albedo_rr, &
                  cosz_rr   , solrad   , dyofyr   , solr_cnst, &
                  inflglw   , iceflglw , liqflglw , &
                  ! cloud parameters
                  zeros     , taucld   , sw_zro   , sw_zro   , sw_zro , &
                  zeros     , zeros    , 10*ones  , 10*ones  , &
                  tauaer    , zro_sw   , zro_sw   , zro_sw   , &
                  ! output
                  swuflx    , swdflx   , swhr     , swuflxc  , swdflxc, swhrc)
          endif
          
          ! make sure we don't have SW radiation at night
          ! there is some optimization possible here: only feed grid points to rrtm_sw where cosz_rr>0
!          do i=1,size(swhr,2)
!             where( cosz_rr <= 0.)
!                swuflx(:,i) = 0.
!                swdflx(:,i) = 0.
!                swhr  (:,i) = 0.
!             endwhere
!          enddo
             
          swijk   = reshape(swhr(:,sk:1:-1),(/ si/lonstep,sj,sk /))*daypersec

          hr = 0.
          dflx = 0.
          uflx = 0.
          if(include_secondary_gases)then
             call rrtmg_lw &
                  (ncols_rrt     , nlay_rrt       , icld           , idrv , &
                  pfull          , phalf          , tfull          , thalf, tsrf  , &
                  h2o            , o3             , co2            , &
                  ! secondary gases
                  ch4_val*ones   , n2o_val*ones   , o2_val*ones    , &
                  cfc11_val*ones , cfc12_val*ones , cfc22_val*ones , ccl4_val*ones , &
                  ! emissivity and cloud composition
                  emis           , inflglw        , iceflglw       , liqflglw      ,  &
                  ! cloud parameters
                  zeros          , taucld         , zeros          , zeros         , 10*ones, 10*ones, &
                  tauaer         , &
                  ! output
                  uflx           , dflx           , hr             , uflxc         , dflxc  , hrc)
          else
             call rrtmg_lw &
                  (ncols_rrt, nlay_rrt, icld    , idrv , &
                  pfull     , phalf   , tfull   , thalf, tsrf , &
                  h2o       , o3      , co2     , zeros, zeros, zeros, &
                  zeros     , zeros   , zeros   , zeros, &
                  ! emissivity and cloud composition
                  emis      , inflglw , iceflglw, liqflglw, &
                  ! cloud parameters
                  zeros     , taucld  , zeros   , zeros, 10*ones, 10*ones, &
                  tauaer    , &
                  ! output
                  uflx      , dflx    , hr      , uflxc, dflxc  , hrc)
          endif

          lwijk   = reshape(hr(:,sk:1:-1),(/ si/lonstep,sj,sk /))*daypersec

!---------------------------------------------------------------------------------------------------------------
          ! get radiation
          if( do_read_radiation ) then
             call interpolator( rad_interp, Time_loc, p_half, tdt_rrtm, trim(radiation_file))
          else
! interpolate back onto GCM grid (latitude is kept the same due to parallelisation)
             dlon=1./lonstep
             do i=1,size(swijk,1)
                i1 = i+1
                ! close toroidally
                if(i1 > size(swijk,1)) i1=1
                do ij=1,lonstep
                   di = (ij-1)*dlon
                   ij1 = (i-1)*lonstep + ij
                   if(do_zm_rad) then
                      tdt_rrtm(ij1,:,:) = sum(swijk+lwijk,1)/max(1,size(swijk,1))
                   else
                      tdt_rrtm(ij1,:,:) =  &
                           +       di *(swijk(i1,:,:) + lwijk(i1,:,:)) &
                           +   (1.-di)*(swijk(i ,:,:) + lwijk(i ,:,:))
                   endif
                   if(id_tdt_sw>0)tdt_sw_rad(ij1,:,:)=di*swijk(i1,:,:)+(1.-di)*swijk(i,:,:)
                   if(id_tdt_lw>0)tdt_lw_rad(ij1,:,:)=di*lwijk(i1,:,:)+(1.-di)*lwijk(i,:,:)
                enddo
             enddo
          endif
          tdt = tdt + tdt_rrtm
          ! store radiation between radiation time steps
          if(store_intermediate_rad .or. id_tdt_rad > 0) tdt_rad = tdt_rrtm


          ! get the surface fluxes
          if(present(flux_sw).and.present(flux_lw))then
             !only surface fluxes are needed
             swflxijk = reshape(swdflx(:,1)-swuflx(:,1),(/ si/lonstep,sj /)) ! net down SW flux
             lwflxijk = reshape(  dflx(:,1)            ,(/ si/lonstep,sj /)) ! down LW flux
             dlon=1./lonstep
             do i=1,size(swijk,1)
                i1 = i+1
                ! close toroidally
                if(i1 > size(swijk,1)) i1=1
                do ij=1,lonstep
                   di = (ij-1)*dlon
                   ij1 = (i-1)*lonstep + ij
                   if(do_zm_rad) then
                      flux_sw(ij1,:) = sum(swflxijk,1)/max(1,size(swflxijk,1))
                      flux_lw(ij1,:) = sum(lwflxijk,1)/max(1,size(lwflxijk,1))
                   else
                      flux_sw(ij1,:) = di*swflxijk(i1,:) + (1.-di)*swflxijk(i ,:)
                      flux_lw(ij1,:) = di*lwflxijk(i1,:) + (1.-di)*lwflxijk(i ,:)
                   endif
                enddo
             enddo
             if ( do_read_sw_flux )then
                call interpolator( fsw_interp, Time_loc, flux_sw, trim(sw_flux_file))
             endif
             if ( do_read_lw_flux )then
                call interpolator( flw_interp, Time_loc, flux_lw, trim(lw_flux_file))
             endif
                
             ! store between radiation steps
             if(store_intermediate_rad)then
                sw_flux = flux_sw
                lw_flux = flux_lw
             else
                if(id_flux_sw > 0)sw_flux = flux_sw
                if(id_flux_lw > 0)lw_flux = flux_lw
             endif
             if(id_coszen  > 0)zencos  = coszen
          endif

          ! get the TOA fluxes (RG)
		  ! OLR:		
          if(id_olr > 0)then
             lwflxijk = reshape(  uflx(:,sk+1)-dflx(:,sk+1),(/ si/lonstep,sj /)) ! OLR
             dlon=1./lonstep
             do i=1,size(swijk,1)
                i1 = i+1
                ! close toroidally
                if(i1 > size(swijk,1)) i1=1
                do ij=1,lonstep
                   di = (ij-1)*dlon
                   ij1 = (i-1)*lonstep + ij
                   olr(ij1,:) = di*lwflxijk(i1,:) + (1.-di)*lwflxijk(i ,:)
                enddo
             enddo
          endif
		  ! TOA SW:
          if(id_toa_sw > 0)then
             swflxijk = reshape(swdflx(:,sk+1)-swuflx(:,sk+1),(/ si/lonstep,sj /)) ! net TOA SW flux, +ve down
             dlon=1./lonstep
             do i=1,size(swijk,1)
                i1 = i+1
                ! close toroidally
                if(i1 > size(swijk,1)) i1=1
                do ij=1,lonstep
                   di = (ij-1)*dlon
                   ij1 = (i-1)*lonstep + ij
                   toa_sw(ij1,:) = di*swflxijk(i1,:) + (1.-di)*swflxijk(i ,:)
                enddo
             enddo
          endif
		  
		  

          ! check if we want surface albedo as a function of precipitation
          !  call diagnostics accordingly
          if(do_precip_albedo)then
             call write_diag_rrtm(Time,is,js,o3f,co2f,fracsun,albedo_loc, t_full=t)
          else
             call write_diag_rrtm(Time,is,js,o3f,co2f,fracsun, t_full=t)
          endif
        end subroutine run_rrtmg

!*****************************************************************************************
!*****************************************************************************************
        subroutine write_diag_rrtm(Time,is,js,ozone,cotwo,fracday,albedo_loc, t_full)
! 
! write out diagnostics fields
!
! Modules
          use rrtm_vars,only:         sw_flux,lw_flux,zencos,tdt_rad,tdt_sw_rad,tdt_lw_rad,t_half,&
                                      &id_tdt_rad,id_tdt_sw,id_tdt_lw,id_coszen,&
                                      &id_flux_sw,id_flux_lw,id_albedo,id_ozone, id_co2, id_fracday,&
									  &id_olr,id_toa_sw,olr,toa_sw, id_half_level_temp, id_full_level_temp
          use diag_manager_mod, only: register_diag_field, send_data
          use time_manager_mod,only:  time_type

! Input variables
          implicit none
          type(time_type)               ,intent(in)          :: Time
          integer                       ,intent(in)          :: is, js
          real(kind=rb),dimension(:,:,:),intent(in),optional :: ozone
          real(kind=rb),dimension(:,:,:),intent(in),optional :: cotwo
          real(kind=rb),dimension(:,:  ),intent(in),optional :: albedo_loc,fracday
          real(kind=rb),dimension(:,:,:),intent(in),optional :: t_full          
! Local variables
          logical :: used

!------- temperature tendency due to radiation ------------
          if ( id_tdt_rad > 0 ) then
!             used = send_data ( id_tdt_rad, tdt_rad, Time, is, js, 1 )
             used = send_data ( id_tdt_rad, tdt_rad, Time)
          endif
!------- temperature tendency due to SW radiation ---------
          if ( id_tdt_sw > 0 ) then
!             used = send_data ( id_tdt_sw, tdt_sw_rad, Time, is, js, 1 )
             used = send_data ( id_tdt_sw, tdt_sw_rad, Time)
          endif
!------- temperature tendency due to LW radiation ---------
          if ( id_tdt_lw > 0 ) then
!             used = send_data ( id_tdt_lw, tdt_lw_rad, Time, is, js, 1 )
             used = send_data ( id_tdt_lw, tdt_lw_rad, Time)
          endif
!------- cosine of zenith angle                ------------
          if ( id_coszen > 0 ) then
             used = send_data ( id_coszen, zencos, Time)
          endif
!------- Net SW surface flux                   ------------
          if ( id_flux_sw > 0 ) then
!             used = send_data ( id_flux_sw, sw_flux, Time, is, js )
             used = send_data ( id_flux_sw, sw_flux, Time)
          endif
!------- Net LW surface flux                   ------------
          if ( id_flux_lw > 0 ) then
!             used = send_data ( id_flux_lw, lw_flux, Time, is, js )
             used = send_data ( id_flux_lw, lw_flux, Time)
          endif
!------- Net LW TOA flux                   ------------
          if ( id_olr > 0 ) then
		     used = send_data ( id_olr, olr, Time)
		  endif
!------- Net SW toa flux                   ------------
		  if ( id_toa_sw > 0 ) then
             used = send_data ( id_toa_sw, toa_sw, Time)
           endif
!------- Interactive albedo                    ------------
          if ( present(albedo_loc)) then
!             used = send_data ( id_albedo, albedo_loc, Time, is, js )
             used = send_data ( id_albedo, albedo_loc, Time)
          endif
!------- Ozone                                 ------------
          if ( present(ozone) .and. id_ozone > 0 ) then
!             used = send_data ( id_ozone, ozone, Time, is, js, 1 )
             used = send_data ( id_ozone, ozone, Time)
          endif
!------- Co2                                   ------------
          if ( present(cotwo) .and. id_co2 > 0 ) then
             used = send_data ( id_co2, cotwo, Time)
          endif

!------- Fracday                                 ------------
          if ( present(fracday) .and. id_fracday > 0 ) then
             used = send_data ( id_fracday, fracday, Time)
          endif
!------- Half-level temperatures                   ------------
          if ( id_half_level_temp > 0 ) then
             used = send_data ( id_half_level_temp, t_half , Time)
          endif          
!------- Full-level temperatures                   ------------
          if (present(t_full) .and. id_full_level_temp > 0 ) then
             used = send_data ( id_full_level_temp, t_full , Time)
          endif          

        end subroutine write_diag_rrtm
!*****************************************************************************************

        subroutine rrtm_radiation_end
          use rrtm_vars, only: do_read_ozone,o3_interp, do_read_co2, co2_interp
          use interpolator_mod, only: interpolator_end
          implicit none

          if(do_read_ozone)call interpolator_end(o3_interp)
          if(do_read_co2)call interpolator_end(co2_interp)

        end subroutine rrtm_radiation_end

!*****************************************************************************************
      end module rrtm_radiation
       


        

