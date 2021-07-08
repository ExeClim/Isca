module column_mod


#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use          constants_mod, only: rdgas, rvgas, pi, grav
use   column_init_cond_mod, only: column_init_cond
use        column_grid_mod, only: column_grid_init, get_deg_lat, get_deg_lon, get_grid_boundaries, get_sin_lat, area_weighted_global_mean
use       diag_manager_mod, only: diag_axis_init, register_diag_field, register_static_field, send_data, diag_manager_end
use      field_manager_mod, only: MODEL_ATMOS
use                fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, NOTE, FATAL, write_version_number, stdlog, &
                                  close_file, open_restart_file, file_exist, set_domain,                      &
                                  read_data, write_data, check_nml_error, lowercase, uppercase, mpp_npes,     &
                                  field_size
use        mpp_mod,         only: NULL_PE, mpp_transmit, mpp_sync, mpp_send, &
                                  mpp_broadcast, mpp_recv, mpp_max
use   press_and_geopot_mod, only: press_and_geopot_init, pressure_variables, &
                                  compute_pressures_and_heights,press_and_geopot_end
use     tracer_manager_mod, only: get_number_tracers, query_method, get_tracer_index, NO_TRACER, get_tracer_names
use           spec_mpp_mod, only: spec_mpp_init, get_grid_domain, grid_domain
use       time_manager_mod, only: time_type, get_time, set_time, get_calendar_type, NO_CALENDAR, &
                                  get_date, interval_alarm, operator( - ), operator( + )
use        tracer_type_mod, only: tracer_type, tracer_type_version, tracer_type_tagname


implicit none
private

public :: column_init, column, column_end, column_diagnostics, get_num_levels, get_surf_geopotential, get_initial_fields, get_axis_id

character(len=128), parameter :: version = '$Id: column.F90,v 0.1 2018/14/11 HH:MM:SS isca Exp $'
character(len=128), parameter :: tagname = '$Name: isca_201811 $'

integer :: id_ps, id_u, id_v, id_t
integer :: id_pres_full, id_pres_half, id_zfull, id_zhalf
integer, allocatable, dimension(:) :: id_tr
character(len=8) :: mod_name = 'column'
integer, dimension(4) :: axis_id

integer, parameter :: num_time_levels = 2
logical :: module_is_initialized = .false.
logical :: dry_model

type(time_type) :: Time_step, Alarm_time, Alarm_interval ! Used to determine when it is time to print global integrals.

real, allocatable, dimension(:) :: sin_lat
real, allocatable, dimension(:,:,:,:  ) :: ug, vg, tg        ! last dimension is for time level
real, allocatable, dimension(:,:,:,:,:) :: grid_tracers      ! 4'th dimension is for time level, last dimension is for tracer number
real, allocatable, dimension(:,:      ) :: surf_geopotential
real, allocatable, dimension(:,:,  :  ) :: psg
real, allocatable, dimension(:) :: pk, bk

! for lon boundaries NTL START HERE !!! 
!real, allocatable, dimension(:) :: lat_boundaries_global, lon_boundaries_global

real :: virtual_factor, dt_real
integer :: pe, npes, num_tracers, nhum, step_number
integer :: is, ie, js, je
integer :: previous, current, future 

!! NAMELIST VARIABLES 

logical :: use_virtual_temperature= .false., &
           graceful_shutdown      = .false.

integer :: lon_max             = 1,  & ! Column 
           lat_max             = 1,  & 
           num_fourier         = 0,  & 
           num_spherical       = 0,  & 
           num_levels          = 31, &
           num_steps           = 1  

integer, dimension(2) ::  print_interval=(/1,0/)


character(len=64) :: vert_coord_option      = 'even_sigma',   &
                     vert_difference_option = 'simmons_and_burridge', &
                     initial_state_option   = 'default'


real              :: scale_heights             =  4., &
                     reference_sea_level_press =  101325., &
                     surf_res                  = .1,  &
                     p_press                   = .1,  &
                     p_sigma                   = .3,  &
                     exponent                  = 2.5, &
                     initial_sphum             = 0.0, &
                     robert_coeff              = 0.0, &
                     raw_filter_coeff          = 1.0 

logical :: json_logging = .false.

real, dimension(2) :: valid_range_t = (/100.,500./)

namelist /column_nml/ use_virtual_temperature, valid_range_t, &
                      lon_max, lat_max, num_levels, &
                      print_interval, vert_coord_option, &
                      vert_difference_option, use_virtual_temperature, &
                      reference_sea_level_press, scale_heights, surf_res, &
                      p_press, p_sigma, exponent, &
                      initial_state_option, initial_sphum, graceful_shutdown, &
                      raw_filter_coeff, robert_coeff, json_logging


contains 

subroutine column_init(Time, Time_step_in, tracer_attributes, dry_model_out, nhum_out)

    type(time_type), intent(in) :: Time, Time_step_in
    type(tracer_type), intent(inout), dimension(:) :: tracer_attributes
    logical, intent(out) :: dry_model_out
    integer, intent(out) :: nhum_out

    integer :: unit, ierr, io, ntr, nsphum, nmix_rat, seconds, days 
    !real :: del_lon, del_lat !!! NTL START HERE
    !real :: longitude_origin_local = 0.0

    integer :: i, j 

    character(len=32) :: params
    character(len=128) :: tname, longname, units
    character(len=8) :: err_msg_1

    if(module_is_initialized) return

#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=column_nml, iostat=io)
        ierr = check_nml_error(io, 'column_nml')
#else
        unit = open_namelist_file()
        ierr=1
        do while (ierr /= 0)
            read(unit, nml=column_nml, iostat=io, end=20)
            ierr = check_nml_error (io, 'column_nml')
        enddo
    20  call close_file (unit)
#endif

    pe   = mpp_pe()
    npes = mpp_npes()
    if (npes .gt. 1) then 
      write(err_msg_1,'(i8)') npes
      call error_mesg('column_init','Can only run column model on one processor but npes = '//err_msg_1, FATAL)
    endif 

    call write_version_number(version, tagname)
    if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=column_nml)
    call write_version_number(tracer_type_version, tracer_type_tagname)

    Time_step  = Time_step_in
    Alarm_interval = set_time(print_interval(2), print_interval(1))
    Alarm_time = Time + Alarm_interval


    call spec_mpp_init( num_fourier, num_spherical, lon_max, lat_max )
    
    !!! MAYBE PUT ALL OF THIS IN A FILE LIKE:
    call column_grid_init(lon_max, lat_max) ! and then get to it with other functions 
    call get_grid_domain(is, ie, js, je)
    call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
    call allocate_fields

    do ntr=1,num_tracers

        call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
        tracer_attributes(ntr)%name = lowercase(tname)
      
    enddo
    nsphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
    nmix_rat = get_tracer_index(MODEL_ATMOS, 'mix_rat')
      
    if(nsphum == NO_TRACER) then
        if(nmix_rat == NO_TRACER) then
          nhum = 0
          dry_model = .true.
        else
          nhum = nmix_rat
          dry_model = .false.
        endif
    else
        if(nmix_rat == NO_TRACER) then
          nhum = nsphum
          dry_model = .false.
        else
          call error_mesg('column_init','sphum and mix_rat cannot both be specified as tracers at the same time', FATAL)
        endif
    endif
    dry_model_out = dry_model
    nhum_out = nhum


    call read_restart_or_do_coldstart(tracer_attributes)

    call press_and_geopot_init(pk, bk, use_virtual_temperature, vert_difference_option)
    call column_diagnostics_init(Time)

    if(use_virtual_temperature) then
      virtual_factor = (rvgas/rdgas) - 1.0
    end if

    allocate(sin_lat(js:je))
    call get_sin_lat(sin_lat)

    call set_domain(grid_domain)
    call get_time(Time_step, seconds, days)
    dt_real = 86400*days + seconds

    module_is_initialized = .true.
    ! NTL: CHECK AGAINST spectral_dynamics_init TO SEE WHAT ELSE NEEDS TO BE INITIALISED 
    return
end subroutine column_init

subroutine column(Time, psg_final, ug_final, vg_final, tg_final, tracer_attributes, grid_tracers_final, &
  time_level_out, dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers, wg_full, p_full, p_half, z_full)

type(time_type),  intent(in) :: Time
real, intent(out), dimension(is:, js:      ) :: psg_final
real, intent(out), dimension(is:, js:, :   ) :: ug_final, vg_final, tg_final
real, intent(out), dimension(is:, js:, :,:,:) :: grid_tracers_final
type(tracer_type),intent(inout), dimension(:) :: tracer_attributes
integer, intent(in)                           :: time_level_out

real, intent(inout), dimension(is:, js:      ) :: dt_psg
real, intent(inout), dimension(is:, js:, :   ) :: dt_ug, dt_vg, dt_tg
real, intent(inout), dimension(is:, js:, :, :) :: dt_tracers
real, intent(out),   dimension(is:, js:, :   ) :: wg_full, p_full
real, intent(out),   dimension(is:, js:, :   ) :: p_half
real, intent(in),    dimension(is:, js:, :   ) :: z_full

type(time_type) :: Time_diag

real, dimension(is:ie, js:je, num_levels, num_tracers) :: dt_tracers_tmp
integer :: p, seconds, days
real :: delta_t
real    :: extrtmp
integer :: ii,jj,kk,i1,j1,k1
integer :: ntr

logical :: pe_is_valid = .true.
logical :: r_pe_is_valid = .true.

! THIS IS WHERE I NEED TO START FROM... 
! TO DO: LOCAL VARIABLES, CHECK INPUTS, HOOK UP TO ATMOSPHERE.F90 AND GLOBAL VARIABLE DEFINTIONS 
! DO something simple like next = current + dt_var * timestep 

if(.not.module_is_initialized) then
  call error_mesg('column','column has not been initialized ', FATAL)
endif

dt_tracers_tmp = dt_tracers

step_loop: do step_number=1,num_steps

if(previous == current) then
  delta_t = dt_real/num_steps
else
  delta_t = 2*dt_real/num_steps
endif
if(num_time_levels == 2) then
  future = 3 - current
else
  call error_mesg('column','Do not know how to set time pointers when num_time_levels does not equal 2',FATAL)
endif


call leapfrog_3d_real(tg, dt_tg, previous, current, future, delta_t, robert_coeff, raw_filter_coeff)




if(minval(tg(:,:,:,future)) < valid_range_t(1) .or. maxval(tg(:,:,:,future)) > valid_range_t(2)) then
  pe_is_valid = .false.
!mj !s This doesn't affect the normal running of the code in any way. It simply allows identification of the point where temp violation has occured.
   if(minval(tg(:,:,:,future)) < valid_range_t(1))then
      extrtmp = minval(tg(:,:,:,future))
   else
      extrtmp = maxval(tg(:,:,:,future))
   endif
   do k1=1,size(tg,3)
      do j1=1,size(tg,2)
         do i1=1,size(tg,1)
            if(tg(i1,j1,k1,future) .eq. extrtmp)then
               ii=i1
               jj=j1
               kk=k1
               exit
            endif
         enddo
      enddo
   enddo
   write(*,'(a,i3,a,3i3,2f10.3)')'PE, location, Textr(curr,future): ',mpp_pe()&
        &,': ',ii,jj,kk&
        &,tg(ii,jj,kk,current)&
        &,tg(ii,jj,kk,future)
   write(*,'(a,i3,a,3i3,2f10.3)')'PE, location, Uextr(curr,future): ',mpp_pe()&
        &,': ',ii,jj,kk&
        &,ug(ii,jj,kk,current)&
        &,ug(ii,jj,kk,future)
!jm
  if (.not.graceful_shutdown) then
    call error_mesg('column','temperatures out of valid range', FATAL)
  endif
endif

! synchronisation between nodes.  THIS WILL SLOW DOWN THE RUN but ensures
! all partially complete diagnostics are written to netcdf output.
if (graceful_shutdown) then
  if (pe == mpp_root_pe()) then
    do p = 0, npes-1
      ! wait for all nodes to report they are error free
      if (p.ne.pe) then
        call mpp_recv(r_pe_is_valid, p)
        if (.not.r_pe_is_valid .or. .not.pe_is_valid) then
          ! node reports error, tell all the others to shutdown
          r_pe_is_valid = .false.
          exit
        end if
      end if
    end do
    ! tell all the nodes to continue, or not
    do p =0, npes-1
      if (p.ne.pe) call mpp_send(r_pe_is_valid, p)
    end do
  else
    call mpp_send(pe_is_valid, mpp_root_pe())
    ! wait to hear back from root that all are ok to continue
    call mpp_recv(r_pe_is_valid, mpp_root_pe())
  endif
  if (.not.r_pe_is_valid) then
    ! one of the nodes has broken the condition.  gracefully shutdown diagnostics
    ! and then raise a fatal error after all have hit the sync point.
    call diag_manager_end(Time)
    call mpp_sync()
    call error_mesg('column','temperatures out of valid range', FATAL)
  endif
endif

! NTL: Need to write a different version of this which uses the leapfrog below...
do ntr = 1, num_tracers  
  call leapfrog_3d_real(grid_tracers(:,:,:,:,ntr),dt_tracers_tmp(:,:,:,ntr),previous,current,future,delta_t,tracer_attributes(ntr)%robert_coeff, raw_filter_coeff)
enddo 


previous = current
current  = future

call get_time(Time, seconds, days)
seconds = seconds + step_number*int(dt_real/2)
Time_diag = set_time(seconds, days)

enddo step_loop

psg_final = psg(:,:,  previous)
ug_final  =  ug(:,:,:,previous)
vg_final  =  vg(:,:,:,previous)
tg_final  =  tg(:,:,:,current)
grid_tracers_final(:,:,:,time_level_out,:) = grid_tracers(:,:,:,current,:)

return 
end subroutine column 

subroutine column_diagnostics_init(Time)

  type(time_type), intent(in) :: Time
  real, dimension(lon_max  ) :: lon
  real, dimension(lon_max+1) :: lonb
  real, dimension(lat_max  ) :: lat
  real, dimension(lat_max+1) :: latb
  real, dimension(num_levels)   :: p_full, ln_p_full
  real, dimension(num_levels+1) :: p_half, ln_p_half
  integer, dimension(3) :: axes_3d_half, axes_3d_full
  integer :: id_lonb, id_latb, id_phalf, id_lon, id_lat, id_pfull
  integer :: id_pk, id_bk, id_zsurf, ntr
  real :: rad_to_deg
  logical :: used
  real,dimension(2) :: vrange
  character(len=128) :: tname, longname, units

  ! NTL: NEED TO DO THIS 

  vrange = (/ -400., 400. /)

  rad_to_deg = 180./pi
  call get_grid_boundaries(lonb,latb,global=.true.)
  call get_deg_lon(lon)
  call get_deg_lat(lat)

  id_lonb=diag_axis_init('lonb', rad_to_deg*lonb, 'degrees_E', 'x', 'longitude edges', set_name=mod_name, Domain2=grid_domain)
  id_latb=diag_axis_init('latb', rad_to_deg*latb, 'degrees_N', 'y', 'latitude edges',  set_name=mod_name, Domain2=grid_domain)
  id_lon =diag_axis_init('lon', lon, 'degrees_E', 'x', 'longitude', set_name=mod_name, Domain2=grid_domain, edges=id_lonb)
  id_lat =diag_axis_init('lat', lat, 'degrees_N', 'y', 'latitude',  set_name=mod_name, Domain2=grid_domain, edges=id_latb)
  
  call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, reference_sea_level_press)
  p_half = .01*p_half
  p_full = .01*p_full
  id_phalf = diag_axis_init('phalf',p_half,'hPa','z','approx half pressure level',direction=-1,set_name=mod_name)
  id_pfull = diag_axis_init('pfull',p_full,'hPa','z','approx full pressure level',direction=-1,set_name=mod_name,edges=id_phalf)
  
  axes_3d_half = (/ id_lon, id_lat, id_phalf /)
  axes_3d_full = (/ id_lon, id_lat, id_pfull /)
  axis_id(1) = id_lon
  axis_id(2) = id_lat
  axis_id(3) = id_pfull
  axis_id(4) = id_phalf

  id_pk = register_static_field(mod_name, 'pk', (/id_phalf/), 'vertical coordinate pressure values', 'pascals')
  id_bk = register_static_field(mod_name, 'bk', (/id_phalf/), 'vertical coordinate sigma values', 'none')
  id_zsurf = register_static_field(mod_name, 'zsurf', (/id_lon,id_lat/), 'geopotential height at the surface', 'm')
  
  if(id_pk    > 0) used = send_data(id_pk, pk, Time)
  if(id_bk    > 0) used = send_data(id_bk, bk, Time)
  if(id_zsurf > 0) used = send_data(id_zsurf, surf_geopotential/grav, Time)

  id_ps  = register_diag_field(mod_name, &
        'ps', (/id_lon,id_lat/),       Time, 'surface pressure',             'pascals')
  
  id_u   = register_diag_field(mod_name, &
        'ucomp',   axes_3d_full,       Time, 'zonal wind component',         'm/sec',      range=vrange)
  
  id_v   = register_diag_field(mod_name, &
        'vcomp',   axes_3d_full,       Time, 'meridional wind component',    'm/sec',      range=vrange)

  id_t   = register_diag_field(mod_name, &
        'temp',    axes_3d_full,       Time, 'temperature',                  'deg_k',      range=valid_range_t)

  id_pres_full = register_diag_field(mod_name, &
        'pres_full',    axes_3d_full,       Time, 'pressure at full model levels', 'pascals')

  id_pres_half = register_diag_field(mod_name, &
        'pres_half',    axes_3d_half,       Time, 'pressure at half model levels', 'pascals')

  id_zfull   = register_diag_field(mod_name, &
        'height',  axes_3d_full,       Time, 'geopotential height at full model levels','m')

  id_zhalf   = register_diag_field(mod_name, &
        'height_half',  axes_3d_half,  Time, 'geopotential height at half model levels','m')

  allocate(id_tr(num_tracers))
  do ntr=1,num_tracers
    call get_tracer_names(MODEL_ATMOS, ntr, tname, longname, units)
    id_tr(ntr) = register_diag_field(mod_name, tname, axes_3d_full, Time, longname, units)
  enddo

  return
end subroutine column_diagnostics_init




subroutine column_diagnostics(Time, p_surf, u_grid, v_grid, t_grid, wg_full, tr_grid, time_level)

  type(time_type), intent(in) :: Time
  real, intent(in), dimension(is:, js:)          :: p_surf
  real, intent(in), dimension(is:, js:, :)       :: u_grid, v_grid, t_grid, wg_full
  real, intent(in), dimension(is:, js:, :, :, :) :: tr_grid
  integer, intent(in) :: time_level
  
  real, dimension(is:ie, js:je, num_levels)    :: ln_p_full, p_full, z_full 
  real, dimension(is:ie, js:je, num_levels+1)  :: ln_p_half, p_half, z_half
  logical :: used
  integer :: ntr, i, j, k
  character(len=8) :: err_msg_1, err_msg_2
  
  if(id_ps  > 0)    used = send_data(id_ps,  p_surf, Time)
  if(id_u   > 0)    used = send_data(id_u,   u_grid, Time)
  if(id_v   > 0)    used = send_data(id_v,   v_grid, Time)
  if(id_t   > 0)    used = send_data(id_t,   t_grid, Time)
  
  if(id_zfull > 0 .or. id_zhalf > 0) then
    call compute_pressures_and_heights(t_grid, p_surf, surf_geopotential, z_full, z_half, p_full, p_half)
  else if(id_pres_half > 0 .or. id_pres_full > 0) then
    call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, p_surf)
  endif
  
  if(id_zfull > 0)   used = send_data(id_zfull,      z_full, Time)
  if(id_zhalf > 0)   used = send_data(id_zhalf,      z_half, Time)
  if(id_pres_full>0) used = send_data(id_pres_full,  p_full, Time)
  if(id_pres_half>0) used = send_data(id_pres_half,  p_half, Time)
  
  if(size(tr_grid,5) /= num_tracers) then
    write(err_msg_1,'(i8)') size(tr_grid,5)
    write(err_msg_2,'(i8)') num_tracers
    call error_mesg('column_diagnostics','size(tracers)='//err_msg_1//' Should be='//err_msg_2, FATAL)
  endif
  do ntr=1,num_tracers
    if(id_tr(ntr) > 0) used = send_data(id_tr(ntr), tr_grid(:,:,:,time_level,ntr), Time)
  enddo
  
  
  if(interval_alarm(Time, Time_step, Alarm_time, Alarm_interval)) then
    call global_integrals(Time, p_surf, u_grid, v_grid, t_grid, wg_full, tr_grid(:,:,:,time_level,:))
  endif
  
  return
  end subroutine column_diagnostics

  subroutine column_end(tracer_attributes, Time)

    type(tracer_type), intent(in), dimension(:) :: tracer_attributes
    type(time_type), intent(in), optional :: Time
    integer :: ntr, nt
    character(len=64) :: file, tr_name
    
    if(.not.module_is_initialized) return
    
    file='RESTART/column_model.res'
    call write_data(trim(file), 'previous', previous, no_domain=.true.)
    call write_data(trim(file), 'current',  current,  no_domain=.true.)
    call write_data(trim(file), 'pk', pk, no_domain=.true.)
    call write_data(trim(file), 'bk', bk, no_domain=.true.)
    do nt=1,num_time_levels
      call write_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain)
      call write_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain)
      call write_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain)
      call write_data(trim(file), 'psg', psg(:,:,  nt), grid_domain)
      do ntr = 1,num_tracers
        tr_name = trim(tracer_attributes(ntr)%name)
        call write_data(trim(file), trim(tr_name), grid_tracers(:,:,:,nt,ntr), grid_domain)
      enddo
    enddo
    call write_data(trim(file), 'surf_geopotential', surf_geopotential, grid_domain)
    
    deallocate(ug, vg, tg, psg)
    deallocate(sin_lat)
    deallocate(pk, bk)
    deallocate(surf_geopotential)
    deallocate(grid_tracers)
    
    call column_diagnostics_end
    call press_and_geopot_end
    call set_domain(grid_domain)
    module_is_initialized = .false.
    
    return
  end subroutine column_end

  subroutine column_diagnostics_end

    if(.not.module_is_initialized) return
    
    deallocate(id_tr)
    
    return
  end subroutine column_diagnostics_end  

  subroutine global_integrals(Time, p_surf, u_grid, v_grid, t_grid, wg_full, tr_grid)
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(is:ie, js:je)                          :: p_surf
    real, intent(in), dimension(is:ie, js:je, num_levels)              :: u_grid, v_grid, t_grid, wg_full
    real, intent(in), dimension(is:ie, js:je, num_levels, num_tracers) :: tr_grid
    integer :: year, month, days, hours, minutes, seconds
    character(len=4), dimension(12) :: month_name
    
    real, dimension(is:ie, js:je, num_levels) :: speed
    real :: max_speed, avgT
    
    month_name=(/' Jan',' Feb',' Mar',' Apr',' May',' Jun',' Jul',' Aug',' Sep',' Oct',' Nov',' Dec'/)
    
    speed = sqrt(u_grid*u_grid + v_grid*v_grid)
    max_speed = maxval(speed)
    call mpp_max(max_speed)
    
    avgT = area_weighted_global_mean(t_grid(:,:, num_levels))
    
    if(mpp_pe() == mpp_root_pe()) then
      if(get_calendar_type() == NO_CALENDAR) then
        call get_time(Time, seconds, days)
        if (json_logging) then
          write(*, 300) days, seconds, max_speed, avgT
        else
          write(*,100) days, seconds
        end if
      else
        call get_date(Time, year, month, days, hours, minutes, seconds)
        if (json_logging) then
          write(*,400) year, month, days, hours, minutes, seconds, max_speed, avgT
        else
          write(*,200) year, month_name(month), days, hours, minutes, seconds
        end if
      endif
    endif
    100 format(' Integration completed through',i6,' days',i6,' seconds')
    200 format(' Integration completed through',i5,a4,i3,2x,i2,':',i2,':',i2)
    300 format(1x, '{"day":',i6,2x,',"second":', i6, &
        2x,',"max_speed":',e13.6,3x,',"avg_T":',e13.6, 3x '}')
    400 format(1x, '{"date": "',i0.4,'-',i0.2,'-',i0.2, &
      '", "time": "', i0.2,':', i0.2,':', i0.2, '", "max_speed":',f6.1,3x,',"avg_T":',f6.1, 3x '}')
    
  end subroutine global_integrals




subroutine allocate_fields

    allocate (psg    (is:ie, js:je,             num_time_levels))
    allocate (ug     (is:ie, js:je, num_levels, num_time_levels))
    allocate (vg     (is:ie, js:je, num_levels, num_time_levels))
    allocate (tg     (is:ie, js:je, num_levels, num_time_levels))
    
    allocate (pk(num_levels+1), bk(num_levels+1))
    
    allocate (surf_geopotential(is:ie, js:je))
    
    allocate (grid_tracers(is:ie, js:je, num_levels, num_time_levels, num_tracers))
    
    ! Filling allocatable arrays with zeros immediately after allocation facilitates code debugging
    psg=0.; ug=0.; vg=0.; tg=0.
    pk=0.; bk=0.; surf_geopotential=0.; grid_tracers=0.
    
    return
end subroutine allocate_fields

subroutine get_num_levels(num_levels_out)
    integer, intent(out) :: num_levels_out
    
    if(.not.module_is_initialized) then
      call error_mesg('get_num_levels', 'column_init has not been called.', FATAL)
    endif
    
    num_levels_out = num_levels
    
    return
end subroutine get_num_levels

subroutine get_surf_geopotential(surf_geopotential_out)
    real, intent(out), dimension(:,:) :: surf_geopotential_out
    character(len=64) :: chtmp='shape(surf_geopotential)=              should be                '
    
    if(.not.module_is_initialized) then
      call error_mesg('get_surf_geopotential', 'column_init has not been called.', FATAL)
    endif
    
    if(any(shape(surf_geopotential_out) /= shape(surf_geopotential))) then
      write(chtmp(26:37),'(3i4)') shape(surf_geopotential_out)
      write(chtmp(50:61),'(3i4)') shape(surf_geopotential)
      call error_mesg('get_surf_geopotential', 'surf_geopotential has wrong shape. '//chtmp, FATAL)
    endif
    
    surf_geopotential_out = surf_geopotential

    return
end subroutine get_surf_geopotential



subroutine read_restart_or_do_coldstart(tracer_attributes)

    ! For backward compatibility, this routine has the capability
    ! to read native data restart files written by inchon code.
    
    type(tracer_type), intent(inout), dimension(:) :: tracer_attributes
    
    integer :: m, n, k, nt, ntr
    integer, dimension(4) :: siz
    character(len=64) :: file, tr_name
    character(len=4) :: ch1,ch2,ch3,ch4,ch5,ch6
    
    file = 'INPUT/column_model.res.nc'
    if(file_exist(trim(file))) then
      call field_size(trim(file), 'ug', siz)
      if(lon_max /= siz(1) .or. lat_max /= siz(2)) then
        write(ch1,'(i4)') siz(1)
        write(ch2,'(i4)') siz(2)
        write(ch3,'(i4)') lon_max
        write(ch4,'(i4)') lat_max
        call error_mesg('column_init','Resolution of restart data does not match resolution specified on namelist.'// &
        ' Restart data: lon_max='//ch1//', lat_max='//ch2//'  Namelist: lon_max='//ch3//', lat_max='//ch4, FATAL)
      endif
      call read_data(trim(file), 'previous', previous, no_domain=.true.)
      call read_data(trim(file), 'current',  current,  no_domain=.true.)
      call read_data(trim(file), 'pk', pk, no_domain=.true.)
      call read_data(trim(file), 'bk', bk, no_domain=.true.)
      do nt=1,num_time_levels
        call read_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain, timelevel=nt)
        call read_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain, timelevel=nt)
        call read_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain, timelevel=nt)
        call read_data(trim(file), 'psg', psg(:,:,  nt), grid_domain, timelevel=nt)
        do ntr = 1,num_tracers
          tr_name = trim(tracer_attributes(ntr)%name)
          call read_data(trim(file), trim(tr_name), grid_tracers(:,:,:,nt,ntr), grid_domain, timelevel=nt)
        enddo ! loop over tracers
      enddo ! loop over time levels
      call read_data(trim(file), 'surf_geopotential', surf_geopotential, grid_domain)

    else
      do ntr = 1,num_tracers
        if(trim(tracer_attributes(ntr)%name) == 'sphum') then
          grid_tracers(:,:,:,:,ntr) = initial_sphum
        else if(trim(tracer_attributes(ntr)%name) == 'mix_rat') then
          grid_tracers(:,:,:,:,ntr) = 0.
        else
          grid_tracers(:,:,:,:,ntr) = 0.
        endif
      enddo
    
      previous = 1
      current  = 1
      call column_init_cond(initial_state_option, tracer_attributes, reference_sea_level_press, use_virtual_temperature,&
                            vert_coord_option, vert_difference_option, scale_heights, surf_res, p_press, p_sigma,  &
                            exponent, pk, bk, ug(:,:,:,1), vg(:,:,:,1), tg(:,:,:,1), psg(:,:,1), &
                            grid_tracers(:,:,:,1,:), surf_geopotential) ! NTL REMOVED LAT AND LON BOUNDARIES 
  
      ug   (:,:,:,2) = ug   (:,:,:,1)
      vg   (:,:,:,2) = vg   (:,:,:,1)
      tg   (:,:,:,2) = tg   (:,:,:,1)
      psg  (:,:,  2) = psg  (:,:,  1)
      grid_tracers(:,:,:,2,:) = grid_tracers(:,:,:,1,:)
    
      
    endif
    
    return
end subroutine read_restart_or_do_coldstart

subroutine get_initial_fields(ug_out, vg_out, tg_out, psg_out, grid_tracers_out)
  real, intent(out), dimension(:,:,:)   :: ug_out, vg_out, tg_out
  real, intent(out), dimension(:,:)     :: psg_out
  real, intent(out), dimension(:,:,:,:) :: grid_tracers_out
  
  if(.not.module_is_initialized) then
    call error_mesg('column, get_initial_fields','column has not been initialized',FATAL)
  endif
  
  if(previous /= 1 .or. current /= 1) then
    call error_mesg('column, get_initial_fields','This routine may be called only to get the&
                    & initial values after a cold_start',FATAL)
  endif
  
  ug_out  =  ug(:,:,:,1)
  vg_out  =  vg(:,:,:,1)
  tg_out  =  tg(:,:,:,1)
  psg_out = psg(:,:,  1)
  grid_tracers_out = grid_tracers(:,:,:,1,:)
  
  end subroutine get_initial_fields

  function get_axis_id()
    integer, dimension(4) :: get_axis_id
    
    if(.not.module_is_initialized) then
      call error_mesg('get_axis_id','column_diagnostics_init has not been called.', FATAL)
    endif
    get_axis_id = axis_id
    return
  end function get_axis_id

  


  subroutine leapfrog_3d_real(a, dt_a, previous, current, future, delta_t, robert_coeff, raw_filter_coeff)

    real, intent(inout), dimension(:,:,:,:) :: a
    real, intent(in),    dimension(:,:,:  ) :: dt_a
    integer, intent(in) :: previous, current, future
    real,    intent(in) :: delta_t, robert_coeff, raw_filter_coeff
    
    real, dimension(size(dt_a,1),size(dt_a,2),size(dt_a,3)) :: prev_curr_part_raw_filter
  
    
    prev_curr_part_raw_filter=a(:,:,:,previous) - 2.0*a(:,:,:,current) !st Defined at the start to get unmodified value of a(:,:,:,current).
    
    if(previous == current) then
      a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
      a(:,:,:,current) = a(:,:,:,current ) + robert_coeff * (prev_curr_part_raw_filter + a(:,:,:,future ))*raw_filter_coeff
    else
      a(:,:,:,current) = a(:,:,:,current ) + robert_coeff * (prev_curr_part_raw_filter                   )*raw_filter_coeff
      a(:,:,:,future ) = a(:,:,:,previous) + delta_t * dt_a
      a(:,:,:,current) = a(:,:,:,current ) + robert_coeff * a(:,:,:,future)*raw_filter_coeff
    endif
    
    a(:,:,:,future ) = a(:,:,:,future ) + robert_coeff * (prev_curr_part_raw_filter + a(:,:,:,future )) * (raw_filter_coeff-1.0) 
    
    !st RAW filter (see e.g. Williams 2011 10.1175/2010MWR3601.1) conserves 3-time-level mean in leap-frog integrations, improving amplitude accuracy of leap-frog scheme from first to third order).
    
    return
  end subroutine leapfrog_3d_real 
  



end module column_mod