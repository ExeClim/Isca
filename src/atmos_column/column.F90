module column_mod


#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use          constants_mod, only: pi
use   column_init_cond_mod, only: column_init_cond
use      field_manager_mod, only: MODEL_ATMOS
use                fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, NOTE, FATAL, write_version_number, stdlog, &
                                  close_file, open_restart_file, file_exist, set_domain,                      &
                                  read_data, write_data, check_nml_error, lowercase, uppercase, mpp_npes,     &
                                  field_size
use     tracer_manager_mod, only: get_number_tracers, query_method, get_tracer_index, NO_TRACER, get_tracer_names
use           spec_mpp_mod, only: spec_mpp_init, get_grid_domain, grid_domain
use                fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, NOTE, FATAL, write_version_number, stdlog, &
                                  close_file, open_restart_file, file_exist, set_domain,                      &
                                  read_data, write_data, check_nml_error, lowercase, uppercase, mpp_npes,     &
                                  field_size
use       time_manager_mod, only: time_type, get_time, set_time, get_calendar_type, NO_CALENDAR, &
                                  get_date, interval_alarm, operator( - ), operator( + )
use        tracer_type_mod, only: tracer_type, tracer_type_version, tracer_type_tagname


implicit none
private

public :: column_init, get_num_levels, get_surf_geopotential

character(len=128), parameter :: version = '$Id: column.F90,v 0.1 2018/14/11 HH:MM:SS isca Exp $'
character(len=128), parameter :: tagname = '$Name: isca_201811 $'

character(len=8) :: mod_name = 'column'

integer, parameter :: num_time_levels = 2
logical :: module_is_initialized = .false.
logical :: dry_model

type(time_type) :: Time_step, Alarm_time, Alarm_interval ! Used to determine when it is time to print global integrals.

real, allocatable, dimension(:,:,:,:  ) :: ug, vg, tg        ! last dimension is for time level
real, allocatable, dimension(:,:,:,:,:) :: grid_tracers      ! 4'th dimension is for time level, last dimension is for tracer number
real, allocatable, dimension(:,:      ) :: surf_geopotential
real, allocatable, dimension(:,:,  :  ) :: psg
real, allocatable, dimension(:) :: pk, bk

! for lon boundaries NTL START HERE !!! 
!real, allocatable, dimension(:) :: lat_boundaries_global, lon_boundaries_global

integer :: num_tracers, nhum
integer :: is, ie, js, je
integer :: previous, current, future 

!! NAMELIST VARIABLES 
integer :: lon_max             = 1,  & ! Column 
           lat_max             = 1,  & 
           num_fourier         = 0,  & 
           num_spherical       = 0,  & 
           num_levels          = 31  

integer, dimension(2) ::  print_interval=(/1,0/)

logical :: use_virtual_temperature = .false.

character(len=64) :: vert_coord_option      = 'even_sigma',   &
                     vert_difference_option = 'simmons_and_burridge', &
                     initial_state_option   = 'default'


real              :: scale_heights             =  4., &
                     reference_sea_level_press =  101325., &
                     surf_res                  = .1,  &
                     p_press                   = .1,  &
                     p_sigma                   = .3,  &
                     exponent                  = 2.5, &
                     ocean_topog_smoothing     = .93, &
                     initial_sphum              = 0.0

namelist /column_nml/ lon_max, lat_max, num_levels, num_fourier, print_interval, vert_coord_option, vert_difference_option, &
                      use_virtual_temperature, reference_sea_level_press, scale_heights, surf_res, p_press, p_sigma, exponent, &
                      ocean_topog_smoothing, initial_state_option, initial_sphum


contains 

subroutine column_init(Time, Time_step_in, tracer_attributes, dry_model_out, nhum_out, ocean_mask)

    type(time_type), intent(in) :: Time, Time_step_in
    type(tracer_type), intent(inout), dimension(:) :: tracer_attributes
    logical, intent(out) :: dry_model_out
    integer, intent(out) :: nhum_out
    logical, optional, intent(in), dimension(:,:) :: ocean_mask

    integer :: unit, ierr, io, ntr, nsphum, nmix_rat
    !real :: del_lon, del_lat !!! NTL START HERE
    !real :: longitude_origin_local = 0.0

    integer :: i, j 

    character(len=32) :: params
    character(len=128) :: tname, longname, units

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

    call write_version_number(version, tagname)
    if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=column_nml)
    call write_version_number(tracer_type_version, tracer_type_tagname)

    Time_step  = Time_step_in
    Alarm_interval = set_time(print_interval(2), print_interval(1))
    Alarm_time = Time + Alarm_interval


    call spec_mpp_init( num_fourier, num_spherical, lon_max, lat_max )
    
    !!! MAYBE PUT ALL OF THIS IN A FILE LIKE:
    call column_grid_init(lon_max, lat_max) ! and then get to it with other functions 
    ! this is done in transforms mod when spectral_dynamics is used 
    ! allocate( lon_boundaries_global(lon_max+1) )   ! NTL START HERE 
    ! allocate( lat_boundaries_global(lat_max+1) )
    ! lat_boundaries_global(1) = -.5 * pi 
    ! if(lat_max == 1)then 
    !   lat_boundaries_global(2) = .5 * pi
    ! else 
    !   del_lat = pi / lat_max 
    !   do j = 1,lat_max - 1
    !     lat_boundaries_global(j+1) = lat_boundaries_global(j) + del_lat
    !   enddo
    !   lat_boundaries_global(lat_max+1) = .5*pi
    ! endif 
    ! del_lon = 2*pi / lon_max
    ! do i=1,lon_max+1
    !   lon_boundaries_global(i) = longitude_origin_local + (i - 1.5)*del_lon
    ! enddo 

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
          call error_mesg('spectral_dynamics_init','sphum and mix_rat cannot both be specified as tracers at the same time', FATAL)
        endif
    endif
    dry_model_out = dry_model
    nhum_out = nhum

    call read_restart_or_do_coldstart(tracer_attributes, ocean_mask)

    module_is_initialized = .true.
    return
end subroutine column_init

















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



subroutine read_restart_or_do_coldstart(tracer_attributes, ocean_mask)

    ! For backward compatibility, this routine has the capability
    ! to read native data restart files written by inchon code.
    
    type(tracer_type), intent(inout), dimension(:) :: tracer_attributes
    logical, optional, intent(in), dimension(:,:) :: ocean_mask
    
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
                            exponent, ocean_topog_smoothing, pk, bk, ug(:,:,:,1), vg(:,:,:,1), tg(:,:,:,1), psg(:,:,1), &
                            grid_tracers(:,:,:,1,:), surf_geopotential, ocean_mask) ! NTL REMOVED LAT AND LON BOUNDARIES 
  
      ug   (:,:,:,2) = ug   (:,:,:,1)
      vg   (:,:,:,2) = vg   (:,:,:,1)
      tg   (:,:,:,2) = tg   (:,:,:,1)
      psg  (:,:,  2) = psg  (:,:,  1)
      grid_tracers(:,:,:,2,:) = grid_tracers(:,:,:,1,:)
    
      
    endif
    
    return
end subroutine read_restart_or_do_coldstart





end module column_mod