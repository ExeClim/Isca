module column_grid_mod 


#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file
#endif

use         fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number, check_nml_error, stdlog
use         mpp_mod, only: mpp_error
use mpp_domains_mod, only: domain1D, mpp_get_compute_domains, mpp_get_domain_components, mpp_get_layout, mpp_global_field
use   constants_mod, only: pi
use    spec_mpp_mod, only: get_grid_domain, grid_domain, get_spec_domain

implicit none 
private 

character(len=128), parameter :: version = '$Id: column_grid.F90,v 0.1 2018/16/11 HH:MM:SS isca Exp $'
character(len=128), parameter :: tagname = '$Name: isca_201811 $'

public :: column_grid_init 
public :: get_sin_lat, area_weighted_global_mean !, got_cos_lat, getcosm_lat, get_cosm2_lat
public :: get_deg_lat, get_deg_lon, get_grid_boundaries 
public :: get_lon_max, get_lat_max!, get_longitude_origin

integer :: num_lon, num_lat, lat_max 
real    :: longitude_origin_local
logical :: module_is_initialized = .false.

real, allocatable, dimension(:) :: deg_lon, deg_lat
real, allocatable, dimension(:) :: sin_lat
real, allocatable, dimension(:) :: cos_lat
real, allocatable, dimension(:) :: cosm_lat
real, allocatable, dimension(:) :: cosm2_lat
real, allocatable, dimension(:) :: wts_lat
real, allocatable, dimension(:) :: sin_hem
real, allocatable, dimension(:) :: wts_hem
real, allocatable, dimension(:) :: lon_boundaries_global
real, allocatable, dimension(:) :: lat_boundaries_global
real :: global_sum_of_wts
real :: sum_wts 

logical :: south_to_north_local

integer :: is, ie, js, je


!! namelist parameters 

real :: lat_value = 0.0
logical :: global_average = .false.

namelist / column_grid_nml / lat_value, global_average

contains 

subroutine column_grid_init(num_lon_in, num_lat_in, longitude_origin, south_to_north)

    integer, intent(in) :: num_lon_in, num_lat_in
    real,    intent(in), optional :: longitude_origin 
    logical, intent(in), optional :: south_to_north

    real, parameter :: total_degrees = 360.
    real :: del_lat, del_lon
    integer :: i, j
    integer :: unit, ierr, io 

    if(module_is_initialized) return 

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=column_grid_nml, iostat=io)
    ierr = check_nml_error(io, 'column_grid_nml')
#else
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=column_grid_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'column_grid_nml')
    enddo
20  call close_file (unit)
#endif

    call write_version_number(version, tagname)
    if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=column_grid_nml)


    if(num_lon_in .eq. 0) then 
      call error_mesg('column_grid_init','num_lon_in cannot be zero', FATAL)
    end if 
    num_lon = num_lon_in

    if(num_lat_in .eq. 0) then 
      call error_mesg('column_grid_init','num_lat_in cannot be zero', FATAL)
    end if 
    num_lat = num_lat_in
    lat_max = num_lat 

    if(present(longitude_origin)) then 
      longitude_origin_local = longitude_origin
    else 
      longitude_origin_local = 0.0
    end if 

    call get_grid_domain(is, ie, js, je)

    allocate(deg_lon(num_lon))
    do i=1, num_lon 
      deg_lon(i) = 180*longitude_origin_local/pi + (i-1) * total_degrees / float(num_lon)
      if(deg_lon(i) .ge. total_degrees) then
        deg_lon(i) = deg_lon(i) - total_degrees
      endif
      if(deg_lon(i) .lt. 0.0) then
        deg_lon(i) = deg_lon(i) + total_degrees
      endif
    end do
    
    allocate(deg_lat(num_lat))
    allocate (sin_lat(lat_max))
    allocate (cos_lat(lat_max))
    allocate (cosm_lat(lat_max))
    allocate (cosm2_lat(lat_max))
    allocate (wts_lat(lat_max))
    if(global_average) then 
      if ((num_lat .ne. 1) .or.(num_lon .ne. 1)) then 
        call error_mesg('column_grid_init', 'cannot set global_average = True with num_lat or num_lon .ne. 1', FATAL)
      endif 
      deg_lat(num_lat) = 180*acos(pi / 4)/pi ! value that yields the fraction needed for insolation at this latitude to be equal to that of the global average (S/4)
      sin_lat = sin(pi / 180 * deg_lat )
      cos_lat = cos(pi / 180 * deg_lat )
      wts_lat = 1. ! need to set wts lat such that it tricks interpolator_mod into doing the right thing.... 
    else
      if(num_lat .eq. 1) then 
        deg_lat(num_lat) = lat_value
        sin_lat = sin(pi / 180 * deg_lat )
        cos_lat = cos(pi / 180 * deg_lat )
        wts_lat = 1. 
      else
        allocate (sin_hem(lat_max/2))
        allocate (wts_hem(lat_max/2))
        ! del_lat = 90. / (num_lat)

        ! deg_lat(1) = -90 + del_lat
        ! do i = 2, num_lat 
        !   deg_lat(i) = deg_lat(i-1) + del_lat
        ! enddo 
        if(present(south_to_north)) then 
          south_to_north_local = south_to_north
        else 
          south_to_north_local = .true.
        endif
        ! if (.not. south_to_north_local) then 
        !   deg_lat(:) = - deg_lat(:)
        ! endif

        call compute_gaussian(sin_hem, wts_hem, lat_max/2)

        if(south_to_north_local) then
          sin_lat(1:lat_max/2)   = - sin_hem
        else
          sin_lat(1:lat_max/2)   =   sin_hem
        end if

        do j=1,lat_max/2
          sin_lat(lat_max+1-j) = - sin_lat(j)
          wts_lat(j)           =   wts_hem(j)
          wts_lat(lat_max+1-j) =   wts_hem(j)
        end do

        cos_lat   = sqrt(1-sin_lat*sin_lat)
        deg_lat   = asin(sin_lat)*180.0/pi
      endif
    endif 
    cosm_lat  = 1./cos_lat
    cosm2_lat = 1./(cos_lat*cos_lat)

    ! this is done in transforms mod when spectral_dynamics is used 
    allocate( lon_boundaries_global(num_lon+1) )
    allocate( lat_boundaries_global(lat_max+1) )
    lat_boundaries_global(1) = 0.0
    if (num_lat .eq. 1) then
      lat_boundaries_global(2) = .5*pi 
    else  
      sum_wts = 0.
      do j=1,lat_max-1
        sum_wts = sum_wts + wts_lat(j)
        lat_boundaries_global(j+1) = asin(sum_wts-1.)
      end do
      lat_boundaries_global(lat_max+1) = .5*pi
      lat_boundaries_global(1) = -0.5*pi
      if (.not. south_to_north_local) then
        lat_boundaries_global(:) = -lat_boundaries_global(:)
      end if
      del_lon = 2*pi/num_lon
      do i=1,num_lon+1
        lon_boundaries_global(i) = longitude_origin_local + (i-1.5)*del_lon
      end do
    endif 

    global_sum_of_wts = sum(wts_lat)

    module_is_initialized = .true.

    return 

end subroutine column_grid_init


subroutine compute_gaussian(sin_hem_lcl, wts_hem_lcl, n_hem_lcl)
  !----------------------------------------------------------------------
  !
  !     reference:
  !       press, h. william, et. al., numerical recipes (fortran version),
  !       cambridge, england: cambridge university press (1990)
  ! 
  !------------------------------------------------------------------------
  
  integer, intent (in) :: n_hem_lcl
  real, intent (out), dimension(n_hem_lcl) :: sin_hem_lcl, wts_hem_lcl 
  
  real :: converg
  integer :: itermax
  integer :: i, iter, j, n, nprec
  real :: pp, p1, p2, p3, z, z1
  
  
  ! must use a more relaxed convergence criteria on the
  ! workstations than that for the cray T90
  ! fez code is commented out
  
  !if(kind(converg).eq.8) then
  !  converg = 1.0E-15
  !else if(kind(converg).eq.4) then
  !  converg = 1.0E-7
  !else
  !  call error_mesg('compute_gaussian','dont know what value to use for converg', FATAL)
  !end if
  
  ! The 2 lines of code below will yeild a different result than the fez code
  ! when kind(converg)=4. converg is 1.0E-6 instead of 1.0E-7
  ! This should be investigated further, but it's OK for now because it yeilds
  ! the same result on the HPCS. -- pjp
  nprec = precision(converg)
  converg = .1**nprec
  
  
  itermax = 10
  
  n=2*n_hem_lcl
  do i=1,n_hem_lcl
    z = cos(pi*(i - 0.25)/(n + 0.5))
    do iter=1,itermax
       p1 = 1.0
       p2 = 0.0
  
       do j=1,n
          p3 = p2
          p2 = p1
          p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
       end do
  
       pp = n*(z*p1 - p2)/(z*z - 1.0E+00)
       z1 = z
       z  = z1 - p1/pp
       if(ABS(z - z1) .LT. converg) go to 10
    end do
    call error_mesg('column_grid, compute_gaussian','abscissas failed to converge in itermax iterations', FATAL)
  
    10  continue

    sin_hem_lcl (i)     = z
    wts_hem_lcl (i)     = 2.0/((1.0 - z*z)*pp*pp)

  end do

end subroutine compute_gaussian


subroutine get_deg_lon(deg_lon_out)

  real, intent (out), dimension(:) :: deg_lon_out
  character(len=8) :: chtmp1, chtmp2
  
  if(.not.module_is_initialized) then
    call error_mesg('column_grid','module column_grid not initialized', FATAL)
  end if
  if(size(deg_lon_out,1).ne.num_lon) then
    write(chtmp1,'(i8)') size(deg_lon_out,1)
    write(chtmp2,'(i8)') num_lon
    call error_mesg('column_grid', &
      'size of deg_lon does not equal num_lon. size(deg_lon)='//chtmp1//' num_lon='//chtmp2, FATAL)
  end if
  
  deg_lon_out(:) = deg_lon(:)
  
  return
end subroutine get_deg_lon

subroutine get_deg_lat(deg_lat_out)
  !-----------------------------------------------------------------------
  
  real, intent (out), dimension(:) :: deg_lat_out
  
  if(.not. module_is_initialized) then
    call error_mesg('column_grid','module column_grid is not initialized', FATAL)
  end if
    
  if(size(deg_lat_out,1).eq.lat_max) then
      deg_lat_out = deg_lat
  else                            !assume grid compute domain
      deg_lat_out = deg_lat(js:je)
  end if
  
  return
end subroutine get_deg_lat

subroutine get_grid_boundaries(lon_boundaries, lat_boundaries,global)
  !-------------------------------------------------------------------------
  
    real, intent(out), dimension(:) :: lon_boundaries, lat_boundaries
    logical,intent(in),optional :: global
  
    logical :: global_tmp
    character(len=3) :: chtmp1, chtmp2
  
    if(.not.module_is_initialized) then
       call error_mesg('get_grid_boundaries','column_grid module is not initialized', FATAL)
    end if
  
    if (present(global)) then
       global_tmp = global
    else
       global_tmp = .false.
    endif
  
    if (.not. global_tmp) then
       if(size(lon_boundaries,1) /= ie-is+2) then
          write(chtmp1,'(i3)') size(lon_boundaries,1)
          write(chtmp2,'(i3)') ie-is+2
          call error_mesg('get_grid_boundaries','size(lon_boundaries) is incorrect. size(lon_boundaries)=' &
                         & //chtmp1//'  Should be'//chtmp2, FATAL)
       endif
  
       if(size(lat_boundaries,1) /= je-js+2) then
          write(chtmp1,'(i3)') size(lat_boundaries,1)
          write(chtmp2,'(i3)') je-js+2
          call error_mesg('get_grid_boundaries','size(lat_boundaries) is incorrect. size(lat_boundaries)=' &
                         & //chtmp1//'  Should be'//chtmp2, FATAL)
       endif
  
    else !global call
       if(size(lon_boundaries,1) /= num_lon+1) then
          write(chtmp1,'(i3)') size(lon_boundaries,1)
          write(chtmp2,'(i3)') num_lon+1
          call error_mesg('get_grid_boundaries','size(lon_boundaries) is incorrect. size(lon_boundaries)=' &
                         & //chtmp1//'  Should be'//chtmp2, FATAL)
       endif
  
       if(size(lat_boundaries,1) /= lat_max+1) then
          write(chtmp1,'(i3)') size(lat_boundaries,1)
          write(chtmp2,'(i3)') lat_max+1
          call error_mesg('get_grid_boundaries','size(lat_boundaries) is incorrect. size(lat_boundaries)=' &
                         & //chtmp1//'  Should be'//chtmp2, FATAL)
       endif
    endif
  
    if (global_tmp) then
       lon_boundaries = lon_boundaries_global
       lat_boundaries = lat_boundaries_global
    else
       lon_boundaries = lon_boundaries_global(is:ie+1)
       lat_boundaries = lat_boundaries_global(js:je+1)
    endif
    return
  end subroutine get_grid_boundaries

  subroutine get_sin_lat(sin_lat_out)
    !-----------------------------------------------------------------------
    
    real, intent (out), dimension(:) :: sin_lat_out
    
    if(.not. module_is_initialized) then
      call error_mesg('get_sin_lat','column_grid is not initialized', FATAL)
    end if
    
    if(size(sin_lat_out,1).eq.lat_max) then
        sin_lat_out = sin_lat
    else                            !assume grid compute domain
        sin_lat_out = sin_lat(js:je)
    end if
    
    return
  end subroutine get_sin_lat

  subroutine get_wts_lat(wts_lat_out)
    !-----------------------------------------------------------------------
    
    real, intent (out), dimension(:) :: wts_lat_out
    
    if(.not. module_is_initialized) then
      call error_mesg('get_wts_lat','column_grid is not initialized', FATAL)
    end if
      
    if(size(wts_lat_out,1).eq.lat_max) then
        wts_lat_out = wts_lat
    else                            !assume grid compute domain
        wts_lat_out = wts_lat(js:je)
    end if
    
    return
    end subroutine get_wts_lat

  function area_weighted_global_mean(field)
    !-------------------------------------------------------------------------
    real :: area_weighted_global_mean
    real, intent(in), dimension(:,:) :: field
    real, dimension(size(field,2))   :: wts_lat
    real, dimension(size(field,1), size(field,2)) :: weighted_field_local
    real, dimension(num_lon, lat_max) :: weighted_field_global
    integer :: j
    
    call get_wts_lat(wts_lat)
    do j=1,size(field,2)
      weighted_field_local(:,j) = wts_lat(j)*field(:,j)
    enddo
    
    call mpp_global_field(grid_domain, weighted_field_local, weighted_field_global)
    area_weighted_global_mean = sum(weighted_field_global)/(global_sum_of_wts*num_lon)
    
    return
    end function area_weighted_global_mean


  subroutine get_lon_max(lon_max_out)

    integer, intent (out) :: lon_max_out
      
    if(.not.module_is_initialized) then
      call error_mesg('get_lon_max','module column_grid not initialized', FATAL)
    end if
      
    lon_max_out = num_lon
      
    return
  end subroutine get_lon_max

  subroutine get_lat_max(lat_max_out)
    !-------------------------------------------------------------------------
        
    integer, intent(out) :: lat_max_out
        
    if(.not.module_is_initialized) then
      call error_mesg('get_lat_max','column_grid module is not initialized', FATAL)
    end if
        
    lat_max_out = lat_max
        
    return
  end subroutine get_lat_max

end module column_grid_mod
