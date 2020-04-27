
module qflux_mod

  use constants_mod, only : pi
  use       fms_mod, only : file_exist, open_namelist_file, check_nml_error, &
                            error_mesg, FATAL, close_file

implicit none

real ::    qflux_amp      = 30.,  &
           qflux_width    = 16.,  &
           warmpool_amp   =  5.,  &
           warmpool_width = 20.,  &
           warmpool_lond  = 20.,  &
           warmpool_latd  = 20.

integer :: warmpool_k     = 1

logical :: qflux_initialized = .false., &
           warmpool_MiMA     = .false., &
           warmpool_gaussian = .false.


namelist /qflux_nml/ qflux_amp,qflux_width,warmpool_amp,&
                     warmpool_MiMA,warmpool_width,warmpool_k,&
                     warmpool_gaussian,warmpool_lond,warmpool_latd

private 

public :: qflux_init,qflux,warmpool

contains

!########################################################
  subroutine qflux_init
    implicit none
    integer :: unit, ierr, io

    if ( file_exist('input.nml') )then
       unit = open_namelist_file()
       ierr=1; 
       do while (ierr /= 0)
          read( unit, nml=qflux_nml, iostat=io, end=10 )
          ierr = check_nml_error(io,'qflux_nml')
       enddo
10     call close_file(unit)
    endif

    qflux_initialized = .true.
    
  end subroutine qflux_init
!########################################################

  subroutine qflux(lat, flux)
! compute Q-flux as in Merlis et al (2013) [Part II]
! modified by Nathanael Zhixin Wong in Jan 2019 (code vectorized)
    implicit none
    real,dimension(:,:),intent(in)    :: lat    !latitude point
    real,dimension(:,:),intent(inout) :: flux   !total ocean heat flux
!
    integer j
    real rad_width

    if( .not. qflux_initialized ) then
       call error_mesg('qflux','qflux module not initialized',FATAL)
    endif

    rad_width = qflux_width * PI/180.
    flux = flux - qflux_amp * (1-2.*lat**2/(rad_width)**2) * &
           exp(- ((lat)**2/(rad_width)**2))/cos(lat)

  end subroutine qflux

!########################################################

  subroutine warmpool(deg_lon, deg_lat, flux)
    implicit none
    real,dimension(:,:),intent(in)    :: deg_lon,deg_lat !lon and lat
    real,dimension(:,:),intent(inout) :: flux            !total ocean heat flux
    integer i,j
    real lat,lon,rad_asym_lon,rad_asym_lat
    real,allocatable,dimension(:) :: vec_lon,vec_lat
!
! warmpool_gaussian option added by Nathanael Zhixin Wong in Jan 2019

    if ( warmpool_MiMA ) then
    
       vec_lon = deg_lon(:,1)
       vec_lat = deg_lat(1,:)

       do j=1,size(vec_lat)
          lat = vec_lat(j)/(warmpool_width*PI/180.)
          if( abs(lat) .le. 1.0 ) then
             do i=1,size(vec_lon)
                lon = vec_lon(i)
                flux(i,j) = flux(i,j) &
                   &+ (1.-lat**2.)*warmpool_amp*cos(warmpool_k*lon)
             enddo
          endif
       enddo

    endif

    if ( warmpool_gaussian ) then
        rad_asym_lon = warmpool_lond * PI/180.
        rad_asym_lat = warmpool_latd * PI/180.
        flux = flux + warmpool_amp * &
                exp(- ((deg_lat)**2/(rad_asym_lat)**2+(deg_lon-PI)**2/(rad_asym_lon)**2))
    endif

  end subroutine warmpool

!########################################################

end module qflux_mod
