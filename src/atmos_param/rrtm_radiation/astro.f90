module rrtm_astro
!
!   Martin Jucker, 2015, https://github.com/mjucker/MiMA.
!
!   Contains all variables needed to
!   run the RRTM code, version for GCMs (hence the 'G'),
!   related to astronomy, i.e. all variables needed
!   for radiation that are not within rrtm_radiation.f90
!
!   Computes zenith angle necessary for SW radiation
!
! Modules
        use parkind, only : im => kind_im, rb => kind_rb
        use fms_mod, only : file_exist, open_namelist_file, check_nml_error, &
                            error_mesg, FATAL, close_file
! Variables
        implicit none
        logical          :: astro_initialized = .false.
!
!---------------------------------------------------------------------------------------------------------------
!                                namelist values
!---------------------------------------------------------------------------------------------------------------
        real(kind=rb)      :: obliq      = 23.439             ! Earth's obliquity
        logical            :: use_dyofyr = .false.            ! use day of the year to compute Earth-Sun distance?
                                                              !  this is done within RRTM, and assumes 365days/year!
        real(kind=rb)      :: solr_cnst= 1368.22              ! solar constant [W/m2]
        real(kind=rb)      :: solrad=1.0                      ! distance Earth-Sun [AU] if use_dyofyr=.false.
        integer(kind=im)   :: solday=0                        ! if >0, do perpetual run corresponding to 
                                                              !  day of the year = solday \in [0,days per year]
        real(kind=rb)      :: equinox_day=0.25                ! fraction of the year defining March equinox \in [0,1]
        

        namelist /astro_nml/ obliq,use_dyofyr,solr_cnst,solrad,solday,equinox_day


        contains
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
          subroutine astro_init
            implicit none
            integer :: unit, ierr, io

            if ( file_exist('input.nml') )then
               unit = open_namelist_file()
               ierr=1; 
               do while (ierr /= 0)
                  read( unit, nml=astro_nml, iostat=io, end=10 )
                  ierr = check_nml_error(io,'astro_nml')
               enddo
10             call close_file(unit)
            endif
    
            astro_initialized = .true.
            
            
          end subroutine astro_init
!--------------------------------------------------------------------------------------
! parts of this are taken from GFDL's astronomy.f90      
          subroutine compute_zenith(Time, equinox_day, dt, lat, lon, cosz, dyofyr)
!
! Computes the zenith angle for RRTM SW radiation
!
! Modules
            use time_manager_mod,only: time_type,get_time,length_of_year
            use constants_mod, only:   PI
            use fms_mod, only:         error_mesg,FATAL
! Local variables
            implicit none
! Inputs
            type(time_type),             intent(in) :: Time        ! time of year, according to calendar
            real(kind=rb),               intent(in) :: equinox_day ! fraction of year for March equinox
            integer(kind=im),            intent(in) :: dt          ! time step over which to average (if > 0)
            real(kind=rb),dimension(:,:),intent(in) :: lat,lon     ! lon/lat grid
            real(kind=rb),dimension(:,:),intent(out):: cosz        ! cosine of zenith angle
            integer(kind=im)            ,intent(out):: dyofyr      ! day of the year to compute cosz at
! Locals
            real(kind=rb),dimension(size(lat,1),size(lat,2)) :: h,cos_h, &
                 lat_h

            real dec_sin,dec_tan,dec,dec_cos,twopi,dt_pi

            integer  :: seconds,sec2,days,daysperyear
            real,dimension(size(lon,1),size(lon,2)) :: time_pi,aa,bb,tt,st,stt,sh
            real     :: radsec,radday

            integer  :: i,j
! Constants
            real     :: radpersec, radperday
            real     :: eps=1.0E-05,deg2rad 
!--------------------------------------------------------------------------------------
            deg2rad = PI/180.
            twopi = 2*PI

            if ( .not. astro_initialized ) then
               call error_mesg('astro','astro_mod not initialized',FATAL)
            endif

            call get_time(length_of_year(),sec2,daysperyear)
            if( daysperyear .ne. 365 .and. use_dyofyr ) then
             print*,' number of days per year: ',daysperyear
             call error_mesg ( 'astro', &
                  ' use_dyofyr is TRUE but the calendar year does not have 365 days. STOPPING', &
                  FATAL)
            endif

            radpersec=2*PI/86400.
            radperday=2*PI/daysperyear

            !get the time for origin
            call get_time(Time,seconds,days)
            !convert into radians
            radsec = seconds*radpersec
            dt_pi  = dt*radpersec
            
            !set local time throughout the globe
            do i=1,size(lon,1)
               !move it into interval [-PI,PI]
               time_pi(i,:) = modulo(radsec + lon(i,:),2*PI) - PI
            enddo
            where(time_pi >= PI) time_pi = time_pi - twopi
            where(time_pi < -PI) time_pi = time_pi + twopi 
            !time_pi now contains local time at each grid point
            !get day of the year relative to equinox. We set equinox at (equinox_day,equinox_day+0.5)*daysperyear
            days = days - int(equinox_day*daysperyear)
            dyofyr   = modulo(days,daysperyear) 
            !convert into radians
            radday = dyofyr*radperday
            !get declination
            dec_sin = sin(obliq*deg2rad)*sin(radday) !check sign ("-" in GFDL's code)
            dec = asin(dec_sin)
            dec_cos = cos(dec)

            !now compute the half day, to determine if it's day or night
            dec_tan = tan(dec)
            lat_h = lat
            where(lat_h == 0.5*PI) lat_h = lat - eps
            where(lat_h ==-0.5*PI) lat_h = lat + eps
            cos_h = -tan(lat_h)*dec_tan
            where(cos_h <= -1.0) h = PI
            where(cos_h >=  1.0) h = 0.0
            where(cos_h > -1.0 .and. cos_h < 1.0)&
                 h = acos(cos_h)
!---------------------------------------------------------------------
!    define terms needed in the cosine zenith angle equation.
!--------------------------------------------------------------------
            aa = sin(lat)*dec_sin
            bb = cos(lat)*dec_cos

            !finally, compute the zenith angle
            if(dt > 0)then !average over some given time interal dt
               tt  = time_pi + dt_pi
               st  = sin(time_pi)
               stt = sin(tt)
               sh  = sin(h)
               cosz = 0.0
!-------------------------------------------------------------------
!    case 1: entire averaging period is before sunrise.
!-------------------------------------------------------------------
               where (time_pi < -h .and. tt < -h) cosz = 0.0

!-------------------------------------------------------------------
!    case 2: averaging period begins before sunrise, ends after sunrise
!    but before sunset
!-------------------------------------------------------------------
               where ( (tt+h) /= 0.0 .and.   time_pi < -h .and. abs(tt) <= h)   &
                    cosz = aa + bb*(stt + sh)/ (tt + h)
!-------------------------------------------------------------------
!    case 3: averaging period begins before sunrise, ends after sunset,
!    but before the next sunrise. modify if averaging period extends 
!    past the next day's sunrise, but if averaging period is less than 
!    a half- day (pi) that circumstance will never occur.
!-------------------------------------------------------------------
               where (time_pi < -h .and. h /= 0.0 .and. h < tt)    &
                    cosz = aa + bb*( sh + sh)/(h+h)
!-------------------------------------------------------------------
!    case 4: averaging period begins after sunrise, ends before sunset.
!-------------------------------------------------------------------
               where ( abs(time_pi) <= h .and. abs(tt) <= h)    &
                    cosz = aa + bb*(stt - st)/ (tt - time_pi)
!-------------------------------------------------------------------
!    case 5: averaging period begins after sunrise, ends after sunset. 
!    modify when averaging period extends past the next day's sunrise.  
!------------------------------------------------------------------- 
               where ((h-time_pi) /= 0.0 .and. abs(time_pi) <= h .and.  h < tt)    &
                    cosz = aa + bb*(sh - st)/(h-time_pi)
!-------------------------------------------------------------------
!    case 6: averaging period begins after sunrise , ends after the
!    next day's sunrise. note that this includes the case when the
!    day length is one day (h = pi).
!-------------------------------------------------------------------
               where (twopi - h < tt .and. (tt+h-twopi) /= 0.0 .and. time_pi <= h ) &
                    cosz = (cosz*(h - time_pi) + (aa*(tt + h - twopi) +     &
                    bb*(stt + sh))) / ((h - time_pi) + (tt + h - twopi))

!-------------------------------------------------------------------
!    case 7: averaging period begins after sunset and ends before the
!    next day's sunrise
!-------------------------------------------------------------------
               where(  h <  time_pi .and. twopi - h >= tt  ) cosz = 0.0

!-------------------------------------------------------------------
!    case 8: averaging period begins after sunset and ends after the
!    next day's sunrise but before the next day's sunset. if the
!    averaging period is less than a half-day (pi) the latter
!    circumstance will never occur.
!-----------------------------------------------------------------
               where(  h <  time_pi .and. twopi - h < tt  ) &
                  cosz = aa + bb*(stt + sh) / (tt + h - twopi)

!mj-----------------------------------------------------------------
!    case 9: averaging period begins after sunset and ends after the
!    next day's sunset but before the day after's sunrise. Typically,
!    this is daily average.
!mj-----------------------------------------------------------------
               where( h < time_pi .and. h + twopi < tt .and. h /= 0. )
                  cosz = aa + bb*(sh + sh) / (h + h)
               end where
    
!----------------------------------------------------------------------
!    if instantaneous values are desired, define cosz at time t.
!----------------------------------------------------------------------
            else !version w/o time averaging
               where(abs(time_pi) <= h)
                  cosz    = aa + bb*cos(time_pi)
               elsewhere
                  cosz    = 0.0
               endwhere
            endif
!----------------------------------------------------------------------
!    be sure that cosz is not negative.
!----------------------------------------------------------------------
          cosz = max(0.0, cosz)

          end subroutine compute_zenith

        end module rrtm_astro
