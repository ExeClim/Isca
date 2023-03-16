module astronomy_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    astronomy_mod provides astronomical variables for use
!    by other modules within fms. the only currently used interface is
!    for determination of astronomical values needed by the shortwave
!    radiation packages.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

!  shared modules:

use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, NOTE, WARNING, close_file
use time_manager_mod,  only: time_type, set_time, get_time, &
                             get_date_julian, set_date_julian, &
                             set_date, length_of_year, length_of_day, &
                             time_manager_init, &
                             operator(-), operator(+), &
                             operator( // ), operator(<), get_calendar_type, NO_CALENDAR
use constants_mod,     only: constants_init, PI, omega, orbital_rate, radius
use mpp_mod,           only: input_nml_file

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    astronomy_mod provides astronomical variables for use
!    by other modules within fms. the only currently used interface is
!    for determination of astronomical values needed by the shortwave
!    radiation packages.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: astronomy.F90,v 17.0.10.1 2010/08/31 14:21:37 z1l Exp $'
character(len=128)  :: tagname =  '$Name: testing $'

! Version Details
! [2016/11/16] <Stephen Thomson>:  Modified v 17.0.10.1 by changing the time averaging of coszen.
!
!   - Problem with original v 17.0.10.1: when `diurnal_solar' was called with the time-averaging option:
!
!     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun, dt)
!
!     diurnal_solar claimed to return coszen time-averaged between the current time (t) and t+dt.
!     However, it actually returned coszen time-averaged over the _hours of daylight_ between t and t+dt.
!
!   - This inconsistency could have been corrected by multiplying coszen by the optional array 'fracday' (an array containing
!     the fraction of daylight experienced in each gridcell), but this was not explained.
!
!   - Therefore, if diurnal_solar was used without multiplying by fracday, it returned answers that were dependent on dt,
!     meaning runs with different radiation timesteps were inconsistent.
!
!   - The time averaging was therefore modified so that diurnal_solar now returns coszen time-averaged between t and t+dt,
!     meaning that runs with different radiation timesteps are now consistent, and _multiplication by fracday is not necessary_.
!
!   - Extra documentation has been added to the time-averaged calculation (shown by !st_doc)

!---------------------------------------------------------------------
!-------  interfaces --------

public       &
              astronomy_init, get_period, set_period, &
              set_orbital_parameters, get_orbital_parameters, &
              set_ref_date_of_ae, get_ref_date_of_ae,  &
              diurnal_solar, daily_mean_solar, annual_mean_solar,  &
              astronomy_end, universal_time, orbital_time, diurnal_exoplanet

interface diurnal_solar
   module procedure diurnal_solar_2d
   module procedure diurnal_solar_1d
   module procedure diurnal_solar_0d
   module procedure diurnal_solar_cal_2d
   module procedure diurnal_solar_cal_1d
   module procedure diurnal_solar_cal_0d
end interface

interface daily_mean_solar
   module procedure daily_mean_solar_2d
   module procedure daily_mean_solar_1d
   module procedure daily_mean_solar_2level
   module procedure daily_mean_solar_0d
   module procedure daily_mean_solar_cal_2d
   module procedure daily_mean_solar_cal_1d
   module procedure daily_mean_solar_cal_2level
   module procedure daily_mean_solar_cal_0d
end interface

interface annual_mean_solar
   module procedure annual_mean_solar_2d
   module procedure annual_mean_solar_1d
   module procedure annual_mean_solar_2level
end interface

interface get_period
   module procedure get_period_time_type, get_period_integer
end interface

interface set_period
   module procedure set_period_time_type, set_period_integer
end interface


private &

! called from astronomy_init and set_orbital_parameters:
              orbit,  &

! called from diurnal_solar, daily_mean_solar and orbit:
              r_inv_squared,   &

! called from  diurnal_solar and daily_mean_solar:
              angle,  declination,  &
              half_day
!             half_day, orbital_time, &

! called from  diurnal_solar:
!             universal_time


interface half_day
   module procedure half_day_2d, half_day_0d
end interface


!---------------------------------------------------------------------
!-------- namelist  ---------


!real   :: ecc   = 0.01671   ! eccentricity of orbital ellipse
real, public    :: ecc   = 0.0       ! eccentricity of orbital ellipse !s changed default eccentricity to zero, as this is how most 'simplified' models are run.
                             ! [ non-dimensional ]
real, public    :: obliq = 23.439    ! obliquity [ degrees ]
real    :: per   = 102.932   ! longitude of perihelion with respect
                             ! to autumnal equinox in NH [ degrees ]
integer :: period = 0        ! specified length of year [ seconds ] ;
                             ! must be specified to override default
                             ! value given by length_of_year in
                             ! time_manager_mod
integer :: day_ae    = 23    ! day of specified autumnal equinox
integer :: month_ae  = 9     ! month of specified autumnal equinox
integer :: year_ae   = 1998  ! year of specified autumnal equinox
integer :: hour_ae   = 5     ! hour of specified autumnal equinox
integer :: minute_ae = 37    ! minute of specified autumnal equinox
integer :: second_ae = 0     ! second of specified autumnal equinox
integer :: num_angles = 3600 ! number of intervals into which the year
                             ! is divided to compute orbital positions


namelist /astronomy_nml/ ecc, obliq, per, period, &
                         year_ae, month_ae,  day_ae,         &
                         hour_ae, minute_ae, second_ae, &
                         num_angles

!--------------------------------------------------------------------
!------   public data ----------


!--------------------------------------------------------------------
!------   private data ----------

type(time_type) :: autumnal_eq_ref  ! time_type variable containing
                                    ! specified time of reference
                                    ! NH autumnal equinox

type(time_type) :: period_time_type ! time_type variable containing
                                    ! period of one orbit

real, dimension(:), allocatable ::      &
                          orb_angle ! table of orbital positions (0 to
                                    ! 2*pi) as a function of time  used
                                    ! to find actual orbital position
                                    ! via interpolation

real    :: seconds_per_day=86400.   ! seconds in a day
real    :: deg_to_rad               ! conversion from degrees to radians
real    :: twopi                    ! 2 *PI
logical :: module_is_initialized=     &
                          .false.   ! has the module been initialized ?

real, dimension(:), allocatable ::       &
                       cosz_ann, &  ! annual mean cos of zenith angle
                       solar_ann, & ! annual mean solar factor
                       fracday_ann  ! annual mean daylight fraction
real    :: rrsun_ann                ! annual mean earth-sun distance
logical :: annual_mean_calculated =      &
                           .false.  ! have the annual mean values been
                                    ! calculated ?
integer :: num_pts = 0              ! count of grid_boxes for which
                                    ! annual mean astronomy values have
                                    ! been calculated
integer :: total_pts                ! number of grid boxes owned by the
                                    ! processor


!--------------------------------------------------------------------
!--------------------------------------------------------------------



                           contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!####################################################################
! <SUBROUTINE NAME="astronomy_init">
!  <OVERVIEW>
!    astronomy_init is the constructor for astronomy_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    astronomy_init is the constructor for astronomy_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call astronomy_init (latb, lonb)
!  </TEMPLATE>
!  <IN NAME="latb" TYPE="real">
!   2d array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   2d array of model longitudes at cell corners [radians]
!  </IN>
! </SUBROUTINE>
!
subroutine astronomy_init (latb, lonb)

!-------------------------------------------------------------------
!    astronomy_init is the constructor for astronomy_mod.
!-------------------------------------------------------------------

real,   dimension(:,:), intent(in), optional   :: latb
real,   dimension(:,:), intent(in), optional   :: lonb

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       latb         2d array of model latitudes at cell corners
!                    [ radians ]
!       lonb         2d array of model longitudes at cell corners
!                    [ radians ]
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      integer                         :: unit, ierr, io, seconds,  &
                                         days, jd, id

!-------------------------------------------------------------------
!  local variables:
!
!      unit
!      ierr
!      io
!      seconds
!      days
!      jd
!      id
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    if module has already been initialized, exit.
!-------------------------------------------------------------------
      if (module_is_initialized) return

!-------------------------------------------------------------------
!    verify that modules used by this module have been initialized.
!-------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call constants_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, astronomy_nml, iostat=io)
      ierr = check_nml_error(io,'astronomy_nml')
#else
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=astronomy_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'astronomy_nml')
        end do
10      call close_file (unit)
      endif
#endif
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) then
        unit = stdlog()
        write (unit, nml=astronomy_nml)
      endif
!--------------------------------------------------------------------
!    be sure input values are within valid ranges.
!    QUESTION : ARE THESE THE RIGHT LIMITS ???
!---------------------------------------------------------------------
      if (ecc < 0.0 .or. ecc > 0.99) &
        call error_mesg ('astronomy_mod', &
                      'ecc must be between 0 and 0.99', FATAL)
      if (obliq < -90. .or. obliq > 90.) &
        call error_mesg ('astronomy_mod', &
                  'obliquity must be between -90 and 90 degrees', FATAL)
      if (per <  0.0 .or. per > 360.0) &
        call error_mesg ('astronomy_mod', &
                 'perihelion must be between 0 and 360 degrees', FATAL)

!----------------------------------------------------------------------
!    set up time-type variable defining specified time of autumnal
!    equinox.
!----------------------------------------------------------------------
    if (get_calendar_type() /=NO_CALENDAR) then
      autumnal_eq_ref = set_date (year_ae,month_ae,day_ae, &
                                  hour_ae,minute_ae,second_ae)
    endif
!---------------------------------------------------------------------
!    set up time-type variable defining length of year.
!----------------------------------------------------------------------
      if (period == 0) then
        period_time_type = length_of_year()
        call get_time (period_time_type, seconds, days)
        period = seconds_per_day*days + seconds
      else
        period_time_type = set_time(period,0)
      endif

!---------------------------------------------------------------------
!    define useful module variables.
!----------------------------------------------------------------------
      twopi = 2.*PI
      deg_to_rad = twopi/360.

!---------------------------------------------------------------------
!    call orbit to define table of orbital angles as function of
!    orbital time.
!----------------------------------------------------------------------
! wfc moved here from orbit
      allocate ( orb_angle(0:num_angles) )
      call orbit

!--------------------------------------------------------------------
!    if annual mean radiation is desired, then latb will be present.
!    allocate arrays to hold the needed astronomical factors. define
!    the total number of points that the processor is responsible for.
!--------------------------------------------------------------------
      if (present(latb)) then
        jd = size(latb,2) - 1
        id = size(lonb,1) - 1
        allocate (cosz_ann(jd))
        allocate (solar_ann(jd))
        allocate (fracday_ann(jd))
        total_pts = jd*id
      endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized=.true.

!---------------------------------------------------------------------

end subroutine astronomy_init




!###################################################################



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    INTERFACE GET_PERIOD
!
!
! call get_period (period)
!
!  separate routines exist within this interface for integer
!  and time_type output:
!
!  integer, intent(out)         :: period
! OR
!  type(time_type), intent(out) :: period
!
!--------------------------------------------------------------------
!
!   intent(out) variable:
!
!    period_out        length of year for calendar in use
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="get_period_integer">
!  <OVERVIEW>
!    get_period_integer returns the length of the year as an integer
!    number of seconds.
!  </OVERVIEW>
!  <DESCRIPTION>
!    get_period_integer returns the length of the year as an integer
!    number of seconds.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_period_integer (period_out)
!  </TEMPLATE>
!  <OUT NAME="period_out" TYPE="integer">
!   number of seconds as the length of the year
!  </OUT>
! </SUBROUTINE>
!
subroutine get_period_integer (period_out)

!--------------------------------------------------------------------
!    get_period_integer returns the length of the year as an integer
!    number of seconds.
!--------------------------------------------------------------------

integer, intent(out) :: period_out

!--------------------------------------------------------------------
!   local variables:

      integer :: seconds, days

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    define length of year in seconds.
!--------------------------------------------------------------------
      call get_time (period_time_type, seconds, days)
      period_out = seconds_per_day*days + seconds


end subroutine get_period_integer

!####################################################################
! <SUBROUTINE NAME="get_period_time_type">
!  <OVERVIEW>
!    get_period_time_type returns the length of the year as a time_type
!    variable.
!  </OVERVIEW>
!  <DESCRIPTION>
!    get_period_time_type returns the length of the year as a time_type
!    variable.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_period_time_type (period_out)
!  </TEMPLATE>
!  <OUT NAME="period_out" TYPE="time_type">
!   the length of the year as a time_type
!  </OUT>
! </SUBROUTINE>
!
subroutine get_period_time_type (period_out)

!--------------------------------------------------------------------
!    get_period_time_type returns the length of the year as a time_type
!    variable.
!--------------------------------------------------------------------

type(time_type), intent(inout) :: period_out

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ( 'astronomy_mod',  &
              ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    define length of year as a time_type variable.
!--------------------------------------------------------------------
      period_out = period_time_type


end subroutine get_period_time_type



!#####################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                END INTERFACE GET_PERIOD
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    INTERFACE SET_PERIOD
!
!
! call set_period (period_in)
!
!  separate routines exist within this interface for integer
!  and time_type output:
!
!  integer, intent(out)         :: period_in
! OR
!  type(time_type), intent(out) :: period_in
!
!--------------------------------------------------------------------
!
!   intent(in) variable:
!
!    period_in        length of year for calendar in use
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="set_period_integer">
!  <OVERVIEW>
!    set_period_integer saves as the input length of the year (an
!    integer) in a time_type module variable.
!  </OVERVIEW>
!  <DESCRIPTION>
!    set_period_integer saves as the input length of the year (an
!    integer) in a time_type module variable.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call set_period_integer (period_in)
!  </TEMPLATE>
!  <IN NAME="period_in" TYPE="time_type">
!   the length of the year as a time_type
!  </IN>
! </SUBROUTINE>
!
subroutine set_period_integer (period_in)

!--------------------------------------------------------------------
!    set_period_integer saves as the input length of the year (an
!    integer) in a time_type module variable.
!--------------------------------------------------------------------

integer, intent(in) :: period_in

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    define time_type variable defining the length of year from the
!    input value (integer seconds).
!---------------------------------------------------------------------
      period_time_type = set_time(period_in, 0)



end subroutine set_period_integer



!#####################################################################

subroutine set_period_time_type(period_in)

!--------------------------------------------------------------------
!    set_period_time_type saves the length of the year (input as a
!    time_type variable) into a time_type module variable.
!--------------------------------------------------------------------

type(time_type), intent(in) :: period_in

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    define time_type variable defining the length of year from the
!    input value (time_type).
!---------------------------------------------------------------------
      period_time_type = period_in


end subroutine set_period_time_type



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                END INTERFACE SET_PERIOD
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!#####################################################################
! <SUBROUTINE NAME="set_orbital_parameters">
!  <OVERVIEW>
!    set_orbital_parameters saves the input values of eccentricity,
!    obliquity and perihelion time as module variables for use by
!    astronomy_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    set_orbital_parameters saves the input values of eccentricity,
!    obliquity and perihelion time as module variables for use by
!    astronomy_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call set_orbital_parameters (ecc_in, obliq_in, per_in)
!  </TEMPLATE>
!  <IN NAME="ecc_in" TYPE="real">
!   eccentricity of orbital ellipse
!  </IN>
!  <IN NAME="obliq_in" TYPE="real">
!   obliquity fof orbital ellipse
!  </IN>
!  <IN NAME="per_in" TYPE="real">
!   longitude of perihelion with respect to autumnal
!                      equinox in northern hemisphere
!  </IN>
! </SUBROUTINE>
!
subroutine set_orbital_parameters (ecc_in, obliq_in, per_in)

!--------------------------------------------------------------------
!    set_orbital_parameters saves the input values of eccentricity,
!    obliquity and perihelion time as module variables for use by
!    astronomy_mod.
!--------------------------------------------------------------------

real, intent(in) :: ecc_in
real, intent(in) :: obliq_in
real, intent(in) :: per_in

!--------------------------------------------------------------------
!
!  intent(in) variables:
!
!     ecc_in           eccentricity of orbital ellipse
!                      [ non-dimensional ]
!     obliq_in         obliquity
!                      [ degrees ]
!     per_in           longitude of perihelion with respect to autumnal
!                      equinox in northern hemisphere
!                      [ degrees ]
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    be sure input values are within valid ranges.
!    QUESTION : ARE THESE THE RIGHT LIMITS ???
!---------------------------------------------------------------------
      if (ecc_in < 0.0 .or. ecc_in > 0.99) &
        call error_mesg ('astronomy_mod', &
                      'ecc must be between 0 and 0.99', FATAL)
      if (obliq_in < -90.0 .or. obliq > 90.0) &
        call error_mesg ('astronomy_mod', &
                'obliquity must be between -90. and 90. degrees', FATAL)
      if (per_in < 0.0 .or. per_in > 360.0) &
        call error_mesg ('astronomy_mod', &
              'perihelion must be between 0.0 and 360. degrees', FATAL)

!---------------------------------------------------------------------
!    save input values into module variables.
!---------------------------------------------------------------------
      ecc   = ecc_in
      obliq = obliq_in
      per   = per_in

!---------------------------------------------------------------------
!    call orbit to define table of orbital angles as function of
!    orbital time using the input values of parameters just supplied.
!----------------------------------------------------------------------
      call orbit

!----------------------------------------------------------------------



end subroutine set_orbital_parameters



!####################################################################
! <SUBROUTINE NAME="get_orbital_parameters">
!  <OVERVIEW>
!    get_orbital_parameters retrieves the orbital parameters for use
!    by another module.
!  </OVERVIEW>
!  <DESCRIPTION>
!    get_orbital_parameters retrieves the orbital parameters for use
!    by another module.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_orbital_parameters (ecc_out, obliq_out, per_out)
!  </TEMPLATE>
!  <OUT NAME="ecc_out" TYPE="real">
!   eccentricity of orbital ellipse
!  </OUT>
!  <OUT NAME="obliq_out" TYPE="real">
!   obliquity fof orbital ellipse
!  </OUT>
!  <OUT NAME="per_out" TYPE="real">
!   longitude of perihelion with respect to autumnal
!                      equinox in northern hemisphere
!  </OUT>
! </SUBROUTINE>
!
subroutine get_orbital_parameters (ecc_out, obliq_out, per_out)

!-------------------------------------------------------------------
!    get_orbital_parameters retrieves the orbital parameters for use
!    by another module.
!--------------------------------------------------------------------

real, intent(out) :: ecc_out, obliq_out, per_out

!--------------------------------------------------------------------
!
!  intent(out) variables:
!
!     ecc_out          eccentricity of orbital ellipse
!                      [ non-dimensional ]
!     obliq_out        obliquity
!                      [ degrees ]
!     per_out          longitude of perihelion with respect to autumnal
!                      equinox in northern hemisphere
!                      [ degrees ]
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    fill the output arguments with the eccentricity, obliquity and
!    perihelion angle.
!--------------------------------------------------------------------
      ecc_out = ecc
      obliq_out = obliq
      per_out = per


end subroutine get_orbital_parameters



!####################################################################
! <SUBROUTINE NAME="set_ref_date_of_ae">
!  <OVERVIEW>
!    set_ref_date_of_ae provides a means of specifying the reference
!    date of the NH autumnal equinox for a particular year.
!  </OVERVIEW>
!  <DESCRIPTION>
!    set_ref_date_of_ae provides a means of specifying the reference
!    date of the NH autumnal equinox for a particular year.  it is only
!    used if calls are made to the calandar versions of the routines
!    diurnal_solar and daily_mean_solar. if the NOLEAP calendar is
!    used, then the date of autumnal equinox will be the same every
!    year. if JULIAN is used, then the date of autumnal equinox will
!    return to the same value every 4th year.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call set_ref_date_of_ae (day_in,month_in,year_in, &
!                               second_in,minute_in,hour_in)
!  </TEMPLATE>
!  <IN NAME="day_in" TYPE="integer">
!   day of reference autumnal equinox
!  </IN>
!  <IN NAME="month_in" TYPE="integer">
!   month of reference autumnal equinox
!  </IN>
!  <IN NAME="year_in" TYPE="integer">
!   year of reference autumnal equinox
!  </IN>
!  <IN NAME="second_in" TYPE="real">
!   OPTIONAL: second of reference autumnal equinox
!  </IN>
!  <IN NAME="minute_in" TYPE="real">
!   OPTIONAL: minute of reference autumnal equinox
!  </IN>
!  <IN NAME="hour_in" TYPE="real">
!   OPTIONAL: hour of reference autumnal equinox
!  </IN>
! </SUBROUTINE>
!
subroutine set_ref_date_of_ae (day_in,month_in,year_in, &
                               second_in,minute_in,hour_in)

!---------------------------------------------------------------------
!    set_ref_date_of_ae provides a means of specifying the reference
!    date of the NH autumnal equinox for a particular year.  it is only
!    used if calls are made to the calandar versions of the routines
!    diurnal_solar and daily_mean_solar. if the NOLEAP calendar is
!    used, then the date of autumnal equinox will be the same every
!    year. if JULIAN is used, then the date of autumnal equinox will
!    return to the same value every 4th year.
!----------------------------------------------------------------------

integer, intent(in)           :: day_in, month_in, year_in
integer, intent(in), optional :: second_in, minute_in, hour_in

!--------------------------------------------------------------------
!
!  intent(in) variables:
!
!     day_in           day of reference autumnal equinox
!                      [ non-dimensional ]
!     month_in         month of reference autumnal equinox
!                      [ non-dimensional ]
!     year_in          year of reference autumnal equinox
!                      [ non-dimensional ]
!
!  intent(in), optional variables:
!
!     second_in        second of reference autumnal equinox
!                      [ non-dimensional ]
!     minute_in        minute of reference autumnal equinox
!                      [ non-dimensional ]
!     hour_in          hour of reference autumnal equinox
!                      [ non-dimensional ]
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    save the input time of ae specification into a time_type module
!    variable autumnal_eq_ref.
!--------------------------------------------------------------------
      day_ae =    day_in
      month_ae =  month_in
      year_ae =   year_in

      if (present(second_in)) then
        second_ae = second_in
        minute_ae = minute_in
        hour_ae =   hour_in
      else
        second_ae = 0
        minute_ae = 0
        hour_ae   = 0
      endif

      autumnal_eq_ref = set_date (year_ae,month_ae,day_ae, &
                                  hour_ae,minute_ae,second_ae)

!---------------------------------------------------------------------


end subroutine set_ref_date_of_ae



!####################################################################
! <SUBROUTINE NAME="get_ref_date_of_ae">
!  <OVERVIEW>
!     get_ref_date_of_ae retrieves the reference date of the autumnal
!     equinox as integer variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!     get_ref_date_of_ae retrieves the reference date of the autumnal
!     equinox as integer variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_ref_date_of_ae (day_out,month_out,year_out, &
!                               second_out, minute_out,hour_out)
!  </TEMPLATE>
!  <OUT NAME="day_out" TYPE="integer">
!   day of reference autumnal equinox
!  </OUT>
!  <OUT NAME="month_out" TYPE="integer">
!   month of reference autumnal equinox
!  </OUT>
!  <OUT NAME="year_out" TYPE="integer">
!   year of reference autumnal equinox
!  </OUT>
!  <OUT NAME="second_out" TYPE="real">
!   second of reference autumnal equinox
!  </OUT>
!  <OUT NAME="minute_out" TYPE="real">
!   minute of reference autumnal equinox
!  </OUT>
!  <OUT NAME="hour_out" TYPE="real">
!   hour of reference autumnal equinox
!  </OUT>
! </SUBROUTINE>
!
subroutine get_ref_date_of_ae (day_out,month_out,year_out,&
                               second_out,minute_out,hour_out)

!---------------------------------------------------------------------
!     get_ref_date_of_ae retrieves the reference date of the autumnal
!     equinox as integer variables.
!---------------------------------------------------------------------

integer, intent(out) :: day_out, month_out, year_out,  &
                        second_out, minute_out, hour_out

!--------------------------------------------------------------------
!
!  intent(out) variables:
!
!     day_out          day of reference autumnal equinox
!                      [ non-dimensional ]
!     month_out        month of reference autumnal equinox
!                      [ non-dimensional ]
!     year_out         year of reference autumnal equinox
!                      [ non-dimensional ]
!     second_out       second of reference autumnal equinox
!                      [ non-dimensional ]
!     minute_out       minute of reference autumnal equinox
!                      [ non-dimensional ]
!     hour_out         hour of reference autumnal equinox
!                      [ non-dimensional ]
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    exit if module has not been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized)   &
        call error_mesg ( 'astronomy_mod',  &
         ' module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    fill the output fields with the proper module data.
!---------------------------------------------------------------------
      day_out    =  day_ae
      month_out  =  month_ae
      year_out   =  year_ae
      second_out =  second_ae
      minute_out =  minute_ae
      hour_out   =  hour_ae


end subroutine get_ref_date_of_ae



!#####################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    INTERFACE DIURNAL_SOLAR
!
! call diurnal_solar (lat, lon, time, cosz, fracday, rrsun, dt_time)
!   or
! call diurnal_solar (lat, lon, gmt, time_since_ae, cosz, fracday,
!                     rrsun, dt)
!
!  the first option (used in conjunction with time_manager_mod)
!  generates the real variables gmt and time_since_ae from the
!  time_type input, and then calls diurnal_solar with these real
!  inputs.
!
!  the time of day is set by
!    real, intent(in) :: gmt
!  the time of year is set by
!    real, intent(in) :: time_since_ae
!  with time_type input, both of these are extracted from
!    type(time_type), intent(in) :: time
!
!
!  separate routines exist within this interface for scalar,
!  1D or 2D input and output fields:
!
!    real, intent(in), dimension(:,:) :: lat, lon
! OR real, intent(in), dimension(:)   :: lat, lon
! OR real, intent(in)                 :: lat, lon
!
!    real, intent(out), dimension(:,:) :: cosz, fracday
! OR real, intent(out), dimension(:)   :: cosz, fracday
! OR real, intent(out)                 :: cosz, fracday
!
!  one may also average the output fields over the time interval
!  between gmt and gmt + dt by including the optional argument dt (or
!  dt_time). dt is measured in radians and must be less than pi
!  (1/2 day). this average is computed analytically, and should be
!  exact except for the fact that changes in earth-sun distance over
!  the time interval dt are ignored. in the context of a diurnal GCM,
!  this option should always be employed to insure that the total flux
!  at the top of the atmosphere is not modified by time truncation
!  error.
!
!    real, intent(in), optional :: dt
!    type(time_type), optional :: dt_time
! (see test.90 for examples of the use of these types)
!
!--------------------------------------------------------------------
!
!  intent(in) variables:
!
!     lat            latitudes of model grid points
!                    [ radians ]
!     lon            longitudes of model grid points
!                    [ radians ]
!     gmt            time of day at longitude 0.0; midnight = 0.0,
!                    one day = 2 * pi
!                    [ radians ]
!     time_since_ae  time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!                    [ radians ]
!     time           time at which astronomical values are desired
!                    time_type variable [ seconds, days]
!
!
!  intent(out) variables:
!
!     cosz           cosine of zenith angle, set to zero when entire
!                    period is in darkness
!                    [ dimensionless ]
!     fracday        daylight fraction of time interval
!                    [ dimensionless ]
!     rrsun          earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!                    [ dimensionless ]
!
!  intent(in), optional variables:
!
!     dt            time interval after gmt over which the astronomical
!                   variables are to be averaged. this produces averaged
!                   output rather than instantaneous.
!                   [ radians ], (1 day = 2 * pi)
!     dt_time       as in dt, but dt_time is a time_type variable
!                   time_type [ days, seconds ]
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="diurnal_solar_2d">
!  <OVERVIEW>
!    diurnal_solar_2d returns 2d fields of cosine of zenith angle,
!    daylight fraction and earth-sun distance at the specified lati-
!    tudes, longitudes and time. these values may be instantaneous
!    or averaged over a specified time interval.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diurnal_solar_2d returns 2d fields of cosine of zenith angle,
!    daylight fraction and earth-sun distance at the specified lati-
!    tudes, longitudes and time. these values may be instantaneous
!    or averaged over a specified time interval.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diurnal_solar_2d (lat, lon, gmt, time_since_ae, cosz, &
!                             fracday, rrsun, dt_time)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   longitude of model grid points
!  </IN>
!  <IN NAME="gmt" TYPE="real">
!   time of day at longitude 0.0; midnight = 0.0,
!                    one day = 2 * pi
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
!  <IN NAME="dt" TYPE="real">
!   OPTIONAL: time interval after gmt over which the astronomical
!                   variables are to be averaged. this produces averaged
!                   output rather than instantaneous.
!  </IN>
! </SUBROUTINE>
!
subroutine diurnal_solar_2d (lat, lon, gmt, time_since_ae, cosz, &
                             fracday, rrsun, dt, allow_negative_cosz, &
                             half_day_out)

!---------------------------------------------------------------------
!    diurnal_solar_2d returns 2d fields of cosine of zenith angle,
!    daylight fraction and earth-sun distance at the specified lati-
!    tudes, longitudes and time. these values may be instantaneous
!    or averaged over a specified time interval.
!---------------------------------------------------------------------

real, dimension(:,:), intent(in)           :: lat, lon
real,                 intent(in)           :: gmt, time_since_ae
real, dimension(:,:), intent(out)          :: cosz, fracday
real,                 intent(out)          :: rrsun
real,                 intent(in), optional :: dt
logical,              intent(in), optional :: allow_negative_cosz
real, dimension(:,:), intent(out), optional :: half_day_out


!---------------------------------------------------------------------
!   local variables

      real, dimension(size(lat,1),size(lat,2)) :: t, tt, h, aa, bb,  &
                                                  st, stt, sh
      real                                     :: ang, dec
      logical :: Lallow_negative

!---------------------------------------------------------------------
!   local variables
!
!    t           time of day with respect to local noon (2 pi = 1 day)
!                [ radians ]
!    tt          end of averaging period [ radians ]
!    h           half of the daily period of daylight, centered at noon
!                [ radians, -pi --> pi ]
!    aa          sin(lat) * sin(declination)
!    bb          cos(lat) * cos(declination)
!    st          sine of time of day
!    stt         sine of time of day at end of averaging period
!    sh          sine of half-day period
!    ang         position of earth in its orbit wrt autumnal equinox
!                [ radians ]
!    dec         earth's declination [ radians ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    be sure the time in the annual cycle is legitimate.
!---------------------------------------------------------------------
      if (time_since_ae < 0.0 .or. time_since_ae > twopi) then
          call error_mesg('astronomy_mod', &
                    'time_since_ae not between 0 and 2pi', FATAL)
      end if
!--------------------------------------------------------------------
!    be sure the time at longitude = 0.0 is legitimate.
!---------------------------------------------------------------------
      if (gmt < 0.0 .or. gmt > twopi) &
         call error_mesg('astronomy_mod', &
                    'gmt not between 0 and 2pi', FATAL)

!---------------------------------------------------------------------
!    define the orbital angle (location in year), solar declination and
!    earth sun distance factor. use functions contained in this module.
!---------------------------------------------------------------------
      ang = angle(time_since_ae)
      dec = declination(ang)
      rrsun  = r_inv_squared(ang)

!---------------------------------------------------------------------
!    define terms needed in the cosine zenith angle equation.
!--------------------------------------------------------------------
      aa = sin(lat)*sin(dec)
      bb = cos(lat)*cos(dec)

!---------------------------------------------------------------------
!    define local time. force it to be between -pi and pi.
!--------------------------------------------------------------------
      t = gmt + lon - PI
      where(t >= PI) t = t - twopi
      where(t < -PI) t = t + twopi

      Lallow_negative = .false.
      if (present(allow_negative_cosz)) then
         if (allow_negative_cosz) Lallow_negative = .true.
      end if

!---------------------------------------------------------------------
!    perform a time integration to obtain cosz and fracday if desired.
!    output is valid over the period from t to t + dt.
!--------------------------------------------------------------------
      h   = half_day   (lat,dec)

      if ( present(half_day_out) ) then
         half_day_out = h
      end if

      if ( present(dt) ) then   ! (perform time averaging)
        tt = t + dt
        st  = sin(t)
        stt = sin(tt)
        sh  = sin(h)
        cosz = 0.0

        if (.not. Lallow_negative) then

!st_doc-------------------------------------------------------------
!     ************** BASIC diurnal SOLAR EXPLANATION ***************
!     Key to understanding code is that coszen is not just coszen = aa + bb*cos(time), but is really max(aa+bb*cos(time),0.).
!     Essentially this is because if it is night-time, the incoming solar = 0, but if it's day time, it's = aa+bb*cos(time).
!     When time-averaging between t and tt=t+dt, this must be accounted for.
!     To illustrate this, I will compare two example cases from the 8 cases used in the calculation of `diurnal_solar':
!     Case 4 below is simple, t is after sunrise, and tt is before sunset, so there is daylight for the entire averaging period.
!     Therefore we time-average by doing:
!
!     \int_t^tt cosz dt' / (tt-t) = (\int_t^tt aa dt') / (tt-t) + bb*(\int_t^tt cos(t') dt')/ (tt - t)
!
!     where t' is a dummy variable and \int_t^tt is the integral between t and tt.
!
!     This comes out to be:
!     coszen_time_av = aa + bb(sin(tt) - sin(t))/(tt-t) = aa + bb(stt-st)/(tt-t)
!     because aa and bb are independent of time (this is not strictly true, as the declination, and hence aa and bb, does change
!     with time. However, we choose to approximate aa and bb as constant in the period t->t+dt, which is a reasonable approximation
!     if dt is a small fraction of the annual cycle).
!
!     Other cases are more complicated, because the averaging periods are between two times that include periods of day and night.
!     For example, case 2 has t that is before sunrise and tt that is before sunset.
!     So now we have to separately consider those times during the day and during the night:

!     \int_t^tt cosz dt' / (tt-t) = \int_t^-h 0 dt' / (tt-t) + \int_-h^tt cosz dt' / (tt-t)
!     where the first integral on the RHS is from t until sunrise (-h). This is night, so we have zero contribution
!     from this period. The second integral is during the day, so is from sunrise (-h) to tt.
!     The final answer is then:
!     coszen_time_av = aa(tt-(-h))/(tt-t) + bb(sin(tt) - sin(-h))/(tt-t) = aa*(tt+h)/(tt-t) + bb*(stt+sh)/(tt-t)
!
!     ******* EXPLANATION OF previous error (/ undocumented feature) in astronomy module *******
!
!     The original `diurnal_solar' subroutine claimed to return the time-averaged value of coszen if the variable dt is passed.
!     However, the time-averaging was only done DURING HOURS OF DAYLIGHT, so in case 2, the answer was:
!     coszen_time_av = aa + bb(stt+sh)/(tt+h)
!     and not:
!     coszen_time_av = aa + bb(stt+sh)/(tt-t)
!     This meant that the original code returned time-averages that were not from t->t+dt, but over the daylight hours between t
!     and t+dt, which was not what the code claimed to return.
!
!     I have therefore changed all the denominators in the code below to be (tt-t), so that the time-averaging takes place
!     over the time-period t->t+dt. This means that the coszen returned by `diurnal_solar' is a true time-average over dt.
!     An equivalent way of doing this with the old code would be to multiply coszen by fracday, where fracday is the fraction of
!     the time period that is daylight.
!end st_doc

!-------------------------------------------------------------------
!    case 1: entire averaging period is before sunrise.
!-------------------------------------------------------------------
        where (t < -h .and. tt < -h) cosz = 0.0

!-------------------------------------------------------------------
!    case 2: averaging period begins before sunrise, ends after sunrise
!    but before sunset
!-------------------------------------------------------------------
! st_doc
! t < -h first time is before sunrise
! abs(tt) <= h for tt>0 this means tt < h, meaning it's before sunset.
! abs(tt) <= h for tt<0 this means tt > -h, meaning it's after sunrise.
! Time average is between t and tt, so the denominator is (tt-t). But in the numerator sin(t) will be < 0 as t is after sunset. So only use up to sin(-h), which gives (stt - sin(-h)) - > stt+sh.
! end st_doc

	where ( t < -h .and. abs(tt) <= h)   &
             cosz = (aa*(tt+h)/(tt-t)) + bb*(stt +sh)/ (tt-t)

!-------------------------------------------------------------------
!    case 3: averaging period begins before sunrise, ends after sunset,
!    but before the next sunrise. modify if averaging period extends
!    past the next day's sunrise, but if averaging period is less than
!    a half- day (pi) that circumstance will never occur.
!-------------------------------------------------------------------
        where (t < -h .and. h /= 0.0 .and. h < tt)    &
              cosz = aa*(2.*h)/(tt-t) + bb*( sh + sh)/(tt-t)

!-------------------------------------------------------------------
!    case 4: averaging period begins after sunrise, ends before sunset.
!-------------------------------------------------------------------
! st_doc
! abs(t) <= h for t>0 this means t < h, meaning it's before sunset.
! abs(t) <= h for t<0 this means t > -h, meaning it's after sunrise.
! end st_doc

        where ( abs(t) <= h .and. abs(tt) <= h)    &
             cosz = aa + bb*(stt - st)/ (tt - t)

!-------------------------------------------------------------------
!    case 5: averaging period begins after sunrise, ends after sunset.
!    modify when averaging period extends past the next day's sunrise.
!-------------------------------------------------------------------
! st_doc
! tt > h final time is after sunset
! abs(t) <= h for t>0 this means t < h, meaning it's before sunset.
! abs(t) <= h for t<0 this means t > -h, meaning it's after sunrise.
! (sh - st) / (tt - t) b/c we're averaging from start time to sunset.
! end st_doc

        where ( abs(t) <= h .and.  h < tt)    &
              cosz = (aa*(h-t)/(tt-t)) + bb*(sh - st)/(tt-t)

!-------------------------------------------------------------------
!    case 6: averaging period begins after sunrise , ends after the
!    next day's sunrise. note that this includes the case when the
!    day length is one day (h = pi).
!-------------------------------------------------------------------
! st_doc
! twopi-h is the time of the next day's sunrise. So tt > twopi - h means end time is after the next day's sunrise.
! t<=h beginning time is after sunrise if t < 0 here(?!)
! end st_doc

        where (twopi - h < tt .and. t <= h ) &
	   cosz = aa*((tt+(2.*h)-t-twopi)/(tt-t)) + bb*(((sh-st)/(tt-t)) + ((stt+sh)/(tt-t)))

!-------------------------------------------------------------------
!    case 7: averaging period begins after sunset and ends before the
!    next day's sunrise
!-------------------------------------------------------------------
        where(  h <  t .and. twopi - h >= tt  ) cosz = 0.0

!-------------------------------------------------------------------
!    case 8: averaging period begins after sunset and ends after the
!    next day's sunrise but before the next day's sunset. if the
!    averaging period is less than a half-day (pi) the latter
!    circumstance will never occur.
!-----------------------------------------------------------------
        where(  h <  t .and. twopi - h < tt .and. tt < twopi+h )  cosz = aa*(tt + h - twopi)/(tt-t) + bb*(stt + sh) / (tt -t )

!-------------------------------------------------------------------
!    case 9: averaging period begins after sunset and ends after the
!    next day's sunset. if the
!    averaging period is less than a half-day (pi) the latter
!    circumstance will never occur.
!-----------------------------------------------------------------
        where(  h <  t .and. twopi - h < tt .and. tt > twopi+h )  cosz = aa*(2.*h)/(tt-t) + bb*(sh + sh) / (tt -t )
        
        else
           cosz = aa + bb*(stt - st)/ (tt - t)
        end if

!-------------------------------------------------------------------
!    day fraction is the fraction of the averaging period contained
!    within the (-h,h) period.
!-------------------------------------------------------------------
        where (t < -h .and.      tt < -h)      fracday = 0.0
        where (t < -h .and. abs(tt) <= h)      fracday = (tt + h )/dt
        where (t < -h .and.       h < tt)      fracday = ( h + h )/dt
        where (abs(t) <= h .and. abs(tt) <= h) fracday = (tt - t )/dt
        where (abs(t) <= h .and.       h < tt) fracday = ( h - t )/dt
        where (      h <  t                 )  fracday = 0.0
        where (twopi - h < tt)                 fracday = fracday +  &
                                                         (tt + h - &
                                                         twopi)/dt
!----------------------------------------------------------------------
!    if instantaneous values are desired, define cosz at time t.
!----------------------------------------------------------------------
      else  ! (no time averaging)
        if (.not. Lallow_negative) then
           where (abs(t) < h)
             cosz = aa + bb*cos(t)
             fracday = 1.0
           elsewhere
             cosz = 0.0
             fracday = 0.0
           end where
        else
           cosz = aa + bb*cos(t)
           where (abs(t) < h)
             fracday = 1.0
           elsewhere
             fracday = 0.0
           end where
        end if
      end if

!----------------------------------------------------------------------
!    be sure that cosz is not negative.
!----------------------------------------------------------------------
      if (.not. Lallow_negative) then
         cosz = max(0.0, cosz)
      end if

!--------------------------------------------------------------------

end subroutine diurnal_solar_2d


!######################################################################
! <SUBROUTINE NAME="diurnal_solar_1d">
!  <OVERVIEW>
!    diurnal_solar_1d takes 1-d input fields, adds a second dimension
!    and calls diurnal_solar_2d. on return, the 2d fields are returned
!    to the original 1d fields.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diurnal_solar_1d takes 1-d input fields, adds a second dimension
!    and calls diurnal_solar_2d. on return, the 2d fields are returned
!    to the original 1d fields.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diurnal_solar_1d (lat, lon, gmt, time_since_ae, cosz, &
!                             fracday, rrsun, dt)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   longitude of model grid points
!  </IN>
!  <IN NAME="gmt" TYPE="real">
!   time of day at longitude 0.0; midnight = 0.0,
!                    one day = 2 * pi
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
!  <IN NAME="dt" TYPE="real">
!   OPTIONAL: time interval after gmt over which the astronomical
!                   variables are to be averaged. this produces averaged
!                   output rather than instantaneous.
!  </IN>
! </SUBROUTINE>
!
subroutine diurnal_solar_1d (lat, lon, gmt, time_since_ae, cosz, &
                             fracday, rrsun, dt, allow_negative_cosz, &
                             half_day_out)

!--------------------------------------------------------------------
!    diurnal_solar_1d takes 1-d input fields, adds a second dimension
!    and calls diurnal_solar_2d. on return, the 2d fields are returned
!    to the original 1d fields.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:),  intent(in)           :: lat, lon
real,                intent(in)           :: gmt, time_since_ae
real, dimension(:),  intent(out)          :: cosz, fracday
real,                intent(out)          :: rrsun
real,                intent(in), optional :: dt
logical,             intent(in), optional :: allow_negative_cosz
real, dimension(:),  intent(out), optional :: half_day_out

!---------------------------------------------------------------------
!  local variables

      real, dimension(size(lat),1) :: lat_2d, lon_2d, cosz_2d,   &
                                      fracday_2d,halfday_2d

!--------------------------------------------------------------------
!    define 2-d versions of input data arrays.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat
      lon_2d(:,1) = lon

!--------------------------------------------------------------------
!    call diurnal_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
!     if (present(dt)) then
        call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae,&
                               cosz_2d, fracday_2d, rrsun, dt=dt, &
                               allow_negative_cosz=allow_negative_cosz, &
                               half_day_out=halfday_2d)
!     else
!       call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae, &
!                              cosz_2d, fracday_2d, rrsun)
!     endif

!-------------------------------------------------------------------
!    place output fields into 1-d arguments for return to
!    calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(:,1)
      cosz  = cosz_2d (:,1)
      if (present(half_day_out)) then
         half_day_out = halfday_2d(:,1)
      end if


end subroutine diurnal_solar_1d


!#####################################################################
! <SUBROUTINE NAME="diurnal_solar_0d">
!  <OVERVIEW>
!    diurnal_solar_0d takes scalar input fields, makes them into 2d
!    arrays dimensioned (1,1), and calls diurnal_solar_2d. on return,
!    the 2d fields are converted back to the desired scalar output.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diurnal_solar_0d takes scalar input fields, makes them into 2d
!    arrays dimensioned (1,1), and calls diurnal_solar_2d. on return,
!    the 2d fields are converted back to the desired scalar output.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diurnal_solar_0d (lat, lon, gmt, time_since_ae, cosz, &
!                             fracday, rrsun, dt)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   longitude of model grid points
!  </IN>
!  <IN NAME="gmt" TYPE="real">
!   time of day at longitude 0.0; midnight = 0.0,
!                    one day = 2 * pi
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
!  <IN NAME="dt" TYPE="real">
!   OPTIONAL: time interval after gmt over which the astronomical
!                   variables are to be averaged. this produces averaged
!                   output rather than instantaneous.
!  </IN>
! </SUBROUTINE>
!
subroutine diurnal_solar_0d (lat, lon, gmt, time_since_ae, cosz,  &
                             fracday, rrsun, dt, allow_negative_cosz, &
                             half_day_out)

!--------------------------------------------------------------------
!    diurnal_solar_0d takes scalar input fields, makes them into 2d
!    arrays dimensioned (1,1), and calls diurnal_solar_2d. on return,
!    the 2d fields are converted back to the desired scalar output.
!----------------------------------------------------------------------

real, intent(in)           :: lat, lon, gmt, time_since_ae
real, intent(out)          :: cosz, fracday, rrsun
real, intent(in), optional :: dt
logical,intent(in), optional :: allow_negative_cosz
real, intent(out), optional :: half_day_out

!--------------------------------------------------------------------
!  local variables:
!
      real, dimension(1,1) :: lat_2d, lon_2d, cosz_2d, fracday_2d, halfday_2d

!---------------------------------------------------------------------
!    create 2d arrays from the scalar input fields.
!---------------------------------------------------------------------
      lat_2d = lat
      lon_2d = lon

!--------------------------------------------------------------------
!    call diurnal_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
!     if (present(dt)) then
        call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae,  &
                               cosz_2d, fracday_2d, rrsun, dt=dt, &
                               allow_negative_cosz=allow_negative_cosz, &
                               half_day_out=halfday_2d)
!     else
!       call diurnal_solar_2d (lat_2d, lon_2d, gmt, time_since_ae, &
!                              cosz_2d, fracday_2d, rrsun)
!     end if

!-------------------------------------------------------------------
!    place output fields into scalars for return to calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(1,1)
      cosz = cosz_2d(1,1)
      if (present(half_day_out)) then
         half_day_out = halfday_2d(1,1)
      end if



end subroutine diurnal_solar_0d



!####################################################################
! <SUBROUTINE NAME="diurnal_solar_cal_2d">
!  <OVERVIEW>
!    diurnal_solar_cal_2d receives time_type inputs, converts
!    them to real variables and then calls diurnal_solar_2d to
!    compute desired astronomical variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diurnal_solar_cal_2d receives time_type inputs, converts
!    them to real variables and then calls diurnal_solar_2d to
!    compute desired astronomical variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diurnal_solar_cal_2d (lat, lon, gmt, time_since_ae, cosz, &
!                             fracday, rrsun, dt)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   longitude of model grid points
!  </IN>
!  <IN NAME="gmt" TYPE="real">
!   time of day at longitude 0.0; midnight = 0.0,
!                    one day = 2 * pi
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
!  <IN NAME="dt_time" TYPE="time_type">
!   OPTIONAL: time interval after gmt over which the astronomical
!                   variables are to be averaged. this produces averaged
!                   output rather than instantaneous.
!  </IN>
! </SUBROUTINE>
!
subroutine diurnal_solar_cal_2d (lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!-------------------------------------------------------------------
!    diurnal_solar_cal_2d receives time_type inputs, converts
!    them to real variables and then calls diurnal_solar_2d to
!    compute desired astronomical variables.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
real, dimension(:,:), intent(in)            :: lat, lon
type(time_type),      intent(in)            :: time
real, dimension(:,:), intent(out)           :: cosz, fracday
real,                 intent(out)           :: rrsun
type(time_type),      intent(in), optional  :: dt_time
logical,              intent(in), optional  :: allow_negative_cosz
real, dimension(:,:), intent(out), optional  :: half_day_out
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      real :: dt
      real :: gmt, time_since_ae

!---------------------------------------------------------------------
!    extract time of day (gmt) from time_type variable time with
!    function universal_time.
!---------------------------------------------------------------------
      gmt = universal_time(time)

!---------------------------------------------------------------------
!    extract the time of year (time_since_ae) from time_type variable
!    time using the function orbital_time.
!---------------------------------------------------------------------
      time_since_ae = orbital_time(time)

!---------------------------------------------------------------------
!    convert optional time_type variable dt_time (length of averaging
!    period) to a real variable dt with the function universal_time.
!---------------------------------------------------------------------
      if (present(dt_time))  then
        dt = universal_time(dt_time)
        if (dt > PI) then
          call error_mesg ( 'astronomy_mod', &
             'radiation time step must be no longer than 12 hrs', &
                                                          FATAL)
        endif
        if (dt == 0.0) then
          call error_mesg ( 'astronomy_mod', &
              'radiation time step must not be an integral &
                                     &number of days', FATAL)
        endif

!--------------------------------------------------------------------
!    call diurnal_solar_2d to calculate astronomy fields, with or
!    without the optional argument dt.
!--------------------------------------------------------------------
        call diurnal_solar_2d (lat, lon, gmt, time_since_ae, cosz, &
               fracday, rrsun, dt=dt, &
               allow_negative_cosz=allow_negative_cosz, &
               half_day_out=half_day_out)
      else
        call diurnal_solar_2d (lat, lon, gmt, time_since_ae, cosz, &
               fracday, rrsun, &
               allow_negative_cosz=allow_negative_cosz, &
               half_day_out=half_day_out)
      end if


end subroutine diurnal_solar_cal_2d


!#####################################################################
! <SUBROUTINE NAME="diurnal_solar_cal_1d">
!  <OVERVIEW>
!    diurnal_solar_cal_1d receives time_type inputs, converts
!    them to real variables and then calls diurnal_solar_2d to
!    compute desired astronomical variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diurnal_solar_cal_1d receives time_type inputs, converts
!    them to real variables and then calls diurnal_solar_2d to
!    compute desired astronomical variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diurnal_solar_cal_1d (lat, lon, gmt, time_since_ae, cosz, &
!                             fracday, rrsun, dt)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   longitude of model grid points
!  </IN>
!  <IN NAME="gmt" TYPE="real">
!   time of day at longitude 0.0; midnight = 0.0,
!                    one day = 2 * pi
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
!  <IN NAME="dt_time" TYPE="time_type">
!   OPTIONAL: time interval after gmt over which the astronomical
!                   variables are to be averaged. this produces averaged
!                   output rather than instantaneous.
!  </IN>
! </SUBROUTINE>
!
subroutine diurnal_solar_cal_1d (lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!--------------------------------------------------------------------
real, dimension(:), intent(in)           :: lat, lon
type(time_type),    intent(in)           :: time
real, dimension(:), intent(out)          :: cosz, fracday
real,               intent(out)          :: rrsun
type(time_type),    intent(in), optional :: dt_time
logical,            intent(in), optional :: allow_negative_cosz
real, dimension(:), intent(out), optional :: half_day_out
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables

      real, dimension(size(lat),1) :: lat_2d, lon_2d, cosz_2d, &
                                      fracday_2d, halfday_2d

!--------------------------------------------------------------------
!    define 2-d versions of input data arrays.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat
      lon_2d(:,1) = lon

!--------------------------------------------------------------------
!    call diurnal_solar_cal_2d to convert the time_types to reals and
!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      if (present(dt_time)) then
        call diurnal_solar_cal_2d (lat_2d, lon_2d, time, cosz_2d,    &
           fracday_2d, rrsun, dt_time=dt_time, &
           allow_negative_cosz=allow_negative_cosz, &
           half_day_out=halfday_2d)
      else
        call diurnal_solar_cal_2d (lat_2d, lon_2d, time, cosz_2d,    &
           fracday_2d, rrsun, &
           allow_negative_cosz=allow_negative_cosz, &
           half_day_out=halfday_2d)
      end if

!-------------------------------------------------------------------
!    place output fields into 1-d arguments for return to
!    calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(:,1)
      cosz  = cosz_2d (:,1)
      if (present(half_day_out)) then
         half_day_out = halfday_2d(:,1)
      end if


end subroutine diurnal_solar_cal_1d



!#####################################################################
! <SUBROUTINE NAME="diurnal_solar_cal_0d">
!  <OVERVIEW>
!    diurnal_solar_cal_0d receives time_type inputs, converts
!    them to real variables and then calls diurnal_solar_2d to
!    compute desired astronomical variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diurnal_solar_cal_0d receives time_type inputs, converts
!    them to real variables and then calls diurnal_solar_2d to
!    compute desired astronomical variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diurnal_solar_cal_0d (lat, lon, gmt, time_since_ae, cosz, &
!                             fracday, rrsun, dt_time)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   longitude of model grid points
!  </IN>
!  <IN NAME="gmt" TYPE="real">
!   time of day at longitude 0.0; midnight = 0.0,
!                    one day = 2 * pi
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
!  <IN NAME="dt_time" TYPE="time_type">
!   OPTIONAL: time interval after gmt over which the astronomical
!                   variables are to be averaged. this produces averaged
!                   output rather than instantaneous.
!  </IN>
! </SUBROUTINE>
!
subroutine diurnal_solar_cal_0d (lat, lon, time, cosz, fracday,   &
                                 rrsun, dt_time, allow_negative_cosz, &
                                 half_day_out)

!---------------------------------------------------------------------
real,            intent(in)           :: lat, lon
type(time_type), intent(in)           :: time
real,            intent(out)          :: cosz, fracday, rrsun
type(time_type), intent(in), optional :: dt_time
logical,         intent(in), optional :: allow_negative_cosz
real,            intent(out), optional :: half_day_out
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      real, dimension(1,1) :: lat_2d, lon_2d, cosz_2d, fracday_2d, halfday_2d

!--------------------------------------------------------------------
!    define 2-d versions of input data arrays.
!--------------------------------------------------------------------
      lat_2d = lat
      lon_2d = lon

!--------------------------------------------------------------------
!    call diurnal_solar_cal_2d to convert the time_types to reals and
!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      if (present(dt_time)) then
        call diurnal_solar_cal_2d (lat_2d, lon_2d, time, cosz_2d,   &
           fracday_2d, rrsun, dt_time=dt_time, &
           allow_negative_cosz=allow_negative_cosz, &
           half_day_out=halfday_2d)
      else
        call diurnal_solar_cal_2d (lat_2d, lon_2d, time, cosz_2d,   &
           fracday_2d, rrsun, &
           allow_negative_cosz=allow_negative_cosz, &
           half_day_out=halfday_2d)
      end if

!-------------------------------------------------------------------
!    place output fields into 1-d arguments for return to
!    calling routine.
!-------------------------------------------------------------------
      fracday= fracday_2d(1,1)
      cosz = cosz_2d(1,1)
      if (present(half_day_out)) then
         half_day_out = halfday_2d(1,1)
      end if


end subroutine diurnal_solar_cal_0d


!####################################################################

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                END INTERFACE DIURNAL_SOLAR
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                 INTERFACE DAILY_MEAN_SOLAR
!
!
! call daily_mean_solar (lat, time, cosz, fracday, rrsun)
!   or
! call daily_mean_solar (lat, time_since_ae, cosz, fracday, rrsun)
!   or
! call daily_mean_solar (lat, time, cosz, solar)
!   or
! call daily_mean_solar (lat, time_since_ae, cosz, solar)
!
!  the first option (used in conjunction with time_manager_mod)
!  generates the real variable time_since_ae from the time_type
!  input time, and then calls daily_mean_solar with this real input
!  (option 2).  the third and fourth options correspond to the first
!  and second and are used with then spectral 2-layer model, where
!  only cosz and solar are desired as output. these routines generate
!  dummy arguments and then call option 2, where the calculation is
!  done.
!
!  the time of year is set by
!    real, intent(in) :: time_since_ae
!  with time_type input, the time of year is extracted from
!    type(time_type), intent(in) :: time
!
!
!  separate routines exist within this interface for scalar,
!  1D or 2D input and output fields:
!
!    real, intent(in), dimension(:,:) :: lat
! OR real, intent(in), dimension(:)   :: lat
! OR real, intent(in)                 :: lat
!
!    real, intent(out), dimension(:,:) :: cosz, fracday
! OR real, intent(out), dimension(:)   :: cosz, fracday
! OR real, intent(out)                 :: cosz, fracday
!
!--------------------------------------------------------------------
!
!  intent(in) variables:
!
!     lat            latitudes of model grid points
!                    [ radians ]
!     time_since_ae  time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!                    [ radians ]
!     time           time at which astronomical values are desired
!                    time_type variable [ seconds, days]
!
!
!  intent(out) variables:
!
!     cosz           cosz is cosine of an effective zenith angle that
!                    would produce the correct daily solar flux if the
!                    sun were fixed at that position for the period of
!                    daylight.
!                    should one also, or instead, compute cosz weighted
!                    by the instantaneous flux, averaged over the day ??
!                    [ dimensionless ]
!     fracday        daylight fraction of time interval
!                    [ dimensionless ]
!     rrsun          earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!                    [ dimensionless ]
!     solar          shortwave flux factor: cosine of zenith angle *
!                    daylight fraction / (earth-sun distance squared)
!                    [ dimensionless ]
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


! <SUBROUTINE NAME="daily_mean_solar_2d">
!  <OVERVIEW>
!    daily_mean_solar_2d computes the daily mean astronomical
!    parameters for the input points at latitude lat and time of year
!    time_since_ae.
!  </OVERVIEW>
!  <DESCRIPTION>
!    daily_mean_solar_2d computes the daily mean astronomical
!    parameters for the input points at latitude lat and time of year
!    time_since_ae.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call daily_mean_solar_2d (lat, time_since_ae, cosz, h_out, rr_out)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="h_out" TYPE="real">
!   2-d array of half-day lengths at the latitudes
!  </OUT>
!  <OUT NAME="rr_out" TYPE="real">
!   the inverse of the square of the earth-sun
!    distance relative to the mean distance at angle ang in the earth's
!    orbit.
!  </OUT>
! </SUBROUTINE>
!
subroutine daily_mean_solar_2d (lat, time_since_ae, cosz, h_out, rr_out)

!----------------------------------------------------------------------
!    daily_mean_solar_2d computes the daily mean astronomical
!    parameters for the input points at latitude lat and time of year
!    time_since_ae.
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
real, dimension(:,:), intent(in)   :: lat
real,                 intent(in)   :: time_since_ae
real, dimension(:,:), intent(out)  :: cosz, h_out
real,                 intent(out)  :: rr_out
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      real, dimension(size(lat,1),size(lat,2)) :: h
      real :: ang, dec, rr

!--------------------------------------------------------------------
!    be sure the time in the annual cycle is legitimate.
!---------------------------------------------------------------------
      if (time_since_ae < 0.0 .or. time_since_ae > twopi) &
        call error_mesg('astronomy_mod', &
                        'time_since_ae not between 0 and 2pi', FATAL)

!---------------------------------------------------------------------
!    define the orbital angle (location in year), solar declination,
!    half-day length and earth sun distance factor. use functions
!    contained in this module.
!---------------------------------------------------------------------
      ang = angle (time_since_ae)
      dec = declination(ang)
      h   = half_day    (lat, dec)
      rr  = r_inv_squared (ang)

!---------------------------------------------------------------------
!    where the entire day is dark, define cosz to be zero. otherwise
!    use the standard formula. define the daylight fraction and earth-
!    sun distance.
!---------------------------------------------------------------------
        cosz = sin(lat)*sin(dec)*(h)/(PI) + cos(lat)*cos(dec)*sin(h)/(PI)

        where (cosz < 0.)
        cosz = 0.
        end where
        
      h_out = h/PI
      rr_out = rr
!--------------------------------------------------------------------

end subroutine daily_mean_solar_2d



!#####################################################################
! <SUBROUTINE NAME="daily_mean_solar_1d">
!  <OVERVIEW>
!    daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!    and calls daily_mean_solar_2d. on return, the 2d fields are
!    returned to the original 1d fields.
!  </OVERVIEW>
!  <DESCRIPTION>
!    daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!    and calls daily_mean_solar_2d. on return, the 2d fields are
!    returned to the original 1d fields.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call daily_mean_solar_1d (lat, time_since_ae, cosz, h_out, rr_out)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="h_out" TYPE="real">
!   2-d array of half-day lengths at the latitudes
!  </OUT>
!  <OUT NAME="rr_out" TYPE="real">
!   the inverse of the square of the earth-sun
!    distance relative to the mean distance at angle ang in the earth's
!    orbit.
!  </OUT>
! </SUBROUTINE>
!
subroutine daily_mean_solar_1d (lat, time_since_ae, cosz, h_out, rr_out)

!--------------------------------------------------------------------
!    daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!    and calls daily_mean_solar_2d. on return, the 2d fields are
!    returned to the original 1d fields.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
real, intent(in), dimension(:) :: lat
real, intent(in) :: time_since_ae
real, intent(out), dimension(size(lat(:))) ::        cosz
real, intent(out), dimension(size(lat(:)))           :: h_out
real, intent(out)           :: rr_out
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables

      real, dimension(size(lat),1) :: lat_2d, cosz_2d, hout_2d

!--------------------------------------------------------------------
!    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat

!--------------------------------------------------------------------
!    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_2d (lat_2d, time_since_ae, cosz_2d,      &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!    place output fields into 1-d arguments for return to
!    calling routine.
!-------------------------------------------------------------------
      h_out = hout_2d(:,1)
      cosz  = cosz_2d(:,1)


end subroutine daily_mean_solar_1d


!######################################################################
! <SUBROUTINE NAME="daily_mean_solar_2level">
!  <OVERVIEW>
!    daily_mean_solar_2level takes 1-d input fields, adds a second
!    dimension and calls daily_mean_solar_2d. on return, the 2d fields
!    are returned to the original 1d fields.
!  </OVERVIEW>
!  <DESCRIPTION>
!    daily_mean_solar_2level takes 1-d input fields, adds a second
!    dimension and calls daily_mean_solar_2d. on return, the 2d fields
!    are returned to the original 1d fields.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call daily_mean_solar_2level (lat, time_since_ae, cosz, solar)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="solar" TYPE="real">
!   shortwave flux factor: cosine of zenith angle *
!                    daylight fraction / (earth-sun distance squared)
!  </OUT>
! </SUBROUTINE>
!
subroutine daily_mean_solar_2level (lat, time_since_ae, cosz, solar)

!--------------------------------------------------------------------
!    daily_mean_solar_2level takes 1-d input fields, adds a second
!    dimension and calls daily_mean_solar_2d. on return, the 2d fields
!    are returned to the original 1d fields.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
real, intent(in), dimension(:)          :: lat
real, intent(in)                        :: time_since_ae
real, intent(out), dimension(size(lat(:))) :: cosz, solar
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables

      real, dimension(size(lat),1) :: lat_2d, cosz_2d, hout_2d
      real                         :: rr_out

!--------------------------------------------------------------------
!    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat

!--------------------------------------------------------------------
!    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_2d (lat_2d, time_since_ae, cosz_2d,      &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!    place output fields into 1-d arguments for return to
!    calling routine.
!-------------------------------------------------------------------
      solar = cosz_2d(:,1)*hout_2d(:,1)*rr_out
      cosz  = cosz_2d(:,1)


end subroutine daily_mean_solar_2level



!####################################################################
! <SUBROUTINE NAME="daily_mean_solar_0d">
!  <OVERVIEW>
!    daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!    and calls daily_mean_solar_2d. on return, the 2d fields are
!    returned to the original 1d fields.
!  </OVERVIEW>
!  <DESCRIPTION>
!    daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!    and calls daily_mean_solar_2d. on return, the 2d fields are
!    returned to the original 1d fields.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call daily_mean_solar_0d (lat, time_since_ae, cosz, h_out, rr_out)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time_since_ae" TYPE="real">
!   time of year; autumnal equinox = 0.0,
!                    one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="h_out" TYPE="real">
!   2-d array of half-day lengths at the latitudes
!  </OUT>
!  <OUT NAME="rr_out" TYPE="real">
!   the inverse of the square of the earth-sun
!    distance relative to the mean distance at angle ang in the earth's
!    orbit.
!  </OUT>
! </SUBROUTINE>
!
subroutine daily_mean_solar_0d (lat, time_since_ae, cosz, h_out, rr_out)

!--------------------------------------------------------------------
!    daily_mean_solar_1d takes 1-d input fields, adds a second dimension
!    and calls daily_mean_solar_2d. on return, the 2d fields are
!    returned to the original 1d fields.
!----------------------------------------------------------------------

real, intent(in)         :: lat, time_since_ae
real, intent(out)        :: cosz, h_out, rr_out

!--------------------------------------------------------------------
!   local variables

      real, dimension(1,1) :: lat_2d, cosz_2d, hout_2d

!--------------------------------------------------------------------
!    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d = lat

!--------------------------------------------------------------------
!    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_2d (lat_2d, time_since_ae, cosz_2d,     &
                                hout_2d, rr_out)

!-------------------------------------------------------------------
!    return output fields to scalars for return to calling routine.
!-------------------------------------------------------------------
      h_out = hout_2d(1,1)
      cosz  = cosz_2d(1,1)


end subroutine daily_mean_solar_0d

!####################################################################
! <SUBROUTINE NAME="daily_mean_solar_cal_2d">
!  <OVERVIEW>
!    daily_mean_solar_cal_2d receives time_type inputs, converts
!    them to real variables and then calls daily_mean_solar_2d to
!    compute desired astronomical variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    daily_mean_solar_cal_2d receives time_type inputs, converts
!    them to real variables and then calls daily_mean_solar_2d to
!    compute desired astronomical variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call daily_mean_solar_cal_2d (lat, time, cosz, fracday, rrsun)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time" TYPE="time_type">
!   time of year; autumnal equinox = 0.0, one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
! </SUBROUTINE>
!
subroutine daily_mean_solar_cal_2d (lat, time, cosz, fracday, rrsun)

!-------------------------------------------------------------------
!    daily_mean_solar_cal_2d receives time_type inputs, converts
!    them to real variables and then calls daily_mean_solar_2d to
!    compute desired astronomical variables.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
real, dimension(:,:), intent(in)  :: lat
type(time_type),      intent(in)  :: time
real, dimension(:,:), intent(out) :: cosz, fracday
real,                 intent(out) :: rrsun
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables

      real :: time_since_ae

!--------------------------------------------------------------------
!    be sure the time in the annual cycle is legitimate.
!---------------------------------------------------------------------
      time_since_ae = orbital_time(time)
      if (time_since_ae < 0.0 .or. time_since_ae > twopi) &
          call error_mesg ('astronomy_mod', &
                         'time_since_ae not between 0 and 2pi', FATAL)

!--------------------------------------------------------------------
!    call daily_mean_solar_2d to calculate astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_2d (lat, time_since_ae, cosz,        &
                                fracday, rrsun)


end subroutine daily_mean_solar_cal_2d



!#####################################################################
! <SUBROUTINE NAME="daily_mean_solar_cal_1d">
!  <OVERVIEW>
!    daily_mean_solar_cal_1d receives time_type inputs, converts
!    them to real, 2d variables and then calls daily_mean_solar_2d to
!    compute desired astronomical variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    daily_mean_solar_cal_1d receives time_type inputs, converts
!    them to real, 2d variables and then calls daily_mean_solar_2d to
!    compute desired astronomical variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call daily_mean_solar_cal_1d (lat, time, cosz, fracday, rrsun)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time" TYPE="time_type">
!   time of year; autumnal equinox = 0.0, one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
! </SUBROUTINE>
!
subroutine daily_mean_solar_cal_1d (lat, time, cosz, fracday, rrsun)

!-------------------------------------------------------------------
!    daily_mean_solar_cal_1d receives time_type inputs, converts
!    them to real, 2d variables and then calls daily_mean_solar_2d to
!    compute desired astronomical variables.
!-------------------------------------------------------------------

real, dimension(:),  intent(in)   :: lat
type(time_type),     intent(in)   :: time
real, dimension(:),  intent(out)  :: cosz, fracday
real,                intent(out)  :: rrsun

!---------------------------------------------------------------------
!  local variables

      real, dimension(size(lat),1) :: lat_2d, cosz_2d, fracday_2d


!--------------------------------------------------------------------
!    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat

!--------------------------------------------------------------------
!    call daily_mean_solar_cal_2d to convert the time_types to reals and
!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_cal_2d (lat_2d, time, cosz_2d,   &
                                    fracday_2d, rrsun)

!-------------------------------------------------------------------
!    place output fields into 1-d arguments for return to
!    calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(:,1)
      cosz  = cosz_2d(:,1)


end subroutine daily_mean_solar_cal_1d


!###################################################################
! <SUBROUTINE NAME="daily_mean_solar_cal_2level">
!  <OVERVIEW>
!    daily_mean_solar_cal_2level receives 1d arrays and time_type input,
!    converts them to real, 2d variables and then calls
!    daily_mean_solar_2d to compute desired astronomical variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    daily_mean_solar_cal_2level receives 1d arrays and time_type input,
!    converts them to real, 2d variables and then calls
!    daily_mean_solar_2d to compute desired astronomical variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call daily_mean_solar_cal_2level (lat, time, cosz, solar)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time" TYPE="time_type">
!   time of year; autumnal equinox = 0.0, one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="solar" TYPE="real">
!   shortwave flux factor: cosine of zenith angle *
!                    daylight fraction / (earth-sun distance squared)
!  </OUT>
! </SUBROUTINE>
!
subroutine daily_mean_solar_cal_2level (lat, time, cosz, solar)

!-------------------------------------------------------------------
!    daily_mean_solar_cal_2level receives 1d arrays and time_type input,
!    converts them to real, 2d variables and then calls
!    daily_mean_solar_2d to compute desired astronomical variables.
!-------------------------------------------------------------------

real, dimension(:),  intent(in)   :: lat
type(time_type),     intent(in)   :: time
real, dimension(:),  intent(out)  :: cosz, solar

!---------------------------------------------------------------------
!  local variables

      real, dimension(size(lat),1) :: lat_2d, cosz_2d, fracday_2d
      real                         :: rrsun


!--------------------------------------------------------------------
!    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d(:,1) = lat

!--------------------------------------------------------------------
!    call daily_mean_solar_cal_2d to convert the time_types to reals and
!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_cal_2d (lat_2d, time, cosz_2d,   &
                                    fracday_2d, rrsun)

!-------------------------------------------------------------------
!    place output fields into 1-d arguments for return to
!    calling routine.
!-------------------------------------------------------------------
      solar = cosz_2d(:,1)*fracday_2d(:,1)*rrsun
      cosz  = cosz_2d(:,1)


end subroutine daily_mean_solar_cal_2level




!###################################################################
! <SUBROUTINE NAME="daily_mean_solar_cal_0d">
!  <OVERVIEW>
!    daily_mean_solar_cal_0d converts scalar input fields to real,
!    2d variables and then calls daily_mean_solar_2d to compute desired
!    astronomical variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    daily_mean_solar_cal_0d converts scalar input fields to real,
!    2d variables and then calls daily_mean_solar_2d to compute desired
!    astronomical variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call daily_mean_solar_cal_0d (lat, time, cosz, fracday, rrsun)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time" TYPE="time_type">
!   time of year; autumnal equinox = 0.0, one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
! </SUBROUTINE>
!
subroutine daily_mean_solar_cal_0d (lat, time, cosz, fracday, rrsun)

!-------------------------------------------------------------------
!    daily_mean_solar_cal_0d converts scalar input fields to real,
!    2d variables and then calls daily_mean_solar_2d to compute desired
!    astronomical variables.
!-------------------------------------------------------------------

!--------------------------------------------------------------------
real,             intent(in)  :: lat
type(time_type),  intent(in)  :: time
real,             intent(out) :: cosz, fracday, rrsun
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables

      real, dimension(1,1) :: lat_2d, cosz_2d, fracday_2d

!--------------------------------------------------------------------
!    define 2-d versions of input data array.
!--------------------------------------------------------------------
      lat_2d = lat

!--------------------------------------------------------------------
!    call daily_mean_solar_cal_2d to convert the time_types to reals and
!    then calculate the astronomy fields.
!--------------------------------------------------------------------
      call daily_mean_solar_cal_2d (lat_2d, time, cosz_2d,           &
                                    fracday_2d, rrsun)

!-------------------------------------------------------------------
!    place output fields into scalar arguments for return to
!    calling routine.
!-------------------------------------------------------------------
      fracday = fracday_2d(1,1)
      cosz  = cosz_2d(1,1)


end subroutine daily_mean_solar_cal_0d



!####################################################################

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                END INTERFACE DAILY_MEAN_SOLAR
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   INTERFACE ANNUAL_MEAN_SOLAR
!
!  call annual_mean_solar (js, je, lat, cosz, solar, fracday, rrsun)
!     or
!  call annual_mean_solar (lat, cosz, solar)
!
!  the second interface above is used by the spectral 2-layer model,
!  which requires only cosz and solar as output arguments, and which
!  makes this call during the initialization phase of the model.

!  separate routines exist within this interface for 1D or 2D input
!  and output fields:
!
!    real, intent(in), dimension(:,:) :: lat
! OR real, intent(in), dimension(:)   :: lat
!
!    real, intent(out), dimension(:,:) :: cosz, solar, fracday
! OR real, intent(out), dimension(:)   :: cosz, solar, fracday
!
!---------------------------------------------------------------------
!  intent(in) variables:
!
!     js, je         starting/ ending subdomain j indices of data in
!                    the physics wiondow being integrated
!     lat            latitudes of model grid points
!                    [ radians ]
!
!
!  intent(out) variables:
!
!     cosz           cosz is the average over the year of the cosine of
!                    an effective zenith angle that would produce the
!                    correct daily solar flux if the sun were fixed at
!                    that single position for the period of daylight on
!                    the given day. in this average, the daily mean
!                    effective cosz is weighted by the daily mean solar
!                    flux.
!                    [ dimensionless ]
!     solar          normalized solar flux, averaged over the year,
!                    equal to the product of fracday*cosz*rrsun
!                    [ dimensionless ]
!     fracday        daylight fraction calculated so as to make the
!                    average flux (solar) equal to the product of the
!                    flux-weighted avg cosz * this fracday * assumed
!                    annual mean avg earth-sun distance of 1.0.
!                    [ dimensionless ]
!     rrsun          annual mean earth-sun distance (r) relative to
!                    semi-major axis of orbital ellipse (a); assumed
!                    to be 1.0
!                    [ dimensionless ]
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





!--------------------------------------------------------------
! <SUBROUTINE NAME="annual_mean_solar_2d">
!  <OVERVIEW>
!    annual_mean_solar_2d returns 2d fields of annual mean values of
!    the cosine of zenith angle, daylight fraction and earth-sun
!    distance at the specified latitude.
!  </OVERVIEW>
!  <DESCRIPTION>
!    annual_mean_solar_2d returns 2d fields of annual mean values of
!    the cosine of zenith angle, daylight fraction and earth-sun
!    distance at the specified latitude.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call annual_mean_solar_2d (js, je, lat, cosz, solar, fracday,  &
!                                 rrsun)
!  </TEMPLATE>
!  <IN NAME="js, je" TYPE="real">
!   Starting/ending index of latitude window
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time" TYPE="time_type">
!   time of year; autumnal equinox = 0.0, one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="solar" TYPE="real">
!   shortwave flux factor: cosine of zenith angle *
!                    daylight fraction / (earth-sun distance squared)
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
! </SUBROUTINE>
!
subroutine annual_mean_solar_2d (js, je, lat, cosz, solar, fracday,  &
                                 rrsun)

!---------------------------------------------------------------------
!    annual_mean_solar_2d returns 2d fields of annual mean values of
!    the cosine of zenith angle, daylight fraction and earth-sun
!    distance at the specified latitude.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                 intent(in)    :: js, je
real, dimension(:,:),    intent(in)    :: lat
real, dimension(:,:),    intent(out)   :: solar, cosz, fracday
real,                    intent(out)   :: rrsun
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables

      real, dimension(size(lat,1),size(lat,2)) :: s,z
      real    :: t
      integer :: n, i

!--------------------------------------------------------------------
!    if the calculation has not yet been done, do it here.
!--------------------------------------------------------------------
      if (.not. annual_mean_calculated) then

!----------------------------------------------------------------------
!    determine annual mean values of solar flux and product of cosz
!    and solar flux by integrating the annual cycle in num_angles
!    orbital increments.
!----------------------------------------------------------------------
        solar = 0.0
        cosz = 0.0
        do n =1, num_angles
          t = float(n-1)*twopi/float(num_angles)
          call daily_mean_solar(lat,t, z, fracday, rrsun)
          s = z*rrsun*fracday
          solar = solar + s
          cosz  = cosz  + z*s
        end do
        solar = solar/float(num_angles)
        cosz  = cosz/float(num_angles)

!--------------------------------------------------------------------
!   define the flux-weighted annual mean cosine of the zenith angle.
!--------------------------------------------------------------------
        where(solar.eq.0.0)
          cosz = 0.0
        elsewhere
          cosz = cosz/solar
        end where

!-------------------------------------------------------------------
!    define avg fracday such as to make the avg flux (solar) equal to
!    the product of the avg cosz * avg fracday * assumed mean avg
!    radius of 1.0. it is unlikely that these avg fracday and avg rr
!    will ever be used.
!--------------------------------------------------------------------
        where(solar  .eq.0.0)
          fracday = 0.0
        elsewhere
          fracday = solar/cosz
        end where
        rrsun = 1.00

!---------------------------------------------------------------------
!    save the values that have been calculated as module variables, if
!    those variables are present; i.e., not the spectral 2-layer model.
!---------------------------------------------------------------------
        if (allocated (cosz_ann)) then
          cosz_ann(js:je) = cosz(1,:)
          solar_ann(js:je)   = solar(1,:)
          fracday_ann(js:je) = fracday(1,:)
          rrsun_ann = rrsun

!--------------------------------------------------------------------
!    increment the points computed counter. set flag to end execution
!    once values have been calculated for all points owned by the
!    processor.
!--------------------------------------------------------------------
          num_pts = num_pts + size(lat,1)*size(lat,2)
          if ( num_pts == total_pts)  annual_mean_calculated = .true.
        endif

!--------------------------------------------------------------------
!    if the calculation has been done, return the appropriate module
!    variables.
!--------------------------------------------------------------------
      else
        if (allocated (cosz_ann)) then
          do i=1, size(lat,1)
            cosz(i,:)    = cosz_ann(js:je)
            solar(i,:)   = solar_ann(js:je)
            fracday(i,:) = fracday_ann(js:je)
          end do
          rrsun = rrsun_ann
        endif
      endif

!----------------------------------------------------------------------


end subroutine annual_mean_solar_2d



!#####################################################################
! <SUBROUTINE NAME="annual_mean_solar_1d">
!  <OVERVIEW>
!    annual_mean_solar_1d creates 2-d input fields from 1-d input fields
!    and then calls annual_mean_solar_2d to obtain 2-d output fields
!    which are then stored into 1-d fields for return to the calling
!    subroutine.
!  </OVERVIEW>
!  <DESCRIPTION>
!    annual_mean_solar_1d creates 2-d input fields from 1-d input fields
!    and then calls annual_mean_solar_2d to obtain 2-d output fields
!    which are then stored into 1-d fields for return to the calling
!    subroutine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call annual_mean_solar_1d (jst, jnd, lat, cosz, solar,  &
!                                 fracday, rrsun_out)
!  </TEMPLATE>
!  <IN NAME="jst, jnd" TYPE="real">
!   Starting/ending index of latitude window
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <IN NAME="time" TYPE="time_type">
!   time of year; autumnal equinox = 0.0, one year = 2 * pi
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="solar" TYPE="real">
!   shortwave flux factor: cosine of zenith angle *
!                    daylight fraction / (earth-sun distance squared)
!  </OUT>
!  <OUT NAME="fracday" TYPE="real">
!   daylight fraction of time interval
!  </OUT>
!  <OUT NAME="rrsun_out" TYPE="real">
!   earth-sun distance (r) relative to semi-major axis
!                    of orbital ellipse (a) : (a/r)**2
!  </OUT>
! </SUBROUTINE>
!
subroutine annual_mean_solar_1d (jst, jnd, lat, cosz, solar,  &
                                 fracday, rrsun_out)

!---------------------------------------------------------------------
!    annual_mean_solar_1d creates 2-d input fields from 1-d input fields
!    and then calls annual_mean_solar_2d to obtain 2-d output fields
!    which are then stored into 1-d fields for return to the calling
!    subroutine.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
integer,            intent(in)     :: jst, jnd
real, dimension(:), intent(in)     :: lat(:)
real, dimension(:), intent(out)    :: cosz, solar,  fracday
real,               intent(out)    :: rrsun_out
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      real, dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d,   &
                                      fracday_2d
      real :: rrsun

!--------------------------------------------------------------------
!    if the calculation has not been done, do it here.
!--------------------------------------------------------------------
      if ( .not. annual_mean_calculated) then

!--------------------------------------------------------------------
!    define 2-d versions of input data array.
!--------------------------------------------------------------------
        lat_2d(:,1) = lat

!--------------------------------------------------------------------
!    call annual_mean_solar_2d to calculate the astronomy fields.
!--------------------------------------------------------------------
        call annual_mean_solar_2d (jst, jnd, lat_2d, cosz_2d,   &
                                   solar_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!    place output fields into 1-D arrays for return to calling routine.
!-------------------------------------------------------------------
        fracday = fracday_2d(:,1)
        rrsun_out = rrsun
        solar = solar_2d(:,1)
        cosz  =  cosz_2d(:,1)

!--------------------------------------------------------------------
!    if the calculation has been done, simply return the module
!    variables contain the results at the desired latitudes.
!--------------------------------------------------------------------
      else
        cosz(:)    = cosz_ann(jst:jnd)
        solar(:)   = solar_ann(jst:jnd)
        fracday(:) = fracday_ann(jst:jnd)
        rrsun      = rrsun_ann
      endif

end subroutine annual_mean_solar_1d



!####################################################################
! <SUBROUTINE NAME="annual_mean_solar_2level">
!  <OVERVIEW>
!    annual_mean_solar_2level creates 2-d input fields from 1-d input
!    fields and then calls annual_mean_solar_2d to obtain 2-d output
!    fields which are then stored into 1-d fields for return to the
!    calling subroutine. this subroutine will be called during model
!    initialization.
!  </OVERVIEW>
!  <DESCRIPTION>
!    annual_mean_solar_2level creates 2-d input fields from 1-d input
!    fields and then calls annual_mean_solar_2d to obtain 2-d output
!    fields which are then stored into 1-d fields for return to the
!    calling subroutine. this subroutine will be called during model
!    initialization.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call annual_mean_solar_2level (lat, cosz, solar)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   latitudes of model grid points
!  </IN>
!  <OUT NAME="cosz" TYPE="real">
!   cosine of solar zenith angle
!  </OUT>
!  <OUT NAME="solar" TYPE="real">
!   shortwave flux factor: cosine of zenith angle *
!                    daylight fraction / (earth-sun distance squared)
!  </OUT>
! </SUBROUTINE>
!
subroutine annual_mean_solar_2level (lat, cosz, solar)

!---------------------------------------------------------------------
!    annual_mean_solar_2level creates 2-d input fields from 1-d input
!    fields and then calls annual_mean_solar_2d to obtain 2-d output
!    fields which are then stored into 1-d fields for return to the
!    calling subroutine. this subroutine will be called during model
!    initialization.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:), intent(in)     :: lat
real, dimension(:), intent(out)    :: cosz, solar
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      real, dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d,   &
                                      fracday_2d
      integer :: jst, jnd
      real    :: rrsun

!--------------------------------------------------------------------
!    if the calculation has not been done, do it here.
!--------------------------------------------------------------------
      if ( .not. annual_mean_calculated) then

!--------------------------------------------------------------------
!    define 2-d versions of input data array.
!--------------------------------------------------------------------
        lat_2d(:,1) = lat
        jst = 1
        jnd = size(lat(:))

!--------------------------------------------------------------------
!    call annual_mean_solar_2d to calculate the astronomy fields.
!--------------------------------------------------------------------
        call annual_mean_solar_2d (jst, jnd, lat_2d, cosz_2d,   &
                                   solar_2d, fracday_2d, rrsun)

!-------------------------------------------------------------------
!    place output fields into 1-D arrays for return to calling routine.
!-------------------------------------------------------------------
        solar = solar_2d(:,1)
        cosz  =  cosz_2d(:,1)

!--------------------------------------------------------------------
!    if the calculation has been done, print an error message since
!    this subroutine should be called only once.
!--------------------------------------------------------------------
      else
        call error_mesg ('astronomy_mod', &
            'annual_mean_solar_2level should be called only once', &
                                                                 FATAL)
      endif
      annual_mean_calculated = .true.

end subroutine annual_mean_solar_2level


!####################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                END INTERFACE ANNUAL_MEAN_SOLAR
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!###################################################################
! <SUBROUTINE NAME="astronomy_end">
!  <OVERVIEW>
!    astronomy_init is the destructor for astronomy_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    astronomy_init is the destructor for astronomy_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call astronomy_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine astronomy_end

!----------------------------------------------------------------------
!    astronomy_end is the destructor for astronomy_mod.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    check if the module has been initialized.
!----------------------------------------------------------------------
      if (.not. module_is_initialized)   &
                call error_mesg ( 'astronomy_mod',  &
                         ' module has not been initialized', FATAL)

!----------------------------------------------------------------------
!    deallocate module variables.
!----------------------------------------------------------------------
      deallocate (orb_angle)
      if (allocated(cosz_ann) ) then
        deallocate (cosz_ann)
        deallocate (fracday_ann)
        deallocate (solar_ann)
      endif

!----------------------------------------------------------------------
!    mark the module as uninitialized.
!----------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------


end subroutine astronomy_end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################
! <SUBROUTINE NAME="orbit">
!  <OVERVIEW>
!    orbit computes and stores a table of value of orbital angles as a
!    function of orbital time (both the angle and time are zero at
!    autumnal equinox in the NH, and range from 0 to 2*pi).
!  </OVERVIEW>
!  <DESCRIPTION>
!    orbit computes and stores a table of value of orbital angles as a
!    function of orbital time (both the angle and time are zero at
!    autumnal equinox in the NH, and range from 0 to 2*pi).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call orbit
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine orbit

!---------------------------------------------------------------------
!    orbit computes and stores a table of value of orbital angles as a
!    function of orbital time (both the angle and time are zero at
!    autumnal equinox in the NH, and range from 0 to 2*pi).
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      integer :: n
      real    :: d1, d2, d3, d4, d5, dt, norm

!--------------------------------------------------------------------
!    allocate the orbital angle array, sized by the namelist parameter
!    num_angles, defining the annual cycle resolution of the earth's
!    orbit. define some constants to be used.
!--------------------------------------------------------------------
! wfc moving to astronomy_init
!     allocate ( orb_angle(0:num_angles) )
      orb_angle(0) = 0.0
      dt = twopi/float(num_angles)
      norm = sqrt(1.0 - ecc**2)
      dt = dt*norm

!---------------------------------------------------------------------
!    define the orbital angle at each of the num_angles locations in
!    the orbit.
!---------------------------------------------------------------------
      do n = 1, num_angles
        d1 = dt*r_inv_squared(orb_angle(n-1))
        d2 = dt*r_inv_squared(orb_angle(n-1)+0.5*d1)
        d3 = dt*r_inv_squared(orb_angle(n-1)+0.5*d2)
        d4 = dt*r_inv_squared(orb_angle(n-1)+d3)
        d5 = d1/6.0 + d2/3.0 +d3/3.0 +d4/6.0
        orb_angle(n) = orb_angle(n-1) + d5
      end do

!-------------------------------------------------------------------



end subroutine orbit



!###################################################################
! <FUNCTION NAME="r_inv_squared">
!  <OVERVIEW>
!    r_inv_squared returns the inverse of the square of the earth-sun
!    distance relative to the mean distance at angle ang in the earth's
!    orbit.
!  </OVERVIEW>
!  <DESCRIPTION>
!    r_inv_squared returns the inverse of the square of the earth-sun
!    distance relative to the mean distance at angle ang in the earth's
!    orbit.
!  </DESCRIPTION>
!  <TEMPLATE>
!    r = r_inv_squared (ang)
!  </TEMPLATE>
!  <IN NAME="ang" TYPE="real">
!    angular position of earth in its orbit, relative to a
!            value of 0.0 at the NH autumnal equinox, value between
!            0.0 and 2 * pi
!  </IN>
! </FUNCTION>
!
function r_inv_squared (ang)

!--------------------------------------------------------------------
!    r_inv_squared returns the inverse of the square of the earth-sun
!    distance relative to the mean distance at angle ang in the earth's
!    orbit.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
real, intent(in) :: ang
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!
!  intent(in) variables:
!
!      ang   angular position of earth in its orbit, relative to a
!            value of 0.0 at the NH autumnal equinox, value between
!            0.0 and 2 * pi
!            [ radians ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      real :: r_inv_squared, r, rad_per

!---------------------------------------------------------------------
!  local variables:
!
!      r_inv_squared    the inverse of the square of the earth-sun
!                       distance relative to the mean distance
!                       [ dimensionless ]
!      r                earth-sun distance relative to mean distance
!                       [ dimensionless ]
!      rad_per          angular position of perihelion
!                       [ radians ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    define the earth-sun distance (r) and then return the inverse of
!    its square (r_inv_squared) to the calling routine.
!--------------------------------------------------------------------
      rad_per       = per*deg_to_rad
      r             = (1. - ecc**2)/(1. + ecc*cos(ang - rad_per))
      r_inv_squared = r**(-2)


end function r_inv_squared




!####################################################################
! <FUNCTION NAME="angle">
!  <OVERVIEW>
!    angle determines the position within the earth's orbit at time t
!    in the year ( t = 0 at NH autumnal equinox.) by interpolating
!    into the orbital position table.
!  </OVERVIEW>
!  <DESCRIPTION>
!    angle determines the position within the earth's orbit at time t
!    in the year ( t = 0 at NH autumnal equinox.) by interpolating
!    into the orbital position table.
!  </DESCRIPTION>
!  <TEMPLATE>
!    r = angle (t)
!  </TEMPLATE>
!  <IN NAME="t" TYPE="real">
!    time of year (between 0 and 2*pi; t=0 at NH autumnal
!                equinox
!  </IN>
! </FUNCTION>
!
function angle (t)

!--------------------------------------------------------------------
!    angle determines the position within the earth's orbit at time t
!    in the year ( t = 0 at NH autumnal equinox.) by interpolating
!    into the orbital position table.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
real, intent(in) :: t
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!
!  intent(in) variables:
!
!         t      time of year (between 0 and 2*pi; t=0 at NH autumnal
!                equinox
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      real :: angle, norm_time, x
      integer :: int, int_1

!--------------------------------------------------------------------
!   local variables:
!
!     angle       orbital position relative to NH autumnal equinox
!                 [ radians ]
!     norm_time   index into orbital table corresponding to input time
!                 [ dimensionless ]
!     x           fractional distance between the orbital table entries
!                 bracketing the input time
!                 [ dimensionless ]
!     int         table index which is lower than actual position, but
!                 closest to it
!                 [ dimensionless ]
!     int_1       next table index just larger than actual orbital
!                 position
!                 [ dimensionless ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    define orbital tables indices bracketing current orbital time
!    (int and int_1). define table index distance between the lower
!    table value (int) and the actual orbital time (x). define orbital
!    position as being  x of the way between int and int_1. renormalize
!    angle to be within the range 0 to 2*pi.
!--------------------------------------------------------------------
      norm_time = t*float(num_angles)/twopi
      int = floor(norm_time)
      int = modulo(int,num_angles)
      int_1 = int+1
      x = norm_time - floor(norm_time)
      angle = (1.0 -x)*orb_angle(int) + x*orb_angle(int_1)
      angle = modulo(angle, twopi)

end function angle

!####################################################################
! <FUNCTION NAME="declination">
!  <OVERVIEW>
!    declination returns the solar declination angle at orbital
!    position ang in earth's orbit.
!  </OVERVIEW>
!  <DESCRIPTION>
!    declination returns the solar declination angle at orbital
!    position ang in earth's orbit.
!  </DESCRIPTION>
!  <TEMPLATE>
!    r =  declination (ang)
!  </TEMPLATE>
!  <IN NAME="ang" TYPE="real">
!     solar orbital position ang in earth's orbit
!  </IN>
! </FUNCTION>
!
function declination (ang)

!--------------------------------------------------------------------
!    declination returns the solar declination angle at orbital
!    position ang in earth's orbit.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
real, intent(in) :: ang
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      real :: declination
      real :: rad_obliq, sin_dec

!--------------------------------------------------------------------
!   local variables:
!
!    declination         solar declination angle
!                        [ radians ]
!    rad_obliq           obliquity of the ecliptic
!                        [ radians ]
!    sin_dec             sine of the solar declination
!                        [ dimensionless ]
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    compute the solar declination.
!---------------------------------------------------------------------
      rad_obliq   =   obliq*deg_to_rad
      sin_dec     = - sin(rad_obliq)*sin(ang)
      declination =   asin(sin_dec)


end function declination


!#####################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                 INTERFACE HALF_DAY
!
! half_day (latitude, dec) result (h)
!
!
!  separate routines exist within this interface for scalar,
!  or 2D input and output fields:
!
!    real, intent(in), dimension(:,:) :: latitude
! OR real, intent(in)                 :: latitude
!
!    real, dimension(size(latitude,1),size(latitude,2))  :: h
! OR real                                                :: h
!
!--------------------------------------------------------------------
!
!  intent(in) variables:
!
!     latitude       latitudes of model grid points
!                    [ radians ]
!     dec            solar declination
!                    [ radians ]
!
!  intent(out) variables:
!
!     h              half of the length of daylight at the given
!                    latitude and orbital position (dec); value
!                    ranges between 0 (all darkness) and pi (all
!                    daylight)
!                    [ dimensionless ]
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! <FUNCTION NAME="half_day_2d">
!  <OVERVIEW>
!    half_day_2d returns a 2-d array of half-day lengths at the
!    latitudes and declination provided.
!  </OVERVIEW>
!  <DESCRIPTION>
!    half_day_2d returns a 2-d array of half-day lengths at the
!    latitudes and declination provided.
!  </DESCRIPTION>
!  <TEMPLATE>
!    h = half_day_2d (latitude, dec)
!  </TEMPLATE>
!  <IN NAME="latitude" TYPE="real">
!   latitutde of view point
!  </IN>
!  <IN NAME="dec" TYPE="real">
!   solar declination angle at view point
!  </IN>
! </FUNCTION>
!
function half_day_2d (latitude, dec) result(h)

!--------------------------------------------------------------------
!    half_day_2d returns a 2-d array of half-day lengths at the
!    latitudes and declination provided.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:,:), intent(in)                     :: latitude
real,                 intent(in)                     :: dec
real, dimension(size(latitude,1),size(latitude,2))   :: h
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!   local variables

      real, dimension (size(latitude,1),size(latitude,2)):: &
                                                  cos_half_day, lat
      real :: tan_dec
      real :: eps = 1.0E-05

!---------------------------------------------------------------------
!   local variables
!
!     cos_half_day       cosine of half-day length
!                        [ dimensionless ]
!     lat                model latitude, adjusted so that it is never
!                        0.5*pi or -0.5*pi
!     tan_dec            tangent of solar declination
!                        [ dimensionless ]
!     eps                small increment
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    define tangent of the declination.
!--------------------------------------------------------------------
      tan_dec = tan(dec)

!--------------------------------------------------------------------
!    adjust latitude so that its tangent will be defined.
!--------------------------------------------------------------------
      lat = latitude
      where (latitude ==  0.5*PI) lat= latitude - eps
      where (latitude == -0.5*PI) lat= latitude + eps

!--------------------------------------------------------------------
!    define the cosine of the half-day length. adjust for cases of
!    all daylight or all night.
!--------------------------------------------------------------------
      cos_half_day = -tan(lat)*tan_dec
      where (cos_half_day <= -1.0)  h = PI
      where (cos_half_day >= +1.0)  h = 0.0
      where(cos_half_day > -1.0 .and. cos_half_day < 1.0) &
                                               h = acos(cos_half_day)


end function half_day_2d

!---------------------------------------------------------------
! <FUNCTION NAME="half_day_0d">
!  <OVERVIEW>
!    half_day_0d takes scalar input fields, makes them into 2-d fields
!    dimensioned (1,1), and calls half_day_2d. on return, the 2-d
!    fields are converted to the desired scalar output.
!  </OVERVIEW>
!  <DESCRIPTION>
!    half_day_0d takes scalar input fields, makes them into 2-d fields
!    dimensioned (1,1), and calls half_day_2d. on return, the 2-d
!    fields are converted to the desired scalar output.
!  </DESCRIPTION>
!  <TEMPLATE>
!    h = half_day_2d (latitude, dec)
!  </TEMPLATE>
!  <IN NAME="latitude" TYPE="real">
!   latitutde of view point
!  </IN>
!  <IN NAME="dec" TYPE="real">
!   solar declination angle at view point
!  </IN>
! </FUNCTION>
!
function half_day_0d(latitude, dec) result(h)

!--------------------------------------------------------------------
!    half_day_0d takes scalar input fields, makes them into 2-d fields
!    dimensioned (1,1), and calls half_day_2d. on return, the 2-d
!    fields are converted to the desired scalar output.
!--------------------------------------------------------------------

real, intent(in) :: latitude, dec
real             :: h

!----------------------------------------------------------------------
!  local variables

      real, dimension(1,1) :: lat_2d, h_2d

!---------------------------------------------------------------------
!    create 2d array from the input latitude field.
!---------------------------------------------------------------------
      lat_2d = latitude

!---------------------------------------------------------------------
!    call half_day with the 2d arguments to calculate half-day length.
!---------------------------------------------------------------------
      h_2d = half_day (lat_2d, dec)

!---------------------------------------------------------------------
!    create scalar from 2d array.
!---------------------------------------------------------------------
      h = h_2d(1,1)



end function half_day_0d


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                 END INTERFACE HALF_DAY
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!####################################################################
! <FUNCTION NAME="orbital_time">
!  <OVERVIEW>
!    orbital time returns the time (1 year = 2*pi) since autumnal
!    equinox
!  </OVERVIEW>
!  <DESCRIPTION>
!    orbital time returns the time (1 year = 2*pi) since autumnal
!    equinox; autumnal_eq_ref is a module variable of time_type and
!    will have been defined by default or by a call to
!    set_ref_date_of_ae; length_of_year is available through the time
!    manager and is set at the value approriate for the calandar being
!    used
!  </DESCRIPTION>
!  <TEMPLATE>
!    t = orbital_time(time)
!  </TEMPLATE>
!  <IN NAME="time" TYPE="time_type">
!   time (1 year = 2*pi) since autumnal equinox
!  </IN>
! </FUNCTION>
!
function orbital_time(time) result(t)

!---------------------------------------------------------------------
!    orbital time returns the time (1 year = 2*pi) since autumnal
!    equinox; autumnal_eq_ref is a module variable of time_type and
!    will have been defined by default or by a call to
!    set_ref_date_of_ae; length_of_year is available through the time
!    manager and is set at the value approriate for the calandar being
!    used
!---------------------------------------------------------------------

type(time_type), intent(in) :: time
real                        :: t

!--------------------------------------------------------------------
      t = real ( (time - autumnal_eq_ref)//period_time_type)
      t = twopi*(t - floor(t))
      if (time < autumnal_eq_ref) t = twopi - t


end function orbital_time


!#####################################################################
! <FUNCTION NAME="universal_time">
!  <OVERVIEW>
!    universal_time returns the time of day at longitude = 0.0
!    (1 day = 2*pi)
!  </OVERVIEW>
!  <DESCRIPTION>
!    universal_time returns the time of day at longitude = 0.0
!    (1 day = 2*pi)
!  </DESCRIPTION>
!  <TEMPLATE>
!    t = universal_time(time)
!  </TEMPLATE>
!  <IN NAME="time" TYPE="time_type">
!   time (1 year = 2*pi) since autumnal equinox
!  </IN>
! </FUNCTION>
!
function universal_time(time) result(t)

!--------------------------------------------------------------------
!    universal_time returns the time of day at longitude = 0.0
!    (1 day = 2*pi)
!--------------------------------------------------------------------

type(time_type), intent(in) :: time
real                        :: t

!--------------------------------------------------------------------
!   local variables

      integer ::  seconds, days

      call get_time (time, seconds, days)
      t = twopi*real(seconds)/seconds_per_day

end function universal_time

subroutine diurnal_exoplanet(lat, lon, Time, coszen, fracsun, rrsun, &
                             dt, allow_negative_cosz, half_day_out)

  ! Returns the diurnal phase on an 'exoplanet'.
  ! An exoplanet is assumed to have no calendar, and may have an arbitrary
  ! length day and year, determined by its rotation rate and orbital period
  ! defined in constants_mod.
  ! Allows for a planet and it's orbital direction to be arbitrary.

  real, dimension(:,:), intent(in)           :: lat, lon
  type(time_type),      intent(in)           :: Time
  real, dimension(:,:), intent(out)          :: coszen, fracsun
  real,                 intent(out)          :: rrsun
  real,                 intent(in), optional :: dt
  logical,              intent(in), optional :: allow_negative_cosz
  real, dimension(:,:), intent(out), optional :: half_day_out

  integer :: seconds, year_in_s
  real :: r_seconds, day_in_s, gmt, time_since_ae, frac_of_day, frac_of_year, substellar_lon

  call get_time(Time, seconds)
  substellar_lon = (orbital_rate - omega)*seconds
  ! negative sign because time relative to gmt = -substellar_lon position
  gmt = modulo(-substellar_lon, 2.*pi)
  ! if ( mpp_pe() == mpp_root_pe() ) then
  !      write (*,*) seconds, substellar_lon, gmt
  ! endif
  frac_of_year = modulo(orbital_rate*seconds, 1.0)
  time_since_ae = frac_of_year * 2.0 * pi
  ! call get_time(length_of_year(), year_in_s)
  ! r_seconds = real(seconds)
  ! day_in_s = length_of_day()


  !frac_of_day = r_seconds / day_in_s
  !frac_of_year = r_seconds / year_in_s

  !gmt = abs(mod(frac_of_day, 1.0)) * 2.0 * pi

  !time_since_ae = modulo(frac_of_year,1.0) * 2.0 * pi
  call diurnal_solar_2d(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun,  &
                     dt, allow_negative_cosz, half_day_out)

end subroutine diurnal_exoplanet

!#####################################################################





                   end module astronomy_mod