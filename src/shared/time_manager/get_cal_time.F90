!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module get_cal_time_mod

!   <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!     fms
!   </CONTACT>
!   <OVERVIEW>
!      Given a time increment as a real number, and base time and calendar
!      as a character strings, returns time as a time_type variable.
!   </OVERVIEW>

use          fms_mod, only: error_mesg, FATAL, write_version_number, lowercase, &
                            open_namelist_file, check_nml_error, stdlog, close_file, &
                            mpp_pe, mpp_root_pe

use time_manager_mod, only: time_type, operator(+), operator(-), set_time, get_time, &
                            NO_CALENDAR, THIRTY_DAY_MONTHS, NOLEAP, JULIAN, GREGORIAN, &
                            set_calendar_type, get_calendar_type, set_date, &
                            get_date, days_in_month, valid_calendar_types
use mpp_mod,          only: input_nml_file

implicit none
private

public :: get_cal_time

logical :: module_is_initialized=.false. ! This module is initialized on
                                         ! the first call to get_cal_time
                                         ! because there is no constructor.
! <NAMELIST NAME="get_cal_time_nml">
! <DATA NAME="allow_calendar_conversion" TYPE="logical"  DEFAULT=".true.">
!   This sets the default value of the optional argument named "permit_calendar_conversion" of get_cal_time.
!   This namelist is deprecated as of the memphis release.
!   If calendar conversion is not desired, then it is recommended that permit_calendar_conversion
!   be present in the call to get_cal_time and that it be set to .false.
! </DATA>

logical :: allow_calendar_conversion=.true.

namelist / get_cal_time_nml / allow_calendar_conversion
! </NAMELIST>

character(len=128) :: version='$Id: get_cal_time.F90,v 19.0 2012/01/06 22:06:10 fms Exp $'
character(len=128) :: tagname='$Name: siena_201211 $'

contains
!------------------------------------------------------------------------
! <FUNCTION NAME="get_cal_time">
!   <TEMPLATE>
!     get_cal_time(time_increment, units, calendar, permit_calendar_conversion)
!   </TEMPLATE>
!   <IN NAME="time_increment" TYPE="real"> A time interval.</IN>
!   <IN NAME="units" TYPE="character">
!
! Examples of acceptable values of units:
!
! 'days since 1980-01-01 00:00:00',
! 'hours since 1980-1-1 0:0:0',
! 'minutes since 0001-4-12'
!
! The first word in the string must be
! 'years', 'months', 'days', 'hours', 'minutes' or 'seconds'.
! The second word must be 'since'
!
! year number must occupy 4 spaces.
! Number of months, days, hours, minutes, seconds may occupy 1 or 2 spaces
! year, month and day must be separated by a '-'
! hour, minute, second must be separated by a ':'
! hour, minute, second are optional. If not present then zero is assumed.
!
! Because months are not equal increments of time, and, for julian calendar,
! neither are years, the 'years since' and 'month since' cases deserve
! further explaination.
!
! When 'years since' is used:
! The year number is increased by floor(time_increment)   to obtain a time T1.
! The year number is increased by floor(time_increment)+1 to obtain a time T2.
! The time returned is T1 + (time_increment-floor(time_increment))*(T2-T1).
!
! When 'months since' is used:
! The month number is increased by floor(time_increment). If it falls outside
! to range 1 to 12 then it is adjusted along with the year number to convert
! to a valid date. The number of days in the month of this date is used to
! compute the time interval of the fraction.
! That is:
! The month number is increased by floor(time_increment) to obtain a time T1.
! delt = the number of days in the month in which T1 falls.
! The time returned is T1 + ((time_increment-floor(time_increment))*delt.
! Two of the consequences of this scheme should be kept in mind.
! -- The time since should not be from the 29'th to 31'st of a month,
!    since an invalid date is likely to result, triggering an error stop.
! -- When time since is from the begining of a month, the fraction of a month
!    will never advance into the month after that which results from only
!    the whole number.
!
! When NO_CALENDAR is in effect, units attribute must specify a starting
! day and second, with day number appearing first
!
! Example: 'days since 100 0' Indicates 100 days 0 seconds
! </IN>
!
! <IN NAME="calendar" TYPE="character">
! Acceptable values of calendar are:
! 'noleap'
! '365_day'
! '360_day'
! 'julian'
! 'thirty_day_months'
! 'no_calendar'
! </IN>
!
! <IN NAME="permit_calendar_conversion" TYPE="logical, optional" DEFAULT="allow_calendar_conversion">
! It is sometimes desirable to allow the value of the intent(in) argument
! "calendar" to be different than the calendar in use by time_manager_mod.
! If this is not desirable, then the optional variable "permit_calendar_conversion"
! should be set to .false. so as to allow an error check.
! When calendar conversion is done, the time returned is the time in the
! time_manager's calendar, but corresponds to the date computed using the input calendar.
! For example, suppose the time_manager is using the julian calendar and
! the values of the input arguments of get_cal_time are:
! time_increment = 59.0
! units = 'days since 1980-1-1 00:00:00'
! calendar = 'noleap'
! Because it will use the noleap calendar to calculate the date, get_cal_time will return
! value of time for midnight March 1 1980, but it will be time in the julian calendar
! rather than the noleap calendar. It will never return a value of time corresponding
! to anytime during the day Feb 29.
!
! Another example:
! Suppose the time_manager is using either the noleap or julian calendars,
! and the values of the input arguments are:
! time_increment = 30.0
! units = 'days since 1980-1-1'
! calendar = 'thirty_day_months'
! In this case get_cal_time will return the value of time for Feb 1 1980 00:00:00,
! but in the time_manager's calendar.

! Calendar conversion may result in a fatal error when the input calendar type is
! a calendar that has more days per year than that of the time_manager's calendar.
! For example, if the input calendar type is julian and the time_manager's calendar
! is thirty_day_months, then get_cal_time will try to convert Jan 31 to a time in
! the thirty_day_months calendar, resulting in a fatal error.

! Note: this option was originally coded to allow noleap calendar as input when
! the julian calendar was in effect by the time_manager.
! </IN>
!
!---------------------------------------------------------------------------------------------

function get_cal_time(time_increment, units, calendar, permit_calendar_conversion)
real, intent(in) :: time_increment
character(len=*), intent(in) :: units
character(len=*), intent(in) :: calendar
logical, intent(in), optional :: permit_calendar_conversion
type(time_type) :: get_cal_time
integer :: year, month, day, hour, minute, second
integer :: i1, i2, i3, i4, i5, i6, increment_seconds, increment_days, increment_years, increment_months
real    :: month_fraction
integer :: calendar_tm_i, calendar_in_i, namelist_unit, ierr, io, logunit
logical :: correct_form
character(len=32) :: calendar_in_c
character(len=64) :: err_msg
character(len=4) :: formt='(i )'
type(time_type) :: base_time, base_time_plus_one_yr, base_time_plus_one_mo
real :: dt
logical :: permit_conversion_local

if(.not.module_is_initialized) then
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, get_cal_time_nml, iostat=io)
      ierr = check_nml_error (io, 'get_cal_time_nml')
#else
  namelist_unit = open_namelist_file()
  ierr=1
  do while (ierr /= 0)
    read(namelist_unit, nml=get_cal_time_nml, iostat=io, end=20)
    ierr = check_nml_error (io, 'get_cal_time_nml')
  enddo
  20 call close_file (namelist_unit)
#endif

  call write_version_number (version, tagname)
  logunit = stdlog()
  if(mpp_pe() == mpp_root_pe()) write (logunit, nml=get_cal_time_nml)
  module_is_initialized = .true.
endif

if(present(permit_calendar_conversion)) then
  permit_conversion_local = permit_calendar_conversion
else
  permit_conversion_local = allow_calendar_conversion
endif

calendar_in_c = lowercase(trim(cut0(calendar)))

correct_form = (trim(calendar_in_c)) == 'noleap'     .or. (trim(calendar_in_c)) == '365_day' .or. &
               (trim(calendar_in_c)) == '360_day'    .or. (trim(calendar_in_c)) == 'julian'  .or. &
               (trim(calendar_in_c)) == 'no_calendar'.or. (trim(calendar_in_c)) == 'thirty_day_months' .or. &
               (trim(calendar_in_c)) == 'gregorian'

if(.not.correct_form) then
  call error_mesg('get_cal_time','"'//trim(calendar_in_c)//'"'// &
   ' is not an acceptable calendar attribute. acceptable calendars are: '// &
   ' noleap, 365_day, 360_day, julian, no_calendar, thirty_day_months, gregorian',FATAL)
endif

calendar_tm_i = get_calendar_type()

if(.not.permit_conversion_local) then
  correct_form = (trim(calendar_in_c) == 'noleap'            .and. calendar_tm_i == NOLEAP)            .or. &
                 (trim(calendar_in_c) == '365_day'           .and. calendar_tm_i == NOLEAP)            .or. &
                 (trim(calendar_in_c) == '360_day'           .and. calendar_tm_i == THIRTY_DAY_MONTHS) .or. &
                 (trim(calendar_in_c) == 'thirty_day_months' .and. calendar_tm_i == THIRTY_DAY_MONTHS) .or. &
                 (trim(calendar_in_c) == 'julian'            .and. calendar_tm_i == JULIAN)            .or. &
                 (trim(calendar_in_c) == 'no_calendar'       .and. calendar_tm_i == NO_CALENDAR)       .or. &
                 (trim(calendar_in_c) == 'gregorian'         .and. calendar_tm_i == GREGORIAN)
  if(.not.correct_form) then
    call error_mesg('get_cal_time','calendar not consistent with calendar type in use by time_manager.'// &
         ' calendar='//trim(calendar_in_c)//'. Type in use by time_manager='//valid_calendar_types(calendar_tm_i),FATAL)
  endif
endif

if (permit_conversion_local) then
    select case (trim(calendar_in_c))
    case ('noleap')
        calendar_in_i = NOLEAP
    case ('365_day')
        calendar_in_i = NOLEAP
    case ('360_day')
        calendar_in_i = THIRTY_DAY_MONTHS
    case ('thirty_day_months')
        calendar_in_i = THIRTY_DAY_MONTHS
    case ('julian')
        calendar_in_i = JULIAN
    case ('no_calendar')
        calendar_in_i = NO_CALENDAR
    case ('gregorian')
        calendar_in_i = GREGORIAN
    case default
        call error_mesg('get_cal_time', &
                 trim(calendar_in_c)//' is an invalid calendar type (specified in call to get_cal_time)',FATAL)
    end select
else
    calendar_in_i = calendar_tm_i
end if

correct_form = lowercase(units(1:10)) == 'days since'    .or. &
               lowercase(units(1:11)) == 'hours since'   .or. &
               lowercase(units(1:13)) == 'minutes since' .or. &
               lowercase(units(1:13)) == 'seconds since'

if(calendar_in_i /= NO_CALENDAR) then
  correct_form = correct_form .or. &
               lowercase(units(1:11)) == 'years since'   .or. &
               lowercase(units(1:12)) == 'months since'
endif

if(.not.correct_form) then
  call error_mesg('get_cal_time',trim(units)//' is an invalid string for units.' // &
        ' units must begin with a time unit then the word "since"' // &
        ' Valid time units are: "seconds" "minutes", "hours", "days", and, ' // &
        ' except when NO_CALENDAR is in effect, "months" and "years"',FATAL)
endif

if(calendar_in_i /= calendar_tm_i) then
! switch to calendar type specified as input argument,
! will switch back before returning.
  call set_calendar_type(calendar_in_i)
endif

! index(string, substring[,back])
! Returns the starting position of substring as a substring of string,
! or zero if it does not occur as a substring. Default value of back is
! .false. If back is .false., the starting position of the first such
! substring is returned. If back is .true., the starting position of the
! last such substring is returned.
! Returns zero if substring is not a substring of string (regardless of value of back)

i1 = index(units,'since') + 5
if(calendar_in_i == NO_CALENDAR) then
  base_time = set_time(units(i1:len_trim(units)))
else
  base_time = set_date(units(i1:len_trim(units)))
endif

if(lowercase(units(1:10)) == 'days since') then
  increment_days = floor(time_increment)
  increment_seconds = 86400*(time_increment - increment_days)
else if(lowercase(units(1:11)) == 'hours since') then
  increment_days = floor(time_increment/24)
  increment_seconds = 86400*(time_increment/24 - increment_days)
else if(lowercase(units(1:13)) == 'minutes since') then
  increment_days = floor(time_increment/1440)
  increment_seconds = 86400*(time_increment/1440 - increment_days)
else if(lowercase(units(1:13)) == 'seconds since') then
  increment_days = floor(time_increment/86400)
  increment_seconds = 86400*(time_increment/86400 - increment_days)
else if(lowercase(units(1:11)) == 'years since') then
! The time period between between (base_time + time_increment) and
! (base_time + time_increment + 1 year) may be 360, 365, or 366 days.
! This must be determined to handle time increments with year fractions.
  call get_date(base_time, year,month,day,hour,minute,second)
  base_time             = set_date(year+floor(time_increment)  ,month,day,hour,minute,second)
  base_time_plus_one_yr = set_date(year+floor(time_increment)+1,month,day,hour,minute,second)
  call get_time(base_time_plus_one_yr - base_time, second, day)
  dt = (day*86400+second)*(time_increment-floor(time_increment))
  increment_days = floor(dt/86400)
  increment_seconds = dt - increment_days*86400
else if(lowercase(units(1:12)) == 'months since') then
  month_fraction = time_increment - floor(time_increment)
  increment_years  = floor(time_increment/12)
  increment_months = floor(time_increment) - 12*increment_years
  call get_date(base_time, year,month,day,hour,minute,second)
  base_time = set_date(year+increment_years,month+increment_months  ,day,hour,minute,second)
  dt = 86400*days_in_month(base_time) * month_fraction
  increment_days = floor(dt/86400)
  increment_seconds = dt - increment_days*86400
else
  call error_mesg('get_cal_time','"'//trim(units)//'"'//' is not an acceptable units attribute of time.'// &
    ' It must begin with: "years since", "months since", "days since", "hours since", "minutes since", or "seconds since"',FATAL)
endif

if (calendar_in_i /= calendar_tm_i) then
    if(calendar_in_i == NO_CALENDAR .or. calendar_tm_i == NO_CALENDAR) then
      call error_mesg('get_cal_time','Cannot do calendar conversion because input calendar is '// &
       trim(valid_calendar_types(calendar_in_i))//' and time_manager is using '//trim(valid_calendar_types(calendar_tm_i))// &
       ' Conversion cannot be done if either is NO_CALENDAR',FATAL)
    endif
    call get_date(base_time,year, month, day, hour, minute, second)
    get_cal_time = set_date(year,month,day,hour,minute,second) + set_time(increment_seconds, increment_days)
    call get_date(get_cal_time,year,month,day,hour,minute,second)
    call set_calendar_type(calendar_tm_i)
    get_cal_time = set_date(year,month,day,hour,minute,second, err_msg=err_msg)
    if(err_msg /= '') then
      call error_mesg('get_cal_time','Error in function get_cal_time: '//trim(err_msg)// &
                      ' Note that the time_manager is using the '//trim(valid_calendar_types(calendar_tm_i))//' calendar '// &
                      'while the calendar type passed to function get_cal_time is '//calendar_in_c,FATAL)
    endif
else
    get_cal_time = base_time + set_time(increment_seconds, increment_days)
endif

end function get_cal_time
! </FUNCTION>
!------------------------------------------------------------------------
function cut0(string)
character(len=256) :: cut0
character(len=*), intent(in) :: string
integer :: i

cut0 = string

do i=1,len(string)
  if(ichar(string(i:i)) == 0 ) then
    cut0(i:i) = ' '
  endif
enddo

return
end function cut0
!------------------------------------------------------------------------
end module get_cal_time_mod
