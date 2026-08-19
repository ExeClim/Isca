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

module time_manager_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   fms
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   A software package that provides a set of simple interfaces for
!   modelers to perform computations related to time and dates.
! </OVERVIEW>

! <DESCRIPTION>
!    The changes between the lima revision and this revision are more
!    extensive that all those between antwerp and lima.
!    A brief description of these changes follows.
!
!    1) Added option to set the smallest time increment to something less than one second.
!       This is controlled by calling the pubic subroutine set_ticks_per_second.
!
!    2) Gregorian calendar fixed.
!
!    3) Optional error flag added to calling arguments of public routines.
!       This allows the using routine to terminate the program. It is likely that more
!       diagnostic information is available from the user than from time_manager alone.
!       If the error flag is present then it is the responsibility of the using
!       routine to test it and add additional information to the error message.
!
!    4) Removed the restriction that time increments be positive in routines that increment or decrement
!       time and date. The option to prohibit negative increments can be turned on via optional argument.
!
!    5) subroutine set_date_c modified to handle strings that include only hours or only hours and minutes.
!       This complies with CF convensions.
!
!    6) Made calendar specific routines private.
!       They are not used, and should not be used, by any using code.
!
!    7) Error messages made more informative.
!
!    The module defines a type that can be used to represent discrete
!    times (accurate to one second) and to map these times into dates
!    using a variety of calendars. A time is mapped to a date by
!    representing the time with respect to an arbitrary base date (refer
!    to <B>NOTES</B> section for the <LINK SRC="#base date">base date</LINK> setting).
!
!    The time_manager provides a single defined type, time_type, which is
!    used to store time and date quantities. A time_type is a positive
!    definite quantity that represents an interval of time. It can be
!    most easily thought of as representing the number of seconds in some
!    time interval. A time interval can be mapped to a date under a given
!    calendar definition by using it to represent the time that has passed
!    since some base date. A number of interfaces are provided to operate
!    on time_type variables and their associated calendars. Time intervals
!    can be as large as n days where n is the largest number represented by
!    the default integer type on a compiler. This is typically considerably
!    greater than 10 million years (assuming 32 bit integer representation)
!    which is likely to be adequate for most applications. The description
!    of the interfaces is separated into two sections. The first deals with
!    operations on time intervals while the second deals with operations
!    that convert time intervals to dates for a given calendar.

!    The smallest increment of time is referred to as a tick.
!    A tick cannot be larger than 1 second, which also is the default.
!    The number of ticks per second is set via pubic subroutine set_ticks_per_second.
!    For example, ticks_per_second = 1000  will set the tick to one millisecond.
! </DESCRIPTION>

! <DATA NAME="time_type" TYPE="derived type">
!    Derived-type data variable used to store time and date quantities. It
!    contains three PRIVATE variables: days, seconds and ticks.
! </DATA>

use constants_mod, only: rseconds_per_day=>seconds_per_day, seconds_per_sol, orbital_period
use fms_mod, only: error_mesg, FATAL, WARNING, write_version_number, stdout

implicit none
private

! Module defines a single type
public time_type

! Operators defined on time_type
public operator(+),  operator(-),   operator(*),   operator(/),  &
       operator(>),  operator(>=),  operator(==),  operator(/=), &
       operator(<),  operator(<=),  operator(//),  assignment(=)

! Subroutines and functions operating on time_type
public set_time, increment_time, decrement_time, get_time, interval_alarm
public repeat_alarm, time_type_to_real, real_to_time_type

! List of available calendar types
public    THIRTY_DAY_MONTHS,    JULIAN,    GREGORIAN,  NOLEAP, NO_CALENDAR, INVALID_CALENDAR

! Subroutines and functions involving relations between time and calendar
public set_calendar_type
public get_calendar_type
public set_ticks_per_second
public get_ticks_per_second
public set_date
public get_date
public increment_date
public decrement_date
public days_in_month
public leap_year
public length_of_year
public length_of_day
public days_in_year
public month_name

public valid_calendar_types

! Subroutines for printing version number and time type
public :: time_manager_init, print_time, print_date

! The following exist only for interpolator.F90
! interpolator.F90 uses them to do a calendar conversion,
! which is also done by get_cal_time. interpolator.F90
! should be modified to use get_cal_time instead.
! After interpolator.F90 is fixed, these can be removed
! and the corresponding private routines can be renamed.
! (e.g., rename set_date_julian_private to be just set_date_julian)
public :: set_date_julian, set_date_no_leap, get_date_julian, get_date_no_leap

public :: date_to_string

!====================================================================

! Global data to define calendar type
integer, parameter :: THIRTY_DAY_MONTHS = 1,      JULIAN = 2, &
                      GREGORIAN = 3,              NOLEAP = 4, &
                      NO_CALENDAR = 0,  INVALID_CALENDAR =-1
integer, private :: calendar_type = NO_CALENDAR
integer, parameter :: max_type = 4

! Define number of days per month
integer, private :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer, parameter :: seconds_per_day = rseconds_per_day  ! This should automatically cast real to integer
integer, parameter :: days_in_400_year_period = 146097    ! Used only for gregorian
integer, dimension(days_in_400_year_period) :: coded_date ! Used only for gregorian
integer, dimension(400,12,31) :: date_to_day              ! Used only for gregorian
integer, parameter :: invalid_date=-1                     ! Used only for gregorian

! time_type is implemented as seconds and days to allow for larger intervals
type time_type
   private
   integer:: seconds
   integer:: days
   integer:: ticks
   integer:: dummy ! added as a workaround bug on IRIX64 (AP)
end type time_type

!======================================================================

interface operator (+);   module procedure time_plus;        end interface
interface operator (-);   module procedure time_minus;       end interface
interface operator (*);   module procedure time_scalar_mult
                          module procedure scalar_time_mult; end interface
interface operator (/);   module procedure time_scalar_divide
                          module procedure time_divide;      end interface
interface operator (>);   module procedure time_gt;          end interface
interface operator (>=);  module procedure time_ge;          end interface
interface operator (<);   module procedure time_lt;          end interface
interface operator (<=);  module procedure time_le;          end interface
interface operator (==);  module procedure time_eq;          end interface
interface operator (/=);  module procedure time_ne;          end interface
interface operator (//);  module procedure time_real_divide; end interface
interface assignment(=);  module procedure time_assignment;  end interface

!======================================================================

interface set_time
  module procedure set_time_i, set_time_c
end interface

interface set_date
  module procedure set_date_i, set_date_c
end interface

!======================================================================

character(len=128) :: version='$Id: time_manager.F90,v 19.0 2012/01/06 22:06:12 fms Exp $'
character(len=128) :: tagname='$Name: siena_201211 $'
logical :: module_is_initialized = .false.

!======================================================================

!  A tick is the smallest increment of time.
!  That is, smallest increment of time = (1/ticks_per_second) seconds

integer :: ticks_per_second = 1

!======================================================================
contains

! First define all operations on time intervals independent of calendar

!=========================================================================
! <FUNCTION NAME="set_time">

!   <OVERVIEW>
!     Given some number of seconds and days, returns the
!     corresponding time_type.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Given some number of seconds and days, returns the
!     corresponding time_type. set_time has two forms;
!     one accepts integer input, the other a character string.
!     For the first form, there are no restrictions on the range of the inputs,
!     except that the result must be positive time.
!     e.g. days=-1, seconds=86401 is acceptable.
!     For the second form, days and seconds must both be positive.
!   </DESCRIPTION>
!   <TEMPLATE>
!     1. set_time(seconds, days, ticks, err_msg)
!   </TEMPLATE>
!   <TEMPLATE>
!     2. set_time(time_string, err_msg, allow_rounding)
!   </TEMPLATE>

!   <IN NAME="seconds" UNITS="" TYPE="integer" DIM="(scalar)">
!     A number of seconds.
!   </IN>
!   <IN NAME="days" UNITS="" TYPE="integer" DIM="(scalar)">
!     A number of days.
!   </IN>
!   <IN NAME="ticks" UNITS="" TYPE="integer, optional" DIM="(scalar)">
!     A number of ticks.
!   </IN>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
!   <IN NAME="time_string" TYPE="character">
!     Contains days and seconds separated by a single blank.
!     days must be integer, seconds may be integer or real.
!     Examples: '100 43200'  '100 43200.50'
!   </IN>
!   <IN NAME="allow_rounding"   TYPE="logical, optional" DEFAULT=".true.">
!     When .true., any fractions of a second will be rounded off to the nearest tick.
!     When .false., it is a fatal error if the second fraction cannot be exactly
!     represented by a number of ticks.
!   </IN>
!   <OUT NAME="set_time" UNITS="" TYPE="" DIM="" DEFAULT="">
!     A time interval corresponding to this number of days and seconds.
!   </OUT>

 function set_time_private(seconds, days, ticks, Time_out, err_msg)

! Returns a time interval corresponding to this number of days, seconds, and ticks.
! days, seconds and ticks may be negative, but resulting time must be positive.

! -- pjp --
! To understand why inputs may be negative,
! one needs to understand the intrinsic function "modulo".
! The expanation below is copied from a web page on fortran 90

! In addition, CEILING, FLOOR  and MODULO  have been added to Fortran 90.
! Only the last one is difficult to explain, which is most easily done with the examples from ISO (1991)

! MOD (8,5)    gives  3     MODULO (8,5)    gives  3
! MOD (-8,5)   gives -3     MODULO (-8,5)   gives  2
! MOD (8,-5)   gives  3     MODULO (8,-5)   gives -2
! MOD (-8,-5)  gives -3     MODULO (-8,-5)  gives -3

! I don't think it is difficult to explain.
! I think that is it sufficient to say this:
! "The result of modulo(n,m) has the sign of m"
! -- pjp --

 logical                       :: set_time_private
 integer, intent(in)           :: seconds, days, ticks
 type(time_type),  intent(out) :: Time_out
 character(len=*), intent(out) :: err_msg
 integer            :: seconds_new, days_new, ticks_new

 seconds_new = seconds + floor(ticks/real(ticks_per_second))
 ticks_new = modulo(ticks,ticks_per_second)
 days_new = days + floor(seconds_new/real(seconds_per_day))
 seconds_new = modulo(seconds_new,seconds_per_day)

 if ( seconds_new < 0 .or. ticks_new < 0) then
   call error_mesg('function set_time_i','Bad result for time. Contact those responsible for maintaining time_manager',FATAL)
 endif

 if(days_new < 0) then
   write(err_msg,'(a,i6,a,i6,a,i6)') 'time is negative. days=',days_new,' seconds=',seconds_new,' ticks=',ticks_new
   set_time_private = .false.
 else
   Time_out%days = days_new
   Time_out%seconds = seconds_new
   Time_out%ticks = ticks_new
   err_msg = ''
   set_time_private = .true.
 endif

 end function set_time_private
!---------------------------------------------------------------------------

 function set_time_i(seconds, days, ticks, err_msg)
 type(time_type)               :: set_time_i
 integer, intent(in)           :: seconds
 integer, intent(in), optional :: days, ticks
 character(len=*), intent(out), optional :: err_msg
 character(len=128) :: err_msg_local
 integer            :: odays, oticks

 if(.not.module_is_initialized) call time_manager_init

 odays  = 0; if(present(days))  odays  = days
 oticks = 0; if(present(ticks)) oticks = ticks
 if(present(err_msg)) err_msg = ''

 if(.not.set_time_private(seconds, odays, oticks, set_time_i, err_msg_local)) then
   if(error_handler('function set_time_i', trim(err_msg_local), err_msg)) return
 endif

 end function set_time_i
!---------------------------------------------------------------------------

 function set_time_c(string, err_msg, allow_rounding)

 type(time_type) :: set_time_c
 character(len=*), intent(in) :: string
 character(len=*), intent(out), optional :: err_msg
 logical, intent(in), optional :: allow_rounding

 character(len=4) :: formt='(i )'
 integer :: i1, i2, i3, day, second, tick, nsps
 character(len=32) :: string_sifted_left
 character(len=128) :: err_msg_local
 logical :: allow_rounding_local

 if(.not.module_is_initialized) call time_manager_init
 if(present(err_msg)) err_msg = ''
 allow_rounding_local=.true.; if(present(allow_rounding)) allow_rounding_local=allow_rounding

 err_msg_local = 'Form of character time stamp is incorrect. The character time stamp is: '//trim(string)

 string_sifted_left = adjustl(string)
 i1 = index(trim(string_sifted_left),' ')
 if(i1 == 0) then
   if(error_handler('function set_time_c', err_msg_local, err_msg)) return
 endif
 if(index(string,'-') /= 0 .or. index(string,':') /= 0) then
   if(error_handler('function set_time_c', err_msg_local, err_msg)) return
 endif

 i2 = index(trim(string_sifted_left),'.')
 i3 = len_trim(cut0(string_sifted_left))

 if(i2 /= 0) then ! There is no decimal point
 ! Check that decimal is on seconds (not days)
   if(i2 < i1) then
     if(error_handler('function set_time_c', err_msg_local, err_msg)) return
   endif
 endif
 write(formt(3:3),'(i1)') i1-1
 read(string_sifted_left(1:i1-1),formt) day

 if(i2 == 0) then ! There is no decimal point
   write(formt(3:3),'(i1)') i3-i1
   read(string_sifted_left(i1+1:i3),formt) second
   tick = 0
 else ! There is a decimal point
 ! nsps = spaces occupied by whole number of seconds
   nsps = i2-i1-1
   if(nsps == 0) then
     second = 0
   else
     write(formt(3:3),'(i1)') nsps
     read(string_sifted_left(i1+1:i2-1),formt) second
   endif

   if(.not.get_tick_from_string(string_sifted_left(i2:i3), err_msg_local, allow_rounding_local, tick)) then
     if(error_handler('function set_time_c', err_msg_local, err_msg)) return
   endif
 ! If tick has been rounded up to ticks_per_second, then bump up second.
   if(tick == ticks_per_second) then
     second = second + 1
     tick = 0
   endif
 endif

 if(.not.set_time_private(second, day, tick, set_time_c, err_msg_local)) then
   if(error_handler('function set_time_c', err_msg_local, err_msg)) return
 endif

 end function set_time_c
!---------------------------------------------------------------------------
! </FUNCTION>

 function get_tick_from_string(string, err_msg, allow_rounding, tick)

 logical :: get_tick_from_string
 character(len=*), intent(in) :: string
 character(len=*), intent(out) :: err_msg
 logical, intent(in) :: allow_rounding
 integer, intent(out) :: tick

 character(len=4) :: formt='(i )'
 integer :: i3, nspf, fraction, magnitude, tpsfrac

 err_msg = ''
 get_tick_from_string = .true.
 i3 = len_trim(string)
 nspf = i3 - 1 ! nspf = spaces occupied by fractional seconds, excluding decimal point
 if(nspf == 0) then
   tick = 0 ! Nothing to the right of the decimal point
 else
   write(formt(3:3),'(i1)') nspf
   read(string(2:i3),formt) fraction
   if(fraction == 0) then
     tick = 0 ! All zeros to the right of the decimal point
   else
     magnitude = 10**nspf
     tpsfrac = ticks_per_second*fraction
     if(allow_rounding) then
       tick = nint((real(tpsfrac)/magnitude))
     else
       if(modulo(tpsfrac,magnitude) == 0) then
         tick = tpsfrac/magnitude
       else
         write(err_msg,'(a,i6)') 'Second fraction cannot be exactly represented with ticks.  '// &
                                 'fraction='//trim(string)//'  ticks_per_second=',ticks_per_second
         get_tick_from_string = .false.
       endif
     endif
   endif
 endif

 end function get_tick_from_string
!---------------------------------------------------------------------------
! <SUBROUTINE NAME="get_time">

!   <OVERVIEW>
!     Given a time interval, returns the corresponding seconds and days.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Given a time interval, returns the corresponding seconds and days.
!   </DESCRIPTION>
!   <TEMPLATE>
!     get_time(time, seconds, days, ticks, err_msg)
!   </TEMPLATE>

!   <IN NAME="time" TYPE="time_type">
!     A time interval.
!   </IN>
!   <OUT NAME="seconds" UNITS="" TYPE="integer" DIM="(scalar)">
!     A number of seconds.
!   </OUT>
!   <OUT NAME="days" UNITS="" TYPE="integer" DIM="(scalar)">
!     A number of days.
!   </OUT>
!   <OUT NAME="ticks" UNITS="" TYPE="integer, optional" DIM="(scalar)">
!     A number of ticks.
!   </OUT>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>

subroutine get_time(Time, seconds, days, ticks, err_msg)

! Returns days and seconds ( < 86400 ) corresponding to a time.

type(time_type), intent(in) :: Time
integer, intent(out) :: seconds
integer, intent(out), optional :: days, ticks
character(len=*), intent(out), optional :: err_msg
character(len=128) :: err_msg_local

if(.not.module_is_initialized) call time_manager_init
if(present(err_msg)) err_msg = ''

seconds = Time%seconds

if(present(ticks)) then
  ticks = Time%ticks
else
  if(Time%ticks /= 0) then
    err_msg_local = 'subroutine get_time: ticks must be present when time has a second fraction'
    if(error_handler('subroutine get_time', err_msg_local, err_msg)) return
  endif
endif

if (present(days)) then
  days = Time%days
else
  if (Time%days > (huge(seconds) - seconds)/seconds_per_day) then
    err_msg_local = 'Integer overflow in seconds. Optional argument days must be present.'
    if(error_handler('subroutine get_time', err_msg_local, err_msg)) return
  endif
  seconds = seconds + Time%days * seconds_per_day
endif

end subroutine get_time
! </SUBROUTINE>

!-------------------------------------------------------------------------
! <FUNCTION NAME="increment_time">

!   <OVERVIEW>
!      Given a time and an increment of days and seconds, returns
!      a time that adds this increment to an input time.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Given a time and an increment of days and seconds, returns
!      a time that adds this increment to an input time.
!      Increments a time by seconds and days.
!   </DESCRIPTION>
!   <TEMPLATE>
!     increment_time(time, seconds, days, ticks, err_msg, allow_neg_inc)
!   </TEMPLATE>

!   <IN NAME="time"  TYPE="time_type" DIM="(scalar)">
!      A time interval.
!   </IN>
!   <IN NAME="seconds"  TYPE="integer" DIM="(scalar)">
!     Increment of seconds.
!   </IN>
!   <IN NAME="days" UNITS="" TYPE="integer, optional" DIM="(scalar)">
!     Increment of days.
!   </IN>
!   <IN NAME="ticks"  TYPE="integer, optional" DIM="(scalar)">
!     Increment of ticks.
!   </IN>
!   <OUT NAME="increment_time"  TYPE="time_type" DIM="(scalar)">
!     A time that adds this increment to the input time.
!     A negative result is a fatal error.
!   </OUT>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
!   <IN NAME="allow_neg_inc" TYPE="logical, optional" DIM="(scalar)" DEFAULT=".true.">
!     When .false., it is a fatal error if any of the input time increments are negative.
!     This mimics the behavior of lima and earlier revisions.
!   </IN>

 function increment_time(Time, seconds, days, ticks, err_msg, allow_neg_inc)

! Increments a time by seconds, days and ticks.

 type(time_type)               :: increment_time
 type(time_type), intent(in)   :: Time
 integer, intent(in)           :: seconds
 integer, intent(in), optional :: days, ticks
 character(len=*), intent(out), optional :: err_msg
 logical, intent(in), optional :: allow_neg_inc

 integer :: odays, oticks
 character(len=128) :: err_msg_local
 logical :: allow_neg_inc_local

 odays  = 0; if(present(days))  odays  = days
 oticks = 0; if(present(ticks)) oticks = ticks
 allow_neg_inc_local=.true.; if(present(allow_neg_inc)) allow_neg_inc_local=allow_neg_inc

 if(.not.allow_neg_inc_local) then
   if(seconds < 0 .or. odays < 0 .or. oticks < 0) then
     write(err_msg_local,10) seconds, odays, oticks
     10 format('One or more time increments are negative: seconds=',i6,'  days=',i6,'  ticks=',i6)
     if(error_handler('function increment_time', err_msg_local, err_msg)) return
   endif
 endif

 if(.not.increment_time_private(Time, seconds, odays, oticks, increment_time, err_msg_local)) then
   if(error_handler('function increment_time', err_msg_local, err_msg)) return
 endif

 end function increment_time
! </FUNCTION>
!--------------------------------------------------------------------------

 function increment_time_private(Time_in, seconds, days, ticks, Time_out, err_msg)

! Increments a time by seconds, days and ticks.

 logical                       :: increment_time_private
 type(time_type),  intent(in)  :: Time_in
 integer,          intent(in)  :: seconds, days, ticks
 type(time_type),  intent(out) :: Time_out
 character(len=*), intent(out) :: err_msg

! Watch for immediate overflow on days or seconds
 if(days >= huge(days) - Time_in%days)  then
   err_msg = 'Integer overflow in days in increment_time'
   increment_time_private = .false.
   return
 endif
 if(seconds >= huge(seconds) - Time_in%seconds) then
   err_msg = 'Integer overflow in seconds in increment_time'
   increment_time_private = .false.
   return
 endif

 increment_time_private = set_time_private(Time_in%seconds+seconds, Time_in%days+days, Time_in%ticks+ticks, Time_out, err_msg)

 end function increment_time_private

!--------------------------------------------------------------------------
! <FUNCTION NAME="decrement_time">

!   <OVERVIEW>
!      Given a time and a decrement of days and seconds, returns
!      a time that subtracts this decrement from an input time.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Decrements a time by seconds and days.
!   </DESCRIPTION>
!   <TEMPLATE>
!     Decrement_time(time, seconds, days, ticks, err_msg, allow_neg_inc)
!   </TEMPLATE>

!   <IN NAME="time"  TYPE="time_type" DIM="(scalar)">
!      A time interval.
!   </IN>
!   <IN NAME="seconds"  TYPE="integer" DIM="(scalar)">
!     Decrement of seconds.
!   </IN>
!   <IN NAME="days"  TYPE="integer, optional" DIM="(scalar)">
!     Decrement of days.
!   </IN>
!   <IN NAME="ticks"  TYPE="integer, optional" DIM="(scalar)">
!     Decrement of ticks.
!   </IN>
!   <OUT NAME="decrement_time"  TYPE="time_type">
!      A time that subtracts this decrement from an input time.
!      A negative result is a fatal error.
!   </OUT>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
!   <IN NAME="allow_neg_inc" TYPE="logical, optional" DIM="(scalar)" DEFAULT=".true.">
!     When .false., it is a fatal error if any of the input time increments are negative.
!     This mimics the behavior of lima and earlier revisions.
!   </IN>

function decrement_time(Time, seconds, days, ticks, err_msg, allow_neg_inc)

! Decrements a time by seconds, days and ticks.

type(time_type)               :: decrement_time
type(time_type), intent(in)   :: Time
integer, intent(in)           :: seconds
integer, intent(in), optional :: days, ticks
character(len=*), intent(out), optional :: err_msg
logical, intent(in), optional :: allow_neg_inc

integer            :: odays, oticks
character(len=128) :: err_msg_local
logical :: allow_neg_inc_local

odays  = 0;  if (present(days))   odays = days
oticks = 0;  if (present(ticks)) oticks = ticks
allow_neg_inc_local=.true.; if(present(allow_neg_inc)) allow_neg_inc_local=allow_neg_inc

if(.not.allow_neg_inc_local) then
  if(seconds < 0 .or. odays < 0 .or. oticks < 0) then
    write(err_msg_local,10) seconds,odays,oticks
    10 format('One or more time increments are negative: seconds=',i6,'  days=',i6,'  ticks=',i6)
    if(error_handler('function decrement_time', err_msg_local, err_msg)) return
  endif
endif

 if(.not.increment_time_private(Time, -seconds, -odays, -oticks, decrement_time, err_msg_local)) then
   if(error_handler('function decrement_time', err_msg_local, err_msg)) return
 endif

end function decrement_time
! </FUNCTION>

!--------------------------------------------------------------------------
! <FUNCTION NAME="time_gt  operator(>)">

!   <OVERVIEW>
!      Returns true if time1 > time2.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns true if time1 > time2.
!   </DESCRIPTION>
!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!       Returns true if time1 > time2
!   </OUT>
!   <TEMPLATE>
!     time_gt(time1, time2)
!   </TEMPLATE>

function time_gt(time1, time2)

! Returns true if time1 > time2

logical :: time_gt
type(time_type), intent(in) :: time1, time2

time_gt = (time1%days > time2%days)
if(time1%days == time2%days) then
   if(time1%seconds == time2%seconds) then
      time_gt = (time1%ticks > time2%ticks)
   else
      time_gt = (time1%seconds > time2%seconds)
   endif
endif

end function time_gt
! </FUNCTION>

!--------------------------------------------------------------------------
! <FUNCTION NAME="time_ge; operator(>=)">

!   <OVERVIEW>
!      Returns true if time1 >= time2.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns true if time1 >= time2.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_ge(time1, time2)
!   </TEMPLATE>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!       Returns true if time1 >= time2
!   </OUT>

function time_ge(time1, time2)

! Returns true if time1 >= time2

logical :: time_ge
type(time_type), intent(in) :: time1, time2

time_ge = (time_gt(time1, time2) .or. time_eq(time1, time2))

end function time_ge
! </FUNCTION>

!--------------------------------------------------------------------------
! <FUNCTION NAME="time_lt; operator(<)">

!   <OVERVIEW>
!      Returns true if time1 < time2.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns true if time1 < time2.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_lt(time1, time2)
!   </TEMPLATE>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!       Returns true if time1 < time2
!   </OUT>

function time_lt(time1, time2)

! Returns true if time1 < time2

logical :: time_lt
type(time_type), intent(in) :: time1, time2

time_lt = (time1%days < time2%days)
if(time1%days == time2%days)then
   if(time1%seconds == time2%seconds) then
      time_lt = (time1%ticks < time2%ticks)
   else
      time_lt = (time1%seconds < time2%seconds)
   endif
endif
end function time_lt
! </FUNCTION>

!--------------------------------------------------------------------------
! <FUNCTION NAME="time_le; operator(<=)">

!   <OVERVIEW>
!      Returns true if time1 <= time2.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns true if time1 <= time2.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_le(time1, time2)
!   </TEMPLATE>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!       Returns true if time1 <= time2
!   </OUT>

function time_le(time1, time2)

! Returns true if time1 <= time2

logical :: time_le
type(time_type), intent(in) :: time1, time2

time_le = (time_lt(time1, time2) .or. time_eq(time1, time2))

end function time_le
! </FUNCTION>

!--------------------------------------------------------------------------
! <FUNCTION NAME="time_eq; operator(==)">

!   <OVERVIEW>
!      Returns true if time1 == time2.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns true if time1 == time2.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_eq(time1, time2)
!   </TEMPLATE>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!       Returns true if time1 == time2
!   </OUT>

function time_eq(time1, time2)

! Returns true if time1 == time2

logical :: time_eq
type(time_type), intent(in) :: time1, time2

if(.not.module_is_initialized) call time_manager_init

time_eq = (time1%seconds == time2%seconds .and. time1%days == time2%days &
     .and. time1%ticks == time2%ticks)

end function time_eq
! </FUNCTION>

!--------------------------------------------------------------------------
! <FUNCTION NAME="time_ne; operator(/=)">

!   <OVERVIEW>
!      Returns true if time1 /= time2.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns true if time1 /= time2.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_ne(time1, time2)
!   </TEMPLATE>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!       Returns true if time1 /= time2
!   </OUT>

function time_ne(time1, time2)

! Returns true if time1 /= time2

logical :: time_ne
type(time_type), intent(in) :: time1, time2

time_ne = (.not. time_eq(time1, time2))

end function time_ne
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="time_plus; operator(+)">

!   <OVERVIEW>
!       Returns sum of two time_types.
!   </OVERVIEW>
!   <TEMPLATE>
!     time1 + time2
!   </TEMPLATE>
!   <DESCRIPTION>
!       Returns sum of two time_types.
!   </DESCRIPTION>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="time_type" DIM="" DEFAULT="">
!       Returns sum of two time_types.
!   </OUT>

function time_plus(time1, time2)

! Returns sum of two time_types

type(time_type) :: time_plus
type(time_type), intent(in) :: time1, time2

if(.not.module_is_initialized) call time_manager_init

time_plus = increment_time(time1, time2%seconds, time2%days, time2%ticks)

end function time_plus
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="time_minus; operator(-)">

!   <OVERVIEW>
!       Returns difference of two time_types.
!   </OVERVIEW>
!   <DESCRIPTION>
!       Returns difference of two time_types. WARNING: a time type is positive
!       so by definition time1 - time2  is the same as time2 - time1.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_minus(time1, time2)
!   </TEMPLATE>
!   <TEMPLATE>
!     time1 - time2
!   </TEMPLATE>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="time_type" DIM="" DEFAULT="">
!       Returns difference of two time_types.
!   </OUT>

function time_minus(time1, time2)

! Returns difference of two time_types. WARNING: a time type is positive
! so by definition time1 - time2  is the same as time2 - time1.

type(time_type) :: time_minus
type(time_type), intent(in) :: time1, time2

if(.not.module_is_initialized) call time_manager_init

if(time1 > time2) then
   time_minus = decrement_time(time1, time2%seconds, time2%days, time2%ticks)
else
   time_minus = decrement_time(time2, time1%seconds, time1%days, time1%ticks)
endif

end function time_minus
! </FUNCTION>

!--------------------------------------------------------------------------
! <FUNCTION NAME="time_scalar_mult; operator(*)">

!   <OVERVIEW>
!       Returns time multiplied by integer factor n.
!   </OVERVIEW>
!   <DESCRIPTION>
!       Returns time multiplied by integer factor n.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_scalar_mult(time, n)
!   </TEMPLATE>

!   <IN NAME="time" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="n" UNITS="" TYPE="integer" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="time_type" DIM="" DEFAULT="">
!       Returns time multiplied by integer factor n.
!   </OUT>

function time_scalar_mult(time, n)

! Returns time multiplied by integer factor n

type(time_type)             :: time_scalar_mult
type(time_type), intent(in) :: time
integer, intent(in)         :: n
integer                     :: days, seconds, ticks, num_sec
double precision            :: sec_prod, tick_prod

if(.not.module_is_initialized) call time_manager_init

! Multiplying here in a reasonable fashion to avoid overflow is tricky
! Could multiply by some large factor n, and seconds could be up to 86399
! Need to avoid overflowing integers and wrapping around to negatives
! ticks could be up to ticks_per_second-1

tick_prod = dble(time%ticks) * dble(n)
num_sec   = tick_prod/dble(ticks_per_second)
sec_prod  = dble(time%seconds) * dble(n) + num_sec
ticks     = tick_prod - num_sec * ticks_per_second

! If sec_prod is large compared to precision of double precision, things
! can go bad.  Need to warn and abort on this.
! The same is true of tick_prod but is is more likely to happen to sec_prod,
! so let's just test sec_prod. (A test of tick_prod would be necessary only
! if ticks_per_second were greater than seconds_per_day)
if(sec_prod /= 0.0) then
   if(log10(sec_prod) > precision(sec_prod) - 3) call error_mesg('time_scalar_mult', &
      'Insufficient precision to handle scalar product in time_scalar_mult; contact developer',FATAL)
end if

days = sec_prod / dble(seconds_per_day)
seconds = sec_prod - dble(days) * dble(seconds_per_day)

time_scalar_mult = set_time(seconds, time%days * n + days, ticks)

end function time_scalar_mult
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="scalar_time_mult; operator(*)">

!   <OVERVIEW>
!       Returns time multiplied by integer factor n.
!   </OVERVIEW>
!   <DESCRIPTION>
!       Returns time multiplied by integer factor n.
!   </DESCRIPTION>
!   <TEMPLATE>
!     n * time
!     scalar_time_mult(n, time)
!   </TEMPLATE>

!   <IN NAME="time" UNITS="" TYPE="time_type" DIM="">A time interval.</IN>
!   <IN NAME="n" UNITS="" TYPE="integer" DIM=""> An integer. </IN>
!   <OUT NAME="" UNITS="" TYPE="time_type" DIM="" DEFAULT="">
!       Returns time multiplied by integer factor n.
!   </OUT>

function scalar_time_mult(n, time)

! Returns time multipled by integer factor n

type(time_type) :: scalar_time_mult
type(time_type), intent(in) :: time
integer, intent(in) :: n

scalar_time_mult = time_scalar_mult(time, n)

end function scalar_time_mult
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="time_divide; operator(/)">

!   <OVERVIEW>
!       Returns the largest integer, n, for which time1 >= time2 * n.
!   </OVERVIEW>
!   <DESCRIPTION>
!       Returns the largest integer, n, for which time1 >= time2 * n.
!   </DESCRIPTION>
!   <TEMPLATE>
!     n = time1 / time2
!     time_divide(time1, time2)
!   </TEMPLATE>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="integer" DIM="" DEFAULT="">
!       Returns the largest integer, n, for which time1 >= time2 * n.
!   </OUT>

function time_divide(time1, time2)

! Returns the largest integer, n, for which time1 >= time2 * n.

integer                     :: time_divide
type(time_type), intent(in) :: time1, time2
double precision            :: d1, d2

if(.not.module_is_initialized) call time_manager_init

! Convert time intervals to floating point days; risky for general performance?
d1 = time1%days * dble(seconds_per_day) + dble(time1%seconds) + time1%ticks/dble(ticks_per_second)
d2 = time2%days * dble(seconds_per_day) + dble(time2%seconds) + time2%ticks/dble(ticks_per_second)

! Get integer quotient of this, check carefully to avoid round-off problems.
time_divide = d1 / d2

! Verify time_divide*time2 is <= time1 and (time_divide + 1)*time2 is > time1
if(time_divide * time2 > time1 .or. (time_divide + 1) * time2 <= time1) &
   call error_mesg('time_divide',' quotient error :: notify developer',FATAL)

end function time_divide
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="time_real_divide; operator(//)">

!   <OVERVIEW>
!       Returns the double precision quotient of two times.
!   </OVERVIEW>
!   <DESCRIPTION>
!       Returns the double precision quotient of two times.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time1 // time2
!     time_real_divide(time1, time2)
!   </TEMPLATE>

!   <IN NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="integer" DIM="double precision" DEFAULT="">
!       Returns the double precision quotient of two times
!   </OUT>

function time_real_divide(time1, time2)

! Returns the double precision quotient of two times

double precision :: time_real_divide
type(time_type), intent(in) :: time1, time2
double precision :: d1, d2

if(.not.module_is_initialized) call time_manager_init

! Convert time intervals to floating point seconds; risky for general performance?
d1 = time1%days * dble(seconds_per_day) + dble(time1%seconds) + dble(time1%ticks)/dble(ticks_per_second)
d2 = time2%days * dble(seconds_per_day) + dble(time2%seconds) + dble(time2%ticks)/dble(ticks_per_second)

time_real_divide = d1 / d2

end function time_real_divide
! </FUNCTION>

!-------------------------------------------------------------------------
! <SUBROUTINE NAME="time_assignment; assignment(=)">

!   <OVERVIEW>
!       Assigns all components of the time_type variable on
!       RHS to same components of time_type variable on LHS.
!   </OVERVIEW>
!   <DESCRIPTION>
!       Assigns all components of the time_type variable on
!       RHS to same components of time_type variable on LHS.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time1 = time2
!   </TEMPLATE>

!   <OUT NAME="time1" UNITS="" TYPE="time_type" DIM="">
!      A time type variable.
!   </OUT>
!   <IN NAME="time2" UNITS="" TYPE="time_type" DIM="">
!      A time type variable.
!   </IN>

subroutine time_assignment(time1, time2)
type(time_type), intent(out) :: time1
type(time_type), intent(in)  :: time2
   time1%seconds = time2%seconds
   time1%days    = time2%days
   time1%ticks   = time2%ticks
end subroutine time_assignment
! </SUBROUTINE>

!-------------------------------------------------------------------------
! <FUNCTION NAME="time_type_to_real">
!   <OVERVIEW>
!       Converts time to seconds and returns it as a real number
!   </OVERVIEW>
!   <DESCRIPTION>
!       Converts time to seconds and returns it as a real number
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_type_to_real(time)
!   </TEMPLATE>
!   <IN NAME="time" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>

function time_type_to_real(time)

double precision            :: time_type_to_real
type(time_type), intent(in) :: time

if(.not.module_is_initialized) call time_manager_init

time_type_to_real = dble(time%days) * 86400.d0 + dble(time%seconds) + &
     dble(time%ticks)/dble(ticks_per_second)

end function time_type_to_real
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="real_to_time_type">
!   <OVERVIEW>
!       Converts a real number of seconds to a time_type variable
!   </OVERVIEW>
!   <DESCRIPTION>
!       Converts a real number of seconds to a time_type variable
!   </DESCRIPTION>
!   <TEMPLATE>
!     real_to_time_type(x, err_msg)
!   </TEMPLATE>
!   <IN NAME="x" UNITS="" TYPE="real" DIM="">
!      A real number of seconds
!   </IN>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
!   <OUT NAME="real_to_time_type"  TYPE="time_type">
!   </OUT>

 function real_to_time_type(x, err_msg)
 type(time_type)  :: real_to_time_type
 real, intent(in) :: x
 character(len=*), intent(out), optional :: err_msg
 integer          :: seconds, days, ticks
 real             :: real_ticks
 character(len=128) :: err_msg_local

 if(.not.module_is_initialized) call time_manager_init

 days = floor(x/86400.)
 seconds = int(x - 86400.*days)
 real_ticks = x - int(x)
 ticks = nint(real_ticks * ticks_per_second)
 if(.not.set_time_private(seconds, days, ticks, real_to_time_type, err_msg_local)) then
   if(error_handler('function real_to_time_type', err_msg_local, err_msg)) return
 endif

 end function real_to_time_type
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="time_scalar_divide; operator(/)">

!   <OVERVIEW>
!       Returns the largest time, t, for which n * t <= time.
!   </OVERVIEW>
!   <DESCRIPTION>
!       Returns the largest time, t, for which n * t <= time.
!   </DESCRIPTION>
!   <TEMPLATE>
!     time_scalar_divide(time, n)
!   </TEMPLATE>

!   <IN NAME="time" UNITS="" TYPE="time_type" DIM="">
!      A time interval.
!   </IN>
!   <IN NAME="n" UNITS="" TYPE="integer" DIM="">
!      An integer factor.
!   </IN>
!   <OUT NAME="" UNITS="" TYPE="integer" DIM="double precision" DEFAULT="">
!       Returns the largest time, t, for which n * t <= time.
!   </OUT>

function time_scalar_divide(time, n)

! Returns the largest time, t, for which n * t <= time

type(time_type) :: time_scalar_divide
type(time_type), intent(in) :: time
integer, intent(in) :: n
double precision :: d, div, dseconds_per_day, dticks_per_second
integer :: days, seconds, ticks
type(time_type) :: prod1, prod2
character(len=128) tmp1,tmp2
logical :: ltmp

! Convert time interval to floating point days; risky for general performance?
dseconds_per_day  = dble(seconds_per_day)
dticks_per_second = dble(ticks_per_second)
d = time%days*dseconds_per_day*dticks_per_second + dble(time%seconds)*dticks_per_second + dble(time%ticks)
div = d/dble(n)

days = div/(dseconds_per_day*dticks_per_second)
seconds = div/dticks_per_second - days*dseconds_per_day
ticks = div - (days*dseconds_per_day + dble(seconds))*dticks_per_second
time_scalar_divide = set_time(seconds, days, ticks)

! Need to make sure that roundoff isn't killing this
prod1 = n * time_scalar_divide
prod2 = n * (increment_time(time_scalar_divide, days=0, seconds=0, ticks=1))
if(prod1 > time .or. prod2 <= time) then
   call get_time(time, seconds, days, ticks)
   write(tmp1,20) days,seconds,ticks
   call get_time(time_scalar_divide, seconds, days, ticks)
   write(tmp2,30) n,days,seconds,ticks
   ltmp = error_handler('time_scalar_divide',' quotient error:'//trim(tmp1)//trim(tmp2))
 20 format('time=',i7,' days, ',i6,' seconds, ',i6,' ticks')
 30 format('   time divided by',i6,'=',i7,' days, ',i6,' seconds, ',i6,' ticks')
endif

end function time_scalar_divide
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="interval_alarm">

!   <OVERVIEW>
!     Given a time, and a time interval, this function returns true
!     if this is the closest time step to the alarm time.
!   </OVERVIEW>
!   <DESCRIPTION>
!      This is a specialized operation that is frequently performed in models.
!      Given a time, and a time interval, this function is true if this is the
!      closest time step to the alarm time. The actual computation is:
!
!             if((alarm_time - time) &#60;&#61; (time_interval / 2))
!
!      If the function is true, the alarm time is incremented by the
!      alarm_interval; WARNING, this is a featured side effect. Otherwise, the
!      function is false and there are no other effects. CAUTION: if the
!      alarm_interval is smaller than the time_interval, the alarm may fail to
!      return true ever again.  Watch
!      for problems if the new alarm time is less than time + time_interval
!   </DESCRIPTION>
!   <TEMPLATE>
!      interval_alarm(time, time_interval, alarm, alarm_interval)
!   </TEMPLATE>

!   <IN NAME="time" TYPE="time_type"> Current time.  </IN>
!   <IN NAME="time_interval" TYPE="time_type"> A time interval.  </IN>
!   <IN NAME="alarm_interval" TYPE="time_type"> A time interval. </IN>
!   <OUT NAME="interval_alarm" TYPE="logical">
!     Returns either True or false.
!   </OUT>
!   <INOUT NAME="alarm" TYPE="time_type">
!     An alarm time, which is incremented by the alarm_interval
!                   if the function is true.
!   </INOUT>

function interval_alarm(time, time_interval, alarm, alarm_interval)

! Supports a commonly used type of test on times for models.  Given the
! current time, and a time for an alarm, determines if this is the closest
! time to the alarm time given a time step of time_interval.  If this
! is the closest time (alarm - time <= time_interval/2), the function
! returns true and the alarm is incremented by the alarm_interval.  Watch
! for problems if the new alarm time is less than time + time_interval

logical :: interval_alarm
type(time_type), intent(in) :: time, time_interval, alarm_interval
type(time_type), intent(inout) :: alarm

if((alarm - time) <= (time_interval / 2)) then
   interval_alarm = .TRUE.
   alarm = alarm + alarm_interval
else
   interval_alarm = .FALSE.
end if

end function interval_alarm
! </FUNCTION>

!--------------------------------------------------------------------------
! <FUNCTION NAME="repeat_alarm">

!   <OVERVIEW>
!      Repeat_alarm supports an alarm that goes off with
!      alarm_frequency and lasts for alarm_length.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Repeat_alarm supports an alarm that goes off with alarm_frequency and
!      lasts for alarm_length.  If the nearest occurence of an alarm time
!      is less than half an alarm_length from the input time, repeat_alarm
!      is true.  For instance, if the alarm_frequency is 1 day, and the
!      alarm_length is 2 hours, then repeat_alarm is true from time 2300 on
!      day n to time 0100 on day n + 1 for all n.
!   </DESCRIPTION>
!   <TEMPLATE>
!      repeat_alarm(time, alarm_frequency, alarm_length)
!   </TEMPLATE>

!   <IN NAME="time" TYPE="time_type"> Current time.  </IN>
!   <IN NAME="alarm_frequency" TYPE="time_type">
!     A time interval for alarm_frequency.
!   </IN>
!   <IN NAME="alarm_length" TYPE="time_type">
!     A time interval for alarm_length.
!   </IN>
!   <OUT NAME="repeat_alarm" TYPE="logical">
!     Returns either True or false.
!   </OUT>

function repeat_alarm(time, alarm_frequency, alarm_length)

! Repeat_alarm supports an alarm that goes off with alarm_frequency and
! lasts for alarm_length.  If the nearest occurence of an alarm time
! is less than half an alarm_length from the input time, repeat_alarm
! is true.  For instance, if the alarm_frequency is 1 day, and the
! alarm_length is 2 hours, then repeat_alarm is true from time 2300 on
! day n to time 0100 on day n + 1 for all n.

logical :: repeat_alarm
type(time_type), intent(in) :: time, alarm_frequency, alarm_length
type(time_type) :: prev, next

prev = (time / alarm_frequency) * alarm_frequency
next = prev + alarm_frequency
if(time - prev <= alarm_length / 2 .or. next - time <= alarm_length / 2) then
   repeat_alarm = .TRUE.
else
   repeat_alarm = .FALSE.
endif

end function repeat_alarm
! </FUNCTION>

!--------------------------------------------------------------------------

!=========================================================================
! CALENDAR OPERATIONS BEGIN HERE
!=========================================================================

! <SUBROUTINE NAME="set_calendar_type">

!   <OVERVIEW>
!     Sets the default calendar type for mapping time intervals to dates.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A constant number for setting the calendar type.
!   </DESCRIPTION>
!   <TEMPLATE> set_calendar_type(type, err_msg) </TEMPLATE>

!   <IN NAME="type" TYPE="integer" DIM="(scalar)" DEFAULT="NO_CALENDAR">
!     A constant number for setting the calendar type.
!   </IN>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>

subroutine set_calendar_type(type, err_msg)

! Selects calendar for default mapping from time to date.

integer, intent(in) :: type
character(len=*), intent(out), optional :: err_msg
integer :: iday, days_this_month, year, month, day
logical :: leap
character(len=256) :: err_msg_local

if(.not.module_is_initialized) call time_manager_init()

if(present(err_msg)) err_msg = ''

if(type <  0 .or. type > max_type) then
  err_msg_local = 'Illegal calendar type'
  if(error_handler('subroutine set_calendar_type', err_msg_local, err_msg)) return
endif

if(seconds_per_day /= 86400 .and. type /= NO_CALENDAR) then
  err_msg_local = 'Only calendar type NO_CALENDAR is allowed when seconds_per_day is not 86400.'// &
                  ' You are using '//trim(valid_calendar_types(type))//' and seconds_per_day='
  write(err_msg_local(len_trim(err_msg_local)+1:len_trim(err_msg_local)+8),'(i8)') seconds_per_day
  if(error_handler('subroutine set_calendar_type', err_msg_local, err_msg)) return
endif

calendar_type = type

if(type == GREGORIAN) then
  date_to_day = invalid_date
  iday = 0
  do year=1,400
    leap = leap_year_gregorian_int(year)
    do month=1,12
      days_this_month = days_per_month(month)
      if(leap .and. month ==2) days_this_month = 29
      do day=1,days_this_month
        date_to_day(year,month,day) = iday
        iday = iday+1
        coded_date(iday) = day + 32*(month + 16*year)
      enddo ! do day
    enddo ! do month
  enddo ! do year
endif

end subroutine set_calendar_type
! </SUBROUTINE>

!------------------------------------------------------------------------
! <FUNCTION NAME="get_calendar_type">

!   <OVERVIEW>
!      Returns the value of the default calendar type for mapping
!      from time to date.
!   </OVERVIEW>
!   <DESCRIPTION>
!     There are no arguments in this function. It returns the value of
!     the default calendar type for mapping from time to date.
!   </DESCRIPTION>
!   <TEMPLATE>
!     get_calendar_type()
!   </TEMPLATE>

function get_calendar_type()

! Returns default calendar type for mapping from time to date.

integer :: get_calendar_type

get_calendar_type = calendar_type

end function get_calendar_type
! </FUNCTION>

!------------------------------------------------------------------------
! <SUBROUTINE NAME="set_ticks_per_second">

!   <OVERVIEW>
!     Sets the number of ticks per second.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Sets the number of ticks per second.
!   </DESCRIPTION>
!   <TEMPLATE> call set_ticks_per_second(ticks_per_second) </TEMPLATE>
!   <IN NAME="type" TYPE="integer" DIM="(scalar)" DEFAULT="1"> </IN>

subroutine set_ticks_per_second(tps)
integer, intent(in) :: tps

ticks_per_second = tps

end subroutine set_ticks_per_second

! </SUBROUTINE>

!------------------------------------------------------------------------
! <FUNCTION NAME="get_ticks_per_second">

!   <OVERVIEW>
!      Returns the number of ticks per second.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns the number of ticks per second.
!   </DESCRIPTION>
!   <TEMPLATE>
!     ticks_per_second = get_ticks_per_second()
!   </TEMPLATE>

function get_ticks_per_second()
integer :: get_ticks_per_second

get_ticks_per_second = ticks_per_second

end function get_ticks_per_second

! </FUNCTION>
!------------------------------------------------------------------------

!========================================================================
! START OF get_date BLOCK
! <SUBROUTINE NAME="get_date">

!   <OVERVIEW>
!      Given a time_interval, returns the corresponding date under
!      the selected calendar.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Given a time_interval, returns the corresponding date under
!      the selected calendar.
!   </DESCRIPTION>
!   <TEMPLATE>
!     get_date(time, year, month, day, hour, minute, second, tick, err_msg)
!   </TEMPLATE>
!   <IN NAME="time"    TYPE="time_type"> A time interval.</IN>
!   <OUT NAME="year"   TYPE="integer"></OUT>
!   <OUT NAME="month"  TYPE="integer"></OUT>
!   <OUT NAME="day"    TYPE="integer"></OUT>
!   <OUT NAME="hour"   TYPE="integer"></OUT>
!   <OUT NAME="minute" TYPE="integer"></OUT>
!   <OUT NAME="second" TYPE="integer"></OUT>
!   <OUT NAME="tick"   TYPE="integer, optional"></OUT>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
 subroutine get_date(time, year, month, day, hour, minute, second, tick, err_msg)

! Given a time, computes the corresponding date given the selected calendar

 type(time_type), intent(in)    :: time
 integer, intent(out)           :: second, minute, hour, day, month, year
 integer, intent(out), optional :: tick
 character(len=*), intent(out), optional :: err_msg
 character(len=128) :: err_msg_local
 integer :: tick1

 if(.not.module_is_initialized) call time_manager_init
 if(present(err_msg)) err_msg = ''

 select case(calendar_type)
 case(THIRTY_DAY_MONTHS)
   call get_date_thirty   (time, year, month, day, hour, minute, second, tick1)
 case(GREGORIAN)
   call get_date_gregorian(time, year, month, day, hour, minute, second, tick1)
 case(JULIAN)
   call get_date_julian_private   (time, year, month, day, hour, minute, second, tick1)
 case(NOLEAP)
   call get_date_no_leap_private  (time, year, month, day, hour, minute, second, tick1)
 case(NO_CALENDAR)
   err_msg_local = 'Cannot produce a date when the calendar type is NO_CALENDAR'
   if(error_handler('subroutine get_date', err_msg_local, err_msg)) return
 case default
   err_msg_local = 'Invalid calendar type'
   if(error_handler('subroutine get_date', err_msg_local, err_msg)) return
 end select

 if(present(tick)) then
   tick = tick1
 else
   if(tick1 /= 0) then
     err_msg_local = 'tick must be present when time has a second fraction'
     if(error_handler('subroutine get_date', err_msg_local, err_msg)) return
   endif
 endif

 end subroutine get_date
! </SUBROUTINE>
!------------------------------------------------------------------------

 subroutine get_date_gregorian(time, year, month, day, hour, minute, second, tick)

! Computes date corresponding to time for gregorian calendar

 type(time_type), intent(in) :: time
 integer, intent(out) :: year, month, day, hour, minute, second
 integer, intent(out) :: tick
 integer :: iday, isec

 if(Time%seconds >= 86400) then ! This check appears to be unecessary.
   call error_mesg('get_date','Time%seconds .ge. 86400 in subroutine get_date_gregorian',FATAL)
 endif

 iday = mod(Time%days+1,days_in_400_year_period)
 if(iday == 0) iday = days_in_400_year_period

 year = coded_date(iday)/512
 day = mod(coded_date(iday),32)
 month = coded_date(iday)/32 - 16*year

 year = year + 400*((Time%days)/days_in_400_year_period)

 hour = Time%seconds / 3600
 isec  = Time%seconds - 3600*hour
 minute = isec / 60
 second = isec - 60*minute
 tick = time%ticks

 end subroutine get_date_gregorian

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

 subroutine get_date_julian_private(time, year, month, day, hour, minute, second, tick)

! Base date for Julian calendar is year 1 with all multiples of 4
! years being leap years.

 type(time_type), intent(in) :: time
 integer, intent(out) :: second, minute, hour, day, month, year
 integer, intent(out) :: tick
 integer :: m, t, nfour, nex, days_this_month
 logical :: leap

! find number of four year periods; also get modulo number of days
 nfour = time%days / (4 * 365 + 1)
 day = modulo(time%days, (4 * 365 + 1))

! Find out what year in four year chunk
 nex = day / 365
 if(nex == 4) then
    nex = 3
    day = 366
 else
    day=modulo(day, 365) + 1
 endif

! Is this a leap year?
 leap = (nex == 3)

 year = 1 + 4 * nfour + nex

! find month and day
 do m = 1, 12
   month = m
   days_this_month = days_per_month(m)
   if(leap .and. m == 2) days_this_month = 29
   if(day <= days_this_month) exit
   day = day - days_this_month
 end do

! find hour,minute and second
 t = time%seconds
 hour = t / (60 * 60)
 t = t - hour * (60 * 60)
 minute = t / 60
 second = t - 60 * minute
 tick = time%ticks
 end subroutine get_date_julian_private

!------------------------------------------------------------------------
 subroutine get_date_julian(time, year, month, day, hour, minute, second)

! No need to include tick in argument list because this routine
! exists only for interpolator.F90, which does not need it.

 type(time_type), intent(in) :: time
 integer, intent(out) :: second, minute, hour, day, month, year
 integer :: tick

 call get_date_julian_private(time, year, month, day, hour, minute, second, tick)

 end subroutine get_date_julian

!------------------------------------------------------------------------

 subroutine get_date_thirty(time, year, month, day, hour, minute, second, tick)

! Computes date corresponding to time interval for 30 day months, 12
! month years.

 type(time_type), intent(in) :: time
 integer, intent(out) :: second, minute, hour, day, month, year
 integer, intent(out) :: tick
 integer :: t, dmonth, dyear

 t = time%days
 dyear = t / (30 * 12)
 year = dyear + 1
 t = t - dyear * (30 * 12)
 dmonth = t / 30
 month = 1 + dmonth
 day = t -dmonth * 30 + 1

 t = time%seconds
 hour = t / (60 * 60)
 t = t - hour * (60 * 60)
 minute = t / 60
 second = t - 60 * minute
 tick = time%ticks

 end subroutine get_date_thirty


!------------------------------------------------------------------------

 subroutine get_date_no_leap_private(time, year, month, day, hour, minute, second, tick)

! Base date for NOLEAP calendar is year 1.

 type(time_type), intent(in) :: time
 integer, intent(out) :: second, minute, hour, day, month, year
 integer, intent(out) :: tick
 integer :: m, t

! get modulo number of days
 year = time%days / 365 + 1
 day = modulo(time%days, 365) + 1

! find month and day
 do m = 1, 12
   month = m
   if(day <= days_per_month(m)) exit
   day = day - days_per_month(m)
 end do

! find hour,minute and second
 t = time%seconds
 hour = t / (60 * 60)
 t = t - hour * (60 * 60)
 minute = t / 60
 second = t - 60 * minute
 tick = time%ticks

 end subroutine get_date_no_leap_private

!------------------------------------------------------------------------
 subroutine get_date_no_leap(time, year, month, day, hour, minute, second)

! No need to include tick in argument list because this routine
! exists only for interpolator.F90, which does not need it.

 type(time_type), intent(in) :: time
 integer, intent(out) :: second, minute, hour, day, month, year
 integer :: tick

 call get_date_no_leap_private(time, year, month, day, hour, minute, second, tick)

 end subroutine get_date_no_leap
!------------------------------------------------------------------------

! END OF get_date BLOCK
!========================================================================
! START OF set_date BLOCK
! <FUNCTION NAME="set_date">

!   <OVERVIEW>
!      Given an input date in year, month, days, etc., creates a
!      time_type that represents this time interval from the
!      internally defined base date.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Given a date, computes the corresponding time given the selected
!      date time mapping algorithm. Note that it is possible to specify
!      any number of illegal dates; these should be checked for and generate
!      errors as appropriate.
!   </DESCRIPTION>
!   <TEMPLATE>
!     1. set_date(year, month, day, hours, minute, second, tick, err_msg)
!   </TEMPLATE>
!   <TEMPLATE>
!     2. set_date_c(time_string, zero_year_warning, err_msg, allow_rounding)
!      time_string is a character string containing a date formatted
!      according to CF conventions. e.g. '1980-12-31 23:59:59.9'
!   </TEMPLATE>
!   <IN NAME="time"   TYPE="time_type"> A time interval.</IN>
!   <IN NAME="year"   TYPE="integer"></IN>
!   <IN NAME="month"  TYPE="integer"></IN>
!   <IN NAME="day"    TYPE="integer"></IN>
!   <IN NAME="hour"   TYPE="integer"></IN>
!   <IN NAME="minute" TYPE="integer"></IN>
!   <IN NAME="second" TYPE="integer"></IN>
!   <IN NAME="tick"   TYPE="integer"></IN>
!   <IN NAME="zero_year_warning"   TYPE="logical">
!     If the year number is zero, it will be silently changed to one,
!     unless zero_year_warning=.true., in which case a WARNING message
!     will also be issued.
!   </IN>
!   <IN NAME="allow_rounding"   TYPE="logical, optional" DEFAULT=".true.">
!     When .true., any fractions of a second will be rounded off to the nearest tick.
!     When .false., it is a fatal error if the second fraction cannot be exactly
!     represented by a number of ticks.
!   </IN>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
!   <OUT NAME="set_date" TYPE="time_type"> A time interval.</OUT>

 function set_date_private(year, month, day, hour, minute, second, tick, Time_out, err_msg)

! Given a date, computes the corresponding time given the selected
! date time mapping algorithm.  Note that it is possible to specify
! any number of illegal dates; these are checked for and generate
! errors as appropriate.

 logical :: set_date_private
 integer, intent(in) :: year, month, day, hour, minute, second, tick
 type(time_type) :: Time_out
 character(len=*), intent(out) :: err_msg

 if(.not.module_is_initialized) call time_manager_init

 err_msg = ''

 select case(calendar_type)
 case(THIRTY_DAY_MONTHS)
   set_date_private = set_date_thirty   (year, month, day, hour, minute, second, tick, Time_out, err_msg)
 case(GREGORIAN)
   set_date_private = set_date_gregorian(year, month, day, hour, minute, second, tick, Time_out, err_msg)
 case(JULIAN)
   set_date_private = set_date_julian_private   (year, month, day, hour, minute, second, tick, Time_out, err_msg)
 case(NOLEAP)
   set_date_private = set_date_no_leap_private  (year, month, day, hour, minute, second, tick, Time_out, err_msg)
 case (NO_CALENDAR)
   err_msg = 'Cannot produce a date when calendar type is NO_CALENDAR'
   set_date_private = .false.
 case default
   err_msg = 'Invalid calendar type'
   set_date_private = .false.
 end select

 end function set_date_private
! </FUNCTION>

!------------------------------------------------------------------------
 function set_date_i(year, month, day, hour, minute, second, tick, err_msg)
 type(time_type) :: set_date_i
 integer, intent(in) :: day, month, year
 integer, intent(in), optional :: second, minute, hour, tick
 character(len=*), intent(out), optional :: err_msg
 integer :: osecond, ominute, ohour, otick
 character(len=128) :: err_msg_local

 if(.not.module_is_initialized) call time_manager_init
 if(present(err_msg)) err_msg = ''

! Missing optionals are set to 0
 osecond = 0; if(present(second)) osecond = second
 ominute = 0; if(present(minute)) ominute = minute
 ohour   = 0; if(present(hour))   ohour   = hour
 otick   = 0; if(present(tick))   otick   = tick

 if(.not.set_date_private(year, month, day, ohour, ominute, osecond, otick, set_date_i, err_msg_local)) then
   if(error_handler('function set_date_i', err_msg_local, err_msg)) return
 endif

 end function set_date_i
!------------------------------------------------------------------------

 function set_date_c(string, zero_year_warning, err_msg, allow_rounding)

 ! Examples of acceptable forms of string:

 ! 1980-01-01 00:00:00
 ! 1980-01-01 00:00:00.50
 ! 1980-1-1 0:0:0
 ! 1980-1-1

 ! year number must occupy 4 spaces.
 ! months, days, hours, minutes, seconds may occupy 1 or 2 spaces
 ! year, month and day must be separated by a '-'
 ! hour, minute, second must be separated by a ':'
 ! hour, minute, second are optional. If not present then zero is assumed.
 ! second may be a real number.

 ! zero_year_warning:
 ! If the year number is zero, it will be silently changed to one,
 ! unless zero_year_warning=.true., in which case a WARNING message
 ! will also be issued

 type(time_type) :: set_date_c
 character(len=*), intent(in) :: string
 logical,          intent(in),  optional :: zero_year_warning
 character(len=*), intent(out), optional :: err_msg
 logical,          intent(in),  optional :: allow_rounding
 character(len=4) :: formt='(i )'
 logical :: correct_form, zero_year_warning_local, allow_rounding_local
 integer :: i1, i2, i3, i4, i5, i6, i7
 character(len=32) :: string_sifted_left
 integer :: year, month, day, hour, minute, second, tick
 character(len=128) :: err_msg_local

 if(.not.module_is_initialized) call time_manager_init()
 if(present(err_msg)) err_msg = ''
 if(present(zero_year_warning)) then
   zero_year_warning_local = zero_year_warning
 else
   zero_year_warning_local = .true.
 endif
 if(present(allow_rounding)) then
   allow_rounding_local = allow_rounding
 else
   allow_rounding_local = .true.
 endif

 string_sifted_left = adjustl(string)
 i1 = index(string_sifted_left,'-')
 i2 = index(string_sifted_left,'-',back=.true.)
 i3 = index(string_sifted_left,':')
 i4 = index(string_sifted_left,':',back=.true.)
 i5 = len_trim(cut0(string_sifted_left))
 i6 = index(string_sifted_left,'.',back=.true.)
 correct_form = (i1 > 1) ! year number must occupy at least 1 space
 correct_form = correct_form .and. (i2-i1 == 2 .or. i2-i1 == 3) ! month number must occupy 1 or 2 spaces
 if(.not.correct_form) then
   err_msg_local = 'Form of character time stamp is incorrect. The character time stamp is: '//trim(string)
   if(error_handler('function set_date_c', err_msg_local, err_msg)) return
 endif
 write(formt(3:3),'(i1)') i1-1
 read(string_sifted_left(1:i1-1),formt) year
 if(year == 0) then
   year = 1
   if(zero_year_warning_local) then
     call error_mesg('set_date_c','Year zero is invalid. Resetting year to 1', WARNING)
   endif
 endif
 write(formt(3:3),'(i1)') i2-i1-1
 read(string_sifted_left(i1+1:i2-1),formt) month
 i7 = min(i2+2,i5)
 read(string_sifted_left(i2+1:i7),'(i2)') day

 if(i3 == 0) then
! There are no minutes or seconds in the string
   minute = 0
   second = 0
   tick   = 0
   if(i5 <= i2+2) then
 !   There is no clocktime in the string at all
     hour = 0
   else
 !   The clocktime includes only hours
     read(string_sifted_left(i5-1:i5),'(i2)') hour
   endif
 else if(i3 == i4) then
 ! The string includes hours and minutes, but no seconds
   read(string_sifted_left(i3-2:i3-1),'(i2)') hour
   write(formt(3:3),'(i1)') i5-i3
   read(string_sifted_left(i3+1:i5),formt) minute
   second = 0
   tick = 0
 else
 ! The string includes hours, minutes, and seconds
   read(string_sifted_left(i3-2:i3-1),'(i2)') hour
   write(formt(3:3),'(i1)') i4-i3-1
   read(string_sifted_left(i3+1:i4-1),formt) minute
   write(formt(3:3),'(i1)') i5-i4
   if(i6 == 0) then
   ! There are no fractional seconds
     read(string_sifted_left(i4+1:i5),formt) second
     tick = 0
   else
     read(string_sifted_left(i4+1:i6-1),formt) second
     if(.not.get_tick_from_string(string_sifted_left(i6:i5), err_msg_local, allow_rounding_local, tick)) then
       if(error_handler('function set_date_c', err_msg_local, err_msg)) return
     endif
 !   If tick has been rounded up to ticks_per_second, then bump up second.
     if(tick == ticks_per_second) then
       second = second + 1
       tick = 0
     endif
   endif
 endif

 if(.not.set_date_private(year, month, day, hour, minute, second, tick, set_date_c, err_msg_local)) then
   if(error_handler('function set_date_c', err_msg_local, err_msg)) return
 endif

 end function set_date_c
!------------------------------------------------------------------------

 function set_date_gregorian(year, month, day, hour, minute, second, tick, Time_out, err_msg)
 logical :: set_date_gregorian

! Computes time corresponding to date for gregorian calendar.

 integer,          intent(in)  :: year, month, day, hour, minute, second, tick
 type(time_type),  intent(out) :: Time_out
 character(len=*), intent(out) :: err_msg
 integer :: yr1, day1

 if( .not.valid_increments(year,month,day,hour,minute,second,tick,err_msg) ) then
   set_date_gregorian = .false.
   return
 endif

 Time_out%seconds = second + 60*(minute + 60*hour)

 yr1 = mod(year,400)
 if(yr1 == 0) yr1 = 400
 day1 = date_to_day(yr1,month,day)
  if(day1 == invalid_date) then
   err_msg = 'Invalid_date. Date='//convert_integer_date_to_char(year,month,day,hour,minute,second)
   set_date_gregorian = .false.
   return
 endif

 Time_out%days = day1 + days_in_400_year_period*((year-1)/400)
 Time_out%ticks = tick
 err_msg = ''
 set_date_gregorian = .true.

 end function set_date_gregorian

!------------------------------------------------------------------------

 function set_date_julian_private(year, month, day, hour, minute, second, tick, Time_out, err_msg)
 logical :: set_date_julian_private

! Returns time corresponding to date for julian calendar.

 integer,          intent(in)  :: year, month, day, hour, minute, second, tick
 type(time_type),  intent(out) :: Time_out
 character(len=*), intent(out) :: err_msg
 integer :: ndays, m, nleapyr
 logical :: leap

 if( .not.valid_increments(year,month,day,hour,minute,second,tick,err_msg) ) then
   set_date_julian_private = .false.
   return
 endif

 if(month /= 2 .and. day > days_per_month(month)) then
   err_msg = 'Invalid date. Date='//convert_integer_date_to_char(year,month,day,hour,minute,second)
   set_date_julian_private = .false.
   return
 endif

! Is this a leap year?
 leap = (modulo(year,4) == 0)
! compute number of complete leap years from year 1
 nleapyr = (year - 1) / 4

! Finish checking for day specication errors
 if(month == 2 .and. (day > 29 .or. ((.not. leap) .and. day > 28))) then
   err_msg = 'Invalid date. Date='//convert_integer_date_to_char(year,month,day,hour,minute,second)
   set_date_julian_private = .false.
   return
 endif

 ndays = 0
 do m = 1, month - 1
   ndays = ndays + days_per_month(m)
   if(leap .and. m == 2) ndays = ndays + 1
 enddo

 Time_out%seconds = second + 60 * (minute + 60 * hour)
 Time_out%days    = day -1 + ndays + 365*(year - nleapyr - 1) + 366*(nleapyr)
 Time_out%ticks   = tick
 err_msg = ''
 set_date_julian_private = .true.

 end function set_date_julian_private

!------------------------------------------------------------------------
 function set_date_julian(year, month, day, hour, minute, second)

! No need to include tick or err_msg in argument list because this
! routine exists only for interpolator.F90, which does not need them.

 type(time_type) :: set_date_julian
 integer, intent(in) :: year, month, day, hour, minute, second
 character(len=128) :: err_msg

 if(.not.set_date_julian_private(year, month, day, hour, minute, second, 0, set_date_julian, err_msg)) then
   call error_mesg('set_date_julian',trim(err_msg),FATAL)
 endif

 end function set_date_julian
!------------------------------------------------------------------------

 function set_date_thirty(year, month, day, hour, minute, second, tick, Time_out, err_msg)
 logical :: set_date_thirty

! Computes time corresponding to date for thirty day months.

 integer,          intent(in)  :: year, month, day, hour, minute, second, tick
 type(time_type),  intent(out) :: Time_out
 character(len=*), intent(out) :: err_msg

 if( .not.valid_increments(year,month,day,hour,minute,second,tick,err_msg) ) then
   set_date_thirty = .false.
   return
 endif

 if(day > 30) then
   err_msg = 'Invalid date. Date='//convert_integer_date_to_char(year,month,day,hour,minute,second)
   set_date_thirty = .false.
   return
 endif

 Time_out%days    = (day - 1) + 30 * ((month - 1) + 12 * (year - 1))
 Time_out%seconds = second + 60 * (minute + 60 * hour)
 Time_out%ticks   = tick
 err_msg = ''
 set_date_thirty = .true.

 end function set_date_thirty

!------------------------------------------------------------------------

 function set_date_no_leap_private(year, month, day, hour, minute, second, tick, Time_out, err_msg)
 logical :: set_date_no_leap_private

! Computes time corresponding to date for fixed 365 day year calendar.

 integer,          intent(in)  :: year, month, day, hour, minute, second, tick
 type(time_type),  intent(out) :: Time_out
 character(len=*), intent(out) :: err_msg
 integer :: ndays, m

 if( .not.valid_increments(year,month,day,hour,minute,second,tick,err_msg) ) then
   set_date_no_leap_private = .false.
   return
 endif

 if(day > days_per_month(month)) then
   err_msg = 'Invalid date. Date='//convert_integer_date_to_char(year,month,day,hour,minute,second)
   set_date_no_leap_private = .false.
   return
 endif

 ndays = 0
 do m = 1, month - 1
   ndays = ndays + days_per_month(m)
 enddo

! No need for err_msg in call to set_time because previous checks ensure positive value of time.
 Time_out = set_time(second + 60 * (minute + 60 * hour), day -1 + ndays + 365 * (year - 1), tick)
 err_msg = ''
 set_date_no_leap_private = .true.

 end function set_date_no_leap_private
!------------------------------------------------------------------------

 function set_date_no_leap(year, month, day, hour, minute, second)

! No need to include tick or err_msg in argument list because this
! routine exists only for interpolator.F90, which does not need them.

 type(time_type) :: set_date_no_leap
 integer, intent(in) :: year, month, day, hour, minute, second
 character(len=128) :: err_msg

 if(.not.set_date_no_leap_private(year, month, day, hour, minute, second, 0, set_date_no_leap, err_msg)) then
   call error_mesg('set_date_no_leap',trim(err_msg),FATAL)
 endif

 end function set_date_no_leap

!=========================================================================

 function valid_increments(year, month, day, hour, minute, second, tick, err_msg)
 logical :: valid_increments
 integer, intent(in) :: year, month, day, hour, minute, second, tick
 character(len=128), intent(out) :: err_msg

!  Check for invalid values

 err_msg = ''
 valid_increments = .true.
 if(second > 59 .or. second < 0 .or. minute > 59 .or. minute < 0 &
   .or. hour > 23 .or. hour < 0 .or. day > 31 .or. day < 1 &
   .or. month > 12 .or. month < 1 .or. year < 1) then
     err_msg = 'Invalid date. Date='//convert_integer_date_to_char(year,month,day,hour,minute,second)
     valid_increments = .false.
     return
 endif
 if(tick < 0 .or. tick >= ticks_per_second) then
   write(err_msg,'(a,i6)') 'Invalid number of ticks. tick=',tick
   valid_increments = .false.
 endif

 end function valid_increments

!=========================================================================

 function convert_integer_date_to_char(year, month, day, hour, minute, second)
 character(len=19) :: convert_integer_date_to_char
 integer, intent(in) :: year, month, day
 integer, intent(in) :: hour, minute, second

 write(convert_integer_date_to_char,10) year,month,day,hour,minute,second
 10 format(i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)

 end function convert_integer_date_to_char

!=========================================================================
! END OF set_date BLOCK
!=========================================================================

! <FUNCTION NAME="increment_date">

!   <OVERVIEW>
!      Increments the date represented by a time interval and the
!      default calendar type by a number of seconds, etc.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Given a time and some date increment, computes a new time.  Depending
!      on the mapping algorithm from date to time, it may be possible to specify
!      undefined increments (i.e. if one increments by 68 days and 3 months in
!      a Julian calendar, it matters which order these operations are done and
!      we don't want to deal with stuff like that, make it an error).
!   </DESCRIPTION>
!   <TEMPLATE>
!      increment_date(time, years, months, days, hours, minutes, seconds, ticks, err_msg)
!   </TEMPLATE>
!   <IN NAME="time"    TYPE="time_type"> A time interval.</IN>
!   <IN NAME="years"   TYPE="integer">An increment of years.</IN>
!   <IN NAME="months"  TYPE="integer">An increment of months.</IN>
!   <IN NAME="days"    TYPE="integer">An increment of days.</IN>
!   <IN NAME="hours"   TYPE="integer">An increment of hours.</IN>
!   <IN NAME="minutes" TYPE="integer">An increment of minutes.</IN>
!   <IN NAME="seconds" TYPE="integer">An increment of seconds.</IN>
!   <IN NAME="ticks"   TYPE="integer">An increment of ticks.</IN>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
!   <OUT NAME="increment_date" TYPE="time_type"> A new time based on the input
!         time interval and the calendar type.
!   </OUT>
!   <IN NAME="allow_neg_inc" TYPE="logical, optional" DIM="(scalar)" DEFAULT=".true.">
!     When .false., it is a fatal error if any of the input time increments are negative.
!     This mimics the behavior of lima and earlier revisions.
!   </IN>
!   <NOTE>
!     For all but the thirty_day_months calendar, increments to months
!     and years must be made separately from other units because of the
!     non-associative nature of addition.
!     If the result is a negative time (i.e. date before the base date)
!     it is considered a fatal error.
!   </NOTE>

 function increment_date(Time, years, months, days, hours, minutes, seconds, ticks, err_msg, allow_neg_inc)

! Given a time and some date increment, computes a new time.  Depending
! on the mapping algorithm from date to time, it may be possible to specify
! undefined increments (i.e. if one increments by 68 days and 3 months in
! a Julian calendar, it matters which order these operations are done and
! we don't want to deal with stuff like that, make it an error).

! This routine operates in one of two modes.
! 1. days, hours, minutes, seconds, ticks are incremented, years and months must be zero or absent arguments.
! 2. years and/or months are incremented, other time increments must be zero or absent arguments.

 type(time_type) :: increment_date
 type(time_type), intent(in) :: Time
 integer, intent(in), optional :: years, months, days, hours, minutes, seconds, ticks
 character(len=*), intent(out), optional :: err_msg
 logical, intent(in), optional :: allow_neg_inc

 integer :: oyears, omonths, odays, ohours, ominutes, oseconds, oticks
 character(len=128) :: err_msg_local
 logical :: allow_neg_inc_local

 if(.not.module_is_initialized) call time_manager_init
 if(present(err_msg)) err_msg = ''

! Missing optionals are set to 0
 oseconds = 0; if(present(seconds)) oseconds = seconds
 ominutes = 0; if(present(minutes)) ominutes = minutes
 ohours   = 0; if(present(hours))   ohours   = hours
 odays    = 0; if(present(days))    odays    = days
 omonths  = 0; if(present(months))  omonths  = months
 oyears   = 0; if(present(years))   oyears   = years
 oticks   = 0; if(present(ticks))   oticks   = ticks
 allow_neg_inc_local=.true.; if(present(allow_neg_inc)) allow_neg_inc_local=allow_neg_inc

 if(.not.allow_neg_inc_local) then
   if(oyears < 0 .or. omonths < 0 .or. odays < 0 .or. ohours < 0 .or. ominutes < 0 .or. oseconds < 0 .or. oticks < 0) then
     write(err_msg_local,10) oyears, omonths, odays, ohours, ominutes, oseconds, oticks
     if(error_handler('function increment_time', err_msg_local, err_msg)) return
   endif
 endif
 10 format('One or more time increments are negative: '// &
   'years=',i6,' months=',i6,' days=',i6,' hours=',i6,' minutes=',i6,' seconds=',i6,' ticks=',i6)

 if(.not.increment_date_private( &
     Time, oyears, omonths, odays, ohours, ominutes, oseconds, oticks, increment_date, err_msg_local)) then
   if(error_handler('function increment_date', err_msg_local, err_msg)) return
 endif

 end function increment_date

! </FUNCTION>

!=======================================================================

 function increment_date_private(Time, years, months, days, hours, minutes, seconds, ticks, Time_out, err_msg)

! Given a time and some date increment, computes a new time.  Depending
! on the mapping algorithm from date to time, it may be possible to specify
! undefined increments (i.e. if one increments by 68 days and 3 months in
! a Julian calendar, it matters which order these operations are done and
! we don't want to deal with stuff like that, make it an error).

! This routine operates in one of two modes.
! 1. days, hours, minutes, seconds, ticks are incremented, years and months must be zero or absent arguments.
! 2. years and/or months are incremented, other time increments must be zero or absent arguments.

! Negative increments are always allowed in the private version of this routine.

 logical :: increment_date_private
 type(time_type),  intent(in)  :: Time
 integer,          intent(in)  :: years, months, days, hours, minutes, seconds, ticks
 type(time_type),  intent(out) :: Time_out
 character(len=*), intent(out) :: err_msg
 integer :: cyear , cmonth , cday , chour , cminute , csecond , ctick
 logical :: mode_1, mode_2

 err_msg = ''
 increment_date_private = .true.

 mode_1 = days /= 0 .or. hours /= 0 .or. minutes /= 0 .or. seconds /= 0 .or. ticks /= 0
 mode_2 = years /= 0 .or. months /= 0

 if(.not.mode_1 .and. .not.mode_2) then
 ! All time increments are zero
   Time_out = Time
   return
 endif

 if(mode_1 .and. mode_2) then
   err_msg = 'years and/or months must not be incremented with other time units'
   increment_date_private = .false.
   return
 endif

 if(mode_1) then
   csecond = seconds + 60 * (minutes + 60 * hours)
   increment_date_private = increment_time_private(Time, csecond, days, ticks, Time_out, err_msg)
 endif

 if(mode_2) then
 ! Convert Time to a date
   select case(calendar_type)
   case(THIRTY_DAY_MONTHS)
     call get_date_thirty   (Time, cyear, cmonth, cday, chour, cminute, csecond, ctick)
   case(NOLEAP)
     call get_date_no_leap_private  (Time, cyear, cmonth, cday, chour, cminute, csecond, ctick)
   case(JULIAN)
     call get_date_julian_private   (Time, cyear, cmonth, cday, chour, cminute, csecond, ctick)
   case(GREGORIAN)
     call get_date_gregorian(Time, cyear, cmonth, cday, chour, cminute, csecond, ctick)
   case(NO_CALENDAR)
     err_msg = 'Cannot increment a date when the calendar type is NO_CALENDAR'
     increment_date_private = .false.
     return
   case default
     err_msg = 'Invalid calendar type'
     increment_date_private = .false.
     return
   end select

 ! Add month increment
   cmonth = cmonth + months

 ! Adjust year and month number when cmonth falls outside the range 1 to 12
   cyear = cyear + floor((cmonth-1)/12.)
   cmonth = modulo((cmonth-1),12) + 1

 ! Add year increment
   cyear = cyear + years

 ! Convert this back into a time.
   select case(calendar_type)
   case(THIRTY_DAY_MONTHS)
     increment_date_private = set_date_thirty   (cyear, cmonth, cday, chour, cminute, csecond, ctick, Time_out, err_msg)
   case(NOLEAP)
     increment_date_private = set_date_no_leap_private  (cyear, cmonth, cday, chour, cminute, csecond, ctick, Time_out, err_msg)
   case(JULIAN)
     increment_date_private = set_date_julian_private   (cyear, cmonth, cday, chour, cminute, csecond, ctick, Time_out, err_msg)
   case(GREGORIAN)
     increment_date_private = set_date_gregorian(cyear, cmonth, cday, chour, cminute, csecond, ctick, Time_out, err_msg)
   end select
 endif ! if(mode_2)

 end function increment_date_private

!=========================================================================
! <FUNCTION NAME="decrement_date">

!   <OVERVIEW>
!      Decrements the date represented by a time interval and the
!      default calendar type by a number of seconds, etc.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Given a time and some date decrement, computes a new time.  Depending
!      on the mapping algorithm from date to time, it may be possible to specify
!      undefined decrements (i.e. if one decrements by 68 days and 3 months in
!      a Julian calendar, it matters which order these operations are done and
!      we don't want to deal with stuff like that, make it an error).
!   </DESCRIPTION>
!   <TEMPLATE>
!      decrement_date(time, years, months, days, hours, minutes, seconds, ticks, err_msg))
!   </TEMPLATE>
!   <IN NAME="time"    TYPE="time_type"> A time interval.</IN>
!   <IN NAME="years"   TYPE="integer">An decrement of years.</IN>
!   <IN NAME="months"  TYPE="integer">An decrement of months.</IN>
!   <IN NAME="days"    TYPE="integer">An decrement of days.</IN>
!   <IN NAME="hours"   TYPE="integer">An decrement of hours.</IN>
!   <IN NAME="minutes" TYPE="integer">An decrement of minutes.</IN>
!   <IN NAME="seconds" TYPE="integer">An decrement of seconds.</IN>
!   <IN NAME="ticks"   TYPE="integer">An decrement of ticks.</IN>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
!   <OUT NAME="decrement_date" TYPE="time_type"> A new time based on the input
!         time interval and the calendar type.
!   </OUT>
!   <IN NAME="allow_neg_inc" TYPE="logical, optional" DIM="(scalar)" DEFAULT=".true.">
!     When .false., it is a fatal error if any of the input time increments are negative.
!     This mimics the behavior of lima and earlier revisions.
!   </IN>
!   <NOTE>
!     For all but the thirty_day_months calendar, decrements to months
!     and years must be made separately from other units because of the
!     non-associative nature of addition.
!     If the result is a negative time (i.e. date before the base date)
!     it is considered a fatal error.
!   </NOTE>

 function decrement_date(Time, years, months, days, hours, minutes, seconds, ticks, err_msg, allow_neg_inc)

 type(time_type) :: decrement_date
 type(time_type), intent(in) :: Time
 integer, intent(in), optional :: seconds, minutes, hours, days, months, years, ticks
 character(len=*), intent(out), optional :: err_msg
 logical, intent(in), optional :: allow_neg_inc

 integer :: oseconds, ominutes, ohours, odays, omonths, oyears, oticks
 character(len=128) :: err_msg_local
 logical :: allow_neg_inc_local

 if(present(err_msg)) err_msg = ''

 ! Missing optionals are set to 0
 oseconds = 0; if(present(seconds)) oseconds = seconds
 ominutes = 0; if(present(minutes)) ominutes = minutes
 ohours   = 0; if(present(hours))   ohours   = hours
 odays    = 0; if(present(days))    odays    = days
 omonths  = 0; if(present(months))  omonths  = months
 oyears   = 0; if(present(years))   oyears   = years
 oticks   = 0; if(present(ticks))   oticks   = ticks
 allow_neg_inc_local=.true.; if(present(allow_neg_inc)) allow_neg_inc_local=allow_neg_inc

 if(.not.allow_neg_inc_local) then
   if(oyears < 0 .or. omonths < 0 .or. odays < 0 .or. ohours < 0 .or. ominutes < 0 .or. oseconds < 0 .or. oticks < 0) then
     write(err_msg_local,10) oyears, omonths, odays, ohours, ominutes, oseconds, oticks
     if(error_handler('function decrement_date', err_msg_local, err_msg)) return
   endif
 endif
 10 format('One or more time increments are negative: '// &
   'years=',i6,' months=',i6,' days=',i6,' hours=',i6,' minutes=',i6,' seconds=',i6,' ticks=',i6)

 if(.not.increment_date_private( &
     Time, -oyears, -omonths, -odays, -ohours, -ominutes, -oseconds, -oticks, decrement_date, err_msg_local)) then
   if(error_handler('function decrement_date', err_msg_local, err_msg)) return
 endif

 end function decrement_date
 ! </FUNCTION>

!=========================================================================
! START days_in_month BLOCK
! <FUNCTION NAME="days_in_month">

!   <OVERVIEW>
!       Given a time interval, gives the number of days in the
!       month corresponding to the default calendar.
!   </OVERVIEW>
!   <DESCRIPTION>
!       Given a time, computes the corresponding date given the selected
!       date time mapping algorithm.
!   </DESCRIPTION>
!   <TEMPLATE> days_in_month(time) </TEMPLATE>

!   <IN NAME="time" UNITS="" TYPE="time_type" DIM="">A time interval.</IN>
!   <OUT NAME="days_in_month" UNITS="" TYPE="integer" DIM="" DEFAULT="">
!       The number of days in the month given the selected time
!       mapping algorithm.
!   </OUT>

function days_in_month(Time, err_msg)

! Given a time, computes the corresponding date given the selected
! date time mapping algorithm

integer :: days_in_month
type(time_type), intent(in) :: Time
character(len=*), intent(out), optional :: err_msg

if(.not.module_is_initialized) call time_manager_init
if(present(err_msg)) err_msg = ''

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   days_in_month = days_in_month_thirty(Time)
case(GREGORIAN)
   days_in_month = days_in_month_gregorian(Time)
case(JULIAN)
   days_in_month = days_in_month_julian(Time)
case(NOLEAP)
   days_in_month = days_in_month_no_leap(Time)
case(NO_CALENDAR)
   if(error_handler('function days_in_month', &
         'days_in_month makes no sense when the calendar type is NO_CALENDAR', err_msg)) return
case default
   if(error_handler('function days_in_month', 'Invalid calendar type', err_msg)) return
end select
end function days_in_month
! </FUNCTION>

!--------------------------------------------------------------------------

function days_in_month_gregorian(Time)

! Returns the number of days in a gregorian month.

integer :: days_in_month_gregorian
type(time_type), intent(in) :: Time
integer :: year, month, day, hour, minute, second, ticks

call get_date_gregorian(Time, year, month, day, hour, minute, second, ticks)
days_in_month_gregorian = days_per_month(month)
if(leap_year_gregorian_int(year) .and. month == 2) days_in_month_gregorian = 29

end function days_in_month_gregorian

!--------------------------------------------------------------------------
function days_in_month_julian(Time)

! Returns the number of days in a julian month.

integer :: days_in_month_julian
type(time_type), intent(in) :: Time
integer :: year, month, day, hour, minute, second, ticks

call get_date_julian_private(Time, year, month, day, hour, minute, second, ticks)
days_in_month_julian = days_per_month(month)
if(leap_year_julian(Time) .and. month == 2) days_in_month_julian = 29

end function days_in_month_julian

!--------------------------------------------------------------------------
function days_in_month_thirty(Time)

! Returns the number of days in a thirty day month (needed for transparent
! changes to calendar type).

integer :: days_in_month_thirty
type(time_type), intent(in) :: Time

days_in_month_thirty = 30

end function days_in_month_thirty

!--------------------------------------------------------------------------
function days_in_month_no_leap(Time)

! Returns the number of days in a 365 day year month.

integer :: days_in_month_no_leap
type(time_type), intent(in) :: Time
integer :: year, month, day, hour, minute, second, ticks

call get_date_no_leap_private(Time, year, month, day, hour, minute, second, ticks)
days_in_month_no_leap= days_per_month(month)

end function days_in_month_no_leap

! END OF days_in_month BLOCK
!==========================================================================
! START OF leap_year BLOCK
! <FUNCTION NAME="leap_year">

!   <OVERVIEW>
!      Returns true if the year corresponding to the input time is
!      a leap year. Always returns false for THIRTY_DAY_MONTHS and NOLEAP.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns true if the year corresponding to the input time is
!      a leap year. Always returns false for THIRTY_DAY_MONTHS and NOLEAP.
!   </DESCRIPTION>
!   <TEMPLATE> leap_year(time) </TEMPLATE>

!   <IN NAME="time" UNITS="" TYPE="time_type" DIM="">A time interval.</IN>
!   <OUT NAME="leap_year" UNITS="" TYPE="calendar_type" DIM="" DEFAULT="">
!      true if the year corresponding to the input time is a leap year.
!   </OUT>

function leap_year(Time, err_msg)

! Is this date in a leap year for default calendar?

logical :: leap_year
type(time_type), intent(in) :: Time
character(len=*), intent(out), optional :: err_msg

if(.not.module_is_initialized) call time_manager_init
if(present(err_msg)) err_msg=''

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   leap_year = leap_year_thirty(Time)
case(GREGORIAN)
   leap_year = leap_year_gregorian(Time)
case(JULIAN)
   leap_year = leap_year_julian(Time)
case(NOLEAP)
   leap_year = leap_year_no_leap(Time)
case default
   if(error_handler('function leap_year', 'Invalid calendar type in leap_year', err_msg)) return
end select
end function leap_year
! </FUNCTION>

!--------------------------------------------------------------------------

function leap_year_gregorian(Time)

! Is this a leap year for gregorian calendar?

logical :: leap_year_gregorian
type(time_type), intent(in) :: Time
integer :: seconds, minutes, hours, day, month, year

call get_date(Time, year, month, day, hours, minutes, seconds)
leap_year_gregorian = leap_year_gregorian_int(year)

end function leap_year_gregorian

!--------------------------------------------------------------------------

function leap_year_gregorian_int(year)
logical :: leap_year_gregorian_int
integer, intent(in) :: year

leap_year_gregorian_int = mod(year,4) == 0
leap_year_gregorian_int = leap_year_gregorian_int .and. .not.mod(year,100) == 0
leap_year_gregorian_int = leap_year_gregorian_int .or. mod(year,400) == 0

end function leap_year_gregorian_int

!--------------------------------------------------------------------------

function leap_year_julian(Time)

! Returns the number of days in a julian month.

logical :: leap_year_julian
type(time_type), intent(in) :: Time
integer :: seconds, minutes, hours, day, month, year

call get_date(Time, year, month, day, hours, minutes, seconds)
leap_year_julian = ((year / 4 * 4) == year)

end function leap_year_julian

!--------------------------------------------------------------------------

function leap_year_thirty(Time)

! No leap years in thirty day months, included for transparency.

logical :: leap_year_thirty
type(time_type), intent(in) :: Time

leap_year_thirty = .FALSE.

end function leap_year_thirty

!--------------------------------------------------------------------------

function leap_year_no_leap(Time)

! Another tough one; no leap year returns false for leap year inquiry.

logical :: leap_year_no_leap
type(time_type), intent(in) :: Time

leap_year_no_leap = .FALSE.

end function leap_year_no_leap

!> length_of_day returns the number of seconds in a "day".
!!
!! On Earth a "day" is considered to be one solar day:
!!    86400.0 seconds +/- a few milliseconds (ignored)
!!
!! On another planet the length depends on the rotation rate
!! and orbital rate of the planet, given by constant SOLS_PER_DAY.
!!
!! * If a specific Earth calendar is used, it is assumed that the length
!! of a day wanted is 86400.0s, regardless of year length.
!! * If NO_CALENDAR is used, the length of a solar day is calculated from
!! the orbital parameters prescribed in constants.f90.
function length_of_day()
  real :: length_of_day

  select case(calendar_type)
  case(NO_CALENDAR)
    length_of_day = SECONDS_PER_SOL
  case default
    length_of_day = SECONDS_PER_DAY
  end select

end function length_of_day

!END OF leap_year BLOCK
!==========================================================================
! START OF length_of_year BLOCK
! <FUNCTION NAME="length_of_year">

!   <OVERVIEW>
!      Returns the mean length of the year in the default calendar setting.
!   </OVERVIEW>
!   <DESCRIPTION>
!      There are no arguments in this function. It returns the mean
!      length of the year in the default calendar setting.
!   </DESCRIPTION>
!   <TEMPLATE> length_of_year() </TEMPLATE>

function length_of_year()

! What is the length of the year for the default calendar type

type(time_type) :: length_of_year

if(.not.module_is_initialized) call time_manager_init

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   length_of_year = length_of_year_thirty()
case(GREGORIAN)
   length_of_year = length_of_year_gregorian()
case(JULIAN)
   length_of_year = length_of_year_julian()
case(NOLEAP)
   length_of_year = length_of_year_no_leap()
case(NO_CALENDAR)
   length_of_year = length_of_year_exoplanet()

case default
   call error_mesg('length_of_year','Invalid calendar type in length_of_year',FATAL)
end select
end function length_of_year
! </FUNCTION>

!--------------------------------------------------------------------------

function length_of_year_exoplanet()

type(time_type) :: length_of_year_exoplanet
integer :: int_period

int_period = abs(orbital_period)
length_of_year_exoplanet = set_time(int_period, 0)

end function length_of_year_exoplanet

!--------------------------------------------------------------------------

function length_of_year_thirty()

type(time_type) :: length_of_year_thirty

length_of_year_thirty = set_time(0, 360)

end function length_of_year_thirty

!---------------------------------------------------------------------------

function length_of_year_gregorian()

type(time_type) :: length_of_year_gregorian
integer :: days, seconds

days = days_in_400_year_period / 400
seconds = 86400*(days_in_400_year_period/400. - days)
length_of_year_gregorian = set_time(seconds, days)

end function length_of_year_gregorian

!--------------------------------------------------------------------------

function length_of_year_julian()

type(time_type) :: length_of_year_julian

length_of_year_julian = set_time((24 / 4) * 60 * 60, 365)

end function length_of_year_julian

!--------------------------------------------------------------------------

function length_of_year_no_leap()

type(time_type) :: length_of_year_no_leap

length_of_year_no_leap = set_time(0, 365)

end function length_of_year_no_leap

!--------------------------------------------------------------------------

! END OF length_of_year BLOCK
!==========================================================================

! START OF days_in_year BLOCK
! <FUNCTION NAME="days_in_year">

!   <OVERVIEW>
!      Returns the number of days in the calendar year corresponding to
!      the date represented by time for the default calendar.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns the number of days in the calendar year corresponding to
!      the date represented by time for the default calendar.
!   </DESCRIPTION>
!   <TEMPLATE> days_in_year(Time) </TEMPLATE>
!   <IN NAME="Time" TYPE="time_type">A time interval.</IN>
!   <OUT>
!      The number of days in this year for the default calendar type.
!   </OUT>


function days_in_year(Time)

! What is the number of days in this year for the default calendar type

integer :: days_in_year
type(time_type), intent(in) :: Time

if(.not.module_is_initialized) call time_manager_init

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   days_in_year = days_in_year_thirty(Time)
case(GREGORIAN)
   days_in_year = days_in_year_gregorian(Time)
case(JULIAN)
   days_in_year = days_in_year_julian(Time)
case(NOLEAP)
   days_in_year = days_in_year_no_leap(Time)
case default
   call error_mesg('days_in_year','Invalid calendar type in days_in_year',FATAL)
end select
end function days_in_year
! </FUNCTION>

!--------------------------------------------------------------------------

function days_in_year_thirty(Time)

integer :: days_in_year_thirty
type(time_type), intent(in) :: Time

days_in_year_thirty = 360

end function days_in_year_thirty

!---------------------------------------------------------------------------

function days_in_year_gregorian(Time)

integer :: days_in_year_gregorian
type(time_type), intent(in) :: Time

if(leap_year_gregorian(Time)) then
  days_in_year_gregorian = 366
else
  days_in_year_gregorian = 365
endif

end function days_in_year_gregorian

!--------------------------------------------------------------------------
function days_in_year_julian(Time)

integer :: days_in_year_julian
type(time_type), intent(in) :: Time

if(leap_year_julian(Time)) then
   days_in_year_julian = 366
else
   days_in_year_julian = 365
endif

end function days_in_year_julian

!--------------------------------------------------------------------------

function days_in_year_no_leap(Time)

integer :: days_in_year_no_leap
type(time_type), intent(in) :: Time

days_in_year_no_leap = 365

end function days_in_year_no_leap

!--------------------------------------------------------------------------

! END OF days_in_year BLOCK

!==========================================================================
! <FUNCTION NAME="month_name">

!   <OVERVIEW>
!      Returns a character string containing the name of the
!      month corresponding to month number n.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns a character string containing the name of the
!      month corresponding to month number n. Definition is the
!      same for all calendar types.
!   </DESCRIPTION>
!   <TEMPLATE> month_name(n) </TEMPLATE>
!   <IN NAME="n" TYPE="integer">Month number.</IN>
!   <OUT NAME="month_name" TYPE="character(len=9)">
!      The character string associated with a month.
!      All calendars have 12 months and return full
!      month names, not abreviations.
!   </OUT>

function month_name(n)

! Returns character string associated with a month, for now, all calendars
! have 12 months and will return standard names.

character (len=9) :: month_name
integer, intent(in) :: n
character (len = 9), dimension(12) :: months = (/'January  ', 'February ', &
          'March    ', 'April    ', 'May      ', 'June     ', 'July     ', &
          'August   ', 'September', 'October  ', 'November ', 'December '/)

if(.not.module_is_initialized) call time_manager_init

if(n < 1 .or. n > 12) call error_mesg('month_name','Illegal month index',FATAL)

month_name = months(n)

end function month_name
! </FUNCTION>

!==========================================================================

 function error_handler(routine, err_msg_local, err_msg)

! The purpose of this routine is to prevent the addition of an excessive amount of code in order to implement
! the error handling scheme involving an optional error flag of type character.
! It allows one line of code to accomplish what would otherwise require 6 lines.
! A value of .true. for this function is a flag to the caller that it should immediately return to it's caller.

 logical :: error_handler
 character(len=*), intent(in) :: routine, err_msg_local
 character(len=*), intent(out), optional :: err_msg

 error_handler = .false.
 if(present(err_msg)) then
   err_msg = err_msg_local
   error_handler = .true.
 else
   call error_mesg(trim(routine),trim(err_msg_local),FATAL)
 endif

 end function error_handler

!==========================================================================
!------------------------------------------------------------------------
! <SUBROUTINE NAME="time_manager_init">

!   <OVERVIEW>
!      Writes the version information to the log file
!   </OVERVIEW>
!   <DESCRIPTION>
!      Initialization routine.
!      Writes the version information to the log file
!   </DESCRIPTION>
!   <TEMPLATE>time_manager_init()</TEMPLATE>

subroutine time_manager_init ( )

  if (module_is_initialized) return  ! silent return if already called

  call write_version_number (version, tagname)
  module_is_initialized = .true.

end subroutine time_manager_init
! </SUBROUTINE>

!------------------------------------------------------------------------
! <SUBROUTINE NAME="print_time">

!   <OVERVIEW>
!      Prints the given time_type argument as a time (using days, seconds and ticks)
!   </OVERVIEW>
!   <DESCRIPTION>
!      Prints the given time_type argument as a time (using days, seconds and ticks)
!      NOTE: there is no check for PE number.
!   </DESCRIPTION>
!   <TEMPLATE>print_time (time,str,unit)</TEMPLATE>
!   <IN NAME="time" TYPE="time_type"> Time that will be printed. </IN>
!   <IN NAME="str" TYPE="character (len=*)" DEFAULT="TIME: or DATE:">
!      Character string that precedes the printed time or date.
!   </IN>
!   <IN NAME="unit" TYPE="integer">
!      Unit number for printed output. The default unit is stdout.
!   </IN>
subroutine print_time (Time,str,unit)
type(time_type)  , intent(in) :: Time
character (len=*), intent(in), optional :: str
integer          , intent(in), optional :: unit
integer :: s,d,ticks, ns,nd,nt, unit_in
character(len=19) :: fmt

! prints the time to standard output (or optional unit) as days and seconds
! NOTE: there is no check for PE number

  unit_in = stdout()
  if (present(unit)) unit_in = unit

  call get_time (Time,s,d,ticks)

! format output
! get number of digits for days and seconds strings
   nd = int(log10(real(max(1,d))))+1
   ns = int(log10(real(max(1,s))))+1
   nt = int(log10(real(max(1,ticks))))+1
   write (fmt,10) nd, ns, nt
10 format ('(a,i',i2.2,',a,i',i2.2,',a,i',i2.2,')')

  if (present(str)) then
     write (unit_in,fmt) trim(str)//' day=', d, ', sec=', s, ', ticks=', ticks
  else
     write (unit_in,fmt)       'TIME: day=', d, ', sec=', s, ', ticks=', ticks
  endif

end subroutine print_time
! </SUBROUTINE>

!------------------------------------------------------------------------
! <SUBROUTINE NAME="print_date">

!   <OVERVIEW>
!      prints the time to standard output (or optional unit) as a date.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Prints the given time_type argument as a date (using year, month, day,
!      hour, minutes, seconds and ticks). NOTE: there is no check for PE number.
!   </DESCRIPTION>
!   <TEMPLATE> print_date (time,str,unit)
!   </TEMPLATE>
!   <IN NAME="time" TYPE="time_type"> Time that will be printed. </IN>
!   <IN NAME="str" TYPE="character (len=*)" DEFAULT="TIME: or DATE:">
!      Character string that precedes the printed time or date.
!   </IN>
!   <IN NAME="unit" TYPE="integer">
!      Unit number for printed output. The default unit is stdout.
!   </IN>

subroutine print_date (Time,str,unit)
type(time_type)  , intent(in) :: Time
character (len=*), intent(in), optional :: str
integer          , intent(in), optional :: unit
integer :: y,mo,d,h,m,s, unit_in
character(len=9) :: mon

! prints the time to standard output (or optional unit) as a date
! NOTE: there is no check for PE number

  unit_in = stdout()
  if (present(unit)) unit_in = unit

  call get_date (Time,y,mo,d,h,m,s)
  mon = month_name(mo)
  if (present(str)) then
     write (unit_in,10) trim(str)//' ', y,mon(1:3),' ',d,' ',h,':',m,':',s
  else
     write (unit_in,10)       'DATE: ', y,mon(1:3),' ',d,' ',h,':',m,':',s
  endif
10 format (a,i4,1x,a3,4(a1,i2.2))

end subroutine print_date
! </SUBROUTINE>

!------------------------------------------------------------------------
! <FUNCTION NAME="valid_calendar_types">

!   <OVERVIEW>
!     Returns a character string that describes the
!     calendar type corresponding to the input integer.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns a character string that describes the
!     calendar type corresponding to the input integer.
!   </DESCRIPTION>
!   <IN NAME="ncal" TYPE="integer">
!     An integer corresponding to a valid calendar type.
!   </IN>
!   <OUT NAME="err_msg" TYPE="character, optional" DIM="(scalar)">
!     When present, and when non-blank, a fatal error condition as been detected.
!     The string itself is an error message.
!     It is recommended that, when err_msg is present in the call
!     to this routine, the next line of code should be something
!     similar to this:
!     if(err_msg /= '') call error_mesg('my_routine','additional info: '//trim(err_msg),FATAL)
!   </OUT>
!   <OUT NAME="valid_calendar_types" TYPE="character(len=24)">
!     A character string describing the calendar type.
!   </OUT>

function valid_calendar_types(ncal, err_msg)
integer, intent(in) :: ncal
character(len=*), intent(out), optional :: err_msg
character(len=24) :: valid_calendar_types
character(len=128) :: err_msg_local

if(.not.module_is_initialized) call time_manager_init

if(present(err_msg)) err_msg = ''

if(ncal == NO_CALENDAR) then
  valid_calendar_types = 'NO_CALENDAR             '
else if(ncal == THIRTY_DAY_MONTHS) then
  valid_calendar_types = 'THIRTY_DAY_MONTHS       '
else if(ncal == JULIAN) then
  valid_calendar_types = 'JULIAN                  '
else if(ncal == GREGORIAN) then
  valid_calendar_types = 'GREGORIAN               '
else if(ncal == NOLEAP) then
  valid_calendar_types = 'NOLEAP                  '
else
  write(err_msg_local,'(a,i4,a)') 'calendar type=',ncal,' is invalid.'
  if(error_handler('function valid_calendar_types', err_msg_local, err_msg)) return
endif
end function valid_calendar_types
! </FUNCTION>
!------------------------------------------------------------------------

!--- get the a character string that represents the time. The format will be
!--- yyyymmdd.hhmmss
function date_to_string(time, err_msg)
  type(time_type),  intent(in)            :: time
  character(len=*), intent(out), optional :: err_msg
  character(len=128)                      :: err_msg_local
  character(len=15)                       :: date_to_string
  integer                                 :: yr,mon,day,hr,min,sec

  if(present(err_msg)) err_msg = ''
  call get_date(time,yr,mon,day,hr,min,sec)
  if (yr <= 9999) then
     write(date_to_string,'(I4.4,I2.2,I2.2,".",I2.2,I2.2,I2.2)') yr, mon, day, hr, min, sec
  else
     write(err_msg_local, '(a,i4.4,a)') 'year = ', yr, ' should be less than 10000'
     if(error_handler('function date_to_string', err_msg_local, err_msg)) return
  endif

end function date_to_string

end module time_manager_mod

! <INFO>

!   <TESTPROGRAM NAME="time_main2">
!    <PRE>
!        use time_manager_mod
!        implicit none
!        type(time_type) :: dt, init_date, astro_base_date, time, final_date
!        type(time_type) :: next_rad_time, mid_date
!        type(time_type) :: repeat_alarm_freq, repeat_alarm_length
!        integer :: num_steps, i, days, months, years, seconds, minutes, hours
!        integer :: months2, length
!        real :: astro_days
!
!   !Set calendar type
!   !    call set_calendar_type(THIRTY_DAY_MONTHS)
!        call set_calendar_type(JULIAN)
!   !    call set_calendar_type(NOLEAP)
!
!   ! Set timestep
!        dt = set_time(1100, 0)
!
!   ! Set initial date
!        init_date = set_date(1992, 1, 1)
!
!   ! Set date for astronomy delta calculation
!        astro_base_date = set_date(1970, 1, 1, 12, 0, 0)
!
!   ! Copy initial time to model current time
!        time = init_date
!
!   ! Determine how many steps to do to run one year
!        final_date = increment_date(init_date, years = 1)
!        num_steps = (final_date - init_date) / dt
!        write(*, *) 'Number of steps is' , num_steps
!
!   ! Want to compute radiation at initial step, then every two hours
!        next_rad_time = time + set_time(7200, 0)
!
!   ! Test repeat alarm
!        repeat_alarm_freq = set_time(0, 1)
!        repeat_alarm_length = set_time(7200, 0)
!
!   ! Loop through a year
!        do i = 1, num_steps
!
!   ! Increment time
!        time = time + dt
!
!   ! Test repeat alarm
!        if(repeat_alarm(time, repeat_alarm_freq, repeat_alarm_length)) &
!        write(*, *) 'REPEAT ALARM IS TRUE'
!
!   ! Should radiation be computed? Three possible tests.
!   ! First test assumes exact interval; just ask if times are equal
!   !     if(time == next_rad_time) then
!   ! Second test computes rad on last time step that is <= radiation time
!   !     if((next_rad_time - time) < dt .and. time < next_rad) then
!   ! Third test computes rad on time step closest to radiation time
!         if(interval_alarm(time, dt, next_rad_time, set_time(7200, 0))) then
!           call get_date(time, years, months, days, hours, minutes, seconds)
!           write(*, *) days, month_name(months), years, hours, minutes, seconds
!
!   ! Need to compute real number of days between current time and astro_base
!           call get_time(time - astro_base_date, seconds, days)
!           astro_days = days + seconds / 86400.
!   !       write(*, *) 'astro offset ', astro_days
!        end if
!
!   ! Can compute daily, monthly, yearly, hourly, etc. diagnostics as for rad
!
!   ! Example: do diagnostics on last time step of this month
!        call get_date(time + dt, years, months2, days, hours, minutes, seconds)
!        call get_date(time, years, months, days, hours, minutes, seconds)
!        if(months /= months2) then
!           write(*, *) 'last timestep of month'
!           write(*, *) days, months, years, hours, minutes, seconds
!        endif
!
!   ! Example: mid-month diagnostics; inefficient to make things clear
!        length = days_in_month(time)
!        call get_date(time, years, months, days, hours, minutes, seconds)
!        mid_date = set_date(years, months, 1) + set_time(0, length) / 2
!
!        if(time < mid_date .and. (mid_date - time) < dt) then
!           write(*, *) 'mid-month time'
!           write(*, *) days, months, years, hours, minutes, seconds
!        endif
!
!        end do
!
!    </PRE>
!   end program time_main2

!   </TESTPROGRAM>
!   <NOTE>
!     The <a name="base date">base date</a> is implicitly defined so users don't
!     need to be concerned with it. For the curious, the base date is defined as
!     0 seconds, 0 minutes, 0 hours, day 1, month 1, year 1
!   </NOTE>
!   <NOTE>
!     Please note that a time is a positive definite quantity.
!   </NOTE>
!   <NOTE>
!     See the <LINK SRC="TEST PROGRAM">Test Program </LINK> for a simple program
!     that shows some of the capabilities of the time manager.
!   </NOTE>
! </INFO>

#ifdef test_time_manager
 program test
 use          mpp_mod, only: input_nml_file
 use          fms_mod, only: fms_init, fms_end, stderr
 use          fms_mod, only: open_namelist_file, check_nml_error, close_file, open_file
 use    constants_mod, only: constants_init, rseconds_per_day=>seconds_per_day
 use       fms_io_mod, only: fms_io_exit
 use time_manager_mod, only: time_type, set_date, get_date, set_time, set_calendar_type, real_to_time_type
 use time_manager_mod, only: length_of_year, leap_year, days_in_month, days_in_year, print_time
 use time_manager_mod, only: set_ticks_per_second, get_ticks_per_second
 use time_manager_mod, only: decrement_date, increment_date, get_time, increment_time, decrement_time
 use time_manager_mod, only: JULIAN, GREGORIAN, THIRTY_DAY_MONTHS, NOLEAP
 use time_manager_mod, only: operator(-), operator(+),  operator(*),  operator(/),  &
                             operator(>), operator(>=), operator(==), operator(/=), &
                             operator(<), operator(<=), operator(//), assignment(=)

 implicit none

 type(time_type) :: Time, time1, time2
 real    :: xx
 integer :: yr, mo, day, hr, min, sec, ticks
 integer :: year, month, dday, days_this_month
 integer :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
 logical :: leap
 integer :: nr, icode, nmlunit, ierr, io, nn, errunit, outunit
 character(len=256) :: err_msg, char_date
 character(len=8),  allocatable, dimension(:) :: test_time
 character(len=23), allocatable, dimension(:) :: test_date
 character(len=8) :: test_name

logical :: test1 =.true.,test2 =.true.,test3 =.true.,test4 =.true.,test5 =.true.,test6 =.true.,test7 =.true.,test8 =.true.
logical :: test9 =.true.,test10=.true.,test11=.true.,test12=.true.,test13=.true.,test14=.true.,test15=.true.,test16=.true.
logical :: test17=.true.,test18=.true.,test19=.true.

 namelist / test_nml / test1 ,test2 ,test3 ,test4 ,test5 ,test6 ,test7 ,test8,  &
                       test9 ,test10,test11,test12,test13,test14,test15,test16, &
                       test17,test18,test19

 call fms_init
 call constants_init

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, test_nml, iostat=io)
      ierr = check_nml_error (io, 'test_nml')
#else
 nmlunit = open_namelist_file()
 ierr=1
 do while (ierr /= 0)
   read(nmlunit, nml=test_nml, iostat=io, end=12)
   ierr = check_nml_error (io, 'test_nml')
 enddo
 12 call close_file (nmlunit)
#endif

 outunit = open_file(file='test_time_manager.out', form='formatted', action='write')
 errunit = stderr()
 call set_ticks_per_second(10)

 !==============================================================================================
 ! Tests of set_time_i and get_time without ticks

 if(test1) then
   write(outunit,'(/,a)') '#################################  test1  #################################'
   Time = set_time(seconds=2, days=1)
   call get_time(Time, sec, day, ticks)
   write(outunit,'(a,i2,a,i8,a,i2)') ' test1.1: days=',day,' seconds=',sec,' ticks=',ticks
   call get_time(Time, sec, day)
   write(outunit,'(a,i2,a,i8)') ' test1.2: days=',day,' seconds=',sec
   call get_time(Time, sec)
   write(outunit,'(a,i8)') ' test1.2: seconds=',sec
 endif
 !==============================================================================================
 ! Tests of set_time_i and get_time with ticks

 if(test2) then
   write(outunit,'(/,a)') '#################################  test2  #################################'
   Time = set_time(seconds=2, days=1, ticks=5)
   call get_time(Time, sec, day, ticks)
   write(outunit,'(a,i2,a,i6,a,i2)') ' test2.1: days=',day,' seconds=',sec,' ticks=',ticks
   call get_time(Time, sec, ticks=ticks)
   write(outunit,'(a,i6,a,i2)') ' test2.2: seconds=',sec,' ticks=',ticks
   call get_time(Time, sec, day, err_msg=err_msg)
   if(err_msg /= '') then
     write(outunit,'(a)') ' test2.3 successful: '//trim(err_msg)
   else
     write(outunit,'(a,i2,a,i8)') ' test2.3 fails. days=',day,' seconds=',sec
   endif
   call get_time(Time, sec, err_msg=err_msg)
   if(err_msg /= '') then
     write(outunit,'(a)') ' test2.4 successful: '//trim(err_msg)
   else
     write(outunit,'(a,i8)') ' test2.4 fails.  seconds=',sec
   endif
 endif
 !==============================================================================================
 ! Tests of time operators
 ! Test of function scalar_time_mult is not necessary, it simply calls time_scalar_mult.
 ! Test of function time_ne is not necessary, it simply calls time_eq.
 ! Test of function time_ge is not necessary, it simply calls time_gt.
 ! Test of function time_le is not necessary, it simply calls time_lt and time_eq.
 ! Test of function time_ne is not necessary, it simply calls time_eq.

  if(test3) then
    write(outunit,'(/,a)') '#################################  test3  #################################'
 !  Test of function time_plus
    call print_time(set_time(seconds=0, days=2, ticks=5) + set_time(seconds=0, days=2, ticks=6), 'test3.1:', unit=outunit)

 !  Test of function time_minus
 !  The minus operator for time ensures a positive result. In effect is does this: abs(time1-time2)
    call print_time(set_time(seconds=0, days=2, ticks=5) - set_time(seconds=0, days=2, ticks=6), 'test3.2:', unit=outunit)

 !  Test of function time_scalar_mult.  Note that 25000*86399 is greater than huge = 2**31 - 1
    call print_time(2*set_time(seconds=0, days=2, ticks=6), 'test3.3:', unit=outunit)
    call print_time(25000*set_time(seconds=86399, days=0, ticks=0), 'test3.4:', unit=outunit)

 !  Test of function time_scalar_divide
    call print_time(set_time(seconds=0, days=60000, ticks=2)/2, 'test3.5:', unit=outunit)

 !  Test of function time_real_divide
    xx = set_time(seconds=0, days=60000, ticks=2)//set_time(seconds=86400)
    write(outunit,'("test3.6: xx=",f15.9)') xx

 !  Test of function time_divide
    nn = set_time(seconds=0, days=60000, ticks=2)//set_time(seconds=86400)
    write(outunit,'("test3.7: nn=",i6)') nn

 !  Test of function time_gt
    if(set_time(seconds=1, days=1, ticks=2) > set_time(seconds=1, days=1, ticks=1)) then
      write(outunit,'("test3.8 successful")')
    else
      write(outunit,'("test3.8 fails")')
    endif
    if(set_time(seconds=1, days=1, ticks=2) > set_time(seconds=1, days=1, ticks=2)) then
      write(outunit,'("test3.9 fails")')
    else
      write(outunit,'("test3.9 successful")')
    endif

 !  Test of function time_lt
    if(set_time(seconds=1, days=1, ticks=1) < set_time(seconds=1, days=1, ticks=2)) then
      write(outunit,'("test3.10 successful")')
    else
      write(outunit,'("test3.10 fails")')
    endif
    if(set_time(seconds=1, days=1, ticks=2) < set_time(seconds=1, days=1, ticks=2)) then
      write(outunit,'("test3.11 fails")')
    else
      write(outunit,'("test3.11 successful")')
    endif

 !  Test of function time_eq
    if(set_time(seconds=1, days=1, ticks=1) == set_time(seconds=1, days=1, ticks=1)) then
      write(outunit,'("test3.12 successful")')
    else
      write(outunit,'("test3.12 fails")')
    endif
    if(set_time(seconds=1, days=1, ticks=1) == set_time(seconds=1, days=1, ticks=2)) then
      write(outunit,'("test3.13 fails")')
    else
      write(outunit,'("test3.13 successful")')
    endif
  endif
 !==============================================================================================
 ! Tests of set_time_c

 if(test4) then
   write(outunit,'(/,a)') '#################################  test4  #################################'
   test_name = 'test4.  '
   allocate(test_time(15))
   test_time( 1: 6) = (/'1 10    ','1 10.   ','1 10.000','1  0.0  ','1   .000','1   .   '/)
   test_time( 7: 9) = (/'1 10.20 ','1 10.300','1  0.40 '/)
   test_time(10:15) = (/'1   .510','2 .50001','1.0 10.2','10.30000','10-0.40 ','10:1.510'/) ! invalid forms
   do nr=1,9
     write(test_name(7:8),'(i2.2)') nr
     Time = set_time(trim(test_time(nr)), err_msg=err_msg, allow_rounding=.false.)
     if(err_msg == '') then
       call print_time(Time, test_name//':', unit=outunit)
     else
       write(outunit,'(a)') test_name//' fails: '//trim(err_msg)
     endif
   enddo

   test_time(1:6) = (/'1   .510','2 .50001','1.0 10.2','10.30000','10-0.40 ','10:1.510'/)
   do nr=10,15
     write(test_name(7:8),'(i2.2)') nr
     Time = set_time(trim(test_time(nr)), err_msg=err_msg, allow_rounding=.false.)
     if(err_msg /= '') then
       write(outunit,'(a)') test_name//' successful: '//trim(err_msg)
     else
       write(outunit,'(a)') test_name//' fails '
     endif
   enddo
 endif

 !==============================================================================================
 ! Tests of set_date_i

 if(test5) then
   write(outunit,'(/,a)') '#################################  test5  #################################'
   call set_calendar_type(JULIAN)
   call print_time(set_date(1980, 1, 1, 0, 0, 0),' test5.1:', unit=outunit)
   call print_time(set_date(1980, 1, 2, 3, 4, 5, 6),' test5.2:', unit=outunit)
   call print_time(set_date(1980, 1, 2, tick=6),' test5.3:', unit=outunit)
   Time = set_date(1980, 1, 2, tick=10, err_msg=err_msg)
   if(err_msg == '') then
     write(outunit,'(a)') ' test5.4 fails'
   else
     write(outunit,'(a)') ' test5.4 successful: '//trim(err_msg)
   endif
 endif
 !==============================================================================================
 ! Tests of set_date_c

 if(test6) then
   write(outunit,'(/,a)') '#################################  test6  #################################'
   test_name = 'test6.  '
   call set_calendar_type(GREGORIAN)
   allocate(test_date(6))
   test_date(1:3) = (/' 1980-12-30 01:01:11   ',' 1980-12-30 01:01:11.50',' 1980-12-30 01:01:11.55'/)
   test_date(4:6) = (/' 1980-12-30 01:01:11.96','   1980-1-3 1:1:11     ','   1980-1-3 1:1:11.99  '/)
   do nr=1,6
     write(test_name(7:8),'(i2.2)') nr
     Time = set_date(trim(test_date(nr)), err_msg=err_msg, allow_rounding=.true., zero_year_warning=.true.)
     if(err_msg == '') then
       call print_time(Time,test_name//' successful:', unit=outunit)
     else
       write(outunit,'(a)') test_name//'fails: '//trim(err_msg)
     endif
   enddo
   call set_calendar_type(THIRTY_DAY_MONTHS)
   call print_time(set_date('1900-02-30 00:00:00'),'test6.7:', unit=outunit)
   Time = set_date('1900-01-31 00:00:00', err_msg=err_msg)
   if(err_msg == '') then
     write(outunit,'(a)') 'test6.8 fails'
   else
     write(outunit,'(a)') 'test6.8 successful '//trim(err_msg)
   endif
   call set_calendar_type(JULIAN)
   Time = set_date('1901-02-29 00:00:00', err_msg=err_msg)
   if(err_msg == '') then
     write(outunit,'(a)') 'test6.9 fails'
   else
     write(outunit,'(a)') 'test6.9 successful '//trim(err_msg)
   endif
 endif
!==============================================================================================
! Tests of decrement_date and increment_date

 if(test7) then
   write(outunit,'(/,a)') '#################################  test7  #################################'
   char_date = '1904-01-01 00:00:00'
   write(outunit,'(a)') ' Initial date='//trim(char_date)//':00'

   do nr=1,4
     write(outunit,'("=================================================================")')
     if(nr == 1) then
       call set_calendar_type(THIRTY_DAY_MONTHS)
       write(outunit,'(" THIRTY_DAY_MONTHS")')
     endif
     if(nr == 2) then
       call set_calendar_type(NOLEAP)
       write(outunit,'(" NOLEAP")')
     endif
     if(nr == 3) then
       call set_calendar_type(JULIAN)
       write(outunit,'(" JULIAN")')
     endif
     if(nr == 4) then
       call set_calendar_type(GREGORIAN)
       write(outunit,'(" GREGORIAN")')
     endif
     time1 = set_date(trim(char_date))
     do year=-1,1
       do month=-1,1
         write(outunit,'(" test of decrement_date increments: year=",i2," month=",i2)') year,month
         time2 = decrement_date(time1, year, month, err_msg=err_msg)
         if(err_msg /= '') then
           write(outunit,'(a)') 'test of decrement_date fails '//trim(err_msg)
         else
           call get_date(time2, yr, mo, day, hr, min, sec, ticks)
           write(outunit,20) yr, mo, day, hr, min, sec, ticks
         endif
       enddo
     enddo
     time1 = set_date(1, 1, 2, 1, 1, 1, 1, err_msg)
     write(outunit,'(" Initial date = 01-01-02 01:01:01:01")')
     do icode=0,242
       day   = modulo(icode/81,3) - 1
       hr    = modulo(icode/27,3) - 1
       min   = modulo(icode/9, 3) - 1
       sec   = modulo(icode/3, 3) - 1
       ticks = modulo(icode   ,3) - 1
       write(outunit,11) day, hr, min, sec, ticks
       time2 = increment_date(time1, 0, 0, day, hr, min, sec, ticks, err_msg)
       call get_date(time2, yr, mo, day, hr, min, sec, ticks)
       write(outunit,20) yr, mo, day, hr, min, sec, ticks
     enddo
   enddo
 endif

  11 format(' test of increment_date increments: day=',i2,' hr=',i2,' min=',i2,' sec=',i2,' ticks=',i2)
  20 format(' time=',i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2, ':', i2.2)
 !==============================================================================================
 ! Tests involving Feb 29

  if(test8) then
    write(outunit,'(/,a)') '#################################  test8  #################################'
    call set_calendar_type(THIRTY_DAY_MONTHS)
    Time = set_date('1904-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      call print_time(Time, 'test8.1 successful', unit=outunit)
    else
      write(outunit,'(a)') 'test8.1 fails: '//trim(err_msg)
    endif

    call set_calendar_type(NOLEAP)
    Time = set_date('1904-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test8.2 fails'
    else
      write(outunit,'(a)') 'test8.2 successful: '//trim(err_msg)
    endif

    call set_calendar_type(GREGORIAN)
    Time = set_date('1900-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test8.3 fails'
    else
      write(outunit,'(a)') 'test8.3 successful: '//trim(err_msg)
    endif
    Time = set_date('2000-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test8.4 successful'
    else
      write(outunit,'(a)') 'test8.4 fails: '//trim(err_msg)
    endif

    call set_calendar_type(JULIAN)
    Time = set_date('1900-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test8.5 successful'
    else
      write(outunit,'(a)') 'test8.5 fails: '//trim(err_msg)
    endif
    Time = set_date('1901-02-29 00:00:00', err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test8.6 fails'
    else
      write(outunit,'(a)') 'test8.6 successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
 ! Tests of days_in_month

  if(test9) then
    write(outunit,'(/,a)') '#################################  test9  #################################'
    day = days_in_month(set_date('1901-02-28 00:00:00'))
    write(outunit,'(a,i4)') ' test9.1: day=',day
    day = days_in_month(set_date('1901-07-01 00:00:00'))
    write(outunit,'(a,i4)') ' test9.2: day=',day
  endif
 !==============================================================================================
 ! Tests of get_time error flag

  if(test10) then
    write(outunit,'(/,a)') '#################################  test10  #################################'
    Time = set_time(seconds=2, days=1, ticks=1)
    call get_time(Time, seconds=sec, days=day, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test10.1 fails'
    else
      write(outunit,'(a)') 'test10.1 successful: '//trim(err_msg)
    endif
    call set_calendar_type(GREGORIAN)
    Time = set_time(seconds=2, days=1, ticks=1)
    call get_date(Time, yr, mo, day, hr, min, sec, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test10.2 fails'
    else
      write(outunit,'(a)') 'test10.2 successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
 ! Tests of increment_time and decrement_time

  if(test11) then
    write(outunit,'(/,a)') '#################################  test11  #################################'
    call print_time(increment_time(set_time(seconds=0, days=2), seconds=0, days=1),'test11.1:', unit=outunit)
    call print_time(decrement_time(set_time(seconds=0, days=2), seconds=0, days=1),'test11.2:', unit=outunit)
    call print_time(increment_time(set_time(seconds=0, days=2, ticks=5), seconds=400, days=1, ticks=14),'test11.3:', unit=outunit)
    call print_time(decrement_time(set_time(seconds=0, days=2, ticks=5), seconds=400, days=1, ticks=14),'test11.4:', unit=outunit)
  endif
 !==============================================================================================
 !  Tests of negative increments in increment_time and decrement_time

  if(test12) then
    write(outunit,'(/,a)') '#################################  test12  #################################'
    call print_time(increment_time(set_time(seconds=0, days=2), seconds=0, days=-1),'test12.1:', unit=outunit)
    call print_time(decrement_time(set_time(seconds=0, days=2), seconds=0, days=-1),'test12.2:', unit=outunit)
    call print_time(increment_time(set_time(seconds=0, days=2, ticks=5),seconds=-400,days=-1,ticks=-14),'test12.3:',unit=outunit)
    call print_time(decrement_time(set_time(seconds=0, days=2, ticks=5),seconds=-400,days=-1,ticks=-14),'test12.4:',unit=outunit)
  endif
 !==============================================================================================
 !  Test of trap for negative time

  if(test13) then
    write(outunit,'(/,a)') '#################################  test13  #################################'
    Time = set_time(seconds= 2, days=0, ticks=-21, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test13.1 fails'
    else
      write(outunit,'(a)') 'test13.1 successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
 !  Tests of negative seconds and/or ticks

  if(test14) then
    write(outunit,'(/,a)') '#################################  test14  #################################'
    call print_time(set_time(seconds=-86399, days=2, ticks=-10),'test14.1:', unit=outunit)
    call print_time(set_time(seconds=-86390, days=2, ticks=-95),'test14.2:', unit=outunit)
    call print_time(set_time(seconds= 86400, days=2, ticks= 95),'test14.3:', unit=outunit)
  endif
 !==============================================================================================
 !  Tests of consistency of day numbering between calendars

  if(test15) then
    write(outunit,'(/,a)') '#################################  test15  #################################'
    call set_calendar_type(GREGORIAN)
    Time = set_date(1, 1, 1)
    call get_time(Time, sec, day)
    write(outunit,10) 'GREGORIAN',day

    call set_calendar_type(JULIAN)
    Time = set_date(1, 1, 1)
    call get_time(Time, sec, day)
    write(outunit,10) 'JULIAN',day

    call set_calendar_type(THIRTY_DAY_MONTHS)
    Time = set_date(1, 1, 1)
    call get_time(Time, sec, day)
    write(outunit,10) 'THIRTY_DAY_MONTHS',day

    call set_calendar_type(NOLEAP)
    Time = set_date(1, 1, 1)
    call get_time(Time, sec, day)
    write(outunit,10) 'NOLEAP',day
  endif

  10 format(a17,' Jan 1 year 1 is day=',i6)

 !==============================================================================================
 ! Tests of error message for invalid dates

  if(test16) then
    write(outunit,'(/,a)') '#################################  test16  #################################'
    call set_calendar_type(GREGORIAN)
    Time = set_date(1900, 1, 32, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.1 fails'
    else
      write(outunit,'(a)') 'test16.1 successful: '//trim(err_msg)
    endif

    Time = set_date(1900, 4, 31, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.2 fails'
    else
      write(outunit,'(a)') 'test16.2 successful: '//trim(err_msg)
    endif

    Time = set_date(1900, 2, 29, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.3 fails'
    else
      write(outunit,'(a)') 'test16.3 successful: '//trim(err_msg)
    endif

    call set_calendar_type(JULIAN)
    Time = set_date(1900, 1, 0, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.4 fails'
    else
      write(outunit,'(a)') 'test16.4 successful: '//trim(err_msg)
    endif

    call set_calendar_type(NOLEAP)
    Time = set_date(1900, 0, 1, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.5 fails'
    else
      write(outunit,'(a)') 'test16.5 successful: '//trim(err_msg)
    endif

    Time = set_date(1900, 1, 1, tick=11, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.6 fails'
    else
      write(outunit,'(a)') 'test16.6 successful: '//trim(err_msg)
    endif

    call set_calendar_type(THIRTY_DAY_MONTHS)
    Time = set_date(1900, 13, 1, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.7 fails'
    else
      write(outunit,'(a)') 'test16.7 successful: '//trim(err_msg)
    endif

    Time = set_date(1900, 12, 31, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.8 fails'
    else
      write(outunit,'(a)') 'test16.8 successful: '//trim(err_msg)
    endif

    call set_calendar_type(JULIAN)
    Time = set_date(1900, 4, 31, err_msg=err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test16.9 fails'
    else
      write(outunit,'(a)') 'test16.9 successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
 !  Tests of Gregorian calendar
 !  This test loops through every day of an 400 year period and writes a line to the output file for each day.

  if(test17) then
    write(outunit,'(/,a)') '#################################  test17  #################################'
    write(errunit,'(/,a)') ' ====================================================='
    write(errunit,'(a)')   '  Warning: test17 produces voluminous output.'
    write(errunit,'(a)')   '  It can be turned off with: &test_nml test17=.false./'
    write(errunit,'(a,/)') ' ====================================================='
    call set_calendar_type(GREGORIAN)
    do year=1801,2200
      leap = mod(year,4) == 0
      leap = leap .and. .not.mod(year,100) == 0
      leap = leap .or. mod(year,400) == 0
      do month=1,12
        days_this_month = days_per_month(month)
        if(leap .and. month == 2) days_this_month = 29
        do dday=1,days_this_month
          Time = set_date(year, month, dday, 0, 0, 0)
          call get_date(Time, yr, mo, day, hr, min, sec)
          write(outunit,100) yr, mo, day, leap_year(Time), days_in_month(Time), days_in_year(Time)
        enddo
      enddo
    enddo
  endif
  100 format('yr=',i4,' mo=',i2,' day=',i2,' leap=',L1,' days_in_month=',i2,' days_in_year=',i3)
 !==============================================================================================
 !  Tests of length_of_year

  if(test18) then
    write(outunit,'(/,a)') '#################################  test18  #################################'
    call set_calendar_type(THIRTY_DAY_MONTHS)
    call print_time(length_of_year(), 'length_of_year for THIRTY_DAY_MONTHS:', unit=outunit)
    call set_calendar_type(NOLEAP)
    call print_time(length_of_year(), 'length_of_year for NOLEAP:', unit=outunit)
    call set_calendar_type(JULIAN)
    call print_time(length_of_year(), 'length_of_year for JULIAN:', unit=outunit)
    call set_calendar_type(GREGORIAN)
    call print_time(length_of_year(), 'length_of_year for GREGORIAN:', unit=outunit)
  endif
 !==============================================================================================
 !  Tests of real_to_time_type

  if(test19) then
    write(outunit,'(/,a)') '#################################  test19  #################################'
    call print_time(real_to_time_type(86401.1), 'real_to_time_type(86401.1):', unit=outunit)
    Time = real_to_time_type(-1.0, err_msg)
    if(err_msg == '') then
      write(outunit,'(a)') 'test of real_to_time_type fails'
    else
      write(outunit,'(a)') 'test successful: '//trim(err_msg)
    endif
  endif
 !==============================================================================================
  write(outunit,'(/,a)') '############################################################################'
  write(outunit,'(a,i6)') ' ticks_per_second=',get_ticks_per_second()

 call fms_io_exit
 call fms_end
 end program test
#endif
