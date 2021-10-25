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

module time_interp_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Computes a weight and dates/indices for linearly interpolating between two dates.
! </OVERVIEW>

! <DESCRIPTION>
!     A time type is converted into two consecutive dates plus
!     a fraction representing the distance between the dates.
!     This information can be used to interpolate between the dates.
!     The dates may be expressed as years, months, or days or
!     as indices in an array.
! </DESCRIPTION>

! <PUBLIC>
!   Description summarizing public interface.
! </PUBLIC>

!-----------------------------------------------------------------------

use time_manager_mod, only: time_type, get_date, set_date, set_time, &
                            days_in_year, days_in_month, leap_year,  &
                            time_type_to_real, real_to_time_type,    &
                            get_calendar_type, JULIAN, GREGORIAN, NO_CALENDAR, &
                            operator(+), operator(-), operator(>),   &
                            operator(<), operator( // ), operator( / ),  &
                            operator(>=), operator(<=), operator( * ), &
                            operator(==), print_date, print_time

use          fms_mod, only: write_version_number, &
                            error_mesg, FATAL, stdout, stdlog, &
                            open_namelist_file, close_file, check_nml_error
use          mpp_mod, only: input_nml_file

implicit none
private

!-----------------------------------------------------------------------

public :: time_interp_init, time_interp, fraction_of_year

! <INTERFACE NAME="time_interp">

!   <OVERVIEW>
!      Returns a weight and dates or indices for interpolating between two dates. The
!      interface fraction_of_year is provided for backward compatibility with the
!      previous version. 
!   </OVERVIEW>
!   <DESCRIPTION>
!      Returns weight by interpolating Time between Time1 and Time2.
!      i.e. weight = (Time-Time1)/(Time2-Time1)
!      Time1 and Time2 may be specified by any of several different ways,
!      which is the reason for multiple interfaces.

!      If Time1 and Time2 are the begining and end of the year in which
!      Time falls, use first interface.

!      If Time1 and Time2 fall on year boundaries, use second interface.

!      If Time1 and Time2 fall on month boundaries, use third.

!      If Time1 and Time2 fall on day boundaries, use fourth.

!      If Time1 and Time2 are consecutive elements of an assending list, use fifth.
!      The fifth also returns the indices of Timelist between which Time falls.

!      The sixth interface is for cyclical data. Time_beg and Time_end specify the
!      begining and end of a repeating period. In this case:
!      weight = (Time_adjusted - Time1) / (Time2 - Time1)
!      Where:
!      Time1 = Timelist(index1)
!      Time2 = Timelist(index2)
!      Time_adjusted = Time - N*Period
!      Period = Time_end-Time_beg
!      N is between (Time-Time_end)/Period and (Time-Time_beg)/Period
!      That is, N is the integer that results in Time_adjusted that is between Time_beg and Time_end.
!      
!   </DESCRIPTION>
!   <TEMPLATE>
!     1. call time_interp( Time, weight )
!   </TEMPLATE>
!   <TEMPLATE>
!     2. call time_interp( Time, weight, year1, year2 )
!   </TEMPLATE>
!   <TEMPLATE>
!     3. call time_interp( Time, weight, year1, year2, month1, month2 )
!   </TEMPLATE>
!   <TEMPLATE>
!     4. call time_interp( Time, weight, year1, year2, month1, month2, day1, day2 )
!   </TEMPLATE>
!   <TEMPLATE>
!     5. call time_interp( Time, Timelist, weight, index1, index2 [, modtime] )
!   </TEMPLATE>
!   <TEMPLATE>
!     6. call time_interp( Time, Time_beg, Time_end, Timelist, weight, index1, index2 [,correct_leap_year_inconsistency])
!   </TEMPLATE>
!   <IN NAME="Time">
!      The time at which the the weight is computed.
!   </IN>
!   <IN NAME="Time_beg">
!      For cyclical interpolation: Time_beg specifies the begining time of a cycle.
!   </IN>
!   <IN NAME="Time_end">
!      For cyclical interpolation: Time_end specifies the ending time of a cycle.
!   </IN>
!   <IN NAME="Timelist">
!      For cyclical interpolation: Timelist is an array of times between Time_beg and Time_end.
!                                  Must be monotonically increasing.
!   </IN>
!   <IN NAME="modtime">
!   </IN>
!   <IN NAME="index1">
!      Timelist(index1) = The largest value of Timelist which is less than mod(Time,Time_end-Time_beg)
!   </IN>
!   <IN NAME="index2">
!      Timelist(index2) = The smallest value of Timelist which is greater than mod(Time,Time_end-Time_beg)
!   </IN>
!   <IN NAME="correct_leap_year_inconsistency">
!       Turns on a kluge for an inconsistency which may occur in a special case.
!       When the modulo time period (i.e. Time_end - Time_beg) is a whole number of years
!       and is not a multiple of 4, and the calendar in use has leap years, then it is
!       likely that the interpolation will involve mapping a common year onto a leap year.
!       In this case it is often desirable, but not absolutely necessary, to use data for
!       Feb 28 of the leap year when it is mapped onto a common year.
!       To turn this on, set correct_leap_year_inconsistency=.true.
!   </IN>
!   <OUT NAME="weight">
!     weight = (mod(Time,Time_end-Time_beg) - Timelist(index1)) / (Timelist(index2) - Timelist(index1))
!   </OUT>
!   <OUT NAME="year1"> </OUT>
!   <OUT NAME="year2"> </OUT>
!   <OUT NAME="month1"> </OUT>
!   <OUT NAME="month2"> </OUT>
!   <OUT NAME="day1"> </OUT>
!   <OUT NAME="day2"> </OUT>
!   <OUT NAME="index1"> </OUT>
!   <OUT NAME="index2"> </OUT>
!   <ERROR MSG="input time list not ascending order" STATUS="ERROR">
!     The list of input time types must have ascending dates.
!   </ERROR>
!   <ERROR MSG="modulo months must have same length" STATUS="ERROR">
!     The length of the current month for input Time and Time_list
!     must be the same when using the modulo month option. The
!     modulo month option is available but not supported. 
!   </ERROR>
!   <ERROR MSG="invalid value for argument modtime" STATUS="ERROR">
!     The optional argument modtime must have a value set by one
!     of the public parameters: NONE, YEAR, MONTH, DAY. The
!     MONTH and DAY options are available but not supported. 
!   </ERROR>
!   <ERROR MSG="period of list exceeds modulo period" STATUS="ERROR">
!     The difference between the last and first values in the input
!     Time list/array exceeds the length of the modulo period.
!   </ERROR>
!   <ERROR MSG="time before range of list or time after range of list" STATUS="ERROR">
!     The difference between the last and first values in the input
!     These errors occur when you are not using a modulo axis and
!     the input Time occurs before the first value in the Time
!     list/array or after the last value in the Time list/array. 
!   </ERROR>
!   <NOTE>
!     Examples: 
!     <PRE>
!       Time: Jan 01 00z    weight = 0.0 
!       Time: Jul 01        weight ~ 0.5 
!       Time: Dec 31 23z    weight ~ 1.0
!     </PRE>
!   </NOTE>

interface time_interp
    module procedure time_interp_frac,  time_interp_year, &
                     time_interp_month, time_interp_day,  &
                     time_interp_list,  time_interp_modulo
end interface
! </INTERFACE>

integer, public, parameter :: NONE=0, YEAR=1, MONTH=2, DAY=3

!-----------------------------------------------------------------------

   integer, parameter ::  secmin = 60, minhour = 60, hourday = 24,  &
                         sechour = secmin*minhour,                  &
                          secday = secmin*minhour*hourday

   integer, parameter :: monyear = 12
   integer, parameter :: halfday = secday/2

   integer :: yrmod, momod, dymod
   logical :: mod_leapyear

   character(len=128) :: version='$Id: time_interp.F90,v 19.0 2012/01/06 22:06:06 fms Exp $'
   character(len=128) :: tagname='$Name: siena_201211 $'

   logical :: module_is_initialized=.FALSE.
   logical :: perthlike_behavior=.FALSE.

   namelist / time_interp_nml / perthlike_behavior

contains


 subroutine time_interp_init()
   integer :: ierr, io, namelist_unit, logunit

   if ( module_is_initialized ) return

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, time_interp_nml, iostat=io)
      ierr = check_nml_error (io, 'time_interp_nml')
#else
   namelist_unit = open_namelist_file()
   ierr=1
   do while (ierr /= 0)
     read(namelist_unit, nml=time_interp_nml, iostat=io, end=20)
     ierr = check_nml_error (io, 'time_interp_nml')
   enddo
   20 call close_file (namelist_unit)
#endif

   call write_version_number( version, tagname )
   logunit = stdlog()
   write(logunit,time_interp_nml)

   module_is_initialized = .TRUE.

 end subroutine time_interp_init

!#######################################################################

! <SUBROUTINE NAME="time_interp_frac" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
! </SUBROUTINE>
!  returns the fractional time into the current year

 subroutine time_interp_frac ( Time, weight )

   type(time_type), intent(in)  :: Time 
   real           , intent(out) :: weight

   integer         :: year, month, day, hour, minute, second
   type(time_type) :: Year_beg, Year_end


   if ( .not. module_is_initialized ) call time_interp_init

!  ---- compute fractional time of year -----

     call get_date (Time, year, month, day, hour, minute, second) 

     Year_beg = set_date(year  , 1, 1) 
     Year_end = set_date(year+1, 1, 1)

     weight = (Time - Year_beg) // (Year_end - Year_beg)

 end subroutine time_interp_frac

!#######################################################################
! <SUBROUTINE NAME="fraction_of_year">
! <OVERVIEW>
!  Wrapper for backward compatibility
! </OVERVIEW>
! </SUBROUTINE>

 function fraction_of_year (Time)
 type(time_type), intent(in)  :: Time
 real :: fraction_of_year

  call time_interp_frac ( Time, fraction_of_year )

 end function fraction_of_year

!#######################################################################
! <SUBROUTINE NAME="time_interp_year" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="year1" TYPE="integer"> </OUT>
!   <OUT NAME="year2" TYPE="integer"> </OUT>
! </SUBROUTINE>
!  returns fractional time between mid points of consecutive years

 subroutine time_interp_year ( Time, weight, year1, year2 )

   type(time_type), intent(in)  :: Time
   real           , intent(out) :: weight
   integer        , intent(out) :: year1, year2

   integer :: year, month, day, hour, minute, second
   type (time_type) :: Mid_year, Mid_year1, Mid_year2


   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, year, month, day, hour, minute, second)

    ! mid point of current year
      Mid_year = year_midpt(year)

      if ( Time >= Mid_year ) then
    ! current time is after mid point of current year
           year1  = year
           year2  = year+1
           Mid_year2 = year_midpt(year2)
           weight = (Time - Mid_year) // (Mid_year2 - Mid_year)
      else
    ! current time is before mid point of current year
           year2  = year
           year1  = year-1
           Mid_year1 = year_midpt(year1)
           weight = (Time - Mid_year1) // (Mid_year - Mid_year1)
      endif

 end subroutine time_interp_year

!#######################################################################
! <SUBROUTINE NAME="time_interp_month" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="year1" TYPE="integer"> </OUT>
!   <OUT NAME="year2" TYPE="integer"> </OUT>
!   <OUT NAME="month1" TYPE="integer"> </OUT>
!   <OUT NAME="month2" TYPE="integer"> </OUT>
! </SUBROUTINE>
!  returns fractional time between mid points of consecutive months

 subroutine time_interp_month ( Time, weight, year1, year2, month1, month2 )

   type(time_type), intent(in)  :: Time
   real           , intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2

   integer :: year, month, day, hour, minute, second,  &
              mid_month, cur_month, mid1, mid2

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, year, month, day, hour, minute, second)

    ! mid point of current month in seconds
      mid_month = days_in_month(Time) * halfday
    ! time into current month in seconds
      cur_month = second + secmin*minute + sechour*hour + secday*(day-1)

      if ( cur_month >= mid_month ) then
    ! current time is after mid point of current month
           year1  = year;  month1 = month
           year2  = year;  month2 = month+1
           if (month2 > monyear)  year2 = year2+1
           if (month2 > monyear) month2 = 1
           mid1 = mid_month
           mid2 = days_in_month(set_date(year2,month2,2)) * halfday
           weight = real(cur_month - mid1) / real(mid1+mid2)
      else
    ! current time is before mid point of current month
           year2  = year;  month2 = month
           year1  = year;  month1 = month-1
           if (month1 < 1)  year1 = year1-1
           if (month1 < 1) month1 = monyear
           mid1 = days_in_month(set_date(year1,month1,2)) * halfday
           mid2 = mid_month
           weight = real(cur_month + mid1) / real(mid1+mid2)
      endif

 end subroutine time_interp_month

!#######################################################################
! <SUBROUTINE NAME="time_interp_day" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="year1" TYPE="integer"> </OUT>
!   <OUT NAME="year2" TYPE="integer"> </OUT>
!   <OUT NAME="month1" TYPE="integer"> </OUT>
!   <OUT NAME="month2" TYPE="integer"> </OUT>
!   <OUT NAME="day1" TYPE="integer"> </OUT>
!   <OUT NAME="day2" TYPE="integer"> </OUT>
! </SUBROUTINE>
!  returns fractional time between mid points of consecutive days

 subroutine time_interp_day ( Time, weight, year1, year2, month1, month2, day1, day2 )

   type(time_type), intent(in)  :: Time
   real           , intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2, day1, day2

   integer :: year, month, day, hour, minute, second, sday

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, year, month, day, hour, minute, second)

    ! time into current day in seconds
      sday = second + secmin*minute + sechour*hour

      if ( sday >= halfday ) then
    ! current time is after mid point of day
           year1 = year;  month1 = month;  day1 = day
           year2 = year;  month2 = month;  day2 = day + 1
           weight  = real(sday - halfday) / real(secday)

           if (day2 > days_in_month(Time)) then
               month2 = month2 + 1
               day2 = 1
               if (month2 > monyear) then
                    month2 = 1;  year2 = year2+1
               endif
           endif
      else
    ! current time is before mid point of day
           year2 = year;  month2 = month;  day2 = day
           year1 = year;  month1 = month;  day1 = day - 1
           weight  = real(sday + halfday) / real(secday)

           if (day1 < 1) then
               month1 = month1 - 1
               if (month1 < 1) then
                   month1 = monyear;  year1 = year1-1
               endif
               day1 = days_in_month(set_date(year1,month1,2))
           endif
      endif

 end subroutine time_interp_day

!#######################################################################
! <SUBROUTINE NAME="time_interp_modulo" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <IN NAME="Time_beg" TYPE="time_type"> </IN>
!   <IN NAME="Time_end" TYPE="time_type"> </IN>
!   <IN NAME="Timelist" TYPE="time_type" DIM="(:)"> </IN>
!   <IN NAME="correct_leap_year_inconsistency" TYPE="logical, optional" DEFAULT=".false.">
!       Turns on a kluge for an inconsistency which may occur in a special case.
!       When the modulo time period (i.e. Time_end - Time_beg) is a whole number of years
!       and is not a multiple of 4, and the calendar in use has leap years, then it is
!       likely that the interpolation will involve mapping a common year onto a leap year.
!       In this case it is often desirable, but not absolutely necessary, to use data for
!       Feb 28 of the leap year when it is mapped onto a common year.
!       To turn this on, set correct_leap_year_inconsistency=.true. </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="index1" TYPE="real"> </OUT>
!   <OUT NAME="index2" TYPE="real"> </OUT>
! </SUBROUTINE>

subroutine time_interp_modulo(Time, Time_beg, Time_end, Timelist, weight, index1, index2, &
                              correct_leap_year_inconsistency, err_msg)
type(time_type), intent(in)  :: Time, Time_beg, Time_end, Timelist(:)
real           , intent(out) :: weight
integer        , intent(out) :: index1, index2
logical, intent(in), optional :: correct_leap_year_inconsistency
character(len=*), intent(out), optional :: err_msg
  
  type(time_type) :: Period, T
  integer :: is, ie,i1,i2
  integer :: ys,ms,ds,hs,mins,ss ! components of the starting date
  integer :: ye,me,de,he,mine,se ! components of the ending date
  integer :: yt,mt,dt,ht,mint,st ! components of the current date
  integer :: dt1                 ! temporary value for day 
  integer :: n                   ! size of Timelist
  integer :: stdoutunit
  logical :: correct_lyr, calendar_has_leap_years, do_the_lyr_correction

  if ( .not. module_is_initialized ) call time_interp_init
  if( present(err_msg) ) err_msg = ''

  stdoutunit = stdout()
  n = size(Timelist)
  
  if (Time_beg>=Time_end) then
     if(present(err_msg)) then
        err_msg = "end of the specified time loop interval must be later than its beginning"
        return
     else
        call error_handler("end of the specified time loop interval must be later than its beginning")
     endif
  endif

  calendar_has_leap_years = (get_calendar_type() == JULIAN .or. get_calendar_type() == GREGORIAN)
  
  Period = Time_end-Time_beg ! period of the time axis

  if(present(correct_leap_year_inconsistency)) then
    correct_lyr = correct_leap_year_inconsistency
  else
    correct_lyr = .false.
  endif
  
  ! bring the requested time inside the specified time period
  T = Time

  do_the_lyr_correction = .false.

  ! Determine if the leap year correction needs to be done.
  ! It never needs to be done unless 3 conditions are met:
  ! 1) We are using a calendar with leap years
  ! 2) optional argument correct_leap_year_inconsistency is present and equals .true.
  ! 3) The modulo time period is an integer number of years
  ! If all of these are true then set do_the_lyr_correction to .true.

  if(calendar_has_leap_years .and. correct_lyr) then
    call get_date(Time_beg,ys,ms,ds,hs,mins,ss)
    call get_date(Time_end,ye,me,de,he,mine,se)
    if(ms==me.and.ds==de.and.hs==he.and.mins==mine.and.ss==se) then
      ! whole number of years
      do_the_lyr_correction = .true.
    endif
  endif

  if(do_the_lyr_correction) then
     call get_date(T,yt,mt,dt,ht,mint,st)
     yt = ys+modulo(yt-ys,ye-ys)
     dt1 = dt
     ! If it is Feb 29, but we map into a common year, use Feb 28
     if(mt==2.and.dt==29.and..not.leap_year(set_date(yt,1,1))) dt1=28
     T = set_date(yt,mt,dt1,ht,mint,st)
     if (T < Time_beg) then
       ! the requested time is within the first year, 
       ! but before the starting date. So we shift it to the last year.
       if(mt==2.and.dt==29.and..not.leap_year(set_date(ye,1,1))) dt=28
       T = set_date(ye,mt,dt,ht,mint,st)
     endif
  else
     do while ( T >= Time_end )
        T = T-Period
     enddo
     do while ( T < Time_beg )
        T = T+Period
     enddo
  endif
  
  ! find indices of the first and last records in the Timelist that are within 
  ! the requested time period.
  if (Time_end<=Timelist(1).or.Time_beg>=Timelist(n)) then
     if(get_calendar_type() == NO_CALENDAR) then
       call print_time(Time_beg,    'Time_beg'    )
       call print_time(Time_end,    'Time_end'    )
       call print_time(Timelist(1), 'Timelist(1)' )
       call print_time(Timelist(n), 'Timelist(n)' )
     else
       call print_date(Time_beg,    'Time_beg'    )
       call print_date(Time_end,    'Time_end'    )
       call print_date(Timelist(1), 'Timelist(1)' )
       call print_date(Timelist(n), 'Timelist(n)' )
     endif
     write(stdoutunit,*)'where n = size(Timelist) =',n
     if(present(err_msg)) then
        err_msg = 'the entire time list is outside the specified time loop interval'
        return
     else
        call error_handler('the entire time list is outside the specified time loop interval')
     endif
  endif
  
  call bisect(Timelist,Time_beg,index1=i1,index2=i2)
  if (i1 < 1) then
     is = 1 ! Time_beg before lower boundary
  else if (Time_beg == Timelist(i1)) then
     is = i1 ! Time_beg right on the lower boundary
  else
     is = i2 ! Time_beg inside the interval or on upper boundary
  endif
  call bisect(Timelist,Time_end,index1=i1,index2=i2)
  if (Time_end > Timelist(i1)) then
    ie = i1
  else if (Time_end == Timelist(i1)) then
    if(Time_beg == Timelist(is)) then
      ! Timelist includes time levels at both the lower and upper ends of the period.
      ! The endpoints of Timelist specify the same point in the cycle.
      ! This ambiguity is resolved by ignoring the last time level.
      ie = i1-1
    else
      ie = i1
    endif
  else
!   This should never happen because bisect does not return i1 such that Time_end < Timelist(i1)
  endif
  if (is>=ie) then
     if(get_calendar_type() == NO_CALENDAR) then
       call print_time(Time_beg,    'Time_beg   =')
       call print_time(Time_end,    'Time_end   =')
       call print_time(Timelist(1), 'Timelist(1)=')
       call print_time(Timelist(n), 'Timelist(n)=')
     else
       call print_date(Time_beg,    'Time_beg   =')
       call print_date(Time_end,    'Time_end   =')
       call print_date(Timelist(1), 'Timelist(1)=')
       call print_date(Timelist(n), 'Timelist(n)=')
     endif
     write(stdoutunit,*)'where n = size(Timelist) =',n
     write(stdoutunit,*)'is =',is,'ie =',ie
     if(present(err_msg)) then
        err_msg = 'error in calculation of time list bounds within the specified time loop interval'
        return
     else
        call error_handler('error in calculation of time list bounds within the specified time loop interval')
     endif
  endif
  
  ! handle special cases:
  if( T>=Timelist(ie) ) then
     ! time is after the end of the portion of the time list within the requested period
     index1 = ie;   index2 = is
     weight = (T-Timelist(ie))//(Period-(Timelist(ie)-Timelist(is)))
  else if (T<Timelist(is)) then
     ! time is before the beginning of the portion of the time list within the requested period
     index1 = ie;   index2 = is
     weight = 1.0-((Timelist(is)-T)//(Period-(Timelist(ie)-Timelist(is))))
  else
     call bisect(Timelist,T,index1,index2)
     weight = (T-Timelist(index1)) // (Timelist(index2)-Timelist(index1))
  endif

end subroutine time_interp_modulo

!#######################################################################
! given an array of times in ascending order and a specific time returns
! values of index1 and index2 such that the Timelist(index1)<=Time and
! Time<=Timelist(index2), and index2=index1+1
! index1=0, index2=1 or index=n, index2=n+1 are returned to indicate that 
! the time is out of range
subroutine bisect(Timelist,Time,index1,index2)
  type(time_type)  , intent(in)  :: Timelist(:)
  type(time_type)  , intent(in)  :: Time
  integer, optional, intent(out) :: index1, index2

  integer :: i,il,iu,n,i1,i2

  n = size(Timelist(:))
  
  if (Time==Timelist(1)) then
     i1 = 1 ; i2 = 2
  else if (Time==Timelist(n)) then
     i1 = n ; i2 = n+1
  else
     il = 0; iu=n+1
     do while(iu-il > 1)
        i = (iu+il)/2
        if(Timelist(i) > Time) then
           iu = i
        else
           il = i
        endif
     enddo
     i1 = il ; i2 = il+1
  endif

  if(PRESENT(index1)) index1 = i1
  if(PRESENT(index2)) index2 = i2
end subroutine bisect


!#######################################################################
! <SUBROUTINE NAME="time_interp_list" INTERFACE="time_interp">
!   <IN NAME="Time" TYPE="time_type" > </IN>
!   <IN NAME="Timelist" TYPE="time_type" DIM="(:)"> </IN>
!   <OUT NAME="weight" TYPE="real"> </OUT>
!   <OUT NAME="index1" TYPE="real"> </OUT>
!   <OUT NAME="index2" TYPE="real"> </OUT>
!   <IN NAME="modtime" TYPE="integer" > </IN>
! </SUBROUTINE>

subroutine time_interp_list ( Time, Timelist, weight, index1, index2, modtime, err_msg )
type(time_type)  , intent(in)  :: Time, Timelist(:)
real             , intent(out) :: weight
integer          , intent(out) :: index1, index2
integer, optional, intent(in)  :: modtime
character(len=*), intent(out), optional :: err_msg

integer :: n, hr, mn, se, mtime
type(time_type) :: T, Ts, Te, Td, Period, Time_mod

  if ( .not. module_is_initialized ) call time_interp_init

  if( present(err_msg) ) err_msg = ''

  weight = 0.; index1 = 0; index2 = 0
  n = size(Timelist(:))

! setup modular time axis?
  mtime = NONE
  if (present(modtime)) then
     mtime = modtime
     Time_mod = (Timelist(1)+Timelist(n))/2
     call get_date (Time_mod, yrmod, momod, dymod, hr, mn, se)
     mod_leapyear = leap_year(Time_mod)
  endif

! set period for modulo axis
  select case (mtime)
     case (NONE)
       ! do nothing
     case (YEAR)
         Period = set_time(0,days_in_year(Time_mod))
     case (MONTH)
       ! month length must be equal
         if (days_in_month(Time_mod) /= days_in_month(Time)) then
            if(present(err_msg)) then
               err_msg = 'modulo months must have same length'
               return
            else
               call error_handler ('modulo months must have same length')
            endif
         endif 
         Period = set_time(0,days_in_month(Time_mod))
     case (DAY)
         Period = set_time(0,1)
     case default
         call error_handler ('invalid value for argument modtime')
  end select

! If modulo time is in effect and Timelist spans a time interval exactly equal to 
! the modulo period, then the endpoints of Timelist specify the same point in the cycle.
! This ambiguity is resolved by ignoring the last time level.
  if (mtime /= NONE .and. Timelist(size(Timelist))-Timelist(1) == Period) then
     n = size(Timelist) - 1
  else
     n = size(Timelist)
  endif

! starting and ending times from list
  Ts = Timelist(1)
  Te = Timelist(n)
  Td = Te-Ts
  T  = set_modtime(Time,mtime)

! Check that Timelist does not span a time interval greater than the modulo period
  if (mtime /= NONE) then
     if (Td > Period) then
        if(present(err_msg)) then
           err_msg = 'period of list exceeds modulo period'
           return
        else
           call error_handler ('period of list exceeds modulo period')
        endif
     endif
  endif

! time falls on start or between start and end list values
  if ( T >= Ts .and. T < Te ) then
     call bisect(Timelist(1:n),T,index1,index2)
     weight = (T-Timelist(index1)) // (Timelist(index2)-Timelist(index1))

! time falls before starting list value
  else if ( T < Ts ) then
     if (mtime == NONE) then
        if(present(err_msg)) then
           err_msg = 'time before range of list'
           return
        else
           call error_handler ('time before range of list')
        endif
     endif
     Td = Te-Ts
     weight = 1. - ((Ts-T) // (Period-Td))
     index1 = n
     index2 = 1

! time falls on ending list value
  else if ( T == Te ) then
    if(perthlike_behavior) then
       weight = 1.0
       index1 = n-1
       index2 = n
    else
       weight = 0.
       index1 = n
       if (mtime == NONE) then
         index2 = n
       else
         index2 = 1
       endif
    endif

! time falls after ending list value
  else if ( T > Te ) then
     if (mtime == NONE) then
        if(present(err_msg)) then
           err_msg = 'time after range of list'
           return
        else
           call error_handler ('time after range of list')
        endif
     endif
     Td = Te-Ts
     weight = (T-Te) // (Period-Td)
     index1 = n
     index2 = 1
  endif

end subroutine time_interp_list

!#######################################################################
!  private routines
!#######################################################################

 function year_midpt (year)

   integer, intent(in) :: year
   type (time_type)    :: year_midpt, year_beg, year_end


   year_beg = set_date(year  , 1, 1)
   year_end = set_date(year+1, 1, 1)

   year_midpt = (year_beg + year_end) / 2
   
 end function year_midpt

!#######################################################################

 function month_midpt (year, month)

   integer, intent(in) :: year, month
   type (time_type)    :: month_midpt, month_beg, month_end

!  --- beginning of this month ---
   month_beg = set_date(year, month, 1)

!  --- start of next month ---
   if (month < 12) then
      month_end = set_date(year, month+1, 1)
   else
      month_end = set_date(year+1, 1, 1)
   endif

   month_midpt = (month_beg + month_end) / 2
   
 end function month_midpt

!#######################################################################

function set_modtime (Tin, modtime) result (Tout)
type(time_type), intent(in) :: Tin
integer, intent(in), optional :: modtime
type(time_type)             :: Tout
integer :: yr, mo, dy, hr, mn, se, mtime

  if(present(modtime)) then
    mtime = modtime
  else
    mtime = NONE
  endif

  select case (mtime)
    case (NONE)
       Tout = Tin
    case (YEAR)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod
        ! correct leap year dates
          if (.not.mod_leapyear .and. mo == 2 .and. dy > 28) then
             mo = 3; dy = dy-28
          endif
       Tout = set_date (yr, mo, dy, hr, mn, se)
    case (MONTH)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod; mo = momod
       Tout = set_date (yr, mo, dy, hr, mn, se)
    case (DAY)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod; mo = momod; dy = dymod
       Tout = set_date (yr, mo, dy, hr, mn, se)
  end select

end function set_modtime

!#######################################################################

subroutine error_handler (string)
character(len=*), intent(in) :: string

  call error_mesg ('time_interp_mod', trim(string), FATAL)

! write (*,'(a)') 'ERROR in time_interp: ' // trim(string)
! stop 111

end subroutine error_handler

!#######################################################################

end module time_interp_mod

! <INFO>

!   <ERROR MSG="input time list not ascending order" STATUS="">
!     The list of input time types must have ascending dates.
!   </ERROR> *
!   <ERROR MSG="modulo months must have same length" STATUS="">
!     The length of the current month for input Time and Time_list
!     must be the same when using the modulo month option.
!     The modulo month option is available but not supported.
!   </ERROR> *
!   <ERROR MSG="invalid value for argument modtime" STATUS="">
!     The optional argument modtime must have a value set by one
!     of the public parameters: NONE, YEAR, MONTH, DAY.
!     The MONTH and DAY options are available but not supported.
!   </ERROR> *
!   <ERROR MSG="period of list exceeds modulo period" STATUS="">
!     The difference between the last and first values in the
!     input Time list/array exceeds the length of the modulo period.
!   </ERROR> *
!   <ERROR MSG="time before range of list or time after range of list" STATUS="">
!     These errors occur when you are not using a modulo axis and the
!     input Time occurs before the first value in the Time list/array
!     or after the last value in the Time list/array.
!   </ERROR> *
!   <NOTE>
!   For all routines in this module the calendar type in module
!   time_manager must be set.
!   </NOTE>
!   <NOTE>
!     The following private parameters are set by this module:
! <PRE>
!           seconds per minute = 60
!           minutes per hour   = 60
!           hours   per day    = 24
!           months  per year   = 12
! </PRE>
!   </NOTE>

! </INFO>

#ifdef test_time_interp_
 program test_time_interp
 use          fms_mod, only: fms_init, fms_end, stdout, stdlog, FATAL, mpp_error
 use time_manager_mod, only: get_date, set_time, set_date, time_manager_init, set_calendar_type, operator(+)
 use time_manager_mod, only: JULIAN, time_type, increment_time, NOLEAP, print_date
 use  time_interp_mod, only: time_interp_init, time_interp, NONE, YEAR, MONTH, DAY

 implicit none

 integer, parameter :: num_Time=6
 type(time_type) :: Time_beg, Time_end, Time(num_Time)
 type(time_type), allocatable, dimension(:) :: Timelist
 integer :: index1, index2, mo, yr, timelist_len, outunit, ntest, nline
 real :: weight

 integer :: nmin, nmax

 namelist / test_time_interp_nml / timelist_len

 call fms_init
 outunit = stdout()
 call set_calendar_type(JULIAN)
 call time_interp_init

 Time_beg = set_date(1, 1, 1)
 Time_end = set_date(2, 1, 1)
 Time(1) = Time_beg
 Time(2) = set_date(1, 1,16)
 Time(3) = set_date(1, 2, 1)
 Time(4) = set_date(1,12, 1)
 Time(5) = set_date(1,12,16)
 Time(6) = Time_end

! Tests with modulo time
 do nline=1,3
   if(nline == 1) then
     allocate(Timelist(12))
     do mo=1,12
       Timelist(mo) = set_date(1, mo, 1)
     enddo
   else if(nline == 2) then
     allocate(Timelist(13))
     do mo=1,12
       Timelist(mo) = set_date(1, mo, 1)
     enddo
     Timelist(13) = set_date(2, 1, 1)
   else if(nline == 3) then
     allocate(Timelist(12))
     do mo=2,12
       Timelist(mo-1) = set_date(1, mo, 1)
     enddo
     Timelist(12) = set_date(2, 1, 1)
   endif

   do ntest=1,num_Time
     call diagram(nline,ntest,modulo_time=.true.)
     call time_interp(Time(ntest), Time_beg, Time_end, Timelist, weight, index1, index2)
     write(outunit,*) 'time_interp_modulo:'
     write(outunit,'()')
     call print_date(Time(ntest),                'Time       =')
     call print_date(Time_beg,                   'Time_beg   =')
     call print_date(Time_end,                   'Time_end   =')
     call print_date(Timelist(1),                'Timelist(1)=')
     call print_date(Timelist(size(Timelist(:))),'Timelist(n)=')
     write(outunit,99) index1,index2,weight
     write(outunit,'()')

     call time_interp(Time(ntest), Timelist, weight, index1, index2, modtime=YEAR)
     write(outunit,*) 'time_interp_list with modtime=YEAR:'
     write(outunit,'()')
     call print_date(Time(ntest),                'Time       =')
     call print_date(Timelist(1),                'Timelist(1)=')
     call print_date(Timelist(size(Timelist(:))),'Timelist(n)=')
     write(outunit,99) index1,index2,weight
   enddo
   deallocate(Timelist)
 enddo

! Tests without modulo time
 do nline=1,3
   if(nline == 1) then
     allocate(Timelist(12))
     do mo=1,12
       Timelist(mo) = set_date(1, mo, 1)
     enddo
   else if(nline == 2) then
     allocate(Timelist(13))
     do mo=1,12
       Timelist(mo) = set_date(1, mo, 1)
     enddo
     Timelist(13) = set_date(2, 1, 1)
   else if(nline == 3) then
     allocate(Timelist(12))
     do mo=2,12
       Timelist(mo-1) = set_date(1, mo, 1)
     enddo
     Timelist(12) = set_date(2, 1, 1)
   endif

   if(nline == 1) then
     nmin = 1; nmax = 4
   else if(nline == 2) then
     nmin = 1; nmax = num_Time
   else if(nline == 3) then
     nmin = 3; nmax = num_Time
   endif
   do ntest=nmin,nmax
     call diagram(nline,ntest,modulo_time=.false.)
     call time_interp(Time(ntest), Timelist, weight, index1, index2, modtime=NONE)
     write(outunit,*) 'time_interp_list with modtime=NONE:'
     write(outunit,'()')
     call print_date(Time(ntest),                'Time       =')
     call print_date(Timelist(1),                'Timelist(1)=')
     call print_date(Timelist(size(Timelist(:))),'Timelist(n)=')
     write(outunit,99) index1,index2,weight
   enddo
   deallocate(Timelist)
 enddo

! More tests with modulo time
 Time_beg = set_date(1999, 1, 1)
 Time_end = set_date(2000, 1, 1)
 Time(1)  = set_date(1998, 1, 1)
 Time(2)  = set_date(1998, 2,28)
 Time(3)  = set_date(1998,12,16)
 Time(4)  = set_date(2000, 1, 1)
 Time(5)  = set_date(2000, 2,28)
 Time(6)  = set_date(2000, 2,29)

 allocate(Timelist(13))
 do mo=1,12
   Timelist(mo) = set_date(1999, mo, 1)
 enddo
 Timelist(13) = set_date(2000, 1, 1)

 write(outunit,'("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>",/)')
 write(outunit,'()')
 write(outunit,*) 'time_interp_modulo with correct_leap_year_inconsistency=.true.'
 write(outunit,'()')
 write(outunit,'(" Jan 1 1999                                     Jan 1 2000")')
 write(outunit,'("    |                                               |")')
 write(outunit,'("    v                                               v")')
 write(outunit,'("    x---x---x---x---x---x---x---x---x---x---x---x---x")')
 write(outunit,'("    ^                                               ^")')
 write(outunit,'("    |                                               |")')
 write(outunit,'(" Time_beg                                        Time_end ")')
 write(outunit,'()')

 do ntest=1,num_Time
   call time_interp(Time(ntest), Time_beg, Time_end, Timelist, weight, index1, index2, correct_leap_year_inconsistency=.true.)
   call print_date(Time(ntest),' Time =')
   write(outunit,99) index1,index2,weight
   write(outunit,'()')
 enddo
 deallocate(Timelist)

! Tests of modulo time and leap year inconsistency
 Time_beg = set_date(1978, 1, 1)
 Time_end = set_date(1981, 1, 1)
 Time(1)  = set_date(1976, 2,28)
 Time(2)  = set_date(1976, 2,29)
 Time(3)  = set_date(1976, 3, 1)
 Time(4)  = set_date(1983, 2,28)
 Time(5)  = set_date(1983, 3, 1)
 Time(6)  = set_date(1981, 1, 1)
 allocate(Timelist(37))
 do yr=1978,1980
   do mo=1,12
     Timelist(12*(yr-1978)+mo) = set_date(yr, mo, 1)
   enddo
 enddo
 Timelist(37) = set_date(1981, 1, 1)

 write(outunit,'("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")')
 write(outunit,'()')
 write(outunit,*) 'time_interp_modulo with correct_leap_year_inconsistency=.true.'
 write(outunit,'()')
 write(outunit,'(" Jan 1 1978              Jan 1 1979              Jan 1 1980              Jan 1 1981")')
 write(outunit,'("     |                       |                       | <---- leap year ----> |")')
 write(outunit,'("     v                       v                       v                       v")')
 write(outunit,'("     x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x")')
 write(outunit,'("     ^                                                                       ^")')
 write(outunit,'("     |                                                                       |")')
 write(outunit,'("  Time_beg                                                               Time_end")')
 write(outunit,'()')

 do ntest=1,num_Time
   call time_interp(Time(ntest), Time_beg, Time_end, Timelist, weight, index1, index2, correct_leap_year_inconsistency=.true.)
   call print_date(Time(ntest),' Time=')
   write(outunit,99) index1,index2,weight
   write(outunit,'()')
 enddo
 deallocate(Timelist)

 allocate(Timelist(12))
 Timelist( 1) = set_date(1,  1, 16, hour=12) ! Jan midmonth
 Timelist( 2) = set_date(1,  2, 15, hour= 0) ! Feb midmonth (common year)
 Timelist( 3) = set_date(1,  3, 16, hour=12) ! Mar midmonth
 Timelist( 4) = set_date(1,  4, 16, hour= 0) ! Apr midmonth
 Timelist( 5) = set_date(1,  5, 16, hour=12) ! May midmonth
 Timelist( 6) = set_date(1,  6, 16, hour= 0) ! Jun midmonth
 Timelist( 7) = set_date(1,  7, 16, hour=12) ! Jul midmonth
 Timelist( 8) = set_date(1,  8, 16, hour=12) ! Aug midmonth
 Timelist( 9) = set_date(1,  9, 16, hour= 0) ! Sep midmonth
 Timelist(10) = set_date(1, 10, 16, hour=12) ! Oct midmonth
 Timelist(11) = set_date(1, 11, 16, hour= 0) ! Nov midmonth
 Timelist(12) = set_date(1, 12, 16, hour=12) ! Dec midmonth
 Time_beg = set_date(1, 1, 1)
 Time_end = set_date(2, 1, 1)
 call diagram(nline=4, ntest=0, modulo_time=.true.)
 do ntest=0,73
   Time(1) = set_date(1996, 1, 1) + set_time(seconds=0, days=5*ntest)
   call print_date(Time(1),' Time=')
   call time_interp(Time(1), Timelist, weight, index1, index2, modtime=YEAR)
   write(outunit,89) 'time_interp_list with modtime=YEAR:   ', index1,index2,weight
   call time_interp(Time(1), Time_beg, Time_end, Timelist, weight, index1, index2, correct_leap_year_inconsistency=.true.)
   write(outunit,89) 'time_interp_modulo: ', index1,index2,weight
   write(outunit,'()')
 enddo 

 99 format(' index1=',i3,'  index2=',i3,'  weight=',f18.15)
 89 format(a20,' index1=',i3,'  index2=',i3,'  weight=',f18.15)
 call fms_end

 contains

 subroutine diagram(nline,ntest,modulo_time)
 integer, intent(in) :: nline,ntest
 logical, intent(in) :: modulo_time

 write(outunit,'("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")')
 write(outunit,'()')
 if(modulo_time) then
   write(outunit,'(" Time_beg                                      Time_end")')
   write(outunit,'("  |                                               |")')
   write(outunit,'("  v                                               v")')
 endif

 if(nline == 1) then
   write(outunit,'("  x---x---x---x---x---x---x---x---x---x---x---x----")')
 else if(nline == 2) then
   write(outunit,'("  x---x---x---x---x---x---x---x---x---x---x---x---x")')
 else if(nline == 3) then
   write(outunit,'("  ----x---x---x---x---x---x---x---x---x---x---x---x")')
 else if(nline == 4) then
   write(outunit,'("  --x---x---x---x---x---x---x---x---x---x---x---x--")')
 endif

 if(ntest == 1) then
   write(outunit,'("  ^")  ')
   write(outunit,'("  |")  ')
   write(outunit,'(" Time")')
 else if(ntest == 2) then
   write(outunit,'("    ^")  ')
   write(outunit,'("    |")  ')
   write(outunit,'("   Time")')
 else if(ntest == 3) then
   write(outunit,'("      ^")  ')
   write(outunit,'("      |")  ')
   write(outunit,'("     Time")')
 else if(ntest == 4) then
   write(outunit,'("                                              ^")  ')
   write(outunit,'("                                              |")  ')
   write(outunit,'("                                             Time")')
 else if(ntest == 5) then
   write(outunit,'("                                                ^")  ')
   write(outunit,'("                                                |")  ')
   write(outunit,'("                                               Time")')
 else if(ntest == 6) then
   write(outunit,'("                                                  ^")  ')
   write(outunit,'("                                                  |")  ')
   write(outunit,'("                                                 Time")')
 endif
 write(outunit,'()')

 end subroutine diagram

 end program test_time_interp
#endif
