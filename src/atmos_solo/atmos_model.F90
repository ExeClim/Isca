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

program atmos_model

!-----------------------------------------------------------------------
!
!  Main program for running a stand-alone atmospheric dynamical core.
!
!-----------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use   atmosphere_mod, only: atmosphere_init, atmosphere_end, atmosphere, atmosphere_domain

use time_manager_mod, only: time_type, set_time, get_time,  &
                            operator(+), operator (<), operator (>), &
                            operator (/=), operator (/), operator (*),&
			    THIRTY_DAY_MONTHS, JULIAN,            &
                            NOLEAP, NO_CALENDAR, set_calendar_type, &
			    set_date, get_date

use          fms_mod, only: file_exist, check_nml_error,                &
                            error_mesg, FATAL, WARNING,                 &
                            mpp_pe, mpp_root_pe, fms_init, fms_end,     &
                            stdlog, stdout, write_version_number,       &
                            open_restart_file,                          &
                            mpp_clock_id, mpp_clock_begin,              &
                            mpp_clock_end, CLOCK_COMPONENT, set_domain, &
                            nullify_domain, uppercase
use       fms_io_mod, only: fms_io_exit

use  mpp_mod,         only: mpp_set_current_pelist
use  mpp_domains_mod, only: domain2d
use       mpp_io_mod, only: mpp_open, mpp_close, MPP_ASCII, MPP_OVERWR, &
                            MPP_SEQUENTIAL, MPP_SINGLE, MPP_RDONLY, MPP_DELETE

use diag_manager_mod, only: diag_manager_init, diag_manager_end, get_base_date

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: register_tracers
use       memutils_mod, only: print_memuse_stats
use   constants_mod,    only: SECONDS_PER_HOUR,  SECONDS_PER_MINUTE, constants_init

implicit none

!-----------------------------------------------------------------------

character(len=128), parameter :: version = &
'$Id: atmos_model.F90,v 17.0.2.1.2.1.4.1.4.1 2012/05/16 13:35:00 z1l Exp $'

character(len=128), parameter :: tag = &
'$Name: siena_201211 $'

!-----------------------------------------------------------------------
!       ----- model time -----
! there is no calendar associated with model of this type
! therefore, year=0, month=0 are assumed
!s Mima 2013 has changed this - a calendar is required to run the rrtm radiation properly.

   type (time_type) :: Time, Time_init, Time_end, Time_step_atmos
   integer :: num_atmos_calls, na

! ----- model initial date -----

   integer :: date_init(6) ! note: year=month=0
   integer :: calendar_type = 0

! ----- timing flags -----

   integer :: id_init, id_loop, id_end
   integer, parameter :: timing_level = 1

!-----------------------------------------------------------------------
   character(len=80) :: text
!-----------------------------------------------------------------------
   type(domain2d), save :: atmos_domain  ! This variable must be treated as read-only
!-----------------------------------------------------------------------

      integer, dimension(4) :: current_time = (/ 0, 0, 0, 0 /)
      integer :: days=0, hours=0, minutes=0, seconds=0
      integer :: dt_atmos = 0
      integer :: memuse_interval = 72
      logical :: print_memuse = .true.
      integer :: atmos_nthreads = 1
      character(len=17) :: calendar = '                 '
      integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)

      namelist /main_nml/ current_date, dt_atmos,  &
                          days, hours, minutes, seconds, memuse_interval, print_memuse, atmos_nthreads, calendar, current_time

!#######################################################################
 call constants_init
 call fms_init ( )
 call atmos_model_init

!   ------ atmosphere integration loop -------

    call mpp_clock_begin (id_loop)

    do na = 1, num_atmos_calls

       call atmosphere (Time)

       Time = Time + Time_step_atmos

       if(modulo(na,memuse_interval) == 0 .and. print_memuse) then
         write( text,'(a,i4)' )'Main loop at timestep=',na
         call print_memuse_stats(text)
       endif

    enddo

    call mpp_clock_end (id_loop)

!   ------ end of atmospheric time step loop -----

 call atmos_model_end
 call fms_io_exit
 call fms_end

contains

!#######################################################################

   subroutine atmos_model_init

!-----------------------------------------------------------------------
    integer :: unit, ierr, io, logunit
    integer :: ntrace, ntprog, ntdiag, ntfamily
    integer :: date(6)
    type (time_type) :: Run_length
    integer :: omp_get_thread_num, get_cpu_affinity, base_cpu
!-----------------------------------------------------------------------
!----- initialization timing identifiers ----

 id_init = mpp_clock_id ('MAIN: initialization', grain=CLOCK_COMPONENT)
 id_loop = mpp_clock_id ('MAIN: time loop'     , grain=CLOCK_COMPONENT)
 id_end  = mpp_clock_id ('MAIN: termination'   , grain=CLOCK_COMPONENT)

 logunit = stdlog()

 call mpp_clock_begin (id_init)

!-------------------------------------------
! how many tracers have been registered?
!  (will print number below)
   call register_tracers ( MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfamily )


!----- read namelist -------

#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=main_nml, iostat=io)
     ierr = check_nml_error(io, 'main_nml')
#else
   unit = open_namelist_file ( )
   ierr=1; do while (ierr /= 0)
          read  (unit, nml=main_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'main_nml')
   enddo
10 call mpp_close (unit)
#endif

!----- write namelist to logfile -----

   call write_version_number (version,tag)
   if ( mpp_pe() == mpp_root_pe() ) write (logunit, nml=main_nml)

   if (dt_atmos == 0) then
     call error_mesg ('program atmos_model', 'dt_atmos has not been specified', FATAL)
   endif

!----- read restart file -----

   if (file_exist('INPUT/atmos_model.res')) then
       call mpp_open (unit, 'INPUT/atmos_model.res', action=MPP_RDONLY, nohdrs=.true.)
       read  (unit,*) date
       read  (unit,*) calendar_type
       call mpp_close (unit)
   else
    ! use namelist time if restart file does not exist
        select case( uppercase(trim(calendar)) )
        case( 'JULIAN' )
            calendar_type = JULIAN
        case( 'NOLEAP' )
            calendar_type = NOLEAP
        case( 'THIRTY_DAY' )
            calendar_type = THIRTY_DAY_MONTHS
        case( 'NO_CALENDAR' )
            calendar_type = NO_CALENDAR
        case default
            call error_mesg( 'program atmos_model', 'main_nml entry calendar must be one of JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.', FATAL )
        end select

       if(calendar_type /= NO_CALENDAR) then !s If using no_calendar we need to set the date using the current time only.
            date = current_date
       else
            date(1:2) = 0
            date(3:6) = current_time
       endif

   endif

   call set_calendar_type (calendar_type)

!----- write current/initial date actually used to logfile file -----

    if ( mpp_pe() == mpp_root_pe() ) then
      write (logunit,16) date(3:6)
    endif

 16 format ('  current time used = day',i5,' hour',i3,2(':',i2.2))

!  print number of tracers to logfile
   if (mpp_pe() == mpp_root_pe()) then
        write (logunit, '(a,i3)') 'Number of tracers =', ntrace
        write (logunit, '(a,i3)') 'Number of prognostic tracers =', ntprog
        write (logunit, '(a,i3)') 'Number of diagnostic tracers =', ntdiag
   endif

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

    call diag_manager_init

!----- always override initial/base date with diag_manager value -----

!----- get the base date in the diag_table from the diag_manager ----
!      this base date is typically the starting date for the
!      experiment and is subtracted from the current date

    call get_base_date ( date_init(1), date_init(2), date_init(3), &
                         date_init(4), date_init(5), date_init(6)  )

  ! make sure base date does not have a year or month specified !s but only with NO_CALENDAR as we DO want a year or month specified for MiMA.
    if ( calendar_type == NO_CALENDAR .and. date_init(1)+date_init(2) /= 0 ) then
         call error_mesg ('program atmos_model', 'invalid base base - &
                          &must have year = month = 0', FATAL)
    endif

!----- set initial and current time types ------
!----- set run length and compute ending time -----
#ifdef MARS_GCM
!               Dont allow minutes in the Mars model
    date_init(5)= 0.0
#endif MARS_GCM

     if ( calendar_type /= NO_CALENDAR) then
!s New way of setting Time_init and Time uses set_date, which subtracts 1 from date_init(3) and date(3) when calculating dates. Leads to correct start date when calendar /= NO_CALENDAR
         Time_init = set_date (date_init(1), date_init(2), date_init(3), &
              date_init(4), date_init(5), date_init(6))

         Time      = set_date (date(1), date(2), date(3),  &
              date(4), date(5), date(6))
     else
!s Works with NO_CALENDAR
         Time_init  = set_time(date_init(4)*int(SECONDS_PER_HOUR)+date_init(5)*int(SECONDS_PER_MINUTE)+date_init(6),date_init(3))
         Time       = set_time(date     (4)*int(SECONDS_PER_HOUR)+date     (5)*int(SECONDS_PER_MINUTE)+date     (6),date     (3))
     endif

    Run_length = set_time(       hours*int(SECONDS_PER_HOUR)+     minutes*int(SECONDS_PER_MINUTE)+     seconds,days        )
    Time_end   = Time + Run_length

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

      call mpp_open (unit, 'time_stamp.out', form=MPP_ASCII, action=MPP_OVERWR, &
                     access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date

!     compute ending time in days,hours,minutes,seconds
      call get_time ( Time_end, date(6), date(3) )  ! gets sec,days
      date(4) = date(6)/int(SECONDS_PER_HOUR); date(6) = date(6) - date(4)*int(SECONDS_PER_HOUR)
#ifdef MARS_GCM
      date(5) = 0                                ; date(6) = date(6) - date(5)*int(SECONDS_PER_MINUTE)
#else
      date(5) = date(6)/int(SECONDS_PER_MINUTE)  ; date(6) = date(6) - date(5)*int(SECONDS_PER_MINUTE)
#endif
      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date

      call mpp_close (unit)

  20  format (6i7,2x,'day')   ! can handle day <= 999999

!-----------------------------------------------------------------------
!--- compute the time steps ---
!    determine number of iterations through the time integration loop
!    must be evenly divisible

      Time_step_atmos = set_time (dt_atmos,0)
      num_atmos_calls = Run_length / Time_step_atmos

!-----------------------------------------------------------------------
!----- initial (base) time must not be greater than current time -----

   if ( Time_init > Time ) call error_mesg ('program atmos_model',  &
                   'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of atmos time step ------

   if ( num_atmos_calls * Time_step_atmos /= Run_length )  &
        call error_mesg ('program atmos_model',  &
           'run length must be multiple of atmosphere time step', FATAL)

!-----------------------------------------------------------------------
!------ initialize atmospheric model ------

!      call omp_set_num_threads(atmos_nthreads)
      if (mpp_pe() .eq. mpp_root_pe()) then
        unit=stdout()
        write(unit,*) ' starting ',atmos_nthreads,' OpenMP threads per MPI-task'
        call flush(unit)
      endif
      base_cpu = get_cpu_affinity()
!$OMP PARALLEL
!      call set_cpu_affinity(base_cpu + omp_get_thread_num())
#ifdef DEBUG
      write(6,*) 'PE: ',mpp_pe(),'  thread_num', omp_get_thread_num(),'  affinity:',get_cpu_affinity()
      call flush(6)
#endif
!$OMP END PARALLEL

      call atmosphere_init (Time_init, Time, Time_step_atmos)
      call atmosphere_domain(atmos_domain)

!-----------------------------------------------------------------------
!   open and close dummy file in restart dir to check if dir exists
      call mpp_set_current_pelist()
      call mpp_open  (unit, 'RESTART/file' )
      call mpp_close (unit, action=MPP_DELETE)

!  ---- terminate timing ----
   call mpp_clock_end (id_init)

!-----------------------------------------------------------------------

   call print_memuse_stats('atmos_model_init')
   end subroutine atmos_model_init

!#######################################################################

   subroutine atmos_model_end

   integer :: unit, date(6)
!-----------------------------------------------------------------------
   call mpp_clock_begin (id_end)

      call atmosphere_end

!----- compute current time in days,hours,minutes,seconds -----

     if( calendar_type /= NO_CALENDAR) then
!s Updated call to get final date in line with updating initial date setting above.
         call get_date (Time, date(1), date(2), date(3),  &
              date(4), date(5), date(6))
     else
         date(1:2) = 0
         call get_time ( Time, date(6), date(3) )
         date(4) = date(6)/int(SECONDS_PER_HOUR); date(6) = date(6) - date(4)*int(SECONDS_PER_HOUR)
         date(5) = date(6)/int(SECONDS_PER_MINUTE); date(6) = date(6) - date(5)*int(SECONDS_PER_MINUTE)
     endif
#ifdef MARS_GCM
      date(5) = 0                              ; date(6) = date(6) - date(5)*int(SECONDS_PER_MINUTE)
!#else
!      date(5) = date(6)/int(SECONDS_PER_MINUTE); date(6) = date(6) - date(5)*int(SECONDS_PER_MINUTE)
#endif MARS_GCM

!----- check time versus expected ending time ----

      if (Time /= Time_end) call error_mesg ('program atmos_model',  &
              'final time does not match expected ending time', WARNING)

!----- write restart file ------

      if ( mpp_pe() == mpp_root_pe() ) then
           call mpp_open (unit, 'RESTART/atmos_model.res', form=MPP_ASCII, action=MPP_OVERWR, &
                          access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )
           write (unit,'(6i6,8x,a)') date, &
                 'Current model time: year, month, day, hour, minute, second'
           write (unit,'(i6,8x,a)') calendar_type, &
                 '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
           call mpp_close (unit)
      endif

!----- final output of diagnostic fields ----
      call set_domain(atmos_domain)  ! This assumes all output fields are on the atmos domain

      call diag_manager_end (Time)

      call nullify_domain()

      call mpp_clock_end (id_end)
!-----------------------------------------------------------------------

   end subroutine atmos_model_end

!#######################################################################
! routines to set/get date when no calendar is set (i.e., yr=0 and mo=0)
!#######################################################################

end program atmos_model
