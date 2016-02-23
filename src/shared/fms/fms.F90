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

module fms_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   The fms module provides routines that are commonly used
!   by most FMS modules.
! </OVERVIEW>

! <DESCRIPTION>
!   Here is a summary of the functions performed by routines
!     in the fms module.
!
!   1. Output module version numbers to a common (<TT>log</TT>) file
!     using a common format.<BR/>
!   2. Open specific types of files common to many FMS modules.
!     These include namelist files, restart files, and 32-bit IEEE
!     data files. There also is a matching interface to close the files.
!     If other file types are needed the <TT>mpp_open</TT> and <TT>mpp_close</TT>
!     interfaces in module <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp_io.html">mpp_io</LINK> must be used.<BR/>
!    3. Read and write distributed data to simple native unformatted files.
!     This type of file (called a restart file) is used to checkpoint
!     model integrations for a subsequent restart of the run.<BR/>
!    4. For convenience there are several routines published from
!     the <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp.html">mpp</LINK> module. These are routines for getting processor
!     numbers, commonly used I/O unit numbers, error handling, and timing sections of code.
! </DESCRIPTION>

!-----------------------------------------------------------------------
!
!         A collection of commonly used routines.
!
!  The routines are primarily I/O related, however, there also
!  exists several simple miscellaneous utility routines.
!
!-----------------------------------------------------------------------
!
!  file_exist         Checks the existence of the given file name.
!
!  check_nml_error    Checks the iostat argument that is returned after
!                     reading a namelist and determines if the error
!                     code is valid.
!
!  write_version_number  Prints to the log file (or a specified unit)
!                        the (cvs) version id string and (cvs) tag name.
!
!  error_mesg          Print notes, warnings and error messages,
!                      terminates program for error messages.
!                      (use error levels NOTE,WARNING,FATAL)
!
!  open_namelist_file  Opens namelist file for reading only.
!
!  open_restart_file   Opens a file that will be used for reading or writing
!                      restart files with native unformatted data.
!
!  open_ieee32_file    Opens a file that will be used for reading or writing
!                      unformatted 32-bit ieee data.
!
!  close_file          Closes a file that was opened using
!                      open_namelist_file, open_restart_file, or
!                      open_ieee32_file.
!
!  set_domain          Call this routine to internally store in fms_mod the
!                      domain2d data type prior to calling the distributed
!                      data I/O routines read_data and write_data.
!
!  read_data           Reads distributed data from a single threaded file.
!
!  write_data          Writes distributed data to a single threaded file.
!
!  fms_init            Initializes the fms module and also the
!                      mpp_io module (which initializes all mpp mods).
!                      Will be called automatically if the user does
!                      not call it.
!
!  fms_end             Calls mpp exit routines.
!
!  lowercase           Convert character strings to all lower case
!
!  uppercase           Convert character strings to all upper case
!
!  monotonic_array     Determines if the real input array has
!                      monotonically increasing or decreasing values.
!
!  string_array_index  Match the input character string to a string
!                      in an array/list of character strings.
!
!-----------------------------------------------------------------------
!---- published routines from mpp_mod ----
!
!   mpp_error, NOTE, WARNING, FATAL
!   mpp_error_state
!   mpp_pe, mpp_npes, mpp_root_pe
!   stdin, stdout, stderr, stdlog
!   mpp_chksum
!
!   mpp_clock_id, mpp_clock_begin , mpp_clock_end
!   MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
!   CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER,
!   CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
!
!-----------------------------------------------------------------------

use          mpp_mod, only:  mpp_error, NOTE, WARNING, FATAL,    &
                             mpp_set_warn_level,                 &
                             mpp_transmit, ALL_PES,              &
                             mpp_pe, mpp_npes, mpp_root_pe,      &
                             mpp_sync, mpp_chksum,               &
                             mpp_clock_begin, mpp_clock_end,     &
                             mpp_clock_id, mpp_init, mpp_exit,   &
                             MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, &
                             CLOCK_COMPONENT, CLOCK_SUBCOMPONENT,&
                             CLOCK_MODULE_DRIVER, CLOCK_MODULE,  &
                             CLOCK_ROUTINE, CLOCK_LOOP,          &
                             CLOCK_INFRA, mpp_clock_set_grain,   &
                             mpp_set_stack_size,                 &
                             stdin, stdout, stderr, stdlog,      &
                             mpp_error_state, lowercase,         &
                             uppercase, mpp_broadcast, input_nml_file

use  mpp_domains_mod, only:  domain2D, mpp_define_domains, &
                             mpp_update_domains, GLOBAL_DATA_DOMAIN, &
                             mpp_domains_init, mpp_domains_exit,     &
                             mpp_global_field, mpp_domains_set_stack_size,  &
                             mpp_get_compute_domain, mpp_get_global_domain, &
                             mpp_get_data_domain

use       mpp_io_mod, only:  mpp_io_init, mpp_open, mpp_close,         &
                       MPP_ASCII, MPP_NATIVE, MPP_IEEE32, MPP_NETCDF,  &
                       MPP_RDONLY, MPP_WRONLY, MPP_APPEND, MPP_OVERWR, &
                       MPP_SEQUENTIAL, MPP_DIRECT,                     &
                       MPP_SINGLE, MPP_MULTI, MPP_DELETE, mpp_io_exit, &
                       fieldtype, mpp_get_atts, mpp_get_info, mpp_get_fields

use fms_io_mod, only : read_data, write_data, fms_io_init, fms_io_exit, field_size, &
                       open_namelist_file, open_restart_file, open_ieee32_file, close_file, &
                       set_domain, get_domain_decomp, nullify_domain, &
                       open_file, open_direct_file, string, get_mosaic_tile_grid, &
                       get_mosaic_tile_file, get_global_att_value, file_exist, field_exist

use memutils_mod, only: print_memuse_stats, memutils_init
!jp constants wants to import fms_mod, this creates a circular dependency
!use constants_mod, only: constants_version=>version, constants_tagname=>tagname !pjp: PI not computed


implicit none
private

! routines for initialization and termination of module
public :: fms_init, fms_end

! routines for opening/closing specific types of file
public :: open_namelist_file, open_restart_file, &
          open_ieee32_file, close_file, &
          open_file, open_direct_file

! routines for reading/writing distributed data
public :: set_domain, read_data, write_data
public :: get_domain_decomp, field_size, nullify_domain
public :: get_global_att_value

! routines for get mosaic information
public :: get_mosaic_tile_grid, get_mosaic_tile_file

! miscellaneous i/o routines
public :: file_exist, check_nml_error, field_exist,     &
          write_version_number, error_mesg, fms_error_handler

! miscellaneous utilities (non i/o)
public :: lowercase, uppercase, string,        &
          string_array_index, monotonic_array

! public mpp interfaces
public :: mpp_error, NOTE, WARNING, FATAL, &
          mpp_error_state,                 &
          mpp_pe, mpp_npes, mpp_root_pe,   &
          stdin, stdout, stderr, stdlog,   &
          mpp_chksum
public :: mpp_clock_id, mpp_clock_begin, mpp_clock_end
public :: MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
public :: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, &
          CLOCK_MODULE_DRIVER, CLOCK_MODULE,   &
          CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA

!Balaji
!this is published by fms and applied to any initialized clocks
!of course you can go and set the flag to SYNC or DETAILED by hand
integer, public :: clock_flag_default

! Namelist read error values
  TYPE nml_errors_type
     INTEGER :: multipleNMLSinFile
     INTEGER :: badType1
     INTEGER :: badType2
     INTEGER :: missingVar
  END TYPE nml_errors_type
  TYPE(nml_errors_type), SAVE :: nml_errors


!------ namelist interface -------
!------ adjustable severity level for warnings ------

  logical           :: read_all_pe   = .true.
  character(len=16) :: clock_grain = 'NONE', clock_flags='NONE'
  character(len=8)  :: warning_level = 'warning'
  character(len=64) :: iospec_ieee32 = '-N ieee_32'
  integer           :: stack_size = 0
  integer           :: domains_stack_size = 0
  logical, public   :: print_memory_usage = .FALSE.

!------ namelist interface -------

! <NAMELIST NAME="fms_nml">
!   <DATA NAME="clock_grain"  TYPE="character"  DEFAULT="'NONE'">
!     The level of clock granularity used for performance timing sections
!     of code. Possible values in order of increasing detail are:
!     'NONE', 'COMPONENT', 'SUBCOMPONENT', 'MODULE_DRIVER', 'MODULE', 'ROUTINE',
!     'LOOP', and 'INFRA'.  Code sections are defined using routines in MPP
!     module: mpp_clock_id, mpp_clock_begin, and mpp_clock_end.
!     The fms module makes these routines public.
!     A list of timed code sections will be printed to STDOUT.
!     See the <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp.html">MPP</LINK>
!     module for more details.
!   </DATA>
!   <DATA NAME="clock_flags"  TYPE="character"  DEFAULT="'NONE'">
!     Possible values are 'NONE', 'SYNC', or 'DETAILED'.
!     SYNC will give accurate information on load balance of the clocked
!     portion of code.
!     DETAILED also turns on detailed message-passing performance diagnosis.
!     Both SYNC and DETAILED will  work correctly on innermost clock nest
!     and distort outer clocks, and possibly the overall code time.
!     See the <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp.html">MPP</LINK>
!     module for more details.
!   </DATA>
!   <DATA NAME="read_all_pe"  TYPE="logical"  DEFAULT="true">
!     Read global data on all processors extracting local part needed (TRUE) or
!     read global data on PE0 and broadcast to all PEs (FALSE).
!   </DATA>
!   <DATA NAME="warning_level"  TYPE="character"  DEFAULT="'warning'">
!     Sets the termination condition for the WARNING flag to interfaces
!     error_mesg/mpp_error. set warning_level = 'fatal' (program crashes for
!     warning messages) or 'warning' (prints warning message and continues).
!   </DATA>
!   <DATA NAME="iospec_ieee32"  TYPE="character"  DEFAULT="'-N ieee_32'">
!     iospec flag used with the open_ieee32_file interface.
!   </DATA>
!   <DATA NAME="stack_size"  TYPE="integer"  DEFAULT="0">
!     The size in words of the MPP user stack. If stack_size > 0, the following
!     MPP routine is called: call mpp_set_stack_size (stack_size). If stack_size
!     = 0 (default) then the default size set by mpp_mod is used.
!   </DATA>
!   <DATA NAME="domains_stack_size" TYPE="integer"  DEFAULT="0">
!     The size in words of the MPP_DOMAINS user stack. If
!     domains_stack_size > 0, the following MPP_DOMAINS routine is called:
!     call mpp_domains_set_stack_size (domains_stack_size). If
!     domains_stack_size = 0 (default) then the default size set by
!     mpp_domains_mod is used.
!   </DATA>
!   <DATA NAME="print_memory_usage"  TYPE="logical"  DEFAULT=".FALSE.">
!     If set to .TRUE., memory usage statistics will be printed at various
!     points in the code. It is used to study memory usage, e.g to detect
!     memory leaks.
!   </DATA>
! </NAMELIST>

  namelist /fms_nml/  read_all_pe, clock_grain, clock_flags,    &
                      warning_level, iospec_ieee32, &
                      stack_size, domains_stack_size, &
                      print_memory_usage

!   ---- private data for check_nml_error ----

   integer, private :: num_nml_error_codes, nml_error_codes(20)
   logical, private :: do_nml_error_init = .true.
   private  nml_error_init


!  ---- version number -----

  character(len=128) :: version = '$Id: fms.F90,v 19.0 2012/01/06 21:57:13 fms Exp $'
  character(len=128) :: tagname = '$Name: siena_201211 $'

  logical :: module_is_initialized = .FALSE.


contains

!#######################################################################

! <SUBROUTINE NAME="fms_init">

!   <OVERVIEW>
!     Initializes the FMS module and also calls the initialization routines for all
!     modules in the MPP package. Will be called automatically if the user does
!     not call it.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Initialization routine for the fms module. It also calls initialization routines
!      for the mpp, mpp_domains, and mpp_io modules. Although this routine
!      will be called automatically by other fms_mod routines, users should
!      explicitly call fms_init. If this routine is called more than once it will
!      return silently. There are no arguments.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call fms_init ( )
!   </TEMPLATE>


!   <ERROR MSG="invalid entry for namelist variable warning_level" STATUS="FATAL">
!     The namelist variable warning_level must be either 'fatal' or 'warning'
!     (case-insensitive).
!   </ERROR>
!   <ERROR MSG="invalid entry for namelist variable clock_grain" STATUS="FATAL">
!     The namelist variable clock_grain must be one of the following values:
!     'NONE', 'COMPONENT', 'SUBCOMPONENT', 'MODULE_DRIVER', 'MODULE', 'ROUTINE',
!     'LOOP', or 'INFRA' (case-insensitive).
!   </ERROR>

! initializes the fms module/package
! also calls mpp initialization routines and reads fms namelist

subroutine fms_init (localcomm )
 integer, intent(in), optional :: localcomm
 integer :: unit, ierr, io

    if (module_is_initialized) return    ! return silently if already called
    module_is_initialized = .true.
!---- initialize mpp routines ----
    if(present(localcomm)) then
       call mpp_init(localcomm=localcomm)
    else
       call mpp_init()
    endif
    call mpp_domains_init
    call fms_io_init

!---- read namelist input ----

    call nml_error_init  ! first initialize namelist iostat error codes

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, fms_nml, iostat=io)
#else
    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
          read  (unit, nml=fms_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'fms_nml')  ! also initializes nml error codes
       enddo
 10    call mpp_close (unit)
    endif
#endif

!---- define mpp stack sizes if non-zero -----

    if (        stack_size > 0) call         mpp_set_stack_size (        stack_size)
    if (domains_stack_size > 0) call mpp_domains_set_stack_size (domains_stack_size)

!---- set severity level for warnings ----

    select case( trim(lowercase(warning_level)) )
    case( 'fatal' )
        call mpp_set_warn_level ( FATAL )
    case( 'warning' )
        call mpp_set_warn_level ( WARNING )
    case default
        call error_mesg ( 'fms_init',  &
             'invalid entry for namelist variable warning_level', FATAL )
    end select

!--- set granularity for timing code sections ---

    select case( trim(uppercase(clock_grain)) )
    case( 'NONE' )
        call mpp_clock_set_grain (0)
    case( 'COMPONENT' )
        call mpp_clock_set_grain (CLOCK_COMPONENT)
    case( 'SUBCOMPONENT' )
        call mpp_clock_set_grain (CLOCK_SUBCOMPONENT)
    case( 'MODULE_DRIVER' )
        call mpp_clock_set_grain (CLOCK_MODULE_DRIVER)
    case( 'MODULE' )
        call mpp_clock_set_grain (CLOCK_MODULE)
    case( 'ROUTINE' )
        call mpp_clock_set_grain (CLOCK_ROUTINE)
    case( 'LOOP' )
        call mpp_clock_set_grain (CLOCK_LOOP)
    case( 'INFRA' )
        call mpp_clock_set_grain (CLOCK_INFRA)
    case default
        call error_mesg ( 'fms_init',  &
             'invalid entry for namelist variable clock_grain', FATAL )
    end select
!Balaji
    select case( trim(uppercase(clock_flags)) )
    case( 'NONE' )
       clock_flag_default = 0
    case( 'SYNC' )
       clock_flag_default = MPP_CLOCK_SYNC
    case( 'DETAILED' )
       clock_flag_default = MPP_CLOCK_DETAILED
    case default
       call error_mesg ( 'fms_init',  &
            'invalid entry for namelist variable clock_flags', FATAL )
   end select

!--- write version info and namelist to logfile ---

    call write_version_number (version, tagname)
    if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog()
      write (unit, nml=fms_nml)
      write (unit,*) 'nml_error_codes=', nml_error_codes(1:num_nml_error_codes)
    endif

    call memutils_init( print_memory_usage )
    call print_memuse_stats('fms_init')

    !call write_version_number (constants_version,constants_tagname)

end subroutine fms_init
! </SUBROUTINE>

!#######################################################################


! <SUBROUTINE NAME="fms_end">

!   <OVERVIEW>
!     Calls the termination routines for all modules in the MPP package.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Termination routine for the fms module. It also calls destructor routines
!      for the mpp, mpp_domains, and mpp_io modules. If this routine is called
!      more than once it will return silently. There are no arguments.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call fms_end ( )
!   </TEMPLATE>

! terminates the fms module/package
! also calls mpp destructor routines

subroutine fms_end ( )

    if (.not.module_is_initialized) return  ! return silently
!    call fms_io_exit  ! now called from coupler_end
    call mpp_io_exit
    call mpp_domains_exit
    call mpp_exit
    module_is_initialized =.FALSE.

end subroutine fms_end
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="error_mesg">

!   <OVERVIEW>
!     Print notes, warnings and error messages; terminates program for warning
!     and error messages. (use error levels NOTE,WARNING,FATAL, see example below)
!   </OVERVIEW>
!   <DESCRIPTION>
!     Print notes, warnings and error messages; and terminates the program for
!     error messages. This routine is a wrapper around mpp_error, and is provided
!     for backward compatibility. This module also publishes mpp_error,
!      <B>users should try to use the mpp_error interface</B>.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call error_mesg ( routine, message, level )
!   </TEMPLATE>

!   <IN NAME="routine"  TYPE="character" >
!     Routine name where the warning or error has occurred.
!   </IN>
!   <IN NAME="message"  TYPE="character" >
!     Warning or error message to be printed.
!   </IN>
!   <IN NAME="level"  TYPE="integer" >
!     Level of severity; set to NOTE, WARNING, or FATAL Termination always occurs
!     for FATAL, never for NOTE, and is settable for WARNING (see namelist).
!   </IN>
!   <NOTE>
!
!     Examples:
!     <PRE>
!        use fms_mod, only: error_mesg, FATAL, NOTE

!        call error_mesg ('fms_mod', 'initialization not called', FATAL)
!        call error_mesg ('fms_mod', 'fms_mod message', NOTE)
!     </PRE>
!   </NOTE>
! wrapper for the mpp error handler
! users should try to use the mpp_error interface

 subroutine error_mesg (routine, message, level)
  character(len=*), intent(in) :: routine, message
  integer,          intent(in) :: level

!  input:
!      routine   name of the calling routine (character string)
!      message   message written to output   (character string)
!      level     set to NOTE, MESSAGE, or FATAL (integer)

    if (.not.module_is_initialized) call fms_init ( )
    call mpp_error ( routine, message, level )

 end subroutine error_mesg
! </SUBROUTINE>

!#######################################################################
! <FUNCTION NAME="fms_error_handler">

!   <OVERVIEW>
!     Facilitates the control of fatal error conditions
!   </OVERVIEW>
!   <DESCRIPTION>
!     When err_msg is present, message is copied into err_msg
!     and the function returns a value of .true.
!     Otherwise calls mpp_error to terminate execution.
!     The intended use is as shown below.
!   </DESCRIPTION>
!   <TEMPLATE>
!     if(fms_error_handler(routine, message, err_msg)) return
!   </TEMPLATE>
!   <IN NAME="routine"  TYPE="character">
!     Routine name where the fatal error has occurred.
!   </IN>
!   <IN NAME="message"  TYPE="character">
!     fatal error message to be printed.
!   </IN>
!   <OUT NAME="fms_error_handler"  TYPE="logical">
!     .true.  when err_msg is present
!     .false. when err_msg is not present
!   </OUT>
!   <OUT NAME="err_msg"  TYPE="character">
!     When err_msg is present: err_msg = message
!   </OUT>

 function fms_error_handler(routine, message, err_msg)

 logical :: fms_error_handler
 character(len=*), intent(in) :: routine, message
 character(len=*), intent(out), optional :: err_msg

 fms_error_handler = .false.
 if(present(err_msg)) then
   err_msg = message
   fms_error_handler = .true.
 else
   call mpp_error(trim(routine),trim(message),FATAL)
 endif

 end function fms_error_handler
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME="check_nml_error">

!   <OVERVIEW>
!     Checks the iostat argument that is returned after reading a namelist
!     and determines if the error code is valid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     The FMS allows multiple namelist records to reside in the same file.
!     Use this interface to check the iostat argument that is returned after
!     reading a record from the namelist file. If an invalid iostat value
!     is detected this routine will produce a fatal error. See the NOTE below.
!   </DESCRIPTION>
!   <TEMPLATE>
!     check_nml_error ( iostat, nml_name )
!   </TEMPLATE>

!   <IN NAME="iostat"  TYPE="integer" >
!     The iostat value returned when reading a namelist record.
!   </IN>
!   <IN NAME="nml_name"  TYPE="character" >
!     The name of the namelist. This name will be printed if an error is
!     encountered, otherwise the name is not used.
!   </IN>
!   <OUT NAME=""  TYPE="integer" >
!     This function returns the input iostat value (integer) if it is an
!     allowable error code. If the iostat error code is not
!     allowable, an error message is printed and the program terminated.
!   </OUT>
!   <NOTE>
!     Some compilers will return non-zero iostat values when reading through
!     files with multiple namelist. This routine
!     will try skip these errors and only terminate for true namelist errors.
!
!     Examples
!
!       The following example checks if a file exists, reads a namelist input
!       from that file, and checks for errors in that
!       namelist. When the correct namelist is read and it has no errors the
!       routine check_nml_error will return zero and the while loop will exit.
!       This code segment should be used to read namelist files.
!       <PRE>
!          integer :: unit, ierr, io
!
!          if ( file_exist('input.nml') ) then
!              unit = open_namelist_file ( )
!              ierr=1
!              do while (ierr > 0)
!                read  (unit, nml=moist_processes_nml, iostat=io)
!                ierr = check_nml_error(io,'moist_processes_nml')
!              enddo
!              call close_file (unit)
!          endif
!       </PRE>
!   </NOTE>

!   <ERROR MSG="Unknown error while reading namelist ...., (IOSTAT = ####)" STATUS="FATAL">
!     There was an error reading the namelist specified. Carefully examine all namelist and variables
!     for anything incorrect (e.g. malformed, hidden characters).
!   </ERROR>
!   <ERROR MSG="Unknown namelist, or mistyped namelist variable in namelist ...., (IOSTAT = ####)" STATUS="FATAL">
!     The name list given doesn't exist in the namelist file, or a variable in the namelist is mistyped or isn't a
!     namelist variable.
!   </ERROR>

! used to check the iostat argument that is
! returned after reading a namelist
! see the online documentation for how this routine might be used
  INTEGER FUNCTION check_nml_error(IOSTAT, NML_NAME)
    INTEGER, INTENT(in) :: IOSTAT
    CHARACTER(len=*), INTENT(in) :: NML_NAME

    CHARACTER(len=256) :: err_str

    IF ( .NOT.module_is_initialized) CALL fms_init()

    check_nml_error = IOSTAT

    ! Return on valid IOSTAT values
    IF ( IOSTAT <= 0 .OR. IOSTAT == nml_errors%multipleNMLSinFile ) RETURN

    ! Everything else is a FATAL
    IF ( mpp_pe() == mpp_root_pe() ) THEN
       IF ( (IOSTAT == nml_errors%badType1 .OR. IOSTAT == nml_errors%badType2) .OR. IOSTAT == nml_errors%missingVar ) THEN
          WRITE (err_str,*) 'Unknown namelist, or mistyped namelist variable in namelist ',TRIM(NML_NAME),', (IOSTAT = ',IOSTAT,')'
          CALL error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
          CALL mpp_sync()
       ELSE
          WRITE (err_str,*) 'Unknown error while reading namelist ',TRIM(NML_NAME),', (IOSTAT = ',IOSTAT,')'
          CALL error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
          CALL mpp_sync()
       END IF
    ELSE
       CALL mpp_sync()
    END IF
  END FUNCTION check_nml_error
! </FUNCTION>

!-----------------------------------------------------------------------
!   private routine for initializing allowable error codes

  SUBROUTINE nml_error_init
    ! Determines the IOSTAT error value for some common Namelist errors.
    ! Also checks if the compiler returns a non-zero status if there are
    ! multiple namelist records in a single file.
    INTEGER, PARAMETER :: unit_begin = 20, unit_end = 1024
    INTEGER :: fileunit, io_stat
    INTEGER, DIMENSION(4) :: nml_iostats
    LOGICAL :: opened

    ! Variables for sample namelists
    INTEGER :: i1, i2
    REAL :: r1, r2
    LOGICAL :: l1
    NAMELIST /a_nml/ i1, r1
    NAMELIST /b_nml/ i2, r2, l1
    NAMELIST /badType1_nml/ i1, r1
    NAMELIST /badType2_nml/ i1, r1
    NAMELIST /missingVar_nml/ i2, r2

    ! Initialize the sample namelist variables
    i1 = 1
    i2 = 2
    r1 = 1.0
    r2 = 2.0
    l1 = .FALSE.

    ! Create a dummy namelist file
    IF ( mpp_pe() == mpp_root_pe() ) THEN
       ! Find a free file unit for a scratch file
       file_opened: DO fileunit = unit_begin, unit_end
          INQUIRE(UNIT=fileunit, OPENED=opened)
          IF ( .NOT.opened ) EXIT file_opened
       END DO file_opened

#if defined __PGI
       OPEN (UNIT=fileunit, FILE='_read_error.nml', IOSTAT=io_stat)
#else
       OPEN (UNIT=fileunit, STATUS='SCRATCH', IOSTAT=io_stat)
#endif

       ! Write sample namelist to the SCRATCH file.
       WRITE (UNIT=fileunit, NML=a_nml, IOSTAT=io_stat)
       WRITE (UNIT=fileunit, NML=b_nml, IOSTAT=io_stat)
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&badType1_nml  i1=1, r1=''bad'' /",/)')
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&badType2_nml  i1=1, r1=.true. /",/)')
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&missingVar_nml  i2=1, r2=1.0e0, l1=.true. /",/)')

       ! Rewind for reading
       REWIND(UNIT=fileunit)

       ! Read the second namelist from the file -- check for namelist bug
       READ (UNIT=fileunit, NML=b_nml, IOSTAT=nml_iostats(1))
       REWIND(UNIT=fileunit)

       ! Read in bad type 1 --- Some compilers treat the string cast differently
       READ (UNIT=fileunit, NML=badType1_nml, IOSTAT=nml_iostats(2))
       REWIND(UNIT=fileunit)

       ! Read in bad type 2
       READ (UNIT=fileunit, NML=badType2_nml, IOSTAT=nml_iostats(3))
       REWIND(UNIT=fileunit)

       ! Read in missing variable/misstyped
       READ (UNIT=fileunit, NML=missingVar_nml, IOSTAT=nml_iostats(4))

       ! Done, close file
       CLOSE (UNIT=fileunit)

       ! Some compilers don't handle the type casting as well as we would like.
       IF ( nml_iostats(2) * nml_iostats(3) .EQ. 0 ) THEN
          IF ( nml_iostats(2) .NE. 0 .AND. nml_iostats(3) .EQ. 0 ) THEN
             nml_iostats(3) = nml_iostats(2)
          ELSE IF ( nml_iostats(2) .EQ. 0 .AND. nml_iostats(3) .NE.0 ) THEN
             nml_iostats(2) = nml_iostats(3)
          ELSE
             nml_iostats(2) = nml_iostats(4)
             nml_iostats(2) = nml_iostats(4)
          END IF
       END IF
    END IF

    ! Broadcast nml_errors
    CALL mpp_broadcast(nml_iostats,4,mpp_root_pe())
    nml_errors%multipleNMLSinFile = nml_iostats(1)
    nml_errors%badType1 = nml_iostats(2)
    nml_errors%badType2 = nml_iostats(3)
    nml_errors%missingVar = nml_iostats(4)

    do_nml_error_init = .FALSE.
  END SUBROUTINE nml_error_init

!#######################################################################
! <SUBROUTINE NAME="write_version_number">

!   <OVERVIEW>
!     Prints to the log file (or a specified unit) the (cvs) version id string and
!     (cvs) tag name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Prints to the log file (stdlog) or a specified unit the (cvs) version id string
!      and (cvs) tag name.
!   </DESCRIPTION>
!   <TEMPLATE>
!    call write_version_number ( version [, tag, unit] )
!   </TEMPLATE>

!   <IN NAME="version" TYPE="character(len=*)">
!    string that contains routine name and version number.
!   </IN>
!   <IN NAME="tag" TYPE="character(len=*)">
!    The tag/name string, this is usually the Name string
!    returned by CVS when checking out the code.
!   </IN>
!   <IN NAME="unit" TYPE="integer">
!    The Fortran unit number of an open formatted file. If this unit number
!    is not supplied the log file unit number is used (stdlog).
!   </IN>
! prints module version number to the log file of specified unit number

 subroutine write_version_number (version, tag, unit)

!   in:  version = string that contains routine name and version number
!
!   optional in:
!        tag = cvs tag name that code was checked out with
!        unit    = alternate unit number to direct output
!                  (default: unit=stdlog)

   character(len=*), intent(in) :: version
   character(len=*), intent(in), optional :: tag
   integer,          intent(in), optional :: unit

   integer :: logunit

   if (.not.module_is_initialized) call fms_init ( )

     logunit = stdlog()
     if (present(unit)) then
         logunit = unit
     else
       ! only allow stdlog messages on root pe
         if ( mpp_pe() /= mpp_root_pe() ) return
     endif

     if (present(tag)) then
         write (logunit,'(/,80("="),/(a))') trim(version), trim(tag)
     else
         write (logunit,'(/,80("="),/(a))') trim(version)
     endif

 end subroutine write_version_number
! </SUBROUTINE>

!#######################################################################


! <FUNCTION NAME="string_array_index">

!   <OVERVIEW>
!     match the input character string to a string
!     in an array/list of character strings
!   </OVERVIEW>
!   <DESCRIPTION>
!      Tries to find a match for a character string in a list of character strings.
!      The match is case sensitive and disregards blank characters to the right of
!      the string.
!   </DESCRIPTION>
!   <TEMPLATE>
!      string_array_index ( string, string_array [, index] )
!   </TEMPLATE>

!   <IN NAME="string"  TYPE="character(len=*), scalar" >
!     Character string of arbitrary length.
!   </IN>
!   <IN NAME="string_array"  TYPE="character(len=*)" DIM="(:)">
!     Array/list of character strings.
!   </IN>
!   <OUT NAME="index"  TYPE="integer" >
!     The index of string_array where the first match was found. If
!            no match was found then index = 0.
!   </OUT>
!   <OUT NAME="string_array_index"  TYPE="logical" >
!     If an exact match was found then TRUE is returned, otherwise FALSE is returned.
!   </OUT>
!   <NOTE>
!     Examples
!      <PRE>
!       string = "def"
!       string_array = (/ "abcd", "def ", "fghi" /)

!       string_array_index ( string, string_array, index )

!       Returns: TRUE, index = 2
!      </PRE>
!   </NOTE>
! match the input character string to a string
! in an array/list of character strings

function string_array_index ( string, string_array, index ) result (found)
character(len=*),  intent(in)  :: string, string_array(:)
integer, optional, intent(out) :: index
logical :: found
integer :: i

! initialize this function to false
! loop thru string_array and exit when a match is found

  found = .false.
  if (present(index)) index = 0

  do i = 1, size(string_array(:))
    ! found a string match ?
    if ( trim(string) == trim(string_array(i)) ) then
         found = .true.
         if (present(index)) index = i
         exit
    endif
  enddo

end function string_array_index
! </FUNCTION>

!#######################################################################

! <FUNCTION NAME="monotonic_array">

!   <OVERVIEW>
!     Determines if a real input array has monotonically increasing or
!     decreasing values.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Determines if the real input array has monotonically increasing or
!     decreasing values.
!   </DESCRIPTION>
!   <TEMPLATE>
!     monotonic_array ( array [, direction] )
!   </TEMPLATE>

!   <IN NAME="array"  TYPE="real" DIM="(:)">
!     An array of real values. If the size(array) < 2 this function
!     assumes the array is not monotonic, no fatal error will occur.
!   </IN>
!   <OUT NAME="direction"  TYPE="integer" >
!     If the input array is:
!                >> monotonic (small to large) then direction = +1.
!                >> monotonic (large to small) then direction = -1.
!                >> not monotonic then direction = 0.
!   </OUT>
!   <OUT NAME="monotonic_array"  TYPE="logical" >
!     If the input array of real values either increases or decreases monotonically
!      then TRUE is returned, otherwise FALSE is returned.
!   </OUT>
! determines if the real input array has
! monotonically increasing or decreasing values

function monotonic_array ( array, direction )
real,    intent(in)            :: array(:)
integer, intent(out), optional :: direction
logical :: monotonic_array
integer :: i

! initialize
  monotonic_array = .false.
  if (present(direction)) direction = 0

! array too short
  if ( size(array(:)) < 2 ) return

! ascending
  if ( array(1) < array(size(array(:))) ) then
     do i = 2, size(array(:))
       if (array(i-1) < array(i)) cycle
       return
     enddo
     monotonic_array = .true.
     if (present(direction)) direction = +1

! descending
  else
     do i = 2, size(array(:))
       if (array(i-1) > array(i)) cycle
       return
     enddo
     monotonic_array = .true.
     if (present(direction)) direction = -1
  endif

end function monotonic_array
! </FUNCTION>

end module fms_mod
! <INFO>
!   <BUG>
!     Namelist error checking may not work correctly with some compilers.
!
!     Users should beware when mixing Fortran reads and read_data calls. If a
!     Fortran read follows read_data and namelist variable read_all_pe = FALSE
!     (not the default), then the code will fail. It is safest if Fortran reads
!     precede calls to read_data.
!   </BUG>
!   <ERROR MSG="unexpected EOF" STATUS="FATAL">
!     An unexpected end-of-file was encountered in a read_data call.
!     You may want to use the optional end argument to detect the EOF.
!   </ERROR>
!   <NOTE>
!     1) If the <B>MPP</B> or <B>MPP_DOMAINS</B> stack size is exceeded the
!     program will terminate after printing the required size.
!
!     2) When running on a very small number of processors or for high
!     resolution models the default domains_stack_size will
!     probably be insufficient.
!
!     3) The following performance routines in the <B>MPP</B> module are published by this module.
!<PRE>
!        mpp_clock_id, mpp_clock_begin, mpp_clock_end
!</PRE>
!        and associated parameters that are published:
!<PRE>
!        MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, CLOCK_COMPONENT, CLOCK_SUBCOMPONENT,
!        CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
!</PRE>
!
!     4) Here is an example of how to time a section of code.<BR/>
!<PRE>
!          use fms_mod, only: mpp_clock_id, mpp_clock_begin, &
!                             mpp_clock_end. MPP_CLOCK_SYNC, &
!                             CLOCK_MODULE_DRIVER
!          integer :: id_mycode
!
!          id_mycode = mpp_clock_id ('mycode loop', flags=MPP_CLOCK_SYNC, grain=CLOCK_MODULE_DRIVER)
!          call mpp_clock_begin (id_mycode)
!                        :
!                        :
!           ~~ this code will be timed ~~
!                        :
!                        :
!          call mpp_clock_end (id_mycode)
! </PRE>
!        Note: <TT>CLOCK_MODULE_DRIVER</TT> can be replaced with
!        CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE,
!        CLOCK_LOOP, or CLOCK_INFRA.
!
!   </NOTE>
!   <FUTURE>
!     NetCDF facilities for reading and writing restart files and (IEEE32)
!       data files.
!    </FUTURE>
!    <FUTURE>
!     May possible split the FMS module into two modules.
!
!      i.general utilities (FMS_MOD) <BR/>
!     ii.I/O utilities (FMS_IO_MOD)
!    </FUTURE>
! </INFO>

