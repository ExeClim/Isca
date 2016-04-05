
! -*-f90-*-
! $Id: mpp_util.inc,v 18.0 2010/03/02 23:58:05 fms Exp $


#ifdef use_libSMA
#include <mpp_util_sma.inc>
#elif defined(use_libMPI)
#include <mpp_util_mpi.inc>
#else
#include <mpp_util_nocomm.inc>
#endif

  !#####################################################################
  ! <FUNCTION NAME="stdin">
  !  <OVERVIEW>
  !    Standard fortran unit numbers.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This function returns the current standard fortran unit numbers for input.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   stdin()
  !  </TEMPLATE>
  ! </FUNCTION>
  function stdin()
    integer :: stdin
    stdin = in_unit
    return
  end function stdin

  !#####################################################################
  ! <FUNCTION NAME="stdout">
  !  <OVERVIEW>
  !    Standard fortran unit numbers.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This function returns the current  standard fortran unit numbers for output.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   stdout()
  !  </TEMPLATE>
  ! </FUNCTION>
  function stdout()
    integer :: stdout
    stdout = out_unit
    if( pe.NE.root_pe )stdout = stdlog()
    return
  end function stdout

  !#####################################################################
  ! <FUNCTION NAME="stderr">
  !  <OVERVIEW>
  !    Standard fortran unit numbers.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This function returns the current standard fortran unit numbers for error messages.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   stderr()
  !  </TEMPLATE>
  ! </FUNCTION>
  function stderr()
    integer :: stderr
    stderr = err_unit
    return
  end function stderr

  !#####################################################################
  ! <FUNCTION NAME="stdlog">
  !  <OVERVIEW>
  !    Standard fortran unit numbers.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This function returns the current  standard fortran unit numbers for log messages.
  !    Log messages, by convention, are written to the file <TT>logfile.out</TT>.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    stdlog()
  !  </TEMPLATE>
  ! </FUNCTION>
  function stdlog()
    integer :: stdlog,istat
    logical :: opened
    character(len=11) :: this_pe
    if( pe.EQ.root_pe )then
       write(this_pe,'(a,i6.6,a)') '.',pe,'.out'
       inquire( file=trim(configfile)//this_pe, opened=opened )
       if( opened )then
          call FLUSH(log_unit)
       else
          log_unit=get_unit()
          open( unit=log_unit, status='UNKNOWN', file=trim(configfile)//this_pe, position='APPEND', err=10 )
       end if
       stdlog = log_unit
    else
       inquire( unit=etc_unit, opened=opened )
       if( opened )then
          call FLUSH(etc_unit)
       else
          open( unit=etc_unit, status='UNKNOWN', file=trim(etcfile), position='APPEND', err=11 )
       end if
       stdlog = etc_unit
    end if
    return
10  call mpp_error( FATAL, 'STDLOG: unable to open '//trim(configfile)//this_pe//'.' )
11  call mpp_error( FATAL, 'STDLOG: unable to open '//trim(etcfile)//'.' )
  end function stdlog

  !#####################################################################
  subroutine mpp_init_logfile()
  integer :: p
  logical :: exist
  character(len=11) :: this_pe
  if( pe.EQ.root_pe )then
     log_unit = get_unit()
     do p=0,npes-1
       write(this_pe,'(a,i6.6,a)') '.',p,'.out'
       inquire( file=trim(configfile)//this_pe, exist=exist )
       if(exist)then
         open( unit=log_unit, file=trim(configfile)//this_pe, status='REPLACE' )
         close(log_unit)
       endif
     end do
  end if
  end subroutine mpp_init_logfile
  !#####################################################################
  subroutine mpp_set_warn_level(flag)
    integer, intent(in) :: flag

    if( flag.EQ.WARNING )then
       warnings_are_fatal = .FALSE.
    else if( flag.EQ.FATAL )then
       warnings_are_fatal = .TRUE.
    else
       call mpp_error( FATAL, 'MPP_SET_WARN_LEVEL: warning flag must be set to WARNING or FATAL.' )
    end if
    return
  end subroutine mpp_set_warn_level

  !#####################################################################
  function mpp_error_state()
    integer :: mpp_error_state
    mpp_error_state = error_state
    return
  end function mpp_error_state

!#####################################################################
!overloads to mpp_error_basic
!support for error_mesg routine in FMS
subroutine mpp_error_mesg( routine, errormsg, errortype )
  character(len=*), intent(in) :: routine, errormsg
  integer,          intent(in) :: errortype

  call mpp_error( errortype, trim(routine)//': '//trim(errormsg) )
  return
end subroutine mpp_error_mesg

!#####################################################################
subroutine mpp_error_noargs()
  call mpp_error(FATAL)
end subroutine mpp_error_noargs

!#####################################################################
subroutine mpp_error_Is(errortype, errormsg1, value, errormsg2)
  integer,          intent(in) :: errortype
  INTEGER,          intent(in) :: value
  character(len=*), intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  call mpp_error( errortype, errormsg1, (/value/), errormsg2)
end subroutine mpp_error_Is
!#####################################################################
subroutine mpp_error_Rs(errortype, errormsg1, value, errormsg2)
  integer,          intent(in) :: errortype
  REAL,             intent(in) :: value
  character(len=*), intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  call mpp_error( errortype, errormsg1, (/value/), errormsg2)
end subroutine mpp_error_Rs
!#####################################################################
subroutine mpp_error_Ia(errortype, errormsg1, array, errormsg2)
  integer,               intent(in) :: errortype
  INTEGER, dimension(:), intent(in) :: array
  character(len=*),      intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array))
  if(present(errormsg2)) string = trim(string)//errormsg2
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_Ia

!#####################################################################
subroutine mpp_error_Ra(errortype, errormsg1, array, errormsg2)
  integer,            intent(in) :: errortype
  REAL, dimension(:), intent(in) :: array
  character(len=*),      intent(in) :: errormsg1
  character(len=*),   intent(in), optional :: errormsg2
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array))
  if(present(errormsg2)) string = trim(string)//errormsg2
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_Ra

!#####################################################################
#define _SUBNAME_ mpp_error_ia_ia
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ integer
#include <mpp_error_a_a.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_ia_ra
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ real
#include <mpp_error_a_a.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_ra_ia
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ integer
#include <mpp_error_a_a.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_ra_ra
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ real
#include <mpp_error_a_a.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_ia_is
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ integer
#include <mpp_error_a_s.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_ia_rs
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ real
#include <mpp_error_a_s.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_ra_is
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ integer
#include <mpp_error_a_s.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_ra_rs
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ real
#include <mpp_error_a_s.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_is_ia
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ integer
#include <mpp_error_s_a.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_is_ra
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ real
#include <mpp_error_s_a.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_rs_ia
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ integer
#include <mpp_error_s_a.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_rs_ra
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ real
#include <mpp_error_s_a.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_is_is
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ integer
#include <mpp_error_s_s.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_is_rs
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ real
#include <mpp_error_s_s.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_rs_is
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ integer
#include <mpp_error_s_s.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
#define _SUBNAME_ mpp_error_rs_rs
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ real
#include <mpp_error_s_s.h>
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
!#####################################################################
function iarray_to_char(iarray) result(string)
integer, intent(in) :: iarray(:)
character(len=256) :: string
character(len=32)  :: chtmp
integer :: i, len_tmp, len_string

 string = ''
 do i=1,size(iarray)
   write(chtmp,'(i16)') iarray(i)
   chtmp = adjustl(chtmp)
   len_tmp = len_trim(chtmp)
   len_string  = len_trim(string)
   string(len_string+1:len_string+len_tmp) = trim(chtmp)
   string(len_string+len_tmp+1:len_string+len_tmp+1) = ','
 enddo
 len_string = len_trim(string)
 string(len_string:len_string) = ' ' ! remove trailing comma

end function iarray_to_char
!#####################################################################
function rarray_to_char(rarray) result(string)
real, intent(in) :: rarray(:)
character(len=256) :: string
character(len=32)  :: chtmp
integer :: i, len_tmp, len_string

 string = ''
 do i=1,size(rarray)
   write(chtmp,'(G16.9)') rarray(i)
   chtmp = adjustl(chtmp)
   len_tmp = len_trim(chtmp)
   len_string  = len_trim(string)
   string(len_string+1:len_string+len_tmp) = trim(chtmp)
   string(len_string+len_tmp+1:len_string+len_tmp+1) = ','
 enddo
 len_string = len_trim(string)
 string(len_string:len_string) = ' ' ! remove trailing comma

end function rarray_to_char

  !#####################################################################
  ! <FUNCTION NAME="mpp_pe">
  !  <OVERVIEW>
  !    Returns processor ID.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This returns the unique ID associated with a PE. This number runs
  !    between 0 and <TT>npes-1</TT>, where <TT>npes</TT> is the total
  !    processor count, returned by <TT>mpp_npes</TT>. For a uniprocessor
  !    application this will always return 0.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    mpp_pe()
  !  </TEMPLATE>
  ! </FUNCTION>
  function mpp_pe()
    integer :: mpp_pe

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_PE: You must first call mpp_init.' )
    mpp_pe = pe
    return
  end function mpp_pe

  !#####################################################################
  function mpp_node()
!calls mld_id from threadloc.c on sgi, which returns the hardware node ID from /hw/nodenum/...
    integer :: mpp_node
    integer :: mld_id

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_NODE: You must first call mpp_init.' )
    mpp_node = mld_id()
    return
  end function mpp_node

  !#####################################################################
  ! <FUNCTION NAME="mpp_npes">
  !  <OVERVIEW>
  !    Returns processor count for current pelist.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This returns the number of PEs in the current pelist. For a
  !    uniprocessor application, this will always return 1.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    mpp_npes()
  !  </TEMPLATE>
  ! </FUNCTION>
  function mpp_npes()
    integer :: mpp_npes

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_NPES: You must first call mpp_init.' )
    mpp_npes = size(peset(current_peset_num)%list(:))
    return
  end function mpp_npes

  !#####################################################################
  function mpp_root_pe()
    integer :: mpp_root_pe

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_ROOT_PE: You must first call mpp_init.' )
    mpp_root_pe = root_pe
    return
  end function mpp_root_pe

  !#####################################################################
  subroutine mpp_set_root_pe(num)
    integer, intent(in) :: num
    logical             :: opened

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SET_ROOT_PE: You must first call mpp_init.' )
    if( .NOT.(ANY(num.EQ.peset(current_peset_num)%list(:))) ) &
         call mpp_error( FATAL, 'MPP_SET_ROOT_PE: you cannot set a root PE outside the current pelist.' )
    !actions to take if root_pe has changed:
    ! open log_unit on new root_pe, close it on old root_pe and point its log_unit to stdout.
    !      if( num.NE.root_pe )then  !root_pe has changed
    !          if( pe.EQ.num )then
    !on the new root_pe
    !              if( log_unit.NE.out_unit )then
    !                  inquire( unit=log_unit, opened=opened )
    !                  if( .NOT.opened )open( unit=log_unit, status='OLD', file=trim(configfile), position='APPEND' )
    !              end if
    !          else if( pe.EQ.root_pe )then
    !on the old root_pe
    !              if( log_unit.NE.out_unit )then
    !                  inquire( unit=log_unit, opened=opened )
    !                  if( opened )close(log_unit)
    !                  log_unit = out_unit
    !              end if
    !          end if
    !      end if
    root_pe = num
    return
  end subroutine mpp_set_root_pe

  !#####################################################################
  ! <SUBROUTINE NAME="mpp_declare_pelist">
  !  <OVERVIEW>
  !    Declare a pelist.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This call is written specifically to accommodate a MPI restriction
  !    that requires a parent communicator to create a child communicator, In
  !    other words: a pelist cannot go off and declare a communicator, but
  !    every PE in the parent, including those not in pelist(:), must get
  !    together for the <TT>MPI_COMM_CREATE</TT> call. The parent is
  !    typically <TT>MPI_COMM_WORLD</TT>, though it could also be a subset
  !    that includes all PEs in <TT>pelist</TT>.
  !  
  !    The restriction does not apply to SMA but to have uniform code, you
  !    may as well call it.
  !  
  !    This call implies synchronization across the PEs in the current
  !    pelist, of which <TT>pelist</TT> is a subset.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call mpp_declare_pelist( pelist,name )
  !  </TEMPLATE>
  !  <IN NAME="pelist" DIM="(:)" TYPE="integer"></IN>
  ! </SUBROUTINE>

  subroutine mpp_declare_pelist( pelist, name )
    integer,                    intent(in) :: pelist(:)
    character(len=*), intent(in), optional :: name
    integer :: i

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DECLARE_PELIST: You must first call mpp_init.' )
    i = get_peset(pelist)
    write( peset(i)%name,'(a,i2.2)' ) 'PElist', i !default name
    if( PRESENT(name) )peset(i)%name = name
    return
  end subroutine mpp_declare_pelist

  !#####################################################################
  ! <SUBROUTINE NAME="mpp_set_current_pelist">
  !  <OVERVIEW>
  !    Set context pelist.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This call sets the value of the current pelist, which is the
  !    context for all subsequent "global" calls where the optional
  !    <TT>pelist</TT> argument is omitted. All the PEs that are to be in the
  !    current pelist must call it.
  !  
  !    In MPI, this call may hang unless <TT>pelist</TT> has been previous
  !    declared using <LINK
  !    SRC="#mpp_declare_pelist"><TT>mpp_declare_pelist</TT></LINK>.
  !  
  !    If the argument <TT>pelist</TT> is absent, the current pelist is
  !    set to the "world" pelist, of all PEs in the job.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_current_pelist( pelist )
  !  </TEMPLATE>
  !  <IN NAME="pliest" TYPE="integer"></IN>
  ! </SUBROUTINE>

  subroutine mpp_set_current_pelist( pelist, no_sync )
    !Once we branch off into a PE subset, we want subsequent "global" calls to
    !sync only across this subset. This is declared as the current pelist (peset(current_peset_num)%list)
    !when current_peset all pelist ops with no pelist should apply the current pelist.
    !also, we set the start PE in this pelist to be the root_pe.
    !unlike mpp_declare_pelist, this is called by the PEs in the pelist only
    !so if the PEset has not been previously declared, this will hang in MPI.
    !if pelist is omitted, we reset pelist to the world pelist.
    integer, intent(in), optional :: pelist(:)
    logical, intent(in), optional :: no_sync
    integer                       :: i

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SET_CURRENT_PELIST: You must first call mpp_init.' )
    if( PRESENT(pelist) )then
       if( .NOT.ANY(pe.EQ.pelist) )call mpp_error( FATAL, 'MPP_SET_CURRENT_PELIST: pe must be in pelist.' )
       current_peset_num = get_peset(pelist)
    else
       current_peset_num = world_peset_num
    end if
    call mpp_set_root_pe( MINVAL(peset(current_peset_num)%list(:)) )
    if(.not.PRESENT(no_sync))call mpp_sync()  !this is called to make sure everyone in the current pelist is here.
    !      npes = mpp_npes()
    return
  end subroutine mpp_set_current_pelist

  !#####################################################################
    !this is created for use by mpp_define_domains within a pelist
    !will be published but not publicized
  subroutine mpp_get_current_pelist( pelist, name, commID )
    integer, intent(out) :: pelist(:)
    character(len=*), intent(out), optional :: name
    integer, intent(out), optional :: commID

    if( size(pelist(:)).NE.size(peset(current_peset_num)%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_CURRENT_PELIST: size(pelist) is wrong.' )
    pelist(:) = peset(current_peset_num)%list(:)
    if( PRESENT(name) )name = peset(current_peset_num)%name
#ifdef use_libMPI
    if( PRESENT(commID) )commID = peset(current_peset_num)%id
#endif

    return
  end subroutine mpp_get_current_pelist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !                        PERFORMANCE PROFILING CALLS                          !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! <SUBROUTINE NAME="mpp_clock_set_grain">
  !  <OVERVIEW>
  !    Set the level of granularity of timing measurements.   
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This routine and three other routines, mpp_clock_id, mpp_clock_begin(id),
  !    and mpp_clock_end(id) may be used to time parallel code sections, and
  !    extract parallel statistics. Clocks are identified by names, which
  !    should be unique in the first 32 characters. The <TT>mpp_clock_id</TT>
  !    call initializes a clock of a given name and returns an integer
  !    <TT>id</TT>. This <TT>id</TT> can be used by subsequent
  !    <TT>mpp_clock_begin</TT> and <TT>mpp_clock_end</TT> calls set around a
  !    code section to be timed. Example:
  !    <PRE>
  !    integer :: id
  !    id = mpp_clock_id( 'Atmosphere' )
  !    call mpp_clock_begin(id)
  !    call atmos_model()
  !    call mpp_clock_end()
  !    </PRE>  
  !     Two flags may be used to alter the behaviour of
  !     <TT>mpp_clock</TT>. If the flag <TT>MPP_CLOCK_SYNC</TT> is turned on
  !     by <TT>mpp_clock_id</TT>, the clock calls <TT>mpp_sync</TT> across all
  !     the PEs in the current pelist at the top of the timed code section,
  !     but allows each PE to complete the code section (and reach
  !     <TT>mpp_clock_end</TT>) at different times. This allows us to measure
  !     load imbalance for a given code section. Statistics are written to
  !     <TT>stdout</TT> by <TT>mpp_exit</TT>.
  !     
  !     The flag <TT>MPP_CLOCK_DETAILED</TT> may be turned on by
  !     <TT>mpp_clock_id</TT> to get detailed communication
  !     profiles. Communication events of the types <TT>SEND, RECV, BROADCAST,
  !     REDUCE</TT> and <TT>WAIT</TT> are separately measured for data volume
  !     and time. Statistics are written to <TT>stdout</TT> by
  !     <TT>mpp_exit</TT>, and individual PE info is also written to the file
  !     <TT>mpp_clock.out.####</TT> where <TT>####</TT> is the PE id given by
  !     <TT>mpp_pe</TT>.
  !     
  !     The flags <TT>MPP_CLOCK_SYNC</TT> and <TT>MPP_CLOCK_DETAILED</TT> are
  !     integer parameters available by use association, and may be summed to
  !     turn them both on.
  !     
  !     While the nesting of clocks is allowed, please note that turning on
  !     the non-optional flags on inner clocks has certain subtle issues.
  !     Turning on <TT>MPP_CLOCK_SYNC</TT> on an inner
  !     clock may distort outer clock measurements of load imbalance. Turning
  !     on <TT>MPP_CLOCK_DETAILED</TT> will stop detailed measurements on its
  !     outer clock, since only one detailed clock may be active at one time.
  !     Also, detailed clocks only time a certain number of events per clock
  !     (currently 40000) to conserve memory. If this array overflows, a
  !     warning message is printed, and subsequent events for this clock are
  !     not timed.
  !     
  !     Timings are done using the <TT>f90</TT> standard
  !     <TT>SYSTEM_CLOCK</TT> intrinsic.
  !     
  !     The resolution of SYSTEM_CLOCK is often too coarse for use except
  !     across large swaths of code. On SGI systems this is transparently
  !     overloaded with a higher resolution clock made available in a
  !     non-portable fortran interface made available by
  !     <TT>nsclock.c</TT>. This approach will eventually be extended to other
  !     platforms.
  !     
  !     New behaviour added at the Havana release allows the user to embed
  !     profiling calls at varying levels of granularity all over the code,
  !     and for any particular run, set a threshold of granularity so that
  !     finer-grained clocks become dormant.
  !     
  !     The threshold granularity is held in the private module variable
  !     <TT>clock_grain</TT>. This value may be modified by the call
  !     <TT>mpp_clock_set_grain</TT>, and affect clocks initiated by
  !     subsequent calls to <TT>mpp_clock_id</TT>. The value of
  !     <TT>clock_grain</TT> is set to an arbitrarily large number initially.
  !     
  !     Clocks initialized by <TT>mpp_clock_id</TT> can set a new optional
  !     argument <TT>grain</TT> setting their granularity level. Clocks check
  !     this level against the current value of <TT>clock_grain</TT>, and are
  !     only triggered if they are <I>at or below ("coarser than")</I> the
  !     threshold. Finer-grained clocks are dormant for that run.
  !
  !The following grain levels are pre-defined:
  !
  !<pre>
  !!predefined clock granularities, but you can use any integer
  !!using CLOCK_LOOP and above may distort coarser-grain measurements
  !  integer, parameter, public :: CLOCK_COMPONENT=1 !component level, e.g model, exchange
  !  integer, parameter, public :: CLOCK_SUBCOMPONENT=11 !top level within a model component, e.g dynamics, physics
  !  integer, parameter, public :: CLOCK_MODULE=21 !module level, e.g main subroutine of a physics module
  !  integer, parameter, public :: CLOCK_ROUTINE=31 !level of individual subroutine or function
  !  integer, parameter, public :: CLOCK_LOOP=41 !loops or blocks within a routine
  !  integer, parameter, public :: CLOCK_INFRA=51 !infrastructure level, e.g halo update
  !</pre>
  !     
  !     Note that subsequent changes to <TT>clock_grain</TT> do not
  !     change the status of already initiated clocks, and that if the
  !     optional <TT>grain</TT> argument is absent, the clock is always
  !     triggered. This guarantees backward compatibility.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !     call mpp_clock_set_grain( grain )
  !  </TEMPLATE>
  !  <IN NAME="grain" TYPE="integer"></IN>
  ! </SUBROUTINE>

  subroutine mpp_clock_set_grain( grain )
    integer, intent(in) :: grain
    !set the granularity of times: only clocks whose grain is lower than
    !clock_grain are triggered, finer-grained clocks are dormant.
    !clock_grain is initialized to CLOCK_LOOP, so all clocks above the loop level
    !are triggered if this is never called.   
    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_SET_GRAIN: You must first call mpp_init.' )

    clock_grain = grain
    return
  end subroutine mpp_clock_set_grain

  !#####################################################################
  subroutine clock_init( id, name, flags, grain )
    integer,           intent(in) :: id
    character(len=*),  intent(in) :: name
    integer, intent(in), optional :: flags, grain
    integer                       :: i

    clocks(id)%name = name
    clocks(id)%tick = 0
    clocks(id)%total_ticks = 0
    clocks(id)%sync_on_begin = .FALSE.
    clocks(id)%detailed      = .FALSE.
    clocks(id)%peset_num = current_peset_num
    if( PRESENT(flags) )then
       if( BTEST(flags,0) )clocks(id)%sync_on_begin = .TRUE.
       if( BTEST(flags,1) )clocks(id)%detailed      = .TRUE.
    end if
    clocks(id)%grain = 0
    if( PRESENT(grain) )clocks(id)%grain = grain
    if( clocks(id)%detailed )then
       allocate( clocks(id)%events(MAX_EVENT_TYPES) )
       clocks(id)%events(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
       clocks(id)%events(EVENT_BROADCAST)%name = 'BROADCAST'
       clocks(id)%events(EVENT_RECV)%name = 'RECV'
       clocks(id)%events(EVENT_SEND)%name = 'SEND'
       clocks(id)%events(EVENT_WAIT)%name = 'WAIT'
       do i=1,MAX_EVENT_TYPES
          clocks(id)%events(i)%ticks(:) = 0
          clocks(id)%events(i)%bytes(:) = 0
          clocks(id)%events(i)%calls = 0
       end do
       clock_summary(id)%name = name
       clock_summary(id)%event(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
       clock_summary(id)%event(EVENT_BROADCAST)%name = 'BROADCAST'
       clock_summary(id)%event(EVENT_RECV)%name = 'RECV'
       clock_summary(id)%event(EVENT_SEND)%name = 'SEND'
       clock_summary(id)%event(EVENT_WAIT)%name = 'WAIT'
       do i=1,MAX_EVENT_TYPES
          clock_summary(id)%event(i)%msg_size_sums(:) = 0.0
          clock_summary(id)%event(i)%msg_time_sums(:) = 0.0
          clock_summary(id)%event(i)%total_data = 0.0
          clock_summary(id)%event(i)%total_time = 0.0
          clock_summary(id)%event(i)%msg_size_cnts(:) = 0
          clock_summary(id)%event(i)%total_cnts = 0
       end do
    end if
    return
  end subroutine clock_init

  !#####################################################################
    !return an ID for a new or existing clock
  function mpp_clock_id( name, flags, grain )
    integer                       :: mpp_clock_id
    character(len=*),  intent(in) :: name
    integer, intent(in), optional :: flags, grain
    integer                       :: i

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_ID: You must first call mpp_init.')

    !if grain is present, the clock is only triggered if it
    !is low ("coarse") enough: compared to clock_grain
    !finer-grained clocks are dormant.
    !if grain is absent, clock is triggered.
    if( PRESENT(grain) )then
       if( grain.GT.clock_grain )then
          mpp_clock_id = 0
          return
       end if
    end if
    mpp_clock_id = 1

    if( clock_num.EQ.0 )then  !first
       clock_num = mpp_clock_id
       call clock_init(mpp_clock_id,name,flags)
    else
       FIND_CLOCK: do while( trim(name).NE.trim(clocks(mpp_clock_id)%name) )
          mpp_clock_id = mpp_clock_id + 1
          if( mpp_clock_id.GT.clock_num )then
             if( mpp_clock_id.GT.MAX_CLOCKS )then
                call mpp_error( FATAL, 'MPP_CLOCK_ID: too many clock requests, ' // &
                      'check your clock id request or increase MAX_CLOCKS.')
             else               !new clock: initialize
                clock_num = mpp_clock_id
                call clock_init(mpp_clock_id,name,flags,grain)
                exit FIND_CLOCK 
             end if
          end if
       end do FIND_CLOCK
    endif
    return
  end function mpp_clock_id

  !#####################################################################
  subroutine mpp_clock_begin(id)
    integer, intent(in) :: id

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: You must first call mpp_init.' )
    if( .not. mpp_record_timing_data)return
    if( id.EQ.0 )return
    if( id.LT.0 .OR. id.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid id.' )

!$OMP MASTER
    if( clocks(id)%peset_num.NE.current_peset_num ) &
         call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: cannot change pelist context of a clock.' )
    if( clocks(id)%is_on) call mpp_error(FATAL, 'MPP_CLOCK_BEGIN: mpp_clock_begin is called again '// &
                'before calling mpp_clock_end for the clock '//trim(clocks(id)%name) )
    if( clocks(id)%sync_on_begin )then
       !do an untimed sync at the beginning of the clock
       !this puts all PEs in the current pelist on par, so that measurements begin together
       !ending time will be different, thus measuring load imbalance for this clock.
       call mpp_sync()
    end if

    num_clock_ids = num_clock_ids+1
    if(num_clock_ids > MAX_CLOCKS)call mpp_error(FATAL,'MPP_CLOCK_BEGIN: max num previous_clock exceeded.' )
    previous_clock(num_clock_ids) = current_clock
    current_clock = id
    call SYSTEM_CLOCK( clocks(id)%tick )
    clocks(id)%is_on = .true.
!$OMP END MASTER
    return
  end subroutine mpp_clock_begin

  !#####################################################################
  subroutine mpp_clock_end(id)
    integer, intent(in) :: id
    integer(LONG_KIND)  :: delta

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_END: You must first call mpp_init.' )
    if( .not. mpp_record_timing_data)return
    if( id.EQ.0 )return
    if( id.LT.0 .OR. id.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid id.' )
!$OMP MASTER
    if( .NOT. clocks(id)%is_on) call mpp_error(FATAL, 'MPP_CLOCK_END: mpp_clock_end is called '// &
                'before calling mpp_clock_begin for the clock '//trim(clocks(id)%name) )

    call SYSTEM_CLOCK(end_tick)
    if( clocks(id)%peset_num.NE.current_peset_num ) &
         call mpp_error( FATAL, 'MPP_CLOCK_END: cannot change pelist context of a clock.' )
    delta = end_tick - clocks(id)%tick
    if( delta.LT.0 )then
       write( stderr(),* )'pe, id, start_tick, end_tick, delta, max_ticks=', pe, id, clocks(id)%tick, end_tick, delta, max_ticks
       delta = delta + max_ticks + 1
       call mpp_error( WARNING, 'MPP_CLOCK_END: Clock rollover, assumed single roll.' )
    end if
    clocks(id)%total_ticks = clocks(id)%total_ticks + delta
    if(num_clock_ids < 1)then
       call mpp_error(FATAL,'MPP_CLOCK_END: min num previous_clock < 1.' )
    endif
    current_clock = previous_clock(num_clock_ids)
    num_clock_ids = num_clock_ids-1
    clocks(id)%is_on = .false.
!$OMP END MASTER
    return
  end subroutine mpp_clock_end

  !#####################################################################
  subroutine increment_current_clock( event_id, bytes )
    integer,           intent(in) :: event_id
    integer, intent(in), optional :: bytes
    integer                       :: n
    integer(LONG_KIND)            :: delta

    if( .not. mpp_record_timing_data )return
    if( current_clock.EQ.0 )return
    if( current_clock.LT.0 .OR. current_clock.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid current_clock.' )
    if( .NOT.clocks(current_clock)%detailed )return
    call SYSTEM_CLOCK(end_tick)
    n = clocks(current_clock)%events(event_id)%calls + 1

    if( n.EQ.MAX_EVENTS )call mpp_error( WARNING, &
         'MPP_CLOCK: events exceed MAX_EVENTS, ignore detailed profiling data for clock '//trim(clocks(current_clock)%name) )
    if( n.GT.MAX_EVENTS )return

    clocks(current_clock)%events(event_id)%calls = n
    delta = end_tick - start_tick
    if( delta.LT.0 )then
       write( stderr(),* )'pe, event_id, start_tick, end_tick, delta, max_ticks=', &
                           pe, event_id, start_tick, end_tick, delta, max_ticks
       delta = delta + max_ticks + 1
       call mpp_error( WARNING, 'MPP_CLOCK_END: Clock rollover, assumed single roll.' )
    end if
    clocks(current_clock)%events(event_id)%ticks(n) = delta
    if( PRESENT(bytes) )clocks(current_clock)%events(event_id)%bytes(n) = bytes
    return
  end subroutine increment_current_clock

  !#####################################################################

  subroutine dump_clock_summary()

    real              :: total_time,total_time_all,total_data
    real              :: msg_size,eff_BW,s
    integer           :: SD_UNIT, total_calls
    integer           :: i,j,k,ct, msg_cnt
    character(len=2)  :: u
    character(len=20) :: filename
    character(len=20),dimension(MAX_BINS),save :: bin

    data bin( 1)  /'  0   -    8    B:  '/
    data bin( 2)  /'  8   -   16    B:  '/
    data bin( 3)  /' 16   -   32    B:  '/
    data bin( 4)  /' 32   -   64    B:  '/
    data bin( 5)  /' 64   -  128    B:  '/
    data bin( 6)  /'128   -  256    B:  '/
    data bin( 7)  /'256   -  512    B:  '/
    data bin( 8)  /'512   - 1024    B:  '/
    data bin( 9)  /'  1.0 -    2.1 KB:  '/
    data bin(10)  /'  2.1 -    4.1 KB:  '/
    data bin(11)  /'  4.1 -    8.2 KB:  '/
    data bin(12)  /'  8.2 -   16.4 KB:  '/
    data bin(13)  /' 16.4 -   32.8 KB:  '/
    data bin(14)  /' 32.8 -   65.5 KB:  '/
    data bin(15)  /' 65.5 -  131.1 KB:  '/
    data bin(16)  /'131.1 -  262.1 KB:  '/
    data bin(17)  /'262.1 -  524.3 KB:  '/
    data bin(18)  /'524.3 - 1048.6 KB:  '/
    data bin(19)  /'  1.0 -    2.1 MB:  '/
    data bin(20)  /' >2.1          MB:  '/

    if( .NOT.ANY(clocks(1:clock_num)%detailed) )return
    write( filename,'(a,i6.6)' )'mpp_clock.out.', pe

    SD_UNIT = get_unit()
    open(SD_UNIT,file=trim(filename),form='formatted')

    COMM_TYPE: do ct = 1,clock_num

       if( .NOT.clocks(ct)%detailed )cycle
       write(SD_UNIT,*) &
            clock_summary(ct)%name(1:15),' Communication Data for PE ',pe

       write(SD_UNIT,*) ' '
       write(SD_UNIT,*) ' '

       total_time_all = 0.0
       EVENT_TYPE: do k = 1,MAX_EVENT_TYPES-1

          if(clock_summary(ct)%event(k)%total_time == 0.0)cycle

          total_time = clock_summary(ct)%event(k)%total_time
          total_time_all = total_time_all + total_time
          total_data = clock_summary(ct)%event(k)%total_data
          total_calls = clock_summary(ct)%event(k)%total_cnts

          write(SD_UNIT,1000) clock_summary(ct)%event(k)%name(1:9) // ':'

          write(SD_UNIT,1001) 'Total Data: ',total_data*1.0e-6, &
               'MB; Total Time: ', total_time, &
               'secs; Total Calls: ',total_calls

          write(SD_UNIT,*) ' '
          write(SD_UNIT,1002) '     Bin            Counts      Avg Size        Eff B/W'
          write(SD_UNIT,*) ' '

          BIN_LOOP: do j=1,MAX_BINS

             if(clock_summary(ct)%event(k)%msg_size_cnts(j)==0)cycle

             if(j<=8)then
                s = 1.0
                u = ' B'
             elseif(j<=18)then
                s = 1.0e-3
                u = 'KB'
             else
                s = 1.0e-6
                u = 'MB'
             endif

             msg_cnt = clock_summary(ct)%event(k)%msg_size_cnts(j)
             msg_size = &
                  s*(clock_summary(ct)%event(k)%msg_size_sums(j)/real(msg_cnt))
             eff_BW = (1.0e-6)*( clock_summary(ct)%event(k)%msg_size_sums(j) / &
                  clock_summary(ct)%event(k)%msg_time_sums(j) )

             write(SD_UNIT,1003) bin(j),msg_cnt,msg_size,u,eff_BW

          end do BIN_LOOP

          write(SD_UNIT,*) ' '
          write(SD_UNIT,*) ' '
       end do EVENT_TYPE

       ! "Data-less" WAIT

       if(clock_summary(ct)%event(MAX_EVENT_TYPES)%total_time>0.0)then

          total_time = clock_summary(ct)%event(MAX_EVENT_TYPES)%total_time
          total_time_all = total_time_all + total_time
          total_calls = clock_summary(ct)%event(MAX_EVENT_TYPES)%total_cnts

          write(SD_UNIT,1000) clock_summary(ct)%event(MAX_EVENT_TYPES)%name(1:9) // ':'

          write(SD_UNIT,1004) 'Total Calls: ',total_calls,'; Total Time: ', &
               total_time,'secs'

       endif

       write(SD_UNIT,*) ' '
       write(SD_UNIT,1005) 'Total communication time spent for ' // &
            clock_summary(ct)%name(1:9) // ': ',total_time_all,'secs'
       write(SD_UNIT,*) ' '
       write(SD_UNIT,*) ' '
       write(SD_UNIT,*) ' '

    end do COMM_TYPE

    close(SD_UNIT)

1000 format(a)
1001 format(a,f8.2,a,f8.2,a,i6)
1002 format(a)
1003 format(a,i6,'    ','  ',f6.1,a,'    ',f7.3,'MB/sec')
1004 format(a,i8,a,f9.2,a)
1005 format(a,f9.2,a)
    return
  end subroutine dump_clock_summary

  !#####################################################################

  integer function get_unit()

    integer,save :: i
    logical      :: l_open

!  9 is reserved for etc_unit
    do i=10,99
       inquire(unit=i,opened=l_open)
       if(.not.l_open)exit
    end do

    if(i==100)then
       call mpp_error(FATAL,'Unable to get I/O unit')
    else
       get_unit = i
    endif

    return
  end function get_unit

  !#####################################################################

  subroutine sum_clock_data()

    integer :: i,j,k,ct,event_size,event_cnt
    real    :: msg_time

    CLOCK_TYPE: do ct=1,clock_num
       if( .NOT.clocks(ct)%detailed )cycle
       EVENT_TYPE: do j=1,MAX_EVENT_TYPES-1
          event_cnt = clocks(ct)%events(j)%calls
          EVENT_SUMMARY: do i=1,event_cnt

             clock_summary(ct)%event(j)%total_cnts = &
                  clock_summary(ct)%event(j)%total_cnts + 1

             event_size = clocks(ct)%events(j)%bytes(i)

             k = find_bin(event_size)

             clock_summary(ct)%event(j)%msg_size_cnts(k) = &
                  clock_summary(ct)%event(j)%msg_size_cnts(k) + 1

             clock_summary(ct)%event(j)%msg_size_sums(k) = &
                  clock_summary(ct)%event(j)%msg_size_sums(k) &
                  + clocks(ct)%events(j)%bytes(i)

             clock_summary(ct)%event(j)%total_data = &
                  clock_summary(ct)%event(j)%total_data &
                  + clocks(ct)%events(j)%bytes(i)

             msg_time = clocks(ct)%events(j)%ticks(i)
             msg_time = tick_rate * real( clocks(ct)%events(j)%ticks(i) )

             clock_summary(ct)%event(j)%msg_time_sums(k) = &
                  clock_summary(ct)%event(j)%msg_time_sums(k) + msg_time

             clock_summary(ct)%event(j)%total_time = &
                  clock_summary(ct)%event(j)%total_time + msg_time

          end do EVENT_SUMMARY
       end do EVENT_TYPE

       j = MAX_EVENT_TYPES ! WAITs
       ! "msg_size_cnts" doesn't really mean anything for WAIT
       ! but position will be used to store number of counts for now.

       event_cnt = clocks(ct)%events(j)%calls
       clock_summary(ct)%event(j)%msg_size_cnts(1) = event_cnt
       clock_summary(ct)%event(j)%total_cnts       = event_cnt

       msg_time = tick_rate * real( sum ( clocks(ct)%events(j)%ticks(1:event_cnt) ) )
       clock_summary(ct)%event(j)%msg_time_sums(1) = &
            clock_summary(ct)%event(j)%msg_time_sums(1) + msg_time

       clock_summary(ct)%event(j)%total_time = clock_summary(ct)%event(j)%msg_time_sums(1)

    end do CLOCK_TYPE

    return
  contains
    integer function find_bin(event_size)

      integer,intent(in) :: event_size
      integer            :: k,msg_size

      msg_size = 8
      k = 1
      do while(event_size>msg_size .and. k<MAX_BINS)
         k = k+1
         msg_size = msg_size*2
      end do
      find_bin = k
      return
    end function find_bin

  end subroutine sum_clock_data

  !#####################################################################

  function uppercase (cs) 
    character(len=*), intent(in) :: cs
    character(len=len(cs)),target       :: uppercase 
    integer                      :: k,tlen
    character, pointer :: ca
    integer, parameter :: co=iachar('A')-iachar('a') ! case offset
    !The transfer function truncates the string with xlf90_r
    tlen = len_trim(cs)
    if(tlen <= 0) then      ! catch IBM compiler bug
       uppercase = cs  ! simply return input blank string
    else
    uppercase = cs(1:tlen)
! #etd
#if defined _CRAYX1
    do k=1, tlen
         if(uppercase(k:k) >= "a" .and. uppercase(k:k) <= 'z') uppercase(k:k) = achar(ichar(uppercase(k:k))+co)
       end do
#else
      do k=1, tlen
       ca => uppercase(k:k)
       if(ca >= "a" .and. ca <= "z") ca = achar(ichar(ca)+co)
    enddo
#endif
    endif
  end function uppercase

!#######################################################################

  function lowercase (cs) 
    character(len=*), intent(in) :: cs
    character(len=len(cs)),target       :: lowercase 
    integer, parameter :: co=iachar('a')-iachar('A') ! case offset
    integer                        :: k,tlen
    character, pointer :: ca
!  The transfer function truncates the string with xlf90_r
    tlen = len_trim(cs)
    if(tlen <= 0) then      ! catch IBM compiler bug
       lowercase = cs  ! simply return input blank string
    else
    lowercase = cs(1:tlen)
! #etd
#if defined _CRAYX1
       do k=1, tlen
         if(lowercase(k:k) >= "A" .and. lowercase(k:k) <= 'Z') lowercase(k:k) = achar(ichar(lowercase(k:k))+co)
       end do
#else
    do k=1, tlen
       ca => lowercase(k:k)
       if(ca >= "A" .and. ca <= "Z") ca = achar(ichar(ca)+co)
    enddo
#endif
    endif
  end function lowercase


  !#######################################################################

