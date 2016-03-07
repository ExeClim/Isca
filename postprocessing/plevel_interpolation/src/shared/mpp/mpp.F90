!-----------------------------------------------------------------------
!                 Communication for message-passing codes
!
! AUTHOR: V. Balaji (V.Balaji@noaa.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------
module mpp_mod
!a generalized communication package for use with shmem and MPI
!will add: co_array_fortran, MPI2
!Balaji (V.Balaji@noaa.gov) 11 May 1998

! <CONTACT EMAIL="V.Balaji@noaa.gov">
!   V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp.html"/>

! <OVERVIEW>
!   <TT>mpp_mod</TT>, is a set of simple calls to provide a uniform interface
!   to different message-passing libraries. It currently can be
!   implemented either in the SGI/Cray native SHMEM library or in the MPI
!   standard. Other libraries (e.g MPI-2, Co-Array Fortran) can be
!   incorporated as the need arises.
! </OVERVIEW>

! <DESCRIPTION>
!   The data transfer between a processor and its own memory is based
!   on <TT>load</TT> and <TT>store</TT> operations upon
!   memory. Shared-memory systems (including distributed shared memory
!   systems) have a single address space and any processor can acquire any
!   data within the memory by <TT>load</TT> and
!   <TT>store</TT>. The situation is different for distributed
!   parallel systems. Specialized MPP systems such as the T3E can simulate
!   shared-memory by direct data acquisition from remote memory. But if
!   the parallel code is distributed across a cluster, or across the Net,
!   messages must be sent and received using the protocols for
!   long-distance communication, such as TCP/IP. This requires a
!   ``handshaking'' between nodes of the distributed system. One can think
!   of the two different methods as involving <TT>put</TT>s or
!   <TT>get</TT>s (e.g the SHMEM library), or in the case of
!   negotiated communication (e.g MPI), <TT>send</TT>s and
!   <TT>recv</TT>s.
!   
!   The difference between SHMEM and MPI is that SHMEM uses one-sided
!   communication, which can have very low-latency high-bandwidth
!   implementations on tightly coupled systems. MPI is a standard
!   developed for distributed computing across loosely-coupled systems,
!   and therefore incurs a software penalty for negotiating the
!   communication. It is however an open industry standard whereas SHMEM
!   is a proprietary interface. Besides, the <TT>put</TT>s or
!   <TT>get</TT>s on which it is based cannot currently be implemented in
!   a cluster environment (there are recent announcements from Compaq that
!   occasion hope).
!   
!   The message-passing requirements of climate and weather codes can be
!   reduced to a fairly simple minimal set, which is easily implemented in
!   any message-passing API. <TT>mpp_mod</TT> provides this API.
!
!    Features of <TT>mpp_mod</TT> include:
!   
!    1) Simple, minimal API, with free access to underlying API for
!       more complicated stuff.<BR/>
!    2) Design toward typical use in climate/weather CFD codes.<BR/>
!    3) Performance to be not significantly lower than any native API.
!   
!   This module is used to develop higher-level calls for <LINK 
!   SRC="mpp_domains.html">domain decomposition</LINK> and <LINK
!   SRC="mpp_io.html">parallel I/O</LINK>.
!   
!   Parallel computing is initially daunting, but it soon becomes
!   second nature, much the way many of us can now write vector code
!   without much effort. The key insight required while reading and
!   writing parallel code is in arriving at a mental grasp of several
!   independent parallel execution streams through the same code (the SPMD
!   model). Each variable you examine may have different values for each
!   stream, the processor ID being an obvious example. Subroutines and
!   function calls are particularly subtle, since it is not always obvious
!   from looking at a call what synchronization between execution streams
!   it implies. An example of erroneous code would be a global barrier
!   call (see <LINK SRC="#mpp_sync">mpp_sync</LINK> below) placed
!   within a code block that not all PEs will execute, e.g:
!   
!   <PRE>
!   if( pe.EQ.0 )call mpp_sync()
!   </PRE>
!   
!   Here only PE 0 reaches the barrier, where it will wait
!   indefinitely. While this is a particularly egregious example to
!   illustrate the coding flaw, more subtle versions of the same are
!   among the most common errors in parallel code.
!   
!   It is therefore important to be conscious of the context of a
!   subroutine or function call, and the implied synchronization. There
!   are certain calls here (e.g <TT>mpp_declare_pelist, mpp_init,
!   mpp_malloc, mpp_set_stack_size</TT>) which must be called by all
!   PEs. There are others which must be called by a subset of PEs (here
!   called a <TT>pelist</TT>) which must be called by all the PEs in the
!   <TT>pelist</TT> (e.g <TT>mpp_max, mpp_sum, mpp_sync</TT>). Still
!   others imply no synchronization at all. I will make every effort to
!   highlight the context of each call in the MPP modules, so that the
!   implicit synchronization is spelt out.  
!   
!   For performance it is necessary to keep synchronization as limited
!   as the algorithm being implemented will allow. For instance, a single
!   message between two PEs should only imply synchronization across the
!   PEs in question. A <I>global</I> synchronization (or <I>barrier</I>)
!   is likely to be slow, and is best avoided. But codes first
!   parallelized on a Cray T3E tend to have many global syncs, as very
!   fast barriers were implemented there in hardware.
!   
!   Another reason to use pelists is to run a single program in MPMD
!   mode, where different PE subsets work on different portions of the
!   code. A typical example is to assign an ocean model and atmosphere
!   model to different PE subsets, and couple them concurrently instead of
!   running them serially. The MPP module provides the notion of a
!   <I>current pelist</I>, which is set when a group of PEs branch off
!   into a subset. Subsequent calls that omit the <TT>pelist</TT> optional
!   argument (seen below in many of the individual calls) assume that the
!   implied synchronization is across the current pelist. The calls
!   <TT>mpp_root_pe</TT> and <TT>mpp_npes</TT> also return the values
!   appropriate to the current pelist. The <TT>mpp_set_current_pelist</TT>
!   call is provided to set the current pelist.

! </DESCRIPTION>
! <PUBLIC>
!  F90 is a strictly-typed language, and the syntax pass of the
!  compiler requires matching of type, kind and rank (TKR). Most calls
!  listed here use a generic type, shown here as <TT>MPP_TYPE_</TT>. This
!  is resolved in the pre-processor stage to any of a variety of
!  types. In general the MPP operations work on 4-byte and 8-byte
!  variants of <TT>integer, real, complex, logical</TT> variables, of
!  rank 0 to 5, leading to 48 specific module procedures under the same
!  generic interface. Any of the variables below shown as
!  <TT>MPP_TYPE_</TT> is treated in this way.
! </PUBLIC>

#include <fms_platform.h>

#if defined(use_libSMA) && defined(sgi_mipspro)
  use shmem_interface
#endif

#if defined(use_libMPI) && defined(sgi_mipspro)
  use mpi
#endif

  use mpp_parameter_mod, only : MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE
  use mpp_parameter_mod, only : NOTE, WARNING, FATAL, MPP_CLOCK_DETAILED,MPP_CLOCK_SYNC
  use mpp_parameter_mod, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER
  use mpp_parameter_mod, only : CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
  use mpp_parameter_mod, only : MAX_EVENTS, MAX_BINS, MAX_EVENT_TYPES, PESET_MAX, MAX_CLOCKS
  use mpp_parameter_mod, only : MAXPES, EVENT_WAIT, EVENT_ALLREDUCE, EVENT_BROADCAST
  use mpp_parameter_mod, only : EVENT_RECV, EVENT_SEND, MPP_READY, MPP_WAIT
  use mpp_parameter_mod, only : mpp_parameter_version=>version, mpp_parameter_tagname=>tagname
  use mpp_data_mod,      only : stat, mpp_stack, ptr_stack, status, ptr_status, sync, ptr_sync  
  use mpp_data_mod,      only : mpp_from_pe, ptr_from, remote_data_loc, ptr_remote
  use mpp_data_mod,      only : mpp_data_version=>version, mpp_data_tagname=>tagname

implicit none
private

#if defined(use_libSMA) 
#include <mpp/shmem.fh>
#endif

#if defined(use_libMPI) && !defined(sgi_mipspro)
#include <mpif.h>   
!sgi_mipspro gets this from 'use mpi'
#endif

  !--- public paramters  -----------------------------------------------
  public :: MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE, NOTE, WARNING, FATAL
  public :: MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
  public :: CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
  public :: MAXPES, EVENT_RECV, EVENT_SEND

  !--- public data from mpp_data_mod ------------------------------
  public :: request

  !--- public interface from mpp_util.h ------------------------------
  public :: stdin, stdout, stderr, stdlog, lowercase, uppercase, mpp_error, mpp_error_state
  public :: mpp_set_warn_level, mpp_sync, mpp_sync_self, mpp_set_stack_size, mpp_pe
  public :: mpp_node, mpp_npes, mpp_root_pe, mpp_set_root_pe, mpp_declare_pelist
  public :: mpp_get_current_pelist, mpp_set_current_pelist, mpp_clock_begin, mpp_clock_end
  public :: mpp_clock_id, mpp_clock_set_grain, mpp_record_timing_data, get_unit

  !--- public interface from mpp_comm.h ------------------------------
  public :: mpp_chksum, mpp_max, mpp_min, mpp_sum, mpp_transmit, mpp_send, mpp_recv
  public :: mpp_broadcast, mpp_malloc, mpp_init, mpp_exit
#ifdef use_MPI_GSM
  public :: mpp_gsm_malloc, mpp_gsm_free
#endif

  !*********************************************************************
  !
  !    public data type
  !
  !*********************************************************************
  !peset hold communicators as SHMEM-compatible triads (start, log2(stride), num)
  type :: communicator
     private
     character(len=32) :: name
     integer, pointer  :: list(:) =>NULL()
     integer           :: count
     integer           :: start, log2stride ! dummy variables when libMPI is defined.
     integer           :: id, group         ! MPI communicator and group id for this PE set.
                                            ! dummy variables when libSMA is defined.
  end type communicator

  type :: event
     private
     character(len=16)                         :: name
     integer(LONG_KIND), dimension(MAX_EVENTS) :: ticks, bytes
     integer                                   :: calls
  end type event

  !a clock contains an array of event profiles for a region
  type :: clock
     private
     character(len=32)    :: name
     integer(LONG_KIND)   :: tick
     integer(LONG_KIND)   :: total_ticks
     integer              :: peset_num
     logical              :: sync_on_begin, detailed
     integer              :: grain
     type(event), pointer :: events(:) =>NULL() !if needed, allocate to MAX_EVENT_TYPES
     logical              :: is_on              !initialize to false. set true when calling mpp_clock_begin
                                                ! set false when calling mpp_clock_end
  end type clock

  type :: Clock_Data_Summary
     private
     character(len=16)  :: name
     real(DOUBLE_KIND)  :: msg_size_sums(MAX_BINS)
     real(DOUBLE_KIND)  :: msg_time_sums(MAX_BINS)
     real(DOUBLE_KIND)  :: total_data
     real(DOUBLE_KIND)  :: total_time
     integer(LONG_KIND) :: msg_size_cnts(MAX_BINS)
     integer(LONG_KIND) :: total_cnts
  end type Clock_Data_Summary

  type :: Summary_Struct
     private
     character(len=16)         :: name
     type (Clock_Data_Summary) :: event(MAX_EVENT_TYPES)
  end type Summary_Struct

!***********************************************************************
!
!     public interface from mpp_util.h
!
!***********************************************************************
  ! <INTERFACE NAME="mpp_error">
  !  <OVERVIEW>
  !    Error handler.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    It is strongly recommended that all error exits pass through
  !    <TT>mpp_error</TT> to assure the program fails cleanly. An individual
  !    PE encountering a <TT>STOP</TT> statement, for instance, can cause the
  !    program to hang. The use of the <TT>STOP</TT> statement is strongly
  !    discouraged.
  !    
  !    Calling mpp_error with no arguments produces an immediate error
  !    exit, i.e:
  !    <PRE>
  !    call mpp_error
  !    call mpp_error(FATAL)
  !    </PRE>
  !    are equivalent.
  !    
  !    The argument order
  !    <PRE>
  !    call mpp_error( routine, errormsg, errortype )
  !    </PRE>
  !    is also provided to support legacy code. In this version of the
  !    call, none of the arguments may be omitted.
  !    
  !    The behaviour of <TT>mpp_error</TT> for a <TT>WARNING</TT> can be
  !    controlled with an additional call <TT>mpp_set_warn_level</TT>.
  !    <PRE>
  !    call mpp_set_warn_level(ERROR)
  !    </PRE>
  !    causes <TT>mpp_error</TT> to treat <TT>WARNING</TT>
  !    exactly like <TT>FATAL</TT>.
  !    <PRE>
  !    call mpp_set_warn_level(WARNING)
  !    </PRE>
  !    resets to the default behaviour described above.
  !    
  !    <TT>mpp_error</TT> also has an internal error state which
  !    maintains knowledge of whether a warning has been issued. This can be
  !    used at startup in a subroutine that checks if the model has been
  !    properly configured. You can generate a series of warnings using
  !    <TT>mpp_error</TT>, and then check at the end if any warnings has been
  !    issued using the function <TT>mpp_error_state()</TT>. If the value of
  !    this is <TT>WARNING</TT>, at least one warning has been issued, and
  !    the user can take appropriate action:
  !    
  !    <PRE>
  !    if( ... )call mpp_error( WARNING, '...' )
  !    if( ... )call mpp_error( WARNING, '...' )
  !    if( ... )call mpp_error( WARNING, '...' )
  !    ...
  !    if( mpp_error_state().EQ.WARNING )call mpp_error( FATAL, '...' )
  !    </PRE>
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_error( errortype, routine, errormsg )
  !  </TEMPLATE>
  !  <IN NAME="errortype">
  !    One of <TT>NOTE</TT>, <TT>WARNING</TT> or <TT>FATAL</TT> 
  !    (these definitions are acquired by use association).
  !    <TT>NOTE</TT> writes <TT>errormsg</TT> to <TT>STDOUT</TT>. 
  !    <TT>WARNING</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>.
  !    <TT>FATAL</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>,
  !    and induces a clean error exit with a call stack traceback.
  !  </IN>
  ! </INTERFACE>
  interface mpp_error
     module procedure mpp_error_basic
     module procedure mpp_error_mesg
     module procedure mpp_error_noargs
     module procedure mpp_error_is
     module procedure mpp_error_rs
     module procedure mpp_error_ia
     module procedure mpp_error_ra
     module procedure mpp_error_ia_ia
     module procedure mpp_error_ia_ra
     module procedure mpp_error_ra_ia
     module procedure mpp_error_ra_ra
     module procedure mpp_error_ia_is
     module procedure mpp_error_ia_rs
     module procedure mpp_error_ra_is
     module procedure mpp_error_ra_rs
     module procedure mpp_error_is_ia
     module procedure mpp_error_is_ra
     module procedure mpp_error_rs_ia
     module procedure mpp_error_rs_ra
     module procedure mpp_error_is_is
     module procedure mpp_error_is_rs
     module procedure mpp_error_rs_is
     module procedure mpp_error_rs_rs
  end interface

  interface array_to_char
     module procedure iarray_to_char
     module procedure rarray_to_char
  end interface

!***********************************************************************
!
!    public interface from mpp_comm.h
!
!***********************************************************************
#ifdef use_libSMA
  !currently SMA contains no generic shmem_wait for different integer kinds:
  !I have inserted one here
  interface shmem_integer_wait
     module procedure shmem_int4_wait_local
     module procedure shmem_int8_wait_local
  end interface
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! <SUBROUTINE NAME="mpp_init">
  !  <OVERVIEW>
  !   Initialize <TT>mpp_mod</TT>.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Called to initialize the <TT>mpp_mod</TT> package. It is recommended
  !   that this call be the first executed line in your program. It sets the
  !   number of PEs assigned to this run (acquired from the command line, or
  !   through the environment variable <TT>NPES</TT>), and associates an ID
  !   number to each PE. These can be accessed by calling <LINK
  !   SRC="#mpp_npes"><TT>mpp_npes</TT></LINK> and <LINK
  !   SRC="#mpp_pe"><TT>mpp_pe</TT></LINK>.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call mpp_init( flags )
  !  </TEMPLATE>
  !  <IN NAME="flags" TYPE="integer">
  !   <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to
  !   have <TT>mpp_mod</TT> keep you informed of what it's up to.
  !  </IN>
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="mpp_exit">
  !  <OVERVIEW>
  !   Exit <TT>mpp_mod</TT>.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Called at the end of the run, or to re-initialize <TT>mpp_mod</TT>,
  !   should you require that for some odd reason.
  !
  !   This call implies synchronization across all PEs.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call mpp_exit()
  !  </TEMPLATE>
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="mpp_malloc">
  !  <OVERVIEW>
  !    Symmetric memory allocation.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This routine is used on SGI systems when <TT>mpp_mod</TT> is
  !    invoked in the SHMEM library. It ensures that dynamically allocated
  !    memory can be used with <TT>shmem_get</TT> and
  !    <TT>shmem_put</TT>. This is called <I>symmetric
  !    allocation</I> and is described in the
  !    <TT>intro_shmem</TT> man page. <TT>ptr</TT> is a <I>Cray
  !    pointer</I> (see the section on <LINK
  !    SRC="#PORTABILITY">portability</LINK>).  The operation can be expensive
  !    (since it requires a global barrier). We therefore attempt to re-use
  !    existing allocation whenever possible. Therefore <TT>len</TT>
  !    and <TT>ptr</TT> must have the <TT>SAVE</TT> attribute
  !    in the calling routine, and retain the information about the last call
  !    to <TT>mpp_malloc</TT>. Additional memory is symmetrically
  !    allocated if and only if <TT>newlen</TT> exceeds
  !    <TT>len</TT>.
  !
  !    This is never required on Cray PVP or MPP systems. While the T3E
  !    manpages do talk about symmetric allocation, <TT>mpp_mod</TT>
  !    is coded to remove this restriction.
  !
  !    It is never required if <TT>mpp_mod</TT> is invoked in MPI.
  !
  !   This call implies synchronization across all PEs.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call mpp_malloc( ptr, newlen, len )
  !  </TEMPLATE>
  !  <IN NAME="ptr">
  !     a cray pointer, points to a dummy argument in this routine.
  !  </IN>
  !  <IN NAME="newlen" TYPE="integer">
  !     the required allocation length for the pointer ptr
  !  </IN>
  !  <IN NAME="len" TYPE="integer">
  !     the current allocation (0 if unallocated).
  !  </IN>
  ! </SUBROUTINE>

  !#####################################################################

  ! <SUBROUTINE NAME="mpp_set_stack_size">
  !  <OVERVIEW>
  !    Allocate module internal workspace.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    <TT>mpp_mod</TT> maintains a private internal array called
  !    <TT>mpp_stack</TT> for private workspace. This call sets the length,
  !    in words, of this array. 
  !
  !    The <TT>mpp_init</TT> call sets this
  !    workspace length to a default of 32768, and this call may be used if a
  !    longer workspace is needed.
  !    
  !    This call implies synchronization across all PEs.
  !    
  !    This workspace is symmetrically allocated, as required for
  !    efficient communication on SGI and Cray MPP systems. Since symmetric
  !    allocation must be performed by <I>all</I> PEs in a job, this call
  !    must also be called by all PEs, using the same value of
  !    <TT>n</TT>. Calling <TT>mpp_set_stack_size</TT> from a subset of PEs,
  !    or with unequal argument <TT>n</TT>, may cause the program to hang.
  !    
  !    If any MPP call using <TT>mpp_stack</TT> overflows the declared
  !    stack array, the program will abort with a message specifying the
  !    stack length that is required. Many users wonder why, if the required
  !    stack length can be computed, it cannot also be specified at that
  !    point. This cannot be automated because there is no way for the
  !    program to know if all PEs are present at that call, and with equal
  !    values of <TT>n</TT>. The program must be rerun by the user with the
  !    correct argument to <TT>mpp_set_stack_size</TT>, called at an
  !    appropriate point in the code where all PEs are known to be present.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_stack_size(n)
  !  </TEMPLATE>
  !  <IN NAME="n" TYPE="integer"></IN>
  ! </SUBROUTINE>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! <INTERFACE NAME="mpp_max">
  !  <OVERVIEW>
  !    Reduction operations.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Find the max of scalar a the PEs in pelist
  !    result is also automatically broadcast to all PEs
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call  mpp_max( a, pelist )
  !  </TEMPLATE>
  !  <IN NAME="a">
  !    <TT>real</TT> or <TT>integer</TT>, of 4-byte of 8-byte kind.
  !  </IN>
  !  <IN NAME="pelist">
  !    If <TT>pelist</TT> is omitted, the context is assumed to be the
  !    current pelist. This call implies synchronization across the PEs in
  !    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !  </IN>
  ! </INTERFACE>

  interface mpp_max
     module procedure mpp_max_real8
#ifndef no_8byte_integers
     module procedure mpp_max_int8
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_max_real4
#endif
     module procedure mpp_max_int4
  end interface

  interface mpp_min
     module procedure mpp_min_real8
#ifndef no_8byte_integers
     module procedure mpp_min_int8
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_min_real4
#endif
     module procedure mpp_min_int4
  end interface


  ! <INTERFACE NAME="mpp_sum">
  !  <OVERVIEW>
  !    Reduction operation.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !    <TT>integer, real, complex</TT> variables, of rank 0 or 1. A
  !    contiguous block from a multi-dimensional array may be passed by its
  !    starting address and its length, as in <TT>f77</TT>.
  !
  !    Library reduction operators are not required or guaranteed to be
  !    bit-reproducible. In any case, changing the processor count changes
  !    the data layout, and thus very likely the order of operations. For
  !    bit-reproducible sums of distributed arrays, consider using the
  !    <TT>mpp_global_sum</TT> routine provided by the <LINK
  !    SRC="mpp_domains.html"><TT>mpp_domains</TT></LINK> module.
  !
  !    The <TT>bit_reproducible</TT> flag provided in earlier versions of
  !    this routine has been removed.
  !
  !
  !    If <TT>pelist</TT> is omitted, the context is assumed to be the
  !    current pelist. This call implies synchronization across the PEs in
  !    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_sum( a, length, pelist )
  !  </TEMPLATE>
  !  <IN NAME="length"></IN>
  !  <IN NAME="pelist"></IN>
  !  <INOUT NAME="a"></INOUT>
  ! </INTERFACE>

  interface mpp_sum
#ifndef no_8byte_integers
     module procedure mpp_sum_int8
     module procedure mpp_sum_int8_scalar
     module procedure mpp_sum_int8_2d
     module procedure mpp_sum_int8_3d
     module procedure mpp_sum_int8_4d
     module procedure mpp_sum_int8_5d
#endif
     module procedure mpp_sum_real8
     module procedure mpp_sum_real8_scalar
     module procedure mpp_sum_real8_2d
     module procedure mpp_sum_real8_3d
     module procedure mpp_sum_real8_4d
     module procedure mpp_sum_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_sum_cmplx8
     module procedure mpp_sum_cmplx8_scalar
     module procedure mpp_sum_cmplx8_2d
     module procedure mpp_sum_cmplx8_3d
     module procedure mpp_sum_cmplx8_4d
     module procedure mpp_sum_cmplx8_5d
#endif
     module procedure mpp_sum_int4
     module procedure mpp_sum_int4_scalar
     module procedure mpp_sum_int4_2d
     module procedure mpp_sum_int4_3d
     module procedure mpp_sum_int4_4d
     module procedure mpp_sum_int4_5d
#ifdef OVERLOAD_R4
     module procedure mpp_sum_real4
     module procedure mpp_sum_real4_scalar
     module procedure mpp_sum_real4_2d
     module procedure mpp_sum_real4_3d
     module procedure mpp_sum_real4_4d
     module procedure mpp_sum_real4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_sum_cmplx4
     module procedure mpp_sum_cmplx4_scalar
     module procedure mpp_sum_cmplx4_2d
     module procedure mpp_sum_cmplx4_3d
     module procedure mpp_sum_cmplx4_4d
     module procedure mpp_sum_cmplx4_5d
#endif
  end interface

  !#####################################################################

  ! <INTERFACE NAME="mpp_transmit">
  !  <OVERVIEW>
  !    Basic message-passing call.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !    <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
  !    contiguous block from a multi-dimensional array may be passed by its
  !    starting address and its length, as in <TT>f77</TT>.
  !    
  !    <TT>mpp_transmit</TT> is currently implemented as asynchronous
  !    outward transmission and synchronous inward transmission. This follows
  !    the behaviour of <TT>shmem_put</TT> and <TT>shmem_get</TT>. In MPI, it
  !    is implemented as <TT>mpi_isend</TT> and <TT>mpi_recv</TT>. For most
  !    applications, transmissions occur in pairs, and are here accomplished
  !    in a single call.
  !    
  !    The special PE designations <TT>NULL_PE</TT>,
  !    <TT>ANY_PE</TT> and <TT>ALL_PES</TT> are provided by use
  !    association.
  !    
  !    <TT>NULL_PE</TT>: is used to disable one of the pair of
  !    transmissions.<BR/>
  !    <TT>ANY_PE</TT>: is used for unspecific remote
  !    destination. (Please note that <TT>put_pe=ANY_PE</TT> has no meaning
  !    in the MPI context, though it is available in the SHMEM invocation. If
  !    portability is a concern, it is best avoided).<BR/>
  !    <TT>ALL_PES</TT>: is used for broadcast operations.
  !    
  !    It is recommended that <LINK
  !    SRC="#mpp_broadcast"><TT>mpp_broadcast</TT></LINK> be used for
  !    broadcasts.
  !    
  !    The following example illustrates the use of
  !    <TT>NULL_PE</TT> and <TT>ALL_PES</TT>:
  !    
  !    <PRE>
  !    real, dimension(n) :: a
  !    if( pe.EQ.0 )then
  !        do p = 1,npes-1
  !           call mpp_transmit( a, n, p, a, n, NULL_PE )
  !        end do
  !    else
  !        call mpp_transmit( a, n, NULL_PE, a, n, 0 )
  !    end if
  !    
  !    call mpp_transmit( a, n, ALL_PES, a, n, 0 )
  !    </PRE>
  !    
  !    The do loop and the broadcast operation above are equivalent.
  !    
  !    Two overloaded calls <TT>mpp_send</TT> and
  !     <TT>mpp_recv</TT> have also been
  !    provided. <TT>mpp_send</TT> calls <TT>mpp_transmit</TT>
  !    with <TT>get_pe=NULL_PE</TT>. <TT>mpp_recv</TT> calls
  !    <TT>mpp_transmit</TT> with <TT>put_pe=NULL_PE</TT>. Thus
  !    the do loop above could be written more succinctly:
  !    
  !    <PRE>
  !    if( pe.EQ.0 )then
  !        do p = 1,npes-1
  !           call mpp_send( a, n, p )
  !        end do
  !    else
  !        call mpp_recv( a, n, 0 )
  !    end if
  !    </PRE>
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_transmit( put_data, put_len, put_pe, get_data, get_len, get_pe )
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_transmit
     module procedure mpp_transmit_real8
     module procedure mpp_transmit_real8_scalar
     module procedure mpp_transmit_real8_2d
     module procedure mpp_transmit_real8_3d
     module procedure mpp_transmit_real8_4d
     module procedure mpp_transmit_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_transmit_cmplx8
     module procedure mpp_transmit_cmplx8_scalar
     module procedure mpp_transmit_cmplx8_2d
     module procedure mpp_transmit_cmplx8_3d
     module procedure mpp_transmit_cmplx8_4d
     module procedure mpp_transmit_cmplx8_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_transmit_int8
     module procedure mpp_transmit_int8_scalar
     module procedure mpp_transmit_int8_2d
     module procedure mpp_transmit_int8_3d
     module procedure mpp_transmit_int8_4d
     module procedure mpp_transmit_int8_5d
     module procedure mpp_transmit_logical8
     module procedure mpp_transmit_logical8_scalar
     module procedure mpp_transmit_logical8_2d
     module procedure mpp_transmit_logical8_3d
     module procedure mpp_transmit_logical8_4d
     module procedure mpp_transmit_logical8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_transmit_real4
     module procedure mpp_transmit_real4_scalar
     module procedure mpp_transmit_real4_2d
     module procedure mpp_transmit_real4_3d
     module procedure mpp_transmit_real4_4d
     module procedure mpp_transmit_real4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_transmit_cmplx4
     module procedure mpp_transmit_cmplx4_scalar
     module procedure mpp_transmit_cmplx4_2d
     module procedure mpp_transmit_cmplx4_3d
     module procedure mpp_transmit_cmplx4_4d
     module procedure mpp_transmit_cmplx4_5d
#endif
     module procedure mpp_transmit_int4
     module procedure mpp_transmit_int4_scalar
     module procedure mpp_transmit_int4_2d
     module procedure mpp_transmit_int4_3d
     module procedure mpp_transmit_int4_4d
     module procedure mpp_transmit_int4_5d
     module procedure mpp_transmit_logical4
     module procedure mpp_transmit_logical4_scalar
     module procedure mpp_transmit_logical4_2d
     module procedure mpp_transmit_logical4_3d
     module procedure mpp_transmit_logical4_4d
     module procedure mpp_transmit_logical4_5d
  end interface
  interface mpp_recv
     module procedure mpp_recv_real8
     module procedure mpp_recv_real8_scalar
     module procedure mpp_recv_real8_2d
     module procedure mpp_recv_real8_3d
     module procedure mpp_recv_real8_4d
     module procedure mpp_recv_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_recv_cmplx8
     module procedure mpp_recv_cmplx8_scalar
     module procedure mpp_recv_cmplx8_2d
     module procedure mpp_recv_cmplx8_3d
     module procedure mpp_recv_cmplx8_4d
     module procedure mpp_recv_cmplx8_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_recv_int8
     module procedure mpp_recv_int8_scalar
     module procedure mpp_recv_int8_2d
     module procedure mpp_recv_int8_3d
     module procedure mpp_recv_int8_4d
     module procedure mpp_recv_int8_5d
     module procedure mpp_recv_logical8
     module procedure mpp_recv_logical8_scalar
     module procedure mpp_recv_logical8_2d
     module procedure mpp_recv_logical8_3d
     module procedure mpp_recv_logical8_4d
     module procedure mpp_recv_logical8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_recv_real4
     module procedure mpp_recv_real4_scalar
     module procedure mpp_recv_real4_2d
     module procedure mpp_recv_real4_3d
     module procedure mpp_recv_real4_4d
     module procedure mpp_recv_real4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_recv_cmplx4
     module procedure mpp_recv_cmplx4_scalar
     module procedure mpp_recv_cmplx4_2d
     module procedure mpp_recv_cmplx4_3d
     module procedure mpp_recv_cmplx4_4d
     module procedure mpp_recv_cmplx4_5d
#endif
     module procedure mpp_recv_int4
     module procedure mpp_recv_int4_scalar
     module procedure mpp_recv_int4_2d
     module procedure mpp_recv_int4_3d
     module procedure mpp_recv_int4_4d
     module procedure mpp_recv_int4_5d
     module procedure mpp_recv_logical4
     module procedure mpp_recv_logical4_scalar
     module procedure mpp_recv_logical4_2d
     module procedure mpp_recv_logical4_3d
     module procedure mpp_recv_logical4_4d
     module procedure mpp_recv_logical4_5d
  end interface
  interface mpp_send
     module procedure mpp_send_real8
     module procedure mpp_send_real8_scalar
     module procedure mpp_send_real8_2d
     module procedure mpp_send_real8_3d
     module procedure mpp_send_real8_4d
     module procedure mpp_send_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_send_cmplx8
     module procedure mpp_send_cmplx8_scalar
     module procedure mpp_send_cmplx8_2d
     module procedure mpp_send_cmplx8_3d
     module procedure mpp_send_cmplx8_4d
     module procedure mpp_send_cmplx8_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_send_int8
     module procedure mpp_send_int8_scalar
     module procedure mpp_send_int8_2d
     module procedure mpp_send_int8_3d
     module procedure mpp_send_int8_4d
     module procedure mpp_send_int8_5d
     module procedure mpp_send_logical8
     module procedure mpp_send_logical8_scalar
     module procedure mpp_send_logical8_2d
     module procedure mpp_send_logical8_3d
     module procedure mpp_send_logical8_4d
     module procedure mpp_send_logical8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_send_real4
     module procedure mpp_send_real4_scalar
     module procedure mpp_send_real4_2d
     module procedure mpp_send_real4_3d
     module procedure mpp_send_real4_4d
     module procedure mpp_send_real4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_send_cmplx4
     module procedure mpp_send_cmplx4_scalar
     module procedure mpp_send_cmplx4_2d
     module procedure mpp_send_cmplx4_3d
     module procedure mpp_send_cmplx4_4d
     module procedure mpp_send_cmplx4_5d
#endif
     module procedure mpp_send_int4
     module procedure mpp_send_int4_scalar
     module procedure mpp_send_int4_2d
     module procedure mpp_send_int4_3d
     module procedure mpp_send_int4_4d
     module procedure mpp_send_int4_5d
     module procedure mpp_send_logical4
     module procedure mpp_send_logical4_scalar
     module procedure mpp_send_logical4_2d
     module procedure mpp_send_logical4_3d
     module procedure mpp_send_logical4_4d
     module procedure mpp_send_logical4_5d
  end interface

  ! <INTERFACE NAME="mpp_broadcast">

  !   <OVERVIEW>
  !     Parallel broadcasts.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     The <TT>mpp_broadcast</TT> call has been added because the original
  !     syntax (using <TT>ALL_PES</TT> in <TT>mpp_transmit</TT>) did not
  !     support a broadcast across a pelist.
  !
  !     <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !     <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
  !     contiguous block from a multi-dimensional array may be passed by its
  !     starting address and its length, as in <TT>f77</TT>.
  !
  !     Global broadcasts through the <TT>ALL_PES</TT> argument to <LINK
  !     SRC="#mpp_transmit"><TT>mpp_transmit</TT></LINK> are still provided for
  !     backward-compatibility.
  !
  !     If <TT>pelist</TT> is omitted, the context is assumed to be the
  !     current pelist. <TT>from_pe</TT> must belong to the current
  !     pelist. This call implies synchronization across the PEs in
  !     <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call mpp_broadcast( data, length, from_pe, pelist )
  !   </TEMPLATE>
  !   <IN NAME="length"> </IN>
  !   <IN NAME="from_pe"> </IN>
  !   <IN NAME="pelist"> </IN>
  !   <INOUT NAME="data(*)"> </INOUT>
  ! </INTERFACE>
  interface mpp_broadcast
     module procedure mpp_broadcast_real8
     module procedure mpp_broadcast_real8_scalar
     module procedure mpp_broadcast_real8_2d
     module procedure mpp_broadcast_real8_3d
     module procedure mpp_broadcast_real8_4d
     module procedure mpp_broadcast_real8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_broadcast_cmplx8
     module procedure mpp_broadcast_cmplx8_scalar
     module procedure mpp_broadcast_cmplx8_2d
     module procedure mpp_broadcast_cmplx8_3d
     module procedure mpp_broadcast_cmplx8_4d
     module procedure mpp_broadcast_cmplx8_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_broadcast_int8
     module procedure mpp_broadcast_int8_scalar
     module procedure mpp_broadcast_int8_2d
     module procedure mpp_broadcast_int8_3d
     module procedure mpp_broadcast_int8_4d
     module procedure mpp_broadcast_int8_5d
     module procedure mpp_broadcast_logical8
     module procedure mpp_broadcast_logical8_scalar
     module procedure mpp_broadcast_logical8_2d
     module procedure mpp_broadcast_logical8_3d
     module procedure mpp_broadcast_logical8_4d
     module procedure mpp_broadcast_logical8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_broadcast_real4
     module procedure mpp_broadcast_real4_scalar
     module procedure mpp_broadcast_real4_2d
     module procedure mpp_broadcast_real4_3d
     module procedure mpp_broadcast_real4_4d
     module procedure mpp_broadcast_real4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_broadcast_cmplx4
     module procedure mpp_broadcast_cmplx4_scalar
     module procedure mpp_broadcast_cmplx4_2d
     module procedure mpp_broadcast_cmplx4_3d
     module procedure mpp_broadcast_cmplx4_4d
     module procedure mpp_broadcast_cmplx4_5d
#endif
     module procedure mpp_broadcast_int4
     module procedure mpp_broadcast_int4_scalar
     module procedure mpp_broadcast_int4_2d
     module procedure mpp_broadcast_int4_3d
     module procedure mpp_broadcast_int4_4d
     module procedure mpp_broadcast_int4_5d
     module procedure mpp_broadcast_logical4
     module procedure mpp_broadcast_logical4_scalar
     module procedure mpp_broadcast_logical4_2d
     module procedure mpp_broadcast_logical4_3d
     module procedure mpp_broadcast_logical4_4d
     module procedure mpp_broadcast_logical4_5d
  end interface

  !#####################################################################
  ! <INTERFACE NAME="mpp_chksum">

  !   <OVERVIEW>
  !     Parallel checksums.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     <TT>mpp_chksum</TT> is a parallel checksum routine that returns an
  !     identical answer for the same array irrespective of how it has been
  !     partitioned across processors. <TT>LONG_KIND</TT>is the <TT>KIND</TT>
  !     parameter corresponding to long integers (see discussion on
  !     OS-dependent preprocessor directives) defined in
  !     the header file <TT>fms_platform.h</TT>. <TT>MPP_TYPE_</TT> corresponds to any
  !     4-byte and 8-byte variant of <TT>integer, real, complex, logical</TT>
  !     variables, of rank 0 to 5.
  !
  !     Integer checksums on FP data use the F90 <TT>TRANSFER()</TT>
  !     intrinsic.
  !
  !     The <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/chksum/chksum.html">serial checksum module</LINK> is superseded
  !     by this function, and is no longer being actively maintained. This
  !     provides identical results on a single-processor job, and to perform
  !     serial checksums on a single processor of a parallel job, you only
  !     need to use the optional <TT>pelist</TT> argument.
  !     <PRE>
  !     use mpp_mod
  !     integer :: pe, chksum
  !     real :: a(:)
  !     pe = mpp_pe()
  !     chksum = mpp_chksum( a, (/pe/) )
  !     </PRE>
  !
  !     The additional functionality of <TT>mpp_chksum</TT> over
  !     serial checksums is to compute the checksum across the PEs in
  !     <TT>pelist</TT>. The answer is guaranteed to be the same for
  !     the same distributed array irrespective of how it has been
  !     partitioned.
  !
  !     If <TT>pelist</TT> is omitted, the context is assumed to be the
  !     current pelist. This call implies synchronization across the PEs in
  !     <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     mpp_chksum( var, pelist )
  !   </TEMPLATE>
  !   <IN NAME="pelist" TYPE="integer" DIM="(:)"> </IN>
  !   <IN NAME="var" TYPE="MPP_TYPE_"> </IN>
  ! </INTERFACE>
  interface mpp_chksum
#ifndef no_8byte_integers
     module procedure mpp_chksum_i8_1d
     module procedure mpp_chksum_i8_2d
     module procedure mpp_chksum_i8_3d
     module procedure mpp_chksum_i8_4d
#endif
     module procedure mpp_chksum_i4_1d
     module procedure mpp_chksum_i4_2d
     module procedure mpp_chksum_i4_3d
     module procedure mpp_chksum_i4_4d
     module procedure mpp_chksum_r8_0d
     module procedure mpp_chksum_r8_1d
     module procedure mpp_chksum_r8_2d
     module procedure mpp_chksum_r8_3d
     module procedure mpp_chksum_r8_4d
     module procedure mpp_chksum_r8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_chksum_c8_0d
     module procedure mpp_chksum_c8_1d
     module procedure mpp_chksum_c8_2d
     module procedure mpp_chksum_c8_3d
     module procedure mpp_chksum_c8_4d
     module procedure mpp_chksum_c8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_chksum_r4_0d
     module procedure mpp_chksum_r4_1d
     module procedure mpp_chksum_r4_2d
     module procedure mpp_chksum_r4_3d
     module procedure mpp_chksum_r4_4d
     module procedure mpp_chksum_r4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_chksum_c4_0d
     module procedure mpp_chksum_c4_1d
     module procedure mpp_chksum_c4_2d
     module procedure mpp_chksum_c4_3d
     module procedure mpp_chksum_c4_4d
     module procedure mpp_chksum_c4_5d
#endif
  end interface

!***********************************************************************
!
!            module variables 
!
!***********************************************************************
  type(communicator),save :: peset(0:PESET_MAX) !0 is a dummy used to hold single-PE "self" communicator
  logical              :: module_is_initialized = .false.
  logical              :: debug = .false.
  integer              :: npes=1, root_pe=0, pe=0
  integer(LONG_KIND)   :: tick, ticks_per_sec, max_ticks, start_tick, end_tick, tick0=0
  integer              :: mpp_comm_private
  logical              :: first_call_system_clock_mpi=.TRUE.
  real(DOUBLE_KIND)    :: mpi_count0=0  ! use to prevent integer overflow
  real(DOUBLE_KIND)    :: mpi_tick_rate=0.d0  ! clock rate for mpi_wtick()
  logical              :: mpp_record_timing_data=.TRUE.
  type(clock),save     :: clocks(MAX_CLOCKS)
  integer              :: log_unit, etc_unit
  character(len=32)    :: configfile='logfile'
  integer              :: peset_num=0, current_peset_num=0
  integer              :: world_peset_num                  !the world communicator
  integer              :: error
  integer              :: clock_num=0, num_clock_ids=0,current_clock=0, previous_clock(MAX_CLOCKS)=0
  real                 :: tick_rate
  integer, allocatable :: request(:)
  integer, allocatable :: request_recv(:)
! if you want to save the non-root PE information uncomment out the following line
! and comment out the assigment of etcfile to '/dev/null'
#ifdef NO_DEV_NULL
  character(len=32)    :: etcfile='._mpp.nonrootpe.msgs'
#else
  character(len=32)    :: etcfile='/dev/null'
#endif

#ifdef SGICRAY
  integer :: in_unit=100, out_unit=101, err_unit=102 !see intro_io(3F): to see why these values are used rather than 5,6,0
#else
  integer :: in_unit=5, out_unit=6, err_unit=0
#endif

  !--- variables used in mpp_util.h
  type(Summary_Struct) :: clock_summary(MAX_CLOCKS)
  logical              :: warnings_are_fatal = .FALSE.
  integer              :: error_state=0
  integer              :: clock_grain=CLOCK_LOOP-1

  !--- variables used in mpp_comm.h
#ifdef use_libMPI
#ifdef _CRAYT3E
  !BWA: mpif.h on t3e currently does not contain MPI_INTEGER8 datatype
  !(O2k and t90 do)
  !(t3e: fixed on 3.3 I believe)
  integer, parameter :: MPI_INTEGER8=MPI_INTEGER
#endif
#endif /* use_libMPI */
#ifdef use_MPI_SMA
#include <mpp/shmem.fh>
  integer :: pSync(SHMEM_BARRIER_SYNC_SIZE)
  pointer( p_pSync, pSync ) !used by SHPALLOC
#endif

  integer            :: clock0    !measures total runtime from mpp_init to mpp_exit
  integer            :: mpp_stack_size=0, mpp_stack_hwm=0
  integer            :: tag=1
  logical            :: verbose=.FALSE.
#ifdef _CRAY
  integer(LONG_KIND) :: word(1)
#endif
#if defined(sgi_mipspro) || defined(__ia64)
  integer(INT_KIND)  :: word(1)
#endif

  character(len=128), public :: version= &
       '$Id mpp.F90 $'
  character(len=128), public :: tagname= &
       '$Name: testing $'

  contains
#include <system_clock.h>
#include <mpp_util.inc>
#include <mpp_comm.inc>

  end module mpp_mod




