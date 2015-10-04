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

! -*-f90-*-
! $Id: mpp_comm_mpi.inc,v 19.0.2.1 2012/05/08 17:45:39 Zhi.Liang Exp $


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    subroutine mpp_init( flags, in, out, err, log )
  !      integer, optional, intent(in) :: flags, in, out, err, log
  subroutine mpp_init( flags, localcomm )
  integer, optional, intent(in) :: flags
  integer, optional, intent(in) :: localcomm
  integer                       :: my_pe, num_pes, len, i, iunit
  logical                       :: opened, existed
  integer                       :: unit_begin, unit_end, unit_nml, io_status
  character(len=5) :: this_pe

  if( module_is_initialized )return

  call MPI_INITIALIZED( opened, error ) !in case called from another MPI package
  if(opened .and. .NOT. PRESENT(localcomm)) call mpp_error( FATAL, 'MPP_INIT: communicator is required' )
  if( .NOT.opened ) then
    call MPI_INIT(error)
    mpp_comm_private = MPI_COMM_WORLD
  else
    mpp_comm_private  = localcomm
  endif
  call MPI_COMM_RANK( mpp_comm_private, pe,   error )
  call MPI_COMM_SIZE( mpp_comm_private, npes, error )
#ifdef use_MPI_SMA
  call SHMEM_BARRIER_ALL()
  call SHPALLOC( p_pSync, SHMEM_BARRIER_SYNC_SIZE, error, -1 )
  call SHMEM_BARRIER_ALL()
#endif

  module_is_initialized = .TRUE.


  !PEsets: make defaults illegal
  allocate(peset(0:current_peset_max))
  peset(:)%count = -1
  peset(:)%id = -1
  peset(:)%group = -1
  peset(:)%start = -1
  peset(:)%log2stride = -1
  peset(:)%name  = " "
  !0=single-PE, initialized so that count returns 1
  peset(0)%count = 1
  allocate( peset(0)%list(1) )
  peset(0)%list = pe
  current_peset_num = 0
  peset(0)%id = mpp_comm_private
  call MPI_COMM_GROUP( mpp_comm_private, peset(0)%group, error )
  world_peset_num = get_peset( (/(i,i=0,npes-1)/) )
  current_peset_num = world_peset_num !initialize current PEset to world

  !initialize clocks
  call SYSTEM_CLOCK( count=tick0, count_rate=ticks_per_sec, count_max=max_ticks )
  tick_rate = 1./ticks_per_sec
  clock0 = mpp_clock_id( 'Total runtime', flags=MPP_CLOCK_SYNC )

  if( PRESENT(flags) )then
     debug   = flags.EQ.MPP_DEBUG
     verbose = flags.EQ.MPP_VERBOSE .OR. debug
  end if

  call mpp_init_logfile()
  call read_input_nml

  !--- read namelist
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, mpp_nml, iostat=io_status)
#else
  unit_begin = 103
  unit_end   = 512
  do unit_nml = unit_begin, unit_end
     inquire( unit_nml,OPENED=opened )
     if( .NOT.opened )exit
  end do

  open(unit_nml,file='input.nml', iostat=io_status)
  read(unit_nml,mpp_nml,iostat=io_status)
  close(unit_nml)  
#endif

! non-root pe messages written to other location than stdout()

  if(etc_unit_is_stderr) then
     etc_unit = stderr()
  else
! 9 is reserved for etc_unit
     etc_unit=9
     inquire(unit=etc_unit,opened=opened)
     if(opened) call mpp_error(FATAL,'Unit 9 is already in use (etc_unit) in mpp_comm_mpi')
     if (trim(etcfile) /= '/dev/null') then
        write( etcfile,'(a,i6.6)' )trim(etcfile)//'.', pe
     endif
     inquire(file=etcfile, exist=existed)
     if(existed) then
        open( unit=etc_unit, file=trim(etcfile), status='REPLACE' )
     else
        open( unit=etc_unit, file=trim(etcfile) )
     endif
  endif

! max_request is set to maximum of npes * REQUEST_MULTIPLY ( default is 20) and MAX_REQUEST_MIN ( default 10000)
  max_request = max(MAX_REQUEST_MIN, mpp_npes()*REQUEST_MULTIPLY)

  allocate( request_send(max_request) )
  allocate( request_recv(max_request) )
  allocate( size_recv(max_request) )
  allocate( type_recv(max_request) )
  request_send(:) = MPI_REQUEST_NULL
  request_recv(:) = MPI_REQUEST_NULL
  size_recv(:)    = 0
  type_recv(:)    = 0

  !if optional argument logunit=stdout, write messages to stdout instead.
  !if specifying non-defaults, you must specify units not yet in use.
  !      if( PRESENT(in) )then
  !          inquire( unit=in, opened=opened )
  !          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stdin.' )
  !          in_unit=in
  !      end if
  !      if( PRESENT(out) )then
  !          inquire( unit=out, opened=opened )
  !          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stdout.' )
  !          out_unit=out
  !      end if
  !      if( PRESENT(err) )then
  !          inquire( unit=err, opened=opened )
  !          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stderr.' )
  !          err_unit=err
  !      end if
  !      log_unit=get_unit()
  !      if( PRESENT(log) )then
  !          inquire( unit=log, opened=opened )
  !          if( opened .AND. log.NE.out_unit )call mpp_error( FATAL, 'MPP_INIT: unable to open stdlog.' )
  !          log_unit=log
  !      end if
  !!log_unit can be written to only from root_pe, all others write to stdout
  !      if( log_unit.NE.out_unit )then
  !          inquire( unit=log_unit, opened=opened )
  !          if( opened )call mpp_error( FATAL, 'MPP_INIT: specified unit for stdlog already in use.' )
  !          if( pe.EQ.root_pe )open( unit=log_unit, file=trim(configfile), status='REPLACE' )
  !          call mpp_sync()
  !          if( pe.NE.root_pe )open( unit=log_unit, file=trim(configfile), status='OLD' )
  !      end if


  !messages
  iunit = stdlog()  ! workaround for lf95.
  if( verbose )call mpp_error( NOTE, 'MPP_INIT: initializing MPP module...' )
  if( pe.EQ.root_pe )then
     write( iunit,'(/a)' )'MPP module '//trim(version)//trim(tagname)
     write( iunit,'(a,i6)' )'MPP started with NPES=', npes
     write( iunit,'(a)' )'Using MPI library for message passing...'
     write( iunit, '(a,es12.4,a,i10,a)' ) &
          'Realtime clock resolution=', tick_rate, ' sec (', ticks_per_sec, ' ticks/sec)'
     write( iunit, '(a,es12.4,a,i20,a)' ) &
          'Clock rolls over after ', max_ticks*tick_rate, ' sec (', max_ticks, ' ticks)'
     write( iunit,'(/a)' )'MPP module '//trim(mpp_parameter_version)//trim(mpp_parameter_tagname)
     write( iunit,'(/a)' )'MPP module '//trim(mpp_data_version)//trim(mpp_data_tagname)
  end if

  call mpp_clock_begin(clock0)

  return
end subroutine mpp_init

!#######################################################################
  !to be called at the end of a run
subroutine mpp_exit()
  integer :: i, j, k, n, nmax, istat, out_unit, log_unit
  real    :: t, tmin, tmax, tavg, tstd
  real    :: m, mmin, mmax, mavg, mstd, t_total
  logical :: opened

  if( .NOT.module_is_initialized )return
  call mpp_set_current_pelist()
  call mpp_clock_end(clock0)
  t_total = clocks(clock0)%total_ticks*tick_rate
  out_unit = stdout()
  log_unit = stdlog()
  if( clock_num.GT.0 )then
     if( ANY(clocks(1:clock_num)%detailed) )then
        call sum_clock_data; call dump_clock_summary
     end if

     call mpp_sync()
     call FLUSH( out_unit )

     if( pe.EQ.root_pe )then
        write( out_unit,'(/a,i6,a)' ) 'Tabulating mpp_clock statistics across ', npes, ' PEs...'
        if( ANY(clocks(1:clock_num)%detailed) ) &
             write( out_unit,'(a)' )'   ... see mpp_clock.out.#### for details on individual PEs.'
        write( out_unit,'(/32x,a)' ) '          tmin          tmax          tavg          tstd  tfrac grain pemin pemax'
     end if
     write( log_unit,'(/37x,a)' ) 'time'

     call FLUSH( out_unit )
     call mpp_sync()
     do i = 1,clock_num
        if( .NOT.ANY(peset(clocks(i)%peset_num)%list(:).EQ.pe) )cycle
        call mpp_set_current_pelist( peset(clocks(i)%peset_num)%list )
        out_unit = stdout()
        log_unit = stdlog()
        !times between mpp_clock ticks
        t = clocks(i)%total_ticks*tick_rate
        tmin = t; call mpp_min(tmin)
        tmax = t; call mpp_max(tmax)
        tavg = t; call mpp_sum(tavg); tavg = tavg/mpp_npes()
        tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt( tstd/mpp_npes() )
        if( pe.EQ.root_pe )write( out_unit,'(a32,4f14.6,f7.3,3i6)' ) &
             clocks(i)%name, tmin, tmax, tavg, tstd, tavg/t_total, &
             clocks(i)%grain, minval(peset(clocks(i)%peset_num)%list), &
             maxval(peset(clocks(i)%peset_num)%list)
        write(log_unit,'(a32,f14.6)') clocks(i)%name, clocks(i)%total_ticks*tick_rate
     end do

     if( ANY(clocks(1:clock_num)%detailed) .AND. pe.EQ.root_pe )write( out_unit,'(/32x,a)' ) &
          '       tmin       tmax       tavg       tstd       mmin       mmax       mavg       mstd  mavg/tavg'

     do i = 1,clock_num
        !messages: bytelengths and times
        if( .NOT.clocks(i)%detailed )cycle
        if( .NOT.ANY(peset(clocks(i)%peset_num)%list(:).EQ.pe) )cycle
        call mpp_set_current_pelist( peset(clocks(i)%peset_num)%list )
        out_unit = stdout()
        do j = 1,MAX_EVENT_TYPES
           n = clocks(i)%events(j)%calls; nmax = n
           call mpp_max(nmax)
           if( nmax.NE.0 )then
              !don't divide by n because n might be 0
              m = 0
              if( n.GT.0 )m = sum(clocks(i)%events(j)%bytes(1:n))
              mmin = m; call mpp_min(mmin)
              mmax = m; call mpp_max(mmax)
              mavg = m; call mpp_sum(mavg); mavg = mavg/mpp_npes()
              mstd = (m-mavg)**2; call mpp_sum(mstd); mstd = sqrt( mstd/mpp_npes() )
              t = 0
              if( n.GT.0 )t = sum(clocks(i)%events(j)%ticks(1:n))*tick_rate
              tmin = t; call mpp_min(tmin)
              tmax = t; call mpp_max(tmax)
              tavg = t; call mpp_sum(tavg); tavg = tavg/mpp_npes()
              tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt( tstd/mpp_npes() )
              if( pe.EQ.root_pe )write( out_unit,'(a32,4f11.3,5es11.3)' ) &
                   trim(clocks(i)%name)//' '//trim(clocks(i)%events(j)%name), &
                   tmin, tmax, tavg, tstd, mmin, mmax, mavg, mstd, mavg/tavg
           end if
        end do
     end do

  end if

  call FLUSH( out_unit )

! close down etc_unit: 9
  inquire(unit=etc_unit, opened=opened)
  if (opened) then
   call FLUSH (etc_unit)
   close(etc_unit)
  endif

  call mpp_set_current_pelist()
  call mpp_sync()
  call mpp_max(mpp_stack_hwm)
  if( pe.EQ.root_pe )write( out_unit,* )'MPP_STACK high water mark=', mpp_stack_hwm
  if(mpp_comm_private == MPI_COMM_WORLD ) call MPI_FINALIZE(error)



  return
end subroutine mpp_exit

!#######################################################################
    subroutine mpp_malloc( ptr, newlen, len )
      integer, intent(in)    :: newlen
      integer, intent(inout) :: len
#ifdef use_MPI_SMA
      real                   :: dummy
      pointer( ptr, dummy )
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_MALLOC: You must first call mpp_init.' )
!use existing allocation if it is enough
      if( newlen.LE.len )return

      call SHMEM_BARRIER_ALL()
!if the pointer is already allocated, deallocate
!      if( len.NE.0 )call SHPDEALLC( ptr, error_8, -1 ) !BWA: error_8 instead of error, see PV 682618 (fixed in mpt.1.3.0.1)
      if( len.NE.0 )call SHPDEALLC( ptr, error, -1 )
!allocate new length: assume that the array is KIND=8
      call SHPALLOC( ptr, newlen, error, -1 )
      len = newlen
      call SHMEM_BARRIER_ALL()
#else
      integer(POINTER_KIND), intent(in) :: ptr
      call mpp_error( FATAL, 'mpp_malloc: requires use_MPI_SMA' )
#endif
      return
    end subroutine mpp_malloc

#ifdef use_MPI_GSM
!#######################################################################
 !--- routine to perform GSM allocations
 !
subroutine mpp_gsm_malloc( ptr, len )
  use mpp_data_mod,      only : peset, current_peset_num
  integer(KIND=MPI_ADDRESS_KIND), intent(in)    :: len
  real                   :: dummy
  integer                :: ierror
  !argument ptr is a cray pointer, points to a dummy argument in this routine 
  pointer( ptr, dummy )
  include "mpi_gsmf.h"
  if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GSM_MALLOC: You must first call mpp_init.' )

  call gsm_alloc(len, GSM_SINGLERANK, 0, peset(current_peset_num)%id, ptr, ierror)
  if (ierror .eq. -1) call mpp_error( FATAL, 'MPP_GSM_MALLOC: gsm_alloc failed.' )
  return
end subroutine mpp_gsm_malloc

!#######################################################################
 !--- routine to free GSM allocations
 !
subroutine mpp_gsm_free( ptr )
  use mpp_data_mod,      only : peset, current_peset_num
  real                   :: dummy
  integer                :: ierror
  !argument ptr is a cray pointer, points to a dummy argument in this routine 
  pointer( ptr, dummy )
  include "mpi_gsmf.h"

  call gsm_free(dummy, peset(current_peset_num)%id, ierror)
  if (ierror .eq. -1) call mpp_error( FATAL, 'MPP_GSM_FREE: gsm_free failed.' )
  return
end subroutine mpp_gsm_free
#endif

!#######################################################################
    !set the mpp_stack variable to be at least n LONG words long
  subroutine mpp_set_stack_size(n)
    integer, intent(in) :: n
    character(len=8)    :: text

    if( n.GT.mpp_stack_size .AND. allocated(mpp_stack) )deallocate(mpp_stack)
    if( .NOT.allocated(mpp_stack) )then
       allocate( mpp_stack(n) )
       mpp_stack_size = n
    end if

    write( text,'(i8)' )n
    if( pe.EQ.root_pe )call mpp_error( NOTE, 'MPP_SET_STACK_SIZE: stack size set to '//text//'.' )

    return
  end subroutine mpp_set_stack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                BASIC MESSAGE PASSING ROUTINE: mpp_transmit                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_broadcast_char(data, length, from_pe, pelist )
      character(len=*), intent(inout) :: data(:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n, i, from_rank, out_unit
      character    :: str1D(length*size(data(:)))
      pointer(lptr, str1D)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'mpp_broadcast_text: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return


      if( debug )then
          call SYSTEM_CLOCK(tick)
          if(mpp_pe() == mpp_root_pe()) then
            out_unit = stdout()
            write( out_unit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, 'mpp_broadcast_text begin: from_pe, length=', from_pe, length
          endif
      end if

      if( .NOT.ANY(from_pe.EQ.peset(current_peset_num)%list) ) &
           call mpp_error( FATAL, 'mpp_broadcast_text: broadcasting from invalid PE.' )

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
 ! find the rank of from_pe in the pelist.     
      do i = 1, mpp_npes()
         if(peset(n)%list(i) == from_pe) then
             from_rank = i - 1
             exit
         endif
      enddo
      lptr = LOC (data)
      if( mpp_npes().GT.1 ) call MPI_BCAST( data, length*size(data(:)), MPI_CHARACTER, from_rank, peset(n)%id, error )
      if( current_clock.NE.0 )call increment_current_clock( EVENT_BROADCAST, length )
      return
    end subroutine mpp_broadcast_char

#undef MPP_TRANSMIT_
#define MPP_TRANSMIT_ mpp_transmit_real8
#undef MPP_TRANSMIT_SCALAR_
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_real8_scalar
#undef MPP_TRANSMIT_2D_
#define MPP_TRANSMIT_2D_ mpp_transmit_real8_2d
#undef MPP_TRANSMIT_3D_
#define MPP_TRANSMIT_3D_ mpp_transmit_real8_3d
#undef MPP_TRANSMIT_4D_
#define MPP_TRANSMIT_4D_ mpp_transmit_real8_4d
#undef MPP_TRANSMIT_5D_
#define MPP_TRANSMIT_5D_ mpp_transmit_real8_5d
#undef MPP_RECV_
#define MPP_RECV_ mpp_recv_real8
#undef MPP_RECV_SCALAR_
#define MPP_RECV_SCALAR_ mpp_recv_real8_scalar
#undef MPP_RECV_2D_
#define MPP_RECV_2D_ mpp_recv_real8_2d
#undef MPP_RECV_3D_
#define MPP_RECV_3D_ mpp_recv_real8_3d
#undef MPP_RECV_4D_
#define MPP_RECV_4D_ mpp_recv_real8_4d
#undef MPP_RECV_5D_
#define MPP_RECV_5D_ mpp_recv_real8_5d
#undef MPP_SEND_
#define MPP_SEND_ mpp_send_real8
#undef MPP_SEND_SCALAR_
#define MPP_SEND_SCALAR_ mpp_send_real8_scalar
#undef MPP_SEND_2D_
#define MPP_SEND_2D_ mpp_send_real8_2d
#undef MPP_SEND_3D_
#define MPP_SEND_3D_ mpp_send_real8_3d
#undef MPP_SEND_4D_
#define MPP_SEND_4D_ mpp_send_real8_4d
#undef MPP_SEND_5D_
#define MPP_SEND_5D_ mpp_send_real8_5d
#undef MPP_BROADCAST_
#define MPP_BROADCAST_ mpp_broadcast_real8
#undef MPP_BROADCAST_SCALAR_
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_real8_scalar
#undef MPP_BROADCAST_2D_
#define MPP_BROADCAST_2D_ mpp_broadcast_real8_2d
#undef MPP_BROADCAST_3D_
#define MPP_BROADCAST_3D_ mpp_broadcast_real8_3d
#undef MPP_BROADCAST_4D_
#define MPP_BROADCAST_4D_ mpp_broadcast_real8_4d
#undef MPP_BROADCAST_5D_
#define MPP_BROADCAST_5D_ mpp_broadcast_real8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_REAL8
#include <mpp_transmit_mpi.h>

#ifdef OVERLOAD_C8
#undef MPP_TRANSMIT_
#define MPP_TRANSMIT_ mpp_transmit_cmplx8
#undef MPP_TRANSMIT_SCALAR_
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_cmplx8_scalar
#undef MPP_TRANSMIT_2D_
#define MPP_TRANSMIT_2D_ mpp_transmit_cmplx8_2d
#undef MPP_TRANSMIT_3D_
#define MPP_TRANSMIT_3D_ mpp_transmit_cmplx8_3d
#undef MPP_TRANSMIT_4D_
#define MPP_TRANSMIT_4D_ mpp_transmit_cmplx8_4d
#undef MPP_TRANSMIT_5D_
#define MPP_TRANSMIT_5D_ mpp_transmit_cmplx8_5d
#undef MPP_RECV_
#define MPP_RECV_ mpp_recv_cmplx8
#undef MPP_RECV_SCALAR_
#define MPP_RECV_SCALAR_ mpp_recv_cmplx8_scalar
#undef MPP_RECV_2D_
#define MPP_RECV_2D_ mpp_recv_cmplx8_2d
#undef MPP_RECV_3D_
#define MPP_RECV_3D_ mpp_recv_cmplx8_3d
#undef MPP_RECV_4D_
#define MPP_RECV_4D_ mpp_recv_cmplx8_4d
#undef MPP_RECV_5D_
#define MPP_RECV_5D_ mpp_recv_cmplx8_5d
#undef MPP_SEND_
#define MPP_SEND_ mpp_send_cmplx8
#undef MPP_SEND_SCALAR_
#define MPP_SEND_SCALAR_ mpp_send_cmplx8_scalar
#undef MPP_SEND_2D_
#define MPP_SEND_2D_ mpp_send_cmplx8_2d
#undef MPP_SEND_3D_
#define MPP_SEND_3D_ mpp_send_cmplx8_3d
#undef MPP_SEND_4D_
#define MPP_SEND_4D_ mpp_send_cmplx8_4d
#undef MPP_SEND_5D_
#define MPP_SEND_5D_ mpp_send_cmplx8_5d
#undef MPP_BROADCAST_
#define MPP_BROADCAST_ mpp_broadcast_cmplx8
#undef MPP_BROADCAST_SCALAR_
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_cmplx8_scalar
#undef MPP_BROADCAST_2D_
#define MPP_BROADCAST_2D_ mpp_broadcast_cmplx8_2d
#undef MPP_BROADCAST_3D_
#define MPP_BROADCAST_3D_ mpp_broadcast_cmplx8_3d
#undef MPP_BROADCAST_4D_
#define MPP_BROADCAST_4D_ mpp_broadcast_cmplx8_4d
#undef MPP_BROADCAST_5D_
#define MPP_BROADCAST_5D_ mpp_broadcast_cmplx8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 16
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_DOUBLE_COMPLEX
#include <mpp_transmit_mpi.h>
#endif

#undef MPP_TRANSMIT_
#define MPP_TRANSMIT_ mpp_transmit_real4
#undef MPP_TRANSMIT_SCALAR_
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_real4_scalar
#undef MPP_TRANSMIT_2D_
#define MPP_TRANSMIT_2D_ mpp_transmit_real4_2d
#undef MPP_TRANSMIT_3D_
#define MPP_TRANSMIT_3D_ mpp_transmit_real4_3d
#undef MPP_TRANSMIT_4D_
#define MPP_TRANSMIT_4D_ mpp_transmit_real4_4d
#undef MPP_TRANSMIT_5D_
#define MPP_TRANSMIT_5D_ mpp_transmit_real4_5d
#undef MPP_RECV_
#define MPP_RECV_ mpp_recv_real4
#undef MPP_RECV_SCALAR_
#define MPP_RECV_SCALAR_ mpp_recv_real4_scalar
#undef MPP_RECV_2D_
#define MPP_RECV_2D_ mpp_recv_real4_2d
#undef MPP_RECV_3D_
#define MPP_RECV_3D_ mpp_recv_real4_3d
#undef MPP_RECV_4D_
#define MPP_RECV_4D_ mpp_recv_real4_4d
#undef MPP_RECV_5D_
#define MPP_RECV_5D_ mpp_recv_real4_5d
#undef MPP_SEND_
#define MPP_SEND_ mpp_send_real4
#undef MPP_SEND_SCALAR_
#define MPP_SEND_SCALAR_ mpp_send_real4_scalar
#undef MPP_SEND_2D_
#define MPP_SEND_2D_ mpp_send_real4_2d
#undef MPP_SEND_3D_
#define MPP_SEND_3D_ mpp_send_real4_3d
#undef MPP_SEND_4D_
#define MPP_SEND_4D_ mpp_send_real4_4d
#undef MPP_SEND_5D_
#define MPP_SEND_5D_ mpp_send_real4_5d
#undef MPP_BROADCAST_
#define MPP_BROADCAST_ mpp_broadcast_real4
#undef MPP_BROADCAST_SCALAR_
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_real4_scalar
#undef MPP_BROADCAST_2D_
#define MPP_BROADCAST_2D_ mpp_broadcast_real4_2d
#undef MPP_BROADCAST_3D_
#define MPP_BROADCAST_3D_ mpp_broadcast_real4_3d
#undef MPP_BROADCAST_4D_
#define MPP_BROADCAST_4D_ mpp_broadcast_real4_4d
#undef MPP_BROADCAST_5D_
#define MPP_BROADCAST_5D_ mpp_broadcast_real4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_REAL4
#include <mpp_transmit_mpi.h>

#ifdef OVERLOAD_C4
#undef MPP_TRANSMIT_
#define MPP_TRANSMIT_ mpp_transmit_cmplx4
#undef MPP_TRANSMIT_SCALAR_
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_cmplx4_scalar
#undef MPP_TRANSMIT_2D_
#define MPP_TRANSMIT_2D_ mpp_transmit_cmplx4_2d
#undef MPP_TRANSMIT_3D_
#define MPP_TRANSMIT_3D_ mpp_transmit_cmplx4_3d
#undef MPP_TRANSMIT_4D_
#define MPP_TRANSMIT_4D_ mpp_transmit_cmplx4_4d
#undef MPP_TRANSMIT_5D_
#define MPP_TRANSMIT_5D_ mpp_transmit_cmplx4_5d
#undef MPP_RECV_
#define MPP_RECV_ mpp_recv_cmplx4
#undef MPP_RECV_SCALAR_
#define MPP_RECV_SCALAR_ mpp_recv_cmplx4_scalar
#undef MPP_RECV_2D_
#define MPP_RECV_2D_ mpp_recv_cmplx4_2d
#undef MPP_RECV_3D_
#define MPP_RECV_3D_ mpp_recv_cmplx4_3d
#undef MPP_RECV_4D_
#define MPP_RECV_4D_ mpp_recv_cmplx4_4d
#undef MPP_RECV_5D_
#define MPP_RECV_5D_ mpp_recv_cmplx4_5d
#undef MPP_SEND_
#define MPP_SEND_ mpp_send_cmplx4
#undef MPP_SEND_SCALAR_
#define MPP_SEND_SCALAR_ mpp_send_cmplx4_scalar
#undef MPP_SEND_2D_
#define MPP_SEND_2D_ mpp_send_cmplx4_2d
#undef MPP_SEND_3D_
#define MPP_SEND_3D_ mpp_send_cmplx4_3d
#undef MPP_SEND_4D_
#define MPP_SEND_4D_ mpp_send_cmplx4_4d
#undef MPP_SEND_5D_
#define MPP_SEND_5D_ mpp_send_cmplx4_5d
#undef MPP_BROADCAST_
#define MPP_BROADCAST_ mpp_broadcast_cmplx4
#undef MPP_BROADCAST_SCALAR_
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_cmplx4_scalar
#undef MPP_BROADCAST_2D_
#define MPP_BROADCAST_2D_ mpp_broadcast_cmplx4_2d
#undef MPP_BROADCAST_3D_
#define MPP_BROADCAST_3D_ mpp_broadcast_cmplx4_3d
#undef MPP_BROADCAST_4D_
#define MPP_BROADCAST_4D_ mpp_broadcast_cmplx4_4d
#undef MPP_BROADCAST_5D_
#define MPP_BROADCAST_5D_ mpp_broadcast_cmplx4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_COMPLEX
#include <mpp_transmit_mpi.h>
#endif 

#ifndef no_8byte_integers
#undef MPP_TRANSMIT_
#define MPP_TRANSMIT_ mpp_transmit_int8
#undef MPP_TRANSMIT_SCALAR_
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_int8_scalar
#undef MPP_TRANSMIT_2D_
#define MPP_TRANSMIT_2D_ mpp_transmit_int8_2d
#undef MPP_TRANSMIT_3D_
#define MPP_TRANSMIT_3D_ mpp_transmit_int8_3d
#undef MPP_TRANSMIT_4D_
#define MPP_TRANSMIT_4D_ mpp_transmit_int8_4d
#undef MPP_TRANSMIT_5D_
#define MPP_TRANSMIT_5D_ mpp_transmit_int8_5d
#undef MPP_RECV_
#define MPP_RECV_ mpp_recv_int8
#undef MPP_RECV_SCALAR_
#define MPP_RECV_SCALAR_ mpp_recv_int8_scalar
#undef MPP_RECV_2D_
#define MPP_RECV_2D_ mpp_recv_int8_2d
#undef MPP_RECV_3D_
#define MPP_RECV_3D_ mpp_recv_int8_3d
#undef MPP_RECV_4D_
#define MPP_RECV_4D_ mpp_recv_int8_4d
#undef MPP_RECV_5D_
#define MPP_RECV_5D_ mpp_recv_int8_5d
#undef MPP_SEND_
#define MPP_SEND_ mpp_send_int8
#undef MPP_SEND_SCALAR_
#define MPP_SEND_SCALAR_ mpp_send_int8_scalar
#undef MPP_SEND_2D_
#define MPP_SEND_2D_ mpp_send_int8_2d
#undef MPP_SEND_3D_
#define MPP_SEND_3D_ mpp_send_int8_3d
#undef MPP_SEND_4D_
#define MPP_SEND_4D_ mpp_send_int8_4d
#undef MPP_SEND_5D_
#define MPP_SEND_5D_ mpp_send_int8_5d
#undef MPP_BROADCAST_
#define MPP_BROADCAST_ mpp_broadcast_int8
#undef MPP_BROADCAST_SCALAR_
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_int8_scalar
#undef MPP_BROADCAST_2D_
#define MPP_BROADCAST_2D_ mpp_broadcast_int8_2d
#undef MPP_BROADCAST_3D_
#define MPP_BROADCAST_3D_ mpp_broadcast_int8_3d
#undef MPP_BROADCAST_4D_
#define MPP_BROADCAST_4D_ mpp_broadcast_int8_4d
#undef MPP_BROADCAST_5D_
#define MPP_BROADCAST_5D_ mpp_broadcast_int8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER8
#include <mpp_transmit_mpi.h>
#endif

#undef MPP_TRANSMIT_
#define MPP_TRANSMIT_ mpp_transmit_int4
#undef MPP_TRANSMIT_SCALAR_
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_int4_scalar
#undef MPP_TRANSMIT_2D_
#define MPP_TRANSMIT_2D_ mpp_transmit_int4_2d
#undef MPP_TRANSMIT_3D_
#define MPP_TRANSMIT_3D_ mpp_transmit_int4_3d
#undef MPP_TRANSMIT_4D_
#define MPP_TRANSMIT_4D_ mpp_transmit_int4_4d
#undef MPP_TRANSMIT_5D_
#define MPP_TRANSMIT_5D_ mpp_transmit_int4_5d
#undef MPP_RECV_
#define MPP_RECV_ mpp_recv_int4
#undef MPP_RECV_SCALAR_
#define MPP_RECV_SCALAR_ mpp_recv_int4_scalar
#undef MPP_RECV_2D_
#define MPP_RECV_2D_ mpp_recv_int4_2d
#undef MPP_RECV_3D_
#define MPP_RECV_3D_ mpp_recv_int4_3d
#undef MPP_RECV_4D_
#define MPP_RECV_4D_ mpp_recv_int4_4d
#undef MPP_RECV_5D_
#define MPP_RECV_5D_ mpp_recv_int4_5d
#undef MPP_SEND_
#define MPP_SEND_ mpp_send_int4
#undef MPP_SEND_SCALAR_
#define MPP_SEND_SCALAR_ mpp_send_int4_scalar
#undef MPP_SEND_2D_
#define MPP_SEND_2D_ mpp_send_int4_2d
#undef MPP_SEND_3D_
#define MPP_SEND_3D_ mpp_send_int4_3d
#undef MPP_SEND_4D_
#define MPP_SEND_4D_ mpp_send_int4_4d
#undef MPP_SEND_5D_
#define MPP_SEND_5D_ mpp_send_int4_5d
#undef MPP_BROADCAST_
#define MPP_BROADCAST_ mpp_broadcast_int4
#undef MPP_BROADCAST_SCALAR_
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_int4_scalar
#undef MPP_BROADCAST_2D_
#define MPP_BROADCAST_2D_ mpp_broadcast_int4_2d
#undef MPP_BROADCAST_3D_
#define MPP_BROADCAST_3D_ mpp_broadcast_int4_3d
#undef MPP_BROADCAST_4D_
#define MPP_BROADCAST_4D_ mpp_broadcast_int4_4d
#undef MPP_BROADCAST_5D_
#define MPP_BROADCAST_5D_ mpp_broadcast_int4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER4
#include <mpp_transmit_mpi.h>

#ifndef no_8byte_integers
#undef MPP_TRANSMIT_
#define MPP_TRANSMIT_ mpp_transmit_logical8
#undef MPP_TRANSMIT_SCALAR_
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_logical8_scalar
#undef MPP_TRANSMIT_2D_
#define MPP_TRANSMIT_2D_ mpp_transmit_logical8_2d
#undef MPP_TRANSMIT_3D_
#define MPP_TRANSMIT_3D_ mpp_transmit_logical8_3d
#undef MPP_TRANSMIT_4D_
#define MPP_TRANSMIT_4D_ mpp_transmit_logical8_4d
#undef MPP_TRANSMIT_5D_
#define MPP_TRANSMIT_5D_ mpp_transmit_logical8_5d
#undef MPP_RECV_
#define MPP_RECV_ mpp_recv_logical8
#undef MPP_RECV_SCALAR_
#define MPP_RECV_SCALAR_ mpp_recv_logical8_scalar
#undef MPP_RECV_2D_
#define MPP_RECV_2D_ mpp_recv_logical8_2d
#undef MPP_RECV_3D_
#define MPP_RECV_3D_ mpp_recv_logical8_3d
#undef MPP_RECV_4D_
#define MPP_RECV_4D_ mpp_recv_logical8_4d
#undef MPP_RECV_5D_
#define MPP_RECV_5D_ mpp_recv_logical8_5d
#undef MPP_SEND_
#define MPP_SEND_ mpp_send_logical8
#undef MPP_SEND_SCALAR_
#define MPP_SEND_SCALAR_ mpp_send_logical8_scalar
#undef MPP_SEND_2D_
#define MPP_SEND_2D_ mpp_send_logical8_2d
#undef MPP_SEND_3D_
#define MPP_SEND_3D_ mpp_send_logical8_3d
#undef MPP_SEND_4D_
#define MPP_SEND_4D_ mpp_send_logical8_4d
#undef MPP_SEND_5D_
#define MPP_SEND_5D_ mpp_send_logical8_5d
#undef MPP_BROADCAST_
#define MPP_BROADCAST_ mpp_broadcast_logical8
#undef MPP_BROADCAST_SCALAR_
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_logical8_scalar
#undef MPP_BROADCAST_2D_
#define MPP_BROADCAST_2D_ mpp_broadcast_logical8_2d
#undef MPP_BROADCAST_3D_
#define MPP_BROADCAST_3D_ mpp_broadcast_logical8_3d
#undef MPP_BROADCAST_4D_
#define MPP_BROADCAST_4D_ mpp_broadcast_logical8_4d
#undef MPP_BROADCAST_5D_
#define MPP_BROADCAST_5D_ mpp_broadcast_logical8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(LONG_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER8
#include <mpp_transmit_mpi.h>
#endif

#undef MPP_TRANSMIT_
#define MPP_TRANSMIT_ mpp_transmit_logical4
#undef MPP_TRANSMIT_SCALAR_
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_logical4_scalar
#undef MPP_TRANSMIT_2D_
#define MPP_TRANSMIT_2D_ mpp_transmit_logical4_2d
#undef MPP_TRANSMIT_3D_
#define MPP_TRANSMIT_3D_ mpp_transmit_logical4_3d
#undef MPP_TRANSMIT_4D_
#define MPP_TRANSMIT_4D_ mpp_transmit_logical4_4d
#undef MPP_TRANSMIT_5D_
#define MPP_TRANSMIT_5D_ mpp_transmit_logical4_5d
#undef MPP_RECV_
#define MPP_RECV_ mpp_recv_logical4
#undef MPP_RECV_SCALAR_
#define MPP_RECV_SCALAR_ mpp_recv_logical4_scalar
#undef MPP_RECV_2D_
#define MPP_RECV_2D_ mpp_recv_logical4_2d
#undef MPP_RECV_3D_
#define MPP_RECV_3D_ mpp_recv_logical4_3d
#undef MPP_RECV_4D_
#define MPP_RECV_4D_ mpp_recv_logical4_4d
#undef MPP_RECV_5D_
#define MPP_RECV_5D_ mpp_recv_logical4_5d
#undef MPP_SEND_
#define MPP_SEND_ mpp_send_logical4
#undef MPP_SEND_SCALAR_
#define MPP_SEND_SCALAR_ mpp_send_logical4_scalar
#undef MPP_SEND_2D_
#define MPP_SEND_2D_ mpp_send_logical4_2d
#undef MPP_SEND_3D_
#define MPP_SEND_3D_ mpp_send_logical4_3d
#undef MPP_SEND_4D_
#define MPP_SEND_4D_ mpp_send_logical4_4d
#undef MPP_SEND_5D_
#define MPP_SEND_5D_ mpp_send_logical4_5d
#undef MPP_BROADCAST_
#define MPP_BROADCAST_ mpp_broadcast_logical4
#undef MPP_BROADCAST_SCALAR_
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_logical4_scalar
#undef MPP_BROADCAST_2D_
#define MPP_BROADCAST_2D_ mpp_broadcast_logical4_2d
#undef MPP_BROADCAST_3D_
#define MPP_BROADCAST_3D_ mpp_broadcast_logical4_3d
#undef MPP_BROADCAST_4D_
#define MPP_BROADCAST_4D_ mpp_broadcast_logical4_4d
#undef MPP_BROADCAST_5D_
#define MPP_BROADCAST_5D_ mpp_broadcast_logical4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(INT_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER4
#include <mpp_transmit_mpi.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max_real8
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_REAL8
#undef MPI_REDUCE_
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce_mpi.h>

#ifdef OVERLOAD_R4
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max_real4
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_REAL4
#undef MPI_REDUCE_
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce_mpi.h>
#endif

#ifndef no_8byte_integers
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max_int8
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER8
#undef MPI_REDUCE_
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce_mpi.h>
#endif

#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max_int4
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER4
#undef MPI_REDUCE_
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce_mpi.h>

#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min_real8
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_REAL8
#undef MPI_REDUCE_
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce_mpi.h>

#ifdef OVERLOAD_R4
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min_real4
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_REAL4
#undef MPI_REDUCE_
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce_mpi.h>
#endif

#ifndef no_8byte_integers
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min_int8
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER8
#undef MPI_REDUCE_
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce_mpi.h>
#endif       
             
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min_int4
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER4
#undef MPI_REDUCE_
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce_mpi.h>

#undef MPP_SUM_
#define MPP_SUM_ mpp_sum_real8
#undef MPP_SUM_SCALAR_
#define MPP_SUM_SCALAR_ mpp_sum_real8_scalar
#undef MPP_SUM_2D_
#define MPP_SUM_2D_ mpp_sum_real8_2d
#undef MPP_SUM_3D_
#define MPP_SUM_3D_ mpp_sum_real8_3d
#undef MPP_SUM_4D_
#define MPP_SUM_4D_ mpp_sum_real8_4d
#undef MPP_SUM_5D_
#define MPP_SUM_5D_ mpp_sum_real8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_REAL8
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum_mpi.h>

#ifdef OVERLOAD_C8
#undef MPP_SUM_
#define MPP_SUM_ mpp_sum_cmplx8
#undef MPP_SUM_SCALAR_
#define MPP_SUM_SCALAR_ mpp_sum_cmplx8_scalar
#undef MPP_SUM_2D_
#define MPP_SUM_2D_ mpp_sum_cmplx8_2d
#undef MPP_SUM_3D_
#define MPP_SUM_3D_ mpp_sum_cmplx8_3d
#undef MPP_SUM_4D_
#define MPP_SUM_4D_ mpp_sum_cmplx8_4d
#undef MPP_SUM_5D_
#define MPP_SUM_5D_ mpp_sum_cmplx8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_DOUBLE_COMPLEX
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 16
#include <mpp_sum_mpi.h>
#endif

#ifdef OVERLOAD_R4
#undef MPP_SUM_
#define MPP_SUM_ mpp_sum_real4
#undef MPP_SUM_SCALAR_
#define MPP_SUM_SCALAR_ mpp_sum_real4_scalar
#undef MPP_SUM_2D_
#define MPP_SUM_2D_ mpp_sum_real4_2d
#undef MPP_SUM_3D_
#define MPP_SUM_3D_ mpp_sum_real4_3d
#undef MPP_SUM_4D_
#define MPP_SUM_4D_ mpp_sum_real4_4d
#undef MPP_SUM_5D_
#define MPP_SUM_5D_ mpp_sum_real4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_REAL4
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#include <mpp_sum_mpi.h>
#endif

#ifdef OVERLOAD_C4
#undef MPP_SUM_
#define MPP_SUM_ mpp_sum_cmplx4
#undef MPP_SUM_SCALAR_
#define MPP_SUM_SCALAR_ mpp_sum_cmplx4_scalar
#undef MPP_SUM_2D_
#define MPP_SUM_2D_ mpp_sum_cmplx4_2d
#undef MPP_SUM_3D_
#define MPP_SUM_3D_ mpp_sum_cmplx4_3d
#undef MPP_SUM_4D_
#define MPP_SUM_4D_ mpp_sum_cmplx4_4d
#undef MPP_SUM_5D_
#define MPP_SUM_5D_ mpp_sum_cmplx4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_COMPLEX
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum_mpi.h>
#endif

#ifndef no_8byte_integers
#undef MPP_SUM_
#define MPP_SUM_ mpp_sum_int8
#undef MPP_SUM_SCALAR_
#define MPP_SUM_SCALAR_ mpp_sum_int8_scalar
#undef MPP_SUM_2D_
#define MPP_SUM_2D_ mpp_sum_int8_2d
#undef MPP_SUM_3D_
#define MPP_SUM_3D_ mpp_sum_int8_3d
#undef MPP_SUM_4D_
#define MPP_SUM_4D_ mpp_sum_int8_4d
#undef MPP_SUM_5D_
#define MPP_SUM_5D_ mpp_sum_int8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER8
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum_mpi.h>
#endif

#undef MPP_SUM_
#define MPP_SUM_ mpp_sum_int4
#undef MPP_SUM_SCALAR_
#define MPP_SUM_SCALAR_ mpp_sum_int4_scalar
#undef MPP_SUM_2D_
#define MPP_SUM_2D_ mpp_sum_int4_2d
#undef MPP_SUM_3D_
#define MPP_SUM_3D_ mpp_sum_int4_3d
#undef MPP_SUM_4D_
#define MPP_SUM_4D_ mpp_sum_int4_4d
#undef MPP_SUM_5D_
#define MPP_SUM_5D_ mpp_sum_int4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPI_TYPE_
#define MPI_TYPE_ MPI_INTEGER4
#undef MPP_TYPE_BYTELEN_
#define MPP_TYPE_BYTELEN_ 4
#include <mpp_sum_mpi.h>
