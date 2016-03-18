
! -*-f90-*-
! $Id: mpp_comm_nocomm.inc,v 18.0 2010/03/02 23:56:55 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    subroutine mpp_init( flags, in, out, err, log )
  !      integer, optional, intent(in) :: flags, in, out, err, log
subroutine mpp_init( flags,localcomm )
  integer, optional, intent(in) :: flags
  integer, optional, intent(in) :: localcomm ! dummy here, used only in MPI
  integer                       :: my_pe, num_pes, len, i, logunit
  logical                       :: opened, existed

  if( module_is_initialized )return

  module_is_initialized = .TRUE.

  !PEsets: make defaults illegal
  peset(:)%count = -1
  peset(:)%id = -1
  peset(:)%group = -1
  !0=single-PE, initialized so that count returns 1
  peset(0)%count = 1
  allocate( peset(0)%list(1) )
  peset(0)%list = pe
  current_peset_num = 0
  peset(0)%id = 0
  world_peset_num = 0
  current_peset_num = world_peset_num !initialize current PEset to world

  !initialize clocks
  call SYSTEM_CLOCK( count=tick0, count_rate=ticks_per_sec, count_max=max_ticks )
  tick_rate = 1./ticks_per_sec
  clock0 = mpp_clock_id( 'Total runtime', flags=MPP_CLOCK_SYNC )

  if( PRESENT(flags) )then
     debug   = flags.EQ.MPP_DEBUG
     verbose = flags.EQ.MPP_VERBOSE .OR. debug
  end if

! non-root pe messages written to other location than stdout()
! 9 is reserved for etc_unit
  etc_unit=9
  inquire(unit=etc_unit,opened=opened)
  if(opened) call mpp_error(FATAL,'Unit 9 is already in use (etc_unit) in mpp_comm_nocomm')
  if (trim(etcfile) /= '/dev/null') then
    write( etcfile,'(a,i6.6)' )trim(etcfile)//'.', pe
  endif
  inquire(file=etcfile, exist=existed)
  if(existed) then
     open( unit=etc_unit, file=trim(etcfile), status='REPLACE' )
  else
     open( unit=etc_unit, file=trim(etcfile) )
  endif
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

  call mpp_init_logfile()

  !messages
  if( verbose )call mpp_error( NOTE, 'MPP_INIT: initializing MPP module...' )
  if( pe.EQ.root_pe )then
     logunit = stdlog()
     write( logunit,'(/a)' )'MPP module '//trim(version)//trim(tagname)
     write( logunit,'(a,i6)' )'MPP started with NPES=', npes
     write( logunit,'(a)' )'Using no library for message passing...'
     write( logunit, '(a26,es12.4,a6,i10,a11)' ) &
          'Realtime clock resolution=', tick_rate, ' sec (', ticks_per_sec, ' ticks/sec)'
     write( logunit, '(a23,es12.4,a6,i20,a7)' ) &
          'Clock rolls over after ', max_ticks*tick_rate, ' sec (', max_ticks, ' ticks)'
  end if
  call mpp_clock_begin(clock0)

  return
end subroutine mpp_init

!#######################################################################
  !to be called at the end of a run
subroutine mpp_exit()
  integer :: i, j, k, n, nmax, istat, out_unit
  real    :: t, tmin, tmax, tavg, tstd
  real    :: m, mmin, mmax, mavg, mstd, t_total
  logical :: opened

  if( .NOT.module_is_initialized )return
  call mpp_set_current_pelist()
  call mpp_clock_end(clock0)
  t_total = clocks(clock0)%total_ticks*tick_rate
  out_unit = stdout()
  if( clock_num.GT.0 )then
     if( ANY(clocks(1:clock_num)%detailed) )then
        call sum_clock_data; call dump_clock_summary
     end if
     if( pe.EQ.root_pe )then
        write( out_unit,'(/a,i6,a)' ) 'Tabulating mpp_clock statistics across ', npes, ' PEs...'
        if( ANY(clocks(1:clock_num)%detailed) ) &
             write( out_unit,'(a)' )'   ... see mpp_clock.out.#### for details on individual PEs.'
        write( out_unit,'(/32x,a)' ) '          tmin          tmax          tavg          tstd  tfrac grain pemin pemax'
     else
        write( out_unit,'(/37x,a)' ) 'time'
     end if
     call FLUSH( out_unit )
     call mpp_sync()
     do i = 1,clock_num
        if( .NOT.ANY(peset(clocks(i)%peset_num)%list(:).EQ.pe) )cycle
        call mpp_set_current_pelist( peset(clocks(i)%peset_num)%list )
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
        if (pe.NE.root_pe) write(out_unit,'(a32,f14.6)') clocks(i)%name, clocks(i)%total_ticks*tick_rate
     end do
     if( ANY(clocks(1:clock_num)%detailed) .AND. pe.EQ.root_pe )write( out_unit,'(/32x,a)' ) &
          '       tmin       tmax       tavg       tstd       mmin       mmax       mavg       mstd  mavg/tavg'

     do i = 1,clock_num
        !messages: bytelengths and times
        if( .NOT.clocks(i)%detailed )cycle
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

  return
end subroutine mpp_exit

!#######################################################################
 !--- dummy routine, should never been called -------------------------
    subroutine mpp_malloc( ptr, newlen, len )
      integer, intent(in)    :: newlen
      integer, intent(inout) :: len
      real                   :: dummy
      pointer( ptr, dummy )

      call mpp_error(FATAL,'mpp_malloc: Should not been called when lib_SMA is not used')

      return
    end subroutine mpp_malloc

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
#include <mpp_transmit_nocomm.h>

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
#include <mpp_transmit_nocomm.h>
#endif

#ifdef OVERLOAD_R4
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
#include <mpp_transmit_nocomm.h>
#endif

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
#include <mpp_transmit_nocomm.h>
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
#include <mpp_transmit_nocomm.h>
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
#include <mpp_transmit_nocomm.h>

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
#include <mpp_transmit_nocomm.h>
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
#include <mpp_transmit_nocomm.h>

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
#include <mpp_reduce_nocomm.h>

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
#include <mpp_reduce_nocomm.h>
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
#include <mpp_reduce_nocomm.h>
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
#include <mpp_reduce_nocomm.h>

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
#include <mpp_reduce_nocomm.h>

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
#include <mpp_reduce_nocomm.h>
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
#include <mpp_reduce_nocomm.h>
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
#include <mpp_reduce_nocomm.h>

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
#include <mpp_sum_nocomm.h>

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
#include <mpp_sum_nocomm.h>
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
#include <mpp_sum_nocomm.h>
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
#include <mpp_sum_nocomm.h>
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
#include <mpp_sum_nocomm.h>
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
#include <mpp_sum_nocomm.h>
