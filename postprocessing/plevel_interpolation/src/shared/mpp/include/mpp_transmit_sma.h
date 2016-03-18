!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_TRANSMIT_( put_data, put_len, to_pe, get_data, get_len, from_pe, block )
!a message-passing routine intended to be reminiscent equally of both MPI and SHMEM

!put_data and get_data are contiguous MPP_TYPE_ arrays

!at each call, your put_data array is put to   to_pe's get_data
!              your get_data array is got from from_pe's put_data
!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get

!special PE designations:
!      NULL_PE: to disable a put or a get (e.g at boundaries)
!      ANY_PE:  if remote PE for the put or get is to be unspecific
!      ALL_PES: broadcast and collect operations (collect not yet implemented)

!ideally we would not pass length, but this f77-style call performs better (arrays passed by address, not descriptor)
!further, this permits <length> contiguous words from an array of any rank to be passed (avoiding f90 rank conformance check)

!caller is responsible for completion checks (mpp_sync_self) before and after

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(*)
      MPP_TYPE_, intent(out) :: get_data(*)
      logical, intent(in), optional :: block
#include <mpp/shmem.fh>
      external shmem_ptr

      integer :: i, outunit
      integer :: np
      integer(LONG_KIND) :: data_loc
!pointer to remote data
      MPP_TYPE_ :: remote_data(get_len)
      pointer( ptr_remote_data, remote_data )
      MPP_TYPE_ :: broadcast_data(get_len)
      pointer( ptr, broadcast_data )
      integer :: words
      character(len=8) :: text
      MPP_TYPE_, allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return
      
      outunit = stdout()
      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
!send data pointer to to_pe
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INT8_WAIT( status(to_pe), MPP_WAIT )
          status(to_pe) = MPP_WAIT !prohibit puts to to_pe until it has retrieved this message
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
#ifdef __ia64
          data_loc = shmem_ptr(put_data,pe)
!          write(0,*)'pe, data_loc, loc(put_data)=', pe, data_loc, loc(put_data)
#else
          data_loc = LOC(put_data)
#endif
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INTEGER_PUT( mpp_from_pe, pe, 1, to_pe )
          call SHMEM_PUT8( remote_data_loc(pe), data_loc, 1, to_pe )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_SEND, put_len*MPP_TYPE_BYTELEN_ )
      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len )call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
              end if
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets
          if( from_pe.LT.0 .OR. from_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe along with to_pe=ANY_PE.' )
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_GET_( get_data, put_data, get_len, from_pe )
          call SHMEM_PUT8( status(pe), MPP_READY, 1, from_pe ) !tell from_pe that you have retrieved this message
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
          return
      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get: for libSMA, a get means do a wait to ensure put on remote PE is complete
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( debug )write( stderr(),* )'pe, from_pe, remote_data_loc(from_pe)=', pe, from_pe, remote_data_loc(from_pe)
          call SHMEM_INT8_WAIT( remote_data_loc(from_pe), MPP_WAIT )
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
          ptr_remote_data = remote_data_loc(from_pe)
          remote_data_loc(from_pe) = MPP_WAIT !reset
!          call SHMEM_PUT8( status(pe), MPP_READY, 1, from_pe ) !tell from_pe we have retrieved the location
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#if defined(CRAYPVP) || defined(sgi_mipspro) || defined(__ia64)
!since we have the pointer to remote data, just retrieve it with a simple copy
          if( LOC(get_data).NE.LOC(remote_data) )then
!dir$ IVDEP
              do i = 1,get_len
                 get_data(i) = remote_data(i)
              end do
          else
              call mpp_error(FATAL)
          end if
#else
          call SHMEM_GET_( get_data, remote_data, get_len, from_pe )
#endif
          call SHMEM_PUT8( status(pe), MPP_READY, 1, from_pe ) !tell from_pe we have retrieved the location
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )

      else if( from_pe.EQ.ANY_PE )then
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
!since we don't know which PE is sending us data, we wait for remote PE to send us its ID
!this is only required for !CRAYPVP  && !sgi_mipspro, but is done there too, so that we can send put_is_done back.
          call shmem_integer_wait( mpp_from_pe, ANY_PE )
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INT8_WAIT( remote_data_loc(mpp_from_pe), MPP_WAIT )
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
          ptr_remote_data = remote_data_loc(mpp_from_pe)
          remote_data_loc(mpp_from_pe) = MPP_WAIT !reset
          call SHMEM_PUT8( status(pe), MPP_READY, 1, mpp_from_pe ) !tell mpp_from_pe we have retrieved the location
#if defined(CRAYPVP) || defined(sgi_mipspro) || defined(__ia64)
!since we have the pointer to remote data, just retrieve it with a simple copy
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( LOC(get_data).NE.LOC(remote_data) )then
!dir$ IVDEP
              do i = 1,get_len
                 get_data(i) = remote_data(i)
              end do
          end if
#else
          call SHMEM_GET_( get_data, remote_data, get_len, mpp_from_pe )
#endif
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
          mpp_from_pe = ANY_PE   !reset
      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning, and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if
      return
    end subroutine MPP_TRANSMIT_


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_BROADCAST_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n
      integer :: np, i, outunit
      integer(LONG_KIND) :: data_loc
!pointer to remote data
      MPP_TYPE_ :: bdata(length)
      pointer( ptr, bdata )
      integer :: words
      character(len=8) :: text

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( debug )then
          call SYSTEM_CLOCK(tick)
          outunit = stdout()
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_BROADCAST begin: from_pe, length=', from_pe, length
      end if

      if( .NOT.ANY(from_pe.EQ.peset(current_peset_num)%list) ) &
           call mpp_error( FATAL, 'MPP_BROADCAST: broadcasting from invalid PE.' )

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
      ptr = LOC(mpp_stack)
      words = size(bdata(:))*size(transfer(bdata(1),word))
      if( words.GT.mpp_stack_size )then
          write( text, '(i8)' )words
          call mpp_error( FATAL, 'MPP_BROADCAST user stack overflow: call mpp_set_stack_size('//text//') from all PEs.' )
      end if
      mpp_stack_hwm = max( words, mpp_stack_hwm )
      if( mpp_npes().GT.1 )then
!dir$ IVDEP
          do i = 1,length
             bdata(i) = data(i)
          end do
          call mpp_sync(pelist) !eliminate?
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          call SHMEM_BROADCAST_( bdata, bdata, length, from_pe, peset(n)%start, peset(n)%log2stride, peset(n)%count, sync )
          call mpp_sync(pelist) !eliminate?
!dir$ IVDEP
          do i = 1,length
             data(i) = bdata(i)
          end do
      end if
      if( current_clock.NE.0 )call increment_current_clock( EVENT_BROADCAST, length*MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_BROADCAST_

!####################################################################################
#include <mpp_transmit.inc>
