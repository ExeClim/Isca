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
      logical                       :: block_comm
      integer                       :: i, out_unit
      MPP_TYPE_, allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return
      
      block_comm = .true.
      if(PRESENT(block)) block_comm = block

      out_unit = stdout()
      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( out_unit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
!use non-blocking sends
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( request(to_pe).NE.MPI_REQUEST_NULL )then !only one message from pe->to_pe in queue 
!              if( debug )write( stderr(),* )'PE waiting for sending', pe, to_pe
              call MPI_WAIT( request(to_pe), stat, error )
          end if
          call MPI_ISEND( put_data, put_len, MPI_TYPE_, to_pe, tag, mpp_comm_private, request(to_pe), error )
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
!...but you cannot have a pure get with MPI
          call mpp_error( FATAL, 'MPP_TRANSMIT: you cannot transmit to ANY_PE using MPI.' )
      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get: for libSMA, a get means do a wait to ensure put on remote PE is complete
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
!receive from from_pe
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( block_comm ) then
             call MPI_RECV( get_data, get_len, MPI_TYPE_, from_pe, MPI_ANY_TAG, mpp_comm_private, stat, error )
          else
             if( request_recv(from_pe).NE.MPI_REQUEST_NULL )then !only one message from from_pe->pe in queue 
                !              if( debug )write( stderr(),* )'PE waiting for receiving', pe, from_pe
                call MPI_WAIT( request_recv(from_pe), stat, error )
             end if
             call MPI_IRECV( get_data, get_len, MPI_TYPE_, from_pe, MPI_ANY_TAG, mpp_comm_private, &
                  request_recv(from_pe), error ) 
          endif
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
      else if( from_pe.EQ.ANY_PE )then
!receive from MPI_ANY_SOURCE
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call MPI_RECV( get_data, get_len, MPI_TYPE_, MPI_ANY_SOURCE, MPI_ANY_TAG, mpp_comm_private, stat, error )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning, and hence is not implemented.' )
      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( out_unit,'(a,i18,a,i5,a,2i5,2i8)' )&
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
      integer :: n, i, from_rank, out_unit

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      out_unit = stdout()
      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( out_unit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_BROADCAST begin: from_pe, length=', from_pe, length
      end if

      if( .NOT.ANY(from_pe.EQ.peset(current_peset_num)%list) ) &
           call mpp_error( FATAL, 'MPP_BROADCAST: broadcasting from invalid PE.' )

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
 ! find the rank of from_pe in the pelist.     
      do i = 1, mpp_npes()
         if(peset(n)%list(i) == from_pe) then
             from_rank = i - 1
             exit
         endif
      enddo
      if( mpp_npes().GT.1 )call MPI_BCAST( data, length, MPI_TYPE_, from_rank, peset(n)%id, error )
      if( current_clock.NE.0 )call increment_current_clock( EVENT_BROADCAST, length*MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_BROADCAST_

!####################################################################################
#include <mpp_transmit.inc>
