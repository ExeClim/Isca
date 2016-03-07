    subroutine MPP_REDUCE_( a, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      MPP_TYPE_ :: work

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
      if( verbose )call mpp_error( NOTE, 'MPP_REDUCE_: using MPI_ALLREDUCE...' )
      call MPI_ALLREDUCE( a, work, 1, MPI_TYPE_, MPI_REDUCE_, peset(n)%id, error )
      a = work
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_REDUCE_
