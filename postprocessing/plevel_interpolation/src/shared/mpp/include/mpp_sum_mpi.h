    subroutine MPP_SUM_( a, length, pelist )
!sums array a over the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast: all PEs have the sum in a at the end
  !we are using f77-style call: array passed by address and not descriptor; further, 
  !the f90 conformance check is avoided.
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_, intent(inout) :: a(*)
      integer :: n
      MPP_TYPE_ :: work(length)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
      if( verbose )call mpp_error( NOTE, 'MPP_SUM: using MPI_ALLREDUCE...' )
      if( debug )write( stderr(),* )'pe, n, peset(n)%id=', pe, n, peset(n)%id
      call MPI_ALLREDUCE( a, work, length, MPI_TYPE_, MPI_SUM, peset(n)%id, error )
      a(1:length) = work(1:length)
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, length*MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_SUM_

!#######################################################################
#include <mpp_sum.inc>
