    subroutine MPP_SUM_( a, length, pelist )
!sums array a over the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast: all PEs have the sum in a at the end
  !we are using f77-style call: array passed by address and not descriptor; further, 
  !the f90 conformance check is avoided.
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_, intent(inout) :: a(*)
      integer :: n
!first <length> words are array, rest are pWrk
      MPP_TYPE_ :: work(length+length/2+1+SHMEM_REDUCE_MIN_WRKDATA_SIZE)
      pointer( ptr, work )
      integer :: words
      character(len=8) :: text

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
!allocate space from the stack for pwrk and b
      ptr = LOC(mpp_stack)
      words = size(work(:))*size(transfer(work(1),word))
      if( words.GT.mpp_stack_size )then
          write( text, '(i8)' )words
          call mpp_error( FATAL, 'MPP_SUM user stack overflow: call mpp_set_stack_size('//text//') from all PEs.' )
      end if
      mpp_stack_hwm = max( words, mpp_stack_hwm )
      work(1:length) = a(1:length)
      call mpp_sync(pelist)
      call SHMEM_SUM_( work, work, length, peset(n)%start, peset(n)%log2stride, peset(n)%count, work(length+1), sync )
      a(1:length) = work(1:length)
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, length*MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_SUM_

!#######################################################################
#include <mpp_sum.inc>
