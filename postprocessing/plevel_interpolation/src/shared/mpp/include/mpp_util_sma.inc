
! -*-f90-*-
! $Id: mpp_util_sma.inc,v 18.0 2010/03/02 23:58:11 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!         MISCELLANEOUS UTILITIES: mpp_error                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mpp_error_basic( errortype, errormsg )
  !a very basic error handler
  !uses ABORT and FLUSH calls, may need to use cpp to rename
  integer, intent(in) :: errortype
  character(len=*), intent(in), optional :: errormsg
  character(len=256) :: text
  logical :: opened
  integer :: istat, errunit, outunit

  if( .NOT.module_is_initialized )call ABORT()

  select case( errortype )
  case(NOTE)
     text = 'NOTE'         !just FYI
  case(WARNING)
     text = 'WARNING'      !probable error
  case(FATAL)
     text = 'FATAL'        !fatal error
  case default
     text = 'WARNING: non-existent errortype (must be NOTE|WARNING|FATAL)'
  end select

  if( npes.GT.1 )write( text,'(a,i5)' )trim(text)//' from PE', pe   !this is the mpp part
  if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)

  outunit = stdout()
  errunit = stderr()
  select case( errortype )
  case(NOTE)
     write( outunit,'(a)' )trim(text)
  case default
#ifdef __SX
     write( errunit, * )trim(text)
#else
     write( errunit,'(/a/)' )trim(text)
#endif
     write( outunit,'(/a/)' )trim(text)
     if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
        call FLUSH(outunit)
#ifdef sgi_mipspro
        call TRACE_BACK_STACK_AND_PRINT()
#endif
        call ABORT() !automatically calls traceback on Cray systems
     end if
  end select

  error_state = errortype
  return
end subroutine mpp_error_basic

!#####################################################################

function get_peset(pelist)
  integer :: get_peset
  !makes a PE set out of a PE list
  !a PE list is an ordered list of PEs
  !a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI
  !if stride is non-uniform or not a power of 2, will return error (not required for MPI but enforced for uniformity)
  integer, intent(in), optional :: pelist(:)
  integer :: group
  integer :: i, n, stride, errunit
  integer, allocatable :: sorted(:)

  if( .NOT.PRESENT(pelist) )then !set it to current_peset_num
     get_peset = current_peset_num; return
  end if
  if( size(pelist(:)).EQ.1 .AND. npes.GT.1 )then    !collective ops on single PEs should return
     get_peset = 0; return
  end if
  errunit = stderr()
  !--- first make sure pelist is monotonically increasing.
  do n = 2, size(pelist(:))
     if(pelist(n) <= pelist(n-1)) call mpp_error(FATAL, "GET_PESET: pelist is not monotonically increasing")
  enddo
  allocate( sorted(size(pelist(:))) )
  sorted = pelist 
  if( debug )write( errunit,* )'pelist=', pelist

  !find if this array matches any existing peset
  do i = 1,peset_num
     if( debug )write( errunit,'(a,3i6)' )'pe, i, peset_num=', pe, i, peset_num
     if( size(sorted(:)).EQ.size(peset(i)%list(:)) )then
        if( ALL(sorted.EQ.peset(i)%list) )then
           deallocate(sorted)
           get_peset = i; return
        end if
     end if
  end do
  !not found, so create new peset
  peset_num = peset_num + 1
  if( peset_num.GE.PESET_MAX )call mpp_error( FATAL, 'GET_PESET: number of PE sets exceeds PESET_MAX.' )
  i = peset_num             !shorthand
  !create list
  allocate( peset(i)%list(size(sorted(:))) )
  peset(i)%list(:) = sorted(:)
  peset(i)%count = size(sorted(:))
  peset(i)%start = sorted(1)
  n = size(sorted(:))
  if( size(sorted(:)).GT.1 )then
     stride = sorted(2)-sorted(1)
     if( ANY(sorted(2:n)-sorted(1:n-1).NE.stride) ) &
          call mpp_error( WARNING, 'GET_PESET: pelist must have constant stride.' )
     peset(i)%log2stride = nint( log(real(stride))/log(2.) )
     if( 2**peset(i)%log2stride.NE.stride )call mpp_error( WARNING, 'GET_PESET: pelist must have power-of-2 stride.' )
  else
     peset(i)%log2stride = 0
  end if

  deallocate(sorted)
  get_peset = i

  return

end function get_peset


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!           SYNCHRONIZATION ROUTINES: mpp_sync, mpp_sync_self                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine mpp_sync( pelist)
  !synchronize PEs in list
  integer, intent(in), optional :: pelist(:)
  integer :: n

  call mpp_sync_self(pelist)

  n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

  if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
  if( n.EQ.world_peset_num )then
     call SHMEM_BARRIER_ALL() !special call is faster
  else
     call SHMEM_BARRIER( peset(n)%start, peset(n)%log2stride, peset(n)%count, sync )
  end if
  if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)

  return
end subroutine mpp_sync

!#####################################################################
! <SUBROUTINE NAME="mpp_sync_self">
!  <OVERVIEW>
!    Local synchronization.
!  </OVERVIEW>
!  <DESCRIPTION>
!   <TT>mpp_transmit</TT> is implemented as asynchronous
!   <TT>put/send</TT> and synchronous
!   <TT>get/recv</TT>. <TT>mpp_sync_self</TT> guarantees that outstanding
!   asynchronous operations from the calling PE are complete. If
!   <TT>pelist</TT> is supplied, <TT>mpp_sync_self</TT> checks only for
!   outstanding puts to the PEs in <TT>pelist</TT>.
!
!   If <TT>pelist</TT> is omitted, the context is assumed to be the
!   current pelist. This call implies synchronization across the PEs in
!   <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
!  </DESCRIPTION>
!  <IN NAME="pelist" TYPE="integer" DIM="(:)"></IN>
! </SUBROUTINE>
subroutine mpp_sync_self( pelist, check )
  !this is to check if current PE's outstanding puts are complete
  !but we can't use shmem_fence because we are actually waiting for
  !a remote PE to complete its get
  integer, intent(in), optional :: pelist(:)
  integer, intent(in), optional :: check

  integer :: i, m, n, stride

  n = get_peset(pelist)

  if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef _CRAYT90
  call SHMEM_UDCFLUSH !invalidate data cache
#endif
  do m = 1,peset(n)%count
     i = peset(n)%list(m)
     call SHMEM_INT8_WAIT( status(i), MPP_WAIT ) !wait for status.NE.MPP_WAIT
  end do
  if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
  return
end subroutine mpp_sync_self


