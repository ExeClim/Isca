#ifdef test_mpp
#ifdef SYSTEM_CLOCK
#undef SYSTEM_CLOCK
#endif

program test   !test various aspects of mpp_mod
#include <fms_platform.h>

#ifdef sgi_mipspro
  use shmem_interface
#endif

  use mpp_mod, only : mpp_init, mpp_exit, mpp_pe, mpp_npes, mpp_root_pe, stdout
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync, mpp_malloc
  use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist, mpp_set_stack_size
  use mpp_mod, only : mpp_broadcast, mpp_transmit, mpp_sum, mpp_max, mpp_chksum, ALL_PES
#ifdef use_MPI_GSM
  use mpp_mod, only : mpp_gsm_malloc, mpp_gsm_free
#endif

  implicit none

  integer, parameter              :: n=1048576
  real, allocatable, dimension(:) :: a, b, c
#ifdef use_MPI_GSM
  real                            :: d(n)
  pointer (locd, d)
#else
  real, allocatable, dimension(:) :: d
  integer(LONG_KIND) :: locd
#endif
  integer                         :: tick, tick0, ticks_per_sec, id
  integer                         :: pe, npes, root, i, j, k, l, m, n2, istat
  real                            :: dt

  call mpp_init()
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()
  call SYSTEM_CLOCK( count_rate=ticks_per_sec )
  allocate( a(n), b(n) )
  id = mpp_clock_id( 'Random number' )
  call mpp_clock_begin(id)
  call random_number(a)
  call mpp_clock_end  (id)
  !---------------------------------------------------------------------!
  !   time transmit, compare against shmem_put and get                  !
  !---------------------------------------------------------------------!
  if( pe.EQ.root )then
     print *, 'Time mpp_transmit for various lengths...'
#ifdef SGICRAY
     print *, 'For comparison, times for shmem_get and shmem_put are also provided.'
#endif
     print *
  end if
  id = mpp_clock_id( 'mpp_transmit' )
  call mpp_clock_begin(id)
  !timing is done for cyclical pass (more useful than ping-pong etc)
  l = n
  do while( l.GT.0 )
     !--- mpp_transmit -------------------------------------------------
     call mpp_sync()
     call SYSTEM_CLOCK(tick0)
     do i = 1,npes
        call mpp_transmit( put_data=a(1), plen=l, to_pe=modulo(pe+npes-i,npes), &
                           get_data=b(1), glen=l, from_pe=modulo(pe+i,npes) )
        !          call mpp_sync_self( (/modulo(pe+npes-i,npes)/) )
     end do
     call mpp_sync()
     call SYSTEM_CLOCK(tick)
     dt = real(tick-tick0)/(npes*ticks_per_sec)
     dt = max( dt, epsilon(dt) )
     if( pe.EQ.root )write( stdout(),'(/a,i8,f13.6,f8.2)' )'MPP_TRANSMIT length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!#ifdef SGICRAY
!     !--- shmem_put ----------------------------------------------------
!     call mpp_sync()
!     call SYSTEM_CLOCK(tick0)
!     do i = 1,npes
!       call shmem_real_put( b, a, l, modulo(pe+1,npes) )
!     end do
!     call mpp_sync()
!     call SYSTEM_CLOCK(tick)
!     dt = real(tick-tick0)/(npes*ticks_per_sec)
!     dt = max( dt, epsilon(dt) )
!     if( pe.EQ.root )write( stdout(),'( a,i8,f13.6,f8.2)' )'SHMEM_PUT    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!     !--- shmem_get ----------------------------------------------------
!     call mpp_sync()
!     call SYSTEM_CLOCK(tick0)
!     do i = 1,npes
!        call shmem_real_get( b, a, l, modulo(pe+1,npes) )
!     end do
!     call SYSTEM_CLOCK(tick)
!     dt = real(tick-tick0)/(npes*ticks_per_sec)
!     dt = max( dt, epsilon(dt) )
!     if( pe.EQ.root )write( stdout(),'( a,i8,f13.6,f8.2)' )'SHMEM_GET    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!#endif
     l = l/2
  end do
  !---------------------------------------------------------------------!
  !                   test mpp_sum                                      !
  !---------------------------------------------------------------------!
  if( pe.EQ.root )then
     print '(/a)', 'Time mpp_sum...'
  end if
  a = real(pe+1)
  call mpp_sync()
  call SYSTEM_CLOCK(tick0)
  call mpp_sum(a(1:1000),1000)
  call SYSTEM_CLOCK(tick)
  dt = real(tick-tick0)/ticks_per_sec
  dt = max( dt, epsilon(dt) )
  if( pe.EQ.root )write( stdout(),'(a,2i6,f9.1,i8,f13.6,f8.2/)' ) &
       'mpp_sum: pe, npes, sum(pe+1), length, time, bw(Mb/s)=', pe, npes, a(1), n, dt, n*8e-6/dt
  call mpp_clock_end(id)
  !---------------------------------------------------------------------!
  !                   test mpp_max                                      !
  !---------------------------------------------------------------------!
  if( pe.EQ.root )then
     print *
     print *, 'Test mpp_max...'
  end if
  a = real(pe+1)
  print *, 'pe,     pe+1 =', pe, a(1)
  call mpp_max( a(1) )
  print *, 'pe, max(pe+1)=', pe, a(1)
  !pelist check
  call mpp_sync()
  call flush(stdout(),istat)
  if( npes.GE.2 )then
     if( pe.EQ.root )print *, 'Test of pelists: bcast, sum and max using PEs 0...npes-2 (excluding last PE)'
     call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
     a = real(pe+1)
     if( pe.NE.npes-1 )call mpp_broadcast( a, n, npes-2, (/(i,i=0,npes-2)/) )
     print *, 'bcast(npes-1) from 0 to npes-2=', pe, a(1)
     a = real(pe+1)
     if( pe.NE.npes-1 )then
        call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
        id = mpp_clock_id( 'Partial mpp_sum' )
        call mpp_clock_begin(id)
        call mpp_sum( a(1:1000), 1000, (/(i,i=0,npes-2)/) )
        call mpp_clock_end  (id)
     end if
     if( pe.EQ.root )print *, 'sum(pe+1) from 0 to npes-2=', a(1)
     a = real(pe+1)
     if( pe.NE.npes-1 )call mpp_max( a(1), (/(i,i=0,npes-2)/) )
     if( pe.EQ.root )print *, 'max(pe+1) from 0 to npes-2=', a(1)
  end if
  call mpp_set_current_pelist()
  
#ifdef use_CRI_pointers
  !---------------------------------------------------------------------!
  !                   test mpp_chksum                                   !
  !---------------------------------------------------------------------!
  if( modulo(n,npes).EQ.0 )then  !only set up for even division
     n2 = 1024
     a = 0.d0
     if( pe.EQ.root )call random_number(a(1:n2))
!    if( pe.EQ.root )call random_number(a)
     call mpp_sync()
     call mpp_transmit( put_data=a(1), plen=n2, to_pe=ALL_PES, &
                        get_data=a(1), glen=n2, from_pe=root )
!    call mpp_transmit( put_data=a(1), plen=n, to_pe=ALL_PES, &
!                       get_data=a(1), glen=n, from_pe=root )
     m= n2/npes
!    m= n/npes
     allocate( c(m) )
     c = a(pe*m+1:pe*m+m)
     
     if( pe.EQ.root )then
        print *
        print *, 'Test mpp_chksum...'
        print *, 'This test shows that a whole array and a distributed array give identical checksums.'
     end if
     print *, 'chksum(a(1:1024))=', mpp_chksum(a(1:n2),(/pe/))
     print *, 'chksum(c(1:1024))=', mpp_chksum(c)
!    print *, 'chksum(a)=', mpp_chksum(a,(/pe/))
!    print *, 'chksum(c)=', mpp_chksum(c)
  end if
!test of pointer sharing
#ifdef use_MPI_GSM
      call mpp_gsm_malloc( locd, sizeof(d) )
#else
  if( pe.EQ.root )then
      allocate( d(n) )
      locd = LOC(d)
  end if
  call mpp_broadcast(locd,root)
#endif
  if( pe.EQ.root )then
      call random_number(d)
  end if
  call mpp_sync()
!  call test_shared_pointers(locd,n)

#ifdef use_MPI_GSM
  call mpp_gsm_free( locd )
#else
  if( pe.EQ.root )then
      deallocate( d )
  end if
#endif
#endif
  call mpp_exit()

contains

  subroutine test_shared_pointers(locd,n)
    integer(LONG_KIND), intent(in) :: locd
    integer :: n
    real :: dd(n)
    pointer( p, dd )

    p = locd
    print *, 'TEST_SHARED_POINTERS: pe, locd=', pe, locd
!    print *, 'TEST_SHARED_POINTERS: pe, chksum(d)=', pe, mpp_chksum(dd,(/pe/))
    print *, 'TEST_SHARED_POINTERS: pe, sum(d)=', pe, sum(dd)
    return
  end subroutine test_shared_pointers
end program test

#else
module null_mpp_test
end module  

#endif /* test_mpp */
