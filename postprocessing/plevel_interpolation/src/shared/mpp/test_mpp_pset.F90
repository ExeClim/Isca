#ifdef test_mpp_pset
#include <fms_platform.h>
program test
  use mpp_mod, only: mpp_init, mpp_exit, mpp_pe, mpp_npes, stderr, stdout, &
       mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_pset_mod, only: mpp_pset_type, mpp_pset_create, mpp_pset_root, &
       mpp_pset_broadcast_ptr, mpp_pset_segment_array, mpp_pset_sync, &
       mpp_pset_stack_push, mpp_pset_print_chksum, mpp_pset_delete
  implicit none
!test program demonstrates how to create PSETs
!  how to distribute allocatable arrays
!  how to distribute automatic arrays
  integer, parameter :: n=96 !divisible by lots of numbers
  real, allocatable, dimension(:,:,:) :: a, b, cc
  real :: c(n,n,n)
#ifdef use_CRI_pointers
  pointer( ptr_c, c )
#endif
  integer(POINTER_KIND) :: ptr !useless declaration, but it will compile
  integer :: i, j, k, ks, ke
!MPP
  integer :: pe, npes
!MPP_PSET
  type(mpp_pset_type) :: pset
  logical :: root
!clocks
  integer :: id_full, id_alloc, id_auto

  call mpp_init()
  pe = mpp_pe()
  npes = mpp_npes()
  write( stdout(),'(a,i6)' )'Starting MPP_PSET unit test, npes=', npes
  call mpp_pset_create( npes, pset )
  root = mpp_pset_root(pset)
  id_full = mpp_clock_id( 'Full array' )
  id_alloc = mpp_clock_id( 'Allocatable array, PSETs' )
  id_auto = mpp_clock_id( 'Automatic array, PSETs' )
!allocate a and b
  allocate( a(n,n,n) )
  allocate( b(n,n,n) )
!allocate shared array c
  if( root )then
      allocate( cc(n,n,n) )
#ifdef use_CRI_pointers
      ptr = LOC(cc)
#endif
  end if
  call mpp_pset_broadcast_ptr( pset, ptr )
#ifdef use_CRI_pointers
  ptr_c = ptr
#endif
!initialize a and b
  call RANDOM_NUMBER(a)
  call mpp_clock_begin(id_full)
  do k = 1,n
     do j = 1,n
        do i = 1,n
           b(i,j,k) = 2*a(i,j,k)
        end do
     end do
  end do
  call mpp_clock_end(id_full)
!divide up among PSETs
  call mpp_pset_segment_array( pset, 1, n, ks, ke )
  write( stderr(),'(a,4i6)' )'pe, n, ks, ke=', pe, n, ks, ke
  call mpp_clock_begin(id_alloc)
  do k = ks,ke
     do j = 1,n
        do i = 1,n
           c(i,j,k) = 2*a(i,j,k)
        end do
     end do
  end do
  call mpp_pset_sync(pset)
  call mpp_clock_end(id_alloc)
  write( stderr(),'(a,i6,2es23.15)' )'b, c should be equal: pe b c=', &
       pe, sum(b), sum(c)
  call mpp_pset_print_chksum( pset, 'test_alloc', c(:,:,ks:ke) )
  call test_auto(n)
  call mpp_pset_delete(pset)
  call mpp_exit()

contains

  subroutine test_auto(m)
!same test as above, on auto array d
!this is how you create shared auto arrays
    integer, intent(in) :: m
    real :: d(m,m,m)
    integer :: js, je
#ifdef use_CRI_pointers
    pointer( pd, d )
    call mpp_pset_stack_push( pset, pd, size(d) )
#endif
    call mpp_pset_segment_array( pset, 1, m, js, je )
    call mpp_clock_begin(id_auto)
    do k = 1,m
       do j = js,je
          do i = 1,m
             d(i,j,k) = 2*a(i,j,k)
          end do
       end do
    end do
    call mpp_pset_sync(pset)
    call mpp_clock_end(id_auto)
    write( stderr(),'(a,i6,2es23.15)' )'b, d should be equal: pe b d=', &
         pe, sum(b), sum(d)
    call mpp_pset_print_chksum( pset, 'test_auto ', d(:,js:je,:) )
  end subroutine test_auto
    
end program test
#else
module null_mpp_pset_test
end module
#endif
