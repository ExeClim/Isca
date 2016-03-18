! module within MPP for handling PSETs:
! PSET: Persistent Shared-memory Execution Thread
!
! AUTHOR: V. Balaji (v.balaji@noaa.gov)
! DATE: 2006-01-15
#include <fms_platform.h>
#ifdef test_mpp_pset
!PSET_DEBUG is always turned on in the test program
#define PSET_DEBUG
#endif

module mpp_pset_mod
  use mpp_mod, only: mpp_pe, mpp_npes, mpp_root_pe, mpp_send, mpp_recv, &
       mpp_sync, mpp_error, FATAL, WARNING, stdout, stderr, mpp_chksum, &
       mpp_declare_pelist, mpp_get_current_pelist, mpp_set_current_pelist, &
       mpp_init
  implicit none
  private

!private variables
  integer :: pe
  integer :: commID !MPI communicator, copy here from pset
  logical :: verbose=.FALSE.
  logical :: module_is_initialized=.FALSE.
  character(len=256) :: text
#ifdef use_SGI_GSM
#include <mpp/shmem.fh>
  integer :: pSync(SHMEM_BARRIER_SYNC_SIZE)
  pointer( p_pSync, pSync ) !used by SHPALLOC
#endif
!generic interfaces
  interface mpp_pset_broadcast_ptr
     module procedure mpp_pset_broadcast_ptr_scalar
     module procedure mpp_pset_broadcast_ptr_array
  end interface
  interface mpp_send_ptr
     module procedure mpp_send_ptr_scalar
     module procedure mpp_send_ptr_array
  end interface
  interface mpp_recv_ptr
     module procedure mpp_recv_ptr_scalar
     module procedure mpp_recv_ptr_array
  end interface
  interface mpp_pset_print_chksum
     module procedure mpp_pset_print_chksum_1D
     module procedure mpp_pset_print_chksum_2D
     module procedure mpp_pset_print_chksum_3D
     module procedure mpp_pset_print_chksum_4D
  end interface
!public type
  type :: mpp_pset_type
     private
     sequence
     integer :: npset !number of PSETs
     integer :: next_in_pset, prev_in_pset !next and prev PE in PSET (cyclic)
     integer :: root_in_pset !PE designated to be the root within PSET
     logical :: root !true if you are the root PSET
     integer :: pos !position of current PE within pset
!stack is allocated by root
!it is then mapped to mpp_pset_stack by mpp_pset_broadcast_ptr
     real, _ALLOCATABLE :: stack(:) _NULL
     integer, _ALLOCATABLE :: pelist(:) _NULL !base PElist
     integer, _ALLOCATABLE :: root_pelist(:) _NULL !a PElist of all the roots
     integer, _ALLOCATABLE :: pset(:) _NULL !PSET IDs
     integer(POINTER_KIND) :: p_stack
     integer :: lstack, maxstack, hiWM !current stack length, max, hiWM
     integer :: commID
     character(len=32) :: name
     logical :: initialized=.FALSE.
  end type mpp_pset_type
!public types
  public :: mpp_pset_type
!public variables
!public member functions
  public :: mpp_pset_create, mpp_pset_sync, mpp_pset_broadcast, &
       mpp_pset_broadcast_ptr, mpp_pset_check_ptr, mpp_pset_segment_array, &
       mpp_pset_stack_push, mpp_pset_stack_reset, mpp_pset_print_chksum, &
       mpp_pset_delete, mpp_pset_root, mpp_pset_numroots, mpp_pset_init, &
       mpp_pset_get_root_pelist, mpp_pset_print_stack_chksum

contains
  subroutine mpp_pset_init
#ifdef use_SGI_GSM
    integer :: err
#ifdef sgi_mipspro
    character(len=8) :: value !won't be read
    integer :: lenname, lenval!won't be read
#endif
    if( module_is_initialized )return
!this part needs to be called _all_ PEs
    call SHMEM_BARRIER_ALL()
    call SHPALLOC( p_pSync, SHMEM_BARRIER_SYNC_SIZE, err, -1 )
    call SHMEM_BARRIER_ALL()
#ifdef sgi_mipspro
    call PXFGETENV( 'SMA_GLOBAL_ALLOC', 0, value, lenval, err )
    if( err.NE.0 )call mpp_error( FATAL, &
         'The environment variable SMA_GLOBAL_ALLOC must be set on Irix.' )
#endif
#endif
    module_is_initialized = .TRUE.
  end subroutine mpp_pset_init
  
  subroutine mpp_pset_create(npset,pset,stacksize,pelist, commID)
!create PSETs
!  called by all PEs in parent pelist
!  mpset must be exact divisor of npes
    integer, intent(in) :: npset !number of PSETs per set
    type(mpp_pset_type), intent(inout) :: pset
    integer, intent(in), optional :: stacksize
    integer, intent(in), optional :: pelist(:)
    integer, intent(in), optional :: commID

    integer :: npes, my_commID
    integer :: i, j, k, out_unit
    integer, allocatable :: my_pelist(:), root_pelist(:)

    call mpp_init()
    call mpp_pset_init()

#ifdef PSET_DEBUG
    verbose=.TRUE.
#endif
    out_unit = stdout()
    pe = mpp_pe()
    if(present(pelist)) then
       npes = size(pelist(:))
    else
       npes = mpp_npes()
    endif
    if( mod(npes,npset).NE.0 )then
        write( text,'(a,2i6)' ) &
             'MPP_PSET_CREATE: PSET size (npset) must divide npes exactly:'// &
             ' npset, npes=', npset, npes
        call mpp_error( FATAL, text )
    end if

    !configure out root_pelist
    allocate(my_pelist(0:npes-1) )
    allocate(root_pelist(0:npes/npset-1) )
    if(present(pelist)) then
       if(.not. present(commID)) call mpp_error(FATAL, &
         'MPP_PSET_CREATE: when pelist is present, commID should also be present')
       my_pelist = pelist
       my_commID = commID
    else
       call mpp_get_current_pelist(my_pelist, commID = my_commID)
    endif
    do i = 0,npes/npset-1
       root_pelist(i) = my_pelist(npset*i)
    enddo
    write( out_unit,'(a,i6)' )'MPP_PSET_CREATE creating PSETs... npset=', npset
    if(ANY(my_pelist == pe) ) then
    if( pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_CREATE: PSET already initialized!' )     
    pset%npset = npset
    allocate( pset%pelist(0:npes-1) )
    allocate( pset%root_pelist(0:npes/npset-1) )
    pset%commID = my_commID
    pset%pelist = my_pelist
!create the root PElist
    pset%root_pelist = root_pelist
    allocate( pset%pset(0:npset-1) )
    do i = 0,npes/npset-1
       k = npset*i
!designate the root PE, next PE, prev PE
       do j = 0,npset-1
          if( pe.EQ.pset%pelist(k+j) )then
              pset%pset(:) =  pset%pelist(k:k+npset-1)
              pset%pos = j
              pset%root_in_pset = pset%root_pelist(i)
              if( j.EQ.0 )then
                  pset%prev_in_pset = pset%pelist(k+npset-1)
              else
                  pset%prev_in_pset = pset%pelist(k+j-1)
              end if
              if( j.EQ.npset-1 )then
                  pset%next_in_pset = pset%pelist(k)
              else
                  pset%next_in_pset = pset%pelist(k+j+1)
              end if
          end if
       end do
    end do

    pset%root = pe.EQ.pset%root_in_pset

!stack
    pset%hiWM = 0 !initialize hi-water-mark
    pset%maxstack = 1000000 !default
    if( PRESENT(stacksize) )pset%maxstack = stacksize
    write( out_unit,'(a,i8)' ) &
         'MPP_PSET_CREATE: setting stacksize=', pset%maxstack
    if( pset%root )then
        allocate( pset%stack(pset%maxstack) )
#ifdef use_CRI_pointers
        pset%p_stack = LOC(pset%stack)
#endif
    end if
    pset%initialized = .TRUE. !must be called before using pset
    call mpp_pset_broadcast_ptr(pset,pset%p_stack)
    endif

    call mpp_declare_pelist(root_pelist)

    if( verbose )then
        write( stderr(),'(a,4i6)' )'MPP_PSET_CREATE: pe, root, next, prev=', &
             pe, pset%root_in_pset, pset%next_in_pset, pset%prev_in_pset
        write( stderr(),* )'PE ', pe, ' pset=', pset%pset(:)
        write( out_unit,* )'root pelist=', pset%root_pelist(:)
    end if
  end subroutine mpp_pset_create

  subroutine mpp_pset_delete(pset)
    type(mpp_pset_type), intent(inout) :: pset
    integer :: out_unit

    out_unit = stdout()
    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_DELETE: called with uninitialized PSET.' )
!deallocate arrays...
    deallocate( pset%pelist )
    deallocate( pset%root_pelist )
    deallocate( pset%pset )
    if( pset%root )deallocate( pset%stack )
    write( out_unit, '(a,i10)' ) &
         'Deleting PSETs... stack high-water-mark=', pset%hiWM
!... and set status flag
    pset%initialized = .FALSE.
  end subroutine mpp_pset_delete

  subroutine mpp_send_ptr_scalar( ptr, pe )
    integer(POINTER_KIND), intent(in) :: ptr
    integer, intent(in) :: pe

!currently only wraps mpp_send
!on some architectures, mangling might occur
    call mpp_send( ptr, pe )
  end subroutine mpp_send_ptr_scalar

  subroutine mpp_send_ptr_array( ptr, pe )
    integer(POINTER_KIND), intent(in) :: ptr(:)
    integer, intent(in) :: pe

!currently only wraps mpp_send
!on some architectures, mangling might occur
    call mpp_send( ptr, size(ptr), pe )
  end subroutine mpp_send_ptr_array

  subroutine mpp_recv_ptr_scalar( ptr, pe )
    integer(POINTER_KIND), intent(inout) :: ptr
    integer, intent(in) :: pe

    call mpp_recv( ptr, pe )
    call mpp_translate_remote_ptr( ptr, pe )
    return
  end subroutine mpp_recv_ptr_scalar

  subroutine mpp_recv_ptr_array( ptr, pe )
    integer(POINTER_KIND), intent(inout) :: ptr(:)
    integer, intent(in) :: pe
    integer :: i

    call mpp_recv( ptr, size(ptr), pe )
    do i = 1, size(ptr)
       call mpp_translate_remote_ptr( ptr(i), pe )
    end do
    return
  end subroutine mpp_recv_ptr_array

  subroutine mpp_translate_remote_ptr( ptr, pe )
!modifies the received pointer to correct numerical address
    integer(POINTER_KIND), intent(inout) :: ptr
    integer, intent(in) :: pe
#ifdef use_SGI_GSM
!from the MPI_SGI_GLOBALPTR manpage
!            POINTER(global_ptr, global_addr)
!            INTEGER rem_rank, comm, ierror
!            INTEGER(KIND=MPI_ADDRESS_KIND) rem_addr, size, global_addr
!
!            CALL MPI_SGI_GLOBALPTR(rem_addr, size, rem_rank, comm, global_ptr, ierror)
    real :: dummy
    pointer( p, dummy )
    integer :: ierror
!length goes in the second argument to MPI_SGI_GLOBALPTR
!    according to Kim Mcmahon, this is only used to ensure the requested array
!    length is within the valid memory-mapped region. We do not have access to
!    the actual array length, so we are only going to set it to 1. This might
!    unexpectedly fail on some large model.
    integer(POINTER_KIND) :: length=1
#ifdef sgi_mipspro
    return !no translation needed on sgi_mipspro if SMA_GLOBAL_ALLOC is set
#endif
#ifdef use_libMPI
!the MPI communicator was stored in pset%commID
!since this routine doesn't take a pset argument, we let the caller store
!it in the module global variable commID (see broadcast_ptr and check_ptr)
    p = ptr
    call MPI_SGI_GLOBALPTR( dummy, length, pe, commID, ptr, ierror )
    if( ierror.EQ.-1 )call mpp_error( FATAL, &
         'MPP_TRANSLATE_REMOTE_PTR: unknown MPI_SGI_GLOBALPTR error.' )
#else
    call mpp_error( FATAL, &
         'MPP_TRANSLATE_REMOTE_PTR now only works under -Duse_libMPI' )
#endif
#endif
    return
  end subroutine mpp_translate_remote_ptr

  subroutine mpp_pset_sync(pset)
!this is a replacement for mpp_sync, doing syncs across
!shared arrays without calling mpp_sync
    type(mpp_pset_type), intent(in) :: pset

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_SYNC: called with uninitialized PSET.' )
#ifdef use_SGI_GSM
!assumes npset contiguous PEs starting with root_in_pset
    call SHMEM_BARRIER( pset%root_in_pset, 0, pset%npset, pSync )
#else
!currently does mpp_sync!!! slow!!!
!try and make a lightweight pset sync
    call mpp_sync
#endif
  end subroutine mpp_pset_sync

  subroutine mpp_pset_broadcast(pset,a)
!broadcast value on the root to its sub-threads
    type(mpp_pset_type), intent(in) :: pset
    real, intent(inout) :: a
    integer :: i

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_BROADCAST: called with uninitialized PSET.' )
    if( pset%root )then
        do i = 1,pset%npset-1
           call mpp_send( a, pset%pset(i) )
        end do
    else
        call mpp_recv( a, pset%root_in_pset )
    end if
    call mpp_pset_sync(pset)
  end subroutine mpp_pset_broadcast

  subroutine mpp_pset_broadcast_ptr_scalar(pset,ptr)
!create a shared array by broadcasting pointer
!root allocates memory and passes pointer in
!on return all other PSETs will have the pointer to a shared object
    type(mpp_pset_type), intent(in) :: pset
    integer(POINTER_KIND), intent(inout) :: ptr
    integer :: i

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_BROADCAST_PTR: called with uninitialized PSET.' )
    commID = pset%commID !pass to mpp_translate_remote_ptr
    if( pset%root )then
        do i = 1,pset%npset-1
           call mpp_send_ptr( ptr, pset%pset(i) )
        end do
    else
        call mpp_recv_ptr( ptr, pset%root_in_pset )
    end if
  end subroutine mpp_pset_broadcast_ptr_scalar

  subroutine mpp_pset_broadcast_ptr_array(pset,ptr)
!create a shared array by broadcasting pointer
!root allocates memory and passes pointer in
!on return all other PSETs will have the pointer to a shared object
    type(mpp_pset_type), intent(in) :: pset
    integer(POINTER_KIND), intent(inout) :: ptr(:)
    integer :: i

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_BROADCAST_PTR: called with uninitialized PSET.' )
    commID = pset%commID !pass to mpp_translate_remote_ptr
    if( pset%root )then
        do i = 1,pset%npset-1
           call mpp_send_ptr( ptr, pset%pset(i) )
        end do
    else
        call mpp_recv_ptr( ptr, pset%root_in_pset )
    end if
  end subroutine mpp_pset_broadcast_ptr_array

  subroutine mpp_pset_check_ptr(pset,ptr)
!checks if the supplied pointer is indeed shared
    type(mpp_pset_type), intent(in) :: pset
#ifdef use_CRI_pointers
    real :: dummy
    pointer( ptr, dummy )
#else
    integer(POINTER_KIND), intent(in) :: ptr
#endif
#ifdef PSET_DEBUG
    integer(POINTER_KIND) :: p
    integer :: i
    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_CHECK_PTR: called with uninitialized PSET.' )
    commID = pset%commID !pass to mpp_translate_remote_ptr
!check if this is a shared pointer
    p = ptr
    if( pset%root )then
        do i = 1,pset%npset-1
           call mpp_send_ptr( p, pset%pset(i) )
        end do
    else
        call mpp_recv_ptr( p, pset%root_in_pset )
    end if
    call mpp_pset_sync(pset)
    if( p.NE.ptr )call mpp_error( FATAL, &
         'MPP_PSET_CHECK_PTR: pointers do not match!' )
#else
!do nothing if the debug CPP flag isn't on
#endif
  end subroutine mpp_pset_check_ptr

  subroutine mpp_pset_segment_array( pset, ls, le, lsp, lep )
!given input indices ls, le, returns indices lsp, lep
!so that segments span the range ls:le with no overlaps.
!attempts load balance: also some PSETs might get lsp>lep
!so that do-loops will be null
    type(mpp_pset_type), intent(in) :: pset
    integer, intent(in) :: ls, le
    integer, intent(out) :: lsp, lep
    integer :: i

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_SEGMENT_ARRAY: called with uninitialized PSET.' )
#ifdef PSET_DEBUG
    if( le-ls+1.LT.pset%npset )then
        write( text,'(3(a,i6))' ) &
             'MPP_PSET_ARRAY_SEGMENT: parallel range (', ls, ',', le, &
             ') is smaller than the number of threads:', pset%npset
        call mpp_error( WARNING, text )
    end if
#endif
    lep = ls-1 !initialize so that lsp is correct on first pass
    do i = 0,pset%pos
       lsp = lep + 1
       lep = lsp + CEILING( REAL(le-lsp+1)/(pset%npset-i) ) - 1
    end do
  end subroutine mpp_pset_segment_array

  subroutine mpp_pset_stack_push( pset, ptr, len )
!mpp_malloc specialized for shared arrays
!len is the length of the required array
!lstack is the stack already in play
!user should zero lstack (call mpp_pset_stack_reset) when the stack is to be cleared
    type(mpp_pset_type), intent(inout) :: pset
    integer, intent(in) :: len
#ifdef use_CRI_pointers
    real :: dummy
    pointer( ptr, dummy )
    real :: stack(pset%maxstack)
    pointer( p, stack )

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_STACK_PUSH: called with uninitialized PSET.' )
    if( pset%lstack+len.GT.pset%maxstack )then
        write( text, '(a,3i12)' ) &
             'MPP_PSET_STACK_PUSH: mpp_pset_stack overflow: '// &
             'len+lstack.GT.maxstack.  len, lstack, maxstack=', &
             len, pset%lstack, pset%maxstack
        call mpp_error( FATAL, text )
    end if
    p = pset%p_stack !point stack to shared stack pointer
    ptr = LOC( stack(pset%lstack+1) )
    call mpp_pset_check_ptr(pset,ptr) !make sure ptr is the same across PSETs
    pset%lstack = pset%lstack + len
    pset%hiWM = max( pset%hiWM, pset%lstack )
#else
    integer(POINTER_KIND), intent(out) :: ptr
    call mpp_error( FATAL, &
         'MPP_PSET_STACK_PUSH only works with Cray pointers.' )
#endif
  end subroutine mpp_pset_stack_push

  subroutine mpp_pset_stack_reset(pset)
    type(mpp_pset_type), intent(inout) :: pset
!reset stack... will reuse any temporary arrays! USE WITH CARE
!next few lines are to zero stack contents...
!but it's better noone tries to use uninitialized stack variables!
!    integer :: l1, l2
!    real :: mpp_pset_stack(maxstack)
!    pointer( p_mpp_pset_stack, mpp_pset_stack )
!    p_mpp_pset_stack = ptr_mpp_pset_stack
!    call mpp_pset_array_segment( 1, lstack, l1, l2 )
!    mpp_pset_stack(l1:l2) = 0.
    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_STACK_RESET: called with uninitialized PSET.' )
    pset%lstack = 0
  end subroutine mpp_pset_stack_reset

  subroutine mpp_pset_print_chksum_1D(pset, caller, array)
!print a checksum of an array
!pass the whole domain seen by root PSET
!add lines to check on shared array?
    type(mpp_pset_type), intent(in) :: pset
    character(len=*), intent(in) :: caller
    real, intent(in) :: array(:)

#ifdef PSET_DEBUG
    logical :: do_print
    integer(LONG_KIND) :: chksum

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_PRINT_CHKSUM: called with uninitialized PSET.' )

    if( pset%root )then
        do_print = pe.EQ.mpp_root_pe() !set to T to print from all PEs
        call mpp_set_current_pelist(pset%root_pelist)
        chksum = mpp_chksum( array )
        if( do_print ) &
             write( stderr(), '(a,z18)' )trim(caller)//' chksum=', chksum
    end if
    call mpp_set_current_pelist(pset%pelist)
#endif
    return
  end subroutine mpp_pset_print_chksum_1D

  subroutine mpp_pset_print_chksum_2D(pset, caller, array)
    type(mpp_pset_type), intent(in) :: pset
    character(len=*), intent(in) :: caller
    real, intent(in) :: array(:,:)
    real :: array1D( size(array) )
#ifdef use_CRI_pointers
    pointer( p, array1D )
    p = LOC(array)
#else
    array1D = TRANSFER( array, array1D )
#endif
    call mpp_pset_print_chksum(pset, caller, array1D)
  end subroutine mpp_pset_print_chksum_2D

  subroutine mpp_pset_print_chksum_3D(pset, caller, array)
    type(mpp_pset_type), intent(in) :: pset
    character(len=*), intent(in) :: caller
    real, intent(in) :: array(:,:,:) !overload for other ranks
    real :: array1D( size(array) )
#ifdef use_CRI_pointers
    pointer( p, array1D )
    p = LOC(array)
#else
    array1D = TRANSFER( array, array1D )
#endif
    call mpp_pset_print_chksum(pset, caller, array1D)
  end subroutine mpp_pset_print_chksum_3D

  subroutine mpp_pset_print_chksum_4D(pset, caller, array)
    type(mpp_pset_type), intent(in) :: pset
    character(len=*), intent(in) :: caller
    real, intent(in) :: array(:,:,:,:)
    real :: array1D( size(array) )
#ifdef use_CRI_pointers
    pointer( p, array1D )
    p = LOC(array)
#else
    array1D = TRANSFER( array, array1D )
#endif
    call mpp_pset_print_chksum(pset, caller, array1D)
  end subroutine mpp_pset_print_chksum_4D

  subroutine mpp_pset_print_stack_chksum( pset, caller )
    type(mpp_pset_type), intent(in) :: pset
    character(len=*), intent(in) :: caller

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_PRINT_STACK_CHKSUM: called with uninitialized PSET.' )
    call mpp_pset_print_chksum( pset, trim(caller)//' stack', &
         pset%stack(1:pset%lstack) )
  end subroutine mpp_pset_print_stack_chksum

!accessor functions
  function mpp_pset_root(pset)
    logical :: mpp_pset_root
    type(mpp_pset_type), intent(in) :: pset

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_ROOT: called with uninitialized PSET.' )
    mpp_pset_root = pset%root
  end function mpp_pset_root
  
  function mpp_pset_numroots(pset)
!necessary to export root_pelist: caller needs to pre-allocate
    integer :: mpp_pset_numroots
    type(mpp_pset_type), intent(in) :: pset

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_NUMROOTS: called with uninitialized PSET.' )
    mpp_pset_numroots = size(pset%root_pelist)
  end function mpp_pset_numroots

  subroutine mpp_pset_get_root_pelist(pset,pelist,commID)
    type(mpp_pset_type), intent(in) :: pset
    integer, intent(out) :: pelist(:)
    integer, intent(out), optional :: commID

    if( .NOT.pset%initialized )call mpp_error( FATAL, &
         'MPP_PSET_GET_ROOT_PELIST: called with uninitialized PSET.' )
    if( size(pelist).NE.size(pset%root_pelist) )then
        write( text,'(a,2i6)' ) &
             'pelist argument has wrong size: requested, actual=', &
             size(pelist), size(pset%root_pelist)
        call mpp_error( FATAL, 'MPP_PSET_GET_ROOT_PELIST: '//text )
    end if
    pelist(:) = pset%root_pelist(:)
    if( PRESENT(commID) )then
#ifdef use_libMPI
        commID = pset%commID
#else
        call mpp_error( WARNING, &
             'MPP_PSET_GET_ROOT_PELIST: commID is only defined under -Duse_libMPI.' )
#endif
    end if
  end subroutine mpp_pset_get_root_pelist
  
end module mpp_pset_mod

