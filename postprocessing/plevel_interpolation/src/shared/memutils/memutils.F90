module memutils_mod
!Author: Balaji (V.Balaji@noaa.gov)
!Various operations for memory management
!these currently include efficient methods for memory-to-memory copy
!including strided data and arbitrary gather-scatter vectors
!also various memory and cache inquiry operators
  implicit none
  private
#ifdef _CRAYT3E
  integer :: pe, shmem_my_pe
#endif

  integer(kind=8) :: l1_cache_line_size, l1_cache_size, l1_associativity
  integer(kind=8) :: l2_cache_line_size, l2_cache_size, l2_associativity

  logical :: memutils_initialized=.FALSE.

  interface memcpy
     module procedure memcpy_r8
     module procedure memcpy_r8_gather
     module procedure memcpy_r8_scatter
     module procedure memcpy_r8_gather_scatter
  end interface

  public :: get_l1_cache_line, get_l2_cache_line, memcpy, memutils_init
  public :: print_memuse_stats
#ifdef _CRAY
  public :: hplen
#endif
#ifdef _CRAYT90
  public :: stklen
#endif
  logical, private :: print_memory_usage=.FALSE.
  contains

    subroutine memutils_init(print_flag)
!initialize memutils module
!currently sets default cache characteristics
!(will provide overrides later)
!also sets pe to my_pe on t3e
      logical, optional :: print_flag
#ifdef _CRAYT3E
!all sizes in bytes
      l1_cache_line_size = 32
      l1_cache_size = 8192
      l1_associativity = 1
      l2_cache_line_size = 64
      l2_cache_size = 98304
      l2_associativity = 3
#else
!defaults
      l1_cache_line_size = 1
      l1_cache_size = 1
      l1_associativity = 1
      l2_cache_line_size = 1
      l2_cache_size = 1
      l2_associativity = 1
#endif
#ifdef _CRAYT3E
      pe = SHMEM_MY_PE()
#endif
      if( PRESENT(print_flag) )print_memory_usage = print_flag
      memutils_initialized = .TRUE.
      return
    end subroutine memutils_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!MEMCPY routines: <nelems> real*8 words are copied from RHS to LHS     !
!  Either side can have constant stride (lhs_stride, rhs_stride)       !
!      or indexed by a gather/scatter array (lhs_indx, rhs_indx)       !
! index arrays are 0-based (i.e C-like not fortran-like: this is       !
! for compatibility with the SHMEM_IXGET/PUT routines)                 !
!                                                                      !
! EXAMPLES:                                                            !
!                                                                      !
!Replace                                                               !
!  a(0:n-1) = b(0:n-1)                                                 !
!with                                                                  !
!  call memcpy(a,b,n)                                                  !
!                                                                      !
!Replace                                                               !
!  a(0:2*n-1:2) = b(0:3*n-1:3)                                         !
!with                                                                  !
!  call memcpy(a,b,dim,n,2,3)    !dim.GE.3*n                           !
!                                                                      !
!Replace                                                               !
!  a(0:n-1) = b(indx(1:n))                                             !
!with                                                                  !
!  call memcpy(a,b,dim,n,1,indx) !dim.GE.max(indx)                     !
!                                                                      !
!Replace                                                               !
!  a(indx(1:n)) = b(0:n-1)                                             !
!with                                                                  !
!  call memcpy(a,b,dim,n,indx,1) !dim.GE.max(indx)                     !
!                                                                      !
!Replace                                                               !
!  a(indxa(1:n)) = b(indxb(1:n))                                       !
!with                                                                  !
!  call memcpy(a,b,dim,n,indx,indxb) !dim.GE.max(indxa,indxb)          !
!                                                                      !
!  There are no error checks!!! (routines are built for speed)         !
!  Specifically there is no bounds-checking: if the stride or          !
!  indexing causes you to exceed <dim> you will have done a            !
!  potentially unsafe memory load                                      !
!                                                                      !
!T3E: we use the shmem routines on-processor to effect the transfer    !
!     via the (faster) E-registers                                     !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine memcpy_r8( lhs, rhs, dim, nelems, lhs_stride, rhs_stride )
!base routine: handles constant stride memcpy
!default strides are of course 1
      integer, intent(in) :: dim
      real(kind=8), dimension(0:dim-1), intent(in)  :: rhs
      real(kind=8), dimension(0:dim-1), intent(out) :: lhs
      integer, intent(in), optional :: nelems, lhs_stride, rhs_stride
      integer :: n, rs, ls

!defaults
      n = dim
      ls = 1
      rs = 1
      if( PRESENT(nelems) )then
          n = nelems
!only check for stride if nelems is present
          if( PRESENT(lhs_stride) )ls = lhs_stride
          if( PRESENT(rhs_stride) )rs = rhs_stride
      endif
      if( ls.EQ.1 .AND. rs.EQ.1 )then
#ifdef _CRAYT3E
          call SHMEM_GET( lhs(0), rhs(0), n, pe )
#else
          lhs(0:n-1) = rhs(0:n-1)
#endif
      else
#ifdef _CRAYT3E
          call SHMEM_IGET( lhs(0), rhs(0), ls, rs, n, pe )
#else
          lhs(0:n*ls-1:ls) = rhs(0:n*rs-1:rs)
#endif
      endif
      return
    end subroutine memcpy_r8

    subroutine memcpy_r8_gather( lhs, rhs, dim, nelems, lhs_stride, rhs_indx )
!memcpy routine with gather: copies nelems words from rhs(indx(:)) to lhs(:)
      integer, intent(in) :: dim, nelems, lhs_stride
      real(kind=8), dimension(0:dim-1), intent(in)  :: rhs
      real(kind=8), dimension(0:dim-1), intent(out) :: lhs
      integer, intent(in), dimension(nelems) :: rhs_indx
#ifdef _CRAYT3E
!dir$ CACHE_BYPASS lhs, rhs, rhs_indx
      real(kind=8), dimension(nelems) :: tmp

      if( lhs_stride.EQ.1 )then
          call SHMEM_IXGET( lhs(0), rhs(0), rhs_indx, nelems, pe )
      else
          call SHMEM_IXGET( tmp, rhs(0), rhs_indx, nelems, pe )
          call SHMEM_IGET( lhs(0), tmp, lhs_stride, 1, nelems, pe )
      endif
#else
      lhs(0:nelems*lhs_stride-1:lhs_stride) = rhs(rhs_indx(1:nelems))
#endif
      return
    end subroutine memcpy_r8_gather

    subroutine memcpy_r8_scatter( lhs, rhs, dim, nelems, lhs_indx, rhs_stride )
!memcpy routine with scatter: copies nelems words from rhs(:) to lhs(indx(:))
      integer, intent(in) :: dim, nelems, rhs_stride
      real(kind=8), dimension(0:dim-1), intent(in)  :: rhs
      real(kind=8), dimension(0:dim-1), intent(out) :: lhs
      integer, intent(in), dimension(nelems) :: lhs_indx
#ifdef _CRAYT3E
!dir$ CACHE_BYPASS lhs, rhs, lhs_indx
      real(kind=8), dimension(nelems) :: tmp

      if( rhs_stride.EQ.1 )then
          call SHMEM_IXPUT( lhs(0), rhs(0), lhs_indx, nelems, pe )
      else
          call SHMEM_IGET( tmp, rhs(0), rhs_stride, 1, nelems, pe )
          call SHMEM_IXPUT( lhs(0), tmp, lhs_indx, nelems, pe )
      endif
      call SHMEM_QUIET          !required to ensure completion of put
#else
      lhs(lhs_indx(1:nelems)) = rhs(0:nelems*rhs_stride-1:rhs_stride)
#endif
      return
    end subroutine memcpy_r8_scatter

    subroutine memcpy_r8_gather_scatter( lhs, rhs, dim, nelems, lhs_indx, rhs_indx )
!memcpy routine with gather/scatter: copies nelems words from rhs(indx(:)) to lhs(indx(:))
      integer, intent(in) :: dim, nelems
      real(kind=8), dimension(0:dim-1), intent(in)  :: rhs
      real(kind=8), dimension(0:dim-1), intent(out) :: lhs
      integer, intent(in), dimension(nelems) :: lhs_indx, rhs_indx
#ifdef _CRAYT3E
!dir$ CACHE_BYPASS lhs, rhs, lhs_indx, rhs_indx
      real(kind=8), dimension(nelems) :: tmp

      call SHMEM_IXGET( tmp, rhs(0), rhs_indx, nelems, pe )
      call SHMEM_IXPUT( lhs(0), tmp, lhs_indx, nelems, pe )
      call SHMEM_QUIET          !required to ensure completion of put
#else
      lhs(lhs_indx(1:nelems)) = rhs(rhs_indx(1:nelems))
#endif
      return
    end subroutine memcpy_r8_gather_scatter

#ifdef _CRAY
  integer function hplen(             hpalloc, hplargest, hpshrink, hpgrow, hpfirst, hplast )
!using IHPSTAT calls from SR-2165 v2.0 p535
!with no arguments returns heap length (in words on PVP, bytes on t3e)
    integer, intent(out), optional :: hpalloc, hplargest, hpshrink, hpgrow, hpfirst, hplast
    integer :: IHPSTAT

    hplen = IHPSTAT(1)	                      !Heap length
    if( present(hpalloc  ) )hpalloc   = IHPSTAT( 4) !Blocks allocated
    if( present(hplargest) )hplargest = IHPSTAT(10) !Largest free block size
    if( present(hpshrink ) )hpshrink  = IHPSTAT(11) !Amount heap can shrink
    if( present(hpgrow   ) )hpgrow    = IHPSTAT(12) !Amount heap can grow
    if( present(hpfirst  ) )hpfirst   = IHPSTAT(13) !First word address
    if( present(hplast   ) )hplast    = IHPSTAT(14) !Last word address
    return
  end function hplen
#endif /* _CRAY */

#ifdef _CRAYT90
  integer function stklen(            stkhiwm, stknumber, stktotal, stkmost, stkgrew, stkgtimes )
!using STKSTAT(3C) struct
    integer, optional, intent(out) :: stkhiwm, stknumber, stktotal, stkmost, stkgrew, stkgtimes
    integer :: istat(20)

    call STKSTAT(istat)
    stklen = istat(1)	!Stack length
    if( present(stkhiwm  ) )stkhiwm   = istat(2) !stack hiwatermark
    if( present(stknumber) )stknumber = istat(3) !current #stacks
    if( present(stktotal ) )stktotal  = istat(4) !total #stacks
    if( present(stkmost  ) )stkmost   = istat(5) !most #stacks at one time
    if( present(stkgrew  ) )stkgrew   = istat(6) !#stacks that grew
    if( present(stkgtimes) )stkgtimes = istat(7) !#times stack grew
    return
  end function stklen
#endif /* _CRAYT90 */

!cache utilities: need to write version for other argument types
  function get_l1_cache_line(a)
    integer(kind=8) :: get_l1_cache_line
    real, intent(in) :: a
    integer(kind=8) :: i
    i = LOC(a)
    get_l1_cache_line = mod(i,l1_cache_size/l1_associativity)/l1_cache_line_size
  end function get_l1_cache_line

  function get_l2_cache_line(a)
    integer(kind=8) :: get_l2_cache_line
    real, intent(in) :: a
    integer(kind=8) :: i
    i = LOC(a)
    get_l2_cache_line = mod(i,l2_cache_size/l2_associativity)/l2_cache_line_size
  end function get_l2_cache_line

  subroutine print_memuse_stats( text, unit, always )
    use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_npes, mpp_min, mpp_max, mpp_sum, stderr
    character(len=*), intent(in) :: text
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: always
    real :: m, mmin, mmax, mavg, mstd
    integer :: mu
!memuse is an external function: works on SGI
!use #ifdef to generate equivalent on other platforms.
    integer :: memuse !default integer OK?

    if( PRESENT(always) )then
        if( .NOT.always )return
    else
        if( .NOT.print_memory_usage )return
    end if
    mu = stderr(); if( PRESENT(unit) )mu = unit
#if defined(__sgi) || defined(__aix) || defined(__SX)
    m = memuse()*1e-3
#else
    call mem_dump(m)
#endif 
    mmin = m; call mpp_min(mmin)
    mmax = m; call mpp_max(mmax)
    mavg = m; call mpp_sum(mavg); mavg = mavg/mpp_npes()
    mstd = (m-mavg)**2; call mpp_sum(mstd); mstd = sqrt( mstd/mpp_npes() )
    if( mpp_pe().EQ.mpp_root_pe() )write( mu,'(a64,4es11.3)' ) &
         'Memuse(MB) at '//trim(text)//'=', mmin, mmax, mstd, mavg

    return
  end subroutine print_memuse_stats

!#######################################################################

subroutine mem_dump ( memuse )
use mpp_mod,    only : stdout
use mpp_io_mod, only : mpp_open, mpp_close, mpp_ascii, mpp_rdonly,     &
                       mpp_sequential, mpp_single

real, intent(out) :: memuse

! This routine returns the memory usage on Linux systems.
! It does this by querying a system file (file_name below).
! It is intended for use by print_memuse_stats above.

character(len=32) :: file_name = '/proc/self/status'
character(len=32) :: string
integer :: mem_unit
real    :: multiplier

  memuse = 0.0
  multiplier = 1.0

  call mpp_open ( mem_unit, file_name,                                 &
                      form=MPP_ASCII,        action=MPP_RDONLY,        &
                      access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
  
  do; read (mem_unit,'(a)', end=10) string
    if ( INDEX ( string, 'VmHWM:' ) == 1 ) then
      read (string(7:LEN_TRIM(string)-2),*) memuse
      exit
    endif
  enddo
  
  if (TRIM(string(LEN_TRIM(string)-1:)) == "kB" ) &
    multiplier = 1.0/1024. ! Convert from kB to MB

10 call mpp_close ( mem_unit )
   memuse = memuse * multiplier

  return
end subroutine mem_dump

end module memutils_mod
