#ifdef test_mpp_domains
program test
  use mpp_mod,         only : FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_node, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync, mpp_sync_self
  use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,         only : mpp_init, mpp_exit, mpp_chksum, stdout, stderr
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, CGRID_NE
  use mpp_domains_mod, only : FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE, FOLD_WEST_EDGE, FOLD_EAST_EDGE
  use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
  use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field, mpp_redistribute, mpp_get_memory_domain
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod, only : mpp_get_neighbor_pe, mpp_define_mosaic, mpp_nullify_domain_list
  use mpp_domains_mod, only : NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
  use mpp_domains_mod, only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST, mpp_define_mosaic_pelist
  use mpp_domains_mod, only : mpp_get_refine_overlap_number, mpp_get_mosaic_refine_overlap
  use mpp_domains_mod, only : mpp_get_global_domain, ZERO, NINETY, MINUS_NINETY
  use mpp_domains_mod, only : mpp_get_boundary
  use mpp_memutils_mod, only : mpp_memuse_begin, mpp_memuse_end

  implicit none
#include <fms_platform.h>
  integer :: pe, npes
  integer :: nx=128, ny=128, nz=40, stackmax=4000000
  integer :: unit=7
  integer :: stdunit = 6
  logical :: debug=.FALSE., opened
  logical :: check_parallel = .FALSE.  ! when check_parallel set to false,
                                       ! mpes should be equal to npes     
  integer :: mpes = 0
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  integer :: x_cyclic_offset = 3   ! to be used in test_cyclic_offset
  integer :: y_cyclic_offset = -4  ! to be used in test_cyclic_offset
  character(len=32) :: warn_level = "fatal"
  namelist / test_mpp_domains_nml / nx, ny, nz, stackmax, debug, mpes, check_parallel, &
                               whalo, ehalo, shalo, nhalo, x_cyclic_offset, y_cyclic_offset, &
                               warn_level
  integer :: i, j, k
  integer :: layout(2)
  integer :: id

  call mpp_memuse_begin()
  call mpp_init()
 
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, status='OLD', file='input.nml', err=10 )
  read( unit,test_mpp_domains_nml )
  close(unit)
10 continue

  select case(trim(warn_level))
  case("fatal")
     call mpp_set_warn_level(FATAL)
  case("warning")
     call mpp_set_warn_level(WARNING)
  case default
     call mpp_error(FATAL, "test_mpp_domains: warn_level should be fatal or warning")
  end select
  
  pe = mpp_pe()
  npes = mpp_npes()
  
  if( debug )then
      call mpp_domains_init(MPP_DEBUG)
  else
      call mpp_domains_init(MPP_DOMAIN_TIME)
  end if
  call mpp_domains_set_stack_size(stackmax)
  
  if( pe.EQ.mpp_root_pe() )print '(a,9i6)', 'npes, mpes, nx, ny, nz, whalo, ehalo, shalo, nhalo =', &
                           npes, mpes, nx, ny, nz, whalo, ehalo, shalo, nhalo
  call mpp_memuse_end("in the begining", stdout())  
  if( .not. check_parallel) then
      call test_modify_domain()
      call test_cyclic_offset('x_cyclic_offset')
      call test_cyclic_offset('y_cyclic_offset')
      call test_cyclic_offset('torus_x_offset')
      call test_cyclic_offset('torus_y_offset')

      call test_get_boundary('Four-Tile')
      call test_get_boundary('Cubic-Grid')
      call test_uniform_mosaic('Single-Tile')
      call test_uniform_mosaic('Folded-north mosaic') ! one-tile tripolar grid
      call test_uniform_mosaic('Folded-north symmetry mosaic') ! one-tile tripolar grid
      call test_uniform_mosaic('Folded-south symmetry mosaic') ! one-tile tripolar grid
      call test_uniform_mosaic('Folded-west symmetry mosaic') ! one-tile tripolar grid
      call test_uniform_mosaic('Folded-east symmetry mosaic') ! one-tile tripolar grid
      call test_uniform_mosaic('Four-Tile')
      call test_uniform_mosaic('Cubic-Grid') ! 6 tiles.
      call test_nonuniform_mosaic('Five-Tile')
      call test_refined_mosaic('Refined-Four-Tile')
      call test_refined_mosaic('Refined-Symmetric-Four-Tile')
      call test_refined_mosaic('Refined-Cubic-Grid')

      call test_halo_update( 'Simple' ) !includes global field, global sum tests
      call test_halo_update( 'Cyclic' )
      call test_halo_update( 'Folded-north' ) !includes vector field test
      call test_halo_update( 'Masked' ) !includes vector field test
      call test_halo_update( 'Folded xy_halo' ) ! 

      call test_halo_update( 'Simple symmetry' ) !includes global field, global sum tests
      call test_halo_update( 'Cyclic symmetry' )
      call test_halo_update( 'Folded-north symmetry' ) !includes vector field test
      call test_halo_update( 'Folded-south symmetry' ) !includes vector field test
      call test_halo_update( 'Folded-west symmetry' ) !includes vector field test
     call test_halo_update( 'Folded-east symmetry' ) !includes vector field test

      !--- z1l: The following will not work due to symmetry and domain%x is cyclic.
      !--- Will solve this problem in the future if needed.
      ! call test_halo_update( 'Masked symmetry' ) !includes vector field test

      call test_global_field( 'Non-symmetry' )
      call test_global_field( 'Symmetry center' )
      call test_global_field( 'Symmetry corner' )
      call test_global_field( 'Symmetry east' )
      call test_global_field( 'Symmetry north' )

      call test_global_reduce( 'Simple')
      call test_global_reduce( 'Simple symmetry center')
      call test_global_reduce( 'Simple symmetry corner')
      call test_global_reduce( 'Simple symmetry east')
      call test_global_reduce( 'Simple symmetry north')
      call test_global_reduce( 'Cyclic symmetry center')
      call test_global_reduce( 'Cyclic symmetry corner')
      call test_global_reduce( 'Cyclic symmetry east')
      call test_global_reduce( 'Cyclic symmetry north')

      call test_redistribute( 'Complete pelist' )
!!$      call test_redistribute( 'Overlap  pelist' )
!!$      call test_redistribute( 'Disjoint pelist' )

      call test_define_mosaic_pelist('One tile', 1)
      call test_define_mosaic_pelist('Two uniform tile', 2)
      call test_define_mosaic_pelist('Two nonuniform tile', 2)
      call test_define_mosaic_pelist('Ten tile', 10)
      call test_define_mosaic_pelist('Ten tile with nonuniform cost', 10)
  else
      call test_parallel( )
  endif

!!$!Balaji adding openMP tests
!!$  call test_openmp()
!!$ 
!!$! Alewxander.Pletzer get_neighbor tests
!!$  call test_get_neighbor_1d
!!$  call test_get_neighbor_non_cyclic
!!$  call test_get_neighbor_cyclic
!!$  call test_get_neighbor_folded_north
!!$  call test_get_neighbor_mask

  call mpp_domains_exit()
  call mpp_exit()
  
contains
  subroutine test_openmp()
#ifdef _OPENMP
    integer :: omp_get_num_thread, omp_get_max_threads, omp_get_thread_num
    real, allocatable :: a(:,:,:)
    type(domain2D) :: domain
    integer :: layout(2)
    integer :: i,j,k, jthr
    integer :: thrnum, maxthr
    integer(LONG_KIND) :: sum1, sum2
    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain )
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( a(isd:ied,jsd:jed,nz) )
    maxthr = omp_get_max_threads()
    write( stdout(),'(a,4i6)' )'pe,js,je,maxthr=', pe, js, je, maxthr
    if( mod(je-js+1,maxthr).NE.0 ) &
         call mpp_error( FATAL, 'maxthr must divide domain (TEMPORARY).' )
    jthr = (je-js+1)/maxthr
!$OMP PARALLEL PRIVATE(i,j,k,thrnum)
    thrnum = omp_get_thread_num()
    write( stdout(),'(a,4i6)' )'pe,thrnum,js,je=', &
         pe, thrnum, js+thrnum*jthr,js+(thrnum+1)*jthr-1
    write( stdout(),'(a,3i6)' )'pe,thrnum,node=', pe, thrnum, mpp_node()
!!$OMP DO
    do k = 1,nz
!when omp DO is commented out, user must compute j loop limits
!with omp DO, let OMP figure it out
       do j = js+thrnum*jthr,js+(thrnum+1)*jthr-1
!       do j = js,je
          do i = is,ie
             a(i,j,k) = global(i,j,k)
          end do
       end do
    end do
!!$OMP END DO
!$OMP END PARALLEL
    sum1 = mpp_chksum( a(is:ie,js:je,:) )
    sum2 = mpp_chksum( global(is:ie,js:je,:) )
    if( sum1.EQ.sum2 )then
        call mpp_error( NOTE, 'OMP parallel test OK.' )
    else
        if( mpp_pe().EQ.mpp_root_pe() )write( stderr(),'(a,2z18)' )'OMP checksums: ', sum1, sum2
        call mpp_error( FATAL, 'OMP parallel test failed.' )
    end if
#endif
    return
  end subroutine test_openmp

  subroutine test_redistribute( type )
!test redistribute between two domains
    character(len=*), intent(in) :: type
    type(domain2D) :: domainx, domainy
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    real, allocatable, dimension(:,:,:)       :: gcheck, global
    real, allocatable, dimension(:,:,:), save :: x, y
    real, allocatable, dimension(:,:,:), save :: x2, y2
    real, allocatable, dimension(:,:,:), save :: x3, y3
    real, allocatable, dimension(:,:,:), save :: x4, y4
    real, allocatable, dimension(:,:,:), save :: x5, y5
    real, allocatable, dimension(:,:,:), save :: x6, y6
    integer, allocatable :: pelist(:)
    integer :: pemax
    integer :: is, ie, js, je, isd, ied, jsd, jed
    
    pemax = npes/2              !the partial pelist will run from 0...pemax
    !--- nullify domain list otherwise it retains memory between calls.
    call mpp_nullify_domain_list(domainx)
    call mpp_nullify_domain_list(domainy)

    allocate( gcheck(nx,ny,nz), global(nx,ny,nz) )
    !fill in global array: with k.iiijjj
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             global(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

!select pelists
    select case(type)
    case( 'Complete pelist' )
!both pelists run from 0...npes-1
        allocate( pelist(0:npes-1) )
        pelist = (/ (i,i=0,npes-1) /)
        call mpp_declare_pelist( pelist )
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)
        call mpp_declare_pelist( pelist )
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1
        if( pemax+1.GE.npes )return
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)

        call mpp_declare_pelist( pelist )
        ! z1l: the follwing will cause deadlock will happen
        ! for npes = 6, x- mpp_global_field will call mpp_sync
        call mpp_declare_pelist( (/ (i,i=pemax+1,npes-1) /))
    case default
        call mpp_error( FATAL, 'TEST_REDISTRIBUTE: no such test: '//type )
    end select
        
!set up x and y arrays
    select case(type)
    case( 'Complete pelist' )
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        allocate( x2(isd:ied,jsd:jed,nz) )
        allocate( x3(isd:ied,jsd:jed,nz) )
        allocate( x4(isd:ied,jsd:jed,nz) )
        allocate( x5(isd:ied,jsd:jed,nz) )
        allocate( x6(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
        x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
!set up y array
        call mpp_define_domains( (/1,nx,1,ny/), (/npes,1/), domainy, name=type )
        call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
        allocate( y(isd:ied,jsd:jed,nz) )
        allocate( y2(isd:ied,jsd:jed,nz) )
        allocate( y3(isd:ied,jsd:jed,nz) )
        allocate( y4(isd:ied,jsd:jed,nz) )
        allocate( y5(isd:ied,jsd:jed,nz) )
        allocate( y6(isd:ied,jsd:jed,nz) )
        y = 0.
        y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        allocate( x2(isd:ied,jsd:jed,nz) )
        allocate( x3(isd:ied,jsd:jed,nz) )
        allocate( x4(isd:ied,jsd:jed,nz) )
        allocate( x5(isd:ied,jsd:jed,nz) )
        allocate( x6(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
        x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            allocate( y2(isd:ied,jsd:jed,nz) )
            allocate( y3(isd:ied,jsd:jed,nz) )
            allocate( y4(isd:ied,jsd:jed,nz) )
            allocate( y5(isd:ied,jsd:jed,nz) )
            allocate( y6(isd:ied,jsd:jed,nz) )
            y = 0.
            y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
        end if
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1
    
!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            allocate( y2(isd:ied,jsd:jed,nz) )
            allocate( y3(isd:ied,jsd:jed,nz) )
            allocate( y4(isd:ied,jsd:jed,nz) )
            allocate( y5(isd:ied,jsd:jed,nz) )
            allocate( y6(isd:ied,jsd:jed,nz) )
            y = 0.
            y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
        else
!set up x array
            call mpp_set_current_pelist( (/ (i,i=pemax+1,npes-1) /) )
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
            call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
            allocate( x(isd:ied,jsd:jed,nz) )
            allocate( x2(isd:ied,jsd:jed,nz) )
            allocate( x3(isd:ied,jsd:jed,nz) )
            allocate( x4(isd:ied,jsd:jed,nz) )
            allocate( x5(isd:ied,jsd:jed,nz) )
            allocate( x6(isd:ied,jsd:jed,nz) )
            x = 0.
            x(is:ie,js:je,:) = global(is:ie,js:je,:)
            x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
         end if
    end select
         
!go global and redistribute
    call mpp_set_current_pelist()
    call mpp_broadcast_domain(domainx)
    call mpp_broadcast_domain(domainy)
    
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_redistribute( domainx, x, domainy, y )
    call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if
        
    call mpp_set_current_pelist()

    call mpp_clock_begin(id)
    if(ALLOCATED(y))y=0.
    call mpp_redistribute( domainx, x,  domainy, y,  complete=.false. )
    call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
    call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
    call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
    call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
    call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
    call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y2, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y3, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y4, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y5, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y6, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if

    call mpp_set_current_pelist()

    if(type == 'Complete pelist')then
      write(stdout(),*) 'Use domain communicator handle'
      call mpp_clock_begin(id)
      if(ALLOCATED(y))then
         y=0.; y2=0.; y3=0.; y4=0.; y5=0.; y6=0.
      endif
      call mpp_redistribute( domainx, x, domainy, y, complete=.false. )
      call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
      call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
      call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
      call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
      call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
      call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y2, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y3, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y4, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y5, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y6, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if
    endif
    dch =>NULL()
        
    call mpp_set_current_pelist()

    call mpp_sync()
    
    deallocate(gcheck, global)
    if(ALLOCATED(pelist)) deallocate(pelist)

    if(ALLOCATED(x))then
      call mpp_redistribute( domainx, x, domainy, y, free=.true.,list_size=6 )
      deallocate(x,x2,x3,x4,x5,x6)
    endif
    if(ALLOCATED(y))deallocate(y,y2,y3,y4,y5,y6)
  end subroutine test_redistribute


  subroutine test_uniform_mosaic( type )
    character(len=*), intent(in) :: type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile, ntile_per_pe, update_flags
    integer        :: i, j, k, l, n, shift, tw, te, ts, tn, tsw, tnw, tse, tne
    integer        :: ism, iem, jsm, jem, wh, eh, sh, nh
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    real           :: gsum, lsum  

    integer, allocatable, dimension(:)       :: tile
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:)     :: global2D
    real,    allocatable, dimension(:,:,:)   :: local1, local2
    real,    allocatable, dimension(:,:,:,:) :: x, y, x1, x2, x3, x4, y1, y2, y3, y4
    real,    allocatable, dimension(:,:,:,:) :: global1, global2, gcheck  
    real,    allocatable, dimension(:,:,:,:) :: global1_all, global2_all, global_all
    character(len=128) :: type2, type3
    logical            :: folded_north, folded_north_sym, folded_north_nonsym
    logical            :: folded_south_sym, folded_west_sym, folded_east_sym
    logical            :: cubic_grid, single_tile, four_tile

    folded_north_nonsym = .false.
    folded_north_sym    = .false.
    folded_north        = .false.
    folded_south_sym    = .false.
    folded_west_sym     = .false.
    folded_east_sym     = .false.
    cubic_grid        = .false.
    single_tile        = .false.
    four_tile          = .false.
    !--- check the type
    select case(type)
    case ( 'Single-Tile' )   !--- single with cyclic along x- and y-direction
       single_tile = .true.
       ntiles = 1
       num_contact = 2
    case ( 'Folded-north mosaic' )
       ntiles = 1
       num_contact = 2
       folded_north_nonsym = .true.
    case ( 'Folded-north symmetry mosaic' )
       ntiles = 1
       num_contact = 2
       folded_north_sym = .true.
    case ( 'Folded-south symmetry mosaic' )
       ntiles = 1
       num_contact = 2
       folded_south_sym = .true.
    case ( 'Folded-west symmetry mosaic' )
       ntiles = 1
       num_contact = 2
       folded_west_sym = .true.
    case ( 'Folded-east symmetry mosaic' )
       ntiles = 1
       num_contact = 2
       folded_east_sym = .true.
    case ( 'Four-Tile' ) !--- cyclic along both x- and y-direction. 
       ntiles = 4
       num_contact = 8
       four_tile = .true.
    case ( 'Cubic-Grid' )
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.
       if( nx .NE. ny) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cubic_grid mosaic, nx should equal ny, '//&
                   'No test is done for Cubic-Grid mosaic. ' )
          return
       end if
    case default
       call mpp_error(FATAL, 'TEST_MPP_DOMAINS: no such test: '//type)
    end select

    folded_north = folded_north_nonsym .OR. folded_north_sym
      
    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    if( mod(npes, ntiles) == 0 ) then
       npes_per_tile = npes/ntiles
       write(stdout(),*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                       '", each tile will be distributed over ', npes_per_tile, ' processors.'
       ntile_per_pe = 1
       allocate(tile(ntile_per_pe))
       tile = pe/npes_per_tile+1
       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       do n = 1, ntiles
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do
    else if ( mod(ntiles, npes) == 0 ) then
       ntile_per_pe = ntiles/npes 
       write(stdout(),*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                        '", there will be ', ntile_per_pe, ' tiles on each processor.'
       allocate(tile(ntile_per_pe))
       do n = 1, ntile_per_pe
          tile(n) = pe*ntile_per_pe + n
       end do
       do n = 1, ntiles
          pe_start(n) = (n-1)/ntile_per_pe
          pe_end(n)   = pe_start(n)
       end do
       layout = 1
    else
       call mpp_error(NOTE,'TEST_MPP_DOMAINS: npes should be multiple of ntiles or ' // &
            'ntiles should be multiple of npes. No test is done for '//trim(type) )       
       return
    end if
 
    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) ) 
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) ) 

    call mpp_memuse_begin()
    !--- define domain
    if(single_tile) then
       !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
       istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                              name = type, symmetry = .false. )
    else if(folded_north) then
       !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (NORTH) and tile 1 (NORTH)  --- folded-north-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = ny;  jend1(2) = ny
       istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = ny;  jend2(2) = ny
       if(folded_north_nonsym) then
          call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                                 istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                                 pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                                 name = type, symmetry = .false.  )
       else
          call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                                 istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                                 pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                                 name = type, symmetry = .true.  )
       endif
    else if(folded_south_sym) then
       !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (SOUTH)  --- folded-south-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = 1;  jend1(2) = 1
       istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = 1;  jend2(2) = 1
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                              name = type, symmetry = .true.  )
    else if(folded_west_sym) then
       !--- Contact line 1, between tile 1 (NORTH) and tile 1 (SOUTH)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = 1; iend1(1) = nx; jstart1(1) = ny;  jend1(1) = ny
       istart2(1) = 1; iend2(1) = nx; jstart2(1) = 1;   jend2(1) = 1
       !--- Contact line 2, between tile 1 (WEST) and tile 1 (WEST)  --- folded-west-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = 1; jstart1(2) = 1;  jend1(2) = ny/2
       istart2(2) = 1;  iend2(2) = 1; jstart2(2) = ny; jend2(2) = ny/2+1
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                              name = type, symmetry = .true.  )
    else if(folded_east_sym) then
       !--- Contact line 1, between tile 1 (NORTH) and tile 1 (SOUTH)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = 1; iend1(1) = nx; jstart1(1) = ny;  jend1(1) = ny
       istart2(1) = 1; iend2(1) = nx; jstart2(1) = 1;   jend2(1) = 1
       !--- Contact line 2, between tile 1 (EAST) and tile 1 (EAST)  --- folded-west-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = nx;  iend1(2) = nx; jstart1(2) = 1;  jend1(2) = ny/2
       istart2(2) = nx;  iend2(2) = nx; jstart2(2) = ny; jend2(2) = ny/2+1
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                              name = type, symmetry = .true.  )
    else if( four_tile ) then
       call define_fourtile_mosaic(type, domain, (/nx,nx,nx,nx/), (/ny,ny,ny,ny/), global_indices, &
                                   layout2D, pe_start, pe_end, symmetry = .false. )
    else if( cubic_grid ) then
       call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                                global_indices, layout2D, pe_start, pe_end )
    endif
    call mpp_memuse_end(trim(type)//" mpp_define_mosaic", stdout() )

    !--- setup data
    allocate(global2(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz, ntile_per_pe) ) 
    allocate(global_all(1:nx,1:ny,nz, ntiles) )    
    global2 = 0
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                global_all(i,j,k,l) = l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global2(1:nx,1:ny,:,n) = global_all(:,:,:,tile(n))
    end do

    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( gcheck(nx, ny, nz, ntile_per_pe) )
    allocate( x (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x1(ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x2(ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x3(ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x4(ism:iem,jsm:jem,nz, ntile_per_pe) )
    x = 0.
    x(isc:iec,jsc:jec,:,:) = global2(isc:iec,jsc:jec,:,:)
    x1 = x; x2 = x; x3 = x; x4 = x;

    !--- test mpp_global_sum
    gsum = 0
    allocate(global2D(nx,ny))
    do n = 1, ntiles
       do j = 1, ny
          do i = 1, nx
             global2D(i,j) = sum(global_all(i,j,:,n))
          end do
       end do
       gsum = gsum + sum(global2D)
    end do

    do n = 1, ntile_per_pe  
       lsum = mpp_global_sum( domain, x(:,:,:,n), tile_count=n )
    end do  
    if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', type//' Fast sum=', lsum, gsum

    !test exact mpp_global_sum
    do n = 1, ntile_per_pe  
       lsum = mpp_global_sum( domain, x(:,:,:,n), BITWISE_EXACT_SUM, tile_count=n)
    end do 
    call compare_data_scalar(lsum, gsum, FATAL, type//' mpp_global_exact_sum')

    !--- test mpp_global_field
    gcheck = 0.    
    id = mpp_clock_id( type//' global field ', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    do n = 1, ntile_per_pe
       call mpp_global_field( domain, x(:,:,:,n), gcheck(:,:,:,n), tile_count=n)
    end do
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    do n = 1, ntile_per_pe
       call compare_checksums( global2(1:nx,1:ny,:,n), gcheck(:,:,:,n), type//' mpp_global_field ' )
    end do

    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    do n = 1, ntile_per_pe
       !--- fill up the value at halo points.
       if(single_tile) then
          call fill_regular_mosaic_halo(global2(:,:,:,n), global_all, 1, 1, 1, 1, 1, 1, 1, 1)
       else if(folded_north) then
          call fill_folded_north_halo(global2(:,:,:,n), 0, 0, 0, 0, 1)
       else if(folded_south_sym) then
          call fill_folded_south_halo(global2(:,:,:,n), 0, 0, 0, 0, 1)
       else if(folded_west_sym) then
          call fill_folded_west_halo(global2(:,:,:,n), 0, 0, 0, 0, 1)
       else if(folded_east_sym) then
          call fill_folded_east_halo(global2(:,:,:,n), 0, 0, 0, 0, 1)
       else if(four_tile) then
          select case ( tile(n) )
          case (1)
             tw = 2; ts = 3; tsw = 4
          case (2)
             tw = 1; ts = 4; tsw = 3
          case (3)
             tw = 4; ts = 1; tsw = 2
          case (4)
             tw = 3; ts = 2; tsw = 1
          end select
          te = tw; tn = ts; tse = tsw; tnw = tsw; tne = tsw
          call fill_regular_mosaic_halo(global2(:,:,:,n), global_all, te, tse, ts, tsw, tw, tnw, tn, tne )
       else if(cubic_grid) then
          call fill_cubic_grid_halo(global2(:,:,:,n), global_all, global_all, tile(n), 0, 0, 1, 1 )
       endif

       !full update
       call mpp_clock_begin(id)
       if(ntile_per_pe == 1) then
          call mpp_update_domains( x(:,:,:,n), domain )
       else
          call mpp_update_domains( x(:,:,:,n), domain, tile_count = n )
       end if
       call mpp_clock_end  (id)
    end do
    type2 = type
    do n = 1, ntile_per_pe  
       if(ntile_per_pe>1)   write(type2, *)type, " at tile_count = ",n
       call compare_checksums( x(ism:ism+ied-isd,jsm:jsm+jed-jsd,:,n), global2(isd:ied,jsd:jed,:,n), trim(type2) )
    end do

    !partial update only be done when there is at most one tile on each pe
    if(ntile_per_pe == 1 ) then
       id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
       call mpp_clock_begin(id)
       call mpp_update_domains( x1, domain, NUPDATE+EUPDATE, complete=.false. )
       call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
       call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
       call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
       call mpp_clock_end  (id)
       call compare_checksums( x1(isc:ied,jsc:jed,:,1), global2(isc:ied,jsc:jed,:,1), type//' partial x1' )
       call compare_checksums( x2(isc:ied,jsc:jed,:,1), global2(isc:ied,jsc:jed,:,1), type//' partial x2' )
       call compare_checksums( x3(isc:ied,jsc:jed,:,1), global2(isc:ied,jsc:jed,:,1), type//' partial x3' )
       call compare_checksums( x4(isc:ied,jsc:jed,:,1), global2(isc:ied,jsc:jed,:,1), type//' partial x4' )

       !arbitrary halo update. not for tripolar grid
       if(single_tile .or. four_tile .or. cubic_grid ) then
          allocate(local2(isd:ied,jsd:jed,nz) )
          do wh = 1-whalo, whalo
             do eh = 1-ehalo, ehalo
                do sh = 1-shalo, shalo
                   do nh = 1-nhalo, nhalo
                      if( wh*eh <= 0 ) cycle
                      if( sh*nh <= 0 ) cycle
                      if( wh*sh <= 0 ) cycle
                      local2(isd:ied,jsd:jed,:) = global2(isd:ied,jsd:jed,:,1)
                      x = 0.
                      x(isc:iec,jsc:jec,:,1) = local2(isc:iec,jsc:jec,:)       
                      call fill_halo_zero(local2, wh, eh, sh, nh, 0, 0, isc, iec, jsc, jec, isd, ied, jsd, jed) 

                      write(type2,'(a,a,i2,a,i2,a,i2,a,i2)') trim(type), ' with whalo = ', wh, &
                           ', ehalo = ',eh, ', shalo = ', sh, ', nhalo = ', nh
                      call mpp_update_domains( x, domain, whalo=wh, ehalo=eh, shalo=sh, nhalo=nh, name = type2  )
                      call compare_checksums( x(isd:ied,jsd:jed,:,1), local2, trim(type2) )
                   end do
                end do
             end do
          end do
          deallocate(local2)
       end if
    end if

    deallocate(global2, global_all, x, x1, x2, x3, x4)
    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    !--- setup data
    shift = 0
    if(single_tile .or. four_tile .or. folded_north_nonsym) then
       shift = 0
    else
       shift = 1
    endif

    allocate(global1(1-whalo:nx+shift+ehalo,1-shalo:ny+shift+nhalo,nz,ntile_per_pe) ) 
    allocate(global2(1-whalo:nx+shift+ehalo,1-shalo:ny+shift+nhalo,nz,ntile_per_pe) ) 
    allocate(global1_all(nx+shift,ny+shift,nz, ntiles),  global2_all(nx+shift,ny+shift,nz, ntiles))    
    global1 = 0; global2 = 0
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny+shift
             do i = 1, nx+shift
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    !-----------------------------------------------------------------------
    !--- make sure consistency on the boundary for cubic grid
    !--- east boundary will take the value of neighbor tile ( west/south),
    !--- north boundary will take the value of neighbor tile ( south/west).
    !--- for the point on the corner, the 12 corner take the following value
    !--- corner between 1, 2, 3 takes the value at 3, 
    !--- corner between 1, 3, 5 takes the value at 3
    !-----------------------------------------------------------------------
    if( cubic_grid ) then
       do l = 1, ntiles
          if(mod(l,2) == 0) then ! tile 2, 4, 6
             te = l + 2
             tn = l + 1
             if(te>6) te = te - 6
             if(tn > 6) tn = tn - 6
             global1_all(nx+shift,1:ny+1,:,l) = global2_all(nx+shift:1:-1,1,:,te)  ! east 
             global2_all(nx+shift,1:ny+1,:,l) = global1_all(nx+shift:1:-1,1,:,te)  ! east 
             global1_all(1:nx,ny+shift,:,l)    = global1_all(1:nx,1,:,tn) ! north
             global2_all(1:nx,ny+shift,:,l)    = global2_all(1:nx,1,:,tn) ! north
          else                   ! tile 1, 3, 5
             te = l + 1
             tn = l + 2
             if(tn > 6) tn = tn - 6
             global1_all(nx+shift,:,:,l)    = global1_all(1,:,:,te)  ! east
             global2_all(nx+shift,:,:,l)    = global2_all(1,:,:,te)  ! east
             global1_all(1:nx+1,ny+shift,:,l) = global2_all(1,ny+shift:1:-1,:,tn) ! north
             global2_all(1:nx+1,ny+shift,:,l) = global1_all(1,ny+shift:1:-1,:,tn) ! north
          end if
       end do
       ! set the corner value to 0 
       global1_all(1,ny+1,:,:) = 0; global1_all(nx+1,1,:,:) = 0; global1_all(1,1,:,:) = 0; global1_all(nx+1,ny+1,:,:) = 0
       global2_all(1,ny+1,:,:) = 0; global2_all(nx+1,1,:,:) = 0; global2_all(1,1,:,:) = 0; global2_all(nx+1,ny+1,:,:) = 0
    end if

    do n = 1, ntile_per_pe
       global1(1:nx+shift,1:ny+shift,:,n) = global1_all(:,:,:,tile(n))
       global2(1:nx+shift,1:ny+shift,:,n) = global2_all(:,:,:,tile(n))
    end do

    if(folded_north) then
       call fill_folded_north_halo(global1(:,:,:,1), 1, 1, shift, shift, -1)    
       call fill_folded_north_halo(global2(:,:,:,1), 1, 1, shift, shift, -1)   
    else if(folded_south_sym) then
       call fill_folded_south_halo(global1(:,:,:,1), 1, 1, shift, shift, -1)    
       call fill_folded_south_halo(global2(:,:,:,1), 1, 1, shift, shift, -1)  
    else if(folded_west_sym) then
       call fill_folded_west_halo(global1(:,:,:,1), 1, 1, shift, shift, -1)    
       call fill_folded_west_halo(global2(:,:,:,1), 1, 1, shift, shift, -1) 
    else if(folded_east_sym) then
       call fill_folded_east_halo(global1(:,:,:,1), 1, 1, shift, shift, -1)    
       call fill_folded_east_halo(global2(:,:,:,1), 1, 1, shift, shift, -1) 
    endif

    allocate( x (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x1(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x2(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x3(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x4(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y1(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y2(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y3(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y4(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )

    x = 0.; y = 0
    x (isc:iec+shift,jsc:jec+shift,:,:) = global1(isc:iec+shift,jsc:jec+shift,:,:)
    y (isc:iec+shift,jsc:jec+shift,:,:) = global2(isc:iec+shift,jsc:jec+shift,:,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points.     
    !-----------------------------------------------------------------------
    if(cubic_grid) then
       type2 = type//' paired-scalar BGRID_NE'
       update_flags = SCALAR_PAIR
    else
       type2 = type//' vector BGRID_NE'
       update_flags = XUPDATE + YUPDATE
    endif

    id = mpp_clock_id( trim(type2), flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    type3 = type2

    do n = 1, ntile_per_pe
       if(single_tile) then
          call fill_regular_mosaic_halo(global1(:,:,:,n), global1_all, 1, 1, 1, 1, 1, 1, 1, 1)       
          call fill_regular_mosaic_halo(global2(:,:,:,n), global2_all, 1, 1, 1, 1, 1, 1, 1, 1)     
       else if(folded_north) then
          !redundant points must be equal and opposite for tripolar grid
          global1(nx/2+shift,                ny+shift,:,:) = 0.  !pole points must have 0 velocity
          global1(nx+shift  ,                ny+shift,:,:) = 0.  !pole points must have 0 velocity
          global1(nx/2+1+shift:nx-1+shift,   ny+shift,:,:) = -global1(nx/2-1+shift:1+shift:-1, ny+shift,:,:)
          global1(1-whalo:shift,             ny+shift,:,:) = -global1(nx-whalo+1:nx+shift,     ny+shift,:,:)
          global1(nx+1+shift:nx+ehalo+shift, ny+shift,:,:) = -global1(1+shift:ehalo+shift,     ny+shift,:,:)
          global2(nx/2+shift,                ny+shift,:,:) = 0.  !pole points must have 0 velocity
          global2(nx+shift  ,                ny+shift,:,:) = 0.  !pole points must have 0 velocity
          global2(nx/2+1+shift:nx-1+shift,   ny+shift,:,:) = -global2(nx/2-1+shift:1+shift:-1, ny+shift,:,:)
          global2(1-whalo:shift,             ny+shift,:,:) = -global2(nx-whalo+1:nx+shift,     ny+shift,:,:)
          global2(nx+1+shift:nx+ehalo+shift, ny+shift,:,:) = -global2(1+shift:ehalo+shift,     ny+shift,:,:)
          !--- the following will fix the +0/-0 problem on altix
          if(nhalo >0) then
             global1(shift,ny+shift,:,:) = 0.  !pole points must have 0 velocity
             global2(shift,ny+shift,:,:) = 0.  !pole points must have 0 velocity
          end if
       else if(folded_south_sym) then
          global1(nx/2+shift,                1,:,:) = 0.  !pole points must have 0 velocity
          global1(nx+shift  ,                1,:,:) = 0.  !pole points must have 0 velocity
          global1(nx/2+1+shift:nx-1+shift,   1,:,:) = -global1(nx/2-1+shift:1+shift:-1, 1,:,:)
          global1(1-whalo:shift,             1,:,:) = -global1(nx-whalo+1:nx+shift,     1,:,:)
          global1(nx+1+shift:nx+ehalo+shift, 1,:,:) = -global1(1+shift:ehalo+shift,     1,:,:)
          global2(nx/2+shift,                1,:,:) = 0.  !pole points must have 0 velocity
          global2(nx+shift  ,                1,:,:) = 0.  !pole points must have 0 velocity
          global2(nx/2+1+shift:nx-1+shift,   1,:,:) = -global2(nx/2-1+shift:1+shift:-1, 1,:,:)
          global2(1-whalo:shift,             1,:,:) = -global2(nx-whalo+1:nx+shift,     1,:,:)
          global2(nx+1+shift:nx+ehalo+shift, 1,:,:) = -global2(1+shift:ehalo+shift,     1,:,:)
          !--- the following will fix the +0/-0 problem on altix
          if(shalo >0) then
             global1(shift,1,:,:) = 0.  !pole points must have 0 velocity
             global2(shift,1,:,:) = 0.  !pole points must have 0 velocity
          endif
       else if(folded_west_sym) then
          global1(1, ny/2+shift, :,:) = 0. !pole points must have 0 velocity
          global1(1, ny+shift,   :,:) = 0. !pole points must have 0 velocity
          global1(1, ny/2+1+shift:ny-1+shift,   :,:) = -global1(1, ny/2-1+shift:1+shift:-1, :,:)
          global1(1, 1-shalo:shift,             :,:) = -global1(1, ny-shalo+1:ny+shift,     :,:)
          global1(1, ny+1+shift:ny+nhalo+shift, :,:) = -global1(1, 1+shift:nhalo+shift,     :,:)
          global2(1, ny/2+shift, :,:) = 0. !pole points must have 0 velocity
          global2(1, ny+shift,   :,:) = 0. !pole points must have 0 velocity
          global2(1, ny/2+1+shift:ny-1+shift,   :,:) = -global2(1, ny/2-1+shift:1+shift:-1, :,:)
          global2(1, 1-shalo:shift,             :,:) = -global2(1, ny-shalo+1:ny+shift,     :,:)
          global2(1, ny+1+shift:ny+nhalo+shift, :,:) = -global2(1, 1+shift:nhalo+shift,     :,:)
       else if(folded_east_sym) then
          global1(nx+shift, ny/2+shift, :,:) = 0. !pole points must have 0 velocity
          global1(nx+shift, ny+shift,   :,:) = 0. !pole points must have 0 velocity
          global1(nx+shift, ny/2+1+shift:ny-1+shift,   :,:) = -global1(nx+shift, ny/2-1+shift:1+shift:-1, :,:)
          global1(nx+shift, 1-shalo:shift,             :,:) = -global1(nx+shift, ny-shalo+1:ny+shift,     :,:)
          global1(nx+shift, ny+1+shift:ny+nhalo+shift, :,:) = -global1(nx+shift, 1+shift:nhalo+shift,     :,:)
          global2(nx+shift, ny/2+shift, :,:) = 0. !pole points must have 0 velocity
          global2(nx+shift, ny+shift,   :,:) = 0. !pole points must have 0 velocity
          global2(nx+shift, ny/2+1+shift:ny-1+shift,   :,:) = -global2(nx+shift, ny/2-1+shift:1+shift:-1, :,:)
          global2(nx+shift, 1-shalo:shift,             :,:) = -global2(nx+shift, ny-shalo+1:ny+shift,     :,:)
          global2(nx+shift, ny+1+shift:ny+nhalo+shift, :,:) = -global2(nx+shift, 1+shift:nhalo+shift,     :,:)
       else if(four_tile) then
          select case ( tile(n) )
          case (1)
             tw = 2; ts = 3; tsw = 4
          case (2)
             tw = 1; ts = 4; tsw = 3
          case (3)
             tw = 4; ts = 1; tsw = 2
          case (4)
             tw = 3; ts = 2; tsw = 1
          end select
          te = tw; tn = ts; tse = tsw; tnw = tsw; tne = tsw
          call fill_regular_mosaic_halo(global1(:,:,:,n), global1_all, te, tse, ts, tsw, tw, tnw, tn, tne )
          call fill_regular_mosaic_halo(global2(:,:,:,n), global2_all, te, tse, ts, tsw, tw, tnw, tn, tne )
       else if(cubic_grid) then
          call fill_cubic_grid_halo(global1(:,:,:,n), global1_all, global2_all, tile(n), 1, 1, 1, 1 )
          call fill_cubic_grid_halo(global2(:,:,:,n), global2_all, global1_all, tile(n), 1, 1, 1, 1 )
       endif

       if(ntile_per_pe > 1) write(type3, *)trim(type2), " at tile_count = ",n
       call mpp_clock_begin(id)
       if(ntile_per_pe == 1) then
          call mpp_update_domains( x(:,:,:,n),  y(:,:,:,n),  domain, flags=update_flags, gridtype=BGRID_NE, name=type3)
       else
          call mpp_update_domains( x(:,:,:,n),  y(:,:,:,n),  domain, flags=update_flags, gridtype=BGRID_NE, &
               name=type3, tile_count = n)
       end if
       call mpp_clock_end  (id)
    end do

    do n = 1, ntile_per_pe
       if(ntile_per_pe > 1) write(type3, *)trim(type2), " at tile_count = ", n
       call compare_checksums( x (isd:ied+shift,jsd:jed+shift,:,n),  global1(isd:ied+shift,jsd:jed+shift,:,n), trim(type3)//' X' )
       call compare_checksums( y (isd:ied+shift,jsd:jed+shift,:,n),  global2(isd:ied+shift,jsd:jed+shift,:,n), trim(type3)//' Y' )
    end do

    if(ntile_per_pe == 1) then
       call mpp_clock_begin(id)
       call mpp_update_domains( x1, y1, domain, flags=update_flags, gridtype=BGRID_NE, complete=.false., name=type2)
       call mpp_update_domains( x2, y2, domain, flags=update_flags, gridtype=BGRID_NE, complete=.false.,  name=type2)
       call mpp_update_domains( x3, y3, domain, flags=update_flags, gridtype=BGRID_NE, complete=.false., name=type2)
       call mpp_update_domains( x4, y4, domain, flags=update_flags, gridtype=BGRID_NE, complete=.true.,  name=type2)
       call mpp_clock_end  (id)

       call compare_checksums( x1(isd:ied+shift,jsd:jed+shift,:,1), global1(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' X1')
       call compare_checksums( x2(isd:ied+shift,jsd:jed+shift,:,1), global1(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' X2')
       call compare_checksums( x3(isd:ied+shift,jsd:jed+shift,:,1), global1(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' X3')
       call compare_checksums( x4(isd:ied+shift,jsd:jed+shift,:,1), global1(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' X4')
       call compare_checksums( y1(isd:ied+shift,jsd:jed+shift,:,1), global2(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' Y1')
       call compare_checksums( y2(isd:ied+shift,jsd:jed+shift,:,1), global2(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' Y2')
       call compare_checksums( y3(isd:ied+shift,jsd:jed+shift,:,1), global2(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' Y3')
       call compare_checksums( y4(isd:ied+shift,jsd:jed+shift,:,1), global2(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' Y4')

       !--- arbitrary halo updates ---------------------------------------
       if(single_tile .or. four_tile .or. cubic_grid ) then
          allocate(local1(isd:ied+shift,jsd:jed+shift,nz) )     
          allocate(local2(isd:ied+shift,jsd:jed+shift,nz) )    
          do wh = 1-whalo, whalo
             do eh = 1-ehalo, ehalo
                do sh = 1-shalo, shalo
                   do nh = 1-nhalo, nhalo
                      if( wh*eh <= 0 ) cycle
                      if( sh*nh <= 0 ) cycle
                      if( wh*sh <= 0 ) cycle

                      local1(isd:ied+shift,jsd:jed+shift,:) = global1(isd:ied+shift,jsd:jed+shift,:,1)
                      local2(isd:ied+shift,jsd:jed+shift,:) = global2(isd:ied+shift,jsd:jed+shift,:,1)
                      x = 0.; y = 0.
                      x(isc:iec+shift,jsc:jec+shift,:,1) = global1_all(isc:iec+shift,jsc:jec+shift,:,tile(1))       
                      y(isc:iec+shift,jsc:jec+shift,:,1) = global2_all(isc:iec+shift,jsc:jec+shift,:,tile(1))    
                      call fill_halo_zero(local1, wh, eh, sh, nh, shift, shift, isc, iec, jsc, jec, isd, ied, jsd, jed)  
                      call fill_halo_zero(local2, wh, eh, sh, nh, shift, shift, isc, iec, jsc, jec, isd, ied, jsd, jed) 

                      write(type3,'(a,a,i2,a,i2,a,i2,a,i2)') trim(type2), ' with whalo = ', wh, &
                           ', ehalo = ',eh, ', shalo = ', sh, ', nhalo = ', nh
                      call mpp_update_domains( x,  y,  domain, flags=update_flags, gridtype=BGRID_NE, &
                           whalo=wh, ehalo=eh, shalo=sh, nhalo=nh, name=type3)
                      call compare_checksums( x(isd:ied+shift,jsd:jed+shift,:,1),  local1, trim(type3)//' X' )
                      call compare_checksums( y(isd:ied+shift,jsd:jed+shift,:,1),  local2, trim(type3)//' Y' )
                   end do
                end do
             end do
          end do
          deallocate(local1, local2)
       end if
    end if
    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- setup data
    if(cubic_grid .or. folded_north .or. folded_south_sym .or. folded_west_sym .or. folded_east_sym ) then
       deallocate(global1_all, global2_all)
       allocate(global1_all(nx+shift,ny,nz, ntiles),  global2_all(nx,ny+shift,nz, ntiles))   
       deallocate(global1, global2, x, y, x1, x2, x3, x4, y1, y2, y3, y4)
       allocate(global1(1-whalo:nx+shift+ehalo,1-shalo:ny  +nhalo,nz,ntile_per_pe) ) 
       allocate( x (ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( y (ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate( x1(ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( x2(ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( x3(ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( x4(ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( y1(ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate( y2(ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate( y3(ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate( y4(ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate(global2(1-whalo:nx  +ehalo,1-shalo:ny+shift+nhalo,nz,ntile_per_pe) ) 
       global1 = 0; global2 = 0
       do l = 1, ntiles
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx+shift
                   global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                end do
             end do
             do j = 1, ny+shift
                do i = 1, nx
                   global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                end do
             end do
          end do
       end do
    endif
    if( folded_north .or. folded_south_sym .or. folded_west_sym .or. folded_east_sym ) then
       do n = 1, ntile_per_pe
          global1(1:nx+shift,1:ny  ,:,n) = global1_all(1:nx+shift,1:ny,  :,tile(n))
          global2(1:nx  ,1:ny+shift,:,n) = global2_all(1:nx  ,1:ny+shift,:,tile(n))
       end do
    endif

    if( cubic_grid ) then
       !-----------------------------------------------------------------------
       !--- make sure consistency on the boundary for cubic grid
       !--- east boundary will take the value of neighbor tile ( west/south),
       !--- north boundary will take the value of neighbor tile ( south/west).
       !-----------------------------------------------------------------------
       do l = 1, ntiles
          if(mod(l,2) == 0) then ! tile 2, 4, 6
             te = l + 2
             tn = l + 1
             if(te>6) te = te - 6
             if(tn > 6) tn = tn - 6
             global1_all(nx+shift,1:ny,:,l) = global2_all(nx:1:-1,1,:,te)  ! east 
             global2_all(1:nx,ny+shift,:,l) = global2_all(1:nx,1,:,tn) ! north
          else                   ! tile 1, 3, 5
             te = l + 1
             tn = l + 2
             if(tn > 6) tn = tn - 6
             global1_all(nx+shift,:,:,l)    = global1_all(1,:,:,te)  ! east
             global2_all(1:nx,ny+shift,:,l) = global1_all(1,ny:1:-1,:,tn) ! north
          end if
       end do
       do n = 1, ntile_per_pe
          global1(1:nx+shift,1:ny  ,:,n) = global1_all(1:nx+shift,1:ny,  :,tile(n))
          global2(1:nx  ,1:ny+shift,:,n) = global2_all(1:nx  ,1:ny+shift,:,tile(n))
       end do
    else if( folded_north ) then
       call fill_folded_north_halo(global1(:,:,:,1), 1, 0, shift, 0, -1)
       call fill_folded_north_halo(global2(:,:,:,1), 0, 1, 0, shift, -1)
    else if(folded_south_sym ) then
       call fill_folded_south_halo(global1(:,:,:,1), 1, 0, shift, 0, -1)
       call fill_folded_south_halo(global2(:,:,:,1), 0, 1, 0, shift, -1)
    else if(folded_west_sym ) then
       call fill_folded_west_halo(global1(:,:,:,1), 1, 0, shift, 0, -1)
       call fill_folded_west_halo(global2(:,:,:,1), 0, 1, 0, shift, -1)
    else if(folded_east_sym ) then
       call fill_folded_east_halo(global1(:,:,:,1), 1, 0, shift, 0, -1)
       call fill_folded_east_halo(global2(:,:,:,1), 0, 1, 0, shift, -1)
    endif
    x = 0.; y = 0.
    x (isc:iec+shift,jsc:jec  ,:,:) = global1(isc:iec+shift,jsc:jec  ,:,:)
    y (isc:iec  ,jsc:jec+shift,:,:) = global2(isc:iec  ,jsc:jec+shift,:,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points for cubic-grid.
    !   On the contact line, the following relation will be used to 
    !   --- fill the value on contact line ( balance send and recv).
    !       2W --> 1E, 1S --> 6N, 3W --> 1N, 4S --> 2E
    !       4W --> 3E, 3S --> 2N, 1W --> 5N, 2S --> 6E
    !       6W --> 5E, 5S --> 4N, 5W --> 3N, 6S --> 4E
    !---------------------------------------------------------------------------
    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    type2 = type
    do n = 1, ntile_per_pe
       if( cubic_grid ) then    
          call fill_cubic_grid_halo(global1(:,:,:,n), global1_all, global2_all, tile(n), 1, 0, 1, -1 )
          call fill_cubic_grid_halo(global2(:,:,:,n), global2_all, global1_all, tile(n), 0, 1, -1, 1 )
       else if( folded_north ) then
          !redundant points must be equal and opposite
          global2(nx/2+1:nx,     ny+shift,:,:) = -global2(nx/2:1:-1, ny+shift,:,:)
          global2(1-whalo:0,     ny+shift,:,:) = -global2(nx-whalo+1:nx, ny+shift,:,:)
          global2(nx+1:nx+ehalo, ny+shift,:,:) = -global2(1:ehalo,       ny+shift,:,:)
       else if( folded_south_sym ) then
          global2(nx/2+1:nx,     1,:,:) = -global2(nx/2:1:-1, 1,:,:)
          global2(1-whalo:0,     1,:,:) = -global2(nx-whalo+1:nx, 1, :,:)
          global2(nx+1:nx+ehalo, 1,:,:) = -global2(1:ehalo,       1, :,:)
       else if( folded_west_sym ) then
          global1(1, ny/2+1:ny,     :,:) = -global1(1, ny/2:1:-1,     :,:)
          global1(1, 1-shalo:0,     :,:) = -global1(1, ny-shalo+1:ny, :,:)
          global1(1, ny+1:ny+nhalo, :,:) = -global1(1, 1:nhalo,       :,:)
       else if( folded_east_sym ) then
          global1(nx+shift, ny/2+1:ny,     :,:) = -global1(nx+shift, ny/2:1:-1,     :,:)
          global1(nx+shift, 1-shalo:0,     :,:) = -global1(nx+shift, ny-shalo+1:ny, :,:)
          global1(nx+shift, ny+1:ny+nhalo, :,:) = -global1(nx+shift, 1:nhalo,       :,:)
       end if

       if(ntile_per_pe > 1) write(type2, *)type, " at tile_count = ",n
       call mpp_clock_begin(id)
       if(ntile_per_pe == 1) then
          call mpp_update_domains( x(:,:,:,n),  y(:,:,:,n),  domain, gridtype=CGRID_NE, name=type2//' vector CGRID_NE')
       else
          call mpp_update_domains( x(:,:,:,n),  y(:,:,:,n),  domain, gridtype=CGRID_NE, &
               name=type2//' vector CGRID_NE', tile_count = n)
       end if
       call mpp_clock_end  (id)
    end do



    do n = 1, ntile_per_pe
       if(ntile_per_pe > 1) write(type2, *)type, " at tile_count = ",n
       call compare_checksums( x(isd:ied+shift,jsd:jed,:,n), global1(isd:ied+shift,jsd:jed,  :,n), &
                               trim(type2)//' CGRID_NE X')
       call compare_checksums( y(isd:ied,jsd:jed+shift,:,n), global2(isd:ied,  jsd:jed+shift,:,n), &
                               trim(type2)//' CGRID_NE Y')       
    end do

    if(ntile_per_pe == 1) then
       call mpp_clock_begin(id)
       call mpp_update_domains( x1, y1, domain, gridtype=CGRID_NE, complete=.false., name=type//' vector CGRID_NE' )
       call mpp_update_domains( x2, y2, domain, gridtype=CGRID_NE, complete=.false., name=type//' vector CGRID_NE')
       call mpp_update_domains( x3, y3, domain, gridtype=CGRID_NE, complete=.false., name=type//' vector CGRID_NE' )
       call mpp_update_domains( x4, y4, domain, gridtype=CGRID_NE, complete=.true. , name=type//' vector CGRID_NE')
       call mpp_clock_end  (id)

       call compare_checksums( x1(isd:ied+shift,jsd:jed,:,1), global1(isd:ied+shift,jsd:jed,:,1), type//' CGRID_NE X1')
       call compare_checksums( x2(isd:ied+shift,jsd:jed,:,1), global1(isd:ied+shift,jsd:jed,:,1), type//' CGRID_NE X2')
       call compare_checksums( x3(isd:ied+shift,jsd:jed,:,1), global1(isd:ied+shift,jsd:jed,:,1), type//' CGRID_NE X3')
       call compare_checksums( x4(isd:ied+shift,jsd:jed,:,1), global1(isd:ied+shift,jsd:jed,:,1), type//' CGRID_NE X4')
       call compare_checksums( y1(isd:ied,jsd:jed+shift,:,1), global2(isd:ied,jsd:jed+shift,:,1), type//' CGRID_NE Y1')
       call compare_checksums( y2(isd:ied,jsd:jed+shift,:,1), global2(isd:ied,jsd:jed+shift,:,1), type//' CGRID_NE Y2')
       call compare_checksums( y3(isd:ied,jsd:jed+shift,:,1), global2(isd:ied,jsd:jed+shift,:,1), type//' CGRID_NE Y3')
       call compare_checksums( y4(isd:ied,jsd:jed+shift,:,1), global2(isd:ied,jsd:jed+shift,:,1), type//' CGRID_NE Y4')

       !--- arbitrary halo updates ---------------------------------------
       if(single_tile .or. four_tile .or. cubic_grid ) then
          allocate(local1(isd:ied+shift,jsd:jed,      nz) )     
          allocate(local2(isd:ied,      jsd:jed+shift,nz) )    

          do wh = 1-whalo, whalo
             do eh = 1-ehalo, ehalo
                do sh = 1-shalo, shalo
                   do nh = 1-nhalo, nhalo
                      if( wh*eh <= 0 ) cycle
                      if( sh*nh <= 0 ) cycle
                      if( wh*sh <= 0 ) cycle
                      local1(isd:ied+shift,jsd:jed,      :) = global1(isd:ied+shift,jsd:jed,      :,1)
                      local2(isd:ied,      jsd:jed+shift,:) = global2(isd:ied,      jsd:jed+shift,:,1)
                      x = 0.; y = 0.
                      x(isc:iec+shift,jsc:jec,      :,1) = global1_all(isc:iec+shift,jsc:jec,      :,tile(1))       
                      y(isc:iec,      jsc:jec+shift,:,1) = global2_all(isc:iec,      jsc:jec+shift,:,tile(1))    
                      call fill_halo_zero(local1, wh, eh, sh, nh, shift, 0, isc, iec, jsc, jec, isd, ied, jsd, jed)  
                      call fill_halo_zero(local2, wh, eh, sh, nh, 0, shift, isc, iec, jsc, jec, isd, ied, jsd, jed) 

                      write(type3,'(a,a,i2,a,i2,a,i2,a,i2)') trim(type), ' vector CGRID_NE with whalo = ', &
                           wh, ', ehalo = ',eh, ', shalo = ', sh, ', nhalo = ', nh
                      !          id = mpp_clock_id( trim(type3), flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
                      !          call mpp_clock_begin(id)
                      call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE, whalo=wh, ehalo=eh, &
                           shalo=sh, nhalo=nh, name=type3)
                      !          call mpp_clock_end  (id)
                      call compare_checksums( x(isd:ied+shift,jsd:jed, :,1),  local1, trim(type3)//' X' )
                      call compare_checksums( y(isd:ied,jsd:jed+shift, :,1),  local2, trim(type3)//' Y' )
                   end do
                end do
             end do
          end do
          deallocate(local1, local2)
       end if
    end if

    deallocate(global1, global2, x, y, x1, x2, x3, x4, y1, y2, y3, y4, global1_all, global2_all)
    deallocate(layout2D, global_indices, pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2 ) 

  end subroutine test_uniform_mosaic

  !#################################################################################

  subroutine fill_halo_zero(data, whalo, ehalo, shalo, nhalo, xshift, yshift, isc, iec, jsc, jec, isd, ied, jsd, jed)
    integer,                         intent(in) :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer,                         intent(in) :: whalo, ehalo, shalo, nhalo, xshift, yshift
    real, dimension(isd:,jsd:,:), intent(inout) :: data

    if(whalo >=0) then
       data(iec+ehalo+1+xshift:ied+xshift,jsd:jed+yshift,:) = 0
       data(isd:isc-whalo-1,jsd:jed+yshift,:) = 0
    else
       data(iec+1+xshift:iec-ehalo+xshift,jsc+shalo:jec-nhalo+yshift,:) = 0
       data(isc+whalo:isc-1,jsc+shalo:jec-nhalo+yshift,:) = 0
    end if

    if(shalo>=0) then
       data(isd:ied+xshift, jec+nhalo+1+yshift:jed+yshift,:) = 0
       data(isd:ied+xshift, jsd:jsc-shalo-1,:) = 0
    else
       data(isc+whalo:iec-ehalo+xshift,jec+1+yshift:jec-nhalo+yshift,:) = 0
       data(isc+whalo:iec-ehalo+xshift,jsc+shalo:jsc-1,:) = 0
    end if

  end subroutine fill_halo_zero

  !##############################################################################
  ! this routine fill the halo points for the regular mosaic. 
  subroutine fill_regular_mosaic_halo(data, data_all, te, tse, ts, tsw, tw, tnw, tn, tne)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data_all
    integer,                              intent(in)    :: te, tse, ts, tsw, tw, tnw, tn, tne

       data(nx+1:nx+ehalo, 1:ny,          :) = data_all(1:ehalo,       1:ny,          :, te) ! east
       data(1:nx,          1-shalo:0,     :) = data_all(1:nx,          ny-shalo+1:ny, :, ts) ! south 
       data(1-whalo:0,     1:ny,          :) = data_all(nx-whalo+1:nx, 1:ny,          :, tw) ! west
       data(1:nx,          ny+1:ny+nhalo, :) = data_all(1:nx,          1:nhalo,       :, tn) ! north  
       data(nx+1:nx+ehalo, 1-shalo:0,     :) = data_all(1:ehalo,       ny-shalo+1:ny, :,tse) ! southeast
       data(1-whalo:0,     1-shalo:0,     :) = data_all(nx-whalo+1:nx, ny-shalo+1:ny, :,tsw) ! southwest
       data(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = data_all(1:ehalo,       1:nhalo,       :,tnw) ! northeast
       data(1-whalo:0,     ny+1:ny+nhalo, :) = data_all(nx-whalo+1:nx, 1:nhalo,       :,tne) ! northwest    



  end subroutine fill_regular_mosaic_halo

  !################################################################################
  subroutine fill_folded_north_halo(data, ioff, joff, ishift, jshift, sign)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer,                              intent(in   ) :: ioff, joff, ishift, jshift, sign    
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,                  1:nyp,:) =      data(nx-whalo+1:nx,        1:ny+jshift,:) ! west
    data(nx+1:nx+ehalo+ishift,       1:nyp,:) =      data(1:ehalo+ishift,       1:ny+jshift,:) ! east
    if(m1 .GE. 1-whalo) data(1-whalo:m1,  nyp+1:nyp+nhalo,:) = sign*data(whalo+m2:1+ishift:-1, nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(m1+1:nx+m2,       nyp+1:nyp+nhalo,:) = sign*data(nx+ishift:1:-1,       nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(nx+m2+1:nxp+ehalo,nyp+1:nyp+nhalo,:) = sign*data(nx:nx-ehalo+m1+1:-1,  nyp-joff:nyp-nhalo-joff+1:-1,:)

  end subroutine fill_folded_north_halo

  !################################################################################
  subroutine fill_folded_south_halo(data, ioff, joff, ishift, jshift, sign)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer,                              intent(in   ) :: ioff, joff, ishift, jshift, sign    
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff


    data(1-whalo:0,                  1:nyp,:) =      data(nx-whalo+1:nx,        1:nyp,:) ! west
    data(nx+1:nx+ehalo+ishift,       1:nyp,:) =      data(1:ehalo+ishift,       1:nyp,:) ! east
    if(m1 .GE. 1-whalo)data(1-whalo:m1, 1-shalo:0,:) = sign*data(whalo+m2:1+ishift:-1, shalo+jshift:1+jshift:-1,:)
    data(m1+1:nx+m2,       1-shalo:0,:) = sign*data(nxp:1:-1,             shalo+jshift:1+jshift:-1,:)
    data(nx+m2+1:nxp+ehalo,1-shalo:0,:) = sign*data(nx:nx-ehalo+m1+1:-1,  shalo+jshift:1+jshift:-1,:)

  end subroutine fill_folded_south_halo

  !################################################################################
  subroutine fill_folded_west_halo(data, ioff, joff, ishift, jshift, sign)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer,                              intent(in   ) :: ioff, joff, ishift, jshift, sign    
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0, :)      = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo, :) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) data(1-whalo:0, 1-shalo:m1, :) = sign*data(whalo+ishift:1+ishift:-1, shalo+m2:1+jshift:-1,:)
    data(1-whalo:0, m1+1:ny+m2, :) = sign*data(whalo+ishift:1+ishift:-1, nyp:1:-1, :)
    data(1-whalo:0, ny+m2+1:nyp+nhalo,:) = sign*data(whalo+ishift:1+ishift:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_west_halo

  !################################################################################
  subroutine fill_folded_east_halo(data, ioff, joff, ishift, jshift, sign)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer,                              intent(in   ) :: ioff, joff, ishift, jshift, sign    
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0, :)      = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo, :) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) data(nxp+1:nxp+ehalo, 1-shalo:m1, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, shalo+m2:1+jshift:-1,:)
    data(nxp+1:nxp+ehalo, m1+1:ny+m2, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, nyp:1:-1, :)
    data(nxp+1:nxp+ehalo, ny+m2+1:nyp+nhalo,:) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_east_halo

  !################################################################################
  subroutine fill_four_tile_bound(data_all, is, ie, js, je, ioff, joff, tile, &
                                   ebound, sbound, wbound, nbound )
    real, dimension(:,:,:,:),       intent(in)    :: data_all
    integer,                        intent(in)    :: is, ie, js, je
    integer,                        intent(in)    :: tile, ioff, joff
    real, dimension(:,:), optional, intent(inout) :: ebound, sbound, wbound, nbound
    integer                                       :: tw, te, ts, tn

    if(tile == 1 .OR. tile == 3) te = tile + 1
    if(tile == 2 .OR. tile == 4) te = tile - 1
    if(tile == 1 .OR. tile == 2) ts = tile + 2
    if(tile == 3 .OR. tile == 4) ts = tile - 2
    tw = te;   tn = ts
    if(present(ebound)) then
       if( ie == nx ) then
          ebound(:,:) = data_all(1, js:je+joff, :, te)
       else
          ebound(:,:) = data_all(ie+ioff, js:je+joff, :, tile)
       end if
    end if

    if(present(wbound)) then
       if( is == 1 ) then
          wbound(:,:) = data_all(nx+ioff, js:je+joff, :, tw)
       else
          wbound(:,:) = data_all(is, js:je+joff, :, tile)
       end if
    end if

    if(present(sbound)) then
       if( js == 1 ) then
          sbound(:,:) = data_all(is:ie+ioff, ny+joff, :, ts)
       else
          sbound(:,:) = data_all(is:ie+ioff, js, :, tile)
       end if
    end if

    if(present(nbound)) then
       if( je == ny ) then
          nbound(:,:) = data_all(is:ie+ioff, 1, :, tn)
       else
          nbound(:,:) = data_all(is:ie+ioff, je+joff, :, tile)
       end if
    end if

    return

  end subroutine fill_four_tile_bound

  !################################################################################
  subroutine fill_cubic_grid_bound(data1_all, data2_all, is, ie, js, je, ioff, joff, tile, sign1, sign2, &
                                   ebound, sbound, wbound, nbound )
    real, dimension(:,:,:,:),       intent(in)    :: data1_all, data2_all
    integer,                        intent(in)    :: is, ie, js, je
    integer,                        intent(in)    :: tile, ioff, joff, sign1, sign2
    real, dimension(:,:), optional, intent(inout) :: ebound, sbound, wbound, nbound
    integer                                       :: tw, te, ts, tn

    if(mod(tile,2) == 0) then ! tile 2, 4, 6
       tw = tile - 1; te = tile + 2; ts = tile - 2; tn = tile + 1
       if(te > 6 ) te = te - 6
       if(ts < 1 ) ts = ts + 6
       if(tn > 6 ) tn = tn - 6
       !--- East bound
       if(present(ebound)) then
          if(ie == nx) then                
             ebound(:,:) = sign1*data2_all(nx+joff-js+1:nx-je+1:-1,1,:,te)
          else
             ebound(:,:) = data1_all(ie+ioff, js:je+joff, :,tile)
          end if
       end if
       !--- South bound
       if(present(sbound)) then
          if(js == 1) then                
             sbound(:,:) = sign2*data2_all(nx+joff, ny+ioff-is+1:ny-ie+1:-1,:,ts)
          else
             sbound(:,:) = data1_all(is:ie+ioff, js, :,tile)
          end if
       end if

       !--- West bound
       if(present(wbound)) then
          if(is == 1) then
             wbound(:,:) = data1_all(nx+ioff, js:je+joff,:,tw)
          else
             wbound(:,:) = data1_all(is, js:je+joff,:,tile)
          end if
       end if

       !--- north bound
       if(present(nbound)) then
          if(je == ny) then                
             nbound(:,:) = data1_all(is:ie+ioff, 1,:,tn)
          else
             nbound(:,:) = data1_all(is:ie+ioff, je+joff, :,tile)
          end if
       end if
    else ! tile 1, 3, 5
       tw = tile - 2; te = tile + 1; ts = tile - 1; tn = tile + 2
       if(tw < 1 ) tw = tw + 6
       if(ts < 1 ) ts = ts + 6
       if(tn > 6 ) tn = tn - 6
       !--- East bound
       if(present(ebound)) then
          if(ie == nx) then                
             ebound(:,:) = data1_all(1, js:je+joff, :,te) 
          else
             ebound(:,:) = data1_all(ie+ioff, js:je+joff, :,tile)
          end if
       end if
       !--- South bound
       if(present(sbound)) then
          if(js == 1) then                
             sbound(:,:) = data1_all(is:ie+ioff,ny+joff,:,ts)
          else
             sbound(:,:) = data1_all(is:ie+ioff, js, :,tile)
          end if
       end if

       !--- West bound
       if(present(wbound)) then
          if(is == 1) then
             wbound(:,:) = sign1*data2_all(nx+joff-js+1:nx-je+1:-1,ny+ioff,:,tw)
          else
             wbound(:,:) = data1_all(is, js:je+joff,:,tile)
          end if
       end if

       !--- north bound
       if(present(nbound)) then
          if(je == ny) then                
             nbound(:,:) = sign2*data2_all(1, ny+ioff-is+1:ny-ie+1:-1,:,tn)
          else
             nbound(:,:) = data1_all(is:ie+ioff, je+joff, :,tile)
          end if
       end if

    end if

  end subroutine fill_cubic_grid_bound

  !##############################################################################
  ! this routine fill the halo points for the cubic grid. ioff and joff is used to distinguish
  ! T, C, E, or N-cell
  subroutine fill_cubic_grid_halo(data, data1_all, data2_all, tile, ioff, joff, sign1, sign2)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data1_all, data2_all
    integer,                              intent(in)    :: tile, ioff, joff, sign1, sign2 
    integer                                             :: lw, le, ls, ln

    if(mod(tile,2) == 0) then ! tile 2, 4, 6
       lw = tile - 1; le = tile + 2; ls = tile - 2; ln = tile + 1
       if(le > 6 ) le = le - 6
       if(ls < 1 ) ls = ls + 6
       if(ln > 6 ) ln = ln - 6
       data(1-whalo:0, 1:ny+joff, :) = data1_all(nx-whalo+1:nx, 1:ny+joff, :, lw) ! west 
       do i = 1, ehalo 
          data(nx+i+ioff, 1:ny+joff, :)    = sign1*data2_all(nx+joff:1:-1, i+ioff, :, le) ! east 
       end do
       do i = 1, shalo 
          data(1:nx+ioff, 1-i, :)     = sign2*data2_all(nx-i+1, ny+ioff:1:-1, :, ls) ! south 
       end do
       data(1:nx+ioff, ny+1+joff:ny+nhalo+joff, :) = data1_all(1:nx+ioff, 1+joff:nhalo+joff, :, ln) ! north
    else ! tile 1, 3, 5
       lw = tile - 2; le = tile + 1; ls = tile - 1; ln = tile + 2
       if(lw < 1 ) lw = lw + 6
       if(ls < 1 ) ls = ls + 6
       if(ln > 6 ) ln = ln - 6
       do i = 1, whalo 
          data(1-i, 1:ny+joff, :)     = sign1*data2_all(nx+joff:1:-1, ny-i+1, :, lw) ! west 
       end do
       data(nx+1+ioff:nx+ehalo+ioff, 1:ny+joff, :) = data1_all(1+ioff:ehalo+ioff, 1:ny+joff, :, le) ! east 
       data(1:nx+ioff, 1-shalo:0, :)     = data1_all(1:nx+ioff, ny-shalo+1:ny, :, ls) ! south 
       do i = 1, nhalo 
          data(1:nx+ioff, ny+i+joff, :)    = sign2*data2_all(i+joff, ny+ioff:1:-1, :, ln) ! north 
       end do
    end if

  end subroutine fill_cubic_grid_halo
    
   !#####################################################################
  subroutine test_nonuniform_mosaic( type )
    character(len=*), intent(in) :: type

    type(domain2D)               :: domain
    integer                      :: num_contact, ntiles, ntile_per_pe
    integer                      :: i, j, k, n, nxm, nym, ni, nj, shift
    integer                      :: ism, iem, jsm, jem, isc, iec, jsc, jec
    integer                      :: isd, ied, jsd, jed
    integer                      :: indices(4), msize(2)
    character(len=128)           :: type2

    integer, allocatable, dimension(:)       :: tile
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:,:,:) :: global1_all, global2_all
    real,    allocatable, dimension(:,:,:,:) :: global1, global2, x, y  

    shift = 0
    select case(type)
    case('Five-Tile') ! one tile will run on pe 0 and other four tiles will run on pe 1
       shift = 1      ! one extra point for symmetry domain
       ntiles = 5     ! tile 1 with resolution 2*nx and 2*ny and the tiles are nx and ny.
       num_contact = 11
       if(npes .NE. 2) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: Five-Tile mosaic will not be tested because npes is not 2')
          return
       end if 
       nxm = 2*nx; nym = 2*ny
       layout = 1
       if( pe == 0) then
          ntile_per_pe = 1
          allocate(tile(ntile_per_pe))
          tile = 1
          indices = (/1,2*nx,1,2*ny/)
          ni = 2*nx; nj = 2*ny
       else
          ntile_per_pe = 4
          allocate(tile(ntile_per_pe))
          do n = 1, ntile_per_pe
             tile(n) = n + 1
          end do
          indices = (/1,nx,1,ny/)
          ni = nx; nj = ny
       end if
       allocate(pe_start(ntiles), pe_end(ntiles) )
       pe_start(1) = 0; pe_start(2:) = 1
       pe_end = pe_start
    case default
       call mpp_error(FATAL, 'TEST_MPP_DOMAINS: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles) )

    do n = 1, ntiles
       if(n==1) then
          global_indices(:,n) = (/1,2*nx,1,2*ny/)
       else
          global_indices(:,n) = (/1,nx,1,ny/)
       endif  
!       global_indices(:,n) = indices
       layout2D(:,n)       = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) ) 
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) ) 

    !--- define domain
    select case(type)
    case( 'Five-Tile' )
       !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
       tile1(1) = 1; tile2(1) = 2
       istart1(1) = 2*nx; iend1(1) = 2*nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;    iend2(1) = 1;    jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (EAST) and tile 4 (WEST)
       tile1(2) = 1; tile2(2) = 4
       istart1(2) = 2*nx; iend1(2) = 2*nx; jstart1(2) = ny+1; jend1(2) = 2*ny
       istart2(2) = 1;    iend2(2) = 1;    jstart2(2) = 1;    jend2(2) = ny
       !--- Contact line 3, between tile 1 (SOUTH) and tile 1 (NORTH)
       tile1(3) = 1; tile2(3) = 1
       istart1(3) = 1; iend1(3) = 2*nx; jstart1(3) = 1;    jend1(3) = 1
       istart2(3) = 1; iend2(3) = 2*nx; jstart2(3) = 2*ny; jend2(3) = 2*ny
       !--- Contact line 4, between tile 1 (WEST) and tile 3 (EAST)
       tile1(4) = 1; tile2(4) = 3
       istart1(4) = 1;  iend1(4) = 1;  jstart1(4) = 1;  jend1(4) = ny
       istart2(4) = nx; iend2(4) = nx; jstart2(4) = 1;  jend2(4) = ny
       !--- Contact line 5, between tile 1 (WEST) and tile 5 (EAST)
       tile1(5) = 1; tile2(5) = 5
       istart1(5) = 1;  iend1(5) = 1;  jstart1(5) = ny+1;  jend1(5) = 2*ny
       istart2(5) = nx; iend2(5) = nx; jstart2(5) = 1;     jend2(5) = ny
       !--- Contact line 6, between tile 2 (EAST) and tile 3 (WEST)
       tile1(6) = 2; tile2(6) = 3
       istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
       istart2(6) = 1;  iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = ny       
       !--- Contact line 7, between tile 2 (SOUTH) and tile 4 (NORTH)  --- cyclic
       tile1(7) = 2; tile2(7) = 4
       istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;   jend1(7) = 1
       istart2(7) = 1;  iend2(7) = nx; jstart2(7) = ny;  jend2(7) = ny
       !--- Contact line 8, between tile 2 (NORTH) and tile 4 (SOUTH) 
       tile1(8) = 2; tile2(8) = 4
       istart1(8) = 1;  iend1(8) = nx; jstart1(8) = ny;  jend1(8) = ny
       istart2(8) = 1;  iend2(8) = nx; jstart2(8) = 1;   jend2(8) = 1
       !--- Contact line 9, between tile 3 (SOUTH) and tile 5 (NORTH)  --- cyclic
       tile1(9) = 3; tile2(9) = 5
       istart1(9) = 1;  iend1(9) = nx; jstart1(9) = 1;   jend1(9) = 1
       istart2(9) = 1;  iend2(9) = nx; jstart2(9) = ny;  jend2(9) = ny
       !--- Contact line 10, between tile 3 (NORTH) and tile 5 (SOUTH) 
       tile1(10) = 3; tile2(10) = 5
       istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny;  jend1(10) = ny
       istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;   jend2(10) = 1
       !--- Contact line 11, between tile 4 (EAST) and tile 5 (WEST)
       tile1(11) = 4; tile2(11) = 5
       istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
       istart2(11) = 1;  iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = ny  
       msize(1) = 2*nx + whalo + ehalo
       msize(2) = 2*ny + shalo + nhalo
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
            istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
            pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo,       &
            name = type, memory_size = msize, symmetry = .true.  )
    end select
    
    !--- setup data
    allocate(global1_all(1:nxm,1:nym,nz, ntiles) )  
    allocate(global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo,nz, ntile_per_pe) )   
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, nym
             do i = 1, nxm
                global1_all(i,j,k,n) = n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global1(1:ni,1:nj,:,n) = global1_all(1:ni,1:nj,:,tile(n))
    end do

    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )

    allocate( x (ism:iem,jsm:jem,nz, ntile_per_pe) )
    x = 0.
    x(isc:iec,jsc:jec,:,:) = global1(isc:iec,jsc:jec,:,:)

    !--- fill up the value at halo points
    do n = 1, ntile_per_pe
       call fill_five_tile_halo(global1(:,:,:,n), global1_all, tile(n), 0, 0 )
    end do

    ! full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    do n = 1, ntile_per_pe
       call mpp_update_domains( x(:,:,:,n), domain, tile_count = n )
    end do
    call mpp_clock_end(id)

   do n = 1, ntile_per_pe
      write(type2, *)type, " at tile_count = ",n
      call compare_checksums( x(isd:ied,jsd:jed,:,n), global1(isd:ied,jsd:jed,:,n), trim(type2) )
   end do

   deallocate(global1_all, global1, x)

    !------------------------------------------------------------------
    !  vector update : BGRID_NE, one extra point in each direction for Five-Tile
    !------------------------------------------------------------------
    !--- setup data
    allocate(global1_all(nxm+shift,nym+shift,nz, ntiles), global2_all(nxm+shift,nym+shift,nz, ntiles) )  
    allocate(global1(1-whalo:ni+ehalo+shift,1-shalo:nj+nhalo+shift,nz, ntile_per_pe) )  
    allocate(global2(1-whalo:ni+ehalo+shift,1-shalo:nj+nhalo+shift,nz, ntile_per_pe) )   
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, nym+shift
             do i = 1, nxm+shift
                global1_all(i,j,k,n) = 1.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,n) = 2.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    !------------------------------------------------------------------------
    ! --- make sure consisency on the boundary for Five-Tile mosaic
    ! --- east boundary will take the value of neighbor tile west,
    ! --- north boundary will take the value of neighbor tile south.
    !------------------------------------------------------------------------
    if(type == 'Five-Tile') then
       global1_all(nxm+1,    1:ny,:,1) = global1_all(1,    1:ny,:,2)  ! east
       global1_all(nxm+1,ny+1:nym,:,1) = global1_all(1,    1:ny,:,4)  ! east
       global1_all(1:nxm+1, nym+1,:,1) = global1_all(1:nxm+1, 1,:,1)  ! north
       global1_all(nx+1,     1:ny,:,2) = global1_all(1,    1:ny,:,3)  ! east
       global1_all(1:nx+1,   ny+1,:,2) = global1_all(1:nx+1,  1,:,4)  ! north
       global1_all(nx+1,     1:ny,:,3) = global1_all(1,    1:ny,:,1)  ! east
       global1_all(1:nx+1,   ny+1,:,3) = global1_all(1:nx+1,  1,:,5)  ! north
       global1_all(nx+1,     1:ny,:,4) = global1_all(1,    1:ny,:,5)  ! east
       global1_all(1:nx+1,   ny+1,:,4) = global1_all(1:nx+1,  1,:,2)  ! north
       global1_all(nx+1,     1:ny,:,5) = global1_all(1,ny+1:nym,:,1)  ! east
       global1_all(1:nx+1,   ny+1,:,5) = global1_all(1:nx+1,  1,:,3)  ! north
       global1_all(nx+1,     ny+1,:,2) = global1_all(1,       1,:,5)  ! northeast 
       global1_all(nx+1,     ny+1,:,3) = global1_all(1,    ny+1,:,1)  ! northeast 
       global2_all(nxm+1,    1:ny,:,1) = global2_all(1,    1:ny,:,2)  ! east
       global2_all(nxm+1,ny+1:nym,:,1) = global2_all(1,    1:ny,:,4)  ! east
       global2_all(1:nxm+1, nym+1,:,1) = global2_all(1:nxm+1, 1,:,1)  ! north
       global2_all(nx+1,     1:ny,:,2) = global2_all(1,    1:ny,:,3)  ! east
       global2_all(1:nx+1,   ny+1,:,2) = global2_all(1:nx+1,  1,:,4)  ! north
       global2_all(nx+1,     1:ny,:,3) = global2_all(1,    1:ny,:,1)  ! east
       global2_all(1:nx+1,   ny+1,:,3) = global2_all(1:nx+1,  1,:,5)  ! north
       global2_all(nx+1,     1:ny,:,4) = global2_all(1,    1:ny,:,5)  ! east
       global2_all(1:nx+1,   ny+1,:,4) = global2_all(1:nx+1,  1,:,2)  ! north
       global2_all(nx+1,     1:ny,:,5) = global2_all(1,ny+1:nym,:,1)  ! east
       global2_all(1:nx+1,   ny+1,:,5) = global2_all(1:nx+1,  1,:,3)  ! north
       global2_all(nx+1,     ny+1,:,2) = global2_all(1,       1,:,5)  ! northeast 
       global2_all(nx+1,     ny+1,:,3) = global2_all(1,    ny+1,:,1)  ! northeast 
    end if

    do n = 1, ntile_per_pe
       global1(1:ni+shift,1:nj+shift,:,n) = global1_all(1:ni+shift,1:nj+shift,:,tile(n))
       global2(1:ni+shift,1:nj+shift,:,n) = global2_all(1:ni+shift,1:nj+shift,:,tile(n))
    end do

    allocate( x (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )

    x = 0.; y = 0
    x (isc:iec+shift,jsc:jec+shift,:,:) = global1(isc:iec+shift,jsc:jec+shift,:,:)
    y (isc:iec+shift,jsc:jec+shift,:,:) = global2(isc:iec+shift,jsc:jec+shift,:,:)

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points.     
    !-----------------------------------------------------------------------
    do n = 1, ntile_per_pe
       call fill_five_tile_halo(global1(:,:,:,n), global1_all, tile(n), shift, shift)
       call fill_five_tile_halo(global2(:,:,:,n), global2_all, tile(n), shift, shift)
    end do

    id = mpp_clock_id( type//' BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    do n = 1, ntile_per_pe
       call mpp_update_domains( x(:,:,:,n), y(:,:,:,n), domain, gridtype=BGRID_NE, tile_count = n )
    end do
    call mpp_clock_end(id)

   do n = 1, ntile_per_pe
      write(type2, *)type, " at tile_count = ",n
      call compare_checksums( x(isd:ied+shift,jsd:jed+shift,:,n), global1(isd:ied+shift,jsd:jed+shift,:,n), &
                              trim(type2)//' BGRID_NE X')
      call compare_checksums( y(isd:ied+shift,jsd:jed+shift,:,n), global2(isd:ied+shift,jsd:jed+shift,:,n), &
                              trim(type2)//' BGRID_NE Y')
   end do

   deallocate(global1_all, global2_all, global1, global2, x, y)

    !------------------------------------------------------------------
    !  vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- setup data
    allocate(global1_all(nxm+shift,nym,nz, ntiles), global2_all(nxm,nym+shift,nz, ntiles) )  
    allocate(global1(1-whalo:ni+ehalo+shift, 1-shalo:nj+nhalo,       nz, ntile_per_pe) )  
    allocate(global2(1-whalo:ni+ehalo,       1-shalo:nj+nhalo+shift, nz, ntile_per_pe) )   
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, nym
             do i = 1, nxm+shift
                global1_all(i,j,k,n) = 1.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
          do j = 1, nym+shift
             do i = 1, nxm
                global2_all(i,j,k,n) = 2.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    !------------------------------------------------------------------------
    ! --- make sure consisency on the boundary for Five-Tile mosaic
    ! --- east boundary will take the value of neighbor tile west,
    ! --- north boundary will take the value of neighbor tile south.
    !------------------------------------------------------------------------
    if(type == 'Five-Tile') then
       global1_all(nxm+1,    1:ny,:,1) = global1_all(1,    1:ny,:,2)  ! east
       global1_all(nxm+1,ny+1:nym,:,1) = global1_all(1,    1:ny,:,4)  ! east
       global1_all(nx+1,     1:ny,:,2) = global1_all(1,    1:ny,:,3)  ! east
       global1_all(nx+1,     1:ny,:,3) = global1_all(1,    1:ny,:,1)  ! east
       global1_all(nx+1,     1:ny,:,4) = global1_all(1,    1:ny,:,5)  ! east
       global1_all(nx+1,     1:ny,:,5) = global1_all(1,ny+1:nym,:,1)  ! east
       global2_all(1:nxm,   nym+1,:,1) = global2_all(1:nxm,   1,:,1)  ! north
       global2_all(1:nx,     ny+1,:,2) = global2_all(1:nx,    1,:,4)  ! north
       global2_all(1:nx,     ny+1,:,3) = global2_all(1:nx,    1,:,5)  ! north
       global2_all(1:nx,     ny+1,:,4) = global2_all(1:nx,    1,:,2)  ! north
       global2_all(1:nx,     ny+1,:,5) = global2_all(1:nx,    1,:,3)  ! north
    end if

    do n = 1, ntile_per_pe
       global1(1:ni+shift,      1:nj,:,n) = global1_all(1:ni+shift,      1:nj,:,tile(n))
       global2(1:ni,      1:nj+shift,:,n) = global2_all(1:ni,      1:nj+shift,:,tile(n))
    end do

    allocate( x (ism:iem+shift,      jsm:jem,nz,ntile_per_pe) )
    allocate( y (ism:iem,      jsm:jem+shift,nz,ntile_per_pe) )

    x = 0.; y = 0
    x (isc:iec+shift,      jsc:jec,:,:) = global1(isc:iec+shift,      jsc:jec,:,:)
    y (isc:iec,      jsc:jec+shift,:,:) = global2(isc:iec,      jsc:jec+shift,:,:)

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points.     
    !-----------------------------------------------------------------------
    do n = 1, ntile_per_pe
       call fill_five_tile_halo(global1(:,:,:,n), global1_all, tile(n), shift, 0)
       call fill_five_tile_halo(global2(:,:,:,n), global2_all, tile(n), 0, shift)
    end do

    id = mpp_clock_id( type//' CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    do n = 1, ntile_per_pe
       call mpp_update_domains( x(:,:,:,n), y(:,:,:,n), domain, gridtype=CGRID_NE, tile_count = n )
    end do
    call mpp_clock_end(id)

   do n = 1, ntile_per_pe
      write(type2, *)type, " at tile_count = ",n
      call compare_checksums( x(isd:ied+shift,jsd:jed,:,n), global1(isd:ied+shift,jsd:jed,:,n), &
                              trim(type2)//' CGRID_NE X')
      call compare_checksums( y(isd:ied,jsd:jed+shift,:,n), global2(isd:ied,jsd:jed+shift,:,n), &
                              trim(type2)//' CGRID_NE Y')
   end do

   deallocate(global1_all, global2_all, global1, global2, x, y)

  end subroutine test_nonuniform_mosaic

  subroutine fill_five_tile_halo(data, data_all, tile, ioff, joff)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data_all
    integer,                              intent(in)    :: tile, ioff, joff
    integer                                             :: nxm, nym

    nxm = 2*nx; nym = 2*ny

    select case(tile)
    case(1)
       data(nxm+1+ioff:nxm+ehalo+ioff,                     1:ny,:) = data_all(1+ioff:ehalo+ioff,              1:ny,:,2) ! east
       data(nxm+1+ioff:nxm+ehalo+ioff,            ny+1:nym+joff,:) = data_all(1+ioff:ehalo+ioff,         1:ny+joff,:,4) ! east
       data(1-whalo:0,                                     1:ny,:) = data_all(nx-whalo+1:nx,                  1:ny,:,3) ! west
       data(1-whalo:0,                            ny+1:nym+joff,:) = data_all(nx-whalo+1:nx,             1:ny+joff,:,5) ! west
       data(1:nxm+ioff,                               1-shalo:0,:) = data_all(1:nxm+ioff,          nym-shalo+1:nym,:,1) ! south
       data(1:nxm+ioff,               nym+1+joff:nym+nhalo+joff,:) = data_all(1:nxm+ioff,        1+joff:nhalo+joff,:,1) ! north
       data(nxm+1+ioff:nxm+ehalo+ioff,                1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,     ny-shalo+1:ny,:,4) ! southeast
       data(1-whalo:0,                                1-shalo:0,:) = data_all(nx-whalo+1:nx,         ny-shalo+1:ny,:,5) ! southwest
       data(nxm+1+ioff:nxm+ehalo+ioff,nym+1+joff:nym+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff, 1+joff:nhalo+joff,:,2) ! northeast
       data(1-whalo:0,                nym+1+joff:nym+nhalo+joff,:) = data_all(nx-whalo+1:nx,     1+joff:nhalo+joff,:,3) ! northwest
    case(2)
       data(nx+1+ioff:nx+ehalo+ioff,              1:ny+joff,:) = data_all(1+ioff:ehalo+ioff,              1:ny+joff,:,3) ! east
       data(1-whalo:0,                            1:ny+joff,:) = data_all(nxm-whalo+1:nxm,                1:ny+joff,:,1) ! west
       data(1:nx+ioff,                            1-shalo:0,:) = data_all(1:nx+ioff,                  ny-shalo+1:ny,:,4) ! south 
       data(1:nx+ioff,              ny+1+joff:ny+nhalo+joff,:) = data_all(1:nx+ioff,              1+joff:nhalo+joff,:,4) ! north
       data(nx+1+ioff:nx+ehalo+ioff,              1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,          ny-shalo+1:ny,:,5) ! southeast
       data(1-whalo:0,                            1-shalo:0,:) = data_all(nxm-whalo+1:nxm,          nym-shalo+1:nym,:,1) ! southwest
       data(nx+1+ioff:nx+ehalo+ioff,ny+1+joff:ny+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff,      1+joff:nhalo+joff,:,5) ! northeast
       data(1-whalo:0,              ny+1+joff:ny+nhalo+joff,:) = data_all(nxm-whalo+1:nxm,  ny+1+joff:ny+nhalo+joff,:,1) ! northwest
    case(3)
       data(nx+1+ioff:nx+ehalo+ioff,              1:ny+joff,:) = data_all(1+ioff:ehalo+ioff,              1:ny+joff,:,1) ! east
       data(1-whalo:0,                            1:ny+joff,:) = data_all(nx-whalo+1:nx,                  1:ny+joff,:,2) ! west
       data(1:nx+ioff,                            1-shalo:0,:) = data_all(1:nx+ioff,                  ny-shalo+1:ny,:,5) ! south 
       data(1:nx+ioff,              ny+1+joff:ny+nhalo+joff,:) = data_all(1:nx+ioff,              1+joff:nhalo+joff,:,5) ! north
       data(nx+1+ioff:nx+ehalo+ioff,              1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,        nym-shalo+1:nym,:,1) ! southeast
       data(1-whalo:0,                            1-shalo:0,:) = data_all(nx-whalo+1:nx,              ny-shalo+1:ny,:,4) ! southwest
       data(nx+1+ioff:nx+ehalo+ioff,ny+1+joff:ny+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff,ny+1+joff:ny+nhalo+joff,:,1) ! northeast
       data(1-whalo:0,              ny+1+joff:ny+nhalo+joff,:) = data_all(nx-whalo+1:nx,          1+joff:nhalo+joff,:,4) ! northwest
    case(4)
       data(nx+1+ioff:nx+ehalo+ioff,              1:ny+joff,:) = data_all(1+ioff:ehalo+ioff,        1:ny+joff,:,5) ! east
       data(1-whalo:0,                            1:ny+joff,:) = data_all(nxm-whalo+1:nxm,     ny+1:2*ny+joff,:,1) ! west
       data(1:nx+ioff,                            1-shalo:0,:) = data_all(1:nx+ioff,            ny-shalo+1:ny,:,2) ! south 
       data(1:nx+ioff,              ny+1+joff:ny+nhalo+joff,:) = data_all(1:nx+ioff,        1+joff:nhalo+joff,:,2) ! north
       data(nx+1+ioff:nx+ehalo+ioff,              1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,    ny-shalo+1:ny,:,3) ! southeast
       data(1-whalo:0,                            1-shalo:0,:) = data_all(nxm-whalo+1:nxm,      ny-shalo+1:ny,:,1) ! southwest
       data(nx+1+ioff:nx+ehalo+ioff,ny+1+joff:ny+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff,1+joff:nhalo+joff,:,3) ! northeast
       data(1-whalo:0,              ny+1+joff:ny+nhalo+joff,:) = data_all(nxm-whalo+1:nxm,  1+joff:nhalo+joff,:,1) ! northwest
    case(5)
       data(nx+1+ioff:nx+ehalo+ioff,            1:  ny+joff,:) = data_all(1+ioff:ehalo+ioff,   ny+1:2*ny+joff,:,1) ! east
       data(1-whalo:0,                            1:ny+joff,:) = data_all(nx-whalo+1:nx,            1:ny+joff,:,4) ! west
       data(1:nx+ioff,                            1-shalo:0,:) = data_all(1:nx+ioff,            ny-shalo+1:ny,:,3) ! south 
       data(1:nx+ioff,              ny+1+joff:ny+nhalo+joff,:) = data_all(1:nx+ioff,        1+joff:nhalo+joff,:,3) ! north
       data(nx+1+ioff:nx+ehalo+ioff,              1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,    ny-shalo+1:ny,:,1) ! southeast
       data(1-whalo:0,                            1-shalo:0,:) = data_all(nx-whalo+1:nx,        ny-shalo+1:ny,:,2) ! southwest
       data(nx+1+ioff:nx+ehalo+ioff,ny+1+joff:ny+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff,1+joff:nhalo+joff,:,1) ! northeast
       data(1-whalo:0,              ny+1+joff:ny+nhalo+joff,:) = data_all(nx-whalo+1:nx,    1+joff:nhalo+joff,:,2) ! northwest
    end select

  end subroutine fill_five_tile_halo

  !##############################################################################
  !--- The following is used to test the refined mosaic. Three cases will be tested, 
  !--- non-symmetric regular mosaic, symmetric regular mosaic cubic grid. The regular mosaic 
  !--- contains 4 tiles. East of tile 1 connected with West of tile 2 (refine = 3)
  !--- and vice verse; East of tile 3 connected with West of tile 4 (refine = 2)
  !--- and vice verse; North of tile 1 connected with South of tile 3 (refine = 2)
  !--- and vice verse; North of tile 2 connected with South of tile 4 (refine = 1)
  !--- and vice verse; So there ar total 8 contacts. 
  subroutine test_refined_mosaic(type)
    character(len=*), intent(in) :: type

    type(domain2D)                           :: domain
    integer,              dimension(4)       :: isMe1, ieMe1, jsMe1, jeMe1
    integer,              dimension(4)       :: isNb1, ieNb1, jsNb1, jeNb1
    integer,              dimension(4)       :: isMe2, ieMe2, jsMe2, jeMe2
    integer,              dimension(4)       :: isNb2, ieNb2, jsNb2, jeNb2
    integer,              dimension(4)       :: rotation1, rotation2, dirMe1, dirMe2
    integer, allocatable, dimension(:,:)     :: from_tile1, from_tile2
    integer                                  :: ntiles, num_contact, npes_on_tile
    integer                                  :: totpoints, maxtotal, pos, ntiles_on_pe
    integer                                  :: tNb, nimax, njmax, avgpoints
    integer                                  :: n, m, l, te, tse, ts, tsw, tw, tnw, tn, tne, nn
    integer                                  :: noverlap1, noverlap2, total1, total2
    integer                                  :: isc, iec, jsc, jec, isg, ieg, jsg, jeg
    integer                                  :: ism, iem, jsm, jem, isd, ied, jsd, jed
    integer, allocatable, dimension(:)       :: tiles, ni, nj
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout, global_indices
    real,    allocatable, dimension(:,:,:,:) :: global_all, global1, global2
    real,    allocatable, dimension(:,:,:,:) :: x, y, x1, y1, x2, y2
    real,    allocatable, dimension(:,:,:,:) :: global1_all, global2_all
    real,    allocatable, dimension(:,:)     :: buffer, buffer1, buffer2, bufferx, buffery
    real,    allocatable, dimension(:,:)     :: bufferx1, buffery1, bufferx2, buffery2
    integer                                  :: shift
    character(len=128)                       :: type2
    logical                                  :: found

    !--- check the type
    select case(type)
    case ("Refined-Four-Tile", "Refined-Symmetric-Four-Tile" )
       ntiles = 4
       allocate(ni(ntiles), nj(ntiles))
       ! "Four-Tile" test case will only run on one pe or multiple 0f 8 ( balanced).
       if( npes .NE. 1 .AND. npes .NE. 8 .AND. npes .NE. 16 .AND. npes .NE. 32) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS(test_refined_mosaic: ' // &
                  type// ' mosaic will not be tested because npes is not 1, 8, 16 or 32')
          return
       end if            
       ni(1) =   nx; nj(1) =   ny
       ni(2) =   nx; nj(2) = 3*ny
       ni(3) = 2*nx; nj(3) =   ny
       ni(4) =   nx; nj(4) = 2*ny
       num_contact = 8   
    case ("Refined-Cubic-Grid")
       ntiles = 6; num_contact = 12
       allocate(ni(ntiles), nj(ntiles))
       ! "Cubic-Grid" will be tested only when nx = ny
       if( nx /= ny ) then
           call mpp_error(NOTE,'TEST_MPP_DOMAINS(test_refined_mosaic: ' // &
           type//' will not be tested because nx is not equal to ny' )
           return
       end if
       ! "Cubic-Grid" test case will only run on one pe or multiple 0f 16 pes ( balanced).
       if( npes .NE. 1 .AND. mod(npes,16) .NE. 0) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS(test_refined_mosaic: ' // &
                  type//' will not be tested because npes is not 1 and can not be divided by 16')
          return
       end if
       ni(1) =   nx; nj(1) =   ny
       ni(2) = 2*nx; nj(2) = 3*ny
       ni(3) = 2*nx; nj(3) =   ny
       ni(4) =   nx; nj(4) = 3*ny
       ni(5) = 2*nx; nj(5) =   ny
       ni(6) =   nx; nj(6) = 2*ny
    case default
       call mpp_error(FATAL, 'TEST_MPP_DOMAINS(test_refined_mosaic): no such test: '//type)
    end select

    allocate(layout(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    totpoints = sum(ni*nj)
    if(mod(totpoints, npes) .NE. 0) call mpp_error(FATAL,    &
        "TEST_MPP_DOMAINS(test_refined_mosaic): totpoints can not be divided by npes")
    avgpoints = totpoints/npes
    layout = 1
    pe_start = 0; pe_end = 0; pos = 0

    do n = 1, ntiles
       global_indices(:,n) = (/1, ni(n), 1, nj(n)/)
       if(npes > 1) then  ! no sharing processor between tiles
          if( mod(ni(n)*nj(n), avgpoints) .NE. 0) call mpp_error(FATAL, &
               'TEST_MPP_DOMAINS(test_refined_mosaic): number of points should be divided by average of points in each pe')
          npes_on_tile = ni(n)*nj(n)/avgpoints
          call mpp_define_layout( (/1,ni(n),1,nj(n)/), npes_on_tile, layout(:,n) )
          pe_start(n) = pos
          pe_end(n)   = pos + npes_on_tile - 1       
          pos         = pos + npes_on_tile
       end if
    end do

    ntiles_on_pe = 1
    if(npes == 1) then
       ntiles_on_pe = ntiles
       allocate(tiles(ntiles_on_pe))
       tiles = (/ (i, i=1,ntiles) /)
    else
       ntiles_on_pe = 1
       allocate(tiles(ntiles_on_pe))
       do n = 1, ntiles
          if( pe .GE. pe_start(n) .AND. pe .LE. pe_end(n) ) tiles = n
       end do
    end if

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) ) 
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) ) 

    !--- define domain
    select case(type)
    case( 'Refined-Four-Tile', 'Refined-Symmetric-Four-Tile' )
       call define_fourtile_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end, &
                                   type == 'Refined-Symmetric-Four-Tile'   )
    case( 'Refined-Cubic-Grid' )
       call define_cubic_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end )
    end select

    !--- first test mpp_get_mosaic_refine_overlap
    maxtotal = 0
    allocate(from_tile1(4, ntiles_on_pe), from_tile2(4, ntiles_on_pe))
    do n = 1, ntiles_on_pe
       rotation2 = ZERO
       noverlap1 = mpp_get_refine_overlap_number(domain, tile_count=n)
       call mpp_get_mosaic_refine_overlap(domain, isMe1, ieMe1, jsMe1, jeMe1, isNb1, ieNb1, jsNb1, jeNb1, &
                                          dirMe1, rotation1, tile_count = n)
       total1 = sum( (ieNb1(1:noverlap1)-isNb1(1:noverlap1)+1) * (jeNb1(1:noverlap1)-jsNb1(1:noverlap1)+1)  )

       !--- the following will figure out the overlapping
       call mpp_get_compute_domain(domain, isc, iec, jsc, jec, tile_count=n)
       call mpp_get_global_domain(domain, isg, ieg, jsg, jeg, tile_count=n)
       noverlap2 = 0; total2 = 0
       select case ( type )
       case ( 'Refined-Four-Tile', 'Refined-Symmetric-Four-Tile' )
          if( iec == ieg ) then   ! --- EAST
             noverlap2 = noverlap2 + 1
             if( mod(tiles(n),2) == 1) then ! tile 1, 3    
                tNb = tiles(n) + 1
             else                           ! tile 2, 4
                tNb = tiles(n) - 1
             endif
             from_tile2(noverlap2,n) = tNb
             dirMe2(noverlap2) = 1
             isMe2(noverlap2) = iec + 1;    ieMe2(noverlap2) = iec + ehalo
             jsMe2(noverlap2) = max(jsg,jsc - shalo); jeMe2(noverlap2) = min(jeg, jec+nhalo)
             isNb2(noverlap2) = 1;          ieNb2(noverlap2) = ehalo
             select case(tiles(n))
             case(1)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)*3+1;  jeNb2(noverlap2) = jeMe2(noverlap2)*3
             case(2)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)/3+1; jeNb2(noverlap2) = ceiling(jeMe2(noverlap2)/3.)
             case(3)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)*2+1;  jeNb2(noverlap2) = jeMe2(noverlap2)*2
             case(4)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)/2+1; jeNb2(noverlap2) = ceiling(jeMe2(noverlap2)/2.)
             end select
             total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
          end if

          if( jsc == jsg .AND. mod(tiles(n),2) == 1 ) then   ! --- SOUTH (only tile 1 and 3)
             noverlap2 = noverlap2 + 1
             tNb = mod(tiles(n)+2, ntiles)
             from_tile2(noverlap2,n) = tNb
             dirMe2(noverlap2) = 3
             isMe2(noverlap2) = max(isg,isc-whalo); ieMe2(noverlap2) = min(ieg, iec+ehalo)
             jsMe2(noverlap2) = jsc - shalo;        jeMe2(noverlap2) = jsc - 1
             if(tiles(n) == 1) then  ! refinement is 2
                isNb2(noverlap2) = (isMe2(noverlap2)-1)*2+1;    ieNb2(noverlap2) = ieMe2(noverlap2)*2
             else                   ! refinement is 2
                isNb2(noverlap2) = (isMe2(noverlap2)-1)/2+1;    ieNb2(noverlap2) = ceiling(ieMe2(noverlap2)/2.)
             end if
             jsNb2(noverlap2) = nj(tNb) - shalo + 1; jeNb2(noverlap2) = nj(tNb)
             total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
          end if

          if( isc == isg ) then   ! --- WEST
             noverlap2 = noverlap2 + 1
             if( mod(tiles(n),2) == 1) then ! tile 1, 3    
                tNb = tiles(n) + 1
             else                           ! tile 2, 4
                tNb = tiles(n) - 1
             endif
             from_tile2(noverlap2,n) = tNb
             dirMe2(noverlap2) = 5
             isMe2(noverlap2) = isc - whalo;          ieMe2(noverlap2) = isc - 1
             jsMe2(noverlap2) = max(jsg,jsc - shalo); jeMe2(noverlap2) = min(jeg, jec+nhalo)
             isNb2(noverlap2) = ni(tNb) - whalo + 1; ieNb2(noverlap2) = ni(tNb)
             select case(tiles(n))
             case(1)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)*3+1;  jeNb2(noverlap2) = jeMe2(noverlap2)*3
             case(2)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)/3+1; jeNb2(noverlap2) = ceiling(jeMe2(noverlap2)/3.)
             case(3)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)*2+1;  jeNb2(noverlap2) = jeMe2(noverlap2)*2
             case(4)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)/2+1; jeNb2(noverlap2) = ceiling(jeMe2(noverlap2)/2.)
             end select
             total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
          end if

          if( jec == jeg .AND. mod(tiles(n),2) == 1 ) then   ! --- NORTH (only tile 1 and 3)
             noverlap2 = noverlap2 + 1
             tNb = mod(tiles(n)+2, ntiles)
             from_tile2(noverlap2,n) = tNb
             dirMe2(noverlap2) = 7
             isMe2(noverlap2) = max(isg,isc-whalo); ieMe2(noverlap2) = min(ieg, iec+ehalo)
             jsMe2(noverlap2) = jec + 1;        jeMe2(noverlap2) = jec + nhalo
             if(tiles(n) == 1) then  ! refinement is 2
                isNb2(noverlap2) = (isMe2(noverlap2)-1)*2+1;    ieNb2(noverlap2) = ieMe2(noverlap2)*2
             else                   ! refinement is 2
                isNb2(noverlap2) = (isMe2(noverlap2)-1)/2+1;    ieNb2(noverlap2) = ceiling(ieMe2(noverlap2)/2.)
             end if
             jsNb2(noverlap2) = 1;          jeNb2(noverlap2) = nhalo
             total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
          end if

       case ( 'Refined-Cubic-Grid' )
          select case( tiles(n) )
          case ( 1 )  ! possible refined overlap will be at EAST, WEST
             if( iec == ieg ) then   ! --- EAST
                noverlap2 = noverlap2 + 1
                tNb = 2
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 1
                isMe2(noverlap2) = iec + 1;                  ieMe2(noverlap2) = iec + ehalo
                jsMe2(noverlap2) = max(jsg,jsc - shalo);     jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = 1;                        ieNb2(noverlap2) = ehalo
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)*3+1; jeNb2(noverlap2) = jeMe2(noverlap2)*3
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1) 
             end if
             if( isc == isg ) then   ! --- WEST
                noverlap2 = noverlap2 + 1
                tNb = 5
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 5
                isMe2(noverlap2) = isc - whalo;                  ieMe2(noverlap2) = isc - 1
                jsMe2(noverlap2) = max(jsg,jsc - shalo);         jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = (nj(1)-jeMe2(noverlap2))*2+1; ieNb2(noverlap2) = (nj(1)-jsMe2(noverlap2)+1)*2
                jsNb2(noverlap2) = nj(tNb) - whalo + 1;          jeNb2(noverlap2) = nj(tNb)
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
                rotation2(noverlap2) = NINETY
             end if
          case ( 2 )  ! possible refined overlap will be at EAST, WEST
             if( iec == ieg ) then   ! --- EAST
                noverlap2 = noverlap2 + 1
                tNb = 4
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 1
                isMe2(noverlap2) = iec + 1;                      ieMe2(noverlap2) = iec + ehalo
                jsMe2(noverlap2) = max(jsg,jsc - shalo);         jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = (nj(2)-jeMe2(noverlap2))/3+1; ieNb2(noverlap2) = ceiling((nj(2)-jsMe2(noverlap2)+1)/3.)
                jsNb2(noverlap2) = 1;                            jeNb2(noverlap2) = ehalo
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
                rotation2(noverlap2) = NINETY       
             end if
             if( isc == isg ) then   ! --- WEST
                noverlap2 = noverlap2 + 1
                tNb = 1
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 5
                isMe2(noverlap2) = isc - whalo;              ieMe2(noverlap2) = isc - 1
                jsMe2(noverlap2) = max(jsg,jsc - shalo);     jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = ni(tNb) - whalo + 1;      ieNb2(noverlap2) = ni(tNb)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)/3+1; jeNb2(noverlap2) = ceiling(jeMe2(noverlap2)/3.)
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
             end if
          case ( 3 )  ! possible refined overlap will be at EAST, NORTH          
             if( iec == ieg ) then   ! --- EAST
                noverlap2 = noverlap2 + 1
                tNb = 4
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 1
                isMe2(noverlap2) = iec + 1;                  ieMe2(noverlap2) = iec + ehalo
                jsMe2(noverlap2) = max(jsg,jsc - shalo);     jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = 1;                        ieNb2(noverlap2) = ehalo
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)*3+1; jeNb2(noverlap2) = jeMe2(noverlap2)*3
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1) 
             end if
             if( jec == jeg ) then   ! --- NORTH
                noverlap2 = noverlap2 + 1
                tNb = 5
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 7
                isMe2(noverlap2) = max(isg,isc - whalo);        ieMe2(noverlap2) = min(ieg, iec+ehalo)
                jsMe2(noverlap2) = jec + 1;                     jeMe2(noverlap2) = jec + nhalo
                isNb2(noverlap2) = 1;                           ieNb2(noverlap2) = nhalo
                jsNb2(noverlap2) = (ni(3)-ieMe2(noverlap2))/2+1; jeNb2(noverlap2) = ceiling((ni(3)-isMe2(noverlap2)+1)/2.)
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1) 
                rotation2(noverlap2) = MINUS_NINETY
             end if
          case ( 4 )  ! possible refined overlap will be at NORTH, EAST, SOUTH, WEST
             if( jec == jeg ) then   ! --- NORTH
                noverlap2 = noverlap2 + 1
                tNb = 5
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 7
                isMe2(noverlap2) = max(isg,isc - whalo);        ieMe2(noverlap2) = min(ieg, iec+ehalo)
                jsMe2(noverlap2) = jec + 1;                     jeMe2(noverlap2) = jec + nhalo
                isNb2(noverlap2) = (isMe2(noverlap2)-1)*2+1;    ieNb2(noverlap2) = ieMe2(noverlap2)*2
                jsNb2(noverlap2) = 1;                           jeNb2(noverlap2) = nhalo
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1) 
             end if
             if( iec == ieg ) then   ! --- EAST
                noverlap2 = noverlap2 + 1
                tNb = 6
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 1
                isMe2(noverlap2) = iec + 1;                        ieMe2(noverlap2) = iec + ehalo
                jsMe2(noverlap2) = max(jsg,jsc - shalo);           jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = (nj(4)-jeMe2(noverlap2)-1)/3+1; ieNb2(noverlap2) = ceiling((nj(4)-jsMe2(noverlap2)+1)/3.)
                jsNb2(noverlap2) = 1;                              jeNb2(noverlap2) = ehalo
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
                rotation2(noverlap2) = NINETY
             end if
             if( jsc == jsg ) then   ! --- SOUTH
                noverlap2 = noverlap2 + 1
                tNb = 2
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 3
                isMe2(noverlap2) = max(isg,isc - whalo);        ieMe2(noverlap2) = min(ieg, iec+ehalo)
                jsMe2(noverlap2) = jsc - shalo;                 jeMe2(noverlap2) = jsc  - 1
                isNb2(noverlap2) = ni(tNb) - shalo + 1;         ieNb2(noverlap2) = ni(tNb)
                jsNb2(noverlap2) = (ni(4)-ieMe2(noverlap2))*3+1; jeNb2(noverlap2) = (ni(4)-isMe2(noverlap2)+1)*3
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1) 
                rotation2(noverlap2) = MINUS_NINETY
             end if
             if( isc == isg ) then   ! --- WEST
                noverlap2 = noverlap2 + 1
                tNb = 3
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 5
                isMe2(noverlap2) = isc - whalo;              ieMe2(noverlap2) = isc - 1
                jsMe2(noverlap2) = max(jsg,jsc - shalo);     jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = ni(tNb) - whalo + 1;      ieNb2(noverlap2) = ni(tNb)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)/3+1; jeNb2(noverlap2) = ceiling(jeMe2(noverlap2)/3.)
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
             end if

         case ( 5 )  ! possible refined overlap will be at EAST, NORTH, WEST, SOUTH       
             if( iec == ieg ) then   ! --- EAST
                noverlap2 = noverlap2 + 1
                tNb = 6
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 1
                isMe2(noverlap2) = iec + 1;                  ieMe2(noverlap2) = iec + ehalo
                jsMe2(noverlap2) = max(jsg,jsc - shalo);     jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = 1;                        ieNb2(noverlap2) = ehalo
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)*2+1; jeNb2(noverlap2) = jeMe2(noverlap2)*2
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1) 
             end if
             if( jec == jeg ) then   ! --- NORTH
                noverlap2 = noverlap2 + 1
                tNb = 1
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 7
                isMe2(noverlap2) = max(isg,isc - whalo);         ieMe2(noverlap2) = min(ieg, iec+ehalo)
                jsMe2(noverlap2) = jec + 1;                      jeMe2(noverlap2) = jec + nhalo
                isNb2(noverlap2) = 1;                            ieNb2(noverlap2) = nhalo
                jsNb2(noverlap2) = (ni(5)-ieMe2(noverlap2))/2+1; jeNb2(noverlap2) = ceiling((ni(5)-isMe2(noverlap2)+1)/2.)
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
                rotation2(noverlap2) = MINUS_NINETY 
             end if
             if( isc == isg ) then   ! --- WEST
                noverlap2 = noverlap2 + 1
                tNb = 3
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 5
                isMe2(noverlap2) = isc - whalo;                  ieMe2(noverlap2) = isc - 1
                jsMe2(noverlap2) = max(jsg,jsc - shalo);         jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = (nj(5)-jeMe2(noverlap2))*2+1; ieNb2(noverlap2) = (nj(5)-jsMe2(noverlap2)+1)*2
                jsNb2(noverlap2) = nj(tNb) - whalo + 1;          jeNb2(noverlap2) = nj(tNb)
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
                rotation2(noverlap2) = NINETY
             end if
             if( jsc == jsg ) then   ! --- SOUTH
                noverlap2 = noverlap2 + 1
                tNb = 4
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 3
                isMe2(noverlap2) = max(isg,isc - whalo);        ieMe2(noverlap2) = min(ieg, iec+ehalo)
                jsMe2(noverlap2) = jsc - shalo;                 jeMe2(noverlap2) = jsc  - 1
                isNb2(noverlap2) = (isMe2(noverlap2)-1)/2+1;    ieNb2(noverlap2) = ceiling(ieMe2(noverlap2)/2.)
                jsNb2(noverlap2) = nj(tNb) - shalo + 1;         jeNb2(noverlap2) = nj(tNb)
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1) 
             end if


          case ( 6 )  ! possible refined overlap will be at SOUTH, WEST
             if( jsc == jsg ) then   ! --- SOUTH
                noverlap2 = noverlap2 + 1
                tNb = 4
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 3
                isMe2(noverlap2) = max(isg,isc - whalo);        ieMe2(noverlap2) = min(ieg, iec+ehalo)
                jsMe2(noverlap2) = jsc - shalo;                 jeMe2(noverlap2) = jsc  - 1
                isNb2(noverlap2) = ni(tNb) - shalo + 1;         ieNb2(noverlap2) = ni(tNb)
                jsNb2(noverlap2) = (ni(6)-ieMe2(noverlap2))*3+1; jeNb2(noverlap2) = (ni(6)-isMe2(noverlap2)+1)*3
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1) 
                rotation2(noverlap2) = MINUS_NINETY
             end if
             if( isc == isg ) then   ! --- WEST
                noverlap2 = noverlap2 + 1
                tNb = 5
                from_tile2(noverlap2,n) = tNb
                dirMe2(noverlap2) = 5
                isMe2(noverlap2) = isc - whalo;              ieMe2(noverlap2) = isc - 1
                jsMe2(noverlap2) = max(jsg,jsc - shalo);     jeMe2(noverlap2) = min(jeg, jec+nhalo)
                isNb2(noverlap2) = ni(tNb) - whalo + 1;      ieNb2(noverlap2) = ni(tNb)
                jsNb2(noverlap2) = (jsMe2(noverlap2)-1)/2+1; jeNb2(noverlap2) = ceiling(jeMe2(noverlap2)/2.)
                total2 = total2 + (ieNb2(noverlap2) - isNb2(noverlap2) + 1) * (jeNb2(noverlap2) - jsNb2(noverlap2) + 1)
             end if
          end select
       end select

       if(total1 .NE. total2) call mpp_error(FATAL, "test_mpp_domains: mismatch on total number of points")
       !--- we add one extra point in each direction for the consideration of symmetric domain.
       total2 = sum( (ieNb2(1:noverlap2) - isNb2(1:noverlap2) + 2) * (jeNb2(1:noverlap2) - jsNb2(1:noverlap2) + 2) )
       maxtotal = max(maxtotal, total2)
       !--- comparing
       if( noverlap1 .NE. noverlap2 ) call mpp_error(FATAL, "test_mpp_domains: mismatch on number of overlapping region")
       do m = 1, noverlap1
          found = .false.
          do l = 1, noverlap2
             if(dirMe1(m) == dirMe2(l)) then
                found = .true.
                exit
             endif
          enddo
          from_tile1(m,n) = from_tile2(l,n)
          if(.not. found) call mpp_error(FATAL, "test_mpp_domains: mismatch on direction")
          if( (isMe1(m) .NE. isMe2(l)) .OR. (ieMe1(m) .NE. ieMe2(l))        &
              .OR. (jsMe1(m) .NE. jsMe2(l)) .OR. (jeMe1(m) .NE. jeMe2(l)) ) &
              call mpp_error(FATAL, "test_mpp_domains: mismatch on myself overlapping index")
          if( (isNb1(m) .NE. isNb2(l)) .OR. (ieNb1(m) .NE. ieNb2(l))        &
              .OR. (jsNb1(m) .NE. jsNb2(l)) .OR. (jeNb1(m) .NE. jeNb2(l)) ) &
              call mpp_error(FATAL, "test_mpp_domains: mismatch on neighbor overlapping index") 
          if(rotation1(m) .NE. rotation2(l)) call mpp_error(FATAL, "test_mpp_domains: mismatch on rotation angle");
       end do
    end do

    !--- setup data
    nimax = maxval(ni); njmax = maxval(nj)
    allocate(global_all(1:nimax,1:njmax,nz,ntiles) )
    allocate(global2(1-whalo:nimax+ehalo,1-shalo:njmax+nhalo,nz, ntiles_on_pe) ) 
    global2 = 0
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, njmax
             do i = 1, nimax
                global_all(i,j,k,n) = n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntiles_on_pe
       nn = tiles(n)
       global2(1:ni(nn),1:nj(nn),:,n) = global_all(1:ni(nn),1:nj(nn),:,nn)
    end do

    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x (ism:iem,jsm:jem,nz, ntiles_on_pe) )
    allocate( x1(ism:iem,jsm:jem,nz, ntiles_on_pe) )
    allocate( x2(ism:iem,jsm:jem,nz, ntiles_on_pe) )
    x = 0

    do n = 1, ntiles_on_pe
       call mpp_get_compute_domain( domain, isc, iec, jsc, jec, tile_count=n  )
       call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, tile_count=n )
       x(isc:iec,jsc:jec,:,:) = global2(isc:iec,jsc:jec,:,:)
    end do
    x1 = x; x2 = x

    allocate(buffer (maxtotal*nz, ntiles_on_pe))
    allocate(buffer1(maxtotal*nz, ntiles_on_pe))
    allocate(buffer2(maxtotal*nz, ntiles_on_pe))

    !--- call mpp_update_domains to update domain 
    do n = 1, ntiles_on_pe
       call mpp_update_domains( x(:,:,:,n), domain, buffer=buffer(:,n), tile_count=n )
    end do

    !--- multiple varaibles update
    do n = 1, ntiles_on_pe
       call mpp_update_domains( x1(:,:,:,n), domain, buffer=buffer1(:,n), complete=.false., tile_count=n )
       call mpp_update_domains( x2(:,:,:,n), domain, buffer=buffer2(:,n), complete=.true., tile_count=n )
    end do

    !--- fill up the value at halo points and compare the value at buffer.
    do n = 1, ntiles_on_pe  

       !--- comparing the buffer.
       noverlap1 = mpp_get_refine_overlap_number(domain, tile_count=n)
       call mpp_get_mosaic_refine_overlap(domain, isMe1, ieMe1, jsMe1, jeMe1, isNb1, ieNb1, jsNb1, jeNb1, &
                                          dirMe1, rotation1, tile_count = n)             
       pos = 0
       do m = 1, noverlap1       
          do k = 1, nz
             do j = jsNb1(m), jeNb1(m)
                do i = isNb1(m), ieNb1(m)
                   pos = pos + 1
                   if(global_all(i,j,k,from_tile1(m,n)) .NE. buffer(pos,n) ) then
                      write(stdunit, 111) 'x', type, mpp_pe(), i, j, k, buffer(pos,n), global_all(i,j,k,from_tile1(m,n)) 
                      call mpp_error(FATAL, "test_refined_mosaic: mismatch between buffer data and actual data for "//type )
                   end if
                   if(global_all(i,j,k,from_tile1(m,n)) .NE. buffer1(pos,n) ) then
                      write(stdunit, 111) 'x1', type, mpp_pe(), i, j, k, buffer1(pos,n), global_all(i,j,k,from_tile1(m,n)) 
                      call mpp_error(FATAL, "test_refined_mosaic: mismatch between buffer data and actual data for "//type )
                   end if
                   if(global_all(i,j,k,from_tile1(m,n)) .NE. buffer2(pos,n) ) then
                      write(stdunit, 111) 'x2', type, mpp_pe(), i, j, k, buffer2(pos,n), global_all(i,j,k,from_tile1(m,n)) 
                      call mpp_error(FATAL, "test_refined_mosaic: mismatch between buffer data and actual data for "//type )
                   end if
                end do
             end do
          end do
       end do

       !--- fill the halo and compare
       select case(type)
       case('Refined-Four-Tile', 'Refined-Symmetric-Four-Tile')
          te = 0; ts = 0; tn = 0; tw = 0
          select case(tiles(n))       
          case(1)
             tsw = 4 
          case(2)
             tsw = 3; tn = 4; ts = 4
          case(3)
             tsw = 2
          case(4)
             tsw = 1; ts = 2; tn = 2
          end select
          tse = tsw; tnw = tsw; tne = tsw
          call fill_regular_refinement_halo( global2(:,:,:,n), global_all, ni, nj, tiles(n), &
                                             te, tse, ts, tsw, tw, tnw, tn, tne, 0, 0)
       case('Refined-Cubic-Grid')       
          call fill_cubicgrid_refined_halo(global2(:,:,:,n), global_all, global_all, ni, nj, tiles(n), 0, 0, 1, 1 )
       end select
       call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
       type2 = type
       if(ntiles_on_pe>1) write(type2, *)trim(type2), " at tile_count = ", tiles(n)
       call compare_checksums( x (isd:ied,jsd:jed,:,n), global2(isd:ied,jsd:jed,:,n), trim(type2)//' X' )
       call compare_checksums( x1(isd:ied,jsd:jed,:,n), global2(isd:ied,jsd:jed,:,n), trim(type2)//' X1' )
       call compare_checksums( x2(isd:ied,jsd:jed,:,n), global2(isd:ied,jsd:jed,:,n), trim(type2)//' X2' )
    end do  

    deallocate(global2, global_all, x, x1, x2)
    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    !--- setup data
    shift = 0
    if( type == 'Refined-Symmetric-Four-Tile' .OR. type == 'Refined-Cubic-Grid' ) shift = 1

    nimax = maxval(ni) + shift; njmax = maxval(nj) + shift
    allocate(global1_all(1:nimax,1:njmax,nz,ntiles) )
    allocate(global2_all(1:nimax,1:njmax,nz,ntiles) )
    allocate(global1(1-whalo:nimax+ehalo,1-shalo:njmax+nhalo,nz, ntiles_on_pe) ) 
    allocate(global2(1-whalo:nimax+ehalo,1-shalo:njmax+nhalo,nz, ntiles_on_pe) ) 
    global1 = 0; global2 = 0
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, njmax
             do i = 1, nimax
                global1_all(i,j,k,n) = 1.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,n) = 2.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    !--- make sure consistency on the no-refinement boundary  for symmetric domain.
    !--- For "Symmetric four tile" mosaic, north of tile 2 and 4 need to be filled.
    !--- For "Cubic-Grid", The following need to be filled: north of tile 1,
    !--- north of tile 2, east and north of tile 6. The velocity at the corner point will be 0.
    select case( type ) 
    case ( 'Refined-Symmetric-Four-Tile' ) 
       global1_all(1:ni(2)+1,nj(2)+1,:,2) = global1_all(1:ni(2)+1,1,:,4)
       global2_all(1:ni(2)+1,nj(2)+1,:,2) = global2_all(1:ni(2)+1,1,:,4)
       global1_all(1:ni(4)+1,nj(4)+1,:,4) = global1_all(1:ni(4)+1,1,:,2)
       global2_all(1:ni(4)+1,nj(4)+1,:,4) = global2_all(1:ni(4)+1,1,:,2)
    case ( 'Refined-Cubic-Grid' )
       global1_all(1:ni(1)+1,nj(1)+1,:,1) = -global2_all(1,nj(3)+1:1:-1,:,3) ! north
       global2_all(1:ni(1)+1,nj(1)+1,:,1) =  global1_all(1,nj(3)+1:1:-1,:,3) ! north
       global1_all(1:ni(2)+1,nj(2)+1,:,2) =  global1_all(1:ni(3)+1,1,   :,3) ! north
       global2_all(1:ni(2)+1,nj(2)+1,:,2) =  global2_all(1:ni(3)+1,1,   :,3) ! north
       global1_all(1:ni(6)+1,nj(6)+1,:,6) =  global1_all(1:ni(1)+1,1,   :,1) ! north
       global2_all(1:ni(6)+1,nj(6)+1,:,6) =  global2_all(1:ni(1)+1,1,   :,1) ! north
       global1_all(ni(6)+1,1:nj(6)+1,:,6) =  global2_all(ni(2)+1:1:-1,1,:,2)  ! east 
       global2_all(ni(6)+1,1:nj(6)+1,:,6) = -global1_all(ni(2)+1:1:-1,1,:,2)  ! east 
       do n = 1, ntiles
          global1_all(1,      1,:,n) = 0;  global1_all(1,      nj(n)+1,:,n) = 0;
          global1_all(ni(n)+1,1,:,n) = 0;  global1_all(ni(n)+1,nj(n)+1,:,n) = 0;
          global2_all(1,      1,:,n) = 0;  global2_all(1,      nj(n)+1,:,n) = 0;
          global2_all(ni(n)+1,1,:,n) = 0;  global2_all(ni(n)+1,nj(n)+1,:,n) = 0;
       end do
    end select

    do n = 1, ntiles_on_pe
       nn = tiles(n)
       global1(1:ni(nn)+shift,1:nj(nn)+shift,:,n) = global1_all(1:ni(nn)+shift,1:nj(nn)+shift,:,nn)
       global2(1:ni(nn)+shift,1:nj(nn)+shift,:,n) = global2_all(1:ni(nn)+shift,1:nj(nn)+shift,:,nn)
    end do

    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x  (ism:iem+shift,jsm:jem+shift,nz, ntiles_on_pe) )
    allocate( y  (ism:iem+shift,jsm:jem+shift,nz, ntiles_on_pe) )
    allocate( x1 (ism:iem+shift,jsm:jem+shift,nz, ntiles_on_pe) )
    allocate( y1 (ism:iem+shift,jsm:jem+shift,nz, ntiles_on_pe) )
    allocate( x2 (ism:iem+shift,jsm:jem+shift,nz, ntiles_on_pe) )
    allocate( y2 (ism:iem+shift,jsm:jem+shift,nz, ntiles_on_pe) )
    x = 0; y = 0

    do n = 1, ntiles_on_pe
       call mpp_get_compute_domain( domain, isc, iec, jsc, jec, tile_count=n )
       call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, tile_count=n )
       x(isc:iec+shift,jsc:jec+shift,:,:) = global1(isc:iec+shift,jsc:jec+shift,:,:)
       y(isc:iec+shift,jsc:jec+shift,:,:) = global2(isc:iec+shift,jsc:jec+shift,:,:)
    end do
    x1 = x; x2 =x; y1 = y; y2 = y

    allocate(bufferx(maxtotal*nz, ntiles_on_pe),  buffery(maxtotal*nz, ntiles_on_pe) )
    allocate(bufferx1(maxtotal*nz, ntiles_on_pe), buffery1(maxtotal*nz, ntiles_on_pe) )
    allocate(bufferx2(maxtotal*nz, ntiles_on_pe), buffery2(maxtotal*nz, ntiles_on_pe) )

    !--- call mpp_update_domains to update domain 
    do n = 1, ntiles_on_pe
       call mpp_update_domains( x(:,:,:,n), y(:,:,:,n), domain, gridtype=BGRID_NE, &
                                bufferx=bufferx(:,n), buffery=buffery(:,n), tile_count=n )
    end do

    !--- multiple update
    do n = 1, ntiles_on_pe
       call mpp_update_domains( x1(:,:,:,n), y1(:,:,:,n), domain, gridtype=BGRID_NE, &
                                bufferx=bufferx1(:,n), buffery=buffery1(:,n), complete=.false., tile_count=n )
       call mpp_update_domains( x2(:,:,:,n), y2(:,:,:,n), domain, gridtype=BGRID_NE, &
                                bufferx=bufferx2(:,n), buffery=buffery2(:,n), complete=.true., tile_count=n )
    end do

    !--- fill up the value at halo points and compare the value at buffer.
    do n = 1, ntiles_on_pe  
       !--- comparing the buffer.
       noverlap1 = mpp_get_refine_overlap_number(domain, tile_count=n, position = CORNER)
       call mpp_get_mosaic_refine_overlap(domain, isMe1, ieMe1, jsMe1, jeMe1, isNb1, ieNb1, jsNb1, jeNb1, &
               dirMe1, rotation1, tile_count = n, position = CORNER)             
       pos = 0

       do m = 1, noverlap1  
          select case( rotation1(m) )
          case (ZERO)     
             do k = 1, nz
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. bufferx(pos,n) ) then
                         write(stdunit,111)'x','BGRID '//type,mpp_pe(),i,j,k,bufferx(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. buffery(pos,n) ) then
                         write(stdunit,111)'y','BGRID '//type,mpp_pe(),i,j,k,buffery(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. bufferx1(pos,n) ) then
                         write(stdunit,111)'x1','BGRID '//type,mpp_pe(),i,j,k,bufferx1(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X1")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. buffery1(pos,n) ) then
                         write(stdunit,111)'y1','BGRID '//type,mpp_pe(),i,j,k,buffery1(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y1")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. bufferx2(pos,n) ) then
                         write(stdunit,111)'x2','BGRID '//type,mpp_pe(),i,j,k,bufferx2(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X2")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. buffery2(pos,n) ) then
                         write(stdunit,111)'y2','BGRID '//type,mpp_pe(),i,j,k,buffery2(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y2")
                      end if
                   end do
                end do
             end do
          case (NINETY) ! S->E, N->W, u->-v, v->u    
             do k = 1, nz
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx(pos,n) ) then
                         write(stdunit,111)'x','BGRID '//type,mpp_pe(),i,j,k,bufferx(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X")
                      end if
                      if(-global1_all(i,j,k,from_tile1(m,n)) .NE. buffery(pos,n) ) then
                         write(stdunit,111)'y','BGRID '//type,mpp_pe(),i,j,k,buffery(pos,n),-global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx1(pos,n) ) then
                         write(stdunit,111)'x1','BGRID '//type,mpp_pe(),i,j,k,bufferx1(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X1")
                      end if
                      if(-global1_all(i,j,k,from_tile1(m,n)) .NE. buffery1(pos,n) ) then
                         write(stdunit,111)'y1','BGRID '//type,mpp_pe(),i,j,k,buffery1(pos,n),-global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y1")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx2(pos,n) ) then
                         write(stdunit,111)'x2','BGRID '//type,mpp_pe(),i,j,k,bufferx2(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X2")
                      end if
                      if(-global1_all(i,j,k,from_tile1(m,n)) .NE. buffery2(pos,n) ) then
                         write(stdunit,111)'y2','BGRID '//type,mpp_pe(),i,j,k,buffery2(pos,n),-global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y2")
                      end if
                   end do
                end do
             end do
          case (MINUS_NINETY) ! S->E, N->W, u->-v, v->u    
             do k = 1, nz
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(-global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx(pos,n) ) then
                         write(stdunit,111)'x','BGRID '//type,mpp_pe(),i,j,k,bufferx(pos,n),-global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. buffery(pos,n) ) then
                         write(stdunit,111)'y','BGRID '//type,mpp_pe(),i,j,k,buffery(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y")
                      end if
                      if(-global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx1(pos,n) ) then
                         write(stdunit,111)'x1','BGRID '//type,mpp_pe(),i,j,k,bufferx1(pos,n),-global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X1")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. buffery1(pos,n) ) then
                         write(stdunit,111)'y1','BGRID '//type,mpp_pe(),i,j,k,buffery1(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y1")
                      end if
                      if(-global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx2(pos,n) ) then
                         write(stdunit,111)'x2','BGRID '//type,mpp_pe(),i,j,k,bufferx2(pos,n),-global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" X2")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. buffery2(pos,n) ) then
                         write(stdunit,111)'y2','BGRID '//type,mpp_pe(),i,j,k,buffery2(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: BGRID "//type//" Y2")
                      end if
                   end do
                end do
             end do
          end select
       end do

       !--- fill the halo data and compare
       select case(type)
       case('Refined-Four-Tile', 'Refined-Symmetric-Four-Tile')
          te = 0; ts = 0; tn = 0; tw = 0
          select case(tiles(n))       
          case(1)
             tsw = 4 
          case(2)
             tsw = 3; tn = 4; ts = 4
          case(3)
             tsw = 2
          case(4)
             tsw = 1; ts = 2; tn = 2
          end select
          tse = tsw; tnw = tsw; tne = tsw
          call fill_regular_refinement_halo( global1(:,:,:,n), global1_all, ni, nj, tiles(n), &
               te, tse, ts, tsw, tw, tnw, tn, tne, shift, shift )
          call fill_regular_refinement_halo( global2(:,:,:,n), global2_all, ni, nj, tiles(n), &
               te, tse, ts, tsw, tw, tnw, tn, tne, shift, shift )
       case('Refined-Cubic-Grid')       
          call fill_cubicgrid_refined_halo(global1(:,:,:,n), global1_all, global2_all, ni, nj, tiles(n), 1, 1, 1, -1 )
          call fill_cubicgrid_refined_halo(global2(:,:,:,n), global2_all, global1_all, ni, nj, tiles(n), 1, 1, -1, 1 )
       end select

       call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
       write(type2, *)"BGRID ", type
       if(ntiles_on_pe>1) write(type2, *)trim(type2), " at tile_count = ", tiles(n)
       call compare_checksums( x (isd:ied+shift,jsd:jed+shift,:,n), global1(isd:ied+shift,jsd:jed+shift,:,n), trim(type2)//' X' )
       call compare_checksums( x1(isd:ied+shift,jsd:jed+shift,:,n), global1(isd:ied+shift,jsd:jed+shift,:,n), trim(type2)//' X1')
       call compare_checksums( x2(isd:ied+shift,jsd:jed+shift,:,n), global1(isd:ied+shift,jsd:jed+shift,:,n), trim(type2)//' X2')
       write(type2, *)"BGRID ", type
       if(ntiles_on_pe>1) write(type2, *)trim(type2), " at tile_count = ", tiles(n)
       call compare_checksums( y (isd:ied+shift,jsd:jed+shift,:,n), global2(isd:ied+shift,jsd:jed+shift,:,n), trim(type2)//' Y' )
       call compare_checksums( y1(isd:ied+shift,jsd:jed+shift,:,n), global2(isd:ied+shift,jsd:jed+shift,:,n), trim(type2)//' Y1')
       call compare_checksums( y2(isd:ied+shift,jsd:jed+shift,:,n), global2(isd:ied+shift,jsd:jed+shift,:,n), trim(type2)//' Y2')
    end do  

    deallocate(global1_all, global2_all, global1, global2, x, y, x1, x2, y1, y2 )
    !------------------------------------------------------------------
    !              vector update : CGRID_NE, one extra point may needed to symmetric domain
    !------------------------------------------------------------------
    !--- setup data

    nimax = maxval(ni); njmax = maxval(nj)
    allocate(global1_all(1:nimax+shift,1:njmax,nz,ntiles) )
    allocate(global2_all(1:nimax,1:njmax+shift,nz,ntiles) )
    allocate(global1(1-whalo:nimax+ehalo+shift,1-shalo:njmax+nhalo,nz, ntiles_on_pe) ) 
    allocate(global2(1-whalo:nimax+ehalo,1-shalo:njmax+nhalo+shift,nz, ntiles_on_pe) ) 
    global1 = 0; global2 = 0
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, njmax
             do i = 1, nimax + shift
                global1_all(i,j,k,n) = 1.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
          do j = 1, njmax + shift
             do i = 1, nimax
                global2_all(i,j,k,n) = 2.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    !--- make sure consistency on the no-refinement boundary  for symmetric domain.
    !--- For "Symmetric four tile" mosaic, north of tile 2 and 4 need to be filled.
    !--- For "Cubic-Grid", The following need to be filled: north of tile 1,
    !--- north of tile 2, east and north of tile 6. The velocity at the corner point will be 0.
    select case( type ) 
    case ('Refined-Symmetric-Four-Tile' ) 
       global2_all(1:ni(2),nj(2)+1,:,2) = global2_all(1:ni(2),1,:,4)
       global2_all(1:ni(4),nj(4)+1,:,4) = global2_all(1:ni(4),1,:,2)
    case ('Refined-Cubic-Grid' )
       global2_all(1:ni(1),nj(1)+1,:,1) =  global1_all(1,nj(3):1:-1,:,3) ! north
       global2_all(1:ni(2),nj(2)+1,:,2) =  global2_all(1:ni(3),1,   :,3) ! north
       global2_all(1:ni(6),nj(6)+1,:,6) =  global2_all(1:ni(1),1,   :,1) ! north
       global1_all(ni(6)+1,1:nj(6),:,6) =  global2_all(ni(2):1:-1,1,:,2)  ! east 
    end select

    do n = 1, ntiles_on_pe
       nn = tiles(n)
       global1(1:ni(nn)+shift,1:nj(nn),:,n) = global1_all(1:ni(nn)+shift,1:nj(nn),:,nn)
       global2(1:ni(nn),1:nj(nn)+shift,:,n) = global2_all(1:ni(nn),1:nj(nn)+shift,:,nn)
    end do

    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x  (ism:iem+shift,jsm:jem,nz, ntiles_on_pe) )
    allocate( y  (ism:iem,jsm:jem+shift,nz, ntiles_on_pe) )
    allocate( x1 (ism:iem+shift,jsm:jem,nz, ntiles_on_pe) )
    allocate( y1 (ism:iem,jsm:jem+shift,nz, ntiles_on_pe) )
    allocate( x2 (ism:iem+shift,jsm:jem,nz, ntiles_on_pe) )
    allocate( y2 (ism:iem,jsm:jem+shift,nz, ntiles_on_pe) )
    x = 0; y = 0
    bufferx  = 0; buffery = 0 
    bufferx1 = 0; buffery1 = 0
    bufferx2 = 0; buffery2 = 0
    do n = 1, ntiles_on_pe
       call mpp_get_compute_domain( domain, isc, iec, jsc, jec, tile_count=n )
       call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, tile_count=n )
       x(isc:iec+shift,jsc:jec,:,:) = global1(isc:iec+shift,jsc:jec,:,:)
       y(isc:iec,jsc:jec+shift,:,:) = global2(isc:iec,jsc:jec+shift,:,:)
    end do
    x1 = x; x2 =x; y1 = y; y2 = y   

    !--- call mpp_update_domains to update domain 
    do n = 1, ntiles_on_pe
       call mpp_update_domains( x(:,:,:,n), y(:,:,:,n), domain, gridtype=CGRID_NE, &
                                bufferx=bufferx(:,n), buffery=buffery(:,n), tile_count=n )
    end do

    !--- multiple update
    do n = 1, ntiles_on_pe
       call mpp_update_domains( x1(:,:,:,n), y1(:,:,:,n), domain, gridtype=CGRID_NE, &
                                bufferx=bufferx1(:,n), buffery=buffery1(:,n), complete=.false., tile_count=n )
       call mpp_update_domains( x2(:,:,:,n), y2(:,:,:,n), domain, gridtype=CGRID_NE, &
                                bufferx=bufferx2(:,n), buffery=buffery2(:,n), complete=.true., tile_count=n )
    end do

    !--- fill up the value at halo points and compare the value at buffer.
    do n = 1, ntiles_on_pe  

       !--- comparing the buffer.
       noverlap1 = mpp_get_refine_overlap_number(domain, tile_count=n)
       call mpp_get_mosaic_refine_overlap(domain, isMe1, ieMe1, jsMe1, jeMe1, isNb1, ieNb1, jsNb1, jeNb1, &
                                          dirMe1, rotation1, tile_count = n, position = EAST)             
       pos = 0
       do m = 1, noverlap1 
          select case( rotation1(m) )
          case (ZERO)     
             do k = 1, nz
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. bufferx(pos,n) ) then
                         write(stdunit,111)'x','CGRID '//type,mpp_pe(),i,j,k,bufferx(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. bufferx1(pos,n) ) then
                         write(stdunit,111)'x1','CGRID '//type,mpp_pe(),i,j,k,bufferx1(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X1")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. bufferx2(pos,n) ) then
                         write(stdunit,111)'x2','CGRID '//type,mpp_pe(),i,j,k,bufferx2(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X2")
                      end if
                   end do
                end do
             end do
          case (NINETY) ! S->E, N->W, u->-v, v->u    
             do k = 1, nz
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx(pos,n) ) then
                         write(stdunit,111)'x','CGRID '//type,mpp_pe(),i,j,k,bufferx(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx1(pos,n) ) then
                         write(stdunit,111)'x1','CGRID '//type,mpp_pe(),i,j,k,bufferx1(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X1")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx2(pos,n) ) then
                         write(stdunit,111)'x2','CGRID '//type,mpp_pe(),i,j,k,bufferx2(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X2")
                      end if
                   end do
                end do
             end do
          case (MINUS_NINETY) ! S->E, N->W, u->-v, v->u    
             do k = 1, nz
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(-global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx(pos,n) ) then
                         write(stdunit,111)'x','CGRID '//type,mpp_pe(),i,j,k,bufferx(pos,n),-global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X")
                      end if
                      if(-global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx1(pos,n) ) then
                         write(stdunit,111)'x1','CGRID '//type,mpp_pe(),i,j,k,bufferx1(pos,n),-global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X1")
                      end if
                      if(-global2_all(i,j,k,from_tile1(m,n)) .NE. bufferx2(pos,n) ) then
                         write(stdunit,111)'x2','CGRID '//type,mpp_pe(),i,j,k,bufferx2(pos,n),-global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" X2")
                      end if
                   end do
                end do
             end do
          end select
       end do

       call mpp_get_mosaic_refine_overlap(domain, isMe1, ieMe1, jsMe1, jeMe1, isNb1, ieNb1, jsNb1, jeNb1, &
            dirMe1, rotation1, tile_count = n, position = NORTH) 
       pos = 0
       do m = 1, noverlap1
          select case( rotation1(m) )
          case (ZERO)               
             do k = 1, nz  
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. buffery(pos,n) ) then
                         write(stdunit,111)'y','CGRID '//type,mpp_pe(),i,j,k,buffery(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. buffery1(pos,n) ) then
                         write(stdunit,111)'y1','CGRID '//type,mpp_pe(),i,j,k,buffery1(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y1")
                      end if
                      if(global2_all(i,j,k,from_tile1(m,n)) .NE. buffery2(pos,n) ) then
                         write(stdunit,111)'y2','CGRID '//type,mpp_pe(),i,j,k,buffery2(pos,n),global2_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y2")
                      end if
                   end do
                end do
             end do
          case (NINETY) ! S->E, N->W, u->-v, v->u    
             do k = 1, nz
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(-global1_all(i,j,k,from_tile1(m,n)) .NE. buffery(pos,n) ) then
                         write(stdunit,111)'y','CGRID '//type,mpp_pe(),i,j,k,buffery(pos,n),-global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y")
                      end if
                      if(-global1_all(i,j,k,from_tile1(m,n)) .NE. buffery1(pos,n) ) then
                         write(stdunit,111)'y1','CGRID '//type,mpp_pe(),i,j,k,buffery1(pos,n),-global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y1")
                      end if
                      if(-global1_all(i,j,k,from_tile1(m,n)) .NE. buffery2(pos,n) ) then
                         write(stdunit,111)'y2','CGRID '//type,mpp_pe(),i,j,k,buffery2(pos,n),-global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y2")
                      end if
                   end do
                end do
             end do
          case (MINUS_NINETY) ! S->E, N->W, u->-v, v->u    
             do k = 1, nz
                do j = jsNb1(m), jeNb1(m)
                   do i = isNb1(m), ieNb1(m)
                      pos = pos + 1
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. buffery(pos,n) ) then
                         write(stdunit,111)'y','CGRID '//type,mpp_pe(),i,j,k,buffery(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. buffery1(pos,n) ) then
                         write(stdunit,111)'y1','CGRID '//type,mpp_pe(),i,j,k,buffery1(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y1")
                      end if
                      if(global1_all(i,j,k,from_tile1(m,n)) .NE. buffery2(pos,n) ) then
                         write(stdunit,111)'y2','CGRID '//type,mpp_pe(),i,j,k,buffery2(pos,n),global1_all(i,j,k,from_tile1(m,n))
                         call mpp_error(FATAL,"test_refined_mosaic: mismatch between buffer and actual data: CGRID "//type//" Y2")
                      end if
                   end do
                end do
             end do
          end select
       end do
       !--- fill the halo data and compare
       select case(type)
       case('Refined-Four-Tile', 'Refined-Symmetric-Four-Tile')
          te = 0; ts = 0; tn = 0; tw = 0
          select case(tiles(n))       
          case(1)
             tsw = 4 
          case(2)
             tsw = 3; tn = 4; ts = 4
          case(3)
             tsw = 2
          case(4)
             tsw = 1; ts = 2; tn = 2
          end select
          tse = tsw; tnw = tsw; tne = tsw
          call fill_regular_refinement_halo( global1(:,:,:,n), global1_all, ni, nj, tiles(n), &
               te, tse, ts, tsw, tw, tnw, tn, tne, shift, 0 )
          call fill_regular_refinement_halo( global2(:,:,:,n), global2_all, ni, nj, tiles(n), &
               te, tse, ts, tsw, tw, tnw, tn, tne, 0, shift )
       case('Refined-Cubic-Grid')       
          call fill_cubicgrid_refined_halo(global1(:,:,:,n), global1_all, global2_all, ni, nj, tiles(n), 1, 0, 1, -1 )
          call fill_cubicgrid_refined_halo(global2(:,:,:,n), global2_all, global1_all, ni, nj, tiles(n), 0, 1, -1, 1 )
       end select

       call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
       write(type2, *)"CGRID ", type
       if(ntiles_on_pe>1) write(type2, *)trim(type2), " at tile_count = ", tiles(n)
       call compare_checksums( x (isd:ied+shift,jsd:jed,:,n), global1(isd:ied+shift,jsd:jed,:,n), trim(type2)//' X' )
       call compare_checksums( x1(isd:ied+shift,jsd:jed,:,n), global1(isd:ied+shift,jsd:jed,:,n), trim(type2)//' X1')
       call compare_checksums( x2(isd:ied+shift,jsd:jed,:,n), global1(isd:ied+shift,jsd:jed,:,n), trim(type2)//' X2')
       write(type2, *)"CGRID ", type
       if(ntiles_on_pe>1) write(type2, *)trim(type2), " at tile_count = ", tiles(n)
       call compare_checksums( y (isd:ied,jsd:jed+shift,:,n), global2(isd:ied,jsd:jed+shift,:,n), trim(type2)//' Y' )
       call compare_checksums( y1(isd:ied,jsd:jed+shift,:,n), global2(isd:ied,jsd:jed+shift,:,n), trim(type2)//' Y1' )
       call compare_checksums( y2(isd:ied,jsd:jed+shift,:,n), global2(isd:ied,jsd:jed+shift,:,n), trim(type2)//' Y2' )
    end do  

  111 format('For variable ', a, ', type = ', a, ', at pe = ', i3, ', at neighbor point (',i3,',',i3,',',i3, &
             '), failed value = ', f14.9, ', but the value should be ', f14.9 )


  end subroutine test_refined_mosaic

  !#######################################################################################
  subroutine test_get_boundary(type)
     character(len=*), intent(in)  :: type

     type(domain2D)       :: domain
     integer              :: ntiles, num_contact, npes_per_tile, ntile_per_pe, layout(2)
     integer              :: n, l, isc, iec, jsc, jec, ism, iem, jsm, jem
     integer, allocatable, dimension(:)       :: tile, ni, nj, pe_start, pe_end
     integer, allocatable, dimension(:,:)     :: layout2D, global_indices
     real,    allocatable, dimension(:,:,:)   :: ebuffer,   sbuffer,   wbuffer,   nbuffer
     real,    allocatable, dimension(:,:,:)   :: ebuffer1,  sbuffer1,  wbuffer1,  nbuffer1
     real,    allocatable, dimension(:,:,:)   :: ebuffer2,  sbuffer2,  wbuffer2,  nbuffer2
     real,    allocatable, dimension(:,:,:)   :: ebound,    sbound,    wbound,    nbound
     real,    allocatable, dimension(:,:,:)   :: ebufferx,  sbufferx,  wbufferx,  nbufferx
     real,    allocatable, dimension(:,:,:)   :: ebufferx1, sbufferx1, wbufferx1, nbufferx1
     real,    allocatable, dimension(:,:,:)   :: ebufferx2, sbufferx2, wbufferx2, nbufferx2
     real,    allocatable, dimension(:,:,:)   :: eboundx,   sboundx,   wboundx,   nboundx
     real,    allocatable, dimension(:,:,:)   :: ebuffery,  sbuffery,  wbuffery,  nbuffery
     real,    allocatable, dimension(:,:,:)   :: ebuffery1, sbuffery1, wbuffery1, nbuffery1
     real,    allocatable, dimension(:,:,:)   :: ebuffery2, sbuffery2, wbuffery2, nbuffery2
     real,    allocatable, dimension(:,:,:)   :: eboundy,   sboundy,   wboundy,   nboundy
     real,    allocatable, dimension(:,:,:,:) :: global_all, global1_all, global2_all
     real,    allocatable, dimension(:,:,:,:) :: global, global1, global2
     real,    allocatable, dimension(:,:,:,:) :: x, x1, x2, y, y1, y2

     !--- check the type
    select case(type)     
    case ( 'Four-Tile' ) !--- cyclic along both x- and y-direction. 
       ntiles = 4
       num_contact = 8
    case ( 'Cubic-Grid' )
       ntiles = 6
       num_contact = 12
       if( nx .NE. ny) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cubic_grid mosaic, nx should equal ny, '//&
                   'No test is done for Cubic-Grid mosaic. ' )
          return
       end if
    case default
       call mpp_error(FATAL, 'TEST_MPP_DOMAINS: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    allocate(ni(ntiles), nj(ntiles))
    ni(:) = nx; nj(:) = ny
    if( mod(npes, ntiles) == 0 ) then
       npes_per_tile = npes/ntiles
       write(stdout(),*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                       '", each tile will be distributed over ', npes_per_tile, ' processors.'
       ntile_per_pe = 1
       allocate(tile(ntile_per_pe))
       tile = pe/npes_per_tile+1
       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       do n = 1, ntiles
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do
    else if ( mod(ntiles, npes) == 0 ) then
       ntile_per_pe = ntiles/npes 
       write(stdout(),*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                        '", there will be ', ntile_per_pe, ' tiles on each processor.'
       allocate(tile(ntile_per_pe))
       do n = 1, ntile_per_pe
          tile(n) = pe*ntile_per_pe + n
       end do
       do n = 1, ntiles
          pe_start(n) = (n-1)/ntile_per_pe
          pe_end(n)   = pe_start(n)
       end do
       layout = 1
    else
       call mpp_error(NOTE,'TEST_MPP_DOMAINS: npes should be multiple of ntiles or ' // &
            'ntiles should be multiple of npes. No test is done for '//trim(type) )       
       return
    end if
 
    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

     select case(type)
     case("Four-Tile")
        call define_fourtile_mosaic(type, domain, (/nx,nx,nx,nx/), (/ny,ny,ny,ny/), global_indices, &
                                    layout2D, pe_start, pe_end, .true. )  
     case("Cubic-Grid")
        call define_cubic_mosaic(type, domain, ni, nj, global_indices, layout2D, pe_start, pe_end )
     end select

    !--- Test the get_boundary of the data at C-cell center. 
    allocate(global_all(1:nx+1,1:ny+1,nz, ntiles) ) 
    allocate(global(1:nx+1,1:ny+1,nz, ntile_per_pe) )      
    global = 0
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny+1
             do i = 1, nx+1
                global_all(i,j,k,l) = l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global(:,:,:,n) = global_all(:,:,:,tile(n))
    end do

    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x (ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( x1(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( x2(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    x = 0.
    x(isc:iec+1,jsc:jec+1,:,:) = global(isc:iec+1,jsc:jec+1,:,:)
    x1 = x; x2 = x*10

    !--- buffer allocation
    allocate(ebuffer(jec-jsc+2, nz, ntile_per_pe), wbuffer(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffer(iec-isc+2, nz, ntile_per_pe), nbuffer(iec-isc+2, nz, ntile_per_pe))
    allocate(ebuffer1(jec-jsc+2, nz, ntile_per_pe), wbuffer1(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffer1(iec-isc+2, nz, ntile_per_pe), nbuffer1(iec-isc+2, nz, ntile_per_pe))
    allocate(ebuffer2(jec-jsc+2, nz, ntile_per_pe), wbuffer2(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffer2(iec-isc+2, nz, ntile_per_pe), nbuffer2(iec-isc+2, nz, ntile_per_pe))
    allocate(ebound(jec-jsc+2, nz, ntile_per_pe), wbound(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbound(iec-isc+2, nz, ntile_per_pe), nbound(iec-isc+2, nz, ntile_per_pe))
    do n = 1, ntile_per_pe 
       call mpp_get_boundary(x(:,:,:,n), domain, ebuffer=ebuffer(:,:,n), sbuffer=sbuffer(:,:,n), wbuffer=wbuffer(:,:,n), &
                             nbuffer=nbuffer(:,:,n), position=CORNER, tile_count=n  )
    end do

    !--- multiple variable 
    do n = 1, ntile_per_pe 
       call mpp_get_boundary(x1(:,:,:,n), domain, ebuffer=ebuffer1(:,:,n), sbuffer=sbuffer1(:,:,n), wbuffer=wbuffer1(:,:,n), &
                             nbuffer=nbuffer1(:,:,n), position=CORNER, tile_count=n, complete = .false.  )
       call mpp_get_boundary(x2(:,:,:,n), domain, ebuffer=ebuffer2(:,:,n), sbuffer=sbuffer2(:,:,n), wbuffer=wbuffer2(:,:,n), &
                             nbuffer=nbuffer2(:,:,n), position=CORNER, tile_count=n, complete = .true.  )
    end do    

    !--- compare the buffer.
    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), ebound(:,:,n), sbound(:,:,n), wbound(:,:,n), nbound(:,:,n) )  
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global_all, global_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, ebound(:,:,n), sbound(:,:,n), wbound(:,:,n), nbound(:,:,n) )  
       end do
    end select

    call compare_checksums( ebound, ebuffer(:,:,:),  "east bound of "//trim(type) )
    call compare_checksums( sbound, sbuffer(:,:,:),  "south bound of "//trim(type) )
    call compare_checksums( wbound, wbuffer(:,:,:),  "west bound of "//trim(type) )
    call compare_checksums( nbound, nbuffer(:,:,:),  "north bound of "//trim(type) )
    call compare_checksums( ebound, ebuffer1(:,:,:),  "east bound of "//trim(type)//" X1" )
    call compare_checksums( sbound, sbuffer1(:,:,:),  "south bound of "//trim(type)//" X1" )
    call compare_checksums( wbound, wbuffer1(:,:,:),  "west bound of "//trim(type)//" X1" )
    call compare_checksums( nbound, nbuffer1(:,:,:),  "north bound of "//trim(type)//" X1" )
    call compare_checksums( ebound*10, ebuffer2(:,:,:),  "east bound of "//trim(type)//" X2" )
    call compare_checksums( sbound*10, sbuffer2(:,:,:),  "south bound of "//trim(type)//" X2" )
    call compare_checksums( wbound*10, wbuffer2(:,:,:),  "west bound of "//trim(type)//" X2" )
    call compare_checksums( nbound*10, nbuffer2(:,:,:),  "north bound of "//trim(type)//" X2" )

    !--- release memory
    deallocate(global, global_all, x, x1, x2)
    deallocate(ebuffer, sbuffer, wbuffer, nbuffer)
    deallocate(ebuffer1, sbuffer1, wbuffer1, nbuffer1)
    deallocate(ebuffer2, sbuffer2, wbuffer2, nbuffer2)
    deallocate(ebound, sbound, wbound, nbound )

    !-------------------------------------------------------------------------------------------
    !
    !             Test SCALAR_PAIR BGRID
    !
    !-------------------------------------------------------------------------------------------
    allocate(global1_all(1:nx+1,1:ny+1,nz, ntiles) ) 
    allocate(global2_all(1:nx+1,1:ny+1,nz, ntiles) ) 
    allocate(global1(1:nx+1,1:ny+1,nz, ntile_per_pe) )   
    allocate(global2(1:nx+1,1:ny+1,nz, ntile_per_pe) )      
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny+1
             do i = 1, nx+1
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global1(:,:,:,n) = global1_all(:,:,:,tile(n))
       global2(:,:,:,n) = global2_all(:,:,:,tile(n))
    end do
    allocate( x (ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( x1(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( x2(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( y (ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( y1(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( y2(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    x = 0.; y = 0
    x(isc:iec+1,jsc:jec+1,:,:) = global1(isc:iec+1,jsc:jec+1,:,:)
    y(isc:iec+1,jsc:jec+1,:,:) = global2(isc:iec+1,jsc:jec+1,:,:)
    x1 = x; x2 = x*10
    y1 = y; y2 = y*10

    !--- buffer allocation
    allocate(ebufferx(jec-jsc+2, nz, ntile_per_pe), wbufferx(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbufferx(iec-isc+2, nz, ntile_per_pe), nbufferx(iec-isc+2, nz, ntile_per_pe))
    allocate(ebufferx1(jec-jsc+2, nz, ntile_per_pe), wbufferx1(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbufferx1(iec-isc+2, nz, ntile_per_pe), nbufferx1(iec-isc+2, nz, ntile_per_pe))
    allocate(ebufferx2(jec-jsc+2, nz, ntile_per_pe), wbufferx2(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbufferx2(iec-isc+2, nz, ntile_per_pe), nbufferx2(iec-isc+2, nz, ntile_per_pe))
    allocate(eboundx(jec-jsc+2, nz, ntile_per_pe), wboundx(jec-jsc+2, nz, ntile_per_pe))
    allocate(sboundx(iec-isc+2, nz, ntile_per_pe), nboundx(iec-isc+2, nz, ntile_per_pe))
    allocate(ebuffery(jec-jsc+2, nz, ntile_per_pe), wbuffery(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery(iec-isc+2, nz, ntile_per_pe), nbuffery(iec-isc+2, nz, ntile_per_pe))
    allocate(ebuffery1(jec-jsc+2, nz, ntile_per_pe), wbuffery1(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery1(iec-isc+2, nz, ntile_per_pe), nbuffery1(iec-isc+2, nz, ntile_per_pe))
    allocate(ebuffery2(jec-jsc+2, nz, ntile_per_pe), wbuffery2(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery2(iec-isc+2, nz, ntile_per_pe), nbuffery2(iec-isc+2, nz, ntile_per_pe))
    allocate(eboundy(jec-jsc+2, nz, ntile_per_pe), wboundy(jec-jsc+2, nz, ntile_per_pe))
    allocate(sboundy(iec-isc+2, nz, ntile_per_pe), nboundy(iec-isc+2, nz, ntile_per_pe))

    do n = 1, ntile_per_pe 
       call mpp_get_boundary(x(:,:,:,n), y(:,:,:,n), domain, ebufferx=ebufferx(:,:,n), sbufferx=sbufferx(:,:,n), &
                             wbufferx=wbufferx(:,:,n), nbufferx=nbufferx(:,:,n), ebuffery=ebuffery(:,:,n),       &
                             sbuffery=sbuffery(:,:,n), wbuffery=wbuffery(:,:,n), nbuffery=nbuffery(:,:,n),       &
                             gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR  )
    end do

    do n = 1, ntile_per_pe 
       call mpp_get_boundary(x1(:,:,:,n), y1(:,:,:,n), domain, ebufferx=ebufferx1(:,:,n), sbufferx=sbufferx1(:,:,n), &
                             wbufferx=wbufferx1(:,:,n), nbufferx=nbufferx1(:,:,n), ebuffery=ebuffery1(:,:,n),       &
                             sbuffery=sbuffery1(:,:,n), wbuffery=wbuffery1(:,:,n), nbuffery=nbuffery1(:,:,n),       &
                             gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR, complete = .false.  )
       call mpp_get_boundary(x2(:,:,:,n), y2(:,:,:,n), domain, ebufferx=ebufferx2(:,:,n), sbufferx=sbufferx2(:,:,n), &
                             wbufferx=wbufferx2(:,:,n), nbufferx=nbufferx2(:,:,n), ebuffery=ebuffery2(:,:,n),       &
                             sbuffery=sbuffery2(:,:,n), wbuffery=wbuffery2(:,:,n), nbuffery=nbuffery2(:,:,n),       &
                             gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR, complete = .true.  )
    end do

    !--- compare the buffer.
    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )  
          call fill_four_tile_bound(global2_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )  
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all, global2_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )  
          call fill_cubic_grid_bound(global2_all, global1_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )  
       end do
    end select

    call compare_checksums( eboundx, ebufferx(:,:,:),   "east bound of SCALAR_PAIR BGRID " //trim(type)//" X" )
    call compare_checksums( sboundx, sbufferx(:,:,:),   "south bound of SCALAR_PAIR BGRID "//trim(type)//" X" )
    call compare_checksums( wboundx, wbufferx(:,:,:),   "west bound of SCALAR_PAIR BGRID " //trim(type)//" X" )
    call compare_checksums( nboundx, nbufferx(:,:,:),   "north bound of SCALAR_PAIR BGRID "//trim(type)//" X" )
    call compare_checksums( eboundy, ebuffery(:,:,:),   "east bound of SCALAR_PAIR BGRID " //trim(type)//" Y" )
    call compare_checksums( sboundy, sbuffery(:,:,:),   "south bound of SCALAR_PAIR BGRID "//trim(type)//" Y" )
    call compare_checksums( wboundy, wbuffery(:,:,:),   "west bound of SCALAR_PAIR BGRID " //trim(type)//" Y" )
    call compare_checksums( nboundy, nbuffery(:,:,:),   "north bound of SCALAR_PAIR BGRID "//trim(type)//" Y" )
    call compare_checksums( eboundx, ebufferx1(:,:,:),  "east bound of SCALAR_PAIR BGRID " //trim(type)//" X1" )
    call compare_checksums( sboundx, sbufferx1(:,:,:),  "south bound of SCALAR_PAIR BGRID "//trim(type)//" X1" )
    call compare_checksums( wboundx, wbufferx1(:,:,:),  "west bound of SCALAR_PAIR BGRID " //trim(type)//" X1" )
    call compare_checksums( nboundx, nbufferx1(:,:,:),  "north bound of SCALAR_PAIR BGRID "//trim(type)//" X1" )
    call compare_checksums( eboundy, ebuffery1(:,:,:),  "east bound of SCALAR_PAIR BGRID " //trim(type)//" Y1" )
    call compare_checksums( sboundy, sbuffery1(:,:,:),  "south bound of SCALAR_PAIR BGRID "//trim(type)//" Y1" )
    call compare_checksums( wboundy, wbuffery1(:,:,:),  "west bound of SCALAR_PAIR BGRID " //trim(type)//" Y1" )
    call compare_checksums( nboundy, nbuffery1(:,:,:),  "north bound of SCALAR_PAIR BGRID "//trim(type)//" Y1" )

    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )  
          call fill_four_tile_bound(global2_all*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )  
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all*10, global2_all*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )  
          call fill_cubic_grid_bound(global2_all*10, global1_all*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )  
       end do
    end select

    call compare_checksums( eboundx, ebufferx2(:,:,:),  "east bound of SCALAR_PAIR BGRID " //trim(type)//" X2" )
    call compare_checksums( sboundx, sbufferx2(:,:,:),  "south bound of SCALAR_PAIR BGRID "//trim(type)//" X2" )
    call compare_checksums( wboundx, wbufferx2(:,:,:),  "west bound of SCALAR_PAIR BGRID " //trim(type)//" X2" )
    call compare_checksums( nboundx, nbufferx2(:,:,:),  "north bound of SCALAR_PAIR BGRID "//trim(type)//" X2" )
    call compare_checksums( eboundy, ebuffery2(:,:,:),  "east bound of SCALAR_PAIR BGRID " //trim(type)//" Y2" )
    call compare_checksums( sboundy, sbuffery2(:,:,:),  "south bound of SCALAR_PAIR BGRID "//trim(type)//" Y2" )
    call compare_checksums( wboundy, wbuffery2(:,:,:),  "west bound of SCALAR_PAIR BGRID " //trim(type)//" Y2" )
    call compare_checksums( nboundy, nbuffery2(:,:,:),  "north bound of SCALAR_PAIR BGRID "//trim(type)//" Y2" )

    !--- release memory
    deallocate(global1, global1_all, global2, global2_all)
    deallocate(x, y, x1, y1, x2, y2)
    deallocate(ebufferx, sbufferx, wbufferx, nbufferx)
    deallocate(ebufferx1, sbufferx1, wbufferx1, nbufferx1)
    deallocate(ebufferx2, sbufferx2, wbufferx2, nbufferx2)
    deallocate(ebuffery, sbuffery, wbuffery, nbuffery)
    deallocate(ebuffery1, sbuffery1, wbuffery1, nbuffery1)
    deallocate(ebuffery2, sbuffery2, wbuffery2, nbuffery2)
    deallocate(eboundx, sboundx, wboundx, nboundx )    
    deallocate(eboundy, sboundy, wboundy, nboundy )  

    !-------------------------------------------------------------------------------------------
    !
    !             Test VECTOR CGRID
    !
    !-------------------------------------------------------------------------------------------
    allocate(global1_all(1:nx+1,1:ny,  nz, ntiles) ) 
    allocate(global2_all(1:nx,  1:ny+1,nz, ntiles) ) 
    allocate(global1(1:nx+1,1:ny,  nz, ntile_per_pe) )   
    allocate(global2(1:nx,  1:ny+1,nz, ntile_per_pe) )      
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx+1
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
          do j = 1, ny+1
             do i = 1, nx
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global1(:,:,:,n) = global1_all(:,:,:,tile(n))
       global2(:,:,:,n) = global2_all(:,:,:,tile(n))
    end do
    allocate( x (ism:iem+1,jsm:jem,  nz, ntile_per_pe) )
    allocate( x1(ism:iem+1,jsm:jem,  nz, ntile_per_pe) )
    allocate( x2(ism:iem+1,jsm:jem,  nz, ntile_per_pe) )
    allocate( y (ism:iem,  jsm:jem+1,nz, ntile_per_pe) )
    allocate( y1(ism:iem,  jsm:jem+1,nz, ntile_per_pe) )
    allocate( y2(ism:iem,  jsm:jem+1,nz, ntile_per_pe) )
    x = 0.; y = 0
    x(isc:iec+1,jsc:jec,  :,:) = global1(isc:iec+1,jsc:jec,  :,:)
    y(isc:iec,  jsc:jec+1,:,:) = global2(isc:iec,  jsc:jec+1,:,:)
    x1 = x; x2 = x*10
    y1 = y; y2 = y*10

    !--- buffer allocation
    allocate(ebufferx(jec-jsc+1, nz, ntile_per_pe), wbufferx(jec-jsc+1, nz, ntile_per_pe))
    allocate(sbufferx(iec-isc+2, nz, ntile_per_pe), nbufferx(iec-isc+2, nz, ntile_per_pe))
    allocate(ebufferx1(jec-jsc+1, nz, ntile_per_pe), wbufferx1(jec-jsc+1, nz, ntile_per_pe))
    allocate(sbufferx1(iec-isc+2, nz, ntile_per_pe), nbufferx1(iec-isc+2, nz, ntile_per_pe))
    allocate(ebufferx2(jec-jsc+1, nz, ntile_per_pe), wbufferx2(jec-jsc+1, nz, ntile_per_pe))
    allocate(sbufferx2(iec-isc+2, nz, ntile_per_pe), nbufferx2(iec-isc+2, nz, ntile_per_pe))
    allocate(ebuffery(jec-jsc+2, nz, ntile_per_pe), wbuffery(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery(iec-isc+1, nz, ntile_per_pe), nbuffery(iec-isc+1, nz, ntile_per_pe))
    allocate(ebuffery1(jec-jsc+2, nz, ntile_per_pe), wbuffery1(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery1(iec-isc+1, nz, ntile_per_pe), nbuffery1(iec-isc+1, nz, ntile_per_pe))
    allocate(ebuffery2(jec-jsc+2, nz, ntile_per_pe), wbuffery2(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery2(iec-isc+1, nz, ntile_per_pe), nbuffery2(iec-isc+1, nz, ntile_per_pe))
    allocate(eboundx(jec-jsc+1, nz, ntile_per_pe), wboundx(jec-jsc+1, nz, ntile_per_pe))
    allocate(sboundy(iec-isc+1, nz, ntile_per_pe), nboundy(iec-isc+1, nz, ntile_per_pe))

    do n = 1, ntile_per_pe 
       call mpp_get_boundary(x(:,:,:,n), y(:,:,:,n), domain, ebufferx=ebufferx(:,:,n), wbufferx=wbufferx(:,:,n), &
                             sbuffery=sbuffery(:,:,n), nbuffery=nbuffery(:,:,n), gridtype=CGRID_NE, tile_count=n  )
    end do

    do n = 1, ntile_per_pe 
       call mpp_get_boundary(x1(:,:,:,n), y1(:,:,:,n), domain, ebufferx=ebufferx1(:,:,n), wbufferx=wbufferx1(:,:,n), &
                             sbuffery=sbuffery1(:,:,n), nbuffery=nbuffery1(:,:,n), gridtype=CGRID_NE, tile_count=n,  &
                             complete = .false.  )
       call mpp_get_boundary(x2(:,:,:,n), y2(:,:,:,n), domain, ebufferx=ebufferx2(:,:,n), wbufferx=wbufferx2(:,:,n), &
                             sbuffery=sbuffery2(:,:,n), nbuffery=nbuffery2(:,:,n), gridtype=CGRID_NE, tile_count=n,  &
                             complete = .true.  )
    end do

    !--- compare the buffer.
    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all, isc, iec, jsc, jec, 1, 0, &
               tile(n), ebound=eboundx(:,:,n), wbound=wboundx(:,:,n) )  
          call fill_four_tile_bound(global2_all, isc, iec, jsc, jec, 0, 1, &
               tile(n), sbound=sboundy(:,:,n), nbound=nboundy(:,:,n) )  
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all, global2_all, isc, iec, jsc, jec, 1, 0, &
               tile(n), 1, -1, ebound=eboundx(:,:,n), wbound=wboundx(:,:,n)  )  
          call fill_cubic_grid_bound(global2_all, global1_all, isc, iec, jsc, jec, 0, 1, &
               tile(n), -1, 1, sbound=sboundy(:,:,n), nbound=nboundy(:,:,n) )  
       end do
    end select

    call compare_checksums( eboundx, ebufferx(:,:,:),   "east bound of CGRID " //trim(type)//" X" )
    call compare_checksums( wboundx, wbufferx(:,:,:),   "west bound of CGRID " //trim(type)//" X" )
    call compare_checksums( sboundy, sbuffery(:,:,:),   "south bound of CGRID "//trim(type)//" Y" )
    call compare_checksums( nboundy, nbuffery(:,:,:),   "north bound of CGRID "//trim(type)//" Y" )
    call compare_checksums( eboundx, ebufferx1(:,:,:),  "east bound of CGRID " //trim(type)//" X1" )
    call compare_checksums( wboundx, wbufferx1(:,:,:),  "west bound of CGRID " //trim(type)//" X1" )
    call compare_checksums( sboundy, sbuffery1(:,:,:),  "south bound of CGRID "//trim(type)//" Y1" )
    call compare_checksums( nboundy, nbuffery1(:,:,:),  "north bound of CGRID "//trim(type)//" Y1" )

    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all*10, isc, iec, jsc, jec, 1, 0, &
               tile(n), ebound=eboundx(:,:,n), wbound=wboundx(:,:,n) )  
          call fill_four_tile_bound(global2_all*10, isc, iec, jsc, jec, 0, 1, &
               tile(n), sbound=sboundy(:,:,n), nbound=nboundy(:,:,n) )  
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all*10, global2_all*10, isc, iec, jsc, jec, 1, 0, &
               tile(n), 1, -1, ebound=eboundx(:,:,n), wbound=wboundx(:,:,n) )  
          call fill_cubic_grid_bound(global2_all*10, global1_all*10, isc, iec, jsc, jec, 0, 1, &
               tile(n), -1, 1, sbound=sboundy(:,:,n), nbound=nboundy(:,:,n) )  
       end do
    end select

    call compare_checksums( eboundx, ebufferx2(:,:,:),  "east bound of CGRID " //trim(type)//" X2" )
    call compare_checksums( wboundx, wbufferx2(:,:,:),  "west bound of CGRID " //trim(type)//" X2" )
    call compare_checksums( sboundy, sbuffery2(:,:,:),  "south bound of CGRID "//trim(type)//" Y2" )
    call compare_checksums( nboundy, nbuffery2(:,:,:),  "north bound of CGRID "//trim(type)//" Y2" )

    !--- release memory
    deallocate(global1, global1_all, global2, global2_all)
    deallocate(x, y, x1, y1, x2, y2)
    deallocate(ebufferx, sbufferx, wbufferx, nbufferx)
    deallocate(ebufferx1, sbufferx1, wbufferx1, nbufferx1)
    deallocate(ebufferx2, sbufferx2, wbufferx2, nbufferx2)
    deallocate(ebuffery, sbuffery, wbuffery, nbuffery)
    deallocate(ebuffery1, sbuffery1, wbuffery1, nbuffery1)
    deallocate(ebuffery2, sbuffery2, wbuffery2, nbuffery2)
    deallocate(eboundx, sboundy, wboundx, nboundy )    

  end subroutine test_get_boundary

  !######################################################################################
  subroutine define_fourtile_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end, symmetry )
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    logical,        intent(in)    :: symmetry
    integer, dimension(8)         :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(8)         :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)

    ntiles = 4
    num_contact = 8
    if(size(pe_start(:)) .NE. 4 .OR. size(pe_end(:)) .NE. 4 ) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of pe_start and pe_end should be 4")
    if(size(global_indices,1) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of first dimension of global_indices should be 4")
    if(size(global_indices,2) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of second dimension of global_indices should be 4")
    if(size(layout,1) .NE. 2) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of first dimension of layout should be 2")
    if(size(layout,2) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of second dimension of layout should be 4")
    if(size(ni(:)) .NE. 4 .OR. size(nj(:)) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of ni and nj should be 4")

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1; tile2(1) = 2
    istart1(1) = ni(1); iend1(1) = ni(1); jstart1(1) = 1;     jend1(1) = nj(1)
    istart2(1) = 1;     iend2(1) = 1;     jstart2(1) = 1;     jend2(1) = nj(2)
    !--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;     iend1(2) = ni(1); jstart1(2) = 1;     jend1(2) = 1    
    istart2(2) = 1;     iend2(2) = ni(3); jstart2(2) = nj(3); jend2(2) = nj(3)
    !--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
    tile1(3) = 1; tile2(3) = 2
    istart1(3) = 1;     iend1(3) = 1;     jstart1(3) = 1;     jend1(3) = nj(1)
    istart2(3) = ni(2); iend2(3) = ni(2); jstart2(3) = 1;     jend2(3) = nj(2)
    !--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH) 
    tile1(4) = 1; tile2(4) = 3
    istart1(4) = 1;     iend1(4) = ni(1); jstart1(4) = nj(1); jend1(4) = nj(1)
    istart2(4) = 1;     iend2(4) = ni(3); jstart2(4) = 1;     jend2(4) = 1    
    !--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
    tile1(5) = 2; tile2(5) = 4
    istart1(5) = 1;     iend1(5) = ni(2); jstart1(5) = 1;     jend1(5) = 1    
    istart2(5) = 1;     iend2(5) = ni(4); jstart2(5) = nj(4); jend2(5) = nj(4)
    !--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = 1;     iend1(6) = ni(2); jstart1(6) = nj(2); jend1(6) = nj(2)
    istart2(6) = 1;     iend2(6) = ni(4); jstart2(6) = 1;     jend2(6) = 1    
    !--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST) 
    tile1(7) = 3; tile2(7) = 4
    istart1(7) = ni(3); iend1(7) = ni(3); jstart1(7) = 1;     jend1(7) = nj(3)
    istart2(7) = 1;     iend2(7) = 1;     jstart2(7) = 1;     jend2(7) = nj(4)
    !--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = 1;     iend1(8) = 1;     jstart1(8) = 1;     jend1(8) = nj(3)
    istart2(8) = ni(4); iend2(8) = ni(4); jstart2(8) = 1;     jend2(8) = nj(4)
    msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1 ! make sure memory domain size is no smaller than
    msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1 ! data domain size       
    call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, tile2,       &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,          &
         pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo,    &
         name = type, memory_size = msize, symmetry = symmetry )
 
    return

  end subroutine define_fourtile_mosaic

  !#######################################################################################
  !--- define mosaic domain for cubic grid
  subroutine define_cubic_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end )
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)


    ntiles = 6
    num_contact = 12
    if(size(pe_start(:)) .NE. 6 .OR. size(pe_end(:)) .NE. 6 ) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of pe_start and pe_end should be 6")
    if(size(global_indices,1) .NE. 4) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of first dimension of global_indices should be 4")
    if(size(global_indices,2) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of second dimension of global_indices should be 6")
    if(size(layout,1) .NE. 2) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of first dimension of layout should be 2")
    if(size(layout,2) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of second dimension of layout should be 6")
    if(size(ni(:)) .NE. 6 .OR. size(nj(:)) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of ni and nj should be 6")

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1; tile2(1) = 2
    istart1(1) = ni(1);  iend1(1) = ni(1);  jstart1(1) = 1;      jend1(1) = nj(1)
    istart2(1) = 1;      iend2(1) = 1;      jstart2(1) = 1;      jend2(1) = nj(2)
    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni(1);  jstart1(2) = nj(1);  jend1(2) = nj(1)
    istart2(2) = 1;      iend2(2) = 1;      jstart2(2) = nj(3);  jend2(2) = 1
    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1; tile2(3) = 5
    istart1(3) = 1;      iend1(3) = 1;      jstart1(3) = 1;      jend1(3) = nj(1)
    istart2(3) = ni(5);  iend2(3) = 1;      jstart2(3) = nj(5);  jend2(3) = nj(5)
    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni(1);  jstart1(4) = 1;      jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni(6);  jstart2(4) = nj(6);  jend2(4) = nj(6)       
    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2; tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni(2);  jstart1(5) = nj(2);  jend1(5) = nj(2)
    istart2(5) = 1;      iend2(5) = ni(3);  jstart2(5) = 1;      jend2(5) = 1
    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni(2);  iend1(6) = ni(2);  jstart1(6) = 1;      jend1(6) = nj(2)
    istart2(6) = ni(4);  iend2(6) = 1;      jstart2(6) = 1;      jend2(6) = 1
    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;      iend1(7) = ni(2);  jstart1(7) = 1;      jend1(7) = 1
    istart2(7) = ni(6);  iend2(7) = ni(6);  jstart2(7) = nj(6);  jend2(7) = 1
    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni(3);  iend1(8) = ni(3);  jstart1(8) = 1;      jend1(8) = nj(3)
    istart2(8) = 1;      iend2(8) = 1;      jstart2(8) = 1;      jend2(8) = nj(4)
    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni(3);  jstart1(9) = nj(3);  jend1(9) = nj(3)
    istart2(9) = 1;      iend2(9) = 1;      jstart2(9) = nj(5);  jend2(9) = 1
    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni(4); jstart1(10) = nj(4); jend1(10) = nj(4)
    istart2(10) = 1;     iend2(10) = ni(5); jstart2(10) = 1;     jend2(10) = 1
    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni(4); iend1(11) = ni(4); jstart1(11) = 1;     jend1(11) = nj(4)
    istart2(11) = ni(6); iend2(11) = 1;     jstart2(11) = 1;     jend2(11) = 1
    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni(5); iend1(12) = ni(5); jstart1(12) = 1;     jend1(12) = nj(5)
    istart2(12) = 1;     iend2(12) = 1;     jstart2(12) = 1;     jend2(12) = nj(6)
    msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1 ! make sure memory domain size is no smaller than
    msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1 ! data domain size       
    call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, tile2, &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
         pe_start, pe_end, symmetry = .true., whalo=whalo, ehalo=ehalo,   &
         shalo=shalo, nhalo=nhalo, name = trim(type), memory_size = msize  )  

    return 

  end subroutine define_cubic_mosaic

  !#######################################################################################
  subroutine fill_regular_refinement_halo( data, data_all, ni, nj, tm, te, tse, ts, tsw, tw, tnw, tn, tne, ioff, joff )
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data_all
    integer, dimension(:),                intent(in)    :: ni, nj
    integer,                              intent(in)    :: tm, te, tse, ts, tsw, tw, tnw, tn, tne
    integer,                              intent(in)    :: ioff, joff


    if(te>0) data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1:nj(tm)+joff,                   :) = &
             data_all(1+ioff:ehalo+ioff,               1:nj(te)+joff,                   :,te)  ! east
    if(ts>0) data    (1:ni(tm)+ioff,                   1-shalo:0,                       :) = &
             data_all(1:ni(ts)+ioff,                   nj(ts)-shalo+1:nj(ts),           :,ts)  ! south 
    if(tw>0) data    (1-whalo:0,                       1:nj(tm)+joff,                   :) = &
             data_all(ni(tw)-whalo+1:ni(tw),           1:nj(tw)+joff,                   :,tw)  ! west
    if(tn>0) data    (1:ni(tm)+ioff,                   nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1:ni(tn)+ioff,                   1+joff:nhalo+joff,               :,tn)  ! north  
    if(tse>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1-shalo:0,                       :) = &
             data_all(1+ioff:ehalo+ioff,               nj(tse)-shalo+1:nj(tse),         :,tse) ! southeast
    if(tsw>0)data    (1-whalo:0,                       1-shalo:0,                       :) = &
             data_all(ni(tsw)-whalo+1:ni(tsw),         nj(tsw)-shalo+1:nj(tsw),         :,tsw) ! southwest
    if(tne>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1+ioff:ehalo+ioff,               1+joff:nhalo+joff,               :,tnw) ! northeast
    if(tnw>0)data    (1-whalo:0,                       nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(ni(tnw)-whalo+1:ni(tnw),         1+joff:nhalo+joff,               :,tne) ! northwest      

  end subroutine fill_regular_refinement_halo

  !##############################################################################
  ! this routine fill the halo points for the refined cubic grid. ioff and joff is used to distinguish
  ! T, C, E, or N-cell
  subroutine fill_cubicgrid_refined_halo(data, data1_all, data2_all, ni, nj, tile, ioff, joff, sign1, sign2)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data1_all, data2_all
    integer, dimension(:),                intent(in)    :: ni, nj
    integer,                              intent(in)    :: tile, ioff, joff, sign1, sign2 
    integer                                             :: lw, le, ls, ln

    if(mod(tile,2) == 0) then ! tile 2, 4, 6
       lw = tile - 1; le = tile + 2; ls = tile - 2; ln = tile + 1
       if(le > 6 ) le = le - 6
       if(ls < 1 ) ls = ls + 6
       if(ln > 6 ) ln = ln - 6
       if( nj(tile) == nj(lw) ) then
          data(1-whalo:0, 1:nj(tile)+joff, :) = data1_all(ni(lw)-whalo+1:ni(lw), 1:nj(lw)+joff, :, lw) ! west 
       end if
       if( nj(tile) == ni(le) ) then
          do i = 1, ehalo 
             data(ni(tile)+i+ioff, 1:nj(tile)+joff, :)    = sign1*data2_all(ni(le)+joff:1:-1, i+ioff, :, le) ! east 
          end do
       end if
       if(ni(tile) == nj(ls) ) then
          do i = 1, shalo 
             data(1:ni(tile)+ioff, 1-i, :)     = sign2*data2_all(ni(ls)-i+1, nj(ls)+ioff:1:-1, :, ls) ! south 
          end do
       end if
       if(ni(tile) == ni(ln) ) then
          data(1:ni(tile)+ioff, nj(tile)+1+joff:nj(tile)+nhalo+joff, :) = data1_all(1:ni(ln)+ioff, 1+joff:nhalo+joff, :, ln) ! north
       end if
    else ! tile 1, 3, 5
       lw = tile - 2; le = tile + 1; ls = tile - 1; ln = tile + 2
       if(lw < 1 ) lw = lw + 6
       if(ls < 1 ) ls = ls + 6
       if(ln > 6 ) ln = ln - 6
       if(nj(tile) == ni(lw) ) then
          do i = 1, whalo 
             data(1-i, 1:nj(tile)+joff, :)     = sign1*data2_all(ni(lw)+joff:1:-1, nj(lw)-i+1, :, lw) ! west 
          end do
       end if
       if(nj(tile) == nj(le) ) then
          data(ni(tile)+1+ioff:ni(tile)+ehalo+ioff, 1:nj(tile)+joff, :) = data1_all(1+ioff:ehalo+ioff, 1:nj(le)+joff, :, le) ! east 
       end if
       if(ni(tile) == ni(ls) ) then
          data(1:ni(tile)+ioff, 1-shalo:0, :)     = data1_all(1:ni(ls)+ioff, nj(ls)-shalo+1:nj(ls), :, ls) ! south 
       end if
       if(ni(tile) == nj(ln) ) then
          do i = 1, nhalo 
             data(1:ni(tile)+ioff, nj(tile)+i+joff, :)    = sign2*data2_all(i+joff, nj(ln)+ioff:1:-1, :, ln) ! north 
          end do
       end if
    end if

  end subroutine fill_cubicgrid_refined_halo
    

  !##################################################################################
  subroutine test_halo_update( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, x1, x2, x3, x4
    real, allocatable, dimension(:,:,:) :: y, y1, y2, y3, y4
    type(domain2D) :: domain
    real,    allocatable :: global1(:,:,:), global2(:,:,:), global(:,:,:)
    logical, allocatable :: maskmap(:,:)
    integer              :: shift, i, xhalo, yhalo
    logical              :: is_symmetry, folded_south, folded_west, folded_east
    integer              :: is, ie, js, je, isd, ied, jsd, jed

    ! when testing maskmap option, nx*ny should be able to be divided by both npes and npes+1
    if(type == 'Masked' .or. type == 'Masked symmetry') then
       if(mod(nx*ny, npes) .NE. 0 .OR. mod(nx*ny, npes+1) .NE. 0 ) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: nx*ny can not be divided by both npes and npes+1, '//&
               'Masked test_halo_update will not be tested')
          return
       end if
    end if

    if(type == 'Folded xy_halo' ) then
       xhalo = max(whalo, ehalo); yhalo = max(shalo, nhalo)
       allocate(global(1-xhalo:nx+xhalo,1-yhalo:ny+yhalo,nz) )
    else
       allocate(global(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )
    end if

    global = 0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             global(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    if(index(type, 'symmetry') == 0) then
       is_symmetry = .false.
    else
       is_symmetry = .true.
    end if
    select case(type)
    case( 'Simple', 'Simple symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, name=type, symmetry = is_symmetry )
    case( 'Cyclic', 'Cyclic symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,        &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=type, symmetry = is_symmetry )
        global(1-whalo:0,                 1:ny,:) = global(nx-whalo+1:nx,             1:ny,:)
        global(nx+1:nx+ehalo,             1:ny,:) = global(1:ehalo,                   1:ny,:)
        global(1-whalo:nx+ehalo,     1-shalo:0,:) = global(1-whalo:nx+ehalo, ny-shalo+1:ny,:)
        global(1-whalo:nx+ehalo, ny+1:ny+nhalo,:) = global(1-whalo:nx+ehalo,       1:nhalo,:)
    case( 'Folded-north', 'Folded-north symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_north_halo(global, 0, 0, 0, 0, 1)
    case( 'Folded-south symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_SOUTH_EDGE, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_south_halo(global, 0, 0, 0, 0, 1)
    case( 'Folded-west symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=FOLD_WEST_EDGE, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_west_halo(global, 0, 0, 0, 0, 1)
    case( 'Folded-east symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=FOLD_EAST_EDGE, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_east_halo(global, 0, 0, 0, 0, 1)
    case( 'Folded xy_halo' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=xhalo, yhalo=yhalo,   &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name=type, symmetry = is_symmetry  )
        global(1-xhalo:0,                1:ny,:) = global(nx-xhalo+1:nx,                   1:ny,:)
        global(nx+1:nx+xhalo,            1:ny,:) = global(1:xhalo,                         1:ny,:)
        global(1-xhalo:nx+xhalo,ny+1:ny+yhalo,:) = global(nx+xhalo:1-xhalo:-1, ny:ny-yhalo+1:-1,:)
    case( 'Masked', 'Masked symmetry' )
!with fold and cyclic, assign to npes+1 and mask out the top-rightdomain
        call mpp_define_layout( (/1,nx,1,ny/), npes+1, layout )
        allocate( maskmap(layout(1),layout(2)) )
        maskmap(:,:) = .TRUE.; maskmap(layout(1),layout(2)) = .FALSE.
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
             maskmap=maskmap, name=type, symmetry = is_symmetry  )
        deallocate(maskmap)
       !we need to zero out the global data on the missing domain.
       !this logic assumes top-right, in an even division
        if( mod(nx,layout(1)).NE.0 .OR. mod(ny,layout(2)).NE.0 )call mpp_error( FATAL, &
             'TEST_MPP_DOMAINS: test for masked domains needs (nx,ny) to divide evenly on npes+1 PEs.' )
        global(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny,:) = 0
        call fill_folded_north_halo(global, 0, 0, 0, 0, 1)        
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select
        
!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( x1(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    x = 0.
    x (is:ie,js:je,:) = global(is:ie,js:je,:)
    x1 = x; x2 = x; x3 = x; x4 = x

!full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain )
    call mpp_clock_end  (id)
    call compare_checksums( x, global(isd:ied,jsd:jed,:), type )

!partial update
    id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x1, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)
    call compare_checksums( x1(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial x1' )
    call compare_checksums( x2(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial x2' )
    call compare_checksums( x3(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial x3' )
    call compare_checksums( x4(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial x4' )
    
    !--- test vector update for FOLDED and MASKED case.
    if(type == 'Simple' .or. type == 'Simple symmetry' .or. type == 'Cyclic' .or. type == 'Cyclic symmetry') then
       deallocate(x,x1,x2,x3,x4)
       return       
    end if

    !------------------------------------------------------------------
    !              vector update : BGRID_NE
    !------------------------------------------------------------------
    shift = 0
    if(is_symmetry) then
       shift = 1
       deallocate(global)
       allocate(global(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
       global = 0.0
       do k = 1,nz
          do j = 1,ny+1
             do i = 1,nx+1
                global(i,j,k) = k + i*1e-3 + j*1e-6
             end do
          end do
       end do
       if(type == 'Masked symmetry') then
           global(nx-nx/layout(1)+1:nx+1,ny-ny/layout(2)+1:ny+1,:) = 0
       endif
       deallocate(x, x1, x2, x3, x4)
       allocate( x (isd:ied+1,jsd:jed+1,nz) )
       allocate( x1(isd:ied+1,jsd:jed+1,nz) )
       allocate( x2(isd:ied+1,jsd:jed+1,nz) )
       allocate( x3(isd:ied+1,jsd:jed+1,nz) )
       allocate( x4(isd:ied+1,jsd:jed+1,nz) )
    endif

    folded_south = .false.
    folded_west  = .false.
    folded_east  = .false.
    select case (type)
    case ('Folded-north', 'Masked')
       !fill in folded north edge, cyclic east and west edge
       call fill_folded_north_halo(global, 1, 1, 0, 0, -1)     
    case ('Folded xy_halo')
       !fill in folded north edge, cyclic east and west edge
       global(1-xhalo:0,                  1:ny,:) =  global(nx-xhalo+1:nx,                     1:ny,:)
       global(nx+1:nx+xhalo,              1:ny,:) =  global(1:xhalo,                           1:ny,:)
       global(1-xhalo:nx+xhalo-1,ny+1:ny+yhalo,:) = -global(nx+xhalo-1:1-xhalo:-1,ny-1:ny-yhalo:-1,:)
       global(nx+xhalo,          ny+1:ny+yhalo,:) = -global(nx-xhalo,             ny-1:ny-yhalo:-1,:)
    case ('Folded-north symmetry', 'Masked symmetry' )
       call fill_folded_north_halo(global, 1, 1, 1, 1, -1)
    case ('Folded-south symmetry' )
       folded_south = .true.
       call fill_folded_south_halo(global, 1, 1, 1, 1, -1)
    case ('Folded-west symmetry' )
       folded_west = .true.
       call fill_folded_west_halo(global, 1, 1, 1, 1, -1)
    case ('Folded-east symmetry' )
       folded_east = .true.
       call fill_folded_east_halo(global, 1, 1, 1, 1, -1)
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

    x = 0.
    x(is:ie+shift,js:je+shift,:) = global(is:ie+shift,js:je+shift,:)
    !set up y array
    allocate( y (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y1(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y2(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y3(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y4(isd:ied+shift,jsd:jed+shift,nz) )
    y = x; x1 = x; x2 = x; x3 = x; x4 = x
    y = x; y1 = x; y2 = x; y3 = x; y4 = x

    id = mpp_clock_id( type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=BGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x3, y3, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=BGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    !redundant points must be equal and opposite

    if(folded_south) then
       global(nx/2+shift,                1,:) = 0.  !pole points must have 0 velocity
       global(nx+shift  ,                1,:) = 0.  !pole points must have 0 velocity
       global(nx/2+1+shift:nx-1+shift,   1,:) = -global(nx/2-1+shift:1+shift:-1, 1,:)
       global(1-whalo:shift,             1,:) = -global(nx-whalo+1:nx+shift,     1,:)
       global(nx+1+shift:nx+ehalo+shift, 1,:) = -global(1+shift:ehalo+shift,     1,:)
       !--- the following will fix the +0/-0 problem on altix
       if(shalo >0) global(shift,1,:) = 0.  !pole points must have 0 velocity
    else if(folded_west) then
       global(1, ny/2+shift, :) = 0. !pole points must have 0 velocity
       global(1, ny+shift,   :) = 0. !pole points must have 0 velocity
       global(1, ny/2+1+shift:ny-1+shift,   :) = -global(1, ny/2-1+shift:1+shift:-1, :)
       global(1, 1-shalo:shift,             :) = -global(1, ny-shalo+1:ny+shift,     :)
       global(1, ny+1+shift:ny+nhalo+shift, :) = -global(1, 1+shift:nhalo+shift,     :)
    else if(folded_east) then
       global(nx+shift, ny/2+shift, :) = 0. !pole points must have 0 velocity
       global(nx+shift, ny+shift,   :) = 0. !pole points must have 0 velocity
       global(nx+shift, ny/2+1+shift:ny-1+shift,   :) = -global(nx+shift, ny/2-1+shift:1+shift:-1, :)
       global(nx+shift, 1-shalo:shift,             :) = -global(nx+shift, ny-shalo+1:ny+shift,     :)
       global(nx+shift, ny+1+shift:ny+nhalo+shift, :) = -global(nx+shift, 1+shift:nhalo+shift,     :)
    else
       global(nx/2+shift,                ny+shift,:) = 0.  !pole points must have 0 velocity
       global(nx+shift  ,                ny+shift,:) = 0.  !pole points must have 0 velocity
       global(nx/2+1+shift:nx-1+shift,   ny+shift,:) = -global(nx/2-1+shift:1+shift:-1, ny+shift,:)
       if(type == 'Folded xy_halo') then
          global(1-xhalo:shift,             ny+shift,:) = -global(nx-xhalo+1:nx+shift,     ny+shift,:)
          global(nx+1+shift:nx+xhalo+shift, ny+shift,:) = -global(1+shift:xhalo+shift,     ny+shift,:)
       else
          global(1-whalo:shift,             ny+shift,:) = -global(nx-whalo+1:nx+shift,     ny+shift,:)
          global(nx+1+shift:nx+ehalo+shift, ny+shift,:) = -global(1+shift:ehalo+shift,     ny+shift,:)
       end if
       !--- the following will fix the +0/-0 problem on altix
       if(nhalo >0) global(shift,ny+shift,:) = 0.  !pole points must have 0 velocity
    endif

    call compare_checksums( x,  global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X' )
    call compare_checksums( y,  global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y' )
    call compare_checksums( x1, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X1' )
    call compare_checksums( x2, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X2' )
    call compare_checksums( x3, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X3' )
    call compare_checksums( x4, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X4' )
    call compare_checksums( y1, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y1' )
    call compare_checksums( y2, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y2' )
    call compare_checksums( y3, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y3' )
    call compare_checksums( y4, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y4' )

    deallocate(global, x, x1, x2, x3, x4, y, y1, y2, y3, y4)

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- global1 is x-component and global2 is y-component
    if(type == 'Folded xy_halo') then
       allocate(global1(1-xhalo:nx+xhalo, 1-yhalo:ny+yhalo, nz))
       allocate(global2(1-xhalo:nx+xhalo, 1-yhalo:ny+yhalo, nz))
    else
       allocate(global1(1-whalo:nx+ehalo+shift, 1-shalo:ny+nhalo, nz))
       allocate(global2(1-whalo:nx+ehalo, 1-shalo:ny+nhalo+shift, nz))
    end if
    allocate(x (isd:ied+shift,jsd:jed,nz), y (isd:ied,jsd:jed+shift,nz) )
    allocate(x1(isd:ied+shift,jsd:jed,nz), y1(isd:ied,jsd:jed+shift,nz) )
    allocate(x2(isd:ied+shift,jsd:jed,nz), y2(isd:ied,jsd:jed+shift,nz) )
    allocate(x3(isd:ied+shift,jsd:jed,nz), y3(isd:ied,jsd:jed+shift,nz) )
    allocate(x4(isd:ied+shift,jsd:jed,nz), y4(isd:ied,jsd:jed+shift,nz) )
    
    global1 = 0.0
    global2 = 0.0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx+shift
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
       do j = 1,ny+shift
          do i = 1,nx
             global2(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    if(type == 'Masked' .or. type == 'Masked symmetry') then
       global1(nx-nx/layout(1)+1:nx+shift,ny-ny/layout(2)+1:ny,:) = 0
       global2(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny+shift,:) = 0
    end if

    select case (type)
    case ('Folded-north', 'Masked')
       !fill in folded north edge, cyclic east and west edge
       call fill_folded_north_halo(global1, 1, 0, 0, 0, -1)
       call fill_folded_north_halo(global2, 0, 1, 0, 0, -1)
    case ('Folded xy_halo')
       global1(1-xhalo:0,                   1:ny,:) =  global1(nx-xhalo+1:nx,                     1:ny,:)
       global1(nx+1:nx+xhalo,               1:ny,:) =  global1(1:xhalo,                           1:ny,:)
       global2(1-xhalo:0,                   1:ny,:) =  global2(nx-xhalo+1:nx,                     1:ny,:)
       global2(nx+1:nx+xhalo,               1:ny,:) =  global2(1:xhalo,                           1:ny,:)
       global1(1-xhalo:nx+xhalo-1, ny+1:ny+yhalo,:) = -global1(nx+xhalo-1:1-xhalo:-1, ny:ny-yhalo+1:-1,:)
       global1(nx+xhalo,           ny+1:ny+yhalo,:) = -global1(nx-xhalo,              ny:ny-yhalo+1:-1,:)
       global2(1-xhalo:nx+xhalo,   ny+1:ny+yhalo,:) = -global2(nx+xhalo:1-xhalo:-1,   ny-1:ny-yhalo:-1,:)
    case ('Folded-north symmetry')
       call fill_folded_north_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_north_halo(global2, 0, 1, 0, 1, -1)
    case ('Folded-south symmetry')
       call fill_folded_south_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_south_halo(global2, 0, 1, 0, 1, -1)
    case ('Folded-west symmetry')
       call fill_folded_west_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_west_halo(global2, 0, 1, 0, 1, -1)
    case ('Folded-east symmetry')
       call fill_folded_east_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_east_halo(global2, 0, 1, 0, 1, -1)
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

    x = 0.; y = 0.
    x(is:ie+shift,js:je,      :) = global1(is:ie+shift,js:je,      :)
    y(is:ie      ,js:je+shift,:) = global2(is:ie,      js:je+shift,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x3, y3, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=CGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    !redundant points must be equal and opposite
    if(folded_south) then
       global2(nx/2+1:nx,     1,:) = -global2(nx/2:1:-1, 1,:)
       global2(1-whalo:0,     1,:) = -global2(nx-whalo+1:nx, 1, :)
       global2(nx+1:nx+ehalo, 1,:) = -global2(1:ehalo,       1, :)
    else if(folded_west) then
       global1(1, ny/2+1:ny,     :) = -global1(1, ny/2:1:-1,     :)
       global1(1, 1-shalo:0,     :) = -global1(1, ny-shalo+1:ny, :)
       global1(1, ny+1:ny+nhalo, :) = -global1(1, 1:nhalo,       :)
    else if(folded_east) then
       global1(nx+shift, ny/2+1:ny,     :) = -global1(nx+shift, ny/2:1:-1,     :)
       global1(nx+shift, 1-shalo:0,     :) = -global1(nx+shift, ny-shalo+1:ny, :)
       global1(nx+shift, ny+1:ny+nhalo, :) = -global1(nx+shift, 1:nhalo,       :)
    else
       global2(nx/2+1:nx,     ny+shift,:) = -global2(nx/2:1:-1, ny+shift,:)
       if(type == 'Folded xy_halo') then
          global2(1-xhalo:0,     ny+shift,:) = -global2(nx-xhalo+1:nx, ny+shift,:)
          global2(nx+1:nx+xhalo, ny+shift,:) = -global2(1:xhalo,       ny+shift,:)
       else
          global2(1-whalo:0,     ny+shift,:) = -global2(nx-whalo+1:nx, ny+shift,:)
          global2(nx+1:nx+ehalo, ny+shift,:) = -global2(1:ehalo,       ny+shift,:)
       end if
    endif

    call compare_checksums( x,  global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y' )
    call compare_checksums( x1, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X1' )
    call compare_checksums( x2, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X2' )
    call compare_checksums( x3, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X3' )
    call compare_checksums( x4, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y4' )

    deallocate(global1, global2, x, x1, x2, x3, x4, y, y1, y2, y3, y4)


  end subroutine test_halo_update

  !##################################################################################
  subroutine test_cyclic_offset( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, x1, x2, x3, x4
    real, allocatable, dimension(:,:,:) :: y, y1, y2, y3, y4
    type(domain2D) :: domain
    real,    allocatable :: global1(:,:,:), global2(:,:,:), global(:,:,:)
    integer              :: i, j, k, jj, ii
    integer              :: is, ie, js, je, isd, ied, jsd, jed
    character(len=128)   :: type2

    allocate(global(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz))

    global = 0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             global(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'x_cyclic_offset' )
        write(type2, *)type, ' x_cyclic=', x_cyclic_offset
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 name=type, x_cyclic_offset = x_cyclic_offset)
        do j = 1, ny
           jj = mod(j + x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(1-whalo:0,j,:) = global(nx-whalo+1:nx, jj,:) ! West
           jj = mod(j - x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(nx+1:nx+ehalo,j,:) = global(1:ehalo,jj,:)    ! East
        end do
    case( 'y_cyclic_offset' )
        write(type2, *)type, ' y_cyclic = ', y_cyclic_offset
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, yflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 name=type, y_cyclic_offset = y_cyclic_offset)
        do i = 1, nx
           ii = mod(i + y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(i, 1-shalo:0,:) = global(ii, ny-shalo+1:ny,:) ! South
           ii = mod(i - y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(i,ny+1:ny+nhalo,:) = global(ii,1:nhalo,:)    ! NORTH
        end do
    case( 'torus_x_offset' )
        write(type2, *)type, ' x_cyclic = ', x_cyclic_offset
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 yflags=CYCLIC_GLOBAL_DOMAIN, name=type,                  &
                                 x_cyclic_offset = x_cyclic_offset)
        do j = 1, ny
           jj = mod(j + x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(1-whalo:0,j,:) = global(nx-whalo+1:nx, jj,:) ! West
           jj = mod(j - x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(nx+1:nx+ehalo,j,:) = global(1:ehalo,jj,:)    ! East
        end do
        global(1:nx,1-shalo:0,:)     = global(1:nx, ny-shalo+1:ny,:) ! South
        global(1:nx,ny+1:ny+nhalo,:) = global(1:nx, 1:nhalo, :)    ! NORTH
        
        do j = 1, shalo
           jj = mod(ny-j+1 + x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(1-whalo:0, 1-j,:) = global(nx-whalo+1:nx, jj, :)  ! Southwest
           jj = mod(ny-j+1-x_cyclic_offset+ny,ny)
           if(jj==0) jj = ny
           global(nx+1:nx+ehalo, 1-j,:) = global(1:ehalo, jj, :)    ! Southeast
        end do
        do j = 1, nhalo
           jj = mod(j + x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(1-whalo:0, ny+j,:) = global(nx-whalo+1:nx, jj, :)  ! northwest
           jj = mod(j - x_cyclic_offset+ny,ny)
           if(jj==0) jj = ny
           global(nx+1:nx+ehalo, ny+j,:) = global(1:ehalo, jj, :)    ! northeast
        end do

    case( 'torus_y_offset' )
        write(type2, *)type, ' y_cyclic = ', y_cyclic_offset
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 yflags=CYCLIC_GLOBAL_DOMAIN, name=type,                  &
                                 y_cyclic_offset = y_cyclic_offset)
        do i = 1, nx
           ii = mod(i + y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(i, 1-shalo:0,:) = global(ii, ny-shalo+1:ny,:) ! South
           ii = mod(i - y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(i,ny+1:ny+nhalo,:) = global(ii,1:nhalo,:)    ! NORTH
        end do
        global(1-whalo:0,1:ny,:)     = global(nx-whalo+1:nx, 1:ny,:) ! West
        global(nx+1:nx+ehalo,1:ny,:) = global(1:ehalo, 1:ny, :)      ! East
        do i = 1, whalo
           ii = mod(nx-i+1 + y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(1-i, 1-shalo:0,:) = global(ii, ny-shalo+1:ny,:) ! southwest
           ii = mod(nx-i+1 - y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(1-i,ny+1:ny+nhalo,:) = global(ii,1:nhalo,:)    ! northwest
        end do
        do i = 1, ehalo
           ii = mod(i + y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(nx+i, 1-shalo:0,:) = global(ii, ny-shalo+1:ny,:) ! southeast
           ii = mod(i - y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(nx+i,ny+1:ny+nhalo,:) = global(ii,1:nhalo,:)    ! northeast
        end do
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select
        
!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( x1(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    x = 0.
    x (is:ie,js:je,:) = global(is:ie,js:je,:)
    x1 = x; x2 = x; x3 = x; x4 = x

!full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain )
    call mpp_clock_end  (id)
    call compare_checksums( x, global(isd:ied,jsd:jed,:), trim(type2) )

!partial update
    id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x1, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)
    call compare_checksums( x1(is:ied,js:jed,:), global(is:ied,js:jed,:), trim(type2)//' partial x1' )
    call compare_checksums( x2(is:ied,js:jed,:), global(is:ied,js:jed,:), trim(type2)//' partial x2' )
    call compare_checksums( x3(is:ied,js:jed,:), global(is:ied,js:jed,:), trim(type2)//' partial x3' )
    call compare_checksums( x4(is:ied,js:jed,:), global(is:ied,js:jed,:), trim(type2)//' partial x4' )
    
    !--- test vector update for FOLDED and MASKED case.
    deallocate(x,x1,x2,x3,x4)


    !------------------------------------------------------------------
    !              vector update : BGRID_NE
    !------------------------------------------------------------------
    !--- global1 is x-component and global2 is y-component
    allocate(global1(1-whalo:nx+ehalo, 1-shalo:ny+nhalo, nz))
    allocate(global2(1-whalo:nx+ehalo, 1-shalo:ny+nhalo, nz))
    allocate(x (isd:ied,jsd:jed,nz), y (isd:ied,jsd:jed,nz) )
    allocate(x1(isd:ied,jsd:jed,nz), y1(isd:ied,jsd:jed,nz) )
    allocate(x2(isd:ied,jsd:jed,nz), y2(isd:ied,jsd:jed,nz) )
    allocate(x3(isd:ied,jsd:jed,nz), y3(isd:ied,jsd:jed,nz) )
    allocate(x4(isd:ied,jsd:jed,nz), y4(isd:ied,jsd:jed,nz) )
    where (global >0)
       global1 = 1000 + global
       global2 = 2000 + global
    elsewhere
       global1 = 0
       global2 = 0
    end where
    x = 0.; y = 0
    x(is:ie,js:je,:) = global1(is:ie,js:je,:)
    y(is:ie,js:je,:) = global2(is:ie,js:je,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    id = mpp_clock_id( type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=BGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x3, y3, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=BGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    !redundant points must be equal and opposite

    call compare_checksums( x,  global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y' )
    call compare_checksums( x1, global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X1' )
    call compare_checksums( x2, global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X2' )
    call compare_checksums( x3, global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X3' )
    call compare_checksums( x4, global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y4' )

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------

    x = 0.; y = 0.
    x(is:ie,js:je,:) = global1(is:ie,js:je,:)
    y(is:ie,js:je,:) = global2(is:ie,js:je,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x3, y3, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=CGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    call compare_checksums( x,  global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y' )
    call compare_checksums( x1, global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X1' )
    call compare_checksums( x2, global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X2' )
    call compare_checksums( x3, global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X3' )
    call compare_checksums( x4, global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y4' )

    deallocate(global1, global2, x, x1, x2, x3, x4, y, y1, y2, y3, y4)


  end subroutine test_cyclic_offset


  subroutine test_global_field( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, gcheck
    type(domain2D) :: domain
    real, allocatable    :: global1(:,:,:)
    integer              :: ishift, jshift, ni, nj, i, j, position
    integer, allocatable :: pelist(:)
    integer              :: is, ie, js, je, isd, ied, jsd, jed

    !--- set up domain    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
        
    !--- determine if an extra point is needed
    ishift = 0; jshift = 0
    position = CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1; jshift = 1; position=CORNER
    case ('Symmetry east')
       ishift = 1; jshift = 0; position=EAST
    case ('Symmetry north')
       ishift = 0; jshift = 1; position=NORTH
    end select

    ie  = ie+ishift;  je  = je+jshift
    ied = ied+ishift; jed = jed+jshift
    ni  = nx+ishift;  nj  = ny+jshift
    allocate(global1(1-whalo:ni+ehalo, 1-shalo:nj+nhalo, nz))
    global1 = 0.0
    do k = 1,nz
       do j = 1,nj
          do i = 1,ni
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    enddo

    allocate( gcheck(ni, nj, nz) )
    allocate( x (isd:ied,jsd:jed,nz) )

    x(:,:,:) = global1(isd:ied,jsd:jed,:)

    !--- test the data on data domain
    gcheck = 0.    
    id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on data domain' )

    !--- Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !--- will be declared. But for the x-direction global field, mpp_sync_self will
    !--- be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !--- in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !--- deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !--- will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !--- on all pe is needed for those partial pelist. But for y-update, it is ok. 
    !--- because the pelist in y-update is not continous.
    allocate(pelist(0:layout(1)-1))    
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !xupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags = XUPDATE, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
                            type//' mpp_global_field xupdate only on data domain' )

    !yupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags = YUPDATE, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
                            type//' mpp_global_field yupdate only on data domain' )

    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )

    call mpp_clock_end  (id)                                          
    !compare checksums between global and x arrays  
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, &
                            type//' mpp_global_field on data domain' )

    !--- test the data on compute domain
    gcheck = 0.    
    id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je, :), gcheck, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on compute domain' )

    !xupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, flags = XUPDATE, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
                            type//' mpp_global_field xupdate only on compute domain' )

    !yupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, flags = YUPDATE, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
                            type//' mpp_global_field yupdate only on compute domain' )


    deallocate(global1, gcheck, x)

  end subroutine test_global_field

    !--- test mpp_global_sum, mpp_global_min and mpp_global_max
  subroutine test_global_reduce (type)
    character(len=*), intent(in) :: type
    real    :: lsum, gsum, lmax, gmax, lmin, gmin
    integer :: ni, nj, ishift, jshift, position
    integer              :: is, ie, js, je, isd, ied, jsd, jed

    type(domain2D) :: domain
    real, allocatable, dimension(:,:,:) :: global1, x
    real, allocatable, dimension(:,:)   :: global2D
    !--- set up domain    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Simple' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type )
    case( 'Simple symmetry center', 'Simple symmetry corner', 'Simple symmetry east', 'Simple symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case( 'Cyclic symmetry center', 'Cyclic symmetry corner', 'Cyclic symmetry east', 'Cyclic symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                                    name=type, symmetry = .true., xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
        
    !--- determine if an extra point is needed
    ishift = 0; jshift = 0; position = CENTER
    select case(type)
    case ('Simple symmetry corner', 'Cyclic symmetry corner')
       ishift = 1; jshift = 1; position = CORNER
    case ('Simple symmetry east', 'Cyclic symmetry east' )
       ishift = 1; jshift = 0; position = EAST
    case ('Simple symmetry north', 'Cyclic symmetry north')
       ishift = 0; jshift = 1; position = NORTH
    end select

    ie  = ie+ishift;  je  = je+jshift
    ied = ied+ishift; jed = jed+jshift
    ni  = nx+ishift;  nj  = ny+jshift
    allocate(global1(1-whalo:ni+ehalo, 1-shalo:nj+nhalo, nz))
    global1 = 0.0
    do k = 1,nz
       do j = 1,nj
          do i = 1,ni
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    enddo

    !--- NOTE: even though the domain is cyclic, no need to apply cyclic condition on the global data

    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( global2D(ni,nj))

    x(:,:,:) = global1(isd:ied,jsd:jed,:)
    do j = 1, nj
       do i = 1, ni
          global2D(i,j) = sum(global1(i,j,:))
       enddo 
    enddo
    !test mpp_global_sum
   
    if(type(1:6) == 'Simple') then
       gsum = sum( global2D(1:ni,1:nj) )
    else
       gsum = sum( global2D(1:nx, 1:ny) )
    endif
    id = mpp_clock_id( type//' sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lsum = mpp_global_sum( domain, x, position = position  )
    call mpp_clock_end  (id)
    if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', type//' Fast sum=', lsum, gsum

    !test exact mpp_global_sum
    id = mpp_clock_id( type//' exact sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lsum = mpp_global_sum( domain, x, BITWISE_EXACT_SUM, position = position )
    call mpp_clock_end  (id)
    !--- The following check will fail on altix in normal mode, but it is ok
    !--- in debugging mode. It is ok on irix.
    call compare_data_scalar(lsum, gsum, FATAL, type//' mpp_global_exact_sum')

    !test mpp_global_min
    gmin = minval(global1(1:ni, 1:nj, :))
    id = mpp_clock_id( type//' min', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lmin = mpp_global_min( domain, x, position = position )
    call mpp_clock_end  (id)
    call compare_data_scalar(lmin, gmin, FATAL, type//' mpp_global_min')

    !test mpp_global_max
    gmax = maxval(global1(1:ni, 1:nj, :))
    id = mpp_clock_id( type//' max', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lmax = mpp_global_max( domain, x, position = position )
    call mpp_clock_end  (id)
    call compare_data_scalar(lmax, gmax, FATAL, type//' mpp_global_max' )

    deallocate(global1, x)

  end subroutine test_global_reduce


  subroutine test_parallel ( )
  
    integer :: npes, layout(2), i, j, k,is, ie, js, je, isd, ied, jsd, jed
    real, dimension(:,:), allocatable :: field, lfield
    real, dimension(:,:,:), allocatable :: field3d, lfield3d
    type(domain2d) :: domain
    integer, dimension(:), allocatable :: pelist1 , pelist2
    logical :: group1, group2
    character(len=128)  :: mesg
    
    npes = mpp_npes()
    allocate(pelist1(npes-mpes), pelist2(mpes))
    pelist1 = (/(i, i = 0, npes-mpes -1)/)
    pelist2 = (/(i, i = npes-mpes, npes - 1)/)
    call mpp_declare_pelist(pelist1)
    call mpp_declare_pelist(pelist2)
    group1 = .FALSE. ; group2 = .FALSE.
    if(any(pelist1==pe)) group1 = .TRUE.
    if(any(pelist2==pe)) group2 = .TRUE.
    mesg = 'parallel checking'
    
    if(group1) then
       call mpp_set_current_pelist(pelist1)
       call mpp_define_layout( (/1,nx,1,ny/), npes-mpes, layout )
    else if(group2) then
       call mpp_set_current_pelist(pelist2)
       call mpp_define_layout( (/1,nx,1,ny/), mpes, layout )
    endif
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)

    call mpp_set_current_pelist() 
     
     call mpp_get_compute_domain(domain, is, ie, js, je)
     call mpp_get_data_domain(domain, isd, ied, jsd, jed)
     allocate(lfield(is:ie,js:je),field(isd:ied,jsd:jed))
     allocate(lfield3d(is:ie,js:je,nz),field3d(isd:ied,jsd:jed,nz))
     
     do i = is, ie
     do j = js, je
        lfield(i,j) = real(i)+real(j)*0.001
     enddo
     enddo
     do i = is, ie
     do j = js, je
     do k = 1, nz
        lfield3d(i,j,k) = real(i)+real(j)*0.001+real(k)*0.00001
     enddo
     enddo
     enddo
     field = 0.0
     field3d = 0.0
     field(is:ie,js:je)= lfield(is:ie,js:je)
     field3d(is:ie,js:je,:) = lfield3d(is:ie,js:je,:)
     call mpp_update_domains(field,domain)
     call mpp_update_domains(field3d,domain)
     
    call mpp_check_field(field, pelist1, pelist2,domain, '2D '//mesg, w_halo = whalo, &
                            s_halo = shalo, e_halo = ehalo, n_halo = nhalo)
    call mpp_check_field(field3d, pelist1, pelist2,domain, '3D '//mesg, w_halo = whalo, &
                            s_halo = shalo, e_halo = ehalo, n_halo = nhalo)
                            
  end subroutine test_parallel
  
  subroutine test_modify_domain( )
  
    type(domain2D) :: domain2d_no_halo, domain2d_with_halo
    integer :: is1, ie1, js1, je1, isd1, ied1, jsd1, jed1
    integer :: is2, ie2, js2, je2, isd2, ied2, jsd2, jed2
    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain2d_no_halo,   &
                            yflags=CYCLIC_GLOBAL_DOMAIN, xhalo=0, yhalo=0)

    call mpp_get_compute_domain(domain2d_no_halo, is1, ie1, js1, je1)
    call mpp_get_data_domain(domain2d_no_halo, isd1, ied1, jsd1, jed1)
    call mpp_modify_domain(domain2d_no_halo, domain2d_with_halo, whalo=whalo,ehalo=ehalo,shalo=shalo,nhalo=nhalo)
    call mpp_get_compute_domain(domain2d_with_halo, is2, ie2, js2, je2)
    call mpp_get_data_domain(domain2d_with_halo, isd2, ied2, jsd2, jed2)
    if( is1 .NE. is2 .OR. ie1 .NE. ie2 .OR. js1 .NE. js2 .OR. je1 .NE. je2 ) then
        print*, "at pe ", pe, " compute domain without halo: ", is1, ie1, js1, je1, &
                " is not equal to the domain with halo ", is2, ie2, js2, je2        
        call mpp_error(FATAL, "compute domain mismatch between domain without halo and domain with halo")
    end if

    if( isd1-whalo .NE. isd2 .OR. ied1+ehalo .NE. ied2 .OR. jsd1-shalo .NE. jsd2 .OR. jed1+nhalo .NE. jed2 ) then
        print*, "at pe ", pe, "halo is w=",whalo,",e=",ehalo,",s=",shalo,"n=",nhalo, &
               ",data domain without halo is ",isd1, ied1, jsd1, jed1,                     &
               ", data domain with halo is ", isd2, ied2, jsd2, jed2 
    else
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, 'test_modify_domain: OK.' )
    end if

    return
    
end subroutine test_modify_domain

  subroutine compare_checksums( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(LONG_KIND) :: sum1, sum2
    integer :: i, j, k

    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
    ! mpp_sync()
    call mpp_sync_self()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
         call mpp_error(FATAL,'compare_chksum: size of a and b does not match')

    do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,k) .ne. b(i,j,k)) then
                write(stdunit,'(a,i3,a,i3,a,i3,a,i3,a,f16.9,a,f16.9)')" at pe ", mpp_pe(), &
                     ", at point (",i,", ", j, ", ", k, "), a = ", a(i,j,k), ", b = ", b(i,j,k)
                call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
             endif
          enddo
       enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
        !--- in some case, even though checksum agree, the two arrays 
        !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
        !--- hence we need to check the value point by point.
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums

  !###########################################################################

  subroutine compare_data_scalar( a, b, action, string )
    real,             intent(in) :: a, b
    integer,          intent(in) :: action
    character(len=*), intent(in) :: string
    if( a .EQ. b)then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': data comparison are OK.' )
    else
        write(stdunit,'(a,i3,a,es12.4,a,es12.4,a,es12.4)')' on pe ', mpp_pe(),' a = ', a, ', b = ', b, ', a - b =', a-b
        call mpp_error( action, trim(string)//': data comparison are not OK.' )
    end if

  end subroutine compare_data_scalar

  subroutine test_get_neighbor_1d
    type(domain1d) :: dmn1d
    integer npes, peN, peS
    npes = mpp_npes()
    call mpp_define_domains((/1,npes/), npes, dmn1d)
    call mpp_get_neighbor_pe(dmn1d, direction=+1, pe=peN)
    call mpp_get_neighbor_pe(dmn1d, direction=-1, pe=peS)
    print '(a,i2,a,2i3)', 'PE: ', mpp_pe(), ' R/L pes: ', peN, peS
  end subroutine test_get_neighbor_1d

  subroutine test_get_neighbor_non_cyclic
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then 
       call mpp_error(NOTE, 'test_mpp_domains: test_get_neighbor_non_cyclic '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (non-cyclic): ', layout,  &
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_non_cyclic

  subroutine test_get_neighbor_cyclic
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then 
       call mpp_error(NOTE, 'test_mpp_domains: test_get_neighbor_cyclic '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (cyclic)    : ', layout, & 
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_cyclic

  subroutine test_get_neighbor_folded_north
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then 
       call mpp_error(NOTE, 'test_mpp_domains: test_get_neighbor_folded_north '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (folded N)  : ', layout, & 
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_folded_north

  subroutine test_get_neighbor_mask
    logical, allocatable ::  mask(:,:)
    integer :: im, jm, n_remove
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    
    n_remove = 2
    if( npes .NE. 8 ) then 
       call mpp_error(NOTE, 'test_mpp_domains: test_get_neighbor_mask '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes+n_remove, layout )
    allocate(mask(layout(1), layout(2)))
    mask = .TRUE.  ! activate domains
    im = min(layout(1), ceiling(layout(1)/2.0))
    jm = min(layout(2), ceiling(layout(2)/2.0))
    mask(im  ,jm  ) = .FALSE. ! deactivate domain
    mask(im  ,jm-1) = .FALSE. ! deactivate domain
    print '(a,2i3,a,2i3)', 'Masked out domains ', im, jm, ' and ', im,jm-1
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         maskmap=mask)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i3,a,2i3,a,8i3)','PE: ', mpp_pe(), ' layout (mask   )  : ', layout, & 
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_mask

  subroutine test_define_mosaic_pelist(type, ntile)
    character(len=*),       intent(in) :: type
    integer,                intent(in) :: ntile
    integer                            :: npes, root_pe, start_pe, n, ntile_per_pe
    integer, dimension(:), allocatable :: pe1_start, pe1_end, pe2_start, pe2_end
    integer, dimension(:), allocatable :: sizes, costpertile

    root_pe = mpp_root_pe()
    npes = mpp_npes()

    allocate(sizes(ntile), pe1_start(ntile), pe1_end(ntile), pe2_start(ntile), pe2_end(ntile),costpertile(ntile) )
    costpertile = 1
    sizes = nx*ny
    if(npes ==1) then
       pe1_start = root_pe; pe1_end = root_pe
    end if
    select case(type)
    case('One tile')
       pe1_start = root_pe; pe1_end = npes+root_pe-1
    case('Two uniform tile')
       if(mod(npes,2) .NE. 0 .AND. npes .NE. 1) then
          call mpp_error(NOTE, 'test_define_mosaic_pelist: npes can not be divided by 2, no test for '//type )
          return
       end if
       if(npes .NE. 1) then
          pe1_start(1) = root_pe;        pe1_end(1) = npes/2+root_pe-1
          pe1_start(2) = npes/2+root_pe; pe1_end(2) = npes+root_pe-1       
       end if
    case('Two nonuniform tile')
       if(mod(npes,3) .NE. 0 .AND. npes .NE. 1) then
          call mpp_error(NOTE, 'test_define_mosaic_pelist: npes can not be divided by 3, no test for '//type )
          return
       end if
       sizes(1) = 2*nx*ny
       if(npes .NE. 1) then
          pe1_start(1) = root_pe;          pe1_end(1) = npes/3*2+root_pe-1
          pe1_start(2) = npes/3*2+root_pe; pe1_end(2) = npes+root_pe-1
       end if
    case('Ten tile')
       if(mod(npes,10) .NE. 0 .AND. npes .NE. 1 .AND. mod(10,npes) .NE. 0) then
          call mpp_error(NOTE, 'test_define_mosaic_pelist: npes can not be divided by 10(or reverse), no test for '//type )
          return
       end if
       if(mod(10, npes)==0) then
          ntile_per_pe = ntile/npes          
          do n = 1, ntile
             pe1_start(n) = root_pe+(n-1)/ntile_per_pe; pe1_end(n) = pe1_start(n)
          end do
       else if(mod(npes,10) == 0) then
          do n = 1, ntile
             pe1_start(n) = npes/10*(n-1)+root_pe; pe1_end(n) = npes/10*n+root_pe-1
          end do
       end if
    case('Ten tile with nonuniform cost')
       if(mod(npes,15) .NE. 0 .AND. npes .NE. 1) then
          call mpp_error(NOTE, 'test_define_mosaic_pelist: npes can not be divided by 15, no test for '//type )
          return
       end if
       costpertile(1:5) = 2; costpertile(6:ntile) = 1
       if(npes .NE. 1) then
          start_pe = root_pe
          do n = 1, ntile
             pe1_start(n) = start_pe
             pe1_end(n)   = start_pe + npes/15*costpertile(n)-1
             start_pe = pe1_end(n) + 1
          end do
       end if
    case default
       call mpp_error(FATAL,"test_define_mosaic_pelist: "//type//" is an invalid type")
    end select

    call mpp_define_mosaic_pelist( sizes, pe2_start, pe2_end, costpertile=costpertile)
    if( ANY(pe1_start .NE. pe2_start) .OR. ANY(pe1_end .NE. pe2_end) ) then
       call mpp_error(FATAL,"test_define_mosaic_pelist: test failed for "//trim(type) )
    else
       call mpp_error(NOTE,"test_define_mosaic_pelist: test successful for "//trim(type) )
    end if

  end subroutine test_define_mosaic_pelist

end program test
#else
module null_mpp_domains_test
end module
#endif
