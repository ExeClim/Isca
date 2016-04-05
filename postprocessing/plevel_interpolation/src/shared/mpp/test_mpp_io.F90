#ifdef test_mpp_io
program test
#include <fms_platform.h>

  use mpp_mod,         only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_sync_self
  use mpp_mod,         only : FATAL, NOTE, mpp_chksum, MPP_DEBUG, mpp_set_stack_size, MPP_CLOCK_SYNC
  use mpp_mod,         only : mpp_sync, mpp_exit, mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_domains_mod, only : mpp_define_domains, mpp_domains_set_stack_size, domain1D, mpp_get_global_domain
  use mpp_domains_mod, only : domain2D, mpp_define_layout, mpp_get_domain_components, mpp_define_mosaic
  use mpp_domains_mod, only : mpp_get_memory_domain, mpp_get_compute_domain, mpp_domains_exit
  use mpp_domains_mod, only : CENTER, EAST, NORTH, CORNER, mpp_get_data_domain
  use mpp_domains_mod, only : mpp_define_io_domain, mpp_deallocate_domain
  use mpp_io_mod,      only : mpp_io_init, mpp_write_meta, axistype, fieldtype, atttype
  use mpp_io_mod,      only : MPP_RDONLY, mpp_open, MPP_OVERWR, MPP_ASCII, MPP_SINGLE
  use mpp_io_mod,      only : MPP_NETCDF, MPP_MULTI, mpp_get_atts, mpp_write, mpp_close
  use mpp_io_mod,      only : mpp_get_info, mpp_get_axes, mpp_get_fields, mpp_get_times
  use mpp_io_mod,      only : mpp_read, mpp_io_exit

  implicit none

#ifdef use_netCDF
#include <netcdf.inc>
#endif

  !--- namelist definition
  integer           :: nx=128, ny=128, nz=40, nt=2
  integer           :: halo=2, stackmax=1500000, stackmaxd=500000
  logical           :: debug=.FALSE.  
  character(len=64) :: file='test', iospec='-F cachea'
  integer           :: layout(2) = (/0,0/)
  integer           :: ntiles_x=1, ntiles_y=1  ! total number of tiles will be ntiles_x*ntiles_y,
                                               ! the grid size for each tile will be (nx/ntiles_x, ny/ntiles_y) 
                                               ! set ntiles > 1 to test the efficiency of mpp_io.
  integer           :: io_layout(2) = (/0,0/)  ! set io_layout to divide each tile into io_layout(1)*io_layout(2)
                                               ! group and write out data from the root pe of each group.

  namelist / test_mpp_io_nml / nx, ny, nz, nt, halo, stackmax, stackmaxd, debug, file, iospec, &
                               ntiles_x, ntiles_y, layout, io_layout

  integer        :: pe, npes
  type(domain2D) :: domain

  integer            :: tks_per_sec
  integer            :: i,j,k, unit=7
  integer            :: id_single_tile_mult_file
  integer            :: id_mult_tile, id_single_tile_with_group, id_mult_tile_with_group
  logical            :: opened
  character(len=64)  :: varname

  real(DOUBLE_KIND)  :: time
  type(axistype)     :: x, y, z, t
  type(fieldtype)    :: f
  type(domain1D)     :: xdom, ydom
  integer(LONG_KIND) :: rchk, chk

  call mpp_init() 
  pe = mpp_pe()
  npes = mpp_npes()

  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, status='OLD', file='input.nml', err=10 )
  read( unit,test_mpp_io_nml )
  close(unit)
10 continue

  call SYSTEM_CLOCK( count_rate=tks_per_sec )
  if( debug )then
      call mpp_io_init(MPP_DEBUG)
  else
      call mpp_io_init()
  end if
  call mpp_set_stack_size(stackmax)
  call mpp_domains_set_stack_size(stackmaxd)

  if( pe.EQ.mpp_root_pe() )then
      print '(a,6i6)', 'npes, nx, ny, nz, nt, halo=', npes, nx, ny, nz, nt, halo
      print *, 'Using NEW domaintypes and calls...'
  end if

  write( file,'(a,i3.3)' )trim(file), npes

  if(ntiles_x == 1 .and. ntiles_y == 1 .and. io_layout(1) == 1 .and. io_layout(2) == 1) then
     call test_netcdf_io('Simple')
     call test_netcdf_io('Symmetry_T_cell')
     call test_netcdf_io('Symmetry_E_cell')
     call test_netcdf_io('Symmetry_N_cell')
     call test_netcdf_io('Symmetry_C_cell')
     call test_netcdf_io('Symmetry_T_cell_memory')
     call test_netcdf_io('Symmetry_E_cell_memory')
     call test_netcdf_io('Symmetry_N_cell_memory')
     call test_netcdf_io('Symmetry_C_cell_memory')
  else
     if(io_layout(1) <1 .OR. io_layout(2) <1) call mpp_error(FATAL, &
            "program test_mpp_io: both elements of test_mpp_io_nml variable io_layout must be positive integer")
     if(ntiles_x <1 .OR. ntiles_y <1) call mpp_error(FATAL, &
            "program test_mpp_io: mpp_io_nml variable ntiles_x and ntiles_y must be positive integer")
     if(mod(nx, ntiles_x) .NE. 0) call mpp_error(FATAL, &
            "program test_mpp_io: nx must be divided by ntiles_x")
     if(mod(ny, ntiles_y) .NE. 0) call mpp_error(FATAL, &
            "program test_mpp_io: ny must be divided by ntiles_y")
     if(mod(npes, ntiles_x*ntiles_y) .NE. 0) call mpp_error(FATAL, &
            "program test_mpp_io: npes should be divided by ntiles = ntiles_x*ntiles_y ")
     if(layout(1) * layout(2) .NE. npes) call mpp_error(FATAL, &
            "program test_mpp_io: npes should equal to layout(1)*layout(2)" )
     if(mod(layout(1), io_layout(1)) .NE. 0 ) call mpp_error(FATAL, &
            "program test_mpp_io: layout(1) must be divided by io_layout(1)")
     if(mod(layout(2), io_layout(2)) .NE. 0 )call mpp_error(FATAL, &
            "program test_mpp_io: layout(2) must be divided by io_layout(2)")

     id_single_tile_mult_file = mpp_clock_id('Single Tile Multiple File', flags=MPP_CLOCK_SYNC)
     call mpp_clock_begin(id_single_tile_mult_file)
     call test_netcdf_io_mosaic('Single_tile_mult_file', layout, 1, 1, (/1,1/) )
     call mpp_clock_end(id_single_tile_mult_file)

     if(io_layout(1) >1 .OR. io_layout(2) > 1) then
        id_single_tile_with_group = mpp_clock_id('Single Tile With Group', flags=MPP_CLOCK_SYNC)
        call mpp_clock_begin(id_single_tile_with_group)
        call test_netcdf_io_mosaic('Single_tile_with_group', layout, 1, 1, io_layout)
        call mpp_clock_end(id_single_tile_with_group)
     endif

     id_mult_tile  = mpp_clock_id('Multiple Tile', flags=MPP_CLOCK_SYNC)
     call mpp_clock_begin(id_mult_tile)
     if(ntiles_x > 1 .OR. ntiles_y > 1) then
        call test_netcdf_io_mosaic('Mult_tile', layout, ntiles_x, ntiles_y, (/1,1/))
     else
        call test_netcdf_io_mosaic('Mult_tile', layout, io_layout(1), io_layout(2), (/1,1/) )
     endif
     call mpp_clock_end(id_mult_tile)

     if( (io_layout(1) >1 .OR. io_layout(2) > 1) .AND. (ntiles_x >1 .OR. ntiles_y > 1) ) then
        id_mult_tile_with_group = mpp_clock_id('Multiple Tile With Group', flags=MPP_CLOCK_SYNC)
        call mpp_clock_begin(id_mult_tile_with_group)
        call test_netcdf_io_mosaic('Mult_tile_with_group', layout, ntiles_x, ntiles_y, io_layout)
        call mpp_clock_end(id_mult_tile_with_group)
     endif
  endif

  call mpp_io_exit()
  call mpp_domains_exit()
  call mpp_exit()

  contains

  !------------------------------------------------------------------

  subroutine test_netcdf_io(type)
  character(len=*), intent(in) :: type
  integer :: ndim, nvar, natt, ntime
  integer :: is, ie, js, je, isd, ied, jsd, jed, ism, iem, jsm, jem
  integer :: position, msize(2), ioff, joff, nxg, nyg
  logical :: symmetry
  type(atttype),          allocatable :: atts(:)
  type(fieldtype),        allocatable :: vars(:)
  type(axistype),         allocatable :: axes(:)
  real(DOUBLE_KIND),      allocatable :: tstamp(:)
  real, dimension(:,:,:), allocatable :: data, gdata, rdata

  !--- determine the shift and symmetry according to type, 
  select case(type)
  case('Simple')
     position = CENTER; symmetry = .false.
  case('Symmetry_T_cell', 'Symmetry_T_cell_memory')
     position = CENTER; symmetry = .true.
  case('Symmetry_E_cell', 'Symmetry_E_cell_memory')
     position = EAST;   symmetry = .true.
  case('Symmetry_N_cell', 'Symmetry_N_cell_memory')
     position = NORTH;  symmetry = .true.
  case('Symmetry_C_cell', 'Symmetry_C_cell_memory')
     position = CORNER; symmetry = .true.
  case default
     call mpp_error(FATAL, "type = "//type//" is not a valid test type")
  end select

!define domain decomposition
  call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
  if(index(type,"memory") == 0) then  
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, symmetry = symmetry )
  else  ! on memory domain
     msize(1) = nx/layout(1) + 2*halo + 2
     msize(2) = ny/layout(2) + 2*halo + 2
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, symmetry = symmetry, &
                              memory_size = msize )
  end if

  call mpp_get_compute_domain( domain, is,  ie,  js,  je, position=position  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, position=position )
  call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )
  call mpp_get_global_domain ( domain, xsize=nxg, ysize=nyg, position=position )
  call mpp_get_domain_components( domain, xdom, ydom )

!define global data array
  allocate( gdata(nxg,nyg,nz) )
  gdata = 0.
  do k = 1,nz
     do j = 1,nyg
        do i = 1,nxg
           gdata(i,j,k) = k + i*1e-3 + j*1e-6
        end do
     end do
  end do

  ioff = ism - isd; joff = jsm - jsd
  allocate( data(ism:iem,jsm:jem,nz) )
  data = 0
  data(is+ioff:ie+ioff,js+joff:je+joff,:) = gdata(is:ie,js:je,:)

!tests

!sequential write: single-threaded formatted: only if small
  if( nx*ny*nz*nt.LT.1000 .AND. index(type,"memory") .NE. 0 )then
      if( pe.EQ.mpp_root_pe() )print *, 'sequential write: single-threaded formatted'
!here the only test is a successful write: please look at test.txt for verification.
      call mpp_open( unit, trim(file)//'s.txt', action=MPP_OVERWR, form=MPP_ASCII, threading=MPP_SINGLE )
      call mpp_write_meta( unit, x, 'X', 'km', 'X distance', domain=xdom, data=(/(i-1.,i=1,nxg)/) )
      call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
      call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance',              data=(/(i-1.,i=1,nz)/) )
      call mpp_write_meta( unit, t, 'T', 'sec', 'Time' )
      call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data' )
      call mpp_write( unit, x )
      call mpp_write( unit, y )
      call mpp_write( unit, z )
      do i = 0,nt-1
         time = i*10.
         call mpp_write( unit, f, domain, data, time )
      end do
      call mpp_close(unit)
  end if

!netCDF distributed write
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF distributed write'
  call mpp_open( unit, trim(type)//"_"//trim(file)//'d', action=MPP_OVERWR, &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nxg)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z', data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=1 )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  do i = 0,nt-1
     time = i*10.
     call mpp_write( unit, f, domain, data, time )
  end do
  call mpp_close(unit)
  
!netCDF single-threaded write
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF single-threaded write'
  call mpp_open( unit, trim(type)//"_"//trim(file)//'s', action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )

  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nxg)/) )

  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z',              data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=1 )

  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )

  do i = 0,nt-1
     time = i*10.
     call mpp_write( unit, f, domain, data, time)
  end do
  call mpp_close(unit)
  allocate( rdata(is:ie,js:je,nz) )

!netCDF multi-threaded read
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF multi-threaded read'
  call mpp_sync()
  call mpp_open( unit, trim(type)//"_"//trim(file)//'s', action=MPP_RDONLY,  &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
  call mpp_get_info( unit, ndim, nvar, natt, ntime )
  allocate( atts(natt) )
  allocate( axes(ndim) )
  allocate( vars(nvar) )
  allocate( tstamp(ntime) )
  call mpp_get_atts ( unit, atts(:) )
  call mpp_get_axes ( unit, axes(:) )
  call mpp_get_fields ( unit, vars(:) )
  call mpp_get_times( unit, tstamp(:) )

  call mpp_get_atts(vars(1),name=varname)

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )
  call mpp_read( unit, vars(1), domain, rdata, 1 )
  rchk = mpp_chksum(rdata(is:ie,js:je,:))
  chk  = mpp_chksum( data(is+ioff:ie+ioff,js+joff:je+joff,:))
  if( pe.EQ.mpp_root_pe() )print '(a,2z18)', trim(type)//' checksum=', rchk, chk
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//': single-fileset: data comparison are OK.' )
  else
      call mpp_error( FATAL, 'Checksum error on multi-threaded/single-fileset netCDF read for type ' &
               //trim(type) )
  end if

  deallocate( atts, axes, vars, tstamp )

!netCDF distributed read
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF multi-threaded read'
  call mpp_sync()               !wait for previous write to complete
  call mpp_open( unit, trim(type)//"_"//trim(file)//'d', action=MPP_RDONLY,  &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  call mpp_get_info( unit, ndim, nvar, natt, ntime )
  allocate( atts(natt) )
  allocate( axes(ndim) )
  allocate( vars(nvar) )
  allocate( tstamp(ntime) )
  call mpp_get_atts ( unit, atts(:) )
  call mpp_get_axes ( unit, axes(:) )
  call mpp_get_fields ( unit, vars(:) )
  call mpp_get_times( unit, tstamp(:) )

  call mpp_get_atts(vars(1),name=varname)
  rdata = 0

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )

  call mpp_read( unit, vars(1), domain, rdata, 1 )

  rchk = mpp_chksum(rdata(is:ie,js:je,:))
  chk  = mpp_chksum( data(is+ioff:ie+ioff,js+joff:je+joff,:))
  if( pe.EQ.mpp_root_pe() )print '(a,2z18)', trim(type)//' checksum=', rchk, chk
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//': multi-fileset: data comparison are OK.' )
  else
      call mpp_error( FATAL, 'Checksum error on multi-threaded/multi-fileset netCDF read for type ' &
           //trim(type) )
  end if

  deallocate( atts, axes, vars, tstamp )

  deallocate( rdata, gdata, data)

  end subroutine test_netcdf_io


  subroutine test_netcdf_io_mosaic(type, layout, ntiles_x, ntiles_y, io_layout)
  character(len=*), intent(in) :: type
  integer,          intent(in) :: layout(:)
  integer,          intent(in) :: io_layout(:)
  integer,          intent(in) :: ntiles_x, ntiles_y

  integer                              :: ndim, nvar, natt, ntime
  integer                              :: isc, iec, jsc, jec, nlon, nlat, n, i, j
  integer                              :: my_tile, ncontacts, npes_per_tile, ntiles
  integer, dimension(:),   allocatable :: tile1, istart1, iend1, jstart1, jend1
  integer, dimension(:),   allocatable :: tile2, istart2, iend2, jstart2, jend2
  integer, dimension(:),   allocatable :: pe_start, pe_end
  integer, dimension(:,:), allocatable :: layout2D, global_indices
  character(len=64)                    :: output_file
  logical                              :: is_root_pe
  real, dimension(:,:,:), allocatable  :: data, rdata
  type(fieldtype), save                :: vars(1)
  integer                              :: id_clock_read, id_clock_write

  ! first get number of tiles of this mosaic. when there is one tile,
  ! the file will be read/write using distributed file.
  ! when there is more than one tile, single fileset will be used
  npes = mpp_npes()
  
  id_clock_read  = mpp_clock_id(trim(type)//" read", flags=MPP_CLOCK_SYNC)
  id_clock_write = mpp_clock_id(trim(type)//" write", flags=MPP_CLOCK_SYNC)

  ncontacts = 0
  ntiles = ntiles_x*ntiles_y

  npes_per_tile = npes/ntiles
  my_tile       = mpp_pe()/npes_per_tile + 1
  is_root_pe = .false.
  if(mpp_pe() == (my_tile-1)*npes_per_tile ) is_root_pe = .true.

  allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
  !--- for simplify purpose, we assume all the tiles have the same size.
  do n = 1, ntiles
     pe_start(n) = (n-1)*npes_per_tile
     pe_end(n)   = n*npes_per_tile-1
  end do
  if(ntiles>1) then
     nlon = nx/ntiles_x
     nlat = ny/ntiles_y
  else
     nlon = nx
     nlat = ny
  endif
  
  do n = 1, ntiles
     global_indices(:,n) = (/1,nlon,1,nlat/)
     layout2D(1,n)         = layout(1)/ntiles_x
     layout2D(2,n)         = layout(2)/ntiles_y
  end do

  call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, ncontacts, tile1, tile2, &
                         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2, pe_start, pe_end, &
                         name = type)
  call mpp_get_compute_domain( domain, isc,  iec,  jsc,  jec  )
  call mpp_get_domain_components(domain, xdom, ydom)
  allocate( data (isc:iec,jsc:jec,nz) )
  allocate( rdata(isc:iec,jsc:jec,nz) )
  do k = 1,nz
     do j = jsc, jec
        do i = isc, iec
           data(i,j,k)  = k + i*1e-3 + j*1e-6
        enddo
     enddo
  enddo

  !--- netcdf distribute write if ntiles = 1, otherwise single-thread write
  output_file = type
  select case(type)
  case("Single_tile_single_file")
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE )
  case("Single_tile_mult_file")
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  case("Mult_tile")
     write(output_file, '(a,I4.4)') type//'.tile', my_tile
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE, is_root_pe=is_root_pe )
  case("Single_tile_with_group")
     call mpp_define_io_domain(domain, io_layout)
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)
  case("Mult_tile_with_group")
     write(output_file, '(a,I4.4)') type//'.tile', my_tile
     call mpp_define_io_domain(domain, io_layout)
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)

  case default
     call mpp_error(FATAL, "program test_mpp_io: invaid value of type="//type)
  end select  

  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nlon)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nlat)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z', data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=1 )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  call mpp_clock_begin(id_clock_write)
  do i = 0,nt-1
     time = i*10.
     call mpp_write( unit, f, domain, data, time )
  end do
  call mpp_clock_end(id_clock_write)
  call mpp_close(unit)
  
  call mpp_sync()               !wait for previous write to complete

  select case(type)
  case("Single_tile_single_file")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
  case("Single_tile_mult_file")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  case("Mult_tile")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, &
         fileset=MPP_SINGLE, is_root_pe=is_root_pe )
  case("Single_tile_with_group", "Mult_tile_with_group")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)
  case default
     call mpp_error(FATAL, "program test_mpp_io: invaid value of type="//type)
  end select  

  call mpp_get_info( unit, ndim, nvar, natt, ntime )
  call mpp_get_fields ( unit, vars(:) )
  call mpp_get_atts(vars(1),name=varname)

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )
  call mpp_clock_begin(id_clock_read)
  do i = 0,nt-1
     call mpp_read( unit, vars(1), domain, rdata, 1 )
  enddo
  call mpp_clock_end(id_clock_read)
  rchk = mpp_chksum(rdata)
  chk  = mpp_chksum( data)
  if( pe.EQ.mpp_root_pe() )print '(a,2z18)', trim(type)//' checksum=', rchk, chk
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//': data comparison are OK.' )
  else
      call mpp_error( FATAL, 'Checksum error on netCDF read for type ' &
               //trim(type) )
  end if

!  deallocate( vars)

  deallocate( rdata, data)
  call mpp_deallocate_domain(domain)

  end subroutine test_netcdf_io_mosaic

end program test

#else
module null_mpp_io_test
end module
#endif
