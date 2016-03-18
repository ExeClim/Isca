#ifdef test_fms_io

 program fms_io_test
#include <fms_platform.h>

 use mpp_mod,         only: mpp_pe, mpp_npes, mpp_root_pe, mpp_init, mpp_exit
 use mpp_mod,         only: stdout, mpp_error, FATAL, NOTE, mpp_chksum
 use mpp_domains_mod, only: domain2D, mpp_define_layout, mpp_define_mosaic
 use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
 use mpp_domains_mod, only: mpp_domains_init, mpp_domains_exit
 use mpp_domains_mod, only: mpp_domains_set_stack_size, mpp_define_io_domain
 use mpp_io_mod,      only: mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY
 use fms_io_mod,      only: read_data, write_data, fms_io_init, fms_io_exit
 use fms_io_mod,      only: file_exist, register_restart_field, save_restart, restore_state
 use fms_io_mod,      only: restart_file_type
 use mpp_io_mod,      only: MAX_FILE_SIZE

 implicit none

 integer :: sizex_latlon_grid = 144
 integer :: sizey_latlon_grid = 90
 integer :: size_cubic_grid = 48
 integer :: nz = 10, nt = 2, halo = 1
 integer :: stackmax =4000000
 integer :: num_step = 4 ! number of time steps to run, this is used for intermediate run.
                         ! set num_step = 0 for no intermediate run.
 logical :: do_write=.true. ! set this to false for high resolution and single file,
                            ! split file capability is not implemented for write_data
 integer :: layout_cubic (2) = (/0,0/)
 integer :: layout_latlon(2) = (/0,0/)  
 integer :: io_layout(2) = (/0,0/) ! set ndivs_x and ndivs_y to divide each tile into io_layout(1)*io_layout(2)
                                   ! group and write out data from the root pe of each group.

 namelist /test_fms_io_nml/ sizex_latlon_grid, sizey_latlon_grid, size_cubic_grid, &
                            nz, nt, halo, num_step, stackmax, do_write, layout_cubic, layout_latlon, io_layout

 integer           :: unit, io_status, step
 character(len=20) :: time_stamp

 type data_storage_type
    real,    allocatable, dimension(:,:,:,:) :: data1_r3d, data2_r3d, data1_r3d_read, data2_r3d_read
    real,    allocatable, dimension(:,:,:)   :: data1_r2d, data2_r2d, data1_r2d_read, data2_r2d_read
    real,    allocatable, dimension(:,:)     :: data1_r1d, data2_r1d, data1_r1d_read, data2_r1d_read
    real,    allocatable, dimension(:)       :: data1_r0d, data2_r0d, data1_r0d_read, data2_r0d_read
    integer, allocatable, dimension(:,:,:,:) :: data1_i3d, data2_i3d, data1_i3d_read, data2_i3d_read
    integer, allocatable, dimension(:,:,:)   :: data1_i2d, data2_i2d, data1_i2d_read, data2_i2d_read
    integer, allocatable, dimension(:,:)     :: data1_i1d, data2_i1d, data1_i1d_read, data2_i1d_read
    integer, allocatable, dimension(:)       :: data1_i0d, data2_i0d, data1_i0d_read, data2_i0d_read
 end type data_storage_type
 
 type(data_storage_type), save :: latlon_data
 type(data_storage_type), save :: cubic_data
 type(domain2d),          save :: domain_latlon
 type(domain2d),          save :: domain_cubic
 type(restart_file_type), save :: restart_latlon
 type(restart_file_type), save :: restart_cubic
 integer                       :: ntile_latlon = 1
 integer                       :: ntile_cubic = 6
 integer                       :: npes

 character(len=128) :: file_latlon, file_cubic

 call mpp_init
 npes = mpp_npes()

 call mpp_domains_init  

 call fms_io_init

 if (file_exist('input.nml') )then
    call mpp_open(unit, 'input.nml',form=MPP_ASCII,action=MPP_RDONLY)
    read(unit,test_fms_io_nml,iostat=io_status)

    if (io_status > 0) then
     call mpp_error(FATAL,'=>test_fms_io: Error reading test_fms_io_nml')
  endif
    call mpp_close (unit)
 end if

 write(stdout(), test_fms_io_nml )
  call mpp_domains_set_stack_size(stackmax)
 !--- currently we assume at most two time level will be written to restart file.
 if(nt > 2) call mpp_error(FATAL, "test_fms_io: test_fms_io_nml variable nt should be no larger than 2")

 file_latlon   = "test.res.latlon_grid.save_restore.nc"
 file_cubic    = "test.res.cubic_grid.save_restore.nc"

 call setup_test_restart(restart_latlon, "latlon_grid", ntile_latlon, latlon_data, file_latlon, layout_latlon, domain_latlon)
 call setup_test_restart(restart_cubic,  "cubic_grid", ntile_cubic, cubic_data, file_cubic, layout_cubic, domain_cubic )

 if(file_exist('INPUT/'//trim(file_latlon), domain_latlon)) then
    call restore_state(restart_latlon)
    call compare_restart("latlon_grid save_restore", latlon_data)
 end if
 if(file_exist('INPUT/'//trim(file_cubic), domain_cubic) ) then
    call restore_state(restart_cubic)
    call compare_restart("cubic_grid save_restore", cubic_data)
 end if
 
 !---copy data
 if(mod(npes,ntile_latlon) == 0) call copy_restart_data(latlon_data)
 if(mod(npes,ntile_cubic) == 0 ) call copy_restart_data(cubic_data)

 do step = 1, num_step
    write(time_stamp, '(a,I4.4)') "step", step
    if(mod(npes,ntile_latlon) == 0) call save_restart(restart_latlon, time_stamp)
    if(mod(npes,ntile_cubic) == 0 ) call save_restart(restart_cubic, time_stamp)
 end do
 if(mod(npes,ntile_latlon) == 0) call save_restart(restart_latlon)
 if(mod(npes,ntile_cubic)  == 0) call save_restart(restart_cubic)

 if(mod(npes,ntile_latlon) == 0) call release_storage_memory(latlon_data)
 if(mod(npes,ntile_cubic) == 0 ) call release_storage_memory(cubic_data)

 if(mod(npes,ntile_cubic) == 0 ) call mpp_error(NOTE, "test_fms_io: restart test is done for latlon_grid")
 if(mod(npes,ntile_cubic) == 0 ) call mpp_error(NOTE, "test_fms_io: restart test is done for cubic_grid")

 call fms_io_exit
 call mpp_domains_exit
 call mpp_exit

contains

  !******************************************************************************
  subroutine setup_test_restart(restart_data, type, ntiles, storage, file, layout_in, domain)
    type(restart_file_type),   intent(inout) :: restart_data
    character(len=*), intent(in)             :: type
    integer,          intent(in)             :: ntiles
    type(data_storage_type), intent(inout)   :: storage
    character(len=*), intent(in)             :: file  
    integer,          intent(in)             :: layout_in(:)
    type(domain2d),   intent(inout)          :: domain
    character(len=128)                       :: file_r
    character(len=128)                       :: file_w
    integer                                  :: pe, npes_per_tile, tile
    integer                                  :: num_contact
    integer                                  :: n, layout(2)
    integer, allocatable, dimension(:,:)     :: global_indices, layout2D
    integer, allocatable, dimension(:)       :: pe_start, pe_end
    integer, dimension(1)                    :: tile1, tile2
    integer, dimension(1)                    :: istart1, iend1, jstart1, jend1
    integer, dimension(1)                    :: istart2, iend2, jstart2, jend2
    integer                                  :: i, j, k, nx, ny
    integer                                  :: isc, iec, jsc, jec
    integer                                  :: isd, ied, jsd, jed
    integer                                  :: id_restart

    file_r = "INPUT/test.res."//trim(type)//".read_write.nc"
    file_w = "RESTART/test.res."//trim(type)//".read_write.nc"

    select case(type)
    case("latlon_grid")
       nx = sizex_latlon_grid
       ny = sizey_latlon_grid
    case("cubic_grid")
       nx = size_cubic_grid
       ny = size_cubic_grid
    case default
       call mpp_error(FATAL, "test_fms_io: "//type//" is not a valid option")
    end select

    pe   = mpp_pe()
    if(mod(npes,ntiles) .NE. 0) then
       call mpp_error(NOTE, "test_fms_io: npes can not be divided by ntiles, no test will be done for "//trim(type))
       return
    end if
    npes_per_tile = npes/ntiles
    tile = pe/npes_per_tile + 1

    if(layout_in(1)*layout_in(2) == npes_per_tile) then
       layout = layout_in
    else
       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
    endif
    if(io_layout(1) <1 .OR. io_layout(2) <1) call mpp_error(FATAL, &
            "program test_fms_io: both elements of test_fms_io_nml variable io_layout must be positive integer")
    if(mod(layout(1), io_layout(1)) .NE. 0 ) call mpp_error(FATAL, &
         "program test_fms_io: layout(1) must be divided by io_layout(1)")
    if(mod(layout(2), io_layout(2)) .NE. 0 ) call mpp_error(FATAL, &
         "program test_fms_io: layout(2) must be divided by io_layout(2)")    
    allocate(global_indices(4,ntiles), layout2D(2,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)       = layout
       pe_start(n)         = (n-1)*npes_per_tile
       pe_end(n)           = n*npes_per_tile-1
    end do
    num_contact = 0
    call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                           istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                           pe_start, pe_end, whalo=halo, ehalo=halo, shalo=halo, nhalo=halo, name = type  )
    if(io_layout(1) .NE. 1 .OR. io_layout(2) .NE. 1) call mpp_define_io_domain(domain, io_layout)

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)

    allocate(storage%data1_r3d(isd:ied, jsd:jed, nz, nt), storage%data1_r3d_read(isd:ied, jsd:jed, nz, nt) )
    allocate(storage%data2_r3d(isd:ied, jsd:jed, nz, nt), storage%data2_r3d_read(isd:ied, jsd:jed, nz, nt) )
    allocate(storage%data1_i3d(isd:ied, jsd:jed, nz, nt), storage%data1_i3d_read(isd:ied, jsd:jed, nz, nt) )
    allocate(storage%data2_i3d(isd:ied, jsd:jed, nz, nt), storage%data2_i3d_read(isd:ied, jsd:jed, nz, nt) )
    allocate(storage%data1_r2d(isd:ied, jsd:jed,     nt), storage%data1_r2d_read(isd:ied, jsd:jed,     nt) )
    allocate(storage%data2_r2d(isd:ied, jsd:jed,     nt), storage%data2_r2d_read(isd:ied, jsd:jed,     nt) )
    allocate(storage%data1_i2d(isd:ied, jsd:jed,     nt), storage%data1_i2d_read(isd:ied, jsd:jed,     nt) )
    allocate(storage%data2_i2d(isd:ied, jsd:jed,     nt), storage%data2_i2d_read(isd:ied, jsd:jed,     nt) )
    allocate(storage%data1_r1d(                  nz, nt), storage%data1_r1d_read(                  nz, nt) )
    allocate(storage%data2_r1d(                  nz, nt), storage%data2_r1d_read(                  nz, nt) )
    allocate(storage%data1_i1d(                  nz, nt), storage%data1_i1d_read(                  nz, nt) )
    allocate(storage%data2_i1d(                  nz, nt), storage%data2_i1d_read(                  nz, nt) )
    allocate(storage%data1_r0d(                      nt), storage%data1_r0d_read(                      nt) )
    allocate(storage%data2_r0d(                      nt), storage%data2_r0d_read(                      nt) )
    allocate(storage%data1_i0d(                      nt), storage%data1_i0d_read(                      nt) )
    allocate(storage%data2_i0d(                      nt), storage%data2_i0d_read(                      nt) )

    storage%data1_r3d = 0; storage%data1_r3d_read = 0; storage%data2_r3d = 0; storage%data2_r3d_read = 0
    storage%data1_i3d = 0; storage%data1_i3d_read = 0; storage%data2_i3d = 0; storage%data2_i3d_read = 0
    storage%data1_r2d = 0; storage%data1_r2d_read = 0; storage%data2_r2d = 0; storage%data2_r2d_read = 0
    storage%data1_i2d = 0; storage%data1_i2d_read = 0; storage%data2_i2d = 0; storage%data2_i2d_read = 0
    storage%data1_r1d = 0; storage%data1_r1d_read = 0; storage%data2_r1d = 0; storage%data2_r1d_read = 0
    storage%data1_i1d = 0; storage%data1_i1d_read = 0; storage%data2_i1d = 0; storage%data2_i1d_read = 0
    storage%data1_r0d = 0; storage%data1_r0d_read = 0; storage%data2_r0d = 0; storage%data2_r0d_read = 0
    storage%data1_i0d = 0; storage%data1_i0d_read = 0; storage%data2_i0d = 0; storage%data2_i0d_read = 0
    do n = 1, nt
       storage%data1_r0d(n) =  tile + n*1e-3
       storage%data2_r0d(n) = -tile - n*1e-3
       storage%data1_i0d(n) =  tile*1e3 + n
       storage%data2_i0d(n) = -tile*1e3 - n
       do k = 1, nz
          storage%data1_r1d(k,n) =   tile*1e3 + n + k*1e-3
          storage%data2_r1d(k,n) =  -tile*1e3 - n - k*1e-3
          storage%data1_i1d(k,n) =   tile*1e6 + n*1e3 + k
          storage%data2_i1d(k,n) =  -tile*1e6 - n*1e3 - k
          do j = jsc, jec
             do i = isc, iec
                storage%data1_r3d(i,j,k,n) =  tile*1e6 + n*1e3 + k + i*1e-3 + j*1e-6; 
                storage%data2_r3d(i,j,k,n) = -tile*1e6 - n*1e3 - k - i*1e-3 - j*1e-6; 
                storage%data1_i3d(i,j,k,n) =  tile*1e9 + n*1e8 + k*1e6 + i*1e3 + j; 
                storage%data2_i3d(i,j,k,n) = -tile*1e9 - n*1e8 - k*1e6 - i*1e3 - j; 
             end do
          end do
       end do

       do j = jsc, jec
          do i = isc, iec
             storage%data1_r2d(i,j,n) =  tile*1e1 + n + i*1e-3 + j*1e-6; 
             storage%data2_r2d(i,j,n) = -tile*1e1 - n - i*1e-3 - j*1e-6; 
             storage%data1_i2d(i,j,n) =  tile*1e7 + n*1e6 + i*1e3 + j; 
             storage%data2_i2d(i,j,n) = -tile*1e7 - n*1e6 - i*1e3 - j; 
          end do
       end do
    end do
    if(file_exist(file_r, domain)) then
       do n = 1, nt
          call read_data(file_r, "data1_r3d", storage%data1_r3d_read(:,:,:,n), domain, timelevel = n )
          call read_data(file_r, "data2_r3d", storage%data2_r3d_read(:,:,:,n), domain, timelevel = n )
          call read_data(file_r, "data1_i3d", storage%data1_i3d_read(:,:,:,n), domain, timelevel = n )
          call read_data(file_r, "data2_i3d", storage%data2_i3d_read(:,:,:,n), domain, timelevel = n )
          call read_data(file_r, "data1_r2d", storage%data1_r2d_read(:,:,  n), domain, timelevel = n )
          call read_data(file_r, "data2_r2d", storage%data2_r2d_read(:,:,  n), domain, timelevel = n )
          call read_data(file_r, "data1_i2d", storage%data1_i2d_read(:,:,  n), domain, timelevel = n )
          call read_data(file_r, "data2_i2d", storage%data2_i2d_read(:,:,  n), domain, timelevel = n )
          call read_data(file_r, "data1_r1d", storage%data1_r1d_read(:,    n), domain, timelevel = n )
          call read_data(file_r, "data2_r1d", storage%data2_r1d_read(:,    n), domain, timelevel = n )
          call read_data(file_r, "data1_i1d", storage%data1_i1d_read(:,    n), domain, timelevel = n )
          call read_data(file_r, "data2_i1d", storage%data2_i1d_read(:,    n), domain, timelevel = n )
          call read_data(file_r, "data1_r0d", storage%data1_r0d_read(      n), domain, timelevel = n )
          call read_data(file_r, "data2_r0d", storage%data2_r0d_read(      n), domain, timelevel = n )
          call read_data(file_r, "data1_i0d", storage%data1_i0d_read(      n), domain, timelevel = n )
          call read_data(file_r, "data2_i0d", storage%data2_i0d_read(      n), domain, timelevel = n )
       end do
       call compare_restart(type//" read_write", storage)
    end if


    !--- high resolution restart is not implemented for write data
    if(do_write ) then 
       do n = 1, nt
          call write_data(file_w, "data1_r3d", storage%data1_r3d(:,:,:,n), domain )
          call write_data(file_w, "data2_r3d", storage%data2_r3d(:,:,:,n), domain )
          call write_data(file_w, "data1_i3d", storage%data1_i3d(:,:,:,n), domain )
          call write_data(file_w, "data2_i3d", storage%data2_i3d(:,:,:,n), domain )
          call write_data(file_w, "data1_r2d", storage%data1_r2d(:,:,  n), domain )
          call write_data(file_w, "data2_r2d", storage%data2_r2d(:,:,  n), domain )
          call write_data(file_w, "data1_i2d", storage%data1_i2d(:,:,  n), domain )
          call write_data(file_w, "data2_i2d", storage%data2_i2d(:,:,  n), domain )
          call write_data(file_w, "data1_r1d", storage%data1_r1d(:,    n), domain )
          call write_data(file_w, "data2_r1d", storage%data2_r1d(:,    n), domain )
          call write_data(file_w, "data1_i1d", storage%data1_i1d(:,    n), domain )
          call write_data(file_w, "data2_i1d", storage%data2_i1d(:,    n), domain )
          call write_data(file_w, "data1_r0d", storage%data1_r0d(      n), domain )
          call write_data(file_w, "data2_r0d", storage%data2_r0d(      n), domain )
          call write_data(file_w, "data1_i0d", storage%data1_i0d(      n), domain )
          call write_data(file_w, "data2_i0d", storage%data2_i0d(      n), domain )
       end do
    end if

    !--- test register_restart_field, save_restart, restore_state

    id_restart = register_restart_field(restart_data, file, "data1_r3d", storage%data1_r3d_read(:,:,:,1), &
                                domain, longname="first data_r3d",units="none")
    id_restart = register_restart_field(restart_data, file, "data1_r3d", storage%data1_r3d_read(:,:,:,2), &
                                domain, longname="first data_r3d",units="none")
    id_restart = register_restart_field(restart_data, file, "data2_r3d", storage%data2_r3d_read(:,:,:,1), &
                                storage%data2_r3d_read(:,:,:,2), &
                                domain, longname="second data_i3d", units="none")

    id_restart = register_restart_field(restart_data, file, "data1_i3d", storage%data1_i3d_read(:,:,:,1), &
                                domain, longname="first data_i3d",units="none")
    id_restart = register_restart_field(restart_data, file, "data1_i3d", storage%data1_i3d_read(:,:,:,2), &
                                domain, longname="first data_i3d",units="none")
    id_restart = register_restart_field(restart_data, file, "data2_i3d", storage%data2_i3d_read(:,:,:,1), &
                                storage%data2_i3d_read(:,:,:,2), &
                                domain, longname="second data_i3d", units="none")

    id_restart = register_restart_field(restart_data, file, "data1_r2d", storage%data1_r2d_read(:,:,  1), &
                                domain, longname="first data_r2d",units="none")
    id_restart = register_restart_field(restart_data, file, "data1_r2d", storage%data1_r2d_read(:,:,  2), &
                                domain, longname="first data_r2d",units="none")
    id_restart = register_restart_field(restart_data, file, "data2_r2d", storage%data2_r2d_read(:,:,  1), &
                                storage%data2_r2d_read(:,:,2), &
                                domain, longname="second data_i2d", units="none")

    id_restart = register_restart_field(restart_data, file, "data1_i2d", storage%data1_i2d_read(:,:,  1), &
                                domain, longname="first data_i2d",units="none")
    id_restart = register_restart_field(restart_data, file, "data1_i2d", storage%data1_i2d_read(:,:,  2), &
                                domain, longname="first data_i2d",units="none")
    id_restart = register_restart_field(restart_data, file, "data2_i2d", storage%data2_i2d_read(:,:,  1), &
                                storage%data2_i2d_read(:,:,2), &
                                domain, longname="second data_i2d", units="none")

    id_restart = register_restart_field(restart_data, file, "data1_r1d", storage%data1_r1d_read(:,    1), &
                                domain, longname="first data_r1d",units="none")
    id_restart = register_restart_field(restart_data, file, "data1_r1d", storage%data1_r1d_read(:,    2), &
                                domain, longname="first data_r1d",units="none")
    id_restart = register_restart_field(restart_data, file, "data2_r1d", storage%data2_r1d_read(:,    1), &
                                storage%data2_r1d_read(:,  2), &
                                domain, longname="second data_i1d", units="none")

    id_restart = register_restart_field(restart_data, file, "data1_i1d", storage%data1_i1d_read(:,    1), &
                                domain, longname="first data_i1d",units="none")
    id_restart = register_restart_field(restart_data, file, "data1_i1d", storage%data1_i1d_read(:,    2), &
                                domain, longname="first data_i1d",units="none")
    id_restart = register_restart_field(restart_data, file, "data2_i1d", storage%data2_i1d_read(:,    1), &
                                storage%data2_i1d_read(:,  2), &
                                domain, longname="second data_i1d", units="none")


    id_restart = register_restart_field(restart_data, file, "data1_r0d", storage%data1_r0d_read(      1), &
                                domain, longname="first data_r0d",units="none")
    id_restart = register_restart_field(restart_data, file, "data1_r0d", storage%data1_r0d_read(      2), &
                                domain, longname="first data_r0d",units="none")
    id_restart = register_restart_field(restart_data, file, "data2_r0d", storage%data2_r0d_read(      1), &
                                storage%data2_r0d_read(    2), &
                                domain, longname="second data_i0d", units="none")

    id_restart = register_restart_field(restart_data, file, "data1_i0d", storage%data1_i0d_read(      1), &
                                domain, longname="first data_i0d",units="none")
    id_restart = register_restart_field(restart_data, file, "data1_i0d", storage%data1_i0d_read(      2), &
                                domain, longname="first data_i0d",units="none")
    id_restart = register_restart_field(restart_data, file, "data2_i0d", storage%data2_i0d_read(      1), &
                                storage%data2_i0d_read(    2), &
                                domain, longname="second data_i0d", units="none")

  end subroutine setup_test_restart

  subroutine compare_restart(type, storage)
    character(len=*), intent(in)             :: type
    type(data_storage_type), intent(inout)   :: storage

       call compare_data_r4d(storage%data1_r3d, storage%data1_r3d_read, type//" data1_r3d")
       call compare_data_r4d(storage%data2_r3d, storage%data2_r3d_read, type//" data2_r3d")
       call compare_data_i4d(storage%data1_i3d, storage%data1_i3d_read, type//" data1_i3d")
       call compare_data_i4d(storage%data2_i3d, storage%data2_i3d_read, type//" data2_i3d")
       call compare_data_r3d(storage%data1_r2d, storage%data1_r2d_read, type//" data1_r2d")
       call compare_data_r3d(storage%data2_r2d, storage%data2_r2d_read, type//" data2_r2d")
       call compare_data_i3d(storage%data1_i2d, storage%data1_i2d_read, type//" data1_i2d")
       call compare_data_i3d(storage%data2_i2d, storage%data2_i2d_read, type//" data2_i2d")
       call compare_data_r2d(storage%data1_r1d, storage%data1_r1d_read, type//" data1_r1d")
       call compare_data_r2d(storage%data2_r1d, storage%data2_r1d_read, type//" data2_r1d")
       call compare_data_i2d(storage%data1_i1d, storage%data1_i1d_read, type//" data1_i1d")
       call compare_data_i2d(storage%data2_i1d, storage%data2_i1d_read, type//" data2_i1d")
       call compare_data_r1d(storage%data1_r0d, storage%data1_r0d_read, type//" data1_r0d")
       call compare_data_r1d(storage%data2_r0d, storage%data2_r0d_read, type//" data2_r0d")
       call compare_data_i1d(storage%data1_i0d, storage%data1_i0d_read, type//" data1_i0d")
       call compare_data_i1d(storage%data2_i0d, storage%data2_i0d_read, type//" data2_i0d")

  end subroutine compare_restart

  subroutine release_storage_memory(storage)
    type(data_storage_type), intent(inout)   :: storage

    deallocate(storage%data1_r3d, storage%data2_r3d, storage%data1_r3d_read, storage%data2_r3d_read)
    deallocate(storage%data1_i3d, storage%data2_i3d, storage%data1_i3d_read, storage%data2_i3d_read)
    deallocate(storage%data1_r2d, storage%data2_r2d, storage%data1_r2d_read, storage%data2_r2d_read)
    deallocate(storage%data1_i2d, storage%data2_i2d, storage%data1_i2d_read, storage%data2_i2d_read)
    deallocate(storage%data1_r1d, storage%data2_r1d, storage%data1_r1d_read, storage%data2_r1d_read)
    deallocate(storage%data1_i1d, storage%data2_i1d, storage%data1_i1d_read, storage%data2_i1d_read)
    deallocate(storage%data1_r0d, storage%data2_r0d, storage%data1_r0d_read, storage%data2_r0d_read)
    deallocate(storage%data1_i0d, storage%data2_i0d, storage%data1_i0d_read, storage%data2_i0d_read)

  end subroutine release_storage_memory

  subroutine copy_restart_data(storage)
    type(data_storage_type), intent(inout)   :: storage

    storage%data1_r3d_read = storage%data1_r3d; storage%data2_r3d_read = storage%data2_r3d
    storage%data1_i3d_read = storage%data1_i3d; storage%data2_i3d_read = storage%data2_i3d
    storage%data1_r2d_read = storage%data1_r2d; storage%data2_r2d_read = storage%data2_r2d
    storage%data1_i2d_read = storage%data1_i2d; storage%data2_i2d_read = storage%data2_i2d
    storage%data1_r1d_read = storage%data1_r1d; storage%data2_r1d_read = storage%data2_r1d
    storage%data1_i1d_read = storage%data1_i1d; storage%data2_i1d_read = storage%data2_i1d
    storage%data1_r0d_read = storage%data1_r0d; storage%data2_r0d_read = storage%data2_r0d
    storage%data1_i0d_read = storage%data1_i0d; storage%data2_i0d_read = storage%data2_i0d

    return

  end subroutine copy_restart_data

  subroutine compare_data_r4d( a, b, string )
    real, intent(in), dimension(:,:,:,:) :: a, b
    character(len=*), intent(in)         :: string
    integer(LONG_KIND)                   :: sum1, sum2
    integer                              :: i, j, k, l
    integer, parameter                   :: stdunit = 6

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) .or. size(a,4) .ne. size(b,4) ) &
         call mpp_error(FATAL,'compare_data_r4d: size of a and b does not match')

    do l = 1, size(a,4)
       do k = 1, size(a,3)
          do j = 1, size(a,2)
             do i = 1, size(a,1)
                if(a(i,j,k,l) .ne. b(i,j,k,l)) then
                   write(stdunit,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,f18.9,a,f18.9)')" at pe ", mpp_pe(), &
                        ", at point (",i,", ", j, ", ", k, ", ", l, "), a = ", a(i,j,k,l), ", b = ", b(i,j,k,l)
                   call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
                endif
             enddo
          enddo
       enddo
    enddo
    sum1 = mpp_chksum( a, (/mpp_pe()/) )
    sum2 = mpp_chksum( b, (/mpp_pe()/) )

    if( sum1.EQ.sum2 )then
       if( mpp_pe() .EQ. mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays 
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_data_r4d

  subroutine compare_data_i4d( a, b, string )
    integer, intent(in), dimension(:,:,:,:) :: a, b
    character(len=*), intent(in)            :: string
    real                                    :: real_a(size(a,1),size(a,2),size(a,3),size(a,4))
    real                                    :: real_b(size(b,1),size(b,2),size(b,3),size(b,4))

    real_a = a 
    real_b = b
    call compare_data_r4d(real_a, real_b, string)

  end subroutine compare_data_i4d


  subroutine compare_data_r3d( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in)       :: string
    integer(LONG_KIND)                 :: sum1, sum2
    integer                            :: i, j, l
    integer, parameter                 :: stdunit = 6

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
         call mpp_error(FATAL,'compare_data_r3d: size of a and b does not match')

    do l = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,l) .ne. b(i,j,l)) then
                write(stdunit,'(a,i3,a,i3,a,i3,a,i3,a,f16.9,a,f16.9)')" at pe ", mpp_pe(), &
                     ", at point (",i,", ", j, ", ", l, "), a = ", a(i,j,l), ", b = ", b(i,j,l)
                call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
             endif
          enddo
       enddo
    enddo
    sum1 = mpp_chksum( a, (/mpp_pe()/) )
    sum2 = mpp_chksum( b, (/mpp_pe()/) )

    if( sum1.EQ.sum2 )then
       if( mpp_pe() .EQ. mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays 
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_data_r3d

  subroutine compare_data_i3d( a, b, string )
    integer, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in)          :: string
    real                                  :: real_a(size(a,1),size(a,2),size(a,3))
    real                                  :: real_b(size(b,1),size(b,2),size(b,3))

    real_a = a 
    real_b = b
    call compare_data_r3d(real_a, real_b, string)

  end subroutine compare_data_i3d


  subroutine compare_data_r2d( a, b, string )
    real, intent(in), dimension(:,:) :: a, b
    character(len=*), intent(in)     :: string
    integer(LONG_KIND)               :: sum1, sum2
    integer                          :: i, l
    integer, parameter               :: stdunit = 6

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) ) &
         call mpp_error(FATAL,'compare_data_r2d: size of a and b does not match')

    do l = 1, size(a,2)
       do i = 1, size(a,1)
          if(a(i,l) .ne. b(i,l)) then
             write(stdunit,'(a,i3,a,i3,a,i3,a,f16.9,a,f16.9)')" at pe ", mpp_pe(), &
                  ", at point (",i, ", ", l, "), a = ", a(i,l), ", b = ", b(i,l)
             call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
          endif
       enddo
    end do
    sum1 = mpp_chksum( a, (/mpp_pe()/) )
    sum2 = mpp_chksum( b, (/mpp_pe()/) )

    if( sum1.EQ.sum2 )then
       if( mpp_pe() .EQ. mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays 
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_data_r2d

  subroutine compare_data_i2d( a, b, string )
    integer, intent(in), dimension(:,:) :: a, b
    character(len=*), intent(in)        :: string
    real                                :: real_a(size(a,1),size(a,2))
    real                                :: real_b(size(b,1),size(b,2))

    real_a = a 
    real_b = b
    call compare_data_r2d(real_a, real_b, string)

  end subroutine compare_data_i2d

  subroutine compare_data_r1d( a, b, string )
    real, intent(in), dimension(:) :: a, b
    character(len=*), intent(in)   :: string
    integer(LONG_KIND)             :: sum1, sum2
    integer                        :: l
    integer, parameter             :: stdunit = 6

    if(size(a,1) .ne. size(b,1) ) &
         call mpp_error(FATAL,'compare_data_r1d: size of a and b does not match')

    do l = 1, size(a(:))
       if(a(l) .ne. b(l)) then
          write(stdunit,'(a,i3,a,i3,a,f16.9,a,f16.9)')" at pe ", mpp_pe(), &
               ", at point (",l, "), a = ", a(l), ", b = ", b(l)
          call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
       endif
    enddo
    sum1 = mpp_chksum( a, (/mpp_pe()/) )
    sum2 = mpp_chksum( b, (/mpp_pe()/) )

    if( sum1.EQ.sum2 )then
       if( mpp_pe() .EQ. mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays 
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_data_r1d

  subroutine compare_data_i1d( a, b, string )
    integer, intent(in), dimension(:) :: a, b
    character(len=*), intent(in)      :: string
    real                              :: real_a(size(a(:)))
    real                              :: real_b(size(b(:)))

    real_a = a 
    real_b = b
    call compare_data_r1d(real_a, real_b, string)

  end subroutine compare_data_i1d

end program fms_io_test

#else
module null_fms_io_test
end module  

#endif  /* test_fms_io */
