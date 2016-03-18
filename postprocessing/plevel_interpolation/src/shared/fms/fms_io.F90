#include <fms_platform.h>

module fms_io_mod

!
!
! <CONTACT EMAIL="Zhi.Liang@noaa.gov">
! Zhi Liang
! </CONTACT>

! <CONTACT EMAIL="Matthew.Harrison@noaa.gov">
! M.J. Harrison 
! </CONTACT>
!
! <REVIEWER EMAIL="Matthew.Harrison@noaa.gov">
! M.J. Harrison 
! </REVIEWER>

! <REVIEWER EMAIL="Bruce.Wyman@noaa.gov">
! B. Wyman
! </REVIEWER>

!<DESCRIPTION>
! This module is for writing and reading restart data in NetCDF format.
! fms_io_init must be called before the first write_data/read_data call
! For writing, fms_io_exit must be called after ALL write calls have
! been made. Typically, fms_io_init and fms_io_exit are placed in the
! main (driver) program while read_data and write_data can be called where needed.
! Presently, two combinations of threading and fileset are supported, users can choose
! one line of the following by setting namelist:
!
! With the introduction of netCDF restart files, there is a need for a global
! switch to turn on/off netCDF restart options in all of the modules that deal with
! restart files. Here two more namelist variables (logical type) are introduced to fms_io
!
! fms_netcdf_override
! fms_netcdf_restart
!
! because default values of both flags are .true., the default behavior of the entire model is 
! to use netCDF IO mode. To turn off netCDF restart, simply set fms_netcdf_restart to .false.
!
! Fei.Liu@noaa.gov
! 05222006
! Read distributed files in NetCDF is available. Details can be found in read_data_3d_new
! <PRE>
!threading_read='multi', threading_write='multi', fileset_write='multi' (default)
!threading_read='multi', threading_write='single', fileset_write='single'
! </PRE>
!</DESCRIPTION>
! <NAMELIST NAME="fms_io_nml">  
! <DATA NAME="threading_read" TYPE="character">
! threading_read can be 'single' or 'multi'
! </DATA>
! <DATA NAME="threading_write" TYPE="character">
! threading_write can be 'single' or 'multi'
! </DATA>
! <DATA NAME="fileset_write" TYPE="character">
! fileset_write can be 'single' or 'multi'
! </DATA>
! <DATA NAME="fms_netcdf_override" TYPE="logical">
!   .true. : fms_netcdf_restart overrides individual do_netcdf_restart value (default behavior)
!   .false.: individual module settings has a precedence over the global setting, therefore fms_netcdf_restart is ignored
! </DATA>
! <DATA NAME="fms_netcdf_restart" TYPE="logical">
!   .true. : all modules deal with restart files will operate under netCDF mode (default behavior)
!   .false.: all modules deal with restart files will operate under binary mode
!   This flag is effective only when fms_netcdf_override is .true. When fms_netcdf_override is .false., individual
!   module setting takes over.
! </DATA>
! <DATA NAME="time_stamped_restart" TYPE="logical">
!   .true. : time_stamp will be added to the restart file name as a prefix when 
!            optional argument time_stamp is passed into routine save_restart. 
!   .false.: time_stmp will not be added to the restart file name even though
!            time_stamp is passed into save_restart.
!    default is true.
! </DATA>
! <DATA NAME="print_chksum" TYPE="logical">
!    set print_chksum (default is false) to true to print out chksum of fields that are
!    read and written through save_restart/restore_state. The chksum is accross all the 
!    processors, so there will be only one chksum even there are multiple-tiles in the 
!    grid. For the multiple case, the filename appeared in the message will contain
!    tile1 because the message is print out from root pe and on root pe the tile id is tile1. 
! </DATA>
!</NAMELIST>

use mpp_io_mod,      only: mpp_open, mpp_close, mpp_io_init, mpp_io_exit, mpp_read, mpp_write
use mpp_io_mod,      only: mpp_write_meta, mpp_get_info, mpp_get_atts, mpp_get_fields
use mpp_io_mod,      only: mpp_get_axes, mpp_get_axis_data, mpp_get_att_char, mpp_get_att_name
use mpp_io_mod,      only: mpp_get_att_real_scalar
use mpp_io_mod,      only: fieldtype, axistype, atttype, default_field, default_axis, default_att
use mpp_io_mod,      only: MPP_NETCDF, MPP_ASCII, MPP_MULTI, MPP_SINGLE, MPP_OVERWR, MPP_RDONLY
use mpp_io_mod,      only: MPP_IEEE32, MPP_NATIVE, MPP_DELETE, MPP_APPEND, MPP_SEQUENTIAL, MPP_DIRECT
use mpp_io_mod,      only: MAX_FILE_SIZE, mpp_get_att_value
use mpp_domains_mod, only: domain2d, domain1d, NULL_DOMAIN1D, NULL_DOMAIN2D, operator( .EQ. ), CENTER
use mpp_domains_mod, only: mpp_get_domain_components, mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only: mpp_get_domain_shift, mpp_get_global_domain, mpp_global_field, mpp_domain_is_tile_root_pe
use mpp_domains_mod, only: mpp_get_ntile_count, mpp_get_current_ntile, mpp_get_tile_id, mpp_mosaic_defined
use mpp_domains_mod, only: mpp_get_io_domain
use mpp_mod,         only: mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe, mpp_npes, stdlog, stdout
use mpp_mod,         only: mpp_broadcast, ALL_PES, mpp_chksum, mpp_get_current_pelist, mpp_npes, lowercase

use platform_mod, only: r8_kind

implicit none
private


integer, parameter, private :: max_split_file = 50
integer, parameter, private :: max_fields=400
integer, parameter, private :: max_axes=40
integer, parameter, private :: max_atts=20
integer, parameter, private :: max_domains = 10
integer, parameter, private :: MAX_TIME_LEVEL_REGISTER = 2
integer, parameter, private :: MAX_TIME_LEVEL_WRITE = 20
integer, parameter          :: max_axis_size=10000

type var_type
   private
   character(len=128)                     :: name
   character(len=128)                     :: longname
   character(len=128)                     :: units
   real, dimension(:,:,:,:), _ALLOCATABLE :: buffer _NULL 
   logical                                :: domain_present
   integer                                :: domain_idx
   logical                                :: is_dimvar
   type(fieldtype)                        :: field
   type(axistype)                         :: axis
   integer                                :: position
   integer                                :: ndim
   integer                                :: siz(4)      ! X/Y/Z/T extent of fields (data domain 
                                                         ! size for distributed writes;global size for reads)
   integer                                :: gsiz(4)     ! global X/Y/Z/T extent of fields 
   integer                                :: csiz(4)     ! actual data size in the file
   integer                                :: id_axes(3)  ! store index for x/y/z axistype.
   logical                                :: initialized ! indicate if the field is read or not in routine save_state.
   logical                                :: mandatory   ! indicate if the field is mandatory to be when restart.
   integer                                :: is, ie, js, je  ! index of the data in compute domain
   real                                   :: default_data 
end type var_type

type Ptr0Dr
   real,                   pointer :: p => NULL()
end type Ptr0Dr

type Ptr1Dr
   real, dimension(:),     pointer :: p => NULL()
end type Ptr1Dr

type Ptr2Dr
   real, dimension(:,:),   pointer :: p => NULL()
end type Ptr2Dr

type Ptr3Dr
   real, dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Dr

type Ptr0Di
   integer,                   pointer :: p => NULL()
end type Ptr0Di

type Ptr1Di
   integer, dimension(:),     pointer :: p => NULL()
end type Ptr1Di

type Ptr2Di
   integer, dimension(:,:),   pointer :: p => NULL()
end type Ptr2Di

type Ptr3Di
   integer, dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Di

type restart_file_type
   private
   integer                                  :: unit ! mpp_io unit for netcdf file
   character(len=128)                       :: name
   integer                                  :: nvar, natt, max_ntime
   logical                                  :: is_root_pe
   integer                                  :: tile_count
   type(var_type), dimension(:),   pointer  :: var  => NULL()
   type(Ptr0Dr),   dimension(:,:), pointer  :: p0dr => NULL()
   type(Ptr1Dr),   dimension(:,:), pointer  :: p1dr => NULL()
   type(Ptr2Dr),   dimension(:,:), pointer  :: p2dr => NULL()
   type(Ptr3Dr),   dimension(:,:), pointer  :: p3dr => NULL()
   type(Ptr0Di),   dimension(:,:), pointer  :: p0di => NULL()
   type(Ptr1Di),   dimension(:,:), pointer  :: p1di => NULL()
   type(Ptr2Di),   dimension(:,:), pointer  :: p2di => NULL()
   type(Ptr3Di),   dimension(:,:), pointer  :: p3di => NULL()
end type restart_file_type

interface read_data
   module procedure read_data_3d_new
   module procedure read_data_2d_new
   module procedure read_data_1d_new
   module procedure read_data_scalar_new
   module procedure read_data_i3d_new
   module procedure read_data_i2d_new
   module procedure read_data_i1d_new
   module procedure read_data_iscalar_new
   module procedure read_data_2d, read_ldata_2d, read_idata_2d
   module procedure read_data_3d, read_data_4d
#ifdef OVERLOAD_C8
   module procedure read_cdata_2d,read_cdata_3d,read_cdata_4d
#endif
   module procedure read_data_text
   module procedure read_data_2d_region
end interface

interface write_data
   module procedure write_data_3d_new
   module procedure write_data_2d_new
   module procedure write_data_1d_new
   module procedure write_data_scalar_new
   module procedure write_data_i3d_new
   module procedure write_data_i2d_new
   module procedure write_data_i1d_new
   module procedure write_data_iscalar_new
   module procedure write_data_2d, write_ldata_2d, write_idata_2d
   module procedure write_data_3d, write_data_4d
#ifdef OVERLOAD_C8
   module procedure write_cdata_2d,write_cdata_3d,write_cdata_4d
#endif
end interface

interface register_restart_field
   module procedure register_restart_field_r0d
   module procedure register_restart_field_r1d
   module procedure register_restart_field_r2d
   module procedure register_restart_field_r3d
   module procedure register_restart_field_i0d
   module procedure register_restart_field_i1d
   module procedure register_restart_field_i2d
   module procedure register_restart_field_i3d
   module procedure register_restart_field_r0d_2level
   module procedure register_restart_field_r1d_2level
   module procedure register_restart_field_r2d_2level
   module procedure register_restart_field_r3d_2level
   module procedure register_restart_field_i0d_2level
   module procedure register_restart_field_i1d_2level
   module procedure register_restart_field_i2d_2level
   module procedure register_restart_field_i3d_2level
end interface

interface reset_field_pointer
   module procedure reset_field_pointer_r0d
   module procedure reset_field_pointer_r1d
   module procedure reset_field_pointer_r2d
   module procedure reset_field_pointer_r3d
   module procedure reset_field_pointer_i0d
   module procedure reset_field_pointer_i1d
   module procedure reset_field_pointer_i2d
   module procedure reset_field_pointer_i3d
   module procedure reset_field_pointer_r0d_2level
   module procedure reset_field_pointer_r1d_2level
   module procedure reset_field_pointer_r2d_2level
   module procedure reset_field_pointer_r3d_2level
   module procedure reset_field_pointer_i0d_2level
   module procedure reset_field_pointer_i1d_2level
   module procedure reset_field_pointer_i2d_2level
   module procedure reset_field_pointer_i3d_2level
end interface

interface restore_state
   module procedure restore_state_all
   module procedure restore_state_one_field
end interface

interface query_initialized
   module procedure query_initialized_id
   module procedure query_initialized_name
   module procedure query_initialized_r2d
end interface

interface get_global_att_value
  module procedure get_global_att_value_text
  module procedure get_global_att_value_real
end interface

interface get_var_att_value
  module procedure get_var_att_value_text
end interface

integer :: num_files_r = 0 ! number of currently opened files for reading
integer :: num_files_w = 0 ! number of currently opened files for writing
integer :: num_domains = 0 ! number of domains in array_domain
integer :: num_registered_files ! mumber of files registered by calling register_restart_file


integer :: thread_r, thread_w, fset_w, form
logical :: module_is_initialized = .FALSE.

character(len=32) :: pelist_name
character(len=7)  :: pe_name
character(len=128):: error_msg  

  
!------ private data, pointer to current 2d domain ------
! entrained from fms_mod.  This will be deprecated in the future.
type(domain2D), pointer, private :: Current_domain =>NULL()

integer, private :: is,ie,js,je      ! compute domain
integer, private :: isd,ied,jsd,jed  ! data domain
integer, private :: isg,ieg,jsg,jeg  ! global domain
character(len=128),      dimension(:), allocatable         :: registered_file ! file names registered through register_restart_file 
type(restart_file_type), dimension(:), allocatable         :: files_read  ! store files that are read through read_data
type(restart_file_type), dimension(:), allocatable, target :: files_write ! store files that are written through write_data
type(domain2d), dimension(max_domains), save       :: array_domain
type(domain1d), dimension(max_domains), save       :: domain_x, domain_y
public  :: read_data, write_data, fms_io_init, fms_io_exit, field_size
public  :: open_namelist_file, open_restart_file, open_ieee32_file, close_file 
public  :: set_domain, nullify_domain, get_domain_decomp, return_domain
public  :: open_file, open_direct_file
public  :: get_restart_io_mode, get_tile_string, string
public  :: get_mosaic_tile_grid, get_mosaic_tile_file
public  :: get_global_att_value, get_var_att_value
public  :: file_exist, field_exist
public  :: register_restart_field, save_restart, restore_state
public  :: restart_file_type, query_initialized
public  :: reset_field_name, reset_field_pointer
private :: lookup_field_r, lookup_axis, unique_axes

public  :: set_filename_appendix, get_instance_filename
character(len=32), save :: filename_appendix = ''

!--- public interface ---
interface string
   module procedure string_from_integer
   module procedure string_from_real
end interface

!--- namelist interface
logical           :: fms_netcdf_override = .true.
logical           :: fms_netcdf_restart  = .true.
character(len=32) :: threading_read      = 'multi'
character(len=32) :: threading_write     = 'multi'
character(len=32) :: fileset_write       = 'multi'
character(len=32) :: format              = 'netcdf' 
logical           :: read_all_pe         = .TRUE.
character(len=64) :: iospec_ieee32       = '-N ieee_32'
integer           :: max_files_w         = 40
integer           :: max_files_r         = 40
logical           :: read_data_bug       = .false.
logical           :: time_stamp_restart  = .true.
logical           :: print_chksum        = .false.
  namelist /fms_io_nml/ fms_netcdf_override, fms_netcdf_restart, &
       threading_read, threading_write, &
       fileset_write, format, read_all_pe, iospec_ieee32,max_files_w,max_files_r, &
       read_data_bug, time_stamp_restart, print_chksum


character(len=128) :: version = '$Id: fms_io.F90,v 18.0 2010/03/02 23:55:56 fms Exp $'
character(len=128) :: tagname = '$Name: testing $'

contains

! <SUBROUTINE NAME="get_restart_io_mode">
! <DESCRIPTION>
! With the introduction of netCDF restart files, there is a need for a global
! switch to turn on/off netCDF restart options in all of the modules that deal with
! restart files. Here two more namelist variables (logical type) are introduced to fms_io
!
! fms_netcdf_override
! fms_netcdf_restart
!
! because default values of both flags are .true., the default behavior of the entire model is 
! to use netCDF IO mode. To turn off netCDF restart, simply set fms_netcdf_restart to .false.
! 
! </DESCRIPTION>
! <TEMPLATE>
!  call get_fms_io_mode(do_netcdf_restart)
! </TEMPLATE>
! <INOUT NAME="do_netcdf_restart" TYPE="logical">
!  This the input argument that contains the individual module setting of restart IO mode.
!  Upon return from this subroutine, this output argument contains the actual setting of restart IO mode
!  the calling module will be using
! </INOUT>
! </SUBROUTINE>
subroutine get_restart_io_mode(do_netcdf_restart)

  logical, intent(inout)  :: do_netcdf_restart

  if(fms_netcdf_override) do_netcdf_restart = fms_netcdf_restart
  
end subroutine get_restart_io_mode
!.....................................................................
! <SUBROUTINE NAME="fms_io_init">
!   <DESCRIPTION>
! Initialize fms_io module
!   </DESCRIPTION>
!   <TEMPLATE>
! call fms_io_init()
!   </TEMPLATE>
subroutine fms_io_init()
    
  integer                            :: i, unit, io_status, logunit
  integer, allocatable, dimension(:) :: pelist


  if (module_is_initialized) return
  call mpp_io_init()

  call mpp_open(unit, 'input.nml',form=MPP_ASCII,action=MPP_RDONLY)
  read(unit,fms_io_nml,iostat=io_status)
  if (io_status > 0) then
     call mpp_error(FATAL,'=>fms_io_init: Error reading input.nml')
  endif
  call mpp_close (unit)
  if (mpp_pe() == mpp_root_pe()) then
    logunit = stdlog() ; write(logunit, fms_io_nml)
    write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)
  end if
! take namelist options if present

  select case (threading_read) 
  case ('multi')
     thread_r = MPP_MULTI
  case ('single')
     thread_r = MPP_SINGLE
  case default
     call mpp_error(FATAL,'fms_io_init: threading_read should be multi/single but you chose'//trim(threading_read))
  end select
! take namelist options if present

  select case (fileset_write) 
  case ('multi')
     fset_w = MPP_MULTI
  case ('single')
     fset_w = MPP_SINGLE
  case default
     call mpp_error(FATAL,'fms_io_init: fileset_write should be multi/single but you chose'//trim(fileset_write))
  end select

  select case (threading_write) 
  case ('multi')
     thread_w = MPP_MULTI
  case ('single')
     thread_w = MPP_SINGLE
  case default
     call mpp_error(FATAL,'fms_io_init: threading_write should be multi/single but you chose'//trim(threading_write))
  end select

  select case(format)
  case ('netcdf')
     form=MPP_NETCDF
  case default
     call mpp_error(FATAL,'fms_io_init: only NetCDF format currently supported in fms_io')
  end select

! Initially allocate  files_write and files_read
  allocate(files_write(max_files_w),files_read(max_files_r))
  allocate(registered_file(max_files_w))

  allocate(pelist(mpp_npes()))        
  call mpp_get_current_pelist(pelist,pelist_name)
  if(mpp_npes()>10000) then
     write(pe_name,'(a,i6.6)' )'.', mpp_pe()    
  else
     write(pe_name,'(a,i4.4)' )'.', mpp_pe()    
  endif
  deallocate(pelist)

  do i = 1, max_domains
     array_domain(i) = NULL_DOMAIN2D
  enddo
  !---- initialize module domain2d pointer ----
  nullify (Current_domain)
  module_is_initialized = .TRUE.
  
end subroutine fms_io_init

! </SUBROUTINE>
! <SUBROUTINE NAME="fms_io_exit">
!   <DESCRIPTION>
! This routine is called after ALL fields have been written to temporary files
! The result NETCDF files are created here.
!   </DESCRIPTION>
!   <TEMPLATE>
! call fms_io_exit
!   </TEMPLATE>

subroutine fms_io_exit()
    integer                             :: num_x_axes, num_y_axes, num_z_axes
    integer                             :: unit
    real, dimension(max_axis_size)      :: axisdata
    real(r8_kind)                       :: tlev  
    integer,        dimension(max_axes) :: id_x_axes, siz_x_axes
    integer,        dimension(max_axes) :: id_y_axes, siz_y_axes
    integer,        dimension(max_axes) :: id_z_axes, siz_z_axes
    type(axistype), dimension(max_axes) :: x_axes, y_axes, z_axes
    type(axistype)                      :: t_axes
    type(var_type), pointer, save       :: cur_var=>NULL()
    integer                             :: i, j, k, kk
    character(len=256)                  :: filename
    character(len=10)                   :: axisname
    logical                             :: domain_present

    if( .NOT.module_is_initialized )return !make sure it's only called once per PE

    do i=1,max_axis_size
       axisdata(i) = i
    enddo

    ! each field has an associated domain type (may be undefined).
    ! each file only needs to write unique axes (i.e. if 2 fields share an identical axis, then only write the axis once)
    ! unique axes are defined by the global size and domain decomposition (i.e. can support identical axis sizes with
    ! different domain decomposition)

    do i = 1, num_files_w
       filename = files_write(i)%name

       !--- check if any field in this file present domain.
       domain_present = .false.
       do j = 1, files_write(i)%nvar
          if (files_write(i)%var(j)%domain_present) then
              domain_present = .true.
              exit
          end if
       end do

       !--- get the unique axes for all the fields.
       num_x_axes = unique_axes(files_write(i), 1, id_x_axes, siz_x_axes, domain_x)
       num_y_axes = unique_axes(files_write(i), 2, id_y_axes, siz_y_axes, domain_y)
       num_z_axes = unique_axes(files_write(i), 3, id_z_axes, siz_z_axes          )

       if( domain_present ) then
          call mpp_open(unit,trim(filename),action=MPP_OVERWR,form=form,threading=thread_w,&
               fileset=fset_w, is_root_pe=files_write(i)%is_root_pe, domain=array_domain(files_write(i)%var(j)%domain_idx))
       else  ! global data
          call mpp_open(unit,trim(filename),action=MPP_OVERWR,form=form,threading=MPP_SINGLE,&
               fileset=MPP_SINGLE, is_root_pe=files_write(i)%is_root_pe)
       end if

       do j = 1, num_x_axes
         if (j < 10) then
             write(axisname,'(a,i1)') 'xaxis_',j
          else
             write(axisname,'(a,i2)') 'xaxis_',j
          endif          
          if(id_x_axes(j) > 0) then
             call mpp_write_meta(unit,x_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_x_axes(j)),domain=domain_x(id_x_axes(j)),cartesian='X')
          else
             call mpp_write_meta(unit,x_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_x_axes(j)),cartesian='X')
          endif             
       end do

       do j = 1, num_y_axes
         if (j < 10) then
             write(axisname,'(a,i1)') 'yaxis_',j
          else
             write(axisname,'(a,i2)') 'yaxis_',j
          endif          
          if(id_y_axes(j) > 0) then
             call mpp_write_meta(unit,y_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_y_axes(j)),domain=domain_y(id_y_axes(j)),cartesian='Y')
          else
             call mpp_write_meta(unit,y_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_y_axes(j)),cartesian='Y')
          endif             
       end do

       do j = 1, num_z_axes
          if (j < 10) then
             write(axisname,'(a,i1)') 'zaxis_',j
          else
             write(axisname,'(a,i2)') 'zaxis_',j
          endif          
          call mpp_write_meta(unit,z_axes(j),axisname,'none',axisname, &
               data=axisdata(1:siz_z_axes(j)),cartesian='Z')
       end do


       ! write time axis  (comment out if no time axis)
       call mpp_write_meta(unit,t_axes,&
            'Time','time level','Time',cartesian='T')

       ! write metadata for fields
       do j = 1, files_write(i)%nvar
          cur_var => files_write(i)%var(j)
          call mpp_write_meta(unit,cur_var%field, (/x_axes(cur_var%id_axes(1)), &
               y_axes(cur_var%id_axes(2)), z_axes(cur_var%id_axes(3)), t_axes/), cur_var%name, &
               'none',cur_var%name,pack=1)
       enddo

       ! write values for ndim of spatial axes
       do j = 1, num_x_axes
          call mpp_write(unit,x_axes(j))
       enddo
       do j = 1, num_y_axes
          call mpp_write(unit,y_axes(j))
       enddo
       do j = 1, num_z_axes
          call mpp_write(unit,z_axes(j))
       enddo

       ! write data of each field
       do k = 1, files_write(i)%max_ntime
          do j = 1, files_write(i)%nvar
             cur_var => files_write(i)%var(j)
             tlev=k
             ! If some fields only have one time level, we do not need to write the second level, just keep
             ! the data missing.
             ! If some fields only have one time level, we just write out 0 to the other level
             if(k > cur_var%siz(4)) then
                cur_var%buffer(:,:,:,1) = 0.0
                kk = 1
             else
                kk = k
             end if
             if(cur_var%domain_present) then
                call mpp_write(unit, cur_var%field,array_domain(cur_var%domain_idx), cur_var%buffer(:,:,:,kk), tlev, &
                               default_data=cur_var%default_data)
             else if (thread_w == MPP_MULTI .or. (files_write(i)%is_root_pe.and.thread_w == MPP_SINGLE)) then
                call mpp_write(unit, cur_var%field, cur_var%buffer(:,:,:,kk), tlev)
             end if
          enddo ! end j loop
       enddo ! end k loop
       call mpp_close(unit)
    enddo ! end i loop

    !--- release the memory 

    do i = 1,  num_files_w
       do j = 1, files_write(i)%nvar
          deallocate(files_write(i)%var(j)%buffer)
       end do
    end do

  cur_var=>NULL()
  module_is_initialized = .false.
  num_files_w = 0
  num_files_r = 0    

end subroutine fms_io_exit
!.....................................................................
! </SUBROUTINE>

! <SUBROUTINE NAME="write_data">
    !<DESCRIPTION>
    ! This subroutine performs writing "fieldname" to file "filename". All values of "fieldname" 
    ! will be written to a temporary file. The final NETCDF file will be created only at a later step
    ! when the user calls fms_io_exit. Therefore, make sure that fms_io_exit is called after all
    ! fields have been written by this subroutine.
    !</DESCRIPTION>
!   <TEMPLATE>
! call write_data(filename, fieldname, data, domain)
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character" DIM="(*)">
!    File name
!   </IN>
!   <IN NAME="fieldname" TYPE="character" DIM="(*)">
!    Field  name
!   </IN>
!   <IN NAME="data"  TYPE="real">
!   array containing data of fieldname
!   </IN>
!   <IN NAME="domain"  TYPE="domain, optional">
!   domain of fieldname
!   </IN>
!=================================================================================
subroutine write_data_i3d_new(filename, fieldname, data, domain,                  &
                              no_domain, position, tile_count, data_default)

  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:,:,:), intent(in) :: data
  type(domain2d), intent(in), optional :: domain
  logical, intent(in), optional :: no_domain
  integer, intent(in), optional :: position, tile_count, data_default
  real :: default_data

  default_data = 0
  if(present(data_default)) default_data = real(data_default)

  call write_data_3d_new(filename, fieldname, real(data), domain,  &
                         no_domain, position, tile_count, data_default=default_data)
end subroutine write_data_i3d_new
!.....................................................................
subroutine write_data_i2d_new(filename, fieldname, data, domain, &
                              no_domain, position, tile_count, data_default)

  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:,:), intent(in) :: data
  type(domain2d), intent(in), optional :: domain
  logical, intent(in), optional :: no_domain
  integer, intent(in), optional :: position, tile_count, data_default
  real :: default_data

  default_data = 0
  if(present(data_default)) default_data = real(data_default)
  call write_data_2d_new(filename, fieldname, real(data), domain, &
                         no_domain, position, tile_count, data_default=default_data)

end subroutine write_data_i2d_new
!.....................................................................
subroutine write_data_i1d_new(filename, fieldname, data, domain, &
                              no_domain, tile_count, data_default)
  type(domain2d), intent(in), optional :: domain
  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:), intent(in) :: data
  logical, intent(in), optional :: no_domain
  integer, intent(in), optional :: tile_count, data_default
  real :: default_data

  default_data = 0
  if(present(data_default)) default_data = real(data_default)
  call write_data_1d_new(filename, fieldname, real(data), domain, &
                         no_domain, tile_count, data_default=default_data)
end subroutine write_data_i1d_new
!.....................................................................
subroutine write_data_iscalar_new(filename, fieldname, data, domain, &
                                  no_domain, tile_count, data_default)
  type(domain2d), intent(in), optional :: domain
  character(len=*), intent(in) :: filename, fieldname 
  integer, intent(in) :: data
  logical, intent(in), optional :: no_domain
  integer, intent(in), optional :: tile_count, data_default
  real :: default_data

  default_data = 0
  if(present(data_default)) default_data = real(data_default)
  call write_data_scalar_new(filename, fieldname, real(data), domain, &
                             no_domain, tile_count, data_default=default_data)

end subroutine write_data_iscalar_new
!.....................................................................
subroutine write_data_3d_new(filename, fieldname, data, domain, no_domain, &
                             position, tile_count, data_default)

  character(len=*),         intent(in)         :: filename, fieldname 
  real, dimension(:,:,:),   intent(in)         :: data
  type(domain2d), optional, intent(in), target :: domain
  real,           optional, intent(in)         :: data_default
  logical,        optional, intent(in)         :: no_domain   
  integer,        optional, intent(in)         :: position, tile_count

  !--- local variables
  real,               allocatable :: tmp_buffer(:,:,:,:)
  integer                         :: index_field ! position of the fieldname in the list of fields
  integer                         :: index_file  ! position of the filename in the list of files_write
  logical                         :: append_pelist, is_no_domain
  character(len=256)              :: fname, filename2,append_string
  real                            :: default_data
  integer                         :: length, i, domain_idx
  integer                         :: ishift, jshift
  integer                         :: gxsize, gysize
  integer                         :: cxsize, cysize
  integer                         :: dxsize, dysize
  type(domain2d), pointer, save   :: d_ptr   =>NULL()
  type(var_type), pointer, save   :: cur_var =>NULL()
  type(restart_file_type), pointer, save :: cur_file =>NULL()

! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_3d_new): need to call fms_io_init')  


  if(PRESENT(data_default))then
     default_data=data_default
  else
     default_data=0.
  endif

  if(present(tile_count) .AND. .not. present(domain)) call mpp_error(FATAL, &
         'fms_io write_data: when tile_count is present, domain must be present')

  is_no_domain = .false.
  if (PRESENT(no_domain)) THEN
     is_no_domain = no_domain
  end if

  if(is_no_domain) then
!     comments the following to allow passing domain when writing 1-D or scalar variables. 
!     if(PRESENT(domain)) &
!       call mpp_error(FATAL, 'fms_io(write_data_3d_new): no_domain cannot be .true. when optional argument domain is present.')  
  else if(PRESENT(domain))then
     d_ptr => domain
  else if (ASSOCIATED(Current_domain) .and. .NOT. is_no_domain ) then
     d_ptr => Current_domain
  endif

  !--- remove .nc from file name
  length = len_trim(filename)
  if(filename(length-2:length) == '.nc') then
     filename2 = filename(1:length-3)
  else
     filename2 = filename(1:length)
  end if

  !Logical append_pelist decides whether to append the pelist_name to file name
  append_pelist = .false.
  !Append a string to the file name
  append_string=''

  !If the filename_appendix  is set override the passed argument. 
  if(len_trim(filename_appendix) > 0)  then
     append_pelist = .true.
     append_string = filename_appendix
  endif

  if(append_pelist) filename2 = trim(filename2)//'.'//trim(append_string)

  !JWD:  This is likely a temporary fix. Since fms_io needs to know tile_count,
  !JWD:  I just don't see how the physics can remain "tile neutral"
  !z1l:  one solution is add one more public interface called set_tile_count
  call get_mosaic_tile_file(filename2, fname, is_no_domain, domain, tile_count)

  ! Check if filename has been open  or not
  index_file = -1
  do i=1,num_files_w
     if (trim(files_write(i)%name) == trim(fname)) then
        index_file = i
        cur_file => files_write(index_file)
        exit
     endif
  enddo

  if (index_file < 0) then 
     if(num_files_w == max_files_w) &  ! need to have bigger max_files_w
          call mpp_error(FATAL,'fms_io(write_data_3d_new): max_files_w exceeded, increase it via fms_io_nml')    
     ! record the file name in array files_write
     num_files_w=num_files_w + 1
     index_file = num_files_w
     cur_file => files_write(index_file)
     cur_file%name = trim(fname)         
     cur_file%tile_count=1
     if(present(tile_count)) cur_file%tile_count = tile_count
     if(ASSOCIATED(d_ptr))then
        cur_file%is_root_pe = mpp_domain_is_tile_root_pe(d_ptr)
     else
        cur_file%is_root_pe = mpp_pe() == mpp_root_pe()
     endif
     cur_file%max_ntime = 1
     !-- allocate memory
     allocate(cur_file%var(max_fields) )
     cur_file%nvar = 0
     do i = 1, max_fields
        cur_file%var(i)%name           = 'none'
        cur_file%var(i)%domain_present = .false.
        cur_file%var(i)%domain_idx     = -1
        cur_file%var(i)%is_dimvar      = .false.
        cur_file%var(i)%position       = CENTER
        cur_file%var(i)%siz(:)         = 0
        cur_file%var(i)%gsiz(:)        = 0
        cur_file%var(i)%id_axes(:)     = -1
     end do
  endif

  ! check if the field is new or not and get position and dimension of the field
  index_field = -1
  do i = 1, cur_file%nvar
     if(trim(cur_file%var(i)%name) == trim(fieldname)) then
        index_field = i
        exit
     end if 
  end do

  if(index_field > 0) then
     cur_var   => cur_file%var(index_field)
     cur_var%siz(4) =  cur_var%siz(4) + 1
     if(cur_file%max_ntime < cur_var%siz(4) ) cur_file%max_ntime = cur_var%siz(4)
     ! the time level should be no larger than MAX_TIME_LEVEL_WRITE ( =20) for write_data.
     if( cur_var%siz(4) > MAX_TIME_LEVEL_WRITE ) call mpp_error(FATAL, 'fms_io(write_data_3d_new): ' // &
          'the time level of field '//trim(cur_var%name)//' in file '//trim(cur_file%name)// &
          ' is greater than MAX_TIME_LEVEL_WRITE(=20), increase MAX_TIME_LEVEL_WRITE or check your code')
  else 
     cur_file%nvar = cur_file%nvar +1
     if(cur_file%nvar>max_fields) then
        write(error_msg,'(I3,"/",I3)') cur_file%nvar, max_fields 
        call  mpp_error(FATAL,'fms_io(write_data_3d_new): max_fields exceeded, needs increasing, nvar/max_fields=' &
             //trim(error_msg))
     endif
     index_field =  cur_file%nvar
     cur_var   => cur_file%var(index_field)
     cur_var%siz(1)  = size(data,1)
     cur_var%siz(2)  = size(data,2)
     cur_var%siz(3)  = size(data,3)
     cur_var%siz(4)  = 1
     cur_var%gsiz(3) = cur_var%siz(3)
     cur_var%name = fieldname
     cur_var%default_data = default_data
     cur_var%ndim = 3
     if(present(position)) cur_var%position = position
     if(ASSOCIATED(d_ptr)) then
        cur_var%domain_present = .true.
        domain_idx = lookup_domain(d_ptr)
        if(domain_idx == -1) then
           num_domains = num_domains + 1
           if(num_domains > max_domains) call  mpp_error(FATAL,'fms_io(write_data_3d_new), 1: max_domains exceeded,' &
                //' needs increasing')
           domain_idx = num_domains
           array_domain(domain_idx) = d_ptr
           call mpp_get_domain_components(array_domain(domain_idx), domain_x(domain_idx), domain_y(domain_idx), &
                tile_count=tile_count)
        endif
        cur_var%domain_idx = domain_idx
        call mpp_get_domain_shift ( array_domain(domain_idx), ishift, jshift, position)
        call mpp_get_global_domain(array_domain(domain_idx), xsize=gxsize,ysize=gysize,tile_count=tile_count)
        call mpp_get_compute_domain(array_domain(domain_idx), xsize = cxsize, ysize = cysize, tile_count=tile_count)
        call mpp_get_data_domain   (array_domain(domain_idx), xsize = dxsize, ysize = dysize, tile_count=tile_count)
        if (ishift .NE. 0) then
           cxsize = cxsize+ishift; dxsize = dxsize+ishift; gxsize = gxsize + ishift
        end if
        if (jshift .NE. 0) then
           cysize = cysize+jshift; dysize = dysize+jshift; gysize = gysize + jshift
        endif
        if( (cur_var%siz(1) .NE. cxsize .AND. cur_var%siz(1) .NE. dxsize ) .OR. &
            (cur_var%siz(2) .NE. cysize .AND. cur_var%siz(2) .NE. dysize ) ) then
            call mpp_error(FATAL, 'fms_io(write_data_3d_new): data should be on either computer domain '//&
              'or data domain when domain is present for field '//trim(fieldname)//' of file '//trim(filename) )
        end if
        cur_var%gsiz(1)   = gxsize
        cur_var%gsiz(2)   = gysize
     else
        cur_var%gsiz(1) = size(data,1)
        cur_var%gsiz(2) = size(data,2)
        cur_var%domain_present=.false.
     endif
  end if

  ! copy the data to the buffer 
  ! if the time level is greater than the size(cur_var%buffer,4), 
  ! need to increase the buffer size 

  if(cur_var%siz(4) == 1) then
     allocate(cur_var%buffer(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3), cur_var%siz(4)) )
  else
     allocate(tmp_buffer(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3), size(cur_var%buffer,4)) )
     tmp_buffer = cur_var%buffer
     deallocate(cur_var%buffer)
     allocate(cur_var%buffer(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3), cur_var%siz(4)) )
     cur_var%buffer(:,:,:,1:size(tmp_buffer,4)) = tmp_buffer 
     deallocate(tmp_buffer)
  endif

  cur_var%buffer(:,:,:,cur_var%siz(4)) = data ! copy current data to buffer for future write out 

  d_ptr =>NULL()
  cur_var =>NULL()
  cur_file =>NULL()

end subroutine write_data_3d_new
! </SUBROUTINE>  

!-------------------------------------------------------------------------------
!
!   The routine will register a scalar real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r0d(fileObj, filename, fieldname, data, domain, mandatory, &
                                      position, tile_count, data_default, longname, units)
  type(restart_file_type),    intent(inout)      :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,                       intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: mandatory
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r0d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_r0d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/1, 1, 1, 1/), index_field, domain,      &
                       mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p0dr(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 0
  register_restart_field_r0d = index_field  

  return    

end function register_restart_field_r0d

!-------------------------------------------------------------------------------
!
!   The routine will register a 1-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r1d(fileObj, filename, fieldname, data, domain, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real, dimension(:),         intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r1d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_r1d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), 1, 1, 1/), index_field, domain, &
                       mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units )

  fileObj%p1dr(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 1
  register_restart_field_r1d = index_field  

  return    
  
end function register_restart_field_r1d

!-------------------------------------------------------------------------------
!
!   The routine will register a 2-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r2d(fileObj, filename, fieldname, data, domain, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:,:),   intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r2d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_r2d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), size(data,2), 1, 1/), &
                       index_field, domain, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p2dr(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 2
  register_restart_field_r2d = index_field    

  return    

end function register_restart_field_r2d


!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r3d(fileObj, filename, fieldname, data, domain, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:,:,:), intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r3d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_r3d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), size(data,2), size(data,3), 1/), &
                       index_field, domain, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p3dr(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 3
  register_restart_field_r3d = index_field   

  return    

end function register_restart_field_r3d

!-------------------------------------------------------------------------------
!
!   The routine will register a scalar integer restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i0d(fileObj, filename, fieldname, data, domain, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,                    intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i0d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_i0d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/1, 1, 1, 1/), index_field, domain, &
                       mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p0di(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 0
  register_restart_field_i0d = index_field
  
  return
   
end function register_restart_field_i0d

!-------------------------------------------------------------------------------
!
!   The routine will register a 1-D integer restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i1d(fileObj, filename, fieldname, data, domain, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer, dimension(:),      intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i1d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_i1d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), 1, 1, 1/), index_field, domain, &
                       mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p1di(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 1
  register_restart_field_i1d = index_field
  
  return
  
end function register_restart_field_i1d


!-------------------------------------------------------------------------------
!
!   The routine will register a 2-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i2d(fileObj, filename, fieldname, data, domain, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:,:),   intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i2d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_i2d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), size(data,2), 1, 1/), &
                       index_field, domain, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p2di(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 2
  register_restart_field_i2d = index_field
  
  return 

end function register_restart_field_i2d

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i3d(fileObj, filename, fieldname, data, domain, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:,:,:), intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i3d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_i3d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), size(data,2), size(data,3), 1/), &
                       index_field, domain, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p3di(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 3
  register_restart_field_i3d = index_field
  
  return  

end function register_restart_field_i3d

!-------------------------------------------------------------------------------
!
!   The routine will register a scalar real restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r0d_2level(fileObj, filename, fieldname, data1, data2, domain, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,                       intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r0d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_r0d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/1, 1, 1, 2/), index_field, domain, &
                       mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p0dr(1, index_field)%p => data1
  fileObj%p0dr(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 0
  register_restart_field_r0d_2level = index_field
  
  return  

end function register_restart_field_r0d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 1-D real restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r1d_2level(fileObj, filename, fieldname, data1, data2, domain, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:),     intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r1d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_r1d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), 1, 1, 2/), index_field, domain, &
                       mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p1dr(1, index_field)%p => data1
  fileObj%p1dr(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 1
  register_restart_field_r1d_2level = index_field
  
  return    

end function register_restart_field_r1d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D real restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r2d_2level(fileObj, filename, fieldname, data1, data2, domain, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:,:),   intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r2d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_r2d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), size(data1,2), 1, 2/), &
                       index_field, domain, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p2dr(1, index_field)%p => data1
  fileObj%p2dr(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 2
  register_restart_field_r2d_2level = index_field
  
  return    

end function register_restart_field_r2d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D real restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r3d_2level(fileObj, filename, fieldname, data1, data2, domain, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:,:,:), intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r3d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_r3d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), size(data1,2), size(data1,3), 2/), &
                       index_field, domain, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p3dr(1, index_field)%p => data1
  fileObj%p3dr(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 3
  register_restart_field_r3d_2level = index_field
  
  return    

end function register_restart_field_r3d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a scalar integer restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i0d_2level(fileObj, filename, fieldname, data1, data2, domain, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,                    intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i0d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_i0d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/1, 1, 1, 2/), index_field, domain, &
                       mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p0di(1, index_field)%p => data1
  fileObj%p0di(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 0
  register_restart_field_i0d_2level = index_field
  
  return   
 
end function register_restart_field_i0d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 1-D integer restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i1d_2level(fileObj, filename, fieldname, data1, data2, domain, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:),     intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i1d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_i1d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), 1, 1, 2/), index_field, domain, &
                       mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p1di(1, index_field)%p => data1
  fileObj%p1di(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 1
  register_restart_field_i1d_2level = index_field
  
  return  
  
end function register_restart_field_i1d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D integer restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i2d_2level(fileObj, filename, fieldname, data1, data2, domain, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:,:),   intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i2d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_i2d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), size(data1,2), 1, 2/), &
                       index_field, domain, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p2di(1, index_field)%p => data1
  fileObj%p2di(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 2
  register_restart_field_i2d_2level = index_field
  
  return    

end function register_restart_field_i2d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D integer restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i3d_2level(fileObj, filename, fieldname, data1, data2, domain, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:,:,:), intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i3d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_i3d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), size(data1,2), size(data1,3), 2/), &
                       index_field, domain, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p3di(1, index_field)%p => data1
  fileObj%p3di(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 3
  register_restart_field_i3d_2level = index_field
  
  return    

end function register_restart_field_i3d_2level

!-------------------------------------------------------------------------------
!
!  saves all registered variables to restart files. Those variables are set 
!  through register_restart_field
!
!-------------------------------------------------------------------------------

subroutine save_restart(fileObj, time_stamp, directory )
  type(restart_file_type), intent(inout) :: fileObj
  character(len=*), intent(in), optional :: directory
  character(len=*), intent(in), optional :: time_stamp
  ! Arguments: 
  !  (in)      directory  - The directory where the restart file goes.
  !  (in)      time_stamp - character format of the time of this restart file.
  character(len=256) :: dir
  character(len=256) :: restartpath          ! The restart file path (dir/file).
  character(len=80)  :: restartname          ! The restart file name (no dir).
  character(len=8)   :: suffix               ! A suffix (like _2) that is appended to the name of files after the first.
  integer            :: var_sz, size_in_file ! The size in bytes of each variable and of the variables already in a file.
  integer            :: start_var, next_var  ! The starting variables of the current and next files.
  integer            :: unit                 ! The mpp unit of the open file.
  real, dimension(max_axis_size)      :: axisdata
  integer,        dimension(max_axes) :: id_x_axes, siz_x_axes
  integer,        dimension(max_axes) :: id_y_axes, siz_y_axes
  integer,        dimension(max_axes) :: id_z_axes, siz_z_axes
  type(axistype), dimension(max_axes) :: x_axes, y_axes, z_axes
  type(axistype)                      :: t_axes            
  integer                             :: num_var_axes
  type(axistype), dimension(4)        :: var_axes
  type(var_type), pointer, save       :: cur_var=>NULL()
  integer                             :: num_x_axes, num_y_axes, num_z_axes
  integer                             :: naxes_x, naxes_y, naxes_z
  integer                             :: nfiles, i, j, k, l, siz, ind_dom
  logical                             :: domain_present
  real(r8_kind)                       :: tlev  
  character(len=10)                   :: axisname  

  real, allocatable, dimension(:,:,:) :: r3d
  real, allocatable, dimension(:,:)   :: r2d, global_r2d
  real, allocatable, dimension(:)     :: r1d  
  real                                :: r0d

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(save_restart): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  do i=1,max_axis_size
     axisdata(i) = i
  enddo

  dir = "RESTART"
  if(present(directory)) dir = directory

  restartname = fileObj%name
  nfiles = 0
  if(time_stamp_restart) then
     if (PRESENT(time_stamp)) then
        restartname = trim(time_stamp)//"."//trim(restartname)
     endif
  end if
  if(len_trim(dir) > 0) then
     restartpath = trim(dir)//"/"// trim(restartname)
  else
     restartpath = trim(restartname)
  end if
  !--- check if any field in this file present domain.
  domain_present = .false.
  do j = 1, fileObj%nvar
     if (fileObj%var(j)%domain_present) then
        domain_present = .true.
        ind_dom = j
        exit
     end if
  end do
  num_x_axes = unique_axes(fileObj, 1, id_x_axes, siz_x_axes, domain_x)
  num_y_axes = unique_axes(fileObj, 2, id_y_axes, siz_y_axes, domain_y)
  num_z_axes = unique_axes(fileObj, 3, id_z_axes, siz_z_axes          )
  next_var = 1
  size_in_file = 0
  do j = 1, num_x_axes
     size_in_file = size_in_file + siz_x_axes(j)
  end do
  do j = 1, num_y_axes
     size_in_file = size_in_file + siz_y_axes(j)
  end do
  do j = 1, num_z_axes
     size_in_file = size_in_file + siz_z_axes(j)
  end do
  size_in_file = 8*(size_in_file*2+1000)

  do while (next_var <= fileObj%nvar )
     start_var = next_var

     do j=start_var,fileObj%nvar
        cur_var => fileObj%var(j)
        var_sz = 8*cur_var%csiz(1)*cur_var%csiz(2)*cur_var%csiz(3)
        if ((j==start_var) .OR. (size_in_file < MAX_FILE_SIZE-var_sz)) then
           size_in_file = size_in_file + var_sz
        else 
           exit
        endif
     enddo
     next_var = j
     ! For distribute write, normally will not over the limit. 
     if( nfiles > 0 ) then
        if(fset_w == MPP_MULTI .AND. domain_present) call mpp_error(FATAL, "fms_io_mod(save_restart): "// &
             "For distribute write(fileset_write='multi'), the file size should not be very large and need to be split")
        if (nfiles < 10) then
           write(suffix,'("_",I1)') nfiles
        else if(nfiles < 100) then
           write(suffix,'("_",I2)') nfiles
        else
           call mpp_error(FATAL, "fms_io(save_restart): num_files should be less than 100")
        endif
        !--- remove .nc from restartpath and attach suffix.
        siz = len_trim(restartpath)
        if(restartpath(siz-2:siz) == ".nc") then
           restartpath = restartpath(1:siz-3)//trim(suffix)
        else      
           restartpath = trim(restartpath) // trim(suffix)
        end if
     end if
     if( domain_present ) then
        call mpp_open(unit,trim(restartpath),action=MPP_OVERWR,form=form,threading=thread_w,&
             fileset=fset_w, is_root_pe=fileObj%is_root_pe, domain=array_domain(fileObj%var(ind_dom)%domain_idx) )
     else  ! global data
        call mpp_open(unit,trim(restartpath),action=MPP_OVERWR,form=form,threading=MPP_SINGLE,&
             fileset=MPP_SINGLE, is_root_pe=fileObj%is_root_pe)
     end if

     ! write_out x_axes
     naxes_x = 0
     do j = 1, num_x_axes
        ! make sure this axis is used by some variable 
        do l=start_var,next_var-1
           if( fileObj%var(l)%id_axes(1) == j ) exit
        end do
        if(l == next_var) cycle  
        naxes_x = naxes_x + 1
        if (naxes_x < 10) then
           write(axisname,'(a,i1)') 'xaxis_',naxes_x
        else
           write(axisname,'(a,i2)') 'xaxis_',naxes_x
        endif
        if(id_x_axes(j) > 0) then
           call mpp_write_meta(unit,x_axes(naxes_x),axisname,'none',axisname, &
                data=axisdata(1:siz_x_axes(j)),domain=domain_x(id_x_axes(j)),cartesian='X')
        else
           call mpp_write_meta(unit,x_axes(naxes_x),axisname,'none',axisname, &
                data=axisdata(1:siz_x_axes(j)),cartesian='X')
        endif
     end do

     ! write out y_axes
     naxes_y = 0
     do j = 1, num_y_axes
        ! make sure this axis is used by some variable 
        do l=start_var,next_var-1
           if( fileObj%var(l)%id_axes(2) == j ) exit
        end do
        if(l == next_var) cycle  
        naxes_y = naxes_y + 1
        if (naxes_y < 10) then
           write(axisname,'(a,i1)') 'yaxis_',naxes_y
        else
           write(axisname,'(a,i2)') 'yaxis_',naxes_y
        endif
        if(id_y_axes(j) > 0) then
           call mpp_write_meta(unit,y_axes(naxes_y),axisname,'none',axisname, &
                data=axisdata(1:siz_y_axes(j)),domain=domain_y(id_y_axes(j)),cartesian='Y')
        else
           call mpp_write_meta(unit,y_axes(naxes_y),axisname,'none',axisname, &
                data=axisdata(1:siz_y_axes(j)),cartesian='Y')
        endif
     end do

     ! write out z_axes
     naxes_z = 0
     do j = 1, num_z_axes
        ! make sure this axis is used by some variable 
        do l=start_var,next_var-1
           if( fileObj%var(l)%id_axes(3) == j ) exit
        end do
        if(l == next_var) cycle  
        naxes_z = naxes_z + 1
        if (naxes_z < 10) then
           write(axisname,'(a,i1)') 'zaxis_',naxes_z
        else
           write(axisname,'(a,i2)') 'zaxis_',naxes_z
        endif
        call mpp_write_meta(unit,z_axes(naxes_z),axisname,'none',axisname, &
             data=axisdata(1:siz_z_axes(j)),cartesian='Z')
     end do

     ! write out time axis  
     call mpp_write_meta(unit,t_axes,&
          'Time','time level','Time',cartesian='T')
     ! write metadata for fields
     do j = start_var,next_var-1
        cur_var => fileObj%var(j)
        if(cur_var%siz(4) > 1 .AND. cur_var%siz(4) .NE. fileObj%max_ntime ) call mpp_error(FATAL, &
         "fms_io(save_restart): "//trim(cur_var%name)//" in file "//trim(fileObj%name)// &
         " has more than one time level, but number of time level is not equal to max_ntime")

        if(cur_var%ndim == 0) then
           num_var_axes = 1
           var_axes(1) = t_axes
        else if(cur_var%ndim == 1) then
           num_var_axes = 1
           var_axes(1) = x_axes(cur_var%id_axes(1))
           if(cur_var%siz(4) == fileObj%max_ntime) then
              num_var_axes = 2
              var_axes(2) = t_axes
           end if
        else if(cur_var%ndim == 2) then
           num_var_axes = 2
           var_axes(1) = x_axes(cur_var%id_axes(1))
           var_axes(2) = y_axes(cur_var%id_axes(2))
           if(cur_var%siz(4) == fileObj%max_ntime) then
              num_var_axes = 3
              var_axes(3) = t_axes
           end if
        else if(cur_var%ndim == 3) then
           num_var_axes = 3
           var_axes(1) = x_axes(cur_var%id_axes(1))
           var_axes(2) = y_axes(cur_var%id_axes(2))
           var_axes(3) = z_axes(cur_var%id_axes(3))
           if(cur_var%siz(4) == fileObj%max_ntime) then
              num_var_axes = 4
              var_axes(4) = t_axes
           end if
        end if
        call mpp_write_meta(unit,cur_var%field, var_axes(1:num_var_axes), cur_var%name, &
                 cur_var%units,cur_var%longname,pack=1)
     enddo

     ! write values for ndim of spatial axes
     do j = 1, naxes_x
        call mpp_write(unit,x_axes(j))
     enddo
     do j = 1, naxes_y
        call mpp_write(unit,y_axes(j))
     enddo
     do j = 1, naxes_z
        call mpp_write(unit,z_axes(j))
     enddo

     ! write data of each field
     do k = 1, fileObj%max_ntime
        do j=start_var,next_var-1
           cur_var => fileObj%var(j)
           tlev=k
           ! If some fields only have one time level, we do not need to write the second level, just keep
           ! the data missing.
           if(k <= cur_var%siz(4)) then
              if(cur_var%domain_present) then  ! one 2-D or 3-D case possible present domain
                 if( Associated(fileObj%p2dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, array_domain(cur_var%domain_idx), fileObj%p2dr(k,j)%p, tlev, &
                                   default_data=cur_var%default_data)
                 else if( Associated(fileObj%p3dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, array_domain(cur_var%domain_idx), fileObj%p3dr(k,j)%p, tlev, &
                                   default_data=cur_var%default_data)
                 else if( Associated(fileObj%p2di(k,j)%p) ) then
                    allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                    r2d = fileObj%p2di(k,j)%p
                    call mpp_write(unit, cur_var%field, array_domain(cur_var%domain_idx), r2d, tlev, &
                                   default_data=cur_var%default_data)
                    deallocate(r2d)
                 else if( Associated(fileObj%p3di(k,j)%p) ) then
                    allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                    r3d = fileObj%p3di(k,j)%p
                    call mpp_write(unit, cur_var%field, array_domain(cur_var%domain_idx), r3d, tlev, &
                                   default_data=cur_var%default_data)
                    deallocate(r3d)
                 else
                    call mpp_error(FATAL, "fms_io(save_restart): domain is present and thread_w  "// &
                         "is MPP_MULTI, field "//trim(cur_var%name)//" of file "//trim(fileObj%name)// &
                         ", but none of p2dr, p3dr, p2di and p3di is associated") 
                 end if
              else if (thread_w == MPP_MULTI .or. (fileObj%is_root_pe.and.thread_w == MPP_SINGLE)) then     
                 if ( Associated(fileObj%p0dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, fileObj%p0dr(k,j)%p, tlev)
                 else if ( Associated(fileObj%p1dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, fileObj%p1dr(k,j)%p, tlev)
                 else if ( Associated(fileObj%p2dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, fileObj%p2dr(k,j)%p, tlev)
                 else if ( Associated(fileObj%p3dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, fileObj%p3dr(k,j)%p, tlev)
                 else if ( Associated(fileObj%p0di(k,j)%p) ) then
                    r0d =  fileObj%p0di(k,j)%p
                    call mpp_write(unit, cur_var%field, r0d,                  tlev)
                 else if ( Associated(fileObj%p1di(k,j)%p) ) then
                    allocate(r1d(cur_var%siz(1)) )
                    r1d = fileObj%p1di(k,j)%p
                    call mpp_write(unit, cur_var%field, r1d,                  tlev)
                    deallocate(r1d)
                 else if ( Associated(fileObj%p2di(k,j)%p) ) then
                    allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                    r2d = fileObj%p2di(k,j)%p
                    call mpp_write(unit, cur_var%field, r2d,                  tlev)
                    deallocate(global_r2d, r2d)
                 else if ( Associated(fileObj%p3di(k,j)%p) ) then
                    allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                    r3d = fileObj%p3di(k,j)%p
                    call mpp_write(unit, cur_var%field, r3d,                  tlev)
                    deallocate(r3d)
                 else
                    call mpp_error(FATAL, "fms_io(save_restart): There is no pointer associated with the data of  field "// &
                         trim(cur_var%name)//" of file "//trim(fileObj%name) )
                 end if
              end if
           end if
        enddo ! end j loop
     enddo ! end k loop
     call mpp_close(unit)
     nfiles = nfiles+1
  enddo

  cur_var =>NULL()

  if(print_chksum) call write_chksum(fileObj, MPP_OVERWR)

end subroutine save_restart

!-------------------------------------------------------------------------------
!    This subroutine will calculate chksum and print out chksum information.
!
subroutine write_chksum(fileObj, action)
  type(restart_file_type), intent(inout) :: fileObj
  integer,                 intent(in)    :: action
  integer(LONG_KIND)                     :: data_chksum
  integer                                :: j, k, outunit
  type(var_type), pointer, save          :: cur_var=>NULL()
  character(len=32)                      :: routine_name

  if(action == MPP_OVERWR) then
     routine_name = "save_restart"
  else if(action == MPP_RDONLY) then
     routine_name = "restore_state"
  else
     call mpp_error(FATAL, "fms_io_mod(write_chksum): action should be MPP_OVERWR or MPP_RDONLY")
  endif

  do j=1,fileObj%nvar
     cur_var => fileObj%var(j)
     if(action == MPP_OVERWR .OR. (action == MPP_RDONLY .AND. cur_var%initialized) ) then 
        do k = 1, cur_var%siz(4)
           if ( Associated(fileObj%p0dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p0dr(k,j)%p, (/mpp_pe()/) )
           else if ( Associated(fileObj%p1dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p1dr(k,j)%p, (/mpp_pe()/) )
           else if ( Associated(fileObj%p2dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p2dr(k,j)%p(cur_var%is:cur_var%ie,cur_var%js:cur_var%je) )
           else if ( Associated(fileObj%p3dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p3dr(k,j)%p(cur_var%is:cur_var%ie,cur_var%js:cur_var%je, :) )
           else if ( Associated(fileObj%p0di(k,j)%p) ) then
              data_chksum = fileObj%p0di(k,j)%p
           else if ( Associated(fileObj%p1di(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p1di(k,j)%p, (/mpp_pe()/) )
           else if ( Associated(fileObj%p2di(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p2di(k,j)%p(cur_var%is:cur_var%ie,cur_var%js:cur_var%je) )
           else if ( Associated(fileObj%p3di(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p3di(k,j)%p(cur_var%is:cur_var%ie,cur_var%js:cur_var%je, :))
           else
              call mpp_error(FATAL, "fms_io(write_chksum): There is no pointer associated with the data of  field "// &
                   trim(cur_var%name)//" of file "//trim(fileObj%name) )
           end if
           outunit = stdout()
           write(outunit,'(a, I1, a, I16)')'fms_io('//trim(routine_name)//'): At time level = ', k, ', chksum for "'// &
                trim(cur_var%name)// '" of "'// trim(fileObj%name)// '" = ', data_chksum    

        enddo
     endif
  enddo
  cur_var =>NULL()

end subroutine write_chksum

!-------------------------------------------------------------------------------
!
!    This subroutine reads the model state from previously
!    generated files.  All restart variables are read from the first
!    file in the input filename list in which they are found.

subroutine restore_state_all(fileObj, directory)
  type(restart_file_type), intent(inout)       :: fileObj
  character(len=*),      intent(in), optional  :: directory

! Arguments: 
!  (in)      directory - The directory where the restart or save
!                        files should be found. The default is 'INPUT'

  character(len=128) :: dir    
  character(len=256) :: restartpath ! The restart file path (dir/file). 
  character(len=200) :: filepath    ! The path (dir/file) to the file being opened.
  character(len=8)   :: suffix      ! A suffix (like "_2") that is added to any
                                    ! additional restart files.
  character(len=80)  :: varname     ! A variable's name.
  character(len=256) :: filename
  integer            :: num_restart ! The number of restart files that have already
                                    ! been opened.
  integer            :: nfile       ! The number of files (restart files and others
                                    ! explicitly in filename) that are open.
  integer   :: unit(max_split_file) ! The mpp unit of all open files.
  type(var_type), pointer, save       :: cur_var=>NULL()
  integer                             :: ndim, nvar, natt, ntime, tlev, siz
  type(fieldtype), allocatable        :: fields(:)
  logical                             :: fexist, domain_present
  integer                             :: j, n, l, k, missing_fields, domain_idx
  integer                             :: tile_id(1)
  real, allocatable, dimension(:,:,:) :: r3d
  real, allocatable, dimension(:,:)   :: r2d
  real, allocatable, dimension(:)     :: r1d  
  real                                :: r0d
  type(domain2d), pointer, save       :: io_domain=>NULL()

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(restore_state_all): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  dir = 'INPUT'
  if(present(directory)) dir = directory

  num_restart = 0
  nfile = 0
  if(len_trim(dir) > 0) then
     restartpath = trim(dir)//"/"// trim(fileObj%name)
  else
     restartpath = trim(fileObj%name)
  end if

  domain_present = .false.
  do j = 1, fileObj%nvar
     if (fileObj%var(j)%domain_present) then
        domain_present = .true.
        domain_idx = fileObj%var(j)%domain_idx
        exit
     end if
  end do

  !--- first open all the restart files
  !--- NOTE: For distributed restart file, we are assuming there is only one file exist.

  inquire (file=trim(restartpath)//trim(pe_name), exist=fexist)
  if(.NOT. fexist .and. domain_present) then
     io_domain => mpp_get_io_domain(array_domain(domain_idx))
     if(associated(io_domain)) then
        tile_id = mpp_get_tile_id(io_domain)
        if(mpp_npes() > 10000) then
           write(filename, '(a,i6.6)' ) trim(restartpath)//'.', tile_id(1)
        else
           write(filename, '(a,i4.4)' ) trim(restartpath)//'.', tile_id(1)
        endif
        inquire (file=trim(filename), exist = fexist)
     endif
     io_domain => NULL()
  endif
  if(fexist) then
     nfile = 1
     if(domain_present) then
        call mpp_open(unit(nfile), trim(restartpath), form=form,action=MPP_RDONLY,threading=thread_r, &
             fileset=MPP_MULTI, domain=array_domain(domain_idx) )
     else
        call mpp_open(unit(nfile), trim(restartpath), form=form,action=MPP_RDONLY,threading=thread_r, &
             fileset=MPP_MULTI)
     endif
  else
     do while(.true.)
        if (num_restart < 10) then
           write(suffix,'("_",I1)') num_restart
        else
           write(suffix,'("_",I2)') num_restart
        endif
        if (num_restart > 0) then
           siz = len_trim(restartpath)
           if(restartpath(siz-2:siz) == ".nc") then
              filepath = restartpath(1:siz-3)//trim(suffix)
           else      
              filepath = trim(restartpath) // trim(suffix)
           end if
        else
           filepath = trim(restartpath)
        end if
        inquire (file=trim(filepath), exist=fexist) 
        if(.not. fexist) inquire(file=trim(filepath)//".nc", exist=fexist)
        if(fexist) then
           nfile = nfile + 1
           if(nfile > max_split_file) call mpp_error(FATAL, &
                "fms_io(restore_state_all): nfile is larger than max_split_file, increase max_split_file")
           call mpp_open(unit(nfile), trim(filepath), form=form,action=MPP_RDONLY,threading=thread_r, &
                fileset=MPP_SINGLE)
        else
           exit
        end if
        num_restart = num_restart + 1
     end do
  end if
  if(nfile == 0) call mpp_error(FATAL, "fms_io(restore_state_all): unable to find any restart files "// &
       "specified by "//trim(restartpath))


  ! Read each variable from the first file in which it is found.
  do n=1,nfile
     call mpp_get_info(unit(n), ndim, nvar, natt, ntime)

     allocate(fields(nvar))
     call mpp_get_fields(unit(n),fields(1:nvar))

     missing_fields = 0

     do j=1,fileObj%nvar
        cur_var => fileObj%var(j)
        domain_present = cur_var%domain_present
        domain_idx = cur_var%domain_idx
        do l=1, nvar
           call mpp_get_atts(fields(l),name=varname)
           if (lowercase(trim(varname)) == lowercase(trim(cur_var%name))) then
              cur_var%initialized = .true.
              do k = 1, cur_var%siz(4)
                 tlev = k
                 if(domain_present) then        
                    if( Associated(fileObj%p2dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), array_domain(domain_idx), fileObj%p2dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p3dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), array_domain(domain_idx), fileObj%p3dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p2di(k,j)%p) ) then
                       allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                       r2d = 0
                       call mpp_read(unit(n), fields(l), array_domain(domain_idx), r2d, tlev)
                       fileObj%p2di(k,j)%p = r2d
                       deallocate(r2d)
                    else if( Associated(fileObj%p3di(k,j)%p) ) then
                       allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                       r3d = 0
                       call mpp_read(unit(n), fields(l), array_domain(domain_idx), r3d, tlev)
                       fileObj%p3di(k,j)%p = r3d
                       deallocate(r3d)
                    else
                       call mpp_error(FATAL, "fms_io(restore_state_all): domain is present for the field "//trim(varname)// &
                            " of file "//trim(fileObj%name)//", but none of p2dr, p3dr, p2di and p3di is associated") 
                    end if
                 else
                    if( Associated(fileObj%p0dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), fileObj%p0dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p1dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), fileObj%p1dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p2dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), fileObj%p2dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p3dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), fileObj%p3dr(k,j)%p, tlev) 
                    else if( Associated(fileObj%p0di(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), r0d, tlev)
                       fileObj%p0di(k,j)%p = r0d
                    else if( Associated(fileObj%p1di(k,j)%p) ) then
                       allocate(r1d(cur_var%siz(1)) )
                       call mpp_read(unit(n), fields(l), r1d, tlev)                
                       fileObj%p1di(k,j)%p = r1d
                       deallocate(r1d)
                    else if( Associated(fileObj%p2di(k,j)%p) ) then
                       allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                       r2d = 0
                       call mpp_read(unit(n), fields(l), r2d, tlev)                
                       fileObj%p2di(k,j)%p = r2d
                       deallocate(r2d)
                    else if( Associated(fileObj%p3di(k,j)%p) ) then
                       allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                       r3d = 0
                       call mpp_read(unit(n), fields(l), r3d, tlev)                
                       fileObj%p3di(k,j)%p = r3d
                       deallocate(r3d)
                    else
                       call mpp_error(FATAL, "fms_io(restore_state_all): There is no pointer "//&
                            "associated with the data of  field "// trim(varname)//" of file "//trim(fileObj%name) )
                    end if
                 end if
              end do
              exit ! Start search for next restart variable.
           endif
        enddo
        if (l>nvar) missing_fields = missing_fields+1
     enddo

     deallocate(fields)
     if (missing_fields == 0) exit
  enddo

  do n=1,nfile
     call close_file(unit(n))
  enddo

  ! check whether all fields have been found
  do j = 1, fileObj%nvar
     if( .NOT. fileObj%var(j)%initialized ) then
        if( fileObj%var(j)%mandatory ) then
           call mpp_error(FATAL, "fms_io(restore_state_all): unable to find mandatory variable "// &
                trim(fileObj%var(j)%name)//" in restart file "//trim(fileObj%name) )
        end if
     end if
  end do
  cur_var =>NULL()

  if(print_chksum) call write_chksum(fileObj, MPP_RDONLY )

end subroutine restore_state_all

!-------------------------------------------------------------------------------
!
!    This subroutine reads the model state from previously
!    generated files.  All restart variables are read from the first
!    file in the input filename list in which they are found.

subroutine restore_state_one_field(fileObj, id_field, directory)
  type(restart_file_type), intent(inout)       :: fileObj
  integer,                 intent(in)          :: id_field
  character(len=*),      intent(in), optional  :: directory

! Arguments: 
!  (in)      directory - The directory where the restart or save
!                        files should be found. The default is 'INPUT'

  character(len=128) :: dir    
  character(len=256) :: restartpath ! The restart file path (dir/file). 
  character(len=200) :: filepath    ! The path (dir/file) to the file being opened.
  character(len=8)   :: suffix      ! A suffix (like "_2") that is added to any
                                    ! additional restart files.
  character(len=80)  :: varname     ! A variable's name.
  character(len=256) :: filename
  integer            :: num_restart ! The number of restart files that have already
                                    ! been opened.
  integer            :: nfile       ! The number of files (restart files and others
                                    ! explicitly in filename) that are open.
  integer   :: unit(max_split_file) ! The mpp unit of all open files.
  type(var_type), pointer, save       :: cur_var=>NULL()
  integer                             :: ndim, nvar, natt, ntime, tlev, siz
  integer                             :: tile_id(1)
  type(fieldtype), allocatable        :: fields(:)
  logical                             :: fexist, domain_present
  integer                             :: j, n, l, k, missing_fields, domain_idx
  real, allocatable, dimension(:,:,:) :: r3d
  real, allocatable, dimension(:,:)   :: r2d
  real, allocatable, dimension(:)     :: r1d  
  real                                :: r0d
  type(domain2d), pointer, save       :: io_domain=>NULL()

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(restore_state_one_field): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  dir = 'INPUT'
  if(present(directory)) dir = directory

  cur_var => fileObj%var(id_field)
  domain_present = cur_var%domain_present
  domain_idx = cur_var%domain_idx

  num_restart = 0
  nfile = 0
  if(len_trim(dir) > 0) then
     restartpath = trim(dir)//"/"// trim(fileObj%name)
  else
     restartpath = trim(fileObj%name)
  end if
  !--- first open all the restart files
  !--- NOTE: For distributed restart file, we are assuming there is only one file exist.
  inquire (file=trim(restartpath)//trim(pe_name), exist=fexist)     
  if(.NOT. fexist .and. domain_present) then
     io_domain => mpp_get_io_domain(array_domain(domain_idx))
     if(associated(io_domain)) then
        tile_id = mpp_get_tile_id(io_domain)
        if(mpp_npes()>10000) then
           write(filename, '(a,i6.6)' ) trim(restartpath)//'.', tile_id(1)
        else
           write(filename, '(a,i4.4)' ) trim(restartpath)//'.', tile_id(1)
        endif
        inquire (file=trim(filename), exist = fexist)
     endif
     io_domain=>NULL()
  endif

  if(fexist) then
     nfile = 1
     if(domain_present) then
        call mpp_open(unit(nfile), trim(restartpath), form=form,action=MPP_RDONLY,threading=thread_r, &
             fileset=MPP_MULTI, domain=array_domain(domain_idx) )
     else
        call mpp_open(unit(nfile), trim(restartpath), form=form,action=MPP_RDONLY,threading=thread_r, &
             fileset=MPP_MULTI)
     endif
  else
     do while(.true.)
        if (num_restart < 10) then
           write(suffix,'("_",I1)') num_restart
        else
           write(suffix,'("_",I2)') num_restart
        endif
        if (num_restart > 0) then
           siz = len_trim(restartpath)
           if(restartpath(siz-2:siz) == ".nc") then
              filepath = restartpath(1:siz-3)//trim(suffix)
           else      
              filepath = trim(restartpath) // trim(suffix)
           end if
        else
           filepath = trim(restartpath)
        end if
        inquire (file=trim(filepath), exist=fexist) 
        if(.not. fexist) inquire(file=trim(filepath)//".nc", exist=fexist)
        if(fexist) then
           nfile = nfile + 1
           if(nfile > max_split_file) call mpp_error(FATAL, &
                "fms_io(restore_state_one_field): nfile is larger than max_split_file, increase max_split_file")
           call mpp_open(unit(nfile), trim(filepath), form=form,action=MPP_RDONLY,threading=thread_r, &
                fileset=MPP_SINGLE)
        else
           exit
        end if
        num_restart = num_restart + 1
     end do
  end if
  if(nfile == 0) call mpp_error(FATAL, "fms_io(restore_state_one_field): unable to find any restart files "// &
       "specified by "//trim(restartpath))


  ! Read each variable from the first file in which it is found.
  do n=1,nfile
     call mpp_get_info(unit(n), ndim, nvar, natt, ntime)

     allocate(fields(nvar))
     call mpp_get_fields(unit(n),fields(1:nvar))

     missing_fields = 0
     j = id_field
     do l=1, nvar
        call mpp_get_atts(fields(l),name=varname)
        if (lowercase(trim(varname)) == lowercase(trim(cur_var%name))) then
           cur_var%initialized = .true.
           do k = 1, cur_var%siz(4)
              tlev = k
              if(domain_present) then        
                 if( Associated(fileObj%p2dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), array_domain(domain_idx), fileObj%p2dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p3dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), array_domain(domain_idx), fileObj%p3dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p2di(k,j)%p) ) then
                    allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                    r2d = 0
                    call mpp_read(unit(n), fields(l), array_domain(domain_idx), r2d, tlev)
                    fileObj%p2di(k,j)%p = r2d
                    deallocate(r2d)
                 else if( Associated(fileObj%p3di(k,j)%p) ) then
                    allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                    r3d = 0
                    call mpp_read(unit(n), fields(l), array_domain(domain_idx), r3d, tlev)
                    fileObj%p3di(k,j)%p = r3d
                    deallocate(r3d)
                 else
                    call mpp_error(FATAL, "fms_io(restore_state_one_field): domain is present for the field "//trim(varname)// &
                         " of file "//trim(fileObj%name)//", but none of p2dr, p3dr, p2di and p3di is associated") 
                 end if
              else
                 if( Associated(fileObj%p0dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), fileObj%p0dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p1dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), fileObj%p1dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p2dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), fileObj%p2dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p3dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), fileObj%p3dr(k,j)%p, tlev) 
                 else if( Associated(fileObj%p0di(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), r0d, tlev)
                    fileObj%p0di(k,j)%p = r0d
                 else if( Associated(fileObj%p1di(k,j)%p) ) then
                    allocate(r1d(cur_var%siz(1)) )
                    call mpp_read(unit(n), fields(l), r1d, tlev)                
                    fileObj%p1di(k,j)%p = r1d
                    deallocate(r1d)
                 else if( Associated(fileObj%p2di(k,j)%p) ) then
                    allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                    r2d = 0
                    call mpp_read(unit(n), fields(l), r2d, tlev)                
                    fileObj%p2di(k,j)%p = r2d
                    deallocate(r2d)
                 else if( Associated(fileObj%p3di(k,j)%p) ) then
                    allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                    r3d = 0
                    call mpp_read(unit(n), fields(l), r3d, tlev)                
                    fileObj%p3di(k,j)%p = r3d
                    deallocate(r3d)
                 else
                    call mpp_error(FATAL, "fms_io(restore_state_one_field): There is no pointer "// &
                         "associated with the data of  field "//trim(varname)//" of file "//trim(fileObj%name) )
                 end if
              end if
           end do
           exit ! Start search for next restart variable.
        endif
     enddo
     if (l>nvar) missing_fields = missing_fields+1
     deallocate(fields)
     if (missing_fields == 0) exit
  enddo

  do n=1,nfile
     call close_file(unit(n))
  enddo

  ! check whether the field have been found
  if( .NOT. fileObj%var(id_field)%initialized ) then
     if( fileObj%var(id_field)%mandatory ) then
        call mpp_error(FATAL, "fms_io(restore_state_one_field): unable to find mandatory variable "// &
             trim(fileObj%var(id_field)%name)//" in restart file "//trim(fileObj%name) )
     end if
  end if
  cur_var =>NULL()

end subroutine restore_state_one_field

!-------------------------------------------------------------------------------
!
!     This routine will setup one entry to be written out 
!
!-------------------------------------------------------------------------------
subroutine setup_one_field(fileObj, filename, fieldname, field_siz, index_field,  domain, mandatory, &
                           no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),         intent(in)           :: filename, fieldname 
  integer, dimension(:),    intent(in)           :: field_siz
  integer,                  intent(out)          :: index_field
  type(domain2d), optional, intent(in), target   :: domain
  real,           optional, intent(in)           :: data_default
  logical,        optional, intent(in)           :: no_domain   
  integer,        optional, intent(in)           :: position, tile_count
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units

  !--- local variables
  integer                         :: i, domain_idx
  integer                         :: ishift, jshift
  integer                         :: gxsize, gysize
  integer                         :: cxsize, cysize
  integer                         :: dxsize, dysize
  real                            :: default_data
  logical                         :: is_no_domain = .false.
  character(len=256)              :: fname, filename2, append_string
  type(domain2d), pointer, save   :: d_ptr   =>NULL()
  type(var_type), pointer, save   :: cur_var =>NULL()
  integer                         :: length

  if(ANY(field_siz < 1)) then
     call mpp_error(FATAL, "fms_io(setup_one_field): each entry of field_size should be a positive integer")
  end if

  if(PRESENT(data_default))then
     default_data=data_default
  else
     default_data=0.
  endif

  if(present(tile_count) .AND. .not. present(domain)) call mpp_error(FATAL, &
         'fms_io(setup_one_field): when tile_count is present, domain must be present')

  is_no_domain = .false.
  if (PRESENT(no_domain)) THEN
     is_no_domain = no_domain
  end if
  if(is_no_domain) then
!     comments the following to allow passing domain when writing 1-D or scalar variables. 
!     if(PRESENT(domain)) &
!       call mpp_error(FATAL, 'fms_io(setup_one_field): no_domain cannot be .true. when optional argument domain is present.')
  else if(PRESENT(domain))then
     d_ptr => domain
  else if (ASSOCIATED(Current_domain)) then
     d_ptr => Current_domain
  endif

  !--- remove .nc from file name
  length = len_trim(filename)
  if(filename(length-2:length) == '.nc') then
     filename2 = filename(1:length-3)
  else
     filename2 = filename(1:length)
  end if

  !Append a string to the file name
  append_string=''
  !If the filename_appendix  is set override the passed argument. 
  if(len_trim(filename_appendix) > 0)   append_string = filename_appendix

  if(len_trim(append_string) > 0) filename2 = trim(filename2)//'.'//trim(append_string)

  !JWD:  This is likely a temporary fix. Since fms_io needs to know tile_count,
  !JWD:  I just don't see how the physics can remain "tile neutral"
  !z1l:  one solution is add one more public interface called set_tile_count
  call get_mosaic_tile_file(filename2, fname, is_no_domain, domain, tile_count)

  if(Associated(fileObj%var) ) then
     ! make sure the consistency of file name
     if(trim(fileObj%name) .NE. trim(fname)) call mpp_error(FATAL, 'fms_io(setup_one_field): filename = '// &
         trim(fname)//' is not consistent with the filename of the restart object = '//trim(fileObj%name) )
  else
     allocate(fileObj%var(max_fields) )
     allocate(fileObj%p0dr(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p1dr(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p2dr(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p3dr(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p0di(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p1di(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p2di(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p3di(MAX_TIME_LEVEL_REGISTER, max_fields))
     !--- make sure fname is not used in other restart_file_type object.
     do i = 1, num_registered_files
        if(trim(fname) == trim(registered_file(i)) ) call mpp_error(FATAL, &
          'fms_io(setup_one_field): '//trim(fname)//' is already registered with other restart_file_type data')
     end do
     num_registered_files = num_registered_files + 1
     registered_file(num_registered_files) = trim(fname)
     fileObj%name = trim(fname)         
     fileObj%tile_count=1
     if(present(tile_count)) fileObj%tile_count = tile_count
     if(ASSOCIATED(d_ptr))then
        fileObj%is_root_pe = mpp_domain_is_tile_root_pe(d_ptr)
     else
        fileObj%is_root_pe = mpp_pe() == mpp_root_pe()
     endif
     fileObj%max_ntime = field_siz(4)
     fileObj%nvar      = 0
     !-- allocate memory
     do i = 1, max_fields
        fileObj%var(i)%name           = 'none'
        fileObj%var(i)%domain_present = .false.
        fileObj%var(i)%domain_idx     = -1
        fileObj%var(i)%is_dimvar      = .false.
        fileObj%var(i)%position       = CENTER
        fileObj%var(i)%siz(:)         = 0
        fileObj%var(i)%gsiz(:)        = 0
        fileObj%var(i)%id_axes(:)     = -1
        fileObj%var(i)%longname       = "";
        fileObj%var(i)%units          = "none";
        fileObj%var(i)%mandatory      = .true.
        fileObj%var(i)%initialized    = .false.
     end do
  endif

  ! check if the field is new or not and get position and dimension of the field
  index_field = -1
  do i = 1, fileObj%nvar
     if(trim(fileObj%var(i)%name) == trim(fieldname)) then
        index_field = i
        exit
     end if 
  end do

  if(index_field > 0) then
     cur_var   => fileObj%var(index_field)
     if(cur_var%siz(1) .NE. field_siz(1) .OR. cur_var%siz(2) .NE. field_siz(2) .OR. cur_var%siz(3) .NE. field_siz(3) ) &
        call mpp_error(FATAL, 'fms_io(setup_one_field): field size mismatch for field '// &
                       trim(fieldname)//' of file '//trim(filename) )

     cur_var%siz(4) =  cur_var%siz(4) + field_siz(4)
     if(fileObj%max_ntime < cur_var%siz(4) ) fileObj%max_ntime = cur_var%siz(4)
     ! the time level should be no larger than MAX_TIME_LEVEL_REGISTER ( = 2) 
     if( cur_var%siz(4) > MAX_TIME_LEVEL_REGISTER ) call mpp_error(FATAL, 'fms_io(setup_one_field): ' // &
          'the time level of field '//trim(cur_var%name)//' in file '//trim(fileObj%name)// &
          ' is greater than MAX_TIME_LEVEL_REGISTER(=2), increase MAX_TIME_LEVEL_REGISTER or check your code')
  else 
     fileObj%nvar = fileObj%nvar +1
     if(fileObj%nvar>max_fields) then
        write(error_msg,'(I3,"/",I3)') fileObj%nvar, max_fields 
        call  mpp_error(FATAL,'fms_io(setup_one_field): max_fields exceeded, needs increasing, nvar/max_fields=' &
             //trim(error_msg))
     endif
     index_field =  fileObj%nvar
     cur_var   => fileObj%var(index_field)
     cur_var%siz(:)  = field_siz(:)
     cur_var%gsiz(3) = field_siz(3)
     cur_var%csiz(3) = field_siz(3)
     cur_var%name = fieldname
     cur_var%default_data = default_data
     if(present(mandatory)) cur_var%mandatory = mandatory
     if(present(longname)) then
        cur_var%longname = longname
     else
        cur_var%longname = fieldname
     end if
     if(present(units))    cur_var%units    = units
     if(present(position)) cur_var%position = position
     cur_var%is = 1; cur_var%ie =  cur_var%siz(1) 
     cur_var%js = 1; cur_var%je =  cur_var%siz(2) 
     if(ASSOCIATED(d_ptr)) then
        cur_var%domain_present = .true.
        domain_idx = lookup_domain(d_ptr)
        if(domain_idx == -1) then
           num_domains = num_domains + 1
           if(num_domains > max_domains) call  mpp_error(FATAL,'fms_io(setup_one_field), 1: max_domains exceeded,' &
                //' needs increasing')
           domain_idx = num_domains
           array_domain(domain_idx) = d_ptr
           call mpp_get_domain_components(array_domain(domain_idx), domain_x(domain_idx), domain_y(domain_idx), &
                tile_count=tile_count)
        endif
        cur_var%domain_idx = domain_idx
        call mpp_get_domain_shift ( array_domain(domain_idx), ishift, jshift, position)
        call mpp_get_global_domain(array_domain(domain_idx), xsize=gxsize,ysize=gysize,tile_count=tile_count)
        call mpp_get_compute_domain(array_domain(domain_idx), xsize = cxsize, ysize = cysize, tile_count=tile_count)
        call mpp_get_data_domain   (array_domain(domain_idx), xsize = dxsize, ysize = dysize, tile_count=tile_count)
        if (ishift .NE. 0) then
           cxsize = cxsize+ishift; dxsize = dxsize+ishift; gxsize = gxsize + ishift
        end if
        if (jshift .NE. 0) then
           cysize = cysize+jshift; dysize = dysize+jshift; gysize = gysize + jshift
        endif
        if( (cur_var%siz(1) .NE. cxsize .AND. cur_var%siz(1) .NE. dxsize ) .OR. &
            (cur_var%siz(2) .NE. cysize .AND. cur_var%siz(2) .NE. dysize ) ) then
            call mpp_error(FATAL, 'fms_io(setup_one_field): data should be on either computer domain '//&
              'or data domain when domain is present for field '//trim(fieldname)//' of file '//trim(filename) )
        end if
        cur_var%is   = 1 + (cur_var%siz(1) - cxsize)/2
        cur_var%ie   = cur_var%is + cxsize - 1;
        cur_var%js   = 1 + (cur_var%siz(2) - cysize)/2
        cur_var%je   = cur_var%js + cysize - 1;
        cur_var%gsiz(1)   = gxsize
        cur_var%gsiz(2)   = gysize
        if(thread_w == MPP_MULTI) then
           call mpp_get_compute_domain(array_domain(domain_idx), xsize=cxsize,ysize=cysize,tile_count=tile_count)
           cur_var%csiz(1)   = cxsize
           cur_var%csiz(2)   = cysize
        else
           cur_var%csiz(1)   = cur_var%gsiz(1)
           cur_var%csiz(2)   = cur_var%gsiz(2)
        end if
     else
        cur_var%gsiz(1:2) = field_siz(1:2)
        cur_var%csiz(1:2) = field_siz(1:2)
        cur_var%domain_present=.false.
     endif
  end if

  d_ptr =>NULL()
  cur_var =>NULL()

end subroutine setup_one_field


!.....................................................................
subroutine write_data_2d_new(filename, fieldname, data, domain,    &
                             no_domain, position,tile_count, data_default)

  character(len=*), intent(in)                 :: filename, fieldname 
  real, dimension(:,:), intent(in)             :: data
  real, dimension(size(data,1),size(data,2),1) :: data_3d
  real, intent(in), optional                   :: data_default
  type(domain2d), intent(in), optional         :: domain
  logical, intent(in), optional                :: no_domain
  integer, intent(in), optional                :: position, tile_count
 
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_2d_new):need to call fms_io_init first')
  data_3d(:,:,1) = data(:,:)

  call write_data_3d_new(filename, fieldname, data_3d, domain, &
                         no_domain, position, tile_count, data_default)

end subroutine write_data_2d_new

! ........................................................
subroutine write_data_1d_new(filename, fieldname, data,domain, &
                             no_domain, tile_count, data_default)
  
  type(domain2d), intent(in), optional   :: domain
  character(len=*), intent(in)           :: filename, fieldname 
  real, dimension(:), intent(in)         :: data
  real, dimension(size(data(:)),1,1)     :: data_3d
  real, intent(in), optional             :: data_default
  logical, intent(in), optional          :: no_domain
  integer, intent(in), optional          :: tile_count
  
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_1d_new): module not initialized')  
  if(present(no_domain)) then
     if(.NOT. no_domain) call mpp_error(FATAL, 'fms_io(write_data_1d_new): no_domain should be true for field ' &
                                 //trim(fieldname)//' of file '//trim(filename) )
  end if
  data_3d(:,1,1) = data(:)  
  call write_data_3d_new(filename, fieldname, data_3d,domain,   &
                         no_domain=.true., tile_count=tile_count, data_default=data_default)  
end subroutine write_data_1d_new

! ..........................................................
subroutine write_data_scalar_new(filename, fieldname, data, domain, &
                                 no_domain, tile_count, data_default)

  type(domain2d), intent(in), optional   :: domain
  character(len=*), intent(in)           :: filename, fieldname 
  real, intent(in)                       :: data
  real, dimension(1,1,1)                 :: data_3d
  real, intent(in), optional             :: data_default
  logical, intent(in), optional          :: no_domain
  integer, intent(in), optional          :: tile_count
    
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_scalar_new):  module not initialized: '//fieldname)  
  if(present(no_domain)) then
     if(.NOT. no_domain) call mpp_error(FATAL, 'fms_io(write_data_scalar_new): no_domain should be true for field ' &
                                 //trim(fieldname)//' of file '//trim(filename) )
  end if

  data_3d(1,1,1) = data
  call write_data_3d_new(filename, fieldname, data_3d,domain, &
                         no_domain=.true., tile_count=tile_count, data_default=data_default)
end subroutine write_data_scalar_new

! ..........................................................

function lookup_field_r(nfile,fieldname)
! Given fieldname, this function returns the field position in the model's fields list

  integer, intent(in)          :: nfile
  character(len=*), intent(in) :: fieldname
  integer                      :: lookup_field_r
  integer                      :: j

  lookup_field_r=-1
  do j = 1, files_read(nfile)%nvar
     if (trim(files_read(nfile)%var(j)%name) == trim(fieldname)) then
        lookup_field_r = j
        exit 
     endif
  enddo
  return
end function lookup_field_r


!..........................................................

function lookup_domain(domain)
! given domain, this function returns the position of domain in array_domain or -1 if not found

  type(domain2d), intent(in) :: domain
  integer                    :: i, lookup_domain
  lookup_domain = -1
  do i =1, num_domains
     if(domain .EQ. array_domain(i)) then
        lookup_domain = i
        exit
     endif
  enddo
end function lookup_domain
!.........................................................
function lookup_axis(axis_sizes,siz,domains,dom)

! Given axis size (global), this function returns the axis id  

  integer, intent(in)      :: axis_sizes(:), siz
  type(domain1d), optional :: domains(:)
  type(domain1d), optional :: dom
  integer :: lookup_axis
  integer :: j


  lookup_axis=-1
  do j=1,size(axis_sizes(:))
     if (siz == axis_sizes(j)) then
        if (PRESENT(domains)) then
           if (dom .EQ. domains(j)) then 
              lookup_axis = j
              exit
           endif
        else
           lookup_axis = j
           exit
        endif
     endif
  enddo
  if (lookup_axis == -1) call mpp_error(FATAL,'fms_io(lookup_axis): could not find axis in set of axes')
end function lookup_axis
!.....................................................................
! <SUBROUTINE NAME="field_size">
!<DESCRIPTION>
! Given filename and fieldname, this subroutine returns the size of field
!</DESCRIPTION>
!   <TEMPLATE>
! call field_size(filename, fieldname, siz) 
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character" DIM="(*)">
!    File name
!   </IN>
!   <IN NAME="fieldname" TYPE="character" DIM="(*)">
!    Field  name
!   </IN>
!   <OUT NAME="siz" TYPE="integer" DIM="(*)">
!    siz must be a dimension(4) array to retrieve the size of the field
!   </OUT>
!   <OUT NAME="field_found" TYPE="logical, optional">
!    if this flag is present, field_size will not abort if
!    called for a non-existent field.
!    Instead it will return T or F depending on
!    whether or not the field was found.
!   </OUT>
subroutine field_size(filename, fieldname, siz, field_found, domain, no_domain )

  character(len=*), intent(in)                 :: filename, fieldname
  integer,       intent(inout)                 :: siz(:)
  logical,       intent(out), optional         :: field_found  
  type(domain2d), intent(in), optional, target :: domain
  logical,       intent(in),  optional         :: no_domain

  integer                              :: nfile, unit
  logical                              :: found, found_file
  character(len=256)                   :: actual_file
  logical                              :: read_dist, io_domain_exist, is_no_domain

  if (size(siz(:)) < 4) call mpp_error(FATAL,'fms_io(field_size): size array must be >=4 to receive field size of ' &
       //trim(fieldname)//' in file '// trim(filename))

  is_no_domain = .false.
  if(present(no_domain)) is_no_domain = no_domain

!--- first need to get the filename, when is_no_domain is true, only check file without tile
!--- if is_no_domain is false, first check no_domain=.false., then check no_domain = .true.
  found_file = get_file_name(filename, actual_file, read_dist, io_domain_exist, no_domain=is_no_domain, &
                             domain=domain)
  !--- when is_no_domain is true and file is not found, send out error message.
  if(is_no_domain .AND. .NOT. found_file) call mpp_error(FATAL, &
         'fms_io_mod(field_size): file '//trim(filename)//' and corresponding distributed file are not found')
  found = .false.
  if(found_file) then
     call get_file_unit(actual_file, unit, nfile, read_dist, io_domain_exist, domain=domain)
     call get_size(unit,fieldname,siz,found)    
  endif

  if(.not.found .AND. .not. is_no_domain) then
     found_file =  get_file_name(filename, actual_file, read_dist, io_domain_exist, no_domain=.true.)
     !--- when is_no_domain is true and file is not found, send out error message.
     if(.NOT. found_file) call mpp_error(FATAL, 'fms_io_mod(field_size): file ' //trim(filename)// &
          '(with the consideration of tile number) and corresponding distributed file are not found')  
     call get_file_unit(actual_file, unit, nfile, read_dist, io_domain_exist, domain=domain)
     call get_size(unit,fieldname,siz,found)
  endif

  if( PRESENT(field_found) )then
     field_found = found
  else if (.not. found )then
     call mpp_error(FATAL, 'fms_io(field_size): field '//trim(fieldname)//' NOT found in file '//trim(actual_file))
  end if

  return
end subroutine field_size
! </SUBROUTINE>

subroutine get_size(unit, fieldname, siz, found)
integer,          intent(in)    :: unit
character(len=*), intent(in)    :: fieldname
integer,          intent(inout) :: siz(:)
logical,          intent(out)   :: found

  character(len=128)             :: name
  character(len=1)               :: cart
  integer                        :: i, ndim, nvar, natt, ntime, siz_in(4), j, len
  type(fieldtype)                :: fields(max_fields)
  type(axistype)                 :: axes(max_fields)
     found = .false.
     call mpp_get_info(unit,ndim,nvar,natt,ntime)
     if (nvar > max_fields) then
        write(error_msg,'(I3,"/",I3)') nvar,max_fields 
        call  mpp_error(FATAL,'fms_io(field_size): max_fields too small, needs increasing, nvar/max_fields=' &
             //trim(error_msg))!//' in file '//trim(filename))
     endif
     call mpp_get_fields(unit,fields(1:nvar))
     do i=1, nvar
        call mpp_get_atts(fields(i),name=name)
        if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
           call mpp_get_atts(fields(i),ndim=ndim)
           call mpp_get_atts(fields(i),axes=axes(1:ndim))
           call mpp_get_atts(fields(i),siz=siz_in)
           siz = siz_in
           siz(4) = ntime
           if(ndim == 1) then
              call mpp_get_atts(axes(1), len=siz(1))
           end if
           do j = 1, ndim
              call mpp_get_atts(axes(j),len=len)
              call get_axis_cart(axes(j),cart)
              select case (cart)
              case ('X')
                 siz(1) = len
              case('Y')
                 siz(2) = len
              case('Z')
                 siz(3) = len
              case('T')
                 siz(4) = len
              end select
           enddo
           found = .true.
           exit
        endif
     enddo

     if(.not. found) then
        call mpp_get_axes(unit,axes(1:ndim))
        do i=1, ndim
           call mpp_get_atts(axes(i),name=name, len= siz_in(1))
           if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
              siz(1)= siz_in(1)
              found = .true.
              exit
           endif
        enddo
     endif
end subroutine get_size

! <SUBROUTINE NAME="read_data">
!<DESCRIPTION>
! This routine performs reading "fieldname" stored in "filename". The data values of fieldname
! will be stored in "data" at the end of this routine. For fieldname with multiple timelevel
! just repeat the routine with explicit timelevel in each call.
!</DESCRIPTION>
!   <TEMPLATE>
! call read_data(filename,fieldname,data,domain,timelevel)
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character" DIM="(*)">
!    File name
!   </IN>
!   <IN NAME="fieldname" TYPE="character" DIM="(*)">
!    Field  name
!   </IN>
!   <IN NAME="domain"  TYPE="domain, optional">
!   domain of fieldname
!   </IN>
!   <IN NAME="timelevel" TYPE="integer, optional">
!     time level of fieldname
!   </IN>
!   <OUT NAME="data"  TYPE="real">
!   array containing data of fieldname
!   </OUT>
!=====================================================================================
subroutine read_data_i3d_new(filename,fieldname,data,domain,timelevel, &
                             no_domain,position, tile_count)
  character(len=*),           intent(in)   :: filename, fieldname
  integer, dimension(:,:,:), intent(inout) :: data ! 3 dimensional data    
  type(domain2d), intent(in),   optional   :: domain
  integer, intent(in),          optional   :: timelevel
  logical, intent(in),          optional   :: no_domain 
  integer, intent(in) ,         optional   :: position, tile_count

  real, dimension(size(data,1),size(data,2),size(data,3)) :: r_data
  r_data = 0
  call read_data_3d_new(filename,fieldname,r_data,domain,timelevel, &
                        no_domain,position, tile_count)
  data = CEILING(r_data)
end subroutine read_data_i3d_new

subroutine read_data_i2d_new(filename,fieldname,data,domain,timelevel, &
                             no_domain,position, tile_count)
  character(len=*),         intent(in)   :: filename, fieldname
  integer, dimension(:,:), intent(inout) :: data ! 2 dimensional data    
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in),        optional   :: timelevel
  logical, intent(in),        optional   :: no_domain
  integer, intent(in) ,       optional   :: position, tile_count
  real, dimension(size(data,1),size(data,2)) :: r_data

  r_data = 0
  call read_data_2d_new(filename,fieldname,r_data,domain,timelevel, &
                        no_domain,position, tile_count)
  data = CEILING(r_data)
end subroutine read_data_i2d_new
!.....................................................................
subroutine read_data_i1d_new(filename,fieldname,data,domain,timelevel, &
                             no_domain, tile_count)
  character(len=*), intent(in)           :: filename, fieldname
  integer, dimension(:), intent(inout)   :: data ! 1 dimensional data    
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in) , optional         :: timelevel
  logical, intent(in), optional          :: no_domain 
  integer, intent(in), optional          :: tile_count

  real, dimension(size(data,1))        :: r_data

  call read_data_1d_new(filename,fieldname,r_data,domain,timelevel, &
                        no_domain, tile_count)
  data = CEILING(r_data)
end subroutine read_data_i1d_new
!.....................................................................
subroutine read_data_iscalar_new(filename,fieldname,data,domain,timelevel, &
                                 no_domain, tile_count)
  character(len=*), intent(in)           :: filename, fieldname
  integer, intent(inout)                 :: data     
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in) , optional         :: timelevel
  logical, intent(in), optional          :: no_domain 
  integer, intent(in), optional          :: tile_count

  real                                 :: r_data
  call read_data_scalar_new(filename,fieldname,r_data,domain,timelevel, &
                            no_domain, tile_count)
  data = CEILING(r_data)
end subroutine read_data_iscalar_new
!=====================================================================================
subroutine read_data_3d_new(filename,fieldname,data,domain,timelevel, &
      no_domain, position, tile_count)
  character(len=*),                  intent(in) :: filename, fieldname
  real, dimension(:,:,:),         intent(inout) :: data ! 3 dimensional data    
  type(domain2d), target, optional,  intent(in) :: domain
  integer,                optional,  intent(in) :: timelevel
  logical,                optional,  intent(in) :: no_domain 
  integer,                optional,  intent(in) :: position, tile_count

  character(len=256)            :: fname
  integer                       :: unit, siz_in(4)
  integer                       :: file_index  ! index of the opened file in array files
  integer                       :: tlev=1
  integer                       :: index_field ! position of the fieldname in the list of variables
  integer                       :: cxsize, cysize
  integer                       :: dxsize, dysize
  integer                       :: gxsize, gysize
  integer                       :: ishift, jshift
  logical                       :: is_no_domain = .false.
  logical                       :: read_dist, io_domain_exist, found_file
  type(domain2d), pointer, save :: d_ptr =>NULL()
  type(domain2d), pointer, save :: io_domain =>NULL()


! read disttributed files is used when reading restart files that are NOT mppnccombined. In this
! case PE 0 will read file_res.nc.0000, PE 1 will read file_res.nc.0001 and so forth.
! 
! namelist to be used with read_dist_files: threading_read=multi,
! threading_write=multi, fileset_write=multi.

! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(read_data_3d_new):  module not initialized')  
  is_no_domain = .false.
  if (PRESENT(no_domain)) THEN
!z1l:  comment out the following error check to allow reading 1-D or scalar variable with passing domain.
!!$     if(PRESENT(domain) .AND. no_domain) &
!!$       call mpp_error(FATAL, 'fms_io(read_data_3d_new): no_domain cannot be .true. when optional argument domain is present.')
     is_no_domain = no_domain
  endif  
 
  if(PRESENT(domain))then
     d_ptr => domain
  elseif (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
     d_ptr => Current_domain
  endif

  if(.not. PRESENT(domain) .and. .not. ASSOCIATED(Current_domain) ) is_no_domain = .true.

  found_file = get_file_name(filename, fname, read_dist, io_domain_exist, is_no_domain, domain,  tile_count)
  if(.not.found_file) call mpp_error(FATAL, 'fms_io_mod(read_data_3d_new): file ' //trim(filename)// &
          '(with the consideration of tile number) and corresponding distributed file are not found')  
  call get_file_unit(fname, unit, file_index, read_dist, io_domain_exist, domain=domain)

  siz_in(3) = size(data,3)
  if(is_no_domain .or. .NOT. associated(d_ptr) ) then
     gxsize = size(data,1)
     gysize = size(data,2)
  else if(read_dist) then
     if(io_domain_exist) then
        io_domain=>mpp_get_io_domain(d_ptr)
        call mpp_get_global_domain(io_domain, xsize = gxsize, ysize = gysize, tile_count=tile_count, position=position)
        io_domain=>NULL()
     else
        call mpp_get_compute_domain(d_ptr, xsize = gxsize, ysize = gysize, tile_count=tile_count, position=position)
     endif
  else
     call mpp_get_compute_domain(d_ptr, xsize = cxsize, ysize = cysize, tile_count=tile_count, position=position)
     call mpp_get_data_domain   (d_ptr, xsize = dxsize, ysize = dysize, tile_count=tile_count, position=position)
     call mpp_get_global_domain (d_ptr, xsize = gxsize, ysize = gysize, tile_count=tile_count, position=position)
     call mpp_get_domain_shift  (d_ptr, ishift, jshift, position)
     if( (size(data,1) .NE. cxsize .AND. size(data,1) .NE. dxsize) .OR. &
         (size(data,2) .NE. cysize .AND. size(data,2) .NE. dysize) )then
       call mpp_error(FATAL,'fms_io(read_data_3d_new): data should be on either computer domain '//&
                            'or data domain when domain is present. '//&
                            'shape(data)=',shape(data),'  cxsize,cysize,dxsize,dysize=',(/cxsize,cysize,dxsize,dysize/))
     end if 
  endif

  if (PRESENT(timelevel)) then
     tlev = timelevel
  else
     tlev = 1
  endif  

  if ((thread_r == MPP_MULTI).or.(mpp_pe()==mpp_root_pe())) then
     call get_field_id(unit, file_index, fieldname, index_field, is_no_domain, .false. )
     siz_in = files_read(file_index)%var(index_field)%siz
     if(files_read(file_index)%var(index_field)%is_dimvar ) then
        if (.not. read_dist) then
           if (siz_in(1) /= gxsize) &
                call mpp_error(FATAL,'fms_io(read_data_3d_new), field '//trim(fieldname)// &
                ' in file '//trim(filename)//' field size mismatch 2')
        endif
     else
        if (siz_in(1) /= gxsize .or. siz_in(2) /= gysize .or. siz_in(3) /= size(data,3)) then
           PRINT *, gxsize, gysize, size(data, 3), siz_in(1), siz_in(2), siz_in(3)
           call mpp_error(FATAL,'fms_io(read_data_3d_new), field '//trim(fieldname)// &
                ' in file '//trim(filename)//': field size mismatch 1')
        endif
     end if
     if ( tlev < 1 .or. files_read(file_index)%max_ntime < tlev)  then
        write(error_msg,'(I5,"/",I5)') tlev, files_read(file_index)%max_ntime
        call mpp_error(FATAL,'fms_io(read_data_3d_new): time level out of range, time level/max_time_level=' &
             //trim(error_msg)//' in field/file: '//trim(fieldname)//'/'//trim(filename))
     endif

     
     if(is_no_domain) then
        if (files_read(file_index)%var(index_field)%is_dimvar) then
           call mpp_get_axis_data(files_read(file_index)%var(index_field)%axis,data(:,1,1))
        else
           call mpp_read(unit,files_read(file_index)%var(index_field)%field,data(:,:,:),tlev)
        endif
     else 
        call mpp_read(unit,files_read(file_index)%var(index_field)%field,d_ptr,data,tlev,tile_count)
     endif
  endif  

  d_ptr =>NULL()

  return
end subroutine read_data_3d_new


!=====================================================================================
subroutine read_data_2d_region(filename,fieldname,data,start,nread,domain, &
                                 no_domain, tile_count)
  character(len=*),                  intent(in) :: filename, fieldname
  real, dimension(:,:),           intent(inout) :: data ! 3 dimensional data    
  integer, dimension(:),             intent(in) :: start, nread
  type(domain2d), target,  optional, intent(in) :: domain
  logical,                 optional, intent(in) :: no_domain 
  integer,                 optional, intent(in) :: tile_count
  character(len=256)            :: fname
  integer                       :: unit, siz_in(4)
  integer                       :: file_index  ! index of the opened file in array files
  integer                       :: tlev=1
  integer                       :: index_field ! position of the fieldname in the list of variables
  integer                       :: cxsize, cysize
  integer                       :: dxsize, dysize
  integer                       :: gxsize, gysize
  integer                       :: ishift, jshift
  logical                       :: is_no_domain = .false.
  logical                       :: read_dist, io_domain_exist, found_file
  type(domain2d), pointer, save :: d_ptr =>NULL()
  type(domain2d), pointer, save :: io_domain =>NULL()


! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(read_data_3d_new):  module not initialized')  
  is_no_domain = .false.
  if (PRESENT(no_domain)) is_no_domain = no_domain
 
  if(PRESENT(domain))then
     d_ptr => domain
  elseif (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
     d_ptr => Current_domain
  endif

  if(.not. PRESENT(domain) .and. .not. ASSOCIATED(Current_domain) ) is_no_domain = .true.

  found_file = get_file_name(filename, fname, read_dist, io_domain_exist, is_no_domain, domain,  tile_count)
  if(.not.found_file) call mpp_error(FATAL, 'fms_io_mod(read_data_2d_region): file ' //trim(filename)// &
          '(with the consideration of tile number) and corresponding distributed file are not found')  
  call get_file_unit(fname, unit, file_index, read_dist, io_domain_exist, domain=domain)


  if ((thread_r == MPP_MULTI).or.(mpp_pe()==mpp_root_pe())) then
     call get_field_id(unit, file_index, fieldname, index_field, is_no_domain, .false. )
     siz_in = files_read(file_index)%var(index_field)%siz
     if(files_read(file_index)%var(index_field)%is_dimvar) then
        call mpp_error(FATAL, 'fms_io_mod(read_data_2d_region): the field should not be a dimension variable')
     endif
     call mpp_read(unit,files_read(file_index)%var(index_field)%field,data,start, nread)
  endif  

  d_ptr =>NULL()

  return
end subroutine read_data_2d_region

!=====================================================================================
!--- we assume any text data are at most 2-dimensional and level is for first dimension
subroutine read_data_text(filename,fieldname,data,level)
  character(len=*), intent(in)   :: filename, fieldname
  character(len=*), intent(out)  :: data
  integer, intent(in) , optional :: level
  logical                        :: file_opened, found_file, read_dist, io_domain_exist
  integer                        :: lev, unit, index_field
  integer                        :: file_index
  character(len=256)             :: fname

! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(read_data_text):  module not initialized')  

  file_opened=.false.
  if (PRESENT(level)) then
     lev = level
  else
     lev = 1
  endif  

  found_file = get_file_name(filename, fname, read_dist, io_domain_exist, no_domain=.true. )
 if(.not.found_file) call mpp_error(FATAL, 'fms_io_mod(read_data_text): file ' //trim(filename)// &
          '(with the consideration of tile number) and corresponding distributed file are not found')  
  call get_file_unit(fname, unit, file_index, read_dist, io_domain_exist )

! Get info of this file and field   
  if ((thread_r == MPP_MULTI).or.(mpp_pe()==mpp_root_pe())) then
     call get_field_id(unit, file_index, fieldname, index_field, .true., .true. )     

     if ( lev < 1 .or. lev > files_read(file_index)%var(index_field)%siz(1) )  then
        write(error_msg,'(I5,"/",I5)') lev, files_read(file_index)%var(index_field)%siz(1)
        call mpp_error(FATAL,'fms_io(read_data_text): text level out of range, level/max_level=' &
             //trim(error_msg)//' in field/file: '//trim(fieldname)//'/'//trim(filename))
     endif

     call mpp_read(unit,files_read(file_index)%var(index_field)%field,data, level=level)
  endif  
  return
end subroutine read_data_text
!.............................................................. 
! </SUBROUTINE>
subroutine read_data_2d_new(filename,fieldname,data,domain,timelevel,&
      no_domain,position,tile_count)
  character(len=*), intent(in)                 :: filename, fieldname
  real, dimension(:,:), intent(inout)          :: data     !2 dimensional data 
  real, dimension(size(data,1),size(data,2),1) :: data_3d
  type(domain2d), intent(in), optional         :: domain
  integer, intent(in) , optional               :: timelevel
  logical, intent(in), optional                ::  no_domain
  integer, intent(in) , optional               :: position, tile_count
  integer                                      :: isc,iec,jsc,jec,isd,ied,jsd,jed
  integer :: isg,ieg,jsg,jeg
  integer                                      :: xsize_c,ysize_c,xsize_d,ysize_d
  integer                                      :: xsize_g,ysize_g, ishift, jshift

!#ifdef use_CRI_pointers
!  pointer( p, data_3d )
!  p = LOC(data)
!#endif
  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
        no_domain,position,tile_count)

  if(PRESENT(domain)) then
     call mpp_get_global_domain( domain,isg,ieg,jsg,jeg,xsize=xsize_g,ysize=ysize_g, tile_count=tile_count, position=position)
     call mpp_get_compute_domain( domain,isc,iec,jsc,jec,xsize=xsize_c,ysize=ysize_c, tile_count=tile_count, position=position)
     call mpp_get_data_domain( domain,isd,ied,jsd,jed,xsize=xsize_d,ysize=ysize_d, tile_count=tile_count, position=position)
     call mpp_get_domain_shift  (domain, ishift, jshift, position)
     if((size(data,1)==xsize_c) .and. (size(data,2)==ysize_c)) then !on_comp_domain
        data(:,:) = data_3d(:,:,1)
     else if((size(data,1)==xsize_d) .and. (size(data,2)==ysize_d)) then !on_data_domain
        data(isc-isd+1:iec-isd+1,jsc-jsd+1:jec-jsd+1) = data_3d(isc-isd+1:iec-isd+1,jsc-jsd+1:jec-jsd+1,1)
     else if((size(data,1)==xsize_g) .and. (size(data,2)==ysize_g)) then !on_global_domain
        data(:,:) = data_3d(:,:,1)
     else 
        call mpp_error(FATAL,'error in read_data_2d_new, field '//trim(fieldname)// &
                      ' in file '//trim(filename)//' data must be in compute or data domain')
     endif
  else     
     data(:,:) = data_3d(:,:,1)
  endif

end subroutine read_data_2d_new
!.....................................................................
subroutine read_data_1d_new(filename,fieldname,data,domain,timelevel,&
      no_domain, tile_count)
  character(len=*), intent(in)           :: filename, fieldname
  real, dimension(:), intent(inout)      :: data     !1 dimensional data 
  real, dimension(size(data,1),1,1)      :: data_3d
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in) , optional         :: timelevel
  logical, intent(in), optional          ::  no_domain
  integer, intent(in), optional          :: tile_count
#ifdef use_CRI_pointers
  pointer( p, data_3d )
  p = LOC(data)
#endif

  if(present(no_domain)) then
     if(.NOT. no_domain) call mpp_error(FATAL, 'fms_io(read_data_1d_new): no_domain should be true for field ' &
                                 //trim(fieldname)//' of file '//trim(filename) )
  end if

  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
        no_domain=.true., tile_count=tile_count)

end subroutine read_data_1d_new
!.....................................................................

subroutine read_data_scalar_new(filename,fieldname,data,domain,timelevel,&
      no_domain, tile_count)

! this subroutine is for reading a single number
  character(len=*), intent(in)           :: filename, fieldname
  real, intent(inout)                    :: data     !zero dimension data 
  real, dimension(1,1,1)                 :: data_3d
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in) , optional         :: timelevel
  logical, intent(in), optional          ::  no_domain
  integer, intent(in), optional          :: tile_count

  if(present(no_domain)) then
     if(.NOT. no_domain) call mpp_error(FATAL, 'fms_io(read_data_scalar_new): no_domain should be true for field ' &
                                 //trim(fieldname)//' of file '//trim(filename) )
  end if

  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
        no_domain=.true., tile_count=tile_count)

  data = data_3d(1,1,1)

end subroutine read_data_scalar_new
!.....................................................................

function unique_axes(file, index, id_axes, siz_axes, dom)
  type(restart_file_type),   intent(inout)           :: file
  integer,                      intent(in)           :: index
  integer, dimension(:),       intent(out)           :: id_axes
  integer, dimension(:),       intent(out)           :: siz_axes
  type(domain1d), dimension(:), intent(in), optional :: dom
  integer                                            :: unique_axes
  type(var_type), pointer, save :: cur_var => NULL()
  integer :: i,j
  logical :: found
  
  unique_axes=0

  if(index <0 .OR. index > 3) call mpp_error(FATAL,"unique_axes(fms_io_mod): index should be 1, 2 or 3")

  do i = 1, file%nvar
     cur_var => file%var(i)
     if(cur_var%ndim < index) cycle
     found = .false.
     do j = 1, unique_axes
        if(siz_axes(j) == cur_var%gsiz(index) ) then
           if(PRESENT(dom)) then
              if(cur_var%domain_idx == id_axes(j) ) then
                 found = .true.
                 exit
              else if(cur_var%domain_idx >0 .AND. id_axes(j) >0) then
                 if(dom(cur_var%domain_idx) .EQ. dom(id_axes(j)) ) then
                    found = .true.
                    exit
                 end if
              end if
           else
              found = .true.
              exit
           end if
        end if  
     end do
     if(found) then
        cur_var%id_axes(index) = j
     else
        unique_axes = unique_axes+1
        if(unique_axes > max_axes) then
           write(error_msg,'(I3,"/",I3)') unique_axes, max_axes
           if(index == 1 ) then
              call mpp_error(FATAL,'# x axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           else if(index == 2 ) then
              call mpp_error(FATAL,'# y axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           else
              call mpp_error(FATAL,'# z axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           end if
        endif
        id_axes(unique_axes)   = cur_var%domain_idx
        siz_axes(unique_axes) = cur_var%gsiz(index)        
        cur_var%id_axes(index) = unique_axes
     end if
  end do       

  cur_var => NULL()
  
  return

end function unique_axes

  !#######################################################################
  !#######################################################################
  !   --------- routines for reading distributed data ---------
  ! before calling these routines the domain decompostion must be set
  ! by calling "set_domain" with the appropriate domain2d data type
  !
  ! reading can be done either by all PEs (default) or by only the root PE
  ! this is controlled by namelist variable "read_all_pe".
  
  ! By default, array data is expected to be declared in data domain and no_halo
  !is NOT needed, however IF data is decalared in COMPUTE domain then optional NO_HALO should be .true.

  !#######################################################################

subroutine read_data_2d ( unit, data, end)

  integer, intent(in)                        :: unit
  real,    intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end  
  real, dimension(isg:ieg,jsg:jeg)           :: gdata
  integer                                    :: len
  logical                                    :: no_halo

  include "read_data_2d.inc"  
end subroutine read_data_2d

!#######################################################################

subroutine read_ldata_2d ( unit, data, end)

  integer, intent(in)                        :: unit
  logical, intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end  
  logical, dimension(isg:ieg,jsg:jeg)        :: gdata
  integer                                    :: len
  logical                                    :: no_halo

  include "read_data_2d.inc"
end subroutine read_ldata_2d
!#######################################################################

subroutine read_idata_2d ( unit, data, end)

  integer, intent(in)                        :: unit
  integer, intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end
  integer, dimension(isg:ieg,jsg:jeg)        :: gdata
  integer                                    :: len
  logical                                    :: no_halo

  include "read_data_2d.inc"
end subroutine read_idata_2d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine read_cdata_2d ( unit, data, end)
  
  integer, intent(in)                           :: unit
  complex,    intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional                :: end
  complex, dimension(isg:ieg,jsg:jeg)           :: gdata
  integer                                       :: len
  logical                                       :: no_halo

  include "read_data_2d.inc"
end subroutine read_cdata_2d
#endif

!#######################################################################

subroutine read_data_3d ( unit, data, end)

  integer, intent(in)                           :: unit
  real,    intent(out), dimension(isd:,jsd:,:)  :: data
  logical, intent(out), optional                :: end  
  real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
  integer                                       :: len
  logical                                       :: no_halo

  include "read_data_3d.inc"
end subroutine read_data_3d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine read_cdata_3d ( unit, data, end)

  integer, intent(in)                              :: unit
  complex, intent(out), dimension(isd:,jsd:,:)     :: data
  logical, intent(out), optional                   :: end
  complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
  integer                                          :: len
  logical                                          :: no_halo

  include "read_data_3d.inc"  
end subroutine read_cdata_3d
#endif

!#######################################################################

subroutine read_data_4d ( unit, data, end)

  integer, intent(in)                                        :: unit
  real,    intent(out), dimension(isd:,jsd:,:,:)             :: data
  logical, intent(out), optional                             :: end
  real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer                                                    :: len
  logical                                                    :: no_halo
! WARNING: memory usage with this routine could be costly
   
  include "read_data_4d.inc"  
end subroutine read_data_4d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine read_cdata_4d ( unit, data, end)

  integer, intent(in)                                           :: unit
  complex, intent(out), dimension(isd:,jsd:,:,:)                :: data
  logical, intent(out), optional                                :: end
  complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer                                                       :: len
  logical                                                       :: no_halo
! WARNING: memory usage with this routine could be costly
 
  include "read_data_4d.inc"
end subroutine read_cdata_4d
#endif

!#######################################################################
!     -------- routines for writing distributed data --------
! before calling these routines the domain decompostion must be set
! by calling "set_domain" with the appropriate domain2d data type
!#######################################################################
subroutine write_data_2d ( unit, data )
  integer, intent(in)                       :: unit
  real,    intent(in), dimension(isd:,jsd:) :: data  
  real, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_data_2d

!#######################################################################

subroutine write_ldata_2d ( unit, data )

  integer, intent(in)                       :: unit
  logical, intent(in), dimension(isd:,jsd:) :: data
  logical, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_ldata_2d

!#######################################################################
subroutine write_idata_2d ( unit, data )

  integer, intent(in)                       :: unit
  integer, intent(in), dimension(isd:,jsd:) :: data
  integer, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_idata_2d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine write_cdata_2d ( unit, data )

  integer, intent(in)                       :: unit
  complex, intent(in), dimension(isd:,jsd:) :: data
  complex, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_cdata_2d
#endif

!#######################################################################

subroutine write_data_3d ( unit, data )

  integer, intent(in) :: unit
  real,    intent(in), dimension(isd:,jsd:,:) :: data
  real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
    
  include "write_data.inc"
end subroutine write_data_3d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine write_cdata_3d ( unit, data )

  integer, intent(in) :: unit
  complex, intent(in), dimension(isd:,jsd:,:) :: data
  complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata

  include "write_data.inc"
end subroutine write_cdata_3d
#endif

!#######################################################################
subroutine write_data_4d ( unit, data )

  integer, intent(in) :: unit
  real,    intent(in), dimension(isd:,jsd:,:,:) :: data
  real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer :: n

  if (.not.associated(Current_domain))  &
       call mpp_error(FATAL,'fms_io(write_data_4d): need to call set_domain ')

! get the global data and write only on root pe
! do this one field at a time to save memory
  do n = 1, size(data,4)
     call mpp_global_field ( Current_domain, data(:,:,:,n), gdata(:,:,:,n) )
  enddo
  if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata
end subroutine write_data_4d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine write_cdata_4d ( unit, data )

  integer, intent(in) :: unit
  complex,    intent(in), dimension(isd:,jsd:,:,:) :: data
  complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer :: n

  if (.not.associated(Current_domain)) call mpp_error(FATAL, 'fms_io(write_cdata_4d): need to call set_domain')

! get the global data and write only on root pe
! do this one field at a time to save memory
  do n = 1, size(data,4)
     call mpp_global_field ( Current_domain, data(:,:,:,n), gdata(:,:,:,n) )
  enddo
  if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata
end subroutine write_cdata_4d
#endif

!#######################################################################
! private routines (read_eof,do_read)
! this routine is called when an EOF is found while
! reading a distributed data file using read_data

subroutine read_eof (end_found)
  logical, intent(out), optional :: end_found
  
  if (present(end_found))then
     end_found = .true.
  else
     call mpp_error(FATAL,'fms_io(read_eof): unexpected EOF')
  endif
end subroutine read_eof

!#######################################################################
! determines if current pe should read data
! checks namelist variable read_all_pe

function do_read ( )
  logical :: do_read
  do_read = mpp_pe() == mpp_root_pe() .or. read_all_pe
end function do_read

!!#######################################################################

subroutine reset_field_name(fileObj, id_field, name)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  character(len=*),        intent(in)         :: name

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_name): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_name): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )

  fileObj%var(id_field)%name = trim(name)

end subroutine reset_field_name 

!#######################################################################

subroutine reset_field_pointer_r0d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real,                    intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r0d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r0d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r0d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p0dr(1, id_field)%p => data

end subroutine reset_field_pointer_r0d

!#######################################################################

subroutine reset_field_pointer_r1d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:),      intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r1d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r1d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r1d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p1dr(1, id_field)%p => data

end subroutine reset_field_pointer_r1d


!#######################################################################
subroutine reset_field_pointer_r2d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:,:),    intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r2d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r2d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r2d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p2dr(1, id_field)%p => data

end subroutine reset_field_pointer_r2d

!#######################################################################

subroutine reset_field_pointer_r3d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:,:,:),  intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r3d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r3d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r3d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p3dr(1, id_field)%p => data

end subroutine reset_field_pointer_r3d

!#######################################################################

subroutine reset_field_pointer_i0d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer,                 intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i0d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i0d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i0d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p0di(1, id_field)%p => data

end subroutine reset_field_pointer_i0d

!#######################################################################

subroutine reset_field_pointer_i1d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer, dimension(:),   intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i1d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i1d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i1d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p1di(1, id_field)%p => data

end subroutine reset_field_pointer_i1d


!#######################################################################
subroutine reset_field_pointer_i2d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer, dimension(:,:), intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i2d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i2d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i2d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p2di(1, id_field)%p => data

end subroutine reset_field_pointer_i2d

!#######################################################################

subroutine reset_field_pointer_i3d(fileObj, id_field, data)
  type(restart_file_type),   intent(inout)      :: fileObj
  integer,                   intent(in)         :: id_field
  integer, dimension(:,:,:), intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i3d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i3d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i3d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p3di(1, id_field)%p => data

end subroutine reset_field_pointer_i3d

!#######################################################################

subroutine reset_field_pointer_r0d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real,                    intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r0d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r0d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r0d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p0dr(1, id_field)%p => data1
  fileObj%p0dr(2, id_field)%p => data2

end subroutine reset_field_pointer_r0d_2level

!#######################################################################

subroutine reset_field_pointer_r1d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:),      intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r1d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r1d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r1d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p1dr(1, id_field)%p => data1
  fileObj%p1dr(2, id_field)%p => data2

end subroutine reset_field_pointer_r1d_2level

!#######################################################################

subroutine reset_field_pointer_r2d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:,:),    intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r2d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r2d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r2d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p2dr(1, id_field)%p => data1
  fileObj%p2dr(2, id_field)%p => data2

end subroutine reset_field_pointer_r2d_2level

!#######################################################################

subroutine reset_field_pointer_r3d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:,:,:),  intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r3d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r3d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r3d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p3dr(1, id_field)%p => data1
  fileObj%p3dr(2, id_field)%p => data2

end subroutine reset_field_pointer_r3d_2level

!#######################################################################

subroutine reset_field_pointer_i0d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer,                 intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i0d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i0d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i0d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p0di(1, id_field)%p => data1
  fileObj%p0di(2, id_field)%p => data2

end subroutine reset_field_pointer_i0d_2level

!#######################################################################

subroutine reset_field_pointer_i1d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer, dimension(:),   intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i1d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i1d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i1d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p1di(1, id_field)%p => data1
  fileObj%p1di(2, id_field)%p => data2

end subroutine reset_field_pointer_i1d_2level

!#######################################################################

subroutine reset_field_pointer_i2d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer, dimension(:,:), intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i2d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i2d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i2d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p2di(1, id_field)%p => data1
  fileObj%p2di(2, id_field)%p => data2

end subroutine reset_field_pointer_i2d_2level

!#######################################################################

subroutine reset_field_pointer_i3d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type),   intent(inout)      :: fileObj
  integer,                   intent(in)         :: id_field
  integer, dimension(:,:,:), intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i3d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i3d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i3d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p3di(1, id_field)%p => data1
  fileObj%p3di(2, id_field)%p => data2

end subroutine reset_field_pointer_i3d_2level

!#########################################################################
!   This function returns .true. if the field referred to by id has
! initialized from a restart file, and .false. otherwise. 
!
! Arguments: id - A integer that is the index of the field in fileObj.
!  (in)  fileObj - The control structure returned by a previous call to
!                  register_restart_field
function query_initialized_id(fileObj, id)
  type(restart_file_type), intent(in) :: fileObj
  integer,                 intent(in) :: id

  logical :: query_initialized_id

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_id): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id < 1 .OR. id > fileObj%nvar) call mpp_error(FATAL, "fms_io(query_initialized_id): " // &
      "argument id must be between 1 and nvar in the restart_file_type object")

  query_initialized_id = fileObj%var(id)%initialized

  return

end function query_initialized_id

!#########################################################################
!   This function returns .true. if the field referred to by name has
! initialized from a restart file, and .false. otherwise. 
!
! Arguments: name - A pointer to the field that is being queried.
!  (in)  fileObj - The control structure returned by a previous call to
!                  register_restart_field
function query_initialized_name(fileObj, name)
  type(restart_file_type), intent(inout) :: fileObj
  character(len=*),           intent(in) :: name

  logical :: query_initialized_name

  integer :: m

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_name): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  query_initialized_name = .false.
  do m=1,fileObj%nvar
    if (trim(name) == fileObj%var(m)%name) then
      if (fileObj%var(m)%initialized) query_initialized_name = .true.
      exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (m<=fileObj%nvar) then
     fileObj%var(m)%initialized = .true.
  else if(mpp_pe() == mpp_root_pe()) then
    call mpp_error(NOTE,"fms_io(query_initialized_name): Unknown restart variable "//name// &
                        " queried for initialization.")
  end if

end function query_initialized_name


!   This function returns 1 if the field pointed to by f_ptr has
! initialized from a restart file, and 0 otherwise.  If f_ptr is
! NULL, it tests whether the entire restart file has been success-
! fully read.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
function query_initialized_r2d(fileObj, f_ptr, name)
  type(restart_file_type),   intent(inout) :: fileObj 
  real, dimension(:,:), target, intent(in) :: f_ptr
  character(len=*),             intent(in) :: name

  logical :: query_initialized_r2d
  integer :: m

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_r2d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  query_initialized_r2d = .false.
  do m=1, fileObj%nvar
     if (ASSOCIATED(fileObj%p2dr(1,m)%p,f_ptr)) then
        if (fileObj%var(m)%initialized) query_initialized_r2d = .true.
        exit
     endif
  enddo
  ! Assume that you are going to initialize it now, so set flag to initialized if
  ! queried again.
  if (m<=fileObj%nvar) then
     fileObj%var(m)%initialized = .true.
  else
     query_initialized_r2d = query_initialized_name(fileObj, name)
     if (mpp_pe() == mpp_root_pe() ) call mpp_error(NOTE, "fms_io(query_initialized_r2d): Unable to find "// &
          trim(name)//" queried by pointer, "//"probably because of the suspect comparison of pointers by ASSOCIATED.")
     query_initialized_r2d = query_initialized_name(fileObj, name)
     if (mpp_pe() == mpp_root_pe() .AND. query_initialized_r2d) call mpp_error(NOTE, &
          "fms_io(query_initialized_r2d): "//trim(name)// " initialization confirmed by name.")
  endif

  return

end function query_initialized_r2d


!#######################################################################
!#######################################################################
!
! routines for opening specific types of files:
!
!                       form        action 
! open_namelist_file  MPP_ASCII   MPP_RDONLY  
! open restart_file   MPP_NATIVE
! open_ieee32_file    MPP_IEEE32
!
! all have: access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true.
! use the close_file interface to close these files
!
! if other types of files need to be opened the mpp_open and
! mpp_close routines in the mpp_io_mod should be used
!
!#######################################################################


! <FUNCTION NAME="open_namelist_file">
!   <DESCRIPTION>
! Opens single namelist file for reading only by all PEs
! the default file opened is called "input.nml".
!   </DESCRIPTION>
! <IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_namelist_file (file) result (unit)
  character(len=*), intent(in), optional :: file
  integer :: unit

  if (.not.module_is_initialized) call fms_io_init ( )    
  if (present(file)) then
     call mpp_open ( unit, file, form=MPP_ASCII, action=MPP_RDONLY, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
  else
     call mpp_open ( unit, 'input.nml', form=MPP_ASCII, action=MPP_RDONLY, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
  endif
end function open_namelist_file
! </FUNCTION>

! <FUNCTION NAME="open_restart_file">
!   <DESCRIPTION> 
! Opens single restart file for reading by all PEs or
! writing by root PE only
! the file has native format and no mpp header records.
!   </DESCRIPTION>
!<IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
!<IN NAME="action" TYPE="character">
! action to be performed: can be 'read' or 'write'
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_restart_file (file, action) result (unit)
  character(len=*), intent(in) :: file, action
  integer :: unit  
  integer :: mpp_action

  if (.not.module_is_initialized) call fms_io_init ( )

!   --- action (read,write) ---

  select case (lowercase(trim(action)))
  case ('read')
     mpp_action = MPP_RDONLY
  case ('write')
     mpp_action = MPP_OVERWR
  case default
     call mpp_error(FATAL,'fms_io(open_restart_file): action should be either read or write in file'//trim(file))
  end select
  
  call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
       access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

end function open_restart_file
! </FUNCTION>


! <FUNCTION NAME="open_direct_file">
!   <DESCRIPTION> 
! Opens single direct access file for reading by all PEs or
! writing by root PE only
! the file has native format and no mpp header records.
!   </DESCRIPTION>

  function open_direct_file (file, action, recl) result (unit)
    character(len=*), intent(in) :: file, action
    integer,          intent(in) :: recl
    integer :: unit

    integer :: mpp_action

    if (.not.module_is_initialized) call fms_io_init ( )

    !   --- action (read,write) ---

    select case (lowercase(trim(action)))
    case ('read')
       mpp_action = MPP_RDONLY
    case ('write')
       mpp_action = MPP_OVERWR
    case default
       call mpp_error(FATAL,'invalid option for argument action')
    end select

    call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
         access=MPP_DIRECT, threading=MPP_SINGLE, nohdrs=.true., recl=recl )

  end function open_direct_file
! </FUNCTION>

! <FUNCTION NAME=" open_ieee32_file">
!   <DESCRIPTION>  
! Opens single 32-bit ieee file for reading by all PEs or 
! writing by root PE only (writing is not recommended)
! the file has no mpp header records.
!   </DESCRIPTION>
!<IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
!<IN NAME="action" TYPE="character">
! action to be performed: can be 'read' or 'write'
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_ieee32_file (file, action) result (unit)
  character(len=*), intent(in) :: file, action
  integer :: unit  
  integer :: mpp_action

  if (.not.module_is_initialized) call fms_io_init ( )

!   --- action (read,write) ---
  select case (lowercase(trim(action)))
  case ('read')
     mpp_action = MPP_RDONLY
  case ('write')
     mpp_action = MPP_OVERWR
  case default
     call mpp_error (FATAL,'fms_io(open_ieee32_file): action should be either read or write in file'//trim(file))
  end select
  
  if (iospec_ieee32(1:1) == ' ') then
     call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
          nohdrs=.true. )
  else
     call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
          nohdrs=.true., iospec=iospec_ieee32 )
  endif
end function open_ieee32_file
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME=" close_file">
!   <DESCRIPTION>
!  Closes files that are opened by: open_namelist_file, open restart_file,
! and open_ieee32_file. Users should use mpp_close for other cases.
!   </DESCRIPTION>
!<IN NAME="unit" TYPE="integer">
! unit number of the file to be closed
! </IN>
!<IN NAME="status" TYPE="character, optional">
! action to be performed: can be 'delete'
! </IN>

subroutine close_file (unit, status)
  integer,          intent(in)           :: unit
  character(len=*), intent(in), optional :: status
  
  if (.not.module_is_initialized) call fms_io_init ( )
  if (unit == stdlog()) return    
  if (present(status)) then
     if (lowercase(trim(status)) == 'delete') then
        call mpp_close (unit, action=MPP_DELETE)
     else
        call mpp_error(FATAL,'fms_io(close_file): status should be DELETE')
     endif
  else
     call mpp_close (unit)
  endif
end subroutine close_file
! </FUNCTION>

!#######################################################################
  

! <SUBROUTINE NAME="set_domain">
!   <DESCRIPTION>
! set_domain is called to save the domain2d data type prior to
! calling the distributed data I/O routines, read_data and write_data.
!   </DESCRIPTION>
! <IN NAME="Domain2" TYPE="domain2D">
! domain to be passed to routines in fms_io_mod, Current_domain will point to
! this Domain2
! </IN>
subroutine set_domain (Domain2)
    
  type(domain2D), intent(in), target :: Domain2

  if (.NOT.module_is_initialized) call fms_io_init ( )

!  --- set_domain must be called before a read_data or write_data ---
  if (associated(Current_domain)) nullify (Current_domain)
  Current_domain => Domain2
  
  !  --- module indexing to shorten read/write routines ---
  
  call mpp_get_compute_domain (Current_domain,is ,ie ,js ,je )
  call mpp_get_data_domain    (Current_domain,isd,ied,jsd,jed)
  call mpp_get_global_domain  (Current_domain,isg,ieg,jsg,jeg)
end subroutine set_domain
!#######################################################################
! </SUBROUTINE>

! <SUBROUTINE NAME="nullify_domain">
subroutine nullify_domain ()
!   <DESCRIPTION>
! Use to nulify domain that has been assigned by set_domain.
!   </DESCRIPTION>
  if (.NOT.module_is_initialized) call fms_io_init ( )

!  --- set_domain must be called before a read_data or write_data ---

  if (associated(Current_domain)) nullify (Current_domain)
  is=0;ie=0;js=0;je=0
  isd=0;ied=0;jsd=0;jed=0
  isg=0;ieg=0;jsg=0;jeg=0
end subroutine nullify_domain
! </SUBROUTINE>

! <SUBROUTINE NAME="return_domain">
!   <DESCRIPTION>
! This routine is the reverse of set_domain above. This routine is called when 
! users want to retrieve the domain2d that is used in fms_io_mod
!   </DESCRIPTION>
! <OUT NAME="domain2" TYPE="domain2D">
! domain returned from  fms_io_mod.
! </OUT>
subroutine return_domain(domain2)
  type(domain2D), intent(inout) :: domain2

  if (associated(Current_domain)) then
     domain2 = Current_domain
  else
     domain2 = NULL_DOMAIN2D       
  endif
end subroutine return_domain
! </SUBROUTINE>

!#######################################################################
! this will be a private routine with the next release
! users should get the domain decomposition from the domain2d data type

!#######################################################################
! <SUBROUTINE NAME="get_domain_decomp">
!   <DESCRIPTION>
! This will be a private routine with the next release.
! Users should get the domain decomposition from the domain2d data type.
!   </DESCRIPTION>
! <OUT NAME="x" TYPE="integer">
! array containing beginning and ending indices of global and compute domain in x direction
! </OUT>
! <OUT NAME="y" TYPE="integer">
! array containing beginning and ending indices of global and compute domain in y direction
! </OUT>
subroutine get_domain_decomp ( x, y )

  integer, intent(out), dimension(4) :: x, y
  
  if (mpp_pe() == mpp_root_pe())  call mpp_error(NOTE, &
       'subroutine get_domain_decomp will be removed with the next release')
  x = (/ isg, ieg, is, ie /)
  y = (/ jsg, jeg, js, je /)

end subroutine get_domain_decomp
! </SUBROUTINE>

subroutine get_axis_cart(axis, cart)      

  type(axistype), intent(in) :: axis
  character(len=1), intent(out) :: cart
  character(len=1) :: axis_cart
  character(len=16), dimension(2) :: lon_names, lat_names
  character(len=16), dimension(3) :: z_names
  character(len=16), dimension(2) :: t_names
  character(len=16), dimension(2) :: lon_units, lat_units
  character(len=8) , dimension(4) :: z_units
  character(len=3) , dimension(4) :: t_units
  character(len=32) :: name
  integer :: i

  lon_names = (/'lon','x  '/)
  lat_names = (/'lat','y  '/)
  z_names = (/'depth ','height','z     '/)
  t_names = (/'time','t   '/)
  lon_units = (/'degrees_e   ', 'degrees_east'/)
  lat_units = (/'degrees_n    ', 'degrees_north'/)
  z_units = (/'cm ','m  ','pa ','hpa'/)
  t_units = (/'sec', 'min','hou','day'/)  
  call mpp_get_atts(axis,cartesian=axis_cart)
  cart = 'N'  
  if (axis_cart == 'x' ) cart = 'X'
  if (axis_cart == 'y' ) cart = 'Y'
  if (axis_cart == 'z' ) cart = 'Z'
  if (axis_cart == 't' ) cart = 'T'
  if (cart /= 'X' .and. cart /= 'Y' .and. cart /= 'Z' .and. cart /= 'T') then
     call mpp_get_atts(axis,name=name)
     name = lowercase(name)
     do i=1,size(lon_names(:))
        if (lowercase(name(1:3)) == trim(lon_names(i))) cart = 'X'
     enddo
     do i=1,size(lat_names(:))
        if (name(1:3) == trim(lat_names(i))) cart = 'Y'
     enddo
     do i=1,size(z_names(:))
        if (name == trim(z_names(i))) cart = 'Z'
     enddo
     do i=1,size(t_names(:))
        if (name(1:3) == t_names(i)) cart = 'T'
     enddo
  end if

  if (cart /= 'X' .and. cart /= 'Y' .and. cart /= 'Z' .and. cart /= 'T') then
     call mpp_get_atts(axis,units=name)
     name = lowercase(name)
     do i=1,size(lon_units(:))
        if (trim(name) == trim(lon_units(i))) cart = 'X'
     enddo
     do i=1,size(lat_units(:))
        if (trim(name) == trim(lat_units(i))) cart = 'Y'
     enddo
     do i=1,size(z_units(:))
        if (trim(name) == trim(z_units(i))) cart = 'Z'
     enddo
     do i=1,size(t_units(:))
        if (name(1:3) == trim(t_units(i))) cart = 'T'
     enddo
  end if
  
  return
end subroutine get_axis_cart


! The following function is here as a last resort.
! This is copied from what was utilities_mod in order that redundant code 
! could be deleted.

 function open_file ( file, form, action, access, threading, recl ) &
             result ( unit )

 character(len=*), intent(in) :: file 
 character(len=*), intent(in), optional :: form, action, access, threading
 integer         , intent(in), optional :: recl 
 integer  :: unit 

 character(len=32) :: form_local, action_local, access_local, thread_local
 character(len=32) :: action_ieee32
 logical :: open, no_headers, do_ieee32
 integer :: mpp_format, mpp_action, mpp_access, mpp_thread
!-----------------------------------------------------------------------

   if ( .not. module_is_initialized ) then
        call fms_io_init ( )
!        do_init = .false.
   endif

!   ---- return stdlog if this is the logfile ----

    if (trim(file) == 'logfile.out') then
       unit = stdlog()
       return
    endif

!   ---- is this file open and connected to a unit ?? ---- 

   inquire (file=trim(file), opened=open, number=unit)

!  cannot open a file that is already open
!  except for the log file

   if ( open .and. unit >= 0 ) then
      call mpp_error (FATAL, 'open_file in fms_mod : '// &
                       'file '//trim(file)//' is already open')
   endif   

!  --- defaults ---

   form_local   = 'formatted';  if (present(form))      form_local   = form
   access_local = 'sequential'; if (present(access))    access_local = access
   thread_local = 'single';     if (present(threading)) thread_local = threading
   no_headers   = .true.
   do_ieee32    = .false.

   if (present(action)) then    ! must be present
      action_local = action
   else
      call mpp_error (FATAL, 'open_file in fms_mod : argument action not present')
   endif


!   --- file format ---

    select case (lowercase(trim(form_local)))
       case ('formatted')
           mpp_format = MPP_ASCII
       case ('ascii')
           mpp_format = MPP_ASCII
       case ('unformatted')
           mpp_format = MPP_NATIVE
       case ('native')
           mpp_format = MPP_NATIVE
       case ('ieee32')
           do_ieee32 = .true.
       case ('netcdf')
           mpp_format = MPP_NETCDF
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument form')
    end select

!   --- action (read,write,append) ---

    select case (lowercase(trim(action_local)))
       case ('read')
           mpp_action = MPP_RDONLY
       case ('write')
           mpp_action = MPP_OVERWR
       case ('append')
           mpp_action = MPP_APPEND
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument action')
    end select

!   --- file access (sequential,direct) ---

    select case (lowercase(trim(access_local)))
       case ('sequential')
           mpp_access = MPP_SEQUENTIAL
       case ('direct')
           mpp_access = MPP_DIRECT
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument access')
    end select

!   --- threading (single,multi) ---

    select case (lowercase(trim(thread_local)))
       case ('single')
           mpp_thread = MPP_SINGLE
       case ('multi')
           mpp_thread = MPP_MULTI
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument thread')
           if (trim(file) /= '_read_error.nml') no_headers = .false.
    end select

!   ---------------- open file -----------------------

    if ( .not.do_ieee32 ) then
       call mpp_open ( unit, file, form=mpp_format, action=mpp_action, &
                       access=mpp_access, threading=mpp_thread,        &
                       nohdrs=no_headers, recl=recl )
    else
     ! special open for ieee32 file
     ! fms_mod has iospec value
     ! pass local action flag to open changing append to write
       action_ieee32 = action_local
       if (lowercase(trim(action_ieee32)) == 'append') action_ieee32 = 'write'
       unit = open_ieee32_file ( file, action_ieee32 )
    endif

!-----------------------------------------------------------------------

 end function open_file

  !#######################################################################

  function string_from_integer(n)
    integer, intent(in) :: n
    character(len=16) :: string_from_integer

    if(n<0) then
       call mpp_error(FATAL, 'fms_io_mod: n should be non-negative integer, contact developer')
    else if( n<10 ) then
       write(string_from_integer,'(i1)') n
    else if( n<100 ) then
       write(string_from_integer,'(i2)') n
    else if( n<1000 ) then
       write(string_from_integer,'(i3)') n
    else if( n<10000 ) then
       write(string_from_integer,'(i4)') n
    else if( n<100000 ) then
       write(string_from_integer,'(i5)') n
    else if( n<1000000 ) then
       write(string_from_integer,'(i6)') n
    else if( n<10000000 ) then
       write(string_from_integer,'(i7)') n
    else if( n<100000000 ) then
       write(string_from_integer,'(i8)') n
    else
       call mpp_error(FATAL, 'fms_io_mod: n is too big, contact developer')
    end if

    return

  end function string_from_integer

  !#######################################################################
  function string_from_real(a)
    real, intent(in) :: a
    character(len=32) :: string_from_real

    write(string_from_real,*) a

    return

  end function string_from_real

  !#######################################################################

 subroutine get_tile_string(str_out, str_in, tile, str2_in)
    character(len=*), intent(inout)        :: str_out
    character(len=*), intent(in)           :: str_in
    integer,          intent(in)           :: tile
    character(len=*), intent(in), optional :: str2_in

    if(tile > 0 .AND. tile < 9) then
       write(str_out,'(a,i1)') trim(str_in), tile
    else if(tile >= 10 .AND. tile < 99) then
       write(str_out,'(a,i2)') trim(str_in), tile
    else
       call mpp_error(FATAL, "FMS_IO: get_tile_string: tile must be a positive number less than 100")
    end if

    if(present(str2_in)) str_out=trim(str_out)//trim(str2_in)

 end subroutine get_tile_string


  !#####################################################################
  subroutine get_mosaic_tile_file(file_in, file_out, is_no_domain, domain, tile_count)
    character(len=*), intent(in)                   :: file_in
    character(len=*), intent(out)                  :: file_out
    logical,          intent(in)                   :: is_no_domain
    type(domain2D),   intent(in), optional, target :: domain
    integer,          intent(in), optional         :: tile_count 
    character(len=256)                             :: basefile, tilename
    integer                                        :: lens, ntiles, ntileMe, tile
    integer, dimension(:), allocatable             :: tile_id
    type(domain2d), pointer, save                  :: d_ptr =>NULL()

    !--- deal with the situation that the file is alreday in the full name.
    lens = len_trim(file_in)
    if(lens > 8) then
       if(file_in(lens-7:lens) == '.nc'//trim(pe_name) ) then
         file_out = file_in
         return
        endif
    endif

    if(index(file_in, '.nc', back=.true.)==0) then
       basefile = trim(file_in)
    else
       lens = len_trim(file_in)
       if(file_in(lens-2:lens) .NE. '.nc') call mpp_error(FATAL, &
            'fms_io_mod: .nc should be at the end of file '//trim(file_in))
       basefile = file_in(1:lens-3)
    end if

    if(mpp_mosaic_defined())then
       !--- get the tile name
       ntiles = 1
       if(PRESENT(domain))then
          ntiles = mpp_get_ntile_count(domain)
          d_ptr => domain
       elseif (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
          ntiles = mpp_get_ntile_count(Current_domain)
          d_ptr => Current_domain
       endif
       if(ntiles > 1 )then
          ntileMe = mpp_get_current_ntile(d_ptr)
          allocate(tile_id(ntileMe))
          tile_id = mpp_get_tile_id(d_ptr)
          tile = 1
          if(present(tile_count)) tile = tile_count
          tilename = 'tile'//string(tile_id(tile))
          deallocate(tile_id)
          if(index(basefile,'.'//trim(tilename),back=.true.) == 0)then
             basefile = trim(basefile)//'.'//trim(tilename);
          end if
       end if
    endif

    file_out = trim(basefile)//'.nc'

    d_ptr =>NULL()

  end subroutine get_mosaic_tile_file

  !#############################################################################
  subroutine get_mosaic_tile_grid(grid_file, mosaic_file, domain, tile_count) 
    character(len=*), intent(out)          :: grid_file
    character(len=*), intent(in)           :: mosaic_file
    type(domain2D),   intent(in)           :: domain
    integer,          intent(in), optional :: tile_count 
    integer                                :: tile, ntileMe
    integer, dimension(:), allocatable     :: tile_id

    tile = 1
    if(present(tile_count)) tile = tile_count
    ntileMe = mpp_get_current_ntile(domain)
    allocate(tile_id(ntileMe))
    tile_id = mpp_get_tile_id(domain)     
    call read_data(mosaic_file, "gridfiles", grid_file, level=tile_id(tile) )
    grid_file = 'INPUT/'//trim(grid_file)
    deallocate(tile_id)

  end subroutine get_mosaic_tile_grid

  subroutine get_var_att_value_text(file, varname, attname, attvalue)
    character(len=*), intent(in)    :: file
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: attname
    character(len=*), intent(inout) :: attvalue
    integer                         :: unit

    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
    call mpp_get_att_value(unit, varname, attname, attvalue)
    call mpp_close(unit)
 
    return

  end subroutine get_var_att_value_text

  !#############################################################################
  ! return false if the attribute is not find in the file.
  function get_global_att_value_text(file, att, attvalue)
    character(len=*), intent(in)    :: file
    character(len=*), intent(in)    :: att
    character(len=*), intent(inout) :: attvalue
    logical                         :: get_global_att_value_text
    integer                         :: unit, ndim, nvar, natt, ntime, i
    type(atttype), allocatable      :: global_atts(:)

    get_global_att_value_text = .false.
    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(global_atts(natt))
    call mpp_get_atts(unit,global_atts)
    do i=1,natt
       if( trim(mpp_get_att_name(global_atts(i))) == trim(att) ) then
          attvalue = trim(mpp_get_att_char(global_atts(i)))
          get_global_att_value_text = .true.
          exit
       end if
    end do
    deallocate(global_atts)

    return

  end function get_global_att_value_text

  !#############################################################################
  ! return false if the attribute is not find in the file.
  function get_global_att_value_real(file, att, attvalue)
    character(len=*), intent(in)    :: file
    character(len=*), intent(in)    :: att
    real,             intent(inout) :: attvalue
    logical                         :: get_global_att_value_real
    integer                         :: unit, ndim, nvar, natt, ntime, i
    type(atttype), allocatable      :: global_atts(:)

    get_global_att_value_real = .false.
    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(global_atts(natt))
    call mpp_get_atts(unit,global_atts)
    do i=1,natt
       if( trim(mpp_get_att_name(global_atts(i))) == trim(att) ) then
          attvalue = mpp_get_att_real_scalar(global_atts(i))
          get_global_att_value_real = .true.
          exit
       end if
    end do
    deallocate(global_atts)

    return

  end function get_global_att_value_real

  !#############################################################################
  ! This routine will get the actual file name, as well as if read_dist is true or false.
  ! return true if such file exist and return false if not.
  function get_file_name(orig_file, actual_file, read_dist, io_domain_exist, no_domain, domain, &
                           tile_count)
    character(len=*),                 intent(in) :: orig_file
    character(len=*),                intent(out) :: actual_file
    logical,                         intent(out) :: read_dist
    logical,                         intent(out) :: io_domain_exist
    logical,                optional, intent(in) :: no_domain
    type(domain2D), target, optional, intent(in) :: domain
    integer,                optional, intent(in) :: tile_count  
    logical                                      :: get_file_name

    type(domain2d), pointer, save :: d_ptr, io_domain
    logical                       :: fexist, is_no_domain
    integer                       :: tile_id(1)
    character(len=256)            :: fname

    is_no_domain=.false.
    if(PRESENT(no_domain)) is_no_domain = no_domain

    if(present(domain)) then
       d_ptr => domain
    elseif (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
       d_ptr => Current_domain
    endif

    fexist          = .false.
    read_dist       = .false.
    get_file_name   = .false.
    io_domain_exist = .false.

  !--- The file maybe not netcdf file, we just check the original file.
    if(index(orig_file, '.nc', back=.true.) == 0) then
       inquire (file=trim(orig_file), exist=fexist)
       if(fexist) then
          actual_file = orig_file
          get_file_name = .true.
          return
       endif
    endif 
    
    !JWD:  This is likely a temporary fix. Since fms_io needs to know tile_count, 
    !JWD:  I just don't see how the physics can remain "tile neutral"
    call get_mosaic_tile_file(orig_file, actual_file, is_no_domain, domain, tile_count)

    !--- check if the file is group redistribution.
    if(ASSOCIATED(d_ptr)) then
       io_domain => mpp_get_io_domain(d_ptr)
       if(associated(io_domain)) then
          tile_id = mpp_get_tile_id(io_domain)       
          if(mpp_npes()>10000) then
             write(fname, '(a,i6.6)' ) trim(actual_file)//'.', tile_id(1)
          else
             write(fname, '(a,i4.4)' ) trim(actual_file)//'.', tile_id(1)
          endif
          inquire (file=trim(fname), exist=fexist)
          if(fexist) io_domain_exist = .true.
       endif   
       io_domain=>NULL()
    endif 

    if(.not. fexist) inquire (file=trim(actual_file)//trim(pe_name), exist=fexist)

    if(fexist) then
       read_dist = .true.
       d_ptr => NULL()
       get_file_name = .true.
       return
    endif

    inquire (file=trim(actual_file), exist=fexist)
    if(fexist) then
       d_ptr => NULL()
       get_file_name = .true.
       return
    endif    

    !Perhaps the file has an ensemble instance appendix
    call get_instance_filename(actual_file, actual_file)
    inquire (file=trim(actual_file)//trim(pe_name), exist=fexist)  
    if(.not. fexist) inquire (file=trim(actual_file)//'.nc'//trim(pe_name), exist=fexist)     
    if(fexist) then
       read_dist = .true. 
       d_ptr => NULL()
       get_file_name = .true.
       return
    endif    
    inquire (file=trim(actual_file), exist=fexist)
          if(.not. fexist) inquire (file=trim(actual_file)//'.nc', exist=fexist)
    
    if(fexist) then
       d_ptr => NULL()
       get_file_name = .true.
       return
    endif

  end function get_file_name
  

  !#############################################################################
  subroutine get_file_unit(filename, unit, index_file, read_dist, io_domain_exist, domain )
    character(len=*),         intent(in) :: filename
    integer,                 intent(out) :: unit, index_file 
    logical,                  intent(in) :: read_dist, io_domain_exist
    type(domain2d), optional, intent(in) :: domain    

    logical  :: file_opened
    integer  :: i

    ! Need to check if filename has been opened or not
    file_opened=.false.
    do i=1,num_files_r
       if (files_read(i)%name == trim(filename))  then
          index_file = i
          unit = files_read(index_file)%unit
          return 
       endif
    enddo

    ! need to open the file now
    ! Increase num_files_r and set file_type 
    if(num_files_r == max_files_r) &  ! need to have bigger max_files_r
         call mpp_error(FATAL,'fms_io(get_file_unit): max_files_r exceeded, increase it via fms_io_nml')
    num_files_r=num_files_r + 1 
    if (read_dist .and. thread_r == MPP_SINGLE) then
       call mpp_error(FATAL,'fms_io(get_file_unit): single-threaded read from distributed fileset not allowed' &
            //'change threading_read to MULTI')
    endif
    if(read_dist) then
       if(io_domain_exist) then
          if(present(domain)) then
             call mpp_open(unit,filename,form=form,action=MPP_RDONLY,threading=MPP_MULTI, &
                fileset=MPP_MULTI, domain=domain)
          else if(ASSOCIATED(current_domain) ) then
             call mpp_open(unit,filename,form=form,action=MPP_RDONLY,threading=MPP_MULTI, &
                fileset=MPP_MULTI, domain=current_domain)
          else
             call mpp_error(FATAL,'fms_io(get_file_unit): when io_domain_exsit = .true., '// &
                   'either domain is present or current_domain is associated')
          endif
       else
          call mpp_open(unit,trim(filename),form=form,action=MPP_RDONLY,threading=thread_r, &
            fileset=MPP_MULTI)
       endif
    else
       call mpp_open(unit,trim(filename),form=form,action=MPP_RDONLY,threading=thread_r, &
            fileset=MPP_SINGLE)
    end if
    files_read(num_files_r)%name = trim(filename)
    allocate(files_read(num_files_r)%var (max_fields) )
    files_read(num_files_r)%nvar = 0
    index_file = num_files_r
    files_read(index_file)%unit = unit   

  end subroutine get_file_unit

  !#############################################################################
  subroutine get_field_id(unit, index_file, fieldname, index_field, is_no_domain, is_not_dim)
    integer,          intent(in) :: unit
    integer,          intent(in) :: index_file
    character(len=*), intent(in) :: fieldname
    integer,         intent(out) :: index_field
    logical,          intent(in) :: is_no_domain
    logical,          intent(in) :: is_not_dim

    character(len=128)                     :: name
    type(axistype),  dimension(max_axes)   :: axes
    type(fieldtype), dimension(max_fields) :: fields
    integer                                :: i, j, ndim, nvar, natt, var_dim
    integer                                :: siz_in(4)

    index_field = -1
    do j = 1, files_read(index_file)%nvar
       if (trim(files_read(index_file)%var(j)%name) == trim(fieldname)) then
          index_field = j
          return
       endif
    enddo

    !--- fieldname is not read, so need to get fieldname from file 
    files_read(index_file)%nvar = files_read(index_file)%nvar + 1
    if(files_read(index_file)%nvar > max_fields) then
       write(error_msg,'(I3,"/",I3)') files_read(index_file)%nvar, max_fields 
       call  mpp_error(FATAL,'fms_io(get_field_id): max_fields exceeded, needs increasing, nvar/max_fields=' &
            //trim(error_msg))
    endif
    call mpp_get_info(unit, ndim, nvar, natt, files_read(index_file)%max_ntime)
    if(files_read(index_file)%max_ntime < 1)  files_read(index_file)%max_ntime = 1
    if(nvar > max_fields) then
       write(error_msg,'(I3,"/",I3)') files_read(index_file)%nvar,max_fields
       call mpp_error(FATAL,'fms_io(get_field_id): max_fields too small needs increasing,nvar/max_fields=' &
            //trim(error_msg)//'in file'//trim(files_read(index_file)%name))
    endif
    call mpp_get_fields(unit, fields(1:nvar))     
    siz_in = 1
    index_field = files_read(index_file)%nvar
    files_read(index_file)%var(index_field)%is_dimvar = .false.

    do i=1, nvar
       call mpp_get_atts(fields(i),name=name,ndim=var_dim,siz=siz_in)
       if (lowercase(trim(name)) == lowercase(trim(fieldname))) then ! found the variable
          if(var_dim .lt.3) then
             do j=var_dim+1,3
                siz_in(j)=1
             enddo
          endif
          files_read(index_file)%var(index_field)%name    = fieldname
          files_read(index_file)%var(index_field)%field   = fields(i)
          files_read(index_file)%var(index_field)%siz(:)  = siz_in
          files_read(index_file)%var(index_field)%gsiz(:) = siz_in
          return
       endif
    enddo

    !--- the fieldname may be a dimension variable.
    if( .not. is_not_dim) then
       if (ndim > max_axes) then
          write(error_msg,'(I3,"/",I3)') ndim, max_axes 
          call  mpp_error(FATAL,'fms_io(get_field_id): max_axes exceeded, needs increasing, ndim/max_fields=' &
               //trim(error_msg)//' in file '//trim(files_read(index_file)%name))             
       endif
       call mpp_get_axes(unit, axes(1:ndim))
       do i=1,ndim
          call mpp_get_atts(axes(i), name=name, len = siz_in(1)) 
          if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
             if(.not. is_no_domain) call mpp_error(FATAL, &
                  'fms_io(get_field_id): the field is a dimension variable, no_domain should be true.')
             files_read(index_file)%var(index_field)%is_dimvar = .true.
             files_read(index_file)%var(index_field)%name      = fieldname
             files_read(index_file)%var(index_field)%axis      = axes(i)
             files_read(index_file)%var(index_field)%siz(:)    = siz_in
             files_read(index_file)%var(index_field)%gsiz(:)   = siz_in
             return
          endif
       enddo
    end if
    !--- the field is not in the file when reaching here.
    call mpp_error(FATAL, 'fms_io(get_field_id): field '//trim(fieldname)// &
                   ' NOT found in file '//trim(files_read(index_file)%name))

  end subroutine get_field_id

!#######################################################################
! check the existence of the given file name
! if the file_name string has zero length or the
! first character is blank return a false result
! <FUNCTION NAME="file_exist">

!   <OVERVIEW>
!     Checks the existence of a given file name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Checks the existence of the given file name.
!     If the file_name string has zero length or the
!     first character is blank return a false result.
!   </DESCRIPTION>
!   <TEMPLATE>
!     file_exist ( file_name )
!   </TEMPLATE>

!   <IN NAME="file_name"  TYPE="character" >
!     A file name (or path name) that is checked for existence.
!   </IN>
!   <OUT NAME=""  TYPE="logical" >
!     This function returns a logical result.  If file_name exists the result 
!     is true, otherwise false is returned.
!     If the length of character string "file_name" is zero or the first
!     character is blank, then the returned value will be false.
!     When reading a file, this function is often used in conjunction with
!     routine open_file.
!   </OUT>
!   <ERROR MSG="set_domain not called" STATUS="FATAL">
!     Before calling write_data you must first call set_domain with domain2d data 
!     type associated with the distributed data you are writing.
!   </ERROR>

 function file_exist (file_name, domain, no_domain)
  character(len=*), intent(in)         :: file_name
  type(domain2d), intent(in), optional :: domain
  logical,        intent(iN), optional :: no_domain  

  logical                              :: file_exist, is_no_domain
  character(len=256)                   :: fname
  logical                              :: read_dist, io_domain_exist

  is_no_domain = .false.
  if(present(no_domain)) is_no_domain = no_domain
   !--- to deal with mosaic file, in this case, the file is assumed to be in netcdf format
   file_exist = get_file_name(file_name, fname, read_dist, io_domain_exist, no_domain=is_no_domain, domain=domain)
   if(is_no_domain) return
   if(.not.file_exist) file_exist=get_file_name(file_name, fname, read_dist, io_domain_exist, no_domain=.true.)   

   return
 
 end function file_exist
! </FUNCTION>


!#######################################################################
! <FUNCTION NAME="field_exist">

!   <OVERVIEW>
!     check if a given field name exists in a given file name. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     check if a given field name exists in a given file name. 
!     If the field_name string has zero length or the
!     first character is blank return a false result.
!     if the file file_name don't exist, return a false result.
!   </DESCRIPTION>
!   <TEMPLATE>
!     field_exist ( file_name, field_name )
!   </TEMPLATE>

!   <IN NAME="file_name"  TYPE="character" >
!     A file name (or path name) that is checked for existence.
!   </IN>
!   <IN NAME="field_name"  TYPE="character" >
!     A field name that is checked for existence.
!   </IN>
!   <OUT NAME=""  TYPE="logical" >
!     This function returns a logical result.  If field exists in the 
!     file file_name, the result is true, otherwise false is returned.
!     If the length of character string "field_name" is zero or the first
!     character is blank, then the returned value will be false.
!     if the file file_name don't exist, return a false result.
!   </OUT>

 function field_exist (file_name, field_name, domain, no_domain)
  character(len=*),                 intent(in) :: file_name
  character(len=*),                 intent(in) :: field_name
  type(domain2d), intent(in), optional, target :: domain
  logical,       intent(in),  optional         :: no_domain
  logical                      :: field_exist, is_no_domain
  integer                      :: unit, ndim, nvar, natt, ntime, i, nfile
  character(len=64)            :: name
  type(fieldtype), allocatable :: fields(:)
  logical                      :: file_exist, read_dist, io_domain_exist
  character(len=256)           :: fname

   field_exist = .false.
   if (len_trim(field_name) == 0) return
   if (field_name(1:1) == ' ')    return

   is_no_domain = .false.
   if(present(no_domain)) is_no_domain = no_domain   

   file_exist=get_file_name(file_name, fname, read_dist, io_domain_exist, no_domain=is_no_domain, domain=domain)
   if(file_exist) then
      call get_file_unit(fname, unit, nfile, read_dist, io_domain_exist, domain=domain)
      call mpp_get_info(unit, ndim, nvar, natt, ntime)
      allocate(fields(nvar))
      call mpp_get_fields(unit,fields)

      do i=1, nvar
         call mpp_get_atts(fields(i),name=name)
         if(lowercase(trim(name)) == lowercase(trim(field_name))) field_exist = .true.
      enddo
      deallocate(fields)
    endif
    if(field_exist .or. is_no_domain) return
    file_exist =  get_file_name(file_name, fname, read_dist, io_domain_exist, no_domain=.true.)
    if(file_exist) then
       call get_file_unit(fname, unit, nfile, read_dist, io_domain_exist)
       call mpp_get_info(unit, ndim, nvar, natt, ntime)
       allocate(fields(nvar))
       call mpp_get_fields(unit,fields)
       do i=1, nvar
          call mpp_get_atts(fields(i),name=name)
          if(lowercase(trim(name)) == lowercase(trim(field_name))) field_exist = .true.
       enddo
       deallocate(fields)
    endif

    return

 end function field_exist
! </FUNCTION>

subroutine set_filename_appendix(string_in)
  character(len=*) , intent(in) :: string_in
  filename_appendix = trim(string_in)
end subroutine set_filename_appendix

subroutine get_instance_filename(name_in,name_out)
  character(len=*)  , intent(in)  :: name_in
  character(len=*), intent(inout) :: name_out
  integer :: length
  
  length = len_trim(name_in)
  name_out = name_in(1:length)
  
  if(len_trim(filename_appendix) > 0) then
     if(name_in(length-2:length) == '.nc') then
        name_out = name_in(1:length-3)//'.'//trim(filename_appendix)//'.nc'
     else
        name_out = name_in(1:length)  //'.'//trim(filename_appendix)
     end if
  end if
  
end subroutine get_instance_filename

end module fms_io_mod


