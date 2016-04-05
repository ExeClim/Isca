!-----------------------------------------------------------------------
!                 Parallel I/O for message-passing codes
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------

! <CONTACT EMAIL="vb@gfdl.noaa.gov">
!   V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp_io.html"/>

! <OVERVIEW>
!   <TT>mpp_io_mod</TT>, is a set of simple calls for parallel I/O on
!   distributed systems. It is geared toward the writing of data in netCDF
!   format. It requires the modules <LINK
!   SRC="mpp_domains.html">mpp_domains_mod</LINK> and <LINK
!   SRC="mpp.html">mpp_mod</LINK>, upon which it is built.
! </OVERVIEW>

! <DESCRIPTION>
!   In massively parallel environments, an often difficult problem is
!   the reading and writing of data to files on disk. MPI-IO and MPI-2 IO
!   are moving toward providing this capability, but are currently not
!   widely implemented. Further, it is a rather abstruse
!   API. <TT>mpp_io_mod</TT> is an attempt at a simple API encompassing a
!   certain variety of the I/O tasks that will be required. It does not
!   attempt to be an all-encompassing standard such as MPI, however, it
!   can be implemented in MPI if so desired. It is equally simple to add
!   parallel I/O capability to <TT>mpp_io_mod</TT> based on vendor-specific
!   APIs while providing a layer of insulation for user codes.
!   
!   The <TT>mpp_io_mod</TT> parallel I/O API built on top of the <LINK
!   SRC="mpp_domains.html">mpp_domains_mod</LINK> and <LINK
!   SRC="mpp.html">mpp_mod</LINK> API for domain decomposition and
!   message passing. Features of <TT>mpp_io_mod</TT> include:
!   
!    1) Simple, minimal API, with free access to underlying API for more
!   complicated stuff.<BR/>
!    2) Self-describing files: comprehensive header information
!   (metadata) in the file itself.<BR/>
!    3) Strong focus on performance of parallel write: the climate models
!   for which it is designed typically read a minimal amount of data
!   (typically at the beginning of the run); but on the other hand, tend
!   to write copious amounts of data during the run. An interface for
!   reading is also supplied, but its performance has not yet been optimized.<BR/>
!    4) Integrated netCDF capability: <LINK SRC
!   ="http://www.unidata.ucar.edu/packages/netcdf/">netCDF</LINK> is a
!   data format widely used in the climate/weather modeling
!   community. netCDF is considered the principal medium of data storage
!   for <TT>mpp_io_mod</TT>. But I provide a raw unformatted
!   fortran I/O capability in case netCDF is not an option, either due to
!   unavailability, inappropriateness, or poor performance.<BR/>
!    5) May require off-line post-processing: a tool for this purpose,
!   <TT>mppnccombine</TT>, is available. GFDL users may use
!   <TT>~hnv/pub/mppnccombine</TT>. Outside users may obtain the
!   source <LINK SRC
!   ="ftp://ftp.gfdl.gov/perm/hnv/mpp/mppnccombine.c">here</LINK>.  It
!   can be compiled on any C compiler and linked with the netCDF
!   library. The program is free and is covered by the <LINK SRC
!   ="ftp://ftp.gfdl.gov/perm/hnv/mpp/LICENSE">GPL license</LINK>.
!   
!   The internal representation of the data being written out is
!   assumed be the default real type, which can be 4 or 8-byte. Time data
!   is always written as 8-bytes to avoid overflow on climatic time scales
!   in units of seconds.
!   
!   <LINK SRC="modes"></LINK><H4>I/O modes in <TT>mpp_io_mod</TT></H4>
!   
!   The I/O activity critical to performance in the models for which
!   <TT>mpp_io_mod</TT> is designed is typically the writing of large
!   datasets on a model grid volume produced at intervals during
!   a run. Consider a 3D grid volume, where model arrays are stored as
!   <TT>(i,j,k)</TT>. The domain decomposition is typically along
!   <TT>i</TT> or <TT>j</TT>: thus to store data to disk as a global
!   volume, the distributed chunks of data have to be seen as
!   non-contiguous. If we attempt to have all PEs write this data into a
!   single file, performance can be seriously compromised because of the
!   data reordering that will be required. Possible options are to have
!   one PE acquire all the data and write it out, or to have all the PEs
!   write independent files, which are recombined offline. These three
!   modes of operation are described in the <TT>mpp_io_mod</TT> terminology
!   in terms of two parameters, <I>threading</I> and <I>fileset</I>,
!   as follows:
!   
!   <I>Single-threaded I/O:</I> a single PE acquires all the data
!   and writes it out.<BR/>
!   <I>Multi-threaded, single-fileset I/O:</I> many PEs write to a
!   single file.<BR/>
!    <I>Multi-threaded, multi-fileset I/O:</I> many PEs write to
!   independent files. This is also called <I>distributed I/O</I>.
!   
!   The middle option is the most difficult to achieve performance. The
!   choice of one of these modes is made when a file is opened for I/O, in
!   <LINK SRC="#mpp_open">mpp_open</LINK>.
!   
!   <LINK name="metadata"></LINK><H4>Metadata in <TT>mpp_io_mod</TT></H4>
!   
!   A requirement of the design of <TT>mpp_io_mod</TT> is that the file must
!   be entirely self-describing: comprehensive header information
!   describing its contents is present in the header of every file. The
!   header information follows the model of netCDF. Variables in the file
!   are divided into <I>axes</I> and <I>fields</I>. An axis describes a
!   co-ordinate variable, e.g <TT>x,y,z,t</TT>. A field consists of data in
!   the space described by the axes. An axis is described in
!   <TT>mpp_io_mod</TT> using the defined type <TT>axistype</TT>:
!   
!   <PRE>
!   type, public :: axistype
!      sequence
!      character(len=128) :: name
!      character(len=128) :: units
!      character(len=256) :: longname
!      character(len=8) :: cartesian
!      integer :: len
!      integer :: sense           !+/-1, depth or height?
!      type(domain1D), pointer :: domain
!      real, dimension(:), pointer :: data
!      integer :: id, did
!      integer :: type  ! external NetCDF type format for axis data
!      integer :: natt
!      type(atttype), pointer :: Att(:) ! axis attributes
!   end type axistype
!   </PRE>
!   
!   A field is described using the type <TT>fieldtype</TT>:
!   
!   <PRE>
!   type, public :: fieldtype
!      sequence
!      character(len=128) :: name
!      character(len=128) :: units
!      character(len=256) :: longname
!      real :: min, max, missing, fill, scale, add
!      integer :: pack
!      type(axistype), dimension(:), pointer :: axes
!      integer, dimension(:), pointer :: size
!      integer :: time_axis_index
!      integer :: id
!      integer :: type ! external NetCDF format for field data
!      integer :: natt, ndim
!      type(atttype), pointer :: Att(:) ! field metadata
!   end type fieldtype
!   </PRE>
!   
!   An attribute (global, field or axis) is described using the <TT>atttype</TT>:
!   
!   <PRE>
!   type, public :: atttype
!      sequence
!      integer :: type, len
!      character(len=128) :: name
!      character(len=256)  :: catt
!      real(FLOAT_KIND), pointer :: fatt(:)
!   end type atttype
!   </PRE>
!   
!   <LINK name="packing"></LINK>This default set of field attributes corresponds
!   closely to various conventions established for netCDF files. The
!   <TT>pack</TT> attribute of a field defines whether or not a
!   field is to be packed on output. Allowed values of
!   <TT>pack</TT> are 1,2,4 and 8. The value of
!   <TT>pack</TT> is the number of variables written into 8
!   bytes. In typical use, we write 4-byte reals to netCDF output; thus
!   the default value of <TT>pack</TT> is 2. For
!   <TT>pack</TT> = 4 or 8, packing uses a simple-minded linear
!   scaling scheme using the <TT>scale</TT> and <TT>add</TT>
!   attributes. There is thus likely to be a significant loss of dynamic
!   range with packing. When a field is declared to be packed, the
!   <TT>missing</TT> and <TT>fill</TT> attributes, if
!   supplied, are packed also.
!   
!   Please note that the pack values are the same even if the default
!   real is 4 bytes, i.e <TT>PACK=1</TT> still follows the definition
!   above and writes out 8 bytes.
!   
!   A set of <I>attributes</I> for each variable is also available. The
!   variable definitions and attribute information is written/read by calling
!   <LINK SRC="#mpp_write_meta">mpp_write_meta</LINK> or <LINK SRC="#mpp_read_meta">mpp_read_meta</LINK>. A typical calling
!   sequence for writing data might be:
!   
!   <PRE>
!   ...
!     type(domain2D), dimension(:), allocatable, target :: domain
!     type(fieldtype) :: field
!     type(axistype) :: x, y, z, t
!   ...
!     call mpp_define_domains( (/1,nx,1,ny/), domain )
!     allocate( a(domain(pe)%x%data%start_index:domain(pe)%x%data%end_index, &
!                 domain(pe)%y%data%start_index:domain(pe)%y%data%end_index,nz) )
!   ...
!     call mpp_write_meta( unit, x, 'X', 'km', 'X distance', &
!          domain=domain(pe)%x, data=(/(float(i),i=1,nx)/) )
!     call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', &
!          domain=domain(pe)%y, data=(/(float(i),i=1,ny)/) )
!     call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', &
!          data=(/(float(i),i=1,nz)/) )
!     call mpp_write_meta( unit, t, 'Time', 'second', 'Time' )
!   
!     call mpp_write_meta( unit, field, (/x,y,z,t/), 'a', '(m/s)', AAA', &
!          missing=-1e36 )
!   ...
!     call mpp_write( unit, x )
!     call mpp_write( unit, y )
!     call mpp_write( unit, z )
!   ...
!   </PRE>
!   
!   In this example, <TT>x</TT> and <TT>y</TT> have been
!   declared as distributed axes, since a domain decomposition has been
!   associated. <TT>z</TT> and <TT>t</TT> are undistributed
!   axes. <TT>t</TT> is known to be a <I>record</I> axis (netCDF
!   terminology) since we do not allocate the <TT>data</TT> element
!   of the <TT>axistype</TT>. <I>Only one record axis may be
!   associated with a file.</I> The call to <LINK
!   SRC="#mpp_write_meta">mpp_write_meta</LINK> initializes
!   the axes, and associates a unique variable ID with each axis. The call
!   to <TT>mpp_write_meta</TT> with argument <TT>field</TT>
!   declared <TT>field</TT> to be a 4D variable that is a function
!   of <TT>(x,y,z,t)</TT>, and a unique variable ID is associated
!   with it. A 3D field will be written at each call to
!   <TT>mpp_write(field)</TT>.
!   
!   The data to any variable, including axes, is written by
!   <TT>mpp_write</TT>.
!   
!   Any additional attributes of variables can be added through
!   subsequent <TT>mpp_write_meta</TT> calls, using the variable ID as a
!   handle. <I>Global</I> attributes, associated with the dataset as a
!   whole, can also be written thus. See the <LINK
!   SRC="#mpp_write_meta">mpp_write_meta</LINK> call syntax below
!   for further details.
!   
!   You cannot interleave calls to <TT>mpp_write</TT> and
!   <TT>mpp_write_meta</TT>: the first call to
!   <TT>mpp_write</TT> implies that metadata specification is
!   complete.
!   
!   A typical calling sequence for reading data might be:
!   
!   <PRE>
!   ...
!     integer :: unit, natt, nvar, ntime
!     type(domain2D), dimension(:), allocatable, target :: domain
!     type(fieldtype), allocatable, dimension(:) :: fields
!     type(atttype), allocatable, dimension(:) :: global_atts
!     real, allocatable, dimension(:) :: times
!   ...
!     call mpp_define_domains( (/1,nx,1,ny/), domain )
!   
!     call mpp_read_meta(unit)
!     call mpp_get_info(unit,natt,nvar,ntime)
!     allocate(global_atts(natt))
!     call mpp_get_atts(unit,global_atts)
!     allocate(fields(nvar))
!     call mpp_get_vars(unit, fields)
!     allocate(times(ntime))
!     call mpp_get_times(unit, times)
!   
!     allocate( a(domain(pe)%x%data%start_index:domain(pe)%x%data%end_index, &
!                 domain(pe)%y%data%start_index:domain(pe)%y%data%end_index,nz) )
!   ...
!     do i=1, nvar
!       if (fields(i)%name == 'a')  call mpp_read(unit,fields(i),domain(pe), a,
!                                                 tindex)
!     enddo
!   ...
!   </PRE>
!   
!   In this example, the data are distributed as in the previous
!   example. The call to <LINK
!   SRC="#mpp_read_meta">mpp_read_meta</LINK> initializes
!   all of the metadata associated with the file, including global
!   attributes, variable attributes and non-record dimension data. The
!   call to <TT>mpp_get_info</TT> returns the number of global
!   attributes (<TT>natt</TT>), variables (<TT>nvar</TT>) and
!   time levels (<TT>ntime</TT>) associated with the file
!   identified by a unique ID (<TT>unit</TT>).
!   <TT>mpp_get_atts</TT> returns all global attributes for
!   the file in the derived type <TT>atttype(natt)</TT>.
!   <TT>mpp_get_vars</TT> returns variable types
!   (<TT>fieldtype(nvar)</TT>).  Since the record dimension data are not allocated for calls to <LINK SRC="#mpp_write">mpp_write</LINK>, a separate call to  <TT>mpp_get_times</TT> is required to access record dimension data.  Subsequent calls to
!   <TT>mpp_read</TT> return the field data arrays corresponding to
!   the fieldtype.  The <TT>domain</TT> type is an optional
!   argument.  If <TT>domain</TT> is omitted, the incoming field
!   array should be dimensioned for the global domain, otherwise, the
!   field data is assigned to the computational domain of a local array.
!   
!   <I>Multi-fileset</I> reads are not supported with <TT>mpp_read</TT>.

! </DESCRIPTION>

module mpp_io_mod
#include <fms_platform.h>

use mpp_parameter_mod,  only : MPP_WRONLY, MPP_RDONLY, MPP_APPEND, MPP_OVERWR, MPP_ASCII
use mpp_parameter_mod,  only : MPP_IEEE32, MPP_NATIVE, MPP_NETCDF, MPP_SEQUENTIAL
use mpp_parameter_mod,  only : MPP_DIRECT, MPP_SINGLE, MPP_MULTI, MPP_DELETE, MPP_COLLECT
use mpp_parameter_mod,  only : MPP_DEBUG, MPP_VERBOSE, NULLUNIT, NULLTIME, ALL_PES
use mpp_parameter_mod,  only : CENTER, EAST, NORTH, CORNER
use mpp_parameter_mod,  only : MAX_FILE_SIZE, GLOBAL_ROOT_ONLY, XUPDATE, YUPDATE
use mpp_mod,            only : mpp_error, FATAL, WARNING, NOTE, stdin, stdout, stderr, stdlog
use mpp_mod,            only : mpp_pe, mpp_root_pe, mpp_npes, lowercase, mpp_transmit
use mpp_mod,            only : mpp_init, mpp_sync, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,            only : MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, CLOCK_ROUTINE
use mpp_domains_mod,    only : domain1d, domain2d, NULL_DOMAIN1D, mpp_domains_init
use mpp_domains_mod,    only : mpp_get_global_domain, mpp_get_compute_domain
use mpp_domains_mod,    only :  mpp_get_data_domain, mpp_get_memory_domain
use mpp_domains_mod,    only : mpp_update_domains, mpp_global_field, mpp_domain_is_symmetry
use mpp_domains_mod,    only : operator( .NE. ), mpp_get_domain_shift
use mpp_domains_mod,    only : mpp_get_io_domain, mpp_domain_is_tile_root_pe, mpp_get_domain_tile_root_pe
use mpp_domains_mod,    only : mpp_get_tile_id, mpp_get_tile_npes, mpp_get_io_domain_layout
use mpp_domains_mod,    only : mpp_get_domain_name

implicit none
private

#ifdef use_netCDF
#include <netcdf.inc>
#endif

  !--- public parameters  -----------------------------------------------
  public :: MPP_WRONLY, MPP_RDONLY, MPP_APPEND, MPP_OVERWR, MPP_ASCII, MPP_IEEE32
  public :: MPP_NATIVE, MPP_NETCDF, MPP_SEQUENTIAL, MPP_DIRECT, MPP_SINGLE
  public :: MPP_MULTI, MPP_DELETE, MPP_COLLECT
  public :: FILE_TYPE_USED
  public :: MAX_FILE_SIZE
  !--- public data type ------------------------------------------------
  public :: axistype, atttype, fieldtype, validtype, filetype

  !--- public data -----------------------------------------------------
  public :: default_field, default_axis, default_att
    
  !--- public interface from mpp_io_util.h ----------------------
  public :: mpp_get_iospec, mpp_get_id, mpp_get_ncid, mpp_get_unit_range, mpp_is_valid
  public :: mpp_set_unit_range, mpp_get_info, mpp_get_atts, mpp_get_fields
  public :: mpp_get_times, mpp_get_axes, mpp_get_recdimid, mpp_get_axis_data
  public :: mpp_io_set_stack_size, mpp_get_field_index, mpp_get_axis_index
  public :: mpp_get_field_name, mpp_get_att_value, mpp_get_att_length
  public :: mpp_get_att_type, mpp_get_att_name, mpp_get_att_real, mpp_get_att_char
  public :: mpp_get_att_real_scalar
  public :: mpp_get_file_name, mpp_file_is_opened 
  public :: mpp_io_clock_on

  !--- public interface from mpp_io_misc.h ----------------------
  public :: mpp_io_init, mpp_io_exit, netcdf_err, mpp_flush

  !--- public interface from mpp_io_write.h ---------------------
  public :: mpp_write, mpp_write_meta, mpp_copy_meta, mpp_modify_meta

  !--- public interface from mpp_io_read.h ---------------------
  public :: mpp_read, mpp_read_meta, mpp_get_tavg_info

  !--- public interface from mpp_io_switch.h ---------------------
  public :: mpp_open, mpp_close

  !-----------------------------------------------------------------------------
  !--- mpp_io data types
  !-----------------------------------------------------------------------------
integer FILE_TYPE_USED  
type :: atttype
     private
     integer             :: type, len
     character(len=128)  :: name
     character(len=1280) :: catt
     real, pointer       :: fatt(:) =>NULL() ! just use type conversion for integers
  end type atttype

  type :: axistype
     private
     character(len=128) :: name
     character(len=128) :: units
     character(len=256) :: longname
     character(len=8)   :: cartesian
     character(len=24)  :: calendar
     integer            :: sense, len          !+/-1, depth or height?
     type(domain1D)     :: domain              !if pointer is associated, it is a distributed data axis
     real, pointer      :: data(:) =>NULL()    !axis values (not used if time axis)
     integer            :: id, did, type, natt !id is the "variable ID", did is the "dimension ID": 
                                               !netCDF requires 2 IDs for axes
     integer            :: shift               !normally is 0. when domain is symmetry, its value maybe 1.
     type(atttype), pointer :: Att(:) =>NULL()
  end type axistype

  type :: validtype
     private
     logical :: is_range ! if true, then the data represent the valid range
     real    :: min,max  ! boundaries of the valid range or missing value
  end type validtype

  type :: fieldtype
     private
     character(len=128)      :: name
     character(len=128)      :: units
     character(len=256)      :: longname
     character(len=128)      :: standard_name   ! CF standard name
     real                    :: min, max, missing, fill, scale, add
     integer                 :: pack
     type(axistype), pointer :: axes(:) =>NULL() !axes associated with field size, time_axis_index redundantly 
                                        !hold info already contained in axes. it's clunky and inelegant, 
                                        !but required so that axes can be shared among multiple files
     integer, pointer        :: size(:) =>NULL()
     integer                 :: time_axis_index
     integer                 :: id, type, natt, ndim
     type(atttype), pointer  :: Att(:) =>NULL()
     integer                 :: position ! indicate the location of the data ( CENTER, NORTH, EAST, CORNER )
  end type fieldtype

  type :: filetype
     private
     character(len=256) :: name
     integer            :: action, format, access, threading, fileset, record, ncid
     logical            :: opened, initialized, nohdrs
     integer            :: time_level
     real(DOUBLE_KIND)  :: time
     logical            :: valid
     logical            :: write_on_this_pe   ! indicate if will write out from this pe
     logical            :: io_domain_exist    ! indicate if io_domain exist or not.
     integer            :: id       !variable ID of time axis associated with file (only one time axis per file)
     integer            :: recdimid !dim ID of time axis associated with file (only one time axis per file)
     real(DOUBLE_KIND), pointer :: time_values(:) =>NULL() ! time axis values are stored here instead of axis%data 
                                                  ! since mpp_write assumes these values are not time values. 
                                                  ! Not used in mpp_write
     ! additional elements of filetype for mpp_read (ignored for mpp_write)
     integer :: ndim, nvar, natt  ! number of dimensions, non-dimension variables and global attributes
                                  ! redundant axis types stored here and in associated fieldtype
                                  ! some axes are not used by any fields, i.e. "edges"
     type(axistype), pointer  :: axis(:) =>NULL()
     type(fieldtype), pointer :: var(:) =>NULL()
     type(atttype), pointer   :: att(:) =>NULL()
     type(domain2d), pointer  :: domain =>NULL()
  end type filetype

!***********************************************************************
!
!     public interface from mpp_io_util.h
!
!***********************************************************************
  interface mpp_get_id
     module procedure mpp_get_axis_id
     module procedure mpp_get_field_id
  end interface

! <INTERFACE NAME="mpp_get_atts">
!   <OVERVIEW>
!     Get file global metdata.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Get file global metdata.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_get_atts( unit, global_atts)
!   </TEMPLATE>
!  <IN NAME="unit"></IN>
!  <IN NAME="global_atts"></IN>
! </INTERFACE>
  interface mpp_get_atts
     module procedure mpp_get_global_atts
     module procedure mpp_get_field_atts
     module procedure mpp_get_axis_atts
  end interface

  interface mpp_get_att_value
     module procedure mpp_get_field_att_text 
  end interface


!***********************************************************************
!
!      public interface from mpp_io_read.h
!
!***********************************************************************
! <INTERFACE NAME="mpp_read">
!   <OVERVIEW>
!     Read from an open file.
!   </OVERVIEW>
!   <DESCRIPTION>
!      <TT>mpp_read</TT> is used to read data to the file on an I/O unit
!      using the file parameters supplied by <LINK
!      SRC="#mpp_open"><TT>mpp_open</TT></LINK>. There are two
!      forms of <TT>mpp_read</TT>, one to read
!      distributed field data, and one to read non-distributed field
!      data. <I>Distributed</I> data refer to arrays whose two
!      fastest-varying indices are domain-decomposed. Distributed data must
!      be 2D or 3D (in space). Non-distributed data can be 0-3D.
!
!      The <TT>data</TT> argument for distributed data is expected by
!      <TT>mpp_read</TT> to contain data specified on the <I>data</I> domain,
!      and will read the data belonging to the <I>compute</I> domain,
!      fetching data as required by the parallel I/O <LINK
!      SRC="#modes">mode</LINK> specified in the <TT>mpp_open</TT> call. This
!      is consistent with our definition of <LINK
!      SRC="http:mpp_domains.html#domains">domains</LINK>, where all arrays are
!      expected to be dimensioned on the data domain, and all operations
!      performed on the compute domain.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_read( unit, field, data, time_index )
!   </TEMPLATE>
!   <TEMPLATE>
!     call mpp_read( unit, field, domain, data, time_index )
!   </TEMPLATE>
!  <IN NAME="unit"></IN>
!  <IN NAME="field"></IN>
!  <INOUT NAME="data"></INOUT>
!  <IN NAME="domain"></IN>
!  <IN NAME="time_index">
!     time_index is an optional argument. It is to be omitted if the
!     field was defined not to be a function of time. Results are
!     unpredictable if the argument is supplied for a time- independent
!     field, or omitted for a time-dependent field.
!  </IN>
!  <NOTE>
!     The type of read performed by <TT>mpp_read</TT> depends on
!     the file characteristics on the I/O unit specified at the <LINK
!     SRC="#mpp_open"><TT>mpp_open</TT></LINK> call. Specifically, the
!     format of the input data (e.g netCDF or IEEE) and the
!     <TT>threading</TT> flags, etc., can be changed there, and
!     require no changes to the <TT>mpp_read</TT>
!     calls. (<TT>fileset</TT> = MPP_MULTI is not supported by
!     <TT>mpp_read</TT>; IEEE is currently not supported).
!
!     Packed variables are unpacked using the <TT>scale</TT> and
!     <TT>add</TT> attributes.
!
!     <TT>mpp_read_meta</TT> must be called prior to calling <TT>mpp_read.</TT>
!  </NOTE>
! </INTERFACE>
  interface mpp_read
     module procedure mpp_read_2ddecomp_r2d
     module procedure mpp_read_2ddecomp_r3d
     module procedure mpp_read_r0D
     module procedure mpp_read_r1D
     module procedure mpp_read_r2D
     module procedure mpp_read_r3D
     module procedure mpp_read_text
     module procedure mpp_read_region_r2D
  end interface

!***********************************************************************
!
!    public interface from mpp_io_write.h
!
!***********************************************************************

! <INTERFACE NAME="mpp_write_meta">
!   <OVERVIEW>
!     Write metadata.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine is used to write the <LINK SRC="#metadata">metadata</LINK>
!     describing the contents of a file being written. Each file can contain
!     any number of fields, which are functions of 0-3 space axes and 0-1
!     time axes. (Only one time axis can be defined per file). The basic
!     metadata defined <LINK SRC="#metadata">above</LINK> for <TT>axistype</TT>
!     and <TT>fieldtype</TT> are written in the first two forms of the call
!     shown below. These calls will associate a unique variable ID with each
!     variable (axis or field). These can be used to attach any other real,
!     integer or character attribute to a variable. The last form is used to
!     define a <I>global</I> real, integer or character attribute that
!     applies to the dataset as a whole.
!   </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_write_meta( unit, axis, name, units, longname,
!      cartesian, sense, domain, data )
!  </TEMPLATE>
!  <NOTE>
!    The first form defines a time or space axis. Metadata corresponding to the type
!    above are written to the file on &lt;unit&gt;. A unique ID for subsequen
!    references to this axis is returned in axis%id. If the &lt;domain&gt;
!    element is present, this is recognized as a distributed data axis
!    and domain decomposition information is also written if required (the
!    domain decomposition info is required for multi-fileset multi-threaded
!    I/O). If the &lt;data&gt; element is allocated, it is considered to be a
!    space axis, otherwise it is a time axis with an unlimited dimension. Only
!    one time axis is allowed per file.
!  </NOTE>
!  <TEMPLATE>
!    call mpp_write_meta( unit, field, axes, name, units, longname,
!                              min, max, missing, fill, scale, add, pack )
!  </TEMPLATE>
!  <NOTE>
!    The second form defines a field. Metadata corresponding to the type
!    above are written to the file on &lt;unit&gt;. A unique ID for subsequen
!    references to this field is returned in field%id. At least one axis
!    must be associated, 0D variables are not considered. mpp_write_meta
!    must previously have been called on all axes associated with this
!    field.
!  </NOTE>
!  <TEMPLATE>
!    call mpp_write_meta( unit, id, name, rval=rval, pack=pack )
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_write_meta( unit, id, name, ival=ival )
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_write_meta( unit, id, name, cval=cval )
!  </TEMPLATE>
!  <NOTE>
!    The third form (3 - 5) defines metadata associated with a previously defined
!    axis or field, identified to mpp_write_meta by its unique ID &lt;id&gt;.
!    The attribute is named &lt;name&gt; and can take on a real, integer
!    or character value. &lt;rval&gt; and &lt;ival&gt; can be scalar or 1D arrays.
!    This need not be called for attributes already contained in
!    the type.
!  </NOTE>
!  <TEMPLATE>
!    call mpp_write_meta( unit, name, rval=rval, pack=pack )
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_write_meta( unit, name, ival=ival )
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_write_meta( unit, name, cval=cval )
!  </TEMPLATE>
!  <NOTE>
!    The last form (6 - 8) defines global metadata associated with the file as a
!    whole. The attribute is named &lt;name&gt; and can take on a real, integer
!    or character value. &lt;rval&gt; and &lt;ival&gt; can be scalar or 1D arrays.
!  </NOTE>
!  <IN NAME="unit"></IN>
!  <OUT NAME="axis"></OUT>
!  <IN NAME="name"></IN>
!  <IN NAME="units"></IN>
!  <IN NAME="longname"></IN>
!  <IN NAME="cartesian"></IN>
!  <IN NAME="sense"></IN>
!  <IN NAME="domain"></IN>
!  <IN NAME="data"></IN>
!  <OUT NAME="field"></OUT>
!  <IN NAME="min, max"></IN>
!  <IN NAME="missing"></IN>
!  <IN NAME="fill"></IN>
!  <IN NAME="scale"></IN>
!  <IN NAME="add"></IN>
!  <IN NAME="pack"></IN>
!  <IN NAME="id"></IN>
!  <IN NAME="cval"></IN>
!  <IN NAME="ival"></IN>
!  <IN NAME="rval"></IN>
! <NOTE>
!    Note that <TT>mpp_write_meta</TT> is expecting axis data on the
!    <I>global</I> domain even if it is a domain-decomposed axis.
!
!    You cannot interleave calls to <TT>mpp_write</TT> and
!    <TT>mpp_write_meta</TT>: the first call to
!    <TT>mpp_write</TT> implies that metadata specification is complete.
! </NOTE>
! </INTERFACE>
  interface mpp_write_meta
     module procedure mpp_write_meta_var
     module procedure mpp_write_meta_scalar_r
     module procedure mpp_write_meta_scalar_i
     module procedure mpp_write_meta_axis
     module procedure mpp_write_meta_field
     module procedure mpp_write_meta_global
     module procedure mpp_write_meta_global_scalar_r
     module procedure mpp_write_meta_global_scalar_i
  end interface
     
  interface mpp_copy_meta
     module procedure mpp_copy_meta_axis
     module procedure mpp_copy_meta_field
     module procedure mpp_copy_meta_global
  end interface

  interface mpp_modify_meta
!     module procedure mpp_modify_att_meta
     module procedure mpp_modify_field_meta
     module procedure mpp_modify_axis_meta
  end interface

! <INTERFACE NAME="mpp_write">
!   <OVERVIEW>
!     Write to an open file.
!   </OVERVIEW>
!   <DESCRIPTION>
!    <TT>mpp_write</TT> is used to write data to the file on an I/O unit
!    using the file parameters supplied by <LINK
!    SRC="#mpp_open"><TT>mpp_open</TT></LINK>. Axis and field definitions must
!    have previously been written to the file using <LINK
!    SRC="#mpp_write_meta"><TT>mpp_write_meta</TT></LINK>.  There are three
!    forms of <TT>mpp_write</TT>, one to write axis data, one to write
!    distributed field data, and one to write non-distributed field
!    data. <I>Distributed</I> data refer to arrays whose two
!    fastest-varying indices are domain-decomposed. Distributed data must
!    be 2D or 3D (in space). Non-distributed data can be 0-3D.
!
!    The <TT>data</TT> argument for distributed data is expected by
!    <TT>mpp_write</TT> to contain data specified on the <I>data</I> domain,
!    and will write the data belonging to the <I>compute</I> domain,
!    fetching or sending data as required by the parallel I/O <LINK
!    SRC="#modes">mode</LINK> specified in the <TT>mpp_open</TT> call. This
!    is consistent with our definition of <LINK
!    SRC="http:mpp_domains.html#domains">domains</LINK>, where all arrays are
!    expected to be dimensioned on the data domain, and all operations
!    performed on the compute domain.
!
!     The type of the <TT>data</TT> argument must be a <I>default
!     real</I>, which can be 4 or 8 byte.
!   </DESCRIPTION>
!  <TEMPLATE>
!    mpp_write( unit, axis )
!  </TEMPLATE>
!  <TEMPLATE>
!    mpp_write( unit, field, data, tstamp )
!  </TEMPLATE>
!  <TEMPLATE>
!    mpp_write( unit, field, domain, data, tstamp )
!  </TEMPLATE>
!  <IN NAME="tstamp">
!    <TT>tstamp</TT> is an optional argument. It is to
!    be omitted if the field was defined not to be a function of time.
!    Results are unpredictable if the argument is supplied for a time-
!    independent field, or omitted for a time-dependent field. Repeated
!    writes of a time-independent field are also not recommended. One
!    time level of one field is written per call. tstamp must be an 8-byte
!    real, even if the default real type is 4-byte.
!  </IN>
!  <NOTE>
!    The type of write performed by <TT>mpp_write</TT> depends on the file
!    characteristics on the I/O unit specified at the <LINK
!    SRC="#mpp_open"><TT>mpp_open</TT></LINK> call. Specifically, the format of
!    the output data (e.g netCDF or IEEE), the <TT>threading</TT> and
!    <TT>fileset</TT> flags, etc., can be changed there, and require no
!    changes to the <TT>mpp_write</TT> calls.
!
!    Packing is currently not implemented for non-netCDF files, and the
!    <TT>pack</TT> attribute is ignored. On netCDF files,
!    <TT>NF_DOUBLE</TT>s (8-byte IEEE floating point numbers) are
!    written for <TT>pack</TT>=1 and <TT>NF_FLOAT</TT>s for
!    <TT>pack</TT>=2. (<TT>pack</TT>=2 gives the customary
!    and default behaviour). We write <TT>NF_SHORT</TT>s (2-byte
!    integers) for <TT>pack=4</TT>, or <TT>NF_BYTE</TT>s
!    (1-byte integers) for <TT>pack=8</TT>. Integer scaling is done
!    using the <TT>scale</TT> and <TT>add</TT> attributes at
!    <TT>pack</TT>=4 or 8, satisfying the relation
!
!    <PRE>
!    data = packed_data*scale + add
!    </PRE>
!
!    <TT>NOTE: mpp_write</TT> does not check to see if the scaled
!    data in fact fits into the dynamic range implied by the specified
!    packing. It is incumbent on the user to supply correct scaling
!    attributes.
!
!    You cannot interleave calls to <TT>mpp_write</TT> and
!    <TT>mpp_write_meta</TT>: the first call to
!    <TT>mpp_write</TT> implies that metadata specification is
!    complete.
! </NOTE>
! </INTERFACE>
  interface mpp_write
     module procedure mpp_write_2ddecomp_r2d
     module procedure mpp_write_2ddecomp_r3d
     module procedure mpp_write_2ddecomp_r4d
     module procedure mpp_write_r0D
     module procedure mpp_write_r1D
     module procedure mpp_write_r2D
     module procedure mpp_write_r3D
     module procedure mpp_write_r4D
     module procedure mpp_write_axis
  end interface

!***********************************************************************
!
!            module variables
!
!***********************************************************************
  logical            :: module_is_initialized = .FALSE.
  logical            :: verbose =.FALSE.
  logical            :: debug = .FALSE.
  integer            :: maxunits, unit_begin, unit_end
  integer            :: mpp_io_stack_size=0, mpp_io_stack_hwm=0
  integer            :: varnum=0
  integer            :: pe, npes
  character(len=256) :: text
  integer            :: error
  integer            :: records_per_pe
  integer            :: mpp_read_clock=0, mpp_write_clock=0
  integer            :: mpp_open_clock=0, mpp_close_clock=0


!initial value of buffer between meta_data and data in .nc file
  integer            :: header_buffer_val = 16384  ! value used in NF__ENDDEF
  logical            :: global_field_on_root_pe = .true.
  logical            :: io_clocks_on = .false.
  namelist /mpp_io_nml/header_buffer_val, global_field_on_root_pe, io_clocks_on

  real(DOUBLE_KIND), allocatable :: mpp_io_stack(:)
  type(axistype),save            :: default_axis      !provided to users with default components
  type(fieldtype),save           :: default_field     !provided to users with default components
  type(atttype),save             :: default_att       !provided to users with default components
  type(filetype), allocatable    :: mpp_file(:)


  character(len=128) :: version= &
       '$Id: mpp_io.F90,v 18.0 2010/03/02 23:56:36 fms Exp $'
  character(len=128) :: tagname= &
       '$Name: testing $'

contains

#include <mpp_io_util.inc>
#include <mpp_io_misc.inc>
#include <mpp_io_connect.inc>
#include <mpp_io_read.inc>
#include <mpp_io_write.inc>

end module mpp_io_mod


