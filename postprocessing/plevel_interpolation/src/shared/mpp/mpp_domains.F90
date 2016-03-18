!-----------------------------------------------------------------------
!   Domain decomposition and domain update for message-passing codes
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

! <CONTACT EMAIL="V.Balaji@noaa.gov">
!   V. Balaji
! </CONTACT>
! <CONTACT EMAIL="Zhi.Liang@noaa.gov">
!   Zhi Liang
! </CONTACT>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp_domains.html"/>

! <OVERVIEW>
!   <TT>mpp_domains_mod</TT> is a set of simple calls for domain
!   decomposition and domain updates on rectilinear grids. It requires the
!   module <LINK SRC="mpp.html">mpp_mod</LINK>, upon which it is built.
! </OVERVIEW>

! <DESCRIPTION>
!   Scalable implementations of finite-difference codes are generally
!   based on decomposing the model domain into subdomains that are
!   distributed among processors. These domains will then be obliged to
!   exchange data at their boundaries if data dependencies are merely
!   nearest-neighbour, or may need to acquire information from the global
!   domain if there are extended data dependencies, as in the spectral
!   transform. The domain decomposition is a key operation in the
!   development of parallel codes.
!   
!   <TT>mpp_domains_mod</TT> provides a domain decomposition and domain
!   update API for <I>rectilinear</I> grids, built on top of the <LINK
!   SRC="mpp.html">mpp_mod</LINK> API for message passing. Features
!   of <TT>mpp_domains_mod</TT> include:
! 
!   Simple, minimal API, with free access to underlying API for more complicated stuff.
!
!   Design toward typical use in climate/weather CFD codes.
!  
!   <H4>Domains</H4>
! 
!   I have assumed that domain decomposition will mainly be in 2
!   horizontal dimensions, which will in general be the two
!   fastest-varying indices. There is a separate implementation of 1D
!   decomposition on the fastest-varying index, and 1D decomposition on
!   the second index, treated as a special case of 2D decomposition, is
!   also possible. We define <I>domain</I> as the grid associated with a <I>task</I>.
!   We define the <I>compute domain</I> as the set of gridpoints that are
!   computed by a task, and the <I>data domain</I> as the set of points
!   that are required by the task for the calculation. There can in
!   general be more than 1 task per PE, though often
!   the number of domains is the same as the processor count. We define
!   the <I>global domain</I> as the global computational domain of the
!   entire model (i.e, the same as the computational domain if run on a
!   single processor). 2D domains are defined using a derived type <TT>domain2D</TT>,
!   constructed as follows (see comments in code for more details):
!   
!   <PRE>
!     type, public :: domain_axis_spec
!        private
!        integer :: begin, end, size, max_size
!        logical :: is_global
!     end type domain_axis_spec
!     type, public :: domain1D
!        private
!        type(domain_axis_spec) :: compute, data, global, active
!        logical :: mustputb, mustgetb, mustputf, mustgetf, folded
!        type(domain1D), pointer, dimension(:) :: list
!        integer :: pe              !PE to which this domain is assigned
!        integer :: pos
!     end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain
!     type, public :: domain2D
!        private
!        type(domain1D) :: x
!        type(domain1D) :: y
!        type(domain2D), pointer, dimension(:) :: list
!        integer :: pe              !PE to which this domain is assigned
!        integer :: pos
!     end type domain2D
!     type(domain1D), public :: NULL_DOMAIN1D
!     type(domain2D), public :: NULL_DOMAIN2D
!   </PRE>

!   The <TT>domain2D</TT> type contains all the necessary information to
!   define the global, compute and data domains of each task, as well as the PE
!   associated with the task. The PEs from which remote data may be
!   acquired to update the data domain are also contained in a linked list
!   of neighbours.
! </DESCRIPTION>

module mpp_domains_mod
!a generalized domain decomposition package for use with mpp_mod
!Balaji (vb@gfdl.gov) 15 March 1999
  use mpp_parameter_mod,      only : MPP_DEBUG, MPP_VERBOSE, MPP_DOMAIN_TIME
  use mpp_parameter_mod,      only : GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, GLOBAL,CYCLIC 
  use mpp_parameter_mod,      only : AGRID, BGRID_SW, BGRID_NE, CGRID_NE, CGRID_SW, DGRID_NE, DGRID_SW
  use mpp_parameter_mod,      only : FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  use mpp_parameter_mod,      only : WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
  use mpp_parameter_mod,      only : NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, MPP_DOMAIN_TIME
  use mpp_parameter_mod,      only : CENTER, CORNER, SCALAR_PAIR, SCALAR_BIT
  use mpp_parameter_mod,      only : NORTH, NORTH_EAST, EAST, SOUTH_EAST
  use mpp_parameter_mod,      only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST
  use mpp_parameter_mod,      only : MAX_DOMAIN_FIELDS, NULL_PE, DOMAIN_ID_BASE
  use mpp_parameter_mod,      only : ZERO, NINETY, MINUS_NINETY, ONE_HUNDRED_EIGHTY, MAX_TILES
  use mpp_parameter_mod,      only : EVENT_SEND, EVENT_RECV, ROOT_GLOBAL
  use mpp_data_mod,           only : mpp_domains_stack, ptr_domains_stack
  use mpp_mod,                only : mpp_pe, mpp_root_pe, mpp_npes, mpp_error, FATAL, WARNING, NOTE
  use mpp_mod,                only : stdout, stderr, stdlog, mpp_send, mpp_recv, mpp_transmit, mpp_sync_self
  use mpp_mod,                only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_mod,                only : mpp_max, mpp_min, mpp_sum, mpp_get_current_pelist, mpp_broadcast
  use mpp_mod,                only : mpp_sync, mpp_init, mpp_malloc, lowercase
  use mpp_memutils_mod,       only : mpp_memuse_begin, mpp_memuse_end
  use mpp_pset_mod, only: mpp_pset_init
  implicit none
  private

#include <fms_platform.h>

  !--- public paramters imported from mpp_domains_parameter_mod
  public :: GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, BGRID_NE, BGRID_SW, CGRID_NE, CGRID_SW
  public :: DGRID_NE, DGRID_SW, FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  public :: WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
  public :: NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, MPP_DOMAIN_TIME
  public :: CENTER, CORNER, SCALAR_PAIR
  public :: NORTH, NORTH_EAST, EAST, SOUTH_EAST
  public :: SOUTH, SOUTH_WEST, WEST, NORTH_WEST
  public :: ZERO, NINETY, MINUS_NINETY, ONE_HUNDRED_EIGHTY 

  !--- public data imported from mpp_data_mod
  public :: NULL_DOMAIN1D, NULL_DOMAIN2D

  public :: domain_axis_spec, domain1D, domain2D, DomainCommunicator2D

  !--- public interface from mpp_domains_util.h
  public :: mpp_domains_set_stack_size, mpp_get_compute_domain, mpp_get_compute_domains
  public :: mpp_get_data_domain, mpp_get_global_domain, mpp_get_domain_components
  public :: mpp_get_layout, mpp_get_pelist, operator(.EQ.), operator(.NE.) 
  public :: mpp_domain_is_symmetry
  public :: mpp_get_neighbor_pe, mpp_nullify_domain_list
  public :: mpp_set_compute_domain, mpp_set_data_domain, mpp_set_global_domain
  public :: mpp_get_memory_domain, mpp_get_domain_shift, mpp_domain_is_tile_root_pe
  public :: mpp_get_tile_id, mpp_get_domain_extents, mpp_get_current_ntile, mpp_get_ntile_count
  public :: mpp_get_refine_overlap_number, mpp_get_mosaic_refine_overlap
  public :: mpp_get_tile_list
  public :: mpp_get_tile_npes
  public :: mpp_get_num_overlap, mpp_get_overlap
  public :: mpp_get_io_domain, mpp_get_domain_pe, mpp_get_domain_tile_root_pe
  public :: mpp_get_domain_name, mpp_get_io_domain_layout
  public :: mpp_copy_domain, mpp_set_domain_symmetry
  public :: mpp_get_update_pelist, mpp_get_update_size

  !--- public interface from mpp_domains_reduce.h
  public :: mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum
!  public :: mpp_global_sum_tl, mpp_global_sum_ad
  !--- public interface from mpp_domains_misc.h
  public :: mpp_broadcast_domain, mpp_domains_init, mpp_domains_exit, mpp_redistribute
  public :: mpp_update_domains, mpp_check_field
!  public :: mpp_update_domains_ad   ! bnc
  public :: mpp_get_boundary
  !--- public interface from mpp_domains_define.h
  public :: mpp_define_layout, mpp_define_domains, mpp_modify_domain, mpp_define_mosaic
  public :: mpp_define_mosaic_pelist, mpp_define_null_domain, mpp_mosaic_defined
  public :: mpp_define_io_domain, mpp_deallocate_domain
  public :: mpp_compute_extent

  integer, parameter :: NAME_LENGTH = 64
  integer, parameter :: MAXLIST = 8
 
  !--- data types used mpp_domains_mod.
  type domain_axis_spec        !type used to specify index limits along an axis of a domain
     private
     integer :: begin, end, size, max_size      !start, end of domain axis, size, max size in set
     logical :: is_global       !TRUE if domain axis extent covers global domain
  end type domain_axis_spec

  type domain1D
     private
     type(domain_axis_spec) :: compute, data, global, memory
     logical :: cyclic
     type(domain1D), pointer :: list(:) =>NULL()
     integer :: pe               !PE to which this domain is assigned
     integer :: pos              !position of this PE within link list, i.e domain%list(pos)%pe = pe
     integer :: goffset, loffset !needed for global sum
  end type domain1D

  type domain1D_spec
     private
     type(domain_axis_spec) :: compute
     integer                :: pos
  end type domain1D_spec
       
  type domain2D_spec
     private
     type(domain1D_spec), pointer :: x(:)       => NULL() ! x-direction domain decomposition
     type(domain1D_spec), pointer :: y(:)       => NULL() ! x-direction domain decomposition
     integer,        pointer :: tile_id(:) => NULL() ! tile id of each tile
     integer                 :: pe                   ! PE to which this domain is assigned
     integer                 :: pos                  ! position of this PE within link list
     integer                 :: tile_root_pe         ! root pe of tile.
  end type domain2D_spec

  type overlap_type
     private
     integer                  :: count = 0                 ! number of ovrelapping
     integer                  :: pe
     integer,         pointer :: tileMe(:)       => NULL() ! my tile id for this overlap
     integer,         pointer :: tileNbr(:)      => NULL() ! neighbor tile id for this overlap
     integer,         pointer :: is(:)           => NULL() ! starting i-index 
     integer,         pointer :: ie(:)           => NULL() ! ending   i-index 
     integer,         pointer :: js(:)           => NULL() ! starting j-index 
     integer,         pointer :: je(:)           => NULL() ! ending   j-index 
     integer,         pointer :: isMe(:)         => NULL() ! starting i-index of my tile on current pe
     integer,         pointer :: ieMe(:)         => NULL() ! ending   i-index of my tile on current pe
     integer,         pointer :: jsMe(:)         => NULL() ! starting j-index of my tile on current pe
     integer,         pointer :: jeMe(:)         => NULL() ! ending   j-index of my tile on current pe
     integer,         pointer :: dir(:)          => NULL() ! direction ( value 1,2,3,4 = E,S,W,N)
     integer,         pointer :: rotation(:)     => NULL() ! rotation angle.
     logical,         pointer :: is_refined(:)   => NULL() ! indicate if the overlap is refined or not.
     integer,         pointer :: index(:)        => NULL() ! for refinement
     logical,         pointer :: from_contact(:) => NULL() ! indicate if the overlap is computed from define_contact_overlap
  end type overlap_type

  type overlapSpec
     private
     integer                     :: whalo, ehalo, shalo, nhalo ! halo size
     integer                     :: xbegin, xend, ybegin, yend
     integer                     :: nsend, nrecv
     type(overlap_type), pointer :: send(:) => NULL()
     type(overlap_type), pointer :: recv(:) => NULL()
     type(refineSpec),   pointer :: rSpec(:)=> NULL()
     type(overlapSpec),  pointer :: next
  end type overlapSpec

  type tile_type
     integer :: xbegin, xend, ybegin, yend
  end type tile_type

  type refineSpec
     private
     integer          :: count                 ! number of ovrelapping
     integer          :: total                 ! total number of points to be saved in buffer.
     integer, pointer :: isMe(:)     => NULL() ! starting i-index on current pe and tile.
     integer, pointer :: ieMe(:)     => NULL() ! ending i-index on current pe and tile.
     integer, pointer :: jsMe(:)     => NULL() ! starting j-index on current pe and tile.
     integer, pointer :: jeMe(:)     => NULL() ! ending j-index on current pe and tile.
     integer, pointer :: isNbr(:)    => NULL() ! starting i-index on neighbor pe or tile
     integer, pointer :: ieNbr(:)    => NULL() ! ending i-index on neighbor pe or tile
     integer, pointer :: jsNbr(:)    => NULL() ! starting j-index on neighbor pe or tile
     integer, pointer :: jeNbr(:)    => NULL() ! ending j-index on neighbor pe or tile
     integer, pointer :: start(:)    => NULL() ! starting index in the buffer
     integer, pointer :: end(:)      => NULL() ! ending index in the buffer
     integer, pointer :: dir(:)      => NULL() ! direction 
     integer, pointer :: rotation(:) => NULL() ! rotation angle.
  end type refineSpec

!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain, not once per axis

  type domain2D
     private
     character(len=NAME_LENGTH)  :: name='unnamed'          ! name of the domain, default is "unspecified"
     integer(LONG_KIND)          :: id 
     integer                     :: pe                      ! PE to which this domain is assigned
     integer                     :: fold          
     integer                     :: pos                     ! position of this PE within link list
     logical                     :: symmetry                ! indicate the domain is symmetric or non-symmetric.
     integer                     :: whalo, ehalo            ! halo size in x-direction
     integer                     :: shalo, nhalo            ! halo size in y-direction
     integer                     :: ntiles                  ! number of tiles within mosaic
     integer                     :: max_ntile_pe            ! maximum value in the pelist of number of tiles on each pe.
     integer                     :: ncontacts               ! number of contact region within mosaic.
     logical                     :: rotated_ninety          ! indicate if any contact rotate NINETY or MINUS_NINETY
     logical                     :: initialized             ! indicate if the overlapping is computed or not.
     integer                     :: tile_root_pe            ! root pe of current tile.
     integer                     :: io_layout(2)            ! io_layout, will be set through mpp_define_io_domain
                                                            ! default = domain layout
     integer,            pointer :: pearray(:,:)  => NULL() ! pe of each layout position 
     integer,            pointer :: tile_id(:)    => NULL() ! tile id of each tile
     type(domain1D),     pointer :: x(:)          => NULL() ! x-direction domain decomposition
     type(domain1D),     pointer :: y(:)          => NULL() ! y-direction domain decomposition
     type(domain2D_spec),pointer :: list(:)       => NULL() ! domain decomposition on pe list
     type(tile_type),    pointer :: tileList(:)   => NULL() ! store tile information
     type(overlapSpec),  pointer :: check_C       => NULL() ! send and recv information for boundary consistency check of C-cell
     type(overlapSpec),  pointer :: check_E       => NULL() ! send and recv information for boundary consistency check of E-cell
     type(overlapSpec),  pointer :: check_N       => NULL() ! send and recv information for boundary consistency check of N-cell
     type(overlapSpec),  pointer :: bound_C       => NULL() ! send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: bound_E       => NULL() ! send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: bound_N       => NULL() ! send information for getting boundary value for symmetry domain.
     type(overlapSpec),  pointer :: update_T      => NULL() ! send and recv information for halo update of T-cell.
     type(overlapSpec),  pointer :: update_E      => NULL() ! send and recv information for halo update of E-cell.
     type(overlapSpec),  pointer :: update_C      => NULL() ! send and recv information for halo update of C-cell.
     type(overlapSpec),  pointer :: update_N      => NULL() ! send and recv information for halo update of N-cell.
     type(domain2d),     pointer :: io_domain     => NULL() ! domain for IO, will be set through calling mpp_set_io_domain ( this will be changed).
  end type domain2D     

  !--- the following type is used to reprsent the contact between tiles.
  !--- this type will only be used in mpp_domains_define.inc
  type contact_type
     private
     integer          :: ncontact                               ! number of neighbor tile.
     integer, pointer :: tile(:) =>NULL()                      ! neighbor tile 
     integer, pointer :: align1(:)=>NULL(), align2(:)=>NULL()   ! alignment of me and neighbor
     real,    pointer :: refine1(:)=>NULL(), refine2(:)=>NULL() !
     integer, pointer :: is1(:)=>NULL(), ie1(:)=>NULL()         ! i-index of current tile repsenting contact
     integer, pointer :: js1(:)=>NULL(), je1(:)=>NULL()         ! j-index of current tile repsenting contact
     integer, pointer :: is2(:)=>NULL(), ie2(:)=>NULL()         ! i-index of neighbor tile repsenting contact
     integer, pointer :: js2(:)=>NULL(), je2(:)=>NULL()         ! j-index of neighbor tile repsenting contact
  end type contact_type


  type DomainCommunicator2D
     private
     logical            :: initialized=.false.
     integer(LONG_KIND) :: id=-9999
     integer(LONG_KIND) :: l_addr  =-9999
     integer(LONG_KIND) :: l_addrx =-9999
     integer(LONG_KIND) :: l_addry =-9999
     type(domain2D), pointer :: domain     =>NULL()
     type(domain2D), pointer :: domain_in  =>NULL()
     type(domain2D), pointer :: domain_out =>NULL()
     type(overlapSpec), pointer :: send(:,:,:,:) => NULL()
     type(overlapSpec), pointer :: recv(:,:,:,:) => NULL()
     integer, dimension(:,:),       _ALLOCATABLE :: sendis _NULL
     integer, dimension(:,:),       _ALLOCATABLE :: sendie _NULL
     integer, dimension(:,:),       _ALLOCATABLE :: sendjs _NULL
     integer, dimension(:,:),       _ALLOCATABLE :: sendje _NULL
     integer, dimension(:,:),       _ALLOCATABLE :: recvis _NULL
     integer, dimension(:,:),       _ALLOCATABLE :: recvie _NULL
     integer, dimension(:,:),       _ALLOCATABLE :: recvjs _NULL
     integer, dimension(:,:),       _ALLOCATABLE :: recvje _NULL
     logical, dimension(:),         _ALLOCATABLE :: S_do_buf _NULL
     logical, dimension(:),         _ALLOCATABLE :: R_do_buf _NULL
     integer, dimension(:),         _ALLOCATABLE :: cto_pe  _NULL
     integer, dimension(:),         _ALLOCATABLE :: cfrom_pe  _NULL
     integer, dimension(:),         _ALLOCATABLE :: S_msize _NULL
     integer, dimension(:),         _ALLOCATABLE :: R_msize _NULL
     integer :: Slist_size=0, Rlist_size=0
     integer :: isize=0, jsize=0, ke=0
     integer :: isize_in=0, jsize_in=0
     integer :: isize_out=0, jsize_out=0
     integer :: isize_max=0, jsize_max=0
     integer :: gf_ioff=0, gf_joff=0
  ! Remote data
     integer, dimension(:)  , _ALLOCATABLE :: isizeR _NULL
     integer, dimension(:)  , _ALLOCATABLE :: jsizeR _NULL
     integer, dimension(:,:), _ALLOCATABLE :: sendisR _NULL
     integer, dimension(:,:), _ALLOCATABLE :: sendjsR _NULL
     integer(LONG_KIND), dimension(:), _ALLOCATABLE :: rem_addr  _NULL
     integer(LONG_KIND), dimension(:), _ALLOCATABLE :: rem_addrx _NULL
     integer(LONG_KIND), dimension(:), _ALLOCATABLE :: rem_addry _NULL
     integer(LONG_KIND), dimension(:,:), _ALLOCATABLE :: rem_addrl  _NULL
     integer(LONG_KIND), dimension(:,:), _ALLOCATABLE :: rem_addrlx  _NULL
     integer(LONG_KIND), dimension(:,:), _ALLOCATABLE :: rem_addrly  _NULL
     integer                             :: position        ! data location. T, E, C, or N.
  end type DomainCommunicator2D

!#######################################################################

!***********************************************************************
!
!     module variables 
!
!***********************************************************************
  integer             :: pe
  logical             :: module_is_initialized = .false.
  logical             :: debug                 = .FALSE.
  logical             :: verbose=.FALSE.
  logical             :: mosaic_defined = .false.
  integer             :: mpp_domains_stack_size=0
  integer             :: mpp_domains_stack_hwm=0
  type(domain1D),save :: NULL_DOMAIN1D
  type(domain2D),save :: NULL_DOMAIN2D

  !-------- The following variables are used in mpp_domains_comm.h
  
  integer, parameter :: MAX_ADDRS=512
  integer(LONG_KIND),dimension(MAX_ADDRS),save :: addrs_sorted=-9999  ! list of sorted local addrs
  integer,           dimension(-1:MAX_ADDRS),save :: addrs_idx=-9999  ! idx of addr assoicated w/ d_comm
  integer,           dimension(MAX_ADDRS),save :: a_salvage=-9999     ! freed idx list of addr
  integer,                                save :: a_sort_len=0        ! len sorted memory list
  integer,                                save :: n_addrs=0           ! num memory addresses used

  integer(LONG_KIND), parameter :: ADDR2_BASE=Z'0000000000010000'
  integer, parameter :: MAX_ADDRS2=128
  integer(LONG_KIND),dimension(MAX_ADDRS2),save :: addrs2_sorted=-9999  ! list of sorted local addrs
  integer,           dimension(-1:MAX_ADDRS2),save :: addrs2_idx=-9999  ! idx of addr2 assoicated w/ d_comm
  integer,           dimension(MAX_ADDRS2),save :: a2_salvage=-9999     ! freed indices of addr2
  integer,                                 save :: a2_sort_len=0        ! len sorted memory list
  integer,                                 save :: n_addrs2=0           ! num memory addresses used

  integer, parameter :: MAX_DOM_IDS=128
  integer(LONG_KIND),dimension(MAX_DOM_IDS),save :: ids_sorted=-9999 ! list of sorted domain identifiers
  integer,           dimension(-1:MAX_DOM_IDS),save :: ids_idx=-9999 ! idx of d_comm associated w/ sorted addr
  integer,                                  save :: i_sort_len=0     ! len sorted domain ids list
  integer,                                  save :: n_ids=0          ! num domain ids used (=i_sort_len; dom ids never removed)

  integer, parameter :: MAX_FIELDS=1024
  integer(LONG_KIND),        dimension(MAX_FIELDS),save           :: dcKey_sorted=-9999  ! list of sorted local addrs
  !     Not sure why static d_comm fails during deallocation of derived type members; allocatable works
  !     type(DomainCommunicator2D),dimension(MAX_FIELDS),save,target    :: d_comm              ! domain communicators
  type(DomainCommunicator2D),dimension(:),allocatable,save,target :: d_comm              ! domain communicators
  integer,                   dimension(-1:MAX_FIELDS),save           :: d_comm_idx=-9999 ! idx of d_comm associated w/ sorted addr
  integer,                   dimension(MAX_FIELDS),save           :: dc_salvage=-9999    ! freed indices of d_comm
  integer,                                         save           :: dc_sort_len=0       ! len sorted comm keys (=num active communicators)
  integer,                                         save           :: n_comm=0            ! num communicators used

  !     integer(LONG_KIND), parameter :: GT_BASE=2**8
  integer(LONG_KIND), parameter :: GT_BASE=Z'0000000000000100'  ! Workaround for 64bit int init problem

  !     integer(LONG_KIND), parameter :: KE_BASE=2**48
  integer(LONG_KIND), parameter :: KE_BASE=Z'0001000000000000'  ! Workaround for 64bit int init problem

  integer, parameter :: MAXOVERLAP = 100 

  integer(LONG_KIND) :: domain_cnt=0

  !--- the following variables are used in mpp_domains_misc.h
  logical :: domain_clocks_on=.FALSE.
  integer :: send_clock=0, recv_clock=0, unpk_clock=0
  integer :: wait_clock=0, pack_clock=0, pack_loop_clock=0

  !--- namelist interface
! <NAMELIST NAME="mpp_domains_nml">
!   <DATA NAME="debug_update_domain" TYPE="character(len=32)"  DEFAULT="none">
!     when debug_update_domain = none, no debug will be done. When debug_update_domain is set to fatal, 
!     the run will be exited with fatal error message. When debug_update_domain is set to 
!     warning, the run will output warning message. when debug update_domain is set to 
!     note, the run will output some note message. Will check the consistency on the boundary between
!     processor/tile when updating doamin for symmetric domain and check the consistency on the north
!     folded edge. 
!   </DATA>
! </NAMELIST>
  character(len=32) :: debug_update_domain = "none"
  logical           :: debug_message_passing = .false.
  namelist /mpp_domains_nml/ debug_update_domain, domain_clocks_on, debug_message_passing

  !***********************************************************************

  integer, parameter :: NO_CHECK = -1
  integer            :: debug_update_level = NO_CHECK
!***********************************************************************
!
!         public interface from mpp_domains_define.h
!
!***********************************************************************

  ! <INTERFACE NAME="mpp_define_layout">
  !  <OVERVIEW>
  !    Retrieve layout associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Given a global 2D domain and the number of divisions in the
  !    decomposition (<TT>ndivs</TT>: usually the PE count unless some
  !    domains are masked) this calls returns a 2D domain layout.
  !    
  !    By default, <TT>mpp_define_layout</TT> will attempt to divide the
  !    2D index space into domains that maintain the aspect ratio of the
  !    global domain. If this cannot be done, the algorithm favours domains
  !    that are longer in <TT>x</TT> than <TT>y</TT>, a preference that could
  !    improve vector performance.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_define_layout( global_indices, ndivs, layout )
  !  </TEMPLATE>
  !  <IN NAME="global_indices"></IN>
  !  <IN NAME="ndivs"></IN>
  !  <OUT NAME="layout"></OUT>
  ! </INTERFACE>

  interface mpp_define_layout
     module procedure mpp_define_layout2D
  end interface


  ! <INTERFACE NAME="mpp_define_domains">

  !   <OVERVIEW>
  !     Set up a domain decomposition.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     There are two forms for the <TT>mpp_define_domains</TT> call. The 2D
  !     version is generally to be used but is built by repeated calls to the
  !     1D version, also provided.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call mpp_define_domains( global_indices, ndivs, domain, &
  !                                   pelist, flags, halo, extent, maskmap )
  !   </TEMPLATE>
  !  <TEMPLATE>
  !    call mpp_define_domains( global_indices, layout, domain, pelist, &
  !                                   xflags, yflags, xhalo, yhalo,           &
  !                                   xextent, yextent, maskmap, name )
  !  </TEMPLATE>
  !   <IN NAME="global_indices" >
  !     Defines the global domain.
  !   </IN>
  !   <IN NAME="ndivs">
  !     Is the number of domain divisions required.
  !   </IN>
  !   <INOUT NAME="domain">
  !     Holds the resulting domain decomposition.
  !   </INOUT>
  !   <IN NAME="pelist">
  !     List of PEs to which the domains are to be assigned.
  !   </IN>
  !   <IN NAME="flags">
  !      An optional flag to pass additional information
  !      about the desired domain topology. Useful flags in a 1D decomposition
  !      include <TT>GLOBAL_DATA_DOMAIN</TT> and
  !      <TT>CYCLIC_GLOBAL_DOMAIN</TT>. Flags are integers: multiple flags may
  !      be added together. The flag values are public parameters available by
  !      use association.
  !   </IN>
  !   <IN NAME="halo">
  !     Width of the halo.
  !   </IN>
  !   <IN NAME="extent">
  !      Normally <TT>mpp_define_domains</TT> attempts
  !      an even division of the global domain across <TT>ndivs</TT>
  !      domains. The <TT>extent</TT> array can be used by the user to pass a
  !      custom domain division. The <TT>extent</TT> array has <TT>ndivs</TT>
  !      elements and holds the compute domain widths, which should add up to
  !      cover the global domain exactly.
  !   </IN>
  !   <IN NAME="maskmap">
  !     Some divisions may be masked
  !     (<TT>maskmap=.FALSE.</TT>) to exclude them from the computation (e.g
  !     for ocean model domains that are all land). The <TT>maskmap</TT> array
  !     is dimensioned <TT>ndivs</TT> and contains <TT>.TRUE.</TT> values for
  !     any domain that must be <I>included</I> in the computation (default
  !     all). The <TT>pelist</TT> array length should match the number of
  !     domains included in the computation.
  !    </IN>   

  !  <IN NAME="layout"></IN>
  !  <IN NAME="xflags, yflags"></IN>
  !  <IN NAME="xhalo, yhalo"></IN>
  !  <IN NAME="xextent, yextent"></IN>
  !  <IN NAME="name" ></IN>

  !  <NOTE>    
  !    For example:
  !    
  !    <PRE>
  !    call mpp_define_domains( (/1,100/), 10, domain, &
  !         flags=GLOBAL_DATA_DOMAIN+CYCLIC_GLOBAL_DOMAIN, halo=2 )
  !    </PRE>
  !    
  !    defines 10 compute domains spanning the range [1,100] of the global
  !    domain. The compute domains are non-overlapping blocks of 10. All the data
  !    domains are global, and with a halo of 2 span the range [-1:102]. And
  !    since the global domain has been declared to be cyclic,
  !    <TT>domain(9)%next => domain(0)</TT> and <TT>domain(0)%prev =>
  !    domain(9)</TT>. A field is allocated on the data domain, and computations proceed on
  !    the compute domain. A call to <LINK
  !    SRC="#mpp_update_domains"><TT>mpp_update_domains</TT></LINK> would fill in
  !    the values in the halo region:

  !    <PRE>
  !    call mpp_get_data_domain( domain, isd, ied ) !returns -1 and 102
  !    call mpp_get_compute_domain( domain, is, ie ) !returns (1,10) on PE 0 ...
  !    allocate( a(isd:ied) )
  !    do i = is,ie
  !       a(i) = &lt;perform computations&gt;
  !    end do
  !    call mpp_update_domains( a, domain )
  !    </PRE>

  !    The call to <TT>mpp_update_domains</TT> fills in the regions outside
  !    the compute domain. Since the global domain is cyclic, the values at
  !    <TT>i=(-1,0)</TT> are the same as at <TT>i=(99,100)</TT>; and
  !    <TT>i=(101,102)</TT> are the same as <TT>i=(1,2)</TT>.
  !    
  !    The 2D version is just an extension of this syntax to two
  !    dimensions.
  !
  !    The 2D version of the above should generally be used in
  !    codes, including 1D-decomposed ones, if there is a possibility of
  !    future evolution toward 2D decomposition. The arguments are similar to
  !    the 1D case, except that now we have optional arguments
  !    <TT>flags</TT>, <TT>halo</TT>, <TT>extent</TT> and <TT>maskmap</TT>
  !    along two axes.
  !    
  !    <TT>flags</TT> can now take an additional possible value to fold
  !    one or more edges. This is done by using flags
  !    <TT>FOLD_WEST_EDGE</TT>, <TT>FOLD_EAST_EDGE</TT>,
  !    <TT>FOLD_SOUTH_EDGE</TT> or <TT>FOLD_NORTH_EDGE</TT>. When a fold
  !    exists (e.g cylindrical domain), vector fields reverse sign upon
  !    crossing the fold. This parity reversal is performed only in the
  !    vector version of <LINK
  !    SRC="#mpp_update_domains"><TT>mpp_update_domains</TT></LINK>. In
  !    addition, shift operations may need to be applied to vector fields on
  !    staggered grids, also described in the vector interface to
  !    <TT>mpp_update_domains</TT>.
  !    
  !    <TT>name</TT> is the name associated with the decomposition,
  !    e.g <TT>'Ocean model'</TT>. If this argument is present,
  !    <TT>mpp_define_domains</TT> will print the domain decomposition
  !    generated to <TT>stdlog</TT>.
  !    
  !    Examples:
  !    
  !    <PRE>
  !    call mpp_define_domains( (/1,100,1,100/), (/2,2/), domain, xhalo=1 )
  !    </PRE>
  !    
  !    will create the following domain layout:
  !    <PRE>
  !                   |---------|-----------|-----------|-------------|
  !                   |domain(1)|domain(2)  |domain(3)  |domain(4)    |
  !    |--------------|---------|-----------|-----------|-------------|
  !    |Compute domain|1,50,1,50|51,100,1,50|1,50,51,100|51,100,51,100|
  !    |--------------|---------|-----------|-----------|-------------|
  !    |Data domain   |0,51,1,50|50,101,1,50|0,51,51,100|50,101,51,100|
  !    |--------------|---------|-----------|-----------|-------------|
  !    </PRE>
  !    
  !    Again, we allocate arrays on the data domain, perform computations
  !    on the compute domain, and call <TT>mpp_update_domains</TT> to update
  !    the halo region.
  !    
  !    If we wished to perfom a 1D decomposition along <TT>Y</TT>
  !    on the same global domain, we could use:

  !    <PRE>
  !    call mpp_define_domains( (/1,100,1,100/), layout=(/4,1/), domain, xhalo=1 )
  !    </PRE>

  !    This will create the following domain layout:
  !    <PRE>
  !                   |----------|-----------|-----------|------------|
  !                   |domain(1) |domain(2)  |domain(3)  |domain(4)   |
  !    |--------------|----------|-----------|-----------|------------|
  !    |Compute domain|1,100,1,25|1,100,26,50|1,100,51,75|1,100,76,100|
  !    |--------------|----------|-----------|-----------|------------|
  !    |Data domain   |0,101,1,25|0,101,26,50|0,101,51,75|1,101,76,100|
  !    |--------------|----------|-----------|-----------|------------|
  !    </PRE>
  !   </NOTE>
  ! </INTERFACE>
  interface mpp_define_domains
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
  end interface

  interface mpp_define_null_domain
     module procedure mpp_define_null_domain1D
     module procedure mpp_define_null_domain2D
  end interface

  interface mpp_copy_domain
     module procedure mpp_copy_domain1D
     module procedure mpp_copy_domain2D
  end interface mpp_copy_domain

  interface mpp_deallocate_domain 
     module procedure mpp_deallocate_domain1D
     module procedure mpp_deallocate_domain2D
  end interface

! <INTERFACE NAME="mpp_modify_domain">
!   <OVERVIEW>
!     modifies the extents (compute, data and global) of domain
!   </OVERVIEW>
!   <IN NAME="domain_in">
!     The source domain.
!   </IN>
!   <IN NAME="halo">
!     Halo size of the returned 1D doamin. Default value is 0.
!   </IN>
!   <IN NAME="cbegin,cend">
!    Axis specifications associated with the compute domain of the returned 1D domain.
!   </IN>
!   <IN NAME="gbegin,gend">
!    Axis specifications associated with the global domain of the returned 1D domain.
!   </IN>
!   <IN NAME="isc,iec">
!    Zonal axis specifications associated with the compute domain of the returned 2D domain.
!   </IN>
!   <IN NAME="jsc,jec">
!    Meridinal axis specifications associated with the compute domain of the returned 2D domain.
!   </IN>
!   <IN NAME="isg,ieg">
!    Zonal axis specifications associated with the global domain of the returned 2D domain.
!   </IN>
!   <IN NAME="jsg,jeg">
!    Meridinal axis specifications associated with the global domain of the returned 2D domain.
!   </IN>
!   <IN NAME="xhalo,yhalo">
!     Halo size of the returned 2D doamin. Default value is 0.
!   </IN>
!   <INOUT NAME="domain_out">
!     The returned domain.
!   </INOUT>

! </INTERFACE>

  interface mpp_modify_domain
     module procedure mpp_modify_domain1D
     module procedure mpp_modify_domain2D
  end interface


!***********************************************************************
!
!        public interface from mpp_domains_misc.h
!
!***********************************************************************

! <INTERFACE NAME="mpp_update_domains">
!  <OVERVIEW>
!     Halo updates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_update_domains</TT> is used to perform a halo update of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!    of 4-byte or 8-byte kind; of rank up to 5. The vector version (with
!    two input data fields) is only present for <TT>real</TT> types.
!    
!    For 2D domain updates, if there are halos present along both
!    <TT>x</TT> and <TT>y</TT>, we can choose to update one only, by
!    specifying <TT>flags=XUPDATE</TT> or <TT>flags=YUPDATE</TT>. In
!    addition, one-sided updates can be performed by setting <TT>flags</TT>
!    to any combination of <TT>WUPDATE</TT>, <TT>EUPDATE</TT>,
!    <TT>SUPDATE</TT> and <TT>NUPDATE</TT>, to update the west, east, north
!    and south halos respectively. Any combination of halos may be used by
!    adding the requisite flags, e.g: <TT>flags=XUPDATE+SUPDATE</TT> or
!    <TT>flags=EUPDATE+WUPDATE+SUPDATE</TT> will update the east, west and
!    south halos.
!    
!    If a call to <TT>mpp_update_domains</TT> involves at least one E-W
!    halo and one N-S halo, the corners involved will also be updated, i.e,
!    in the example above, the SE and SW corners will be updated.
!    
!    If <TT>flags</TT> is not supplied, that is
!    equivalent to <TT>flags=XUPDATE+YUPDATE</TT>.
!    
!    The vector version is passed the <TT>x</TT> and <TT>y</TT>
!    components of a vector field in tandem, and both are updated upon
!    return. They are passed together to treat parity issues on various
!    grids. For example, on a cubic sphere projection, the <TT>x</TT> and
!    <TT>y</TT> components may be interchanged when passing from an
!    equatorial cube face to a polar face. For grids with folds, vector
!    components change sign on crossing the fold.  Paired scalar quantities
!    can also be passed with the vector version if flags=SCALAR_PAIR, in which
!    case components are appropriately interchanged, but signs are not.
!    
!    Special treatment at boundaries such as folds is also required for
!    staggered grids. The following types of staggered grids are
!    recognized:
!    
!    1) <TT>AGRID</TT>: values are at grid centers.<BR/>
!    2) <TT>BGRID_NE</TT>: vector fields are at the NE vertex of a grid
!    cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT> are
!    actually at (i+&#189;,j+&#189;) with respect to the grid centers.<BR/>
!    3) <TT>BGRID_SW</TT>: vector fields are at the SW vertex of a grid
!    cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT> are
!    actually at (i-&#189;,j-&#189;) with respect to the grid centers.<BR/>
!    4) <TT>CGRID_NE</TT>: vector fields are at the N and E faces of a
!    grid cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT>
!    are actually at (i+&#189;,j) and (i,j+&#189;) with respect to the
!    grid centers.<BR/>
!    5) <TT>CGRID_SW</TT>: vector fields are at the S and W faces of a
!    grid cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT>
!    are actually at (i-&#189;,j) and (i,j-&#189;) with respect to the
!    grid centers.
!
!    The gridtypes listed above are all available by use association as
!    integer parameters. The scalar version of <TT>mpp_update_domains</TT>
!    assumes that the values of a scalar field are always at <TT>AGRID</TT>
!    locations, and no special boundary treatment is required. If vector
!    fields are at staggered locations, the optional argument
!    <TT>gridtype</TT> must be appropriately set for correct treatment at
!    boundaries.
!    
!    It is safe to apply vector field updates to the appropriate arrays
!    irrespective of the domain topology: if the topology requires no
!    special treatment of vector fields, specifying <TT>gridtype</TT> will
!    do no harm.
!
!    <TT>mpp_update_domains</TT> internally buffers the date being sent
!    and received into single messages for efficiency. A turnable internal
!    buffer area in memory is provided for this purpose by
!    <TT>mpp_domains_mod</TT>. The size of this buffer area can be set by
!    the user by calling <LINK SRC="mpp_domains.html#mpp_domains_set_stack_size">
!    <TT>mpp_domains_set_stack_size</TT></LINK>.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_update_domains( field, domain, flags )
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_update_domains( fieldx, fieldy, domain, flags, gridtype )
!  </TEMPLATE>
! </INTERFACE>
  interface mpp_update_domains
     module procedure mpp_update_domain2D_r8_2d
     module procedure mpp_update_domain2D_r8_3d
     module procedure mpp_update_domain2D_r8_4d
     module procedure mpp_update_domain2D_r8_5d
     module procedure mpp_update_domain2D_r8_2dv
     module procedure mpp_update_domain2D_r8_3dv
     module procedure mpp_update_domain2D_r8_4dv
     module procedure mpp_update_domain2D_r8_5dv
#ifdef OVERLOAD_C8
     module procedure mpp_update_domain2D_c8_2d
     module procedure mpp_update_domain2D_c8_3d
     module procedure mpp_update_domain2D_c8_4d
     module procedure mpp_update_domain2D_c8_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_update_domain2D_i8_2d
     module procedure mpp_update_domain2D_i8_3d
     module procedure mpp_update_domain2D_i8_4d
     module procedure mpp_update_domain2D_i8_5d
!!$     module procedure mpp_update_domain2D_l8_2d
!!$     module procedure mpp_update_domain2D_l8_3d
!!$     module procedure mpp_update_domain2D_l8_4d
!!$     module procedure mpp_update_domain2D_l8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_update_domain2D_r4_2d
     module procedure mpp_update_domain2D_r4_3d
     module procedure mpp_update_domain2D_r4_4d
     module procedure mpp_update_domain2D_r4_5d
     module procedure mpp_update_domain2D_r4_2dv
     module procedure mpp_update_domain2D_r4_3dv
     module procedure mpp_update_domain2D_r4_4dv
     module procedure mpp_update_domain2D_r4_5dv
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_update_domain2D_c4_2d
     module procedure mpp_update_domain2D_c4_3d
     module procedure mpp_update_domain2D_c4_4d
     module procedure mpp_update_domain2D_c4_5d
#endif
     module procedure mpp_update_domain2D_i4_2d
     module procedure mpp_update_domain2D_i4_3d
     module procedure mpp_update_domain2D_i4_4d
     module procedure mpp_update_domain2D_i4_5d
!!$     module procedure mpp_update_domain2D_l4_2d
!!$     module procedure mpp_update_domain2D_l4_3d
!!$     module procedure mpp_update_domain2D_l4_4d
!!$     module procedure mpp_update_domain2D_l4_5d
  end interface

!--------------------------------------------------------------
!bnc: for adjoint update
!--------------------------------------------------------------
!!$  interface mpp_update_domains_ad
!!$     module procedure mpp_update_domain2D_ad_r8_2d
!!$     module procedure mpp_update_domain2D_ad_r8_3d
!!$     module procedure mpp_update_domain2D_ad_r8_4d
!!$     module procedure mpp_update_domain2D_ad_r8_5d
!!$     module procedure mpp_update_domain2D_ad_r8_2dv
!!$     module procedure mpp_update_domain2D_ad_r8_3dv
!!$     module procedure mpp_update_domain2D_ad_r8_4dv
!!$     module procedure mpp_update_domain2D_ad_r8_5dv
!!$#ifdef OVERLOAD_C8
!!$     module procedure mpp_update_domain2D_ad_c8_2d
!!$     module procedure mpp_update_domain2D_ad_c8_3d
!!$     module procedure mpp_update_domain2D_ad_c8_4d
!!$     module procedure mpp_update_domain2D_ad_c8_5d
!!$#endif
!!$#ifndef no_8byte_integers
!!$     module procedure mpp_update_domain2D_ad_i8_2d
!!$     module procedure mpp_update_domain2D_ad_i8_3d
!!$     module procedure mpp_update_domain2D_ad_i8_4d
!!$     module procedure mpp_update_domain2D_ad_i8_5d
!!$     module procedure mpp_update_domain2D_ad_l8_2d
!!$     module procedure mpp_update_domain2D_ad_l8_3d
!!$     module procedure mpp_update_domain2D_ad_l8_4d
!!$     module procedure mpp_update_domain2D_ad_l8_5d
!!$#endif
!!$#ifdef OVERLOAD_R4
!!$     module procedure mpp_update_domain2D_ad_r4_2d
!!$     module procedure mpp_update_domain2D_ad_r4_3d
!!$     module procedure mpp_update_domain2D_ad_r4_4d
!!$     module procedure mpp_update_domain2D_ad_r4_5d
!!$     module procedure mpp_update_domain2D_ad_r4_2dv
!!$     module procedure mpp_update_domain2D_ad_r4_3dv
!!$     module procedure mpp_update_domain2D_ad_r4_4dv
!!$     module procedure mpp_update_domain2D_ad_r4_5dv
!!$#endif
!!$#ifdef OVERLOAD_C4
!!$     module procedure mpp_update_domain2D_ad_c4_2d
!!$     module procedure mpp_update_domain2D_ad_c4_3d
!!$     module procedure mpp_update_domain2D_ad_c4_4d
!!$     module procedure mpp_update_domain2D_ad_c4_5d
!!$#endif
!!$     module procedure mpp_update_domain2D_ad_i4_2d
!!$     module procedure mpp_update_domain2D_ad_i4_3d
!!$     module procedure mpp_update_domain2D_ad_i4_4d
!!$     module procedure mpp_update_domain2D_ad_i4_5d
!!$  end interface
!bnc


  interface mpp_do_update
     module procedure mpp_do_update_r8_3d
     module procedure mpp_do_update_r8_3dv
#ifdef OVERLOAD_C8
     module procedure mpp_do_update_c8_3d
#endif
#ifndef no_8byte_integers
     module procedure mpp_do_update_i8_3d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_do_update_r4_3d
     module procedure mpp_do_update_r4_3dv
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_do_update_c4_3d
#endif
     module procedure mpp_do_update_i4_3d
  end interface

  interface mpp_do_check
     module procedure mpp_do_check_r8_3d
     module procedure mpp_do_check_r8_3dv
#ifdef OVERLOAD_C8
     module procedure mpp_do_check_c8_3d
#endif
#ifndef no_8byte_integers
     module procedure mpp_do_check_i8_3d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_do_check_r4_3d
     module procedure mpp_do_check_r4_3dv
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_do_check_c4_3d
#endif
     module procedure mpp_do_check_i4_3d
  end interface


!-------------------------------------------------------
!bnc  for adjoint do_update
!-------------------------------------------------------
!!$  interface mpp_do_update_ad
!!$     module procedure mpp_do_update_ad_r8_3d
!!$     module procedure mpp_do_update_ad_r8_3dv
!!$#ifdef OVERLOAD_C8
!!$     module procedure mpp_do_update_ad_c8_3d
!!$#endif
!!$#ifndef no_8byte_integers
!!$     module procedure mpp_do_update_ad_i8_3d
!!$#endif
!!$#ifdef OVERLOAD_R4
!!$     module procedure mpp_do_update_ad_r4_3d
!!$     module procedure mpp_do_update_ad_r4_3dv
!!$#endif
!!$#ifdef OVERLOAD_C4
!!$     module procedure mpp_do_update_ad_c4_3d
!!$#endif
!!$     module procedure mpp_do_update_ad_i4_3d
!!$  end interface
!bnc

! <INTERFACE NAME="mpp_get_boundary">
! <OVERVIEW>
!    Get the boundary data for symmetric domain when the data is at C, E, or N-cell center
! </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_get_boundary</TT> is used to get the boundary data for symmetric domain 
!        when the data is at C, E, or N-cell center. For cubic grid, the data should 
!        always at C-cell center. 
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_get_boundary
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_get_boundary
!  </TEMPLATE>
! </INTERFACE>
  interface mpp_get_boundary
     module procedure mpp_get_boundary_r8_2d
     module procedure mpp_get_boundary_r8_3d
     module procedure mpp_get_boundary_r8_4d
     module procedure mpp_get_boundary_r8_5d
     module procedure mpp_get_boundary_r8_2dv
     module procedure mpp_get_boundary_r8_3dv
     module procedure mpp_get_boundary_r8_4dv
     module procedure mpp_get_boundary_r8_5dv
#ifdef OVERLOAD_R4
     module procedure mpp_get_boundary_r4_2d
     module procedure mpp_get_boundary_r4_3d
     module procedure mpp_get_boundary_r4_4d
     module procedure mpp_get_boundary_r4_5d
     module procedure mpp_get_boundary_r4_2dv
     module procedure mpp_get_boundary_r4_3dv
     module procedure mpp_get_boundary_r4_4dv
     module procedure mpp_get_boundary_r4_5dv
#endif
  end interface

  interface mpp_do_get_boundary
     module procedure mpp_do_get_boundary_r8_3d
     module procedure mpp_do_get_boundary_r8_3dv
#ifdef OVERLOAD_R4
     module procedure mpp_do_get_boundary_r4_3d
     module procedure mpp_do_get_boundary_r4_3dv
#endif
  end interface

! <INTERFACE NAME="mpp_redistribute">
!  <OVERVIEW>
!    Reorganization of distributed global arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_redistribute</TT> is used to reorganize a distributed
!    array.  <TT>MPP_TYPE_</TT> can be of type <TT>integer</TT>,
!    <TT>complex</TT>, or <TT>real</TT>; of 4-byte or 8-byte kind; of rank
!    up to 5.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_redistribute( domain_in, field_in, domain_out, field_out )
!  </TEMPLATE>
!  <IN NAME="field_in" TYPE="MPP_TYPE_">
!    <TT>field_in</TT> is dimensioned on the data domain of <TT>domain_in</TT>.
!  </IN>
!  <OUT NAME="field_out" TYPE="MPP_TYPE_">
!    <TT>field_out</TT> on the data domain of <TT>domain_out</TT>.
!  </OUT>
! </INTERFACE>
  interface mpp_redistribute
     module procedure mpp_redistribute_r8_2D
     module procedure mpp_redistribute_r8_3D
     module procedure mpp_redistribute_r8_4D
     module procedure mpp_redistribute_r8_5D
#ifdef OVERLOAD_C8
     module procedure mpp_redistribute_c8_2D
     module procedure mpp_redistribute_c8_3D
     module procedure mpp_redistribute_c8_4D
     module procedure mpp_redistribute_c8_5D
#endif
#ifndef no_8byte_integers
     module procedure mpp_redistribute_i8_2D
     module procedure mpp_redistribute_i8_3D
     module procedure mpp_redistribute_i8_4D
     module procedure mpp_redistribute_i8_5D
!!$     module procedure mpp_redistribute_l8_2D
!!$     module procedure mpp_redistribute_l8_3D
!!$     module procedure mpp_redistribute_l8_4D
!!$     module procedure mpp_redistribute_l8_5D
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_redistribute_r4_2D
     module procedure mpp_redistribute_r4_3D
     module procedure mpp_redistribute_r4_4D
     module procedure mpp_redistribute_r4_5D
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_redistribute_c4_2D
     module procedure mpp_redistribute_c4_3D
     module procedure mpp_redistribute_c4_4D
     module procedure mpp_redistribute_c4_5D
#endif
     module procedure mpp_redistribute_i4_2D
     module procedure mpp_redistribute_i4_3D
     module procedure mpp_redistribute_i4_4D
     module procedure mpp_redistribute_i4_5D
!!$     module procedure mpp_redistribute_l4_2D
!!$     module procedure mpp_redistribute_l4_3D
!!$     module procedure mpp_redistribute_l4_4D
!!$     module procedure mpp_redistribute_l4_5D
  end interface

  interface mpp_do_redistribute
     module procedure mpp_do_redistribute_r8_3D
#ifdef OVERLOAD_C8
     module procedure mpp_do_redistribute_c8_3D
#endif
#ifndef no_8byte_integers
     module procedure mpp_do_redistribute_i8_3D
     module procedure mpp_do_redistribute_l8_3D
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_do_redistribute_r4_3D
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_do_redistribute_c4_3D
#endif
     module procedure mpp_do_redistribute_i4_3D
     module procedure mpp_do_redistribute_l4_3D
  end interface


! <INTERFACE NAME="mpp_check_field">
!   <OVERVIEW>
!     Parallel checking between two ensembles which run
!     on different set pes at the same time.
!   </OVERVIEW>
!   <DESCRIPTION>
!     There are two forms for the <TT>mpp_check_field</TT> call. The 2D
!     version is generally to be used and 3D version is  built by repeated calls to the
!     2D version.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_check_field(field_in, pelist1, pelist2, domain, mesg, &
!                                w_halo, s_halo, e_halo, n_halo, force_abort  )
!   </TEMPLATE>
!   <IN NAME="field_in" >
!     Field to be checked
!   </IN>
!   <IN NAME="pelist1, pelist2">
!     Pelist of the two ensembles to be compared
!   </IN>
!   <IN NAME="domain">
!     Domain of current pe
!   </IN>
!   <IN NAME="mesg" >
!     Message to be printed out
!   </IN>
!   <IN NAME="w_halo, s_halo, e_halo, n_halo">
!     Halo size to be checked. Default value is 0.
!   </IN>
!   <IN NAME="force_abort">
!     When true, abort program when any difference found. Default value is false.
!   </IN>
! </INTERFACE>

  interface mpp_check_field
     module procedure mpp_check_field_2D
     module procedure mpp_check_field_3D
  end interface

!***********************************************************************
!
!         public interface from mpp_domains_reduce.h
!
!***********************************************************************

! <INTERFACE NAME="mpp_global_field">
!  <OVERVIEW>
!    Fill in a global array from domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_field</TT> is used to get an entire
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!    of 4-byte or 8-byte kind; of rank up to 5.
!    
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_field</TT>, and each will have a complete global field
!    at the end. Please note that a global array of rank 3 or higher could
!    occupy a lot of memory.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_global_field( domain, local, global, flags )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="local" TYPE="MPP_TYPE_">
!    <TT>local</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <OUT NAME="global" TYPE="MPP_TYPE_">
!    <TT>global</TT> is dimensioned on the corresponding global domain.
!  </OUT>
!  <IN NAME="flags" TYPE="integer">
!    <TT>flags</TT> can be given the value <TT>XONLY</TT> or
!    <TT>YONLY</TT>, to specify a globalization on one axis only.
!  </IN>
! </INTERFACE>
  interface mpp_global_field
     module procedure mpp_global_field2D_r8_2d
     module procedure mpp_global_field2D_r8_3d
     module procedure mpp_global_field2D_r8_4d
     module procedure mpp_global_field2D_r8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_global_field2D_c8_2d
     module procedure mpp_global_field2D_c8_3d
     module procedure mpp_global_field2D_c8_4d
     module procedure mpp_global_field2D_c8_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_field2D_i8_2d
     module procedure mpp_global_field2D_i8_3d
     module procedure mpp_global_field2D_i8_4d
     module procedure mpp_global_field2D_i8_5d
     module procedure mpp_global_field2D_l8_2d
     module procedure mpp_global_field2D_l8_3d
     module procedure mpp_global_field2D_l8_4d
     module procedure mpp_global_field2D_l8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_global_field2D_r4_2d
     module procedure mpp_global_field2D_r4_3d
     module procedure mpp_global_field2D_r4_4d
     module procedure mpp_global_field2D_r4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_global_field2D_c4_2d
     module procedure mpp_global_field2D_c4_3d
     module procedure mpp_global_field2D_c4_4d
     module procedure mpp_global_field2D_c4_5d
#endif
     module procedure mpp_global_field2D_i4_2d
     module procedure mpp_global_field2D_i4_3d
     module procedure mpp_global_field2D_i4_4d
     module procedure mpp_global_field2D_i4_5d
     module procedure mpp_global_field2D_l4_2d
     module procedure mpp_global_field2D_l4_3d
     module procedure mpp_global_field2D_l4_4d
     module procedure mpp_global_field2D_l4_5d
  end interface

  interface mpp_do_global_field
     module procedure mpp_do_global_field2D_r8_3d
#ifdef OVERLOAD_C8
     module procedure mpp_do_global_field2D_c8_3d
#endif
#ifndef no_8byte_integers
     module procedure mpp_do_global_field2D_i8_3d
     module procedure mpp_do_global_field2D_l8_3d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_do_global_field2D_r4_3d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_do_global_field2D_c4_3d
#endif
     module procedure mpp_do_global_field2D_i4_3d
     module procedure mpp_do_global_field2D_l4_3d
  end interface

! <INTERFACE NAME="mpp_global_max">
!  <OVERVIEW>
!    Global max/min of domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_max</TT> is used to get the maximum value of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>integer</TT> or <TT>real</TT>; of 4-byte or 8-byte kind; of rank
!    up to 5. The dimension of <TT>locus</TT> must equal the rank of
!    <TT>field</TT>.
!    
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_max</TT>, and each will have the result upon exit.
!    
!    The function <TT>mpp_global_min</TT>, with an identical syntax. is
!    also available.
!  </DESCRIPTION>
!  <TEMPLATE>
!    mpp_global_max( domain, field, locus )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="field" TYPE="MPP_TYPE_">  
!    <TT>field</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <OUT NAME="locus" TYPE="integer" DIM="(:)">
!    <TT>locus</TT>, if present, can be used to retrieve the location of
!    the maximum (as in the <TT>MAXLOC</TT> intrinsic of f90).
!  </OUT>
! </INTERFACE>

  interface mpp_global_max
     module procedure mpp_global_max_r8_2d
     module procedure mpp_global_max_r8_3d
     module procedure mpp_global_max_r8_4d
     module procedure mpp_global_max_r8_5d
#ifdef OVERLOAD_R4
     module procedure mpp_global_max_r4_2d
     module procedure mpp_global_max_r4_3d
     module procedure mpp_global_max_r4_4d
     module procedure mpp_global_max_r4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_max_i8_2d
     module procedure mpp_global_max_i8_3d
     module procedure mpp_global_max_i8_4d
     module procedure mpp_global_max_i8_5d
#endif
     module procedure mpp_global_max_i4_2d
     module procedure mpp_global_max_i4_3d
     module procedure mpp_global_max_i4_4d
     module procedure mpp_global_max_i4_5d
  end interface

  interface mpp_global_min
     module procedure mpp_global_min_r8_2d
     module procedure mpp_global_min_r8_3d
     module procedure mpp_global_min_r8_4d
     module procedure mpp_global_min_r8_5d
#ifdef OVERLOAD_R4
     module procedure mpp_global_min_r4_2d
     module procedure mpp_global_min_r4_3d
     module procedure mpp_global_min_r4_4d
     module procedure mpp_global_min_r4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_min_i8_2d
     module procedure mpp_global_min_i8_3d
     module procedure mpp_global_min_i8_4d
     module procedure mpp_global_min_i8_5d
#endif
     module procedure mpp_global_min_i4_2d
     module procedure mpp_global_min_i4_3d
     module procedure mpp_global_min_i4_4d
     module procedure mpp_global_min_i4_5d
  end interface

! <INTERFACE NAME="mpp_global_sum">
!  <OVERVIEW>
!    Global sum of domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_sum</TT> is used to get the sum of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>integer</TT>, <TT>complex</TT>, or <TT>real</TT>; of 4-byte or
!    8-byte kind; of rank up to 5.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_global_sum( domain, field, flags )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="field" TYPE="MPP_TYPE_">
!    <TT>field</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <IN NAME="flags" TYPE="integer">
!    <TT>flags</TT>, if present, must have the value
!    <TT>BITWISE_EXACT_SUM</TT>. This produces a sum that is guaranteed to
!    produce the identical result irrespective of how the domain is
!    decomposed. This method does the sum first along the ranks beyond 2,
!    and then calls <LINK
!    SRC="#mpp_global_field"><TT>mpp_global_field</TT></LINK> to produce a
!    global 2D array which is then summed. The default method, which is
!    considerably faster, does a local sum followed by <LINK
!    SRC="mpp.html#mpp_sum"><TT>mpp_sum</TT></LINK> across the domain
!    decomposition.
!  </IN>
!  <NOTE>
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_sum</TT>, and each will have the result upon exit.
!  </NOTE>
! </INTERFACE>

  interface mpp_global_sum
     module procedure mpp_global_sum_r8_2d
     module procedure mpp_global_sum_r8_3d
     module procedure mpp_global_sum_r8_4d
     module procedure mpp_global_sum_r8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_global_sum_c8_2d
     module procedure mpp_global_sum_c8_3d
     module procedure mpp_global_sum_c8_4d
     module procedure mpp_global_sum_c8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_global_sum_r4_2d
     module procedure mpp_global_sum_r4_3d
     module procedure mpp_global_sum_r4_4d
     module procedure mpp_global_sum_r4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_global_sum_c4_2d
     module procedure mpp_global_sum_c4_3d
     module procedure mpp_global_sum_c4_4d
     module procedure mpp_global_sum_c4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_sum_i8_2d
     module procedure mpp_global_sum_i8_3d
     module procedure mpp_global_sum_i8_4d
     module procedure mpp_global_sum_i8_5d
#endif
     module procedure mpp_global_sum_i4_2d
     module procedure mpp_global_sum_i4_3d
     module procedure mpp_global_sum_i4_4d
     module procedure mpp_global_sum_i4_5d
  end interface

!gag
  interface mpp_global_sum_tl
     module procedure mpp_global_sum_tl_r8_2d
     module procedure mpp_global_sum_tl_r8_3d
     module procedure mpp_global_sum_tl_r8_4d
     module procedure mpp_global_sum_tl_r8_5d
#ifdef OVERLOAD_C8
     module procedure mpp_global_sum_tl_c8_2d
     module procedure mpp_global_sum_tl_c8_3d
     module procedure mpp_global_sum_tl_c8_4d
     module procedure mpp_global_sum_tl_c8_5d
#endif
#ifdef OVERLOAD_R4
     module procedure mpp_global_sum_tl_r4_2d
     module procedure mpp_global_sum_tl_r4_3d
     module procedure mpp_global_sum_tl_r4_4d
     module procedure mpp_global_sum_tl_r4_5d
#endif
#ifdef OVERLOAD_C4
     module procedure mpp_global_sum_tl_c4_2d
     module procedure mpp_global_sum_tl_c4_3d
     module procedure mpp_global_sum_tl_c4_4d
     module procedure mpp_global_sum_tl_c4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_sum_tl_i8_2d
     module procedure mpp_global_sum_tl_i8_3d
     module procedure mpp_global_sum_tl_i8_4d
     module procedure mpp_global_sum_tl_i8_5d
#endif
     module procedure mpp_global_sum_tl_i4_2d
     module procedure mpp_global_sum_tl_i4_3d
     module procedure mpp_global_sum_tl_i4_4d
     module procedure mpp_global_sum_tl_i4_5d
  end interface
!gag

!bnc
!!$  interface mpp_global_sum_ad
!!$     module procedure mpp_global_sum_ad_r8_2d
!!$     module procedure mpp_global_sum_ad_r8_3d
!!$     module procedure mpp_global_sum_ad_r8_4d
!!$     module procedure mpp_global_sum_ad_r8_5d
!!$#ifdef OVERLOAD_C8
!!$     module procedure mpp_global_sum_ad_c8_2d
!!$     module procedure mpp_global_sum_ad_c8_3d
!!$     module procedure mpp_global_sum_ad_c8_4d
!!$     module procedure mpp_global_sum_ad_c8_5d
!!$#endif
!!$#ifdef OVERLOAD_R4
!!$     module procedure mpp_global_sum_ad_r4_2d
!!$     module procedure mpp_global_sum_ad_r4_3d
!!$     module procedure mpp_global_sum_ad_r4_4d
!!$     module procedure mpp_global_sum_ad_r4_5d
!!$#endif
!!$#ifdef OVERLOAD_C4
!!$     module procedure mpp_global_sum_ad_c4_2d
!!$     module procedure mpp_global_sum_ad_c4_3d
!!$     module procedure mpp_global_sum_ad_c4_4d
!!$     module procedure mpp_global_sum_ad_c4_5d
!!$#endif
!!$#ifndef no_8byte_integers
!!$     module procedure mpp_global_sum_ad_i8_2d
!!$     module procedure mpp_global_sum_ad_i8_3d
!!$     module procedure mpp_global_sum_ad_i8_4d
!!$     module procedure mpp_global_sum_ad_i8_5d
!!$#endif
!!$     module procedure mpp_global_sum_ad_i4_2d
!!$     module procedure mpp_global_sum_ad_i4_3d
!!$     module procedure mpp_global_sum_ad_i4_4d
!!$     module procedure mpp_global_sum_ad_i4_5d
!!$  end interface
!bnc

!***********************************************************************
!
!            public interface from mpp_domain_util.h
!
!***********************************************************************

  ! <INTERFACE NAME="mpp_get_neighbor_pe">
  !  <OVERVIEW>
  !    Retrieve PE number of a neighboring domain.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Given a 1-D or 2-D domain decomposition, this call allows users to retrieve 
  !    the PE number of an adjacent PE-domain while taking into account that the 
  !    domain may have holes (masked) and/or have cyclic boundary conditions and/or a 
  !    folded edge. Which PE-domain will be retrived will depend on "direction": 
  !    +1 (right) or -1 (left) for a 1-D domain decomposition and either NORTH, SOUTH, 
  !    EAST, WEST, NORTH_EAST, SOUTH_EAST, SOUTH_WEST, or NORTH_WEST for a 2-D 
  !    decomposition. If no neighboring domain exists (masked domain), then the 
  !    returned "pe" value will be set to NULL_PE.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_neighbor_pe( domain1d, direction=+1   , pe)
  !    call mpp_get_neighbor_pe( domain2d, direction=NORTH, pe)
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_neighbor_pe
     module procedure mpp_get_neighbor_pe_1d
     module procedure mpp_get_neighbor_pe_2d
  end interface 
  ! <INTERFACE NAME="operator">
  !  <OVERVIEW>
  !    Equality/inequality operators for domaintypes.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The module provides public operators to check for
  !    equality/inequality of domaintypes, e.g:
  !    
  !    <PRE>
  !    type(domain1D) :: a, b
  !    type(domain2D) :: c, d
  !    ...
  !    if( a.NE.b )then
  !        ...
  !    end if
  !    if( c==d )then
  !        ...
  !    end if
  !    </PRE>
  !    
  !    Domains are considered equal if and only if the start and end
  !    indices of each of their component global, data and compute domains
  !    are equal.
  !  </DESCRIPTION>
  ! </INTERFACE>
  interface operator(.EQ.)
     module procedure mpp_domain1D_eq
     module procedure mpp_domain2D_eq
  end interface

  interface operator(.NE.)
     module procedure mpp_domain1D_ne
     module procedure mpp_domain2D_ne
  end interface

  ! <INTERFACE NAME="mpp_get_compute_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the compute domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the compute domains
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_compute_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_compute_domain
     module procedure mpp_get_compute_domain1D
     module procedure mpp_get_compute_domain2D
  end interface

  ! <INTERFACE NAME="mpp_get_compute_domains">
  !  <OVERVIEW>
  !    Retrieve the entire array of compute domain extents associated with a decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Retrieve the entire array of compute domain extents associated with a decomposition.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_compute_domains( domain, xbegin, xend, xsize, &
  !                                                ybegin, yend, ysize )
  !  </TEMPLATE>
  !  <IN NAME="domain" TYPE="type(domain2D)"></IN>
  !  <OUT NAME="xbegin,ybegin" TYPE="integer" DIM="(:)"></OUT>
  !  <OUT NAME="xend,yend" TYPE="integer" DIM="(:)"></OUT>
  !  <OUT NAME="xsize,ysize" TYPE="integer" DIM="(:)"></OUT>
  ! </INTERFACE>
  interface mpp_get_compute_domains
     module procedure mpp_get_compute_domains1D
     module procedure mpp_get_compute_domains2D
  end interface

  ! <INTERFACE NAME="mpp_get_data_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the data domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the data domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_data_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_data_domain
     module procedure mpp_get_data_domain1D
     module procedure mpp_get_data_domain2D
  end interface

  ! <INTERFACE NAME="mpp_get_global_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the global domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the global domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_global_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_global_domain
     module procedure mpp_get_global_domain1D
     module procedure mpp_get_global_domain2D
  end interface

  ! <INTERFACE NAME="mpp_get_memory_domain">
  !  <OVERVIEW>
  !    These routines retrieve the axis specifications associated with the memory domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    retrieve the axis specifications associated with the memory domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_memory_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_get_memory_domain
     module procedure mpp_get_memory_domain1D
     module procedure mpp_get_memory_domain2D
  end interface

  interface mpp_get_domain_extents
     module procedure mpp_get_domain_extents1D
     module procedure mpp_get_domain_extents2D
  end interface

  ! <INTERFACE NAME="mpp_set_compute_domain">
  !  <OVERVIEW>
  !    These routines set the axis specifications associated with the compute domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    set the axis specifications associated with the compute domains
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_compute_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_set_compute_domain
     module procedure mpp_set_compute_domain1D
     module procedure mpp_set_compute_domain2D
  end interface

  ! <INTERFACE NAME="mpp_set_data_domain">
  !  <OVERVIEW>
  !    These routines set the axis specifications associated with the data domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    set the axis specifications associated with the data domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_data_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_set_data_domain
     module procedure mpp_set_data_domain1D
     module procedure mpp_set_data_domain2D
  end interface

  ! <INTERFACE NAME="mpp_set_global_domain">
  !  <OVERVIEW>
  !    These routines set the axis specifications associated with the global domains.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The domain is a derived type with private elements. These routines 
  !    set the axis specifications associated with the global domains.
  !    The 2D version of these is a simple extension of 1D.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_set_global_domain
  !  </TEMPLATE>
  ! </INTERFACE>
  interface mpp_set_global_domain
     module procedure mpp_set_global_domain1D
     module procedure mpp_set_global_domain2D
  end interface


  ! <INTERFACE NAME="mpp_get_pelist">
  !  <OVERVIEW>
  !    Retrieve list of PEs associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The 1D version of this call returns an array of the PEs assigned to this 1D domain
  !    decomposition. In addition the optional argument <TT>pos</TT> may be
  !    used to retrieve the 0-based position of the domain local to the
  !    calling PE, i.e <TT>domain%list(pos)%pe</TT> is the local PE,
  !    as returned by <LINK SRC="mpp.html#mpp_pe"><TT>mpp_pe()</TT></LINK>.
  !    The 2D version of this call is identical to 1D version.
  !  </DESCRIPTION>
  !  <IN NAME="domain"></IN>
  !  <OUT NAME="pelist"></OUT>
  !  <OUT NAME="pos"></OUT>
  ! </INTERFACE>
  interface mpp_get_pelist
     module procedure mpp_get_pelist1D
     module procedure mpp_get_pelist2D
  end interface

  ! <INTERFACE NAME="mpp_get_layout">
  !  <OVERVIEW>
  !    Retrieve layout associated with a domain decomposition.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    The 1D version of this call returns the number of divisions that was assigned to this
  !    decomposition axis. The 2D version of this call returns an array of
  !    dimension 2 holding the results on two axes.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_get_layout( domain, layout )
  !  </TEMPLATE>
  !  <IN NAME="domain"></IN>
  !  <OUT NAME="layout"></OUT>
  ! </INTERFACE>
  interface mpp_get_layout
     module procedure mpp_get_layout1D
     module procedure mpp_get_layout2D
  end interface

  ! <INTERFACE NAME="mpp_nullify_domain_list">
  !  <OVERVIEW>
  !    nullify domain list.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    Nullify domain list. This interface is needed in mpp_domains_test.
  !    1-D case can be added in if needed.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_nullify_domain_list( domain)
  !  </TEMPLATE>
  !  <INOUT NAME="domain"></INOUT>
  ! </INTERFACE>
  interface mpp_nullify_domain_list
     module procedure nullify_domain2d_list
  end interface  

  !--- version information variables
  character(len=128), public :: version= &
       '$Id: mpp_domains.F90,v 18.0 2010/03/02 23:56:33 fms Exp $'
  character(len=128), public :: tagname= &
       '$Name: testing $'


contains

#include <mpp_domains_util.inc>
#include <mpp_domains_comm.inc>
#include <mpp_domains_define.inc>
#include <mpp_domains_misc.inc>
#include <mpp_domains_reduce.inc>


end module mpp_domains_mod

! <INFO>

!   <COMPILER NAME="">     
!     Any module or program unit using <TT>mpp_domains_mod</TT>
!     must contain the line

!     <PRE>
!     use mpp_domains_mod
!     </PRE>

!     <TT>mpp_domains_mod</TT> <TT>use</TT>s <LINK
!     SRC="mpp.html">mpp_mod</LINK>, and therefore is subject to the <LINK
!     SRC="mpp.html#COMPILING AND LINKING SOURCE">compiling and linking requirements of that module.</LINK>
!   </COMPILER>
!   <PRECOMP FLAG="">      
!     <TT>mpp_domains_mod</TT> uses standard f90, and has no special
!     requirements. There are some OS-dependent
!     pre-processor directives that you might need to modify on
!     non-SGI/Cray systems and compilers. The <LINK
!     SRC="mpp.html#PORTABILITY">portability of mpp_mod</LINK>
!     obviously is a constraint, since this module is built on top of
!     it. Contact me, Balaji, SGI/GFDL, with questions.
!   </PRECOMP> 
!   <LOADER FLAG="">       
!     The <TT>mpp_domains</TT> source consists of the main source file
!     <TT>mpp_domains.F90</TT> and also requires the following include files:
!    <PRE>
!     <TT>fms_platform.h</TT>
!     <TT>mpp_update_domains2D.h</TT>
!     <TT>mpp_global_reduce.h</TT>
!     <TT>mpp_global_sum.h</TT>
!     <TT>mpp_global_field.h</TT>
!    </PRE>
!    GFDL users can check it out of the main CVS repository as part of
!    the <TT>mpp</TT> CVS module. The current public tag is <TT>galway</TT>.
!    External users can download the latest <TT>mpp</TT> package <LINK SRC=
!    "ftp://ftp.gfdl.gov/pub/vb/mpp/mpp.tar.Z">here</LINK>. Public access
!    to the GFDL CVS repository will soon be made available.

!   </LOADER>

! </INFO>
