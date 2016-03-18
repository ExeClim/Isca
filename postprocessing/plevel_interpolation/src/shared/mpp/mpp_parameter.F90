module mpp_parameter_mod
#include <fms_platform.h>

  implicit none
  private

  character(len=128), public :: version= &
       '$Id mpp_parameter.F90 $'
  character(len=128), public :: tagname= &
       '$Name: testing $'

  !--- public paramters which is used by mpp_mod and its components. 
  !--- All othere modules should import these parameters from mpp_mod. 
  public :: MAXPES, MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE, NOTE, WARNING, FATAL
  public :: MPP_WAIT, MPP_READY, MAX_CLOCKS, MAX_EVENT_TYPES, MAX_EVENTS, MPP_CLOCK_SYNC
  public :: MPP_CLOCK_DETAILED, CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER
  public :: CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA, MAX_BINS, PESET_MAX
  public :: EVENT_ALLREDUCE, EVENT_BROADCAST, EVENT_RECV, EVENT_SEND, EVENT_WAIT

  !--- public paramters which is used by mpp_domains_mod and its components. 
  !--- All othere modules should import these parameters from mpp_domains_mod. 
  public :: GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, BGRID_NE, BGRID_SW, CGRID_NE, CGRID_SW
  public :: DGRID_NE, DGRID_SW, FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  public :: WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
  public :: MPP_DOMAIN_TIME, WEST, EAST, SOUTH, NORTH, SCALAR_BIT, SCALAR_PAIR
  public :: NORTH_EAST, SOUTH_EAST, SOUTH_WEST, NORTH_WEST
  public :: AGRID, GLOBAL, CYCLIC, DOMAIN_ID_BASE, CENTER, CORNER
  public :: MAX_DOMAIN_FIELDS, MAX_TILES
  public :: ZERO, NINETY, MINUS_NINETY, ONE_HUNDRED_EIGHTY

  !--- public paramters which is used by mpp_domains_mod and its components. 
  !--- All othere modules should import these parameters from mpp_io_mod. 
  public :: MPP_WRONLY, MPP_RDONLY, MPP_APPEND, MPP_OVERWR, MPP_ASCII, MPP_IEEE32
  public :: MPP_NATIVE, MPP_NETCDF, MPP_SEQUENTIAL, MPP_DIRECT, MPP_SINGLE, MPP_MULTI
  public :: MPP_DELETE, MPP_COLLECT, NULLUNIT, NULLTIME
  public :: MAX_FILE_SIZE, ROOT_GLOBAL, GLOBAL_ROOT_ONLY

  !--- The following paramters are used by mpp_mod and its components.
  integer, parameter :: MAXPES=2048            !used for dimensioning stuff that might be indexed by pe
  integer, parameter :: MPP_VERBOSE=1, MPP_DEBUG=2
  integer, parameter :: ALL_PES=-1, ANY_PE=-2, NULL_PE=-3
  integer, parameter :: NOTE=0, WARNING=1, FATAL=2
  integer, parameter :: MAX_CLOCKS=400, MAX_EVENT_TYPES=5, MAX_EVENTS=40000
  integer, parameter :: EVENT_ALLREDUCE=1, EVENT_BROADCAST=2, EVENT_RECV=3, EVENT_SEND=4, EVENT_WAIT=5
  integer, parameter :: MPP_CLOCK_SYNC=1, MPP_CLOCK_DETAILED=2
  !--- predefined clock granularities, but you can use any integer
  !--- using CLOCK_LOOP and above may distort coarser-grain measurements
  integer, parameter :: CLOCK_COMPONENT=1      !component level, e.g model, exchange
  integer, parameter :: CLOCK_SUBCOMPONENT=11  !top level within a model component, e.g dynamics, physics
  integer, parameter :: CLOCK_MODULE_DRIVER=21 !module driver level, e.g adriver that calls multiple 
                                               !related physics routines
  integer, parameter :: CLOCK_MODULE=31        !module level, e.g main subroutine of a physics module
  integer, parameter :: CLOCK_ROUTINE=41       !level of individual subroutine or function
  integer, parameter :: CLOCK_LOOP=51          !loops or blocks within a routine
  integer, parameter :: CLOCK_INFRA=61         !infrastructure level, e.g halo update
  integer, parameter :: MAX_BINS=20
  integer, parameter :: PESET_MAX=32           !should be .LE. max num of MPI communicators
  integer(LONG_KIND), parameter :: MPP_WAIT=-1, MPP_READY=-2

  !--- The following paramters are used by mpp_domains_mod and its components.
  integer, parameter :: GLOBAL=0, CYCLIC=1
  integer, parameter :: WEST=2, EAST=3, SOUTH=4, NORTH=5, SCALAR_BIT=6, CENTER=7, CORNER=8
  integer, parameter :: SOUTH_WEST=7, SOUTH_EAST=8, NORTH_WEST=9, NORTH_EAST=10
  integer, parameter :: SEND=1, RECV=2
  integer, parameter :: GLOBAL_DATA_DOMAIN=2**GLOBAL, CYCLIC_GLOBAL_DOMAIN=2**CYCLIC
  integer, parameter :: AGRID=0, BGRID=1, CGRID=2, DGRID=3
  integer, parameter :: BGRID_NE=BGRID+2**NORTH+2**EAST
  integer, parameter :: BGRID_SW=BGRID+2**SOUTH+2**WEST
  integer, parameter :: CGRID_NE=CGRID+2**NORTH+2**EAST
  integer, parameter :: CGRID_SW=CGRID+2**SOUTH+2**WEST
  integer, parameter :: DGRID_NE=DGRID+2**NORTH+2**EAST
  integer, parameter :: DGRID_SW=DGRID+2**SOUTH+2**WEST
  integer, parameter :: FOLD_WEST_EDGE = 2**WEST, FOLD_EAST_EDGE = 2**EAST
  integer, parameter :: FOLD_SOUTH_EDGE=2**SOUTH, FOLD_NORTH_EDGE=2**NORTH
  integer, parameter :: WUPDATE=2**WEST, EUPDATE=2**EAST, SUPDATE=2**SOUTH, NUPDATE=2**NORTH
  integer, parameter :: XUPDATE=WUPDATE+EUPDATE, YUPDATE=SUPDATE+NUPDATE, SCALAR_PAIR=2**SCALAR_BIT
  integer, parameter :: ZERO=0, NINETY=90, MINUS_NINETY=-90, ONE_HUNDRED_EIGHTY=180

! DOMAIN_ID_BASE acts as a counter increment for domains as they are defined. It's used in
! combination with the flag parameter defined above to create a unique identifier for
! each Domain+flags combination. Therefore, the value of any flag must not exceed DOMAIN_ID_BASE.
! integer(LONG_KIND), parameter :: DOMAIN_ID_BASE=INT( 2**(4*LONG_KIND),KIND=LONG_KIND )
  integer(LONG_KIND), parameter :: DOMAIN_ID_BASE=Z'0000000100000000' ! Workaround for 64bit init problem
  integer, parameter :: NON_BITWISE_EXACT_SUM=0
  integer, parameter :: BITWISE_EXACT_SUM=1
  integer, parameter :: MPP_DOMAIN_TIME=MPP_DEBUG+1
  integer, parameter :: MAX_DOMAIN_FIELDS=100
  integer, parameter :: MAX_TILES=100

  !--- The following paramters are used by mpp_io_mod and its components.
  integer, parameter :: MPP_WRONLY=100, MPP_RDONLY=101, MPP_APPEND=102, MPP_OVERWR=103 !action on open
  integer, parameter :: MPP_ASCII=200,  MPP_IEEE32=201, MPP_NATIVE=202, MPP_NETCDF=203 !format
  integer, parameter :: MPP_SEQUENTIAL=300, MPP_DIRECT=301 !access
  integer, parameter :: MPP_SINGLE=400, MPP_MULTI=401      !threading, fileset
  integer, parameter :: MPP_DELETE=501, MPP_COLLECT=502    !action on close
  integer, parameter :: NULLUNIT=-1                        !returned by PEs not participating in 
                                                           !IO after a collective call with threading
                                                           !equal to MPP_SINGLE
  integer, parameter :: ROOT_GLOBAL = 9
  integer, parameter :: GLOBAL_ROOT_ONLY = 2**ROOT_GLOBAL 
  real(DOUBLE_KIND), parameter :: NULLTIME=-1.
#ifdef LARGE_FILE
  integer(LONG_KIND), parameter :: MAX_FILE_SIZE = 4294967295
#else
  integer(LONG_KIND), parameter :: MAX_FILE_SIZE = 2147483647
#endif

  !#####################################################################

end module mpp_parameter_mod
