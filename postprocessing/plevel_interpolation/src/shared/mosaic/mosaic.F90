module mosaic_mod

! <CONTACT EMAIL="Zhi.Liang@noaa.gov">
!   Zhi Liang
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!    <TT>mosaic_mod</TT> implements some utility routines to read mosaic information.
! </OVERVIEW>

! <DESCRIPTION>
!    <TT>mosaic_mod</TT> implements some utility routines to read mosaic information.
!    The information includes number of tiles and contacts in the mosaic, 
!    mosaic grid resolution of each tile, mosaic contact information, mosaic exchange
!    grid information. Each routine will call a C-version routine to get these information.
! </DESCRIPTION>

use mpp_mod, only : mpp_error, FATAL

implicit none
private


! --- public interface


public :: get_mosaic_ntiles
public :: get_mosaic_ncontacts
public :: get_mosaic_grid_sizes
public :: get_mosaic_contact
public :: get_mosaic_xgrid_size
public :: get_mosaic_xgrid
public :: calc_mosaic_grid_area

logical :: module_is_initialized = .true.
! version information varaible
 character(len=128) :: version = '$Id: mosaic.F90,v 15.0 2007/08/14 04:14:22 fms Exp $'
 character(len=128) :: tagname = '$Name: testing $'

contains

!#######################################################################

! <SUBROUTINE NAME="mosaic_init">
!   <OVERVIEW>
!     Initialize the mosaic_mod. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initialization routine for the mosaic module. It writes the 
!     version information to the log file.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mosaic_init ( )
!   </TEMPLATE>
subroutine mosaic_init() 

  if (module_is_initialized) return
  module_is_initialized = .TRUE.

!--------- write version number and namelist ------------------
!  call write_version_number (version, tagname)

end subroutine mosaic_init
! </SUBROUTINE>

!#######################################################################
! <FUNCTION NAME="get_mosaic_xgrid_size">
!   <OVERVIEW>
!     return exchange grid size of mosaic xgrid file.
!   </OVERVIEW>
!   <DESCRIPTION>
!     return exchange grid size of mosaic xgrid file.
!   </DESCRIPTION>
!   <TEMPLATE>
!    nxgrid = get_mosaic_xgrid_size(xgrid_file)
!   </TEMPLATE>
!   <IN NAME="xgrid_file" TYPE="character(len=*)">
!     The file that contains exchange grid information.
!   </IN>
  function get_mosaic_xgrid_size(xgrid_file)
    character(len=*), intent(in)          :: xgrid_file
    integer                               :: get_mosaic_xgrid_size
    character(len=len_trim(xgrid_file)+1) :: xfile    
    integer                               :: read_mosaic_xgrid_size
    integer                               :: strlen

    !---- transfer to C-stype string
    strlen = len_trim(xgrid_file)
    xfile(1:strlen) = xgrid_file(1:strlen)
    strlen = strlen+1
    xfile(strlen:strlen) = CHAR(0)

    get_mosaic_xgrid_size = read_mosaic_xgrid_size(xfile)

    return   

  end function get_mosaic_xgrid_size
! </FUNCTION>
!#######################################################################
! <SUBROUTINE NAME="get_mosaic_xgrid">
!   <OVERVIEW>
!     get exchange grid information from mosaic xgrid file.
!   </OVERVIEW>
!   <DESCRIPTION>
!     get exchange grid information from mosaic xgrid file.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_mosaic_xgrid(xgrid_file, nxgrid, i1, j1, i2, j2, area)
!   </TEMPLATE>
!   <IN NAME="xgrid_file" TYPE="character(len=*)">
!     The file that contains exchange grid information.
!   </IN>
!   <INOUT NAME="nxgrid" TYPE="integer">
!     number of exchange grid in xgrid_file
!   </INOUT>
!   <INOUT NAME="i1, j1" TYPE="integer, dimension(:)">
!     i and j-index in grid 1 of exchange grid.
!   </INOUT>
!   <INOUT NAME="i2, j2" TYPE="integer, dimension(:)">
!     i and j-index in grid 2 of exchange grid.
!   </INOUT>
!   <INOUT NAME="area" TYPE="real, dimension(:)">
!     area of the exchange grid. The area is scaled to represent unit earth area.
!   </INOUT>
  subroutine get_mosaic_xgrid(xgrid_file, i1, j1, i2, j2, area, di, dj)
    character(len=*), intent(in) :: xgrid_file
    integer,       intent(inout) :: i1(:), j1(:), i2(:), j2(:)
    real,          intent(inout) :: area(:)
    real, optional,intent(inout) :: di(:), dj(:)

    character(len=len_trim(xgrid_file)+1) :: xfile
    integer :: n, strlen, nxgrid

    !---- transfer to C-stype string
    strlen = len_trim(xgrid_file)
    xfile(1:strlen) = xgrid_file(1:strlen)
    strlen = strlen+1
    xfile(strlen:strlen) = CHAR(0)

    !--- order 2 xgrid will be implemented later 
    nxgrid = size(i1(:))

    if(PRESENT(di)) then
       if(.NOT. PRESENT(dj) ) call mpp_error(FATAL, "mosaic_mod: when di is present, dj should be present")
       call read_mosaic_xgrid_order2(xfile, i1, j1, i2, j2, area, di, dj)
    else
       call read_mosaic_xgrid_order1(xfile, i1, j1, i2, j2, area)
    end if

    ! in C, programming, the starting index is 0, so need add 1 to the index.
    do n = 1, nxgrid
       i1(n) = i1(n) + 1
       j1(n) = j1(n) + 1
       i2(n) = i2(n) + 1
       j2(n) = j2(n) + 1
    end do
  end subroutine get_mosaic_xgrid
! </SUBROUTINE>

  !###############################################################################
  ! <SUBROUTINE NAME="get_mosaic_ntiles">
  !   <OVERVIEW>
  !     get number of tiles in the mosaic_file.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     get number of tiles in the mosaic_file.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     ntiles = get_mosaic_ntiles( mosaic_file)
  !   </TEMPLATE>
  !   <IN NAME="mosaic_file" TYPE="character(len=*)">
  !     The file that contains mosaic information.
  !   </IN>
  function get_mosaic_ntiles(mosaic_file)
    character(len=*), intent(in) :: mosaic_file
    integer                      :: get_mosaic_ntiles 

    character(len=len_trim(mosaic_file)+1) :: mfile    
    integer                                :: strlen
    integer                                :: read_mosaic_ntiles

    !---- transfer to C-stype string
    strlen = len_trim(mosaic_file)
    mfile(1:strlen) = mosaic_file(1:strlen)
    strlen = strlen+1
    mfile(strlen:strlen) = CHAR(0)

    get_mosaic_ntiles = read_mosaic_ntiles(mfile)

  end function get_mosaic_ntiles
! </SUBROUTINE>

  !###############################################################################
  ! <SUBROUTINE NAME="get_mosaic_ncontacts">
  !   <OVERVIEW>
  !     get number of contacts in the mosaic_file.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     get number of contacts in the mosaic_file.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     ntiles = get_mosaic_ncontacts( mosaic_file)
  !   </TEMPLATE>
  !   <IN NAME="mosaic_file" TYPE="character(len=*)">
  !     The file that contains mosaic information.
  !   </IN>
  function get_mosaic_ncontacts( mosaic_file)
    character(len=*), intent(in) :: mosaic_file
    integer                      :: get_mosaic_ncontacts 

    character(len=len_trim(mosaic_file)+1) :: mfile    
    integer                                :: strlen
    integer                                :: read_mosaic_ncontacts

    !---- transfer to C-stype string
    strlen = len_trim(mosaic_file)
    mfile(1:strlen) = mosaic_file(1:strlen)
    strlen = strlen+1
    mfile(strlen:strlen) = CHAR(0)

    get_mosaic_ncontacts = read_mosaic_ncontacts(mfile)

  end function get_mosaic_ncontacts
! </SUBROUTINE>


  !###############################################################################
  ! <SUBROUTINE NAME="get_mosaic_grid_sizes">
  !   <OVERVIEW>
  !     get grid size of each tile from mosaic_file
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     get grid size of each tile from mosaic_file
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call get_mosaic_grid_sizes(mosaic_file, nx, ny)
  !   </TEMPLATE>
  !   <IN NAME="mosaic_file" TYPE="character(len=*)">
  !     The file that contains mosaic information.
  !   </IN>
  !   <INOUT NAME="nx" TYPE="integer, dimension(:)">
  !     List of grid size in x-direction of each tile.
  !   </INOUT>
  !   <INOUT NAME="ny" TYPE="integer, dimension(:)">
  !     List of grid size in y-direction of each tile.
  !   </INOUT>
  subroutine get_mosaic_grid_sizes( mosaic_file, nx, ny)
    character(len=*),         intent(in) :: mosaic_file
    integer, dimension(:), intent(inout) :: nx, ny

    character(len=len_trim(mosaic_file)+1) :: mfile    
    integer                                :: strlen

    !---- transfer to C-stype string
    strlen = len_trim(mosaic_file)
    mfile(1:strlen) = mosaic_file(1:strlen)
    strlen = strlen+1
    mfile(strlen:strlen) = CHAR(0)

    call read_mosaic_grid_sizes(mfile, nx, ny)

  end subroutine get_mosaic_grid_sizes
! </SUBROUTINE>

  !###############################################################################
  ! <SUBROUTINE NAME="get_mosaic_contact">
  !   <OVERVIEW>
  !     get contact information from mosaic_file
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     get contact information from mosaic_file
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call get_mosaic_contact(mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1,
  !                             istart2, iend2, jstart2, jend2)
  !   </TEMPLATE>
  !   <IN NAME="mosaic_file" TYPE="character(len=*)">
  !     The file that contains mosaic information.
  !   </IN>
  !   <INOUT NAME="tile1" TYPE="integer, dimension(:)">
  !     list tile number in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="tile1" TYPE="integer, dimension(:)">
  !     list tile number in tile 2 of each contact.
  !   </INOUT>
  !   <INOUT NAME="istart1" TYPE="integer, dimension(:)">
  !     list starting i-index in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="iend1" TYPE="integer, dimension(:)">
  !     list ending i-index in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="jstart1" TYPE="integer, dimension(:)">
  !     list starting j-index in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="jend1" TYPE="integer, dimension(:)">
  !     list ending j-index in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="istart2" TYPE="integer, dimension(:)">
  !     list starting i-index in tile 2 of each contact.
  !   </INOUT>
  !   <INOUT NAME="iend2" TYPE="integer, dimension(:)">
  !     list ending i-index in tile 2 of each contact.
  !   </INOUT>
  !   <INOUT NAME="jstart2" TYPE="integer, dimension(:)">
  !     list starting j-index in tile 2 of each contact.
  !   </INOUT>
  !   <INOUT NAME="jend2" TYPE="integer, dimension(:)">
  !     list ending j-index in tile 2 of each contact.
  !   </INOUT>
  subroutine get_mosaic_contact( mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, &
                                   istart2, iend2, jstart2, jend2)
    character(len=*),         intent(in) :: mosaic_file
    integer, dimension(:), intent(inout) :: tile1, tile2
    integer, dimension(:), intent(inout) :: istart1, iend1, jstart1, jend1
    integer, dimension(:), intent(inout) :: istart2, iend2, jstart2, jend2
    character(len=len_trim(mosaic_file)+1) :: mfile    
    integer                                :: strlen

    !---- transfer to C-stype string
    strlen = len_trim(mosaic_file)
    mfile(1:strlen) = mosaic_file(1:strlen)
    strlen = strlen+1
    mfile(strlen:strlen) = CHAR(0)

    call read_mosaic_contact(mfile, tile1, tile2, istart1, iend1, jstart1, jend1, &
                            istart2, iend2, jstart2, jend2)
    !--- transfer C-index to Fortran-index.
    istart1 = istart1 + 1
    iend1   = iend1   + 1
    jstart1 = jstart1 + 1
    jend1   = jend1   + 1
    istart2 = istart2 + 1
    iend2   = iend2   + 1
    jstart2 = jstart2 + 1
    jend2   = jend2   + 1

  end subroutine get_mosaic_contact
! </SUBROUTINE>

  !###############################################################################
  ! <SUBROUTINE NAME="calc_mosaic_grid_area">
  !   <OVERVIEW>
  !     calculate grid cell area.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     calculate the grid cell area. The purpose of this routine is to make 
  !     sure the consistency between model grid area and exchange grid area.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call calc_mosaic_grid_area(lon, lat, area)
  !   </TEMPLATE>
  !   <IN NAME="lon" TYPE="real, dimension(:,:)">
  !     geographical longitude of grid cell vertices.
  !   </IN>
  !   <IN NAME="lat" TYPE="real, dimension(:,:)">
  !     geographical latitude of grid cell vertices.
  !   </IN>
  !   <INOUT NAME="area" TYPE="real, dimension(:,:)">
  !     grid cell area.
  !   </INOUT>
  subroutine calc_mosaic_grid_area(lon, lat, area)
     real, dimension(:,:), intent(in)    :: lon
     real, dimension(:,:), intent(in)    :: lat
     real, dimension(:,:), intent(inout) :: area
     integer                             :: nlon, nlat

     nlon = size(area,1)
     nlat = size(area,2)
     ! make sure size of lon, lat and area are consitency
     if( size(lon,1) .NE. nlon+1 .OR. size(lat,1) .NE. nlon+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,1) and size(lat,1) should equal to size(area,1)+1")
     if( size(lon,2) .NE. nlat+1 .OR. size(lat,2) .NE. nlat+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,2) and size(lat,2) should equal to size(area,2)+1")

     call get_grid_area( nlon, nlat, lon, lat, area)

  end subroutine calc_mosaic_grid_area
  ! </SUBROUTINE>

end module mosaic_mod


#ifdef TEST_MOSAIC
program test_mosaic

use mosaic_mod, only : get_mosaic_ntiles, get_mosaic_ncontacts
use mosaic_mod, only : get_mosaic_grid_sizes, get_mosaic_contact

implicit none

integer              :: ntiles, ncontacts, n
integer, allocatable :: tile1(:), tile2(:), nx(:), ny(:)
integer, allocatable :: istart1(:), iend1(:), jstart1(:), jend1(:)
integer, allocatable :: istart2(:), iend2(:), jstart2(:), jend2(:)
character(len=128)   :: mosaic_file = "INPUT/mosaic.nc"

ntiles = get_mosaic_ntiles(mosaic_file)
ncontacts = get_mosaic_ncontacts(mosaic_file)
allocate(nx(ntiles), ny(ntiles))
allocate(tile1(ncontacts), tile2(ncontacts) )
allocate(istart1(ncontacts), iend1(ncontacts), jstart1(ncontacts), jend1(ncontacts) )
allocate(istart2(ncontacts), iend2(ncontacts), jstart2(ncontacts), jend2(ncontacts) )

call get_mosaic_grid_sizes(mosaic_file, nx, ny )
call get_mosaic_contact(mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2)

! print out information

print '(a,i3,a,a)', "****** There is ", ntiles, " tiles in ", trim(mosaic_file)
do n = 1, ntiles
   print '(a,i3,a,i3,a,i3)', " tile = ", n, ", nx = ", nx(n), ", ny = ", ny(n)
end do

print '(a,i3,a,a)', "****** There is ", ncontacts, " contacts in ", trim(mosaic_file)
do n = 1, ncontacts
   print '(a,i3,a,i3,a,i3,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)', &
           "contact=", n, ": tile1=", tile1(n), " tile2=", tile2(n),   &
           " is1=", istart1(n), " ie1=", iend1(n),                   &
           " js1=", jstart1(n), " je1=", jend1(n),                   &
           " is2=", istart2(n), " ie2=", iend2(n),                   &
           " js2=", jstart2(n), " je2=", jend2(n)
end do

deallocate(tile1, tile2, nx, ny)
deallocate(istart1, iend1, jstart1, jend1)
deallocate(istart2, iend2, jstart2, jend2)


end program test_mosaic
#endif
