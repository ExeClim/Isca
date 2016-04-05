module grid_mod

use constants_mod, only : PI, radius
use fms_mod, only : uppercase, lowercase, field_exist, field_size, read_data, &
     error_mesg, string, FATAL
use mosaic_mod, only : get_mosaic_ntiles, get_mosaic_xgrid_size, get_mosaic_grid_sizes, &
     get_mosaic_xgrid, calc_mosaic_grid_area

! the following two use statement are only needed for define_cube_mosaic
use mpp_domains_mod, only : domain2d, mpp_define_mosaic
use mosaic_mod, only : get_mosaic_ncontacts, get_mosaic_contact

implicit none;private

! ==== public interfaces =====================================================
! grid dimension inquiry subroutines
public :: get_grid_ntiles ! returns number of tiles
public :: get_grid_size   ! returns horizontal sizes of the grid
! grid geometry inquiry subroutines
public :: get_grid_cell_centers 
public :: get_grid_cell_vertices
! grid area inquiry subroutines
public :: get_grid_cell_area
public :: get_grid_comp_area
! decompose cubed sphere domains -- probably does not belong here, but it should 
! be in some place available for component models
public :: define_cube_mosaic
! ==== end of public interfaces ==============================================

interface get_grid_size
   module procedure get_grid_size_for_all_tiles
   module procedure get_grid_size_for_one_tile
end interface

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'grid_mod', &
     version     = '$Id: grid.F90,v 18.0 2010/03/02 23:56:19 fms Exp $', &
     tagname     = '$Name: testing $'

character(len=*), parameter :: &
     grid_dir  = 'INPUT/',     &      ! root directory for all grid files
     grid_file = 'INPUT/grid_spec.nc' ! name of the grid spec file

integer, parameter :: &
     MAX_NAME = 256,  & ! max length of the variable names
     MAX_FILE = 1024, & ! max length of the file names
     VERSION_0 = 0,   &
     VERSION_1 = 1,   &
     VERSION_2 = 2

! ==== module variables ======================================================
integer :: grid_version = -1

contains 

function get_grid_version()
  integer :: get_grid_version
  if(grid_version<0) then
    if(field_exist(grid_file, 'geolon_t')) then
       grid_version = VERSION_0 
    else if(field_exist(grid_file, 'x_T')) then
       grid_version = VERSION_1
    else if(field_exist(grid_file, 'ocn_mosaic_file') ) then
       grid_version = VERSION_2
    else
       call error_mesg(module_name,&
            'Can''t determine the version of the grid spec: none of "x_T", "geolon_t", or "ocn_mosaic_file" exist in file "'//trim(grid_file)//'"', &
            FATAL )
    endif
  endif
  get_grid_version = grid_version
end function get_grid_version


! ============================================================================
! returns number of tiles for a given component
! ============================================================================
subroutine get_grid_ntiles(component,ntiles)
  character(len=*)     :: component
  integer, intent(out) :: ntiles

  ! local vars
  character(len=MAX_FILE) :: component_mosaic

  select case (get_grid_version())
  case(VERSION_0,VERSION_1)
     ntiles = 1
  case(VERSION_2)
     call read_data(grid_file,trim(lowercase(component))//'_mosaic_file',component_mosaic)
     ntiles = get_mosaic_ntiles(grid_dir//trim(component_mosaic))
  end select
end subroutine get_grid_ntiles


! ============================================================================
! returns size of the grid for each of the tiles
! ============================================================================
subroutine get_grid_size_for_all_tiles(component,nx,ny)
  character(len=*)     :: component
  integer, intent(inout) :: nx(:),ny(:)

  ! local vars
  integer :: siz(4) ! for the size of external fields
  character(len=MAX_NAME) :: varname1, varname2
  character(len=MAX_FILE) :: component_mosaic
  
  varname1 = 'AREA_'//trim(uppercase(component))
  varname2 = trim(lowercase(component))//'_mosaic_file'

  select case (get_grid_version())
  case(VERSION_0,VERSION_1)
     call field_size(grid_file, varname1, siz)
     nx(1) = siz(1); ny(1)=siz(2)
  case(VERSION_2) ! mosaic file
     call read_data(grid_file,varname2, component_mosaic)
     call get_mosaic_grid_sizes(grid_dir//trim(component_mosaic),nx,ny)
  end select
end subroutine get_grid_size_for_all_tiles


! ============================================================================
! returns size of the grid for one of the tiles
! ============================================================================
subroutine get_grid_size_for_one_tile(component,tile,nx,ny)
  character(len=*)       :: component
  integer, intent(in)    :: tile
  integer, intent(inout) :: nx,ny
  
  ! local vars
  integer, allocatable :: nnx(:), nny(:)
  integer :: ntiles

  call get_grid_ntiles(component, ntiles)
  if(tile>0.and.tile<=ntiles) then
     allocate(nnx(ntiles),nny(ntiles))
     call get_grid_size_for_all_tiles(component,nnx,nny)
     nx = nnx(tile); ny = nny(tile)
     deallocate(nnx,nny)
  else
     call error_mesg('get_grid_size',&
          'requested tile index '//trim(string(tile))//' is out of bounds (1:'//trim(string(ntiles))//')',&
          FATAL)
  endif
end subroutine get_grid_size_for_one_tile

! ============================================================================
! return grid cell area for the specified model component and tile
! ============================================================================
subroutine get_grid_cell_area(component, tile, cellarea)
  character(len=*), intent(in)    :: component
  integer         , intent(in)    :: tile
  real            , intent(inout) :: cellarea(:,:)

  ! local vars
  integer :: nlon, nlat
  real, allocatable :: glonb(:,:), glatb(:,:)

  select case(get_grid_version())
  case(VERSION_0,VERSION_1)
     select case(trim(component))
     case('LND')
        call read_data(grid_file, 'AREA_LND_CELL', cellarea, no_domain=.TRUE.)
     case('ATM','OCN')
        call read_data(grid_file, 'AREA_'//trim(uppercase(component)),cellarea,no_domain=.TRUE.)
     case default
        call error_mesg(module_name,&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select
     ! convert area to m2
     cellarea = cellarea*4*PI*radius**2
  case(VERSION_2)
     call get_grid_size(component,tile,nlon,nlat)
     allocate(glonb(nlon+1,nlat+1),glatb(nlon+1,nlat+1))
     call get_grid_cell_vertices(component, tile, glonb, glatb)
     call calc_mosaic_grid_area(glonb*pi/180.0, glatb*pi/180.0, cellarea)
     deallocate(glonb,glatb)
  end select

end subroutine get_grid_cell_area


! ============================================================================
! get the area of the component per grid cell
! ============================================================================
subroutine get_grid_comp_area(component,tile,area)
  character(len=*) :: component
  integer, intent(in) :: tile
  real, intent(inout) :: area(:,:)
  ! locala vars
  integer :: n_xgrid_files ! number of exchange grid files in the mosaic
  integer :: siz(4), nxgrid
  integer :: i,j,m,n
  integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
  real, allocatable :: xgrid_area(:)
  real, allocatable :: rmask(:,:)
  character(len=MAX_NAME) :: &
     xgrid_name, & ! name of the variable holding xgrid names
     tile_name,  & ! name of the tile
     xgrid_file    ! name of the current xgrid file
  character(len=MAX_NAME) :: varname1, varname2

  select case (get_grid_version())
  case(VERSION_0,VERSION_1)
     select case(component)
     case('ATM')
        call read_data(grid_file,'AREA_ATM',area, no_domain=.TRUE.)
     case('OCN')
        allocate(rmask(size(area,1),size(area,2)))
        call read_data(grid_file,'AREA_OCN',area, no_domain=.TRUE.)
        call read_data(grid_file,'wet',     rmask,no_domain=.TRUE.)
        area = area*rmask
        deallocate(rmask)
     case('LND')
        call read_data(grid_file,'AREA_LND',area,no_domain=.TRUE.)
     case default
        call error_mesg(module_name,&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select
  case(VERSION_2) ! mosaic gridspec
     select case (component)
     case ('ATM')
        ! just read the grid cell area and return
        call get_grid_cell_area(component,tile,area)
        return
     case ('LND')
        xgrid_name = 'aXl_file'
        tile_name  = 'land_mosaic_tile'//char(tile+ichar('0'))
     case ('OCN')
        xgrid_name = 'aXo_file'
        tile_name  = 'ocean_mosaic_tile'//char(tile+ichar('0'))
     case default
        call error_mesg(module_name,&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select

     area(:,:) = 0
     ! get the number of the exchange-grid files
     if(field_exist(grid_file,xgrid_name)) then

        call field_size(grid_file,xgrid_name,siz)
        ! loop through all exchange grid files

        n_xgrid_files = siz(2)
        do n = 1, n_xgrid_files
           ! get the name of the current exchange grid file
           call read_data(grid_file,xgrid_name,xgrid_file,level=n)
           ! skip the rest of the loop if the name of the current tile isn't found 
           ! in the file name, but check this only if there is more than 1 tile
           if(n_xgrid_files>1) then
              if(index(xgrid_file,trim(tile_name))==0) cycle
           endif
           ! finally read the exchange grid
           nxgrid = get_mosaic_xgrid_size(grid_dir//xgrid_file)
           allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), xgrid_area(nxgrid))
           call get_mosaic_xgrid(grid_dir//xgrid_file, i1, j1, i2, j2, xgrid_area)
           ! and sum the exchange grid areas
           do m = 1, nxgrid
              i = i2(m); j = j2(m)
              area(i,j) = area(i,j) + xgrid_area(m)
           end do
           deallocate(i1, j1, i2, j2, xgrid_area)
        enddo
     endif
  end select ! version
  ! convert area to m2
  area = area*4*PI*radius**2
end subroutine


! ============================================================================
! returns cell vertices for the specified model component and mosaic tile number
! ============================================================================
subroutine get_grid_cell_vertices(component, tile, glonb, glatb)
  character(len=*) :: component
  integer, intent(in) :: tile
  real, intent(inout) :: glonb(:,:),glatb(:,:)
  ! local vars
  character(len=MAX_FILE) :: filename1, filename2
  integer :: nlon, nlat
  integer :: i,j
  real, allocatable :: buffer(:), tmpx(:,:), tmpy(:,:), x_vert_t(:,:,:), y_vert_t(:,:,:)

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)

  select case(get_grid_version())
  case(VERSION_0)
     select case(component)
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)+1))
        ! read coordinates of grid cell vertices
        call read_data(grid_file, 'xb'//lowercase(component(1:1)), buffer(1:nlon+1), no_domain=.true.)
        do j = 1, nlat+1
           do i = 1, nlon+1
              glonb(i,j) = buffer(i)
           enddo
        enddo
        call read_data(grid_file, 'yb'//lowercase(component(1:1)), buffer(1:nlat+1), no_domain=.true.)
        do j = 1, nlat+1
           do i = 1, nlon+1
              glatb(i,j) = buffer(j)
           enddo
        enddo
        deallocate(buffer)
     case('OCN')
        call read_data(grid_file, 'geolon_vert_t', glonb, no_domain=.TRUE. )
        call read_data(grid_file, 'geolat_vert_t', glatb, no_domain=.TRUE. )
     case default
        call error_mesg(module_name,&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select
  case(VERSION_1)
     select case(component)
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)+1))
        ! read coordinates of grid cell vertices
        call read_data(grid_file, 'xb'//lowercase(component(1:1)), buffer(1:nlon+1), no_domain=.true.)
        do j = 1, nlat+1
           do i = 1, nlon+1
              glonb(i,j) = buffer(i)
           enddo
        enddo
        call read_data(grid_file, 'yb'//lowercase(component(1:1)), buffer(1:nlat+1), no_domain=.true.)
        do j = 1, nlat+1
           do i = 1, nlon+1
              glatb(i,j) = buffer(j)
           enddo
        enddo
        deallocate(buffer)
     case('OCN')
        allocate (x_vert_t(nlon,nlat,4), y_vert_t(nlon,nlat,4) ) 
        call read_data(grid_file, 'x_vert_T', x_vert_t, no_domain=.TRUE.)
        call read_data(grid_file, 'y_vert_T', y_vert_t, no_domain=.TRUE. )
        glonb(1:nlon,1:nlat) = x_vert_t(1:nlon,1:nlat,1)
        glonb(nlon+1,1:nlat) = x_vert_t(nlon,1:nlat,2)
        glonb(1:nlon,nlat+1) = x_vert_t(1:nlon,nlat,4)
        glonb(nlon+1,nlat+1) = x_vert_t(nlon,nlat,3)
        glatb(1:nlon,1:nlat) = y_vert_t(1:nlon,1:nlat,1)
        glatb(nlon+1,1:nlat) = y_vert_t(nlon,1:nlat,2)
        glatb(1:nlon,nlat+1) = y_vert_t(1:nlon,nlat,4)
        glatb(nlon+1,nlat+1) = y_vert_t(nlon,nlat,3)
        deallocate(x_vert_t, y_vert_t)
     case default
        call error_mesg(module_name,&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select
  case(VERSION_2)
     ! get the name of the mosaic file for the component
     call read_data(grid_file, trim(lowercase(component))//'_mosaic_file', filename1)
     filename1=grid_dir//trim(filename1)
     ! get the name of the grid file for the component and tile
     call read_data(filename1, 'gridfiles', filename2, level=tile)
     filename2 = grid_dir//trim(filename2)

     allocate(tmpx(2*nlon+1,2*nlat+1),tmpy(2*nlon+1,2*nlat+1))
     call read_data(filename2, 'x', tmpx, no_domain=.TRUE.)
     call read_data(filename2, 'y', tmpy, no_domain=.TRUE.)
     do j = 1, nlat+1
        do i = 1, nlon+1
           glonb(i,j) = tmpx(2*i-1,2*j-1)
           glatb(i,j) = tmpy(2*i-1,2*j-1)
        end do
     end do
     deallocate(tmpx,tmpy)
  end select

end subroutine get_grid_cell_vertices


! ============================================================================
! returns grid cell centers for specified model component and mosaic tile number
! ============================================================================
subroutine get_grid_cell_centers(component, tile, glon, glat)
  character(len=*), intent(in) :: component
  integer, intent(in) :: tile
  real, intent(inout) :: glon(:,:),glat(:,:)
  ! local vars
  character(len=MAX_NAME) :: varname
  character(len=MAX_FILE) :: filename1, filename2
  integer :: nlon, nlat
  integer :: i,j
  real, allocatable :: buffer(:),tmpx(:,:),tmpy(:,:)

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)

  select case(get_grid_version())
  case(VERSION_0)
     select case (trim(component))
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)))
        ! read coordinates of grid cell vertices
        call read_data(grid_file, 'xt'//lowercase(component(1:1)), buffer(1:nlon), no_domain=.true.)
        do j = 1, nlat
           do i = 1, nlon
              glon(i,j) = buffer(i)
           enddo
        enddo
        call read_data(grid_file, 'yt'//lowercase(component(1:1)), buffer(1:nlat), no_domain=.true.)
        do j = 1, nlat
           do i = 1, nlon
              glat(i,j) = buffer(j)
           enddo
        enddo
        deallocate(buffer)
     case('OCN')
        call read_data(grid_file, 'geolon_t', glon, no_domain=.TRUE. )
        call read_data(grid_file, 'geolat_t', glat, no_domain=.TRUE. )
     case default
        ! error processing here
     end select
  case(VERSION_1)
     select case(trim(component))
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)))
        ! read coordinates of grid cell vertices
        call read_data(grid_file, 'xt'//lowercase(component(1:1)), buffer(1:nlon), no_domain=.true.)
        do j = 1, nlat
           do i = 1, nlon
              glon(i,j) = buffer(i)
           enddo
        enddo
        call read_data(grid_file, 'yt'//lowercase(component(1:1)), buffer(1:nlat), no_domain=.true.)
        do j = 1, nlat
           do i = 1, nlon
              glat(i,j) = buffer(j)
           enddo
        enddo
        deallocate(buffer)
     case('OCN')
        call read_data(grid_file, 'x_T', glon, no_domain=.TRUE. )
        call read_data(grid_file, 'y_T', glat, no_domain=.TRUE. )
     case default
        call error_mesg(module_name,&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select
  case(VERSION_2) ! mosaic grid file
     ! get the name of the mosaic file for the component
     call read_data(grid_file, trim(lowercase(component))//'_mosaic_file', filename1)
     filename1=grid_dir//trim(filename1)
     ! get the name of the grid file for the component and tile
     call read_data(filename1, 'gridfiles', filename2, level=tile)
     filename2 = grid_dir//trim(filename2)

     allocate(tmpx(2*nlon+1,2*nlat+1),tmpy(2*nlon+1,2*nlat+1))
     call read_data(filename2, 'x', tmpx, no_domain=.TRUE.)
     call read_data(filename2, 'y', tmpy, no_domain=.TRUE.)
     do j = 1, nlat
        do i = 1, nlon
           glon(i,j) = tmpx(2*i,2*j)
           glat(i,j) = tmpy(2*i,2*j)
        end do
     end do
     deallocate(tmpx,tmpy)
  end select
end subroutine get_grid_cell_centers


! ============================================================================
! given a model component, a layout, and (optionally) a halo size, returns a 
! domain for current processor
! ============================================================================
! this subroutine probably does not belong in the grid_mod 
subroutine define_cube_mosaic ( component, domain, layout, halo )
  character(len=*) , intent(in)    :: component
  type(domain2d)   , intent(inout) :: domain
  integer          , intent(in)    :: layout(2)
  integer, optional, intent(in)    :: halo 

  ! ---- local constants
  
  ! ---- local vars
  character(len=MAX_NAME) :: varname
  character(len=MAX_FILE) :: mosaic_file
  integer :: ntiles     ! number of tiles
  integer :: ncontacts  ! number of contacts between mosaic tiles
  integer :: n
  integer :: ng         ! halo size
  integer, allocatable :: nlon(:), nlat(:), global_indices(:,:)
  integer, allocatable :: pe_start(:), pe_end(:), layout_2d(:,:)
  integer, allocatable :: tile1(:),tile2(:)
  integer, allocatable :: is1(:),ie1(:),js1(:),je1(:)
  integer, allocatable :: is2(:),ie2(:),js2(:),je2(:)

  call get_grid_ntiles(component,ntiles)
  allocate(nlon(ntiles), nlat(ntiles))
  allocate(global_indices(4,ntiles))
  allocate(pe_start(ntiles),pe_end(ntiles))
  allocate(layout_2d(2,ntiles))
  call get_grid_size(component,nlon,nlat)

  do n = 1, ntiles
     global_indices(:,n) = (/ 1, nlon(n), 1, nlat(n) /)
     layout_2d     (:,n) = layout
     pe_start        (n) = (n-1)*layout(1)*layout(2)
     pe_end          (n) =     n*layout(1)*layout(2) - 1
  enddo

  varname=trim(lowercase(component))//'_mosaic_file'
  call read_data(grid_file,varname,mosaic_file)
  mosaic_file = grid_dir//mosaic_file

  ! get the contact information from mosaic file
  ncontacts = get_mosaic_ncontacts(mosaic_file)
  allocate(tile1(ncontacts),tile2(ncontacts))
  allocate(is1(ncontacts),ie1(ncontacts),js1(ncontacts),je1(ncontacts))
  allocate(is2(ncontacts),ie2(ncontacts),js2(ncontacts),je2(ncontacts))
  call get_mosaic_contact(mosaic_file, tile1, tile2, &
       is1, ie1, js1, je1, is2, ie2, js2, je2)

  ng = 0
  if(present(halo)) ng = halo
  ! create the domain2d variable
  call mpp_define_mosaic ( global_indices, layout_2d, domain, &
       ntiles, ncontacts, tile1, tile2,                  &
       is1, ie1, js1, je1, &
       is2, ie2, js2, je2, &
       pe_start=pe_start, pe_end=pe_end, symmetry=.true.,  &
       shalo = ng, nhalo = ng, whalo = ng, ehalo = ng,     &
       name = trim(component)//'Cubic-Sphere Grid' )

  deallocate(nlon,nlat,global_indices,pe_start,pe_end,layout_2d)
  deallocate(tile1,tile2)
  deallocate(is1,ie1,js1,je1)
  deallocate(is2,ie2,js2,je2)

end subroutine define_cube_mosaic

end module grid_mod
