!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module interpolator_mod
!
! Purpose: Module to interpolate climatology data to model grid.
!
! author: William Cooke GFDL.Climate.Model.Info@noaa.gov
!
#include <fms_platform.h>

use mpp_mod,           only : mpp_error, &
                              FATAL,     &
                              mpp_pe,    &
                              mpp_init,  &
                              mpp_exit,  &
                              mpp_npes,  &
                              WARNING,   &
                              NOTE
use mpp_io_mod,        only : mpp_open,          &
                              mpp_close,         &
                              mpp_get_times,     &
                              mpp_get_atts,      &
                              mpp_get_info,      &
                              mpp_read,          &
                              mpp_get_axes,      &
                              mpp_get_axis_data, &
                              mpp_get_fields,    &
                              fieldtype,         &
                              atttype,           &
                              axistype,          &
                              MPP_RDONLY,        &
                              MPP_NETCDF,        &
                              MPP_MULTI,         &
                              MPP_APPEND,        &
                              MPP_SINGLE
use mpp_domains_mod,   only : mpp_domains_init,      &
                              mpp_update_domains,    &
                              mpp_define_domains,    &
                              mpp_global_field,      &
                              domain2d,              &
                              mpp_define_layout,     &
                              mpp_get_compute_domain
use diag_manager_mod,  only : diag_manager_init, get_base_time, &
                              register_diag_field, send_data, &
                              diag_axis_init
use fms_mod,           only : lowercase, write_version_number, &
                              fms_init, &
                              file_exist, mpp_root_pe, stdlog
use horiz_interp_mod,  only : horiz_interp_type, &
                              horiz_interp_new,  &
                              horiz_interp_init, &
                              assignment(=), &
                              horiz_interp,      &
                              horiz_interp_del
use time_manager_mod,  only : time_type,   &
                              set_time,    &
                              set_date,    &
                              get_date,    &
                              get_calendar_type, &
                              JULIAN, NOLEAP, &
                              THIRTY_DAY_MONTHS, & !mj
                              get_date_julian, set_date_no_leap, &
                              set_date_julian, get_date_no_leap, &
                              print_date, &
                              operator(+), &
                              operator(-), &
                              operator(*), &
                              operator(>), &
                              operator(<), &
                              assignment(=), &
                              decrement_time
use time_interp_mod,   only : time_interp, YEAR
use constants_mod,     only : grav, PI

implicit none
private 

public interpolator_init, &
       interpolator,      &
       interpolate_type_eq, &
       obtain_interpolator_time_slices, &
       unset_interpolator_time_flag, &
       interpolator_end,  &
       init_clim_diag,    &
       query_interpolator,&
       read_data

interface interpolator
   module procedure interpolator_4D
   module procedure interpolator_3D
   module procedure interpolator_2D
   module procedure interpolator_4D_no_time_axis
   module procedure interpolator_3D_no_time_axis
   module procedure interpolator_2D_no_time_axis
end interface 

interface assignment(=)
   module procedure interpolate_type_eq
end interface

interface interp_weighted_scalar
   module procedure interp_weighted_scalar_1D
   module procedure interp_weighted_scalar_2D
end interface interp_weighted_scalar
character(len=128) :: version = &
'$Id: interpolator.F90,v 19.0.8.2.2.1 2012/08/23 18:29:52 Zhi.Liang Exp $'
character(len=128) :: tagname = '$Name: siena_201211 $'
logical            :: module_is_initialized = .false.
logical            :: clim_diag_initialized = .false.

type, public  :: interpolate_type
private
!Redundant data between fields
!All climatology data
real, pointer            :: lat(:) =>NULL()
real, pointer            :: lon(:) =>NULL()
real, pointer            :: latb(:) =>NULL()
real, pointer            :: lonb(:) =>NULL()
real, pointer            :: levs(:) =>NULL()
real, pointer            :: halflevs(:) =>NULL()
type(horiz_interp_type)  :: interph
type(time_type), pointer :: time_slice(:) =>NULL() ! An array of the times within the climatology.
integer                  :: unit          ! Unit number on which file is being read.
character(len=64)        :: file_name     ! Climatology filename
integer                  :: TIME_FLAG     ! Linear or seaonal interpolation?
integer                  :: level_type    ! Pressure or Sigma level
integer                  :: is,ie,js,je
integer                  :: vertical_indices ! direction of vertical 
                                              ! data axis
logical                  :: climatological_year ! Is data for year = 0000?

!Field specific data  for nfields
type(fieldtype),   pointer :: field_type(:) =>NULL()   ! NetCDF field type
character(len=64), pointer :: field_name(:) =>NULL()   ! name of this field
integer,           pointer :: time_init(:,:) =>NULL()  ! second index is the number of time_slices being kept. 2 or ntime.
integer,           pointer :: mr(:) =>NULL()           ! Flag for conversion of climatology to mixing ratio. 
integer,           pointer :: out_of_bounds(:) =>NULL()! Flag for when surface pressure is out of bounds.
!++lwh
integer,           pointer :: vert_interp(:) =>NULL()  ! Flag for type of vertical interpolation.
!--lwh
real,              pointer :: data(:,:,:,:,:) =>NULL() ! (nlatmod,nlonmod,nlevclim,size(time_init,2),nfields)

real,              pointer :: pmon_pyear(:,:,:,:) =>NULL()
real,              pointer :: pmon_nyear(:,:,:,:) =>NULL()
real,              pointer :: nmon_nyear(:,:,:,:) =>NULL()
real,              pointer :: nmon_pyear(:,:,:,:) =>NULL()
!integer                    :: indexm, indexp, climatology
integer,dimension(:),  pointer :: indexm =>NULL() 
integer,dimension(:),  pointer :: indexp =>NULL()
integer,dimension(:),  pointer :: climatology =>NULL() 

type(time_type), pointer :: clim_times(:,:) => NULL()
logical :: separate_time_vary_calc
real :: tweight, tweight1, tweight2, tweight3
integer :: itaum, itaup
end type interpolate_type


integer :: ndim, nvar,natt,ntime
integer :: nlat,nlatb,nlon,nlonb,nlev,nlevh
integer ::          len, ntime_in, num_fields
type(axistype), allocatable :: axes(:)
type(axistype),save          :: time_axis
type(fieldtype), allocatable :: varfields(:)

! pletzer real, allocatable :: time_in(:)
! sjs real, allocatable :: climdata(:,:,:), climdata2(:,:,:)

character(len=32) :: name, units       
integer           :: sense

integer, parameter :: max_diag_fields = 30

! flags to indicate direction of vertical axis in  data file
integer, parameter :: INCREASING_DOWNWARD = 1, INCREASING_UPWARD = -1
!++lwh
! Flags to indicate whether the time interpolation should be linear or some other scheme for seasonal data.
! NOTIME indicates that data file has no time axis.
integer, parameter :: LINEAR = 1, SEASONAL = 2, BILINEAR = 3, NOTIME = 4

! Flags to indicate where climatology pressure levels are pressure or sigma levels
integer, parameter :: PRESSURE = 1, SIGMA = 2 

! Flags to indicate whether the climatology units are mixing ratio (kg/kg) or column integral (kg/m2).
! Vertical interpolation scheme requires mixing ratio at this time.
integer, parameter :: NO_CONV = 1, KG_M2 = 2 

! Flags to indicate what to do when the model surface pressure exceeds the  climatology surface pressure level.
integer, parameter, public :: CONSTANT = 1, ZERO = 2 

! Flags to indicate the type of vertical interpolation
integer, parameter, public :: INTERP_WEIGHTED_P = 10, INTERP_LINEAR_P = 20, INTERP_LOG_P = 30
!--lwh

integer :: num_clim_diag = 0
character(len=64) :: climo_diag_name(max_diag_fields)
integer :: climo_diag_id(max_diag_fields), hinterp_id(max_diag_fields)
real ::  missing_value = -1.e10
! sjs integer :: itaum, itaup

#ifdef NO_QUAD_PRECISION
! 64-bit precision (kind=8)
 integer, parameter:: f_p = selected_real_kind(15)
#else
! Higher precision (kind=16) for grid geometrical factors:
 integer, parameter:: f_p = selected_real_kind(20)
#endif

logical :: read_all_on_init = .false.
integer :: verbose = 0
character(len=64) :: interp_method = "conserve_latlon"  

namelist /interpolator_nml/    &
                             read_all_on_init, verbose, interp_method

contains

!#####################################################################

subroutine interpolate_type_eq (Out, In)

type(interpolate_type), intent(in) :: In
type(interpolate_type), intent(inout) :: Out


     if (associated(In%lat))      Out%lat      =>  In%lat
     if (associated(In%lon))      Out%lon      =>  In%lon
     if (associated(In%latb))     Out%latb     =>  In%latb
     if (associated(In%lonb))     Out%lonb     =>  In%lonb
     if (associated(In%levs))     Out%levs     =>  In%levs
     if (associated(In%halflevs)) Out%halflevs =>  In%halflevs

     Out%interph = In%interph
     if (associated(In%time_slice)) Out%time_slice =>  In%time_slice
     Out%unit = In%unit
     Out%file_name = In%file_name
     Out%time_flag = In%time_flag
     Out%level_type = In%level_type
     Out%is = In%is
     Out%ie = In%ie
     Out%js = In%js
     Out%je = In%je
     Out%vertical_indices = In%vertical_indices
     Out%climatological_year = In%climatological_year
     Out%field_type => In%field_type
     if (associated(In%field_name   )) Out%field_name    =>  In%field_name
     if (associated(In%time_init    )) Out%time_init     =>  In%time_init 
     if (associated(In%mr           )) Out%mr            =>  In%mr         
     if (associated(In%out_of_bounds)) Out%out_of_bounds =>  In%out_of_bounds
     if (associated(In%vert_interp  )) Out%vert_interp   =>  In%vert_interp  
     if (associated(In%data         )) Out%data          =>  In%data  
     if (associated(In%pmon_pyear   )) Out%pmon_pyear    =>  In%pmon_pyear
     if (associated(In%pmon_nyear   )) Out%pmon_nyear    =>  In%pmon_nyear
     if (associated(In%nmon_nyear   )) Out%nmon_nyear    =>  In%nmon_nyear
     if (associated(In%nmon_pyear   )) Out%nmon_pyear    =>  In%nmon_pyear
     if (associated(In%indexm       )) Out%indexm        =>  In%indexm    
     if (associated(In%indexp       )) Out%indexp        =>  In%indexp    
     if (associated(In%climatology  )) Out%climatology   =>  In%climatology
     if (associated(In%clim_times   )) Out%clim_times    =>  In%clim_times
      Out%separate_time_vary_calc = In%separate_time_vary_calc
      Out%tweight = In%tweight
      Out%tweight1 = In%tweight1
      Out%tweight2 = In%tweight2
      Out%tweight3 = In%tweight3
      Out%itaum = In%itaum
      Out%itaup = In%itaup



end subroutine interpolate_type_eq 



 
!#######################################################################
!
subroutine interpolator_init( clim_type, file_name, lonb_mod, latb_mod, &
                              data_names, data_out_of_bounds,           &
                              vert_interp, clim_units, single_year_file)
type(interpolate_type), intent(inout) :: clim_type
character(len=*), intent(in)            :: file_name
real            , intent(in)            :: lonb_mod(:,:), latb_mod(:,:)
character(len=*), intent(in) , optional :: data_names(:)
!++lwh
integer         , intent(in)            :: data_out_of_bounds(:) 
integer         , intent(in), optional  :: vert_interp(:) 
!--lwh
character(len=*), intent(out), optional :: clim_units(:)
logical,          intent(out), optional :: single_year_file
!
! INTENT IN
!  file_name  :: Climatology filename
!  lonb_mod   :: The corners of the model grid-box longitudes.
!  latb_mod   :: The corners of the model grid_box latitudes.
!  data_names :: A list of the names of components within the climatology file which you wish to read.
!  data_out_of_bounds :: A list of the flags that are to be used in determining what to do if the pressure levels in the model
!                        go out of bounds from those of the climatology.
!  vert_interp:: Flag to determine type of vertical interpolation
!
! INTENT OUT
!  clim_type  :: An interpolate type containing the necessary file and field data to be passed to the interpolator routine.
!  clim_units :: A list of the units for the components listed in data_names.
!

integer                      :: unit
character(len=64)            :: src_file
!++lwh
real                         :: dlat, dlon
!--lwh
type(time_type)              :: base_time
logical                      :: NAME_PRESENT
real                         :: dtr,tpi
integer                      :: fileday, filemon, fileyr, filehr, filemin,filesec, m,m1
character(len= 20)           :: fileunits
real, dimension(:), allocatable  :: alpha
integer   :: j, i
logical :: non_monthly
character(len=24) :: file_calendar
integer :: model_calendar
integer :: yr, mo, dy, hr, mn, sc
integer :: n
type(time_type) :: Julian_time, Noleap_time
real, allocatable :: time_in(:)
real, allocatable, save :: agrid_mod(:,:,:)
integer :: nx, ny
integer :: io, ierr

if (.not. module_is_initialized) then
  call fms_init
  call diag_manager_init
  call horiz_interp_init
endif

clim_type%separate_time_vary_calc = .false.

tpi = 2.0*PI ! 4.*acos(0.)
dtr = tpi/360.

num_fields = 0

!--------------------------------------------------------------------
! open source file containing fields to be interpolated
!--------------------------------------------------------------------
src_file = 'INPUT/'//trim(file_name)

if(file_exist(trim(src_file))) then
   call mpp_open( unit, trim(src_file), action=MPP_RDONLY, &
                  form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
else
!Climatology file doesn't exist, so exit
   call mpp_error(FATAL,'Interpolator_init : Data file '//trim(src_file)//' does not exist')
endif

!Find the number of variables (nvar) in this file
call mpp_get_info(unit, ndim, nvar, natt, ntime)
clim_type%unit      = unit
clim_type%file_name = trim(file_name)


num_fields = nvar
if(present(data_names)) num_fields= size(data_names(:))

! -------------------------------------------------------------------
! Allocate space for the number of axes in the data file.
! -------------------------------------------------------------------
allocate(axes(ndim))
call mpp_get_axes(unit, axes, time_axis)

nlon=0 ! Number of longitudes (center-points) in the climatology.
nlat=0 ! Number of latitudes (center-points) in the climatology.
nlev=0 ! Number of levels (center-points) in the climatology.
nlatb=0 ! Number of longitudes (boundaries) in the climatology.
nlonb=0 ! Number of latitudes (boundaries) in the climatology.
nlevh=0 ! Number of levels (boundaries) in the climatology.

clim_type%level_type = 0 ! Default value

!++lwh
! -------------------------------------------------------------------
! For 2-D fields, set a default value of nlev=nlevh=1
! -------------------------------------------------------------------
nlev = 1
nlevh = 1
!--lwh
        clim_type%vertical_indices = 0  ! initial value

do i = 1, ndim
  call mpp_get_atts(axes(i), name=name,len=len,units=units,  &
                    calendar=file_calendar, sense=sense)
  select case(name)
    case('lat')
      nlat=len
      allocate(clim_type%lat(nlat))
      call mpp_get_axis_data(axes(i),clim_type%lat)
      select case(units(1:6))
        case('degree')
          clim_type%lat = clim_type%lat*dtr
        case('radian')
        case default  
          call mpp_error(FATAL, "interpolator_init : Units for lat not recognised in file "//file_name)
      end select
    case('lon')
      nlon=len
      allocate(clim_type%lon(nlon))
      call mpp_get_axis_data(axes(i),clim_type%lon)
      select case(units(1:6))
        case('degree')
          clim_type%lon = clim_type%lon*dtr
        case('radian')
        case default  
          call mpp_error(FATAL, "interpolator_init : Units for lon not recognised in file "//file_name)
      end select
    case('latb')
      nlatb=len
      allocate(clim_type%latb(nlatb))
      call mpp_get_axis_data(axes(i),clim_type%latb)
      select case(units(1:6))
        case('degree')
          clim_type%latb = clim_type%latb*dtr
        case('radian')
        case default  
          call mpp_error(FATAL, "interpolator_init : Units for latb not recognised in file "//file_name)
      end select
    case('lonb')
      nlonb=len
      allocate(clim_type%lonb(nlonb))
      call mpp_get_axis_data(axes(i),clim_type%lonb)
      select case(units(1:6))
        case('degree')
          clim_type%lonb = clim_type%lonb*dtr
        case('radian')
        case default  
          call mpp_error(FATAL, "interpolator_init : Units for lonb not recognised in file "//file_name)
      end select
    case('pfull')
      nlev=len
      allocate(clim_type%levs(nlev))
      call mpp_get_axis_data(axes(i),clim_type%levs)
      clim_type%level_type = PRESSURE
  ! Convert to Pa
      if( chomp(units) == "mb" .or. chomp(units) == "hPa") then
         clim_type%levs = clim_type%levs * 100.
      end if
! define the direction of the vertical data axis
! switch index order if necessary so that indx 1 is at lowest pressure,
! index nlev at highest pressure.
      if( sense == 1 ) then
        clim_type%vertical_indices = INCREASING_UPWARD
          allocate (alpha(nlev))
          do n = 1, nlev
          alpha(n) = clim_type%levs(nlev-n+1)
          end do
          do n = 1, nlev
          clim_type%levs(n) = alpha(n)
          end do
          deallocate (alpha)
      else 
        clim_type%vertical_indices = INCREASING_DOWNWARD
      endif
      
    case('phalf')
      nlevh=len
      allocate(clim_type%halflevs(nlevh))
      call mpp_get_axis_data(axes(i),clim_type%halflevs)
      clim_type%level_type = PRESSURE
  ! Convert to Pa
      if( chomp(units) == "mb" .or. chomp(units) == "hPa") then
         clim_type%halflevs = clim_type%halflevs * 100.
      end if
! define the direction of the vertical data axis
! switch index order if necessary so that indx 1 is at lowest pressure,
! index nlev at highest pressure.
      if( sense == 1 ) then
        clim_type%vertical_indices = INCREASING_UPWARD
          allocate (alpha(nlevh))
          do n = 1, nlevh
          alpha(n) = clim_type%halflevs(nlevh-n+1)
          end do
          do n = 1, nlevh
          clim_type%halflevs(n) = alpha(n)
          end do
          deallocate (alpha)
      else 
        clim_type%vertical_indices = INCREASING_DOWNWARD
      endif
    case('sigma_full')
      nlev=len
      allocate(clim_type%levs(nlev))
      call mpp_get_axis_data(axes(i),clim_type%levs)
      clim_type%level_type = SIGMA
    case('sigma_half')
      nlevh=len
      allocate(clim_type%halflevs(nlevh))
      call mpp_get_axis_data(axes(i),clim_type%halflevs)
      clim_type%level_type = SIGMA
    
    case('time')
      model_calendar = get_calendar_type() 
      fileday = 0
      filemon = 0
      fileyr = 0
      filehr = 0
      filemin= 0
      filesec = 0
      select case(units(:3))
        case('day')
          fileunits = units(12:) !Assuming "days since YYYY-MM-DD HH:MM:SS"
          read(fileunits(1:4)  , *)  fileyr
          read(fileunits(6:7)  , *)  filemon
          read(fileunits(9:10) , *)  fileday
          read(fileunits(12:13), *)  filehr
          read(fileunits(15:16), *)  filemin
          read(fileunits(18:19), *)  filesec
        case('mon')
          fileunits = units(14:) !Assuming "months since YYYY-MM-DD HH:MM:SS"
          read(fileunits(1:4)  , *)  fileyr
          read(fileunits(6:7)  , *)  filemon
          read(fileunits(9:10) , *)  fileday
          read(fileunits(12:13), *)  filehr
          read(fileunits(15:16), *)  filemin
          read(fileunits(18:19), *)  filesec
        case default
          call mpp_error(FATAL,'Interpolator_init : Time units not recognised in file '//file_name)
      end select

       clim_type%climatological_year = (fileyr == 0)
      if (.not. clim_type%climatological_year) then

!----------------------------------------------------------------------
!    if file date has a non-zero year in the base time, determine that
!    base_time based on the netcdf info.
!----------------------------------------------------------------------
        if ( (model_calendar == JULIAN .and.   &
              trim(file_calendar) == 'julian')  .or. &
              (model_calendar == NOLEAP .and.   &
               trim(file_calendar) == 'noleap')  .or. &
              (model_calendar == THIRTY_DAY_MONTHS .and. & !mj
               trim(file_calendar) == 'thirty_day_months'))  then
          call mpp_error (NOTE, 'interpolator_mod: Model and file&
                    & calendars are the same for file ' //   &
                    & trim(file_name) // '; no calendar conversion  &
                    &needed')
          base_time = set_date (fileyr, filemon, fileday, filehr, &
                                filemin,filesec)
        else if ( (model_calendar == JULIAN .and.   &
                   trim(file_calendar) == 'noleap')) then  
          call mpp_error (NOTE, 'interpolator_mod: Using julian &
                            &model calendar and noleap file calendar&
                            & for file ' // trim(file_name) //   &
                            &'; calendar conversion needed')
          base_time = set_date_no_leap (fileyr, filemon, fileday,  &
                                        filehr, filemin, filesec)
        else if ( (model_calendar == NOLEAP .and.   &
                   trim(file_calendar) == 'julian')) then  
          call mpp_error (NOTE, 'interpolator_mod: Using noleap &
                            &model calendar and julian file calendar&
                            & for file ' // trim(file_name) //  &
                            &'; calendar conversion needed')
          base_time = set_date_julian (fileyr, filemon, fileday,  &
                                       filehr, filemin, filesec)
        else
          call mpp_error (FATAL , 'interpolator_mod: Model and file&
               & calendars ( ' // trim(file_calendar) // ' ) differ  &
               &for file ' // trim(file_name) // ';  this calendar  &
               &conversion not currently available')
        endif

      else

!! if the year is specified as '0000', then the file is intended to 
!! apply to all years -- the time variables within the file refer to
!! the displacement from the start of each year to the time of the
!! associated data. Time interpolation is to be done with interface 
!! time_interp_list, with the optional argument modtime=YEAR. base_time
!! is set to an arbitrary value here; it's only use will be as a 
!! timestamp for optionally generated diagnostics.

        base_time = get_base_time ()
      endif


      ntime_in = 1
      if (ntime > 0) then
        allocate(time_in(ntime), clim_type%time_slice(ntime))
        allocate(clim_type%clim_times(12,(ntime+11)/12))
        time_in = 0.0
        clim_type%time_slice = set_time(0,0) + base_time
        clim_type%clim_times = set_time(0,0) + base_time
        call mpp_get_times(clim_type%unit, time_in)
        ntime_in = ntime
! determine whether the data is a continuous set of monthly values or
! a series of annual cycles spread throughout the period of data
        non_monthly = .false.
        do n = 1, ntime-1
!  Assume that the times in the data file correspond to days only.
          if (time_in(n+1) > (time_in(n) + 32)) then
            non_monthly = .true.
            exit
          endif
        end do
        if (clim_type%climatological_year) then
          call mpp_error (NOTE, 'interpolator_mod :'  // &
          trim(file_name) // ' is a year-independent climatology file') 
        else
          call mpp_error (NOTE, 'interpolator_mod :' // &
            trim(file_name) // ' is a timeseries file') 
        endif

        do n = 1, ntime
!Assume that the times in the data file correspond to days only.
            

          if (clim_type%climatological_year) then
!! RSH NOTE:
!! for this case, do not add base_time. time_slice will be sent to
!! time_interp_list with the optional argument modtime=YEAR, so that
!! the time that is needed in time_slice is the displacement into the
!! year, not the displacement from a base_time.
            clim_type%time_slice(n) = &
                set_time(INT( ( time_in(n) - INT(time_in(n)) ) * 86400 ),INT(time_in(n)))
          else

!--------------------------------------------------------------------
!    if fileyr /= 0 (i.e., climatological_year=F),
!    then define the times associated with each time-
!    slice. if calendar conversion between data file and model calendar
!    is needed, do it so that data from the file is associated with the
!    same calendar time in the model. here the time_slice needs to 
!    include the base_time; values will be generated relative to the 
!    "real" time.
!--------------------------------------------------------------------
            if ( (model_calendar == JULIAN .and.   &
                  trim(file_calendar) == 'julian')  .or. &
                 (model_calendar == NOLEAP .and.   &
                  trim(file_calendar) == 'noleap') .or. &
                 (model_calendar == THIRTY_DAY_MONTHS .and. & !mj
                  trim(file_calendar) == 'thirty_day_months') )  then

!---------------------------------------------------------------------
!    no calendar conversion needed.
!---------------------------------------------------------------------
              clim_type%time_slice(n) = &
                 set_time(INT( ( time_in(n) - INT(time_in(n)) ) * 86400 ),INT(time_in(n)))  &
                  + base_time

!---------------------------------------------------------------------
!    convert file times from noleap to julian.
!---------------------------------------------------------------------
            else if ( (model_calendar == JULIAN .and.   &
                       trim(file_calendar) == 'noleap')) then  
              Noleap_time = set_time (0, INT(time_in(n))) + base_time
              call get_date_no_leap (Noleap_time, yr, mo, dy, hr,  &
                                     mn, sc)
              clim_type%time_slice(n) = set_date_julian (yr, mo, dy,  &
                                                         hr, mn, sc)
              if (n == 1) then
                call print_date (clim_type%time_slice(1), &
                        str= 'for file ' // trim(file_name) // ', the &
                              &first time slice is mapped to :')
              endif
              if (n == ntime) then
                call print_date (clim_type%time_slice(ntime), &
                         str= 'for file ' // trim(file_name) // ', the &
                               &last time slice is mapped to:')
              endif
  

!---------------------------------------------------------------------
!    convert file times from julian to noleap.
!---------------------------------------------------------------------
            else if ( (model_calendar == NOLEAP .and.   &
                       trim(file_calendar) == 'julian')) then  
              Julian_time = set_time (0, INT(time_in(n))) + base_time
              call get_date_julian (Julian_time, yr, mo, dy, hr, mn, sc)
              clim_type%time_slice(n) = set_date_no_leap (yr, mo, dy, &
                                                          hr, mn, sc)
              if (n == 1) then
                call print_date (clim_type%time_slice(1), &
                         str= 'for file ' // trim(file_name) // ', the &
                               &first time slice is mapped to :')
              endif
              if (n == ntime) then
                call print_date (clim_type%time_slice(ntime), &
                         str= 'for file ' // trim(file_name) // ', the &
                               &last time slice is mapped to:')
              endif

!---------------------------------------------------------------------
!    any other calendar combinations would have caused a fatal error 
!    above.
!---------------------------------------------------------------------
            endif
          endif

          m = (n-1)/12 +1 ; m1 = n- (m-1)*12
          clim_type%clim_times(m1,m) = clim_type%time_slice(n)
        enddo
      else
        allocate(time_in(1), clim_type%time_slice(1))
        allocate(clim_type%clim_times(1,1))
        time_in = 0.0
        clim_type%time_slice = set_time(0,0) + base_time
        clim_type%clim_times(1,1) = set_time(0,0) + base_time
      endif
      deallocate(time_in)
  end select ! case(name)
enddo


! -------------------------------------------------------------------
! For 2-D fields, allocate levs and halflevs here
!  code is still needed for case when only halflevs are in data file.
! -------------------------------------------------------------------
    if( .not. associated(clim_type%levs) ) then
        allocate( clim_type%levs(nlev) )
        clim_type%levs = 0.0        
    endif  
    if( .not. associated(clim_type%halflevs) )  then
        allocate( clim_type%halflevs(nlev+1) )
        clim_type%halflevs(1) = 0.0
        if (clim_type%level_type == PRESSURE) then
          clim_type%halflevs(nlev+1) = 1013.25* 100.0   ! MKS
        else if (clim_type%level_type == SIGMA   ) then
          clim_type%halflevs(nlev+1) = 1.0
        endif
        do n=2,nlev
           clim_type%halflevs(n) = 0.5*(clim_type%levs(n) + &
                                         clim_type%levs(n-1))
        end do
    endif
deallocate(axes)


! In the case where only the midpoints of the longitudes are defined we force the definition
! of the boundaries to be half-way between the midpoints.
if (.not. associated(clim_type%lon) .and. .not. associated(clim_type%lonb)) &
   call mpp_error(FATAL,'Interpolator_init : There appears to be no longitude axis in file '//file_name)

if (.not. associated(clim_type%lonb) ) then

  if (size(clim_type%lon(:)) /= 1) then
    allocate(clim_type%lonb(size(clim_type%lon(:))+1))
    dlon = (clim_type%lon(2)-clim_type%lon(1))/2.0
    clim_type%lonb(1) = clim_type%lon(1) - dlon
    clim_type%lonb(2:) = clim_type%lon(1:) + dlon
  else

!! this is the case for zonal mean data, lon = 1, lonb not present 
!! in file.

    allocate(clim_type%lonb(2))
    clim_type%lonb(1) = -360.*dtr
    clim_type%lonb(2) = 360.0*dtr
    clim_type%lon(1) = 0.0
  endif    
endif

!clim_type%lonb=clim_type%lonb*dtr 
! This assumes the lonb are in degrees in the NetCDF file!

if (.not. associated(clim_type%lat) .and. .not. associated(clim_type%latb)) &
   call mpp_error(FATAL,'Interpolator_init : There appears to be no latitude axis in file '//file_name)
! In the case where only the grid midpoints of the latitudes are defined we force the 
! definition of the boundaries to be half-way between the midpoints.
if (.not. associated(clim_type%latb) ) then
   allocate(clim_type%latb(nlat+1))
   dlat = (clim_type%lat(2)-clim_type%lat(1)) * 0.5
!  clim_type%latb(1) = min( 90., max(-90., clim_type%lat(1) - dlat) )
   clim_type%latb(1) = min( PI/2., max(-PI/2., clim_type%lat(1) - dlat) )
   clim_type%latb(2:nlat) = ( clim_type%lat(1:nlat-1) + clim_type%lat(2:nlat) ) * 0.5
   dlat = ( clim_type%lat(nlat) - clim_type%lat(nlat-1) ) * 0.5
!  clim_type%latb(nlat+1) = min( 90., max(-90., clim_type%lat(nlat) + dlat) )
   clim_type%latb(nlat+1) = min( PI/2., max(-PI/2., clim_type%lat(nlat) + dlat) )
endif
!clim_type%latb=clim_type%latb*dtr

!Assume that the horizontal interpolation within a file is the same for each variable.

 if (trim(interp_method) == "conserve_latlon" ) then
    call horiz_interp_new (clim_type%interph, &
                        clim_type%lonb, clim_type%latb, &
                        lonb_mod, latb_mod)
 else if (trim(interp_method) == "conserve_great_circle") then
    call horiz_interp_new(clim_type%interph, &
                        clim_type%lonb, clim_type%latb, &
                        lonb_mod, latb_mod, interp_method="conserve_great_circle")
 else if( trim(interp_method) == "bilinear" ) then
    
    call mpp_error(NOTE, "Using Bilinear interpolation")

    !!! DEBUG CODE
    if (.not. allocated(agrid_mod)) then
       nx = size(lonb_mod,1)-1
       ny = size(latb_mod,2)-1
       allocate(agrid_mod(nx,ny,2))
       do j=1,ny
       do i=1,nx
          call cell_center2((/lonb_mod(i,j),latb_mod(i,j)/), & 
               (/lonb_mod(i+1,j),latb_mod(i+1,j)/), & 
               (/lonb_mod(i,j+1),latb_mod(i,j+1)/), & 
               (/lonb_mod(i+1,j+1),latb_mod(i+1,j+1)/),  agrid_mod(i,j,:))          
       enddo
       enddo
    endif

    !!! END DEBUG CODE

    call horiz_interp_new (clim_type%interph, &
                        clim_type%lonb, clim_type%latb, &
                        agrid_mod(:,:,1), agrid_mod(:,:,2), interp_method="bilinear")
 else
   call mpp_error(FATAL,'Interpolator_init : interp_methos should be "conserve_latlon", "conserve_great_circle" or "bilinear"')
    
 endif

!--------------------------------------------------------------------
!  allocate the variable clim_type%data . This will be the climatology 
!  data horizontally interpolated, so it will be on the model horizontal
!  grid, but it will still be on the climatology vertical grid.
!--------------------------------------------------------------------

select case(ntime)
 case (13:)
! This may  be data that does not have a continous time-line
! i.e. IPCC data where decadal data is present but we wish to retain 
! the seasonal nature of the data.
!! RSH: the following test will not always work; instead use the
!! RSH: non-monthly variable to test on.
!RSHlast_time = clim_type%time_slice(1) + ( ntime -1 ) * &
!RSH        ( clim_type%time_slice(2) - clim_type%time_slice(1) )

!RSHif ( last_time < clim_type%time_slice(ntime)) then

 if (non_monthly) then
! We have a broken time-line. e.g. We have monthly data but only for years ending in 0. 1960,1970 etc.
!   allocate(clim_type%data(size(lonb_mod(:))-1, size(latb_mod(:))-1, nlev, 2, num_fields))
   allocate(clim_type%pmon_pyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%pmon_nyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%nmon_nyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%nmon_pyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   clim_type%pmon_pyear = 0.0
   clim_type%pmon_nyear = 0.0
   clim_type%nmon_nyear = 0.0
   clim_type%nmon_pyear = 0.0
   clim_type%TIME_FLAG = BILINEAR
else
! We have a continuous time-line so treat as for 5-12 timelevels as below.
   if ( .not. read_all_on_init) then
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 2, num_fields))
   else
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, &
               ntime, num_fields))
   endif
   clim_type%data = 0.0
   clim_type%TIME_FLAG = LINEAR
endif


!++lwh
 case (1:12)
!--lwh
! We have more than 4 timelevels 
! Assume we have monthly or higher time resolution datasets (climatology or time series)
! So we only need to read 2 datasets and apply linear temporal interpolation.
   if ( .not. read_all_on_init) then
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 2, num_fields))
   else
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, &
               ntime, num_fields))
   endif
   clim_type%data = 0.0
   clim_type%TIME_FLAG = LINEAR
!++lwh
!case (1:4) 
! Assume we have seasonal data and read in all the data.
! We can apply sine curves to these data.
 
!  allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, ntime, num_fields))
!  clim_type%data = 0.0
!  clim_type%TIME_FLAG = SEASONAL
!--lwh
! case (default)
 case(:0)
   clim_type%TIME_FLAG = NOTIME
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 1, num_fields))

   call mpp_error(NOTE,"Interpolator: Using input file with no time dimension - "&
                             // trim(file_name))
end select


!------------------------------------------------------------------
!    Allocate space for the single time level of the climatology on its 
!    grid size.
!----------------------------------------------------------------------

   if(clim_type%TIME_FLAG .eq. LINEAR ) then
   allocate(clim_type%time_init(num_fields,2))
   else
   allocate(clim_type%time_init(num_fields,ntime))
   endif
   allocate (clim_type%indexm(num_fields),   &
             clim_type%indexp(num_fields),   &
             clim_type%climatology(num_fields))
   clim_type%time_init(:,:) = 0
   clim_type%indexm(:)      = 0
   clim_type%indexp(:)      = 0
   clim_type%climatology(:) = 0
   

allocate(clim_type%field_name(num_fields))
allocate(clim_type%field_type(num_fields))
allocate(clim_type%mr(num_fields))
allocate(clim_type%out_of_bounds(num_fields))
clim_type%out_of_bounds(:)=0
allocate(clim_type%vert_interp(num_fields))
clim_type%vert_interp(:)=0
!--------------------------------------------------------------------
!Allocate the space for the fields within the climatology data file.
allocate(varfields(nvar))
!--------------------------------------------------------------------
! Get the variable names out of the file.
call mpp_get_fields(clim_type%unit, varfields)

if(present(data_names)) then

!++lwh
   if ( size(data_out_of_bounds(:)) /= size(data_names(:)) .and. size(data_out_of_bounds(:)) /= 1 ) &
      call mpp_error(FATAL,'interpolator_init : The size of the data_out_of_bounds array must be 1&
                            & or size(data_names)')
   if (present(vert_interp)) then
      if( size(vert_interp(:)) /= size(data_names(:)) .and. size(vert_interp(:)) /= 1 ) &
      call mpp_error(FATAL,'interpolator_init : The size of the vert_interp array must be 1&
                            & or size(data_names)')
   endif
! Only read the fields named in data_names
   do j=1,size(data_names(:))
      NAME_PRESENT = .FALSE.
      do i=1,nvar
         call mpp_get_atts(varfields(i),name=name,ndim=ndim,units=units)
         if( name == data_names(j) ) then
            units=chomp(units)
            if (mpp_pe() == 0 ) write(*,*) 'Initializing src field : ',trim(name)
            clim_type%field_name(j) = name
            clim_type%field_type(j) = varfields(i)
            clim_type%mr(j)         = check_climo_units(units)
            NAME_PRESENT = .TRUE.
            if (present(clim_units)) clim_units(j) = units
            clim_type%out_of_bounds(j) = data_out_of_bounds( MIN(j,SIZE(data_out_of_bounds(:))) )
            if( clim_type%out_of_bounds(j) /= CONSTANT .and. &
                clim_type%out_of_bounds(j) /= ZERO ) &
               call mpp_error(FATAL,"Interpolator_init: data_out_of_bounds must be&
                                    & set to ZERO or CONSTANT")               
            if( present(vert_interp) ) then
               clim_type%vert_interp(j) = vert_interp( MIN(j,SIZE(vert_interp(:))) )
               if( clim_type%vert_interp(j) /= INTERP_WEIGHTED_P .and. &
                   clim_type%vert_interp(j) /= INTERP_LINEAR_P ) &
                  call mpp_error(FATAL,"Interpolator_init: vert_interp must be&
                                       & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
            else
               clim_type%vert_interp(j) = INTERP_WEIGHTED_P
            end if
         endif
      enddo
      if(.not. NAME_PRESENT) &
         call mpp_error(FATAL,'interpolator_init : Check names of fields being passed. ' &
                              //trim(data_names(j))//' does not exist.')
   enddo
else

   if ( size(data_out_of_bounds(:)) /= nvar .and. size(data_out_of_bounds(:)) /= 1 ) &
      call mpp_error(FATAL,'interpolator_init : The size of the out of bounds array must be 1&
                           & or the number of fields in the climatology dataset')
   if ( present(vert_interp) ) then
      if (size(vert_interp(:)) /= nvar .and. size(vert_interp(:)) /= 1 ) & 
      call mpp_error(FATAL,'interpolator_init : The size of the vert_interp array must be 1&
                           & or the number of fields in the climatology dataset')
   endif

! Read all the fields within the climatology data file.
   do i=1,nvar
      call mpp_get_atts(varfields(i),name=name,ndim=ndim,units=units)
         if (mpp_pe() ==0 ) write(*,*) 'Initializing src field : ',trim(name)
         clim_type%field_name(i) = lowercase(trim(name))
         clim_type%field_type(i) = varfields(i)
         clim_type%mr(i)         = check_climo_units(units)
         if (present(clim_units)) clim_units(i) = units
         clim_type%out_of_bounds(i) = data_out_of_bounds( MIN(i,SIZE(data_out_of_bounds(:))) )
         if( clim_type%out_of_bounds(i) /= CONSTANT .and. &
             clim_type%out_of_bounds(i) /= ZERO ) &
            call mpp_error(FATAL,"Interpolator_init: data_out_of_bounds must be&
                                 & set to ZERO or CONSTANT")
         if( present(vert_interp) ) then
            clim_type%vert_interp(i) = vert_interp( MIN(i,SIZE(vert_interp(:))) )
            if( clim_type%vert_interp(i) /= INTERP_WEIGHTED_P .and. &
                clim_type%vert_interp(i) /= INTERP_LINEAR_P ) &
               call mpp_error(FATAL,"Interpolator_init: vert_interp must be&
                                    & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
         else
            clim_type%vert_interp(i) = INTERP_WEIGHTED_P
         end if
   end do
!--lwh
endif

deallocate(varfields)


if( clim_type%TIME_FLAG .eq. SEASONAL ) then
! Read all the data at this point.
   do i=1,num_fields
      do n = 1, ntime
         call read_data( clim_type, clim_type%field_type(i), &
                         clim_type%data(:,:,:,n,i), n, i, base_time )
      enddo
   enddo
endif

if( clim_type%TIME_FLAG .eq. LINEAR  .and. read_all_on_init) then
! Read all the data at this point.
   do i=1,num_fields
      do n = 1, ntime
         call read_data( clim_type, clim_type%field_type(i), &
                         clim_type%data(:,:,:,n,i), n, i, base_time )
      enddo
   enddo

   call mpp_close (unit)
endif

if( clim_type%TIME_FLAG .eq. NOTIME ) then
! Read all the data at this point.
   do i=1,num_fields
     call read_data_no_time_axis( clim_type, clim_type%field_type(i), &
                                  clim_type%data(:,:,:,1,i), i )
   enddo
   call mpp_close (unit)
endif

if (present (single_year_file)) then
  single_year_file = clim_type%climatological_year
endif

module_is_initialized = .true.

call write_version_number (version, tagname)

end subroutine interpolator_init

 subroutine cell_center2(q1, q2, q3, q4, e2)
      real , intent(in ) :: q1(2), q2(2), q3(2), q4(2)
      real , intent(out) :: e2(2)
! Local
      real p1(3), p2(3), p3(3), p4(3)
      real ec(3)
      real dd
      integer k

      call latlon2xyz(q1, p1)
      call latlon2xyz(q2, p2)
      call latlon2xyz(q3, p3)
      call latlon2xyz(q4, p4)

      do k=1,3
         ec(k) = p1(k) + p2(k) + p3(k) + p4(k)
      enddo
      dd = sqrt( ec(1)**2 + ec(2)**2 + ec(3)**2 )

      do k=1,3
         ec(k) = ec(k) / dd
      enddo

      call cart_to_latlon(1, ec, e2(1), e2(2))

 end subroutine cell_center2

 subroutine cart_to_latlon(np, q, xs, ys)
! vector version of cart_to_latlon1
  integer, intent(in):: np
  real, intent(inout):: q(3,np)
  real, intent(inout):: xs(np), ys(np)
! local
  real, parameter:: esl=1.e-10
  real (f_p):: p(3)
  real (f_p):: dist, lat, lon
  integer i,k

  do i=1,np
     do k=1,3
        p(k) = q(k,i)
     enddo
     dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
     do k=1,3
        p(k) = p(k) / dist
     enddo

     if ( (abs(p(1))+abs(p(2)))  < esl ) then
          lon = 0.
     else
          lon = atan2( p(2), p(1) )   ! range [-pi,pi]
     endif

     if ( lon < 0.) lon = 2.*pi + lon
     lat = asin(p(3))
     
     xs(i) = lon
     ys(i) = lat
! q Normalized:
     do k=1,3
        q(k,i) = p(k)
     enddo
  enddo

 end  subroutine cart_to_latlon

 subroutine latlon2xyz(p, e)
!
! Routine to map (lon, lat) to (x,y,z)
!
 real, intent(in) :: p(2)
 real, intent(out):: e(3)

 integer n
 real (f_p):: q(2)
 real (f_p):: e1, e2, e3

    do n=1,2
       q(n) = p(n)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3

 end subroutine latlon2xyz

!
!#######################################################################
!
function check_climo_units(units)
! Function to check the units that the climatology data is using. 
! This is needed to allow for conversion of datasets to mixing ratios which is what the 
! vertical interpolation scheme requires
! The default is to assume no conversion is needed.
! If the units are those of a column burden (kg/m2) then conversion to mixing ratio is flagged.
!
character(len=*), intent(in) :: units

integer :: check_climo_units

check_climo_units = NO_CONV
select case(chomp(units))
  case('kg/m2', 'kg/m^2', 'kg/m**2', 'kg m^-2', 'kg m**-2')
     check_climo_units = KG_M2  
  case('molecules/cm2/s', 'molecule/cm2/s', 'molec/cm2/s')
     check_climo_units = KG_M2  
  case('kg/m2/s')
     check_climo_units = KG_M2  
end select

end function check_climo_units
!
!#######################################################################
!
subroutine init_clim_diag(clim_type, mod_axes, init_time)
!
! Routine to register diagnostic fields for the climatology file. 
! This routine calculates the domain decompostion of the climatology fields 
! for later export through send_data.
! The ids created here are for column burdens that will diagnose the vertical interpolation routine.
! climo_diag_id : 'module_name = climo' is intended for use with the model vertical resolution.
! hinterp_id    : 'module_name = 'hinterp' is intended for use with the climatology vertical resolution.

! INTENT INOUT :
!    clim_type : The interpolate type containing the names of the fields in the climatology file.
!
! INTENT IN    :
!   mod_axes   : The axes of the model.
!   init_time  : The model initialization time.
!
type(interpolate_type), intent(inout)  :: clim_type
integer               , intent(in)     :: mod_axes(:)
type(time_type)       , intent(in)     :: init_time

integer :: axes(2),nxd,nyd,ndivs,i
type(domain2d) :: domain
integer :: domain_layout(2), iscomp, iecomp,jscomp,jecomp


if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "init_clim_diag : You must call interpolator_init before calling init_clim_diag")


ndivs = mpp_npes()
nxd = size(clim_type%lon(:))
nyd = size(clim_type%lat(:))

! Define the domain decomposition of the climatology file. This may be (probably is) different from the model domain.
call mpp_define_layout ((/1,nxd,1,nyd/), ndivs, domain_layout)
call mpp_define_domains((/1,nxd,1,nyd/),domain_layout, domain,xhalo=0,yhalo=0)  
call mpp_get_compute_domain (domain, iscomp, iecomp, jscomp, jecomp)
   axes(1) = diag_axis_init(clim_type%file_name(1:5)//'x',clim_type%lon,units='degrees',cart_name='x',domain2=domain)
   axes(2) = diag_axis_init(clim_type%file_name(1:5)//'y',clim_type%lat,units='degrees',cart_name='y',domain2=domain)
clim_type%is = iscomp
clim_type%ie = iecomp
clim_type%js = jscomp
clim_type%je = jecomp

!init_time = set_date(1980,1,1,0,0,0)

if ((num_clim_diag + size(clim_type%field_name(:))) .gt. max_diag_fields )  &
   call mpp_error(FATAL, "init_clim_diag : Trying to set up too many diagnostic fields for the climatology data")
do i=1,size(clim_type%field_name(:))
climo_diag_name(i+num_clim_diag) = clim_type%field_name(i)
climo_diag_id(i+num_clim_diag) =  register_diag_field('climo',clim_type%field_name(i),axes(1:2),init_time,&
                                'climo_'//clim_type%field_name(i), 'kg/kg', missing_value)
hinterp_id(i+num_clim_diag) =  register_diag_field('hinterp',clim_type%field_name(i),mod_axes(1:2),init_time,&
                                'interp_'//clim_type%field_name(i),'kg/kg' , missing_value)
enddo
! Total number of climatology diagnostics (num_clim_diag). This can be from multiple climatology fields with different spatial axes. 
! It is simply a holder for the diagnostic indices.
num_clim_diag = num_clim_diag+size(clim_type%field_name(:))

clim_diag_initialized = .true.

end subroutine init_clim_diag



!----------------------------------------------------------------------------

subroutine obtain_interpolator_time_slices (clim_type, Time)

!  Makes sure that appropriate time slices are available for interpolation 
!  on this time step
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   Time        : The model time that you wish to interpolate to.

type(interpolate_type), intent(inout)  :: clim_type
type(time_type)       , intent(in)  :: Time

integer :: taum, taup
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: i, n


    if (clim_type%climatological_year) then
!++lwh
       if (size(clim_type%time_slice) > 1) then
          call time_interp(Time, clim_type%time_slice, clim_type%tweight, taum, taup, modtime=YEAR )
       else
          taum = 1
          taup = 1
          clim_type%tweight = 0.
       end if
!--lwh
    else
       call time_interp(Time, clim_type%time_slice, clim_type%tweight, taum, taup )
    endif


    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
      ! The difference between the model time and the last climatology time-slice previous to the model time.
      ! We need 2 time levels.
        clim_type%itaum=0
        clim_type%itaup=0
      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute 
      ! the climatology year into the appropriate place.

     
      ! We need to get the previous months data for the climatology year before 
      ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:,:),2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1,m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m 
        enddo
        do m = 1,12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m,climatology), year1, month1, day, hour, minute, second)
          if ( month1 == modmonth ) then
!RSHBUGFX   if ( modday <= day ) then 
            if ( modday <  day ) then 
              indexm = m-1 ; indexp = m
            else
              indexm = m ; indexp = m+1
            endif
          endif
        
        enddo
        if ( indexm == 0 ) then 
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        endif
          call get_date(clim_type%time_slice(indexm+(climatology-1)*12), &
                        climyear, climmonth, climday, climhour, climminute, climsecond)
          month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if ( indexp == 13 ) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        endif
          call get_date(clim_type%time_slice(indexp+(climatology-1)*12), &
                        climyear, climmonth, climday, climhour, climminute, climsecond)
          month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)
        
        call time_interp(Time, month, clim_type%tweight3, taum, taup ) ! tweight3 is the time weight between the months.

        month(1) = clim_type%time_slice(indexm+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexm+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, clim_type%tweight1, taum, taup ) !tweight1 is the time weight between the climatology years.
        month(1) = clim_type%time_slice(indexp+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexp+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, clim_type%tweight2, taum, taup ) !tweight1 is the time weight between the climatology years.

        if (indexm == clim_type%indexm(1) .and.  &
            indexp == clim_type%indexp(1) .and. &
            climatology == clim_type%climatology(1)) then
        else
          clim_type%indexm(:) = indexm
          clim_type%indexp(:) = indexp
          clim_type%climatology(:) = climatology
          do i=1, size(clim_type%field_name(:))
            call read_data(clim_type,clim_type%field_type(i),  &
             clim_type%pmon_pyear(:,:,:,i),   &
             clim_type%indexm(i)+(clim_type%climatology(i)-1)*12,i,Time)
! Read the data for the next month in the previous climatology.
            call read_data(clim_type,clim_type%field_type(i),  &
             clim_type%nmon_pyear(:,:,:,i),   &
             clim_type%indexp(i)+(clim_type%climatology(i)-1)*12,i,Time)
            call read_data(clim_type,clim_type%field_type(i),  &
              clim_type%pmon_nyear(:,:,:,i),  &
              clim_type%indexm(i)+clim_type%climatology(i)*12,i,Time)
            call read_data(clim_type,clim_type%field_type(i),  &
              clim_type%nmon_nyear(:,:,:,i),  &
              clim_type%indexp(i)+clim_type%climatology(i)*12,i,Time)
          end do
        endif



      else ! We are within a climatology data set
        

        do i=1, size(clim_type%field_name(:))
          if (taum /= clim_type%time_init(i,1) .or. &
              taup /= clim_type%time_init(i,2) ) then
 
     
            call read_data(clim_type,clim_type%field_type(i),   &
                           clim_type%pmon_pyear(:,:,:,i), taum,i,Time)
! Read the data for the next month in the previous climatology.
            call read_data(clim_type,clim_type%field_type(i),   &
                           clim_type%nmon_pyear(:,:,:,i), taup,i,Time)
            clim_type%time_init(i,1) = taum
            clim_type%time_init(i,2) = taup
          endif
        end do

!       clim_type%pmon_nyear = 0.0
!       clim_type%nmon_nyear = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
          clim_type%indexm(:) = 0       
          clim_type%indexp(:) = 0        
          clim_type%climatology(:) = 0             


!       clim_type%tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        clim_type%tweight1 = 0.0 
        clim_type%tweight2 = 0.0 
        clim_type%tweight3 = clim_type%tweight                                          
      endif
    endif   !(BILINEAR)

    if(clim_type%TIME_FLAG .eq. LINEAR  .and.   &
        (.not. read_all_on_init) ) then
! We need 2 time levels. Check we have the correct data.
      clim_type%itaum=0
      clim_type%itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(1,n) .eq. taum ) clim_type%itaum = n
        if (clim_type%time_init(1,n) .eq. taup ) clim_type%itaup = n
      enddo

      if (clim_type%itaum.eq.0 .and. clim_type%itaup.eq.0) then
!Neither time is set so we need to read 2 time slices.
!Set up 
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
    do i=1, size(clim_type%field_name(:))
    call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,1,i), taum,i,Time)
          clim_type%time_init(i,1) = taum
          clim_type%itaum = 1
    call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,2,i), taup,i,Time)
          clim_type%time_init(i,2) = taup
          clim_type%itaup = 2
    end do
      endif ! clim_type%itaum.eq.clim_type%itaup.eq.0
      if (clim_type%itaum.eq.0 .and. clim_type%itaup.ne.0) then
! Can't think of a situation where we would have the next time level but not the previous.
 call mpp_error(FATAL,'interpolator_timeslice : No data from the previous climatology time &
                         & but we have the next time. How did this happen?')
      endif
      if (clim_type%itaum.ne.0 .and. clim_type%itaup.eq.0) then
!We have the previous time step but not the next time step data
        clim_type%itaup = 1
        if (clim_type%itaum .eq. 1 ) clim_type%itaup = 2
    do i=1, size(clim_type%field_name(:))
        call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,clim_type%itaup,i), taup,i, Time)
        clim_type%time_init(i,clim_type%itaup)=taup
     end do
      endif


    endif! TIME_FLAG

    clim_type%separate_time_vary_calc = .true.

!-------------------------------------------------------------------


end subroutine obtain_interpolator_time_slices 


!#####################################################################

subroutine unset_interpolator_time_flag (clim_type)

type(interpolate_type), intent(inout) :: clim_type


      clim_type%separate_time_vary_calc = .false.


end subroutine unset_interpolator_time_flag 



!#####################################################################

!---------------------------------------------------------------------

subroutine interpolator_4D(clim_type, Time, phalf, interp_data,  &
                           field_name, is,js, clim_units)
!
! Return 4-D field interpolated to model grid and time
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   field_name  : The name of a field that you wish to interpolate.
!                 all variables within this interpolate_type variable
!                 will be interpolated on this call. field_name may
!                 be any one of the variables.
!   Time        : The model time that you wish to interpolate to.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.
!
! INTENT OUT
!   interp_data : The model fields with the interpolated climatology data.
!   clim_units  : The units of field_name
!
type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
type(time_type)       , intent(in)  :: Time
real, dimension(:,:,:), intent(in)  :: phalf
real, dimension(:,:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: taum, taup, ilon
real :: hinterp_data(size(interp_data,1),size(interp_data,2),size(clim_type%levs(:)),size(clim_type%field_name(:)))
real :: p_fact(size(interp_data,1),size(interp_data,2))
real :: col_data(size(interp_data,1),size(interp_data,2),   &
                           size(clim_type%field_name(:)))
real :: pclim(size(clim_type%halflevs(:)))
integer :: istart,iend,jstart,jend
logical :: result, found
logical :: found_field=.false.
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: i, j, k, n

!RG add option to go to 'no time axis' version if clim_type%TIME_FLAG==NOTIME i.e. data has no time, despite input received
if (clim_type%TIME_FLAG .eq. NOTIME) then
	call interpolator_4D_no_time_axis(clim_type, phalf, interp_data, field_name, is,js, clim_units)
	return
endif

if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "interpolator_4D : You must call interpolator_init before calling interpolator")

   do n=2,size(clim_type%field_name(:))
     if (clim_type%vert_interp(n) /= clim_type%vert_interp(n-1) .or. &
      clim_type%out_of_bounds(n) /= clim_type%out_of_bounds(n-1)) then
       if (mpp_pe() == mpp_root_pe() ) then
         print *, 'processing file ' // trim(clim_type%file_name)
       endif
       call mpp_error (FATAL, 'interpolator_mod: &
               &cannot use 4D interface to interpolator for this file')
     endif
   end do
     



istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

  do i= 1,size(clim_type%field_name(:))
!!++lwh
   if ( field_name == clim_type%field_name(i) ) then
!--lwh
    found_field=.true.
    exit 
 endif
end do
   i = 1

    if(present(clim_units)) then
      call mpp_get_atts(clim_type%field_type(i),units=clim_units)
      clim_units = chomp(clim_units)
    endif




!----------------------------------------------------------------------
!   skip the time interpolation portion of this routine if subroutine
!   obtain_interpolator_time_slices has already been called on this
!   stewp for this interpolate_type variable.
!----------------------------------------------------------------------

if ( .not. clim_type%separate_time_vary_calc) then     
!   print *, 'TIME INTERPOLATION NOT SEPARATED 4d--',  &
!                                trim(clim_type%file_name), mpp_pe()

    if (clim_type%climatological_year) then
!++lwh
       if (size(clim_type%time_slice) > 1) then
          call time_interp(Time, clim_type%time_slice, clim_type%tweight, taum, taup, modtime=YEAR )
       else
          taum = 1
          taup = 1
          clim_type%tweight = 0.
       end if
!--lwh
    else
       call time_interp(Time, clim_type%time_slice, clim_type%tweight, taum, taup )
    endif


    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
      ! The difference between the model time and the last climatology time-slice previous to the model time.
      ! We need 2 time levels.
        clim_type%itaum=0
        clim_type%itaup=0
      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute 
      ! the climatology year into the appropriate place.

     
      ! We need to get the previous months data for the climatology year before 
      ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:,:),2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1,m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m 
        enddo
        do m = 1,12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m,climatology), year1, month1, day, hour, minute, second)
          if ( month1 == modmonth ) then
!RSHBUGFX   if ( modday <= day ) then 
            if ( modday <  day ) then 
              indexm = m-1 ; indexp = m
            else
              indexm = m ; indexp = m+1
            endif
          endif
        
        enddo
        if ( indexm == 0 ) then 
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        endif
          call get_date(clim_type%time_slice(indexm+(climatology-1)*12), &
                        climyear, climmonth, climday, climhour, climminute, climsecond)
          month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if ( indexp == 13 ) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        endif
          call get_date(clim_type%time_slice(indexp+(climatology-1)*12), &
                        climyear, climmonth, climday, climhour, climminute, climsecond)
          month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)
        
        call time_interp(Time, month, clim_type%tweight3, taum, taup ) ! tweight3 is the time weight between the months.

        month(1) = clim_type%time_slice(indexm+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexm+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, clim_type%tweight1, taum, taup ) !tweight1 is the time weight between the climatology years.
        month(1) = clim_type%time_slice(indexp+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexp+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, clim_type%tweight2, taum, taup ) !tweight1 is the time weight between the climatology years.

        if (indexm == clim_type%indexm(1) .and.  &
            indexp == clim_type%indexp(1) .and. &
            climatology == clim_type%climatology(1)) then
        else
          clim_type%indexm(:) = indexm
          clim_type%indexp(:) = indexp
          clim_type%climatology(:) = climatology
          do i=1, size(clim_type%field_name(:))
            call read_data(clim_type,clim_type%field_type(i),  &
             clim_type%pmon_pyear(:,:,:,i),   &
             clim_type%indexm(i)+(clim_type%climatology(i)-1)*12,i,Time)
! Read the data for the next month in the previous climatology.
            call read_data(clim_type,clim_type%field_type(i),  &
             clim_type%nmon_pyear(:,:,:,i),   &
             clim_type%indexp(i)+(clim_type%climatology(i)-1)*12,i,Time)
            call read_data(clim_type,clim_type%field_type(i),  &
              clim_type%pmon_nyear(:,:,:,i),  &
              clim_type%indexm(i)+clim_type%climatology(i)*12,i,Time)
            call read_data(clim_type,clim_type%field_type(i),  &
              clim_type%nmon_nyear(:,:,:,i),  &
              clim_type%indexp(i)+clim_type%climatology(i)*12,i,Time)
          end do
        endif



      else ! We are within a climatology data set
        

        do i=1, size(clim_type%field_name(:))
          if (taum /= clim_type%time_init(i,1) .or. &
              taup /= clim_type%time_init(i,2) ) then
 
     
            call read_data(clim_type,clim_type%field_type(i),   &
                           clim_type%pmon_pyear(:,:,:,i), taum,i,Time)
! Read the data for the next month in the previous climatology.
            call read_data(clim_type,clim_type%field_type(i),   &
                           clim_type%nmon_pyear(:,:,:,i), taup,i,Time)
            clim_type%time_init(i,1) = taum
            clim_type%time_init(i,2) = taup
          endif
        end do

!       clim_type%pmon_nyear = 0.0
!       clim_type%nmon_nyear = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
          clim_type%indexm(:) = 0       
          clim_type%indexp(:) = 0        
          clim_type%climatology(:) = 0             


!       clim_type%tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        clim_type%tweight1 = 0.0 
        clim_type%tweight2 = 0.0 
        clim_type%tweight3 = clim_type%tweight                                          
      endif
    endif   !(BILINEAR)

    if(clim_type%TIME_FLAG .eq. LINEAR  .and.   &
        (.not. read_all_on_init) ) then
! We need 2 time levels. Check we have the correct data.
      clim_type%itaum=0
      clim_type%itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(1,n) .eq. taum ) clim_type%itaum = n
        if (clim_type%time_init(1,n) .eq. taup ) clim_type%itaup = n
      enddo

      if (clim_type%itaum.eq.0 .and. clim_type%itaup.eq.0) then
!Neither time is set so we need to read 2 time slices.
!Set up 
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
    do i=1, size(clim_type%field_name(:))
    call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,1,i), taum,i,Time)
          clim_type%time_init(i,1) = taum
          clim_type%itaum = 1
    call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,2,i), taup,i,Time)
          clim_type%time_init(i,2) = taup
          clim_type%itaup = 2
    end do
      endif ! clim_type%itaum.eq.clim_type%itaup.eq.0
      if (clim_type%itaum.eq.0 .and. clim_type%itaup.ne.0) then
! Can't think of a situation where we would have the next time level but not the previous.
 call mpp_error(FATAL,'interpolator_3D : No data from the previous climatology time &
                         & but we have the next time. How did this happen?')
      endif
      if (clim_type%itaum.ne.0 .and. clim_type%itaup.eq.0) then
!We have the previous time step but not the next time step data
        clim_type%itaup = 1
        if (clim_type%itaum .eq. 1 ) clim_type%itaup = 2
    do i=1, size(clim_type%field_name(:))
        call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,clim_type%itaup,i), taup,i, Time)
        clim_type%time_init(i,clim_type%itaup)=taup
     end do
      endif


    endif! TIME_FLAG


endif ! (.not. separate_time_vary_calc)


select case(clim_type%TIME_FLAG)
  case (LINEAR)
    do n=1, size(clim_type%field_name(:))
      hinterp_data(:,:,:,n) = (1-clim_type%tweight)*  &
                clim_type%data(istart:iend,jstart:jend,:,clim_type%itaum,n)  +  &
                                 clim_type%tweight*   &
                clim_type%data(istart:iend,jstart:jend,:,clim_type%itaup,n)
    end do
! case (SEASONAL)
! Do sine fit to data at this point
  case (BILINEAR)
    do n=1, size(clim_type%field_name(:))
      hinterp_data(:,:,:,n) = (1-clim_type%tweight1)*(1-clim_type%tweight3)*   &
                   clim_type%pmon_pyear(istart:iend,jstart:jend,:,n) + &
                              (1-clim_type%tweight2)*clim_type%tweight3*    &
                   clim_type%nmon_pyear(istart:iend,jstart:jend,:,n) + &
                               clim_type%tweight1* (1-clim_type%tweight3)*  &
                   clim_type%pmon_nyear(istart:iend,jstart:jend,:,n) + &
                               clim_type%tweight2* clim_type%tweight3*   &
                   clim_type%nmon_nyear(istart:iend,jstart:jend,:,n)
    
    end do

end select
    
select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select

col_data(:,:,:)=0.0
     do i= 1, size(clim_type%field_name(:))

select case(clim_type%mr(i))
  case(NO_CONV)
    do k = 1,size(hinterp_data,3)
   col_data(:,:,i) = col_data(:,:,i) + hinterp_data(:,:,k,i)* &
      (clim_type%halflevs(k+1)-clim_type%halflevs(k))/grav
    enddo
    
  case(KG_M2)
    do k = 1,size(hinterp_data,3)
       col_data(:,:,i) = col_data(:,:,i) + hinterp_data(:,:,k,i)
       hinterp_data(:,:,k,i) = hinterp_data(:,:,k,i)/ &
         ((clim_type%halflevs(k+1)-clim_type%halflevs(k))*p_fact)
    enddo
end select
    enddo

     do i= 1, size(clim_type%field_name(:))
found = .false.
do j = 1,size(climo_diag_name(:))
  if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
    found = .true.
    exit
  endif
enddo

if (found) then
  if (hinterp_id(j) > 0 ) then
       result = send_data(hinterp_id(j),col_data(:,:,i),Time)
  endif
endif

  end do

   i = 1

!++lwh
do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model surface pressure&
                             & is greater than climatology surface pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model top pressure&
                             & is less than climatology top pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,:),interp_data(ilon,j,:,:))
         case(INTERP_LINEAR_P)
          do n=1, size(clim_type%field_name(:))
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,n),interp_data(ilon,j,:,n))
          end do
!        case(INTERP_LOG)
      end select
   enddo
enddo

!--lwh
     do i= 1, size(clim_type%field_name(:))

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k,i) = interp_data(:,:,k,i)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select

     end do

if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine interpolator_4D
!
!#######################################################################
!#######################################################################
!
subroutine interpolator_3D(clim_type, Time, phalf, interp_data,field_name, is,js, clim_units)
!
! Return 3-D field interpolated to model grid and time
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   Time        : The model time that you wish to interpolate to.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.
!
! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name
!
type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
type(time_type)       , intent(in)  :: Time
real, dimension(:,:,:), intent(in)  :: phalf
real, dimension(:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: taum, taup, ilon
real :: hinterp_data(size(interp_data,1),size(interp_data,2),size(clim_type%levs(:)))
real :: p_fact(size(interp_data,1),size(interp_data,2))
real :: col_data(size(interp_data,1),size(interp_data,2))
real :: pclim(size(clim_type%halflevs(:)))
integer :: istart,iend,jstart,jend
logical :: result, found
logical :: found_field=.false.
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: i, j, k, n

!RG add option to go to 'no time axis' version if clim_type%TIME_FLAG==NOTIME i.e. data has no time, despite input received
if (clim_type%TIME_FLAG .eq. NOTIME) then
	call interpolator_3D_no_time_axis(clim_type, phalf, interp_data, field_name, is,js, clim_units)
	return
endif

if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "interpolator_3D : You must call interpolator_init before calling interpolator")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
!++lwh
  if ( field_name == clim_type%field_name(i) ) then
!--lwh
    found_field=.true.
    if(present(clim_units)) then
      call mpp_get_atts(clim_type%field_type(i),units=clim_units)
      clim_units = chomp(clim_units)
    endif

!----------------------------------------------------------------------
!   skip the time interpolation portion of this routine if subroutine
!   obtain_interpolator_time_slices has already been called on this
!   stewp for this interpolate_type variable.
!----------------------------------------------------------------------


if ( .not. clim_type%separate_time_vary_calc) then     
!   print *, 'TIME INTERPOLATION NOT SEPARATED 3d--',  &
!                                trim(clim_type%file_name), mpp_pe()
    if (clim_type%climatological_year) then
!++lwh
       if (size(clim_type%time_slice) > 1) then
          call time_interp(Time, clim_type%time_slice, clim_type%tweight, taum, taup, modtime=YEAR )
       else
          taum = 1
          taup = 1
          clim_type%tweight = 0.
       end if
!--lwh
    else
       call time_interp(Time, clim_type%time_slice, clim_type%tweight, taum, taup )
    endif

!   if(clim_type%TIME_FLAG .ne. LINEAR ) then
    if(clim_type%TIME_FLAG .ne. LINEAR .or. read_all_on_init ) then
      clim_type%itaum=taum
      clim_type%itaup=taup
    endif

    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
      ! The difference between the model time and the last climatology time-slice previous to the model time.
      ! We need 2 time levels.
        clim_type%itaum=0
        clim_type%itaup=0
      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute 
      ! the climatology year into the appropriate place.

     
      ! We need to get the previous months data for the climatology year before 
      ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:,:),2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1,m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m 
        enddo
        do m = 1,12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m,climatology), year1, month1, day, hour, minute, second)
          if ( month1 == modmonth ) then
!RSHBUGFX   if ( modday <= day ) then 
            if ( modday <  day ) then 
              indexm = m-1 ; indexp = m
            else
              indexm = m ; indexp = m+1
            endif
          endif
        
        enddo
        if ( indexm == 0 ) then 
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        endif
        call get_date(clim_type%time_slice(indexm+(climatology-1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if ( indexp == 13 ) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        endif
        call get_date(clim_type%time_slice(indexp+(climatology-1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)
        
        call time_interp(Time, month, clim_type%tweight3, taum, taup ) ! tweight3 is the time weight between the months.

        month(1) = clim_type%time_slice(indexm+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexm+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, clim_type%tweight1, taum, taup ) !tweight1 is the time weight between the climatology years.

        month(1) = clim_type%time_slice(indexp+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexp+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, clim_type%tweight2, taum, taup ) !tweight1 is the time weight between the climatology years.



        if (indexm == clim_type%indexm(i) .and.  &
          indexp == clim_type%indexp(i) .and. &
          climatology == clim_type%climatology(i)) then
        else
          clim_type%indexm(i) = indexm
          clim_type%indexp(i) = indexp
          clim_type%climatology(i) = climatology
          call read_data(clim_type,clim_type%field_type(i),  &
            clim_type%pmon_pyear(:,:,:,i),  &
            clim_type%indexm(i)+(clim_type%climatology(i)-1)*12,i,Time)
! Read the data for the next month in the previous climatology.
          call read_data(clim_type,clim_type%field_type(i),  &
            clim_type%nmon_pyear(:,:,:,i),   &
            clim_type%indexp(i)+(clim_type%climatology(i)-1)*12,i,Time)
          call read_data(clim_type,clim_type%field_type(i),   &
            clim_type%pmon_nyear(:,:,:,i),  &
            clim_type%indexm(i)+clim_type%climatology(i)*12,i,Time)
          call read_data(clim_type,clim_type%field_type(i),  &
            clim_type%nmon_nyear(:,:,:,i),  &
            clim_type%indexp(i)+clim_type%climatology(i)*12,i,Time)
        endif




      else ! We are within a climatology data set
        
        if (taum /= clim_type%time_init(i,1) .or. &
            taup /= clim_type%time_init(i,2) ) then
 
          call read_data(clim_type,clim_type%field_type(i), clim_type%pmon_pyear(:,:,:,i), taum,i,Time)
! Read the data for the next month in the previous climatology.
          call read_data(clim_type,clim_type%field_type(i), clim_type%nmon_pyear(:,:,:,i), taup,i,Time)
!RSHbug   clim_type%pmon_nyear = 0.0
!RSHbug   clim_type%nmon_nyear = 0.0

!         clim_type%pmon_nyear(:,:,:,i) = 0.0
!         clim_type%nmon_nyear(:,:,:,i) = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
          clim_type%indexm(i) = 0       
          clim_type%indexp(i) = 0        
          clim_type%climatology(i) = 0             


          clim_type%time_init(i,1) = taum
          clim_type%time_init(i,2) = taup
        endif
!       clim_type%tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        clim_type%tweight1 = 0.0 ; clim_type%tweight2 = 0.0
        clim_type%tweight3 = clim_type%tweight                                          
      endif

    endif ! (BILINEAR)


    if(clim_type%TIME_FLAG .eq. LINEAR  .and.   &
        (.not. read_all_on_init) ) then
! We need 2 time levels. Check we have the correct data.
      clim_type%itaum=0
      clim_type%itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(i,n) .eq. taum ) clim_type%itaum = n
        if (clim_type%time_init(i,n) .eq. taup ) clim_type%itaup = n
      enddo

      if (clim_type%itaum.eq.0 .and. clim_type%itaup.eq.0) then
!Neither time is set so we need to read 2 time slices.
!Set up 
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
    call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,1,i), taum,i,Time)
          clim_type%time_init(i,1) = taum
          clim_type%itaum = 1
    call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,2,i), taup,i,Time)
          clim_type%time_init(i,2) = taup
          clim_type%itaup = 2
      endif ! clim_type%itaum.eq.clim_type%itaup.eq.0
      if (clim_type%itaum.eq.0 .and. clim_type%itaup.ne.0) then
! Can't think of a situation where we would have the next time level but not the previous.
 call mpp_error(FATAL,'interpolator_3D : No data from the previous climatology time &
                         & but we have the next time. How did this happen?')
      endif
      if (clim_type%itaum.ne.0 .and. clim_type%itaup.eq.0) then
!We have the previous time step but not the next time step data
        clim_type%itaup = 1
        if (clim_type%itaum .eq. 1 ) clim_type%itaup = 2
        call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,clim_type%itaup,i), taup,i, Time)
        clim_type%time_init(i,clim_type%itaup)=taup
      endif


    endif! TIME_FLAG

    endif   !( .not. clim_type%separate_time_vary_calc) 
select case(clim_type%TIME_FLAG)
  case (LINEAR)
    hinterp_data = (1-clim_type%tweight) * clim_type%data(istart:iend,jstart:jend,:,clim_type%itaum,i) + &
                       clim_type%tweight * clim_type%data(istart:iend,jstart:jend,:,clim_type%itaup,i)
! case (SEASONAL)
! Do sine fit to data at this point
  case (BILINEAR)
    hinterp_data = &
    (1-clim_type%tweight1)  * (1-clim_type%tweight3) * clim_type%pmon_pyear(istart:iend,jstart:jend,:,i) + &
    (1-clim_type%tweight2)  *    clim_type%tweight3  * clim_type%nmon_pyear(istart:iend,jstart:jend,:,i) + &
         clim_type%tweight1 * (1-clim_type%tweight3) * clim_type%pmon_nyear(istart:iend,jstart:jend,:,i) + &
         clim_type%tweight2 *     clim_type%tweight3 * clim_type%nmon_nyear(istart:iend,jstart:jend,:,i)
    


end select

select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select

col_data(:,:)=0.0
select case(clim_type%mr(i))
  case(NO_CONV)
    do k = 1,size(hinterp_data,3)
   col_data(:,:) = col_data(:,:) + hinterp_data(:,:,k)* &
      (clim_type%halflevs(k+1)-clim_type%halflevs(k))/grav
    enddo
    
  case(KG_M2)
    do k = 1,size(hinterp_data,3)
       col_data(:,:) = col_data(:,:) + hinterp_data(:,:,k)
       hinterp_data(:,:,k) = hinterp_data(:,:,k)/ &
         ((clim_type%halflevs(k+1)-clim_type%halflevs(k))*p_fact)
    enddo
end select

found = .false.
do j = 1,size(climo_diag_name(:))
  if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
    found = .true.
    exit
  endif
enddo

if (found) then
  if (hinterp_id(j) > 0 ) then
       result = send_data(hinterp_id(j),col_data,Time)
  endif
endif


!++lwh
do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model surface pressure&
                             & is greater than climatology surface pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model top pressure&
                             & is less than climatology top pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
         case(INTERP_LINEAR_P)
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
!        case(INTERP_LOG)
      end select
   enddo
enddo

!--lwh

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k) = interp_data(:,:,k)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select

  endif !field_name
enddo !End of i loop
if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine interpolator_3D
!
!#######################################################################
!
!++lwh
subroutine interpolator_2D(clim_type, Time, interp_data, field_name, is, js, clim_units)
!
! Return 2-D field interpolated to model grid and time
!
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   Time        : The model time that you wish to interpolate to.
!   is, js      : The indices of the physics window.
!
! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name
!
type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)     :: field_name
type(time_type)       , intent(in)     :: Time
real, dimension(:,:),   intent(out)    :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: taum, taup
real :: hinterp_data(size(interp_data,1),size(interp_data,2),size(clim_type%levs(:)))
integer :: istart,iend,jstart,jend
logical :: result, found
logical :: found_field=.false.
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: j, i, n

!RG add option to go to 'no time axis' version if clim_type%TIME_FLAG==NOTIME i.e. data has no time, despite input received
if (clim_type%TIME_FLAG .eq. NOTIME) then
	call interpolator_2D_no_time_axis(clim_type, interp_data, field_name, is,js, clim_units)
	return
endif
	
if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "interpolator_2D : You must call interpolator_init before calling interpolator")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
!++lwh
  if ( field_name == clim_type%field_name(i) ) then
!--lwh

    found_field=.true.

    if(present(clim_units)) then
      call mpp_get_atts(clim_type%field_type(i),units=clim_units)
      clim_units = chomp(clim_units)
    endif


!----------------------------------------------------------------------
!   skip the time interpolation portion of this routine if subroutine
!   obtain_interpolator_time_slices has already been called on this
!   stewp for this interpolate_type variable.
!----------------------------------------------------------------------

if ( .not. clim_type%separate_time_vary_calc) then     
!   print *, 'TIME INTERPOLATION NOT SEPARATED 2d--',  &
!                                   trim(clim_type%file_name), mpp_pe()
    if (clim_type%climatological_year) then
!++lwh
       if (size(clim_type%time_slice) > 1) then
          call time_interp(Time, clim_type%time_slice, clim_type%tweight, taum, taup, modtime=YEAR )
       else
          taum = 1
          taup = 1
          clim_type%tweight = 0.
       end if
!--lwh
    else
       call time_interp(Time, clim_type%time_slice, clim_type%tweight, taum, taup )
    endif

! If the climatology file has seasonal, a split time-line or has all the data 
! read in then enter this loop.
! 
    if(clim_type%TIME_FLAG .ne. LINEAR .or. read_all_on_init) then
      clim_type%itaum=taum
      clim_type%itaup=taup
    endif

!    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
!      ! Check if delta-time is greater than delta of first two climatology time-slices.
!      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
!           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
!      ! The difference between the model time and the last climatology time-slice previous to the model time.
!      ! We need 2 time levels. Check we have the correct data.
!        itaum=0
!        itaup=0
!      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute 
!      ! the climatology year into the appropriate place.
!      
!        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
!        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)
!        clim_datem = set_date(climyear, modmonth, modday, modhour, modminute, modsecond)
!        call time_interp(clim_datem, clim_type%time_slice, tweight1, taum1, taup1 )
!
!
!        call get_date(clim_type%time_slice(taup), climyear, climmonth, climday, climhour, climminute, climsecond)
!        clim_datep = set_date(climyear, modmonth, modday, modhour, modminute, modsecond)
!
!
!      endif
!
!    endif
    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
      ! The difference between the model time and the last climatology time-slice previous to the model time.
      ! We need 2 time levels.
        clim_type%itaum=0
        clim_type%itaup=0
      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute 
      ! the climatology year into the appropriate place.

     
      ! We need to get the previous months data for the climatology year before 
      ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:,:),2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1,m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m 
        enddo
        do m = 1,12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m,climatology), year1, month1, day, hour, minute, second)
          if ( month1 == modmonth ) then
!RSHBUGFX   if ( modday <= day ) then 
            if ( modday <  day ) then 
              indexm = m-1 ; indexp = m
            else
              indexm = m ; indexp = m+1
            endif
          endif
        
        enddo
        if ( indexm == 0 ) then 
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        endif
        call get_date(clim_type%time_slice(indexm+(climatology-1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if ( indexp == 13 ) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        endif
        call get_date(clim_type%time_slice(indexp+(climatology-1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)
        
        call time_interp(Time, month, clim_type%tweight3, taum, taup ) ! tweight3 is the time weight between the months.

        month(1) = clim_type%time_slice(indexm+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexm+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, clim_type%tweight1, taum, taup ) !tweight1 is the time weight between the climatology years.

        month(1) = clim_type%time_slice(indexp+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexp+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, clim_type%tweight2, taum, taup ) !tweight1 is the time weight between the climatology years.



        if (indexm == clim_type%indexm(i) .and.  &
          indexp == clim_type%indexp(i) .and. &
          climatology == clim_type%climatology(i)) then
        else
          clim_type%indexm(i) = indexm
          clim_type%indexp(i) = indexp
          clim_type%climatology(i) = climatology
          call read_data(clim_type,clim_type%field_type(i),  &
            clim_type%pmon_pyear(:,:,:,i),  &
            clim_type%indexm(i)+(clim_type%climatology(i)-1)*12,i,Time)
! Read the data for the next month in the previous climatology.
          call read_data(clim_type,clim_type%field_type(i),  &
            clim_type%nmon_pyear(:,:,:,i),   &
            clim_type%indexp(i)+(clim_type%climatology(i)-1)*12,i,Time)
          call read_data(clim_type,clim_type%field_type(i),   &
            clim_type%pmon_nyear(:,:,:,i),  &
            clim_type%indexm(i)+clim_type%climatology(i)*12,i,Time)
          call read_data(clim_type,clim_type%field_type(i),  &
            clim_type%nmon_nyear(:,:,:,i),  &
            clim_type%indexp(i)+clim_type%climatology(i)*12,i,Time)
        endif




      else ! We are within a climatology data set
        
        if (taum /= clim_type%time_init(i,1) .or. &
            taup /= clim_type%time_init(i,2) ) then
 
          call read_data(clim_type,clim_type%field_type(i), clim_type%pmon_pyear(:,:,:,i), taum,i,Time)
! Read the data for the next month in the previous climatology.
          call read_data(clim_type,clim_type%field_type(i), clim_type%nmon_pyear(:,:,:,i), taup,i,Time)
!RSHbug   clim_type%pmon_nyear = 0.0
!RSHbug   clim_type%nmon_nyear = 0.0

!         clim_type%pmon_nyear(:,:,:,i) = 0.0
!         clim_type%nmon_nyear(:,:,:,i) = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
          clim_type%indexm(i) = 0       
          clim_type%indexp(i) = 0        
          clim_type%climatology(i) = 0             


          clim_type%time_init(i,1) = taum
          clim_type%time_init(i,2) = taup
        endif
!       clim_type%tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        clim_type%tweight1 = 0.0 ; clim_type%tweight2 = 0.0
        clim_type%tweight3 = clim_type%tweight                                          
      endif

    endif ! (BILINEAR)

    if(clim_type%TIME_FLAG .eq. LINEAR .and. &
        (.not. read_all_on_init) ) then
! We need 2 time levels. Check we have the correct data.
      clim_type%itaum=0
      clim_type%itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(i,n) .eq. taum ) clim_type%itaum = n
        if (clim_type%time_init(i,n) .eq. taup ) clim_type%itaup = n
      enddo

      if (clim_type%itaum.eq.0 .and. clim_type%itaup.eq.0) then
      !Neither time is set so we need to read 2 time slices.
      !Set up 
      ! field(:,:,:,1) as the previous time slice.
      ! field(:,:,:,2) as the next time slice.
        call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,1,i), taum,i,Time)
          clim_type%time_init(i,1) = taum
          clim_type%itaum = 1
        call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,2,i), taup,i,Time)
          clim_type%time_init(i,2) = taup
          clim_type%itaup = 2
      endif ! clim_type%itaum.eq.clim_type%itaup.eq.0
      if (clim_type%itaum.eq.0 .and. clim_type%itaup.ne.0) then
      ! Can't think of a situation where we would have the next time level but not the previous.
        call mpp_error(FATAL,'interpolator_2D : No data from the previous climatology time but we have&
                            & the next time. How did this happen?')
      endif
      if (clim_type%itaum.ne.0 .and. clim_type%itaup.eq.0) then
      !We have the previous time step but not the next time step data
        clim_type%itaup = 1
        if (clim_type%itaum .eq. 1 ) clim_type%itaup = 2
        call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,clim_type%itaup,i), taup,i, Time)
        clim_type%time_init(i,clim_type%itaup)=taup
      endif
    endif! TIME_FLAG .eq. LINEAR .and. (.not. read_all_on_init)

  endif ! (.not. separate_time_vary_calc)



select case(clim_type%TIME_FLAG)
  case (LINEAR)
    hinterp_data = (1-clim_type%tweight)*clim_type%data(istart:iend,jstart:jend,:,clim_type%itaum,i) &
                     + clim_type%tweight*clim_type%data(istart:iend,jstart:jend,:,clim_type%itaup,i)
! case (SEASONAL)
! Do sine fit to data at this point
  case (BILINEAR)
    hinterp_data = &
    (1-clim_type%tweight1)  * (1-clim_type%tweight3) * clim_type%pmon_pyear(istart:iend,jstart:jend,:,i) + &
    (1-clim_type%tweight2)  *    clim_type%tweight3  * clim_type%nmon_pyear(istart:iend,jstart:jend,:,i) + &
         clim_type%tweight1 * (1-clim_type%tweight3) * clim_type%pmon_nyear(istart:iend,jstart:jend,:,i) + &
         clim_type%tweight2 *     clim_type%tweight3 * clim_type%nmon_nyear(istart:iend,jstart:jend,:,i)

end select

found = .false.
do j = 1,size(climo_diag_name(:))
  if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
    found = .true.
    exit
  endif
enddo

if (found) then
  if (hinterp_id(j) > 0 ) then
       result = send_data(hinterp_id(j),hinterp_data,Time)
  endif
endif

  interp_data(:,:) = hinterp_data(:,:,1)

  endif !field_name
enddo !End of i loop

if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine interpolator_2D
!--lwh
!
!#######################################################################

subroutine interpolator_4D_no_time_axis(clim_type, phalf, interp_data, field_name, is,js, clim_units)

! Return 4-D field interpolated to model grid

! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init

! INTENT IN
!   field_name  : The name of a field that you wish to interpolate.
!                 all variables within this interpolate_type variable
!                 will be interpolated on this call. field_name may
!                 be any one of the variables.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.

! INTENT OUT
!   interp_data : The model fields
!   clim_units  : The units of field_name

type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
real, dimension(:,:,:), intent(in)  :: phalf
real, dimension(:,:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: ilon
real :: hinterp_data(size(interp_data,1),size(interp_data,2),size(clim_type%levs(:)),size(clim_type%field_name(:)))
real :: p_fact(size(interp_data,1),size(interp_data,2))
real :: pclim(size(clim_type%halflevs(:)))
integer :: istart,iend,jstart,jend
logical :: result
logical :: found_field=.false.
integer :: i, j, k, n

if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "interpolator_4D_no_time_axis : You must call interpolator_init before calling interpolator")

do n=2,size(clim_type%field_name(:))
  if (clim_type%vert_interp(n) /= clim_type%vert_interp(n-1) .or. &
   clim_type%out_of_bounds(n) /= clim_type%out_of_bounds(n-1)) then
    if (mpp_pe() == mpp_root_pe() ) then
      print *, 'processing file ' // trim(clim_type%file_name)
    endif
    call mpp_error (FATAL, 'interpolator_mod: &
            &cannot use 4D interface to interpolator for this file')
  endif
end do

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
  if ( field_name == clim_type%field_name(i) ) then
    found_field=.true.
    exit 
  endif
end do
i = 1

if(present(clim_units)) then
  call mpp_get_atts(clim_type%field_type(i),units=clim_units)
  clim_units = chomp(clim_units)
endif

do n=1, size(clim_type%field_name(:))
  hinterp_data(:,:,:,n) = clim_type%data(istart:iend,jstart:jend,:,1,n)
end do
    
select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select

    do i= 1, size(clim_type%field_name(:))
      select case(clim_type%mr(i))
      case(KG_M2)
        do k = 1,size(hinterp_data,3)
          hinterp_data(:,:,k,i) = hinterp_data(:,:,k,i)/((clim_type%halflevs(k+1)-clim_type%halflevs(k))*p_fact)
        enddo
      end select
    enddo

   i = 1

do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model surface pressure&
                             & is greater than surface pressure of input data for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model top pressure&
                             & is less than top pressure of input data for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,:),interp_data(ilon,j,:,:))
         case(INTERP_LINEAR_P)
          do n=1, size(clim_type%field_name(:))
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,n),interp_data(ilon,j,:,n))
          end do
      end select
   enddo
enddo

     do i= 1, size(clim_type%field_name(:))

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k,i) = interp_data(:,:,k,i)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select

     end do

if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine interpolator_4D_no_time_axis

!#######################################################################

subroutine interpolator_3D_no_time_axis(clim_type, phalf, interp_data, field_name, is,js, clim_units)

! Return 3-D field interpolated to model grid

! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init

! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.

! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name

type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
real, dimension(:,:,:), intent(in)  :: phalf
real, dimension(:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
real :: tweight, tweight1, tweight2, tweight3
integer :: taum, taup, ilon
real :: hinterp_data(size(interp_data,1),size(interp_data,2),size(clim_type%levs(:)))
real :: p_fact(size(interp_data,1),size(interp_data,2))
real :: pclim(size(clim_type%halflevs(:)))
integer :: istart,iend,jstart,jend
logical :: result
logical :: found_field=.false.
integer :: i, j, k, n

if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "interpolator_3D_no_time_axis : You must call interpolator_init before calling interpolator")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
  if ( field_name == clim_type%field_name(i) ) then
    found_field=.true.
    if(present(clim_units)) then
      call mpp_get_atts(clim_type%field_type(i),units=clim_units)
      clim_units = chomp(clim_units)
    endif

    hinterp_data = clim_type%data(istart:iend,jstart:jend,:,1,i)

select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(hinterp_data,3)
       hinterp_data(:,:,k) = hinterp_data(:,:,k)/((clim_type%halflevs(k+1)-clim_type%halflevs(k))*p_fact)
    enddo
end select

do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model surface pressure&
                             & is greater than climatology surface pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model top pressure&
                             & is less than climatology top pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
         case(INTERP_LINEAR_P)
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
      end select
   enddo
enddo

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k) = interp_data(:,:,k)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select

  endif !field_name
enddo !End of i loop
if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine interpolator_3D_no_time_axis

!#######################################################################

subroutine interpolator_2D_no_time_axis(clim_type, interp_data, field_name, is, js, clim_units)

! Return 2-D field interpolated to model grid

! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init

! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   is, js      : The indices of the physics window.

! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name

type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)     :: field_name
real, dimension(:,:),   intent(out)    :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
real :: tweight, tweight1, tweight2, tweight3
integer :: taum, taup, ilon
real :: hinterp_data(size(interp_data,1),size(interp_data,2),size(clim_type%levs(:)))
real :: p_fact(size(interp_data,1),size(interp_data,2))
integer :: istart,iend,jstart,jend
logical :: result
logical :: found_field=.false.
integer :: j, k, i, n

if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "interpolator_2D_no_time_axis : You must call interpolator_init before calling interpolator")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
  if ( field_name == clim_type%field_name(i) ) then

    found_field=.true.

    if(present(clim_units)) then
      call mpp_get_atts(clim_type%field_type(i),units=clim_units)
      clim_units = chomp(clim_units)
    endif

    hinterp_data = clim_type%data(istart:iend,jstart:jend,:,1,i)

    interp_data(:,:) = hinterp_data(:,:,1)

  endif !field_name
enddo !End of i loop

if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif

end subroutine interpolator_2D_no_time_axis

!#######################################################################
!
subroutine interpolator_end(clim_type)
! Subroutine to deallocate the interpolate type clim_type.
!
! INTENT INOUT
!  clim_type : allocate type whose components will be deallocated.
!
type(interpolate_type), intent(inout) :: clim_type
integer :: logunit

logunit=stdlog()
if ( mpp_pe() == mpp_root_pe() ) then
   write (logunit,'(/,(a))') 'Exiting interpolator, have a nice day ...'
end if

if (associated (clim_type%lat     )) deallocate(clim_type%lat)
if (associated (clim_type%lon     )) deallocate(clim_type%lon)
if (associated (clim_type%latb    )) deallocate(clim_type%latb)
if (associated (clim_type%lonb    )) deallocate(clim_type%lonb)
if (associated (clim_type%levs    )) deallocate(clim_type%levs)
if (associated (clim_type%halflevs)) deallocate(clim_type%halflevs) 
call horiz_interp_del(clim_type%interph)
if (associated (clim_type%time_slice)) deallocate(clim_type%time_slice)
if (associated (clim_type%field_type)) deallocate(clim_type%field_type)
if (associated (clim_type%field_name)) deallocate(clim_type%field_name)
if (associated (clim_type%time_init )) deallocate(clim_type%time_init)
if (associated (clim_type%mr        )) deallocate(clim_type%mr)
if (associated (clim_type%data)) then
  deallocate(clim_type%data)
endif
if (associated (clim_type%pmon_pyear)) then
  deallocate(clim_type%pmon_pyear)
  deallocate(clim_type%pmon_nyear)
  deallocate(clim_type%nmon_nyear)
  deallocate(clim_type%nmon_pyear)
endif

!! RSH mod   
if(  .not. ((clim_type%TIME_FLAG .eq. LINEAR  .and.    &
!     read_all_on_init)) .or. clim_type%TIME_FLAG .eq. BILINEAR  ) then
      read_all_on_init).or. ( clim_type%TIME_FLAG .eq. NOTIME ))  ) then
 call mpp_close(clim_type%unit)
endif


module_is_initialized = .false.

end subroutine interpolator_end
!
!#######################################################################
!
subroutine read_data(clim_type,src_field, hdata, nt,i, Time)
!
!  INTENT IN
!    clim_type : The interpolate type which contains the data 
!    src_field : The field type 
!    nt        : The index of the time slice of the climatology that you wish to read.
!    i         : The index of the field name that you are trying to read. (optional)
!    Time      : The model time. Used for diagnostic purposes only. (optional)
!
!  INTENT OUT
!
!    hdata     : The horizontally interpolated climatology field. This 
!                field will still be on the climatology vertical grid.
!
type(interpolate_type)   , intent(in)  :: clim_type
type(fieldtype)          , intent(in)  :: src_field
integer                  , intent(in)  :: nt
real                     , intent(out) :: hdata(:,:,:)
integer        , optional, intent(in)  :: i
type(time_type), optional, intent(in)  :: Time

integer   :: k, km
! sjs
real, allocatable :: climdata(:,:,:), climdata2(:,:,:)

      allocate(climdata(size(clim_type%lon(:)),size(clim_type%lat(:)), &
                        size(clim_type%levs(:))))

      call mpp_read(clim_type%unit,src_field, climdata,nt)

!  if vertical index increases upward, flip the data so that lowest
!  pressure level data is at index 1, rather than the highest pressure
!  level data. the indices themselves were previously flipped.
      if (clim_type%vertical_indices == INCREASING_UPWARD) then
        allocate(climdata2(size(clim_type%lon(:)),   &
                           size(clim_type%lat(:)), &
                           size(clim_type%levs(:))))
        km = size(clim_type%levs(:))
        do k=1, km                      
          climdata2(:,:,k) = climdata(:,:,km+1-k)
        end do
        climdata = climdata2
        deallocate (climdata2)
      endif

      call horiz_interp(clim_type%interph, climdata, hdata)
      if (clim_diag_initialized) &
        call diag_read_data(clim_type,climdata,i, Time)
      deallocate(climdata)


end subroutine read_data

!#######################################################################

subroutine read_data_no_time_axis(clim_type,src_field, hdata, i)

!  INTENT IN
!    clim_type : The interpolate type which contains the data 
!    src_field : The field type 
!    i         : The index of the field name that you are trying to read. (optional)

!  INTENT OUT

!    hdata     : The horizontally interpolated climatology field. This 
!                field will still be on the climatology vertical grid.

type(interpolate_type)   , intent(in)  :: clim_type
type(fieldtype)          , intent(in)  :: src_field
real                     , intent(out) :: hdata(:,:,:)
integer        , optional, intent(in)  :: i

integer   :: k, km
! sjs
real, allocatable :: climdata(:,:,:), climdata2(:,:,:)

      allocate(climdata(size(clim_type%lon(:)),size(clim_type%lat(:)), size(clim_type%levs(:))))

      call mpp_read(clim_type%unit,src_field, climdata)

!  if vertical index increases upward, flip the data so that lowest
!  pressure level data is at index 1, rather than the highest pressure
!  level data. the indices themselves were previously flipped.
      if (clim_type%vertical_indices == INCREASING_UPWARD) then
        allocate(climdata2(size(clim_type%lon(:)),   &
                           size(clim_type%lat(:)), &
                           size(clim_type%levs(:))))
        km = size(clim_type%levs(:))
        do k=1, km                      
          climdata2(:,:,k) = climdata(:,:,km+1-k)
        end do
        climdata = climdata2
        deallocate (climdata2)
      endif

      call horiz_interp(clim_type%interph, climdata, hdata)
      deallocate(climdata)

end subroutine read_data_no_time_axis

!#######################################################################
!
subroutine diag_read_data(clim_type,model_data, i, Time)
!
! A routine to diagnose the data read in by read_data
!
!  INTENT IN
!    clim_type  : The interpolate type.
!    model_data : The data read in from file that is being diagnosed.
!    i          : The index of the field name that you are diagnosing.
!    Time       : The model time
!
type(interpolate_type), intent(in) :: clim_type
real                  , intent(in) :: model_data(:,:,:)
integer               , intent(in) :: i
type(time_type)       , intent(in) :: Time

integer :: j,k
real :: col_data(size(model_data,1),size(model_data,2))
logical :: result, found


found = .false.
do j = 1,size(climo_diag_name(:))
  if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
      found = .true.
      exit
  endif
enddo

if(found) then
  if(climo_diag_id(j)>0) then
  col_data(:,:)=0.0
    do k=1,size(model_data,3)
      col_data(:,:) = col_data(:,:) + &
        model_data(:,:,k)* &
        (clim_type%halflevs(k+1)-clim_type%halflevs(k))/grav
    enddo
    result = send_data(climo_diag_id(j),col_data(clim_type%is:clim_type%ie,clim_type%js:clim_type%je),Time)
  endif
endif

end subroutine diag_read_data
!
!#######################################################################
!
!++lwh
subroutine query_interpolator( clim_type, nfields, field_names )
!
! Query an interpolate_type variable to find the number of fields and field names. 
!
type(interpolate_type), intent(in)                    :: clim_type
integer, intent(out), optional                        :: nfields
character(len=*), dimension(:), intent(out), optional :: field_names

if( present( nfields ) )     nfields     = SIZE( clim_type%field_name(:) )
if( present( field_names ) ) field_names = clim_type%field_name

end subroutine query_interpolator
!--lwh
!
!#######################################################################
!
function chomp(string)
!
! A function to remove CHAR(0) from the end of strings read from NetCDF files.
!
character(len=*), intent(in) :: string
character(len=64) :: chomp

integer :: len

len = len_trim(string)
if (string(len:len) == CHAR(0)) len = len -1

chomp = string(:len)

end function chomp
!
!#################################################################
!
 subroutine interp_weighted_scalar_2D (grdin, grdout, datin, datout )
real, intent(in),  dimension(:) :: grdin, grdout
real, intent(in),  dimension(:,:) :: datin
real, intent(out), dimension(:,:) :: datout

integer :: j, k, n

if (size(grdin(:)).ne. (size(datin,1)+1)) &
 call mpp_error(FATAL,'interp_weighted_scalar : input data and pressure do not have the same number of levels')
if (size(grdout(:)).ne. (size(datout,1 )+1)) &
 call mpp_error(FATAL,'interp_weighted_scalar : output data and pressure do not have the same number of levels')

  do k = 1, size(datout,1 )
   datout(k,:) = 0.0

     do j = 1, size(datin,1 )

        if ( grdin(j)   <= grdout(k) .and. &
             grdin(j+1) >= grdout(k) .and. &
             grdin(j+1) <= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdin(j+1)-grdout(k))
          end do

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j)   <= grdout(k+1) .and. &
                  grdin(j+1) >= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdout(k+1)-grdin(j))
          end do

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j+1) <= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdin(j+1)-grdin(j))
          end do

        else if ( grdin(j)   <= grdout(k)   .and. &
                  grdin(j+1) >= grdout(k+1) ) then

          do n= 1, size(datin,2)
          datout(k,n) = datout(k,n) + datin(j,n)*(grdout(k+1)-grdout(k))

          end do
        endif

     enddo

     do n= 1, size(datin,2)
       datout(k,n) = datout(k,n)/(grdout(k+1)-grdout(k))
     end do

  enddo

end subroutine interp_weighted_scalar_2D


!---------------------------------------------------------------------
 
 subroutine interp_weighted_scalar_1D (grdin, grdout, datin, datout )
real, intent(in),  dimension(:) :: grdin, grdout, datin
real, intent(out), dimension(:) :: datout

integer :: j, k

if (size(grdin(:)).ne. (size(datin(:))+1)) &
 call mpp_error(FATAL,'interp_weighted_scalar : input data and pressure do not have the same number of levels')
if (size(grdout(:)).ne. (size(datout(:))+1)) &
 call  mpp_error(FATAL,'interp_weighted_scalar : output data and pressure do not have the same number of levels')

  do k = 1, size(datout(:))
   datout(k) = 0.0

     do j = 1, size(datin(:))

        if ( grdin(j)   <= grdout(k) .and. &
             grdin(j+1) >= grdout(k) .and. &
             grdin(j+1) <= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdin(j+1)-grdout(k))

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j)   <= grdout(k+1) .and. &
                  grdin(j+1) >= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdout(k+1)-grdin(j))

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j+1) <= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdin(j+1)-grdin(j))

        else if ( grdin(j)   <= grdout(k)   .and. &
                  grdin(j+1) >= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdout(k+1)-grdout(k))

        endif

     enddo

     datout(k) = datout(k)/(grdout(k+1)-grdout(k))

  enddo

end subroutine interp_weighted_scalar_1D
!
!#################################################################
!
subroutine interp_linear ( grdin, grdout, datin, datout )
real, intent(in),  dimension(:) :: grdin, grdout, datin
real, intent(out), dimension(:) :: datout

integer :: j, k, n
real    :: wt


if (size(grdin(:)).ne. (size(datin(:))+1)) &
 call mpp_error(FATAL,'interp_linear : input data and pressure do not have the same number of levels')
if (size(grdout(:)).ne. (size(datout(:))+1)) &
 call mpp_error(FATAL,'interp_linear : output data and pressure do not have the same number of levels')


  n = size(grdin(:))

  do k= 1, size(datout(:))

   ! ascending grid values
     if (grdin(1) < grdin(n)) then
         do j = 2, size(grdin(:))-1
           if (grdout(k) <= grdin(j)) exit
         enddo
   ! descending grid values
     else
         do j = size(grdin(:)), 3, -1
           if (grdout(k) <= grdin(j-1)) exit
         enddo
     endif

   ! linear interpolation
     wt = (grdout(k)-grdin(j-1)) / (grdin(j)-grdin(j-1))
!print '(a,2i3,4f6.1)', 'k,j=',k,j,grdout(k),grdin(j-1),grdin(j),wt
   ! constant value extrapolation
   ! wt = min(max(wt,0.),1.)

     datout(k) = (1.-wt)*datin(j-1) + wt*datin(j)
     
  enddo

end subroutine interp_linear
!
!########################################################################

end module interpolator_mod
!
!#######################################################################
!
#ifdef test_interp
program test

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#endif

use mpp_mod
use mpp_io_mod
use mpp_domains_mod
use fms_mod
use time_manager_mod
use diag_manager_mod!, only : diag_axis_init, file_exist, MPP_NPES, &
                    !  MPP_PE, REGISTER_DIAG_FIELD, SEND_DATA, SET_DATE,&
                    !  SET_TIME

use interpolator_mod
!use sulfate_mod
!use ozone_mod
use constants_mod, only : grav, constants_init, PI
use time_interp_mod, only : time_interp_init

implicit none
integer, parameter :: nsteps_per_day = 8, ndays = 16
real, parameter :: delt = 1.0/nsteps_per_day
! integer, parameter :: nxd = 144, nyd = 90, ntsteps = 240, two_delt = 2*delt
integer, parameter :: nxd = 20, nyd = 40, ntsteps = nsteps_per_day*ndays, two_delt = 2*delt
integer :: delt_days, delt_secs
integer, parameter :: max_fields = 20 ! maximum number of fields to be interpolated

integer :: i,k,n,level
integer :: unit, io_status
integer :: ndivs
integer :: jscomp, jecomp, iscomp, iecomp, isd,ied,jsd,jed
integer :: numfields, domain_layout(2)
integer :: num_nbrs, nbins,axes(3), interp_diagnostic_id
integer :: column_diagnostic_id1, column_diagnostic_id(max_fields)

real ::  missing_value = -1.e10

character(len=1) :: dest_grid
character(len=128) :: src_file, file_out, title, units, colaer
logical :: vector_field=.false., result

type(axistype), allocatable, dimension(:)  :: axes_out, axes_src
type(axistype) :: time_axis
type(fieldtype), allocatable, dimension(:) :: fields
type(fieldtype) :: dest_field(max_fields), src_field(max_fields), field_geolon_t, &
     field_geolat_t, field_geolon_c, field_geolat_c
type(atttype), allocatable, dimension(:) :: global_atts
type(domain2d) :: domain
type(time_type) :: model_time

type(interpolate_type) :: o3, aerosol

real, dimension(:,:), allocatable :: col_data
real, dimension(:,:,:), allocatable :: model_data, p_half, p_full
real, dimension(:), allocatable :: latb_mod(:,:),lonb_mod(:,:),lon_mod(:),lat_mod(:)
real :: dx,dy
real :: dtr,tpi
real :: p_bot,p_top,lambda
character(len=64) :: names(13)
data names(:) /"so4_anthro","so4_natural","organic_carbon","black_carbon","sea_salt",&
"anthro_dust_0.2","anthro_dust_0.8","anthro_dust_2.0","anthro_dust_8.0",&
"natural_dust_0.2","natural_dust_0.8","natural_dust_2.0","natural_dust_8.0"/

integer :: out_of_bounds(1)
data out_of_bounds / CONSTANT/!, CONSTANT/!, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, &
!ZERO, ZERO, ZERO, ZERO /

namelist /interpolator_nml/ src_file

! initialize communication modules

delt_days = INT(delt)
delt_secs = INT(delt*86400.0) - delt_days*86400.0

write(*,*) delt, delt_days,delt_secs

call mpp_init
call mpp_io_init
call mpp_domains_init
call set_calendar_type(JULIAN)
call diag_manager_init
call constants_init
call time_interp_init

level = 18
tpi = 2.0*PI !4.*acos(0.)
dtr = tpi/360.

src_file = 'src_file'  ! input file containing fields to be interpolated


model_time = set_date(1979,12,1,0,0,0)

!if (numfields.ne.2.and.vector_field) call mpp_error(FATAL,'2 components of vector field not specified')
!if (numfields.gt.1.and..not.vector_field) call mpp_error(FATAL,'only 1 scalar at a time')
!if (numfields .gt. max_fields) call mpp_error(FATAL,'max num fields exceeded')

!--------------------------------------------------------------------
! namelist input
!--------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=interpolator_nml, iostat=io)
  ierr = check_nml_error(io, 'interpolator_nml')
#else
  call mpp_open(unit, 'input.nml',  action=MPP_RDONLY, form=MPP_ASCII)
  read  (unit, interpolator_nml,iostat=io_status)
  if (io_status .gt. 0) then
    call mpp_error(FATAL,'=>Error reading interpolator_nml')
  endif
call mpp_close(unit)
#endif

! decompose model grid points
! mapping can get expensive so we distribute the task at this level

ndivs = mpp_npes()

call mpp_define_layout ((/1,nxd,1,nyd/), ndivs, domain_layout)
call mpp_define_domains((/1,nxd,1,nyd/),domain_layout, domain,xhalo=0,yhalo=0)  
call mpp_get_data_domain(domain,isd,ied,jsd,jed)
call mpp_get_compute_domain (domain, iscomp, iecomp, jscomp, jecomp)

allocate(lonb_mod(nxd+1,nyd+1),lon_mod(nxd))
allocate(latb_mod(nxd+1,nyd+1),lat_mod(nyd))
allocate(col_data(isd:ied,jsd:jed)) ; col_data = 0.0
allocate(p_half(isd:ied,jsd:jed,level+1),p_full(isd:ied,jsd:jed,level))
p_top = 1.0
p_bot = 101325.0 !Model level in Pa
lambda = -1.0*log(p_top/p_bot)/(level+1)

p_half(:,:,level+1) = p_bot
do i=level,1,-1
  p_half(:,:,i)=p_half(:,:,i+1)*exp(-1.0*lambda)
enddo
do i=1,level
  p_full(:,:,i)=(p_half(:,:,i+1)+p_half(:,:,i))/2.0
enddo

allocate(model_data(isd:ied,jsd:jed,level))

dx = 360./nxd
dy = 180./nyd
do i = 1,nxd+1
  lonb_mod(i,:) = (i-1)*dx 
enddo
do i = 1,nyd+1
  latb_mod(:,i) = -90. + (i-1)*dy 
enddo
do i=1,nxd
  lon_mod(i)=(lonb_mod(i+1,1)+lonb_mod(i,1))/2.0
enddo
do i=1,nyd
  lat_mod(i)=(latb_mod(1,i+1)+latb_mod(1,i))/2.0
enddo

lonb_mod = lonb_mod * dtr
latb_mod = latb_mod * dtr

   axes(1) = diag_axis_init('x',lon_mod,units='degrees',cart_name='x',domain2=domain)
   axes(2) = diag_axis_init('y',lat_mod,units='degrees',cart_name='y',domain2=domain)
   axes(3) = diag_axis_init('z',p_full(isd,jsd,:),units='mb',cart_name='z')

interp_diagnostic_id =  register_diag_field('interp','ozone',axes(1:3),model_time,&
                                'interpolated_ozone_clim', 'kg/kg', missing_value)      
column_diagnostic_id1 =  register_diag_field('interp','colozone',axes(1:2),model_time,&
                                'column_ozone_clim', 'kg/m2', missing_value)      

do i=1,size(names(:))
colaer = 'col'//trim(names(i))
column_diagnostic_id(i) =  register_diag_field('interp',colaer,axes(1:2),model_time,&
                                'column_aerosol_clim', 'kg/m2', missing_value)      
enddo


call ozone_init(o3,lonb_mod(isd:ied+1,jsd:jed+1), latb_mod(isd:ied+1,jsd:jed+1), axes, model_time, &
                data_out_of_bounds=out_of_bounds)
call init_clim_diag(o3, axes, model_time)
call sulfate_init(aerosol,lonb_mod(isd:ied+1,jsd:jed+1), latb_mod(isd:ied+1,jsd:jed+1), names, &
                data_out_of_bounds=(/CONSTANT/) )
call init_clim_diag(aerosol, axes, model_time)

do n=1,ntsteps
  if( mpp_pe() == mpp_root_pe() ) write(*,*) n

  call get_ozone(o3,model_time,p_half,model_data)

  if(interp_diagnostic_id>0) &
       result = send_data(interp_diagnostic_id,&
            model_data(iscomp:iecomp,jscomp:jecomp,:),model_time)

  if(column_diagnostic_id1>0) then

    col_data(iscomp:iecomp,jscomp:jecomp)=0.0
    do k=1,level
       col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
          model_data(iscomp:iecomp,jscomp:jecomp,k)* &
          (p_half(iscomp:iecomp,jscomp:jecomp,k+1)-p_half(iscomp:iecomp,jscomp:jecomp,k))/grav
    enddo
       result = send_data(column_diagnostic_id1,col_data(:,:),model_time)
  endif



  do i=1,size(names(:))

call get_anthro_sulfate(aerosol,model_time,p_half,names(i),model_data,clim_units=units)

    if(column_diagnostic_id(i)>0) then

      col_data(iscomp:iecomp,jscomp:jecomp)=0.0
      do k=1,level
        if (trim(units) .eq. 'kg/m^2') then
           col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
              model_data(iscomp:iecomp,jscomp:jecomp,k)
        else
           col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
              model_data(iscomp:iecomp,jscomp:jecomp,k)* &
              (p_half(iscomp:iecomp,jscomp:jecomp,k+1)-p_half(iscomp:iecomp,jscomp:jecomp,k))/grav
        endif
      enddo
      result = send_data(column_diagnostic_id(i),&
      col_data(iscomp:iecomp,jscomp:jecomp),model_time)
    endif

  enddo

   model_time = model_time + set_time(delt_secs,delt_days)      

   if (n.eq. ntsteps) call diag_manager_end(model_time)

enddo

call interpolator_end(aerosol)
call interpolator_end(o3)

deallocate(lonb_mod, lon_mod, latb_mod,lat_mod, col_data, p_half, p_full, model_data)

call mpp_exit

contains
!
!#######################################################################
!
subroutine sulfate_init(aerosol,lonb, latb, names, data_out_of_bounds, vert_interp, units)
type(interpolate_type), intent(inout)         :: aerosol
real,                   intent(in)            :: lonb(:,:),latb(:,:)
character(len=64),      intent(in)            :: names(:)
integer,                intent(in)            :: data_out_of_bounds(:) 
integer,                intent(in), optional  :: vert_interp(:)
character(len=*),       intent(out),optional  :: units(:)

if (.not. file_exist("INPUT/aerosol.climatology.nc") ) return
call interpolator_init( aerosol, "aerosol.climatology.nc", lonb, latb, &
                        data_names=names, data_out_of_bounds=data_out_of_bounds, &
                        vert_interp=vert_interp, clim_units=units )

end subroutine sulfate_init
!
!#######################################################################
!
subroutine get_anthro_sulfate( sulfate, model_time, p_half, name, model_data, is, js, clim_units )
type(interpolate_type), intent(inout) :: sulfate
type(time_type), intent(in) :: model_time
real, intent(in)           :: p_half(:,:,:)
character(len=*), intent(in) :: name
character(len=*), intent(out), optional :: clim_units
real, intent(out) :: model_data(:,:,:)
integer, intent(in), optional :: is,js

call interpolator( sulfate, model_time, p_half, model_data, name, is, js, clim_units)

end subroutine get_anthro_sulfate
!
!#######################################################################
!
subroutine ozone_init( o3, lonb, latb, axes, model_time, data_out_of_bounds, vert_interp )
real,                  intent(in)           :: lonb(:,:),latb(:,:)
integer,               intent(in)           :: axes(:)
type(time_type),       intent(in)           :: model_time
type(interpolate_type),intent(inout)        :: o3
integer,               intent(in)           :: data_out_of_bounds(:)
integer,               intent(in), optional :: vert_interp(:)

if (.not. file_exist("INPUT/o3.climatology.nc") ) return
call interpolator_init( o3, "o3.climatology.nc", lonb, latb, &
                        data_out_of_bounds=data_out_of_bounds, vert_interp=vert_interp )

end subroutine ozone_init
!
!#######################################################################
!
subroutine get_ozone( o3, model_time, p_half, model_data, is, js )
type(interpolate_type),intent(inout) :: o3
type(time_type), intent(in) :: model_time
real, intent(in)           :: p_half(:,:,:)
real, intent(out) :: model_data(:,:,:)
integer, intent(in), optional :: is,js

call interpolator( o3, model_time, p_half, model_data, "ozone", is, js)

end subroutine get_ozone

end program test

#endif
