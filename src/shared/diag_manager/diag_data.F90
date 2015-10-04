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

#include <fms_platform.h>

MODULE diag_data_mod
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
  !   Seth Underwood
  ! </CONTACT>
  
  ! <OVERVIEW>
  !   Type descriptions and global variables for the diag_manager modules.
  ! </OVERVIEW>

  ! <DESCRIPTION>
  !   Notation: 
  !   <DL>
  !     <DT>input field</DT>
  !     <DD>The data structure describing the field as
  !       registered by the model code.</DD>
  !
  !     <DT>output field</DT>
  !     <DD>The data structure describing the actual
  !       diagnostic output with requested frequency and
  !       other options.</DD>
  !   </DL>
  !
  !   Input fields, output fields, and output files are gathered in arrays called
  !   "input_fields", "output_fields", and "files", respectively. Indices in these
  !   arrays are used as pointers to create associations between various data
  !   structures.
  !
  !   Each input field associated with one or several output fields via array of
  !   indices output_fields; each output field points to the single "parent" input
  !   field with the input_field index, and to the output file with the output_file 
  !   index
  ! </DESCRIPTION>

  USE time_manager_mod, ONLY: time_type
  USE mpp_domains_mod,  ONLY: domain1d, domain2d
  USE mpp_io_mod,       ONLY: fieldtype
  USE fms_mod, ONLY: WARNING
#ifdef use_netCDF
  ! NF90_FILL_REAL has value of 9.9692099683868690e+36.
  USE netcdf, ONLY: NF_FILL_REAL => NF90_FILL_REAL
#endif

  IMPLICIT NONE

  PUBLIC


  ! <!-- PARAMETERS for diag_data.F90 -->
  ! <DATA NAME="MAX_FIELDS_PER_FILE" TYPE="INTEGER, PARAMETER" DEFAULT="300">
  !   Maximum number of fields per file.
  ! </DATA>
  ! <DATA NAME="MAX_OUT_PER_IN_FIELD" TYPE="INTEGER, PARAMETER" DEFAULT="150">
  !   Maximum number of output_fields per input_field.
  ! </DATA>
  ! <DATA NAME="DIAG_OTHER" TYPE="INTEGER, PARAMETER" DEFAULT="0" />
  ! <DATA NAME="DIAG_OCEAN" TYPE="INTEGER, PARAMETER" DEFAULT="1" />
  ! <DATA NAME="DIAG_ALL" TYPE="INTEGER, PARAMETER" DEFAULT="2" />
  ! <DATA NAME="VERY_LARGE_FILE_FREQ" TYPE="INTEGER, PARAMETER" DEFAULT="100000" />
  ! <DATA NAME="VERY_LARGE_AXIS_LENGTH" TYPE="INTEGER, PARAMETER" DEFAUTL="10000" />
  ! <DATA NAME="EVERY_TIME" TYPE="INTEGER, PARAMETER" DEFAULT="0" />
  ! <DATA NAME="END_OF_RUN" TYPE="INTEGER, PARAMETER" DEFAULT="-1" />
  ! <DATA NAME="DIAG_SECONDS" TYPE="INTEGER, PARAMETER" DEFAULT="1" />
  ! <DATA NAME="DIAG_MINUTES" TYPE="INTEGER, PARAMETER" DEFAULT="2" />
  ! <DATA NAME="DIAG_HOURS" TYPE="INTEGER, PARAMETER" DEFAULT="3" />
  ! <DATA NAME="DIAG_DAYS" TYPE="INTEGER, PARAMETER" DEFAULT="4" />
  ! <DATA NAME="DIAG_MONTHS" TYPE="INTEGER, PARAMETER" DEFAULT="5" />
  ! <DATA NAME="DIAG_YEARS" TYPE="INTEGER, PARAMETER" DEFAULT="6" />
  ! <DATA NAME="MAX_SUBAXES" TYPE="INTEGER, PARAMETER" DEFAULT="10" />
  ! <DATA NAME="CMOR_MISSING_VALUE" TYPE="REAL, PARAMETER" DEFAULT="1.0e20" />

  ! Specify storage limits for fixed size tables used for pointers, etc.
  INTEGER, PARAMETER :: MAX_FIELDS_PER_FILE = 300 !< Maximum number of fields per file.
  INTEGER, PARAMETER :: MAX_OUT_PER_IN_FIELD = 150 !< Maximum number of output_fields per input_field
  INTEGER, PARAMETER :: DIAG_OTHER = 0
  INTEGER, PARAMETER :: DIAG_OCEAN = 1
  INTEGER, PARAMETER :: DIAG_ALL   = 2
  INTEGER, PARAMETER :: VERY_LARGE_FILE_FREQ = 100000
  INTEGER, PARAMETER :: VERY_LARGE_AXIS_LENGTH = 10000
  INTEGER, PARAMETER :: EVERY_TIME =  0
  INTEGER, PARAMETER :: END_OF_RUN = -1
  INTEGER, PARAMETER :: DIAG_SECONDS = 1, DIAG_MINUTES = 2, DIAG_HOURS = 3
  INTEGER, PARAMETER :: DIAG_DAYS = 4, DIAG_MONTHS = 5, DIAG_YEARS = 6
  INTEGER, PARAMETER :: MAX_SUBAXES = 10
  REAL, PARAMETER :: CMOR_MISSING_VALUE = 1.0e20 !< CMOR standard missing value

  ! <TYPE NAME="diag_grid">
  !   <DESCRIPTION>
  !     Contains the coordinates of the local domain to output.
  !   </DESCRIPTION>
  !   <DATA NAME="start" TYPE="REAL, DIMENSION(3)">
  !     Start coordinates (Lat, Lon, Depth) of the local domain to output.
  !   </DATA>
  !   <DATA NAME="end" TYPE="REAL, DIMENSION(3)">
  !     End coordinates (Lat, Lon, Depth) of the local domain to output.
  !   </DATA>
  !   <DATA NAME="l_start_indx" TYPE="INTEGER, DIMENSION(3)">
  !     Start indices at each local PE.
  !   </DATA>
  !   <DATA NAME="l_end_indx" TYPE="INTEGER, DIMENSION(3)">
  !     End indices at each local PE.
  !   </DATA>
  !   <DATA NAME="subaxes" TYPE="INTEGER, DIMENSION(3)">
  !     ID returned from diag_subaxes_init of 3 subaces.
  !   </DATA>
  TYPE diag_grid
     REAL, DIMENSION(3) :: start, END ! start and end coordinates (lat,lon,depth) of local domain to output   
     INTEGER, DIMENSION(3) :: l_start_indx, l_end_indx ! start and end indices at each LOCAL PE
     INTEGER, DIMENSION(3) :: subaxes ! id returned from diag_subaxes_init of 3 subaxes
  END TYPE diag_grid
  ! </TYPE>
  
  ! <TYPE NAME="diag_fieldtype">
  !   <DESCRIPTION>
  !     Diagnostic field type
  !   </DESCRIPTION>
  !   <DATA NAME="Field" TYPE="TYPE(fieldtype)">
  !   </DATA>
  !   <DATA NAME="Domain" TYPE="TYPE(domain2d)">
  !   </DATA>
  !   <DATA NAME="miss" TYPE="REAL">
  !   </DATA>
  !   <DATA NAME="miss_pack" TYPE="REAL">
  !   </DATA>
  !   <DATA NAME="miss_present" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="miss_pack_present" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="tile_count" TYPE="INTEGER">
  !   </DATA>
  TYPE diag_fieldtype
     TYPE(fieldtype) :: Field
     TYPE(domain2d) :: Domain
     REAL :: miss, miss_pack
     LOGICAL :: miss_present, miss_pack_present
     INTEGER :: tile_count
  END TYPE diag_fieldtype
  ! </TYPE>

  ! <TYPE NAME="coord_type">
  !   <DESCRIPTION>
  !     Define the region for field output.
  !   </DESCRIPTION>
  !   <DATA NAME="xbegin" TYPE="REAL">
  !   </DATA>
  !   <DATA NAME="xend" TYPE="REAL">
  !   </DATA>
  !   <DATA NAME="ybegin" TYPE="REAL">
  !   </DATA>
  !   <DATA NAME="yend" TYPE="REAL">
  !   </DATA>
  !   <DATA NAME="zbegin" TYPE="REAL">
  !   </DATA>
  !   <DATA NAME="zend" TYPE="REAL">
  !   </DATA>
  TYPE coord_type
     REAL :: xbegin
     REAL :: xend
     REAL :: ybegin
     REAL :: yend
     REAL :: zbegin
     REAL :: zend
  END TYPE coord_type
  ! </TYPE>
  
  ! <TYPE NAME="file_type">
  !   <DESCRIPTION>
  !     Type to define the diagnostic files that will be written as defined by the diagnostic table.
  !   </DESCRIPTION>
  !   <DATA NAME="name" TYPE="CHARACTER(len=128)">
  !     Name of the output file.
  !   </DATA>
  !   <DATA NAME="long_name" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="fields" TYPE="INTEGER, dimension(max_fields_per_file)">
  !   </DATA>
  !   <DATA NAME="num_fields" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="output_freq" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="output_units" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="format" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="time_units" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="file_unit" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="bytes_written" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="time_axis_id" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="time_bounds_id" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="new_file_freq" TYPE="INTEGER">
  !     Frequency to create a new file.
  !   </DATA>
  !   <DATA NAME="new_file_freq_units" TYPE="INTEGER">
  !     Time units of new_file_freq ( days, hours, years, ...)
  !   </DATA>
  !   <DATA NAME="duration" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="duration_units" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="tile_count" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="local" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="last_flush" TYPE="TYPE(time_type)">
  !   </DATA>
  !   <DATA NAME="next_open" TYPE="TYPE(time_type)">
  !     Time to open next file.
  !   </DATA>
  !   <DATA NAME="start_time" TYPE="TYPE(time_type)">
  !     Time file opened
  !   </DATA>
  !   <DATA NAME="close_time" TYPE="TYPE(time_type)">
  !     Time file closed.  File does not allow data after close time
  !   </DATA>
  !   <DATA NAME="f_avg_start" TYPE="TYPE(diag_fieldtype)">
  !   </DATA>
  !   <DATA NAME="f_avg_end" TYPE="TYPE(diag_fieldtype)">
  !   </DATA>
  !   <DATA NAME="f_avg_nitems" TYPE="TYPE(diag_fieldtype)">
  !   </DATA>
  !   <DATA NAME="f_bounds" TYPE="TYPE(diag_fieldtype)">
  !   </DATA>
  TYPE file_type
     CHARACTER(len=128) :: name !< Name of the output file.
     CHARACTER(len=128) :: long_name
     INTEGER, DIMENSION(max_fields_per_file) :: fields
     INTEGER :: num_fields
     INTEGER :: output_freq
     INTEGER :: output_units
     INTEGER :: FORMAT
     INTEGER :: time_units
     INTEGER :: file_unit
     INTEGER :: bytes_written
     INTEGER :: time_axis_id, time_bounds_id
     INTEGER :: new_file_freq !< frequency to create new file
     INTEGER :: new_file_freq_units !< time units of new_file_freq (days, hours, years, ...)
     INTEGER :: duration
     INTEGER :: duration_units
     INTEGER :: tile_count
     LOGICAL :: local !< .TRUE. if fields are output in a region instead of global.
     TYPE(time_type) :: last_flush
     TYPE(time_type) :: next_open !< Time to open a new file.
     TYPE(time_type) :: start_time !< Time file opened.
     TYPE(time_type) :: close_time !< Time file closed.  File does not allow data after close time
     TYPE(diag_fieldtype):: f_avg_start, f_avg_end, f_avg_nitems, f_bounds
  END TYPE file_type
  ! </TYPE>  
  
  ! <TYPE NAME="input_field_type">
  !   <DESCRIPTION>
  !     Type to hold the input field description
  !   </DESCRIPTION>
  !   <DATA NAME="module_name" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="field_name" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="long_name" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="units" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="standard_name" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="interp_method" TYPE="CHARACTER(len=64)">
  !   </DATA>
  !   <DATA NAME="axes" TYPE="INTEGER, DIMENSION(3)">
  !   </DATA>
  !   <DATA NAME="num_axes" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="missing_value_present" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="range_present" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="missing_value" TYPE="REAL">
  !   </DATA>
  !   <DATA NAME="range" TYPE="REAL, DIMENSION(2)">
  !   </DATA>
  !   <DATA NAME="output_fields" TYPE="INTEGER, DIMENSION(max_out_per_in_field)">
  !   </DATA>
  !   <DATA NAME="num_output_fields" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="size" TYPE="INTEGER, DIMENSION(3)">
  !   </DATA>
  !   <DATA NAME="static" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="register" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="mask_variant" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="local" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="tile_count" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="local_coord" TYPE="TYPE(coord_type)">
  !   </DATA>
  TYPE input_field_type
     CHARACTER(len=128) :: module_name, field_name, long_name, units, standard_name
     CHARACTER(len=64) :: interp_method
     INTEGER, DIMENSION(3) :: axes
     INTEGER :: num_axes 
     LOGICAL :: missing_value_present, range_present
     REAL :: missing_value
     REAL, DIMENSION(2) :: range
     INTEGER, DIMENSION(max_out_per_in_field) :: output_fields
     INTEGER :: num_output_fields
     INTEGER, DIMENSION(3) :: size
     LOGICAL :: static, register, mask_variant, local
     INTEGER :: numthreads
     INTEGER :: tile_count
     TYPE(coord_type) :: local_coord
     TYPE(time_type)  :: time
     LOGICAL :: issued_mask_ignore_warning
  END TYPE input_field_type
  ! </TYPE>

  ! <TYPE NAME="output_field_type">
  !   <DESCRIPTION>
  !     Type to hold the output field description.
  !   </DESCRIPTION>
  !   <DATA NAME="input_field" TYPE="INTEGER">
  !     Index of the corresponding input field in the table
  !   </DATA>
  !   <DATA NAME="output_file" TYPE="INTEGER">
  !     Index of the output file in the table
  !   </DATA>
  !   <DATA NAME="output_name" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="static" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="time_max" TYPE="LOGICAL">
  !     .TRUE. if the output field is maximum over time interval
  !   </DATA>
  !   <DATA NAME="time_min" TYPE="LOGICAL">
  !     .TRUE. if the output field is minimum over time interval
  !   </DATA>
  !   <DATA NAME="time_average" TYPE="LOGICAL">
  !     .TRUE. if the output field is averaged over time interval.
  !   </DATA>
  !   <DATA NAME="time_ops" TYPE="LOGICAL">
  !     .TRUE. if any of time_min, time_max, or time_average is true
  !   </DATA>
  !   <DATA NAME="pack" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="time_method" TYPE="CHARACTER(len=50)">
  !     Time method field from the input file
  !   </DATA>
  !   <DATA NAME="buffer" TYPE="REAL, _ALLOCATABLE, DIMENSION(:,:,:,:)" DEFAULT="_NULL">
  !     Coordinates of buffer are (x, y, z, time-of-day)
  !   </DATA>
  !   <DATA NAME="counter" TYPE="REAL, _ALLOCATABLE, DIMENSION(:,:,:,:)" DEFAULT="_NULL">
  !     Coordinates of buffer are (x, y, z, time-of-day)
  !   </DATA>
  !   <DATA NAME="count_0d" TYPE="REAL, _ALLOCATABLE, DIMENSION(:)">
  !   </DATA>
  !   <DATA NAME="num_elements" TYPE="REAL, _ALLOCATABLE, DIMENSION(:)">
  !   </DATA>
  !   <DATA NAME="last_output" TYPE="TYPE(time_type)">
  !   </DATA>
  !   <DATA NAME="next_output" TYPE="TYPE(time_type)">
  !   </DATA>
  !   <DATA NAME="next_next_output" TYPE="TYPE(time_type)">
  !   </DATA>
  !   <DATA NAME="f_type" TYPE="TYPE(diag_fieldtype)">
  !   </DATA>
  !   <DATA NAME="axes" TYPE="INTEGER, DIMENSION(4)">
  !   </DATA>
  !   <DATA NAME="num_axes" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="total_elements" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="region_elements" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="n_diurnal_samples" TYPE="INTEGER">
  !     Number of diurnal sample intervals, 1 or more
  !   </DATA>
  !   <DATA NAME="output_grid" TYPE="TYPE(diag_grid)">
  !   </DATA>
  !   <DATA NAME="local_output" TYPE="LOGICAL">
  !     .TRUE. if this field is written out on a region and not globally.
  !   </DATA>
  !   <DATA NAME="need_compute" TYPE="LOGICAL">
  !     .TRUE. if this field is written out on a region, not global.
  !   </DATA>
  !   <DATA NAME="phys_window" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="written_once" TYPE="LOGICAL">
  !   </DATA>
  !   <DATA NAME="reduced_k_range" TYPE="LOGICAL">
  !     .TRUE. if dealing with vertical sub-level output.
  !   </DATA>
  !   <DATA NAME="imin" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="imax" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="jmin" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="jmax" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="kmin" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="kmax" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="Time_of_prev_field_data" TYPE="TYPE(time_type)">
  !   </DATA>
  TYPE output_field_type
     INTEGER :: input_field ! index of the corresponding input field in the table
     INTEGER :: output_file ! index of the output file in the table
     CHARACTER(len=128) :: output_name
     LOGICAL :: time_average ! true if the output field is averaged over time interval
     LOGICAL :: static
     LOGICAL :: time_max ! true if the output field is maximum over time interval
     LOGICAL :: time_min ! true if the output field is minimum over time interval
     LOGICAL :: time_ops ! true if any of time_min, time_max, or time_average is true
     INTEGER  :: pack
     CHARACTER(len=50) :: time_method ! time method field from the input file 
     ! coordianes of the buffer and counter are (x, y, z, time-of-day)
     REAL, _ALLOCATABLE, DIMENSION(:,:,:,:) :: buffer _NULL
     REAL, _ALLOCATABLE, DIMENSION(:,:,:,:) :: counter _NULL
     ! the following two counters are used in time-averaging for some 
     ! combination of the field options. Their size is the length of the 
     ! diurnal axis; the counters must be tracked separately for each of
     ! the diurnal interval, becaus the number of time slices accumulated
     ! in each can be different, depending on time step and the number of
     ! diurnal samples.
     REAL, _ALLOCATABLE, DIMENSION(:)  :: count_0d
     INTEGER, _ALLOCATABLE, dimension(:) :: num_elements
     
     TYPE(time_type) :: last_output, next_output, next_next_output
     TYPE(diag_fieldtype) :: f_type
     INTEGER, DIMENSION(4) :: axes
     INTEGER :: num_axes, total_elements, region_elements
     INTEGER :: n_diurnal_samples ! number of diurnal sample intervals, 1 or more
     TYPE(diag_grid) :: output_grid
     LOGICAL :: local_output, need_compute, phys_window, written_once
     LOGICAL :: reduced_k_range
     INTEGER :: imin, imax, jmin, jmax, kmin, kmax
     TYPE(time_type) :: Time_of_prev_field_data
  END TYPE output_field_type
  ! </TYPE>

  ! <TYPE NAME="diag_axis_type">
  !   <DESCRIPTION>
  !     Type to hold the diagnostic axis description.
  !   </DESCRIPTION>
  !   <DATA NAME="name" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="units" TYPE="CHARACTER(len=256)">
  !   </DATA>
  !   <DATA NAME="long_name" TYPE="CHARACTER(len=256)">
  !   </DATA>
  !   <DATA NAME="cart_name" TYPE="CHARACTER(len=1)">
  !   </DATA>
  !   <DATA NAME="data" TYPE="REAL, DIMENSION(:), POINTER">
  !   </DATA>
  !   <DATA NAME="start" TYPE="INTEGER, DIMENSION(MAX_SUBAXES)">
  !   </DATA>
  !   <DATA NAME="end" TYPE="INTEGER, DIMENSION(MAX_SUBAXES)">
  !   </DATA>
  !   <DATA NAME="subaxis_name" TYPE="CHARACTER(len=128), DIMENSION(MAX_SUBAXES)">
  !   </DATA>
  !   <DATA NAME="length" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="direction" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="edges" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="set" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="shift" TYPE="INTEGER">
  !   </DATA>
  !   <DATA NAME="Domain" TYPE="TYPE(domain1d)">
  !   </DATA>
  !   <DATA NAME="Domain2" TYPE="TYPE(domain2d)">
  !   </DATA>
  !   <DATA NAME="subaxis_domain2" TYPE="TYPE(domain2d), dimension(MAX_SUBAXES)">
  !   </DATA>
  !   <DATA NAME="aux" TYPE="CHARACTER(len=128)">
  !   </DATA>
  !   <DATA NAME="tile_count" TYPE="INTEGER">
  !   </DATA>
  TYPE diag_axis_type
     CHARACTER(len=128) :: name
     CHARACTER(len=256) :: units, long_name
     CHARACTER(len=1) :: cart_name
     REAL, DIMENSION(:), POINTER :: data
     INTEGER, DIMENSION(MAX_SUBAXES) :: start
     INTEGER, DIMENSION(MAX_SUBAXES) :: end
     CHARACTER(len=128), DIMENSION(MAX_SUBAXES) :: subaxis_name
     INTEGER :: length, direction, edges, set, shift
     TYPE(domain1d) :: Domain
     TYPE(domain2d) :: Domain2
     TYPE(domain2d), dimension(MAX_SUBAXES) :: subaxis_domain2
     CHARACTER(len=128) :: aux
     INTEGER :: tile_count
  END TYPE diag_axis_type
  ! </TYPE>

  ! <TYPE NAME="diag_global_att_type">
  !   <DESCRIPTION>
  !   </DESCRIPTION>
  !   <DATA NAME="grid_type" TYPE="CHARACTER(len=128)" DEFAULT="regular">
  !   </DATA>
  !   <DATA NAME="tile_name" TYPE="CHARACTER(len=128)" DEFAULT="N/A">
  !   </DATA>
  TYPE diag_global_att_type
     CHARACTER(len=128)   :: grid_type='regular'
     CHARACTER(len=128)   :: tile_name='N/A'
  END TYPE diag_global_att_type
  ! </TYPE>
  
  ! Private CHARACTER Arrays for the CVS version and tagname.
  CHARACTER(len=128),PRIVATE  :: version =&
       & '$Id: diag_data.F90,v 19.0.2.3 2012/05/14 18:40:11 Seth.Underwood Exp $'
  CHARACTER(len=128),PRIVATE  :: tagname =&
       & '$Name: siena_201211 $'

  ! <!-- Other public variables -->
  ! <DATA NAME="num_files" TYPE="INTEGER" DEFAULT="0">
  !   Number of output files currenly in use by the diag_manager.
  ! </DATA>
  ! <DATA NAME="num_input_fields" TYPE="INTEGER" DEFAULT="0">
  !   Number of input fields in use.
  ! </DATA>
  ! <DATA NAME="num_output_fields" TYPE="INTEGER" DEFAULT="0">
  !   Number of output fields in use.
  ! </DATA>
  ! <DATA NAME="null_axis_id" TYPE="INTEGER" />
  INTEGER :: num_files = 0
  INTEGER :: num_input_fields = 0
  INTEGER :: num_output_fields = 0
  INTEGER :: null_axis_id

  ! <!-- Namelist variables -->
  ! <DATA NAME="append_pelist_name" TYPE="LOGICAL" DEFAULT=".FALSE." />
  ! <DATA NAME="mix_snapshot_average_fields" TYPE="LOGICAL" DEFAULT=".FALSE." />
  ! <DATA NAME="max_files" TYPE="INTEGER" DEFAULT="31">
  !   Maximum number of output files allowed.  Increase via the diag_manager_nml namelist.
  ! </DATA>
  ! <DATA NAME="max_output_fields" TYPE="INTEGER" DEFAULT="300">
  !   Maximum number of output fields.  Increase via the diag_manager_nml namelist.
  ! </DATA>
  ! <DATA NAME="max_input_fields" TYPE="INTEGER" DEFAULT="300">
  !   Maximum number of input fields.  Increase via the diag_manager_nml namelist.
  ! </DATA>
  ! <DATA NAME="max_axes" TYPE="INTEGER" DEFAULT="60">
  !   Maximum number of independent axes.
  ! </DATA>
  ! <DATA NAME="do_diag_field_log" TYPE="LOGICAL" DEFAULT=".FALSE." />
  ! <DATA NAME="write_bytes_in_file" TYPE="LOGICAL" DEFAULT=".FALSE." />
  ! <DATA NAME="debug_diag_manager" TYPE="LOGICAL" DEFAULT=".FALSE." />
  ! <DATA NAME="max_num_axis_sets" TYPE="INTEGER" DEFAULT="25" />
  ! <DATA NAME="use_cmor" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   Indicates if we should overwrite the MISSING_VALUE to use the CMOR missing value.
  ! </DATA>
  ! <DATA NAME="ISSUE_OOR_WARNINGS" TYPE="LOGICAL" DEFAULT=".TRUE.">
  !   Issue warnings if the output field has values outside the given
  !   range for a variable.
  ! </DATA>
  ! <DATA NAME="OOR_WARNINGS_FATAL" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   Cause a fatal error if the output field has a value outside the
  !   given range for a variable.
  ! </DATA>
  LOGICAL :: append_pelist_name = .FALSE.
  LOGICAL :: mix_snapshot_average_fields =.FALSE.
  INTEGER :: max_files = 31 !< Maximum number of output files allowed.  Increase via diag_manager_nml.
  INTEGER :: max_output_fields = 300 !< Maximum number of output fields.  Increase via diag_manager_nml.
  INTEGER :: max_input_fields = 300 !< Maximum number of input fields.  Increase via diag_manager_nml.
  INTEGER :: max_axes = 60 !< Maximum number of independent axes.
  LOGICAL :: do_diag_field_log = .FALSE.
  LOGICAL :: write_bytes_in_file = .FALSE.
  LOGICAL :: debug_diag_manager = .FALSE.
  LOGICAL :: conserve_water = .TRUE. ! Undocumented namelist to control flushing of output files.
  INTEGER :: max_num_axis_sets = 25
  LOGICAL :: use_cmor = .FALSE.
  LOGICAL :: issue_oor_warnings = .TRUE.
  LOGICAL :: oor_warnings_fatal = .FALSE.

  ! <!-- netCDF variable -->
  ! <DATA NAME="FILL_VALUE" TYPE="REAL" DEFAULT="NF90_FILL_REAL">
  !   Fill value used.  Value will be <TT>NF90_FILL_REAL</TT> if using the
  !   netCDF module, otherwise will be 9.9692099683868690e+36.
  ! </DATA>
#ifdef use_netCDF
  REAL :: FILL_VALUE = NF_FILL_REAL  ! from file /usr/local/include/netcdf.inc
#else
  REAL :: FILL_VALUE = 9.9692099683868690e+36 
#endif

  INTEGER :: pack_size = 1 ! 1 for double and 2 for float

  ! <!-- REAL public variables -->
  ! <DATA NAME="EMPTY" TYPE="REAL" DEFAULT="0.0" />
  ! <DATA NAME="MAX_VALUE" TYPE="REAL" />
  ! <DATA NAME="MIN_VALUE" TYPE="REAL" />
  REAL :: EMPTY = 0.0
  REAL :: MAX_VALUE, MIN_VALUE

  ! <!-- Global data for all files -->
  ! <DATA NAME="base_time" TYPE="TYPE(time_type)" />
  ! <DATA NAME="base_year" TYPE="INTEGER" />
  ! <DATA NAME="base_month" TYPE="INTEGER" />
  ! <DATA NAME="base_day" TYPE="INTEGER" />
  ! <DATA NAME="base_hour" TYPE="INTEGER" />
  ! <DATA NAME="base_minute" TYPE="INTEGER" />
  ! <DATA NAME="base_second" TYPE="INTEGER" />
  ! <DATA NAME="global_descriptor" TYPE="CHARACTER(len=256)" />
  TYPE(time_type) :: base_time
  INTEGER :: base_year, base_month, base_day, base_hour, base_minute, base_second
  CHARACTER(len = 256):: global_descriptor

  ! <!-- ALLOCATABLE variables -->
  ! <DATA NAME="files" TYPE="TYPE(file_type), DIMENSION(:), SAVE, ALLOCATABLE" />
  ! <DATA NAME="input_fields" TYPE="TYPE(input_field_type), DIMENSION(:), ALLOCATABLE" />
  ! <DATA NAME="output_fields" TYPE="TYPE(output_field_type), DIMENSION(:), ALLOCATABLE" />
  TYPE(file_type), SAVE, ALLOCATABLE :: files(:)
  TYPE(input_field_type), ALLOCATABLE :: input_fields(:)
  TYPE(output_field_type), ALLOCATABLE :: output_fields(:)

  ! <!-- Even More Variables -->
  ! <DATA NAME="time_zero" TYPE="TYPE(time_type)" />
  ! <DATA NAME="first_send_data_call" TYPE="LOGICAL" DEFAULT=".TRUE." />
  ! <DATA NAME="module_is_initialized" TYPE="LOGICAL" DEFAULT=".FALSE." />
  ! <DATA NAME="diag_log_unit" TYPE="INTEGER" />
  ! <DATA NAME="time_unit_list" TYPE="CHARACTER(len=10), DIMENSION(6)"
  !       DEFAULT="(/'seconds   ', 'minutes   ', 'hours     ', 'days      ', 'months    ', 'years     '/)" />
  ! <DATA NAME="filename_appendix" TYPE="CHARACTER(len=32)" DEFAULT="" />
  ! <DATA NAME="pelist_name" TYPE="CHARACTER(len=32)" />
  TYPE(time_type) :: time_zero
  LOGICAL :: first_send_data_call = .TRUE.
  LOGICAL :: module_is_initialized = .FALSE.
  INTEGER :: diag_log_unit
  CHARACTER(len=10), DIMENSION(6) :: time_unit_list = (/'seconds   ', 'minutes   ',&
       & 'hours     ', 'days      ', 'months    ', 'years     '/)
  CHARACTER(len=32), SAVE :: filename_appendix = ''
  CHARACTER(len=32) :: pelist_name
  INTEGER :: oor_warning = WARNING
  
END MODULE diag_data_mod
