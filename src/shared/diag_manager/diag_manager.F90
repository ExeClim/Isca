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

MODULE diag_manager_mod
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
  !   Matt Harrison
  ! </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
  !   Giang Nong
  ! </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
  !   Seth Underwood
  ! </CONTACT>
  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/" />
  ! <OVERVIEW>
  !   <TT>diag_manager_mod</TT> is a set of simple calls for parallel diagnostics
  !   on distributed systems. It is geared toward the writing of data in netCDF
  !   format.
  ! </OVERVIEW>
  ! <DESCRIPTION>
  !   <TT>diag_manager_mod</TT> provides a convenient set of interfaces for
  !   writing data to disk.  It is built upon the parallel I/O interface of FMS
  !   code <TT>/shared/mpp/mpp_io.F90</TT>.
  !
  !   A single group of calls to the <TT>diag_manager_mod</TT> interfaces
  !   provides data to disk at any number of sampling and/or averaging intervals
  !   specified at run-time. Run-time specification of diagnostics are input
  !   through the diagnostics table.
  !
  !   <H4>Usage</H4>
  !   Use of <TT>diag_manager</TT> includes the following steps:
  !   <OL>
  !     <LI> Create diag_table as described in the
  !          <LINK SRC="diag_table.html">diag_table.F90</LINK>
  !          documentation.</LI>
  !     <LI> Call <LINK SRC="#diag_manager_init"><TT>diag_manager_init</TT></LINK> to initialize
  !          diag_manager_mod.</LI>
  !     <LI> Call <LINK SRC="#register_diag_field"><TT>register_diag_field</TT></LINK> to register the field to be
  !          output.
  !          <B>NOTE:</B> ALL fields in diag_table should be registered <I>BEFORE</I>
  !          the first send_data call</LI>
  !     <LI> Call <LINK SRC="#send_data"><TT>send_data</TT></LINK> to send data to output fields </LI>
  !     <LI> Call <LINK SRC="#diag_manager_end"><TT>diag_manager_end</TT></LINK> to exit diag_manager </LI>
  !   </OL>
  !
  !   <H4>Features</H4>
  !   Features of <TT>diag_manager_mod</TT>:
  !   <OL>
  !     <LI> Ability to output from 0D arrays (scalars) to 3D arrays.</LI>
  !     <LI> Ability to output time average of fields that have time dependent
  !          mask.</LI>
  !     <LI> Give optional warning if <TT>register_diag_field</TT> fails due to
  !          misspelled module name or field name.</LI>
  !     <LI> Check if a field is registered twice.</LI>
  !     <LI> Check for duplicate lines in diag_table. </LI>
  !     <LI> <LINK SRC="diag_table.html">diag_table</LINK> can contain fields
  !          that are NOT written to any files. The file name in diag_table of
  !          these fields is <TT>null</TT>.</LI>
  !     <LI> By default, a field is output in its global grid.  The user can now 
  !          output a field in a specified region.  See
  !          <LINK SRC="#send_data"><TT>send_data</TT></LINK> for more details.</LI>
  !     <LI> To check if the diag table is set up correctly, user should set
  !          <TT>debug_diag_manager=.true.</TT> in diag_manager namelist, then
  !          the the content of diag_table is printed in stdout.</LI>
  !     <LI> New optional format of file information in <LINK SRC="diag_table.html">diag_table</LINK>.It is possible to have just
  !          one file name and reuse it many times. A time string will be appended to the base file name each time a new file is
  !          opened. The time string can be any combination from year to second of current model time.
  !
  !          Here is an example file line: <BR />
  !          <PRE>"file2_yr_dy%1yr%3dy",2,"hours",1,"hours","Time", 10, "days", "1 1 7 0 0 0", 6, "hours"</PRE>
  !          <BR />
  !
  !          From left to right we have: 
  !          <UL>
  !            <LI>file name</LI> 
  !            <LI>output frequency</LI>
  !            <LI>output frequency unit</LI>
  !            <LI>Format (should always be 1)</LI>
  !            <LI>time axis unit</LI>
  !            <LI>time axis name</LI>
  !            <LI>frequency for creating new file</LI>
  !            <LI>unit for creating new file</LI>
  !            <LI>start time of the new file</LI>
  !            <LI>file duration</LI>
  !            <LI>file duration unit.</LI>
  !          </UL>
  !          The 'file duration', if absent, will be equal to frequency for creating a new file.
  !
  !          Thus, the above means: create a new file every 10 days, each file will last 6 hours from creation time, no files will
  !          be created before time "1 1 7 0 0 0".
  !
  !          In this example the string
  !          <TT>10, "days", "1 1 7 0 0 0", 6, "hours"</TT> is optional.
  !
  !          Keywords for the time string suffix is
  !          <TT>%xyr,%xmo,%xdy,%xhr,%xmi,%xsc</TT> where <TT>x</TT> is a
  !          mandatory 1 digit number specifying the width of field used in
  !          writing the string</LI>
  !     <LI> New time axis for time averaged fields.  Users can use a namelist option to handle the time value written
  !          to time axis for time averaged fields.
  !
  !          If <TT>mix_snapshot_average_fields=.true.</TT> then a time averaged file will have time values corresponding to
  !          ending time_bound e.g. January monthly average is labeled Feb01. Users can have both snapshot and averaged fields in
  !          one file.
  !
  !          If <TT>mix_snapshot_average_fields=.false.</TT> The time value written to time axis for time averaged fields is the
  !          middle on the averaging time. For example, January monthly mean will be written at Jan 16 not Feb 01 as
  !          before. However, to use this new feature users should <B>separate</B> snapshot fields and time averaged fields in
  !          <B>different</B> files or a fatal error will occur.
  !
  !          The namelist <B>default</B> value is <TT>mix_snapshot_average_fields=.false.</TT></LI>
  !     <LI> Time average, Max and Min, and diurnal. In addition to time average users can also get then Max or Min value
  !          during the same interval of time as time average. For this purpose, in the diag table users must replace
  !          <TT>.true.</TT> or <TT>.false.</TT> by "<TT>max</TT>" or "<TT>min</TT>".  <B><I>Note:</I></B> Currently, max
  !          and min are not available for regional output.
  !
  !          A diurnal average can also be requested using <TT>diurnal##</TT> where <TT>##</TT> are the number of diurnal
  !          sections to average.</LI>
  !     <LI> <TT>standard_name</TT> is added as optional argument in <LINK SRC="#register_diag_field"><TT>register_diag_field</TT>
  !          </LINK>.</LI>
  !     <LI>When namelist variable <TT>debug_diag_manager = .true.</TT> array
  !         bounds are checked in <LINK SRC="#send_data"><TT>send_data</TT></LINK>.</LI>
  !     <LI>Coordinate attributes can be written in the output file if the
  !         argument "<TT>aux</TT>" is given in <LINK SRC="diag_axis.html#diag_axis_init"><TT>diag_axis_init</TT></LINK>. The
  !         corresponding fields (geolat/geolon) should also be written to the
  !         same file.</LI>
  !   </OL>
  !
  ! </DESCRIPTION>

  ! <NAMELIST NAME="diag_manager_nml">
  !   <DATA NAME="append_pelist_name" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   </DATA>
  !   <DATA NAME="mix_snapshot_average_fields" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !     Set to .TRUE. to allow both time average and instantaneous fields in the same output file.
  !   </DATA>
  !   <DATA NAME="max_files" TYPE="INTEGER" DEFULT="31">
  !   </DATA>
  !   <DATA NAME="max_output_fields" TYPE="INTEGER" DEFAULT="300">
  !   </DATA>
  !   <DATA NAME="max_input_fields" TYPE="INTEGER" DEFAULT="300">
  !   </DATA>
  !   <DATA NAME="max_axes" TYPE="INTEGER" DEFAULT="60">
  !   </DATA>
  !   <DATA NAME="do_diag_field_log" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   </DATA>
  !   <DATA NAME="write_bytes_in_files" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   </DATA>
  !   <DATA NAME="debug_diag_manager" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   </DATA>
  !   <DATA NAME="max_num_axis_sets" TYPE="INTEGER" DEFAULT="25">
  !   </DATA>
  !   <DATA NAME="use_cmor" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !     Let the <TT>diag_manager</TT> know if the missing value (if supplied) should be overridden to be the
  !     CMOR standard value of -1.0e20.
  !   </DATA>
  !   <DATA NAME="issue_oor_warnings" TYPE="LOGICAL" DEFAULT=".TRUE.">
  !     If <TT>.TRUE.</TT>, then the <TT>diag_manager</TT> will check for values outside the valid range.  This range is defined in
  !     the model, and passed to the <TT>diag_manager_mod</TT> via the OPTIONAL variable range in the <TT>register_diag_field</TT>
  !     function.
  !   </DATA>
  !   <DATA NAME="oor_warnings_fatal" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !     If <TT>.TRUE.</TT> then <TT>diag_manager_mod</TT> will issue a <TT>FATAL</TT> error if any values for the output field are
  !     outside the given range.
  !   </DATA>
  ! </NAMELIST>

  USE time_manager_mod, ONLY: set_time, set_date, get_time, time_type, OPERATOR(>=), OPERATOR(>),&
       & OPERATOR(<), OPERATOR(==), OPERATOR(/=), OPERATOR(/), OPERATOR(+)
  USE mpp_io_mod, ONLY: mpp_open, mpp_close
  USE mpp_mod, ONLY: mpp_get_current_pelist, mpp_pe, mpp_npes, mpp_root_pe, mpp_sum

#ifdef INTERNAL_FILE_NML
  USE mpp_mod, ONLY: input_nml_file
#else
  USE fms_mod, ONLY: open_namelist_file, close_file
#endif

  USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE, stdout, stdlog, write_version_number,&
       & file_exist, fms_error_handler, check_nml_error
  USE diag_axis_mod, ONLY: diag_axis_init, get_axis_length, max_axes, get_axis_num
  USE diag_util_mod, ONLY: get_subfield_size, log_diag_field_info, update_bounds,&
       & check_out_of_bounds, check_bounds_are_exact_dynamic, check_bounds_are_exact_static,&
       & diag_time_inc, find_input_field, init_input_field, init_output_field,&
       & diag_data_out, write_static, get_date_dif, get_subfield_vert_size, sync_file_times
  USE diag_data_mod, ONLY: max_files, CMOR_MISSING_VALUE, DIAG_OTHER, DIAG_OCEAN, DIAG_ALL, EVERY_TIME,&
       & END_OF_RUN, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, DIAG_YEARS, num_files,&
       & max_input_fields, max_output_fields, num_output_fields, EMPTY, FILL_VALUE, null_axis_id,&
       & MAX_VALUE, MIN_VALUE, base_time, base_year, base_month, base_day,&
       & base_hour, base_minute, base_second, global_descriptor, coord_type, files, input_fields,&
       & output_fields, Time_zero, append_pelist_name, mix_snapshot_average_fields,&
       & first_send_data_call, do_diag_field_log, write_bytes_in_file, debug_diag_manager,&
       & diag_log_unit, time_unit_list, pelist_name, module_is_initialized, max_num_axis_sets,&
       & use_cmor, issue_oor_warnings, oor_warnings_fatal, oor_warning, filename_appendix, pack_size,&
       & conserve_water
  USE diag_table_mod, ONLY: parse_diag_table
  USE diag_output_mod, ONLY: get_diag_global_att, set_diag_global_att
  USE diag_grid_mod, ONLY: diag_grid_init, diag_grid_end
  USE constants_mod, ONLY: SECONDS_PER_HOUR, SECONDS_PER_MINUTE

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: diag_manager_init, send_data, send_tile_averaged_data, diag_manager_end,&
       & register_diag_field, register_static_field, diag_axis_init, get_base_time, get_base_date,&
       & need_data, average_tiles, DIAG_ALL, DIAG_OCEAN, DIAG_OTHER, get_date_dif, DIAG_SECONDS,&
       & DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, DIAG_YEARS, get_diag_global_att,&
       & set_diag_global_att
  ! Public interfaces from diag_grid_mod
  PUBLIC :: diag_grid_init, diag_grid_end
  PUBLIC :: set_diag_filename_appendix
  PUBLIC :: diag_manager_set_time_end, diag_send_complete

  ! version number of this module
  CHARACTER(len=128), PARAMETER :: version =&
       & '$Id: diag_manager.F90,v 19.0.4.8 2012/05/21 14:06:11 Seth.Underwood Exp $'
  CHARACTER(len=128), PARAMETER :: tagname =&
       & '$Name: siena_201211 $'  

  type(time_type) :: Time_end

  ! <INTERFACE NAME="send_data">
  !   <TEMPLATE>
  !     send_data(diag_field_id, field, time, is_in, js_in, ks_in,
  !             mask, rmask, ie_in, je_in, ke_in, weight)
  !   </TEMPLATE>
  !   <OVERVIEW>
  !     Send data over to output fields. 
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     <TT>send_data</TT> is overloaded for fields having zero dimension
  !     (scalars) to 3 dimension.  <TT>diag_field_id</TT> corresponds to the id
  !     returned from a previous call to <TT>register_diag_field</TT>. The field
  !     array is restricted to the computational range of the array. Optional
  !     argument <TT>is_in</TT> can be used to update sub-arrays of the entire
  !     field. Additionally, an optional logical or real mask can be used to
  !     apply missing values to the array.
  !
  !     If a field is declared to be <TT>mask_variant</TT> in
  !     <TT>register_diag_field</TT> logical mask should be mandatory.
  !
  !     For the real  mask, the mask is applied if the mask value is less than
  !     0.5.
  !
  !     By default, a field will be written out entirely in its global grid.
  !     Users can also specify regions in which the field will be output. The
  !     region is specified in diag-table just before the end of output_field
  !     replacing "none". 
  !
  !     For example, by default:
  !
  !     "ocean_mod","Vorticity","vorticity","file1","all",.false.,"none",2 
  !
  !     for regional output:
  !
  !     "ocean_mod","Vorticity","vorticity_local","file2","all",.false.,"0.5 53.5 -89.5 -28.5 -1 -1",2
  !
  !     The format of region is "<TT>xbegin xend ybegin yend zbegin zend</TT>".
  !     If it is a 2D field use (-1 -1) for (zbegin zend) as in the example
  !     above. For a 3D field use (-1 -1) for (zbegin zend) when you want to
  !     write the entire vertical extent, otherwise specify real coordinates.
  !     The units used for region are the actual units used in grid_spec.nc
  !     (for example degrees for lat, lon). a FATAL error will occur if the
  !     region's boundaries are not found in grid_spec.nc.
  !
  !     Regional output on the cubed sphere is also supported.  To use regional output on the cubed sphere, first the grid
  !     information needs to be sent to <TT>diag_manager_mod</TT> using the <LINK SRC="diag_grid.html#diag_grid_init"><TT>
  !     diag_grid_init</TT></LINK> subroutine.  <B><I>NOTE:</I></B> Regions must be confined to a single tile.  Regions spanning
  !     tiles will be ignored.  A future release will allow multi-tile regions.
  ! 
  !     <B><I>NOTE:</I></B> When using regional output the files containing regional 
  !     outputs should be different from files containing global (default) output. 
  !     It is a FATAL error to have one file containing both regional and global 
  !     results. For maximum flexibility and independence from PE counts one file 
  !     should contain just one region.
  ! 
  !     Time averaging is supported in regional output.
  ! 
  !     Physical fields (written in "physics windows" of atmospheric code) are 
  !     currently fully supported for regional outputs.
  !
  !     Note of dimension of field in send_data
  !
  !     Most fields are defined in data_domain but used in compute domain. In 
  !     <TT>send_data</TT> users can pass EITHER field in data domain OR field in 
  !     compute domain. If data domain is used, users should also pass the starting and 
  !     ending indices of compute domain (isc, iec ...). If compute domain is used no 
  !     indices are needed. These indices are for determining halo exclusively. If 
  !     users want to ouput the field partially they should use regional output as 
  !     mentioned above.
  !
  !     Weight in Time averaging is now supported, each time level may have a
  !     different weight. The default of weight is 1.
  !   </DESCRIPTION>
  !   <IN NAME="diag_field_id" TYPE="INTEGER"> </IN>
  !   <IN NAME="field" TYPE="REAL, DIMENSION(:,:,:)"> </IN>
  !   <IN NAME="time" TYPE="TYPE(time_type)"> </IN>
  !   <IN NAME="is_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="js_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="ks_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="mask" TYPE="LOGICAL, DIMENSION(:,:,:), OPTIONAL"></IN>
  !   <IN NAME="rmask" TYPE="REAL, DIMENSION(:,:,:), OPTIONAL"></IN>
  !   <IN NAME="ie_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="je_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="ke_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="weight" TYPE="REAL, OPTIONAL"></IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL"></OUT>
  INTERFACE send_data
     MODULE PROCEDURE send_data_0d
     MODULE PROCEDURE send_data_1d
     MODULE PROCEDURE send_data_2d
     MODULE PROCEDURE send_data_3d
  END INTERFACE
  ! </INTERFACE>

  ! <INTERFACE NAME="register_diag_field">
  !   <OVERVIEW>
  !      Register Diagnostic Field.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION register_diag_field (module_name, field_name, axes, init_time,
  !           long_name, units, missing_value, range, mask_variant, standard_name,
  !           verbose)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !      Return field index for subsequent calls to
  !      <LINK SRC="#send_data">send_data</LINK>.
  !
  !      <TT>axes</TT> are the axis ID returned from <TT>diag_axis_init</TT>,
  !      <TT>axes</TT> are required for fields of 1-3 dimension and NOT required
  !      for scalars.
  !
  !      For a static scalar (constant) <TT>init_time</TT> is not needed.
  !
  !      Optional <TT>mask_variant</TT> is for fields that have a time-dependent
  !      mask. If <TT>mask_variant</TT> is true then <TT>mask</TT> must be
  !      present in argument list of <TT>send_data</TT>.
  !
  !      The pair (<TT>module_name</TT>, <TT>fieldname</TT>) should be registered
  !      only once or a FATAL error will occur.
  !    </DESCRIPTION>
  !    <IN NAME="module_name" TYPE="CHARACTER(len=*)" />
  !    <IN NAME="field_name" TYPE="CHARACTER(len=*)" />
  !    <IN NAME="axes" TYPE="INTEGER, DIMENSION(:)" />
  !    <IN NAME="init_time" TYPE="TYPE(time_type)" />
  !    <IN NAME="long_name" TYPE="CHARACTER(len=*)" />
  !    <IN NAME="units" TYPE="CHARACTER(len=*)" />
  !    <IN NAME="missing_value" TYPE="REAL" />
  !    <IN NAME="range" TYPE="REAL, DIMENSION(2)" />
  !    <IN NAME="mask_variant" TYPE="LOGICAL" /> 
  !    <IN NAME="standard_name" TYPE="CHARACTER(len=*)" />
  INTERFACE register_diag_field
     MODULE PROCEDURE register_diag_field_scalar
     MODULE PROCEDURE register_diag_field_array
  END INTERFACE
  ! </INTERFACE>

  !  <INTERFACE NAME="send_tile_averaged_data">
  !    <OVERVIEW>
  !      Send tile-averaged data over to output fields. 
  !    </OVERVIEW>
  !    <TEMPLATE>
  !      LOGICAL send_tile_averaged_data(diag_field_id, field, area, time, mask)
  !    </TEMPLATE>
  !    <DESCRIPTION>
  !      <TT>send_tile_averaged_data</TT> is overloaded for 3D and 4D arrays. 
  !      <TT>diag_field_id</TT> corresponds to the ID returned by previous call
  !      to <TT>register_diag_field</TT>. Logical masks can be used to mask out
  !      undefined and/or unused values.  Note that the dimension of output field
  !      is smaller by one than the dimension of the data, since averaging over
  !      tiles (3D dimension) is performed.
  !    </DESCRIPTION>
  !    <IN NAME="diag_field_id" TYPE="INTEGER" />
  !    <IN NAME="field" TYPE="REAL" DIM="(:,:,:)" />
  !    <IN NAME="area" TYPE="REAL" DIM="(:,:,:)" />
  !    <IN NAME="time" TYPE="TYPE(time_type)" DIM="(:,:,:)" />
  !    <IN NAME="mask" TYPE="LOGICAL" DIM="(:,:,:)" />
  INTERFACE send_tile_averaged_data
     MODULE PROCEDURE send_tile_averaged_data2d
     MODULE PROCEDURE send_tile_averaged_data3d
  END INTERFACE
  ! </INTERFACE>

CONTAINS

  ! <FUNCTION NAME="register_diag_field_scalar" INTERFACE="register_diag_field">
  !   <IN NAME="module_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="field_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="axes" TYPE="Not Applicable" />
  !   <IN NAME="init_time" TYPE="TYPE(time_type), OPTIONAL" />
  !   <IN NAME="long_name" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="units" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="missing_value" TYPE="REAL, OPTIONAL" />
  !   <IN NAME="range" TYPE="REAL, DIMENSION(2), OPTIONAL" />
  !   <IN NAME="mask_variant" TYPE="Not Applicable" />
  !   <IN NAME="standard_name" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="do_not_log" TYPE="LOGICAL, OPTIONAL" />
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL" />
  INTEGER FUNCTION register_diag_field_scalar(module_name, field_name, init_time, &
       & long_name, units, missing_value, range, standard_name, do_not_log, err_msg)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name
    TYPE(time_type), OPTIONAL, INTENT(in) :: init_time
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name, units, standard_name
    REAL, OPTIONAL, INTENT(in) :: missing_value
    REAL,  DIMENSION(2), OPTIONAL, INTENT(in) :: RANGE
    LOGICAL, OPTIONAL, INTENT(in) :: do_not_log ! if TRUE, field information is not logged
    CHARACTER(len=*), OPTIONAL, INTENT(out):: err_msg
 
    IF ( PRESENT(init_time) ) THEN
       register_diag_field_scalar = register_diag_field_array(module_name, field_name,&
            & (/null_axis_id/), init_time,long_name, units, missing_value, range, &
            & standard_name=standard_name, do_not_log=do_not_log, err_msg=err_msg)
    ELSE
       IF ( PRESENT(err_msg) ) err_msg = ''
       register_diag_field_scalar = register_static_field(module_name, field_name,&
            & (/null_axis_id/),long_name, units, missing_value, range,&
            & standard_name=standard_name, do_not_log=do_not_log)
    END IF
  END FUNCTION register_diag_field_scalar
  ! </FUNCTION>

  ! <FUNCTION NAME="register_diag_field_array" INTERFACE="register_diag_field">
  !   <IN NAME="module_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="field_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="axes" TYPE="INTEGER, DIMENSION(:)" />
  !   <IN NAME="init_time" TYPE="TYPE(time_type)" />
  !   <IN NAME="long_name" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="units" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="missing_value" TYPE="REAL, OPTIONAL" />
  !   <IN NAME="range" TYPE="REAL, DIMENSION(2), OPTIONAL" />
  !   <IN NAME="mask_variant" TYPE="LOGICAL, OPTIONAL" />
  !   <IN NAME="standard_name" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="do_not_log" TYPE="LOGICAL, OPTIONAL" />
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="interp_method" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="tile_count" TYPE="INTEGER, OPTIONAL" />
  INTEGER FUNCTION register_diag_field_array(module_name, field_name, axes, init_time, &
       & long_name, units, missing_value, range, mask_variant, standard_name, verbose,&
       & do_not_log, err_msg, interp_method, tile_count)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name
    INTEGER, INTENT(in) :: axes(:)
    TYPE(time_type), INTENT(in) :: init_time
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name, units, standard_name
    REAL, OPTIONAL, INTENT(in) :: missing_value, RANGE(2)
    LOGICAL, OPTIONAL, INTENT(in) :: mask_variant,verbose
    LOGICAL, OPTIONAL, INTENT(in) :: do_not_log ! if TRUE, field info is not logged
    CHARACTER(len=*), OPTIONAL, INTENT(out):: err_msg
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: interp_method
    INTEGER, OPTIONAL, INTENT(in) :: tile_count

    INTEGER :: field, j, ind, file_num, freq
    INTEGER :: output_units
    INTEGER :: stdout_unit
    LOGICAL :: mask_variant1, verbose1
    CHARACTER(len=128) :: msg

    ! get stdout unit number
    stdout_unit = stdout()

    IF ( PRESENT(mask_variant) ) THEN
       mask_variant1 = mask_variant
    ELSE
       mask_variant1 = .FALSE.
    END IF

    IF ( PRESENT(verbose) ) THEN 
       verbose1 = verbose
    ELSE 
       verbose1 = .FALSE.
    END IF

    IF ( PRESENT(err_msg) ) err_msg = ''
    
    ! Call register static, then set static back to false
    register_diag_field_array = register_static_field(module_name, field_name, axes,&
         & long_name, units, missing_value, range, mask_variant1, standard_name=standard_name,&
         & DYNAMIC=.TRUE., do_not_log=do_not_log, interp_method=interp_method, tile_count=tile_count)

    IF ( .NOT.first_send_data_call ) THEN 
       ! <ERROR STATUS="WARNING">
       !   module/output_field <module_name>/<field_name> registered AFTER first
       !   send_data call, TOO LATE
       ! </ERROR>
       IF ( mpp_pe() == mpp_root_pe() ) &
            & CALL  error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '&
            &//TRIM(module_name)//'/'// TRIM(field_name)//&
            &' registered AFTER first send_data call, TOO LATE', WARNING)  
    END IF

    IF ( register_diag_field_array < 0 ) THEN
       ! <ERROR STATUS="WARNING">
       !   module/output_field <modul_name>/<field_name> NOT found in diag_table
       ! </ERROR>
       IF ( debug_diag_manager .OR. verbose1 ) THEN 
          IF ( mpp_pe() == mpp_root_pe() ) &
               & CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '&
               &//TRIM(module_name)//'/'// TRIM(field_name)//' NOT found in diag_table',&
               & WARNING) 
       END IF
    ELSE 
       input_fields(register_diag_field_array)%static = .FALSE.
       field = register_diag_field_array
       IF ( PRESENT(standard_name) ) input_fields(field)%standard_name = standard_name  

       DO j = 1, input_fields(field)%num_output_fields
          ind = input_fields(field)%output_fields(j)
          output_fields(ind)%static = .FALSE.
          ! Set up times in output_fields
          output_fields(ind)%last_output = init_time
          ! Get output frequency from for the appropriate output file
          file_num = output_fields(ind)%output_file
          IF ( file_num == max_files ) CYCLE
          IF ( output_fields(ind)%local_output ) THEN
             IF ( output_fields(ind)%need_compute) THEN         
                files(file_num)%local = .TRUE.
             END IF
          END IF

          ! Need to sync start_time of file with init time of model
          ! and close_time calculated with the duration of the file.
          ! Also, increase next_open until it is greater than init_time.
          CALL sync_file_times(file_num, init_time, err_msg=msg)
          IF ( msg /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::register_diag_field', TRIM(msg), err_msg) ) RETURN
          END IF

          freq = files(file_num)%output_freq
          output_units = files(file_num)%output_units
          output_fields(ind)%next_output = diag_time_inc(init_time, freq, output_units, err_msg=msg)
          IF ( msg /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::register_diag_field',&
                  & ' file='//TRIM(files(file_num)%name)//': '//TRIM(msg),err_msg)) RETURN
          END IF
          output_fields(ind)%next_next_output = &
               & diag_time_inc(output_fields(ind)%next_output, freq, output_units, err_msg=msg)
          IF ( msg /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::register_diag_field',&
                  &' file='//TRIM(files(file_num)%name)//': '//TRIM(msg),err_msg) ) RETURN
          END IF
          IF ( debug_diag_manager .AND. mpp_pe() == mpp_root_pe() .AND. output_fields(ind)%local_output ) THEN
             WRITE (msg,'(" lon(",F5.1,", ",F5.1,"), lat(",F5.1,", ",F5.1,"), dep(",F5.1,", ",F5.1,")")') &
                  & output_fields(ind)%output_grid%start(1),output_fields(ind)%output_grid%end(1),&
                  & output_fields(ind)%output_grid%start(2),output_fields(ind)%output_grid%end(2),&
                  & output_fields(ind)%output_grid%start(3),output_fields(ind)%output_grid%end(3)
             WRITE(stdout_unit,* ) 'module/output_field '//TRIM(module_name)//'/'//TRIM(field_name)// &
                  & ' will be output in region:'//TRIM(msg)
          END IF
       END DO
    END IF
  END FUNCTION register_diag_field_array
  ! </FUNCTION>

  ! <FUNCTION NAME="register_static_field">
  !   <OVERVIEW>
  !     Register Static Field.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION register_static_field(module_name, field_name, axes,
  !       long_name, units, missing_value, range, mask_variant, standard_name,
  !       dynamic, do_not_log, interp_method, tile_count)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return field index for subsequent call to send_data.
  !   </DESCRIPTION>
  !   <IN NAME="module_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="field_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="axes" TYPE="INTEGER, DIMENSION(:)" />
  !   <IN NAME="long_name" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="units" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="missing_value" TYPE="REAL, OPTIONAL" />
  !   <IN NAME="range" TYPE="REAL, DIMENSION(2), OPTIONAL" />
  !   <IN NAME="mask_variang" TYPE="LOGICAL, OPTIONAL" DEFAULT=".FALSE."/>
  !   <IN NAME="standard_name" TYPE="CHARACTER(len=*), OPTIONAL" />
  !   <IN NAME="dynamic" TYPE="LOGICAL, OPTIONAL" DEFAULT=".FALSE."/>
  !   <IN NAME="do_not_log" TYPE="LOGICAL, OPTIONAL" DEFAULT=".TRUE."/>
  !   <IN NAME="interp_method" TYPE="CHARACTER(len=*), OPTIOANL" />
  !   <IN NAME="tile_count" TYPE="INTEGER, OPTIONAL" />
  INTEGER FUNCTION register_static_field(module_name, field_name, axes, long_name, units,&
       & missing_value, range, mask_variant, standard_name, DYNAMIC, do_not_log, interp_method,&
       & tile_count)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name
    INTEGER, DIMENSION(:), INTENT(in) :: axes
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name, units, standard_name
    REAL, OPTIONAL, INTENT(in) :: missing_value
    REAL, DIMENSION(2), OPTIONAL, INTENT(in) :: range
    LOGICAL, OPTIONAL, INTENT(in) :: mask_variant
    LOGICAL, OPTIONAL, INTENT(in) :: DYNAMIC
    LOGICAL, OPTIONAL, INTENT(in) :: do_not_log ! if TRUE, field information is not logged
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: interp_method
    INTEGER,          OPTIONAL, INTENT(in) :: tile_count

    REAL :: missing_value_use
    INTEGER :: field, num_axes, j, out_num, k
    INTEGER, DIMENSION(3) :: siz, local_siz, local_start, local_end ! indices of local domain of global axes
    INTEGER :: tile, file_num
    LOGICAL :: mask_variant1, dynamic1, allow_log
    CHARACTER(len=128) :: msg

    ! Fatal error if the module has not been initialized.
    IF ( .NOT.module_is_initialized ) THEN 
       ! <ERROR STATUS="FATAL">diag_manager has NOT been initialized</ERROR>
       CALL error_mesg ('diag_manager_mod::register_static_field', 'diag_manager has NOT been initialized', FATAL)
    END IF

    ! Check if OPTIONAL parameters were passed in.
    IF ( PRESENT(missing_value) ) THEN
       IF ( use_cmor ) THEN 
          missing_value_use = CMOR_MISSING_VALUE
       ELSE
          missing_value_use = missing_value
       END IF
    END IF
    
    IF ( PRESENT(mask_variant) ) THEN 
       mask_variant1 = mask_variant
    ELSE 
       mask_variant1 = .FALSE.
    END IF
    
    IF ( PRESENT(DYNAMIC) ) THEN
       dynamic1 = DYNAMIC
    ELSE
       dynamic1 = .FALSE.
    END IF

    IF ( PRESENT(tile_count) ) THEN 
       tile = tile_count
    ELSE
       tile = 1
    END IF

    IF ( PRESENT(do_not_log) ) THEN
       allow_log = .NOT.do_not_log
    ELSE
       allow_log = .TRUE. 
    END IF

    ! Namelist do_diag_field_log is by default false.  Thus to log the
    ! registration of the data field, but the OPTIONAL parameter
    ! do_not_log == .FALSE. and the namelist variable
    ! do_diag_field_log == .TRUE..
    IF ( do_diag_field_log.AND.allow_log ) THEN
       CALL log_diag_field_info (module_name, field_name, axes, &
            & long_name, units, missing_value=missing_value, range=range, &
            & DYNAMIC=dynamic1)
    END IF

    register_static_field = find_input_field(module_name, field_name, 1)
    field = register_static_field
    ! Negative index returned if this field was not found in the diag_table.
    IF ( register_static_field < 0 ) RETURN

    IF ( tile > 1 ) THEN
       IF ( .NOT.input_fields(field)%register ) THEN
          ! <ERROR STATUS="FATAL">
          !   module/output_field <module_name>/<field_name> is not registered for tile_count = 1,
          !   should not register for tile_count > 1
          ! </ERROR>
          CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '//trim(module_name)//'/'//&
           & TRIM(field_name)//' is not registered for tile_count = 1, should not register for tile_count > 1',&
           & FATAL)    
       END IF

       CALL init_input_field(module_name, field_name, tile)
       register_static_field = find_input_field(module_name, field_name, tile)
       DO j = 1, input_fields(field)%num_output_fields
          out_num = input_fields(field)%output_fields(j)
          file_num = output_fields(out_num)%output_file
          IF(input_fields(field)%local) THEN
             CALL init_output_field(module_name, field_name,output_fields(out_num)%output_name,&
                  & files(file_num)%name,output_fields(out_num)%time_method, output_fields(out_num)%pack,&
                  & tile, input_fields(field)%local_coord)
          ELSE
             CALL init_output_field(module_name, field_name,output_fields(out_num)%output_name,&
                  & files(file_num)%name,output_fields(out_num)%time_method, output_fields(out_num)%pack, tile)
          END IF
       END DO
       field = register_static_field       
    END IF

    ! Store information for this input field into input field table

    ! Set static to true, if called by register_diag_field this is
    ! flipped back to false
    input_fields(field)%static = .TRUE.
    ! check if the field is registered twice
    IF ( input_fields(field)%register .AND. mpp_pe() == mpp_root_pe() ) THEN
       ! <ERROR STATUS="FATAL">
       !   module/output_field <module_name>/<field_name> ALREADY Registered, should
       !   not register twice
       ! </ERROR>
       CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '//trim(module_name)//'/'//&
            & TRIM(field_name)//' ALREADY registered, should not register twice', FATAL) 
    END IF

    ! Set flag that this field was registered
    input_fields(field)%register = .TRUE.
    ! set flag for mask: does it change with time?
    input_fields(field)%mask_variant = mask_variant1
    ! Set flag for mask warning
    input_fields(field)%issued_mask_ignore_warning = .FALSE.

    ! Check for more OPTIONAL parameters.
    IF ( PRESENT(long_name) ) THEN
       input_fields(field)%long_name = TRIM(long_name)
    ELSE
       input_fields(field)%long_name = input_fields(field)%field_name
    END IF
    
    IF ( PRESENT(standard_name) ) input_fields(field)%standard_name = standard_name 
    
    IF ( PRESENT(units) ) THEN
       input_fields(field)%units = TRIM(units)
    ELSE
       input_fields(field)%units = 'none'
    END IF
    
    IF ( PRESENT(missing_value) ) THEN
       input_fields(field)%missing_value = missing_value_use
       input_fields(field)%missing_value_present = .TRUE.
    ELSE
       input_fields(field)%missing_value_present = .FALSE.
    END IF
    
    IF ( PRESENT(range) ) THEN
       input_fields(field)%range = range
       input_fields(field)%range_present = .TRUE.
    ELSE
       input_fields(field)%range = (/ 1., 0. /)
       input_fields(field)%range_present = .FALSE.
    END IF

    IF ( PRESENT(interp_method) ) THEN
       IF ( TRIM(interp_method) .NE. 'conserve_order1' ) THEN 
          ! <ERROR STATUS="FATAL">
          !   when registering module/output_field <module_name>/<field_name> then optional
          !   argument interp_method = <interp_method>, but it should be "conserve_order1"
          ! </ERROR>
          CALL error_mesg ('diag_manager_mod::register_diag_field',&
               & 'when registering module/output_field '//TRIM(module_name)//'/'//&
               & TRIM(field_name)//', the optional argument interp_method = '//TRIM(interp_method)//&
               & ', but it should be "conserve_order1"', FATAL)
       END IF
       input_fields(field)%interp_method = TRIM(interp_method)
    ELSE 
       input_fields(field)%interp_method = ''
    END IF

    ! Store the axis info
    num_axes = SIZE(axes(:)) ! num_axes should be <= 3.
    input_fields(field)%axes(1:num_axes) = axes
    input_fields(field)%num_axes = num_axes
    
    siz = 1
    DO j = 1, num_axes
       IF ( axes(j) .LE. 0 ) THEN
          ! <ERROR STATUS="FATAL">
          !   module/output_field <module_name>/<field_name> has non-positive axis_id
          ! </ERROR>
          CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '//trim(module_name)//'/'//&
               & TRIM(field_name)//' has non-positive axis_id', FATAL) 
       END IF
       siz(j) = get_axis_length(axes(j))
    END DO

    ! Default length for axes is 1
    DO j = 1, 3
       input_fields(field)%size(j) = siz(j)
    END DO

    local_siz = 1
    local_start = 1
    local_end= 1
    ! Need to loop through all output_fields associated and allocate their buffers
    DO j = 1, input_fields(field)%num_output_fields
       out_num = input_fields(field)%output_fields(j)
       ! Range is required when pack >= 4 
       IF ( output_fields(out_num)%pack>=4 .AND. .NOT.input_fields(field)%range_present ) THEN
          IF(mpp_pe() .EQ. mpp_root_pe()) THEN
             ! <ERROR STATUS="FATAL">
             !   output_field <field_name> has pack >= 4, range is REQUIRED in register_diag_field
             ! </ERROR>
             CALL error_mesg ('diag_manager_mod::register_diag_field ', 'output_field '//TRIM(field_name)// &
                  ' has pack >=4, range is REQUIRED in register_diag_field', FATAL)
          END IF
       END IF
       ! reset the number of diurnal samples to 1 if the field is static (and, therefore,
       ! doesn't vary diurnally)
       IF ( .NOT.dynamic1 ) output_fields(out_num)%n_diurnal_samples = 1
       !  if local_output (size of output_fields does NOT equal size of input_fields)
       IF ( output_fields(out_num)%reduced_k_range ) THEN
          CALL get_subfield_vert_size(axes, out_num)
   
          local_start(3) = output_fields(out_num)%output_grid%l_start_indx(3)
          local_end(3)   = output_fields(out_num)%output_grid%l_end_indx(3)
          local_siz(3)   = local_end(3) - local_start(3) +1                         

          ALLOCATE(output_fields(out_num)%buffer(siz(1), siz(2), local_siz(3),&
               & output_fields(out_num)%n_diurnal_samples))

          IF ( output_fields(out_num)%time_max ) THEN
             output_fields(out_num)%buffer = MAX_VALUE
          ELSE IF ( output_fields(out_num)%time_min ) THEN
             output_fields(out_num)%buffer = MIN_VALUE
          ELSE
             output_fields(out_num)%buffer = EMPTY
          END IF
          output_fields(out_num)%region_elements = siz(1)*siz(2)*local_siz(3)
          output_fields(out_num)%total_elements  = siz(1)*siz(2)*siz(3)
       ELSE IF ( output_fields(out_num)%local_output ) THEN
          IF ( SIZE(axes(:)) .LE. 1 ) THEN
             ! <ERROR STATUS="FATAL">axes of <field_name> must >= 2 for local output</ERROR>
             CALL error_mesg ('diag_manager_mod::register_diag_field', 'axes of '//TRIM(field_name)//&
                  & ' must >= 2 for local output', FATAL)
          END IF
          CALL get_subfield_size(axes, out_num)
          IF ( output_fields(out_num)%need_compute ) THEN
             DO k = 1, num_axes
                local_start(k) = output_fields(out_num)%output_grid%l_start_indx(k)
                local_end(k) = output_fields(out_num)%output_grid%l_end_indx(k)
                local_siz(k) = local_end(k) - local_start(k) +1                         
             END DO
             ALLOCATE(output_fields(out_num)%buffer(local_siz(1), local_siz(2), local_siz(3),&
                  & output_fields(out_num)%n_diurnal_samples))
             IF(output_fields(out_num)%time_max) THEN
                output_fields(out_num)%buffer = MAX_VALUE
             ELSE IF(output_fields(out_num)%time_min) THEN
                output_fields(out_num)%buffer = MIN_VALUE
             ELSE
                output_fields(out_num)%buffer = EMPTY
             END IF
             output_fields(out_num)%region_elements = local_siz(1)*local_siz(2)*local_siz(3)
             output_fields(out_num)%total_elements = siz(1)*siz(2)*siz(3)
             files(output_fields(out_num)%output_file)%local = .true.
          END IF
       ELSE ! the field is output globally
          ! size of output_fields equal size of input_fields 
          ALLOCATE(output_fields(out_num)%buffer(siz(1), siz(2), siz(3),&
               & output_fields(out_num)%n_diurnal_samples))
          IF(output_fields(out_num)%time_max) THEN
             output_fields(out_num)%buffer = MAX_VALUE
          ELSE IF(output_fields(out_num)%time_min) THEN
             output_fields(out_num)%buffer = MIN_VALUE
          ELSE
             output_fields(out_num)%buffer = EMPTY
          END IF
          output_fields(out_num)%total_elements = siz(1)*siz(2)*siz(3)
       END IF
  
       ! Reset to false in register_field if this is not static
       output_fields(out_num)%static = .TRUE.
       ! check if time average is true for static field
       IF ( .NOT.dynamic1 .AND. output_fields(out_num)%time_ops ) THEN
          WRITE (msg,'(a,"/",a)') TRIM(module_name), TRIM(field_name)
          IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
             ! <ERROR STATUS="WARNING">
             !   module/field <module_name>/<field_name> is STATIC.
             !   Cannot perform time operations average, maximum or
             !   minimum on static fields.  Setting the time operation to 'NONE'
             !   for this field.
             ! </ERROR>
             CALL error_mesg ('diag_manager_mod::register_static_field',&
                  & 'module/field '//TRIM(msg)//' is STATIC.  Cannot perform time operations&
                  & average, maximum, or minimum on static fields.  Setting the time operation&
                  & to "NONE" for this field.', WARNING)
          END IF
          output_fields(out_num)%time_ops = .FALSE.
          output_fields(out_num)%time_average = .FALSE.
          output_fields(out_num)%time_method = 'point'
       END IF

       ! assume that the number of axes of output_fields = that of input_fields
       ! this should be changed later to take into account time-of-day axis
       output_fields(out_num)%num_axes = input_fields(field)%num_axes
       ! Axes are copied from input_fields if output globally or from subaxes if output locally
       IF ( .NOT.output_fields(out_num)%local_output ) THEN 
          output_fields(out_num)%axes(1:input_fields(field)%num_axes) =&
               & input_fields(field)%axes(1:input_fields(field)%num_axes)
       ELSE
          output_fields(out_num)%axes(1:input_fields(field)%num_axes) =&
               & output_fields(out_num)%output_grid%subaxes(1:input_fields(field)%num_axes)
       END IF

       ! if necessary, initialize the diurnal time axis and append its index in the 
       ! output field axes array
       IF ( output_fields(out_num)%n_diurnal_samples > 1 ) THEN
          output_fields(out_num)%axes(output_fields(out_num)%num_axes+1) =&
               & init_diurnal_axis(output_fields(out_num)%n_diurnal_samples)
          output_fields(out_num)%num_axes = output_fields(out_num)%num_axes+1
       END IF

       IF ( output_fields(out_num)%reduced_k_range ) THEN 
          output_fields(out_num)%axes(3) = output_fields(out_num)%output_grid%subaxes(3)
       END IF

       ! Initialize a time variable used in an error check
       output_fields(out_num)%Time_of_prev_field_data = Time_zero
    END DO

    IF ( input_fields(field)%mask_variant ) THEN
       DO j = 1, input_fields(field)%num_output_fields
          out_num = input_fields(field)%output_fields(j)
          IF(output_fields(out_num)%time_average) THEN
             ALLOCATE(output_fields(out_num)%counter(siz(1), siz(2), siz(3),&
                  & output_fields(out_num)%n_diurnal_samples))
             output_fields(out_num)%counter = 0.0
          END IF
       END DO
    END IF
  END FUNCTION register_static_field
  ! </FUNCTION>

  ! <FUNCTION NAME="send_data_0d" INTERFACE="send_data">
  !   <IN NAME="diag_field_id" TYPE="INTEGER"> </IN>
  !   <IN NAME="field" TYPE="REAL"> </IN>
  !   <IN NAME="time" TYPE="TYPE(time_type), OPTIONAL"> </IN>
  !   <IN NAME="is_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="js_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="ks_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="mask" TYPE="Not Applicable"></IN>
  !   <IN NAME="rmask" TYPE="Not Applicable"></IN>
  !   <IN NAME="ie_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="je_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="ke_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="weight" TYPE="Not Applicable"></IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL"></OUT>
  LOGICAL FUNCTION send_data_0d(diag_field_id, field, time, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, INTENT(in) :: field
    TYPE(time_type), INTENT(in), OPTIONAL :: time
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL :: field_out(1, 1, 1)

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_0d = .FALSE.
       RETURN
    END IF
    ! First copy the data to a three d array with last element 1
    field_out(1, 1, 1) = field
    send_data_0d = send_data_3d(diag_field_id, field_out, time, err_msg=err_msg)
  END FUNCTION send_data_0d
  ! </FUNCTION>

  ! <FUNCTION NAME="send_data_1d" INTERFACE="send_data">
  !   <IN NAME="diag_field_id" TYPE="INTEGER"> </IN>
  !   <IN NAME="field" TYPE="REAL, DIMENSION(:)"> </IN>
  !   <IN NAME="time" TYPE="TYPE(time_type)"> </IN>
  !   <IN NAME="is_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="js_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="ks_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="mask" TYPE="LOGICAL, DIMENSION(:), OPTIONAL"></IN>
  !   <IN NAME="rmask" TYPE="REAL, DIMENSION(:), OPTIONAL"></IN>
  !   <IN NAME="ie_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="je_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="ke_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="weight" TYPE="REAL, OPTIONAL"></IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL"></OUT>
  LOGICAL FUNCTION send_data_1d(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:), INTENT(in) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    REAL, INTENT(in), DIMENSION(:), OPTIONAL :: rmask
    TYPE (time_type), INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, ie_in
    LOGICAL, INTENT(in), DIMENSION(:), OPTIONAL :: mask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL, DIMENSION(SIZE(field(:)), 1, 1) :: field_out
    LOGICAL, DIMENSION(SIZE(field(:)), 1, 1) ::  mask_out

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_1d = .FALSE.
       RETURN
    END IF

    ! First copy the data to a three d array with last element 1
    field_out(:, 1, 1) = field

    ! Default values for mask
    IF ( PRESENT(mask) ) THEN 
       mask_out(:, 1, 1) = mask
    ELSE
       mask_out = .TRUE.
    END IF

    IF ( PRESENT(rmask) ) WHERE (rmask < 0.5) mask_out(:, 1, 1) = .FALSE.
    IF ( PRESENT(mask) .OR. PRESENT(rmask) ) THEN
       IF ( PRESENT(is_in) .OR. PRESENT(ie_in) ) THEN
          send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=1, ks_in=1,&
               & mask=mask_out, ie_in=ie_in, je_in=1, ke_in=1, weight=weight, err_msg=err_msg)
       ELSE
          send_data_1d = send_data_3d(diag_field_id, field_out, time, mask=mask_out,&
               & weight=weight, err_msg=err_msg)
       END IF
    ELSE
       IF ( PRESENT(is_in) .OR. PRESENT(ie_in) ) THEN
          send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=1, ks_in=1,&
               & ie_in=ie_in, je_in=1, ke_in=1, weight=weight, err_msg=err_msg)
       ELSE
          send_data_1d = send_data_3d(diag_field_id, field_out, time, weight=weight, err_msg=err_msg)
       END IF
    END IF
  END FUNCTION send_data_1d
  ! </FUNCTION>

  ! <FUNCTION NAME="send_data_2d" INTERFACE="send_data">
  !   <IN NAME="diag_field_id" TYPE="INTEGER"> </IN>
  !   <IN NAME="field" TYPE="REAL, DIMENSION(:,:)"> </IN>
  !   <IN NAME="time" TYPE="TYPE(time_type)"> </IN>
  !   <IN NAME="is_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="js_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="ks_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="mask" TYPE="LOGICAL, DIMENSION(:,:), OPTIONAL"></IN>
  !   <IN NAME="rmask" TYPE="REAL, DIMENSION(:,:), OPTIONAL"></IN>
  !   <IN NAME="ie_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="je_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="ke_in" TYPE="Not Applicable"></IN>
  !   <IN NAME="weight" TYPE="REAL, OPTIONAL"></IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL"></OUT>
  LOGICAL FUNCTION send_data_2d(diag_field_id, field, time, is_in, js_in, &
       & mask, rmask, ie_in, je_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, INTENT(in), DIMENSION(:,:) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    TYPE (time_type), INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ie_in, je_in
    LOGICAL, INTENT(in), DIMENSION(:,:), OPTIONAL :: mask
    REAL, INTENT(in), DIMENSION(:,:),OPTIONAL :: rmask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2),1) :: field_out
    LOGICAL, DIMENSION(SIZE(field,1),SIZE(field,2),1) ::  mask_out

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_2d = .FALSE.
       RETURN
    END IF

    ! First copy the data to a three d array with last element 1
    field_out(:, :, 1) = field

    ! Default values for mask
    IF ( PRESENT(mask) ) THEN 
       mask_out(:, :, 1) = mask
    ELSE 
       mask_out = .TRUE.
    END IF
    
    IF ( PRESENT(rmask) ) WHERE ( rmask < 0.5 ) mask_out(:, :, 1) = .FALSE.
    IF ( PRESENT(mask) .OR. PRESENT(rmask) ) THEN
       send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=js_in, ks_in=1, mask=mask_out,&
            & ie_in=ie_in, je_in=je_in, ke_in=1, weight=weight, err_msg=err_msg)
    ELSE
       send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=js_in, ks_in=1,&
            & ie_in=ie_in, je_in=je_in, ke_in=1, weight=weight, err_msg=err_msg)
    END IF
  END FUNCTION send_data_2d
  ! </FUNCTION>

  ! <FUNCTION NAME="send_data_3d" INTERFACE="send_data">
  !   <IN NAME="diag_field_id" TYPE="INTEGER"> </IN>
  !   <IN NAME="field" TYPE="REAL, DIMENSION(:,:,:)"> </IN>
  !   <IN NAME="time" TYPE="TYPE(time_type)"> </IN>
  !   <IN NAME="is_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="js_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="ks_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="mask" TYPE="LOGICAL, DIMENSION(:,:,:), OPTIONAL"></IN>
  !   <IN NAME="rmask" TYPE="REAL, DIMENSION(:,:,:), OPTIONAL"></IN>
  !   <IN NAME="ie_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="je_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="ke_in" TYPE="INTEGER, OPTIONAL"></IN>
  !   <IN NAME="weight" TYPE="REAL, OPTIONAL"></IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL"></OUT>
  LOGICAL FUNCTION send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
             & mask, rmask, ie_in, je_in, ke_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:,:,:), INTENT(in) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    TYPE (time_type), INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ks_in,ie_in,je_in, ke_in 
    LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
    REAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: rmask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL :: weight1
    REAL :: missvalue
    INTEGER :: ksr, ker
    INTEGER :: i, out_num, file_num, n1, n2, n3, number_of_outputs, ii,f1,f2,f3,f4
    INTEGER :: freq, units, is, js, ks, ie, je, ke, i1, j1,k1, j, k
    INTEGER, DIMENSION(3) :: l_start, l_end ! local start and end indices on 3 axes for regional output
    INTEGER   :: hi, hj, twohi, twohj  ! halo size in x and y direction
    INTEGER :: sample ! index along the diurnal time axis
    INTEGER :: day,second,tick ! components of the current date
    INTEGER :: status
    INTEGER :: numthreads
#if defined(_OPENMP)
    INTEGER :: omp_get_num_threads !< OMP function
#endif
    LOGICAL :: average, phys_window, need_compute
    LOGICAL :: reduced_k_range, local_output
    LOGICAL :: time_max, time_min
    LOGICAL :: missvalue_present
    LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask
    CHARACTER(len=256) :: err_msg_local
    CHARACTER(len=128) :: error_string, error_string1
    TYPE(time_type) :: dt ! time interval for diurnal output

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_3d = .FALSE.
       RETURN
    ELSE
       send_data_3d = .TRUE.
    END IF

    IF ( PRESENT(err_msg) ) err_msg = ''
    IF ( .NOT.module_is_initialized ) THEN
       IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'diag_manager NOT initialized', err_msg) ) RETURN
    END IF
    err_msg_local = ''

    ! oor_mask is only used for checking out of range values.
    ALLOCATE(oor_mask(SIZE(field,1),SIZE(field,2),SIZE(field,3)), STAT=status)
    IF ( status .NE. 0 ) THEN
       WRITE (err_msg_local, FMT='("Unable to allocate oor_mask(",I5,",",I5,",",I5"). (STAT: ",I5,")")')&
            & SIZE(field,1), SIZE(field,2), SIZE(field,3), status
       IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) RETURN
    END IF

    IF ( PRESENT(mask) ) THEN 
       oor_mask = mask
    ELSE 
       oor_mask = .TRUE.
    END IF
    IF ( PRESENT(rmask) ) WHERE ( rmask < 0.5 ) oor_mask = .FALSE.

    ! send_data works in either one or another of two modes.
    ! 1. Input field is a window (e.g. FMS physics)
    ! 2. Input field includes halo data
    ! It cannot handle a window of data that has halos.
    ! (A field with no windows or halos can be thought of as a special case of either mode.)
    ! The logic for indexing is quite different for these two modes, but is not clearly separated.
    ! If both the beggining and ending indices are present, then field is assumed to have halos.
    ! If only beggining indices are present, then field is assumed to be a window.

    ! There are a number of ways a user could mess up this logic, depending on the combination
    ! of presence/absence of is,ie,js,je. The checks below should catch improper combinations.
    IF ( PRESENT(ie_in) ) THEN
       IF ( .NOT.PRESENT(is_in) ) THEN
          IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'ie_in present without is_in', err_msg) ) THEN 
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF
       IF ( PRESENT(js_in) .AND. .NOT.PRESENT(je_in) ) THEN
          IF ( fms_error_handler('diag_manager_modsend_data_3d',&
               & 'is_in and ie_in present, but js_in present without je_in', err_msg) ) THEN 
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF
    END IF
    IF ( PRESENT(je_in) ) THEN
       IF ( .NOT.PRESENT(js_in) ) THEN
          IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'je_in present without js_in', err_msg) ) THEN 
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF
       IF ( PRESENT(is_in) .AND. .NOT.PRESENT(ie_in) ) THEN
          IF ( fms_error_handler('diag_manager_mod::send_data_3d',&
               & 'js_in and je_in present, but is_in present without ie_in', err_msg)) THEN 
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF
    END IF

    ! If is, js, or ks not present default them to 1
    is = 1
    js = 1
    ks = 1
    IF ( PRESENT(is_in) ) is = is_in
    IF ( PRESENT(js_in) ) js = js_in
    IF ( PRESENT(ks_in) ) ks = ks_in
    n1 = SIZE(field, 1)
    n2 = SIZE(field, 2)
    n3 = SIZE(field, 3)
    ie = is+n1-1
    je = js+n2-1
    ke = ks+n3-1
    IF ( PRESENT(ie_in) ) ie = ie_in
    IF ( PRESENT(je_in) ) je = je_in
    IF ( PRESENT(ke_in) ) ke = ke_in
    twohi = n1-(ie-is+1)
    IF ( MOD(twohi,2) /= 0 ) THEN
       IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'non-symmetric halos in first dimension', err_msg) ) THEN 
          DEALLOCATE(oor_mask)
          RETURN
       END IF
    END IF
    twohj = n2-(je-js+1)
    IF ( MOD(twohj,2) /= 0 ) THEN
       IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'non-symmetric halos in second dimension', err_msg) ) THEN 
          DEALLOCATE(oor_mask)
          RETURN
       END IF
    END IF
    hi = twohi/2
    hj = twohj/2

    ! The next line is necessary to ensure that is,ie,js,ie are relative to field(1:,1:)
    ! But this works only when there is no windowing.
    IF ( PRESENT(ie_in) .AND. PRESENT(je_in) ) THEN
       is=1+hi
       ie=n1-hi
       js=1+hj
       je=n2-hj
    END IF

    ! used for field, mask and rmask bounds
    f1=1+hi
    f2=n1-hi
    f3=1+hj
    f4=n2-hj

    ! weight is for time averaging where each time level may has a different weight
    IF ( PRESENT(weight) ) THEN 
       weight1 = weight
    ELSE
       weight1 = 1.
    END IF
    
    ! Is there a missing_value?
    missvalue_present = input_fields(diag_field_id)%missing_value_present
    IF ( missvalue_present ) missvalue = input_fields(diag_field_id)%missing_value

    number_of_outputs = input_fields(diag_field_id)%num_output_fields
!$OMP CRITICAL
    input_fields(diag_field_id)%numthreads = 1
#if defined(_OPENMP)
    input_fields(diag_field_id)%numthreads = omp_get_num_threads()
#endif
    numthreads = input_fields(diag_field_id)%numthreads
!$OMP END CRITICAL

    if(present(time)) input_fields(diag_field_id)%time = time

    ! Issue a warning if any value in field is outside the valid range
    IF ( input_fields(diag_field_id)%range_present ) THEN
       IF ( ISSUE_OOR_WARNINGS .OR. OOR_WARNINGS_FATAL ) THEN
          WRITE (error_string, '("[",ES14.5E3,",",ES14.5E3,"]")')&
               & input_fields(diag_field_id)%range(1:2)
          WRITE (error_string1, '("(Min: ",ES14.5E3,", Max: ",ES14.5E3, ")")')&
                  & MINVAL(field(f1:f2,f3:f4,ks:ke),MASK=oor_mask(f1:f2,f3:f4,ks:ke)),&
                  & MAXVAL(field(f1:f2,f3:f4,ks:ke),MASK=oor_mask(f1:f2,f3:f4,ks:ke))
          IF ( missvalue_present ) THEN
             IF ( ANY(oor_mask(f1:f2,f3:f4,ks:ke) .AND.&
                  &   ((field(f1:f2,f3:f4,ks:ke) < input_fields(diag_field_id)%range(1) .OR.&
                  &     field(f1:f2,f3:f4,ks:ke) > input_fields(diag_field_id)%range(2)).AND.&
                  &     field(f1:f2,f3:f4,ks:ke) .NE. missvalue)) ) THEN
                ! <ERROR STATUS="WARNING/FATAL">
                !   A value for <module_name> in field <field_name> (Min: <min_val>, Max: <max_val>)
                !   is outside the range [<lower_val>,<upper_val>] and not equal to the missing
                !   value.
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::send_data_3d',&
                     & 'A value for '//&
                     &TRIM(input_fields(diag_field_id)%module_name)//' in field '//&
                     &TRIM(input_fields(diag_field_id)%field_name)//' '&
                     &//TRIM(error_string1)//&
                     &' is outside the range '//TRIM(error_string)//',&
                     & and not equal to the missing value.',&
                     &OOR_WARNING)
             END IF
          ELSE
             IF ( ANY(oor_mask(f1:f2,f3:f4,ks:ke) .AND.&
                  &   (field(f1:f2,f3:f4,ks:ke) < input_fields(diag_field_id)%range(1) .OR.&
                  &    field(f1:f2,f3:f4,ks:ke) > input_fields(diag_field_id)%range(2))) ) THEN
                ! <ERROR STATUS="WARNING/FATAL">
                !   A value for <module_name> in field <field_name> (Min: <min_val>, Max: <max_val>)
                !   is outside the range [<lower_val>,<upper_val>].
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::send_data_3d',&
                     & 'A value for '//&
                     &TRIM(input_fields(diag_field_id)%module_name)//' in field '//&
                     &TRIM(input_fields(diag_field_id)%field_name)//' '&
                     &//TRIM(error_string1)//&
                     &' is outside the range '//TRIM(error_string)//'.',&
                     &OOR_WARNING)
             END IF
          END IF
       END IF
    END IF

    ! Loop through each output field that depends on this input field
    num_out_fields: DO ii = 1, number_of_outputs
       ! Get index to an output field
       out_num = input_fields(diag_field_id)%output_fields(ii)

      ! is this field output on a local domain only?
       local_output = output_fields(out_num)%local_output
       ! if local_output, does the current PE take part in send_data?
       need_compute = output_fields(out_num)%need_compute

       reduced_k_range = output_fields(out_num)%reduced_k_range

      ! skip all PEs not participating in outputting this field
       IF ( local_output .AND. (.NOT.need_compute) ) CYCLE

       ! Get index to output file for this field
       file_num = output_fields(out_num)%output_file
       IF(file_num == max_files) CYCLE
       ! Output frequency and units for this file is
       freq = files(file_num)%output_freq
       units = files(file_num)%output_units
       ! Is this output field being time averaged?
       average = output_fields(out_num)%time_average
       ! Looking for max and min value of this field over the sampling interval?
       time_max = output_fields(out_num)%time_max
       time_min = output_fields(out_num)%time_min   
       IF ( output_fields(out_num)%total_elements > SIZE(field(f1:f2,f3:f4,ks:ke)) ) THEN
          output_fields(out_num)%phys_window = .TRUE.
       ELSE
          output_fields(out_num)%phys_window = .FALSE.
       END IF
       phys_window = output_fields(out_num)%phys_window
       IF ( need_compute ) THEN        
          l_start = output_fields(out_num)%output_grid%l_start_indx
          l_end = output_fields(out_num)%output_grid%l_end_indx
       END IF

       ! compute the diurnal index
       sample = 1
       IF ( PRESENT(time) ) THEN
          dt = set_time(0,1)/output_fields(out_num)%n_diurnal_samples ! our time interval
          CALL get_time(time,second,day,tick) ! current date
          sample = set_time(second,0,tick)/dt + 1
       END IF
   
       ! Get the vertical layer start and end index.
       IF ( reduced_k_range ) THEN
          l_start(3) = output_fields(out_num)%output_grid%l_start_indx(3)
          l_end(3) = output_fields(out_num)%output_grid%l_end_indx(3)
       END IF
       ksr= l_start(3)
       ker= l_end(3)

       ! Initialize output time for fields output every time step
       IF ( freq == EVERY_TIME .AND. .NOT.output_fields(out_num)%static ) THEN
          IF (output_fields(out_num)%next_output == output_fields(out_num)%last_output) THEN
             IF(PRESENT(time)) THEN
                output_fields(out_num)%next_output = time
             ELSE
                WRITE (error_string,'(a,"/",a)')&
                     & TRIM(input_fields(diag_field_id)%module_name),&
                     & TRIM(output_fields(out_num)%output_name)
                IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                     & ', time must be present when output frequency = EVERY_TIME', err_msg)) THEN 
                   DEALLOCATE(oor_mask)
                   RETURN
                END IF
             END IF
          END IF
       END IF 
       IF ( .NOT.output_fields(out_num)%static .AND. .NOT.PRESENT(time) ) THEN
          WRITE (error_string,'(a,"/",a)')&
               & TRIM(input_fields(diag_field_id)%module_name), &
               & TRIM(output_fields(out_num)%output_name)
          IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
               & ', time must be present for nonstatic field', err_msg)) THEN 
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF

       ! Is it time to output for this field; CAREFUL ABOUT > vs >= HERE
       !--- The fields send out within openmp parallel region will be written out in
       !--- diag_send_complete. 
       IF ( numthreads == 1) then
          IF ( .NOT.output_fields(out_num)%static .AND. freq /= END_OF_RUN ) THEN
             IF ( time > output_fields(out_num)%next_output ) THEN
                ! A non-static field that has skipped a time level is an error
                IF ( time > output_fields(out_num)%next_next_output .AND. freq > 0 ) THEN
                   IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN 
                      WRITE (error_string,'(a,"/",a)')&
                           & TRIM(input_fields(diag_field_id)%module_name), &
                           & TRIM(output_fields(out_num)%output_name)
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                           & ' is skipped one time level in output data', err_msg)) THEN 
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF

                status = writing_field(out_num, .FALSE., error_string, time)
                IF(status == -1) THEN
                   IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
                      IF(fms_error_handler('diag_manager_mod::send_data_3d','module/output_field '//TRIM(error_string)//&
                           & ', write EMPTY buffer', err_msg)) THEN 
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
             END IF  !time > output_fields(out_num)%next_output
          END IF  !.not.output_fields(out_num)%static .and. freq /= END_OF_RUN
          ! Finished output of previously buffered data, now deal with buffering new data   
       END IF

       IF ( .NOT.output_fields(out_num)%static .AND. .NOT.need_compute .AND. debug_diag_manager ) THEN
          CALL check_bounds_are_exact_dynamic(out_num, diag_field_id, Time, err_msg=err_msg_local)
          IF ( err_msg_local /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                DEALLOCATE(oor_mask)
                RETURN
             END IF
          END IF
       END IF
 
       ! Take care of submitted field data
       IF ( average ) THEN
          IF ( input_fields(diag_field_id)%mask_variant ) THEN
             IF ( need_compute ) THEN
                WRITE (error_string,'(a,"/",a)')  &
                     & TRIM(input_fields(diag_field_id)%module_name), &
                     & TRIM(output_fields(out_num)%output_name)   
                IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                     & ', regional output NOT supported with mask_variant', err_msg)) THEN 
                   DEALLOCATE(oor_mask)
                   RETURN
                END IF
             END IF

             ! Should reduced_k_range data be supported with the mask_variant option   ?????
             ! If not, error message should be produced and the reduced_k_range loop below eliminated 
             IF ( PRESENT(mask) ) THEN
                IF ( missvalue_present ) THEN              
                   IF ( debug_diag_manager ) THEN
                      CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                      CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                      IF ( err_msg_local /= '' ) THEN
                         IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                            DEALLOCATE(oor_mask)
                            RETURN
                         END IF
                      END IF
                   END IF
                   IF( numthreads>1 .AND. phys_window ) then
                      IF ( reduced_k_range ) THEN 
                         DO k= ksr, ker
                            k1= k - ksr + 1 
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                          & field(i-is+1+hi, j-js+1+hj, k) * weight1  
                                     output_fields(out_num)%counter(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%counter(i-hi,j-hj,k1,sample) + weight1
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
                         DO k=ks, ke 
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k)*weight1
                                     output_fields(out_num)%counter(i-hi,j-hj,k,sample) =&
                                          &output_fields(out_num)%counter(i-hi,j-hj,k,sample) + weight1
                                  END IF
                               END DO
                            END DO
                         END DO
                      END IF
                   ELSE
!$OMP CRITICAL
                      IF ( reduced_k_range ) THEN 
                         DO k= ksr, ker
                            k1= k - ksr + 1 
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                          & field(i-is+1+hi, j-js+1+hj, k) * weight1  
                                     output_fields(out_num)%counter(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%counter(i-hi,j-hj,k1,sample) + weight1
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
                         DO k=ks, ke 
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k)*weight1
                                     output_fields(out_num)%counter(i-hi,j-hj,k,sample) =&
                                          &output_fields(out_num)%counter(i-hi,j-hj,k,sample) + weight1
                                  END IF
                               END DO
                            END DO
                         END DO
                      END IF
!$OMP END CRITICAL
                   END IF
                ELSE
                   WRITE (error_string,'(a,"/",a)')&
                        & TRIM(input_fields(diag_field_id)%module_name), &
                        & TRIM(output_fields(out_num)%output_name)
                   IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                        & ', variable mask but no missing value defined', err_msg)) THEN 
                      DEALLOCATE(oor_mask)
                      RETURN
                   END IF
                END IF
             ELSE  ! no mask present
                WRITE (error_string,'(a,"/",a)')&
                     & TRIM(input_fields(diag_field_id)%module_name), &
                     & TRIM(output_fields(out_num)%output_name)
                IF(fms_error_handler('diag_manager_mod::send_data_3d','module/output_field '//TRIM(error_string)//&
                     & ', variable mask but no mask given', err_msg)) THEN 
                   DEALLOCATE(oor_mask)
                   RETURN
                END IF
             END IF
          ELSE ! mask_variant=false
             IF ( PRESENT(mask) ) THEN
                IF ( missvalue_present ) THEN
                   IF ( need_compute ) THEN
                      IF (numthreads>1 .AND. phys_window) then
                         DO k = l_start(3), l_end(3)
                            k1 = k-l_start(3)+1
                            DO j = js, je 
                               DO i = is, ie
                                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                     i1 = i-l_start(1)-hi+1 
                                     j1=  j-l_start(2)-hj+1
                                     IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                             & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1                              
                                     ELSE
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue   
                                     END IF
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k = l_start(3), l_end(3)
                            k1 = k-l_start(3)+1
                            DO j = js, je 
                               DO i = is, ie
                                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                     i1 = i-l_start(1)-hi+1 
                                     j1=  j-l_start(2)-hj+1
                                     IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                             & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1                              
                                     ELSE
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue   
                                     END IF
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      ENDIF
!$OMP CRITICAL
                      DO j = js, je 
                         DO i = is, ie
                            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                               output_fields(out_num)%num_elements(sample) = &
                                    output_fields(out_num)%num_elements(sample) + l_end(3) - l_start(3) + 1
                            END IF
                         END DO
                      END DO
!$OMP END CRITICAL
                   ELSE IF ( reduced_k_range ) THEN 
                      IF (numthreads>1 .AND. phys_window) then
                         DO k=ksr, ker
                            k1 = k - ksr + 1 
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k) * weight1  
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k=ksr, ker
                            k1 = k - ksr + 1 
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k) * weight1  
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
                   ELSE
                      IF ( debug_diag_manager ) THEN
                         CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                         CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                         IF ( err_msg_local /= '' ) THEN
                            IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                               DEALLOCATE(oor_mask)
                               RETURN
                            END IF
                         END IF
                      END IF
                      IF (numthreads>1 .AND. phys_window) then
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k) * weight1  
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k) * weight1  
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                                  END IF
                               END DO
                            END DO
                         END DO 
!$OMP END CRITICAL                     
                      END IF
                   END IF
!$OMP CRITICAL
                   IF ( need_compute .AND. .NOT.phys_window ) THEN
                      IF ( ANY(mask(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3))) ) &
                           & output_fields(out_num)%count_0d(sample) =&
                           & output_fields(out_num)%count_0d(sample) + weight1
                   ELSE
                      IF ( ANY(mask(f1:f2,f3:f4,ks:ke)) ) output_fields(out_num)%count_0d(sample) =&
                           & output_fields(out_num)%count_0d(sample)+weight1
                   END IF
!$OMP END CRITICAL 

                ELSE ! missing value NOT present
                   IF (   (.NOT.ALL(mask(f1:f2,f3:f4,ks:ke)) .AND. mpp_pe() .EQ. mpp_root_pe()).AND.&
                        &  .NOT.input_fields(diag_field_id)%issued_mask_ignore_warning ) THEN
                      ! <ERROR STATUS="WARNING">
                      !   Mask will be ignored since missing values were not specified for field <field_name>
                      !   in module <module_name>
                      ! </ERROR>
                      CALL error_mesg('diag_manager_mod::send_data_3d',&
                           & 'Mask will be ignored since missing values were not specified for field '//&
                           & trim(input_fields(diag_field_id)%field_name)//' in module '//&
                           & trim(input_fields(diag_field_id)%module_name), WARNING)
                      input_fields(diag_field_id)%issued_mask_ignore_warning = .TRUE.
                   END IF
                   IF ( need_compute ) THEN                 
                      IF (numthreads>1 .AND. phys_window) then
                         DO j = js, je 
                            DO i = is, ie
                               IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                  i1 = i-l_start(1)-hi+1 
                                  j1 =  j-l_start(2)-hj+1
                                  output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                       & field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                               END IF
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO j = js, je 
                            DO i = is, ie
                               IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                  i1 = i-l_start(1)-hi+1 
                                  j1 =  j-l_start(2)-hj+1
                                  output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                       & field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                               END IF
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
!$OMP CRITICAL
                      DO j = js, je 
                         DO i = is, ie
                            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                               output_fields(out_num)%num_elements(sample)=&
                                    & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1

                            END IF
                         END DO
                      END DO
!$OMP END CRITICAL
                   ELSE IF ( reduced_k_range ) THEN 
                      IF (numthreads>1 .AND. phys_window) then
                         ksr= l_start(3)
                         ker= l_end(3)
                         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                              & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample)&
                              & + field(f1:f2,f3:f4,ksr:ker)*weight1
                      ELSE 
!$OMP CRITICAL  
                         ksr= l_start(3)
                         ker= l_end(3)
                         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                              & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample)&
                              & + field(f1:f2,f3:f4,ksr:ker)*weight1
!$OMP END CRITICAL
                      END IF
                   ELSE
                      IF ( debug_diag_manager ) THEN
                         CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                         CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                         IF ( err_msg_local /= '') THEN
                            IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                               DEALLOCATE(oor_mask)
                               RETURN
                            END IF
                         END IF
                      END IF
                      IF (numthreads>1 .AND. phys_window) then
                         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                              & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample)&
                              & + field(f1:f2,f3:f4,ks:ke)*weight1 
                      ELSE
!$OMP CRITICAL  
                         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                              & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample)&
                              & + field(f1:f2,f3:f4,ks:ke)*weight1 
!$OMP END CRITICAL 
                      END IF
                   END IF
!$OMP CRITICAL 
                   IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
                        & output_fields(out_num)%count_0d(sample) + weight1
!$OMP END CRITICAL 
                END IF
             ELSE ! mask NOT present
                IF ( missvalue_present ) THEN
                   IF ( need_compute ) THEN 
                      if( numthreads>1 .AND. phys_window ) then
                         DO k = l_start(3), l_end(3)
                            k1 = k - l_start(3) + 1                    
                            DO j = js, je 
                               DO i = is, ie
                                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj) THEN
                                     i1 = i-l_start(1)-hi+1 
                                     j1=  j-l_start(2)-hj+1 
                                     IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                             & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     ELSE
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                                     END IF
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k = l_start(3), l_end(3)
                            k1 = k - l_start(3) + 1                    
                            DO j = js, je 
                               DO i = is, ie
                                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj) THEN
                                     i1 = i-l_start(1)-hi+1 
                                     j1=  j-l_start(2)-hj+1 
                                     IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                             & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     ELSE
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                                     END IF
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
!$OMP CRITICAL
                      DO j = js, je 
                         DO i = is, ie
                            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj) THEN
                               output_fields(out_num)%num_elements(sample) =&
                                    & output_fields(out_num)%num_elements(sample) + l_end(3) - l_start(3) + 1
                            END IF
                         END DO
                      END DO
                      IF ( .NOT.phys_window ) THEN
                         !rab if(any(field(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3)) /= &
                         !rab        & missvalue)) &
                         !rab        & output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1 
                         outer0: DO k = l_start(3), l_end(3)
                            DO j=l_start(2)+hj, l_end(2)+hj
                               DO i=l_start(1)+hi, l_end(1)+hi
                                  IF ( field(i,j,k) /= missvalue ) THEN
                                     output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) + weight1    
                                     EXIT outer0
                                  END IF
                               END DO
                            END DO
                         END DO outer0
                      END IF
!$OMP END CRITICAL
                   ELSE IF ( reduced_k_range ) THEN 
                      if( numthreads>1 .AND. phys_window ) then
                         ksr= l_start(3)
                         ker= l_end(3)
                         DO k = ksr, ker
                            k1 = k - ksr + 1
                            DO j=js, je
                               DO i=is, ie
                                  IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k) * weight1  
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) = missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
                      else
!$OMP CRITICAL
                         ksr= l_start(3)
                         ker= l_end(3)
                         DO k = ksr, ker
                            k1 = k - ksr + 1
                            DO j=js, je
                               DO i=is, ie
                                  IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k) * weight1  
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) = missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
                      !rab
                      !rab if(any(field(f1:f2,f3:f4,ks:ke) /= missvalue)) &
                      !rab       & output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1    
!$OMP CRITICAL
                      outer3: DO k = ksr, ker
                         k1=k-ksr+1
                         DO j=f3, f4 
                            DO i=f1, f2 
                               IF ( field(i,j,k) /= missvalue ) THEN
                                  output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1    
                                  EXIT outer3
                               END IF
                            END DO
                         END DO
                      END DO outer3
!$OMP END CRITICAL
                   ELSE
                      IF ( debug_diag_manager ) THEN
                         CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                         CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                         IF ( err_msg_local /= '' ) THEN
                            IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                               DEALLOCATE(oor_mask)
                               RETURN
                            END IF
                         END IF
                      END IF
                      IF( numthreads > 1 .AND. phys_window ) then
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k) * weight1  
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,k) * weight1  
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
!$OMP CRITICAL
                      outer1: DO k=ks, ke 
                         DO j=f3, f4 
                            DO i=f1, f2 
                               IF ( field(i,j,k) /= missvalue ) THEN
                                  output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) + weight1    
                                  EXIT outer1
                               END IF
                            END DO
                         END DO
                      END DO outer1
!$OMP END CRITICAL
                   END IF
                ELSE ! no missing value defined, No mask
                   IF ( need_compute ) THEN
                      IF( numthreads > 1 .AND. phys_window ) then
                         DO j = js, je  
                            DO i = is, ie                                         
                               IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                  i1 = i-l_start(1)-hi+1 
                                  j1=  j-l_start(2)-hj+1
                                  output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                       & field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                               END IF
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO j = js, je  
                            DO i = is, ie                                         
                               IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                  i1 = i-l_start(1)-hi+1 
                                  j1=  j-l_start(2)-hj+1
                                  output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                       & field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                               END IF
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF

!$OMP CRITICAL
                      DO j = js, je  
                         DO i = is, ie                                         
                            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                               output_fields(out_num)%num_elements(sample) =&
                                    & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1
                            END IF
                         END DO
                      END DO
!$OMP END CRITICAL
                      ! Accumulate time average 
                   ELSE IF ( reduced_k_range ) THEN 
                      ksr= l_start(3)
                      ker= l_end(3)
                      IF( numthreads > 1 .AND. phys_window ) then
                         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                              & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                              & field(f1:f2,f3:f4,ksr:ker)*weight1
                      ELSE
!$OMP CRITICAL
                         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                              & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                              & field(f1:f2,f3:f4,ksr:ker)*weight1
!$OMP END CRITICAL
                      END IF
                   ELSE 
                      IF ( debug_diag_manager ) THEN
                         CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                         CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                         IF ( err_msg_local /= '' ) THEN
                            IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                               DEALLOCATE(oor_mask)
                               RETURN
                            END IF
                         END IF
                      END IF
                      IF( numthreads > 1 .AND. phys_window ) then
                         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                              & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                              & field(f1:f2,f3:f4,ks:ke)*weight1
                      ELSE
!$OMP CRITICAL
                         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                              & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                              & field(f1:f2,f3:f4,ks:ke)*weight1
!$OMP END CRITICAL
                      END IF
                   END IF
!$OMP CRITICAL
                   IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
                        & output_fields(out_num)%count_0d(sample) + weight1
!$OMP END CRITICAL
                END IF
             END IF ! if mask present
          END IF  !if mask_variant
!$OMP CRITICAL
          IF ( .NOT.need_compute )&
               & output_fields(out_num)%num_elements(sample) =&
               & output_fields(out_num)%num_elements(sample) + (ie-is+1)*(je-js+1)*(ke-ks+1)
          IF ( reduced_k_range ) &
               & output_fields(out_num)%num_elements(sample) = output_fields(out_num)%num_elements(sample) +&
               & (ie-is+1)*(je-js+1)*(ker-ksr+1)
!$OMP END CRITICAL
             ! Add processing for Max and Min
       ELSE IF ( time_max ) THEN
          IF ( PRESENT(mask) ) THEN
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1 
                            j1=  j-l_start(2)-hj+1
                            IF ( mask(i-is+1+hi,j-js+1+hj,k) .AND.&
                                 & field(i-is+1+hi,j-js+1+hj,k)>output_fields(out_num)%buffer(i1,j1,k1,sample)) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Maximum time value with masking 
             ELSE IF ( reduced_k_range ) THEN 
                ksr = l_start(3)
                ker = l_end(3)
                WHERE ( mask(f1:f2,f3:f4,ksr:ker) .AND. &
                     & field(f1:f2,f3:f4,ksr:ker) > output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample))&
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( mask(f1:f2,f3:f4,ks:ke) .AND.&
                     & field(f1:f2,f3:f4,ks:ke)>output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample))&
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
             END IF
          ELSE
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF(l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1 
                            j1 =  j-l_start(2)-hj+1
                            IF ( field(i-is+1+hi,j-js+1+hj,k) > output_fields(out_num)%buffer(i1,j1,k1,sample) ) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Maximum time value 
             ELSE IF ( reduced_k_range ) THEN 
                ksr = l_start(3)
                ker = l_end(3)
                WHERE ( field(f1:f2,f3:f4,ksr:ker) > output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) ) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( field(f1:f2,f3:f4,ks:ke) > output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) ) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
             END IF
          END IF
          output_fields(out_num)%count_0d(sample) = 1
       ELSE IF ( time_min ) THEN
          IF ( PRESENT(mask) ) THEN
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1 
                            j1 =  j-l_start(2)-hj+1
                            IF ( mask(i-is+1+hi,j-js+1+hj,k) .AND.&
                                 & field(i-is+1+hi,j-js+1+hj,k) < output_fields(out_num)%buffer(i1,j1,k1,sample) ) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Minimum time value with masking 
             ELSE IF ( reduced_k_range ) THEN 
                ksr= l_start(3)
                ker= l_end(3)
                WHERE ( mask(f1:f2,f3:f4,ksr:ker) .AND.&
                     & field(f1:f2,f3:f4,ksr:ker) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample)) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( mask(f1:f2,f3:f4,ks:ke) .AND.&
                     & field(f1:f2,f3:f4,ks:ke) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) ) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke) 
             END IF
          ELSE
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                            i1 = i-l_start(1)-hi+1 
                            j1=  j-l_start(2)-hj+1
                            IF ( field(i-is+1+hi,j-js+1+hj,k) < output_fields(out_num)%buffer(i1,j1,k1,sample) ) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Minimum time value 
             ELSE IF ( reduced_k_range ) THEN 
                ksr= l_start(3)
                ker= l_end(3)
                WHERE ( field(f1:f2,f3:f4,ksr:ker) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) ) &
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( field(f1:f2,f3:f4,ks:ke) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) )&
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
             END IF
          END IF
          output_fields(out_num)%count_0d(sample) = 1
       ELSE  ! ( not average, not min, max)
          output_fields(out_num)%count_0d(sample) = 1
          IF ( need_compute ) THEN
             DO j = js, je 
                DO i = is, ie           
                   IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                      i1 = i-l_start(1)-hi+1 
                      j1 = j-l_start(2)-hj+1
                      output_fields(out_num)%buffer(i1,j1,:,sample) = field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))
                   END IF
                END DO
             END DO
             ! instantaneous output
          ELSE IF ( reduced_k_range ) THEN 
             ksr = l_start(3)
             ker = l_end(3)
             output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
          ELSE
             IF ( debug_diag_manager ) THEN
                CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                IF ( err_msg_local /= '' ) THEN
                   IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN 
                      DEALLOCATE(oor_mask)
                      RETURN
                   END IF
                END IF
             END IF
             output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
          END IF
               
          IF ( PRESENT(mask) .AND. missvalue_present ) THEN
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie                                                       
                         IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1 
                            j1 =  j-l_start(2)-hj+1                    
                            IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) )&
                                 & output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue                     
                         END IF
                      END DO
                   END DO
                END DO
             ELSE IF ( reduced_k_range ) THEN 
                ksr= l_start(3)
                ker= l_end(3)
                DO k=ksr, ker 
                   k1= k - ksr + 1
                   DO j=js, je
                      DO i=is, ie
                         IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) ) &
                              & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                      END DO
                   END DO
                END DO
             ELSE
                DO k=ks, ke
                   DO j=js, je
                      DO i=is, ie
                         IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) )&
                              & output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                      END DO
                   END DO
                END DO
             END IF
          END IF
       END IF !average

       IF ( output_fields(out_num)%static .AND. .NOT.need_compute .AND. debug_diag_manager ) THEN
          CALL check_bounds_are_exact_static(out_num, diag_field_id, err_msg=err_msg_local)
          IF ( err_msg_local /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg)) THEN 
                DEALLOCATE(oor_mask)
                RETURN
             END IF
          END IF
       END IF
 
       ! If rmask and missing value present, then insert missing value     
       IF ( PRESENT(rmask) .AND. missvalue_present ) THEN
          IF ( need_compute ) THEN
             DO k = l_start(3), l_end(3)
                k1 = k - l_start(3) + 1
                DO j = js, je 
                   DO i = is, ie                                               
                      IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                         i1 = i-l_start(1)-hi+1 
                         j1 =  j-l_start(2)-hj+1
                         IF ( rmask(i-is+1+hi,j-js+1+hj,k) < 0.5 ) &
                              & output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue                                   
                      END IF
                   END DO
                END DO
             END DO
          ELSE IF ( reduced_k_range ) THEN 
             ksr= l_start(3)
             ker= l_end(3)
             DO k= ksr, ker 
                k1 = k - ksr + 1
                DO j=js, je
                   DO i=is, ie
                      IF ( rmask(i-is+1+hi,j-js+1+hj,k) < 0.5 ) &
                           & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                   END DO
                END DO
             END DO
          ELSE
             DO k=ks, ke
                DO j=js, je
                   DO i=is, ie
                      IF ( rmask(i-is+1+hi,j-js+1+hj,k) < 0.5 ) &
                           & output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                   END DO
                END DO
             END DO
          END IF
       END IF

    END DO num_out_fields
    
    DEALLOCATE(oor_mask)
  END FUNCTION send_data_3d
  ! </FUNCTION>

  ! <FUNCTION NAME="send_tile_averaged_data2d" INTERFACE="send_tile_averaged_data">
  !   <IN NAME="diag_field_id" TYPE="INTEGER"></IN>
  !   <IN NAME="field" TYPE="REAL, DIMENSION(:,:,:)"></IN>
  !   <IN NAME="area" TYPE="REAL, DIMENSION(:,:,:)">  </IN>
  !   <IN NAME="time" TYPE="TYPE(time_type)">  </IN>
  !   <IN NAME="mask" TYPE="LOGICAL, DIMENSION(:,:,:), OPTIONAL"></IN>
  LOGICAL FUNCTION send_tile_averaged_data2d ( id, field, area, time, mask )
    INTEGER, INTENT(in) :: id  ! id od the diagnostic field 
    REAL, INTENT(in) :: field(:,:,:) ! field to average and send
    REAL, INTENT(in) :: area (:,:,:) ! area of tiles (== averaging weights), arbitrary units
    TYPE(time_type), INTENT(in)  :: time ! current time
    LOGICAL, INTENT(in),OPTIONAL :: mask (:,:,:) ! land mask

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2)) :: out(SIZE(field,1), SIZE(field,2))

    CALL average_tiles(id, field, area, mask, out)
    send_tile_averaged_data2d = send_data(id, out, time, mask=ANY(mask,DIM=3))
  END FUNCTION send_tile_averaged_data2d
  ! </FUNCTION>

  ! <FUNCTION NAME="send_tile_averaged_data3d" INTERFACE="send_tile_averaged_data">
  !   <IN NAME="diag_field_id" TYPE="INTEGER"></IN>
  !   <IN NAME="field" TYPE="REAL, DIMENSION(:,:,:,:)"></IN>
  !   <IN NAME="area" TYPE="REAL, DIMENSION(:,:,:)"></IN>
  !   <IN NAME="time" TYPE="TYPE(time_type)"></IN>
  !   <IN NAME="mask" TYPE="LOGICAL, DIMENSION(:,:,:), OPTIONAL">  </IN>
  LOGICAL FUNCTION send_tile_averaged_data3d( id, field, area, time, mask )
    INTEGER, INTENT(in) :: id ! id of the diagnostic field
    REAL, DIMENSION(:,:,:,:), INTENT(in) :: field ! (lon, lat, tile, lev) field to average and send
    REAL, DIMENSION(:,:,:), INTENT(in) :: area (:,:,:) ! (lon, lat, tile) tile areas ( == averaging weights), arbitrary units
    TYPE(time_type), INTENT(in)  :: time ! current time
    LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask ! (lon, lat, tile) land mask

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2),SIZE(field,4)) :: out
    LOGICAL, DIMENSION(SIZE(field,1),SIZE(field,2),SIZE(field,4)) :: mask3
    INTEGER :: it

    DO it=1, SIZE(field,4)
       CALL average_tiles(id, field(:,:,:,it), area, mask, out(:,:,it) )
    END DO

    mask3(:,:,1) = ANY(mask,DIM=3)
    DO it = 2, SIZE(field,4)
       mask3(:,:,it) = mask3(:,:,1)
    END DO

    send_tile_averaged_data3d = send_data( id, out, time, mask=mask3 )
  END FUNCTION send_tile_averaged_data3d
  ! </FUNCTION>

  ! <SUBROUTINE NAME="average_tiles">
  !   <OVERVIEW>
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE average_tiles(diag_field_id, x, area, mask, out)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !   </DESCRIPTION>
  !   <IN NAME="diag_field_id" TYPE="INTEGER"></IN>
  !   <IN NAME="x" TYPE="REAL, DIMENSION(:,:,:)">(lon, lat, tile) field to average</IN>
  !   <IN NAME="area" TYPE="REAL, DIMENSION(:,:,:)">(lon, lat, tile) fractional area</IN>
  !   <IN NAME="mask" TYPE="LOGICAL, DIMENSION(:,:,:)">(lon, lat, tile) land mask</IN>
  !   <OUT NAME="out" TYPE="REAL, DIMENSION(:,:)">(lon, lat) result of averaging</OUT>
  SUBROUTINE average_tiles(diag_field_id, x, area, mask, out)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:,:,:), INTENT(in) :: x 
    REAL, DIMENSION(:,:,:), INTENT(in) :: area
    LOGICAL, DIMENSION(:,:,:), INTENT(in) :: mask
    REAL, DIMENSION(:,:), INTENT(out) :: out 

    INTEGER  :: it ! iterator over tile number
    REAL, DIMENSION(SIZE(x,1),SIZE(x,2)) :: s ! area accumulator
    REAL :: local_missing_value

    ! Initialize local_missing_value
    IF ( input_fields(diag_field_id)%missing_value_present ) THEN
       local_missing_value = input_fields(diag_field_id)%missing_value
    ELSE 
       local_missing_value = 0.0
    END IF
    
    ! Initialize s and out to zero.
    s(:,:) = 0.0
    out(:,:) = 0.0
    
    DO it = 1, SIZE(area,3)
       WHERE ( mask(:,:,it) ) 
          out(:,:) = out(:,:) + x(:,:,it)*area(:,:,it)
          s(:,:) = s(:,:) + area(:,:,it)
       END WHERE
    END DO

    WHERE ( s(:,:) > 0 ) 
       out(:,:) = out(:,:)/s(:,:)
    ELSEWHERE
       out(:,:) = local_missing_value
    END WHERE
  END SUBROUTINE average_tiles
  ! </SUBROUTINE>

  INTEGER FUNCTION writing_field(out_num, at_diag_end, error_string, time)
    INTEGER, INTENT(in) :: out_num
    LOGICAL, INTENT(in) :: at_diag_end
    CHARACTER(len=*), INTENT(out) :: error_string
    TYPE(time_type), INTENT(in) :: time

    TYPE(time_type) :: middle_time
    LOGICAL :: time_max, time_min, reduced_k_range, missvalue_present
    LOGICAL :: average, need_compute, phys_window
    INTEGER :: in_num, file_num, freq, units
    INTEGER :: b1,b2,b3,b4 ! size of buffer along x,y,z,and diurnal axes
    INTEGER :: i, j, k, m
    REAL    :: missvalue, num

    writing_field = 0
    
    need_compute = output_fields(out_num)%need_compute

    in_num = output_fields(out_num)%input_field
    IF ( input_fields(in_num)%static ) RETURN

    missvalue = input_fields(in_num)%missing_value
    missvalue_present = input_fields(in_num)%missing_value_present
    reduced_k_range = output_fields(out_num)%reduced_k_range
    phys_window = output_fields(out_num)%phys_window
    ! Is this output field being time averaged?
    average = output_fields(out_num)%time_average
    ! Looking for max and min value of this field over the sampling interval?
    time_max = output_fields(out_num)%time_max
    time_min = output_fields(out_num)%time_min   
    file_num = output_fields(out_num)%output_file
    freq = files(file_num)%output_freq
    units = files(file_num)%output_units

    ! If average get size: Average intervals are last_output, next_output
    IF ( average ) THEN
       b1=SIZE(output_fields(out_num)%buffer,1)
       b2=SIZE(output_fields(out_num)%buffer,2) 
       b3=SIZE(output_fields(out_num)%buffer,3)
       b4=SIZE(output_fields(out_num)%buffer,4)
       IF ( input_fields(in_num)%mask_variant ) THEN           
          DO m=1, b4 
             DO k=1, b3
                DO j=1, b2
                   DO i=1, b1 
                      IF ( output_fields(out_num)%counter(i,j,k,m) > 0. )THEN
                         output_fields(out_num)%buffer(i,j,k,m) = &
                              & output_fields(out_num)%buffer(i,j,k,m)/output_fields(out_num)%counter(i,j,k,m)
                      ELSE
                         output_fields(out_num)%buffer(i,j,k,m) =  missvalue
                      END IF
                   END DO
                END DO
             END DO
          END DO
       ELSE  !not mask variant
          DO m = 1, b4
             IF ( phys_window ) THEN
                IF ( need_compute .OR. reduced_k_range ) THEN
                   num = REAL(output_fields(out_num)%num_elements(m)/output_fields(out_num)%region_elements)
                ELSE
                   num = REAL(output_fields(out_num)%num_elements(m)/output_fields(out_num)%total_elements)
                END IF
             ELSE
                num = output_fields(out_num)%count_0d(m)
             END IF
             IF ( num > 0. ) THEN
                IF ( missvalue_present ) THEN
                   DO k=1, b3
                      DO j=1, b2
                         DO i=1, b1
                            IF ( output_fields(out_num)%buffer(i,j,k,m) /= missvalue ) &
                                 & output_fields(out_num)%buffer(i,j,k,m) = output_fields(out_num)%buffer(i,j,k,m)/num  
                         END DO
                      END DO
                   END DO
                ELSE
                   output_fields(out_num)%buffer(:,:,:,m) = output_fields(out_num)%buffer(:,:,:,m)/num
                END IF
             ELSE IF ( .NOT. at_diag_end ) THEN
                IF ( missvalue_present ) THEN
                   IF(ANY(output_fields(out_num)%buffer /= missvalue)) THEN
                      WRITE (error_string,'(a,"/",a)')&
                           & TRIM(input_fields(in_num)%module_name), &
                           & TRIM(output_fields(out_num)%output_name)
                      writing_field = -1
                      RETURN
                   END IF
                END IF
             END IF
          END DO
       END IF ! mask_variant
    ELSE IF ( time_min .OR. time_max ) THEN
       IF ( missvalue_present ) THEN
          WHERE ( ABS(output_fields(out_num)%buffer) == MIN_VALUE ) 
             output_fields(out_num)%buffer = missvalue
          END WHERE
       END IF ! if missvalue is NOT present buffer retains max_value or min_value
    END IF !average

    ! Output field
    IF ( at_diag_end .AND. freq == END_OF_RUN ) output_fields(out_num)%next_output = time
    IF ( (output_fields(out_num)%time_ops) .AND. (.NOT. mix_snapshot_average_fields) ) THEN
       middle_time = (output_fields(out_num)%last_output+output_fields(out_num)%next_output)/2
       CALL diag_data_out(file_num, out_num, output_fields(out_num)%buffer, middle_time)
    ELSE
       CALL diag_data_out(file_num, out_num, &
            & output_fields(out_num)%buffer, output_fields(out_num)%next_output)
    END IF

    IF ( at_diag_end ) RETURN

    ! Take care of cleaning up the time counters and the storeage size
    output_fields(out_num)%last_output = output_fields(out_num)%next_output
    IF ( freq == END_OF_RUN ) THEN
       output_fields(out_num)%next_output = time
    ELSE
       IF ( freq == EVERY_TIME ) THEN
          output_fields(out_num)%next_output = time
       ELSE
          output_fields(out_num)%next_output = output_fields(out_num)%next_next_output
          output_fields(out_num)%next_next_output = &
               & diag_time_inc(output_fields(out_num)%next_next_output, freq, units)
       END IF
       output_fields(out_num)%count_0d(:) = 0.0
       output_fields(out_num)%num_elements(:) = 0
       IF ( time_max ) THEN 
          output_fields(out_num)%buffer = MAX_VALUE
       ELSE IF ( time_min ) THEN
          output_fields(out_num)%buffer = MIN_VALUE
       ELSE
          output_fields(out_num)%buffer = EMPTY
       END IF
       IF ( input_fields(in_num)%mask_variant .AND. average ) output_fields(out_num)%counter = 0.0
    END IF

  END FUNCTION writing_field

  SUBROUTINE diag_manager_set_time_end(Time_end_in)
    TYPE (time_type), INTENT(in) :: Time_end_in

    Time_end = Time_end_in

  END SUBROUTINE diag_manager_set_time_end

  !-----------------------------------------------------------------------
  SUBROUTINE diag_send_complete(time_step, err_msg)
    TYPE (time_type), INTENT(in)           :: time_step
    character(len=*), INTENT(out), optional :: err_msg
     
    type(time_type)    :: next_time, time
    integer            :: file, j, out_num, in_num, freq, status
    logical            :: local_output, need_compute
    CHARACTER(len=128) :: error_string

    IF ( Time_end == Time_zero ) THEN
       ! <ERROR STATUS="FATAL">
       !   diag_manager_set_time_end must be called before diag_send_complete
       ! </ERROR>
       CALL error_mesg('diag_manager_mod::diag_send_complete',&
            & "diag_manager_set_time_end must be called before diag_send_complete", FATAL)
    END IF

    DO file = 1, num_files
       freq = files(file)%output_freq
       DO j = 1, files(file)%num_fields
          out_num = files(file)%fields(j) !this is position of output_field in array output_fields
          in_num = output_fields(out_num)%input_field

          IF ( input_fields(in_num)%numthreads == 1 ) CYCLE
          IF ( output_fields(out_num)%static .OR. freq == END_OF_RUN ) CYCLE
          time = input_fields(in_num)%time
          IF ( time >= time_end ) CYCLE

          ! is this field output on a local domain only?
          local_output = output_fields(out_num)%local_output
          ! if local_output, does the current PE take part in send_data?
          need_compute = output_fields(out_num)%need_compute
          ! skip all PEs not participating in outputting this field
          IF ( local_output .AND. (.NOT.need_compute) ) CYCLE
          next_time = time + time_step

          IF ( next_time > output_fields(out_num)%next_output ) THEN
             ! A non-static field that has skipped a time level is an error
             IF ( next_time > output_fields(out_num)%next_next_output .AND. freq > 0 ) THEN
                IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN 
                   WRITE (error_string,'(a,"/",a)')&
                        & TRIM(input_fields(in_num)%module_name), &
                        & TRIM(output_fields(out_num)%output_name)
                   IF ( fms_error_handler('diag_manager_mod::diag_send_complete', 'module/output_field '//TRIM(error_string)//&
                        & ' is skipped one time level in output data', err_msg)) RETURN
                END IF
             END IF

             status = writing_field(out_num, .FALSE., error_string, next_time)
             IF ( status == -1 ) THEN
                IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
                   IF(fms_error_handler('diag_manager_mod::diag_send_complete','module/output_field '//TRIM(error_string)//&
                        & ', write EMPTY buffer', err_msg)) RETURN
                END IF
             END IF
          END IF  !time > output_fields(out_num)%next_output
       END DO
    END DO

  END SUBROUTINE diag_send_complete

  ! <SUBROUTINE NAME="diag_manager_end">
  !   <OVERVIEW>
  !     Exit Diagnostics Manager.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Flushes diagnostic buffers where necessary. Close diagnostics files.
  !
  !     A warning will be issued here if a field in diag_table is not registered
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     SUBROUTINE diag_manager_end(time)
  !   </TEMPLATE>
  !   <IN NAME="TIME" TYPE="time_type"></IN>
  SUBROUTINE diag_manager_end(time)
    TYPE(time_type), INTENT(in) :: time

    INTEGER :: file

    IF ( do_diag_field_log ) THEN
       CALL mpp_close (diag_log_unit)
    END IF
    DO file = 1, num_files
       CALL closing_file(file, time)   
    END DO
  END SUBROUTINE diag_manager_end
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="closing_file">
  !   <OVERVIEW>
  !     Replaces diag_manager_end; close just one file: files(file)
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE closing_file(file, time)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !   </DESCRIPTION>
  !   <IN NAME="file" TYPE="INTEGER"></IN>
  !   <IN NAME="tile" TYPE="TYPE(time_type)"></IN>
  SUBROUTINE closing_file(file, time)
    INTEGER, INTENT(in) :: file
    TYPE(time_type), INTENT(in) :: time

    INTEGER :: j, i, input_num, freq, status
    INTEGER :: stdout_unit
    LOGICAL :: reduced_k_range, need_compute, local_output
    CHARACTER(len=128) :: message

    stdout_unit = stdout()

    ! Output all registered, non_static output_fields
    DO j = 1, files(file)%num_fields
       i = files(file)%fields(j) !this is position of output_field in array output_fields

       ! is this field output on a local domain only?
       local_output = output_fields(i)%local_output
       ! if local_output, does the current PE take part in send_data?
       need_compute = output_fields(i)%need_compute

       reduced_k_range = output_fields(i)%reduced_k_range

       ! skip all PEs not participating in outputting this field
       IF ( local_output .AND. (.NOT. need_compute) ) CYCLE
       ! skip fields that were not registered or non-static   
       input_num = output_fields(i)%input_field
       IF ( input_fields(input_num)%static ) CYCLE
       IF ( .NOT.input_fields(input_num)%register ) CYCLE 
       freq = files(file)%output_freq
       IF ( freq /= END_OF_RUN .AND. files(file)%file_unit < 0 &
            & .AND. ALL(output_fields(i)%num_elements(:) == 0) .AND. ALL(output_fields(i)%count_0d(:) == 0) ) CYCLE
       ! Is it time to output for this field; CAREFUL ABOUT >= vs > HERE
       ! For end should be >= because no more data is coming 
       IF ( time >= output_fields(i)%next_output .OR. freq == END_OF_RUN ) THEN
          IF ( time >= output_fields(i)%next_next_output .AND. freq > 0 ) THEN
             WRITE (message,'(a,"/",a)') TRIM(input_fields(input_num)%module_name), &
                  & TRIM(output_fields(i)%output_name)
             ! <ERROR STATUS="WARNING">
             !   <input_fields(input_num)%module_name>/<output_fields(i)%output_name> skip one time
             !   level, maybe send_data never called
             ! </ERROR>
             IF ( mpp_pe() .EQ. mpp_root_pe() ) & 
                  & CALL error_mesg('diag_manager_mod::closing_file', 'module/output_field ' //&
                  & TRIM(message)//', skip one time level, maybe send_data never called', WARNING)
          END IF
          
          status = writing_field(i, .TRUE., message, time)

       ELSEIF ( .NOT.output_fields(i)%written_once ) THEN
          ! <ERROR STATUS="NOTE">
          !   <output_fields(i)%output_name) NOT available, check if output interval > runlength.
          !   NetCDF fill_values are written
          ! </ERROR>
          CALL error_mesg('Potential error in diag_manager_end ',TRIM(output_fields(i)%output_name)//' NOT available,'//&
               & ' check if output interval > runlength. Netcdf fill_values are written', NOTE)
          output_fields(i)%buffer = FILL_VALUE
          CALL diag_data_out(file, i, output_fields(i)%buffer, time, .TRUE.)   
       END IF
    END DO
    ! Now it's time to output static fields
    CALL write_static(file)

    ! Write out the number of bytes of data saved to this file
    IF ( write_bytes_in_file ) THEN
       CALL mpp_sum (files(file)%bytes_written)
       IF ( mpp_pe() == mpp_root_pe() ) WRITE (stdout_unit,'(a,i12,a,a)') 'Diag_Manager: ',files(file)%bytes_written, &
            & ' bytes of data written to file ',TRIM(files(file)%name)
    END IF
  END SUBROUTINE closing_file
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="diag_manager_init">
  !   <OVERVIEW>
  !     Initialize Diagnostics Manager.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE diag_manager_init(diag_model_subset, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Open and read diag_table. Select fields and files for diagnostic output.
  !   </DESCRIPTION>
  !   <IN NAME="diag_model_subset" TYPE="INTEGER, OPTIONAL"></IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL"></OUT>
  SUBROUTINE diag_manager_init(diag_model_subset, err_msg)
    INTEGER, OPTIONAL, INTENT(IN) :: diag_model_subset
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    CHARACTER(len=*), PARAMETER :: SEP = '|'

    INTEGER, PARAMETER :: FltKind = FLOAT_KIND
    INTEGER, PARAMETER :: DblKind = DOUBLE_KIND
    INTEGER :: diag_subset_output
    INTEGER :: mystat
    INTEGER, ALLOCATABLE, DIMENSION(:) :: pelist
    INTEGER :: stdlog_unit, stdout_unit
#ifndef INTERNAL_FILE_NML
    INTEGER :: nml_unit
#endif
    CHARACTER(len=256) :: err_msg_local

    NAMELIST /diag_manager_nml/ append_pelist_name, mix_snapshot_average_fields, max_output_fields, &
         & max_input_fields, max_axes, do_diag_field_log, write_bytes_in_file, debug_diag_manager,&
         & max_num_axis_sets, max_files, use_cmor, issue_oor_warnings,&
         & oor_warnings_fatal, conserve_water

    ! If the module was already initialized do nothing
    IF ( module_is_initialized ) RETURN

    ! Clear the err_msg variable if contains any residual information
    IF ( PRESENT(err_msg) ) err_msg = ''

    ! Determine pack_size from how many bytes a real value has (how compiled)
    pack_size = SIZE(TRANSFER(0.0_DblKind, (/0.0, 0.0, 0.0, 0.0/)))
    IF ( pack_size.NE.1 .AND. pack_size.NE.2 ) THEN
       IF ( fms_error_handler('diag_manager_mod::diag_manager_init', 'unknown pack_size.  Must be 1, or 2.', err_msg) ) RETURN
    END IF

    ! Get min and max values for real(kind=FLOAT_KIND)
    min_value = HUGE(0.0_FltKind)
    max_value = -min_value

    ! get stdlog and stdout unit number
    stdlog_unit = stdlog()
    stdout_unit = stdout()

    ! version number to logfile
    CALL write_version_number(version, tagname)

    Time_zero = set_time(0,0)
    !--- initialize time_end to time_zero
    Time_end  = Time_zero
    diag_subset_output = DIAG_ALL
    IF ( PRESENT(diag_model_subset) ) THEN
       IF ( diag_model_subset >= DIAG_OTHER .AND. diag_model_subset <= DIAG_ALL ) THEN
          diag_subset_output = diag_model_subset
       ELSE
          IF ( fms_error_handler('diag_manager_mod::diag_manager_init', 'invalid value of diag_model_subset',err_msg) ) RETURN
       END IF
    END IF

#ifdef INTERNAL_FILE_NML
    READ (input_nml_file, NML=diag_manager_nml, IOSTAT=mystat)
#else
    IF ( file_exist('input.nml') ) THEN
       nml_unit = open_namelist_file()
       READ (nml_unit, diag_manager_nml, iostat=mystat)
       CALL close_file(nml_unit)
    ELSE
       ! Set mystat to an arbitrary positive number if input.nml does not exist.
       mystat = 100
    END IF
#endif
    ! Check the status of reading the diag_manager_nml
    
    IF ( check_nml_error(IOSTAT=mystat, NML_NAME='DIAG_MANAGER_NML') < 0 ) THEN
       IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_manager_mod::diag_manager_init', 'DIAG_MANAGER_NML not found in input.nml.  Using defaults.',&
               & WARNING)
       END IF
    END IF

    IF ( mpp_pe() == mpp_root_pe() ) THEN 
       WRITE (stdlog_unit, diag_manager_nml)
    END IF

    ! Issue note about using the CMOR missing value.
    IF ( use_cmor ) THEN 
       err_msg_local = ''
       WRITE (err_msg_local,'(ES8.1E2)') CMOR_MISSING_VALUE
       CALL error_mesg('diag_manager_mod::diag_manager_init', 'Using CMOR missing value ('//TRIM(err_msg_local)//').', NOTE)
    END IF

    ! How to handle Out of Range Warnings.
    IF ( oor_warnings_fatal ) THEN
       oor_warning = FATAL
       CALL error_mesg('diag_manager_mod::diag_manager_init', 'Out &
            &of Range warnings are fatal.', NOTE)
    ELSEIF ( .NOT.issue_oor_warnings ) THEN
       CALL error_mesg('diag_manager_mod::diag_manager_init', 'Out &
            &of Range warnings will be ignored.', NOTE)
    END IF

    IF ( mix_snapshot_average_fields ) THEN
       IF ( mpp_pe() == mpp_root_pe() ) THEN 
          CALL error_mesg('diag_manager_mod::diag_manager_init', 'Setting diag_manager_nml variable '//&
               & 'mix_snapshot_average_fields = .TRUE. will cause ERRORS in the time coordinates '//&
               & 'of all time averaged fields.  Strongly recommend setting mix_snapshot_average_fields '//&
               & '= .FALSE.', WARNING)
       END IF
    END IF
    ALLOCATE(output_fields(max_output_fields))
    ALLOCATE(input_fields(max_input_fields))
    ALLOCATE(files(max_files))
    ALLOCATE(pelist(mpp_npes()))
    CALL mpp_get_current_pelist(pelist, pelist_name)

    CALL parse_diag_table(DIAG_SUBSET=diag_subset_output, ISTAT=mystat, ERR_MSG=err_msg_local)
    IF ( mystat /= 0 ) THEN
       IF ( fms_error_handler('diag_manager_mod::diag_manager_init',&
            & 'Error parsing diag_table. '//TRIM(err_msg_local), err_msg) ) RETURN
    END IF
    
    !initialize files%bytes_written to zero
    files(:)%bytes_written = 0

    ! open diag field log file
    IF ( do_diag_field_log ) THEN
       CALL mpp_open(diag_log_unit, 'diag_field_log.out', nohdrs=.TRUE.)
       WRITE (diag_log_unit,'(777a)') &
            & 'Module',        SEP, 'Field',          SEP, 'Long Name',    SEP,&
            & 'Units',         SEP, 'Number of Axis', SEP, 'Time Axis',    SEP,&
            & 'Missing Value', SEP, 'Min Value',      SEP, 'Max Value',    SEP,&
            & 'AXES LIST'
    END IF

    module_is_initialized = .TRUE.
    ! create axis_id for scalars here
    null_axis_id = diag_axis_init('scalar_axis', (/0./), 'none', 'X', 'none')
    RETURN
  END SUBROUTINE diag_manager_init
  ! </SUBROUTINE>


  ! <FUNCTION NAME="get_base_time">
  !   <OVERVIEW>
  !     Return base time for diagnostics. 
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     TYPE(time_type) FUNCTION get_base_time()
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return base time for diagnostics (note: base time must be >= model time).
  !   </DESCRIPTION>
  TYPE(time_type) FUNCTION get_base_time ()
    ! <ERROR STATUS="FATAL">
    !   MODULE has not been initialized
    ! </ERROR>
    IF ( .NOT.module_is_initialized ) CALL error_mesg('diag_manager_mod::get_base_time', &
         & 'module has not been initialized', FATAL)
    get_base_time = base_time
  END FUNCTION get_base_time
  ! </FUNCTION>

  ! <SUBROUTINE NAME="get_base_date">
  !   <OVERVIEW>
  !     Return base date for diagnostics.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_base_date(year, month, day, hour, minute, second)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return date information for diagnostic reference time.
  !   </DESCRIPTION>
  !   <OUT NAME="year" TYPE="INTEGER"></OUT>
  !   <OUT NAME="month" TYPE="INTEGER"></OUT>
  !   <OUT NAME="day" TYPE="INTEGER"></OUT>
  !   <OUT NAME="hour" TYPE="INTEGER"></OUT>
  !   <OUT NAME="minute" TYPE="INTEGER"></OUT>
  !   <OUT NAME="second" TYPE="INTEGER"></OUT>
  SUBROUTINE get_base_date(year, month, day, hour, minute, second)
    INTEGER, INTENT(out) :: year, month, day, hour, minute, second

    ! <ERROR STATUS="FATAL">module has not been initialized</ERROR>
    IF (.NOT.module_is_initialized) CALL error_mesg ('diag_manager_mod::get_base_date', &
         & 'module has not been initialized', FATAL)
    year   = base_year
    month  = base_month
    day    = base_day
    hour   = base_hour
    minute = base_minute
    second = base_second
  END SUBROUTINE get_base_date
  ! </SUBROUTINE>

  ! <FUNCTION NAME="need_data">
  !   <OVERVIEW>
  !     Determine whether data is needed for the current model time step.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     LOGICAL need_data(diag_field_id, next_model_time)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Determine whether data is needed for the current model time step.
  !     Since diagnostic data are buffered, the "next" model time is passed
  !     instead of the current model time. This call can be used to minimize
  !     overhead for complicated diagnostics.
  !   </DESCRIPTION>
  !   <IN NAME="next_model_time" TYPE="TYPE(time_type)">
  !     next_model_time = current model time + model time_step
  !   </IN>
  !   <IN NAME="diag_field_id" TYPE="INTEGER"></IN>
  LOGICAL FUNCTION need_data(diag_field_id, next_model_time)
    TYPE(time_type), INTENT(in) :: next_model_time
    INTEGER, INTENT(in) :: diag_field_id

    INTEGER :: i, out_num 

    need_data = .FALSE.
    IF ( diag_field_id < 0 ) RETURN ! this field is unused
    DO i = 1, input_fields(diag_field_id)%num_output_fields
       ! Get index to an output field
       out_num = input_fields(diag_field_id)%output_fields(i)
       IF ( .NOT.output_fields(out_num)%static ) THEN
          IF ( next_model_time > output_fields(out_num)%next_output ) need_data=.TRUE.
          ! Is this output field being time averaged?
          ! assume average data based on every timestep
          ! needs to be changed when different forms of averaging are implemented 
          IF ( output_fields(out_num)%time_average) need_data = .TRUE. 
       END IF
    END DO
    RETURN
  END FUNCTION need_data
  ! </FUNCTION>

  ! <SUBROUTINE NAME="set_diag_filename_appendix">
  !   <OVERVIEW>
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE set_diag_filename_appendix(string_in)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !   </DESCRIPTION>
  !   <IN NAME="string_in" TYPE="CHARACTER(len=*)"></IN>
  SUBROUTINE set_diag_filename_appendix(string_in)
    CHARACTER(len=*) , INTENT(in) :: string_in
    
    filename_appendix = TRIM(string_in)
  END SUBROUTINE set_diag_filename_appendix
  ! </SUBROUTINE>

  ! <FUNCTION NAME="init_diurnal_axis">
  !   <OVERVIEW>
  !     Finds or initializes a diurnal time axis and returns its' ID.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION init_diurnal_axis(n_samples)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Given number of time intervals in the day, finds or initializes a diurnal time axis
  !     and returns its ID. It uses get_base_date, so should be in the file where it's accessible.
  !     The units are 'days since BASE_DATE', all diurnal axes belong to the set 'diurnal'
  !   </DESCRIPTION>
  !   <IN NAME="n_samples" TYPE="INTEGER">Number of intervals during the day</IN>
  INTEGER FUNCTION init_diurnal_axis(n_samples)
    INTEGER, INTENT(in) :: n_samples ! number of intervals during the day

    REAL :: DATA  (n_samples)   ! central points of time intervals
    REAL :: edges (n_samples+1) ! boundaries of time intervals
    INTEGER :: edges_id ! id of the corresponding edges
    INTEGER :: i
    INTEGER :: year, month, day, hour, minute, second ! components of the base date
    CHARACTER(32)  :: name  ! name of the axis
    CHARACTER(128) :: units ! units of time

    CALL get_base_date(year, month, day, hour, minute, second)
    WRITE (units,11) 'hours', year, month, day, hour, minute, second
11  FORMAT(a,' since ',i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2)
    ! compute central points and units
    edges(1) = 0.0
    DO i = 1, n_samples
       DATA (i) = 24.0*(REAL(i)-0.5)/n_samples
       edges(i+1) = 24.0* REAL(i)/n_samples
    END DO

    ! define edges
    name = ''
    WRITE (name,'(a,i2.2)') 'time_of_day_edges_', n_samples
    edges_id = get_axis_num(name, 'diurnal')
    IF ( edges_id <= 0 ) THEN
       edges_id =  diag_axis_init(name,edges,units,'N','time of day edges', set_name='diurnal')
    END IF
  
    ! define axis itself
    name = ''
    WRITE (name,'(a,i2.2)') 'time_of_day_', n_samples
    init_diurnal_axis = get_axis_num(name, 'diurnal')
    IF ( init_diurnal_axis <= 0 ) THEN
       init_diurnal_axis = diag_axis_init(name, DATA, units, 'N', 'time of day', set_name='diurnal', edges=edges_id)
    END IF
  END FUNCTION init_diurnal_axis
  ! </FUNCTION>
END MODULE diag_manager_mod

! <INFO>
!   <COMPILER NAME="PORTABILITY">
!     <TT>diag_manager_mod</TT> uses standard Fortran 90.
!   </COMPILER>
!   <COMPILER NAME="ACQUIRING SOURCE">
!     Use the following commands to check out the source at GFDL.
!     <PRE>
!       setenv CVSROOT '/home/fms/cvs'
!       cvs co diag_manager
!     </PRE>
!   </COMPILER>
!   <COMPILER NAME="COMPILING AND LINKING SOURCE">
!     Any module or program unit using <TT>diag_manager_mod</TT> must contain the line
!     <PRE>
!     use diag_manager_mod
!     </PRE>
!     If netCDF output is desired, the cpp flag <TT>-Duse_netCDF</TT>
!     must be turned on. 
!   </COMPILER>
!   <PRECOMP FLAG="-Duse_netCDF"> 
!     Used to write out <LINK SRC="http://www.unidata.ucar.edu/software/netcdf">NetCDF</LINK> files.
!   </PRECOMP>
!   <PRECOMP FLAG="-Dtest_diag_manager">
!     Used to build the unit test suite for the <TT>diag_manager_mod</TT>.
!   </PRECOMP>
!   <LOADER FLAG="-lnetcdf">
!     Link in the NetCDF libraries.
!   </LOADER>
!   <TESTPROGRAM NAME="test">
!     Unit test for the <TT>diag_manager_mod</TT>.  Each test must be run separately, and ends with an intentional fatal error.
!     Each test has its own <TT>diag_table</TT>, see the source of <TT>diag_manager.F90</TT> for the list of <TT>diag_tables</TT>
!     for the unit tests.
!   </TESTPROGRAM>
!   <FUTURE>
!     Regional output for the cubed-sphere grid.
!   </FUTURE>
! </INFO>

! ********** Test Program **********
#ifdef test_diag_manager
! This program runs only one of many possible tests with each execution.
! Each test ends with an intentional fatal error.
! diag_manager_mod is not a stateless module, and there are situations
! where a fatal error leaves the module in a state that does not allow
! it to function properly if used again. Therefore, the program must
! be terminated after each intentional fatal error.

! Each test is dependent on the diag_table, and different diag_tables
! exist for each test. Depending on the test, an intentional fatal error
! may be triggered upon the call to diag_manager_init, register_diag_field or send_data.
! Because of this, the calls to all of those routines differ depending on the test.

! The diag_table for each test is included below.

!--------------------------------------------------------------------------------------------------
! diag_table for test 1

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 2

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 3

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 4

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 5

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 6

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 7

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 8

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 9

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "bk",   "bk",   "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 10

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "bk",   "bk",   "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 11

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 12

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
! # Test of the error check that duplicate field names do not appear in same file,
!  "test_mod",              "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 13

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "months", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
! # Test of WARNING message that no data is written when run length is less than output interval  
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 14

! test_diag_manager
! 1990 1 29 0 0 0
! #output files
!  "diag_test2", 1, "months", 1, "days", "time",
! #output variables
! # Test of check for invalid date. (Jan 29 1990 + one month = Feb 29 1990)
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------

PROGRAM test
  ! This program runs only one of many possible tests with each execution.
  ! Each test ends with an intentional fatal error.
  ! diag_manager_mod is not a stateless module, and there are situations
  ! where a fatal error leaves the module in a state that does not allow
  ! it to function properly if used again. Therefore, the program must
  ! be terminated after each intentional fatal error.

  ! Each test is dependent on the diag_table, and different diag_tables
  ! exist for each test. Depending on the test, an intentional fatal error
  ! may be triggered upon the call to diag_manager_init, register_diag_field or send_data.
  ! Because of this, the calls to all of those routines differ depending on the test.

  USE mpp_mod, ONLY: mpp_pe, mpp_error, FATAL
  USE mpp_domains_mod, ONLY: domain2d, mpp_define_domains, mpp_get_compute_domain
  USE mpp_domains_mod, ONLY: mpp_define_io_domain, mpp_define_layout
  USE fms_mod, ONLY: fms_init, fms_end, mpp_npes, file_exist, check_nml_error, close_file, open_file
  USE fms_mod, ONLY: error_mesg, FATAL, stdlog
#ifdef INTERNAL_FILE_NML
  USE mpp_mod, ONLY: input_nml_file
#else
  USE fms_mod, ONLY:  open_namelist_file
#endif
  USE fms_io_mod, ONLY: fms_io_exit
  USE constants_mod, ONLY: constants_init, PI, RAD_TO_DEG

  USE time_manager_mod, ONLY: time_type, set_calendar_type, set_date, decrement_date, OPERATOR(+), set_time
  USE time_manager_mod, ONLY: NOLEAP, JULIAN, GREGORIAN, THIRTY_DAY_MONTHS, OPERATOR(*), assignment(=)
  use time_manager_mod, only: operator(+), operator(-), operator(/), days_in_month

  USE diag_manager_mod, ONLY: diag_manager_init, send_data, diag_axis_init, diag_manager_end
  USE diag_manager_mod, ONLY: register_static_field, register_diag_field, diag_send_complete
  USE diag_manager_mod, ONLY: diag_manager_set_time_end

  IMPLICIT NONE

  TYPE(domain2d) :: Domain1
  TYPE(domain2d) :: Domain2

  REAL, ALLOCATABLE, DIMENSION(:) :: lon_global1, lonb_global1
  REAL, ALLOCATABLE, DIMENSION(:) :: lat_global1, latb_global1
  REAL, ALLOCATABLE, DIMENSION(:) :: lon_global2, lonb_global2
  REAL, ALLOCATABLE, DIMENSION(:) :: lat_global2, latb_global2
  REAL, ALLOCATABLE, DIMENSION(:) :: pfull, bk, phalf
  REAL, ALLOCATABLE, DIMENSION(:) :: lon1, lat1, lonb1, latb1
  REAL, ALLOCATABLE, DIMENSION(:) :: lon2, lat2, lonb2, latb2
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dat1, dat1h
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dat2, dat2h
  REAL, ALLOCATABLE, DIMENSION(:,:) :: dat2_2d
  real    :: solar_constant = 1600
  REAL :: dp, surf_press=1.e5
  INTEGER :: id_phalf, id_pfull, id_bk
  INTEGER :: id_lon1, id_lonb1, id_latb1, id_lat1, id_dat1
  INTEGER :: id_lon2, id_lat2, id_dat2, id_dat2_2d, id_sol_con
  INTEGER :: i, j, k, is1, ie1, js1, je1, nml_unit, io, ierr, log_unit, out_unit, m
  INTEGER :: is_in, ie_in, js_in, je_in
  INTEGER :: is2, ie2, js2, je2, hi=1, hj=1
  INTEGER :: nlon1, nlat1, nlon2, nlat2
  INTEGER, DIMENSION(2) :: layout = (/0,0/)
  INTEGER :: test_number=1
  INTEGER :: nlon=18, nlat=18, nlev=2
  INTEGER :: io_layout(2) = (/0,0/) 
  INTEGER :: nstep = 2
  TYPE(time_type) :: Time, Time_step, Time_end, Time_start, Run_length
  LOGICAL :: used, test_successful
  CHARACTER(len=256) :: err_msg
  integer :: omp_get_num_threads

  integer :: nyc1, n, jsw, jew, isw, iew
  integer :: numthreads=1, ny_per_thread, idthread
  integer :: months=0, days=0, dt_step=0


  NAMELIST /test_diag_manager_nml/ layout, test_number, nlon, nlat, nlev, io_layout, numthreads, &
                                   dt_step, months, days

  CALL fms_init
  nml_unit = open_namelist_file()
  log_unit = stdlog()
  out_unit = open_file(file='test_diag_manager.out', form='formatted', threading='multi', action='write')
  CALL constants_init
  CALL set_calendar_type(JULIAN)

#ifdef INTERNAL_FILE_NML
  READ (input_nml_file, NML=test_diag_manager_nml, IOSTAT=io)
  ierr = check_nml_error(io, 'test_diag_manager_nml')
#else
  IF ( file_exist('input.nml') ) THEN
     ierr=1
     DO WHILE (ierr > 0)
        READ(nml_unit, nml=test_diag_manager_nml, iostat=io)
        ierr = check_nml_error(io, 'test_diag_manager_nml')
     END DO
     CALL close_file(nml_unit)
  END IF
#endif
  WRITE (log_unit,test_diag_manager_nml)

  IF ( test_number == 12 ) THEN
     CALL diag_manager_init(err_msg=err_msg)
     IF ( err_msg /= '' ) THEN
        WRITE (out_unit,'(a)') 'test12 successful: err_msg='//TRIM(err_msg)
        CALL error_mesg('test_diag_manager','test12 successful.',FATAL)
     ELSE
        WRITE (out_unit,'(a)') 'test12 fails'
        CALL error_mesg('test_diag_manager','test12 fails',FATAL)
     END IF
  ELSE
     CALL diag_manager_init
  END IF

  IF ( layout(1)*layout(2) .NE. mpp_npes() ) THEN
     CALL mpp_define_layout((/1,nlon,1,nlat/), mpp_npes(), layout )
  END IF

  nlon1 = nlon
  nlat1 = nlat
  nlon2 = nlon * 2
  nlat2 = nlat * 2

  CALL mpp_define_domains((/1,nlon1,1,nlat1/), layout, Domain1, name='test_diag_manager')
  CALL mpp_get_compute_domain(Domain1, is1, ie1, js1, je1)
  ALLOCATE(lon_global1(nlon1), lonb_global1(nlon1+1))
  ALLOCATE(lat_global1(nlat1), latb_global1(nlat1+1))
  ALLOCATE(lon_global2(nlon2), lonb_global2(nlon2+1))
  ALLOCATE(lat_global2(nlat2), latb_global2(nlat2+1))
  ALLOCATE(pfull(nlev), bk(nlev), phalf(nlev+1))

  ALLOCATE(lon1(is1:ie1), lat1(js1:je1), lonb1(is1:ie1+1), latb1(js1:je1+1))
  CALL compute_grid(nlon1, nlat1, is1, ie1, js1, je1, lon_global1, lat_global1, lonb_global1, latb_global1, lon1, lat1, lonb1, latb1)
  CALL mpp_define_domains((/1,nlon2,1,nlat2/), layout, Domain2, name='test_diag_manager')
  CALL mpp_get_compute_domain(Domain2, is2, ie2, js2, je2)
  CALL mpp_define_io_domain(Domain1, io_layout)
  CALL mpp_define_io_domain(Domain2, io_layout)    

  ALLOCATE(lon2(is2:ie2), lat2(js2:je2), lonb2(is2:ie2+1), latb2(js2:je2+1))
  CALL compute_grid(nlon2, nlat2, is2, ie2, js2, je2, lon_global2, lat_global2, lonb_global2, latb_global2, lon2, lat2, lonb2, latb2)
  dp = surf_press/nlev
  DO k=1, nlev+1
     phalf(k) = dp*(k-1)
  END DO
  DO k=1, nlev
     pfull(k) = .5*(phalf(k) + phalf(k+1))
     bk(k) = pfull(k)/surf_press
  END DO

  ALLOCATE(dat1(is1:ie1,js1:je1,nlev))
  ALLOCATE(dat1h(is1-hi:ie1+hi,js1-hj:je1+hj,nlev))
  dat1h = 0.
  DO j=js1, je1
     DO i=is1, ie1
        dat1(i,j,1) = SIN(lon1(i))*COS(lat1(j))
     END DO
  END DO
  dat1h(is1:ie1,js1:je1,1) = dat1(:,:,1)
  dat1(:,:,2) = -dat1(:,:,1)
  dat1h(:,:,2) = -dat1h(:,:,1)

  ALLOCATE(dat2(is2:ie2,js2:je2,nlev))
  ALLOCATE(dat2_2d(is2:ie2,js2:je2))
  ALLOCATE(dat2h(is2-hi:ie2+hi,js2-hj:je2+hj,nlev))
  dat2h = 0.
  DO j=js2, je2
     DO i=is2, ie2
        dat2(i,j,1) = SIN(lon2(i))*COS(lat2(j))
     END DO
  END DO
  dat2h(is2:ie2,js2:je2,1) = dat2(:,:,1)
  dat2(:,:,2) = -dat2(:,:,1)
  dat2h(:,:,2) = -dat2h(:,:,1)
  dat2_2d = dat2(:,:,1)

  id_lonb1 = diag_axis_init('lonb1', RAD_TO_DEG*lonb_global1, 'degrees_E', 'x', long_name='longitude edges', Domain2=Domain1)
  id_latb1 = diag_axis_init('latb1', RAD_TO_DEG*latb_global1, 'degrees_N', 'y', long_name='latitude edges',  Domain2=Domain1)

  id_lon1  = diag_axis_init('lon1',  RAD_TO_DEG*lon_global1, 'degrees_E','x',long_name='longitude',Domain2=Domain1,edges=id_lonb1)
  id_lat1  = diag_axis_init('lat1',  RAD_TO_DEG*lat_global1, 'degrees_N','y',long_name='latitude', Domain2=Domain1,edges=id_latb1)

  id_phalf= diag_axis_init('phalf', phalf, 'Pa', 'z', long_name='half pressure level', direction=-1)
  id_pfull= diag_axis_init('pfull', pfull, 'Pa', 'z', long_name='full pressure level', direction=-1, edges=id_phalf)

  id_lon2 = diag_axis_init('lon2',  RAD_TO_DEG*lon_global2,  'degrees_E', 'x', long_name='longitude', Domain2=Domain2)
  id_lat2 = diag_axis_init('lat2',  RAD_TO_DEG*lat_global2,  'degrees_N', 'y', long_name='latitude',  Domain2=Domain2)

  IF ( test_number == 14 ) THEN
     Time = set_date(1990,1,29,0,0,0)
  ELSE
     Time = set_date(1990,1,1,0,0,0)
  END IF

  id_dat1 = register_diag_field('test_diag_manager_mod', 'dat1', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K')
  id_dat2 = register_diag_field('test_diag_manager_mod', 'dat2', (/id_lon2,id_lat2,id_pfull/), Time, 'sample data', 'K')

  id_sol_con = register_diag_field ('test_diag_manager_mod', 'solar_constant', Time, &
                  'solar constant', 'watts/m2') 

  if( dt_step == 0 ) call error_mesg ('test_diag_manager',&
       & 'dt_step is not set', FATAL)

  Time_step = set_time(dt_step,0)
  Time_start = Time
  Time_end = Time
  do m = 1,months
       Time_end = Time_end + set_time(0,days_in_month(Time_end))
    end do
    Time_end   = Time_end + set_time(0, days)
  Run_length = Time_end - Time_start
  nstep = Run_length / Time_step

  !-- The following is used to test openMP
  IF ( test_number == 15 ) THEN
!$      call omp_set_num_threads(numthreads)
      nyc1 = je1 - js1 + 1
      IF (MOD(nyc1, numthreads ) /= 0) THEN
         CALL error_mesg ('test_diag_manager',&
              & 'The number of OpenMP threads must be an integral multiple &
              &of the number of rows in the compute domain', FATAL)
     END IF
     ny_per_thread = nyc1/numthreads

     dat1 = 1
     CALL diag_manager_set_time_end(Time_end)
     DO n = 1, nstep

        Time = Time + Time_step
        !$OMP parallel do default(shared) private(isw, iew, jsw, jew )
        
        DO jsw = js1, je1, ny_per_thread
           jew = jsw + ny_per_thread -1
           isw = is1 
           iew = ie1
           if(id_dat1>0) used = send_data(id_dat1, dat1(isw:iew, jsw:jew,:), Time, &
                                is_in=isw-is1+1, js_in=jsw-js1+1,err_msg=err_msg)
           if(id_sol_con>0) used = send_data(id_sol_con, solar_constant, Time ) 
        END DO
        !$OMP END parallel do
        CALL diag_send_complete(Time_step) 
     END DO
  END IF


  IF ( test_number == 14 ) THEN
     id_dat2_2d = register_diag_field('test_mod', 'dat2', (/id_lon2,id_lat2/), Time, 'sample data', 'K', err_msg=err_msg)
     IF ( err_msg /= '' ) THEN
        WRITE (out_unit,'(a)') 'test14 successful. err_msg='//TRIM(err_msg)
     ELSE
        WRITE (out_unit,'(a)') 'test14 fails.'
     END IF
  ELSE
     id_dat2_2d = register_diag_field('test_mod', 'dat2', (/id_lon2,id_lat2/), Time, 'sample data', 'K')
  END IF

  id_bk = register_static_field('test_diag_manager_mod', 'bk', (/id_pfull/), 'half level sigma', 'none')

  IF ( test_number == 13 ) THEN
     IF ( id_dat2_2d > 0 ) used=send_data(id_dat2_2d, dat2(:,:,1), Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test13: successful if a WARNING message appears that refers to output interval greater than runlength'
     ELSE
        WRITE (out_unit,'(a)') 'test13 fails: err_msg='//TRIM(err_msg)
     END IF
  END IF

  ! Note: test12 involves diag_manager_init, it does not require a call to send_data.
  !       See call to diag_manager_init above.

  IF ( test_number == 11 ) THEN
     is_in = 1+hi
     js_in = 1+hj
     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj

     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test11.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test11.1 fails. err_msg='//TRIM(err_msg)
     END IF

     ! intentional_error: je_in is missing
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test11.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test11.2 successful. err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 10 ) THEN
     !  1 window, no halos, static, 1 dimension, global data.

     IF ( id_bk > 0 ) used = send_data(id_bk, bk, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test10.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test10.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too large.
     IF ( id_bk > 0 ) used = send_data(id_bk, phalf, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE(out_unit,'(a)') 'test10.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test10.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 9 ) THEN
     !  1 window, no halos, static, 1 dimension, global data
     IF ( id_bk > 0 ) used = send_data(id_bk, bk, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test9.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test9.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small
     IF ( id_bk > 0 ) used = send_data(id_bk, bk(1:nlev-1), err_msg=err_msg) ! intentional_error
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test9.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test9.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 8 ) THEN
     !  1 window with halos
     is_in = 1+hi
     js_in = 1+hj

     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test8.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test8.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small in both x and y directions
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     ie_in = ie1-is1+1+hi
     je_in = je1-js1+1+hj
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat1h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat1h, Time, is_in=is_in, js_in=js_in, &
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test8.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test8.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 7 ) THEN
     !  1 window with halos
     is_in = 1+hi
     js_in = 1+hj

     ie_in = ie1-is1+1+hi
     je_in = je1-js1+1+hj
     IF ( id_dat1 > 0 ) used=send_data(id_dat1, dat1h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test7.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test7.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too large in both x and y directions
     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj
     IF ( id_dat1 > 0 ) used=send_data(id_dat1, dat2h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test7.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test7.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 6 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within either do loop for test6.1
     test_successful = .TRUE.
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,:,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test6.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,:,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test6.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test6.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test6.1 fails.'
     END IF

     !  intentional_error: data array too small in y direction
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,js2:je2-1,:), Time, i-is2+1, 1)
     END DO
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,js2:je2-1,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test6.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test6.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 5 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within either do loop for test5.1
     test_successful = .TRUE.
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(:,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test5.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(:,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test5.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test5.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test5.1 fails.'
     END IF

     !  intentional_error: data array too small in x direction.
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(is2:ie2-1,j:j,:), Time, 1, j-js2+1)
     END DO
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(is2:ie2-1,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test5.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test5.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 4 ) THEN
     !  1 window, no halos
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2, Time, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test4.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test4.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small in both x and y directions
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test4.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test4.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 3 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within do loop for test3.1
     test_successful = .TRUE.
     DO i=is1, ie1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1(i:i,:,:), Time, i-is1+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test3.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test3.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test3.1 fails.'
     END IF

     !  intentional_error: data array too large in y direction
     DO i=is1, ie1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2(i:i,:,:), Time, i-is1+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test3.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test3.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 2 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within do loop for test2.1
     test_successful = .TRUE.
     DO j=js1, je1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1(:,j:j,:), Time, 1, j-js1+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test2.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test2.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test2.1 fails.'
     END IF

     !  intentional_error: data array too large in x direction
     DO j=js1, je1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2(:,j:j,:), Time, 1, j-js1+1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test2.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test2.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 1 ) THEN
     !  1 window, no halos
     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test1.1 fails: Intentional error not detected'
     ELSE
        WRITE (out_unit,'(a)') 'test1.1 successful: '//TRIM(err_msg)
     END IF

     !  intentional_error: data array too large in both x and y directions
     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test1.2 successful'
     ELSE
        WRITE (out_unit,'(a)') 'test1.2 fails: '//TRIM(err_msg)
     END IF
  END IF

  CALL diag_manager_end(Time)
  CALL fms_io_exit
  CALL fms_end

CONTAINS

  SUBROUTINE compute_grid(nlon, nlat, is, ie, js, je, lon_global, lat_global, lonb_global, latb_global, lon, lat, lonb, latb)
    INTEGER, INTENT(in) :: nlon, nlat, is, ie, js, je
    REAL, INTENT(out), DIMENSION(:) :: lon_global, lat_global, lonb_global, latb_global, lon, lat, lonb, latb

    REAL :: dlon, dlat
    INTEGER :: i, j

    dlon = 2*PI/nlon
    dlat = PI/nlat

    DO i=1, nlon+1
       lonb_global(i) = dlon*(i-1)
    END DO
    DO j=1,nlat+1
       latb_global(j) = dlat*(j-1) - .5*PI
    END DO
    DO i=1,nlon
       lon_global(i) = .5*(lonb_global(i) + lonb_global(i+1))
    END DO
    DO j=1,nlat
       lat_global(j) = .5*(latb_global(j) + latb_global(j+1))
    END DO
    lon  = lon_global(is:ie)
    lat  = lat_global(js:je)
    lonb = lonb_global(is:ie+1)
    latb = latb_global(js:je+1)
  END SUBROUTINE compute_grid
END PROGRAM test
#endif
