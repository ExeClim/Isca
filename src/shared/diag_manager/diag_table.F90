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

MODULE diag_table_mod
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
  !   Seth Underwood
  ! </CONTACT>
  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/" />
  ! <OVERVIEW>
  !   <TT>diag_table_mod</TT> is a set of subroutines use to parse out the data from a <TT>diag_table</TT>.  This module
  !   will also setup the arrays required to store the information by counting the number of input fields, output files, and
  !   files.
  ! </OVERVIEW>
  ! <DESCRIPTION>
  !   <TT>diag_table_mod</TT> parses the <TT>diag_table</TT> file, and sets up the required arrays to hold the information
  !   needed for the <TT>diag_manager_mod</TT> to correctly write out the model history files.
  !
  !   The <I>diagnostics table</I> allows users to specify sampling rates and the choice of fields at run time.  The
  !   <TT>diag_table</TT> file consists of comma-separated ASCII values.  The <TT>diag_table</TT> essentially has three sections:
  !   <B>Global</B>, <B>File</B>, and <B>Field</B> sections.  The <B>Global</B> section must be the first two lines of the file,
  !   whereas the <B>File</B> and <B>Field</B> sections can be inter mixed to allow the file to be organized as desired.
  !   Comments can be added to the <TT>diag_table</TT> file by using the hash symbol (#) as the first character in the line.
  !
  !   All errors in the <TT>diag_table</TT> will throw a <TT>FATAL</TT> error.  A simple utility <TT>diag_table_chk</TT>has been
  !   added to the FRE tools suite to check a <TT>diag_table</TT> for errors.  A brief usage statement can be obtained by running
  !   <TT>diag_table_chk --help</TT>, and a man page like description can views by running <TT>perldoc diag_table_chk</TT>.
  !
  !   Below is a description of the three sections.
  !   <OL>
  !     <LI>
  !       <B>Global Section:</B>  The first two lines of the <TT>diag_table</TT> must contain the <I>title</I> and the <I>base
  !       date</I> of the experiment respectively.  The <I>title</I> must be a Fortran CHARACTER string.  The <I>base date</I>
  !       is the reference time used for the time units, and must be greater than or equal to the model start time.
  !       The <I>base date</I> consists of six space-separated integer in the following format.<BR />
  !       <TT><NOBR>year month day hour minute second</NOBR></TT><BR />
  !     </LI>
  !     <LI>
  !       <B>File Section:</B> File lines contain 6 required and 5 optional fields (optional fields are surrounded with
  !       square brackets ([]).  File lines can be intermixed with the field lines, but the file must be defined before any
  !       fields that are to be written to the file.  File lines have the following format:<BR />
  !       <PRELN>
  !         "file_name", output_freq, "output_freq_units", file_format, "time_axis_units", "time_axis_name"
  !         [, new_file_freq, "new_file_freq_units"[, "start_time"[, file_duration, "file_duration_units"]]]
  !       </PRELN>
  !       <BR />
  !       with the following descriptions.
  !       <DL>
  !         <DT><TT>CHARACTER(len=128) :: file_name</TT></DT>
  !         <DD>
  !           Output file name without the trailing "<TT>.nc</TT>".
  !
  !           A single file description can produce multiple files using special time string suffix keywords.  This time string
  !           will append the time strings to the base file name each time a new file is opened.  They syntax for the time string
  !           suffix keywords are <TT>%#tt</TT> Where <TT>#</TT> is a mandatory single digit number specifying the width of the
  !           field, and <TT>tt</TT> can be as follows:
  !           <NL>
  !             <LI><TT>yr</TT> <EN /> Years</LI>
  !             <LI><TT>mo</TT> <EN /> Months</LI>
  !             <LI><TT>dy</TT> <EN /> Days</LI>
  !             <LI><TT>hr</TT> <EN /> Hours</LI>
  !             <LI><TT>mi</TT> <EN /> Minutes</LI>
  !             <LI><TT>sc</TT> <EN /> Seconds</LI>
  !           </NL>
  !           Thus, a file name of <TT>file2_yr_dy%1yr%3dy</TT> will have a base file name of <TT>file2_yr_dy_1_001</TT> if the
  !           file is created on year 1 day 1 of the model run.  <B><I>NOTE:</I></B> The time suffix keywords must be used if the
  !           optional fields <TT>new_file_freq</TT> and <TT>new_file_freq_units</TT> are used, otherwise a <TT>FATAL</TT> error
  !           will occur.
  !         </DD>
  !
  !         <DT><TT>INTEGER :: output_freq</TT></DT>
  !         <DD>How often to write fields to file.
  !           <NL>
  !             <LI><TT>> 0</TT> <EN /> Output frequency in <TT>output_freq_units</TT>.</LI>
  !             <LI><TT>= 0</TT> <EN /> Output frequency every time set. (<TT>output_freq_units</TT> is ignored.)</LI>
  !             <LI><TT>=-1</TT> <EN /> Output at end of run only. (<TT>output_freq_units</TT> is ignored.)</LI>
  !           </NL>
  !         </DD>
  !         <DT><TT>CHARACTER(len=10) :: output_freq_units</TT></DT>
  !         <DD>
  !           Time units for output.  Can be either <TT>years</TT>, <TT>months</TT>, <TT>days</TT>, <TT>minutes</TT>,
  !           <TT>hours</TT>, or <TT>seconds</TT>.
  !         </DD>
  !         <DT><TT>INTEGER :: file_format</TT></DT>
  !         <DD>
  !           Output file format.  Currently only the <I>netCDF</I> file format is supported.
  !           <NL>
  !             <LI><TT>= 1</TT> <EN /> netCDF</LI>
  !           </NL>
  !         </DD>
  !         <DT><TT>CHARACTER(len=10) :: time_axis_units</TT></DT>
  !         <DD>
  !           Time units for the output file time axis.  Can be either <TT>years</TT>, <TT>months</TT>, <TT>days</TT>,
  !           <TT>minutes</TT>, <TT>hours</TT>, or <TT>seconds</TT>.
  !         </DD>
  !         <DT><TT>CHARACTER(len=128) :: time_axis_name</TT></DT>
  !         <DD>
  !           Axis name for the output file time axis.  The character sting must contain the string 'time'. (mixed upper and
  !           lowercase allowed.)
  !         </DD>
  !         <DT><TT>INTEGER, OPTIONAL :: new_file_freq</TT></DT>
  !         <DD>
  !           Frequency for closing the existing file, and creating a new file in <TT>new_file_freq_units</TT>.
  !         </DD>
  !         <DT><TT>CHARACTER(len=10), OPTIONAL :: new_file_freq_units</TT></DT>
  !         <DD>
  !           Time units for creating a new file.  Can be either <TT>years</TT>, <TT>months</TT>, <TT>days</TT>,
  !           <TT>minutes</TT>, <TT>hours</TT>, or <TT>seconds</TT>.  <B><I>NOTE:</I></B> If the <TT>new_file_freq</TT> field is
  !           present, then this field must also be present.
  !         </DD>
  !         <DT><TT>CHARACTER(len=25), OPTIONAL :: start_time</TT></DT>
  !         <DD>
  !           Time to start the file for the first time.  The format of this string is the same as the <I>global date</I>.  <B><I>
  !           NOTE:</I></B> The <TT>new_file_freq</TT> and the <TT>new_file_freq_units</TT> fields must be present to use this field.
  !         </DD>
  !         <DT><TT>INTEGER, OPTIONAL :: file_duration</TT></DT>
  !         <DD>
  !           How long file should receive data after start time in <TT>file_duration_units</TT>.  This optional field can only
  !           be used if the <TT>start_time</TT> field is present.  If this field is absent, then the file duration will be equal
  !           to the frequency for creating new files.  <B><I>NOTE:</I></B> The <TT>file_duration_units</TT> field must also be
  !           present if this field is present.
  !         </DD>
  !         <DT><TT>CHARACTER(len=10), OPTIONAL :: file_duration_units</TT></DT>
  !         <DD>
  !           File duration units. Can be either <TT>years</TT>, <TT>months</TT>, <TT>days</TT>,
  !           <TT>minutes</TT>, <TT>hours</TT>, or <TT>seconds</TT>.  <B><I>NOTE:</I></B> If the <TT>file_duration</TT> field is
  !           present, then this field must also be present.
  !         </DD>
  !       </DL>
  !     </LI>
  !     <LI>
  !       <B>Field Section:</B> Field lines contain 8 fields.  Field lines can be intermixed with file lines, but the file must
  !       be defined before any fields that are to be written to the file.  Fields line can contain fields that are not written
  !       to any files.  The file name for these fields is <TT>null</TT>.
  !
  !       Field lines have the following format:<BR />
  !       <PRE>
  ! "module_name", "field_name", "output_name", "file_name", "time_sampling", "reduction_method", "regional_section", packing
  !       </PRE>
  !       with the following descriptions.
  !       <DL>
  !         <DT><TT>CHARACTER(len=128) :: module_name</TT></DT>
  !         <DD>Module that contains the <TT>field_name</TT> variable.  (e.g. <TT>atmos_mod</TT>, <TT>land_mod</TT>)</DD>
  !         <DT><TT>CHARACTER(len=128) :: field_name</TT></DT>
  !         <DD>Module variable name that has data to be written to file.</DD>
  !         <DT><TT>CHARACTER(len=128) :: output_name</TT></DT>
  !         <DD>Name of the field as written in <TT>file_name</TT>.</DD>
  !         <DT><TT>CHARACTER(len=128) :: file_name</TT></DT>
  !         <DD>
  !           Name of the file where the field is to be written. <B><I>NOTE:</I></B> The file <TT>file_name</TT> must be
  !           defined first.
  !         </DD>
  !         <DT><TT>CHARACTER(len=50) :: time_sampling</TT></DT>
  !         <DD>Currently not used.  Please use the string "all".</DD>
  !         <DT><TT>CHARACTER(len=50) :: reduction_method</TT></DT>
  !         <DD>
  !           The data reduction method to perform prior to writing data to disk.  Valid options are (redundant names are
  !           separated with commas):
  !           <DL>
  !             <DT><TT>.TRUE.</TT>, average</DT>
  !             <DD>Average from the last time written to the current time.</DD>
  !             <DT><TT>.FALSE.</TT>, none</DT>
  !             <DD>No reduction performed.  Write current time step value only.</DD>
  !             <DT>min</DT> <DD>Minimum value from last write to current time.</DD>
  !             <DT>max</DT> <DD>Maximum value from last write to current time.</DD>
  !             <DT>diurnal##</DT> <DD>## diurnal averages</DD>
  !           </DL>
  !         </DD>
  !         <DT><TT>CHARACTER(len=50) :: regional_section</TT></DT>
  !         <DD>
  !           Bounds of the regional section to capture.  A value of <TT>none</TT> indicates a global region.  The regional
  !           section has the following format:<BR />
  !           <TT>lat_min, lat_max, lon_min, lon_max, vert_min, vert_max</TT><BR />
  !           Use <TT>vert_min = -1</TT> and <TT>vert_max = -1</TT> to get the entire vertical axis.  <B><I>NOTE:</I></B>
  !           Currently, the defined region <I>MUST</I> be confined to a single tile.
  !         </DD>
  !         <DT><TT>INTEGER :: packing</TT></DT>
  !         <DD>
  !           Fortran number <TT>KIND</TT> of the data written.  Valid values:
  !           <NL>
  !             <LI><TT>= 1</TT> <EN /> double precision</LI>
  !             <LI><TT>= 2</TT> <EN /> float</LI>
  !             <LI><TT>= 4</TT> <EN /> packed 16-bit integers</LI>
  !             <LI><TT>= 8</TT> <EN /> packed 1-byte (not tested).</LI>
  !           </NL>
  !         </DD>
  !       </DL>
  !     </LI>
  !   </OL>
  !
  !   <H4><B>Sample <TT>diag_table</TT></B></H4>
  !   <NL>
  !     <LI>
  !       <PRE>
  ! "diag manager test"
  ! 1999 1 1 0 0 0
  !
  ! #output files
  ! 10_days,               10, "days", 1, "hours", "Time"
  ! "file1_hr%hr3",         5, "days", 1, "hours", "Time", 15, "days"
  ! "file2_yr_dy%yr1%dy3",  5, "days", 1, "hours", "Time", 10, "days", "1 1 7 0 0 0"
  ! "file3_yr_dy%yr1%dy3",  5, "days", 1, "hours", "Time", 20, "days", "1 1 7 0 0 0", 5, "years"
  !
  ! #output variables
  ! "ice_mod", "ice", "ice", "10_days", "all", .false., "none", 2
  !
  ! # temp_local file and fields.
  ! temp_local, 1, "days", 1, "hours", "Time"
  ! "ocean_mod", "temp", "temp", "temp_local", "all", .FALSE., "5 259.5 -59.5 59.5 1 1", 2
  !       </PRE>
  !     </LI>
  !   </NL>
  !
  !   <H4>Useful Additional Utility</H4>
  !   A simple utility has been created to help discover
  ! </DESCRIPTION>
  USE mpp_io_mod, ONLY: mpp_open, MPP_RDONLY
  USE mpp_mod, ONLY: read_ascii_file
  USE fms_mod, ONLY: fms_error_handler, error_mesg, file_exist, stdlog, mpp_pe, mpp_root_pe, FATAL, WARNING, lowercase, close_file
  USE time_manager_mod, ONLY: get_calendar_type, NO_CALENDAR, set_date, set_time, month_name, time_type
  USE constants_mod, ONLY: SECONDS_PER_HOUR, SECONDS_PER_MINUTE
  
  USE diag_data_mod, ONLY: global_descriptor, base_time, base_year, base_month, base_day, base_hour, base_minute, base_second,&
       & DIAG_OTHER, DIAG_OCEAN, DIAG_ALL, coord_type, append_pelist_name, pelist_name, filename_appendix
  USE diag_util_mod, ONLY: init_file, check_duplicate_output_fields, init_input_field, init_output_field

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: parse_diag_table
  
  TYPE field_description_type
     CHARACTER(len=128) :: module_name, field_name, output_name, file_name
     CHARACTER(len=50) :: time_sampling
     CHARACTER(len=50) :: time_method   
     CHARACTER(len=50) :: spatial_ops
     TYPE(coord_type) :: regional_coords
     INTEGER :: pack
  END TYPE field_description_type

  TYPE file_description_type
     INTEGER :: output_freq 
     INTEGER :: file_format
     INTEGER :: new_file_freq
     INTEGER :: file_duration
     INTEGER :: iTime_units
     INTEGER :: iOutput_freq_units
     INTEGER :: iNew_file_freq_units
     INTEGER :: iFile_duration_units
     CHARACTER(len=128) :: file_name
     CHARACTER(len=10) :: output_freq_units
     CHARACTER(len=10) :: time_units
     CHARACTER(len=128) :: long_name
     CHARACTER(len=10) :: new_file_freq_units
     CHARACTER(len=25) :: start_time_s
     CHARACTER(len=10) :: file_duration_units
     TYPE(time_type) :: start_time
  END TYPE file_description_type

  CHARACTER(len=*), PARAMETER :: UNALLOWED_QTE = "'"//'"'
  CHARACTER(len=*), PARAMETER :: UNALLOWED_ALL = UNALLOWED_QTE//","

CONTAINS
  
  ! <SUBROUTINE NAME="parse_diag_table">
  !   <OVERVIEW>
  !     Parse the <TT>diag_table</TT> in preparation for diagnostic output.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE parse_diag_table(diag_subset, istat, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>parse_diag_table</TT> is the public interface to parse the diag_table, and setup the arrays needed to store the
  !     requested diagnostics from the <TT>diag_table</TT>.  <TT>parse_diag_table</TT> will return a non-zero <TT>istat</TT> if
  !     a problem parsing the <TT>diag_table</TT>.
  !
  !     NOT YET IMPLEMENTED: <TT>parse_diag_table</TT> will parse through the <TT>diag_table</TT> twice.  The first pass, will be
  !     to get a good "guess" of array sizes.  These arrays, that will hold the requested diagnostic fields and files, will then be
  !     allocated to the size of the "guess" plus a slight increase.
  !   </DESCRIPTION>
  !   <IN NAME="diag_subset" TYPE="INTEGER, OPTIONAL">
  !     Diagnostic sampling subset.
  !   </IN>
  !   <OUT NAME="iunit" TYPE="INTEGER, OPTIONAL">
  !     Status of parsing the <TT>diag_table</TT>.  A non-zero status indicates a problem parsing the table.
  !   </OUT>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL">
  !     Error message corresponding to the <TT>istat</TT> return value.
  !   </OUT>
  !   <ERROR STATUS="FATAL">
  !     diag_table file does not exist.
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Error reading the global descriptor from the diagnostic table.
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Error reading the base date from the diagnostic table.
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     The base_year/month/day can not equal zero
  !   </ERROR>
  !   <ERROR STATUS="WARNING">
  !     Problem reading diag_table, line numbers in errors may be incorrect.
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Problem reading the diag_table (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Incorrect file description FORMAT in diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Invalid file FORMAT for file description in the diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Invalid time axis units in diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Invalid output frequency units in diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Invalid NEW file frequency units in diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Invalid file duration units in diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Invalid start time in the file description in diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Field description FORMAT is incorrect in diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Packing is out of range for the field description in diag_table. (line: <line_number>)
  !   </ERROR>
  !   <ERROR STATUS="FATAL">
  !     Error in regional output description for field description in diag_table. (line: <line_number>)
  !   </ERROR>
  SUBROUTINE parse_diag_table(diag_subset, istat, err_msg)
    INTEGER, INTENT(in), OPTIONAL :: diag_subset
    INTEGER, INTENT(out), OPTIONAL, TARGET :: istat
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    INTEGER, PARAMETER :: DT_LINE_LENGTH = 512

    INTEGER :: stdlog_unit !< Fortran file unit number for the stdlog file.
    INTEGER :: record_len !< String length of the diag_table line read in.
    INTEGER :: num_lines !< Number of lines in diag_table
    INTEGER :: line_num !< Integer representation of the line number.
    INTEGER :: commentStart !< Index location of first '#' on line
    INTEGER :: diag_subset_output !< local value of diag_subset
    INTEGER :: nfields, nfiles !< Number of fields and files.  Not used yet.
    INTEGER, TARGET :: mystat !< variable to hold return status of function/subroutine calls.
    INTEGER, POINTER :: pstat !< pointer that points to istat if preset, otherwise, points to mystat.

    CHARACTER(len=5) :: line_number !< String representation of the line number.
    CHARACTER(len=9) :: amonth !< Month name
    CHARACTER(len=256) :: record_line !< Current line from the diag_table.
    CHARACTER(len=256) :: local_err_msg !< Sting to hold local error messages.
    CHARACTER(len=DT_LINE_LENGTH), DIMENSION(:), ALLOCATABLE :: diag_table

    TYPE(file_description_type) :: temp_file
    TYPE(field_description_type) :: temp_field

    ! set up the pstat pointer
    IF ( PRESENT(istat) ) THEN
       pstat => istat
    ELSE
       pstat => mystat
    END IF
    ! Default return value (success)
    pstat = 0

    IF ( PRESENT(diag_subset) ) THEN
       diag_subset_output = diag_subset
    ELSE 
       diag_subset_output = DIAG_ALL
    END IF
    
    ! get the stdlog unit number
    stdlog_unit = stdlog()

    call read_ascii_file('diag_table', DT_LINE_LENGTH, diag_table, num_lines)
    
    ! Read in the global file labeling string
    READ (UNIT=diag_table(1), FMT=*, IOSTAT=mystat) global_descriptor
    IF ( mystat /= 0 ) THEN
       pstat = mystat
       IF ( fms_error_handler('diag_table_mod::parse_diag_table', 'Error reading the global descriptor from the diagnostic table.',&
            & err_msg) ) RETURN
    END IF
    
    ! Read in the base date
    READ (UNIT=diag_table(2), FMT=*, IOSTAT=mystat) base_year, base_month, base_day, base_hour, base_minute, base_second
    IF ( mystat /= 0 ) THEN
       pstat = mystat
       IF ( fms_error_handler('diag_manager_init', 'Error reading the base date from the diagnostic table.', err_msg) ) RETURN
    END IF
    
    ! Set up the time type for base time
    IF ( get_calendar_type() /= NO_CALENDAR ) THEN
       IF ( base_year==0 .OR. base_month==0 .OR. base_day==0 ) THEN
          pstat = 101
          IF ( fms_error_handler('diag_table_mod::parse_diag_table', 'The base_year/month/day can not equal zero', err_msg) ) RETURN
       END IF
       base_time = set_date(base_year, base_month, base_day, base_hour, base_minute, base_second)
       amonth = month_name(base_month)
    ELSE
       ! No calendar - ignore year and month
       base_time = set_time(NINT(base_hour*SECONDS_PER_HOUR)+NINT(base_minute*SECONDS_PER_MINUTE)+base_second, base_day)
       base_year = 0
       base_month = 0
       amonth = 'day'
    END IF

    IF ( mpp_pe() == mpp_root_pe() ) THEN
       WRITE (stdlog_unit,'("base date used = ",I4,1X,A,2I3,2(":",I2.2)," gmt")') base_year, TRIM(amonth), base_day, &
            & base_hour, base_minute, base_second
    END IF

    nfiles=0
    nfields=0
    parser: DO line_num=3, num_lines
       ! Read in the entire line from the file.
       ! If there is a read error, give a warning, and
       ! cycle the parser loop.
       READ (diag_table(line_num), FMT='(A)', IOSTAT=mystat) record_line
       ! Increase line counter, and put in string for use in warning/error messages.
       WRITE (line_number, '(I5)') line_num

       IF ( mystat > 0 ) THEN
          IF ( mpp_pe() == mpp_root_pe() ) &
               & CALL error_mesg("diag_table_mod::parse_diag_table",&
               & "Problem reading the diag_table (line:" //line_number//").", FATAL)
          CYCLE parser
       ELSE IF ( mystat < 0 ) THEN
          EXIT parser
       END IF
       
       ! How long is the read in string?
       record_len = LEN_TRIM(record_line)

       ! ignore blank lines and  lines with comments only (comment marker '#')
       commentStart = INDEX(record_line,'#')
       IF ( commentStart .NE. 0 ) record_line = record_line(1:commentStart-1)
       IF ( LEN_TRIM(record_line) == 0 .OR. record_len == 0 ) CYCLE parser

       init: IF ( is_a_file(TRIM(record_line)) ) THEN
          temp_file = parse_file_line(LINE=record_line, ISTAT=mystat, ERR_MSG=local_err_msg)
          
          IF ( mystat > 0 ) THEN 
             CALL error_mesg("diag_table_mod::parse_diag_table",&
                  & TRIM(local_err_msg)//" (line:" //TRIM(line_number)//").", FATAL)
          ELSE IF ( mystat < 0 ) THEN
             IF ( mpp_pe() == mpp_root_pe() )&
                  & CALL error_mesg("diag_table_mod::parse_diag_table",&
                  & TRIM(local_err_msg)//" (line: "//TRIM(line_number)//").", WARNING)
             CYCLE parser
          ELSE IF ( (diag_subset_output == DIAG_OTHER .AND. VERIFY('ocean', lowercase(temp_file%file_name)) == 0).OR.&
               &    (diag_subset_output == DIAG_OCEAN .AND. VERIFY('ocean', lowercase(temp_file%file_name)) /= 0) ) THEN
             CYCLE parser
          ELSE IF ( temp_file%new_file_freq > 0 ) THEN ! Call the init_file subroutine.  The '1' is for the tile_count
             CALL init_file(temp_file%file_name, temp_file%output_freq, temp_file%iOutput_freq_units, temp_file%file_format,&
                  & temp_file%iTime_units, temp_file%long_name, 1, temp_file%new_file_freq, temp_file%iNew_file_freq_units,&
                  & temp_file%start_time, temp_file%file_duration, temp_file%iFile_duration_units)
          ELSE
             CALL init_file(temp_file%file_name, temp_file%output_freq, temp_file%iOutput_freq_units, temp_file%file_format,&
                  & temp_file%iTime_units, temp_file%long_name, 1)
          END IF
          
          ! Increment number of files
          nfiles = nfiles + 1
       ELSE ! We have a field.
          temp_field = parse_field_line(LINE=record_line, ISTAT=mystat, ERR_MSG=local_err_msg)

          ! Check for errors, then initialize the input and output field
          IF (  mystat > 0 ) THEN
             CALL error_mesg("diag_table_mod::parse_diag_table",&
                  & TRIM(local_err_msg)//" (line: "//TRIM(line_number)//").",FATAL)
          ELSE IF ( mystat < 0 ) THEN
             IF ( mpp_pe() == mpp_root_pe() )&
                  & CALL error_mesg("diag_table_mod::Parse_diag_table",&
                  & TRIM(local_err_msg)//" (line: "//TRIM(line_number)//").",WARNING)
             CYCLE parser
          ELSE IF ( (diag_subset_output == DIAG_OTHER .AND. VERIFY('ocean', lowercase(temp_field%file_name)) == 0).OR.&
               &    (diag_subset_output == DIAG_OCEAN .AND. VERIFY('ocean', lowercase(temp_field%file_name)) /= 0) ) THEN 
             CYCLE parser
          ELSE IF ( lowercase(TRIM(temp_field%spatial_ops)) == 'none' ) THEN
             CALL init_input_field(temp_field%module_name, temp_field%field_name, 1)
             CALL init_output_field(temp_field%module_name, temp_field%field_name, temp_field%output_name, temp_field%file_name,&
                  & temp_field%time_method, temp_field%pack, 1)
          ELSE 
             CALL init_input_field(temp_field%module_name, temp_field%field_name, 1)
             CALL init_output_field(temp_field%module_name, temp_field%field_name, temp_field%output_name, temp_field%file_name,&
                  & temp_field%time_method, temp_field%pack, 1, temp_field%regional_coords)
          END IF

          ! Increment number of fields
          nfields = nfields + 1
       END IF init
    END DO parser

    ! Close the diag_table file.
    deallocate(diag_table)

    ! check duplicate output_fields in the diag_table
    CALL check_duplicate_output_fields(ERR_MSG=local_err_msg)
    IF ( local_err_msg /= '' ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_diag_table', TRIM(local_err_msg), err_msg) ) RETURN
    END IF

  END SUBROUTINE parse_diag_table
  ! </SUBROUTINE>

  ! <PRIVATE>
  ! <SUBROUTINE NAME="open_diag_table">
  !   <OVERVIEW>
  !     Open the diag_table file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE open_diag_table(iunit, iostat)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Open the <TT>diag_table</TT> file, and return the Fortran file unit number.
  !   </DESCRIPTION>
  !   <OUT NAME="iunit" TYPE="INTEGER">Fortran file unit number of the <TT>diag_table</TT>.</OUT>
  !   <IN NAME="iostat" TYPE="INTEGER, OPTIONAL">
  !     Status of opening file.  If iostat == 0, file exists.  If iostat > 0, the diag_table file does not exist.
  !   </IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL">
  !     String to hold the return error message.
  !   </OUT>
  SUBROUTINE open_diag_table(iunit, iostat, err_msg)
    INTEGER, INTENT(out) :: iunit
    INTEGER, INTENT(out), OPTIONAL, TARGET :: iostat
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    INTEGER, TARGET :: mystat
    INTEGER, POINTER :: pstat

    IF ( PRESENT(iostat) ) THEN
       pstat => iostat
    ELSE 
       pstat => mystat
    END IF
    
    IF ( .NOT.file_exist('diag_table') ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::open_diag_table',&
            & 'diag_table file does not exist.', err_msg) ) RETURN
    ELSE 
       pstat = 0
    END IF

    CALL mpp_open(iunit, 'diag_table', action=MPP_RDONLY)
  END SUBROUTINE open_diag_table
  ! </SUBROUTINE>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <SUBROUTINE NAME="close_diag_table">
  !   <OVERVIEW>
  !     Close the diag_table file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE close_diag_table(iunit)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Closes the diag_table file.
  !   </DESCRIPTION>
  !   <IN NAME="iunit" TYPE="INTEGER">Fortran file unit number of the <TT>diag_table</TT>.</IN>
  SUBROUTINE close_diag_table(iunit)
    INTEGER, INTENT(in) :: iunit

    CALL close_file(iunit)
  END SUBROUTINE close_diag_table
  ! </SUBROUTINE>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="parse_file_line">
  !   <OVERVIEW>
  !     Parse a file description line from the <TT>diag_table</TT> file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     TYPE(file_description_type) FUNCTION parse_file_line(line, istat, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>parse_file_line</TT> parses a file description line from the <TT>diag_table</TT> file, and returns a
  !     <TT>TYPE(file_description_type)</TT>.  The calling function, would then need to call the <TT>init_file</TT> to initialize
  !     the diagnostic output file.
  !   </DESCRIPTION>
  !   <IN NAME="line" TYPE="CHARACTER(len=*)">Line to parse from the <TT>diag_table</TT> file.</IN>
  !   <OUT NAME="istat" TYPE="INTEGER, OPTIONAL">
  !     Return state of the function.  A value of 0 indicates success.  A positive value indicates a <TT>FATAL</TT> error occurred,
  !     and a negative value indicates a <TT>WARNING</TT> should be issued.
  !   </OUT>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL">
  !     Error string to include in the <TT>FATAL</TT> or <TT>WARNING</TT> message.
  !   </OUT>
  TYPE(file_description_type) FUNCTION parse_file_line(line, istat, err_msg)
    CHARACTER(len=*), INTENT(in) :: line
    INTEGER, INTENT(out), OPTIONAL, TARGET :: istat
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    INTEGER, TARGET :: mystat
    INTEGER, POINTER :: pstat
    INTEGER :: year, month, day, hour, minute, second
    CHARACTER(len=256) :: local_err_msg !< Hold the return error message from routine calls.

    IF ( PRESENT(istat) ) THEN
       pstat => istat
    ELSE
       pstat => mystat
    END IF
    pstat = 0 ! default success return value

    ! Initialize the optional file description fields.
    parse_file_line%new_file_freq = 0
    parse_file_line%new_file_freq_units = ''
    parse_file_line%start_time_s = ''
    parse_file_line%file_duration = 0
    parse_file_line%file_duration_units = ''

    ! Read in the file description line..
    READ (line, FMT=*, IOSTAT=mystat) parse_file_line%file_name, parse_file_line%output_freq, parse_file_line%output_freq_units,&
         & parse_file_line%file_format, parse_file_line%time_units, parse_file_line%long_name,&
         & parse_file_line%new_file_freq, parse_file_line%new_file_freq_units, parse_file_line%start_time_s,&
         & parse_file_line%file_duration, parse_file_line%file_duration_units
    IF ( mystat > 0 ) THEN
       pstat = mystat
       IF ( fms_error_handler('diag_table_mod::parse_file_line', 'Incorrect file description format in diag_table.', err_msg) )&
            & RETURN
    END IF

    ! Check for unallowed characters in strings
    IF ( SCAN(parse_file_line%file_name, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line',&
            & 'Unallowed character in file_name in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_file_line%output_freq_units, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line',&
            & 'Unallowed character in output_freq_units in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_file_line%time_units, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line',&
            & 'Unallowed character in time_units in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_file_line%long_name, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line',&
            & 'Unallowed character in long_name in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_file_line%new_file_freq_units, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line',&
            & 'Unallowed character in new_file_freq_units in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_file_line%start_time_s, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line',&
            & 'Unallowed character in start_time_s in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_file_line%file_duration_units, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line',&
            & 'Unallowed character in file_duration_units in the diag_table.', err_msg) ) RETURN
    END IF

            
    ! Fix the file name
    parse_file_line%file_name = fix_file_name(TRIM(parse_file_line%file_name))

    ! Verify values / formats are correct
    IF ( parse_file_line%file_format > 2 .OR. parse_file_line%file_format < 1 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line', 'Invalid file format for file description in the diag_table.',&
            & err_msg) ) RETURN
    END IF
    
    ! check for known units
    parse_file_line%iTime_units = find_unit_ivalue(parse_file_line%time_units)
    parse_file_line%iOutput_freq_units = find_unit_ivalue(parse_file_line%output_freq_units)
    parse_file_line%iNew_file_freq_units = find_unit_ivalue(parse_file_line%new_file_freq_units)
    parse_file_line%iFile_duration_units = find_unit_ivalue(parse_file_line%file_duration_units)
    ! Verify the units are valid
    IF ( parse_file_line%iTime_units < 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line', 'Invalid time axis units in diag_table.', err_msg) )&
            & RETURN
    END IF
    IF ( parse_file_line%iOutput_freq_units < 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line', 'Invalid output frequency units in diag_table.', err_msg) )&
            & RETURN
    END IF
    IF ( parse_file_line%iNew_file_freq_units < 0 .AND. parse_file_line%new_file_freq > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line', 'Invalid new file frequency units in diag_table.', err_msg) )&
            & RETURN
    END IF
    IF ( parse_file_line%iFile_duration_units < 0 .AND. parse_file_line%file_duration > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_file_line', 'Invalid file duration units in diag_table.', err_msg) )&
            & RETURN
    END IF

    !::sdu::
    !::sdu:: Here is where we would want to parse the regional/global string
    !::sdu::

    ! Check for file frequency, start time and duration presence.
    ! This will determine how the init subroutine is called.
    new_file_freq_present: IF ( parse_file_line%new_file_freq > 0 ) THEN ! New file frequency present.
       IF ( LEN_TRIM(parse_file_line%start_time_s) > 0 ) THEN ! start time present
          READ (parse_file_line%start_time_s, FMT=*, IOSTAT=mystat) year, month, day, hour, minute, second
          IF ( mystat /= 0 ) THEN 
             pstat = 1
             IF ( fms_error_handler('diag_table_mod::parse_file_line',&
                  & 'Invalid start time in the file description in diag_table.', err_msg) ) RETURN
          END IF
          parse_file_line%start_time = set_date(year, month, day, hour, minute, second, err_msg=local_err_msg)
          IF ( local_err_msg /= '' ) THEN
             pstat = 1
             IF ( fms_error_handler('diag_table_mod::parse_file_line', local_err_msg, err_msg) ) RETURN
          END IF
          IF ( parse_file_line%file_duration <= 0 ) THEN ! file_duration not present
             parse_file_line%file_duration = parse_file_line%new_file_freq
             parse_file_line%iFile_duration_units = parse_file_line%iNew_file_freq_units
          END IF
       ELSE
          parse_file_line%start_time = base_time
          parse_file_line%file_duration = parse_file_line%new_file_freq
          parse_file_line%iFile_duration_units = parse_file_line%iNew_file_freq_units
       END IF
    END IF new_file_freq_present

  END FUNCTION parse_file_line
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="parse_field_line">
  !   <OVERVIEW>
  !     Parse a field description line from the <TT>diag_table</TT> file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     TYPE(field_description_type) FUNCTION parse_field_line(line, istat, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>parse_field_line</TT> parses a field description line from the <TT>diag_table</TT> file, and returns a
  !     <TT>TYPE(field_description_type)</TT>.  The calling function, would then need to call the <TT>init_input_field</TT> and
  !     <TT>init_output_field</TT> to initialize the diagnostic output field.
  !   </DESCRIPTION>
  !   <IN NAME="line" TYPE="CHARACTER(len=*)">Line to parse from the <TT>diag_table</TT> file.</IN>
  !   <OUT NAME="istat" TYPE="INTEGER, OPTIONAL">
  !     Return state of the function.  A value of 0 indicates success.  A positive value indicates a <TT>FATAL</TT> error occurred,
  !     and a negative value indicates a <TT>WARNING</TT> should be issued.
  !   </OUT>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL">
  !     Error string to include in the <TT>FATAL</TT> or <TT>WARNING</TT> message.
  !   </OUT>
  TYPE(field_description_type) FUNCTION parse_field_line(line, istat, err_msg)
    CHARACTER(len=*), INTENT(in) :: line
    INTEGER, INTENT(out), OPTIONAL, TARGET :: istat
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: err_msg

    INTEGER, TARGET :: mystat
    INTEGER, POINTER :: pstat

    IF ( PRESENT(istat) ) THEN
       pstat => istat
    ELSE
       pstat => mystat
    END IF
    pstat = 0 ! default success return value

    READ (line, FMT=*, IOSTAT=mystat) parse_field_line%module_name, parse_field_line%field_name, parse_field_line%output_name,&
         & parse_field_line%file_name, parse_field_line%time_sampling, parse_field_line%time_method, parse_field_line%spatial_ops,&
         & parse_field_line%pack
    IF ( mystat /= 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Field description format is incorrect in diag_table.', err_msg) ) RETURN
    END IF

    ! Check for unallowed characters in the string
    IF ( SCAN(parse_field_line%module_name, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Unallowed Unallowed character in module_name in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_field_line%field_name, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Unallowed Unallowed character in field_name in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_field_line%output_name, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Unallowed Unallowed character in output_name in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_field_line%file_name, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Unallowed Unallowed character in file_name in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_field_line%time_sampling, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Unallowed Unallowed character in time_sampling in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_field_line%time_method, UNALLOWED_ALL) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Unallowed Unallowed character in time_method in the diag_table.', err_msg) ) RETURN
    END IF
    IF ( SCAN(parse_field_line%spatial_ops, UNALLOWED_QTE) > 0 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Unallowed Unallowed character in spatial_ops in the diag_table.', err_msg) ) RETURN
    END IF

    ! Fix the file name
    ! Removes any added '.nc' and appends additional information.
    parse_field_line%file_name = fix_file_name(TRIM(parse_field_line%file_name))

    IF ( parse_field_line%pack > 8 .OR. parse_field_line%pack < 1 ) THEN
       pstat = 1
       IF ( fms_error_handler('diag_table_mod::parse_field_line',&
            & 'Packing is out of range for the field description in diag_table.', err_msg) ) RETURN
    END IF
    
    IF ( lowercase(TRIM(parse_field_line%spatial_ops)) /= 'none' ) THEN
       READ (parse_field_line%spatial_ops, FMT=*, IOSTAT=mystat) parse_field_line%regional_coords
       IF ( mystat /= 0 ) THEN
          IF ( fms_error_handler('diag_table_mod::parse_field_line',&
               & 'Error in regional output description for field description in diag_table.', err_msg) ) RETURN
       END IF
    END IF
  END FUNCTION parse_field_line
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="is_a_file">
  !   <OVERVIEW>
  !     Determines if a line from the diag_table file is a file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE LOGICAL FUNCTION is_a_file(line)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>is_a_file</TT> checks a diag_table line to determine if the line describes a file.  If the line describes a file, the
  !     <TT>is_a_file</TT> will return <TT>.TRUE.</TT>.  Otherwise, it will return <TT>.FALSE.</TT>
  !   </DESCRIPTION>
  !   <IN NAME="line" TYPE="CARACTER(len=*)">String containing the <TT>diag_table</TT> line.</IN>
  PURE LOGICAL FUNCTION is_a_file(line)
    CHARACTER(len=*), INTENT(in) :: line

    CHARACTER(len=5) :: first
    INTEGER :: second 
    INTEGER :: mystat !< IO status from read

#if defined __PATHSCALE__ || defined _CRAYFTN
    ! This portion is to 'fix' pathscale's and Cray's Fortran compilers inability to handle the FMT=* correctly in the read
    ! statement.
    CHARACTER(len=10) :: secondString
    INTEGER :: comma1, comma2, linelen

    linelen = LEN(line)
    comma1 = INDEX(line,',') + 1 ! +1 to go past the comma
    comma2 = INDEX(line(comma1:linelen),',') + comma1 - 2 ! -2 to get rid of +1 in comma1 and to get 1 character before the comma

    secondString = ADJUSTL(line(comma1:comma2))
    READ (UNIT=secondString, FMT='(I)', IOSTAT=mystat) second
#else
    READ (UNIT=line, FMT=*, IOSTAT=mystat) first, second
#endif

    ! The line is a file if my status is zero after the read.
    is_a_file = mystat == 0
  END FUNCTION is_a_file
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="fix_file_name(file_name_string)">
  !   <OVERVIEW>
  !     Fixes the file name for use with diagnostic file and field initializations.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE CHARACTER(len=128) FUNCTION fix_file_name(file_name_string)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Removes any trailing '.nc' and appends to the file name additional information
  !     depending on if we are running an ensemble, or requesting append_pelist_name.
  !     
  !     Presently, the ensemble appendix will override the append_pelist_name variable.
  !   </DESCRIPTION>
  !   <IN NAME="file_name_string" TYPE="CHARACTER(len=*)">String containing the file name from the <TT>diag_table</TT>.</IN>
  PURE CHARACTER(len=128) FUNCTION fix_file_name(file_name_string)
    CHARACTER(len=*), INTENT(IN) :: file_name_string

    INTEGER :: file_name_len

    fix_file_name = file_name_string ! Default return value

    file_name_len = LEN_TRIM(file_name_string)

    ! Remove trailing '.nc' from the file_name, and append suffixes
    IF ( file_name_len > 2 ) THEN 
       IF ( file_name_string(file_name_len-2:file_name_len) == '.nc' ) THEN
          fix_file_name = file_name_string(1:file_name_len-3)
          file_name_len = file_name_len - 3
       END IF
    END IF
       
    ! If using ensembles, then append the ensemble information
    ! Or add the optional suffix based on the pe list name if the
    ! append_pelist_name == .TRUE.
    IF ( LEN_TRIM(filename_appendix) > 0 ) THEN 
       fix_file_name(file_name_len+1:) = TRIM(filename_appendix)    
    ELSE IF ( append_pelist_name ) THEN
       fix_file_name(file_name_len+1:) = TRIM(pelist_name)
    END IF
  END FUNCTION fix_file_name
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="find_unit_ivalue">
  !   <OVERVIEW>
  !     Return the integer value for the given time unit.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     PURE INTEGER FUNCTION find_unit_ivalue(unit_string)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Returns the corresponding integer value for the given time unit.
  !     <UL>
  !       <LI> seconds = 1 </LI>
  !       <LI> minutes = 2 </LI>
  !       <LI> hours = 3 </LI>
  !       <LI> days = 4 </LI>
  !       <LI> months = 5 </LI>
  !       <LI> years = 6 </LI>
  !       <LI> unknown = -1 </LI>
  !     </UL>
  !   </DESCRIPTION>
  !   <IN NAME="unit_string" TYPE="CHARACTER(len=*)">String containing the unit.</IN>
  PURE INTEGER FUNCTION find_unit_ivalue(unit_string)
       CHARACTER(len=*), INTENT(IN) :: unit_string !< Input string, containing the unit.

    SELECT CASE (TRIM(unit_string))
    CASE ('seconds')
       find_unit_ivalue = 1
    CASE ('minutes')
       find_unit_ivalue = 2
    CASE ('hours')
       find_unit_ivalue = 3
    CASE ('days')
       find_unit_ivalue = 4
    CASE ('months') 
       find_unit_ivalue = 5
    CASE ('years')
       find_unit_ivalue = 6
    CASE DEFAULT
       find_unit_ivalue = -1 ! Return statement if an incorrect / unknown unit used.
    END SELECT
  END FUNCTION find_unit_ivalue
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <SUBROUTINE NAME="initialize_output_arrays">
  !   <OVERVIEW>
  !     Allocate the file, in and out field arrays after reading the <TT>diag_table</TT> file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE initialize_output_arrays()
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     After reading in the <TT>diag_table</TT> file, the arrays that will hold the file, in, and out field data need to be
  !     allocated.  This routine will determine the size of the arrays, and then allocate the arrays.
  !   </DESCRIPTION>
  SUBROUTINE initialize_output_arrays()
    ! Place Holder
  END SUBROUTINE initialize_output_arrays
  ! </SUBROUTINE>
  ! </PRIVATE>
END MODULE diag_table_mod

