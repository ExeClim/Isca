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

MODULE diag_output_mod
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
  !   Seth Underwood
  ! </CONTACT>

  ! <OVERVIEW> <TT>diag_output_mod</TT> is an integral part of 
  !   <TT>diag_manager_mod</TT>. Its function is to write axis-meta-data, 
  !   field-meta-data and field data
  ! </OVERVIEW>

  USE mpp_io_mod, ONLY: axistype, fieldtype, mpp_io_init, mpp_open,  mpp_write_meta,&
       & mpp_write, mpp_flush, mpp_close, mpp_get_id, MPP_WRONLY, MPP_OVERWR,&
       & MPP_NETCDF, MPP_MULTI, MPP_SINGLE
  USE mpp_domains_mod, ONLY: domain1d, domain2d, mpp_define_domains, mpp_get_pelist,&
       &  mpp_get_global_domain, mpp_get_compute_domains, null_domain1d, null_domain2d,&
       & OPERATOR(.NE.), mpp_get_layout, OPERATOR(.EQ.)
  USE mpp_mod, ONLY: mpp_npes, mpp_pe
  USE diag_axis_mod, ONLY: diag_axis_init, get_diag_axis, get_axis_length,&
       & get_axis_global_length, get_domain1d, get_domain2d, get_axis_aux, get_tile_count
  USE diag_data_mod, ONLY: diag_fieldtype, diag_global_att_type, CMOR_MISSING_VALUE
  USE time_manager_mod, ONLY: get_calendar_type, valid_calendar_types
  USE fms_mod, ONLY: error_mesg, mpp_pe, write_version_number, FATAL

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: diag_output_init, write_axis_meta_data, write_field_meta_data, done_meta_data,&
       & diag_field_out, diag_flush, diag_fieldtype, get_diag_global_att, set_diag_global_att

  TYPE(diag_global_att_type), SAVE :: diag_global_att

  INTEGER, PARAMETER      :: NETCDF1 = 1
  INTEGER, PARAMETER      :: mxch  = 128
  INTEGER, PARAMETER      :: mxchl = 256
  INTEGER                 :: current_file_unit = -1
  INTEGER, DIMENSION(2,2) :: max_range = RESHAPE((/ -32767, 32767, -127,   127 /),(/2,2/))
!  DATA max_range / -32767, 32767, -127,   127 /
  INTEGER, DIMENSION(2)   :: missval = (/ -32768, -128 /)
  
  INTEGER, PARAMETER      :: max_axis_num = 20
  INTEGER                 :: num_axis_in_file = 0
  INTEGER, DIMENSION(max_axis_num) :: axis_in_file   
  LOGICAL, DIMENSION(max_axis_num) :: time_axis_flag, edge_axis_flag
  TYPE(axistype), DIMENSION(max_axis_num), SAVE :: Axis_types

  LOGICAL :: module_is_initialized = .FALSE.

  CHARACTER(len=128), PRIVATE :: version= &
       '$Id: diag_output.F90,v 19.0.10.3 2012/10/19 19:01:31 Zhi.Liang Exp $'
  CHARACTER(len=128), PRIVATE :: tagname= &
       '$Name: siena_201211 $'

CONTAINS

  ! <SUBROUTINE NAME="diag_output_init">
  !   <OVERVIEW>
  !     Registers the time axis and opens the output file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE diag_output_init (file_name, format, file_title, file_unit,
  !      all_scalar_or_1d, domain)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Registers the time axis, and opens the file for output.
  !   </DESCRIPTION>
  !   <IN NAME="file_name" TYPE="CHARACTER(len=*)">Output file name</IN>
  !   <IN NAME="format" TYPE="INTEGER">File format (Currently only 'NETCDF' is valid)</IN>
  !   <IN NAME="file_title" TYPE="CHARACTER(len=*)">Descriptive title for the file</IN>
  !   <OUT NAME="file_unit" TYPE="INTEGER">
  !     File unit number assigned to the output file.  Needed for subsuquent calls to
  !     <TT>diag_output_mod</TT>
  !   </OUT>
  !   <IN NAME="all_scalar_or_1d" TYPE="LOGICAL" />
  !   <IN NAME="domain" TYPE="TYPE(domain2d)" />
  SUBROUTINE diag_output_init(file_name, FORMAT, file_title, file_unit,&
       & all_scalar_or_1d, domain)
    CHARACTER(len=*), INTENT(in)  :: file_name, file_title
    INTEGER         , INTENT(in)  :: FORMAT
    INTEGER         , INTENT(out) :: file_unit
    LOGICAL         , INTENT(in)  :: all_scalar_or_1d
    TYPE(domain2d)  , INTENT(in)  :: domain

    INTEGER :: form, threading, fileset
    TYPE(diag_global_att_type) :: gAtt

    !---- initialize mpp_io ----
    IF ( .NOT.module_is_initialized ) THEN
       CALL mpp_io_init ()
       module_is_initialized = .TRUE.
    END IF
    CALL write_version_number( version, tagname )
   
    !---- set up output file ----
    SELECT CASE (FORMAT)
    CASE (NETCDF1)
       form      = MPP_NETCDF
       threading = MPP_MULTI
       fileset   = MPP_MULTI
    CASE default
       ! <ERROR STATUS="FATAL">invalid format</ERROR>
       CALL error_mesg('diag_output_init', 'invalid format', FATAL)
    END SELECT

    IF(all_scalar_or_1d) THEN
       threading = MPP_SINGLE
       fileset   = MPP_SINGLE
    END IF

    !---- open output file (return file_unit id) -----
    IF ( domain .EQ. NULL_DOMAIN2D ) THEN
       CALL mpp_open(file_unit, file_name, action=MPP_OVERWR, form=form,&
            & threading=threading, fileset=fileset)
    ELSE
       CALL mpp_open(file_unit, file_name, action=MPP_OVERWR, form=form,&
            & threading=threading, fileset=fileset, domain=domain) 
    END IF

    !---- write global attributes ----
    IF ( file_title(1:1) /= ' ' ) THEN
       CALL mpp_write_meta(file_unit, 'title', cval=TRIM(file_title))
    END IF

    !---- write grid type (mosaic or regular)
    CALL get_diag_global_att(gAtt)
    CALL mpp_write_meta(file_unit, 'grid_type', cval=TRIM(gAtt%grid_type))
    CALL mpp_write_meta(file_unit, 'grid_tile', cval=TRIM(gAtt%tile_name))

  END SUBROUTINE diag_output_init
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="write_axis_meta_data">
  !   <OVERVIEW>
  !     Write the axes meta data to file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE write_axis_meta_data(file_unit, axes, time_ops)
  !   </TEMPLATE>
  !   <IN NAME="file_unit" TYPE="INTEGER">File unit number</IN>
  !   <IN NAME="axes" TYPE="INTEGER, DIMENSION(:)">Array of axis ID's, including the time axis</IN>
  !   <IN NAME="time_ops" TYPE="LOGICAL, OPTIONAL">
  !     .TRUE. if this file contains any min, max, or time_average
  !   </IN>
  SUBROUTINE write_axis_meta_data(file_unit, axes, time_ops)
    INTEGER, INTENT(in) :: file_unit, axes(:)
    LOGICAL, INTENT(in), OPTIONAL :: time_ops

    TYPE(domain1d)       :: Domain
    TYPE(domain1d)       :: Edge_Domain

    CHARACTER(len=mxch)  :: axis_name, axis_units
    CHARACTER(len=mxchl) :: axis_long_name
    CHARACTER(len=1)     :: axis_cart_name
    CHARACTER(len=17)    :: calendar_out_name
    INTEGER              :: axis_direction, axis_edges
    REAL, ALLOCATABLE    :: axis_data(:)
    INTEGER, ALLOCATABLE :: axis_extent(:), pelist(:)

    INTEGER              :: calendar, id_axis, id_time_axis
    INTEGER              :: i, index, num, length, edges_index
    INTEGER              :: gbegin, gend, gsize, ndivs
    LOGICAL              :: time_ops1

    IF ( PRESENT(time_ops) ) THEN 
       time_ops1 = time_ops
    ELSE
       time_ops1 = .FALSE.
    END IF

    !---- save the current file_unit ----
    IF ( num_axis_in_file == 0 ) current_file_unit = file_unit

    !---- dummy checks ----
    num = SIZE(axes(:))
    ! <ERROR STATUS="FATAL">number of axes < 1 </ERROR>
    IF ( num < 1 ) CALL error_mesg('write_axis_meta_data', 'number of axes < 1.', FATAL)

    ! <ERROR STATUS="FATAL">writing meta data out-of-order to different files.</ERROR>
    IF ( file_unit /= current_file_unit ) CALL error_mesg('write_axis_meta_data',&
         & 'writing meta data out-of-order to different files.', FATAL)

    !---- check all axes ----
    !---- write axis meta data for new axes ----
    DO i = 1, num
       id_axis = axes(i)
       index = get_axis_index ( id_axis )

       !---- skip axes already written -----
       IF ( index > 0 ) CYCLE

       !---- create new axistype (then point to) -----
       num_axis_in_file = num_axis_in_file + 1
       axis_in_file(num_axis_in_file) = id_axis
       edge_axis_flag(num_axis_in_file) = .FALSE.
       length = get_axis_global_length(id_axis)
       ALLOCATE(axis_data(length))

       CALL get_diag_axis(id_axis, axis_name, axis_units, axis_long_name,&
            & axis_cart_name, axis_direction, axis_edges, Domain, axis_data)

       IF ( Domain .NE. null_domain1d ) THEN
          IF ( length > 0 ) THEN
             CALL mpp_write_meta(file_unit, Axis_types(num_axis_in_file),&
                  & axis_name, axis_units, axis_long_name, axis_cart_name,&
                  & axis_direction, Domain, axis_data )
          ELSE
             CALL mpp_write_meta(file_unit, Axis_types(num_axis_in_file), axis_name,&
                  & axis_units, axis_long_name, axis_cart_name, axis_direction, Domain)
          END IF
       ELSE
          IF ( length > 0 ) THEN
             CALL mpp_write_meta(file_unit, Axis_types(num_axis_in_file), axis_name,&
                  & axis_units, axis_long_name, axis_cart_name, axis_direction, DATA=axis_data)
          ELSE
             CALL mpp_write_meta(file_unit, Axis_types(num_axis_in_file), axis_name,&
                  & axis_units, axis_long_name, axis_cart_name, axis_direction)
          END IF
       END IF

       !---- write additional attribute (calendar_type) for time axis ----
       !---- NOTE: calendar attribute is compliant with CF convention 
       !---- http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-current.htm#cal
       IF ( axis_cart_name == 'T' ) THEN
          time_axis_flag(num_axis_in_file) = .TRUE.
          id_time_axis = mpp_get_id(Axis_types(num_axis_in_file))
          calendar = get_calendar_type()
          
             if(TRIM(valid_calendar_types(calendar)).eq.'THIRTY_DAY_MONTHS') THEN
                calendar_out_name='360_day'
             else
                calendar_out_name=TRIM(valid_calendar_types(calendar))
             endif
          
          CALL mpp_write_meta(file_unit, id_time_axis, 'calendar_type', cval=calendar_out_name)
          CALL mpp_write_meta(file_unit, id_time_axis, 'calendar', cval=calendar_out_name)
          IF ( time_ops1 ) THEN 
             CALL mpp_write_meta( file_unit, id_time_axis, 'bounds', cval = TRIM(axis_name)//'_bounds')        
          END IF
       ELSE
          time_axis_flag(num_axis_in_file) = .FALSE.
       END IF
    
       DEALLOCATE(axis_data)

       !------------- write axis containing edge information ---------------

       !  --- this axis has no edges -----
       IF ( axis_edges <= 0 ) CYCLE

       !  --- was this axis edge previously defined? ---
       id_axis = axis_edges
       edges_index = get_axis_index(id_axis)
       IF ( edges_index > 0 ) CYCLE
    
       !  ---- get data for axis edges ----
       length = get_axis_global_length ( id_axis )
       ALLOCATE(axis_data(length))
       CALL get_diag_axis(id_axis, axis_name, axis_units, axis_long_name, axis_cart_name,&
            & axis_direction, axis_edges, Domain, axis_data )

       !  ---- write edges attribute to original axis ----
       CALL mpp_write_meta(file_unit, mpp_get_id(Axis_types(num_axis_in_file)),&
            & 'edges', cval=axis_name )

       !  ---- add edges index to axis list ----
       !  ---- assume this is not a time axis ----
       num_axis_in_file = num_axis_in_file + 1
       axis_in_file(num_axis_in_file) = id_axis
       edge_axis_flag(num_axis_in_file) = .TRUE.
       time_axis_flag (num_axis_in_file) = .FALSE.

       !  ---- write edges axis to file ----
       IF ( Domain .NE. null_domain1d ) THEN
          ! assume domain decomposition is irregular and loop through all prev and next
          ! domain pointers extracting domain extents.  Assume all pes are used in
          ! decomposition
          CALL mpp_get_global_domain(Domain, begin=gbegin, END=gend, size=gsize)
          CALL mpp_get_layout(Domain, ndivs)
          IF ( ndivs .EQ. 1 ) THEN
             CALL mpp_write_meta(file_unit, Axis_types(num_axis_in_file), axis_name,&
                  & axis_units, axis_long_name, axis_cart_name, axis_direction, DATA=axis_data )
          ELSE
             IF ( ALLOCATED(axis_extent) ) DEALLOCATE(axis_extent)
             ALLOCATE(axis_extent(0:ndivs-1))
             CALL mpp_get_compute_domains(Domain,size=axis_extent(0:ndivs-1))
             gend=gend+1
             axis_extent(ndivs-1)= axis_extent(ndivs-1)+1
             IF ( ALLOCATED(pelist) ) DEALLOCATE(pelist)      
             ALLOCATE(pelist(0:ndivs-1))
             CALL mpp_get_pelist(Domain,pelist)
             CALL mpp_define_domains((/gbegin,gend/),ndivs,Edge_Domain,&
                  & pelist=pelist(0:ndivs-1), extent=axis_extent(0:ndivs-1))
             CALL mpp_write_meta(file_unit, Axis_types(num_axis_in_file),&
                  & axis_name, axis_units, axis_long_name, axis_cart_name,&
                  & axis_direction, Edge_Domain,  DATA=axis_data)
          END IF
       ELSE
          CALL mpp_write_meta(file_unit, Axis_types(num_axis_in_file), axis_name, axis_units,&
               & axis_long_name, axis_cart_name, axis_direction, DATA=axis_data)
       END IF
       DEALLOCATE (axis_data)
    END DO
  END SUBROUTINE write_axis_meta_data
  ! </SUBROUTINE>

  ! <FUNCTION NAME="write_field_meta_data">
  !   <OVERVIEW>
  !     Write the field meta data to file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     TYPE(diag_fieldtype) FUNCTION write_field_meta_data(file_unit, name, axes, units,
  !     long_name, rnage, pack, mval, avg_name, time_method, standard_name, interp_method)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     The meta data for the field is written to the file indicated by file_unit
  !   </DESCRIPTION>
  !   <IN NAME="file_unit" TYPE="INTEGER">Output file unit number</IN>
  !   <IN NAME="name" TYPE="CHARACTER(len=*)">Field name</IN>
  !   <IN NAME="axes" TYPE="INTEGER, DIMENSION(:)">Array of axis IDs</IN>
  !   <IN NAME="units" TYPE="CHARACTER(len=*)">Field units</IN>
  !   <IN NAME="long_name" TYPE="CHARACTER(len=*)">Field's long name</IN>
  !   <IN NAME="range" TYPE="REAL, DIMENSION(2), OPTIONAL">
  !     Valid range (min, max).  If min > max, the range will be ignored
  !   </IN>
  !   <IN NAME="pack" TYPE="INTEGER, OPTIONAL" DEFAULT="2">
  !     Packing flag.  Only valid when range specified.  Valid values:
  !     <UL>
  !       <LI> 1 = 64bit </LI>
  !       <LI> 2 = 32bit </LI>
  !       <LI> 4 = 16bit </LI>
  !       <LI> 8 =  8bit </LI>
  !     </UL>
  !   </IN>
  !   <IN NAME="mval" TYPE="REAL, OPTIONAL">Missing value, must be within valid range</IN>
  !   <IN NAME="avg_name" TYPE="CHARACTER(len=*), OPTIONAL">
  !     Name of variable containing time averaging info
  !   </IN>
  !   <IN NAME="time_method" TYPE="CHARACTER(len=*), OPTIONAL">
  !     Name of transformation applied to the time-varying data, i.e. "avg", "min", "max"
  !   </IN>
  !   <IN NAME="standard_name" TYPE="CHARACTER(len=*), OPTIONAL">Standard name of field</IN>
  !   <IN NAME="interp_method" TYPE="CHARACTER(len=*), OPTIONAL" />
  FUNCTION write_field_meta_data ( file_unit, name, axes, units, long_name, range, pack,&
       & mval, avg_name, time_method,standard_name,interp_method) result ( Field )
    INTEGER, INTENT(in) :: file_unit, axes(:)
    CHARACTER(len=*), INTENT(in) :: name, units, long_name
    REAL, OPTIONAL, INTENT(in) :: RANGE(2), mval
    INTEGER, OPTIONAL, INTENT(in) :: pack
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: avg_name, time_method,standard_name
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: interp_method

    CHARACTER(len=128) :: standard_name2
    TYPE(diag_fieldtype) :: Field
    LOGICAL :: coord_present
    CHARACTER(len=40) :: aux_axes(SIZE(axes))
    CHARACTER(len=160) :: coord_att

    REAL :: scale, add
    INTEGER :: i, indexx, num, ipack, np
    LOGICAL :: use_range
    INTEGER :: axis_indices(SIZE(axes))

    !---- dummy checks ----
    coord_present = .FALSE.
    IF( PRESENT(standard_name) ) THEN 
       standard_name2 = standard_name
    ELSE
       standard_name2 = 'none'
    END IF
    
    num = SIZE(axes(:))
    ! <ERROR STATUS="FATAL">number of axes < 1</ERROR>
    IF ( num < 1 ) CALL error_mesg ( 'write_meta_data', 'number of axes < 1', FATAL)
    ! <ERROR STATUS="FATAL">writing meta data out-of-order to different files</ERROR>
    IF ( file_unit /= current_file_unit ) CALL error_mesg ( 'write_meta_data',  &
         & 'writing meta data out-of-order to different files', FATAL)


    !---- check all axes for this field ----
    !---- set up indexing to axistypes ----
    DO i = 1, num
       indexx = get_axis_index(axes(i))
       !---- point to existing axistype -----
       IF ( indexx > 0 ) THEN
          axis_indices(i) = indexx
       ELSE
          ! <ERROR STATUS="FATAL">axis data not written for field</ERROR>
          CALL error_mesg ('write_field_meta_data',&
               & 'axis data not written for field '//TRIM(name), FATAL)
       END IF
    END DO

    !  Create coordinate attribute
    IF ( num >= 2 ) THEN     
       coord_att = ' '
       DO i = 1, num
          aux_axes(i) = get_axis_aux(axes(i))
          IF( TRIM(aux_axes(i)) /= 'none' ) THEN
             IF(LEN_TRIM(coord_att) == 0) THEN
                coord_att = TRIM(aux_axes(i))
             ELSE
                coord_att = TRIM(coord_att)// ' '//TRIM(aux_axes(i))
             ENDIF
             coord_present = .TRUE.
          END IF
       END DO
    END IF

    !--------------------- write field meta data ---------------------------

    !---- select packing? ----
    !(packing option only valid with range option)
    IF ( PRESENT(pack) ) THEN
       ipack = pack
    ELSE
       ipack = 2
    END IF
    
    !---- check range ----
    use_range = .FALSE.
    add = 0.0
    scale = 1.0
    IF ( PRESENT(range) ) THEN
       IF ( RANGE(2) > RANGE(1) ) THEN
          use_range = .TRUE.
          !---- set packing parameters ----
          IF ( ipack > 2 ) THEN
             np = ipack/4
             add = 0.5*(RANGE(1)+RANGE(2))
             scale = (RANGE(2)-RANGE(1)) / real(max_range(2,np)-max_range(1,np))
          END IF
       END IF
    END IF

    !---- select packing? ----
    IF ( PRESENT(mval) ) THEN
       Field%miss = mval
       Field%miss_present = .TRUE.
       IF ( ipack > 2 ) THEN
          np = ipack/4
          Field%miss_pack = REAL(missval(np))*scale+add
          Field%miss_pack_present = .TRUE.
       ELSE
          Field%miss_pack = mval
          Field%miss_pack_present = .FALSE.
       END IF
    ELSE
       Field%miss_present = .FALSE.
       Field%miss_pack_present = .FALSE.
    END IF

    !------ write meta data and return fieldtype -------
    IF ( use_range ) THEN
       IF ( Field%miss_present ) THEN
          CALL mpp_write_meta(file_unit, Field%Field,&
               & Axis_types(axis_indices(1:num)),&
               & name, units, long_name,&
               & RANGE(1), RANGE(2),&
               & missing=Field%miss_pack,&
               & fill=Field%miss_pack,&
               & scale=scale, add=add, pack=ipack,&
               & time_method=time_method)
       ELSE
          CALL mpp_write_meta(file_unit, Field%Field,&
               & Axis_types(axis_indices(1:num)),&
               & name, units,  long_name,&
               & RANGE(1), RANGE(2),&
               & missing=CMOR_MISSING_VALUE,&
               & fill=CMOR_MISSING_VALUE,&
               & scale=scale, add=add, pack=ipack,&
               & time_method=time_method)
       END IF
    ELSE
       IF ( Field%miss_present ) THEN
          CALL mpp_write_meta(file_unit, Field%Field,&
               & Axis_types(axis_indices(1:num)),&
               & name, units, long_name,&
               & missing=Field%miss_pack,&
               & fill=Field%miss_pack,&
               & pack=ipack, time_method=time_method)
       ELSE
          CALL mpp_write_meta(file_unit, Field%Field,&
               & Axis_types(axis_indices(1:num)),&
               & name, units, long_name,&
               & missing=CMOR_MISSING_VALUE,&
               & fill=CMOR_MISSING_VALUE,&
               & pack=ipack, time_method=time_method)
       END IF
    END IF

    !---- write additional attribute for time averaging -----
    IF ( PRESENT(avg_name) ) THEN
       IF ( avg_name(1:1) /= ' ' ) THEN
          CALL mpp_write_meta(file_unit, mpp_get_id(Field%Field),&
             & 'time_avg_info',&
             & cval=trim(avg_name)//'_T1,'//trim(avg_name)//'_T2,'//trim(avg_name)//'_DT')
       END IF
    END IF

    ! write coordinates attribute for CF compliance
    IF ( coord_present ) &
         CALL mpp_write_meta(file_unit, mpp_get_id(Field%Field),&
         & 'coordinates', cval=TRIM(coord_att))
    IF ( TRIM(standard_name2) /= 'none' ) CALL mpp_write_meta(file_unit, mpp_get_id(Field%Field),&
         & 'standard_name', cval=TRIM(standard_name2))

    !---- write attribute for interp_method ----
    IF( PRESENT(interp_method) ) THEN
       CALL mpp_write_meta ( file_unit, mpp_get_id(Field%Field),&
            & 'interp_method', cval=TRIM(interp_method)) 
    END IF

    !---- get axis domain ----
    Field%Domain = get_domain2d ( axes )
    Field%tile_count = get_tile_count ( axes )

  END FUNCTION write_field_meta_data
  ! </FUNCTION>

  ! <SUBROUTINE NAME="done_meta_data">
  !   <OVERVIEW>
  !     Writes axis data to file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE done_meta_data(file_unit)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Writes axis data to file.  This subroutine is to be called once per file
  !     after all <TT>write_meta_data</TT> calls, and before the first 
  !     <TT>diag_field_out</TT> call.
  !   </DESCRIPTION>
  !   <IN NAME="file_unit" TYPE="INTEGER">Output file unit number</IN>
  SUBROUTINE done_meta_data(file_unit)
    INTEGER,  INTENT(in)  :: file_unit  

    INTEGER               :: i

    !---- write data for all non-time axes ----
    DO i = 1, num_axis_in_file
       IF ( time_axis_flag(i) ) CYCLE
       CALL mpp_write(file_unit, Axis_types(i))
    END DO

    num_axis_in_file = 0
  END SUBROUTINE done_meta_data
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="diag_field_out">
  !   <OVERVIEW>
  !     Writes field data to an output file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE diag_field_out(file_unit, field, data, time)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Writes field data to an output file.
  !   </DESCRIPTION>
  !   <IN NAME="file_unit" TYPE="INTEGER">Output file unit number</IN>
  !   <INOUT NAME="field" TYPE="TYPE(diag_fieldtype)"></INOUT>
  !   <INOUT NAME="data" TYPE="REAL, DIMENSIONS(:,:,:,:)"></INOUT>
  !   <IN NAME="time" TYPE="REAL, OPTIONAL"></IN>
  SUBROUTINE diag_field_out(file_unit, Field, DATA, time)
    INTEGER, INTENT(in) :: file_unit
    TYPE(diag_fieldtype), INTENT(inout) :: Field
    REAL , INTENT(inout) :: data(:,:,:,:)
    REAL, OPTIONAL, INTENT(in) :: time

    !---- replace original missing value with (un)packed missing value ----
    !print *, 'PE,name,miss_pack_present=',mpp_pe(), &
    !  trim(Field%Field%name),Field%miss_pack_present
    IF ( Field%miss_pack_present ) THEN
       WHERE ( DATA == Field%miss ) DATA = Field%miss_pack
    END IF

    !---- output data ----
    IF ( Field%Domain .NE. null_domain2d ) THEN
       IF( Field%miss_present ) THEN
          CALL mpp_write(file_unit, Field%Field, Field%Domain, DATA, time, &
                      tile_count=Field%tile_count, default_data=Field%miss_pack)
       ELSE
          CALL mpp_write(file_unit, Field%Field, Field%Domain, DATA, time, &
                      tile_count=Field%tile_count, default_data=CMOR_MISSING_VALUE)
       END IF
    ELSE
       CALL mpp_write(file_unit, Field%Field, DATA, time)
    END IF
  END SUBROUTINE diag_field_out
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="diag_flush">
  !   <OVERVIEW>
  !     Flush buffer and insure data is not lost.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     CALL diag_flush(file_unit)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     This subroutine can be called periodically to flush the buffer, and
  !     insure that data is not lost if the execution fails.
  !   </DESCRIPTION>
  !   <IN NAME="file_unit" TYPE="INTEGER">Output file unit number to flush</IN>
  SUBROUTINE diag_flush(file_unit)
    INTEGER, INTENT(in) :: file_unit

    CALL mpp_flush (file_unit)
  END SUBROUTINE diag_flush
  ! </SUBROUTINE>


  ! <FUNCTION NAME="get_axis_index">
  !   <OVERVIEW>
  !     Return the axis index number.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION get_axis_index(num)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the axis index number.
  !   </DESCRIPTION>
  !   <IN NAME="num" TYPE="INTEGER"></IN>
  FUNCTION get_axis_index(num) RESULT ( index )
    INTEGER, INTENT(in) :: num

    INTEGER :: index
    INTEGER :: i

    !---- get the array index for this axis type ----
    !---- set up pointers to axistypes ----
    !---- write axis meta data for new axes ----
    index = 0
    DO i = 1, num_axis_in_file
       IF ( num == axis_in_file(i) ) THEN
          index = i
          EXIT
       END IF
    END DO
  END FUNCTION get_axis_index
  ! </FUNCTION>

  ! <SUBROUTINE NAME="get_diag_global_att">
  !   <OVERVIEW>
  !     Return the global attribute type.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     CALL get_diag_global_att(gAtt)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the global attribute type.
  !   </DESCRIPTION>
  !   <OUT NAME="gAtt" TYPE="TYPE(diag_global_att_type"></OUT>
  SUBROUTINE get_diag_global_att(gAtt)
    TYPE(diag_global_att_type), INTENT(out) :: gAtt

    gAtt=diag_global_att
  END SUBROUTINE get_diag_global_att
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="set_diag_global_att">
  !   <OVERVIEW>
  !     Set the global attribute type.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     CALL set_diag_global_att(component, gridType, timeName)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Set the global attribute type.
  !   </DESCRIPTION>
  !   <IN NAME="component" TYPE="CHARACTER(len=*)"></IN>
  !   <IN NAME="gridType" TYPE="CHARACTER(len=*)"></IN>
  !   <IN NAME="tileName" TYPE="CHARACTER(len=*)"></IN>
  SUBROUTINE set_diag_global_att(component, gridType, tileName)
    CHARACTER(len=*),INTENT(in) :: component, gridType, tileName 

    ! The following two lines are set to remove compile time warnings
    ! about 'only used once'.
    CHARACTER(len=64) :: component_tmp
    component_tmp = component
    ! Don't know how to set these for specific component
    ! Want to be able to say 
    ! if(output_file has component) then
    diag_global_att%grid_type = gridType
    diag_global_att%tile_name = tileName
    ! endif
  END SUBROUTINE set_diag_global_att
  ! </SUBROUTINE>

END MODULE diag_output_mod

