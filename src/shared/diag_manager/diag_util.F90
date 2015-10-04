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

MODULE diag_util_mod
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
  !   Seth Underwood
  ! </CONTACT>
  ! <HISTORY SRC="http://cobweb.gfdl.noaa.gov/fms-cgi-bin/viewcvs/FMS/shared/diag_manager/"/>

  ! <OVERVIEW>
  !   Functions and subroutines necessary for the <TT>diag_manager_mod</TT>.
  ! </OVERVIEW>

  ! <DESCRIPTION>
  !   <TT>diag_util_mod</TT> is a set of Fortran functions and subroutines used by the <TT>diag_manager_mod</TT>.
  ! </DESCRIPTION>

  ! <INFO>
  !   <FUTURE>
  !     Make an interface <TT>check_bounds_are_exact</TT> for the subroutines <TT>check_bounds_are_exact_static</TT> and
  !     <TT>check_bounds_are_exact_dynamic</TT>. 
  !     <PRE>
  !       INTERFACE check_bounds_are_exact
  !         MODULE PROCEDURE check_bounds_are_exact_static
  !         MODULE PROCEDURE check_bounds_are_exact_dynamic
  !       END INTERFACE check_bounds_are_exact
  !     </PRE>
  !   </FUTURE>
  ! </INFO>
  USE diag_data_mod, ONLY  : output_fields, input_fields, files, do_diag_field_log, diag_log_unit,&
       & VERY_LARGE_AXIS_LENGTH, time_zero, VERY_LARGE_FILE_FREQ, END_OF_RUN, EVERY_TIME,&
       & DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, DIAG_YEARS, base_time,&
       & time_unit_list, max_files, base_year, base_month, base_day, base_hour, base_minute,&
       & base_second, num_files, max_files, max_fields_per_file, max_out_per_in_field,&
       & max_input_fields,num_input_fields, max_output_fields, num_output_fields, coord_type,&
       & mix_snapshot_average_fields, global_descriptor, CMOR_MISSING_VALUE, use_cmor, pack_size,&
       & debug_diag_manager, conserve_water
  USE diag_axis_mod, ONLY  : get_diag_axis_data, get_axis_global_length, get_diag_axis_cart,&
       & get_domain1d, get_domain2d, diag_subaxes_init, diag_axis_init, get_diag_axis, get_axis_aux,&
       & get_axes_shift, get_diag_axis_name, get_diag_axis_domain_name
  USE diag_output_mod, ONLY: diag_flush, diag_field_out, diag_output_init, write_axis_meta_data,&
       & write_field_meta_data, done_meta_data
  USE diag_grid_mod, ONLY: get_local_indexes
  USE fms_mod, ONLY        : error_mesg, FATAL, WARNING, mpp_pe, mpp_root_pe, lowercase, fms_error_handler
  USE fms_io_mod, ONLY     : get_tile_string, return_domain, string
  USE mpp_domains_mod,ONLY : domain1d, domain2d, mpp_get_compute_domain, null_domain1d, null_domain2d,&
       & OPERATOR(.NE.), OPERATOR(.EQ.), mpp_modify_domain, mpp_get_domain_components,&
       & mpp_get_ntile_count, mpp_get_current_ntile, mpp_get_tile_id, mpp_mosaic_defined, mpp_get_tile_npes
  USE time_manager_mod,ONLY: time_type, OPERATOR(==), OPERATOR(>), NO_CALENDAR, increment_date,&
       & increment_time, get_calendar_type, get_date, get_time, leap_year, OPERATOR(-),&
       & OPERATOR(<), OPERATOR(>=), OPERATOR(<=)
  USE mpp_io_mod, ONLY : mpp_close
  USE mpp_mod, ONLY : mpp_npes
  USE constants_mod, ONLY : SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_MINUTE

  IMPLICIT NONE
  PRIVATE
  PUBLIC get_subfield_size, log_diag_field_info, update_bounds, check_out_of_bounds,&
       & check_bounds_are_exact_dynamic, check_bounds_are_exact_static, init_file, diag_time_inc,&
       & find_input_field, init_input_field, init_output_field, diag_data_out, write_static,&
       & check_duplicate_output_fields, get_date_dif, get_subfield_vert_size, sync_file_times

  CHARACTER(len=128),PRIVATE  :: version =&
       & '$Id: diag_util.F90,v 19.0.2.2 2012/04/03 18:41:44 sdu Exp $'
  CHARACTER(len=128),PRIVATE  :: tagname =&
       & '$Name: siena_201211 $'

CONTAINS

  ! <SUBROUTINE NAME="get_subfield_size">
  !   <OVERVIEW>
  !     Get the size, start, and end indices for output fields.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_subfield_size(axes, outnum)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Get the size, start and end indices for <TT>output_fields(outnum)</TT>, then  
  !     fill in <TT>output_fields(outnum)%output_grid%(start_indx, end_indx)</TT>
  !   </DESCRIPTION>
  !   <IN NAME="axes" TYPE="INTEGER, DIMENSION(:)">Axes of the <TT>input_field</TT>.</IN>
  !   <IN NAME="outnum" TYPE="INTEGER">Position in array <TT>output_fields</TT>.</IN>
  SUBROUTINE get_subfield_size(axes, outnum)
    INTEGER, INTENT(in) :: axes(:) ! axes of the input_field
    INTEGER, INTENT(in) :: outnum  ! position in array output_fields

    REAL, ALLOCATABLE   :: global_lat(:), global_lon(:), global_depth(:)
    INTEGER :: global_axis_size
    INTEGER :: i,xbegin,xend,ybegin,yend,xbegin_l,xend_l,ybegin_l,yend_l 
    CHARACTER(len=1) :: cart
    TYPE(domain2d) :: Domain2, Domain2_new
    TYPE(domain1d) :: Domain1, Domain1x, Domain1y
    REAL :: start(3), end(3) ! start and end coordinates in 3 axes
    INTEGER :: gstart_indx(3), gend_indx(3) ! global start and end indices of output domain in 3 axes 
    REAL, ALLOCATABLE :: subaxis_x(:), subaxis_y(:), subaxis_z(:) !containing local coordinates in x,y,z axes
    CHARACTER(len=128) :: msg
    INTEGER :: ishift, jshift
    CHARACTER(len=128), DIMENSION(2) :: axis_domain_name

    !initilization for local output
    ! initially out of (lat/lon/depth) range
    start = -1.e10
    end = -1.e10 
    gstart_indx = -1
    gend_indx=-1

    ! get axis data (lat, lon, depth) and indices
    start = output_fields(outnum)%output_grid%start
    end = output_fields(outnum)%output_grid%end

    CALL get_diag_axis_domain_name(axes(1), axis_domain_name(1))
    CALL get_diag_axis_domain_name(axes(2), axis_domain_name(2))

    IF (   INDEX(lowercase(axis_domain_name(1)), 'cubed') == 0 .AND. &
         & INDEX(lowercase(axis_domain_name(2)), 'cubed') == 0 ) THEN
       DO i = 1, SIZE(axes(:))
          global_axis_size = get_axis_global_length(axes(i))
          output_fields(outnum)%output_grid%subaxes(i) = -1
          CALL get_diag_axis_cart(axes(i), cart)
          SELECT CASE(cart)
          CASE ('X')
             ! <ERROR STATUS="FATAL">wrong order of axes.  X should come first.</ERROR>
             IF( i.NE.1 ) CALL error_mesg('diag_util_mod::get_subfield_size',&
                  & 'wrong order of axes, X should come first',FATAL)
             ALLOCATE(global_lon(global_axis_size))
             CALL get_diag_axis_data(axes(i),global_lon)
             IF( INT( start(i)*END(i) ) == 1 ) THEN 
                gstart_indx(i) = 1
                gend_indx(i) = global_axis_size
                output_fields(outnum)%output_grid%subaxes(i) = axes(i)
             ELSE 
                gstart_indx(i) = get_index(start(i),global_lon)
                gend_indx(i) = get_index(END(i),global_lon)
             END IF
             ALLOCATE(subaxis_x(gstart_indx(i):gend_indx(i)))
             subaxis_x=global_lon(gstart_indx(i):gend_indx(i))   
          CASE ('Y')
             ! <ERROR STATUS="FATAL">wrong order of axes, Y should come second.</ERROR>
             IF( i.NE.2 ) CALL error_mesg('diag_util_mod::get_subfield_size',&
                  & 'wrong order of axes, Y should come second',FATAL)
             ALLOCATE(global_lat(global_axis_size))
             CALL get_diag_axis_data(axes(i),global_lat)
             IF( INT( start(i)*END(i) ) == 1 ) THEN 
                gstart_indx(i) = 1
                gend_indx(i) = global_axis_size
                output_fields(outnum)%output_grid%subaxes(i) = axes(i)
             ELSE
                gstart_indx(i) = get_index(start(i),global_lat)
                gend_indx(i) = get_index(END(i),global_lat)
             END IF
             ALLOCATE(subaxis_y(gstart_indx(i):gend_indx(i)))
             subaxis_y=global_lat(gstart_indx(i):gend_indx(i))
          CASE ('Z')
             ! <ERROR STATUS="FATAL">wrong values in vertical axis of region</ERROR>
             IF ( start(i)*END(i)<0 ) CALL error_mesg('diag_util_mod::get_subfield_size',&
                  & 'wrong values in vertical axis of region',FATAL)
             IF ( start(i)>=0 .AND. END(i)>0 ) THEN 
                ALLOCATE(global_depth(global_axis_size))
                CALL get_diag_axis_data(axes(i),global_depth)
                gstart_indx(i) = get_index(start(i),global_depth)
                gend_indx(i) = get_index(END(i),global_depth)
                ALLOCATE(subaxis_z(gstart_indx(i):gend_indx(i)))
                subaxis_z=global_depth(gstart_indx(i):gend_indx(i))
                output_fields(outnum)%output_grid%subaxes(i) =&
                     & diag_subaxes_init(axes(i),subaxis_z, gstart_indx(i),gend_indx(i))
                DEALLOCATE(subaxis_z,global_depth)
             ELSE ! regional vertical axis is the same as global vertical axis
                gstart_indx(i) = 1
                gend_indx(i) = global_axis_size
                output_fields(outnum)%output_grid%subaxes(i) = axes(i)
                ! <ERROR STATUS="FATAL">i should equal 3 for z axis</ERROR>
                IF( i /= 3 ) CALL error_mesg('diag_util_mod::get_subfield_size',&
                     & 'i should equal 3 for z axis', FATAL)
             END IF
          CASE default
             ! <ERROR STATUS="FATAL">Wrong axis_cart</ERROR>
             CALL error_mesg('diag_util_mod::get_subfield_size', 'Wrong axis_cart', FATAL)
          END SELECT
       END DO

       DO i = 1, SIZE(axes(:))
          IF( gstart_indx(i) == -1 .OR. gend_indx(i) == -1 ) THEN
             ! <ERROR STATUS="FATAL">
             !   can not find gstart_indx/gend_indx for <output_fields(outnum)%output_name>,
             !   check region bounds for axis <i>.
             ! </ERROR>
             WRITE(msg,'(A,I2)') ' check region bounds for axis ', i
             CALL error_mesg('diag_util_mod::get_subfield_size', 'can not find gstart_indx/gend_indx for '&
                  & //TRIM(output_fields(outnum)%output_name)//','//TRIM(msg), FATAL)
          END IF
       END DO
    ELSE ! cubed sphere
       ! get the i and j start and end indexes
       CALL get_local_indexes(LONSTART=start(1), LONEND=END(1), &
            &                 LATSTART=start(2), LATEND=END(2), &
            &                 ISTART=gstart_indx(1), IEND=gend_indx(1), &
            &                 JSTART=gstart_indx(2), JEND=gend_indx(2))
       global_axis_size =  get_axis_global_length(axes(1))
       ALLOCATE(global_lon(global_axis_size))
       global_axis_size = get_axis_global_length(axes(2))
       ALLOCATE(global_lat(global_axis_size))
       CALL get_diag_axis_data(axes(1),global_lon)
       CALL get_diag_axis_data(axes(2),global_lat)
       IF (   (gstart_indx(1) > 0 .AND. gstart_indx(2) > 0) .AND. &
            & (gend_indx(1) > 0 .AND. gend_indx(2) > 0) ) THEN
          ALLOCATE(subaxis_x(gstart_indx(1):gend_indx(1)))
          ALLOCATE(subaxis_y(gstart_indx(2):gend_indx(2)))
          subaxis_x=global_lon(gstart_indx(1):gend_indx(1))
          subaxis_y=global_lat(gstart_indx(2):gend_indx(2))
       END IF

       ! Now deal with the Z component
       IF ( SIZE(axes(:)) > 2 ) THEN
          global_axis_size = get_axis_global_length(axes(3))
          output_fields(outnum)%output_grid%subaxes(3) = -1
          CALL get_diag_axis_cart(axes(3), cart)
          ! <ERROR STATUS="FATAL">
          !   axis(3) should be Z-axis
          ! </ERROR>
          IF ( lowercase(cart) /= 'z' ) CALL error_mesg('diag_util_mod::get_subfield_size', &
               &'axis(3) should be Z-axis', FATAL)
          ! <ERROR STATUS="FATAL">
          !   wrong values in vertical axis of region
          ! </ERROR>
          IF ( start(3)*END(3)<0 ) CALL error_mesg('diag_util_mod::get_subfield_size',&
               & 'wrong values in vertical axis of region',FATAL)
          IF ( start(3)>=0 .AND. END(3)>0 ) THEN 
             ALLOCATE(global_depth(global_axis_size))
             CALL get_diag_axis_data(axes(3),global_depth)
             gstart_indx(3) = get_index(start(3),global_depth)
             IF( start(3) == 0.0 )  gstart_indx(3) = 1
             gend_indx(3) = get_index(END(3),global_depth)
             IF( start(3) >= MAXVAL(global_depth) ) gstart_indx(3)= global_axis_size
             IF( END(3)   >= MAXVAL(global_depth) ) gend_indx(3)  = global_axis_size
             
             ALLOCATE(subaxis_z(gstart_indx(3):gend_indx(3)))
             subaxis_z=global_depth(gstart_indx(3):gend_indx(3))
             output_fields(outnum)%output_grid%subaxes(3) =&
                  & diag_subaxes_init(axes(3),subaxis_z, gstart_indx(3),gend_indx(3))
             DEALLOCATE(subaxis_z,global_depth)
          ELSE ! regional vertical axis is the same as global vertical axis
             gstart_indx(3) = 1
             gend_indx(3) = global_axis_size
             output_fields(outnum)%output_grid%subaxes(3) = axes(3)
          END IF
       END IF
    END IF
    
    ! get domain and compute_domain(xbegin,xend,ybegin,yend)
    xbegin=-1
    xend=-1
    ybegin=-1
    yend=-1

    Domain2 = get_domain2d(axes)
    IF ( Domain2 .NE. NULL_DOMAIN2D ) THEN
       CALL mpp_get_compute_domain(Domain2,xbegin,xend,ybegin,yend)
       CALL mpp_get_domain_components(Domain2, Domain1x, Domain1y)
    ELSE
       DO i = 1, MIN(SIZE(axes(:)),2)    
          Domain1 = get_domain1d(axes(i))
          IF ( Domain1 .NE. NULL_DOMAIN1D ) THEN
             CALL get_diag_axis_cart(axes(i),cart)
             SELECT CASE(cart)
             CASE ('X')
                Domain1x = get_domain1d(axes(i))
                CALL mpp_get_compute_domain(Domain1x, xbegin, xend)
             CASE ('Y')
                Domain1y = get_domain1d(axes(i))
                CALL mpp_get_compute_domain(Domain1y, ybegin, yend)
             CASE default ! do nothing here
             END SELECT
          ELSE
             ! <ERROR STATUS="FATAL">No domain available</ERROR>
             CALL error_mesg('diag_util_mod::get_subfield_size', 'NO domain available', FATAL)
          END IF
       END DO
    END IF

    CALL get_axes_shift(axes, ishift, jshift)
    xend = xend+ishift
    yend = yend+jshift

    IF ( xbegin== -1 .OR. xend==-1 .OR. ybegin==-1 .OR. yend==-1 ) THEN
       ! <ERROR STATUS="FATAL">wrong compute domain indices</ERROR>
       CALL error_mesg('diag_util_mod::get_subfield_size', 'wrong compute domain indices',FATAL)  
    END IF
      
    ! get the area containing BOTH compute domain AND local output area
    IF(gstart_indx(1)> xend .OR. xbegin > gend_indx(1)) THEN
       output_fields(outnum)%output_grid%l_start_indx(1) = -1
       output_fields(outnum)%output_grid%l_end_indx(1) = -1
       output_fields(outnum)%need_compute = .FALSE. ! not involved
    ELSEIF (gstart_indx(2)> yend .OR. ybegin > gend_indx(2)) THEN
       output_fields(outnum)%output_grid%l_start_indx(2) = -1
       output_fields(outnum)%output_grid%l_end_indx(2) = -1
       output_fields(outnum)%need_compute = .FALSE. ! not involved
    ELSE
       output_fields(outnum)%output_grid%l_start_indx(1) = MAX(xbegin, gstart_indx(1))
       output_fields(outnum)%output_grid%l_start_indx(2) = MAX(ybegin, gstart_indx(2))
       output_fields(outnum)%output_grid%l_end_indx(1) = MIN(xend, gend_indx(1))
       output_fields(outnum)%output_grid%l_end_indx(2) = MIN(yend, gend_indx(2))
       output_fields(outnum)%need_compute = .TRUE.  ! involved in local output
    END IF

    IF ( output_fields(outnum)%need_compute ) THEN
       ! need to modify domain1d and domain2d for subaxes
       xbegin_l = output_fields(outnum)%output_grid%l_start_indx(1)
       xend_l = output_fields(outnum)%output_grid%l_end_indx(1)
       ybegin_l = output_fields(outnum)%output_grid%l_start_indx(2)
       yend_l = output_fields(outnum)%output_grid%l_end_indx(2)
       CALL mpp_modify_domain(Domain2, Domain2_new, xbegin_l,xend_l, ybegin_l,yend_l,&
            & gstart_indx(1),gend_indx(1), gstart_indx(2),gend_indx(2))

       output_fields(outnum)%output_grid%subaxes(1) =&
            & diag_subaxes_init(axes(1),subaxis_x, gstart_indx(1),gend_indx(1),Domain2_new)
       output_fields(outnum)%output_grid%subaxes(2) =&
            & diag_subaxes_init(axes(2),subaxis_y, gstart_indx(2),gend_indx(2),Domain2_new)
       DO i = 1, SIZE(axes(:))
          IF(output_fields(outnum)%output_grid%subaxes(i) == -1) THEN  
             ! <ERROR STATUS="FATAL">
             !   <output_fields(outnum)%output_name> error at i = <i>
             ! </ERROR>
             WRITE(msg,'(a,"/",I4)') 'at i = ',i
             CALL error_mesg('diag_util_mod::get_subfield_size '//TRIM(output_fields(outnum)%output_name),&
                  'error '//TRIM(msg), FATAL)   
          END IF
       END DO

       ! local start index should start from 1
       output_fields(outnum)%output_grid%l_start_indx(1) = MAX(xbegin, gstart_indx(1)) - xbegin + 1   
       output_fields(outnum)%output_grid%l_start_indx(2) = MAX(ybegin, gstart_indx(2)) - ybegin + 1
       output_fields(outnum)%output_grid%l_end_indx(1) = MIN(xend, gend_indx(1)) - xbegin + 1 
       output_fields(outnum)%output_grid%l_end_indx(2) = MIN(yend, gend_indx(2)) - ybegin + 1
       IF ( SIZE(axes(:))>2 ) THEN
          output_fields(outnum)%output_grid%l_start_indx(3) = gstart_indx(3)
          output_fields(outnum)%output_grid%l_end_indx(3) = gend_indx(3)
       ELSE
          output_fields(outnum)%output_grid%l_start_indx(3) = 1
          output_fields(outnum)%output_grid%l_end_indx(3) = 1
       END IF
    END IF
    IF ( ALLOCATED(subaxis_x) ) DEALLOCATE(subaxis_x, global_lon)
    IF ( ALLOCATED(subaxis_y) ) DEALLOCATE(subaxis_y, global_lat)

  END SUBROUTINE get_subfield_size
  ! </SUBROUTINE>
  
  ! <SUBROUTINE NAME="get_subfield_vert_size">
  !   <OVERVIEW>
  !     Get size, start and end indices for output fields.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_subfield_vert_size(axes, outnum)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Get size, start and end indices for <TT>output_fields(outnum)</TT>, fill in
  !     <TT>output_fields(outnum)%output_grid%(start_indx, end_indx)</TT>.
  !   </DESCRIPTION>
  !   <IN NAME="axes" TYPE="INTEGER, DIMENSION(:)">Axes of the <TT>input_field</TT></IN>
  !   <IN NAME="outnum" TYPE="INTEGER">Position in array <TT>output_fields</TT>.</IN>
  SUBROUTINE get_subfield_vert_size(axes, outnum)
    INTEGER, DIMENSION(:), INTENT(in) :: axes ! axes of the input_field
    INTEGER, INTENT(in) :: outnum  ! position in array output_fields

    REAL, DIMENSION(3) :: start, end ! start and end coordinates in 3 axes
    REAL, ALLOCATABLE, DIMENSION(:) :: global_depth
    REAL, ALLOCATABLE, DIMENSION(:) :: subaxis_z !containing local coordinates in x,y,z axes
    INTEGER :: i, global_axis_size
    INTEGER, DIMENSION(3) :: gstart_indx, gend_indx ! global start and end indices of output domain in 3 axes 
    CHARACTER(len=1) :: cart
    CHARACTER(len=128) :: msg

    !initilization for local output
    start = -1.e10
    end = -1.e10 ! initially out of (lat/lon/depth) range
    gstart_indx = -1 
    gend_indx=-1

    ! get axis data (lat, lon, depth) and indices
    start= output_fields(outnum)%output_grid%start
    end = output_fields(outnum)%output_grid%end

    DO i = 1, SIZE(axes(:))   
       global_axis_size = get_axis_global_length(axes(i))
       output_fields(outnum)%output_grid%subaxes(i) = -1
       CALL get_diag_axis_cart(axes(i), cart)
       SELECT CASE(cart)
       CASE ('X')
          ! <ERROR STATUS="FATAL">wrong order of axes, X should come first</ERROR>
          IF ( i.NE.1 ) CALL error_mesg('diag_util_mod::get_subfield_vert_size',&
               & 'wrong order of axes, X should come first',FATAL)
          gstart_indx(i) = 1
          gend_indx(i) = global_axis_size
          output_fields(outnum)%output_grid%subaxes(i) = axes(i)
       CASE ('Y')
          ! <ERROR STATUS="FATAL">wrong order of axes, Y should come second</ERROR>
          IF( i.NE.2 ) CALL error_mesg('diag_util_mod::get_subfield_vert_size',&
               & 'wrong order of axes, Y should come second',FATAL)
          gstart_indx(i) = 1
          gend_indx(i) = global_axis_size
          output_fields(outnum)%output_grid%subaxes(i) = axes(i)
       CASE ('Z')
          ! <ERROR STATUS="FATAL">wrong values in vertical axis of region</ERROR>
          IF( start(i)*END(i) < 0 ) CALL error_mesg('diag_util_mod::get_subfield_vert_size',&
               & 'wrong values in vertical axis of region',FATAL)
          IF( start(i) >= 0 .AND. END(i) > 0 ) THEN 
             ALLOCATE(global_depth(global_axis_size))
             CALL get_diag_axis_data(axes(i),global_depth)
             gstart_indx(i) = get_index(start(i),global_depth)
             IF( start(i) == 0.0 )  gstart_indx(i) = 1

             gend_indx(i) = get_index(END(i),global_depth)
             IF( start(i) >= MAXVAL(global_depth) ) gstart_indx(i)= global_axis_size
             IF( END(i)   >= MAXVAL(global_depth) ) gend_indx(i)  = global_axis_size

             ALLOCATE(subaxis_z(gstart_indx(i):gend_indx(i)))
             subaxis_z=global_depth(gstart_indx(i):gend_indx(i))
             output_fields(outnum)%output_grid%subaxes(i) =&
                  & diag_subaxes_init(axes(i),subaxis_z, gstart_indx(i),gend_indx(i))
             DEALLOCATE(subaxis_z,global_depth)
          ELSE !   vertical axis is the same as global vertical axis
             gstart_indx(i) = 1
             gend_indx(i) = global_axis_size
             output_fields(outnum)%output_grid%subaxes(i) = axes(i)
             ! <ERROR STATUS="FATAL">i should equal 3 for z axis</ERROR>
             IF( i /= 3 ) CALL error_mesg('diag_util_mod::get_subfield_vert_size',&
                  & 'i should equal 3 for z axis', FATAL)
          END IF
       CASE default
          ! <ERROR STATUS="FATAL">Wrong axis_cart</ERROR>
          CALL error_mesg('diag_util_mod::get_subfield_vert_size', 'Wrong axis_cart', FATAL)
       END SELECT
    END DO

    DO i = 1,SIZE(axes(:))
       IF ( gstart_indx(i)== -1 .OR. gend_indx(i)== -1 ) THEN
          ! <ERROR STATUS="FATAL">
          !   can not find gstart_indx/gend_indx for <output_fields(outnum)%output_name>
          !   check region bounds for axis
          ! </ERROR>
          WRITE(msg,'(A,I2)') ' check region bounds for axis ', i
          CALL error_mesg('diag_util_mod::get_subfield_vert_size', 'can not find gstart_indx/gend_indx for '&
               & //TRIM(output_fields(outnum)%output_name)//','//TRIM(msg), FATAL)
       END IF
    END DO

    DO i= 1, 2
       output_fields(outnum)%output_grid%l_start_indx(i) = gstart_indx(i)
       output_fields(outnum)%output_grid%l_end_indx(i)   = gend_indx(i)
    END DO

    IF( SIZE(axes(:)) > 2 ) THEN
       output_fields(outnum)%output_grid%l_start_indx(3) = gstart_indx(3)
       output_fields(outnum)%output_grid%l_end_indx(3)   = gend_indx(3)
    ELSE
       output_fields(outnum)%output_grid%l_start_indx(3) = 1
       output_fields(outnum)%output_grid%l_end_indx(3)   = 1
    END IF
  END SUBROUTINE get_subfield_vert_size
  ! </SUBROUTINE>
  
  ! <PRIVATE>
  ! <FUNCTION NAME="get_index">
  !   <OVERVIEW>
  !     Find index <TT>i</TT> of array such that <TT>array(i)</TT> is closest to number.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION get_index(number, array)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Find index <TT>i</TT> of array such that <TT>array(i)</TT> is closest to number.
  !     Array must be  monotonouslly ordered.
  !   </DESCRIPTION>
  !   <IN NAME="number" TYPE="REAL"></IN>
  !   <IN NAME="array" TYPE="REAL, DIMENSION(:)"></IN>
  INTEGER FUNCTION get_index(number, array)
    REAL, INTENT(in) :: number
    REAL, INTENT(in), DIMENSION(:) :: array

    INTEGER :: i, n
    LOGICAL :: found

    n = SIZE(array(:))
    ! check if array is monotonous
    DO i = 2, n-1
       IF( (array(i-1)<array(i).AND.array(i)>array(i+1)) .OR. (array(i-1)>array(i).AND.array(i)<array(i+1))) THEN
          ! <ERROR STATUS="FATAL">array NOT monotonously ordered</ERROR>
          CALL error_mesg('diag_util_mod::get_index', 'array NOT monotonously ordered',FATAL) 
       END IF
    END DO
    get_index = -1
    found = .FALSE.
    ! search in increasing array 
    DO i = 1, n-1                
       IF ( (array(i)<=number).AND.(array(i+1)>= number) ) THEN
          IF( number - array(i) <= array(i+1) - number ) THEN
             get_index = i
             found=.TRUE.
          ELSE
             get_index = i+1
             found=.TRUE.
          ENDIF
          EXIT
       END IF
    END DO
    ! if not found, search in decreasing array
    IF( .NOT.found ) THEN
       DO i = 1, n-1
          IF ( (array(i)>=number).AND.(array(i+1)<= number) ) THEN
             IF ( array(i)-number <= number-array(i+1) ) THEN
                get_index = i 
                found = .TRUE.
             ELSE
                get_index = i+1
                found = .TRUE.
             END IF
             EXIT
          END IF
       END DO
    END IF
    ! if still not found, is it less than the first element
    ! or greater than last element? (Increasing Array)
    IF ( .NOT. found ) THEN
       IF ( array(1).GT.number ) THEN
          get_index = 1
          found = .TRUE.
       ELSE IF ( array(n).LT.number ) THEN
          get_index = n
          found = .TRUE.
       ELSE
          found = .FALSE.
       END IF
    END IF
   
   ! if still not found, is it greater than the first element
   ! or less than the last element? (Decreasing Array)
    IF ( .NOT. found ) THEN
       IF ( array(1).LT.number ) THEN
          get_index = 1
          found = .TRUE.
       ELSE IF ( array(n).GT.number ) THEN
          get_index = n
          found = .TRUE.
       ELSE
          found = .FALSE.
       END IF
    END IF
  END FUNCTION get_index
  ! </FUNCTION>
  ! </PRIVATE>

  ! <SUBROUTINE NAME="log_diag_field_info">
  !   <OVERVIEW>
  !     Writes brief diagnostic field info to the log file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE log_diag_field_info(module_name, field_name, axes, long_name, units,
  !     missing_value, range, dynamic)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     If the <TT>do_diag_field_log</TT> namelist parameter is .TRUE.,
  !     then a line briefly describing diagnostic field is added to
  !     the log file.  Normally users should not call this subroutine
  !     directly, since it is called by register_static_field and
  !     register_diag_field if do_not_log is not set to .TRUE..  It is
  !     used, however, in LM3 to avoid excessive logs due to the
  !     number of fields registered for each of the tile types.  LM3
  !     code uses a do_not_log parameter in the registration calls,
  !     and subsequently calls this subroutine to log field information
  !     under a generic name.
  !   </DESCRIPTION>
  !   <IN NAME="module_name" TYPE="CHARACTER(len=*)">Module name.</IN>
  !   <IN NAME="field_name" TYPE="CHARACTER(len=*)">Field name.</IN>
  !   <IN NAME="axes" TYPE="INTEGER, DIMENSION(:)">Axis IDs.</IN>
  !   <IN NAME="long_name" TYPE="CHARACTER(len=*), OPTIONAL">Long name for field.</IN>
  !   <IN NAME="units" TYPE="CHARACTER(len=*), OPTIONAL">Unit of field.</IN>
  !   <IN NAME="missing_value" TYPE="REAL, OPTIONAL">Missing value value.</IN>
  !   <IN NAME="range" TYPE="REAL, DIMENSION(2), OPTIONAL">Valid range of values for field.</IN>
  !   <IN NAME="dynamic" TYPE="LOGICAL, OPTIONAL"><TT>.TRUE.</TT> if field is not static.</IN>
  SUBROUTINE log_diag_field_info(module_name, field_name, axes, long_name, units,&
       & missing_value, range, dynamic)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name
    INTEGER, DIMENSION(:), INTENT(in) :: axes
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name, units
    REAL, OPTIONAL, INTENT(in) :: missing_value
    REAL, DIMENSION(2), OPTIONAL, INTENT(IN) :: range
    LOGICAL, OPTIONAL, INTENT(in) :: dynamic

    ! ---- local vars
    CHARACTER(len=256) :: lmodule, lfield, lname, lunits
    CHARACTER(len=64)  :: lmissval, lmin, lmax
    CHARACTER(len=8)   :: numaxis, timeaxis
    CHARACTER(len=1)   :: sep = '|'
    CHARACTER(len=256) :: axis_name, axes_list
    INTEGER :: i

    IF ( .NOT.do_diag_field_log ) RETURN
    IF ( mpp_pe().NE.mpp_root_pe() ) RETURN

    lmodule = TRIM(module_name)
    lfield = TRIM(field_name)

    IF ( PRESENT(long_name) ) THEN
       lname  = TRIM(long_name)
    ELSE 
       lname  = ''
    END IF
    
    IF ( PRESENT(units) ) THEN
       lunits = TRIM(units)
    ELSE
       lunits = ''
    END IF
 
    WRITE (numaxis,'(i1)') SIZE(axes)

    IF (PRESENT(missing_value)) THEN
       IF ( use_cmor ) THEN
          WRITE (lmissval,*) CMOR_MISSING_VALUE
       ELSE
          WRITE (lmissval,*) missing_value
       END IF
    ELSE
       lmissval = ''
    ENDIF

    IF ( PRESENT(range) ) THEN
       WRITE (lmin,*) range(1)
       WRITE (lmax,*) range(2)
    ELSE
       lmin = ''
       lmax = ''
    END IF

    IF ( PRESENT(dynamic) ) THEN
       IF (dynamic) THEN
          timeaxis = 'T'
       ELSE
          timeaxis = 'F'
       END IF
    ELSE
       timeaxis = ''
    END IF

    axes_list=''
    DO i = 1, SIZE(axes)
       CALL get_diag_axis_name(axes(i),axis_name)
       IF ( TRIM(axes_list) /= '' ) axes_list = TRIM(axes_list)//','
       axes_list = TRIM(axes_list)//TRIM(axis_name)
    END DO

    !write (diag_log_unit,'(8(a,a),a)') &
    WRITE (diag_log_unit,'(777a)') &
         & TRIM(lmodule),  sep, TRIM(lfield),  sep, TRIM(lname),    sep,&
         & TRIM(lunits),   sep, TRIM(numaxis), sep, TRIM(timeaxis), sep,&
         & TRIM(lmissval), sep, TRIM(lmin),    sep, TRIM(lmax),     sep,&
         & TRIM(axes_list)
  END SUBROUTINE log_diag_field_info
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="update_bounds">
  !   <OVERVIEW>
  !     Update the <TT>output_fields</TT> min and max boundaries.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE update_bounds(out_num, lower_i, upper_i, lower_j, upper_j, lower_k, upper_k)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Update the <TT>output_fields</TT> x, y, and z min and max boundaries (array indices).
  !   </DESCRIPTION>
  !   <IN NAME="out_num" TYPE="INTEGER"><TT>output_field</TT> ID.</IN>
  !   <IN NAME="lower_i" TYPE="INTEGER">Lower <TT>i</TT> bound.</IN>
  !   <IN NAME="upper_i" TYPE="INTEGER">Upper <TT>i</TT> bound.</IN>
  !   <IN NAME="lower_j" TYPE="INTEGER">Lower <TT>j</TT> bound.</IN>
  !   <IN NAME="upper_j" TYPE="INTEGER">Upper <TT>j</TT> bound.</IN>
  !   <IN NAME="lower_k" TYPE="INTEGER">Lower <TT>k</TT> bound.</IN>
  !   <IN NAME="upper_k" TYPE="INTEGER">Upper <TT>k</TT> bound.</IN>
  SUBROUTINE update_bounds(out_num, lower_i, upper_i, lower_j, upper_j, lower_k, upper_k)
    INTEGER, INTENT(in) :: out_num, lower_i, upper_i, lower_j, upper_j, lower_k, upper_k
    
    output_fields(out_num)%imin = MIN(output_fields(out_num)%imin, lower_i)
    output_fields(out_num)%imax = MAX(output_fields(out_num)%imax, upper_i)
    output_fields(out_num)%jmin = MIN(output_fields(out_num)%jmin, lower_j)
    output_fields(out_num)%jmax = MAX(output_fields(out_num)%jmax, upper_j)
    output_fields(out_num)%kmin = MIN(output_fields(out_num)%kmin, lower_k)
    output_fields(out_num)%kmax = MAX(output_fields(out_num)%kmax, upper_k)
  END SUBROUTINE update_bounds
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="check_out_of_bounds">
  !   <OVERVIEW>
  !     Checks if the array indices for <TT>output_fields(out_num)</TT> are outside the <TT>output_fields(out_num)%buffer</TT> upper
  !     and lower bounds.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE check_out_of_bounds(out_num, diag_field_id, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>check_out_of_bounds</TT> verifies the array min and max indices in the x, y, and z directions of <TT>
  !     output_fields(out_num)</TT> are not outside the upper and lower array boundaries of
  !     <TT>output_fields(out_num)%buffer</TT>.  If the min and max indices are outside the upper and lower bounds of the buffer
  !     array, then <TT>check_out_of_bounds</TT> returns an error string.
  !   </DESCRIPTION>
  !   <IN NAME="out_num" TYPE="INTEGER">
  !     Output field ID number.
  !   </IN>
  !   <IN NAME="diag_field_id" TYPE="INTEGER">
  !     Input field ID number.
  !   </IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*)">
  !     Return status of <TT>check_out_of_bounds</TT>.  An empty error string indicates the x, y, and z indices are not outside the
  !     buffer array boundaries.
  !   </OUT>
  SUBROUTINE check_out_of_bounds(out_num, diag_field_id, err_msg)
    INTEGER, INTENT(in) :: out_num, diag_field_id
    CHARACTER(len=*), INTENT(out) :: err_msg

    CHARACTER(len=128) :: error_string1, error_string2

    IF (   output_fields(out_num)%imin < LBOUND(output_fields(out_num)%buffer,1) .OR.&
         & output_fields(out_num)%imax > UBOUND(output_fields(out_num)%buffer,1) .OR.&
         & output_fields(out_num)%jmin < LBOUND(output_fields(out_num)%buffer,2) .OR.&
         & output_fields(out_num)%jmax > UBOUND(output_fields(out_num)%buffer,2) .OR.&
         & output_fields(out_num)%kmin < LBOUND(output_fields(out_num)%buffer,3) .OR.&
         & output_fields(out_num)%kmax > UBOUND(output_fields(out_num)%buffer,3) ) THEN
       WRITE(error_string1,'(a,"/",a)') TRIM(input_fields(diag_field_id)%module_name),&
            & TRIM(output_fields(out_num)%output_name)
       error_string2 ='Buffer bounds=   :   ,   :   ,   :     Actual bounds=   :   ,   :   ,   :   '
       WRITE(error_string2(15:17),'(i3)') LBOUND(output_fields(out_num)%buffer,1)
       WRITE(error_string2(19:21),'(i3)') UBOUND(output_fields(out_num)%buffer,1)
       WRITE(error_string2(23:25),'(i3)') LBOUND(output_fields(out_num)%buffer,2)
       WRITE(error_string2(27:29),'(i3)') UBOUND(output_fields(out_num)%buffer,2)
       WRITE(error_string2(31:33),'(i3)') LBOUND(output_fields(out_num)%buffer,3)
       WRITE(error_string2(35:37),'(i3)') UBOUND(output_fields(out_num)%buffer,3)
       WRITE(error_string2(54:56),'(i3)') output_fields(out_num)%imin
       WRITE(error_string2(58:60),'(i3)') output_fields(out_num)%imax
       WRITE(error_string2(62:64),'(i3)') output_fields(out_num)%jmin
       WRITE(error_string2(66:68),'(i3)') output_fields(out_num)%jmax
       WRITE(error_string2(70:72),'(i3)') output_fields(out_num)%kmin
       WRITE(error_string2(74:76),'(i3)') output_fields(out_num)%kmax
       err_msg = 'module/output_field='//TRIM(error_string1)//&
            & '  Bounds of buffer exceeded.  '//TRIM(error_string2)
       !   imax, imin, etc need to be reset in case the program is not terminated.
       output_fields(out_num)%imax = 0
       output_fields(out_num)%imin = VERY_LARGE_AXIS_LENGTH
       output_fields(out_num)%jmax = 0
       output_fields(out_num)%jmin = VERY_LARGE_AXIS_LENGTH
       output_fields(out_num)%kmax = 0
       output_fields(out_num)%kmin = VERY_LARGE_AXIS_LENGTH
    ELSE
       err_msg = ''
    END IF

  END SUBROUTINE check_out_of_bounds
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="check_bounds_are_exact_dynamic">
  !   <OVERVIEW>
  !     Check if the array indices for <TT>output_fields(out_num)</TT> are equal to the <TT>output_fields(out_num)%buffer</TT>
  !     upper and lower bounds.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE check_bounds_are_exact_dynamic(out_num, diag_field_id, Time, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>check_bounds_are_exact_dynamic</TT> checks if the min and max array indices for <TT>output_fields(out_num)</TT> are
  !     equal to the upper and lower bounds of <TT>output_fields(out_num)%buffer</TT>.  This check is only performed if
  !     <TT>output_fields(out_num)%Time_of_prev_field_data</TT> doesn't equal <TT>Time</TT> or <TT>Time_zero</TT>.
  !     <TT>check_bounds_are_exact_dynamic</TT> returns an error string if the array indices do not match the buffer bounds.
  !   </DESCRIPTION>
  !   <IN NAME="out_num" TYPE="INTEGER">
  !     Output field ID number.
  !   </IN>
  !   <IN NAME="diag_field_id" TYPE="INTEGER">
  !     Input field ID number.
  !   </IN>
  !   <IN NAME="Time" TYPE="TYPE(time_type)">
  !     Time to use in check.  The check is only performed if <TT>output_fields(out_num)%Time_of_prev_field_data</TT> is not
  !     equal to <TT>Time</TT> or <TT>Time_zero</TT>.
  !   </IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*)">
  !     Return status of <TT>check_bounds_are_exact_dynamic</TT>.  An empty error string indicates the x, y, and z indices are
  !     equal to the buffer array boundaries.
  !   </OUT>
  SUBROUTINE check_bounds_are_exact_dynamic(out_num, diag_field_id, Time, err_msg)
    INTEGER, INTENT(in) :: out_num, diag_field_id
    TYPE(time_type), INTENT(in) :: Time
    CHARACTER(len=*), INTENT(out) :: err_msg

    CHARACTER(len=128) :: error_string1, error_string2
    LOGICAL :: do_check

    err_msg = ''

    ! Check bounds only when the value of Time changes. When windows are used,
    ! a change in Time indicates that a new loop through the windows has begun,
    !  so a check of the previous loop can be done.
    IF ( Time == output_fields(out_num)%Time_of_prev_field_data ) THEN
       do_check = .FALSE.
    ELSE
       IF ( output_fields(out_num)%Time_of_prev_field_data == Time_zero ) THEN
          ! It may or may not be OK to check, I don't know how to tell.
          ! Check will be done on subsequent calls anyway.
          do_check = .FALSE.
       ELSE
          do_check = .TRUE.
       END IF
       output_fields(out_num)%Time_of_prev_field_data = Time
    END IF

    IF ( do_check ) THEN
       IF (   output_fields(out_num)%imin /= LBOUND(output_fields(out_num)%buffer,1) .OR.&
            & output_fields(out_num)%imax /= UBOUND(output_fields(out_num)%buffer,1) .OR.&
            & output_fields(out_num)%jmin /= LBOUND(output_fields(out_num)%buffer,2) .OR.&
            & output_fields(out_num)%jmax /= UBOUND(output_fields(out_num)%buffer,2) .OR.&
            & output_fields(out_num)%kmin /= LBOUND(output_fields(out_num)%buffer,3) .OR.&
            & output_fields(out_num)%kmax /= UBOUND(output_fields(out_num)%buffer,3) ) THEN
          WRITE(error_string1,'(a,"/",a)') TRIM(input_fields(diag_field_id)%module_name),&
               & TRIM(output_fields(out_num)%output_name)
          error_string2 ='Buffer bounds=   :   ,   :   ,   :     Actual bounds=   :   ,   :   ,   :   '
          WRITE(error_string2(15:17),'(i3)') LBOUND(output_fields(out_num)%buffer,1)
          WRITE(error_string2(19:21),'(i3)') UBOUND(output_fields(out_num)%buffer,1)
          WRITE(error_string2(23:25),'(i3)') LBOUND(output_fields(out_num)%buffer,2)
          WRITE(error_string2(27:29),'(i3)') UBOUND(output_fields(out_num)%buffer,2)
          WRITE(error_string2(31:33),'(i3)') LBOUND(output_fields(out_num)%buffer,3)
          WRITE(error_string2(35:37),'(i3)') UBOUND(output_fields(out_num)%buffer,3)
          WRITE(error_string2(54:56),'(i3)') output_fields(out_num)%imin
          WRITE(error_string2(58:60),'(i3)') output_fields(out_num)%imax
          WRITE(error_string2(62:64),'(i3)') output_fields(out_num)%jmin
          WRITE(error_string2(66:68),'(i3)') output_fields(out_num)%jmax
          WRITE(error_string2(70:72),'(i3)') output_fields(out_num)%kmin
          WRITE(error_string2(74:76),'(i3)') output_fields(out_num)%kmax
          err_msg = TRIM(error_string1)//' Bounds of data do not match those of buffer. '//TRIM(error_string2)
       END IF
       output_fields(out_num)%imax = 0
       output_fields(out_num)%imin = VERY_LARGE_AXIS_LENGTH
       output_fields(out_num)%jmax = 0
       output_fields(out_num)%jmin = VERY_LARGE_AXIS_LENGTH
       output_fields(out_num)%kmax = 0
       output_fields(out_num)%kmin = VERY_LARGE_AXIS_LENGTH
    END IF
  END SUBROUTINE check_bounds_are_exact_dynamic
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="check_bounds_are_exact_static">
  !   <OVERVIEW>
  !     Check if the array indices for <TT>output_fields(out_num)</TT> are equal to the <TT>output_fields(out_num)%buffer</TT>
  !     upper and lower bounds.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE check_bounds_are_exact_static(out_num, diag_field_id, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !   </DESCRIPTION>
  !   <IN NAME="out_num" TYPE="INTEGER">Output field ID</IN>
  !   <IN NAME="diag_field_id" TYPE="INTEGER">Input field ID.</IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*)"></OUT>
  SUBROUTINE check_bounds_are_exact_static(out_num, diag_field_id, err_msg)
    INTEGER, INTENT(in) :: out_num, diag_field_id
    CHARACTER(len=*), INTENT(out) :: err_msg

    CHARACTER(len=128)  :: error_string1, error_string2

    err_msg = ''

    IF (   output_fields(out_num)%imin /= LBOUND(output_fields(out_num)%buffer,1) .OR.&
         & output_fields(out_num)%imax /= UBOUND(output_fields(out_num)%buffer,1) .OR.&
         & output_fields(out_num)%jmin /= LBOUND(output_fields(out_num)%buffer,2) .OR.&
         & output_fields(out_num)%jmax /= UBOUND(output_fields(out_num)%buffer,2) .OR.&
         & output_fields(out_num)%kmin /= LBOUND(output_fields(out_num)%buffer,3) .OR.&
         & output_fields(out_num)%kmax /= UBOUND(output_fields(out_num)%buffer,3) ) THEN
       WRITE(error_string1,'(a,"/",a)') TRIM(input_fields(diag_field_id)%module_name),&
            & TRIM(output_fields(out_num)%output_name)
       error_string2 ='Buffer bounds=   :   ,   :   ,   :     Actual bounds=   :   ,   :   ,   :   '
       WRITE(error_string2(15:17),'(i3)') LBOUND(output_fields(out_num)%buffer,1)
       WRITE(error_string2(19:21),'(i3)') UBOUND(output_fields(out_num)%buffer,1)
       WRITE(error_string2(23:25),'(i3)') LBOUND(output_fields(out_num)%buffer,2)
       WRITE(error_string2(27:29),'(i3)') UBOUND(output_fields(out_num)%buffer,2)
       WRITE(error_string2(31:33),'(i3)') LBOUND(output_fields(out_num)%buffer,3)
       WRITE(error_string2(35:37),'(i3)') UBOUND(output_fields(out_num)%buffer,3)
       WRITE(error_string2(54:56),'(i3)') output_fields(out_num)%imin
       WRITE(error_string2(58:60),'(i3)') output_fields(out_num)%imax
       WRITE(error_string2(62:64),'(i3)') output_fields(out_num)%jmin
       WRITE(error_string2(66:68),'(i3)') output_fields(out_num)%jmax
       WRITE(error_string2(70:72),'(i3)') output_fields(out_num)%kmin
       WRITE(error_string2(74:76),'(i3)') output_fields(out_num)%kmax
       err_msg = TRIM(error_string1)//' Bounds of data do not match those of buffer. '//TRIM(error_string2)
    END IF
    output_fields(out_num)%imax = 0
    output_fields(out_num)%imin = VERY_LARGE_AXIS_LENGTH
    output_fields(out_num)%jmax = 0
    output_fields(out_num)%jmin = VERY_LARGE_AXIS_LENGTH
    output_fields(out_num)%kmax = 0
    output_fields(out_num)%kmin = VERY_LARGE_AXIS_LENGTH
    
  END SUBROUTINE check_bounds_are_exact_static
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="init_file">
  !   <OVERVIEW>
  !     Initialize the output file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE init_file(name, output_freq, output_units, format, time_units
  !     long_name, tile_count, new_file_freq, new_file_freq_units, start_time,
  !     file_duration, file_duration_units)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Initialize the output file.
  !   </DESCRIPTION>
  !   <IN NAME="name" TYPE="CHARACTER(len=*)">File name.</IN>
  !   <IN NAME="output_freq" TYPE="INTEGER">How often data is to be written to the file.</IN>
  !   <IN NAME="output_units" TYPE="INTEGER">The output frequency unit.  (MIN, HOURS, DAYS, etc.)</IN>
  !   <IN NAME="format" TYPE="INTEGER">Number type/kind the data is to be written out to the file.</IN>
  !   <IN NAME="time_units" TYPE="INTEGER">Time axis units.</IN>
  !   <IN NAME="log_name" TYPE="CHARACTER(len=*)">Long name for time axis.</IN>
  !   <IN NAME="tile_count" TYPE="INTEGER">Tile number.</IN>
  !   <IN NAME="new_file_freq" TYPE="INTEGER, OPTIONAL">How often a new file is to be created.</IN>
  !   <IN NAME="new_file_freq_units" TYPE="INTEGER, OPTIONAL">The new file frequency unit.  (MIN, HOURS, DAYS, etc.)</IN>
  !   <IN NAME="start_time" TYPE="TYPE(time_type), OPTIONAL">Time when the file is to start </IN>
  !   <IN NAME="file_duration" TYPE="INTEGER, OPTIONAL">How long file is to be used.</IN>
  !   <IN NAME="file_duration_units" TYPE="INTEGER, OPTIONAL">File duration unit.  (MIN, HOURS, DAYS, etc.)</IN>
  SUBROUTINE init_file(name, output_freq, output_units, FORMAT, time_units, long_name, tile_count,&
       & new_file_freq, new_file_freq_units, start_time, file_duration, file_duration_units)
    CHARACTER(len=*), INTENT(in) :: name, long_name
    INTEGER, INTENT(in) :: output_freq, output_units, FORMAT, time_units
    INTEGER, INTENT(in) :: tile_count
    INTEGER, INTENT(in), OPTIONAL :: new_file_freq, new_file_freq_units
    INTEGER, INTENT(in), OPTIONAL :: file_duration, file_duration_units
    TYPE(time_type), INTENT(in), OPTIONAL :: start_time

    INTEGER :: new_file_freq1, new_file_freq_units1
    INTEGER :: file_duration1, file_duration_units1
    REAL, DIMENSION(1) :: tdata
    CHARACTER(len=128) :: time_units_str

    ! Get a number for this file
    num_files = num_files + 1
    IF ( num_files >= max_files ) THEN
       ! <ERROR STATUS="FATAL">
       !   max_files exceeded, increase max_files via the max_files variable
       !   in the namelist diag_manager_nml.
       ! </ERROR>
       CALL error_mesg('diag_util_mod::init_file',&
            & ' max_files exceeded, increase max_files via the max_files variable&
            & in the namelist diag_manager_nml.', FATAL)
    END IF

    IF ( PRESENT(new_file_freq) ) THEN 
       new_file_freq1 = new_file_freq
    ELSE 
       new_file_freq1 = VERY_LARGE_FILE_FREQ
    END IF
    
    IF ( PRESENT(new_file_freq_units) ) THEN 
       new_file_freq_units1 = new_file_freq_units 
    ELSE IF ( get_calendar_type() == NO_CALENDAR ) THEN
       new_file_freq_units1 = DIAG_DAYS
    ELSE 
       new_file_freq_units1 = DIAG_YEARS
    END IF
    
    IF ( PRESENT(file_duration) ) THEN
       file_duration1 = file_duration 
    ELSE
       file_duration1 = new_file_freq1
    END IF
    
    IF ( PRESENT(file_duration_units) ) THEN 
       file_duration_units1 = file_duration_units
    ELSE 
       file_duration_units1 = new_file_freq_units1
    END IF
    
    files(num_files)%tile_count = tile_count
    files(num_files)%name = TRIM(name)
    files(num_files)%output_freq = output_freq
    files(num_files)%output_units = output_units
    files(num_files)%format = FORMAT
    files(num_files)%time_units = time_units
    files(num_files)%long_name = TRIM(long_name)
    files(num_files)%num_fields = 0
    files(num_files)%local = .FALSE.
    files(num_files)%last_flush = base_time
    files(num_files)%file_unit = -1
    files(num_files)%new_file_freq = new_file_freq1
    files(num_files)%new_file_freq_units = new_file_freq_units1
    files(num_files)%duration = file_duration1
    files(num_files)%duration_units = file_duration_units1
    IF ( PRESENT(start_time) ) THEN 
       files(num_files)%start_time = start_time
    ELSE
       files(num_files)%start_time = base_time
    END IF
    files(num_files)%next_open=diag_time_inc(files(num_files)%start_time,new_file_freq1,new_file_freq_units1)
    files(num_files)%close_time = diag_time_inc(files(num_files)%start_time,file_duration1, file_duration_units1)
    IF ( files(num_files)%close_time>files(num_files)%next_open ) THEN
       ! <ERROR STATUS="FATAL">
       !   close time GREATER than next_open time, check file duration,
       !   file frequency in <files(num_files)%name>
       ! </ERROR>
       CALL error_mesg('diag_util_mod::init_file', 'close time GREATER than next_open time, check file duration,&
            & file frequency in '//files(num_files)%name, FATAL)
    END IF
    
    ! add time_axis_id and time_bounds_id here
    WRITE(time_units_str, 11) TRIM(time_unit_list(files(num_files)%time_units)), base_year,&
         & base_month, base_day, base_hour, base_minute, base_second
11  FORMAT(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
    files(num_files)%time_axis_id = diag_axis_init (TRIM(long_name), tdata, time_units_str, 'T',&
         & TRIM(long_name) , set_name=TRIM(name) )
    !---- register axis for storing time boundaries
    files(num_files)%time_bounds_id = diag_axis_init( 'nv',(/1.,2./),'none','N','vertex number',&
         & set_name=TRIM(name))
  END SUBROUTINE init_file
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="sync_file_times">
  !   <OVERVIEW>
  !     Synchronize the file's start and close times with the model start and end times.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE sync_file_times(init_time)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>sync_file_times</TT> checks to see if the file start time is less than the
  !     model's init time (passed in as the only argument).  If it is less, then the
  !     both the file start time and end time are synchronized using the passed in initial time
  !     and the duration as calculated by the <TT>diag_time_inc</TT> function.  <TT>sync_file_times</TT>
  !     will also increase the <TT>next_open</TT> until it is greater than the init_time.
  !   </DESCRIPTION>
  !   <IN NAME="file_id" TYPE="INTEGER">The file ID</IN>
  !   <IN NAME="init_time" TYPE="TYPE(time_type)">Initial time use for the synchronization.</IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL">Return error message</OUT>
  SUBROUTINE sync_file_times(file_id, init_time, err_msg)
    INTEGER, INTENT(in) :: file_id
    TYPE(time_type), INTENT(in) :: init_time
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: err_msg

    CHARACTER(len=128) :: msg

    IF ( PRESENT(err_msg) ) err_msg = ''

    IF ( files(file_id)%start_time < init_time ) THEN
       ! Sync the start_time of the file with the initial time of the model
       files(file_id)%start_time = init_time
       ! Sync the file's close time also
       files(file_id)%close_time = diag_time_inc(files(file_id)%start_time,&
            & files(file_id)%duration, files(file_id)%duration_units)
    END IF

    ! Need to increase next_open until it is greate than init_time
    DO WHILE ( files(file_id)%next_open <= init_time )
       files(file_id)%next_open = diag_time_inc(files(file_id)%next_open,&
            & files(file_id)%new_file_freq, files(file_id)%new_file_freq_units, err_msg=msg)
       IF ( msg /= '' ) THEN
          IF ( fms_error_handler('diag_util_mod::sync_file_times',&
               & ' file='//TRIM(files(file_id)%name)//': '//TRIM(msg), err_msg) ) RETURN
       END IF
    END DO
  END SUBROUTINE sync_file_times
  ! </SUBROUTINE>

  ! <FUNCTION NAME="diag_time_inc">
  !   <OVERVIEW>
  !     Return the next time data/file is to be written based on the frequency and units.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     TYPE(time_type) FUNCTION diag_time_inc(time, output_freq, output_units, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the next time data/file is to be written.  This value is based on the current time and the frequency and units.
  !     Function completed successful if the optional <TT>err_msg</TT> is empty.
  !   </DESCRIPTION>
  !   <IN NAME="time" TYPE="TYPE(time_type)">Current model time.</IN>
  !   <IN NAME="output_freq" TYPE="INTEGER">Output frequency number value.</IN>
  !   <IN NAME="output_units" TYPE="INTEGER">Output frequency unit.</IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER, OPTIONAL">
  !     Function error message.  An empty string indicates the next output time was found successfully.
  !   </OUT>
  TYPE(time_type) FUNCTION diag_time_inc(time, output_freq, output_units, err_msg)
    TYPE(time_type), INTENT(in) :: time
    INTEGER, INTENT(in):: output_freq, output_units
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    CHARACTER(len=128) :: error_message_local

    IF ( PRESENT(err_msg) ) err_msg = ''
    error_message_local = ''

    ! special values for output frequency are -1 for output at end of run
    ! and 0 for every timestep.  Need to check for these here?
    ! Return zero time increment, hopefully this value is never used
    IF ( output_freq == END_OF_RUN .OR. output_freq == EVERY_TIME ) THEN
       diag_time_inc = time
       RETURN
    END IF

    ! Make sure calendar was not set after initialization
    IF ( output_units == DIAG_SECONDS ) THEN
       IF ( get_calendar_type() == NO_CALENDAR ) THEN
          diag_time_inc = increment_time(time, output_freq, 0, err_msg=error_message_local)
       ELSE
          diag_time_inc = increment_date(time, 0, 0, 0, 0, 0, output_freq, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_MINUTES ) THEN
       IF ( get_calendar_type() == NO_CALENDAR ) THEN
          diag_time_inc = increment_time(time, NINT(output_freq*SECONDS_PER_MINUTE), 0, &
               &err_msg=error_message_local)
       ELSE
          diag_time_inc = increment_date(time, 0, 0, 0, 0, output_freq, 0, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_HOURS ) THEN
       IF ( get_calendar_type() == NO_CALENDAR ) THEN
          diag_time_inc = increment_time(time, NINT(output_freq*SECONDS_PER_HOUR), 0, err_msg=error_message_local)
       ELSE
          diag_time_inc = increment_date(time, 0, 0, 0, output_freq, 0, 0, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_DAYS ) THEN
       IF (get_calendar_type() == NO_CALENDAR) THEN
          diag_time_inc = increment_time(time, 0, output_freq, err_msg=error_message_local)
       ELSE
          diag_time_inc = increment_date(time, 0, 0, output_freq, 0, 0, 0, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_MONTHS ) THEN
       IF (get_calendar_type() == NO_CALENDAR) THEN
          error_message_local = 'output units of months NOT allowed with no calendar'
       ELSE
          diag_time_inc = increment_date(time, 0, output_freq, 0, 0, 0, 0, err_msg=error_message_local)
       END IF
    ELSE IF ( output_units == DIAG_YEARS ) THEN
       IF ( get_calendar_type() == NO_CALENDAR ) THEN
          error_message_local = 'output units of years NOT allowed with no calendar'
       ELSE
          diag_time_inc = increment_date(time, output_freq, 0, 0, 0, 0, 0, err_msg=error_message_local)
       END IF
    ELSE 
       error_message_local = 'illegal output units'
    END IF

    IF ( error_message_local /= '' ) THEN
       IF ( fms_error_handler('diag_time_inc',error_message_local,err_msg) ) RETURN
    END IF
  END FUNCTION diag_time_inc
  ! </FUNCTION>

  ! <PRIVATE>
  ! <FUNCTION NAME="find_file">
  !   <OVERVIEW>
  !     Return the file number for file name and tile.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION fild_file(name, time_count)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Find the file number for the file name and tile number given.  A return value of <TT>-1</TT> indicates the file was not found.
  !   </DESCRIPTION>
  !   <IN NAME="name=" TYPE="CHARACTER(len=*)">File name.</IN>
  !   <IN NAME="tile_count" TYPE="INTEGER">Tile number.</IN>
  INTEGER FUNCTION find_file(name, tile_count)
    INTEGER, INTENT(in) :: tile_count
    CHARACTER(len=*), INTENT(in) :: name

    INTEGER :: i

    find_file = -1
    DO i = 1, num_files
       IF( TRIM(files(i)%name) == TRIM(name) .AND. tile_count == files(i)%tile_count ) THEN
          find_file = i
          RETURN
       END IF
    END DO
  END FUNCTION find_file
  ! </FUNCTION>
  ! </PRIVATE>

  ! <FUNCTION NAME="find_input_field">
  !   <OVERVIEW>
  !     Return the field number for the given module name, field name, and tile number.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION find_input_field(module_name, field_name, tile_count)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the field number for the given module name, field name and tile number.  A return value of <TT>-1</TT> indicates
  !     the field was not found.
  !   </DESCRIPTION>
  !   <IN NAME="module_name" TYPE="CHARACTER(len=*)">Module name.</IN>
  !   <IN NAME="field_name" TYPE="CHARACTER(len=*)">field name.</IN>
  !   <IN NAME="tile_count" TYPE="INTEGER">Tile number.</IN>
  INTEGER FUNCTION find_input_field(module_name, field_name, tile_count)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name
    INTEGER, INTENT(in) :: tile_count

    INTEGER :: i

    find_input_field = -1 ! Default return value if not found.
    DO i = 1, num_input_fields
       IF(tile_count == input_fields(i)%tile_count .AND.&
            & TRIM(input_fields(i)%module_name) == TRIM(module_name) .AND.&
            & lowercase(TRIM(input_fields(i)%field_name)) == lowercase(TRIM(field_name))) THEN 
          find_input_field = i
          RETURN
       END IF
    END DO
  END FUNCTION find_input_field
  ! </FUNCTION>

  ! <SUBROUTINE NAME="init_input_field">
  !   <OVERVIEW>
  !     Initialize the input field.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE init_input_field(module_name, field_name, tile_count)
  !   </TEMPLATE>
  !     Initialize the input field.
  !   <DESCRIPTION>
  !   </DESCRIPTION>
  !   <IN NAME="module_name" TYPE="CHARACTER(len=*)">Module name.</IN>
  !   <IN NAME="field_name" TYPE="CHARACTER(len=*)">Input field name.</IN>
  !   <IN NAME="tile_count" TYPE="INTEGER">Tile number.</IN>
  SUBROUTINE init_input_field(module_name, field_name, tile_count)
    CHARACTER(len=*),  INTENT(in) :: module_name, field_name
    INTEGER, INTENT(in) :: tile_count

    ! Get a number for this input_field if not already set up
    IF ( find_input_field(module_name, field_name, tile_count) < 0 ) THEN
       num_input_fields = num_input_fields + 1
       IF ( num_input_fields > max_input_fields ) THEN
          ! <ERROR STATUS="FATAL">max_input_fields exceeded, increase it via diag_manager_nml</ERROR>
          CALL error_mesg('diag_util_mod::init_input_field',&
               & 'max_input_fields exceeded, increase it via diag_manager_nml', FATAL)
       END IF
    ELSE
       ! If this is already initialized do not need to do anything
       RETURN
    END IF

    input_fields(num_input_fields)%module_name = TRIM(module_name)
    input_fields(num_input_fields)%field_name = TRIM(field_name)
    input_fields(num_input_fields)%num_output_fields = 0
    ! Set flag that this field has not been registered
    input_fields(num_input_fields)%register = .FALSE.
    input_fields(num_input_fields)%local = .FALSE.
    input_fields(num_input_fields)%standard_name = 'none'
    input_fields(num_input_fields)%tile_count = tile_count
    input_fields(num_input_fields)%numthreads = 1
    input_fields(num_input_fields)%time = time_zero
  END SUBROUTINE init_input_field
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="init_output_field">
  !   <OVERVIEW>
  !     Initialize the output field.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE init_output_field(module_name, field_name, output_name, output_file
  !     time_method, pack, tile_count, local_coord)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Initialize the output field.
  !   </DESCRIPTION>
  !   <IN NAME="module_name" TYPE="CHARACTER(len=*)">Module name.</IN>
  !   <IN NAME="field_name" TYPE="CHARACTER(len=*)">Output field name.</IN>
  !   <IN NAME="output_name" TYPE="CHARACTER(len=*)">Output name written to file.</IN>
  !   <IN NAME="output_file" TYPE="CHARACTER(len=*)">File where field should be written.</IN>
  !   <IN NAME="time_method" TYPE="CHARACTER(len=*)">
  !     Data reduction method.  See <LINK SRC="diag_manager.html">diag_manager_mod</LINK> for valid methods.</IN>
  !   <IN NAME="pack" TYPE="INTEGER">Packing method.</IN>
  !   <IN NAME="tile_count" TYPE="INTEGER">Tile number.</IN>
  !   <IN NAME="local_coord" TYPE="INTEGER, OPTIONAL">Region to be written.  If missing, then all data to be written.</IN>
  SUBROUTINE init_output_field(module_name, field_name, output_name, output_file,&
       & time_method, pack, tile_count, local_coord)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name, output_name, output_file
    CHARACTER(len=*), INTENT(in) :: time_method
    INTEGER, INTENT(in) :: pack
    INTEGER, INTENT(in) :: tile_count
    TYPE(coord_type), INTENT(in), OPTIONAL :: local_coord
    INTEGER :: out_num, in_num, file_num, file_num_tile1
    INTEGER :: num_fields, i, method_selected, l1
    INTEGER :: ioerror
    CHARACTER(len=128) :: error_msg
    CHARACTER(len=50) :: t_method

    ! Get a number for this output field
    num_output_fields = num_output_fields + 1
    IF ( num_output_fields > max_output_fields ) THEN
       ! <ERROR STATUS="FATAL">max_output_fields = <max_output_fields> exceeded.  Increase via diag_manager_nml</ERROR>
       WRITE (UNIT=error_msg,FMT=*) max_output_fields
       CALL error_mesg('diag_util_mod::init_output_field', 'max_output_fields = '//TRIM(error_msg)//' exceeded.&
            &  Increase via diag_manager_nml', FATAL)
    END IF
    out_num = num_output_fields

    ! First, find the index to the associated input field
    in_num = find_input_field(module_name, field_name, tile_count)
    IF ( in_num < 0 ) THEN
       IF ( tile_count > 1 ) THEN
          WRITE (error_msg,'(A,"/",A,"/",A)') TRIM(module_name),TRIM(field_name),&
               & "tile_count="//TRIM(string(tile_count))
       ELSE
          WRITE (error_msg,'(A,"/",A)') TRIM(module_name),TRIM(field_name)
       END IF
       ! <ERROR STATUS="FATAL">module_name/field_name <module_name>/<field_name>[/tile_count=<tile_count>] NOT registered</ERROR>
       CALL error_mesg('diag_util_mod::init_output_field',&
            & 'module_name/field_name '//TRIM(error_msg)//' NOT registered', FATAL)
    END IF

    ! Add this output field into the list for this input field
    input_fields(in_num)%num_output_fields =&
         & input_fields(in_num)%num_output_fields + 1
    IF ( input_fields(in_num)%num_output_fields > max_out_per_in_field ) THEN
       ! <ERROR STATUS="FATAL">
       !   MAX_OUT_PER_IN_FIELD = <MAX_OUT_PER_IN_FIELD> exceeded for <module_name>/<field_name>, increase MAX_OUT_PER_IN_FIELD
       !   in diag_data.F90.
       ! </ERROR>
       WRITE (UNIT=error_msg,FMT=*) MAX_OUT_PER_IN_FIELD
       CALL error_mesg('diag_util_mod::init_output_field',&
        & 'MAX_OUT_PER_IN_FIELD exceeded for '//TRIM(module_name)//"/"//TRIM(field_name)//&
        &', increase MAX_OUT_PER_IN_FIELD in diag_data.F90', FATAL)
    END IF
    input_fields(in_num)%output_fields(input_fields(in_num)%num_output_fields) = out_num

    ! Also put pointer to input field in this output field
    output_fields(out_num)%input_field = in_num

    ! Next, find the number for the corresponding file
    IF ( TRIM(output_file).EQ.'null' ) THEN
       file_num = max_files
    ELSE
       file_num = find_file(output_file, 1)
       IF ( file_num < 0 ) THEN
          ! <ERROR STATUS="FATAL">
          !   file <file_name> is NOT found in the diag_table.
          ! </ERROR>
          CALL error_mesg('diag_util_mod::init_output_field', 'file '&
               & //TRIM(output_file)//' is NOT found in the diag_table', FATAL)
       END IF
       IF ( tile_count > 1 ) THEN
          file_num_tile1 = file_num
          file_num = find_file(output_file, tile_count)
          IF(file_num < 0) THEN
             CALL init_file(files(file_num_tile1)%name, files(file_num_tile1)%output_freq,&
                  & files(file_num_tile1)%output_units, files(file_num_tile1)%format,&
                  & files(file_num_tile1)%time_units, files(file_num_tile1)%long_name,&
                  & tile_count, files(file_num_tile1)%new_file_freq,&
                  & files(file_num_tile1)%new_file_freq_units, files(file_num_tile1)%start_time,&
                  & files(file_num_tile1)%duration, files(file_num_tile1)%duration_units  )
             file_num = find_file(output_file, tile_count)
             IF ( file_num < 0 ) THEN
                ! <ERROR STATUS="FATAL">
                !   file <output_file> is not initialized for tile_count = <tile_count>
                ! </ERROR>
                CALL error_mesg('diag_util_mod::init_output_field', 'file '//TRIM(output_file)//&
                     & ' is not initialized for tile_count = '//TRIM(string(tile_count)), FATAL)
             END IF
          END IF
       END IF
    END IF

    ! Insert this field into list for this file
    files(file_num)%num_fields = files(file_num)%num_fields + 1
    IF ( files(file_num)%num_fields > MAX_FIELDS_PER_FILE ) THEN
       WRITE (UNIT=error_msg, FMT=*) MAX_FIELDS_PER_FILE
       ! <ERROR STATUS="FATAL">
       !   MAX_FIELDS_PER_FILE = <MAX_FIELDS_PER_FILE> exceeded.  Increase MAX_FIELDS_PER_FILE in diag_data.F90.
       ! </ERROR>
       CALL error_mesg('diag_util_mod::init_output_field',&
            & 'MAX_FIELDS_PER_FILE = '//TRIM(error_msg)//' exceeded.  Increase MAX_FIELDS_PER_FILE in diag_data.F90.', FATAL)
    END IF
    num_fields = files(file_num)%num_fields
    files(file_num)%fields(num_fields) = out_num

    ! Set the file for this output field
    output_fields(out_num)%output_file = file_num

    ! Enter the other data for this output field
    output_fields(out_num)%output_name = TRIM(output_name)
    output_fields(out_num)%pack = pack
    output_fields(out_num)%num_axes = 0
    output_fields(out_num)%total_elements = 0
    output_fields(out_num)%region_elements = 0
    output_fields(out_num)%imax = 0
    output_fields(out_num)%jmax = 0
    output_fields(out_num)%kmax = 0
    output_fields(out_num)%imin = VERY_LARGE_AXIS_LENGTH
    output_fields(out_num)%jmin = VERY_LARGE_AXIS_LENGTH
    output_fields(out_num)%kmin = VERY_LARGE_AXIS_LENGTH

    ! initialize the size of the diurnal axis to 1
    output_fields(out_num)%n_diurnal_samples = 1

    ! Initialize all time method to false
    method_selected = 0
    output_fields(out_num)%time_average = .FALSE.
    output_fields(out_num)%time_min = .FALSE.
    output_fields(out_num)%time_max = .FALSE. 
    output_fields(out_num)%time_ops = .FALSE.
    output_fields(out_num)%written_once = .FALSE.

    t_method = lowercase(time_method)
    ! cannot time average fields output every time
    IF ( files(file_num)%output_freq == EVERY_TIME ) THEN
       output_fields(out_num)%time_average = .FALSE.
       method_selected = method_selected+1
       t_method = 'point'
    ELSEIF ( INDEX(t_method,'diurnal') == 1 ) THEN
       ! get the integer number from the t_method
       READ (UNIT=t_method(8:LEN_TRIM(t_method)), FMT=*, IOSTAT=ioerror) output_fields(out_num)%n_diurnal_samples
       IF ( ioerror /= 0 ) THEN
          ! <ERROR STATUS="FATAL">
          !   could not find integer number of diurnal samples in string "<t_method>"
          ! </ERROR>
          CALL error_mesg('diag_util_mod::init_output_field',&
               & 'could not find integer number of diurnal samples in string "' //TRIM(t_method)//'"', FATAL)
       ELSE IF ( output_fields(out_num)%n_diurnal_samples <= 0 ) THEN
          ! <ERROR STATUS="FATAL">
          !   The integer value of diurnal samples must be greater than zero.
          ! </ERROR>
          CALL error_mesg('diag_util_mod::init_output_field',&
               & 'The integer value of diurnal samples must be greater than zero.', FATAL)
       END IF
       output_fields(out_num)%time_average = .TRUE.
       method_selected = method_selected+1
       t_method='mean'
    ELSE
       SELECT CASE(TRIM(t_method))
       CASE ( '.true.', 'mean', 'average', 'avg' )
          output_fields(out_num)%time_average = .TRUE.
          method_selected = method_selected+1
          t_method = 'mean'
       CASE ( '.false.', 'none', 'point' )
          output_fields(out_num)%time_average = .FALSE.
          method_selected = method_selected+1
          t_method = 'point'
       CASE ( 'maximum', 'max' )
          output_fields(out_num)%time_max = .TRUE.
          l1 = LEN_TRIM(output_fields(out_num)%output_name)
          IF ( output_fields(out_num)%output_name(l1-2:l1) /= 'max' ) &
               output_fields(out_num)%output_name = TRIM(output_name)//'_max'
          method_selected = method_selected+1
          t_method = 'max'        
       CASE ( 'minimum', 'min' )
          output_fields(out_num)%time_min = .TRUE.
          l1 = LEN_TRIM(output_fields(out_num)%output_name)
          IF ( output_fields(out_num)%output_name(l1-2:l1) /= 'min' )&
               & output_fields(out_num)%output_name = TRIM(output_name)//'_min'
          method_selected = method_selected+1
          t_method = 'min'        
       END SELECT
    END IF
    
    ! reconcile logical flags
    output_fields(out_num)%time_ops = output_fields(out_num)%time_min.OR.output_fields(out_num)%time_max&
         & .OR.output_fields(out_num)%time_average

    output_fields(out_num)%phys_window = .FALSE.
    ! need to initialize grid_type = -1(start, end, l_start_indx,l_end_indx etc...)
    IF ( PRESENT(local_coord) ) THEN
       input_fields(in_num)%local = .TRUE.
       input_fields(in_num)%local_coord = local_coord
       IF ( INT(local_coord%xbegin * local_coord%xbegin) == 1 .AND.&
            & INT(local_coord%ybegin * local_coord%ybegin) ==1 ) THEN
          output_fields(out_num)%local_output = .FALSE.
          output_fields(out_num)%need_compute = .FALSE.
          output_fields(out_num)%reduced_k_range = .TRUE.
       ELSE
          output_fields(out_num)%local_output = .TRUE.
          output_fields(out_num)%need_compute = .FALSE.
          output_fields(out_num)%reduced_k_range = .FALSE.
       END IF

       output_fields(out_num)%output_grid%start(1) = local_coord%xbegin
       output_fields(out_num)%output_grid%start(2) = local_coord%ybegin
       output_fields(out_num)%output_grid%start(3) = local_coord%zbegin
       output_fields(out_num)%output_grid%end(1) = local_coord%xend
       output_fields(out_num)%output_grid%end(2) = local_coord%yend
       output_fields(out_num)%output_grid%end(3) = local_coord%zend
       DO i = 1, 3
          output_fields(out_num)%output_grid%l_start_indx(i) = -1
          output_fields(out_num)%output_grid%l_end_indx(i) = -1
          output_fields(out_num)%output_grid%subaxes(i) = -1
       END DO
    ELSE
       output_fields(out_num)%local_output = .FALSE.
       output_fields(out_num)%need_compute = .FALSE.
       output_fields(out_num)%reduced_k_range = .FALSE.
    END IF

    ! <ERROR STATUS="FATAL">
    !   improper time method in diag_table for output field <output_name>
    ! </ERROR>
    IF ( method_selected /= 1 ) CALL error_mesg('diag_util_mod::init_output_field',&
         &'improper time method in diag_table for output field:'//TRIM(output_name),FATAL)

    output_fields(out_num)%time_method = TRIM(t_method)

    ! allocate counters: NOTE that for simplicity we always allocate them, even 
    ! if they are superceeded by 4D "counter" array. This isn't most memory 
    ! efficient, approach, but probably tolerable since they are so small anyway
    ALLOCATE(output_fields(out_num)%count_0d(output_fields(out_num)%n_diurnal_samples))
    ALLOCATE(output_fields(out_num)%num_elements(output_fields(out_num)%n_diurnal_samples))
    output_fields(out_num)%count_0d(:) = 0
    output_fields(out_num)%num_elements(:) = 0
  END SUBROUTINE init_output_field
  ! </SUBROUTINE>
  
  ! <PRIVATE>
  ! <SUBROUTINE NAME="opening_file">
  !   <OVERVIEW>
  !     Open file for output.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE opening_file(file, time)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Open file for output, and write the meta data.  <BB>Warning:</BB> Assumes all data structures have been fully initialized.
  !   </DESCRIPTION>
  !   <IN NAME="file" TYPE="INTEGER">File ID.</IN>
  !   <IN NAME="tile" TYPE="TYPE(time_type)">Tile number.</IN>
  SUBROUTINE opening_file(file, time)
    ! WARNING: Assumes that all data structures are fully initialized
    INTEGER, INTENT(in) :: file
    TYPE(time_type), INTENT(in) :: time  

    REAL, DIMENSION(2) :: DATA
    INTEGER :: j, field_num, input_field_num, num_axes, k 
    INTEGER :: field_num1
    INTEGER :: position
    INTEGER :: dir, edges
    INTEGER :: ntileMe
    INTEGER, ALLOCATABLE :: tile_id(:)
    INTEGER, DIMENSION(1) :: time_axis_id, time_bounds_id
    ! size of this axes array must be at least max num. of
    ! axes per field + 2; the last two elements are for time
    ! and time bounds dimensions
    INTEGER, DIMENSION(6) :: axes 
    LOGICAL :: time_ops, aux_present, match_aux_name
    LOGICAL :: all_scalar_or_1d
    CHARACTER(len=7) :: prefix
    CHARACTER (len = 7) :: avg_name = 'average'
    CHARACTER(len=128) :: time_units, timeb_units, avg, error_string, filename, aux_name,fieldname
    CHARACTER(len=128) :: suffix, base_name
    CHARACTER(len=32) :: time_name, timeb_name,time_longname, timeb_longname, cart_name
    CHARACTER(len=256) :: fname
    TYPE(domain1d) :: domain
    TYPE(domain2d) :: domain2

    aux_present = .FALSE.
    match_aux_name = .FALSE.
    ! it's unlikely that a file starts with word "rregion", need to check anyway.
    IF ( LEN(files(file)%name) >=7 .AND. .NOT.files(file)%local ) THEN
       prefix = files(file)%name(1:7)
       IF ( lowercase(prefix) == 'rregion' ) THEN 
          ! <ERROR STATUS="WARNING">
          !   file name should not start with word "rregion"
          ! </ERROR>
          IF ( mpp_pe() == mpp_root_pe() ) CALL error_mesg('diag_util_mod::opening_file',&
               & 'file name should not start with word "rregion"', WARNING)
       END IF
    END IF
    
    ! Here is where time_units string must be set up; time since base date
    WRITE (time_units, 11) TRIM(time_unit_list(files(file)%time_units)), base_year,&
         & base_month, base_day, base_hour, base_minute, base_second
11  FORMAT(A, ' since ', I4.4, '-', I2.2, '-', I2.2, ' ', I2.2, ':', I2.2, ':', I2.2)
    base_name = files(file)%name
    IF ( files(file)%new_file_freq < VERY_LARGE_FILE_FREQ ) THEN
       position = INDEX(files(file)%name, '%')
       IF ( position > 0 )  THEN
          base_name = base_name(1:position-1)
       ELSE
          ! <ERROR STATUS="FATAL">
          !   filename <files(file)%name> does not contain % for time stamp string
          ! </ERROR>
          CALL error_mesg('diag_util_mod::opening_file',&
               & 'file name '//TRIM(files(file)%name)//' does not contain % for time stamp string', FATAL) 
       END IF
       suffix = get_time_string(files(file)%name, time)
    ELSE
       suffix = ' '
    END IF
    ! Add CVS tag as prefix of filename  (currently not implemented)
    !  i1 = INDEX(tagname,':') + 2
    !  i2 = len_trim(tagname) - 2
    !  if(i2 <=i1)  call error_mesg('diag_util opening_file','error in CVS tagname index',FATAL)
    !  prefix2 = tagname(i1:i2)//'_'
    IF ( files(file)%local ) THEN      
       ! prepend "rregion" to all local files for post processing, the prefix will be removed in postprocessing
       filename = 'rregion'//TRIM(base_name)//TRIM(suffix)
    ELSE
       ! filename = trim(prefix2)//trim(base_name)//trim(suffix)
       filename = TRIM(base_name)//TRIM(suffix)
    END IF

    ! Loop through all fields with this file to output axes
    ! JWD: This is a klooge; need something more robust
    domain2 = NULL_DOMAIN2D
    all_scalar_or_1d = .TRUE.
    DO j = 1, files(file)%num_fields
       field_num = files(file)%fields(j)
       if (output_fields(field_num)%local_output .AND. .NOT. output_fields(field_num)%need_compute) CYCLE
       num_axes = output_fields(field_num)%num_axes
       IF ( num_axes > 1 ) THEN
          all_scalar_or_1d = .FALSE.
          domain2 = get_domain2d ( output_fields(field_num)%axes(1:num_axes) )
          IF ( domain2 .NE. NULL_DOMAIN2D ) EXIT
       END IF
    END DO
    IF( .NOT.all_scalar_or_1d ) THEN
       IF ( domain2 .EQ. NULL_DOMAIN2D ) CALL return_domain(domain2)
       IF ( domain2 .EQ. NULL_DOMAIN2D ) THEN
          ! <ERROR STATUS="FATAL">
          !   Domain not defined through set_domain interface; cannot retrieve tile info
          ! </ERROR>
          CALL error_mesg('diag_util_mod::opening_file',&
               & 'Domain not defined through set_domain interface; cannot retrieve tile info', FATAL)
       END IF
       IF ( mpp_get_ntile_count(domain2) > 1 ) THEN
          ntileMe = mpp_get_current_ntile(domain2)
          ALLOCATE(tile_id(ntileMe))
          tile_id = mpp_get_tile_id(domain2)
          fname = TRIM(filename)
          CALL get_tile_string(filename, TRIM(fname)//'.tile' , tile_id(files(file)%tile_count))
          DEALLOCATE(tile_id)
       END IF
    END IF

    CALL diag_output_init(filename, files(file)%format, global_descriptor,&
         & files(file)%file_unit, all_scalar_or_1d, domain2) 
    files(file)%bytes_written = 0 
    ! Does this file contain time_average fields?
    time_ops = .FALSE.
    DO j = 1, files(file)%num_fields
       field_num = files(file)%fields(j)
       IF ( output_fields(field_num)%time_ops ) THEN
          time_ops = .TRUE.
          EXIT
       END IF
    END DO
    ! Loop through all fields with this file to output axes
    DO j = 1, files(file)%num_fields
       field_num = files(file)%fields(j)
       input_field_num = output_fields(field_num)%input_field
       IF (.NOT.input_fields(input_field_num)%register) THEN
          WRITE (error_string,'(A,"/",A)') TRIM(input_fields(input_field_num)%module_name),&
               & TRIM(input_fields(input_field_num)%field_name)
          IF(mpp_pe() .EQ. mpp_root_pe()) THEN
             ! <ERROR STATUS="WARNING">
             !   module/field_name (<input_fields(input_field_num)%module_name>/<input_fields(input_field_num)%field_name>)
             !   NOT registered
             ! </ERROR>
             CALL error_mesg('diag_util_mod::opening_file',&
                  & 'module/field_name ('//TRIM(error_string)//') NOT registered', WARNING)  
          END IF
          CYCLE
       END IF
       if (output_fields(field_num)%local_output .AND. .NOT. output_fields(field_num)%need_compute) CYCLE

       ! Put the time axis in the axis field
       num_axes = output_fields(field_num)%num_axes
       axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
       ! make sure that axis_id are not -1
       DO k = 1, num_axes
          IF ( axes(k) < 0 ) THEN
             WRITE(error_string,'(a)') output_fields(field_num)%output_name
             ! <ERROR STATUS="FATAL">
             !   ouptut_name <output_fields(field_num)%output_name> has axis_id = -1
             ! </ERROR>
             CALL error_mesg('diag_util_mod::opening_file','output_name '//TRIM(error_string)//&
                  & ' has axis_id = -1', FATAL)
          END IF
       END DO
       ! check if aux is present in any axes
       IF ( .NOT.aux_present ) THEN
          DO k = 1, num_axes
             aux_name = get_axis_aux(axes(k))
             IF ( TRIM(aux_name) /= 'none' ) THEN
                aux_present = .TRUE.
                EXIT
             END IF
          END DO
       END IF

       axes(num_axes + 1) = files(file)%time_axis_id
       CALL write_axis_meta_data(files(file)%file_unit, axes(1:num_axes + 1), time_ops)
       IF ( time_ops ) THEN
          axes(num_axes + 2) = files(file)%time_bounds_id
          CALL write_axis_meta_data(files(file)%file_unit, axes(1:num_axes + 2))     
       END IF
    END DO

    ! Looking for the first NON-static field in a file
    field_num1 = files(file)%fields(1)
    DO j = 1, files(file)%num_fields
       field_num = files(file)%fields(j)
       IF ( output_fields(field_num)%time_ops ) THEN
          field_num1 = field_num
          EXIT
       END IF
    END DO
    DO j = 1, files(file)%num_fields
       field_num = files(file)%fields(j)
       input_field_num = output_fields(field_num)%input_field
       IF (.NOT.input_fields(input_field_num)%register) CYCLE
       IF (output_fields(field_num)%local_output .AND. .NOT. output_fields(field_num)%need_compute) CYCLE
       ! Make sure that 1 file contains either time_average or instantaneous fields
       ! cannot have both time_average and instantaneous in 1 file
       IF ( .NOT.mix_snapshot_average_fields ) THEN
          IF ( (output_fields(field_num)%time_ops.NEQV.output_fields(field_num1)%time_ops) .AND.&
               & .NOT.output_fields(field_num1)%static .AND. .NOT.output_fields(field_num)%static) THEN
             IF ( mpp_pe() == mpp_root_pe() ) THEN
                ! <ERROR STATUS="FATAL">
                !   <files(file)%name> can NOT have BOTH time average AND instantaneous fields.
                !   Create a new file or set mix_snapshot_average_fields=.TRUE. in the namelist diag_manager_nml.
                ! </ERROR>
                CALL error_mesg('diag_util_mod::opening_file','file '//&
                     & TRIM(files(file)%name)//' can NOT have BOTH time average AND instantaneous fields.'//&
                     & ' Create a new file or set mix_snapshot_average_fields=.TRUE. in the namelist diag_manager_nml.' , FATAL)
             END IF
          END IF
       END IF
       ! check if any field has the same name as aux_name
       IF ( aux_present .AND. .NOT.match_aux_name ) THEN
          fieldname = output_fields(field_num)%output_name
          IF ( INDEX(aux_name, TRIM(fieldname)) > 0 ) match_aux_name = .TRUE.   
       END IF

       ! Put the time axis in the axis field
       num_axes = output_fields(field_num)%num_axes
       axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
       IF ( .NOT.output_fields(field_num)%static ) THEN
          num_axes=num_axes+1
          axes(num_axes) = files(file)%time_axis_id
       END IF
       IF(output_fields(field_num)%time_average) THEN
          avg = avg_name
       ELSE IF(output_fields(field_num)%time_max) THEN
          avg = avg_name
       ELSE IF(output_fields(field_num)%time_min) THEN
          avg = avg_name
       ELSE
          avg = " "
       END IF
       IF ( input_fields(input_field_num)%missing_value_present ) THEN
          IF ( LEN_TRIM(input_fields(input_field_num)%interp_method) > 0 ) THEN
             output_fields(field_num)%f_type = write_field_meta_data(files(file)%file_unit,&
                  & output_fields(field_num)%output_name, axes(1:num_axes),&
                  & input_fields(input_field_num)%units,&
                  & input_fields(input_field_num)%long_name,&
                  & input_fields(input_field_num)%range, output_fields(field_num)%pack,&
                  & input_fields(input_field_num)%missing_value, avg_name = avg,&
                  & time_method=output_fields(field_num)%time_method,&
                  & standard_name = input_fields(input_field_num)%standard_name,&
                  & interp_method = input_fields(input_field_num)%interp_method)
          ELSE
             output_fields(field_num)%f_type = write_field_meta_data(files(file)%file_unit,&
                  & output_fields(field_num)%output_name, axes(1:num_axes),&
                  & input_fields(input_field_num)%units,&
                  & input_fields(input_field_num)%long_name,&
                  & input_fields(input_field_num)%range, output_fields(field_num)%pack,&
                  & input_fields(input_field_num)%missing_value, avg_name = avg,&
                  & time_method=output_fields(field_num)%time_method,&
                  & standard_name = input_fields(input_field_num)%standard_name)
          END IF
          ! NEED TO TAKE CARE OF TIME AVERAGING INFO TOO BOTH CASES
       ELSE
          IF ( LEN_TRIM(input_fields(input_field_num)%interp_method) > 0 ) THEN
             output_fields(field_num)%f_type = write_field_meta_data(files(file)%file_unit,&
                  & output_fields(field_num)%output_name, axes(1:num_axes),&
                  & input_fields(input_field_num)%units,&
                  & input_fields(input_field_num)%long_name,&
                  & input_fields(input_field_num)%range, output_fields(field_num)%pack,&
                  & avg_name = avg,&
                  & time_method=output_fields(field_num)%time_method,&
                  & standard_name = input_fields(input_field_num)%standard_name,&
                  & interp_method = input_fields(input_field_num)%interp_method)
          ELSE
             output_fields(field_num)%f_type = write_field_meta_data(files(file)%file_unit,&
                  & output_fields(field_num)%output_name, axes(1:num_axes),&
                  & input_fields(input_field_num)%units,&
                  & input_fields(input_field_num)%long_name,&
                  & input_fields(input_field_num)%range, output_fields(field_num)%pack,&
                  & avg_name = avg,&
                  & time_method=output_fields(field_num)%time_method,&
                  & standard_name = input_fields(input_field_num)%standard_name)
          END IF
       END IF
    END DO

    ! If any of the fields in the file are time averaged, need to output the axes
    ! Use double precision since time axis is double precision
    IF ( time_ops ) THEN
       time_axis_id(1) = files(file)%time_axis_id
       files(file)%f_avg_start = write_field_meta_data(files(file)%file_unit,&
            & avg_name // '_T1', time_axis_id, time_units,&
            & "Start time for average period", pack=pack_size)
       files(file)%f_avg_end = write_field_meta_data(files(file)%file_unit,&
            & avg_name // '_T2', time_axis_id, time_units,&
            & "End time for average period", pack=pack_size)
       files(file)%f_avg_nitems = write_field_meta_data(files(file)%file_unit,&
            & avg_name // '_DT', time_axis_id,&
            & TRIM(time_unit_list(files(file)%time_units)),& 
            & "Length of average period", pack=pack_size)
    END IF

    IF ( time_ops ) THEN
       time_axis_id(1) = files(file)%time_axis_id
       time_bounds_id(1) = files(file)%time_bounds_id
       CALL get_diag_axis( time_axis_id(1), time_name, time_units, time_longname,&
            & cart_name, dir, edges, Domain, DATA)
       CALL get_diag_axis( time_bounds_id(1), timeb_name, timeb_units, timeb_longname,&
            & cart_name, dir, edges, Domain, DATA)     
       files(file)%f_bounds =  write_field_meta_data(files(file)%file_unit,&
            & TRIM(time_name)//'_bounds', (/time_bounds_id,time_axis_id/),&
            & TRIM(time_unit_list(files(file)%time_units)),&
            & TRIM(time_name)//' axis boundaries', pack=pack_size)      
    END IF
    ! Let lower levels know that all meta data has been sent
    CALL done_meta_data(files(file)%file_unit)
    IF( aux_present .AND. .NOT.match_aux_name ) THEN
       ! <ERROR STATUS="WARNING">
       !   one axis has auxiliary but the corresponding field is NOT
       !   found in file <file_name>
       ! </ERROR>
       IF ( mpp_pe() == mpp_root_pe() ) CALL error_mesg('diag_util_mod::opening_file',&
            &'one axis has auxiliary but the corresponding field is NOT found in file '//TRIM(files(file)%name), WARNING)
    END IF
  END SUBROUTINE opening_file
  ! </SUBROUTINE>
  ! </PRIVATE>
  
  ! <PRIVATE>
  ! <FUNCTION NAME="get_time_string">
  !   <OVERVIEW>
  !     This function determines a string based on current time.
  !     This string is used as suffix in output file name
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     CHARACTER(len=128) FUNCTION get_time_string(filename, current_time)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     This function determines a string based on current time.
  !     This string is used as suffix in output file name
  !   </DESCRIPTION>
  !   <IN NAME="filename" TYPE="CHARACTER(len=128)">File name.</IN>
  !   <IN NAME="current_time" TYPE="TYPE(time_type)">Current model time.</IN>
  CHARACTER(len=128) FUNCTION get_time_string(filename, current_time)
    CHARACTER(len=128), INTENT(in) :: filename
    TYPE(time_type), INTENT(in) :: current_time

    INTEGER :: yr1, mo1, dy1, hr1, mi1, sc1  ! get from current time
    INTEGER :: yr2, dy2, hr2, mi2            ! for computing next_level time unit
    INTEGER :: yr1_s, mo1_s, dy1_s, hr1_s, mi1_s, sc1_s ! actual values to write string
    INTEGER :: abs_sec, abs_day              ! component of current_time
    INTEGER :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    INTEGER :: julian_day, i, position, len, first_percent
    CHARACTER(len=1) :: width  ! width of the field in format write
    CHARACTER(len=10) :: format
    CHARACTER(len=20) :: yr, mo, dy, hr, mi, sc        ! string of current time (output)
    CHARACTER(len=128) :: filetail

    format = '("_",i*.*)'
    CALL get_date(current_time, yr1, mo1, dy1, hr1, mi1, sc1)
    len = LEN_TRIM(filename)
    first_percent = INDEX(filename, '%')
    filetail = filename(first_percent:len)
    ! compute year string 
    position = INDEX(filetail, 'yr')
    IF ( position > 0 ) THEN
       width = filetail(position-1:position-1)
       yr1_s = yr1
       format(7:9) = width//'.'//width
       WRITE(yr, format) yr1_s   
       yr2 = 0
    ELSE  
       yr = ' '
       yr2 = yr1 - 1
    END IF
    ! compute month string 
    position = INDEX(filetail, 'mo')
    IF ( position > 0 ) THEN   
       width = filetail(position-1:position-1)
       mo1_s = yr2*12 + mo1  
       format(7:9) = width//'.'//width
       WRITE(mo, format) mo1_s
    ELSE
       mo = ' '
    END IF
    ! compute day string        
    IF ( LEN_TRIM(mo) > 0 ) THEN ! month present
       dy1_s = dy1 
       dy2 = dy1_s - 1
    ELSE IF ( LEN_TRIM(yr) >0 )  THEN ! no month, year present
       ! compute julian day
       IF ( mo1 == 1 ) THEN
          dy1_s = dy1
       ELSE
          julian_day = 0
          DO i = 1, mo1-1
             julian_day = julian_day + days_per_month(i)
          END DO
          IF ( leap_year(current_time) .AND. mo1 > 2 ) julian_day = julian_day + 1
          julian_day = julian_day + dy1
          dy1_s = julian_day
       END IF
       dy2 = dy1_s - 1
    ELSE ! no month, no year
       CALL get_time(current_time, abs_sec, abs_day)
       dy1_s = abs_day  
       dy2 = dy1_s 
    END IF
    position = INDEX(filetail, 'dy')
    IF ( position > 0 ) THEN 
       width = filetail(position-1:position-1)
       FORMAT(7:9) = width//'.'//width
       WRITE(dy, FORMAT) dy1_s
    ELSE
       dy = ' '
    END IF
    ! compute hour string
    IF ( LEN_TRIM(dy) > 0 ) THEN
       hr1_s = hr1
    ELSE
       hr1_s = dy2*24 + hr1
    END IF
    hr2 = hr1_s
    position = INDEX(filetail, 'hr')
    IF ( position > 0 ) THEN
       width = filetail(position-1:position-1)
       format(7:9) = width//'.'//width
       WRITE(hr, format) hr1_s
    ELSE
       hr = ' '
    END IF
    ! compute minute string
    IF ( LEN_TRIM(hr) > 0 ) THEN
       mi1_s = mi1
    ELSE
       mi1_s = hr2*60 + mi1
    END IF
    mi2 = mi1_s
    position = INDEX(filetail, 'mi')
    IF(position>0) THEN
       width = filetail(position-1:position-1)
       format(7:9) = width//'.'//width
       WRITE(mi, format) mi1_s
    ELSE
       mi = ' '
    END IF
    ! compute second string
    IF ( LEN_TRIM(mi) > 0 ) THEN
       sc1_s = sc1
    ELSE
       sc1_s = NINT(mi2*SECONDS_PER_MINUTE) + sc1
    END IF
    position = INDEX(filetail, 'sc')
    IF ( position > 0 ) THEN
       width = filetail(position-1:position-1)
       format(7:9) = width//'.'//width
       WRITE(sc, format) sc1_s
    ELSE
       sc = ' '
    ENDIF
    get_time_string = TRIM(yr)//TRIM(mo)//TRIM(dy)//TRIM(hr)//TRIM(mi)//TRIM(sc)
  END FUNCTION get_time_string
  ! </FUNCTION>
  ! </PRIVATE>

  ! <FUNCTION NAME="get_date_dif">
  !   <OVERVIEW>
  !     Return the difference between two times in units.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     REAL FUNCTION get_date_dif(t2, t1, units)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Calculate and return the difference between the two times given in the unit given using the function <TT>t2 - t1</TT>.
  !   </DESCRIPTION>
  !   <IN NAME="t2" TYPE="TYPE(time_type)">Most recent time.</IN>
  !   <IN NAME="t1" TYPE="TYPE(time_type)">Most distant time.</IN>
  !   <IN NAME="units" TYPE="INTEGER">Unit of return value.</IN>
  REAL FUNCTION get_date_dif(t2, t1, units)
    TYPE(time_type), INTENT(in) :: t2, t1
    INTEGER, INTENT(in) :: units

    INTEGER :: dif_seconds, dif_days
    TYPE(time_type) :: dif_time

    ! Compute time axis label value
    ! <ERROR STATUS="FATAL">
    !   variable t2 is less than in variable t1
    ! </ERROR>
    IF ( t2 < t1 ) CALL error_mesg('diag_util_mod::get_date_dif', &
         & 'in variable t2 is less than in variable t1', FATAL)

    dif_time = t2 - t1

    CALL get_time(dif_time, dif_seconds, dif_days)

    IF ( units == DIAG_SECONDS ) THEN
       get_date_dif = dif_seconds + SECONDS_PER_DAY * dif_days
    ELSE IF ( units == DIAG_MINUTES ) THEN
       get_date_dif = 1440 * dif_days + dif_seconds / SECONDS_PER_MINUTE
    ELSE IF ( units == DIAG_HOURS ) THEN
       get_date_dif = 24 * dif_days + dif_seconds / SECONDS_PER_HOUR
    ELSE IF ( units == DIAG_DAYS ) THEN
       get_date_dif = dif_days + dif_seconds / SECONDS_PER_DAY
    ELSE IF ( units == DIAG_MONTHS ) THEN
       ! <ERROR STATUS="FATAL">
       !   months not supported as output units
       ! </ERROR>
       CALL error_mesg('diag_util_mod::get_date_dif', 'months not supported as output units', FATAL)
    ELSE IF ( units == DIAG_YEARS ) THEN
       ! <ERROR STATUS="FATAL">
       !   years not suppored as output units
       ! </ERROR>
       CALL error_mesg('diag_util_mod::get_date_dif', 'years not supported as output units', FATAL)
    ELSE
       ! <ERROR STATUS="FATAL">
       !   illegal time units
       ! </ERROR>
       CALL error_mesg('diag_util_mod::diag_date_dif', 'illegal time units', FATAL)
    END IF
  END FUNCTION get_date_dif
  ! </FUNCTION>

  ! <SUBROUTINE NAME="diag_data_out">
  !   <OVERVIEW>
  !     Write data out to file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE diag_data_out(file, field, dat, time, fianl_call_in, static_write_in)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Write data out to file, and if necessary flush the buffers.
  !   </DESCRIPTION>
  !   <IN NAME="file" TYPE="INTEGER">File ID.</IN>
  !   <IN NAME="field" TYPE="INTEGER">Field ID.</IN>
  !   <INOUT NAME="dat" TYPE="REAL, DIMENSION(:,:,:,:)">Data to write out.</INOUT>
  !   <IN NAME="time" TYPE="TYPE(time_type)">Current model time.</IN>
  !   <IN NAME="final_call_in" TYPE="LOGICAL, OPTIONAL"><TT>.TRUE.</TT> if this is the last write for file.</IN>
  !   <IN NAME="static_write_in" TYPE="LOGICAL, OPTIONAL"><TT>.TRUE.</TT> if static fields are to be written to file.</IN>
  SUBROUTINE diag_data_out(file, field, dat, time, final_call_in, static_write_in)
    INTEGER, INTENT(in) :: file, field
    REAL, DIMENSION(:,:,:,:), INTENT(inout) :: dat
    TYPE(time_type), INTENT(in) :: time
    LOGICAL, OPTIONAL, INTENT(in):: final_call_in, static_write_in

    LOGICAL :: final_call, do_write, static_write
    INTEGER :: i, num
    REAL :: dif, time_data(2, 1, 1, 1), dt_time(1, 1, 1, 1), start_dif, end_dif

    do_write = .TRUE.
    final_call = .FALSE.
    IF ( PRESENT(final_call_in) ) final_call = final_call_in
    static_write = .FALSE.
    IF ( PRESENT(static_write_in) ) static_write = static_write_in
    dif = get_date_dif(time, base_time, files(file)%time_units)
    ! get file_unit, open new file and close curent file if necessary
    IF ( .NOT.static_write .OR. files(file)%file_unit < 0 ) CALL check_and_open(file, time, do_write)
    IF ( .NOT.do_write ) RETURN  ! no need to write data
    CALL diag_field_out(files(file)%file_unit,output_fields(field)%f_type, dat, dif)
    ! record number of bytes written to this file
    files(file)%bytes_written = files(file)%bytes_written +&
         & (SIZE(dat,1)*SIZE(dat,2)*SIZE(dat,3))*(8/output_fields(field)%pack)
    IF ( .NOT.output_fields(field)%written_once ) output_fields(field)%written_once = .TRUE.
    ! *** inserted this line because start_dif < 0 for static fields ***
    IF ( .NOT.output_fields(field)%static ) THEN 
       start_dif = get_date_dif(output_fields(field)%last_output, base_time,files(file)%time_units)
       IF ( .NOT.mix_snapshot_average_fields ) THEN
          end_dif = get_date_dif(output_fields(field)%next_output, base_time, files(file)%time_units)
       ELSE
          end_dif = dif
       END IF
    END IF

    ! Need to write average axes out;
    DO i = 1, files(file)%num_fields
       num = files(file)%fields(i)
       IF ( output_fields(num)%time_ops .AND. &
            input_fields(output_fields(num)%input_field)%register) THEN
          IF ( num == field ) THEN
             ! Output the axes if this is first time-averaged field
             time_data(1, 1, 1, 1) = start_dif
             CALL diag_field_out(files(file)%file_unit, files(file)%f_avg_start, time_data(1:1,:,:,:), dif)
             time_data(2, 1, 1, 1) = end_dif
             CALL diag_field_out(files(file)%file_unit, files(file)%f_avg_end, time_data(2:2,:,:,:), dif)
             ! Compute the length of the average
             dt_time(1, 1, 1, 1) = end_dif - start_dif
             CALL diag_field_out(files(file)%file_unit, files(file)%f_avg_nitems, dt_time(1:1,:,:,:), dif)

             ! Include boundary variable for CF compliance
             CALL diag_field_out(files(file)%file_unit, files(file)%f_bounds, time_data(1:2,:,:,:), dif)         
             EXIT
          END IF
       END IF
    END DO

    ! If write time is greater (equal for the last call) than last_flush for this file, flush it
    IF ( final_call ) THEN
       IF ( time >= files(file)%last_flush ) THEN
          CALL diag_flush(files(file)%file_unit)
          files(file)%last_flush = time
       END IF
    ELSE
       IF ( time > files(file)%last_flush .AND. (.NOT.conserve_water.OR.debug_diag_manager) ) THEN
          CALL diag_flush(files(file)%file_unit)
          files(file)%last_flush = time
       END IF
    END IF
  END SUBROUTINE diag_data_out
  ! </SUBROUTINE>

  ! <PRIVATE>
  ! <SUBROUTINE NAME="check_and_open">
  !   <OVERVIEW>
  !     Checks if it is time to open a new file.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE check_and_open(file, time, do_write)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Checks if it is time to open a new file. If yes, it first closes the
  !     current file, opens a new file and returns file_unit
  !     previous diag_manager_end is replaced by closing_file and output_setup by opening_file.
  !   </DESCRIPTION>
  !   <IN NAME="file" TYPE="INTEGER">File ID.</IN>
  !   <IN NAME="time" TYPE="TYPE(time_type)">Current model time.</IN>
  !   <OUT NAME="do_write" TYPE="LOGICAL"><TT>.TRUE.</TT> if file is expecting more data to write, <TT>.FALSE.</TT> otherwise.</OUT>
  SUBROUTINE check_and_open(file, time, do_write)
    INTEGER, INTENT(in) :: file
    TYPE(time_type), INTENT(in) :: time
    LOGICAL, INTENT(out) :: do_write

    IF ( time >= files(file)%start_time ) THEN 
       IF ( files(file)%file_unit < 0 ) THEN ! need to open a new file
          CALL opening_file(file, time)
          do_write = .TRUE.
       ELSE
          do_write = .TRUE.
          IF ( time > files(file)%close_time .AND. time < files(file)%next_open ) THEN
             do_write = .FALSE. ! file still open but receives NO MORE data
          ELSE IF ( time > files(file)%next_open ) THEN ! need to close current file and open a new one 
             CALL write_static(file)  ! write all static fields and close this file
             CALL opening_file(file, time)        
             files(file)%start_time = files(file)%next_open
             files(file)%close_time =&
                  & diag_time_inc(files(file)%start_time,files(file)%duration, files(file)%duration_units)  
             files(file)%next_open =&
                  & diag_time_inc(files(file)%next_open, files(file)%new_file_freq,&
                  & files(file)%new_file_freq_units)
             IF ( files(file)%close_time > files(file)%next_open ) THEN 
                ! <ERROR STATUS="FATAL">
                !   <file_name> has close time GREATER than next_open time,
                !   check file duration and frequency
                ! </ERROR>
                CALL error_mesg('diag_util_mod::check_and_open',&
                     & files(file)%name//' has close time GREATER than next_open time, check file duration and frequency',FATAL)
             END IF
          END IF ! no need to open new file, simply return file_unit
       END IF
    ELSE
       do_write = .FALSE.
    END IF
  END SUBROUTINE check_and_open
  ! </SUBROUTINE>
  ! </PRIVATE>

  ! <SUBROUTINE NAME="write_static">
  !   <OVERVIEW>
  !     Output all static fields in this file
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE write_static(file)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Write the static data to the file.
  !   </DESCRIPTION>
  !   <IN NAME="file" TYPE="INTEGER">File ID.</IN>
  SUBROUTINE write_static(file)
    INTEGER, INTENT(in) :: file

    INTEGER :: j, i, input_num

    DO j = 1, files(file)%num_fields
       i = files(file)%fields(j)
       input_num = output_fields(i)%input_field
       ! skip fields that were not registered
       IF ( .NOT.input_fields(input_num)%register ) CYCLE
       if( output_fields(i)%local_output .AND. .NOT. output_fields(i)%need_compute) CYCLE
       ! only output static fields here
       IF ( .NOT.output_fields(i)%static ) CYCLE
       CALL diag_data_out(file, i, output_fields(i)%buffer, files(file)%last_flush, .TRUE., .TRUE.)
    END DO
    ! Close up this file   
    CALL mpp_close(files(file)%file_unit)
    files(file)%file_unit = -1
  END SUBROUTINE write_static
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="check_duplicate_output_fields">
  !   <OVERVIEW>
  !     Checks to see if <TT>output_name</TT> and <TT>output_file</TT> are unique in <TT>output_fields</TT>.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE check_duplicate_output_fields(err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Check to see if <TT>output_name</TT> and <TT>output_file</TT> are unique in <TT>output_fields</TT>.  An empty
  !     <TT>err_msg</TT> indicates no duplicates found.
  !   </DESCRIPTION>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL">Error message.  If empty, then no duplicates found.</OUT>
  SUBROUTINE check_duplicate_output_fields(err_msg)
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    INTEGER :: i, j, tmp_file
    CHARACTER(len=128) :: tmp_name
    CHARACTER(len=256) :: err_msg_local

    IF ( PRESENT(err_msg) ) err_msg=''
    ! Do the checking when more than 1 output_fileds present
    IF ( num_output_fields <= 1 ) RETURN 
    err_msg_local = ''

    i_loop: DO i = 1, num_output_fields-1
       tmp_name = TRIM(output_fields(i)%output_name)
       tmp_file =  output_fields(i)%output_file
       DO j = i+1, num_output_fields
          IF ( (tmp_name == TRIM(output_fields(j)%output_name)) .AND. &
               &(tmp_file == output_fields(j)%output_file)) THEN
             err_msg_local = ' output_field "'//TRIM(tmp_name)//&
                  &'" duplicated in file "'//TRIM(files(tmp_file)%name)//'"'
             EXIT i_loop
          END IF
       END DO
    END DO i_loop
    IF ( err_msg_local /= '' ) THEN
       IF ( fms_error_handler(' ERROR in diag_table',err_msg_local,err_msg) ) RETURN
    END IF
  END SUBROUTINE check_duplicate_output_fields
  ! </SUBROUTINE>
END MODULE diag_util_mod
