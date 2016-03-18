
! -*-f90-*-
! $Id: mpp_io_util.inc,v 18.0 2010/03/02 23:57:42 fms Exp $


  !#####################################################################
! <SUBROUTINE NAME="mpp_get_info">
!   <OVERVIEW>
!     Get some general information about a file.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Get some general information about a file.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_get_info( unit, ndim, nvar, natt, ntime )
!   </TEMPLATE>
!   <IN NAME="unit" TYPE="integer"> </IN>
!   <OUT NAME="ndim" TYPE="integer"> </OUT>
!   <OUT NAME="nvar" TYPE="integer"> </OUT>
!   <OUT NAME="natt" TYPE="integer"> </OUT>
!   <OUT NAME="ntime" TYPE="integer"> </OUT>
! </SUBROUTINE>

    subroutine mpp_get_info( unit, ndim, nvar, natt, ntime )

      integer, intent(in) :: unit
      integer, intent(out) :: ndim, nvar, natt, ntime


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_INFO: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )&
           call mpp_error(FATAL, 'MPP_GET_INFO: invalid unit number, file '//trim(mpp_file(unit)%name))

      ndim = mpp_file(unit)%ndim
      nvar = mpp_file(unit)%nvar
      natt = mpp_file(unit)%natt
      ntime = mpp_file(unit)%time_level

      return

    end subroutine mpp_get_info

  !#####################################################################
! <SUBROUTINE NAME="mpp_get_global_atts" INTERFACE="mpp_get_atts">
!  <IN NAME="unit" TYPE="integer"></IN>
!  <IN NAME="global_atts" TYPE="atttype" DIM="(:)"></IN>
! </SUBROUTINE>
    subroutine mpp_get_global_atts( unit, global_atts )
!
!  copy global file attributes for use by user
!
!  global_atts is an attribute type which is allocated from the
!  calling routine

      integer,       intent(in)    :: unit
      type(atttype), intent(inout) :: global_atts(:)
      integer :: natt,i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_INFO: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )&
           call mpp_error( FATAL, 'MPP_GET_INFO: invalid unit number,file '//trim(mpp_file(unit)%name))

      if (size(global_atts(:)).lt.mpp_file(unit)%natt) &
           call mpp_error(FATAL, 'MPP_GET_ATTS: atttype not dimensioned properly in calling routine, file '// &
           trim(mpp_file(unit)%name))

      natt = mpp_file(unit)%natt
      global_atts = default_att

      do i=1,natt
         global_atts(i) = mpp_file(unit)%Att(i)
      enddo

      return
   end subroutine mpp_get_global_atts

  !#####################################################################
   subroutine mpp_get_field_atts(field, name, units, longname, min, max, missing, ndim, siz, axes, atts, valid, scale, add)

     type(fieldtype), intent(in) :: field
     character(len=*), intent(out) , optional :: name, units
     character(len=*), intent(out), optional :: longname
     real,intent(out), optional :: min,max,missing
     integer, intent(out), optional :: ndim
     integer, intent(out), dimension(:), optional :: siz
     type(validtype), intent(out), optional :: valid
     real,            intent(out), optional :: scale
     real,            intent(out), optional :: add

     type(atttype), intent(inout), optional, dimension(:) :: atts
     type(axistype), intent(inout), optional, dimension(:) :: axes

     integer :: n,m

     if (PRESENT(name)) name = field%name
     if (PRESENT(units)) units = field%units
     if (PRESENT(longname)) longname = field%longname
     if (PRESENT(min)) min = field%min
     if (PRESENT(max)) max = field%max
     if (PRESENT(missing)) missing = field%missing
     if (PRESENT(ndim)) ndim = field%ndim
     if (PRESENT(atts)) then
        atts = default_att
        n = size(atts(:));m=size(field%Att(:))
        if (n.LT.m)&
             call mpp_error(FATAL,'attribute array not large enough in mpp_get_field_atts, field '//&
             trim(field%name))
        do n=1,m
          atts(n) = field%Att(n)
        end do
     end if
     if (PRESENT(axes)) then
        axes = default_axis
        n = size(axes(:));m=field%ndim
        if (n.LT.m) &
             call mpp_error(FATAL,'axis array not large enough in mpp_get_field_atts, field '//&
             trim(field%name))
        do n=1,m
          axes(n) = field%axes(n)
        end do
     end if
     if (PRESENT(siz)) then
        siz = -1
        n = size(siz(:));m=field%ndim
        if (n.LT.m) &
             call mpp_error(FATAL,'size array not large enough in mpp_get_field_atts, field '//&
             trim(field%name))
        do n=1,m
          siz(n) = field%size(n)
        end do
     end if

     if(PRESENT(valid)) then
       call mpp_get_valid(field,valid)
     endif

     if(PRESENT(scale)) scale = field%scale
     if(present(add))   add   = field%add

     return
   end subroutine mpp_get_field_atts

  !#####################################################################
   subroutine mpp_get_axis_atts( axis, name, units, longname, cartesian, calendar, sense, len, natts, atts )

     type(axistype), intent(in) :: axis
     character(len=*), intent(out) , optional :: name, units
     character(len=*), intent(out), optional :: longname, cartesian,calendar
     integer,intent(out), optional :: sense, len , natts
     type(atttype), intent(inout), optional, dimension(:) :: atts

     integer :: n,m

     if (PRESENT(name)) name = axis%name
     if (PRESENT(units)) units = axis%units
     if (PRESENT(longname)) longname = axis%longname
     if (PRESENT(cartesian)) cartesian = axis%cartesian
     if (PRESENT(calendar)) calendar = axis%calendar
     if (PRESENT(sense)) sense = axis%sense
     if (PRESENT(len)) len = axis%len
     if (PRESENT(atts)) then
        atts = default_att
        n = size(atts(:));m=size(axis%Att(:))
        if (n.LT.m) &
             call mpp_error(FATAL,'attribute array not large enough in mpp_get_field_atts, axis '//&
             trim(axis%name))
        do n=1,m
          atts(n) = axis%Att(n)
        end do
     end if
     if (PRESENT(natts)) natts = size(axis%Att(:))

     return
   end subroutine mpp_get_axis_atts


  !#####################################################################
    subroutine mpp_get_fields( unit, variables )
!
!  copy variable information from file (excluding data)
!  global_atts is an attribute type which is allocated from the
!  calling routine
!
      integer,         intent(in)    :: unit
      type(fieldtype), intent(inout) :: variables(:)

      integer :: nvar,i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_FIELDS: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_GET_FIELDS: invalid unit number.' )

      if (size(variables(:)).ne.mpp_file(unit)%nvar) &
          call mpp_error(FATAL,'MPP_GET_FIELDS: fieldtype not dimensioned properly in calling routine, file '//&
          trim(mpp_file(unit)%name))

      nvar = mpp_file(unit)%nvar

      do i=1,nvar
         variables(i) = mpp_file(unit)%Var(i)
      enddo

      return
   end subroutine mpp_get_fields

  !#####################################################################
    subroutine mpp_get_axes( unit, axes, time_axis )
!
!  copy variable information from file (excluding data)
!  global_atts is an attribute type which is allocated from the
!  calling routine
!
      integer, intent(in) :: unit
      type(axistype), intent(inout) :: axes(:)
      type(axistype), intent(inout), optional :: time_axis
      integer :: ndim,i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_AXES: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )&
           call mpp_error( FATAL, 'MPP_GET_AXES: invalid unit number, file '//trim(mpp_file(unit)%name))

      if (size(axes(:)).ne.mpp_file(unit)%ndim) &
           call mpp_error(FATAL, 'MPP_GET_AXES: axistype not dimensioned properly in calling routine, file '//&
           trim(mpp_file(unit)%name))


      if (PRESENT(time_axis)) time_axis = default_axis
      ndim = mpp_file(unit)%ndim

      do i=1,ndim
         axes(i)=mpp_file(unit)%Axis(i)

         if (PRESENT(time_axis) &
             .AND. .NOT. ASSOCIATED(mpp_file(unit)%Axis(i)%data) &
             .AND. mpp_file(unit)%Axis(i)%type /= -1) then
            time_axis = mpp_file(unit)%Axis(i)
         endif
      enddo

      return
   end subroutine mpp_get_axes

  !#####################################################################
! <SUBROUTINE NAME="mpp_get_times">
!   <OVERVIEW>
!     Get file time data.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Get file time data.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_get_times( unit, time_values )
!   </TEMPLATE>
!   <IN NAME="unit" TYPE="integer"> </IN>
!   <INOUT NAME="time_values" TYPE="real(DOUBLE_KIND)" DIM="(:)"> </INOUT>
! </SUBROUTINE>

   subroutine mpp_get_times( unit, time_values )
!
!  copy time information from file and convert to time_type
!
      integer, intent(in) :: unit
      real(DOUBLE_KIND), intent(inout) :: time_values(:)

      integer :: ntime,i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_TIMES: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )&
           call mpp_error(FATAL, 'MPP_GET_TIMES: invalid unit number, file '//trim(mpp_file(unit)%name))

! NF_INQ_DIM returns -1 for the length of a record dimension if
! it does not exist

      if (mpp_file(unit)%time_level == -1) then
          time_values = 0.0
          return
      endif

      if (size(time_values(:)).ne.mpp_file(unit)%time_level) &
         call mpp_error(FATAL,'MPP_GET_TIMES: time_values not dimensioned properly in calling routine, file '//&
         trim(mpp_file(unit)%name))

      ntime = mpp_file(unit)%time_level

      do i=1,ntime
         time_values(i) = mpp_file(unit)%time_values(i)
      enddo

      return
    end subroutine mpp_get_times

  !#####################################################################
   function mpp_get_field_index(fields,fieldname)

     type(fieldtype), dimension(:) :: fields
     character(len=*) :: fieldname
     integer :: mpp_get_field_index

     integer :: n

     mpp_get_field_index = -1

     do n=1,size(fields(:))
        if (lowercase(fields(n)%name) == lowercase(fieldname)) then
           mpp_get_field_index = n
           exit
        endif
     enddo

     return
   end function mpp_get_field_index

  !#####################################################################
   function mpp_get_axis_index(axes,axisname)

     type(axistype), dimension(:) :: axes
     character(len=*) :: axisname
     integer :: mpp_get_axis_index

     integer :: n

     mpp_get_axis_index = -1

     do n=1,size(axes(:))
        if (lowercase(axes(n)%name) == lowercase(axisname)) then
           mpp_get_axis_index = n
           exit
        endif
     enddo

     return
   end function mpp_get_axis_index

  !#####################################################################
   function mpp_get_field_size(field)

     type(fieldtype) :: field
     integer :: mpp_get_field_size(4)

     mpp_get_field_size = -1

     mpp_get_field_size(1) = field%size(1)
     mpp_get_field_size(2) = field%size(2)
     mpp_get_field_size(3) = field%size(3)
     mpp_get_field_size(4) = field%size(4)

     return
   end function mpp_get_field_size


  !#####################################################################
   subroutine mpp_get_axis_data( axis, data )

     type(axistype), intent(in) :: axis
     real, dimension(:), intent(out) :: data


     if (size(data(:)).lt.axis%len)&
          call mpp_error(FATAL,'MPP_GET_AXIS_DATA: data array not large enough, axis '//trim(axis%name))
     if (.NOT.ASSOCIATED(axis%data)) then
        call mpp_error(NOTE,'MPP_GET_AXIS_DATA: use mpp_get_times for record dims')
        data = 0.
     else
        data(1:axis%len) = axis%data
     endif
 
     return
   end subroutine mpp_get_axis_data


  !#####################################################################
   function mpp_get_recdimid(unit)
!
      integer, intent(in) :: unit
      integer  :: mpp_get_recdimid


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_RECDIMID: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_GET_RECDIMID: invalid unit number.' )

      mpp_get_recdimid = mpp_file(unit)%recdimid

      return
   end function mpp_get_recdimid

    subroutine mpp_get_iospec( unit, iospec )
      integer, intent(in) :: unit
      character(len=*), intent(out) :: iospec

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_IOSPEC: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_GET_IOSPEC: invalid unit number.' )
#ifdef SGICRAY
!currently will write to stdout: don't know how to trap and return as string to iospec
      call ASSIGN( 'assign -V f:'//trim(mpp_file(unit)%name), error )
#endif
      return
    end subroutine mpp_get_iospec


  !#####################################################################
! <FUNCTION NAME="mpp_get_ncid">
!   <OVERVIEW>
!     Get netCDF ID of an open file.
!   </OVERVIEW>
!   <DESCRIPTION>
!    This returns the <TT>ncid</TT> associated with the open file on
!    <TT>unit</TT>. It is used in the instance that the user desires to
!    perform netCDF calls upon the file that are not provided by the
!    <TT>mpp_io_mod</TT> API itself.
!   </DESCRIPTION>
!   <TEMPLATE>
!     mpp_get_ncid(unit)
!   </TEMPLATE>
!   <IN NAME="unit" TYPE="integer"> </IN>
! </FUNCTION>

    function mpp_get_ncid(unit)
      integer :: mpp_get_ncid
      integer, intent(in) :: unit

      mpp_get_ncid = mpp_file(unit)%ncid
      return
    end function mpp_get_ncid

  !#####################################################################
    function mpp_get_axis_id(axis)
      integer mpp_get_axis_id
      type(axistype), intent(in) :: axis
      mpp_get_axis_id = axis%id
      return
    end function mpp_get_axis_id

  !#####################################################################
    function mpp_get_field_id(field)
      integer mpp_get_field_id
      type(fieldtype), intent(in) :: field
      mpp_get_field_id = field%id
      return
    end function mpp_get_field_id

  !#####################################################################
    subroutine mpp_get_unit_range( unit_begin_out, unit_end_out )
      integer, intent(out) ::      unit_begin_out, unit_end_out

      unit_begin_out = unit_begin; unit_end_out = unit_end
      return
    end subroutine mpp_get_unit_range

  !#####################################################################
    subroutine mpp_set_unit_range( unit_begin_in, unit_end_in )
      integer, intent(in) ::       unit_begin_in, unit_end_in

      if( unit_begin_in.GT.unit_end_in )call mpp_error( FATAL, 'MPP_SET_UNIT_RANGE: unit_begin_in.GT.unit_end_in.' )
      if( unit_begin_in.LT.0           )call mpp_error( FATAL, 'MPP_SET_UNIT_RANGE: unit_begin_in.LT.0.' )
      if( unit_end_in  .GT.maxunits    )call mpp_error( FATAL, 'MPP_SET_UNIT_RANGE: unit_end_in.GT.maxunits.' )
      unit_begin = unit_begin_in; unit_end = unit_end_in
      return
    end subroutine mpp_set_unit_range

  !#####################################################################
    subroutine mpp_io_set_stack_size(n)
!set the mpp_io_stack variable to be at least n LONG words long
      integer, intent(in) :: n
      character(len=10) :: text

      if( n.GT.mpp_io_stack_size .AND. allocated(mpp_io_stack) )deallocate(mpp_io_stack)
      if( .NOT.allocated(mpp_io_stack) )then
          allocate( mpp_io_stack(n) )
          mpp_io_stack_size = n
          write( text,'(i10)' )n
          if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, 'MPP_IO_SET_STACK_SIZE: stack size set to '//text//'.' )
      end if

      return
    end subroutine mpp_io_set_stack_size

  !#####################################################################
  ! based on presence/absence of attributes, defines valid range or missing
  ! value. For details, see section 8.1 of NetCDF User Guide
  subroutine mpp_get_valid(f,v)
     type(fieldtype),intent(in)  :: f ! field
     type(validtype),intent(out) :: v ! validator

     integer :: irange,imin,imax,ifill,imissing,iscale
     integer :: valid_T, scale_T ! types of attributes

     v%is_range = .true.
     v%min = -HUGE(v%min); v%max = HUGE(v%max)
     ! find indices of relevant attributes
     irange   = mpp_find_att(f%att,'valid_range')
     imin     = mpp_find_att(f%att,'valid_min')
     imax     = mpp_find_att(f%att,'valid_max')
     ifill    = mpp_find_att(f%att,'_FillValue')
     imissing = mpp_find_att(f%att,'missing_value')

     ! find the widest type of scale and offset; note that the code
     ! uses assumption that NetCDF types are arranged in th order of rank,
     ! that is NF_BYTE < NF_SHORT < NF_INT < NF_FLOAT < NF_DOUBLE
     scale_T = 0
     iscale   = mpp_find_att(f%att,'scale_factor')
     if(iscale>0) scale_T = f%att(iscale)%type
     iscale = mpp_find_att(f%att,'add_offest')
     if(iscale>0) scale_T = max(scale_T,f%att(iscale)%type)


     ! examine possible range attributes
     valid_T = 0
     if (irange>0) then
        v%min = f%att(irange)%fatt(1)
        v%max = f%att(irange)%fatt(2)
        valid_T = f%att(irange)%type
     else if (imax>0.or.imin>0) then
        if(imax>0) then
           v%max = f%att(imax)%fatt(1)
           valid_T = max(valid_T,f%att(imax)%type)
        endif
        if(imin>0) then
           v%min = f%att(imin)%fatt(1)
           valid_T = max(valid_T,f%att(imin)%type)
        endif
     else if (imissing > 0) then
        v%is_range = .false.
        ! here we always scale, since missing_value is supposed to be in
        ! external representation
        v%min = f%att(imissing)%fatt(1)*f%scale + f%add
     else if (ifill>0) then
     !z1l ifdef is added in to be able to compile without using use_netCDF.
#ifdef use_netCDF  
        ! define min and max according to _FillValue
        if(f%att(ifill)%fatt(1)>0) then
            ! if _FillValue is positive, then it defines valid maximum
            v%max = f%att(ifill)%fatt(1)
            select case(f%type)
            case (NF_BYTE,NF_SHORT,NF_INT)
               v%max = v%max-1
            case (NF_FLOAT)
               v%max = nearest(nearest(real(v%max,4),-1.0),-1.0)
            case (NF_DOUBLE)
               v%max = nearest(nearest(real(v%max,8),-1.0),-1.0)
            end select
            ! always do the scaling, as the _FillValue is in external
            ! representation
            v%max = v%max*f%scale + f%add
        else
            ! if _FillValue is negative or zero, then it defines valid minimum
            v%min = f%att(ifill)%fatt(1)
            select case(f%type)
            case (NF_BYTE,NF_SHORT,NF_INT)
               v%min = v%min+1
            case (NF_FLOAT)
               v%min = nearest(nearest(real(v%min,4),+1.0),+1.0)
            case (NF_DOUBLE)
               v%min = nearest(nearest(real(v%min,8),+1.0),+1.0)
            end select
            ! always do the scaling, as the _FillValue is in external
            ! representation
            v%min = v%min*f%scale + f%add
        endif
#endif
    endif
	! If valid_range is the same type as scale_factor (actually the wider of
	! scale_factor and add_offset) and this is wider than the external data, then it
	! will be interpreted as being in the units of the internal (unpacked) data.
	! Otherwise it is in the units of the external (packed) data.
    ! Note that it is not relevant if we went through the missing_data of _FillValue
    ! brances, because in this case all irange, imin, and imax are less then 0
    if(.not.((valid_T == scale_T).and.(scale_T>f%type))) then
       if(irange>0 .or. imin>0) then
          v%min = v%min*f%scale + f%add
       endif
       if(irange>0 .or. imax>0) then
          v%max = v%max*f%scale + f%add
       endif
    endif

   end subroutine mpp_get_valid

   !#####################################################################
   logical elemental function mpp_is_valid(x, v)
      real           , intent(in) :: x ! real value to be eaxmined
      type(validtype), intent(in) :: v ! validator

      if (v%is_range) then
         mpp_is_valid = (v%min<=x).and.(x<=v%max)
      else
         mpp_is_valid = x/=v%min
      endif
    end function mpp_is_valid

   !#####################################################################
   ! finds an attribute by name in the array; returns -1 if it is not
   ! found
   function mpp_find_att(atts, name)
     integer                   :: mpp_find_att
     type(atttype), intent(in) :: atts(:) ! array of attributes
     character(len=*)          :: name ! name of the attributes

     integer :: i

     mpp_find_att = -1
     do i = 1, size(atts)
        if (trim(name)==trim(atts(i)%name)) then
           mpp_find_att=i
           exit
        endif
     enddo
   end function mpp_find_att
   !#####################################################################

   ! return the name of an attribute.
   function mpp_get_att_name(att)
      type(atttype),    intent(in) :: att
      character(len=len(att%name)) :: mpp_get_att_name
     
      mpp_get_att_name = att%name
      return

   end function mpp_get_att_name

   !#####################################################################

   ! return the type of an attribute.
   function mpp_get_att_type(att)
      type(atttype), intent(in) :: att
      integer                   :: mpp_get_att_type
     
      mpp_get_att_type = att%type
      return

   end function mpp_get_att_type

   !#####################################################################

   ! return the length of an attribute.
   function mpp_get_att_length(att)
      type(atttype), intent(in) :: att
      integer                   :: mpp_get_att_length
     
      mpp_get_att_length = att%len

      return

   end function mpp_get_att_length

   !#####################################################################

   ! return the char value of an attribute.
   function mpp_get_att_char(att)
      type(atttype), intent(in) :: att
      character(len=att%len)    :: mpp_get_att_char
     
      mpp_get_att_char = att%catt
      return

   end function mpp_get_att_char

   !#####################################################################

   ! return the real array value of an attribute.
   function mpp_get_att_real(att)
      type(atttype), intent(in)          :: att
      real, dimension(size(att%fatt(:))) :: mpp_get_att_real
     
      mpp_get_att_real = att%fatt
      return

   end function mpp_get_att_real

  !#####################################################################

   ! return the real array value of an attribute.
   function mpp_get_att_real_scalar(att)
      type(atttype), intent(in)          :: att
      real                               :: mpp_get_att_real_scalar

      mpp_get_att_real_scalar = att%fatt(1)
      return

   end function mpp_get_att_real_scalar

   !#####################################################################
   ! return the name of an field
   function mpp_get_field_name(field)
      type(fieldtype), intent(in) :: field
      character(len=len(field%name)) :: mpp_get_field_name

      mpp_get_field_name = field%name
      return
   end function mpp_get_field_name

   !#####################################################################
   ! return the  file name of corresponding unit
   function mpp_get_file_name(unit)
      integer,                  intent(in) :: unit
      character(len=len(mpp_file(1)%name)) :: mpp_get_file_name      

      mpp_get_file_name = mpp_file(unit)%name
      return

   end function mpp_get_file_name

   !####################################################################
   ! return if certain file with unit is opened or not
   function mpp_file_is_opened(unit)
      integer,  intent(in) :: unit
      logical              :: mpp_file_is_opened      

      mpp_file_is_opened = mpp_file(unit)%opened
      return

   end function mpp_file_is_opened

   !####################################################################
   ! return the attribute value of given field name
   subroutine mpp_get_field_att_text(unit, fieldname, attname, attvalue)
     integer,           intent(in) :: unit  
     character(len=*),  intent(in) :: fieldname, attname
     character(len=*), intent(out) :: attvalue 
     logical                       :: found_field,  found_att
     integer                       :: i, j

     found_field = .false.
     found_att = .false.
     do i=1,mpp_file(unit)%nvar
        if( trim(mpp_file(unit)%Var(i)%name) == trim(fieldname)) then
           found_field = .true.
           do j=1, size(mpp_file(unit)%Var(i)%Att(:))
              if( trim(mpp_file(unit)%Var(i)%Att(j)%name) == trim(attname) ) then
                 found_att = .true.
                 if(len(attvalue) .LE. len_trim(mpp_file(unit)%Var(i)%Att(j)%catt) ) call mpp_error(FATAL, &
                      'mpp_io_util.inc: length of attvalue is less than the length of catt')
                 attvalue = trim(mpp_file(unit)%Var(i)%Att(j)%catt)
                 exit
              end if
           end do
           exit
        end if
     end do

     if(.NOT. found_field) call mpp_error(FATAL,"mpp_io_util.inc: field "//trim(fieldname)// &
               " does not exist in the file "//trim(mpp_file(unit)%name) )
     if(.NOT. found_att) call mpp_error(FATAL,"mpp_io_util.inc: attribute "//trim(attname)//" of field "&
               //trim(fieldname)// " does not exist in the file "//trim(mpp_file(unit)%name) )

     return

   end subroutine mpp_get_field_att_text


   !####################################################################
   ! return mpp_io_nml variable io_clock_on
   function mpp_io_clock_on()
      logical :: mpp_io_clock_on

      mpp_io_clock_on = io_clocks_on
      return

   end function mpp_io_clock_on   
