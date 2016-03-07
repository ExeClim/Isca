
! -*-f90-*-
! $Id: mpp_io_write.inc,v 18.0 2010/03/02 23:57:44 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!                             MPP_WRITE_META                                 !
!                                                                            !
!This series of routines is used to describe the contents of the file        !
!being written on <unit>. Each file can contain any number of fields,        !
!which can be functions of 0-3 spatial axes and 0-1 time axes. Axis          !
!descriptors are stored in the <axistype> structure and field                !
!descriptors in the <fieldtype> structure.                                   !
!                                                                            !
!  type, public :: axistype                                                  !
!     sequence                                                               !
!     character(len=128) :: name                                             !
!     character(len=128) :: units                                            !
!     character(len=256) :: longname                                         !
!     integer :: sense           !+/-1, depth or height?                     !
!     type(domain1D) :: domain                                      !
!     real, pointer :: data(:) !axis values (not used if time axis)          !
!     integer :: id                                                          !
!  end type axistype                                                         !
!                                                                            !
!  type, public :: fieldtype                                                 !
!     sequence                                                               !
!     character(len=128) :: name                                             !
!     character(len=128) :: units                                            !
!     character(len=256) :: longname                                         !
!     character(len=128) :: standard_name  !CF standard name                 !
!     real :: min, max, missing, fill, scale, add                            !
!     type(axistype), pointer :: axis(:)                                     !
!     integer :: id                                                          !
!  end type fieldtype                                                        !
!                                                                            !
!The metadata contained in the type is always written for each axis and      !
!field. Any other metadata one wishes to attach to an axis or field          !
!can subsequently be passed to mpp_write_meta using the ID, as shown below.  !
!                                                                            !
!mpp_write_meta can take several forms:                                      !
!                                                                            !
!  mpp_write_meta( unit, name, rval=rval, pack=pack )                        !
!  mpp_write_meta( unit, name, ival=ival )                                   !
!  mpp_write_meta( unit, name, cval=cval )                                   !
!      integer, intent(in) :: unit                                           !
!      character(len=*), intent(in) :: name                                  !
!      real, intent(in), optional :: rval(:)                                 !
!      integer, intent(in), optional :: ival(:)                              !
!      character(len=*), intent(in), optional :: cval                        !
!                                                                            !
!    This form defines global metadata associated with the file as a         !
!    whole. The attribute is named <name> and can take on a real, integer    !
!    or character value. <rval> and <ival> can be scalar or 1D arrays.       !
!                                                                            !
!  mpp_write_meta( unit, id, name, rval=rval, pack=pack )                    !
!  mpp_write_meta( unit, id, name, ival=ival )                               !
!  mpp_write_meta( unit, id, name, cval=cval )                               !
!      integer, intent(in) :: unit, id                                       !
!      character(len=*), intent(in) :: name                                  !
!      real, intent(in), optional :: rval(:)                                 !
!      integer, intent(in), optional :: ival(:)                              !
!      character(len=*), intent(in), optional :: cval                        !
!                                                                            !
!    This form defines metadata associated with a previously defined         !
!    axis or field, identified to mpp_write_meta by its unique ID <id>.      !
!    The attribute is named <name> and can take on a real, integer           !
!    or character value. <rval> and <ival> can be scalar or 1D arrays.       !
!    This need not be called for attributes already contained in             !
!    the type.                                                               !
!                                                                            !
!    PACK can take values 1,2,4,8. This only has meaning when writing        !
!    floating point numbers. The value of PACK defines the number of words   !
!    written into 8 bytes. For pack=4 and pack=8, an integer value is        !
!    written: rval is assumed to have been scaled to the appropriate dynamic !
!    range.                                                                  !
!    PACK currently only works for netCDF files, and is ignored otherwise.   !
!                                                                            !
!   subroutine mpp_write_meta_axis( unit, axis, name, units, longname, &     !
!        cartesian, sense, domain, data )                                    !
!     integer, intent(in) :: unit                                            !
!     type(axistype), intent(inout) :: axis                                  !
!     character(len=*), intent(in) :: name, units, longname                  !
!     character(len=*), intent(in), optional :: cartesian                    !
!     integer, intent(in), optional :: sense                                 !
!     type(domain1D), intent(in), optional :: domain                 !
!     real, intent(in), optional :: data(:)                                  !
!                                                                            !
!    This form defines a time or space axis. Metadata corresponding to the   !
!    type above are written to the file on <unit>. A unique ID for subsequent!
!    references to this axis is returned in axis%id. If the <domain>         !
!    element is present, this is recognized as a distributed data axis       !
!    and domain decomposition information is also written if required (the   !
!    domain decomposition info is required for multi-fileset multi-threaded  !
!    I/O). If the <datLINK> element is allocated, it is considered to be a space!
!    axis, otherwise it is a time axis with an unlimited dimension. Only one !
!    time axis is allowed per file.                                          !
!                                                                            !
!   subroutine mpp_write_meta_field( unit, field, axes, name, units, longname!
!       stanadard_name, min, max, missing, fill, scale, add, pack)           !
!     integer, intent(in) :: unit                                            !
!     type(fieldtype), intent(out) :: field                                  !
!     type(axistype), intent(in) :: axes(:)                                  !
!     character(len=*), intent(in) :: name, units, longname,standard_name    !
!     real, intent(in), optional :: min, max, missing, fill, scale, add      !
!     integer, intent(in), optional :: pack                                  !
!                                                                            !
!    This form defines a field. Metadata corresponding to the type           !
!    above are written to the file on <unit>. A unique ID for subsequent     !
!    references to this field is returned in field%id. At least one axis     !
!    must be associated, 0D variables are not considered. mpp_write_meta     !
!    must previously have been called on all axes associated with this       !
!    field.                                                                  !
!                                                                            !
! The mpp_write_meta package also includes subroutines write_attribute and   !
! write_attribute_netcdf, that are private to this module.                   !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_write_meta_global( unit, name, rval, ival, cval, pack)
!writes a global metadata attribute to unit <unit>
!attribute <name> can be an real, integer or character
!one and only one of rval, ival, and cval should be present
!the first found will be used
!for a non-netCDF file, it is encoded into a string "GLOBAL <name> <val>"
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real,             intent(in), optional :: rval(:)
      integer,          intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
      integer, intent(in), optional :: pack

!      call mpp_clock_begin(mpp_write_clock)
      if( .NOT.module_is_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT. mpp_file(unit)%write_on_this_pe) then
!         call mpp_clock_end(mpp_write_clock)
         return
      endif
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )
           
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          call write_attribute_netcdf( unit, NF_GLOBAL, name, rval, ival, cval, pack )
#endif    
      else
          call write_attribute( unit, 'GLOBAL '//trim(name), rval, ival, cval, pack )
      end if
!      call mpp_clock_end(mpp_write_clock)
          
      return
    end subroutine mpp_write_meta_global
      
!versions of above to support <rval> and <ival> as scalars (because of f90 strict rank matching)
    subroutine mpp_write_meta_global_scalar_r( unit, name, rval, pack )
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real, intent(in) :: rval
      integer, intent(in), optional :: pack
      
      call mpp_write_meta_global( unit, name, rval=(/rval/), pack=pack )
      return
    end subroutine mpp_write_meta_global_scalar_r
      
    subroutine mpp_write_meta_global_scalar_i( unit, name, ival )
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      integer, intent(in) :: ival
      
      call mpp_write_meta_global( unit, name, ival=(/ival/) )
      return
    end subroutine mpp_write_meta_global_scalar_i
      
    subroutine mpp_write_meta_var( unit, id, name, rval, ival, cval, pack)
!writes a metadata attribute for variable <id> to unit <unit>
!attribute <name> can be an real, integer or character
!one and only one of rval, ival, and cval should be present
!the first found will be used
!for a non-netCDF file, it is encoded into a string "<id> <name> <val>"
      integer, intent(in) :: unit, id
      character(len=*), intent(in) :: name
      real,             intent(in), optional :: rval(:)
      integer,          intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
      integer,          intent(in), optional :: pack

      if( .NOT.module_is_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT. mpp_file(unit)%write_on_this_pe) then
         return
      endif
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )
           
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
          call write_attribute_netcdf( unit, id, name, rval, ival, cval, pack )
      else
          write( text, '(a,i4,a)' )'VARIABLE ', id, ' '//name
          call write_attribute( unit, trim(text), rval, ival, cval, pack )
      end if

      return
    end subroutine mpp_write_meta_var
      
!versions of above to support <rval> and <ival> as scalar (because of f90 strict rank matching)
    subroutine mpp_write_meta_scalar_r( unit, id, name, rval, pack )
      integer, intent(in) :: unit, id
      character(len=*), intent(in) :: name
      real, intent(in) :: rval
      integer, intent(in), optional :: pack
      
      call mpp_write_meta( unit, id, name, rval=(/rval/), pack=pack )
      return
    end subroutine mpp_write_meta_scalar_r
      
    subroutine mpp_write_meta_scalar_i( unit, id, name, ival )
      integer, intent(in) :: unit, id
      character(len=*), intent(in) :: name
      integer, intent(in) :: ival
      
      call mpp_write_meta( unit, id, name, ival=(/ival/) )
      return
    end subroutine mpp_write_meta_scalar_i
      
    subroutine mpp_write_meta_axis( unit, axis, name, units, longname, cartesian, sense, domain, data)
!load the values in an axistype (still need to call mpp_write)
!write metadata attributes for axis
!it is declared intent(inout) so you can nullify pointers in the incoming object if needed
!the f90 standard doesn't guarantee that intent(out) on a type guarantees that its pointer components will be unassociated
      integer, intent(in) :: unit
      type(axistype), intent(inout) :: axis
      character(len=*), intent(in) :: name, units, longname
      character(len=*), intent(in), optional :: cartesian
      integer, intent(in), optional :: sense
      type(domain1D), intent(in), optional :: domain
      real, intent(in), optional :: data(:)
      integer :: is, ie, isg, ieg
      integer :: istat      
      logical :: domain_exist
      type(domain2d), pointer :: io_domain => NULL()

!      call mpp_clock_begin(mpp_write_clock)
      !--- the shift and cartesian information is needed in mpp_write_meta_field from all the pe.
      !--- we may revise this in the future.
      axis%cartesian = 'N'
      if( PRESENT(cartesian) )axis%cartesian = cartesian
           
      domain_exist = .false.

     if( PRESENT(domain) ) then
         domain_exist = .true.
         call mpp_get_global_domain( domain, isg, ieg )
         if(mpp_file(unit)%io_domain_exist) then
            io_domain => mpp_get_io_domain(mpp_file(unit)%domain)
            if(axis%cartesian=='X') then
               call mpp_get_global_domain( io_domain, xbegin=is, xend=ie)
            else if(axis%cartesian=='Y') then                        
               call mpp_get_global_domain( io_domain, ybegin=is, yend=ie)
            endif
         else
            call mpp_get_compute_domain( domain, is, ie )
         endif
      else if( PRESENT(data) )then
         isg=1; ieg=size(data(:)); is=isg; ie=ieg
      endif

      axis%shift = 0         
      if( PRESENT(data) .AND. domain_exist ) then
         if( size(data(:)) == ieg-isg+2 ) then
            axis%shift = 1
            ie = ie   + 1
            ieg = ieg + 1
         endif
      endif

      if( .NOT.module_is_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT. mpp_file(unit)%write_on_this_pe) then
!         call mpp_clock_end(mpp_write_clock)
         return
      endif
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )

!pre-existing pointers need to be nullified
      if( ASSOCIATED(axis%data) ) then
         DEALLOCATE(axis%data, stat=istat)
      endif
!load axistype
      axis%name     = name
      axis%units    = units
      axis%longname = longname
      if( PRESENT(sense)     )axis%sense     = sense
      if( PRESENT(data) )then
         if( mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. domain_exist ) then
            axis%len = ie - is + 1
            allocate(axis%data(axis%len))
            axis%data = data(is-isg+1:ie-isg+1)
         else
            axis%len = size(data(:))
            allocate(axis%data(axis%len))
            axis%data = data
         endif
      endif
!write metadata
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
!write axis def
!space axes are always floats, time axis is always double
          if( ASSOCIATED(axis%data) )then !space axis
              error = NF_DEF_DIM( mpp_file(unit)%ncid, axis%name, axis%len, axis%did )
              call netcdf_err( error, mpp_file(unit), axis )
              error = NF_DEF_VAR( mpp_file(unit)%ncid, axis%name, NF_DOUBLE, 1, axis%did, axis%id )
              call netcdf_err( error, mpp_file(unit), axis )
          else                            !time axis
              if( mpp_file(unit)%id.NE.-1 ) &
                   call mpp_error( FATAL, 'MPP_WRITE_META_AXIS: There is already a time axis for this file.' )
              error = NF_DEF_DIM( mpp_file(unit)%ncid, axis%name, NF_UNLIMITED, axis%did )
              call netcdf_err( error, mpp_file(unit), axis )
              error = NF_DEF_VAR( mpp_file(unit)%ncid, axis%name, NF_DOUBLE, 1, axis%did, axis%id )
              call netcdf_err( error, mpp_file(unit), axis )
              mpp_file(unit)%id = axis%id !file ID is the same as time axis varID
          end if
#endif        
      else
          varnum = varnum + 1
          axis%id = varnum
          axis%did = varnum
!write axis def
          write( text, '(a,i4,a)' )'AXIS ', axis%id, ' name'
          call write_attribute( unit, trim(text), cval=axis%name )
          write( text, '(a,i4,a)' )'AXIS ', axis%id, ' size'
          if( ASSOCIATED(axis%data) )then !space axis
!              if( mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. axis%domain.NE.NULL_DOMAIN1D )then
!                  call write_attribute( unit, trim(text), ival=(/ie-is+1/) )
!              else
                  call write_attribute( unit, trim(text), ival=(/size(axis%data(:))/) )
!              end if
          else                            !time axis
              if( mpp_file(unit)%id.NE.-1 ) &
                   call mpp_error( FATAL, 'MPP_WRITE_META_AXIS: There is already a time axis for this file.' )
              call write_attribute( unit, trim(text), ival=(/0/) ) !a size of 0 indicates time axis
              mpp_file(unit)%id = axis%id
          end if
      end if  
!write axis attributes
      call mpp_write_meta( unit, axis%id, 'long_name', cval=axis%longname)
      call mpp_write_meta( unit, axis%id, 'units',     cval=axis%units)
      if( PRESENT(cartesian) )call mpp_write_meta( unit, axis%id, 'cartesian_axis', cval=axis%cartesian)
      if( PRESENT(sense) )then
          if( sense.EQ.-1 )then
              call mpp_write_meta( unit, axis%id, 'positive', cval='down')
          else if( sense.EQ.1 )then
              call mpp_write_meta( unit, axis%id, 'positive', cval='up')
          end if
!silently ignore values of sense other than +/-1.
      end if
      if( mpp_file(unit)%threading.EQ.MPP_MULTI .AND. mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. domain_exist )then
          call mpp_write_meta( unit, axis%id, 'domain_decomposition', ival=(/isg,ieg,is,ie/))
      end if
      if( verbose )print '(a,2i3,x,a,2i3)', 'MPP_WRITE_META: Wrote axis metadata, pe, unit, axis%name, axis%id, axis%did=', &
           pe, unit, trim(axis%name), axis%id, axis%did

!      call mpp_clock_end(mpp_write_clock)           
      return
    end subroutine mpp_write_meta_axis
      
    subroutine mpp_write_meta_field( unit, field, axes, name, units, longname,&
         min, max, missing, fill, scale, add, pack, time_method,standard_name)
!define field: must have already called mpp_write_meta(axis) for each axis
      integer, intent(in) :: unit
      type(fieldtype), intent(inout) :: field
      type(axistype), intent(in) :: axes(:)
      character(len=*), intent(in) :: name, units, longname
      real, intent(in), optional :: min, max, missing, fill, scale, add
      integer, intent(in), optional :: pack
      character(len=*), intent(in), optional :: time_method
      character(len=*), intent(in), optional :: standard_name
!this array is required because of f77 binding on netCDF interface
      integer, allocatable :: axis_id(:)
      real :: a, b
      integer :: i, istat, ishift, jshift

!      call mpp_clock_begin(mpp_write_clock)

      !--- figure out the location of data, this is needed in mpp_write.
      !--- for NON-symmetry domain, the position is not an issue.
      !--- we may need to rethink how to address the symmetric issue.
      ishift = 0; jshift = 0
      do i = 1, size(axes(:))
         select case ( lowercase( axes(i)%cartesian ) )
         case ( 'x' )
            ishift = axes(i)%shift
         case ( 'y' )
            jshift = axes(i)%shift 
         end select
      end do

      field%position = CENTER
      if(ishift == 1 .AND. jshift == 1) then
         field%position = CORNER
      else if(ishift == 1) then
         field%position = EAST
      else if(jshift == 1) then
         field%position = NORTH
      endif
      
      if( .NOT.module_is_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )

      if( .NOT.mpp_file(unit)%write_on_this_pe) then
         if( .NOT. ASSOCIATED(field%axes) )allocate(field%axes(1)) !temporary fix
!         call mpp_clock_end(mpp_write_clock)
         return
      endif
      if( .NOT.mpp_file(unit)%opened ) call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )
           
!pre-existing pointers need to be nullified
      if( ASSOCIATED(field%axes) ) then
         DEALLOCATE(field%axes, stat=istat)
      endif
!fill in field metadata
      field%name = name
      field%units = units
      field%longname = longname
      allocate( field%axes(size(axes(:))) )
      field%axes = axes
      field%ndim = size(axes(:))
      field%time_axis_index = -1 !this value will never match any axis index
!size is buffer area for the corresponding axis info: it is required to buffer this info in the fieldtype
!because axis might be reused in different files
      allocate( field%size(size(axes(:))) )
      do i = 1,size(axes(:))
         if( ASSOCIATED(axes(i)%data) )then !space axis
             field%size(i) = size(axes(i)%data(:))
         else               !time
             field%size(i) = 1
             field%time_axis_index = i
         end if
      end do 
!attributes
      if( PRESENT(min) )field%min = min
      if( PRESENT(max) )field%max = max
      if( PRESENT(missing) )field%missing = missing
      if( PRESENT(fill) )field%fill = fill
      if( PRESENT(scale) )field%scale = scale
      if( PRESENT(add) )field%add = add
      if( PRESENT(standard_name)) field%standard_name = standard_name
      
!pack is currently used only for netCDF
      field%pack = 2        !default write 32-bit floats
      if( PRESENT(pack) )field%pack = pack
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          allocate( axis_id(size(field%axes(:))) )
          do i = 1,size(field%axes(:))
             axis_id(i) = field%axes(i)%did
          end do
!write field def
          select case (field%pack)
              case(1)
                  error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_DOUBLE, size(field%axes(:)), axis_id, field%id )
              case(2)
                  error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_FLOAT,  size(field%axes(:)), axis_id, field%id )
              case(4)
                  if( .NOT.PRESENT(scale) .OR. .NOT.PRESENT(add) ) &
                       call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: scale and add must be supplied when pack=4.' )
                  error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_SHORT,  size(field%axes(:)), axis_id, field%id )
              case(8)  
                  if( .NOT.PRESENT(scale) .OR. .NOT.PRESENT(add) ) &
                       call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: scale and add must be supplied when pack=8.' )
                  error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_BYTE,   size(field%axes(:)), axis_id, field%id )
              case default
                  call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: only legal packing values are 1,2,4,8.' )
          end select
          call netcdf_err( error, mpp_file(unit), field=field )
#endif    
      else
          varnum = varnum + 1
          field%id = varnum
          if( PRESENT(pack) )call mpp_error( WARNING, 'MPP_WRITE_META: Packing is currently available only on netCDF files.' )
!write field def
          write( text, '(a,i4,a)' )'FIELD ', field%id, ' name'
          call write_attribute( unit, trim(text), cval=field%name )
          write( text, '(a,i4,a)' )'FIELD ', field%id, ' axes'
          call write_attribute( unit, trim(text), ival=field%axes(:)%did )
      end if
!write field attributes: these names follow netCDF conventions
      call mpp_write_meta( unit, field%id, 'long_name', cval=field%longname)
      call mpp_write_meta( unit, field%id, 'units',     cval=field%units)
!all real attributes must be written as packed
      if( PRESENT(min) .AND. PRESENT(max) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_range', rval=(/min,max/), pack=pack )
          else
              a = nint((min-add)/scale)
              b = nint((max-add)/scale)
              call mpp_write_meta( unit, field%id, 'valid_range', rval=(/a,  b  /), pack=pack )
          end if
      else if( PRESENT(min) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_min', rval=field%min, pack=pack )
          else
              a = nint((min-add)/scale)
              call mpp_write_meta( unit, field%id, 'valid_min', rval=a, pack=pack )
          end if
      else if( PRESENT(max) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_max', rval=field%max, pack=pack )
          else
              a = nint((max-add)/scale)
              call mpp_write_meta( unit, field%id, 'valid_max', rval=a, pack=pack )
          end if
      end if  
      if( PRESENT(missing) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'missing_value', rval=field%missing, pack=pack )
          else
              a = nint((missing-add)/scale)
              call mpp_write_meta( unit, field%id, 'missing_value', rval=a, pack=pack )
          end if
      end if  
      if( PRESENT(fill) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, '_FillValue', rval=field%missing, pack=pack )
          else
              a = nint((fill-add)/scale)
              call mpp_write_meta( unit, field%id, '_FillValue', rval=a, pack=pack )
          end if
      end if  
      if( field%pack.NE.1 .AND. field%pack.NE.2 )then
          call mpp_write_meta( unit, field%id, 'packing', ival=field%pack )
          if( PRESENT(scale) )call mpp_write_meta( unit, field%id, 'scale_factor',  rval=field%scale )
          if( PRESENT(add)   )call mpp_write_meta( unit, field%id, 'add_offset',    rval=field%add   )
      end if

      if ( PRESENT(time_method) ) then
          call mpp_write_meta(unit,field%id, 'cell_methods',cval='time: '//trim(time_method))
      endif
      if ( PRESENT(standard_name)) &
           call mpp_write_meta(unit,field%id,'standard_name ', cval=field%standard_name) 

      if( verbose )print '(a,2i3,x,a,i3)', 'MPP_WRITE_META: Wrote field metadata: pe, unit, field%name, field%id=', &
           pe, unit, trim(field%name), field%id
          
!      call mpp_clock_end(mpp_write_clock) 
      return
    end subroutine mpp_write_meta_field
      
    subroutine write_attribute( unit, name, rval, ival, cval, pack )
!called to write metadata for non-netCDF I/O
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real, intent(in), optional :: rval(:)
      integer, intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
!pack is currently ignored in this routine: only used by netCDF I/O
      integer, intent(in), optional :: pack
      
      if( mpp_file(unit)%nohdrs )return
!encode text string
      if( PRESENT(rval) )then
          write( text,* )trim(name)//'=', rval
      else if( PRESENT(ival) )then
          write( text,* )trim(name)//'=', ival
      else if( PRESENT(cval) )then
          text = ' '//trim(name)//'='//trim(cval)
      else
          call mpp_error( FATAL, 'WRITE_ATTRIBUTE: one of rval, ival, cval must be present.' )
      end if
      if( mpp_file(unit)%format.EQ.MPP_ASCII )then
!implies sequential access
          write( unit,fmt='(a)' )trim(text)//char(10)
      else                      !MPP_IEEE32 or MPP_NATIVE
          if( mpp_file(unit)%access.EQ.MPP_SEQUENTIAL )then
              write(unit)trim(text)//char(10)
          else                  !MPP_DIRECT
              write( unit,rec=mpp_file(unit)%record )trim(text)//char(10)
              if( verbose )print '(a,i3,a,i3)', 'WRITE_ATTRIBUTE: PE=', pe, ' wrote record ', mpp_file(unit)%record
              mpp_file(unit)%record = mpp_file(unit)%record + 1
          end if
      end if  
      return
    end subroutine write_attribute
      
    subroutine write_attribute_netcdf( unit, id, name, rval, ival, cval, pack )
!called to write metadata for netCDF I/O
      integer, intent(in) :: unit
      integer, intent(in) :: id
      character(len=*), intent(in) :: name
      real,             intent(in), optional :: rval(:)
      integer,          intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
      integer, intent(in), optional :: pack
      integer, allocatable :: rval_i(:)
#ifdef use_netCDF
      if( PRESENT(rval) )then
!pack is only meaningful for FP numbers
          if( PRESENT(pack) )then
              if( pack.EQ.1 )then
                  if( KIND(rval).EQ.DOUBLE_KIND )then
                      error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_DOUBLE, size(rval(:)), rval )
                  else if( KIND(rval).EQ.FLOAT_KIND )then
                      call mpp_error( WARNING, &
                           'WRITE_ATTRIBUTE_NETCDF: attempting to write internal 32-bit real as external 64-bit.' )
                      error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_DOUBLE, size(rval(:)), rval )
                  end if   
                  call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
              else if( pack.EQ.2 )then
                  if( KIND(rval).EQ.DOUBLE_KIND )then
                      error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_FLOAT,  size(rval(:)), rval )
                  else if( KIND(rval).EQ.FLOAT_KIND )then
                      error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_FLOAT,  size(rval(:)), rval )
                  end if
                  call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
              else if( pack.EQ.4 )then
                  allocate( rval_i(size(rval(:))) )
                  rval_i = rval
                  if( KIND(rval).EQ.DOUBLE_KIND )then
                      error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_SHORT,  size(rval_i(:)), rval )
                  else if( KIND(rval).EQ.FLOAT_KIND )then
                      error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_SHORT,  size(rval_i(:)), rval )
                  end if
                  call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
                  deallocate(rval_i)
              else if( pack.EQ.8 )then
                  allocate( rval_i(size(rval(:))) )
                  rval_i = rval
                  if( KIND(rval).EQ.DOUBLE_KIND )then
                      error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_BYTE,   size(rval_i(:)), rval )
                  else if( KIND(rval).EQ.FLOAT_KIND )then
                      error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_BYTE,   size(rval_i(:)), rval )
                  end if
                  call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
                  deallocate(rval_i)
              else
                  call mpp_error( FATAL, 'WRITE_ATTRIBUTE_NETCDF: only legal packing values are 1,2,4,8.' )
              end if
          else    
!default is to write FLOATs (32-bit)
              if( KIND(rval).EQ.DOUBLE_KIND )then
                  error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_FLOAT,  size(rval(:)), rval )
              else if( KIND(rval).EQ.FLOAT_KIND )then
                  error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_FLOAT,  size(rval(:)), rval )
              end if
              call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
          end if
      else if( PRESENT(ival) )then
          error = NF_PUT_ATT_INT ( mpp_file(unit)%ncid, id, name, NF_INT, size(ival(:)), ival )
          call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
      else if( present(cval) )then
          error = NF_PUT_ATT_TEXT( mpp_file(unit)%ncid, id, name, len_trim(cval), cval )
          call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
      else
          call mpp_error( FATAL, 'WRITE_ATTRIBUTE_NETCDF: one of rval, ival, cval must be present.' )
      end if
#endif /* use_netCDF */
      return
    end subroutine write_attribute_netcdf
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                             MPP_WRITE                                !
!                                                                      !
! mpp_write is used to write data to the file on <unit> using the      !
! file parameters supplied by mpp_open(). Axis and field definitions   !
! must have previously been written to the file using mpp_write_meta.  !
!                                                                      !
! mpp_write can take 2 forms, one for distributed data and one for     !
! non-distributed data. Distributed data refer to arrays whose two     !
! fastest-varying indices are domain-decomposed. Distributed data      !
! must be 2D or 3D (in space). Non-distributed data can be 0-3D.       !
!                                                                      !
! In all calls to mpp_write, tstamp is an optional argument. It is to  !
! be omitted if the field was defined not to be a function of time.    !
! Results are unpredictable if the argument is supplied for a time-    !
! independent field, or omitted for a time-dependent field. Repeated   !
! writes of a time-independent field are also not recommended. One     !
! time level of one field is written per call.                         !
!                                                                      !
!                                                                      !
! For non-distributed data, use                                        !
!                                                                      !
!  mpp_write( unit, field, data, tstamp )                              !
!     integer, intent(in) :: unit                                      !
!     type(fieldtype), intent(in) :: field                             !
!     real(DOUBLE_KIND), optional :: tstamp                                         !
!     data is real and can be scalar or of rank 1-3.                   !
!                                                                      !
! For distributed data, use                                            !
!                                                                      !
!  mpp_write( unit, field, domain, data, tstamp )                      !
!     integer, intent(in) :: unit                                      !
!     type(fieldtype), intent(in) :: field                             !
!     type(domain2D), intent(in) :: domain                             !
!     real(DOUBLE_KIND), optional :: tstamp                                         !
!     data is real and can be of rank 2 or 3.                          !
!                                                                      !
!  mpp_write( unit, axis )                                             !
!     integer, intent(in) :: unit                                      !
!     type(axistype), intent(in) :: axis                               !
!                                                                      !
! This call writes the actual co-ordinate values along each space      !
! axis. It must be called once for each space axis after all other     !
! metadata has been written.                                           !
!                                                                      !
! The mpp_write package also includes the routine write_record which   !
! performs the actual write. This routine is private to this module.   !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef MPP_WRITE_2DDECOMP_2D_
#define MPP_WRITE_2DDECOMP_2D_ mpp_write_2ddecomp_r2d
#undef MPP_WRITE_2DDECOMP_3D_
#define MPP_WRITE_2DDECOMP_3D_ mpp_write_2ddecomp_r3d
#undef MPP_WRITE_2DDECOMP_4D_
#define MPP_WRITE_2DDECOMP_4D_ mpp_write_2ddecomp_r4d
#undef MPP_TYPE_
#define MPP_TYPE_ real
#include <mpp_write_2Ddecomp.h>

#undef MPP_WRITE_
#define MPP_WRITE_ mpp_write_r0D
#undef MPP_TYPE_
#define MPP_TYPE_ real
#undef MPP_RANK_
#define MPP_RANK_ !
#undef MPP_WRITE_RECORD_
#define MPP_WRITE_RECORD_ call write_record( unit, field, 1, (/data/), tstamp)
#include <mpp_write.h>

#undef MPP_WRITE_
#define MPP_WRITE_ mpp_write_r1D
#undef MPP_TYPE_
#define MPP_TYPE_ real
#undef MPP_WRITE_RECORD_
#define MPP_WRITE_RECORD_ call write_record( unit, field, size(data(:)), data, tstamp)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_write.h>

#undef MPP_WRITE_
#define MPP_WRITE_ mpp_write_r2D
#undef MPP_TYPE_
#define MPP_TYPE_ real
#undef MPP_WRITE_RECORD_
#define MPP_WRITE_RECORD_ call write_record( unit, field, size(data(:,:)), data, tstamp )
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_write.h>

#undef MPP_WRITE_
#define MPP_WRITE_ mpp_write_r3D
#undef MPP_TYPE_
#define MPP_TYPE_ real
#undef MPP_WRITE_RECORD_
#define MPP_WRITE_RECORD_ call write_record( unit, field, size(data(:,:,:)), data, tstamp)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_write.h>

#undef MPP_WRITE_
#define MPP_WRITE_ mpp_write_r4D
#undef MPP_TYPE_
#define MPP_TYPE_ real
#undef MPP_WRITE_RECORD_
#define MPP_WRITE_RECORD_ call write_record( unit, field, size(data(:,:,:,:)), data, tstamp)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_write.h>

    subroutine mpp_write_axis( unit, axis )
      integer, intent(in) :: unit
      type(axistype), intent(in) :: axis
      type(fieldtype) :: field

      call mpp_clock_begin(mpp_write_clock)
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT. mpp_file(unit)%write_on_this_pe ) then
         call mpp_clock_end(mpp_write_clock)
         return
      endif
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )
!we convert axis to type(fieldtype) in order to call write_record
      field = default_field
      allocate( field%axes(1) )
      field%axes(1) = axis
      allocate( field%size(1) )
      field%size(1) = axis%len
      field%id = axis%id

      call write_record( unit, field, axis%len, axis%data )
      call mpp_clock_end(mpp_write_clock)
      return
    end subroutine mpp_write_axis
      
    subroutine write_record( unit, field, nwords, data, time_in, domain, tile_count)
!routine that is finally called by all mpp_write routines to perform the write
!a non-netCDF record contains:
!      field ID
!      a set of 4 coordinates (is:ie,js:je) giving the data subdomain
!      a timelevel and a timestamp (=NULLTIME if field is static)
!      3D real data (stored as 1D)
!if you are using direct access I/O, the RECL argument to OPEN must be large enough for the above
!in a global direct access file, record position on PE is given by %record.

!Treatment of timestamp:
!   We assume that static fields have been passed without a timestamp.
!   Here that is converted into a timestamp of NULLTIME.
!   For non-netCDF fields, field is treated no differently, but is written
!   with a timestamp of NULLTIME. There is no check in the code to prevent
!   the user from repeatedly writing a static field.

      integer,           intent(in)           :: unit, nwords
      type(fieldtype),   intent(in)           :: field
      real,              intent(in)           :: data(nwords)
      real(DOUBLE_KIND), intent(in), optional :: time_in
      type(domain2D),    intent(in), optional :: domain
      integer,           intent(in), optional :: tile_count
      integer, dimension(size(field%axes(:))) :: start, axsiz
      real(DOUBLE_KIND) :: time
      integer :: time_level
      logical :: newtime
      integer :: subdomain(4)
      integer :: packed_data(nwords)
      integer :: i, is, ie, js, je

      real(FLOAT_KIND) :: data_r4(nwords)
      pointer( ptr1, data_r4)
      pointer( ptr2, packed_data)
      
      if (mpp_io_stack_size < nwords) call mpp_io_set_stack_size(nwords)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%write_on_this_pe) return
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )
      if( .NOT.mpp_file(unit)%initialized )then
!this is the first call to mpp_write
!we now declare the file to be initialized
!if this is netCDF we switch file from DEFINE mode to DATA mode
          if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
!NOFILL is probably required for parallel: any circumstances in which not advisable?
              error = NF_SET_FILL( mpp_file(unit)%ncid, NF_NOFILL, i ); call netcdf_err( error, mpp_file(unit) )
              if( mpp_file(unit)%action.EQ.MPP_WRONLY )then
                 if(header_buffer_val>0) then
                    error = NF__ENDDEF(mpp_file(unit)%ncid,header_buffer_val,4,0,4)
                 else
                    error = NF_ENDDEF(mpp_file(unit)%ncid)
                 endif
              endif
              call netcdf_err( error, mpp_file(unit) )
#endif        
          else
              call mpp_write_meta( unit, 'END', cval='metadata' )
          end if
          mpp_file(unit)%initialized = .TRUE.
          if( verbose )print '(a,i3,a)', 'MPP_WRITE: PE=', pe, ' initialized file '//trim(mpp_file(unit)%name)//'.'
      end if
          
!initialize time: by default assume NULLTIME
      time = NULLTIME
      time_level = -1
      newtime = .FALSE.
      if( PRESENT(time_in) )time = time_in
!increment time level if new time
      if( time.GT.mpp_file(unit)%time+EPSILON(time) )then !new time
          mpp_file(unit)%time_level = mpp_file(unit)%time_level + 1
          mpp_file(unit)%time = time
          newtime = .TRUE.
      end if
      if( verbose )print '(a,2i3,2i5,es13.5)', 'MPP_WRITE: PE, unit, %id, %time_level, %time=',&
           pe, unit, mpp_file(unit)%id, mpp_file(unit)%time_level, mpp_file(unit)%time
           
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
          ptr2 = LOC(mpp_io_stack(1))
!define netCDF data block to be written:
!  time axis: START = time level
!             AXSIZ = 1
!  space axis: if there is no domain info
!              START = 1
!              AXSIZ = field%size(axis)
!          if there IS domain info:
!              start of domain is compute%start_index for multi-file I/O
!                                 global%start_index for all other cases
!              this number must be converted to 1 for NF_PUT_VAR
!                  (netCDF fortran calls are with reference to 1),
!          So, START = compute%start_index - <start of domain> + 1
!              AXSIZ = usually compute%size
!          However, if compute%start_index-compute%end_index+1.NE.compute%size,
!              we assume that the call is passing a subdomain.
!              To pass a subdomain, you must pass a domain2D object that satisfies the following:
!                  global%start_index must contain the <start of domain> as defined above;
!                  the data domain and compute domain must refer to the subdomain being passed.
!              In this case, START = compute%start_index - <start of domain> + 1
!                            AXSIZ = compute%start_index - compute%end_index + 1
! NOTE: passing of subdomains will fail for multi-PE single-threaded I/O,
!       since that attempts to gather all data on PE 0.
          start = 1
          do i = 1,size(field%axes(:))
             axsiz(i) = field%size(i)
             if( i.EQ.field%time_axis_index )start(i) = mpp_file(unit)%time_level
             start(i) = max(start(i),1)
          end do

          if( debug )print '(a,2i6,12i6)', 'WRITE_RECORD: PE, unit, start, axsiz=', pe, unit, start, axsiz
#ifdef use_netCDF
!write time information if new time
          if( newtime )then
              if( KIND(time).EQ.DOUBLE_KIND )then
                  error = NF_PUT_VAR1_DOUBLE( mpp_file(unit)%ncid, mpp_file(unit)%id, mpp_file(unit)%time_level, time )
              else if( KIND(time).EQ.FLOAT_KIND )then
                  error = NF_PUT_VAR1_REAL  ( mpp_file(unit)%ncid, mpp_file(unit)%id, mpp_file(unit)%time_level, time )
              end if
          end if  
          if( field%pack.LE.2 )then
              if( KIND(data).EQ.DOUBLE_KIND )then
                  error = NF_PUT_VARA_DOUBLE( mpp_file(unit)%ncid, field%id, start, axsiz, data )
              else if( KIND(data).EQ.FLOAT_KIND )then
                  error = NF_PUT_VARA_REAL  ( mpp_file(unit)%ncid, field%id, start, axsiz, data )
              end if
          else              !convert to integer using scale and add: no error check on packed data representation
              packed_data = nint((data-field%add)/field%scale)
              error = NF_PUT_VARA_INT   ( mpp_file(unit)%ncid, field%id, start, axsiz, packed_data )
          end if
          call netcdf_err( error, mpp_file(unit), field=field )
#endif    
      else                      !non-netCDF
          ptr1 = LOC(mpp_io_stack(1))
!subdomain contains (/is,ie,js,je/)  
          if( PRESENT(domain) )then
              subdomain(:) = (/ is, ie, js, je /)  ! ??? is, ie, js, je are never initialized.
          else
              subdomain(:) = -1    ! -1 means use global value from axis metadata
          end if
          if( mpp_file(unit)%format.EQ.MPP_ASCII )then
!implies sequential access
              write( unit,* )field%id, subdomain, time_level, time, data
          else                      !MPP_IEEE32 or MPP_NATIVE
              if( mpp_file(unit)%access.EQ.MPP_SEQUENTIAL )then
#ifdef __sgi  
                  if( mpp_file(unit)%format.EQ.MPP_IEEE32 )then
                      data_r4 = data !IEEE conversion layer on SGI until assign -N ieee_32 is supported
                      write(unit)field%id, subdomain, time_level, time, data_r4
                  else
                      write(unit)field%id, subdomain, time_level, time, data
                  end if
#else                 
                  write(unit)field%id, subdomain, time_level, time, data
#endif            
              else                  !MPP_DIRECT
#ifdef __sgi  
                  if( mpp_file(unit)%format.EQ.MPP_IEEE32 )then
                      data_r4 = data !IEEE conversion layer on SGI until assign -N ieee_32 is supported
                      write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data_r4
                  else
                      write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data
                  end if
#else                 
                  write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data
#endif            
                  if( debug )print '(a,i3,a,i3)', 'MPP_WRITE: PE=', pe, ' wrote record ', mpp_file(unit)%record
              end if
          end if  
      end if  
          
!recompute current record for direct access I/O
      if( mpp_file(unit)%access.EQ.MPP_DIRECT )then
          if( mpp_file(unit)%fileset.EQ.MPP_SINGLE )then
!assumes all PEs participate in I/O: modify later
              mpp_file(unit)%record = mpp_file(unit)%record + records_per_pe*npes
          else
              mpp_file(unit)%record = mpp_file(unit)%record + records_per_pe
          end if
      end if  
          
      return
    end subroutine write_record
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                          MPP_COPY_META                               !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_copy_meta_global( unit, gatt )
!writes a global metadata attribute to unit <unit>
!attribute <name> can be an real, integer or character
!one and only one of rval, ival, and cval should be present
!the first found will be used
!for a non-netCDF file, it is encoded into a string "GLOBAL <name> <val>"
      integer, intent(in) :: unit
      type(atttype), intent(in) :: gatt
      integer :: len
      
      if( .NOT.module_is_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT. mpp_file(unit)%write_on_this_pe )return
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )
#ifdef use_netCDF
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
         if( gatt%type.EQ.NF_CHAR )then
            len = gatt%len
            call write_attribute_netcdf( unit, NF_GLOBAL, gatt%name, cval=gatt%catt(1:len) )
         else
            call write_attribute_netcdf( unit, NF_GLOBAL, gatt%name, rval=gatt%fatt )
         endif
      else  
         if( gatt%type.EQ.NF_CHAR )then
            len=gatt%len
            call write_attribute( unit, 'GLOBAL '//trim(gatt%name), cval=gatt%catt(1:len) )
         else
            call write_attribute( unit, 'GLOBAL '//trim(gatt%name), rval=gatt%fatt )
         endif
     end if 
#else    
     call mpp_error( FATAL, 'MPP_READ currently requires use_netCDF option' )
#endif
      return
    end subroutine mpp_copy_meta_global
      
    subroutine mpp_copy_meta_axis( unit, axis, domain )
!load the values in an axistype (still need to call mpp_write)
!write metadata attributes for axis.  axis is declared inout
!because the variable and dimension ids are altered

      integer, intent(in) :: unit
      type(axistype), intent(inout) :: axis
      type(domain1D), intent(in), optional :: domain
      character(len=512) :: text
      integer :: i, len, is, ie, isg, ieg

!      call mpp_clock_begin(mpp_write_clock)      
      if( .NOT.module_is_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT. mpp_file(unit)%write_on_this_pe ) then 
!         call mpp_clock_end(mpp_write_clock)
         return
      endif
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )
           
! redefine domain if present
      if( PRESENT(domain) )then
          axis%domain = domain
      else
          axis%domain = NULL_DOMAIN1D
      end if
          
#ifdef use_netCDF
!write metadata
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
      
!write axis def
          if( ASSOCIATED(axis%data) )then !space axis
              if( mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. axis%domain.NE.NULL_DOMAIN1D )then
                  call mpp_get_compute_domain( axis%domain, is, ie )
                  call mpp_get_global_domain( axis%domain, isg, ieg )
                  ie  = ie + axis%shift
                  ieg = ieg + axis%shift
                  error = NF_DEF_DIM( mpp_file(unit)%ncid, axis%name, ie-is+1, axis%did )
              else
                  error = NF_DEF_DIM( mpp_file(unit)%ncid, axis%name, size(axis%data(:)),          axis%did )
              end if
              call netcdf_err( error, mpp_file(unit), axis )
              error = NF_DEF_VAR( mpp_file(unit)%ncid, axis%name, NF_FLOAT, 1, axis%did, axis%id )
              call netcdf_err( error, mpp_file(unit), axis )
          else                            !time axis
              error = NF_DEF_DIM( mpp_file(unit)%ncid, axis%name, NF_UNLIMITED, axis%did )
              call netcdf_err( error, mpp_file(unit), axis )
              error = NF_DEF_VAR( mpp_file(unit)%ncid, axis%name, NF_DOUBLE, 1, axis%did, axis%id )
              call netcdf_err( error, mpp_file(unit), axis )
              mpp_file(unit)%id = axis%id !file ID is the same as time axis varID
              mpp_file(unit)%recdimid = axis%did ! record dimension id
          end if
      else    
          varnum = varnum + 1
          axis%id = varnum
          axis%did = varnum
!write axis def
          write( text, '(a,i4,a)' )'AXIS ', axis%id, ' name'
          call write_attribute( unit, trim(text), cval=axis%name )
          write( text, '(a,i4,a)' )'AXIS ', axis%id, ' size'
          if( ASSOCIATED(axis%data) )then !space axis
              if( mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. axis%domain.NE.NULL_DOMAIN1D )then
                  call write_attribute( unit, trim(text), ival=(/ie-is+1/) ) ! ??? is, ie is not initialized 
              else
                  call write_attribute( unit, trim(text), ival=(/size(axis%data(:))/) )
              end if
          else                            !time axis
              if( mpp_file(unit)%id.NE.-1 ) &
                   call mpp_error( FATAL, 'MPP_WRITE_META_AXIS: There is already a time axis for this file.' )
              call write_attribute( unit, trim(text), ival=(/0/) ) !a size of 0 indicates time axis
              mpp_file(unit)%id = axis%id
          end if
      end if  
!write axis attributes
      
      do i=1,axis%natt
         if( axis%Att(i)%name.NE.default_att%name )then
            if( axis%Att(i)%type.EQ.NF_CHAR )then
               len = axis%Att(i)%len
               call mpp_write_meta( unit, axis%id, axis%Att(i)%name, cval=axis%Att(i)%catt(1:len) )
            else
               call mpp_write_meta( unit, axis%id, axis%Att(i)%name, rval=axis%Att(i)%fatt)
            endif
         endif 
      enddo 
         
      if( mpp_file(unit)%threading.EQ.MPP_MULTI .AND. mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. axis%domain.NE.NULL_DOMAIN1D )then
          call mpp_write_meta( unit, axis%id, 'domain_decomposition', ival=(/isg,ieg,is,ie/) )
      end if
      if( verbose )print '(a,2i3,x,a,2i3)', 'MPP_WRITE_META: Wrote axis metadata, pe, unit, axis%name, axis%id, axis%did=', &
           pe, unit, trim(axis%name), axis%id, axis%did
#else      
      call mpp_error( FATAL, 'MPP_READ currently requires use_netCDF option' )
#endif
!      call mpp_clock_end(mpp_write_clock)
      return
    end subroutine mpp_copy_meta_axis
      
    subroutine mpp_copy_meta_field( unit, field, axes )
!useful for copying field metadata from a previous call to mpp_read_meta
!define field: must have already called mpp_write_meta(axis) for each axis
      integer, intent(in) :: unit
      type(fieldtype), intent(inout) :: field
      type(axistype), intent(in), optional :: axes(:)
!this array is required because of f77 binding on netCDF interface
      integer, allocatable :: axis_id(:)
      real :: a, b
      integer :: i
      
!      call mpp_clock_begin(mpp_write_clock)
      if( .NOT.module_is_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT. mpp_file(unit)%write_on_this_pe ) then
!         call mpp_clock_end(mpp_write_clock)
         return
      endif
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )
           
       if( field%pack.NE.1 .AND. field%pack.NE.2 )then
            if( field%pack.NE.4 .AND. field%pack.NE.8 ) &
               call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: only legal packing values are 1,2,4,8.' )
      end if   
               
      if (PRESENT(axes)) then
         deallocate(field%axes)
         deallocate(field%size)
         allocate(field%axes(size(axes(:))))
         allocate(field%size(size(axes(:))))
         field%axes = axes
         do i=1,size(axes(:))
            if (ASSOCIATED(axes(i)%data)) then
               field%size(i) = size(axes(i)%data(:))
            else
               field%size(i) = 1
               field%time_axis_index = i
            endif
         enddo 
      endif 
         
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          allocate( axis_id(size(field%axes(:))) )
          do i = 1,size(field%axes(:))
             axis_id(i) = field%axes(i)%did
          end do
!write field def
          select case (field%pack)
              case(1)
                  error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_DOUBLE, size(field%axes(:)), axis_id, field%id )
              case(2)
                  error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_FLOAT,  size(field%axes(:)), axis_id, field%id )
              case(4)
!                 if( field%scale.EQ.default_field%scale .OR. field%add.EQ.default_field%add ) &
!                      call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: scale and add must be supplied when pack=4.' )
                  error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_SHORT,  size(field%axes(:)), axis_id, field%id )
              case(8)  
!                 if( field%scale.EQ.default_field%scale .OR. field%add.EQ.default_field%add ) &
!                      call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: scale and add must be supplied when pack=8.' )
                  error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_BYTE,   size(field%axes(:)), axis_id, field%id )
              case default
                  call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: only legal packing values are 1,2,4,8.' )
          end select
#endif            
      else
          varnum = varnum + 1
          field%id = varnum
          if( field%pack.NE.default_field%pack ) &
           call mpp_error( WARNING, 'MPP_WRITE_META: Packing is currently available only on netCDF files.' )
!write field def
          write( text, '(a,i4,a)' )'FIELD ', field%id, ' name'
          call write_attribute( unit, trim(text), cval=field%name )
          write( text, '(a,i4,a)' )'FIELD ', field%id, ' axes'
          call write_attribute( unit, trim(text), ival=field%axes(:)%did )
      end if
!write field attributes: these names follow netCDF conventions
      call mpp_write_meta( unit, field%id, 'long_name', cval=field%longname )
      call mpp_write_meta( unit, field%id, 'units',     cval=field%units    )
!all real attributes must be written as packed
      if( (field%min.NE.default_field%min) .AND. (field%max.NE.default_field%max) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_range', rval=(/field%min,field%max/), pack=field%pack )
          else
              a = nint((field%min-field%add)/field%scale)
              b = nint((field%max-field%add)/field%scale)
              call mpp_write_meta( unit, field%id, 'valid_range', rval=(/a,  b  /), pack=field%pack )
          end if
      else if( field%min.NE.default_field%min )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_min', rval=field%min, pack=field%pack )
          else
              a = nint((field%min-field%add)/field%scale)
              call mpp_write_meta( unit, field%id, 'valid_min', rval=a, pack=field%pack )
          end if
      else if( field%max.NE.default_field%max )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_max', rval=field%max, pack=field%pack )
          else
              a = nint((field%max-field%add)/field%scale)
              call mpp_write_meta( unit, field%id, 'valid_max', rval=a, pack=field%pack )
          end if
      end if  
      if( field%missing.NE.default_field%missing )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'missing_value', rval=field%missing, pack=field%pack )
          else
              a = nint((field%missing-field%add)/field%scale)
              call mpp_write_meta( unit, field%id, 'missing_value', rval=a, pack=field%pack )
          end if
      end if  
      if( field%fill.NE.default_field%fill )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, '_FillValue', rval=field%missing, pack=field%pack )
          else
              a = nint((field%fill-field%add)/field%scale)
              call mpp_write_meta( unit, field%id, '_FillValue', rval=a, pack=field%pack )
          end if
      end if  
      if( field%pack.NE.1 .AND. field%pack.NE.2 )then
          call mpp_write_meta( unit, field%id, 'packing', ival=field%pack )
          if( field%scale.NE.default_field%scale )call mpp_write_meta( unit, field%id, 'scale_factor',  rval=field%scale )
          if( field%add.NE.default_field%add   )call mpp_write_meta( unit, field%id, 'add_offset',    rval=field%add   )
      end if
      if( verbose )print '(a,2i3,x,a,i3)', 'MPP_WRITE_META: Wrote field metadata: pe, unit, field%name, field%id=', &
           pe, unit, trim(field%name), field%id

!      call mpp_clock_end(mpp_write_clock)           
      return
    end subroutine mpp_copy_meta_field

    subroutine mpp_modify_axis_meta( axis, name, units, longname, cartesian, data )
    
      type(axistype), intent(inout) :: axis
      character(len=*), intent(in), optional :: name, units, longname, cartesian
      real, dimension(:), intent(in), optional :: data
      
      if (PRESENT(name)) axis%name = trim(name)
      if (PRESENT(units)) axis%units = trim(units)
      if (PRESENT(longname)) axis%longname = trim(longname)
      if (PRESENT(cartesian)) axis%cartesian = trim(cartesian)
      if (PRESENT(data)) then
         axis%len = size(data(:))
         if (ASSOCIATED(axis%data)) deallocate(axis%data)
         allocate(axis%data(axis%len))
         axis%data = data
      endif
         
      return
    end subroutine mpp_modify_axis_meta
      
    subroutine mpp_modify_field_meta( field, name, units, longname, min, max, missing, axes )
    
      type(fieldtype), intent(inout) :: field
      character(len=*), intent(in), optional :: name, units, longname
      real, intent(in), optional :: min, max, missing
      type(axistype), dimension(:), intent(inout), optional :: axes
      
      if (PRESENT(name)) field%name = trim(name)
      if (PRESENT(units)) field%units = trim(units)
      if (PRESENT(longname)) field%longname = trim(longname)
      if (PRESENT(min)) field%min = min
      if (PRESENT(max)) field%max = max
      if (PRESENT(missing)) field%missing = missing
!      if (PRESENT(axes)) then
!         axis%len = size(data(:))
!         deallocate(axis%data)
!         allocate(axis%data(axis%len))
!         axis%data = data
!      endif
        
      return
    end subroutine mpp_modify_field_meta

