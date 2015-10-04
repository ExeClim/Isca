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

! -*-f90-*-
! $Id: mpp_io_misc.inc,v 19.0 2012/01/06 22:03:22 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!               mpp_io_init: initialize parallel I/O                   !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <SUBROUTINE NAME="mpp_io_init">
!   <OVERVIEW>
!    Initialize <TT>mpp_io_mod</TT>.
!   </OVERVIEW>
!   <DESCRIPTION>
!    Called to initialize the <TT>mpp_io_mod</TT> package. Sets the range
!    of valid fortran units and initializes the <TT>mpp_file</TT> array of
!    <TT>type(filetype)</TT>.  <TT>mpp_io_init</TT> will call <TT>mpp_init</TT> and
!    <TT>mpp_domains_init</TT>, to make sure its parent modules have been
!    initialized. (Repeated calls to the <TT>init</TT> routines do no harm,
!    so don't worry if you already called it).
!   </DESCRIPTION>
!   <TEMPLATE>
!    call mpp_io_init( flags, maxunit )
!   </TEMPLATE>
!   <IN NAME="flags" TYPE="integer"></IN>
!   <IN NAME="maxunit" TYPE="integer"></IN>
! </SUBROUTINE>

    subroutine mpp_io_init( flags, maxunit )
      integer, intent(in), optional :: flags, maxunit
      integer                       :: unit_nml, io_status, iunit
      integer                       :: logunit, outunit, inunit, errunit
      logical                       :: opened
      real(DOUBLE_KIND)             :: doubledata = 0
      real                          :: realarray(4)

      if( module_is_initialized )return

!initialize IO package: initialize mpp_file array, set valid range of units for fortran IO

      call mpp_init(flags)           !if mpp_init has been called, this call will merely return
      pe = mpp_pe()
      npes = mpp_npes()
      call mpp_domains_init(flags)

      maxunits = 1024
      if( PRESENT(maxunit) )maxunits = maxunit
      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
      end if

!set range of allowed fortran unit numbers: could be compiler-dependent (should not overlap stdin/out/err)
      call mpp_set_unit_range( 103, maxunits )

      !--- namelist 
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, mpp_io_nml, iostat=io_status)
#else
      do unit_nml = unit_begin, unit_end
         inquire( unit_nml,OPENED=opened )
         if( .NOT.opened )exit
      end do
      open(unit_nml,file='input.nml')
      read(unit_nml,mpp_io_nml,iostat=io_status)
      close(unit_nml)
#endif
      outunit = stdout(); logunit=stdlog()
      write(outunit, mpp_io_nml)
      write(logunit, mpp_io_nml)

!--- check the deflate level, set deflate = 1 if deflate_level is greater than equal to 0
      if(deflate_level .GE. 0) deflate = 1
      if(deflate .NE. 0) then
         if(deflate_level <0 .OR. deflate > 9) then
            call mpp_error(FATAL, "mpp_io_mod(mpp_io_init): mpp_io_nml variable must be between 0 and 9 when set")
         endif
      endif

! determine the pack_size
      pack_size = size(transfer(doubledata, realarray))
      if( pack_size .NE. 1 .AND. pack_size .NE. 2) call mpp_error(FATAL,'mpp_io_mod(mpp_io_init): pack_size should be 1 or 2')  

!initialize default_field
      default_field%name = 'noname'
      default_field%units = 'nounits'
      default_field%longname = 'noname'
      default_field%id = -1
      default_field%type = -1
      default_field%natt = -1
      default_field%ndim = -1
!largest possible 4-byte reals
      default_field%min = -huge(1._4)
      default_field%max =  huge(1._4)
      default_field%missing = -1e36
      default_field%fill = -1e36
      default_field%scale = 1.0
      default_field%add = 0.0
      default_field%pack = 1
      default_field%time_axis_index = -1 !this value will never match any index
! Initialize default axis
      default_axis%name = 'noname'
      default_axis%units = 'nounits'
      default_axis%longname = 'noname'
      default_axis%cartesian = 'none'
      default_axis%calendar = 'unspecified'
      default_axis%sense = 0
      default_axis%len = -1
      default_axis%id = -1
      default_axis%did = -1
      default_axis%type = -1
      default_axis%natt = -1
! Initialize default attribute
      default_att%name = 'noname'
      default_att%type = -1
      default_att%len = -1
      default_att%catt = 'none'
      
!up to MAXUNITS fortran units and MAXUNITS netCDF units are supported
!file attributes (opened, format, access, threading, fileset) are saved against the unit number
!external handles to netCDF units are saved from maxunits+1:2*maxunits
      allocate( mpp_file(NULLUNIT:2*maxunits) ) !starts at NULLUNIT=-1, used by non-participant PEs in single-threaded I/O
      mpp_file(:)%name   = ' '
      mpp_file(:)%action    = -1
      mpp_file(:)%format    = -1
      mpp_file(:)%threading = -1
      mpp_file(:)%fileset   = -1
      mpp_file(:)%record    = -1
      mpp_file(:)%ncid      = -1
      mpp_file(:)%opened = .FALSE.
      mpp_file(:)%initialized = .FALSE.
      mpp_file(:)%write_on_this_pe = .FALSE.
      mpp_file(:)%io_domain_exist = .FALSE.
      mpp_file(:)%time_level = 0
      mpp_file(:)%time = NULLTIME
      mpp_file(:)%id = -1
      mpp_file(:)%valid = .FALSE.     
      mpp_file(:)%ndim = -1
      mpp_file(:)%nvar = -1
!NULLUNIT "file" is always single-threaded, open and initialized (to pass checks in mpp_write)
      mpp_file(NULLUNIT)%threading = MPP_SINGLE
      mpp_file(NULLUNIT)%opened = .TRUE.
      mpp_file(NULLUNIT)%valid  = .TRUE.
      mpp_file(NULLUNIT)%initialized = .TRUE.
!declare the stdunits to be open
      mpp_file(outunit)%opened = .TRUE.
      mpp_file(logunit)%opened = .TRUE.
      inunit  = stdin()  ; mpp_file(inunit)%opened  = .TRUE.
      errunit = stderr() ; mpp_file(errunit)%opened = .TRUE.
      
      if( pe.EQ.mpp_root_pe() )then
          iunit = stdlog()  ! PGI compiler does not like stdlog() doing I/O within write call
          write( iunit,'(/a)' )'MPP_IO module '//trim(version)
          write( iunit,'( a)' )'MPP_IO module '//trim(tagname)
#ifdef use_netCDF
          text = NF_INQ_LIBVERS()
          write( iunit,'(/a)' )'Using netCDF library version '//trim(text)
#endif    
      endif
      
#ifdef CRAYPVP
!we require every file to be assigned threadwise: PVPs default to global, and are reset here
      call ASSIGN( 'assign -P thread p:%', error )
#endif
      
      call mpp_io_set_stack_size(131072) ! default initial value
      call mpp_sync()
      if( io_clocks_on )then
          mpp_read_clock  = mpp_clock_id( 'mpp_read')
          mpp_write_clock  = mpp_clock_id( 'mpp_write')
          mpp_open_clock  = mpp_clock_id( 'mpp_open')
          mpp_close_clock  = mpp_clock_id( 'mpp_close')
      endif
      module_is_initialized = .TRUE.
      return
    end subroutine mpp_io_init


! <SUBROUTINE NAME="mpp_io_exit">
!   <OVERVIEW>
!    Exit <TT>mpp_io_mod</TT>.
!   </OVERVIEW>
!   <DESCRIPTION>
!    It is recommended, though not at present required, that you call this
!    near the end of a run. This will close all open files that were opened
!    with <LINK SRC="#mpp_open"><TT>mpp_open</TT></LINK>. Files opened otherwise
!    are not affected.
!   </DESCRIPTION>
!   <TEMPLATE>
!    call mpp_io_exit()
!   </TEMPLATE>
! </SUBROUTINE>

    subroutine mpp_io_exit(string)
      character(len=*), optional :: string
      integer :: unit,istat
      logical :: dosync
      
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_IO_EXIT: must first call mpp_io_init.' )
      dosync = .TRUE.
      if( PRESENT(string) )then
          dosync = .NOT.( trim(string).EQ.'NOSYNC' )
      end if
!close all open fortran units
      do unit = unit_begin,unit_end
         if( mpp_file(unit)%opened )call FLUSH(unit)
      end do
      if( dosync )call mpp_sync()
      do unit = unit_begin,unit_end
         if( mpp_file(unit)%opened )close(unit)
      end do
#ifdef use_netCDF
!close all open netCDF units
      do unit = maxunits+1,2*maxunits
         if( mpp_file(unit)%opened )error = NF_CLOSE(mpp_file(unit)%ncid)
      end do
#endif   
      
!      call mpp_max(mpp_io_stack_hwm)
      
      if( pe.EQ.mpp_root_pe() )then
!          write( stdout,'(/a)' )'Exiting MPP_IO module...'
!          write( stdout,* )'MPP_IO_STACK high water mark=', mpp_io_stack_hwm
      end if
      deallocate(mpp_file)
      module_is_initialized = .FALSE.
      return
    end subroutine mpp_io_exit


    subroutine netcdf_err( err, file, axis, field, attr, string )
      integer, intent(in) :: err
      type(filetype), optional :: file
      type(axistype), optional :: axis
      type(fieldtype), optional :: field
      type(atttype), optional :: attr
      character(len=*), optional :: string
      character(len=256) :: errmsg

#ifdef use_netCDF
      if( err.EQ.NF_NOERR )return
      errmsg = NF_STRERROR(err)
      if( PRESENT(file) )errmsg = trim(errmsg)//' File='//file%name
      if( PRESENT(axis) )errmsg = trim(errmsg)//' Axis='//axis%name
      if( PRESENT(field) )errmsg = trim(errmsg)//' Field='//field%name
      if( PRESENT(attr) )errmsg = trim(errmsg)//' Attribute='//attr%name
      if( PRESENT(string) )errmsg = trim(errmsg)//string
      call mpp_io_exit('NOSYNC')        !make sure you close all open files
      call mpp_error( FATAL, 'NETCDF ERROR: '//trim(errmsg) )
#endif
      return
    end subroutine netcdf_err


    subroutine mpp_flush(unit)
!flush the output on a unit, syncing with disk
      integer, intent(in) :: unit

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_FLUSH: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%write_on_this_pe) return
      if( .NOT.mpp_file(unit)%opened ) call mpp_error( FATAL, 'MPP_FLUSH: invalid unit number.' )
      if( .NOT.mpp_file(unit)%initialized )call mpp_error( FATAL, 'MPP_FLUSH: cannot flush a file during writing of metadata.' )

      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          error = NF_SYNC(mpp_file(unit)%ncid); call netcdf_err( error, mpp_file(unit) )
#endif
      else
          call FLUSH(unit)
      end if
      return
    end subroutine mpp_flush


