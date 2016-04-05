
! -*-f90-*-
! $Id: mpp_io_connect.inc,v 18.0 2010/03/02 23:57:36 fms Exp $

! <SUBROUTINE NAME="mpp_open">

!   <OVERVIEW>
!     Open a file for parallel I/O.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Open a file for parallel I/O.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_open( unit, file, action, form, access, threading, fileset,
!             iospec, nohdrs, recl, pelist )
!   </TEMPLATE>

!   <OUT NAME="unit" TYPE="integer">
!     unit is intent(OUT): always _returned_by_ mpp_open().
!   </OUT>
!   <IN NAME="file" TYPE="character(len=*)">
!     file is the filename: REQUIRED
!    we append .nc to filename if it is a netCDF file
!    we append .<pppp> to filename if fileset is private (pppp is PE number)
!   </IN>
!   <IN NAME="action" TYPE="integer">
!     action is one of MPP_RDONLY, MPP_APPEND, MPP_WRONLY or MPP_OVERWR.
!   </IN>
!   <IN NAME="form" TYPE="integer">
!     form is one of MPP_ASCII:  formatted read/write
!                   MPP_NATIVE: unformatted read/write with no conversion
!                   MPP_IEEE32: unformatted read/write with conversion to IEEE32
!                   MPP_NETCDF: unformatted read/write with conversion to netCDF
!   </IN>
!   <IN NAME="access" TYPE="integer">
!     access is one of MPP_SEQUENTIAL or MPP_DIRECT (ignored for netCDF).
!     RECL argument is REQUIRED for direct access IO.
!   </IN>
!   <IN NAME="threading" TYPE="integer">
!     threading is one of MPP_SINGLE or MPP_MULTI
!      single-threaded IO in a multi-PE run is done by PE0.
!   </IN>
!   <IN NAME="fileset" TYPE="integer">
!     fileset is one of MPP_MULTI and MPP_SINGLE
!     fileset is only used for multi-threaded I/O
!     if all I/O PEs in <pelist> use a single fileset, they write to the same file
!     if all I/O PEs in <pelist> use a multi  fileset, they each write an independent file
!   </IN>
!   <IN NAME="pelist" TYPE="integer">
!     pelist is the list of I/O PEs (currently ALL).
!   </IN>
!   <IN NAME="recl" TYPE="integer">
!     recl is the record length in bytes.
!   </IN>
!   <IN NAME="iospec" TYPE="character(len=*)">
!     iospec is a system hint for I/O organization, e.g assign(1) on SGI/Cray systems.
!   </IN>
!   <IN NAME="nohdrs" TYPE="logical">
!     nohdrs has no effect when action=MPP_RDONLY|MPP_APPEND or when form=MPP_NETCDF
!   </IN>
!   <NOTE>
!     The integer parameters to be passed as flags (<TT>MPP_RDONLY</TT>,
!   etc) are all made available by use association. The <TT>unit</TT>
!   returned by <TT>mpp_open</TT> is guaranteed unique. For non-netCDF I/O
!   it is a valid fortran unit number and fortran I/O can be directly called
!   on the file.
!
!   <TT>MPP_WRONLY</TT> will guarantee that existing files named
!   <TT>file</TT> will not be clobbered. <TT>MPP_OVERWR</TT>
!   allows overwriting of files.
!
!   Files opened read-only by many processors will give each processor
!   an independent pointer into the file, i.e:
!
!   <PRE>
!      namelist / nml / ...
!   ...
!      call mpp_open( unit, 'input.nml', action=MPP_RDONLY )
!      read(unit,nml)
!   </PRE>
!
!   will result in each PE independently reading the same namelist.
!
!   Metadata identifying the file and the version of
!   <TT>mpp_io_mod</TT> are written to a file that is opened
!   <TT>MPP_WRONLY</TT> or <TT>MPP_OVERWR</TT>. If this is a
!   multi-file set, and an additional global attribute
!   <TT>NumFilesInSet</TT> is written to be used by post-processing
!   software.
!
!   If <TT>nohdrs=.TRUE.</TT> all calls to write attributes will
!   return successfully <I>without</I> performing any writes to the
!   file. The default is <TT>.FALSE.</TT>.
!
!   For netCDF files, headers are always written even if
!   <TT>nohdrs=.TRUE.</TT>
!
!   The string <TT>iospec</TT> is passed to the OS to
!   characterize the I/O to be performed on the file opened on
!   <TT>unit</TT>. This is typically used for I/O optimization. For
!   example, the FFIO layer on SGI/Cray systems can be used for
!   controlling synchronicity of reads and writes, buffering of data
!   between user space and disk for I/O optimization, striping across
!   multiple disk partitions, automatic data conversion and the like
!   (<TT>man intro_ffio</TT>). All these actions are controlled through
!   the <TT>assign</TT> command. For example, to specify asynchronous
!   caching of data going to a file open on <TT>unit</TT>, one would do:
!
!   <PRE>
!   call mpp_open( unit, ... iospec='-F cachea' )
!   </PRE>
!
!   on an SGI/Cray system, which would pass the supplied
!   <TT>iospec</TT> to the <TT>assign(3F)</TT> system call.
!
!   Currently <TT>iospec </TT>performs no action on non-SGI/Cray
!   systems. The interface is still provided, however: users are cordially
!   invited to add the requisite system calls for other systems.
!   </NOTE>
! </SUBROUTINE>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!           OPENING AND CLOSING FILES: mpp_open() and mpp_close()            !
!                                                                            !
! mpp_open( unit, file, action, form, access, threading, &                   !
!           fileset, iospec, nohdrs, recl, pelist )                          !
!      integer, intent(out) :: unit                                          !
!      character(len=*), intent(in) :: file                                  !
!      integer, intent(in), optional :: action, form, access, threading,     !
!                                       fileset, recl                        !
!      character(len=*), intent(in), optional :: iospec                      !
!      logical, intent(in), optional :: nohdrs                               !
!      integer, optional, intent(in) :: pelist(:) !default ALL               !
!                                                                            !
!  unit is intent(OUT): always _returned_by_ mpp_open()                      !
!  file is the filename: REQUIRED                                            !
!    we append .nc to filename if it is a netCDF file                        !
!    we append .<pppp> to filename if fileset is private (pppp is PE number) !
!  iospec is a system hint for I/O organization                              !
!         e.g assign(1) on SGI/Cray systems.                                 !
!  if nohdrs is .TRUE. headers are not written on non-netCDF writes.         !
!  nohdrs has no effect when action=MPP_RDONLY|MPP_APPEND                    !
!                    or when form=MPP_NETCDF                                 !
! FLAGS:                                                                     !
!    action is one of MPP_RDONLY, MPP_APPEND or MPP_WRONLY                   !
!    form is one of MPP_ASCII:  formatted read/write                         !
!                   MPP_NATIVE: unformatted read/write, no conversion        !
!                   MPP_IEEE32: unformatted read/write, conversion to IEEE32 !
!                   MPP_NETCDF: unformatted read/write, conversion to netCDF !
!    access is one of MPP_SEQUENTIAL or MPP_DIRECT (ignored for netCDF)      !
!      RECL argument is REQUIRED for direct access IO                        !
!    threading is one of MPP_SINGLE or MPP_MULTI                             !
!      single-threaded IO in a multi-PE run is done by PE0                   !
!    fileset is one of MPP_MULTI and MPP_SINGLE                              !
!      fileset is only used for multi-threaded I/O                           !
!      if all I/O PEs in <pelist> use a single fileset,                      !
!              they write to the same file                                   !
!      if all I/O PEs in <pelist> use a multi  fileset,                      !
!              they each write an independent file                           !
!  recl is the record length in bytes                                        !
!  pelist is the list of I/O PEs (currently ALL)                             !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_open( unit, file, action, form, access, threading, &
                         fileset, iospec, nohdrs, recl,               &
                         iostat, is_root_pe, domain )
      integer,         intent(out)           :: unit   
      character(len=*), intent(in)           :: file
      integer,          intent(in), optional :: action, form, access 
      integer,          intent(in), optional :: threading, fileset, recl
      character(len=*), intent(in), optional :: iospec
      logical,          intent(in), optional :: nohdrs
      integer,         intent(out), optional :: iostat
      logical,          intent(in), optional :: is_root_pe      
      type(domain2d),   intent(in), optional :: domain

      character(len=16)  :: act, acc, for, pos
      character(len=128) :: mesg
      integer :: action_flag, form_flag, access_flag, threading_flag, fileset_flag, length
      integer :: nfiles, tile_id(1), io_layout(2)
      logical :: exists, on_root_pe, dist_file
      logical :: write_on_this_pe, io_domain_exist
      integer :: ios, nc_pos     !position of .nc in file name
      type(axistype) :: unlim    !used by netCDF with mpp_append
      type(domain2d), pointer :: io_domain=>NULL()
      integer*8				     :: lenp
      integer			             :: comm
      integer				     :: info, ierror
      integer,dimension(:), allocatable     :: glist(:)  
      integer				    ::lena, lenb
      character(len=12)		            ::ncblk
      character(len=128)                    ::nc_name 
      integer				    ::f_size, f_stat
      integer				    ::fsize, inital = 0
      character(len=128)		    :: f_test
      
      
      
      call mpp_clock_begin(mpp_open_clock)
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_OPEN: must first call mpp_io_init.' )
      on_root_pe = mpp_pe() == mpp_root_pe()
      if(present(is_root_pe)) on_root_pe = is_root_pe

      dist_file = .false.
!set flags
      action_flag = MPP_WRONLY        !default
      if( PRESENT(action) )action_flag = action
      form_flag = MPP_ASCII
      if( PRESENT(form) )form_flag = form
#ifndef use_netCDF
      if( form_flag.EQ.MPP_NETCDF ) &
           call mpp_error( FATAL, 'MPP_OPEN: To open a file with form=MPP_NETCDF, you must compile mpp_io with -Duse_netCDF.' )
#endif     
      access_flag = MPP_SEQUENTIAL
      if( PRESENT(access) )access_flag = access
      threading_flag = MPP_SINGLE
      if( npes.GT.1 .AND. PRESENT(threading) )threading_flag = threading
      fileset_flag = MPP_MULTI
      if( PRESENT(fileset) )fileset_flag = fileset
      if( threading_flag.EQ.MPP_SINGLE )fileset_flag = MPP_SINGLE

      io_domain_exist = .false.
      if( PRESENT(domain) ) then
         io_domain => mpp_get_io_domain(domain)
         if(associated(io_domain)) io_domain_exist = .true.
      endif

      !--- when io_domain_exist, will ignore the threading_flag and fileset_flag.
      if(io_domain_exist) then
         call mpp_error(NOTE, "mpp_io_connect.inc(mpp_open): io_domain exists for domain "// &
            trim(mpp_get_domain_name(domain))//", optional argument fileset will be ignored")
      endif

      write_on_this_pe = .true.
      if( threading_flag.EQ.MPP_SINGLE .AND. .NOT. on_root_pe ) write_on_this_pe = .false.
      if(form_flag == MPP_NETCDF .AND. action_flag .NE. MPP_RDONLY) then
         if(fileset_flag .EQ.MPP_SINGLE .AND.  threading_flag.EQ.MPP_MULTI) then
            call mpp_error(FATAL, "mpp_io_connect.inc(mpp_open): multiple thread and single "// &
                "file writing/appending is not supported for netCDF file")
         endif
         if( fileset_flag  .EQ.MPP_SINGLE .AND. .NOT. on_root_pe ) write_on_this_pe = .false.
      endif   
   
      if( io_domain_exist) then
          ! in this case, only write out from the root_pe of io_domain.
          write_on_this_pe = mpp_domain_is_tile_root_pe(io_domain)
      endif 
      if( action_flag == MPP_RDONLY) write_on_this_pe = .false.
      !get a unit number
      if( .NOT. write_on_this_pe .AND. action_flag.NE.MPP_RDONLY .AND. .NOT. io_domain_exist)then
            unit = NULLUNIT           !PEs not participating in IO from this mpp_open() will return this value for unit
         call mpp_clock_end(mpp_open_clock)
         return
      end if
      if( form_flag.EQ.MPP_NETCDF )then
          do unit = maxunits+1,2*maxunits
             if( .NOT.mpp_file(unit)%valid )exit
          end do
          if( unit.GT.2*maxunits ) then
            write(mesg,*) 'all the units between ',maxunits+1,' and ',2*maxunits,' are used'
            call mpp_error( FATAL, 'MPP_OPEN: too many open netCDF files.'//trim(mesg) )
          endif
      else
          do unit = unit_begin, unit_end
             inquire( unit,OPENED=mpp_file(unit)%opened )
             if( .NOT.mpp_file(unit)%opened )exit
          end do
          if( unit.GT.unit_end ) then
            write(mesg,*) 'all the units between ',unit_begin,' and ',unit_end,' are used'
            call mpp_error( FATAL, 'MPP_OPEN: no available units.'//trim(mesg) )
          endif
      end if
      mpp_file(unit)%valid            = .true.
      mpp_file(unit)%write_on_this_pe = write_on_this_pe
      mpp_file(unit)%io_domain_exist  = io_domain_exist
      if( PRESENT(domain) ) then
         allocate(mpp_file(unit)%domain)
         mpp_file(unit)%domain = domain
      endif

!get a filename
      nc_pos = index(file,'.nc.')
      dist_file = nc_pos>0  ! this is a distributed file ending with filename.nc.0???
      text = file
      length = len_trim(file)      
      if(form_flag.EQ.MPP_NETCDF.AND. file(length-2:length) /= '.nc' .AND. .NOT.dist_file) &
         text = trim(file)//'.nc'
      if(present(domain)) io_layout = mpp_get_io_domain_layout(domain) 
      if( io_domain_exist) then
         if( io_layout(1)*io_layout(2) > 1) then
            fileset_flag = MPP_MULTI
            threading_flag = MPP_MULTI
            tile_id = mpp_get_tile_id(io_domain)
            if(mpp_npes() > 10000) then 
               write( text,'(a,i6.6)' )trim(text)//'.', tile_id(1)
            else
               write( text,'(a,i4.4)' )trim(text)//'.', tile_id(1)
            endif
         else
            fileset_flag   = MPP_SINGLE
            threading_flag = MPP_SINGLE
         endif
      else if( fileset_flag.EQ.MPP_MULTI ) then
         if(mpp_npes() > 10000) then
            write( text,'(a,i6.6)' )trim(text)//'.', pe
         else
            write( text,'(a,i4.4)' )trim(text)//'.', pe
         endif
      endif
      mpp_file(unit)%name = text
      if( verbose )print '(a,2i6,x,a,5i5)', 'MPP_OPEN: PE, unit, filename, action, format, access, threading, fileset=', &
           pe, unit, trim(mpp_file(unit)%name), action_flag, form_flag, access_flag, threading_flag, fileset_flag
           
!action: read, write, overwrite, append: act and pos are ignored by netCDF
      if( action_flag.EQ.MPP_RDONLY )then
          act = 'READ'
          pos = 'REWIND'
      else if( action_flag.EQ.MPP_WRONLY .OR. action_flag.EQ.MPP_OVERWR )then
          act = 'WRITE'
          pos = 'REWIND'
      else if( action_flag.EQ.MPP_APPEND )then
          act = 'WRITE'
          pos = 'APPEND'
      else
          call mpp_error( FATAL, 'MPP_OPEN: action must be one of MPP_WRONLY, MPP_APPEND or MPP_RDONLY.' )
      end if

      mpp_file(unit)%threading = threading_flag
      mpp_file(unit)%fileset = fileset_flag

      if( .NOT. write_on_this_pe .AND. action_flag.NE.MPP_RDONLY ) then
         call mpp_clock_end(mpp_open_clock)
         return
      endif

!access: sequential or direct: ignored by netCDF
      if( form_flag.NE.MPP_NETCDF )then
          if( access_flag.EQ.MPP_SEQUENTIAL )then
              acc = 'SEQUENTIAL'
          else if( access_flag.EQ.MPP_DIRECT )then
              acc = 'DIRECT'
              if( form_flag.EQ.MPP_ASCII )call mpp_error( FATAL, 'MPP_OPEN: formatted direct access I/O is prohibited.' )
              if( .NOT.PRESENT(recl) ) &
                   call mpp_error( FATAL, 'MPP_OPEN: recl (record length in bytes) must be specified with access=MPP_DIRECT.' )
              mpp_file(unit)%record = 1
              records_per_pe = 1 !each PE writes 1 record per mpp_write
          else
              call mpp_error( FATAL, 'MPP_OPEN: access must be one of MPP_SEQUENTIAL or MPP_DIRECT.' )
          end if
      end if  
          
!threading: SINGLE or MULTI
      if( threading_flag.EQ.MPP_MULTI )then
!fileset: MULTI or SINGLE (only for multi-threaded I/O
          if( fileset_flag.EQ.MPP_SINGLE )then
              if( form_flag.EQ.MPP_NETCDF .AND. act.EQ.'WRITE' ) &
                   call mpp_error( FATAL, 'MPP_OPEN: netCDF currently does not support single-file multi-threaded output.' )
                   
#ifdef _CRAYT3E    
              call ASSIGN( 'assign -I -F global.privpos f:'//trim(mpp_file(unit)%name), error )
#endif        
          else if( fileset_flag.NE.MPP_MULTI )then
              call mpp_error( FATAL, 'MPP_OPEN: fileset must be one of MPP_MULTI or MPP_SINGLE.' )
          end if
      else if( threading_flag.NE.MPP_SINGLE )then
          call mpp_error( FATAL, 'MPP_OPEN: threading must be one of MPP_SINGLE or MPP_MULTI.' )
      end if
          
!apply I/O specs before opening the file
!note that -P refers to the scope of a fortran unit, which is always thread-private even if file is shared
#ifdef CRAYPVP
#ifndef _CRAYX1
      call ASSIGN( 'assign -I -P thread  f:'//trim(mpp_file(unit)%name), error )
#endif
#endif
#ifdef _CRAYT3E
      call ASSIGN( 'assign -I -P private f:'//trim(mpp_file(unit)%name), error )
#endif
#ifdef _CRAYX1
      if (file(length-3:length) == '.nml') then
        call ASSIGN( 'assign -I -f77 f:'//trim(mpp_file(unit)%name), error )
!       call ASSIGN( 'assign -I -F global f:'//trim(mpp_file(unit)%name), error )
      endif
#endif
      if( PRESENT(iospec) )then
!iospec provides hints to the system on how to organize I/O
!on Cray systems this is done through 'assign', see assign(1) and assign(3F)
!on other systems this will be expanded as needed
!no error checks here on whether the supplied iospec is valid
#if defined(SGICRAY) || defined(_CRAYX1)
          call ASSIGN( 'assign -I '//trim(iospec)//' f:'//trim(mpp_file(unit)%name), error )
          if( form_flag.EQ.MPP_NETCDF )then
!for netCDF on SGI/Cray systems we pass it to the environment variable NETCDF_XFFIOSPEC
!ideally we should parse iospec, pass the argument of -F to NETCDF_FFIOSPEC, and the rest to NETCDF_XFFIOSPEC
!maybe I'll get around to it someday
!PXFSETENV is a POSIX-standard routine for setting environment variables from fortran
              !if we ever use this again....F2003 non-intel flavor of getenv
	      !call get_enviornment_variable( 'NETCDF_XFFIOSPEC', trim(iospec))
	      call PXFSETENV( 'NETCDF_XFFIOSPEC', 0, trim(iospec), 0, 1, error )
          end if
#endif        
      end if
      
!open the file as specified above for various formats
      if( form_flag.EQ.MPP_NETCDF )then
#ifdef PAR_ACCESS  
      !adding some items for netcdf-4.... --fmi
      if( .NOT.PRESENT(pelist)) then
	allocate (glist(0:npes-1))
	call mpp_get_current_pelist(glist, name, comm) 
      end if
      
      !create info parts...
      call MPI_INFO_CREATE(info, ierror)
      
     !F2003 convention
     !call get_enviornment_variable('NC_BLKSZ', ncblk)
     call GETENV( 'NC_BLKSZ', ncblk)
     ncblk = trim(ncblk)
     
  
     if (ncblk /= "") then
      


	call MPI_INFO_SET(info, "cb_buffer_size", ncblk, ierror)
     	call MPI_INFO_SET(info, "ind_rd_buffer_size", ncblk, ierror)
     	call MPI_INFO_SET(info, "ind_wr_buffer_size", ncblk, ierror)

     else
     	call MPI_INFO_SET(info, "ind_rd_buffer_size", "16777216", ierror)
	call MPI_INFO_SET(info, "ind_wr_buffer_size", "16777216", ierror)
     endif

#else
!added by fmi to read NC_BLKSZ and NC_BLKSZ_filename...

		
		!get regular nc_blksz...
     		 !build env var for check
		 !write (*,*) 'hello', trim(mpp_file(unit)%name)
     		 nc_name = 'NC_BLKSZ_'//trim(mpp_file(unit)%name)
		 !write (*,*) 'nc_name: ', nc_name,  ' bcblk: ', ncblk
		 
      
		  
    		  !make the call.....
   		  !f2003 replaces GETENV with get_enviornment_variable so the guts are here if we need to switch
		  !call get_enviornment_variable(trim(nc_name),ncblk )
		  call GETENV( trim(nc_name),ncblk )
		   
   		   !might not be there...use the general setting
    	  
   		   if (ncblk .EQ. '') then
		   	!call get_enviornment_variable( 'NC_BLKSZ', ncblk)
    		  	call GETENV( 'NC_BLKSZ', ncblk)
			
			
  		   endif
      
   		   !if no general setting then use default
  		    if (ncblk .EQ. '') then
  		   	ncblk = '64k'  !change for platform...perhaps we should set an ifdef for this....
 		    endif
      
  		   !set or convert the chunksize
		   
		    call file_size(ncblk, mpp_file(unit)%name, fsize)
		    !write (*,*) 'this is fsize after: ', fsize
		 
     	
    		  if(debug) write(*,*) 'Blocksize for ', trim(mpp_file(unit)%name),' is ', fsize
		  !ends addition from fmi - oct.22.2008
#endif

#ifdef use_netCDF
#ifdef use_netCDF3
          if( action_flag.EQ.MPP_WRONLY )then
	      if(debug) write(*,*) 'Blocksize for create of ', trim(mpp_file(unit)%name),' is ', fsize
	      error = NF__CREATE( trim(mpp_file(unit)%name), NF_NOCLOBBER, inital, fsize, mpp_file(unit)%ncid )
              call netcdf_err( error, mpp_file(unit) )
              if( verbose )print '(a,i6,i16)', 'MPP_OPEN: new netCDF file: pe, ncid=', pe, mpp_file(unit)%ncid
          else if( action_flag.EQ.MPP_OVERWR )then
	      if(debug) write(*,*) 'Blocksize for create of ', trim(mpp_file(unit)%name),' is ', fsize
	      error = NF__CREATE( trim(mpp_file(unit)%name),NF_CLOBBER, inital, fsize, mpp_file(unit)%ncid )
              call netcdf_err( error, mpp_file(unit) )
              action_flag = MPP_WRONLY !after setting clobber, there is no further distinction btwn MPP_WRONLY and MPP_OVERWR
              if( verbose )print '(a,i6,i16)', 'MPP_OPEN: overwrite netCDF file: pe, ncid=', pe, mpp_file(unit)%ncid
          else if( action_flag.EQ.MPP_APPEND )then
              inquire(file=trim(mpp_file(unit)%name),EXIST=exists)
              if (.NOT.exists) call mpp_error(FATAL,'MPP_OPEN:'&
                   &//trim(mpp_file(unit)%name)//' does not exist.')
              error = NF__OPEN( trim(mpp_file(unit)%name), NF_WRITE, fsize, mpp_file(unit)%ncid ); call netcdf_err( error, mpp_file(unit))
!get the current time level of the file: writes to this file will be at next time level
              error = NF_INQ_UNLIMDIM( mpp_file(unit)%ncid, unlim%did )
              if( error.EQ.NF_NOERR )then
                  error = NF_INQ_DIM( mpp_file(unit)%ncid, unlim%did, unlim%name, mpp_file(unit)%time_level )
                  call netcdf_err( error, mpp_file(unit) )
                  error = NF_INQ_VARID( mpp_file(unit)%ncid, unlim%name, mpp_file(unit)%id )
                  call netcdf_err( error, mpp_file(unit), unlim )
              end if
              if( verbose )print '(a,i6,i16,i4)', 'MPP_OPEN: append to existing netCDF file: pe, ncid, time_axis_id=',&
                   pe, mpp_file(unit)%ncid, mpp_file(unit)%id
              mpp_file(unit)%format=form_flag ! need this for mpp_read
              call mpp_read_meta(unit)
          else if( action_flag.EQ.MPP_RDONLY )then
               inquire(file=trim(mpp_file(unit)%name),EXIST=exists)
              if (.NOT.exists) call mpp_error(FATAL,'MPP_OPEN:'&
                   &//trim(mpp_file(unit)%name)//' does not exist.')
              error = NF__OPEN( trim(mpp_file(unit)%name), NF_NOWRITE, fsize, mpp_file(unit)%ncid ); call netcdf_err( error, mpp_file(unit))
              if( verbose )print '(a,i6,i16,i4)', 'MPP_OPEN: opening existing netCDF file: pe, ncid, time_axis_id=',&
                   pe, mpp_file(unit)%ncid, mpp_file(unit)%id
              mpp_file(unit)%format=form_flag ! need this for mpp_read
              call mpp_read_meta(unit)
          end if
          mpp_file(unit)%opened = .TRUE.
#elif use_LARGEFILE
          if( action_flag.EQ.MPP_WRONLY )then
	      if(debug) write(*,*) 'Blocksize for create of ', trim(mpp_file(unit)%name),' is ', fsize
	      error = NF__CREATE( trim(mpp_file(unit)%name), IOR(NF_64BIT_OFFSET, NF_NOCLOBBER), inital, fsize, mpp_file(unit)%ncid )
              call netcdf_err( error, mpp_file(unit) )
              if( verbose )print '(a,i6,i16)', 'MPP_OPEN: new netCDF file: pe, ncid=', pe, mpp_file(unit)%ncid
          else if( action_flag.EQ.MPP_OVERWR )then
	      if(debug) write(*,*) 'Blocksize for create of ', trim(mpp_file(unit)%name),' is ', fsize
	      error = NF__CREATE( trim(mpp_file(unit)%name),IOR(NF_64BIT_OFFSET,NF_CLOBBER), inital, fsize, mpp_file(unit)%ncid )
              call netcdf_err( error, mpp_file(unit) )
              action_flag = MPP_WRONLY !after setting clobber, there is no further distinction btwn MPP_WRONLY and MPP_OVERWR
              if( verbose )print '(a,i6,i16)', 'MPP_OPEN: overwrite netCDF file: pe, ncid=', pe, mpp_file(unit)%ncid
          else if( action_flag.EQ.MPP_APPEND )then
              inquire(file=trim(mpp_file(unit)%name),EXIST=exists)
              if (.NOT.exists) call mpp_error(FATAL,'MPP_OPEN:'&
                   &//trim(mpp_file(unit)%name)//' does not exist.')
              error = NF__OPEN( trim(mpp_file(unit)%name), NF_WRITE, fsize, mpp_file(unit)%ncid ); call netcdf_err( error, mpp_file(unit))
!get the current time level of the file: writes to this file will be at next time level
              error = NF_INQ_UNLIMDIM( mpp_file(unit)%ncid, unlim%did )
              if( error.EQ.NF_NOERR )then
                  error = NF_INQ_DIM( mpp_file(unit)%ncid, unlim%did, unlim%name, mpp_file(unit)%time_level )
                  call netcdf_err( error, mpp_file(unit) )
                  error = NF_INQ_VARID( mpp_file(unit)%ncid, unlim%name, mpp_file(unit)%id )
                  call netcdf_err( error, mpp_file(unit), unlim )
              end if
              if( verbose )print '(a,i6,i16,i4)', 'MPP_OPEN: append to existing netCDF file: pe, ncid, time_axis_id=',&
                   pe, mpp_file(unit)%ncid, mpp_file(unit)%id
              mpp_file(unit)%format=form_flag ! need this for mpp_read
              call mpp_read_meta(unit)
          else if( action_flag.EQ.MPP_RDONLY )then
               inquire(file=trim(mpp_file(unit)%name),EXIST=exists)
              if (.NOT.exists) call mpp_error(FATAL,'MPP_OPEN:'&
                   &//trim(mpp_file(unit)%name)//' does not exist.')
              error = NF__OPEN( trim(mpp_file(unit)%name), NF_NOWRITE, fsize, mpp_file(unit)%ncid ); call netcdf_err( error, mpp_file(unit))
              if( verbose )print '(a,i6,i16,i4)', 'MPP_OPEN: opening existing netCDF file: pe, ncid, time_axis_id=',&
                   pe, mpp_file(unit)%ncid, mpp_file(unit)%id
              mpp_file(unit)%format=form_flag ! need this for mpp_read
              call mpp_read_meta(unit)
          end if
          mpp_file(unit)%opened = .TRUE.
#else
          if( action_flag.EQ.MPP_WRONLY )then
	      if(debug) write(*,*) 'Blocksize for create of ', trim(mpp_file(unit)%name),' is ', fsize
	      error = NF__CREATE( trim(mpp_file(unit)%name), IOR(NF_NETCDF4,NF_CLASSIC_MODEL), inital, fsize, mpp_file(unit)%ncid )
              call netcdf_err( error, mpp_file(unit) )
              if( verbose )print '(a,i6,i16)', 'MPP_OPEN: new netCDF file: pe, ncid=', pe, mpp_file(unit)%ncid
          else if( action_flag.EQ.MPP_OVERWR )then
	      if(debug) write(*,*) 'Blocksize for create of ', trim(mpp_file(unit)%name),' is ', fsize
	      error = NF__CREATE( trim(mpp_file(unit)%name), IOR(NF_NETCDF4,NF_CLASSIC_MODEL), inital, fsize, mpp_file(unit)%ncid )
              call netcdf_err( error, mpp_file(unit) )
              action_flag = MPP_WRONLY !after setting clobber, there is no further distinction btwn MPP_WRONLY and MPP_OVERWR
              if( verbose )print '(a,i6,i16)', 'MPP_OPEN: overwrite netCDF file: pe, ncid=', pe, mpp_file(unit)%ncid
          else if( action_flag.EQ.MPP_APPEND )then
              inquire(file=trim(mpp_file(unit)%name),EXIST=exists)
              if (.NOT.exists) call mpp_error(FATAL,'MPP_OPEN:'&
                   &//trim(mpp_file(unit)%name)//' does not exist.')
              error = NF__OPEN( trim(mpp_file(unit)%name), NF_WRITE, fsize, mpp_file(unit)%ncid ); call netcdf_err( error, mpp_file(unit))
!get the current time level of the file: writes to this file will be at next time level
              error = NF_INQ_UNLIMDIM( mpp_file(unit)%ncid, unlim%did )
              if( error.EQ.NF_NOERR )then
                  error = NF_INQ_DIM( mpp_file(unit)%ncid, unlim%did, unlim%name, mpp_file(unit)%time_level )
                  call netcdf_err( error, mpp_file(unit) )
                  error = NF_INQ_VARID( mpp_file(unit)%ncid, unlim%name, mpp_file(unit)%id )
                  call netcdf_err( error, mpp_file(unit), unlim )
              end if
              if( verbose )print '(a,i6,i16,i4)', 'MPP_OPEN: append to existing netCDF file: pe, ncid, time_axis_id=',&
                   pe, mpp_file(unit)%ncid, mpp_file(unit)%id
              mpp_file(unit)%format=form_flag ! need this for mpp_read
              call mpp_read_meta(unit)
          else if( action_flag.EQ.MPP_RDONLY )then
               inquire(file=trim(mpp_file(unit)%name),EXIST=exists)
              if (.NOT.exists) call mpp_error(FATAL,'MPP_OPEN:'&
                   &//trim(mpp_file(unit)%name)//' does not exist.')
              error = NF__OPEN( trim(mpp_file(unit)%name), NF_NOWRITE, fsize, mpp_file(unit)%ncid ); call netcdf_err( error, mpp_file(unit))
              if( verbose )print '(a,i6,i16,i4)', 'MPP_OPEN: opening existing netCDF file: pe, ncid, time_axis_id=',&
                   pe, mpp_file(unit)%ncid, mpp_file(unit)%id
              mpp_file(unit)%format=form_flag ! need this for mpp_read
              call mpp_read_meta(unit)
          end if
          mpp_file(unit)%opened = .TRUE.

#endif    
#endif
      else
!format: ascii, native, or IEEE 32 bit
          if( form_flag.EQ.MPP_ASCII )then
              for = 'FORMATTED'
          else if( form_flag.EQ.MPP_IEEE32 )then
              for = 'UNFORMATTED'
!assign -N is currently unsupported on SGI
#ifdef _CRAY  
#ifndef _CRAYX1  
              call ASSIGN( 'assign -I -N ieee_32 f:'//trim(mpp_file(unit)%name), error )
#endif        
#endif        
          else if( form_flag.EQ.MPP_NATIVE )then
              for = 'UNFORMATTED'
          else
              call mpp_error( FATAL, 'MPP_OPEN: form must be one of MPP_ASCII, MPP_NATIVE, MPP_IEEE32 or MPP_NETCDF.' )
          end if
          inquire( file=trim(mpp_file(unit)%name), EXIST=exists )
          if( exists .AND. action_flag.EQ.MPP_WRONLY ) &
               call mpp_error( WARNING, 'MPP_OPEN: File '//trim(mpp_file(unit)%name)//' opened WRONLY already exists!' )
          if( action_flag.EQ.MPP_OVERWR )action_flag = MPP_WRONLY
!perform the OPEN here
          ios = 0
          if( PRESENT(recl) )then
              if( verbose )print '(2(x,a,i6),5(x,a),a,i8)', 'MPP_OPEN: PE=', pe, &
                   'unit=', unit, trim(mpp_file(unit)%name), 'attributes=', trim(acc), trim(for), trim(act), ' RECL=', recl
              open( unit, file=trim(mpp_file(unit)%name), access=acc, form=for, action=act, recl=recl,iostat=ios )
          else     
              if( verbose )print '(2(x,a,i6),6(x,a))',      'MPP_OPEN: PE=', pe, &
                   'unit=', unit, trim(mpp_file(unit)%name), 'attributes=', trim(acc), trim(for), trim(pos), trim(act)
              open( unit, file=trim(mpp_file(unit)%name), access=acc, form=for, action=act, position=pos, iostat=ios)
          end if   
!check if OPEN worked
          inquire( unit,OPENED=mpp_file(unit)%opened )
          if (ios/=0) then
              if (PRESENT(iostat)) then
                  iostat=ios
                  call mpp_error( WARNING, 'MPP_OPEN: error in OPEN for '//trim(mpp_file(unit)%name)//'.' )
                  call mpp_clock_end(mpp_open_clock)
                  return
              else
                  call mpp_error( FATAL, 'MPP_OPEN: error in OPEN for '//trim(mpp_file(unit)%name)//'.' )
              endif
          endif   
          if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_OPEN: error in OPEN() statement.' )
      end if
      mpp_file(unit)%action = action_flag
      mpp_file(unit)%format = form_flag
      mpp_file(unit)%access = access_flag
      if( PRESENT(nohdrs) )mpp_file(unit)%nohdrs = nohdrs
      
      if( action_flag.EQ.MPP_WRONLY )then
          if( form_flag.NE.MPP_NETCDF .AND. access_flag.EQ.MPP_DIRECT )call mpp_write_meta( unit, 'record_length', ival=recl )
!actual file name
          call mpp_write_meta( unit, 'filename', cval=mpp_file(unit)%name)
!MPP_IO package version
!          call mpp_write_meta( unit, 'MPP_IO_VERSION', cval=trim(version) )
!filecount for multifileset.
          if( threading_flag.EQ.MPP_MULTI .AND. fileset_flag.EQ.MPP_MULTI ) then
             if(present(domain)) then
                nfiles = io_layout(1)*io_layout(2)
             else
                nfiles = mpp_npes()
             endif
             call mpp_write_meta( unit, 'NumFilesInSet', ival=nfiles)
          end if
      end if   

      if(ASSOCIATED(io_domain) )io_domain=>NULL()
      call mpp_clock_end(mpp_open_clock)              
      return
    end subroutine mpp_open


! <SUBROUTINE NAME="mpp_close">
!   <OVERVIEW>
!     Close an open file.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Closes the open file on <TT>unit</TT>. Clears the
!     <TT>type(filetype)</TT> object <TT>mpp_file(unit)</TT> making it
!     available for reuse.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_close( unit, action )
!   </TEMPLATE>
!   <IN NAME="unit" TYPE="integer"> </IN>
!   <IN NAME="action" TYPE="integer"> </IN>
! </SUBROUTINE>

    subroutine mpp_close( unit, action )
      integer, intent(in) :: unit
      integer, intent(in), optional :: action
      character(len=8) :: status
      logical :: collect

      call mpp_clock_begin(mpp_close_clock)      
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOSE: must first call mpp_io_init.' )
      if( unit.EQ.NULLUNIT ) then
         call mpp_clock_end(mpp_close_clock)
         return !nothing was actually opened on this unit
      endif
!action on close
      status = 'KEEP'
!collect is supposed to launch the post-processing collector tool for multi-fileset
      collect = .FALSE.
      if( PRESENT(action) )then
          if( action.EQ.MPP_DELETE )then
              if( pe.EQ.mpp_root_pe() .OR. mpp_file(unit)%fileset.EQ.MPP_MULTI )status = 'DELETE'
          else if( action.EQ.MPP_COLLECT )then
              collect = .FALSE.         !should be TRUE but this is not yet ready
              call mpp_error( WARNING, 'MPP_CLOSE: the COLLECT operation is not yet implemented.' )
          else
              call mpp_error( FATAL, 'MPP_CLOSE: action must be one of MPP_DELETE or MPP_COLLECT.' )
          end if
      end if  
      if( mpp_file(unit)%fileset.NE.MPP_MULTI )collect = .FALSE.
      if( mpp_file(unit)%opened) then
         if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
            error = NF_CLOSE(mpp_file(unit)%ncid); call netcdf_err( error, mpp_file(unit) )
#endif    
         else
            close(unit,status=status)
         end if
      endif
#ifdef SGICRAY
!this line deleted: since the FILENV is a shared file, this might cause a problem in
! multi-threaded I/O if one PE does assign -R before another one has opened it.
!      call ASSIGN( 'assign -R f:'//trim(mpp_file(unit)%name), error )
#endif
      mpp_file(unit)%name = ' '
      mpp_file(unit)%action    = -1
      mpp_file(unit)%format    = -1
      mpp_file(unit)%access    = -1
      mpp_file(unit)%threading = -1
      mpp_file(unit)%fileset   = -1
      mpp_file(unit)%record    = -1
      mpp_file(unit)%ncid      = -1
      mpp_file(unit)%opened = .FALSE.
      mpp_file(unit)%initialized = .FALSE.
      mpp_file(unit)%id = -1
      mpp_file(unit)%ndim = -1
      mpp_file(unit)%nvar = -1
      mpp_file(unit)%time_level = 0
      mpp_file(unit)%time = NULLTIME
      mpp_file(unit)%valid = .false.
      mpp_file(unit)%io_domain_exist = .false.
      mpp_file(unit)%write_on_this_pe = .false.
      if(Associated(mpp_file(unit)%domain)) mpp_file(unit)%domain=>NULL() 

      call mpp_clock_end(mpp_close_clock)
      return
    end subroutine mpp_close
    
    
    subroutine file_size(fsize, fname, size)
    
    character(len=12), intent(in)    ::fsize
    character(len=128)		      ::filesize
    character(len=128), intent(in),optional :: fname
    character(len=128) :: filename
    integer*4 :: fstat(13)
    integer :: length
    character(len=16)			::number
    integer,intent(OUT) :: size 
    integer*4 ::ierr, stat
    integer	:: tend
    logical :: there
    
   	size = 0 
  
    	filesize = fsize
   	
    	length = len(trim(fsize))
	tend = length - 1
	
	
    	 
	if (filesize .EQ. 'file') then
		filename = trim(fname)
		INQUIRE( FILE=filename, EXIST=THERE )
		if (THERE) then
			ierr = stat(filename, fstat)
			if (ierr .EQ. 0) then
				size = fstat(8)
			else
				size = 0
			end if
		end if
	elseif ((filesize(length:length) >= 'a' .AND. fsize(length:length) <= 'z') .OR. (filesize(length:length) >= 'A' .AND. fsize(length:length) <= 'Z')) then
		number = filesize(1:tend)
		READ(number, FMT='(I5)') size
		if (filesize(length:length) >= 'a' .AND. fsize(length:length) <= 'z') then
			filesize(length:length) = ACHAR ( ICHAR (filesize(length:length)) - 32)
		end if
		if ( filesize(length:length) .EQ. 'K') then
			size = size*1024
		elseif ( filesize(length:length) .EQ. 'M') then
			size = (size*1024)*1024
		elseif ( filesize(length:length) .EQ. 'G') then
			size = (((size*1024)*1024)*1024)
		else
			size = size
		end if
	else   
		READ(number, FMT='(I5)') size
	endif
	
	if (size .eq. 0) then
		size = 65536
	endif
	
	return
    
    end subroutine file_size



