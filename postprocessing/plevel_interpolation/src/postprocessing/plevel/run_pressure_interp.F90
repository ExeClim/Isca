
program run_pressure_interp

use netcdf

use    plev_constants_mod, only: GRAV, RDGAS, RVGAS

use pinterp_utilities_mod, only: copy_axis, copy_variable,  &
                                 set_verbose_level,         &
                                 define_new_variable,       &
                                 error_handler
use pressure_interp_mod


implicit none

!-----------------------------------------------------------------------

integer, parameter :: MAX_PLEVS = 200
integer, parameter :: MAX_FIELD_NAMES = 300
integer, parameter :: MAX_AXIS_NAMES  =  50
integer, parameter  :: header_buffer_val = 16384

 character(len=256)           ::  in_file_name, out_file_name
!character(len=NF90_MAX_NAME) ::  in_file_name, out_file_name
  real, dimension(MAX_PLEVS) :: pout
character(len=NF90_MAX_NAME) ::    ps_name = 'ps',     &
                                zsurf_name = 'zsurf',  &
                                 temp_name = 'temp',   &
                                sphum_name = 'sphum',  &
                                   bk_name = 'bk',     &
                                   pk_name = 'pk',     &
                                  res_name = 'res',    &
                                 hght_name = 'hght',   &
                                  slp_name = 'slp'

 character(len=NF90_MAX_NAME) :: pres_name_out = 'pfull'
!character(len=NF90_MAX_NAME) :: time_name_out = 'time'   ! not used

 character(len=NF90_MAX_NAME), dimension(MAX_FIELD_NAMES) :: field_names

type field_type
   character(len=NF90_MAX_NAME) :: name
   integer                      :: varid_in, varid_out, ndim, xtype
   integer                      :: dimids(NF90_MAX_VAR_DIMS), &
                                   dimlen(NF90_MAX_VAR_DIMS)
   integer                      :: y_axis, z_axis
   real                         :: pack(2), missval_in, unpkd_missval_in, missval_out, unpkd_missval_out
   logical                      :: static, skip, do_avg, do_pack, do_miss_in, do_miss_out
end type

type axis_type
   character(len=NF90_MAX_NAME) :: name_in, name_out
   integer                      :: varid_in, varid_out
   integer                      :: dimlen, start(1)
end type

type(field_type) :: Fields(MAX_FIELD_NAMES)
type(axis_type)  :: Axes(MAX_AXIS_NAMES)

logical :: mask_extrap = .false.
logical :: unpack = .false.  ! not implemented
integer :: verbose     = 1

logical :: do_all_3d_fields = .false.
logical :: do_all_fields = .false.
logical :: use_default_missing_value = .false.

logical :: allow_zero_topog = .false.
logical :: allow_zero_sphum = .false.

integer :: time_beg=1, time_end=0, time_inc=1
integer :: blksz

!#######################################################################
!         -------- NAMELIST -------
!
!  in_file_name = input file name/path
!               must contain atmospheric dynamics fields (required)
!
!  out_file_name = output file name/path
!
!  field_names = additional fields (up to 40) for output
!                (dynamics files are searched for these names)
!
!         pout = output pressure levels (max 200), in pascals
!
!--- specify these names for the following output fields ---
!         sphum_name = specific humidity
!--------------------------------------------------------
!
!  mask_extrap = set points below the bottom level to the missing value
!
!  use_default_missing_value = use the NetCDF default for missing/fill value
!                              for all fields
!
!   verbose   = integer verbose flag
!
!-----------------------------------------------------------------------
!
!  *** NOTES ***
!
!    1)  currently only ONE input atmospheric file and ONE output file
!        are allowed
!
!    2)  the namelist file must be passed to the program by
!        redirection: e.g., a.out < namelist
!
!#######################################################################

namelist /names/ temp_name, sphum_name, hght_name, slp_name, &
                 ps_name, zsurf_name, bk_name, pk_name 

namelist /input/ in_file_name, out_file_name,  field_names,  pout, &
                 do_all_fields,  do_all_3d_fields,                 &
                 mask_extrap,                          &
                 time_beg, time_end, time_inc,                     &
                 allow_zero_topog, allow_zero_sphum,               &
                 use_default_missing_value, verbose

!----- axis names -----

!character(len=NF90_MAX_NAME) :: xmass_in = 'lon',  &
!                                ymass_in = 'lat',  &
!character(len=NF90_MAX_NAME) :: pfull_in = 'pfull', &
!                                phalf_in = 'phalf', &
!                                time_in  = 'time'


character(len=NF90_MAX_NAME) :: xmass_in, ymass_in
character(len=NF90_MAX_NAME) :: pfull_in, time_in

character(len=NF90_MAX_NAME) :: xmass_out, ymass_out
character(len=NF90_MAX_NAME) :: edge_name
character(len=NF90_MAX_NAME) :: attname
character(len=NF90_MAX_NAME) :: file_title

integer :: dimids_mass(4)
!-----------------------------------------------------------------------

integer :: ncid_in, ncid_out
type (pres_interp_type) :: Interp(2)


character(len=NF90_MAX_NAME) :: version
character(len=NF90_MAX_NAME) :: time_avg_info
character(len=NF90_MAX_NAME) :: cell_methods
character(len=NF90_MAX_NAME) :: name_time_bounds, name_time_average

logical    :: do_time_average, do_time_bounds
character(len=NF90_MAX_NAME) :: avg_name = 'average'
character(len=8) :: avg_subnames(3) = (/ '_T1', '_T2', '_DT' /)

logical    :: temp_present, sphum_present
logical    :: exists

logical :: do_temp  = .false.
logical :: do_sphum = .false.
logical :: do_zsurf = .false.
logical :: do_hght  = .false.
logical :: do_slp   = .false.

logical :: need_zsurf, need_temp, need_sphum

integer :: varid_res,      varid_pout
integer :: varid_bk,       varid_pk,     &
           varid_temp_in,  varid_sphum_in, &
           varid_zsurf_in, varid_ps_in

integer :: varid_hght_out, varid_slp_out

integer :: varid_time_in,        varid_time_out
integer :: varid_time_bounds_in, varid_time_bounds_out
integer :: varid_tavg_in(3),     varid_tavg_out(3)

integer :: index_ps, index_zsurf, index_temp, index_sphum

integer :: istat, ivar, itime, otime, i, n, kv, attnum, numatts, recdim,  &
           nlon, nlat, nlev, ndim, dimids(NF90_MAX_VAR_DIMS), dimid, naxes, dimid_tbnds
integer :: nlev_out, start(4), kount(4), ntimes, num_field_names
integer :: itime_beg, itime_end, itime_inc
integer :: isg, ieg, jsg, jeg, is, ie, js, je
real    :: mval_hght, mval_slp
real(8) :: time, time_bnds(2)
logical :: done, do_phalf
integer :: idim, nlen
integer :: unit, ierr, io
real, allocatable :: axisdata(:)

real, allocatable, dimension(:)     :: bk, pk
real, allocatable, dimension(:,:)   :: res
real, allocatable, dimension(:,:,:) :: zsurf, ps

real, parameter :: tlapse = 6.5e-3
real, parameter :: tref   = 288.
real, parameter :: ginv = 1./GRAV
real, parameter :: gorg = GRAV / (RDGAS*tlapse)
real, parameter :: mrgog = -1./gorg
real, parameter :: d608 = (RVGAS-RDGAS)/RDGAS
real, parameter :: pref = 101325.
real, parameter :: rgog = RDGAS*tlapse/GRAV

!#######################################################################

!--------- Initalize blksz
      blksz = 65536

!--------- initialize namelist variable and then read namelist ---------

      do i=1,len(out_file_name);  out_file_name(i:i) = ' '; enddo
      do i=1,len( in_file_name);   in_file_name(i:i) = ' '; enddo
      pout = 0.0

      do ivar = 1, MAX_FIELD_NAMES
         do i=1,NF90_MAX_NAME; field_names(ivar)(i:i) = ' '; enddo
      enddo

      inquire (file='plev.input.nml',exist=exists)
      if (exists) then
         open (10, file='plev.input.nml')
         read (10, nml=input, end=1, err=1)
      1  close (10)
         !write (*,input)
      endif

         do i = 1, MAX_PLEVS
            if (pout(i) < 1.e-6) exit
            nlev_out = i
         enddo
 
!     ---- count the number of field names ----
         num_field_names = 0
      do ivar = 1, MAX_FIELD_NAMES
         if (field_names(ivar)(1:1) == ' ') exit
         num_field_names = ivar
      enddo

      call set_verbose_level (verbose)

  ! create version string (may replace with CVS $Id: run_pressure_interp.F90,v 18.0 2010/03/03 00:00:45 fms Exp $)
    version = 'pressure level interpolator, version 3.0'
    if (precision(ps) == precision(time)) then
        version = trim(version)//', precision=double'
    else
        version = trim(version)//', precision=float'
    endif

!-----------------------------------------------------------------------
!----------------------- loop through files ----------------------------
!-----------------------------------------------------------------------
!   --- if any file name has no length then terminate the program ---
      
     if (  in_file_name(1:1) == ' ' .or.   &
          out_file_name(1:1) == ' ' )  call error_handler  &
           ('must specifiy atmospheric and output file names')

!-----------------------------------------------------------------------
!-----------------  get input file information -------------------------

    !--- opening input files ---
     istat = NF90_OPEN (trim(in_file_name), NF90_NOWRITE, ncid_in, blksz)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

    !--- get the time axis identifier ---
     istat = NF90_INQUIRE ( ncid_in, nattributes=numatts, unlimiteddimid=recdim )
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

!    --------- set flags for special computed fields --------

     if (verbose > 1) print *, 'num_field_names=',num_field_names
     do i = 1, num_field_names
        if (verbose > 1) print *, 'field_name=',trim(field_names(i))
        if (field_names(i) ==  hght_name) then
            do_hght  = .true.; cycle; endif
        if (field_names(i) ==   slp_name) then
            do_slp   = .true.; cycle; endif
     enddo

!    ------- check to see if required fields are in atmos file -------

     call check_fields ( ncid_in )

!-----------------------------------------------------------------------
!--- get the size of the vertical axis ---
     istat = NF90_INQUIRE_VARIABLE (ncid_in, varid_bk, ndims=ndim, dimids=dimids)
     if (istat /= NF90_NOERR) call error_handler ('getting dimids for '//trim(bk_name), ncode=istat)
     istat = NF90_INQUIRE_DIMENSION (ncid_in, dimids(1), len=nlev)
     if (istat /= NF90_NOERR) call error_handler ('getting dim length for '//trim(bk_name), ncode=istat)
     nlev = nlev-1

    !--- get the names of the input horizontal mass axes ---
     istat = NF90_INQUIRE_VARIABLE (ncid_in, varid_ps_in, ndims=ndim, dimids=dimids)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

     if (ndim /= 3) call error_handler ('variable '//trim(ps_name)//&
                                    &' must have three dimensions: x,y,t')

    !--- size and name of input horizontal mass axes ---
     istat = NF90_INQUIRE_DIMENSION (ncid_in, dimids(1), name=xmass_in, len=nlon)
     if (istat /= NF90_NOERR) call error_handler ('getting name and length of first mass dimension', ncode=istat)
     istat = NF90_INQUIRE_DIMENSION (ncid_in, dimids(2), name=ymass_in, len=nlat)
     if (istat /= NF90_NOERR) call error_handler ('getting name and length of second mass dimension', ncode=istat)
     istat = NF90_INQUIRE_DIMENSION (ncid_in, dimids(3), name=time_in)
     if (istat /= NF90_NOERR) call error_handler ('getting name of last mass dimension (time)', ncode=istat)

!-----------------------------------------------------------------------
!--- set output fields using all input fields ---

     if (do_all_fields .or. do_all_3d_fields) then
         call set_all_field_names ( ncid_in, field_names, num_field_names )
     endif

!--- error: if there were no requested fields ---
     if (num_field_names == 0) call error_handler &
                               ('no field names have been selected')

!--- set flags for special non-computed fields ---

     if (verbose > 1) print *, 'num_field_names=',num_field_names
     do i = 1, num_field_names
        if (verbose > 1) print *, 'field_name=',trim(field_names(i))
        if (field_names(i) == sphum_name) then
            if (sphum_present) do_sphum = .true.; cycle; endif
        if (field_names(i) ==  temp_name) then
            if (temp_present)  do_temp  = .true.; cycle; endif
        if (field_names(i) == zsurf_name) then
            do_zsurf = .true.; cycle; endif
     enddo

!-----------------------------------------------------------------------

     isg = 1;   ieg = nlon; jsg = 1;   jeg = nlat
     is  = isg; ie  = ieg;  js  = jsg; je  = jeg

!    ---- allocate space ----

     allocate (bk(nlev+1), pk(nlev+1))
     allocate (res  (is:ie,js:je),    &
               zsurf(is:ie,js:je,1),  &
               ps   (is:ie,js:je,1)   )

!    ---- read vert coord variables and set up vert coord ----

     istat = NF90_GET_VAR (ncid_in, varid_bk, bk)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     if (varid_pk > 0) then
        istat = NF90_GET_VAR (ncid_in, varid_pk, pk)
        if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     else
        pk = 0.0
     endif

!-----------------------------------------------------------------------
!----------------- output file header info -----------------------------

     istat = NF90_CREATE (trim(out_file_name), NF90_CLOBBER, ncid_out, blksz)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

   ! copy global attributes
     do attnum = 1, numatts
        ! get name
        istat = NF90_INQ_ATTNAME (ncid_in, NF90_GLOBAL, attnum, attname)
        if (istat /= NF90_NOERR) call error_handler (ncode=istat)
        ! change global filename attribute
        if (trim(attname) == 'filename') then
            istat = NF90_PUT_ATT (ncid_out, NF90_GLOBAL, 'filename', trim(out_file_name))
            if (istat /= NF90_NOERR) call error_handler (ncode=istat)
        ! append to history attribute
        else if (trim(attname) == 'history') then
           !istat = NF90_INQUIRE_ATTRIBUTE (ncid_in, NF90_GLOBAL, trim(attname), len=nc)
           !istat = NF90_GET_ATT (ncid_in, NF90_GLOBAL, 'history', history)
        else    
            ! copy all others
            istat = NF90_COPY_ATT (ncid_in, NF90_GLOBAL, attname, ncid_out, NF90_GLOBAL)
            if (istat /= NF90_NOERR) call error_handler ('copying global attributes', &
                                                          ncode=istat)
        endif   
     enddo
     ! add comment attribute
     istat = NF90_PUT_ATT (ncid_out, NF90_GLOBAL, 'comment', trim(version))
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

!-----------------------------------------------------------------------
!------------------ generate axes for output file ----------------------

    xmass_out = xmass_in;  ymass_out = ymass_in

    !---- copy mass axes ----
     naxes = 1
     call copy_axis_information (ncid_in, xmass_in, ncid_out, (/isg,ieg,is,ie/), Axes(naxes), &
                                 dimid_out=dimids_mass(1), name_out=xmass_out, name_edges=edge_name)
         ! copy edge axis if available
         if (edge_name(1:1) .ne. ' ') then
             naxes = naxes + 1
             call copy_axis_information (ncid_in, edge_name, ncid_out, (/isg,ieg+1,is,ie+1/), Axes(naxes))
         endif

     naxes = naxes + 1
     call copy_axis_information (ncid_in, ymass_in, ncid_out, (/jsg,jeg,js,je/), Axes(naxes), &
                                 dimid_out=dimids_mass(2), name_out=ymass_out, name_edges=edge_name)
         ! copy edge axis if available
         if (edge_name(1:1) .ne. ' ') then
             naxes = naxes + 1
             call copy_axis_information (ncid_in, edge_name, ncid_out, (/jsg,jeg+1,js,je+1/), Axes(naxes))
         endif

!    -------- vertical axis (convert to mb) -------

   ! define new output pressure dimension and variable
     istat = NF90_DEF_DIM (ncid_out, trim(pres_name_out), nlev_out, dimids_mass(3))
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     istat = NF90_DEF_VAR (ncid_out, trim(pres_name_out), NF90_REAL4, (/dimids_mass(3)/), varid_pout)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
   ! attributes for this axis
     istat = NF90_PUT_ATT (ncid_out, varid_pout, 'units', 'hPa')
     if (istat /= NF90_NOERR) call error_handler ('putting attributes for pout dimension', ncode=istat)
     istat = NF90_PUT_ATT (ncid_out, varid_pout, 'long_name', 'pressure')
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     istat = NF90_PUT_ATT (ncid_out, varid_pout, 'cartesian_axis', 'Z')
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     istat = NF90_PUT_ATT (ncid_out, varid_pout, 'axis', 'Z') ! for CF compliance
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     istat = NF90_PUT_ATT (ncid_out, varid_pout, 'positive', 'down')
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

   ! copy the time axis
   ! get variable identifiers and length
     call copy_axis (ncid_in, time_in, ncid_out, dimids_mass(4), &
                     varid_time_in, varid_time_out, ntimes)

   ! set up time axis indexing
     itime_beg = time_beg
     itime_end = time_end
     itime_inc = time_inc
     if (itime_end == 0 .or. itime_end > ntimes) itime_end = ntimes
     if (itime_beg <= 0) itime_beg = 1
     if (itime_inc <= 0) itime_inc = 1
    !if (itime_beg < itime_end) itime_end = itime_beg  ! process at least one time level


   ! check for time averaging "bounds"
   ! make sure axis is copied (time axis already copied)
     do i=1,NF90_MAX_NAME;  name_time_bounds(i:i) = ' '; enddo
     do i=1,NF90_MAX_NAME;  cell_methods(i:i) = ' '    ; enddo
     istat = NF90_GET_ATT (ncid_in, varid_time_in, 'bounds', name_time_bounds)
     if (istat /= NF90_NOERR) then
         istat = NF90_GET_ATT (ncid_in, varid_time_in, 'climatology', name_time_bounds)
     endif
     if (istat == NF90_NOERR) then
         if (verbose > 1) print *, 'name_time_bounds=',trim(name_time_bounds)
         istat = NF90_INQ_VARID (ncid_in, trim(name_time_bounds), varid_time_bounds_in)
         if (istat /= NF90_NOERR) call error_handler ('getting varid for time bounds', ncode=istat)
         istat = NF90_INQUIRE_VARIABLE (ncid_in, varid_time_bounds_in, ndims=ndim, dimids=dimids)
         if (istat /= NF90_NOERR) call error_handler ('getting dimensions for time bounds', ncode=istat)
         if (ndim == 1 .or. ndim == 2) then
             istat = NF90_INQUIRE_DIMENSION (ncid_in, dimids(1), name=edge_name)
             if (istat == NF90_NOERR) then
                 naxes = naxes + 1
                 call copy_axis_information (ncid_in, edge_name, ncid_out, (/1,1,1,1/), &
                                             Axes(naxes), dimid_out=dimid_tbnds)
             endif
         else
             call error_handler ('number of dimensions not consistent for '//trim(name_time_bounds))
         endif
         ! get cell methods for new variables
         istat = NF90_GET_ATT (ncid_in, varid_ps_in, 'cell_methods', cell_methods)
         do_time_bounds = .true.
     else
         do_time_bounds = .false.
     endif

!-----------------------------------------------------------------------
!-------------- define new variables for output file -----------------

!  ---- geopotential height ----

if (do_hght) then
   mval_hght = -1000.
   if (use_default_missing_value) mval_hght = NF90_FILL_REAL4

   call define_new_variable ( ncid_out, 'hght', NF90_REAL4, dimids_mass, varid_hght_out, &
                              units='m', long_name='height', missing_value=mval_hght,    &
                              time_avg_info=time_avg_info, cell_methods=cell_methods)
   
endif

!  ---- sea level pressure ----

if (do_slp) then
   mval_slp = -1.
   if (use_default_missing_value) mval_slp = NF90_FILL_REAL4

   call define_new_variable ( ncid_out, 'slp', NF90_REAL4,                      &
                              (/dimids_mass(1),dimids_mass(2),dimids_mass(4)/), &
                              varid_slp_out, units='hPa', long_name='sea level pressure', &
                              missing_value=mval_slp, time_avg_info=time_avg_info,   &
                              cell_methods=cell_methods)
endif

!-----------------------------------------------------------------------
!---- copy variable info for other fields ------

    Fields(:)%varid_in = 0
    Fields(:)%y_axis = 0
    Fields(:)%z_axis = 0
    Fields(:)%static = .false.   ! only allowing 2d static fields
    do_phalf = .false.

    do i = 1, num_field_names

       call copy_variable_information ( field_names(i), ncid_in, ncid_out, Fields(i) )

       ! save index for special fields
       if (trim(field_names(i)) == trim(   ps_name)) index_ps    = i
       if (trim(field_names(i)) == trim(zsurf_name)) index_zsurf = i
       if (trim(field_names(i)) == trim( temp_name)) index_temp  = i
       if (trim(field_names(i)) == trim(sphum_name)) index_sphum = i

       ! check valid fields for half levels then go to end of loop
       if (.not.Fields(i)%skip) then
           if (Fields(i)%z_axis == 2) do_phalf = .true.
           cycle
       endif

       if (trim(field_names(i)) == trim( hght_name)) cycle
       if (trim(field_names(i)) == trim(  slp_name)) cycle
       print *, 'WARNING: unable to process field '//trim(field_names(i))
    enddo

!-----------------------------------------------------------------------
!---- input variable only are needed to process some field ----
!---- check to see if these variables are in the field_names list ----

      if (.not.field_name_present(ps_name,field_names,num_field_names)) then
          num_field_names = num_field_names + 1
          index_ps = num_field_names
          call get_input_variable_info ( ps_name, ncid_in, Fields(num_field_names) )
      endif

      if (need_zsurf) then
          if (.not.field_name_present(zsurf_name,field_names,num_field_names)) then
              num_field_names = num_field_names + 1
              index_zsurf = num_field_names
              call get_input_variable_info ( zsurf_name, ncid_in, Fields(num_field_names) )
          endif
      endif

      if (need_temp) then
          if (.not.field_name_present(temp_name,field_names,num_field_names)) then
              num_field_names = num_field_names + 1
              index_temp = num_field_names
              call get_input_variable_info ( temp_name, ncid_in, Fields(num_field_names) )
          endif
      endif

      if (need_sphum) then
          if (.not.field_name_present(sphum_name,field_names,num_field_names)) then
              num_field_names = num_field_names + 1
              index_sphum = num_field_names
              call get_input_variable_info ( sphum_name, ncid_in, Fields(num_field_names) )
          endif
      endif

!-----------------------------------------------------------------------
! ---- copy variable info for time averaging -----

      if (do_time_bounds) then
          call copy_variable ( ncid_in, trim(name_time_bounds), ncid_out,    &
                               (/dimid_tbnds,dimids_mass(4)/), varid_time_bounds_out )
      endif

      if (do_time_average) then
          do i = 1, 3
             name_time_average = trim(avg_name)//trim(avg_subnames(i))
             istat = NF90_INQ_VARID (ncid_in, trim(name_time_average), varid_tavg_in(i))
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)
             call copy_variable ( ncid_in, trim(name_time_average), ncid_out,    &
                                  (/dimids_mass(4)/), varid_tavg_out(i) )
          enddo
       endif
          
!-----------------------------------------------------------------------
! ---- write header info for output file ----

      istat = NF90_ENDDEF (ncid_out,header_buffer_val,4,0,4)
      if (istat /= NF90_NOERR) call error_handler (ncode=istat)

!-----------------------------------------------------------------------
!------------ read and write static axis data -------------

       do idim = 1, naxes
          if (Axes(idim)%varid_in < 0) cycle
          nlen = Axes(idim)%dimlen
          allocate (axisdata(nlen))
          istat = NF90_GET_VAR (ncid_in, Axes(idim)%varid_in, axisdata, start=Axes(idim)%start)
          if (istat /= NF90_NOERR) call error_handler ('getting axis data for var='//trim(Axes(idim)%name_in), ncode=istat)
          istat = NF90_PUT_VAR (ncid_out, Axes(idim)%varid_out, axisdata)
          if (istat /= NF90_NOERR) call error_handler ('putting axis data for var='//trim(Axes(idim)%name_out), ncode=istat)
          deallocate (axisdata)
       enddo
        ! output pressure levels in millibars
          istat = NF90_PUT_VAR (ncid_out, varid_pout, pout(1:nlev_out)*0.01)
          if (istat /= NF90_NOERR) call error_handler ('putting axis data for var='//trim(pres_name_out), ncode=istat)

!-----------------------------------------------------------------------
!------------ first time: read constant fields -------------------------

       start(1:3) = (/ is, js, 1 /)
       istat = NF90_INQ_VARID (ncid_in, trim(res_name), varid_res)
       if ( istat == NF90_NOERR ) then
          istat = NF90_GET_VAR (ncid_in, varid_res, res, start=start(1:2))
          if (istat /= NF90_NOERR) call error_handler (ncode=istat)
       else
          res = 1.0
       endif
       if (need_zsurf .or. do_zsurf) then
          call read_variable (ncid_in, Fields(index_zsurf), zsurf, start=start(1:2))
       else
          zsurf = 0.0
       endif
       if (do_zsurf) then
          call write_variable (ncid_out, Fields(index_zsurf), zsurf)
       endif
!      --- eta coord fix ----
       if (count(res > 1.0001) > 0) zsurf = 0.0

!-----------------------------------------------------------------------
!      ---- read/write static fields ----

    do i = 1, num_field_names
       if ( Fields(i)%skip )         cycle
       if ( .not.Fields(i)%static )  cycle
       if (trim(Fields(i)%name) == trim(zsurf_name)) cycle
       call read_write_static (ncid_in, ncid_out, is, ie, js, je, Fields(i))
    enddo

!-----------------------------------------------------------------------
!-------------------- loop through time axis ---------------------------
!-----------------------------------------------------------------------

      otime = 0
      do itime = itime_beg, itime_end, itime_inc
          otime = otime+1
          if ( verbose > 0 ) then
              if (itime == otime) then
                 print '(a,i8)', 'itime =', itime
              else
                 print '(a,2i8)', 'itime, otime =', itime, otime
              endif
          endif

!---------------- read/write time level --------------

          istat =  NF90_GET_VAR (ncid_in, varid_time_in, time, start=(/itime/))
          if (istat /= NF90_NOERR) call error_handler ('getting time coord value', &
                                                       ncode=istat)

          istat =  NF90_PUT_VAR (ncid_out, varid_time_out, time, start=(/otime/))
          if (istat /= NF90_NOERR) call error_handler ('putting time coord value', &
                                                       ncode=istat)

!-----------------------------------------------------------------------
!-------------- interpolation setup ------------------------------------

     !---- surface pressure ----

          start(1:3) = (/ is, js, itime /)
          call read_variable (ncid_in, Fields(index_ps), ps, start=start(1:3))
          ps(is:ie,js:je,1) = ps(is:ie,js:je,1) * res

!-----------------------------------------------------------------------
!----------------- mass field interpolation ---------------------

!     ---- interpolation initialization for mass fields ------

     call interp_field_setup ( ncid_in, ncid_out, bk, pk,   &
                               is, js, itime, otime,        &
                               ps(is:ie,js:je,1),           &
                               zsurf(:,:,1), pout(1:nlev_out) )

!     ---- interpolation of mass fields ----

     do i = 1, num_field_names
        if ( Fields(i)%skip   )    cycle
        if ( Fields(i)%static )    cycle
        if ( Fields(i)%y_axis /= nlat ) cycle

        kv = max(1,Fields(i)%z_axis)
        if ( verbose > 1 ) then
           print *, 'Interp field = ',trim(Fields(i)%name)
        endif

!     ---- interpolation of mass fields ----
        call interp_field ( ncid_in, ncid_out, Interp(kv),  &
                            is, js, itime, otime, Fields(i) )

     enddo

!     ---- free up allocated space ----
                    call pres_interp_free (Interp(1))
      if (do_phalf) call pres_interp_free (Interp(2))

!-----------------------------------------------------------------------
!   ---- copy time average data ----

    if (do_time_bounds) then
        istat = NF90_GET_VAR (ncid_in,  varid_time_bounds_in,  time_bnds(:), start=(/1,itime/))
        if (istat /= NF90_NOERR) call error_handler ('getting time_bounds', ncode=istat)
        istat = NF90_PUT_VAR (ncid_out, varid_time_bounds_out, time_bnds(:), start=(/1,otime/))
        if (istat /= NF90_NOERR) call error_handler ('putting time_bounds', ncode=istat)
    endif

    if (do_time_average) then
        do i = 1, 3
           istat = NF90_GET_VAR (ncid_in,  varid_tavg_in(i),  time_bnds(1), start=(/itime/))
           if (istat /= NF90_NOERR) call error_handler ('getting tavg_info', ncode=istat)
           istat = NF90_PUT_VAR (ncid_out, varid_tavg_out(i), time_bnds(1), start=(/otime/))
           if (istat /= NF90_NOERR) call error_handler ('putting tavg_info', ncode=istat)
        enddo
    endif

!-----------------------------------------------------------------------

   enddo  ! end of time loop

!-----------------------------------------------------------------------

       istat = NF90_CLOSE (ncid_in)
       if (istat /= NF90_NOERR) call error_handler (ncode=istat)
       istat = NF90_CLOSE (ncid_out)
       if (istat /= NF90_NOERR) call error_handler (ncode=istat)

!-----------------------------------------------------------------------

contains

!#######################################################################

   subroutine interp_field_setup ( ncid_in, ncid_out, bk, pk, &
                                   is, js, itime, otime, ps,  &
                                   zsurf, pout )

   integer, intent(in) :: ncid_in, ncid_out
   real,    intent(in) :: bk(:), pk(:)
   integer, intent(in) :: is, js, itime, otime
   real   , intent(in) :: ps(:,:), zsurf(:,:), pout(:)

 real, dimension(size(ps,1),size(ps,2),size(bk(:))-1) :: &
            temp_in, sphum_in, mask_in, pfull, hght_in
 real, dimension(size(ps,1),size(ps,2),size(bk(:)))   :: phalf
 real, dimension(size(ps,1),size(ps,2),size(pout)) ::  &
                         temp_out, hght_out, sphum_out
 real   , dimension(size(ps,1),size(ps,2)) :: zbot, tbot, pbot, slp
 integer, dimension(size(ps,1),size(ps,2)) :: kbot
   character(len=NF90_MAX_NAME) :: name
   integer :: start(4), kount(4), ndim, ivar, nlev, kb, kr, i, j, k
   real    :: sig


      do k = 1, size(bk(:))
         phalf(:,:,k) = pk(k) + bk(k)*ps(:,:)
      enddo
      call compute_pres_full (pk, bk, phalf, pfull)

  !---- read in temperature data ----
   if (do_temp .or. need_temp) then
       start = (/ is, js, 1, itime /)
       call read_variable (ncid_in, Fields(index_temp), temp_in, start)
      !---- mask ----
       call set_up_input_mask (ncid_in, Fields(index_temp), temp_in, mask_in, kbot)
   endif

  !---- read in specific humidity data ----
   if (do_sphum .or. need_sphum) then
       if (sphum_present) then
           start = (/ is, js, 1, itime /)
           call read_variable (ncid_in, Fields(index_sphum), sphum_in, start)
          !---- mask ----
           if (.not.do_temp) &
           call set_up_input_mask (ncid_in, Fields(index_sphum), sphum_in, mask_in, kbot)
       else
           sphum_in = 0.0
       endif
   endif

  !---- interpolation setup for standard full pressure levels -----
   if (do_hght) then
       if ( verbose > 1 ) print *, 'Computing field = hght'
      !---- geopotential height (in meters) -----
       call compute_height ( zsurf, temp_in, sphum_in, pfull, phalf, hght_in, mask_in)
      !---- interpolate temp, sphum, hght ----
       if ( verbose > 1 ) then
           print *, 'Interp field = ',trim(Fields(index_temp)%name)
           print *, 'Interp field = ',trim(Fields(index_sphum)%name)
       endif
       Interp(1) = pres_interp_init ( pfull, pout, kbot,   &
                                      tin=temp_in,  tout=temp_out,  &
                                      zin=hght_in,  zout=hght_out,  &
                                      qin=sphum_in, qout=sphum_out, &
                                      use_extrap=.not.mask_extrap   )
   else if ((do_temp.and.do_sphum).and.(.not.do_hght)) then
       if ( verbose > 1 ) then
           print *, 'Interp field = ',trim(Fields(index_temp)%name)
           print *, 'Interp field = ',trim(Fields(index_sphum)%name)
       endif
      !---- interpolate temp, sphum ----
       Interp(1) = pres_interp_init ( pfull, pout, kbot,   &
                                      tin=temp_in,  tout=temp_out,  &
                                      qin=sphum_in, qout=sphum_out, &
                                      use_extrap=.not.mask_extrap   )
   else if ((do_temp).and.(.not.do_hght)) then
       if ( verbose > 1 ) print *, 'Interp field = ',trim(Fields(index_temp)%name)
      !---- interpolate temp ----
       Interp(1) = pres_interp_init ( pfull, pout, kbot,   &
                                      tin=temp_in,  tout=temp_out,  &
                                      use_extrap=.not.mask_extrap   )
   else if ((do_sphum).and.(.not.do_hght)) then
       if ( verbose > 1 ) print *, 'Interp field = ',trim(Fields(index_sphum)%name)
      !---- interpolate sphum ----
       Interp(1) = pres_interp_init ( pfull, pout, kbot,   &
                                      qin=sphum_in, qout=sphum_out, &
                                      use_extrap=.not.mask_extrap   )
   else
      !---- initialize only (need kbot) ----
       start = (/ is, js, 1, itime /)
       call read_and_set_up_mask (ncid_in, start, mask_in, kbot)
       Interp(1) = pres_interp_init ( pfull, pout, kbot,   &
                                      use_extrap=.not.mask_extrap   )
   endif
   
  !---- interpolation setup for standard full pressure levels -----
   if (do_phalf) then
       Interp(2) = pres_interp_init ( phalf, pout, kbot )
   endif

! ---- write out data ----

     !start = (/ is, js,  1, otime /)
      start = (/  1,  1,  1, otime /)

!     ---- temperature ----

   if (do_temp) then
      if ( verbose > 1 ) print *, 'Writing field = ',trim(Fields(index_temp)%name)
      if (Fields(index_temp)%do_miss_out) then
          call mask_data (Fields(index_temp)%unpkd_missval_out, Interp(1)%mask, temp_out)
      endif
      call write_variable (ncid_out, Fields(index_temp), temp_out, start)
   endif

!        ---- specific humidity ----

   if (do_sphum) then
      if ( verbose > 1 ) print *, 'Writing field = ',trim(Fields(index_sphum)%name)
      if (Fields(index_sphum)%do_miss_out) then
          call mask_data (Fields(index_sphum)%unpkd_missval_out, Interp(1)%mask, sphum_out)
      endif
      call write_variable (ncid_out, Fields(index_sphum), sphum_out, start)
   endif

!        ---- height ----

   if (do_hght) then
      if ( verbose > 1 ) print *, 'Writing field = hght'
      call mask_data (mval_hght, Interp(1)%mask, hght_out)
      istat = NF90_PUT_VAR (ncid_out, varid_hght_out, hght_out, start)
      if (istat /= NF90_NOERR) call error_handler ('putting hght data', ncode=istat)
   endif

!        ---- sea level pressure ----

   if (do_slp) then
     if ( verbose > 1 ) print *, 'Computing field = slp'
     do j = 1, size(ps,2)
     do i = 1, size(ps,1)
        kb = kbot(i,j)
        pbot(i,j) = phalf(i,j,kb+1)

        if ( abs(zsurf(i,j)) > 0.0001 ) then

!            ---- get ref level for temp (as in spectral model) ----
             do k = 1, kb
                sig = pfull(i,j,k)/pbot(i,j)
                if ( sig > 0.8 ) then
                     kr = k
                     exit
                endif
             enddo

            !zbot(i,j) = zhalf(i,j,kb+1)
!!!          tbot(i,j) = tref - tlapse*zbot(i,j)           ! original
             tbot(i,j) = temp_in(i,j,kr) * sig ** mrgog    ! spectral
             slp(i,j) = 0.01 * pbot(i,j) *  &
                       ( 1.0 + tlapse * zsurf(i,j) / tbot(i,j) ) ** gorg
                      !( 1.0 + tlapse * zbot(i,j) / tbot(i,j) ) ** gorg
        else
             slp(i,j) = pbot(i,j) * 0.01
        endif
     enddo
     enddo

     if ( verbose > 1 ) print *, 'Writing field = slp'
     start(1:3) = (/ is, js, otime /)
     istat = NF90_PUT_VAR (ncid_out, varid_slp_out, slp, start)
     if (istat /= NF90_NOERR) call error_handler ('putting slp data', ncode=istat)
   endif


   end subroutine interp_field_setup

!#######################################################################

   subroutine interp_field ( ncid_in, ncid_out, Interp, is, js, itime, otime, Field )

   integer               , intent(in) :: ncid_in, ncid_out
   type(pres_interp_type), intent(in) :: Interp
   integer               , intent(in) :: is, js, itime, otime
   type(field_type)      , intent(in) :: Field

   real, dimension(Interp%nx,Interp%ny,Interp%kx) :: data_in
   real, dimension(Interp%nx,Interp%ny,Interp%nz) :: data_out

   integer :: start_in(4), start_out(4), istat

!     ---- skip fields have been processed in the initialization routine
       if (trim(Field%name) == trim(temp_name))  return
       if (trim(Field%name) == trim(sphum_name)) return

!     ---- read in data ----

       start_in (1:3) = (/ is, js, 1 /)
      !start_out(1:3) = (/ is, js, 1 /)
       start_out(1:3) = (/  1,  1, 1 /)
       start_in (Field%ndim) = itime
       start_out(Field%ndim) = otime

!     ---- check horiz axis lengths ----

      !if ( Field%dimlen(1) /= Interp%nx .or. Field%dimlen(2) /= Interp%ny ) then
      !    print *, 'x,y=',Field%dimlen(1),Interp%nx,Field%dimlen(2),Interp%ny
      !    call error_handler (   &
      !      'horizontal coordinate resolution inconsistent, var = '//trim(Field%name))
      !endif


!  ------ 2-d or 3-d data allowed -----

 select case ( Field%ndim )

   case (4)
!     ---- check vert axis length ----
     if ( Field%dimlen(3) /= Interp%kx ) print *, trim(Field%name)
     if ( Field%dimlen(3) /= Interp%kx ) call error_handler (    &
                'vertical coordinate resolution inconsistent')
     call read_variable (ncid_in, Field, data_in, start_in(1:4))
!    ---- interpolation ----
     call pres_interp (Interp, data_in, data_out)
     if (Field%do_miss_out) &
     call mask_data (Field%unpkd_missval_out, Interp%mask, data_out)
!    ---- write out data ----
     call write_variable (ncid_out, Field, data_out, start_out(1:4))

   case (3)
     call read_variable (ncid_in, Field, data_in(:,:,1:1), start_in(1:3))
     call write_variable (ncid_out, Field, data_in(:,:,1:1), start_out(1:3))

   case default
     print *, 'ERROR: only 2-d or 3-d fields allowed'
     stop

 end select



 end subroutine interp_field

!#######################################################################

  subroutine mask_data (mval, mask, data)

   real            , intent(in)    :: mval
   logical         , intent(in)    :: mask(:,:,:)
   real            , intent(inout) :: data(:,:,:)

     if (.not.mask_extrap) return
     where (mask) data = mval

  end subroutine mask_data

!#######################################################################

 subroutine check_fields ( ncid )

   integer, intent(in) :: ncid

   integer :: id(6), ierr, iwarn, i
   integer :: avg(6), istat(6)
   character(len=NF90_MAX_NAME), dimension(6) :: names

!-----------------------------------------------------------------------

     names = (/ bk_name, pk_name, zsurf_name, ps_name, temp_name, sphum_name /)

     need_zsurf = .false.
     need_temp  = .false.
     need_sphum = .false.

     id = 0

   ! get variable index for all important fields
     do i = 1, 6
       istat(i) = NF90_INQ_VARID (ncid, trim(names(i)), id(i))
     enddo

!------ perform checks ------

       ierr  = 0
       iwarn = 0

    ! minimum: bk, ps (, pk)
      if (istat(1) /= NF90_NOERR) then
          print *, 'ERROR: required field does not exist: ',trim(names(1))
          ierr = 1
      endif
      if (istat(2) /= NF90_NOERR) then
          print *, 'WARNING: required field does not exist: ',trim(names(2))
          iwarn = 1
      endif
      if (istat(4) /= NF90_NOERR) then
          print *, 'ERROR: required field does not exist: ',trim(names(4))
          ierr = 1
      else
          ! check if time average to be done
          avg(4) = NF90_INQUIRE_ATTRIBUTE (ncid, id(4), 'time_avg_info')
      endif
    ! time averaging for all variable must be consistent
    ! initialize all variables to surf pres value
      avg(5:6) = avg(4)

    ! surf height & temperature
      temp_present = .false.;  if (istat(5) == NF90_NOERR ) temp_present = .true.
      if (do_hght .or. do_slp) then
          if (istat(3) /= NF90_NOERR) then
              if (allow_zero_topog) then
                  print *, 'NOTE: field does not exist: ',trim(names(3)),' = 0 assumed'
              else
                  print *, 'ERROR: required field does not exist: ',trim(names(3))
                  ierr = 1
                  do_hght = .false.
                  do_slp  = .false.
              endif
          else
              need_zsurf = .true.
          endif
          if (temp_present) then
              avg(5) = NF90_INQUIRE_ATTRIBUTE (ncid, id(5), 'time_avg_info')
              need_temp = .true.
          else
              print *, 'ERROR: required field does not exist: ',trim(names(5))
              ierr = 1
              do_hght = .false.
              do_slp  = .false.
          endif
      endif

    ! specific humidity
      sphum_present = .false.;  if (istat(6) == NF90_NOERR ) sphum_present = .true.
      if (do_hght) then
          if (sphum_present) then
              avg(6) = NF90_INQUIRE_ATTRIBUTE (ncid, id(6), 'time_avg_info')
              need_sphum = .true.
          else
              if (allow_zero_sphum) then
                  print *, 'NOTE: field does not exist: ',trim(names(6)),' = 0 assumed'
              else
                  print *, 'ERROR: required field does not exist: ',trim(names(6))
                  ierr = 1
                  do_hght = .false.
              endif
          endif
      
      endif


   ! set variable identifiers in the input file
      varid_bk       = id(1)
      varid_pk       = id(2)
      varid_zsurf_in = id(3)
      varid_ps_in    = id(4)
      varid_temp_in  = id(5)
      varid_sphum_in = id(6)

      if ( ierr  /= 0 ) call error_handler ('run_pressure_interp')
      if ( iwarn /= 0 ) print *, 'WARNING: run_pressure_interp'
             
!-----------------------------------------------------------------------
   ! check time averaging flag for consistency
     do i = 5, 6
        if (avg(4) /= avg(i)) call error_handler ('inconsistent time averages')
     enddo

    ! set averaging flag for program ----
      do_time_average = avg(4) == NF90_NOERR

    ! set metadata time average string
      do i=1,NF90_MAX_NAME;  time_avg_info(i:i) = ' '; enddo
      if (do_time_average) then
         time_avg_info = trim(avg_name)//trim(avg_subnames(1))//','// &
                         trim(avg_name)//trim(avg_subnames(2))//','// &
                         trim(avg_name)//trim(avg_subnames(3))
      endif
      
 end subroutine check_fields

!#######################################################################

 subroutine set_all_field_names ( ncid, field_names, num_field_names )

   integer,          intent(in)    :: ncid
   character(len=*), intent(inout) :: field_names(:)
   integer,          intent(inout) :: num_field_names

   character(len=NF90_MAX_NAME) :: name
   integer :: ivar, idim, rdim, nvar, ndim, i, istat
   integer :: dimids(NF90_MAX_VAR_DIMS), axlen(NF90_MAX_VAR_DIMS)
   logical :: static, threed

   if ( num_field_names >= MAX_FIELD_NAMES ) return

   istat = NF90_INQUIRE (ncid, nvariables=nvar, unlimiteddimid=rdim)
   if (istat /= NF90_NOERR) call error_handler (ncode=istat)

     do ivar = 1, nvar

        ! get number and names of dimensions
        istat = NF90_INQUIRE_VARIABLE (ncid, ivar, name=name, ndims=ndim, dimids=dimids)
        if (istat /= NF90_NOERR) call error_handler (ncode=istat)
        ! get dimension lengths
        do idim = 1, ndim
           istat = NF90_INQUIRE_DIMENSION (ncid, dimids(idim), len=axlen(idim))
           if (istat /= NF90_NOERR) call error_handler (ncode=istat)
        enddo

!     ----  has this name been specified already? ----
       if (field_name_present(name,field_names,num_field_names)) cycle

!     ---- skip surface pressure and topography fields ----
        if ( trim(name) == 'res') cycle    ! for eta coordinate

!     ---- must be 2-d, 3-d, or 4-d field ----
        static = dimids(ndim) .ne. rdim
        threed = .false.
!!!del  if (static) cycle   ! will not allow static for now (see below)

!     ---- skip 3-d fields if flag set ----

        select case (ndim)
          case (4)
            if (static) cycle   ! not possible
            threed = .true. 
          case (3)
           !if (.not.static .and. do_all_3d_fields) cycle
            if (static) threed = .true.
            if (.not.static .and. .not.do_all_fields) cycle
           !if (static) cycle           ! no 3d static fields allowed
          case (2)
            if (.not.static) cycle
            if (.not.do_all_fields) cycle
          case default
            cycle
        end select

!     ---- axes must have correct size ----
        if (axlen(1) /= nlon .or.  axlen(2) /= nlat)   cycle
        if (threed) then
        if (axlen(3) /= nlev .and. axlen(3) /= nlev+1) cycle
        endif

!     ---- output this field ----
        num_field_names = num_field_names + 1
        field_names(num_field_names) = name
        if ( num_field_names == MAX_FIELD_NAMES ) return

     enddo

 end subroutine set_all_field_names

!#######################################################################

 function field_name_present ( name, field_names, num_field_names ) 
 character(len=*), intent(in) :: name, field_names(:)
 integer,          intent(in) :: num_field_names
 logical :: field_name_present
 integer :: i

     field_name_present = .false.

    !----  has this name been specified already? ----
     do i = 1, num_field_names
        if (trim(name) == trim(field_names(i))) then
           field_name_present = .true.
           return
        endif
     enddo

 end function field_name_present

!#######################################################################

 function skip_field_name ( field_name )

   character(len=NF90_MAX_NAME), intent(in) :: field_name
   logical :: skip_field_name
   integer :: i

!  ------ skip fields that are written separately ------

       !skip_field_name = .false.

       !if (trim(field_name) == trim(sphum_name)) then
       !    skip_field_name = .true.; return; endif
       !if (trim(field_name) ==  trim(temp_name)) then
       !    skip_field_name = .true.; return; endif
       !if (trim(field_name) == trim(zsurf_name)) then
       !    skip_field_name = .true.; return; endif

        if (trim(field_name) ==  trim(hght_name)) then
            skip_field_name = .true.; return; endif
        if (trim(field_name) ==   trim(slp_name)) then
            skip_field_name = .true.; return; endif

 end function skip_field_name

!#######################################################################

 subroutine get_input_variable_info ( name, ncid_in, Field )

   character(len=NF90_MAX_NAME), intent(in)    :: name
   integer          ,       intent(in)    :: ncid_in
   type(field_type) ,       intent(inout) :: Field
   logical :: answer
   integer :: rdim, zcoord, idim
   integer :: i, istat, ipack, ioffs, imiss, ifill
   real    :: mval, fval

 ! reset skip flag
 ! skip all fields unless properties are acceptable
   Field%skip = .true.
   Field%varid_in = 0

 ! skip certain field names
   if (trim(name) == trim(hght_name))  return
   if (trim(name) == trim(slp_name))   return

 ! get time axis id
   istat = NF90_INQUIRE (ncid_in, unlimiteddimid=rdim)
   if (istat /= NF90_NOERR) call error_handler (ncode=istat)
 ! get variable id 
   istat = NF90_INQ_VARID (ncid_in, trim(name), Field%varid_in)

 ! proceed if this variable exists
   if (istat /= NF90_NOERR) return

       if (verbose > 2) print *, 'getting input variable info, name=',trim(name)

       ! get number and names of dimensions
       istat = NF90_INQUIRE_VARIABLE (ncid_in, Field%varid_in, name=Field%name, &
                        xtype=Field%xtype, ndims=Field%ndim, dimids=Field%dimids)
       if (istat /= NF90_NOERR) call error_handler (ncode=istat)
       ! get dimension lengths
       do idim = 1, Field%ndim
          istat = NF90_INQUIRE_DIMENSION (ncid_in, Field%dimids(idim), len=Field%dimlen(idim))
          if (istat /= NF90_NOERR) call error_handler (ncode=istat)
       enddo
       Field%static = Field%dimids(Field%ndim) .ne. rdim
       Field%z_axis = 0
       Field%y_axis = 0

     ! check time averaging
       istat = NF90_INQUIRE_ATTRIBUTE (ncid_in, Field%varid_in, 'time_avg_info')
       Field%do_avg = istat == NF90_NOERR

     ! get packing attributes of input data
       ipack = NF90_GET_ATT (ncid_in, Field%varid_in, 'scale_factor', Field%pack(1))
       ioffs = NF90_GET_ATT (ncid_in, Field%varid_in, 'add_offset',   Field%pack(2))
       if (ipack == NF90_NOERR .and. ioffs == NF90_NOERR) then
           Field%do_pack = .true.
       else
           Field%do_pack = .false.
           Field%pack = (/1.,0./)
       endif

     ! get missing value attribute or fill value attribute
     ! NOTE: conversion of missing value data type to default real type
       imiss = NF90_GET_ATT (ncid_in, Field%varid_in, 'missing_value', mval)
       ifill = NF90_GET_ATT (ncid_in, Field%varid_in, '_FillValue',    fval)
       if (imiss == NF90_NOERR) then
          Field%do_miss_in = .true.
          Field%missval_in = mval
       else if (ifill == NF90_NOERR) then
          Field%do_miss_in = .true.
          Field%missval_in = fval
       else
          Field%do_miss_in = .false.
          Field%missval_in = 0.
       endif

     ! unpack input missing value (for comparison with unpacked data)
       Field%unpkd_missval_in = Field%missval_in
       if (Field%do_miss_in .and. Field%do_pack) then
           Field%unpkd_missval_in = Field%unpkd_missval_in*Field%pack(1) + Field%pack(2)
       endif

       Field%skip = .false.

 end subroutine get_input_variable_info

!#######################################################################

 subroutine copy_variable_information ( name, ncid_in, ncid_out, Field )

   character(len=NF90_MAX_NAME), intent(in)    :: name
   integer          ,       intent(in)    :: ncid_in, ncid_out
   type(field_type) ,       intent(inout) :: Field
   logical :: answer
   integer :: rdim, zcoord, idim
   integer :: i, istat, ipack, ioffs, imiss, ifill
   integer :: dimids(4)
   real    :: mval, fval

   call get_input_variable_info ( name, ncid_in, Field )

   if (Field%skip) return
   Field%skip = .true.
   zcoord = 0

   ! set output missing value

     Field%do_miss_out = Field%do_miss_in
     Field%missval_out = Field%missval_in

     if (use_default_missing_value) then
         if (Field%xtype == NF90_REAL8) Field%missval_out = NF90_FILL_REAL8
         if (Field%xtype == NF90_REAL4) Field%missval_out = NF90_FILL_REAL4
         if (Field%xtype == NF90_INT)   Field%missval_out = NF90_FILL_INT
         if (Field%xtype == NF90_INT2)  Field%missval_out = NF90_FILL_INT2
     endif

     ! unpack output missing value (for comparison with unpacked data)
       Field%unpkd_missval_out = Field%missval_out
       if (Field%do_miss_out .and. Field%do_pack) then
           Field%unpkd_missval_out = Field%unpkd_missval_out*Field%pack(1) + Field%pack(2)
       endif

   ! do not process fields with one dimension
     if (Field%ndim == 1) return

   !---- check horizontal dimensions (again) ---
     if (Field%ndim >= 2) then
         if (Field%dimlen(1) /= nlon .or. Field%dimlen(2) /= nlat) then
               if (Field%ndim == 2 .and. Field%static) return
               if (Field%ndim >= 3)                    return
         endif
     endif

   ! set dimension identifiers for output file

     if (Field%dimlen(2) == nlat)   dimids = dimids_mass

   !---- copy input variable metadata to output file ----

     select case (Field%ndim)
       case (4)
            if (Field%dimlen(3) == nlev)   zcoord = 1
            if (Field%dimlen(3) == nlev+1) zcoord = 2
            if (zcoord == 0) return
            if (Field%static) return   ! 4d static fields not allowed
            call copy_variable ( ncid_in, trim(name), ncid_out, &
                                 dimids, Field%varid_out,       &
                                 missvalue=Field%missval_out  )
       case (3)
          if (Field%static) then
            if (Field%dimlen(3) == nlev)   zcoord = 1
            if (Field%dimlen(3) == nlev+1) zcoord = 2
            if (zcoord == 0) return
            return               ! no 3d static fields allowed
            call copy_variable ( ncid_in, trim(name), ncid_out, &
                                 dimids(1:3), Field%varid_out,  &
                                 missvalue=Field%missval_out )
          else
             !---- 2d non-static field ----
             if (Field%do_miss_in) then
                call copy_variable ( ncid_in, trim(name), ncid_out,      &
                      (/dimids(1),dimids(2),dimids(4)/), Field%varid_out,&
                      missvalue=Field%missval_out )
             else
                call copy_variable ( ncid_in, trim(name), ncid_out,      &
                      (/dimids(1),dimids(2),dimids(4)/), Field%varid_out )
             endif
          endif
       case(2)
          !---- 2d static field ----
          if (Field%static) then
             if (Field%do_miss_in) then
                call copy_variable ( ncid_in, trim(name), ncid_out, &
                                     dimids(1:2), Field%varid_out,  &
                                     missvalue=Field%missval_out    )
             else
                call copy_variable ( ncid_in, trim(name), ncid_out, &
                                     dimids(1:2), Field%varid_out   )
             endif
          endif
       case default
            return
     end select

     Field%skip = .false.
     Field%y_axis = Field%dimlen(2)
     Field%z_axis = zcoord

 end subroutine copy_variable_information

!#######################################################################

 subroutine copy_axis_information (ncid_in, name_in, ncid_out, decomp, &
                                   Axis, dimid_out, name_out, name_edges)

 integer,          intent(in)    :: ncid_in, ncid_out, decomp(4)
 character(len=*), intent(in)    :: name_in
 type(axis_type),  intent(inout) :: Axis
 integer,          intent(out), optional :: dimid_out
 character(len=*), intent(out), optional :: name_out, name_edges

 integer :: varid_in, varid_out, dimlen_out, dimlen_in, did_out

   if (verbose > 2) print *, 'COPYING, name_in=',trim(name_in)

   if (decomp(1) == decomp(3) .and. decomp(2) == decomp(4)) then
      ! global domain indices == compute domain indices
      call copy_axis ( ncid_in, name_in, ncid_out, did_out, varid_in, varid_out, &
                       dimlen_out, name_out=name_out, name_edges=name_edges)
   else
      dimlen_in = decomp(4)-decomp(3)+1
      call copy_axis ( ncid_in, name_in, ncid_out, did_out, varid_in, varid_out, dimlen_out, &
                       dimlen_in=dimlen_in, name_out=name_out, name_edges=name_edges)
   endif

   Axis%name_in   = name_in
   Axis%name_out  = name_in
   if (present(name_out)) Axis%name_out  = name_out
   Axis%varid_in  = varid_in 
   Axis%varid_out = varid_out
   Axis%dimlen    = dimlen_out
   Axis%start     = (/decomp(3)/)

   if (present(dimid_out)) dimid_out = did_out

 end subroutine copy_axis_information

!#######################################################################

 subroutine set_up_input_mask (ncid, Field, data, mask, kbot)

   integer,          intent(in)  :: ncid
   type(field_type), intent(in)  :: Field
   real,             intent(in)  :: data(:,:,:)
   real,             intent(out) :: mask(:,:,:)
   integer,          intent(out) :: kbot(:,:)

   integer :: k, nlev

!     ----- set up mask and bottom level index -----

      nlev = size(data,3)

      kbot = nlev
      mask = 1.0

      if (Field%do_miss_in) then
         where ( data == Field%unpkd_missval_in )
                 mask = 0.0
         elsewhere
                 mask = 1.0
         endwhere

         do k = nlev, 1, -1
            if (count(mask(:,:,k) < 0.5) == 0) exit
            where (mask(:,:,k) < 0.5) kbot(:,:) = k-1
         enddo
      endif

 end subroutine set_up_input_mask

!#######################################################################

 subroutine read_and_set_up_mask (ncid, start, mask, kbot)

   integer,           intent(in)  :: ncid, start(4)
   real,              intent(out) :: mask(:,:,:)
   integer,           intent(out) :: kbot(:,:)

   character(len=NF90_MAX_NAME) :: name
   real,  dimension(size(mask,1),size(mask,2),size(mask,3)) :: values
   integer :: ivar, nvar, istat
   logical :: done
   type(field_type) :: Field

     done = .false.

   ! get number of variables and time axis id
     istat = NF90_INQUIRE (ncid, nvariables=nvar)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

     do ivar = 1, nvar
        istat = NF90_INQUIRE_VARIABLE ( ncid, ivar, name=name )
        if (istat /= NF90_NOERR) call error_handler (ncode=istat)
        call get_input_variable_info ( name, ncid, Field )

      ! find first 3d mass field
        if (Field%ndim <= 2) cycle
        if (Field%ndim <= 3 .and. .not.Field%static) cycle
        if (Field%ndim == 4 .and.      Field%static) cycle

      ! must have mass dimension lengths and full vertical levels
        if (Field%dimlen(1) /= nlon .or. Field%dimlen(2) /= nlat .or. Field%dimlen(3) /= nlev) cycle

      ! must have input missing value
        if (.not.Field%do_miss_in) cycle

        call read_variable     (ncid, Field, values, start )
        call set_up_input_mask (ncid, Field, values, mask, kbot)
        done = .true.
        exit

     enddo

     if (.not.done) then
        mask = 1.0
        kbot = nlev
     endif

 end subroutine read_and_set_up_mask

!#######################################################################

 subroutine read_write_static (ncid_in, ncid_out, is, ie, js, je, Field)
   integer         , intent(in) :: ncid_in, ncid_out, is, ie, js, je
   type(field_type), intent(in) :: Field

   real :: values(is:ie,js:je,1)

 ! read and write static fields that are on the mass grid (nlon x nlat)

 ! read data
   call read_variable ( ncid_in, Field, values, start=(/is, js/))

 ! write data
   call write_variable ( ncid_out, Field, values)

 end subroutine read_write_static

!#######################################################################

 subroutine read_variable ( ncid, Field, values, start )
 integer,           intent(in)  :: ncid
 type(field_type),  intent(in)  :: Field
 real,              intent(out) :: values(:,:,:)
 integer, optional, intent(in)  :: start(:)

 integer :: ndim, nd, begin(4)

 ! number of dimensions (less the time dimension)
   ndim = Field%ndim
   if (.not.Field%static) ndim = ndim-1

 ! start indices
   nd    = Field%ndim
   begin = 1
   if (present(start)) then
       nd = size(start)
       begin(1:nd) = start
   endif

   select case (ndim)
      case(1)
          istat = NF90_GET_VAR (ncid, Field%varid_in, values(:,1,1), start=begin(1:nd))
          if (istat /= NF90_NOERR) call error_handler ('getting 1d var='//trim(Field%name), ncode=istat)
      case(2)
          istat = NF90_GET_VAR (ncid, Field%varid_in, values(:,:,1), start=begin(1:nd))
          if (istat /= NF90_NOERR) call error_handler ('getting 2d var='//trim(Field%name), ncode=istat)
      case(3)
          istat = NF90_GET_VAR (ncid, Field%varid_in, values(:,:,:), start=begin(1:nd))
          if (istat /= NF90_NOERR) call error_handler ('getting 3d var='//trim(Field%name), ncode=istat)
      case default
          call error_handler ('invalid number of dimensions in read_variable, var='//trim(Field%name))
   end select

 ! unpack values
   if (Field%do_pack) then
       if (Field%do_miss_in) then
           where (values .ne. Field%missval_in)
                 values = values*Field%pack(1) + Field%pack(2)
           elsewhere
                 values = Field%unpkd_missval_in
           endwhere
       else
           values = values*Field%pack(1) + Field%pack(2)
       endif
   endif

 end subroutine read_variable

!#######################################################################

 subroutine write_variable ( ncid, Field, values, start )
 integer,           intent(in)    :: ncid
 type(field_type),  intent(in)    :: Field
 real,              intent(inout) :: values(:,:,:)
 integer, optional, intent(in)    :: start(:)

 integer :: ndim, nd, begin(4)

 ! pack values (overwrite input data)
   if (Field%do_pack) then
       if (Field%do_miss_out) then
           where (values .ne. Field%unpkd_missval_out)
               values = nint((values-Field%pack(2))/Field%pack(1))
           elsewhere
               values = Field%missval_out
           endwhere
       else
           values = nint((values-Field%pack(2))/Field%pack(1))
       endif
   endif

 ! number of dimensions (less the time dimension)
   ndim = Field%ndim
   if (.not.Field%static) ndim = ndim-1

 ! start indices
   nd    = Field%ndim
   begin = 1
   if (present(start)) then
       nd = size(start)
       begin(1:nd) = start
   endif

   select case (ndim)
      case(1)
          istat = NF90_PUT_VAR (ncid, Field%varid_out, values(:,1,1), start=begin(1:nd))
          if (istat /= NF90_NOERR) call error_handler ('putting 1d var='//trim(Field%name), ncode=istat)
      case(2)
          istat = NF90_PUT_VAR (ncid, Field%varid_out, values(:,:,1), start=begin(1:nd))
          if (istat /= NF90_NOERR) call error_handler ('putting 2d var='//trim(Field%name), ncode=istat)
      case(3)
          istat = NF90_PUT_VAR (ncid, Field%varid_out, values(:,:,:), start=begin(1:nd))
          if (istat /= NF90_NOERR) call error_handler ('putting 3d var='//trim(Field%name), ncode=istat)
      case default
          call error_handler ('invalid number of dimensions in write_variable, var='//trim(Field%name))
   end select

 end subroutine write_variable

!#######################################################################

 subroutine compute_pres_full (pk, bk, ph, pf)
 real, intent(in)  :: pk(:), bk(:), ph(:,:,:)
 real, intent(out) :: pf(:,:,:)
 real, dimension(size(pf,1),size(pf,2))            :: dp, lpf
 real, dimension(size(pf,1),size(pf,2),size(ph,3)) :: lph
 integer :: k, nlev

 nlev = size(pf,3)
 
 ! compute p*logp at half levels
   where (ph(:,:,1) > 0.)
      lph(:,:,1) = ph(:,:,1) * log(ph(:,:,1))
   elsewhere
      lph(:,:,1) = 0.0
   endwhere
   lph(:,:,2:nlev+1) = ph(:,:,2:nlev+1) * log(ph(:,:,2:nlev+1))

! compute pressure at full levels
  do k = 1, nlev
      dp(:,:) = ph(:,:,k+1)-ph(:,:,k)
      lpf(:,:) = (lph(:,:,k+1)-lph(:,:,k))/dp(:,:) - 1.0
      pf(:,:,k) = exp((lph(:,:,k+1)-lph(:,:,k))/dp(:,:) - 1.0)
  enddo

 end subroutine compute_pres_full

!#######################################################################

 subroutine compute_height ( zsurf, temp, sphum, pfull, phalf, zfull, mask)
 real, intent(in),  dimension(:,:)   :: zsurf
 real, intent(in),  dimension(:,:,:) :: temp, sphum, pfull, phalf
 real, intent(out), dimension(:,:,:) :: zfull
 real, intent(in),  optional         :: mask(:,:,:)
 real, dimension(size(temp,1),size(temp,2)) ::  wta, wtb, vt, zb, zt, &
                                                lpb, lpt, lpf
 integer :: k
 logical :: flag

! flag for zero ptop
flag = count(phalf(:,:,1)<=0.) > 0

! vertical integration (bottom to top)
 zb(:,:) = zsurf(:,:)*GRAV
 lpb(:,:) = log(phalf(:,:,size(zfull,3)+1))
 do k = size(zfull,3), 1, -1
    lpf = log(pfull(:,:,k))
    wtb = lpb - lpf
    ! check for ptop=0
    if (k == 1 .and. flag) then
       wta = wtb
    else
       lpt = log(phalf(:,:,k))
       wta = lpf - lpt
    endif

    vt = temp(:,:,k)*(1.0+d608*sphum(:,:,k))*RDGAS ! virtual temperature * Rgas
    zt = zb + vt*(wta+wtb)
    zfull(:,:,k) = (zb + vt*wtb)/GRAV
    if (present(mask)) then
       where (mask(:,:,k) < 0.01)
          zt = GRAV*tref*(1.0-(pfull(:,:,k)/pref)**rgog)/tlapse
          zfull(:,:,k) = 0.5*(zt+zb)/GRAV
       endwhere
    endif
    zb  = zt
    lpb = lpt
 enddo

 end subroutine compute_height

!#######################################################################

end program run_pressure_interp

