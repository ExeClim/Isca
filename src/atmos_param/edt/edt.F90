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

module edt_mod

!=======================================================================
!
!
!
!      EDT (Entrainment and Diagnostic Turbulence) MODULE
!
!
!      February 2002
!      Contact person: Steve Klein
!
!
!      These routines calculate the diffusivity coefficients for
!      momentum and temperature-moisture-scalars using the moist
!      thermodynamcs modules based on:
!
!      H. Grenier and C. Bretherton, 2001: A moist PBL parameterization
!      for large-scale models and its application to subtropical
!      cloud-topped marine boundary layers. Mon. Wea. Rev., 129,
!      357-377.
!
!      The actual routine is not described in this paper but is
!      a simplified extension of the parameterization discussed
!      here.  The original code, given to Steve Klein from 
!      Chris Bretherton in May 2001, was tested in the NCAR 
!      atmospheric model, formerly known as CCM. The code has 
!      been adapted for the FMS system by Steve Klein and Paul
!      Kushner.
!
!
!      To quote the Bretherton and Grenier description:
!
!      Driver routine to compute eddy diffusion coefficients for 
!      momentum, moisture, trace constituents and static energy.  Uses 
!      first order closure for stable turbulent layers. For convective 
!      layers, an entrainment closure is used, coupled to a diagnosis 
!      of layer-average TKE from the instantaneous thermodynamic and 
!      velocity profiles. Convective layers are diagnosed by extending 
!      layers of moist static instability into adjacent weakly stably 
!      stratified interfaces, stopping if the stability is too strong.  
!      This allows a realistic depiction of dry convective boundary 
!      layers with a downgradient approach."
! 
!      Authors:  Herve Grenier, 06/2000, Chris Bretherton 09/2000
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
! outside modules used
!

use      constants_mod, only: grav,vonkarm,cp_air,rdgas,rvgas,hlv,hls, &
                              tfreeze,radian

use            mpp_mod, only: input_nml_file
use            fms_mod, only: file_exist, open_namelist_file, error_mesg, FATAL,&
                              NOTE, mpp_pe, mpp_root_pe, close_file, read_data, &
                              write_data, write_version_number, stdlog, &
                              open_file, check_nml_error
use fms_io_mod,         only: register_restart_field, restart_file_type, &
                              save_restart, restore_state

use   diag_manager_mod, only: register_diag_field, send_data
        
use   time_manager_mod, only: time_type, get_date, month_name
 
use  monin_obukhov_mod, only: mo_diff

use sat_vapor_pres_mod, only: compute_qs

implicit none
private

!-----------------------------------------------------------------------
!
!      public interfaces

public edt, edt_init, edt_end, edt_on, qaturb, qcturb,tblyrtau

!-----------------------------------------------------------------------
!
!      global storage variable
!

real, allocatable, dimension(:,:,:) :: qaturb ! cloud fraction diagnosed
                                              ! from turbulence model
      ! (fraction)
real, allocatable, dimension(:,:,:) :: qcturb ! cloud condensate 
                                              ! diagnosed from turb.
      ! model (kg liq/kg air)
real, allocatable, dimension(:,:,:) :: tblyrtau  ! turbulent layer
                                                 ! time scale
 ! (seconds)
real, allocatable, dimension(:,:,:) :: sigmas ! standard deviation of 
                                              ! water perturbation
                                              ! (kg water/kg air)     

type(restart_file_type), pointer, save :: edt_restart => NULL()
!-----------------------------------------------------------------------
!
!      set default values to namelist parameters       
!

character :: sftype  = "z"           ! method for calculating sat frac
real      :: qcminfrac = 1.e-3       ! Min condensate counted as cloud
                                     ! as a fraction of qsat.
logical   :: use_qcmin = .true.      ! Use qcminfrac as indicator of
                                     ! of cloud top.
logical   :: use_extrapolated_ql  = .false.  ! should the layer top
                                             ! liquid water be used to
     ! estimate the evaporative
     ! enhancement
integer   :: n_print_levels = 14     ! how many of the lowest levels 
                                     ! should be printed out
integer, dimension(2) :: edt_pts = 0 ! the global indices for i,j
                                     ! at which diagnostics will 
                                     ! print out
logical   :: do_print = .false.      ! should selected variables 
                                     ! be sent to logfile
logical   :: column_match = .false.  ! should this column be printed 
                                     ! out?
integer   :: dpu = 0                 ! unit # for do_print output
 
real      :: min_adj_time = 0.25     ! minimum adjustment time of
                                     ! turbulent layer for thermo-
                                     ! dynamics, as a fraction of
                                     ! the physics time step. 
! 
! the following quantities only deal with the gaussian cloud model 
!
 
logical   :: do_gaussian_cloud = .false.  

real      :: kappa   = 0.5       ! absolute value of the correlation
                                 ! coefficient between sli and qt 
 ! fluctions
real      :: mesovar = 0.02      ! amplitude to mesoscale fluctuations
                                 ! in sigma-s as a fraction of the 
 ! saturation specific humidity
         

integer, parameter                 :: MAX_PTS = 20
integer, dimension (MAX_PTS)       :: i_edtprt_gl=0, j_edtprt_gl=0
real, dimension(MAX_PTS)           :: lat_edtprt=999., lon_edtprt=999.
integer                            :: num_pts_ij = 0
integer                            :: num_pts_latlon = 0



namelist /edt_nml/sftype, qcminfrac, use_qcmin,kappa, mesovar,edt_pts,    &
                  do_gaussian_cloud, n_print_levels, use_extrapolated_ql, &
                  min_adj_time, i_edtprt_gl, j_edtprt_gl,                 &
                  num_pts_ij, num_pts_latlon, lat_edtprt, lon_edtprt
  
 
integer     :: num_pts           !  total number of columns in which
                                 !  diagnostics are desired

!---------------------------------------------------------------------
!    deglon1 and deglat1 are the longitude and latitude of the columns
!    at which diagnostics will be calculated (degrees).
!---------------------------------------------------------------------
real,    dimension(:), allocatable  :: deglon1, deglat1
 
!---------------------------------------------------------------------
!    iradprt and jradprt are the processor-based i and j coordinates 
!    of the desired diagnostics columns.
!---------------------------------------------------------------------
integer, dimension(:), allocatable  :: j_edtprt, i_edtprt

!---------------------------------------------------------------------
!    do_raddg is an array of logicals indicating which latitude rows
!    belonging to the processor contain diagnostics columns.
!---------------------------------------------------------------------
logical, dimension(:), allocatable  :: do_edt_dg

!-----------------------------------------------------------------------
!
!      diagnostic fields       
!

character(len=10) :: mod_name = 'edt'
real              :: missing_value = -999.
integer           :: id_fq_cv_int,    id_fq_cv_top,      id_fq_cv_bot, &
                     id_fq_st,   id_n2,      id_s2,      id_ri,        &
                     id_leng,    id_bprod,   id_sprod,   id_trans,     &
                     id_diss,    id_qaedt,   id_qcedt,   id_tauinv,    &
                     id_eddytau, id_fq_turb, id_sigmas,  id_radf,      &
                     id_sh,      id_sm,      id_gh,      id_evhc
       
      
!-----------------------------------------------------------------------
!
!      set default values to parameters       
!

logical         :: edt_on = .false.
logical         :: init = .false.
real, parameter :: small  = 1.e-8      
real, parameter :: fpi  = 3.14159     ! pi
real :: d608 = 0.
real :: d622 = 0.
real :: d378 = 0.
real, parameter :: frac_sfclyr = 0.1  ! height of surface layer top as a 
                                      ! fraction of the pbl height

!-----------------------------------------------------------------------
!
!      the following parameters are those defined only in the routines
!      provided by Chris Bretherton and Herve Grenier
!

real, parameter :: ntzero  = 1.e-10   ! not zero, used to set min value
                                      ! to s2, shear-squared.
real :: zvir    = 0.

real, parameter :: b1      =   5.8    ! TKE dissipation = e^3/(b1*leng)
real            :: b123    =   3.2281 ! b1**(2/3)
real, parameter :: tunl    =   0.085  ! Asympt leng = 
                                      !       tunl*(turb lyr depth)
real, parameter :: alph1   =   0.5562 ! Galperin stability fn params
real, parameter :: alph2   =  -4.3640
real, parameter :: alph3   = -34.6764
real, parameter :: alph4   =  -6.1272
real, parameter :: alph5   =   0.6986
real, parameter :: ricrit  =   0.19   ! Critical Richardson # for turb.
real, parameter :: mu      =  70.     ! used in finding e/ebar in 
                                      ! convective layers (CLs)
real, parameter :: rinc    =  -0.5    ! Min W/<W> for incorp into CL

! params governing entr. effic. A=a1l*evhc, evhc=1+a2l*a3l*L*ql/jt2slv
! where ql is cloud-top liq. water and jt2slv is the jump in slv across
! the cloud-top entrainment zone.

real, parameter :: a1l       =  0.10  ! Dry entrainment efficiency param
                                      ! Herve set to 0.05 due to excess 
      ! TKE but a1l = 0.10 = 0.2*tunl*
      ! erat^-1.5 should be the "real" 
      ! value, where erat = <e>/wstar^2 
      ! for dry CBL = 0.3
real, parameter :: a2l       = 15.    ! Moist entrainment enhancement 
                                      ! param Herve's SCCM value was 15
real, parameter :: a3l       =  0.8   ! 
real, parameter :: jbumin    =  0.001 ! Min buoyancy jump at an entrain-
                                      ! ment interface (m/s2)
      ! (~ jump in K/30)
real, parameter :: evhcmax   = 10.    ! Max entrainment efficiency
real, parameter :: rimaxentr =  0.    ! Limiting Ri for entraining turb
                                      ! layer
      
!  parameters affecting TKE

real, parameter :: rmin    =   0.1    ! Min allowable e/<e> in a CL
real, parameter :: rmax    =   2.0    ! Max allowable e/<e> in a CL
real, parameter :: tkemax  =  20.     ! tke capped at tkemax (m2/s2)
real, parameter :: tkemin  =   1.e-6  ! tke minimum (m2/s2)

!-----------------------------------------------------------------------
!
! declare version number 
!

character(len=128) :: Version = '$Id: edt.F90,v 19.0 2012/01/06 20:09:18 fms Exp $'
character(len=128) :: Tagname = '$Name: siena_201211 $'
logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------
!
! fms module subroutines include:
!
!      edt         main driver program of the module
!
!      edt_init    initialization routine       
!
!      edt_tend    adds in the longwave heating rate to the 
!                         global storage variable
!
!      edt_end     ending routine
!
!
!      Grenier-Bretherton subroutines are described after the 
!      subroutines listed above.
!


      integer, dimension(1) :: restart_versions = (/ 1 /)

contains




!======================================================================= 
!
!      subroutine edt_init 
!        
!
!      this subroutine reads the namelist file and restart data
!      and initializes some constants.
!        

subroutine edt_init(lonb, latb, axes,time,idim,jdim,kdim)

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
! 
!      idim,jdim,kdim    size of the first 3 dimensions 
!      axes, time        variables needed for netcdf diagnostics
!      latb, lonb        latitudes and longitudes at grid box corners
!
!
!      --------
!      internal
!      --------
! 
!      unit              unit number for namelist and restart file
!      io                internal variable for reading of namelist file
!      full              indices for full level axes coordinates
!      half              indices for half level axes coordinates
!
!-----------------------------------------------------------------------

integer,              intent(in) :: idim,jdim,kdim,axes(4)
type(time_type),      intent(in) :: time
real, dimension(:,:), intent(in) :: lonb, latb

integer               :: unit, io, logunit, ierr, vers, vers2, nn, i, j, id_restart
real                  :: dellat, dellon
character(len=4)      :: chvers
integer, dimension(3) :: full = (/1,2,3/), half = (/1,2,4/)

!-----------------------------------------------------------------------
!
!      namelist functions

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=edt_nml, iostat=io)
   ierr = check_nml_error(io,'edt_nml')
#else   
   if (file_exist('input.nml')) then
      unit = open_namelist_file ()
      ierr=1 ; Do While (ierr .ne. 0)
        Read  (unit, nml=edt_nml, iostat=io, End=10)
        ierr = check_nml_error(io,'edt_nml')
      enddo
10    Call Close_File (unit)
   endif
#endif

!------- write version number and namelist ---------

   if ( mpp_pe() == mpp_root_pe() ) then
       call write_version_number(version, tagname)
       logunit = stdlog()
       write (logunit, nml=edt_nml)
   endif

!s initialised here as rdgas now an namelist parameter
d608 = (rvgas-rdgas)/rdgas
d622 = rdgas/rvgas
d378 = 1. - d622
zvir    = d608      

!---------------------------------------------------------------------
!    allocate and initialize a flag array which indicates the latitudes
!    containing columns where radiation diagnostics are desired.
!---------------------------------------------------------------------
   allocate (do_edt_dg (size(latb,2)-1) )
   do_edt_dg(:) = .false.

!-------------------------------------------------------------------
!    define the total number of points at which diagnostics are desired.
!    points may be specified either by lat-lon pairs or by global index
!    pairs. 
!-------------------------------------------------------------------
   num_pts = num_pts_latlon + num_pts_ij

!-------------------------------------------------------------------
!    continue on only if diagnostics are desired in at least one column.
!-------------------------------------------------------------------
   if (num_pts > 0) then

!-------------------------------------------------------------------
!    if more points are desired than space has been reserved for, print 
!    a message.
!-------------------------------------------------------------------
     if (num_pts > MAX_PTS) then
       call error_mesg ( 'edt_mod', &
      'must reset MAX_PTS or reduce number of diagnostics points', &
                                                          FATAL)
     endif

!-------------------------------------------------------------------
!    allocate space for arrays which will contain the lat and lon and
!    processor-local i and j indices.
!-------------------------------------------------------------------
     allocate ( deglon1 (num_pts))
     allocate ( deglat1 (num_pts))
     allocate ( j_edtprt (num_pts))
     allocate ( i_edtprt (num_pts))

!---------------------------------------------------------------------
!    if any points for diagnostics are specified by (i,j) global 
!    indices, determine their lat-lon coordinates. assumption is made 
!    that the deltas of latitude and longitude are uniform over 
!    the globe.
!---------------------------------------------------------------------
     do nn=1,num_pts_ij
       dellat = latb(1,2) - latb(1,1)
       dellon = lonb(2,1) - lonb(1,1)
       lat_edtprt(nn + num_pts_latlon) =     &
                   (-0.5*acos(-1.0) + (j_edtprt_gl(nn) - 0.5)*  &
                                        dellat) * radian
       lon_edtprt(nn + num_pts_latlon) =                & 
                    (i_edtprt_gl(nn) - 0.5)*dellon*radian
     end do

!--------------------------------------------------------------------
!    determine if the lat/lon values are within the global grid,
!    latitude between -90 and 90 degrees and longitude between 0 and
!    360 degrees.
!--------------------------------------------------------------------
     do nn=1,num_pts
       j_edtprt(nn) = 0
       i_edtprt(nn) = 0
       deglat1(nn) = 0.0
       deglon1(nn) = 0.0
       if (lat_edtprt(nn) .ge. -90. .and. &
           lat_edtprt(nn) .le.  90.) then
       else
         call error_mesg ('edt_mod', &
             ' invalid latitude for edt diagnostics ', FATAL)
       endif

       if (lon_edtprt(nn) .ge. 0. .and. &
           lon_edtprt(nn) .le. 360.) then
       else
         call error_mesg ('edt_mod', &
             ' invalid longitude for edt diagnostics ', FATAL)
       endif

!--------------------------------------------------------------------
!    determine if the diagnostics column is within the current 
!    processor's domain. if so, set a logical flag indicating the
!    presence of a diagnostic column on the particular row, define the 
!    i and j processor-coordinates and the latitude and longitude of 
!    the diagnostics column.
!--------------------------------------------------------------------
       do j=1,size(latb,2) - 1
         if (lat_edtprt(nn) .ge. latb(1,j)*radian .and.  &
             lat_edtprt(nn) .lt. latb(1,j+1)*radian) then
           do i=1,size(lonb,1) - 1
             if (lon_edtprt(nn) .ge. lonb(i,1)*radian   &
                               .and.&
                 lon_edtprt(nn) .lt. lonb(i+1,1)*radian)  &
                                then
               do_edt_dg(j) = .true.
               j_edtprt(nn) = j
               i_edtprt(nn) = i
               deglon1(nn) = 0.5*(lonb(i,1) + lonb(i+1,1))*radian
               deglat1(nn) = 0.5*(latb(1,j) + latb(1,j+1))*radian
               exit
             endif
           end do
           exit
         endif
       end do
     end do

!----------------------------------------------------------------------
!    open a unit for the radiation diagnostics output.
!---------------------------------------------------------------------
     dpu = open_file ('edt.out', action='write', &
                                 threading='multi', form='formatted')
     do_print = .true.

     if ( mpp_pe() == mpp_root_pe() ) then
        call write_version_number(version, tagname)
        write (dpu ,nml=edt_nml)
     endif
   endif     ! (num_pts > 0)

!-----------------------------------------------------------------------
!
!      initialize edt_on

   edt_on = .TRUE.
   module_is_initialized = .true.
!-----------------------------------------------------------------------
!
!      initialize b123 = b1**(2/3)

   b123 = b1**(2./3.)
       
!-----------------------------------------------------------------------
!
!      handle global storage
        
   if (allocated(qaturb)) deallocate (qaturb)
   allocate(qaturb(idim,jdim,kdim))
   if (allocated(qcturb)) deallocate (qcturb)
   allocate(qcturb(idim,jdim,kdim))
   if (allocated(tblyrtau)) deallocate (tblyrtau)
   allocate(tblyrtau(idim,jdim,kdim))
   if (allocated(sigmas)) deallocate (sigmas)
   allocate(sigmas(idim,jdim,kdim+1))


   allocate(edt_restart)
   id_restart = register_restart_field(edt_restart, 'edt.res.nc', 'qaturb'  , qaturb )
   id_restart = register_restart_field(edt_restart, 'edt.res.nc', 'qcturb'  , qcturb )
   id_restart = register_restart_field(edt_restart, 'edt.res.nc', 'tblyrtau', tblyrtau )
   id_restart = register_restart_field(edt_restart, 'edt.res.nc', 'sigmas'  , sigmas )

   if (File_Exist('INPUT/edt.res.nc')) then

     if (mpp_pe() == mpp_root_pe() ) then
       call error_mesg ('edt_mod',  'Reading netCDF formatted restart file: INPUT/edt.res.nc', &
                         NOTE)
     endif
     call restore_state(edt_restart)

   elseif (File_Exist('INPUT/edt.res')) then
     call error_mesg ('edt_mod', 'Native format restart file read no longer supported.',&
                       FATAL)
!     unit = Open_restart_File (FILE='INPUT/edt.res', ACTION='read')
!     read (unit, iostat=io, err=142) vers, vers2
!142  continue
! 
!!--------------------------------------------------------------------
!!    if eor is not encountered, then the file includes tdtlw as the
!!    first record (which this read statement read). that data is not 
!!    needed; note this and continue by reading next record.
!!--------------------------------------------------------------------
!     if (io == 0) then
!       call error_mesg ('edt_mod',  &
!           'reading pre-version number edt.res file, '//&
!           'ignoring tdtlw', NOTE)
!
!!--------------------------------------------------------------------
!!    if the first record was only one word long, then the file is a 
!!    newer one, and that record was the version number, read into vers. 
!!    if it is not a valid version, stop execution with a message.
!!--------------------------------------------------------------------
!     else
!       if (.not. any(vers == restart_versions) ) then
!         write (chvers, '(i4)') vers
!         call error_mesg ('edt_mod',  &
!              'restart version ' // chvers//' cannot be read '//&
!              'by this version of edt_mod.', FATAL)
!       endif
!     endif
!
!!---------------------------------------------------------------------
!     call read_data (unit, qaturb)
!     call read_data (unit, qcturb)
!     call read_data (unit, tblyrtau)
!     call read_data (unit, sigmas)
!     call Close_File (unit)
   else
     qaturb  (:,:,:) = 0.
     qcturb  (:,:,:) = 0.
     tblyrtau(:,:,:) = 0.
     sigmas  (:,:,:) = 0.
   endif

!-----------------------------------------------------------------------
!
! register diagnostic fields       

   id_fq_cv_int = register_diag_field (mod_name, 'fq_cv_int', axes(half), &
                  time, 'Frequency that the interface is in '//           &
                  'the interior of a convective layer from EDT',          &
                  'none', missing_value=missing_value )

   id_fq_cv_top = register_diag_field (mod_name, 'fq_cv_top', axes(half), &
                  time, 'Frequency that the interface is at '//           &
                  'the top of a convective layer from EDT',               &
                  'none', missing_value=missing_value )

   id_fq_cv_bot = register_diag_field (mod_name, 'fq_cv_bot', axes(half), &
                  time, 'Frequency that the interface is at '//           &
                  'the bottom of a convective layer from EDT',            &
                  'none', missing_value=missing_value )

   id_fq_st = register_diag_field (mod_name, 'fq_st', axes(half),  &
              time, 'Frequency of stable turbulence from EDT',     &
              'none', missing_value=missing_value )

   id_fq_turb = register_diag_field (mod_name, 'fq_turb', axes(full), &
                time, 'Frequency that the layer is fully '//          &
                'turbulent from EDT','none', missing_value=missing_value )

   id_n2 = register_diag_field (mod_name, 'n2', axes(half),        &
           time,'Moist Vaisala Frequency from EDT',                &
           '(1/sec)**2', missing_value=missing_value )
       
   id_s2 = register_diag_field (mod_name, 's2', axes(half),        &
           time,'Shear vector magnitude squared from EDT',         &
           '(1/sec)**2', missing_value=missing_value )
       
   id_ri = register_diag_field (mod_name, 'ri', axes(half),        &
           time, 'Moist Richardson number from EDT',               &
           'none', missing_value=missing_value )

   id_leng = register_diag_field (mod_name, 'leng', axes(half),    &
             time, 'Turbulent Length Scale from EDT',              &
             'meters', missing_value=missing_value )
       
   id_bprod = register_diag_field (mod_name, 'bprod', axes(half),  &
              time, 'Buoyancy production of TKE from EDT',         &
              '(meters**2)/(sec**3)', missing_value=missing_value )

   id_sprod = register_diag_field (mod_name, 'sprod', axes(half),  &
              time, 'Shear production of TKE from EDT',            &
              '(meters**2)/(sec**3)', missing_value=missing_value )
       
   id_trans = register_diag_field (mod_name, 'trans', axes(half),  &
              time, 'TKE transport from EDT',                      &
              '(meters**2)/(sec**3)', missing_value=missing_value )
       
   id_diss = register_diag_field (mod_name, 'diss', axes(half),    &
             time, 'TKE dissipation from EDT',                     &
             '(meters**2)/(sec**3)', missing_value=missing_value )

   id_radf  = register_diag_field (mod_name, 'radf', axes(half),   &
              time, 'TKE radiative forcing from EDT',              &
              '(meters**2)/(sec**3)', missing_value=missing_value )
       
   id_sh  = register_diag_field (mod_name, 'sh', axes(half),       &
            time, 'Galperin heat stability coefficient from EDT',  &
            'none', missing_value=missing_value )
 
   id_sm  = register_diag_field (mod_name, 'sm', axes(half),          &
            time, 'Galperin momentum stability coefficient from EDT', &
            'none', missing_value=missing_value )
 
   id_gh  = register_diag_field (mod_name, 'gh', axes(half),          &
            time, 'Galperin stability ratio from EDT',                &
            'none', missing_value=missing_value )
 
   id_evhc  = register_diag_field (mod_name, 'evhc', axes(half),       &
              time, 'Evaporative cooling enhancement factor from EDT', &
              'none', missing_value=missing_value )
 
   id_sigmas = register_diag_field (mod_name, 'sigmas', axes(half), &
               time, 'Std. dev. of water perturbation from EDT',    &
               '(kg water)/(kg air)', missing_value=missing_value )
       
   id_qaedt = register_diag_field (mod_name, 'qaedt', axes(full),  &
              time, 'statistical cloud fraction from EDT',         &
              'fraction', missing_value=missing_value )
       
   id_qcedt = register_diag_field (mod_name, 'qcedt', axes(full),  &
              time, 'statistical cloud condensate from EDT',       &
              'kg condensate/kg air', missing_value=missing_value )
       
   id_tauinv = register_diag_field (mod_name, 'tauinv', axes(half),&
               time, 'inverse large-eddy turnover time from EDT',  &
               '1/second', missing_value=missing_value )
       
   id_eddytau = register_diag_field (mod_name, 'eddytau', axes(full), &
                time, 'large-eddy turnover time from EDT',      &
                'seconds', missing_value=missing_value )
                          
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine edt_init

!
!======================================================================= 




!======================================================================= 
!
!      subroutine edt
!        
!
!       this subroutine is the main driver program to the routines
!       provided by Chris Bretherton
!        

subroutine edt(is,ie,js,je,dt,time,tdtlw_in, u_star,b_star,q_star,t,qv,ql,qi,qa, &
               u,v,z_full,p_full,z_half,p_half,stbltop,k_m,k_t,pblh,kbot,tke)

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
!
!      is,ie,js,je  i,j indices marking the slab of model working on
!      dt        physics time step (seconds)
!      time      variable needed for netcdf diagnostics
!      u_star    friction velocity (m/s)
!      b_star    buoyancy scale (m/(s**2))
!      q_star    moisture scale (kg vapor/kg air)
!
!      three dimensional fields on model full levels, reals dimensioned
!      (:,:,pressure), third index running from top of atmosphere to 
!      bottom
!          
!      t         temperature (K)
!      qv        water vapor specific humidity (kg vapor/kg air)
!      ql        liquid water specific humidity (kg cond/kg air)
!      qi        ice water specific humidity (kg cond/kg air)
!      qa        cloud fraction 
!      u         zonal wind (m/s)
!      v         meridional wind (m/s) 
!      z_full    height of full levels (m)
!      p_full    pressure (Pa)
!
!      the following two fields are on the model half levels, with
!      size(z_half,3) = size(t,3) +1, z_half(:,:,size(z_half,3)) 
!      must be height of surface (if you are not using eta-model)
!
!      z_half    height at half levels (m)
!      p_half    pressure at half levels (Pa)
!        
!      ------
!      output
!      ------
!
!      stbltop   maximum altitude the very stable boundary layer
!                is permitted to operate
!
!      The following variables are defined at half levels and are
!      dimensions 1:nlev+1.
!
!      k_m       diffusivity for momentum (m**2/s)
!      k_t       diffusivity for temperature and scalars (m**2/s)
!
!      k_m and k_t are defined at half-levels so that size(k_m,3) 
!      should be at least as large as size(t,3). Note, however, that 
!      only the returned values at levels 2 to size(t,3) are 
!      meaningful; other values will be returned as zero.
!
!      --------------
!      optional input
!      --------------
!
!      kbot      integer indicating the lowest true layer of atmosphere
!
!      ---------------
!      optional output
!      ---------------
!
!      pblh      depth of planetary boundary layer (m)
!      tke       turbulent kinetic energy (m*m)/(s*s)
!
!      --------
!      internal
!      --------
!
!      z_surf    height of surface (m)
!      z_full_ag height of full model levels above the surface (m)
!      qt        total water specific humidity (kg water/kg air)
!      qc        cloud condensate spec. hum. (kg condensate/kg air)
!      qsl       saturation specific humidity at the midpoint pressure
!                and liquid-ice water temperature (kg water/kg air)
!      dqsldtl   temperature derivative of qsl (kg H20/kg air/Kelvin)
!      hleff     effective latent heat of vaporization (J/kg)
!      sli       ice-liq water static energy (J/kg) 
!      sliv      ice-liq water virtual static energy (J/kg)                                      
!      khfs      surface kinematic heat flux (K*(m/s)) 
!      kqfs      surface kinematic vapor flux (kg/kg)*(m/s)
!      slislope  sli slope wrt pressure in thermo layer (J/kg/Pa)
!      qtslope   qt slope wrt pressure in thermo layer (kg/kg/Pa)
!      qxtop     saturation excess at the top of the layer 
!                (kg wat/kg air)
!      qxbot     saturation excess at the bottom of the layer 
!                (kg wat/kg air)
!      sfuh      saturated fraction in upper half-layer
!      sflh      sflh saturated fraction in lower half-layer
!      sfclyr_h  height of the surface layer top (m)
!      ql_new    new value of cloud liquid  (kg cond/kg air)
!      qi_new    new value of cloud ice     (kg cond/kg air)
!      qa_new    new value of cloud fraction
!      mask      real array indicating the point is above the surface  &
!                if equal to 1.0 and indicating the point is below the &
!                surface if equal to 0.
!
!      mineddytau  minimum value to adjustment time (1/sec)
!
!
!      the following variables are defined on model half levels
!      (1:kdim+1)
!
!      z_half_ag height of half model levels above the surface (m)
!      chu       heat var. coef for dry states (1/m)
!      chs       heat var. coef for sat states (1/m)
!      cmu       moisture var. coef for dry states (kg/kg)*(m/s*s)
!      cms       moisture var. coef for sat states (kg/kg)*(m/s*s)
!      n2        moist squared buoyancy freq (1/s*s)
!      s2        squared deformation, or shear vector mag. (1/s*s)
!      ri        gradient Richardson number
!
!      formulas for selected internal variables
!      ----------------------------------------
!
!      qt   = qv + ql + qi                 qc   = ql + qi
!      sli  = cp*T + g*z - hleff*qc        sliv = sli*(1.+d608*qt)
!      khfs = mean of (w'T')               kqfs = mean of (w'q')
!
!-----------------------------------------------------------------------

integer,         intent(in)                            :: is,ie,js,je
real,            intent(in)                            :: dt
type(time_type), intent(in)                            :: time
real,            intent(in),  dimension(:,:)           :: u_star, b_star, q_star
real,            intent(in),  dimension(:,:,:)         :: tdtlw_in, &
                                                          t,qv,ql,qi,qa, &
                                                          u, v, &
                                                          z_full, p_full, &
                                                          z_half, p_half
real,            intent(out), dimension(:,:)           :: stbltop, pblh
real,            intent(out), dimension(:,:,:)         :: k_m,k_t
integer,         intent(in),  dimension(:,:),  optional:: kbot
real,            intent(out), dimension(:,:,:),optional:: tke

integer                                         :: i,j,k,kk,ibot
integer                                         :: ipt,jpt
integer                                         :: nlev,nlat,nlon
integer, dimension(4,size(t,1),size(t,2),size(t,3)+1):: turbtype

real                                            :: khfs,kqfs
real                                            :: sfclyr_h_max
real                                            :: mineddytau
real, dimension(size(t,1),size(t,2))            :: z_surf, sfclyr_h
real, dimension(size(t,1),size(t,2),size(t,3))  :: z_full_ag
real, dimension(size(t,1),size(t,2),size(t,3))  :: qt,sli,sliv
real, dimension(size(t,1),size(t,2),size(t,3))  :: esl, qsl, dqsldtl
real, dimension(size(t,1),size(t,2),size(t,3))  :: hleff,density
real, dimension(size(t,1),size(t,2),size(t,3))  :: slislope, qtslope
real, dimension(size(t,1),size(t,2),size(t,3))  :: qxtop, qxbot
real, dimension(size(t,1),size(t,2),size(t,3))  :: qa_new, qc_new
real, dimension(size(t,1),size(t,2),size(t,3))  :: mask, isturb
real, dimension(size(t,1),size(t,2),size(t,3))  :: eddytau,tauinvtmp

real, dimension(size(t,1),size(t,2),size(t,3)+1):: z_half_ag
real, dimension(size(t,1),size(t,2),size(t,3)+1):: n2, s2, ri, tauinv
real, dimension(size(t,1),size(t,2),size(t,3)+1):: leng, bprod, sprod
real, dimension(size(t,1),size(t,2),size(t,3)+1):: trans, diss, evhc
real, dimension(size(t,1),size(t,2),size(t,3)+1):: radf, sh, sm, gh
real, dimension(size(t,1),size(t,2),size(t,3)+1):: mask3, tmpdat
real, dimension(size(t,1),size(t,2),size(t,3)+1):: k_t_mo, k_m_mo

integer, dimension(4,size(t,3)+1) :: gb_turbtype

real                        :: gb_pblh,     gb_u_star
real, dimension(size(t,3))  :: gb_t,        gb_u,        gb_v
real, dimension(size(t,3))  :: gb_qv
real, dimension(size(t,3))  :: gb_dqsldtl,  gb_qa_new,   gb_qc_new
real, dimension(size(t,3))  :: gb_sli,      gb_qt
real, dimension(size(t,3))  :: gb_qc,       gb_sliv,     gb_sflh
real, dimension(size(t,3))  :: gb_slisl,    gb_qtsl,     gb_sfuh
real, dimension(size(t,3))  :: gb_tdtlw,    gb_hleff
real, dimension(size(t,3))  :: gb_qsl,      gb_esl,      gb_qxtop
real, dimension(size(t,3))  :: gb_qxbot,    gb_density,  gb_z_full
real, dimension(size(t,3))  :: gb_p_full,   gb_isturb,   gb_qltop
real, dimension(size(t,3)+1):: gb_cmu,      gb_chu
real, dimension(size(t,3)+1):: gb_chs,      gb_cms
real, dimension(size(t,3)+1):: gb_k_m,      gb_k_t,      gb_tke
real, dimension(size(t,3)+1):: gb_n2,       gb_s2,       gb_ri
real, dimension(size(t,3)+1):: gb_leng,     gb_bprod,    gb_sprod
real, dimension(size(t,3)+1):: gb_diss,     gb_trans,    gb_tauinv
real, dimension(size(t,3)+1):: gb_z_half,   gb_sigmas,   gb_p_half
real, dimension(size(t,3)+1):: gb_radf,     gb_sh,       gb_sm
real, dimension(size(t,3)+1):: gb_gh,       gb_evhc
logical :: used, topfound
 
integer, dimension(MAX_PTS) :: nsave
integer :: iloc(MAX_PTS), jloc(MAX_PTS), nn, npts, nnsave
integer :: year, month, day, hour, minute, second
character(len=16) :: mon


                    
!-----------------------------------------------------------------------
!
!      initialize variables

   pblh   = 0.0       
   
   k_t    = 0.0
   k_m    = 0.0
   tke    = 0.0
   
   turbtype = 0
   isturb = 0.0
   
   n2     = 0.0
   s2     = 0.0
   ri     = 0.0
   tauinv = 0.0
   leng   = 0.0
   bprod  = 0.0
   sprod  = 0.0
   trans  = 0.0
   diss   = 0.0
   radf   = 0.0
   sh     = 0.0
   sm     = 0.0
   evhc   = 0.0
   gh     = 0.0
   
   slislope = 0.0
   qtslope  = 0.0
   
   qxbot  = 0.0
   qxtop  = 0.0
   
   qc_new = 0.0
   qa_new = 0.0

   mineddytau =  min_adj_time * dt

!-----------------------------------------------------------------------
!
!      compute height above surface
!

   nlev = size(t,3)
   nlat = size(t,2)
   nlon = size(t,1)

   mask = 1.0
               
   if (present(kbot)) then
        do j=1,nlat
        do i=1,nlon
             z_surf(i,j) = z_half(i,j,kbot(i,j)+1)
        enddo
        enddo
   else
        z_surf(:,:) = z_half(:,:,nlev+1)
   end if

   do k = 1, nlev
        z_full_ag(:,:,k) = z_full(:,:,k) - z_surf(:,:)
        z_half_ag(:,:,k) = z_half(:,:,k) - z_surf(:,:)
   end do
   z_half_ag(:,:,nlev+1) = z_half(:,:,nlev+1) - z_surf(:,:)
   
   if (present(kbot)) then
        where (z_full_ag < 0)
           mask = 0.0
        endwhere
   end if
      
!-----------------------------------------------------------------------
!
!      Calculate saturation specific humidity and its temperature 
!      derivative, the effective latent heat
!
!      These are calculated according to the formulas:
!
!      qsl  = d622*esl/ [p_full  -  (1.-d622)*esl]
!
!      dqsldtl = d622*p_full*(desat/dT)/[p_full-(1.-d622)*esl]**2.
!
!       
!      where d622 = rdgas/rvgas; esl = saturation vapor pressure;
!      and desat/dT is the temperature derivative of esl. Note that: 
!
!              {             hlv          for t > tfreeze             }
!      hleff = { 0.05*(t-tfreeze+20.)*hlv + 0.05*(tfreeze-t)*hls      }
!              {                          for tfreeze-20.< t < tfreeze}
!              {             hls          for t < tfreeze-20.         }
!
!      This linear form is chosen because at Tfreeze-20. es = esi, and
!      at Tfreeze, es = esl, with linear interpolation in between.
!
!

   !calculate effective latent heat
   hleff = (min(1.,max(0.,0.05*(t       -tfreeze+20.)))*hlv + &
            min(1.,max(0.,0.05*(tfreeze -t          )))*hls)
   
   !calculate qsl and dqsldtl. return es for diagnostic use.
   call compute_qs (t-(hleff*(ql+qi)/cp_air), p_full, qsl, &
                                         esat = esl, dqsdT=dqsldtl)

!-----------------------------------------------------------------------
!
!      set up specific humidities and static energies  
!      compute airdensity

   qt      = qv + ql + qi
   
   sli     = cp_air*t + grav*z_full_ag - hleff*(ql + qi)
   sliv    = sli*(1+zvir*qt)
   density = p_full/rdgas/(t *(1.+d608*qv-ql-qi))              

!-----------------------------------------------------------------------
! 
!      big loop over points
!

   ibot = nlev

   do j=1,nlat
     npts = 0
     if (do_edt_dg(j+js-1) ) then       
       do nn=1,num_pts
         if ( js == j_edtprt(nn) .and.  &
              i_edtprt(nn) >= is .and. i_edtprt(nn) <= ie) then
           iloc(npts+1) = i_edtprt(nn) - is + 1
           jloc(npts+1) = j_edtprt(nn) - js + 1
           nsave(npts+1) = nn
           npts = npts + 1
         endif
       end do    ! (num_points)
     else
       ipt = 0
       jpt = 0
       column_match = .false.
     endif
     do i=1,nlon
       if (npts > 0) then
         do nn=1,npts

           ipt = iloc(nn)
           jpt = jloc(nn)
           if (i == ipt ) then
             column_match = .true.
             nnsave = nsave(nn)
             exit
           else
             column_match = .false.
           endif
         end do
         nn = nnsave
       else 
         column_match = .false.
         nn = 0
       endif
       !-----------------------------------------------------------
       ! should diagnostics be printed out for this column
       !
       !-----------------------------------------------------------
       !
       ! extract column data to pass to caleddy
       !
    
       if (present(kbot)) ibot = kbot(i,j)
    
       gb_t      (1:ibot  ) = t         (i,j,1:ibot  )
       gb_qv     (1:ibot  ) = max(0., qv(i,j,1:ibot  ) )
       gb_qc     (1:ibot  ) = max(0., ql(i,j,1:ibot  ) +          &
                                     qi(i,j,1:ibot  ) )
       gb_qt     (1:ibot  ) = max(0., qt(i,j,1:ibot  ) )
       gb_sli    (1:ibot  ) = sli       (i,j,1:ibot  )
       gb_sliv   (1:ibot  ) = sliv      (i,j,1:ibot  )
       gb_u      (1:ibot  ) = u         (i,j,1:ibot  )
       gb_v      (1:ibot  ) = v         (i,j,1:ibot  )
       gb_z_half (1:ibot+1) = z_half_ag (i,j,1:ibot+1)
       gb_z_full (1:ibot  ) = z_full_ag (i,j,1:ibot  )
       gb_p_half (1:ibot+1) = p_half    (i,j,1:ibot+1)
       gb_p_full (1:ibot  ) = p_full    (i,j,1:ibot  )
       gb_qsl    (1:ibot  ) = qsl       (i,j,1:ibot  )
       gb_esl    (1:ibot  ) = esl       (i,j,1:ibot  )
       gb_hleff  (1:ibot  ) = hleff     (i,j,1:ibot  )
       gb_density(1:ibot  ) = density   (i,j,1:ibot  )
       gb_dqsldtl(1:ibot  ) = dqsldtl   (i,j,1:ibot  )
       gb_tdtlw  (1:ibot  ) = tdtlw_in (i,j,1:ibot)
       gb_sigmas (1:ibot+1) = sigmas(is-1+i,js-1+j,1:ibot+1)
       kqfs                 = u_star(i,j)*q_star(i,j)
       khfs                 = u_star(i,j)*b_star(i,j)*gb_t(ibot)/ &
                              grav
       gb_u_star            = u_star(i,j)          

      !compute sigmas at surface
       gb_sigmas(ibot+1) =  ((mesovar*gb_qsl(ibot))**2.0) +     &
         ((q_star(i,j)-gb_dqsldtl(ibot)*b_star(i,j)*gb_t(ibot)/ &
         grav)**2.0)
       gb_sigmas(ibot+1) = sqrt( gb_sigmas(ibot+1) ) / ( 1. + &
                                 gb_hleff(ibot)*gb_dqsldtl(ibot)/cp_air )

       if (column_match) then
         call get_date(Time, year, month, day, hour, minute, second)
         mon = month_name (month)
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  '===================================='//&
                            '=================='
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  '               ENTERING EDT    '
         write (dpu,'(a)')  ' '
         write (dpu,'(a, i6,a,i4,i4,i4,i4)')  ' time stamp:',   &
                                       year, trim(mon), day, &
                                       hour, minute, second
         write (dpu,'(a)')  '  DIAGNOSTIC POINT COORDINATES :'
         write (dpu,'(a)')  ' '
         write (dpu,'(a,f8.3,a,f8.3)') ' longitude = ', deglon1(nn),&
                                       ' latitude  = ', deglat1(nn)
         write (dpu,'(a,i6,a,i6)')     ' global i =', i_edtprt_gl(nn), &
                                       ' global j = ', j_edtprt_gl(nn)
         write (dpu,'(a,i6,a,i6)')     ' processor i =', i_edtprt(nn),     &
                                       ' processor j = ',j_edtprt(nn)
         write (dpu,'(a,i6,a,i6)')     ' window    i =', ipt,          &
                                       ' window    j = ',jpt
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' sigmas at the surface .... '
         write (dpu,'(a)')  ' ' 
         write (dpu,'(a,f14.7,a)')  ' sigmas = ',1000.*             &
                                     gb_sigmas(ibot+1), ' g/kg'
         write (dpu,'(a,f14.7)  ')  ' acoef = ', 1./( 1.+           &
                                     gb_hleff(ibot)* gb_dqsldtl(ibot)/cp_air )            
         write (dpu,'(a,f14.7,a)')  ' sigmas/a = ', 1000.*          &
                                     gb_sigmas(ibot+1)*( 1. + gb_hleff(ibot)* &
                                     gb_dqsldtl(ibot)/cp_air ), ' g/kg'
         write (dpu,'(a,f14.7)'  )  ' mesovar = ', mesovar
         write (dpu,'(a,f14.7,a)')  ' mesovar*qsl = ',mesovar*      &
                                    gb_qsl(ibot)*1000.,' g/kg'
         write (dpu,'(a,f14.7,a)')  ' turb.fluct = ',1000.*         &
              (q_star(i,j)-gb_dqsldtl(ibot)*b_star(i,j)*gb_t(ibot)/ &
               grav),' g/kg'
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' k      T         u         v       '//&
                            '  qv        qt ' 
         write (dpu,'(a)')  '       (K)      (m/s)     (m/s)     '//&
                            '(g/kg)    (g/kg)'
         write (dpu,'(a)')  '------------------------------------'//&
                            '-----------------'
         write (dpu,'(a)')  ' '
         do kk = nlev-n_print_levels,nlev
           write(dpu,18) kk,gb_t(kk),gb_u(kk),gb_v(kk),1000.*    &
                         gb_qv(kk), 1000.*gb_qt(kk)
         end do
18       format(1X,i2,1X,5(f9.4,1X))
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' k      qa        qc      sli/cp_air    '//&
                            'sliv/cp_air    tdtlw'
         write (dpu,'(a)')  '                (g/kg)     (K)      '//&
                            ' (K)      (K/day)'
         write (dpu,'(a)')  '------------------------------------'//&
                            '-----------------'
         write (dpu,'(a)')  ' '
         do kk = nlev-n_print_levels,nlev
           write(dpu,19) kk,qa(i,j,kk),1000.*gb_qc(kk),          &
                         gb_sli(kk)/cp_air,gb_sliv(kk)/cp_air,gb_tdtlw(kk)*86400.
         enddo    
19       format(1X,i2,1X,5(f9.4,1X))
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' k   z_full    z_half    p_full    p'//&
                            '_half    sigmas'
         write (dpu,'(a)')  '      (m)      (m)        (mb)      '//&
                            '(mb)     (g/kg)'
         write (dpu,'(a)')  '------------------------------------'//&
                            '-----------------'
         write (dpu,'(a)')  ' '
         do kk = nlev-n_print_levels,nlev
           write(dpu,19) kk,gb_z_full(kk),gb_z_half(kk+1),       &
                         gb_p_full(kk)/100.,gb_p_half(kk+1)/100.,1000.*   &
                         gb_sigmas(kk+1)
         enddo
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' k     esl       qsl      dqsldtl   '//&
                            ' hleff '
         write (dpu,'(a)')  '       (mb)    (g/kg)     (g/kg/K)  '//&
                            '(MJ/kg)' 
         write (dpu,'(a)')  '------------------------------------'//&
                            '-------'
         write (dpu,'(a)')  ' '
         do kk = nlev-n_print_levels,nlev
           write(dpu,19) kk, gb_esl(kk)/100.,gb_qsl(kk)*1000.,   &
                         gb_dqsldtl(kk)*1000.,gb_hleff(kk)/1.0e+06
         enddo
         write (dpu,'(a)')  ' '
       end if
       
       ! call the initialization routine for variables needed 
       ! by caleddy
         
       call trbintd(gb_t,       gb_qv,     gb_qt,     gb_qc,      &
                    gb_sli,     gb_sliv,   gb_u,      gb_v,       &
                    gb_z_full,  gb_z_half, gb_p_full, gb_p_half,  &
                    gb_sigmas,  gb_qsl,    gb_esl,    gb_hleff,   &
                    gb_dqsldtl, gb_slisl,  gb_qtsl,   gb_qxtop,   &
                    gb_qxbot,   gb_qltop,  gb_sfuh,   gb_sflh,    &
                    gb_qc_new,  gb_qa_new, gb_chu,    gb_chs,     &
                    gb_cmu,     gb_cms,    gb_n2,     gb_s2,      &
                    gb_ri    )
            
       !-----------------------------------------------------------
       !
       ! call caleddy
       !

            
       call caleddy(gb_u_star, kqfs,      khfs,      gb_sli,      &
                    gb_qt,     gb_qc_new, gb_qa_new, gb_sliv,     &
                    gb_u,      gb_v,      gb_p_full, gb_z_full,   &
                    gb_sfuh,   gb_sflh,   gb_slisl,  gb_qtsl,     &
                    gb_tdtlw,  gb_hleff,  gb_density,gb_qsl,      &
                    gb_dqsldtl,gb_qltop,  gb_p_half, gb_z_half,   &
                    gb_chu,    gb_chs,    gb_cmu,    gb_cms,      &
                    gb_n2,     gb_s2,     gb_ri,     gb_pblh,     &
                    gb_turbtype,gb_k_t,   gb_k_m,    gb_tke,      &
                    gb_leng,   gb_bprod,  gb_sprod,  gb_trans,    &
                    gb_diss,   gb_isturb, gb_tauinv, gb_sigmas,   &
                    gb_radf,   gb_sh,     gb_sm,     gb_gh,       &
                    gb_evhc)
 
       !-----------------------------------------------------------
       !
       ! paste back outputs
       !

       pblh(i,j)              = gb_pblh
       k_t(i,j,1:ibot)        = gb_k_t(1:ibot)
       k_m(i,j,1:ibot)        = gb_k_m(1:ibot)
       tke(i,j,1:ibot+1)      = gb_tke(1:ibot+1)
    
       turbtype(:,i,j,1:ibot+1) = gb_turbtype(:,1:ibot+1)
       isturb(i,j,1:ibot)       = gb_isturb(1:ibot)
       n2(i,j,1:ibot+1)         = gb_n2(1:ibot+1)
       s2(i,j,1:ibot+1)         = gb_s2(1:ibot+1)
       ri(i,j,1:ibot+1)         = gb_ri(1:ibot+1)    
       tauinv(i,j,1:ibot+1)     = gb_tauinv(1:ibot+1)
       leng(i,j,1:ibot+1)       = gb_leng(1:ibot+1)
       bprod(i,j,1:ibot+1)      = gb_bprod(1:ibot+1)
       sprod(i,j,1:ibot+1)      = gb_sprod(1:ibot+1)
       trans(i,j,1:ibot+1)      = gb_trans(1:ibot+1)
       diss(i,j,1:ibot+1)       = gb_diss(1:ibot+1)
       radf(i,j,1:ibot+1)       = gb_radf(1:ibot+1)
       sh(i,j,1:ibot+1)         = gb_sh(1:ibot+1)
       sm(i,j,1:ibot+1)         = gb_sm(1:ibot+1)
       gh(i,j,1:ibot+1)         = gb_gh(1:ibot+1)
       evhc(i,j,1:ibot+1)       = gb_evhc(1:ibot+1)
                                                    
       slislope(i,j,1:ibot  )   = gb_slisl(1:ibot  )
       qtslope (i,j,1:ibot  )   = gb_qtsl (1:ibot  )
       qxtop   (i,j,1:ibot  )   = gb_qxtop(1:ibot  )
       qxbot   (i,j,1:ibot  )   = gb_qxbot(1:ibot  )
       qc_new  (i,j,1:ibot  )   = gb_qc_new(1:ibot  )
       qa_new  (i,j,1:ibot  )   = gb_qa_new(1:ibot  )
    
       sigmas(is-1+i,js-1+j,1:ibot+1) = gb_sigmas (1:ibot+1) 
    
       !determine maximum altitude that the very stable pbl
       !is permitted to operate.  this is set to the half
       !level altitude just beneath the first turbulent
       !level
       stbltop(i,j) = 0.
       topfound = .false.
     
       kk = ibot + 1
       do while (.not.topfound.and.kk.gt.1)
         kk = kk - 1
         if (gb_isturb(kk).gt.0.5) topfound = .true.
       enddo 
       stbltop(i,j) = gb_z_half(kk+1)
       
       if (column_match) then
         write (dpu,'(a)') ' '
         write (dpu,'(a,f14.7,a)') ' stbltop = ',stbltop(i,j), ' meters'
         write (dpu,'(a)') ' '
       end if
    
     enddo
   enddo

!-----------------------------------------------------------------------
! 
!      diagnose cloud fraction and condensate tendencies
!
      
!----------------------------------------------------------------
! If the interface at the top or the bottom of the layer is part
! of a turbulent layer then set the eddytau to the minimum 
! timescale of the two.
!
       
   do k = 1, nlev
     tauinvtmp(:,:,k) = max(tauinv(:,:,k),tauinv(:,:,k+1))
   enddo
       
   where (tauinvtmp .gt. 1.e-10 .and. isturb .gt. 0.5)
     eddytau = max ( 1./tauinvtmp, mineddytau )
   elsewhere
     eddytau = missing_value
   endwhere
       
   qaturb(is:ie,js:je,:)   = qa_new
   qcturb(is:ie,js:je,:)   = qc_new
   tblyrtau(is:ie,js:je,:) = eddytau
             
!-----------------------------------------------------------------------
! 
!      blend in monin-obukhov similarity theory mixing coefficients at 
!      interface levels which are inside the surface layer
!

   sfclyr_h     = frac_sfclyr*pblh
   sfclyr_h_max = maxval(sfclyr_h)

   kk = nlev
   do k = 2, nlev
     if (minval(z_half_ag(:,:,k)) < sfclyr_h_max) then
       kk = k
       exit
     end if
   end do
   k_m_mo = 0.0
   k_t_mo = 0.0
   call mo_diff(z_half_ag(:,:,kk:nlev), u_star, b_star, &
                k_m_mo(:,:,kk:nlev), k_t_mo(:,:,kk:nlev))

   do k = kk, nlev
     where(z_half_ag(:,:,k) < sfclyr_h(:,:) .and. &
           z_half_ag(:,:,k) > 0.)
       k_t(:,:,k) = k_t_mo(:,:,k)
       k_m(:,:,k) = k_m_mo(:,:,k)
     endwhere
   enddo

!-----------------------------------------------------------------------
! 
!      Diagnostics
!

   if ( id_fq_cv_int > 0 .or. id_fq_cv_top > 0 .or. id_n2 > 0 .or. & 
        id_fq_turb   > 0 .or. id_diss      > 0 .or. id_s2 > 0 .or. &
        id_fq_cv_bot > 0 .or. id_fq_st     > 0 .or. id_ri > 0 .or. &
        id_eddytau   > 0 .or. id_tauinv    > 0 .or.                &
        id_leng      > 0 .or. id_qaedt     > 0 .or.                &
        id_bprod     > 0 .or. id_sprod     > 0 .or.                &
        id_trans     > 0 .or. id_qcedt     > 0 .or.                &
        id_sigmas    > 0 ) then  
      
     mask3(:,:,1:(nlev+1)) = 1.
     if (present(kbot)) then
       where (z_half_ag < 0.)
         mask3(:,:,:) = 0.
       end where
     endif
     
     if ( id_fq_cv_int > 0 ) then
       where (turbtype(2,:,:,:) .eq. 1) 
         tmpdat = 1.
       elsewhere
         tmpdat = 0.
       end where
       used = send_data (id_fq_cv_int,tmpdat,time, is, js, 1,&
                         rmask=mask3 )
     end if

     if ( id_fq_cv_top > 0 ) then
       where (turbtype(4,:,:,:) .eq. 1) 
         tmpdat = 1.
       elsewhere
         tmpdat = 0.
       end where
       used = send_data (id_fq_cv_top,tmpdat,time, is, js, 1,&
                         rmask=mask3 )
     end if
            
     if ( id_fq_cv_bot > 0 ) then
       where (turbtype(3,:,:,:) .eq. 1) 
         tmpdat = 1.
       elsewhere
         tmpdat = 0.
       end where
       used = send_data (id_fq_cv_bot,tmpdat,time, is, js, 1,&
                         rmask=mask3 )
     end if
         
     if ( id_fq_st > 0 ) then
       where (turbtype(1,:,:,:) .eq. 1) 
         tmpdat = 1.
       elsewhere
         tmpdat = 0.
       end where                 
       used = send_data ( id_fq_st, tmpdat, time, is, js, 1, &
                          rmask=mask3 )
     end if

     if ( id_fq_turb > 0 ) then
       used = send_data ( id_fq_turb, isturb, time, is, js,1,&
                          rmask=mask )
     end if
                     
     if ( id_n2 > 0 ) then
       used = send_data ( id_n2, n2, time, is, js, 1,        &
                          rmask=mask3 )
     end if
         
     if ( id_s2 > 0 ) then
       used = send_data ( id_s2, s2, time, is, js, 1,        &
                          rmask=mask3 )
     end if
         
     if ( id_ri > 0 ) then
       used = send_data ( id_ri, ri, time, is, js, 1,        &
                          rmask=mask3 )
     end if
         
     if ( id_leng > 0 ) then
       used = send_data ( id_leng, leng, time, is, js, 1,    &
                          rmask=mask3 )
     end if
            
     if ( id_bprod > 0 ) then
       used = send_data ( id_bprod, bprod, time, is, js, 1,  &
                          rmask=mask3 )
     end if
         
     if ( id_sprod > 0 ) then
       used = send_data ( id_sprod, sprod, time, is, js, 1,  &
                          rmask=mask3 )            
     end if
    
     if ( id_trans > 0 ) then
          used = send_data ( id_trans, trans, time, is, js, 1,  &
                             rmask=mask3 )            
     end if
    
     if ( id_diss > 0 ) then
          used = send_data ( id_diss, diss, time, is, js, 1,    &
                             rmask=mask3 )            
     end if
    
     if ( id_radf > 0 ) then
          used = send_data ( id_radf, radf, time, is, js, 1,    &
                             rmask=mask3 )            
     end if
    
     if ( id_sh > 0 ) then
          used = send_data ( id_sh, sh, time, is, js, 1,        &
                             rmask=mask3 )            
     end if
    
     if ( id_sm > 0 ) then
          used = send_data ( id_sm, sm, time, is, js, 1,        & 
                             rmask=mask3 )            
     end if
    
     if ( id_gh > 0 ) then
          used = send_data ( id_gh, gh, time, is, js, 1,        &
                             rmask=mask3 )            
     end if
    
     if ( id_evhc > 0 ) then
          used = send_data ( id_evhc, evhc, time, is, js, 1,    &
                             rmask=mask3 )            
     end if
    
     if ( id_tauinv > 0 ) then
          used = send_data ( id_tauinv, tauinv, time, is, js, 1,&
                             rmask=mask3 )
     end if
     
     if ( id_eddytau > 0 ) then
          used = send_data ( id_eddytau, eddytau, time, is, js, &
                             1, rmask=mask )
     end if
     
     if ( id_sigmas > 0 ) then
          used = send_data ( id_sigmas, sigmas(is:ie,js:je,:),  &
                             time, is, js, 1, rmask=mask3 )
     end if
            
     if (id_qaedt > 0) then
       tmpdat = 0.0
       do k = 1, nlev   
         tmpdat(:,:,k) = isturb(:,:,k)*qa_new(:,:,k)
       enddo
       used = send_data ( id_qaedt, tmpdat(:,:,1:nlev),      &
                          time, is, js, 1, rmask=mask )                 
     end if

     if (id_qcedt > 0) then
       tmpdat = 0.0
       do k = 1, nlev   
         tmpdat(:,:,k) = isturb(:,:,k)*qc_new(:,:,k)
       enddo
       used = send_data ( id_qcedt, tmpdat(:,:,1:nlev),      &
                          time, is, js, 1, rmask=mask )      
     end if
 
   end if  ! do diagnostics if


!-----------------------------------------------------------------------
! 
!      close edt output file if data was written for this window

!!RSH  if (ipt .gt. 0 .and. jpt .gt. 0) call Close_File (dpu)
       
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine edt

!
!======================================================================= 





!======================================================================= 
!
!      subroutine edt_end
!        
!
!      this subroutine writes out the restart field
!        

subroutine edt_end()

!-----------------------------------------------------------------------
!
!      variables
!
!      --------
!      internal
!      --------
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!      write out restart file
!
       call save_restart(edt_restart)

!-----------------------------------------------------------------------
! 
!      close edt output file if data was written for this window

       if (do_print ) call Close_File (dpu)
       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 
!      subroutine end
!
       module_is_initialized = .false.
end subroutine edt_end

!
!======================================================================= 




!======================================================================= 
!======================================================================= 
!
!
!      Grenier-Bretherton turbulence subroutines follow. These include:
!
!
!      sfdiag             diagnoses the fraction of each interval
!                         between layer midpoints which is saturated;
!                         this fraction is used to calculate the 
!                         proper buoyancy frequency, n2.
!
!      trbintd            calculates the subgrid vertical variability
!                         to sl, the liquid water static energy, and
!                         qt, the total water specific humidity.
!
!                         also calculates the richardson number ri 
!                         from n2:
!
!                         n2 = ch*dslidz +  cm*dqtdz
!
!                         where
!
!                         dslidz and dqtdz are the vertical slopes
!                         of sl and qt in the interior of the grid
!                         box and ch and cm are thermodynamic 
!                         coefficients which are solely functions of 
!                         temperature and pressure, but have very 
!                         different values in saturated and unsaturated 
!                         air.
!                          
!      caleddy            main subroutine which calculates diffusivity
!                         coefficients
!
!      exacol          determines whether the column has adjacent 
!                         regions where Ri < 0 (unstable layers or ULs) 
!                         and determine the indices kbase, ktop which 
!                         delimit these unstable layers : 
!
!                         ri(kbase) > 0 and ri(ktop) > 0, 
!                         but ri(k)     < 0 for ktop < k < kbase. 
!
!      zisocl             solves for mean tke <e>, W, Sh, and Sm for 
!                         each convective layer.  here,
!
!                         W = leng**2 * (-Sh*n2 + Sm*s2), where
!
!                         leng = mixing length,
!                         s2   = shear vector magnitude
!                         Sh   = Galperin stab. function for buoyancy
!                         Sm   = Galperin stab. function for momentum
!
!                         in addition, it merges adjacent convective 
!                         layers when the intervening stable layers 
!                         would not consume too much tke for the 
!                         combined convective layer
!
!
!======================================================================= 
!======================================================================= 




!======================================================================= 

subroutine sfdiag (qsl, esl, dqsldtl, hleff, qt, qtslope,sli, slislope,&
                   p_full, z_full, p_half, z_half, sdevs, qxtop, qxbot,&
                   qltop, sfuh, sflh, qc_new, qa_new, sfi)
   
!----------------------------------------------------------------------- 
! 
!      Purpose: 
!      Interface for computing cloud variables for use by turb. scheme
! 
!      Authors: B. Stevens and C. Bretherton (August 2000)
! 
!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
!      field 1-d arrays on model full levels, reals dimensioned
!      (1:nlev), index running from top of atmosphere to bottom
!      
!      qsl         saturation spec. humidity at the liquid-ice
!                  water temperature (kg water/kg air)
!      esl         saturation vapor pressure at the liquid-ice 
!                  water temperature (Pa)
!      dqsldtl     temperature derivative of qsl (kg water/kg air/K)
!      hleff       effective latent heat of vaporization (J/kg)        
!      qt          total water specific humidity (kg water/kg air)
!      qtslope     qt slope wrt pressure in thermo layer (kg/kg/Pa)
!      sli         ice-liq water static energy (J/kg) 
!      slislope    sli slope wrt pressure in thermo layer (J/kg/Pa)
!      p_full      pressure (Pa)
!      z_full      height of full level above the surface (m)
!
!      the following fields are on the model half levels, 
!      dimension(1:nlev+1)
!
!      p_half      pressure at half levels (Pa)
!      z_half      height of half model levels above the surface (m)
!      sdevs       standard deviation of water perturbation
!                  (kg water/kg air)
!
!      ------
!      output
!      ------
!
!      qxtop       saturation excess at top of layer (kg wat/kg air)
!      qxbot       saturation excess at bottom of layer (kg wat/kg air)
!      qltop       cloud liquid at top of layer (kg wat/kg air)
!      sfuh        saturated fraction in upper half-layer
!      sflh        sflh saturated fraction in lower half-layer
!      qc_new      thermodynamically diagnosed condensate value
!      qa_new      thermodynamically diagnosed cloud fraction
!
!      the following fields are on the model half levels, 
!      dimension(1:nlev+1)
!
!      sfi         interfacial saturated fraction
!


real, intent(in) , dimension(:) :: qsl, esl, dqsldtl, hleff
real, intent(in) , dimension(:) :: qt, qtslope, sli, slislope
real, intent(in) , dimension(:) :: p_full, z_full
real, intent(in) , dimension(:) :: p_half, z_half, sdevs 
real, intent(out), dimension(:) :: qxtop, qxbot, qltop, sfuh, sflh
real, intent(out), dimension(:) :: qc_new, qa_new
real, intent(out), dimension(:) :: sfi

integer :: k                      ! vertical index
integer :: nlev                   ! number of vertical levels
real    :: slitop, slibot         ! sli at top/bot of layer
real    :: qttop , qtbot          ! qt  at top/bot of layer
real    :: qsltop, qslbot         ! qsl at top/bot of layer
real    :: tlitop, tlibot         ! liq wat temp at top/bot of layer          
real    :: qxm                    ! sat excess at midpoint
real    :: qlm, qlbot             ! liq wat at midpoint, and bottom
real    :: tlim                   ! tli at midpoint
real    :: dqsldp                 ! pressure derivative of qsl
real    :: sigmasf                ! sdevs on full model levels
real    :: acoef                  ! 1./(1+L*dqsldT/cp_air)

!----------------------------------------------------------------------- 
!
!      code
!

   nlev   = size(p_full,1)
   sfi    = 0.0
   sfuh   = 0.0
   sflh   = 0.0    
   qc_new = 0.0
   qa_new = 0.0
   qltop  = 0.0
   
   if (column_match) then
     write (dpu,'(a)')  '-----------------------------------------'//&
                        '-------------'
     write (dpu,'(a)')  '        ENTERING SFDIAG                     '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k  tlim   dqsldp slitop tlitop  qttop  q'//&
                        'sltop  qxtop  qltop   qxm    qlm'
     write (dpu,'(a)')  '    (K)(g/kg/100mb)(K)    (K)    (g/kg) ('//&
                        'g/kg)  (g/kg) (g/kg) (g/kg) (g/kg)'
     write (dpu,'(a)')  '-----------------------------------------'//&
                        '----------------------------------'
     write (dpu,'(a)')  ' '       
21 format(1X,i2,1X,f6.2,1X,f7.2,1X,8(f6.2,1X))
   end if
         
   do k = 2,nlev

!-----------------------------------------------------------
! calculate midpoint liquid-ice water temperature
     tlim = (sli(k) - grav*z_full(k))/cp_air

!-----------------------------------------------------------
! calculate dqsl/dp
     dqsldp = -1.*qsl(k)*qsl(k)/d622/esl(k)
            
!-----------------------------------------------------------
! Compute saturation excess at top and bottom of layer k

!extrapolation to layer top
     slitop = sli(k) + slislope(k) * (p_half(k) - p_full(k))
     qttop  = qt (k) + qtslope (k) * (p_half(k) - p_full(k))
     slitop = min ( max ( 0.25*sli(k), slitop), 4.*sli(k))
     qttop  = min ( max ( 0.25* qt(k), qttop ), 4.* qt(k))
     tlitop = (slitop - grav*z_half(k))/cp_air   
     qsltop  = qsl(k) + dqsldp      * (p_half(k) - p_full(k)) + &
                        dqsldtl (k) * (tlitop    - tlim     )
     qsltop  = min ( max ( 0.25*qsl(k), qsltop ) , 4.*qsl(k) )
     qxtop(k) = qttop  - qsltop
     qltop(k)  = max( 0., qxtop(k)/( 1.+ hleff(k)*dqsldtl(k)/cp_air))

    !extrapolation to layer bottom
     slibot = sli(k) + slislope(k) * (p_half(k+1) - p_full(k))
     qtbot  = qt (k) + qtslope (k) * (p_half(k+1) - p_full(k))
     slibot = min ( max ( 0.25*sli(k), slibot), 4.*sli(k))
     qtbot  = min ( max ( 0.25* qt(k), qtbot ), 4.* qt(k))
     tlibot = (slibot - grav*z_half(k+1))/cp_air
     qslbot = qsl(k) + dqsldp      * (p_half(k+1) - p_full(k))+ &
                       dqsldtl (k) * (tlibot      - tlim     ) 
     qslbot  = min ( max ( 0.25*qsl(k), qslbot ) ,4.*qsl(k) )     
     qxbot(k) = qtbot  - qslbot
     qlbot  = max( 0., qxbot(k)/( 1.+hleff(k)*dqsldtl(k)/cp_air ) )

!-----------------------------------------------------------
! Compute saturation excess at midpoint of layer k           
     qxm    = qxtop(k) + (qxbot(k)-qxtop(k))*                   &
              (p_full(k)-p_half(k))/(p_half(k+1) - p_half(k))
     qlm  = max( 0., qxm / ( 1. +  hleff(k)*dqsldtl(k)/cp_air ) )

     if (column_match) then
       write(dpu,21) k,tlim,1000.*dqsldp*100.*100.,slitop/cp_air,    &
                     tlitop,1000.*qttop,1000.*qsltop,1000.*qxtop(k), &
                     1000.*qltop(k),1000.*qxm,   1000.*qlm
       write(dpu,21) k,tlim,1000.*dqsldp*100.*100.,slibot/cp_air,    &
                     tlibot,1000.*qtbot,1000.*qslbot,1000.*qxbot(k), &
                     1000.*qlbot,1000.*qxm,   1000.*qlm
       write (dpu,'(a)')  ' '
     end if

!-----------------------------------------------------------
!
! TWO WAYS OF CALCULATING SATURATION FRACTION AND QA AND QC
!
!
!          FIRST WAY:  USE GAUSSIAN CLOUD MODEL
!
!

     if (do_gaussian_cloud) then
   
!calculate sigmas on model full levels
       sigmasf = 0.5 * (sdevs(k) + sdevs(k+1))
       acoef = 1. / ( 1. + hleff(k)*dqsldtl(k)/cp_air )
!sigmasf = max ( sigmasf, acoef*mesovar*qsl(k) )
       sigmasf = acoef*mesovar*qsl(k)
       call gaussian_cloud(qxtop(k), qxm,      qxbot(k),     &
                           acoef,    sigmasf,  qa_new(k),    &
                           qc_new(k),sfuh(k),  sflh(k))

!------------------------------------------------------
! Combine with sflh (still for layer k-1) to get 
! interface layer saturation fraction
!
! N.B.:
!
! if sfuh(k)>sflh(k-1),sfi(k) = sflh(k-1)
! if sfuh(k)<sflh(k-1),sfi(k) = mean(sfuh(k),sflh(k-1))
!
       sfi(k) =  0.5 * ( sflh(k-1) + min(sflh(k-1),sfuh(k)) )
 
     else
    
!
!
!          SECOND WAY:  ORIGINAL BRETHERTON-GRENIER WAY
!
!

!------------------------------------------------------
! Compute saturation fraction sfuh(k) of the upper half 
! of layer k.

       if      ( (qxtop(k).lt.0.) .and. (qxm.lt.0.) ) then
         sfuh(k) = 0.  ! Upper half-layer unsaturated
       else if ( (qxtop(k).gt.0.) .and. (qxm.gt.0.) ) then
         sfuh(k) = 1.  ! Upper half-layer fully saturated
       else               ! Either qxm < 0 and qxtop > 0 
                          ! or vice versa
         sfuh(k) = max(qxtop(k),qxm) / abs(qxtop(k) - qxm)
       end if

!------------------------------------------------------
! Combine with sflh (still for layer k-1) to get      
! interfac layer sat frac
!
! N.B.:
!
! if sfuh(k)>sflh(k-1),sfi(k) = sflh(k-1)
! if sfuh(k)<sflh(k-1),sfi(k) = mean(sfuh(k),sflh(k-1))
!
       sfi(k) =  0.5 * ( sflh(k-1) + min( sflh(k-1), sfuh(k)))

!------------------------------------------------------
! Update sflh to be for the lower half of layer k.             

       if      ( (qxbot(k).lt.0.) .and. (qxm.lt.0.) ) then
         sflh(k) = 0.  ! Upper half-layer unsaturated
       else if ( (qxbot(k).gt.0.) .and. (qxm.gt.0.) ) then
         sflh(k) = 1.  ! Upper half-layer fully saturated
       else            ! Either qxm < 0 and qxbot > 0 or vice versa
         sflh(k) = max(qxbot(k),qxm) / abs(qxbot(k) - qxm)
       end if
         
         !------------------------------------------------------
           !Compute grid volume mean condensate and cloud fraction
       qc_new(k) = 0.5 * ( sfuh(k) * 0.5 * (qltop(k) + qlm) )&
                 + 0.5 * ( sflh(k) * 0.5 * (qlbot    + qlm) )
       qa_new(k) = 0.5 * ( sfuh(k) + sflh(k) )
   
     end if
    
   end do

       
       !set surface saturated fraction equal to sflh(nlev)
       sfi(nlev+1)  = sflh(nlev)! Sat frac in lowest half-layer. 

   if (column_match) then
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k     sfuh       sflh      sfi      qc_n'//&
                        'ew    qa_new'
     write (dpu,'(a)')  '                                     (g/kg)'
     write (dpu,'(a)')  '-----------------------------------------'//&
                        '------------'
     write (dpu,'(a)')  ' '
     do k = nlev-n_print_levels,nlev
       write(dpu,22) k,sfuh(k),sflh(k),sfi(k),1000.*qc_new(k),    &
                     qa_new(k)
     enddo    
22   format(1X,i2,1X,5(f9.4,1X))
     write (dpu,'(a)')  ' '
   end if
     
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine sfdiag

!
!======================================================================= 




!======================================================================= 
!
!      subroutine trbintd
!

 subroutine trbintd (t, qv, qt, qc, sli, sliv, u, v, z_full, z_half,   &
                     p_full, p_half, sdevs, qsl, esl, hleff, dqsldtl,  &
                     slislope, qtslope, qxtop, qxbot, qltop, sfuh,sflh,&
                     qc_new, qa_new, chu, chs, cmu, cms, n2, s2, ri)

!----------------------------------------------------------------------- 
! 
! Purpose: 
!  time dependent initialization
! 
! Method: 
!  Diagnosis of variables that do not depend on mixing assumptions or
!  PBL depth.
!
! Authors: B. Stevens (extracted from pbldiff, August, 2000)
!          C. Bretherton (Dec. 2000)
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!       variables
!
!       -----
!       input
!       -----
!       field 1-d arrays on model full levels, reals dimensioned
!       (1:nlev), index running from top of atmosphere to 
!       bottom
!           
!       t           temperature (K)
!       qv          water vapor spec. humidity (kg vapor/kg air)
!       qt          total water specific humidity (kg water/kg air)
!       qc          condensed water spec. humidity (kg cond/kg air)
!       sli         ice-liq water static energy (J/kg) 
!       sliv        ice-liq water virtual static energy (J/kg) 
!       u           zonal wind (m/s)
!       v           meridional wind (m/s) 
!       z_full      height of full level above the surface (m)
!       p_full      pressure (Pa)
!       sdevs       standard deviation of water perturbation s
!                   (kg water/kg air)
!       qsl         saturation spec. humidity at the liquid-ice
!                   water temperature (kg water/kg air)
!       esl         saturation vapor pressure at the liquid-ice 
!                   water temperature (Pa)
!       hleff       effective latent heat (J/kg condensate)
!       dqsldtl     temperature derivative of qsl (kg water/kg air/K)
!
!       the following fields are on the model half levels, 
!       dimension(1:nlev+1)
!       p_half      pressure at half levels (Pa)
!       z_half      height of half model levels above the surface (m)
!
!       ------
!       output
!       ------
!
!       slislope    sli slope wrt pressure in thermo layer (J/kg/Pa)
!       qtslope     qt slope wrt pressure in thermo layer (kg/kg/Pa)
!       qxtop       saturation excess at the top of the layer 
!                   (kg wat/kg air)
!       qxbot       saturation excess at the bottom of the layer 
!                   (kg wat/kg air)
!       qltop       liquid water at top of thermo layer (kg condensate/
!                   kg air)
!       sfuh        saturated fraction in upper half-layer
!       sflh        sflh saturated fraction in lower half-layer
!       qc_new      thermodynamically defined cloud condensate (kg/kg)
!       qa_new      thermodynamically defined cloud fraction 
!
!       the following fields are defined on model half levels,
!       dimension(1:nlev+1)
!
!       chu         heat var. coef for dry states (1/m)
!       chs         heat var. coef for sat states (1/m)
!       cmu         moisture var. coef for dry states (kg/kg)*(m/s*s)
!       cms         moisture var. coef for sat states (kg/kg)*(m/s*s)
!       n2          moist squared buoyancy freq (1/s*s)
!       s2          squared deformation, or shear vector mag. (1/s*s)
!       ri          gradient Richardson number
!
!

real, intent(in),    dimension (:) :: t, qv, qt, qc, sli, sliv, u, v,  &
                                      z_full, p_full, qsl, esl, hleff, &
                                      dqsldtl, sdevs, z_half, p_half
real, intent(out),   dimension (:) :: slislope, qtslope, qxtop, qxbot,   &
                                      sfuh, sflh, qc_new, qa_new, qltop, &
                                      chu, chs, cmu, cms, n2, s2, ri
  
! internal variables

integer            :: kdim        ! # of levels in the vertical
integer            :: k, km1, kp  ! level indexes
real               :: rdz         ! 1 / (delta z) between midpoints
real               :: dslidz      ! delta sli / delta z at interface
real               :: dqtdz       ! delta qt  / delta z at interface
real               :: ch          ! sfi weighted ch at the interface
real               :: cm          ! sfi weighted cm at the interface
real               :: product     ! temporary variable
real               :: dslidp_a    ! sli slope across interface above
real               :: dqtdp_a     ! qt  slope across interface above
real               :: dslidp_b    ! sli slope across interface below
real               :: dqtdp_b     ! qt  slope across interface below
real, dimension(size(t,1)) :: bfact ! buoyancy factor in n2 calculation
real, dimension(size(t,1)+1) :: sfi ! saturated fraction at interfaces

!----------------------------------------------------------------------- 
!
!  code
!

   kdim = size(t,1)

!-----------------------------------------------------------------------
! 
!  Thermodynamic coefficients for buoyancy flux - these
!  are calculated at midpoints; they will be averaged to interfaces,
!  where they will ultimately be used. At the surface, the coeff-
!  icients are taken from the lowest midpoint.
!
!  These formulas come from the following expression
!
!  grav* tv' / tv = ch * sli'   + cm * qt'
!
!  chu and cmu are the values for unsaturated air, whereas
!  chs and cms are the values for   saturated air.

   bfact       = grav/(t*(1.+zvir*qv - qc))
   chu(1:kdim) = (1. + zvir*qt)*bfact/cp_air
   chs(1:kdim) = ( (1. + (1. + zvir)*dqsldtl*t) / &
                   (1. + (hleff*dqsldtl/cp_air)   ) ) * bfact/cp_air
   cmu(1:kdim) = zvir  * bfact * t
   cms(1:kdim) = hleff * chs(1:kdim)  -  bfact * t

   chu(kdim+1) = chu(kdim)
   chs(kdim+1) = chs(kdim)
   cmu(kdim+1) = cmu(kdim)
   cms(kdim+1) = cms(kdim)
   
   if (column_match) then
     write (dpu,'(a)')  '--------------------------------------------'
     write (dpu,'(a)')  '        ENTERING TRBINTD                    '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k     Tv       bfact '
     write (dpu,'(a)')  '       (K)    (m/K/s**2)' 
     write (dpu,'(a)')  '------------------------'
     write (dpu,'(a)')  ' '
     do k = kdim-n_print_levels,kdim
       write(dpu,17) k,t(k)*(1.+zvir*qv(k) - qc(k)),bfact(k)
     end do
17   format(1X,i2,1X,2(f9.4,1X))
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k     cp_air*chu    cp_air*chs      cmu      cms'
     write (dpu,'(a)')  '     (m/K/s**2)(m/K/s**2) (m/s**2) (m/s**2)'
     write (dpu,'(a)')  '-------------------------------------------'
     write (dpu,'(a)')  ' '
     do k = kdim+1-n_print_levels,kdim+1
       write(dpu,20) k,cp_air*chu(k),cp_air*chs(k),cmu(k),cms(k)
     enddo    
20   format(1X,i2,1X,4(f9.4,1X))
     write (dpu,'(a)')  ' '
   end if
       
!-----------------------------------------------------------------------
!     
!  Compute slopes in conserved variab. sl, qt within thermo layer k. 
!  a indicates the 'above' gradient from layer k-1 to layer k and 
!  b indicates the 'below' gradient from layer k   to layer k+1.
!  We take the smaller (in absolute value) of these gradients
!  as the slope within layer k. If they have opposite sign, gradient
!  in layer k is taken to be zero.
!
!  Slopes at endpoints determined by extrapolation

   slislope(kdim) = (sli   (kdim) - sli   (kdim-1))/ &
                    (p_full(kdim) - p_full(kdim-1))
   qtslope (kdim) = (qt    (kdim) - qt    (kdim-1))/ &
                    (p_full(kdim) - p_full(kdim-1))
   slislope(1)    = (sli   (2)    - sli   (2)     )/ &
                    (p_full(2)    - p_full(1)     )     
   qtslope (1)    = (qt    (2)    - qt    (2)     )/ &
                    (p_full(2)    - p_full(1)     ) 
   dslidp_b        = slislope(1)
   dqtdp_b         = qtslope (1)
    
   do k = 2, kdim-1
   
     kp = k + 1
     dslidp_a  = dslidp_b
     dqtdp_a   = dqtdp_b
     dslidp_b  = (sli(kp)-sli(k))/(p_full(kp)-p_full(k))
     dqtdp_b   = (qt (kp)-qt (k))/(p_full(kp)-p_full(k))
     product   = dslidp_a*dslidp_b
     if (product .le. 0.) then 
       slislope(k) = 0.
     else if (product.gt.0. .and. dslidp_a.lt.0.) then 
       slislope(k) = max(dslidp_a,dslidp_b)
     else if (product.gt.0. .and. dslidp_a.gt.0.) then 
       slislope(k) = min(dslidp_a,dslidp_b)
     end if
 
     product   = dqtdp_a*dqtdp_b
     if (product .le. 0.) then 
       qtslope (k) = 0.
     else if (product.gt.0. .and. dqtdp_a.lt.0.) then 
       qtslope (k) = max(dqtdp_a,dqtdp_b)
     else if (product.gt.0. .and. dqtdp_a.gt.0.) then 
       qtslope (k) = min(dqtdp_a,dqtdp_b)
     end if
   
   end do ! loop over k


   if (column_match) then
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k   slislope   qtslope  '
     write (dpu,'(a)')  '    (K/100mb)(g/kg/100mb)'
     write (dpu,'(a)')  '-------------------------'
     write (dpu,'(a)')  ' '
     do k = kdim-n_print_levels,kdim
       write(dpu,17) k,slislope(k)*100.*100./cp_air,                  &
                     1000.*qtslope(k)*100.*100.
     enddo
     write (dpu,'(a)')  ' '
   end if

!-----------------------------------------------------------------------
!     
!      Compute saturation fraction in the interfacial layers for use in
!      buoyancy flux computation.

   call sfdiag(qsl, esl, dqsldtl, hleff, qt, qtslope,sli, slislope,&
               p_full, z_full, p_half, z_half, sdevs, qxtop, qxbot,&
               qltop, sfuh, sflh, qc_new, qa_new, sfi)

!-----------------------------------------------------------------------
!     
!  Compute shear squared (s2), squared buoyancy frequency (n2) and 
!  Ri.  For the n2 calculation use gradients of sl and qt, weighted 
!  according to sfi, the fraction of the interfacial layer that is 
!  saturated.
!
!  This loop has to be done in increasing levels to avoid over-
!  writing chu, etc. arrays to interface values before we are done 
!  with their midpoint values.
!
!  Note that n2,s2,ri are set to zero at interfaces k = 1 and
!  k = kdim + 1, where they are not used.

   n2(kdim+1) = 0.
   s2(kdim+1) = 0.
   ri(kdim+1) = 0.

   if (column_match) then
     write (dpu,'(a)')  '---------------------------------------------'
     write (dpu,'(a)')  '                 IN TRBINTD'
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k    rdz   dslidz   dqtdz    chu     chs'//&
                        '     ch       cmu    cms      cm'
     write (dpu,'(a)')  '     (1/km) (K/km) (g/kg/km)       (m/K/s'//&
                        '**2)              (m/s**2)'
     write (dpu,'(a)')  '-----------------------------------------'//&
                        '----------------------------------'
     write (dpu,'(a)')  ' '       
22   format(1X,i2,1X,9(f7.3,1X))
   end if
      
   do k = kdim, 2, -1
     km1     = k - 1
     rdz     = 1. / (z_full(km1) - z_full(k))
     dslidz  = (sli(km1) - sli(k)) * rdz
     dqtdz   = (qt (km1) - qt (k)) * rdz 
     chu(k)  = (chu(km1) + chu(k))*0.5
     chs(k)  = (chs(km1) + chs(k))*0.5
     cmu(k)  = (cmu(km1) + cmu(k))*0.5
     cms(k)  = (cms(km1) + cms(k))*0.5
     ch      = chu(k)*(1.-sfi(k)) + chs(k)*sfi(k)
     cm      = cmu(k)*(1.-sfi(k)) + cms(k)*sfi(k)
     n2(k)   = ch*dslidz +  cm*dqtdz
     s2(k)   = ((u(km1)-u(k))**2 + (v(km1)-v(k))**2)*(rdz**2)
     s2(k)   = max(ntzero,s2(k))
     ri(k)   = n2(k) / s2(k)
    
     if (column_match) then
       write(dpu,22) k,1000.*rdz,1000.*dslidz/cp_air,            &
                     1000.*1000.*dqtdz,cp_air*chu(k),cp_air*chs(k),cp_air*ch, &
                     cmu(k),cms(k),cm
     end if
    
   end do
   n2(1) = 0.
   s2(1) = 0.
   ri(1) = 0.
      
   if (column_match) then
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k    sqrt(n2)  sqrt(s2)    ri'
     write (dpu,'(a)')  '      (1/hr)    (1/hr) '
     write (dpu,'(a)')  '---------------------------------'
     write (dpu,'(a)')  ' '
     do k = kdim-n_print_levels,kdim
       write(dpu,23) k,n2(k)*sqrt(abs(n2(k)))*3600./max(small,    &
                     abs(n2(k))),sqrt(s2(k))*3600.,ri(k)
     enddo
23   format(1X,i2,1X,3(f9.3,1X))
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  '--------------------------------------------'
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
   end if

!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine trbintd

!
!======================================================================= 




!======================================================================= 
!
!      subroutine exacol
!        
     
subroutine exacol(ri, bflxs, ktop, kbase, ncvfin) 

!----------------------------------------------------------------------- 
!
! object : determine whether the column has adjacent regions where 
!          Ri < 0 (unstable layers or ULs) and determine the indices 
!          kbase, ktop which delimit these unstable layers : 
!          ri(kbase) > 0 and ri(ktop) > 0, but 
!          ri(k) < 0 for ktop < k < kbase. 
!
! author : H. Grenier    05/2000, 
!          C. Bretherton 08/2000
!
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!
!

real,    intent(in),  dimension(:) :: ri     ! Moist gradient Ri. #
real,    intent(in)                :: bflxs  ! Surface buoyancy flux 
                                             ! (m*m)/(s*s*s)
integer, intent(out), dimension(:) :: kbase  ! vertical index of UL base
integer, intent(out), dimension(:) :: ktop   ! vertical index of UL top
integer, intent(out)               :: ncvfin ! number of ULs

!
! internal variables
!
      
integer :: k, kdim, ncv         
real    :: riex(size(ri,1)) ! Column Ri profile extended to surface
                            ! by taking riex > rimaxentr for bflxs < 0
    !           riex < rimaxentr for bflxs > 0.


!-----------------------------------------------------------------------
!
!  Initialize variables

   ncvfin    = 0
   kdim      = size(ri,1) - 1
   ktop(:)   = 0
   kbase(:)  = 0
   
   riex(1:kdim) = ri(1:kdim)
   riex(kdim+1) = rimaxentr-bflxs ! Allows consistent treatment
                                  ! of surface with other intrfcs.
   ncv = 0

!-----------------------------------------------------------------------
!
!  Work upward from surf interface

   k = kdim+1
                   
   do while ( k.gt.2 )
   
     if (riex(k) .lt. rimaxentr) then 
             
 !------------------------------------------------------
 !
 ! A new convective layer has been found.
 ! Define kbase as interface below first unstable one
 ! then decrement k until top unstable level is found.
 ! Set ktop to the first interface above unstable layer. 
 
       ncv = ncv + 1
       kbase(ncv) = min(k+1,kdim+1)
       do while (riex(k) .lt. rimaxentr .and. k.gt.2)
         k = k-1
       end do  
       ktop(ncv) = k
     else
                 
!------------------------------------------------------
!
! Keep on looking for a CL.
              
       k = k-1
       
     end if
       
   end do
       
!--------------------------
!
! Set total number of CLs

   ncvfin = ncv   

   if (column_match) then
     write (dpu,'(a)')  '--------------------------------------------'
     write (dpu,'(a)')  '                 IN EXACOL'
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  '  k     riex'
     write (dpu,'(a)')  '-------------'
     write (dpu,'(a)')  ' '
     do k = kdim+1-n_print_levels,kdim+1
       write(dpu,28) k,riex(k)
     enddo
28 format(1X,i2,1X,f9.3)
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a,i3)')  ' ncvfin = ',ncvfin
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  'ncv  ktop  kbase'
     write (dpu,'(a)')  '--------------------'
     write (dpu,'(a)')  ' '
     do ncv = 1, ncvfin
       write(dpu,29) ncv, ktop(ncv), kbase(ncv)
     enddo
29 format(1X,i2,3X,i2,5X,i2)
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
   end if

!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine exacol

!
!======================================================================= 
       

!=======================================================================
!
!      subroutine zisocl
!        
     
subroutine zisocl(u_star, bflxs, tkes, zm, ql, zi, n2, s2, ri, ncvfin, &
                  kbase, ktop, belongcv, ebrk, wbrk, ghcl, shcl, smcl, &
  lbrk)

!-----------------------------------------------------------------------
!
! object : find <e>, <W>, <Sh>, <Sm>, ktop(ncv), accounting for the 
!          presence of stably stratified layers inside the convective 
!          layer(CL, to be defined) but with r2 > ratinv.
!
!          Re-arrange the indexing of arrays kbase/ktop if some CLs are 
!          found to be coupled such that ncv defines the index of each 
!          CL increasing with height.
!
! author : H. Grenier 05/08/2000
!
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!
!
   
real, intent(in)               :: u_star ! friction velocity (m/s)   
real, intent(in)               :: bflxs  ! Surface buoyancy flux 
                                         ! (m*m*m)/(s*s)
real, intent(in)               :: tkes   ! TKE at the surface (m2/s2)
real, intent(in), dimension(:) :: zm ! Layer midpoint height (m)
real, intent(in), dimension(:) :: ql ! condensate spec. hum. (kg/kg)
real, intent(in), dimension(:) :: zi ! Interface height (m)
real, intent(in), dimension(:) :: n2 ! Moist squared buoy freq (s-2)
real, intent(in), dimension(:) :: s2 ! Shear deformation (s-2)
real, intent(in), dimension(:) :: ri ! Gradient Richardson number
    
integer, intent(inout)               :: ncvfin ! Total number of CLs
integer, intent(inout), dimension(:) :: kbase  ! Vert index of CL base
integer, intent(inout), dimension(:) :: ktop   ! Vert index of CL top

logical, intent(out), dimension(:) :: belongcv ! = T if flux level in CL
real,    intent(out), dimension(:) :: ebrk  ! vert ave. of TKE  in CL
real,    intent(out), dimension(:) :: wbrk  !   "       of W^2  "
real,    intent(out), dimension(:) :: ghcl  !   "       of Gh   "
real,    intent(out), dimension(:) :: shcl  !   "       of Sh   "
real,    intent(out), dimension(:) :: smcl  !   "       of Sm   "
real,    intent(out), dimension(:) :: lbrk  ! CL depth not within entr
                                            ! layers
  
! internal variables

logical :: extend       ! True if CL is extended in zisocl
logical :: bottom       ! True if CL base at surface(kb = kdim+1)
integer :: ncv          ! Index enumerating convective layers in col
integer :: incv
integer :: k
integer :: kdim         ! number of full vertical levels
integer :: kb           ! Local index for kbase
integer :: kt           ! Local index for ktop
integer :: ncvinit      ! Value of ncv at routine entrance 
integer :: cntu         ! counts upward no. of merged CLs
integer :: cntd         ! counts downward  "          "
integer :: kbinc        ! Index for incorporating underlying CL
integer :: ktinc        ! Index for incorporating  overlying CL
real    :: ebar
real    :: wint
real    :: dwinc
real    :: dzinc
real    :: dwsurf
real    :: gh
real    :: sh
real    :: sm
real    :: l2n2         ! Vert. integral of l^2N^2 over CL
real    :: l2s2         ! Vert. integral of l^2S^2 over CL
real    :: dl2n2        ! Vert. int. of l^2N^2 over incorp. layer
real    :: dl2s2        ! Vert. int. of l^2S^2 over incorp. layer
real    :: lint         ! CL depth excluding entrainment layers
real    :: lbulk        ! Depth of the convective layer
real    :: lz           ! Turbulent length scale
real    :: ricl         ! Ri Number for the whole convective layer
real    :: zbot         ! Height of CL base
real    :: l2rat        ! Square of ratio of actual to initial CL depth
real    :: tmpr

!-----------------------------------------------------------------------
!
!  Initialize variables

   kdim = size(ql,1)
   ncv = 1 

   if (column_match) then
     write (dpu,'(a)')  '-------------------------------------------'
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  '        ENTERING ZISOCL               '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  '  k    lint      lz      l2n2         l2s'//&
                        '2' 
     write (dpu,'(a)')  '       (m)       (m)    (m3/s2)      (m3/'//&
                        's2)'
     write (dpu,'(a)')  '-----------------------------------------'//&
                        '------'
   end if       
      
!-----------------------------------------------------------------------
!
!      Loop over convective layers to see if they need to be extended

   do while ( ncv .le. ncvfin )
       
     ncvinit = ncv
     cntu    = 0
     cntd    = 0
     kb      = kbase(ncv) 
     kt      = ktop(ncv)
     lbulk   = zi(kt)-zi(kb)
            
     if (column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a,i3)')  ' ncv    = ', ncv
       write (dpu,'(a,i3)')  ' kb     = ', kb
       write (dpu,'(a,i3)')  ' kt     = ', kt
       write (dpu,'(a,f14.7,a)')  ' zi(kt) = ', zi(kt), ' meters'
       write (dpu,'(a,f14.7,a)')  ' zi(kb) = ', zi(kb), ' meters'
       write (dpu,'(a,f14.7,a)')  ' lbulk  = ', lbulk, ' meters'
       write (dpu,'(a)')  ' '
     end if
    
!-----------------------------------------------------------
!            
! Add contribution (if any) from surface interfacial layer 
! to turbulent production and lengthscales.  If there is 
! positive surface buoyancy flux, the CL extends to the 
! surface and there is a surface interfacial layer contri-
! bution to W and the CL interior depth. If there is neg-
! ative buoyancy flux, the surface interfacial layer is 
! treated as energetically isolated from the rest of the CL
! and does not contribute to the layer-interior W or depth.
! This case also requires a redefinition of lbulk.

     bottom = kb .eq. kdim+1
     if (bottom .and. (bflxs .ge. 0.)) then
       lint = zm(kdim)
       dwsurf = (tkes/b1)*zm(kdim)
     else if (bottom .and. (bflxs .lt. 0.)) then
       lint = 0.
       dwsurf = 0.
       lbulk = zi(kt)-zm(kdim)
     else
       lint = 0.
       dwsurf = 0.
     end if
     l2n2 = 0.
     l2s2 = 0.
            
     if (column_match .and. bottom) then
       write(dpu,30) kdim+1,lint,lz,l2n2,l2s2
30   format(1X,i2,1X,2(f8.4,1X),2(f12.9,1X))
     end if
    
!-----------------------------------------------------------
!            
! Turbulence contribution from conv layer (CL) interior 
! kt < k < kb, which at this point contains only unstable 
! interfaces. Based on the CL interior stratification, 
! initial guesses at the stability functions are made. If 
! there is no CL interior interface, neutral stability is 
! assumed for now.

     if (kt .lt. kb-1) then 
            
       do k = kb-1, kt+1, -1
         lz   = lengthscale(zi(k),lbulk)
         l2n2 = l2n2 + lz*lz*n2(k)*(zm(k-1)-zm(k))
         l2s2 = l2s2 + lz*lz*s2(k)*(zm(k-1)-zm(k))
         lint = lint + (zm(k-1)-zm(k))      
         if (column_match) write(dpu,30) k,lint,lz,l2n2,  &
           l2s2
       enddo

!------------------------------------------------------
!
! Solve for bulk Sh, Sm, and wint over the CL interior
 
       ricl = min(l2n2/l2s2,ricrit) ! actually we should have 
                              ! ricl < 0 
       call galperin(ricl,gh,sh,sm)
       wint = -sh*l2n2 + sm*l2s2 + dwsurf 
       ebar = b1*wint/lint

     else

!-----------------------------------------------------------
!            
! There is no CL interior interface. The only way that 
! should happen at this point is if there is upward surface 
! buoy flux but no unstable interior interfaces. In that 
! case, the surface interface turbulent production terms are
! used as its CL 'interior'.

       if (bottom) then
         wint = dwsurf
         ebar = tkes     
      
!use neutral stability fns for layer extension
         call galperin(0.,gh,sh,sm) 
      
       else
         call error_mesg ('edt_mod', &
              'no convective layers found although ncv <= ncvfin', &
                           FATAL)
       endif
 
     endif
    
     if(column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a,f14.7)')  ' ricrit = ',ricrit
       write (dpu,'(a,f14.7)')  ' ricl   = ', ricl
       write (dpu,'(a,f14.7)')  ' gh     = ', gh
       write (dpu,'(a,f14.7)')  ' sh     = ', sh
       write (dpu,'(a,f14.7)')  ' sm     = ', sm
       write (dpu,'(a,f14.7,a)')  ' lbulk  = ', lbulk,  ' meters'
       write (dpu,'(a,f14.7,a)')  ' wint   = ', wint,   ' m3/s2'
       write (dpu,'(a,f14.7,a)')  ' dwsurf = ', dwsurf, ' m3/s2'
       write (dpu,'(a,f14.7,a)')  ' ebar   = ', ebar,   ' m2/s2' 
       write (dpu,'(a)')  ' '
     end if
    
!-----------------------------------------------------------
!            
! Try to extend the top of the convective layer. Compute 
! possible contributions to TKE production and lengthscale 
! were the CL top interfacial layer found by exacol incor-
! porated into the CL interior.

     extend = .false.    ! will become true if CL top is extended
     dzinc  = zm(kt-1)-zm(kt)
     lz     = lengthscale(zi(kt),lbulk)
     dl2n2  = lz*lz*n2(kt)*dzinc
     dl2s2  = lz*lz*s2(kt)*dzinc
     dwinc  = -sh*dl2n2 + sm*dl2s2

     if (column_match) then 
       write (dpu,'(a)')  '------------------------------------'//&
                          '-------'
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  '  trying to extend a layer upwards     '
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' k    dzinc      lz      dl2n2      '//&
                          '  dl2s2        dwinc' 
       write (dpu,'(a)')  '       (m)       (m)    (m3/s2)     '//&
                          ' (m3/s2)      (m3/s2)'
       write (dpu,'(a)')  '------------------------------------'//&
                          '----------------------'
       write (dpu,'(a)')  ' '
       write(dpu,31) kt,dzinc,lz,dl2n2,dl2s2,dwinc
31     format(1X,i2,1X,2(f8.4,1X),3(f12.9,1X))
     end if
    
!-----------------------------------------------------------
!            
! Test for incorporation of top layer kt into CL interior. 
! If true, extend CL by incorporating top layers until test 
! fails.
      
     l2n2 = -max(min(-l2n2,tkemax*lint/(b1*sh)),tkemin*lint/ (b1*sh))
     tmpr = -rinc*dzinc*l2n2/(lint+(1-rinc)*dzinc)
    
     if (column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a,f14.7,a)')  '-dl2n2 = ', -1.*dl2n2, ' m3/s2'
       write (dpu,'(a,f14.7,a)')  '  l2n2 = ',      l2n2, ' m3/s2'
       write (dpu,'(a,f14.7)')    '  rinc = ',      rinc
       write (dpu,'(a,f14.7,a)')  '  tmpr = ',      tmpr, ' m3/s2'
       !write (dpu,'(a,f14.7)')  ' will layer be extended (-dl2'//&
                        !    'n2 .gt. tmpr)? ', real( -dl2n2 .gt. tmpr )
       write (dpu,'(a)')  ' '
     end if
    
     do while (-dl2n2 .gt. tmpr)

!------------------------------------------------------
! Add contributions from layer kt to interior length-
! scale/TKE prod

       lint = lint + dzinc
       wint = wint + dwinc
       l2n2 = l2n2 + dl2n2
       l2n2 = -max(min(-l2n2,tkemax*lint/(b1*sh)),tkemin*lint/(b1*sh))
       l2s2 = l2s2 + dl2s2
       kt = kt-1
       extend = .true.
       if (kt .eq. 1) then
         call error_mesg ('edt_mod', &
              'trying to extend convective layer at model top',&
                           FATAL)
       end if

!------------------------------------------------------
! Check for existence of an overlying CL which might 
! be merged. If such exists (ktinc > 1), check for 
! merging by testing for incorporation of its top 
! interior interface into current CL. If no such layer 
! exists, ktop(ncv+cntu+1) will equal its default value
! of zero, so ktinc will be 1 and the test kt=ktinc
! will fail.

       ktinc = ktop(ncv+cntu+1)+1
       if (kt .eq. ktinc) then
         ncvfin = ncvfin - 1
         cntu   = cntu   + 1 
       end if

!------------------------------------------------------
! Compute possible lengthscale and TKE production
! contributions were layer kt incorporated into CL 
! interior. Then go back to top of loop to test for 
! incorporation.
            
       dzinc = zm(kt-1)-zm(kt)
       lz    = lengthscale(zi(kt),lbulk)
       dl2n2 = lz*lz*n2(kt)*dzinc
       dl2s2 = lz*lz*s2(kt)*dzinc
       dwinc = -sh*dl2n2 + sm*dl2s2
    
       if (column_match) write(dpu,31) kt,dzinc,lz,dl2n2, dl2s2,dwinc
    
 !------------------------------------------------------
 ! Recalculate tmpr
      
       tmpr  = -rinc*dzinc*l2n2/(lint+(1-rinc)*dzinc)
            
       if (column_match) then
         write (dpu,'(a)')  ' '
         write (dpu,'(a,f14.7,a)')  '-dl2n2 = ', -1.*dl2n2,    &
                                    ' m3/s2'
         write (dpu,'(a,f14.7,a)')  '  l2n2 = ',      l2n2,    &
                                    ' m3/s2'
         write (dpu,'(a,f14.7)')    '  rinc = ',      rinc
         write (dpu,'(a,f14.7,a)')  '  tmpr = ',      tmpr,    &
                                    ' m3/s2'
         !write (dpu,'(a,f14.7,a)')  ' will layer be extende'//&
                 !     'd (-dl2n2 .gt. tmpr)? ', real(-dl2n2 .gt. tmpr)
         write (dpu,'(a)')  ' '
       end if
    
     end do   ! Done with top extension of CL

     !-----------------------------------------------------------
     ! Shift indices appropriately if layers have been merged

     if (cntu .gt. 0) then
       do incv = 1, ncvfin - ncv
         kbase(ncv+incv) = kbase(ncv+cntu+incv)
         ktop(ncv+incv) = ktop(ncv+cntu+incv)
       end do
     end if

     !-----------------------------------------------------------
     !            
     ! Extend the CL base if possible.

     if (column_match .and. .not. bottom) then 
       write (dpu,'(a)')  '------------------------------------'//&
                          '-------'
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  '  trying to extend a layer downwards'
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  'k     dzinc       lz       dl2n2    '//&
                          '       dl2s2          dwinc' 
       write (dpu,'(a)')  '       (m)        (m)     (m3/s2)   '//&
                          '     (m3/s2)        (m3/s2)'
       write (dpu,'(a)')  '------------------------------------'//&
                          '---------------------------'
       write (dpu,'(a)')  ' '
     end if
    
     if (.not. bottom) then

!------------------------------------------------------
! Compute possible contributions to TKE production and 
! lengthscale, were the CL base interfacial layer found 
! by exacol incorporated into the CL interior.

       dzinc = zm(kb-1)-zm(kb)
       lz    = lengthscale(zi(kb),lbulk)
       dl2n2 = lz*lz*n2(kb)*dzinc
       dl2s2 = lz*lz*s2(kb)*dzinc
       dwinc = -sh*dl2n2 + sm*dl2s2

       if (column_match) write(dpu,31) kb,dzinc,lz,dl2n2,dl2s2,dwinc

!------------------------------------------------------
! Test for incorporation of base layer kb into CL 
! interior. If true, extend CL by incorporating base 
! layers until test fails.

       tmpr = -rinc*dzinc*l2n2/(lint+(1-rinc)*dzinc)
 
       if (column_match) then
         write (dpu,'(a)')  ' '
         write (dpu,'(a,f14.7,a)')  '-dl2n2 = ', -1.*dl2n2,    &
                                    ' m3/s2'
         write (dpu,'(a,f14.7,a)')  '  l2n2 = ',      l2n2,    &
                                    ' m3/s2'
         write (dpu,'(a,f14.7)')    '  rinc = ',      rinc
         write (dpu,'(a,f14.7,a)')  '  tmpr = ',      tmpr,    &
                                    ' m3/s2'
        !write (dpu,'(a,f14.7,a)')  ' will layer be extende'//&
                !     'd (-dl2n2 .gt. tmpr)? ', real(-dl2n2 .gt. tmpr)
         write (dpu,'(a)')  ' '
       end if

       do while((-dl2n2.gt. tmpr) .and. (.not. bottom) )

!-------------------------------------------------
! Add contributions from layer kb to interior 
! lengthscale/TKE prod
 
         lint = lint + dzinc
         wint = wint + dwinc
         l2n2 = l2n2 + dl2n2
         l2n2 = -max(min(-l2n2,tkemax*lint/(b1*sh)),tkemin&
                        *lint/(b1*sh))
         l2s2 = l2s2 + dl2s2

!-------------------------------------------------
! Extend base of CL downward a layer

         kb = kb+1
         extend = .true.

!-------------------------------------------------
! Check for existence of an underlying CL which 
! might be merged. If such exists (kbinc > 1), 
! check for merging by testing for incorporation 
! of its top interior interface into current CL.
! Note that this top 'interior' interface could be
! the surface.

         kbinc = 0
         if (ncv .gt. 1) kbinc = ktop(ncv-1)+1
         if (kb .eq. kbinc) then

!--------------------------------------------
! We are incorporating interior of CL ncv-1, 
! so merge this CL into the current CL.

           ncv    = ncv    - 1
           ncvfin = ncvfin - 1
           cntd   = cntd   + 1 
         end if

!-------------------------------------------------
! If CL would now reach the surface, check sign of
! surface buoyancy flux. If positive, add contri-
! butions of surface interfacial layer to TKE 
! production and lengthscale. If negative, we 
! regard the surface layer as stable and do not
! add surface interfacial layer contributions to 
! the CL. In either case the surface interface is 
! classified as part of the CL for bookkeeping 
! purposes (to ensure no base entrainment calcula-
! tion is done). If we are merging with a surface-
! driven CL with no interior unstable interfaces, 
! the above code will already have handled the 
! merging book-keeping.

         bottom = kb .eq. kdim+1
         if (bottom) then 
           if (bflxs .gt. 0.) then 
             dwsurf = (tkes/b1)*zm(kdim)
             lint = lint + zm(kdim)
           end if
         else

!--------------------------------------------
! Compute possible lengthscale and TKE prod-
! uction contributions were layer kb incor-
! porated into CL interior,then go back to 
! top of loop to test for incorporation
            
           dzinc = zm(kb-1) - zm(kb)
           lz    = lengthscale(zi(kb),lbulk)
           dl2n2 = lz*lz*n2(kb)*dzinc
           dl2s2 = lz*lz*s2(kb)*dzinc
           dwinc = -sh*dl2n2 + sm*dl2s2
                        
         end if

         if (column_match) write(dpu,31) kb,dzinc,lz,dl2n2,dl2s2,dwinc
    
!-------------------------------------------------
! Recalculate tmpr
      
         tmpr = -rinc*dzinc*l2n2/(lint+(1-rinc)*dzinc)

         if (column_match) then
           write (dpu,'(a)')  ' '
           write (dpu,'(a,f14.7,a)')  '-dl2n2 = ',-1.*dl2n2,' m3/s2'
           write (dpu,'(a,f14.7,a)')  '  l2n2 = ',     l2n2,' m3/s2'
           write (dpu,'(a,f14.7)')    '  rinc = ',     rinc
           write (dpu,'(a,f14.7,a)')  '  tmpr = ',     tmpr,' m3/s2'
           !write (dpu,'(a,f14.7,a)')  ' will layer be ex'//&
                     !     'tended (-dl2n2 .gt. tmpr)? ', real(-dl2n2 &
      !     .gt. tmpr)
           write (dpu,'(a)')  ' '
         end if
    
       end do ! for downward extension
 
       if (bottom .and. ncv .ne. 1) then 
         call error_mesg ('edt_mod', &
                          'bottom convective layer not indexed 1',&
                           FATAL)
       end if

     end if   ! Done with bottom extension of CL 

!-----------------------------------------------------------
! Shift indices if some layers with N2 < 0 have been found

     if (cntd .gt. 0) then
       do incv = 1, ncvfin - ncv
         kbase(ncv+incv) = kbase(ncvinit+incv)
         ktop(ncv+incv) = ktop(ncvinit+incv)
       end do
     end if

!-----------------------------------------------------------
! Sanity check for positive wint.
     if (wint .lt. 0.) then
       call error_mesg ('edt_mod', &
                        'interior avg TKE < 0', FATAL)
     end if

!-----------------------------------------------------------
! Recompute base and top indices, Ri_cl, Sh, Sm, and <W> 
! after layer extension if necessary. Ideally, we would 
! recompute l2n2 and l2s2 to account for the incorrect lbulk
! used in the computation of lz, but we take the simpler 
! approach of simply multiplying the lz's by the ratio of 
! the actual PBL depth to lbulk.
    
     if (extend) then

       ktop (ncv) = kt
       kbase(ncv) = kb
       zbot       = zi(kb)
       if (bottom .and. (bflxs.lt.0)) zbot = zm(kdim)
       l2rat      = ((zi(kt) - zbot)/lbulk)**2
       l2n2       = l2n2*l2rat
       l2s2       = l2s2*l2rat
       ricl = min(l2n2/l2s2,ricrit)
       call galperin(ricl,gh,sh,sm)
                 
 !------------------------------------------------------
 ! It is conceivable that even though the original wint 
 ! was positive, it will be negative after correction. 
 ! In this case, correct wint to be a small positive 
 ! number
       wint = max(dwsurf + (-sh*l2n2 + sm*l2s2),0.01*wint)

     end if  ! for extend if

     lbrk(ncv) = lint
     wbrk(ncv) = wint/lint
     ebrk(ncv) = b1*wbrk(ncv)
     ebrk(ncv) = max(min(ebrk(ncv),tkemax),tkemin)
     ghcl(ncv) = gh 
     shcl(ncv) = sh
     smcl(ncv) = sm
           
     if (column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a,i4)') ' FINAL RESULTS FOR CONVECTIVE LAYER', ncv
       write (dpu,'(a)')  ' '
       write (dpu,'(a,i4,i4)')   ' ktop(ncv), kbase(ncv) = ',     &
                                   ktop(ncv), kbase(ncv)
       write (dpu,'(a,f14.7,a)') ' lbrk(ncv) = ', lbrk(ncv),' m'
       write (dpu,'(a,f14.7,a)') ' wbrk(ncv) = ', wbrk(ncv),' m2/s2'
       write (dpu,'(a,f14.7,a)') ' sqrt(ebrk(ncv)) = ',           &
                                   sqrt(ebrk(ncv)), ' m/s'
       write (dpu,'(a,f14.7)')   ' ghcl(ncv) = ', ghcl(ncv)
       write (dpu,'(a,f14.7)')   ' shcl(ncv) = ', shcl(ncv)
       write (dpu,'(a,f14.7)')   ' smcl(ncv) = ', smcl(ncv)
     end if

            !-----------------------------------------------------------
    ! Increment counter for next CL

     ncv = ncv + 1

   end do     ! Loop over convective layers

       !----------------------------------------------------------------
       ! 
       ! set belongcv 

   belongcv(:) = .false.
   do ncv = 1, ncvfin
     do k = ktop(ncv), kbase(ncv)
       belongcv(k) = .true.
     enddo
   enddo

   if (column_match) then
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  '                  END OF ZISOCL '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  '-------------------------------------------'
   end if
              
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine zisocl

!
!=======================================================================

!======================================================================= 
!
!      subroutine caleddy
!        
!
!      this is the main program provided by Chris Bretherton and Herve
!      Grenier
!        

subroutine caleddy( u_star,       kqfs,         khfs,         sl,      &
                    qt,           ql,           qa,           slv,     &
                    u,            v,            pm,           zm,      &
                    sfuh,         sflh,         slslope,      qtslope, &
                    qrl,          hleff,        density,      qsl,     &
                    dqsldtl,      qltop,        pi,           zi,      &
                    chu,          chs,          cmu,          cms,     &
                    n2,           s2,           ri,           pblh,    &
                    turbtype,     kvh,          kvm,          tke,     &
                    leng,         bprod,        sprod,        trans,   &
                    diss,         isturb,       adj_time_inv, sdevs,   &
                    radfvec,      shvec,        smvec,        ghvec,   &
                    evhcvec)


!-----------------------------------------------------------------------
!
!      Driver routine to compute eddy diffusion coefficients for 
!      momentum, moisture, trace constituents and static energy.  Uses 
!      first order closure for stable turbulent layers. For convective 
!      layers, an entrainment closure is used, coupled to a diagnosis 
!      of layer-average TKE from the instantaneous thermodynamic and 
!      velocity profiles. Convective layers are diagnosed by extending 
!      layers of moist static instability into adjacent weakly stably 
!      stratified interfaces, stopping if the stability is too strong.  
!      This allows a realistic depiction of dry convective boundary 
!      layers with a downgradient approach.
! 
!      Authors:  Herve Grenier, 06/2000, Chris Bretherton 09/2000
!
!-----------------------------------------------------------------------
! 
! Variable declarations
!
! Inputs
!

real, intent(in)               :: u_star  ! surface friction velocity
real, intent(in)               :: kqfs    ! kinematic surf constituent
                                          ! flux (kg/kg)*(m/s)
real, intent(in)               :: khfs    ! kinematic surface heat flux
                                          ! (K*m/s)

!
!       the following fields are defined on model full
!       levels with dimension (1:nlev)
!

real, intent(in), dimension(:) :: sl      ! liq water static energy 
                                          !  cp*T+g*z-L*ql (J/kg)
real, intent(in), dimension(:) :: qt      ! total water spec. hum.
                                          ! (kg water/kg air) 
real, intent(in), dimension(:) :: ql      ! liq water spec. hum.
                                          ! (kg condensate/kg air)
real, intent(in), dimension(:) :: qa      ! cloud fraction (fraction)
real, intent(in), dimension(:) :: slv     ! liq water virtual static
                                          ! energy (J/kg)
                                          !  sl*(1 + .608*qt)
real, intent(in), dimension(:) :: u       ! u wind input (m/s)
real, intent(in), dimension(:) :: v       ! v wind input (m/s)
real, intent(in), dimension(:) :: pm      ! midpoint pressures (Pa)
real, intent(in), dimension(:) :: zm      ! layer midpoint height 
                                          ! above sfc (m)
real, intent(in), dimension(:) :: sfuh    ! sat frac in upper half-lyr
real, intent(in), dimension(:) :: sflh    ! sat frac in lower half-lyr    
real, intent(in), dimension(:) :: slslope ! sl slope with respect to
                                          ! pressure in thermo lyr
  ! (J/kg/Pa)
real, intent(in), dimension(:) :: qtslope ! qt slope with respect to
                                          ! pressure in thermo lyr
                                          ! (kg water/kg air/Pa)
real, intent(in), dimension(:) :: qrl     ! LW heating rate (K/s)
real, intent(in), dimension(:) :: hleff   ! effective latent heat of 
                                          ! condensation (J/kg)
real, intent(in), dimension(:) :: density ! air density (kg/m3)
real, intent(in), dimension(:) :: qsl     ! saturation specific hum.
                                          ! (kg water/kg air)
real, intent(in), dimension(:) :: dqsldtl ! temperature derivative of
                                          ! qsl (kg water/kg air/K)
real, intent(in), dimension(:) :: qltop   ! cloud liquid at the top
                                          ! of the thermo layer
  
!
!       the following fields are defined on model half levels,
!       dimension(1:nlev+1)
!

real, intent(in), dimension(:) :: pi      ! interface pressures (Pa)
real, intent(in), dimension(:) :: zi      ! interface height above
                                          ! sfc (m)
real, intent(in), dimension(:) :: chu     ! Unsat sl (heat) coef (1/m)
real, intent(in), dimension(:) :: chs     ! Sat sl (heat) coef (1/m)
real, intent(in), dimension(:) :: cmu     ! Unsat qt (moisture) coef 
                                          ! (kg/kg)*(m/s*s)
real, intent(in), dimension(:) :: cms     ! Sat qt (moisture) coef 
                                          ! (kg/kg)*(m/s*s)
real, intent(in), dimension(:) :: n2      ! Moist squared buoy freq 
                                          ! (1/sec**2)
real, intent(in), dimension(:) :: s2      ! Squared deformation (s-2)
                                          ! (1/sec**2)
real, intent(in), dimension(:) :: ri      ! gradient Richardson number

!
! Outputs
!

real,    intent(out)   :: pblh     ! planetary boundary layer height (m)


!
!       the following field is defined on model full levels,
!       dimension(1:nlev)
!

real,    intent(out), dimension(:) :: isturb ! is full layer part of a 
     ! turbulent layer?

!
!       the following fields are defined on model half levels,
!       dimension(1:nlev+1)
!

integer, intent(out), dimension(:,:) :: turbtype! Interface turb. type
                                                ! 1 = stable turb 
          ! 2 = CL interior
                                                ! 3 = bottom entr intfc 
                                                ! 4 = upper entr intfc 
real,    intent(out), dimension(:) :: kvh       ! diffusivity for heat 
                                                ! and tracers (m*m/s)
real,    intent(out), dimension(:) :: kvm       ! diffusivity for mom.
                                                ! (m*m/s)
real,    intent(out), dimension(:) :: tke       ! turb. kin. energy 
                                                ! (m*m)/(s*s)
real,    intent(out), dimension(:) :: leng      ! turbulent length scale
                                                ! (m)
real,    intent(out), dimension(:) :: bprod     ! Buoyancy production
                                                ! (m*m)/(s*s*s)
real,    intent(out), dimension(:) :: sprod     ! shear production
                                                ! (m*m)/(s*s*s)
real,    intent(out), dimension(:) :: trans     ! TKE transport
                                                ! (m*m)/(s*s*s)
real,    intent(out), dimension(:) :: diss      ! TKE dissipation
                                                ! (m*m)/(s*s*s)
real,    intent(out), dimension(:) :: adj_time_inv  ! inverse adjustment 
                                                    ! time for turbulent
    ! layer (1/sec)
real,    intent(out), dimension(:) :: sdevs     ! std. dev. of water
                                                ! perturbation
        ! (kg water/kg air)
        ! defined on model half
        ! levels
real,    intent(out), dimension(:) :: radfvec   ! Buoyancy production
                                                ! from lw radiation
                                                ! (m*m)/(s*s*s)
real,    intent(out), dimension(:) :: shvec     ! Galperin heat stab. 
                                                ! fn. (none)
real,    intent(out), dimension(:) :: smvec     ! Galperin mom. stab. 
                                                ! fn. (none)
real,    intent(out), dimension(:) :: ghvec     ! Galperin stability
                                                ! ratio (none)
real,    intent(out), dimension(:) :: evhcvec   ! Evaporative cooling
                                                ! entrainment factor
! (none)


!
! Internal variables
!

logical :: in_CL                             ! True if interfaces k,k+1
                                             ! both in same CL
logical :: any_stable                        ! Are there any stable 
                                             ! turbulent layers?
logical :: cloudtop                          ! Is the interface at 
                                             ! cloudtop?
logical, dimension(size(sl,1)+1) :: belongcv ! True for interfaces 
                                             ! interior to convective
     ! layer (CL)
logical, dimension(size(sl,1)+1) :: belongst ! True for interfaces 
                                             ! interior to a stable
     ! turbulent layer
integer :: k                    ! vertical index
integer :: ks                   ! vertical index
integer :: kk                   ! vertical index
integer :: kdim                 ! number of full vertical levels
integer :: ncvfin               ! Total number of CL in column
integer :: ncvf                 ! Total number of CL in column prior to
                                ! addition of one layer rad-driven CLs
integer :: ncv                  ! index of current CL
integer :: ncvnew               ! index of added one layer rad-driven CL
integer :: ncvsurf              ! if nonzero, index of CL including 
                                ! surface
integer :: kb, kt               ! kbase and ktop for current CL


integer, dimension(size(sl,1)+1) :: kbase     ! vert. index for base
                                              ! interface of CL
integer, dimension(size(sl,1)+1) :: ktop      ! vert. index for top 
                                              ! interface of CL

real    :: bflxs    ! Surface buoyancy flux (m2/s3)
real    :: tkes     ! Surface TKE
real    :: jtzm     ! Interface layer thickness atop conv layer (CL) ncv
real    :: jtsl     ! Jump in s_l               atop       "
real    :: jtqt     ! Jump in q_t               atop       "
real    :: jtbu     ! Jump in buoyancy          atop       "
real    :: jtu      ! Jump in zonal wind        atop       "
real    :: jtv      ! Jump in meridional wind   atop       "
real    :: jt2slv   ! 2-layer Jump in s_lv              atop       "
real    :: radf     ! buoy flx jump at cloudtop from lw rad flx div
real    :: jbzm     ! Interface layer thickness at base of CL ncv
real    :: jbsl     ! Jump in s_l               at base    "
real    :: jbqt     ! Jump in qt                at base    "
real    :: jbbu     ! Jump in buoyancy          at base    "
real    :: jbu      ! Jump in zonal wind        at base    "
real    :: jbv      ! Jump in merid. wind       at base    "
real    :: ch       ! buoy flux coefs for sl, qt in a half-layer 
real    :: cm       ! 
real    :: n2h      ! Moist squared buoy freq for a half-layer (s-2)
real    :: ckh      ! Galperin stability function for heat
real    :: ckm      ! Galperin stability function for momentum
real    :: gh       ! Normalised buoyancy production (m*m*m)/(s*s)
real    :: lbulk    ! Depth of turbulent layer (m)
real    :: trma     ! intermediate variables
real    :: vus
real    :: vub
real    :: trmp
real    :: trmq
real    :: angle
real    :: qq
real    :: rootp             
real    :: evhc          ! (1+E) with E = evap. cool. efficiency [nd]
real    :: vys           ! n2h/n2 at upper inversion [nd]
real    :: vyb           ! Same at lower inversion [nd]
real    :: kentr         ! effective entrainment diffusivity we*dz (m*m/s)
real    :: lwp           ! liquid water path in layer kt (kg cond/m2)
real    :: opt_depth     ! optical depth of layer kt
real    :: radinvfrac    ! frac of lw cooling in layer kt put at inv.
real    :: qtsltmp       ! dqt/dz (kg/kg/m)
real    :: slsltmp       ! dsl/dz (J/kg/m)
real    :: qsltmp        ! qsl (kg water/kg air)
real    :: dqsldtltmp    ! dqsldtl (kg water/kg air/K)
real    :: hlefftmp      ! effective latent heat (J/kg water)
real    :: qstartmp      ! q* at upper/lower inversion (kg water/kg air)
real    :: bstartmp      ! b* at upper/lower inversion (m/s2)
real    :: temperature   ! actual temperature (K)
real    :: tmpfrac       ! fraction of cloud at top of convective layer
                         ! exposed to clear air above using maximum 
 ! overlap assumption
real, dimension(size(sl,1)+1) :: ebrk,wbrk,lbrk,ghcl,shcl,smcl
real, dimension(size(sl,1)+1) :: wcap          ! W (m2/s2)
real, dimension(size(sl,1)+1) :: rcap          ! e/<e> 

!-----------------------------------------------------------------------
!
!  Initialize to zero outputs needed at all interfaces, but
!  calculated only at turbulent interfaces. Only kvh and kvm are 
!  outputs, the other arrays are zeroed for plotting or diagnostic 
!  purposes. 

   kdim             = size(sl,1)
   kvh(:)           = 0.0
   kvm(:)           = 0.0
   wcap(:)          = 0.0
   leng(:)          = 0.0
   rcap(:)          = 0.0
   bprod(:)         = 0.0
   sprod(:)         = 0.0
   trans(:)         = 0.0
   diss(:)          = 0.0
   tke(:)           = 0.0
   turbtype(:,:)    = 0
   adj_time_inv(:)  = missing_value
   isturb(:)        = 0.0
   sdevs(:)         = missing_value
   radfvec(:)       = 0.0
   shvec(:)         = missing_value
   smvec(:)         = missing_value
   ghvec(:)         = missing_value
   evhcvec(:)       = missing_value
       
!-----------------------------------------------------------------------
!
!      Optional printout
        
   if (column_match) then
     write (dpu,'(a)')  '--------------------------------------------'
     write (dpu,'(a)')  '        ENTERING CALEDDY                '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  'checking inputs: '
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7,a)') ' u_star             = ',u_star, ' m/s'
     write (dpu,'(a,f14.7,a)') ' latent   heat flux = ',density(kdim)&
                               *hleff(kdim)*kqfs,' W/m2'
     write (dpu,'(a,f14.7,a)') ' sensible heat flux = ',density(kdim)&
                               *cp_air*khfs,' W/m2'
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k     sl/cp_air        qt        qc        q'//&
                        'a     slv/cp_air' 
     write (dpu,'(a)')  '        (K)       (g/kg)    (g/kg)       '//&
                        '       (K)  '
     write (dpu,'(a)')  '-----------------------------------------'//&
                        '------------'
     write (dpu,'(a)')  ' '
     do kk = kdim-n_print_levels,kdim
       write(dpu,224) kk,sl (kk)/cp_air,1000.*qt(kk),1000.*ql(kk),    &
                      qa(kk),slv(kk)/cp_air
     end do
224  format(1X,i2,1X,5(f9.4,1X))
24   format(1X,i2,1X,4(f9.4,1X))
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k      u          v      z_full     p_full'
     write (dpu,'(a)')  '      (m/s)      (m/s)     (m)        (mb) '
     write (dpu,'(a)')  '-------------------------------------------'
     write (dpu,'(a)')  ' '
     do kk = kdim-n_print_levels,kdim
       write(dpu,24) kk,u(kk),v(kk),zm(kk),pm(kk)/100.
     enddo    
     write (dpu,'(a)') ' '
     write (dpu,'(a)') ' '
     write (dpu,'(a)') ' k     sfuh      sflh     slslope   qtslope'
     write (dpu,'(a)') '                        (K/100mb)(g/kg/100mb)'              
     write (dpu,'(a)') '---------------------------------------------'
     write (dpu,'(a)') ' '
     do kk = kdim-n_print_levels,kdim
       write(dpu,24) kk,sfuh(kk),sflh(kk),slslope(kk)*100.*100./cp_air&
                     ,1000.*qtslope(kk)*100.*100.
     enddo
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  'k       rho       qrl      hleff'
     write (dpu,'(a)')  '      (kg/m3)   (K/day)   (MJ/kg)'
     write (dpu,'(a)')  '---------------------------------'
     write (dpu,'(a)')  ' '
     do kk = kdim-n_print_levels,kdim
       write(dpu,124) kk,density(kk),qrl(kk)*86400,hleff(kk)/1.e+06
124  format(1X,i2,1X,3(f9.4,1X))
     enddo
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k   p_half  z_half  chu    chs     cmu  '//&
                        '   cms   sqrt(n2)  sqrt(s2)     ri'
     write (dpu,'(a)')  '     (mb)    (m)    (m/K/s**2)      (m/s*'//&
                        '*2)      (1/hr)    (1/hr) '
     write (dpu,'(a)')  '-----------------------------------------'//&
                        '--------------------------------------'
     write (dpu,'(a)')  ' '
     do kk = kdim+1-n_print_levels,kdim+1
       write(dpu,26) kk,pi(kk)/100.,zi(kk),cp_air*chu(kk),cp_air*chs(kk), &
                     cmu(kk),cms(kk),n2(kk)*sqrt(abs(n2(kk)))*3600./&
                     max(small,abs(n2(kk))),sqrt(s2(kk))*3600.,ri(kk)
     enddo
26   format(1X,i2,1X,f7.2,1X,f7.1,1X,2(f7.4,1X),2(f6.3,1X),3(f9.3,1X))
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  '-----------------------------------------'//&
                        '-------------------------------------'
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
   end if
       
!-----------------------------------------------------------------------
!
!  calculate 'surface' (actually lowest half-layer) buoyancy flux.
!
!  Note that sensible heat flux in W/m2 = density_air * cp * khfs
!              latent heat flux in W/m2 = density_air * L  * kqfs
!
!  units of bflxs below is m2/s3
!
!  Buoyancy flux in W/m2 = rho * bflxs / ch
!
   ch = chu(kdim+1)*(1-sflh(kdim)) + chs(kdim+1)*sflh(kdim)   
   cm = cmu(kdim+1)*(1-sflh(kdim)) + cms(kdim+1)*sflh(kdim)   
   bflxs  = ch*cp_air*khfs + cm*kqfs
   bprod(kdim+1) = bflxs       
             
!-----------------------------------------------------------------------
!
!  Diagnostic surface TKE used in calculating layer-average TKE.
!
!  Chris B. note:
!
!  Originally tkes was computed from the commented out line below.
!  This commented-out line can make tkes small or even negative. 
!  The replacement has at least equal justifiability, trying to 
!  represent tke at the surface rather than at zm(kdim).  Near the
!  surface is where averaging W to get <e> is not really just-
!  ifiable.
!
!  tkes = (b1*(ustar**3+vonkarm*zm(kdim)*bflxs))**2./3.
!  
!  [Steve Klein note:
!
!  The replacement line appears to create an inconsistency in the 
!  sense that the buoyancy flux in the lower half of level one does
!  not effectively contribute to dwsurf, the contribution to W,
!  in the lower half of level one. Hence diagnosing shear production
!  and dissipation at the surface is a bit of a ruse.]

   tkes          = max(min(b123*u_star**2,tkemax),tkemin)
   tke(kdim+1)   = tkes
   diss(kdim+1)  = tkes**(3/2)/b1/zm(kdim)              
   sprod(kdim+1) = tkes**(3/2)/b1/zm(kdim)

   if (column_match) then       
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7,a)')  ' ch (surface)  = ',cp_air*ch,' m/K/s**2'
     write (dpu,'(a,f14.7,a)')  ' cm (surface)  = ',cm,' m/s**2'
     write (dpu,'(a,f14.7,a)')  ' bflxs         = ',bflxs,' m2/s3'
     write (dpu,'(a,f14.7,a)')  ' buoyancy flux = ',density(kdim)*   &
                                bflxs/ch,' W/m2'
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7,a)')  ' u_star           = ', u_star,' m/s'
     write (dpu,'(a,f14.7,a)')  ' surface tke sqrt = ',sqrt(tkes),   &
                                ' m/s'
     write (dpu,'(a)')  ' '
   end if        

!-----------------------------------------------------------------------
!
!  Examine each column and determine whether it is convective.

    call exacol(ri, bflxs, ktop, kbase, ncvfin)

!-----------------------------------------------------------------------
!
!  CONVECTIVE LAYER (CL) computations
!
!  If some convective layers have been found, determine their bulk 
!  properties (<e>, <Sh>, <Sm>, and indices ktop/bkase for upper/
!  lower inversion ).

   ncvsurf = 0  
   
   if (ncvfin .gt. 0) then
   
      call zisocl(u_star, bflxs, tkes, zm, ql, zi, n2, s2, ri,     &
                  ncvfin, kbase, ktop, belongcv, ebrk, wbrk, ghcl, &
                  shcl, smcl, lbrk)
 
!-------------------------------------------------------------
!  CLs found by zisocl are in order of height, so if any CL 
!  contains the surface, it will be CL1.
         
     if (kbase(1) .eq. kdim+1) ncvsurf = 1
  
   else
       
     belongcv(:) = .false.
       
   end if

   if (column_match) write (dpu,'(a,i3)')  ' ncvsurf = ', ncvsurf

!-----------------------------------------------------------------------
!
!  Find single-level radiatively-driven cloud-topped convective 
!  layers (SRCLs). SRCLs extend through a single thermo layer k, 
!  with entrainment at interfaces k and k+1 (unless k+1 is the 
!  surface, in which case surface shear generation contributes to 
!  the layer-averaged energy). The conditions for an SRCL are:
!      1. cloud at level k
!      2. no cloud at level k+1 (else assuming that some fraction 
!         of the longwave flux div in layer k is concentrated at 
!         the top interface is invalid.
!      3. Longwave radiative cooling (shortwave heating is assumed
!         uniformly distributed through layer k, so not relevant to
!         buoyancy production of TKE)
!      4. Internal stratification n2h of half-layer from level k to
!         interface k is unstable using similar method as in sfdiag,
!         but applied to internal slopes of sl, qt in layer k.
!      5. Interfaces k, k+1 not both in the same existing convective
!         layer.
!      6. k >= 2 
!      7. Ri at interface k > ricrit, otherwise stable turb mixing 
!         will broadly distribute the cloud top in the vertical, 
!         preventing localized radiative radiative destabilization 
!         at the interface height.

   ncv = 1
   ncvf = ncvfin
   
   do k = kdim, 2, -1
             
     cloudtop  = .false. 
     if (ql(k)  .gt. qcminfrac*qsl(k)   .and. use_qcmin .and.   &
         ql(k-1).lt. qcminfrac*qsl(k-1))        cloudtop = .true.
     if (.not.use_qcmin .and. qa(k).gt.qa(k-1)) cloudtop = .true.
     
     if (qrl(k).lt. 0. .and. ri(k).gt.ricrit .and. cloudtop) then

       ch  = (1 -sfuh(k))*chu(k) + sfuh(k)*chs(k)
       cm  = (1 -sfuh(k))*cmu(k) + sfuh(k)*cms(k)
       n2h = ch*slslope(k) + cm*qtslope(k)

       if (n2h.le.0.) then

       !-------------------------------------------------
       ! Test if k and k+1 are part of the same preexist-
       ! ing CL. If not, find appropriate index for new 
       ! SRCL. Note that this calculation makes use of 
       ! ncv set from prior passes through the k do loop
       
         in_CL = .false.

         do while (ncv .le. ncvf)

           if (ktop(ncv) .le. k) then

           !---------------------------------------
           ! If kbase > k, k and k+1 are part of 
           ! same prior CL
             if (kbase(ncv) .gt. k) in_CL = .true.

             !---------------------------------------
             ! exit from do-loop once CL top at/above
             ! intfc k.
             exit  

           else

           !---------------------------------------
           !  Go up one CL
             ncv = ncv + 1  

           end if
         end do ! ncv

         !-------------------------------------------------
         ! Add a new SRCL

         if (.not.in_CL) then
           ncvfin        = ncvfin+1
           ncvnew        = ncvfin
           ktop(ncvnew)  = k
           kbase(ncvnew) = k+1
           belongcv(k)   = .true.
           belongcv(k+1) = .true.

           if (k.lt.kdim) then
             ebrk(ncvnew) = 0.
             lbrk(ncvnew) = 0.
             shcl(ncvnew) = 0.
             smcl(ncvnew) = 0.

           else 

           !---------------------------------------
           ! surface radiatively driven fog
             if (bflxs.gt.0.) then 
               !----------------------------------
               ! unstable surface layer 
               ! incorporate surface TKE into
               ebrk(ncvnew) = tkes
               lbrk(ncvnew) = zm(k)
               adj_time_inv(kdim+1) = sqrt(tkes)/zm(k)
             else   
               !----------------------------------
               ! stable surface layer 
               ! don't incorporate surface TKE 
               ebrk(ncvnew) = 0.
               lbrk(ncvnew) = 0.
             end if
             shcl(ncvnew) = 0.
             smcl(ncvnew) = 0.
             ncvsurf = ncvnew

           end if    ! k < kdim  

         end if    ! new SRCL

       end if    ! n2h < 0 

     end if ! qrl < 0, ri(k) > ricrit and cloudtop

   end do ! k do loop, end of SRCL section

!-----------------------------------------------------------------------
!
!      For each CL, compute length scale, r^2, e, Kh and Km

   if (column_match) then
     write (dpu,'(a)') ' '
     write (dpu,'(a)') ' '
     write (dpu,'(a)') ' IN CALEDDY CALCULATION OF L, R2, e, Kh, a'//&
                       'nd Km'
     write (dpu,'(a)') ' '       
   end if

   do ncv = 1, ncvfin

     kt    = ktop(ncv)
     kb    = kbase(ncv)
     lbulk = zi(kt)-zi(kb)          

     if (column_match) then
       write (dpu,'(a,i3)')  ' convective layer #',ncv
       write (dpu,'(a,i4,i4)')  ' kt, kb = ', kt,kb
       write (dpu,'(a,f14.7,a)')  ' lbulk = ', lbulk, ' m'
       write (dpu,'(a)')  ' '
     end if
    
     do k = min(kb,kdim), kt, -1             
       leng(k) = lengthscale(zi(k),lbulk)
       wcap(k) = (leng(k)**2)*(-shcl(ncv)*n2(k)+ &
                  smcl(ncv)*s2(k))
     end do    ! k do loop

     if (column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' k     leng          wcap '
       write (dpu,'(a)')  '       (m)          (m2/s2)'
       write (dpu,'(a)')  '----------------------------'
       write (dpu,'(a)')  ' '
       do k = min(kb,kdim),kt,-1
         write(dpu,33) k,leng(k),wcap(k)
       enddo
33     format(1X,i2,1X,f9.4,1X,f14.9)    
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' '
     end if
    
     !-----------------------------------------------------------
     ! Calculate jumps at the lower inversion 

     if (kb .lt. kdim+1) then 

       jbzm = zm(kb-1) - zm(kb)
       jbsl = sl(kb-1) - sl(kb)
       jbqt = qt(kb-1) - qt(kb)
       jbbu = n2(kb)   * jbzm
       jbbu = max(jbbu,jbumin)
       jbu  = u(kb-1)  - u(kb)
       jbv  = v(kb-1)  - v(kb)
       ch   = (1 -sflh(kb-1))*chu(kb) + sflh(kb-1)*chs(kb)
       cm   = (1 -sflh(kb-1))*cmu(kb) + sflh(kb-1)*cms(kb)
       n2h  = (ch*jbsl + cm*jbqt)/jbzm
       vyb  = n2h*jbzm/jbbu
       vub  = min(1.,(jbu**2+jbv**2)/(jbbu*jbzm) )

     else  
 !------------------------------------------------------
 ! Zero bottom entrainment contribution for CL extending
 ! down to sfc
            
       vyb = 0.
       vub = 0.
          
     end if
            
     if (column_match .and. kb .lt. kdim+1) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' jumps at lower inversion '
       write (dpu,'(a)')  ' ------------------------ '
       write (dpu,'(a)')  ' '
       write (dpu,'(a,i3)')  ' inversion half level = ',kb
       write (dpu,'(a,f14.7,a)')  ' jbzm = ',jbzm,' m'
       write (dpu,'(a,f14.7,a)')  ' jbsl = ',jbsl/cp_air,' K'
       write (dpu,'(a,f14.7,a)')  ' jbqt = ',jbqt*1000.,' g/kg'
       write (dpu,'(a,f14.7,a)')  ' jbbu = ',jbbu,' m/s2'
       write (dpu,'(a,f14.7,a)')  ' jbumin = ',jbumin,' m/s2'
       write (dpu,'(a,f14.7,a)')  ' jbu = ',jbu,' m/s'
       write (dpu,'(a,f14.7,a)')  ' jbv = ',jbv,' m/s'
       write (dpu,'(a,f14.7,a)')  ' ch  = ',ch*cp_air, ' m/K/s2'
       write (dpu,'(a,f14.7,a)')  ' cm  = ',cm, ' m/s2'
       write (dpu,'(a,f14.7,a)')  ' sqrt(n2h) = ',n2h*            &
                                    sqrt(abs(n2h))*3600./abs(n2h),' 1/sec'
       write (dpu,'(a,f14.7,a)')  ' vyb = ', vyb
       write (dpu,'(a,f14.7,a)')  ' vub = ', vub
       write (dpu,'(a)')  ' '
     end if
    
     !-----------------------------------------------------------
     ! Calculate jumps at the upper inversion
     ! Note the check to force jtbu to be greater than or equal
     ! to jbumin.

     jtzm = zm(kt-1) - zm(kt)
     jtsl = sl(kt-1) - sl(kt)
     jtqt = qt(kt-1) - qt(kt)
     jtbu = n2(kt)   * jtzm 
     jtbu = max(jtbu,jbumin)
     jtu  = u(kt-1) - u(kt)
     jtv  = v(kt-1) - v(kt)
     ch   = (1 -sfuh(kt))*chu(kt) + sfuh(kt)*chs(kt)
     cm   = (1 -sfuh(kt))*cmu(kt) + sfuh(kt)*cms(kt)
     n2h  = (ch*jtsl + cm*jtqt)/jtzm
            
     ! Ratio of buoy flux to w'(b_l)'
     vys  = n2h*jtzm/jtbu 
            
     ! Inverse of shear prodution divided by buoyancy production
     vus  = min(1.,(jtu**2+jtv**2)/(jtbu*jtzm)) 
            
     if (column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' jumps at upper inversion '
       write (dpu,'(a)')  ' ------------------------ '
       write (dpu,'(a)')  ' '
       write (dpu,'(a,i3)')  ' inversion half level = ',kt
       write (dpu,'(a,f14.7,a)')  ' jtzm = ',jtzm,' m'
       write (dpu,'(a,f14.7,a)')  ' jtsl = ',jtsl/cp_air,' K'
       write (dpu,'(a,f14.7,a)')  ' jtqt = ',jtqt*1000.,' g/kg'
       write (dpu,'(a,f14.7,a)')  ' jtbu = ',jtbu,' m/s2'
       write (dpu,'(a,f14.7,a)')  ' jbumin = ',jbumin,' m/s2'
       write (dpu,'(a,f14.7,a)')  ' jtu = ',jtu,' m/s'
       write (dpu,'(a,f14.7,a)')  ' jtv = ',jtv,' m/s'
       write (dpu,'(a,f14.7,a)')  ' ch  = ',ch*cp_air, ' m/K/s2'
       write (dpu,'(a,f14.7,a)')  ' cm  = ',cm, ' m/s2'
       write (dpu,'(a,f14.7,a)')  ' sqrt(n2h) = ',n2h*            &
                                    sqrt(abs(n2h))*3600./abs(n2h),' 1/sec'
       write (dpu,'(a,f14.7,a)')  ' vys = ', vys
       write (dpu,'(a,f14.7,a)')  ' vus = ', vus
       write (dpu,'(a)')  ' '
     end if
        
     !-----------------------------------------------------------
     ! 
     ! Calculate evaporative entrainment enhancement factor evhc. 
     ! We take the full inversion strength to be jt2slv, where
     ! jt2slv = slv(kt-2)  - slv(kt), and kt - 1 is in the 
     ! ambiguous layer.  However, for a cloud-topped CL overlain
     ! by another convective layer, it is possible that 
     ! slv(kt-2) < slv(kt). To avoid negative or excessive evhc, 
     ! we lower-bound jt2slv and upper-bound evhc.

     evhc = 1.
     cloudtop  = .false. 
     if (ql(kt)  .gt. qcminfrac*qsl(kt)   .and. use_qcmin .and. &
         ql(kt-1).lt. qcminfrac*qsl(kt-1))       cloudtop= .true.
     if (.not.use_qcmin .and. qa(kt).gt.qa(kt-1))cloudtop= .true.
    
     if (cloudtop) then 
       if (use_qcmin) then
         tmpfrac = 1.
       else
         tmpfrac = 1. - (qa(kt-1)/qa(kt))
       end if     
       jt2slv = slv(max(kt-2,1)) - slv(kt)
       jt2slv = max(jt2slv, jbumin*slv(kt-1)/grav)
       if (use_extrapolated_ql) then
         evhc = 1.+tmpfrac*a2l*a3l*hleff(kt)*qltop(kt)/ &
                jt2slv
       else
         evhc = 1.+tmpfrac*a2l*a3l*hleff(kt)*ql(kt)   / &
                jt2slv
       end if
       evhc   = min(evhc,evhcmax)
       evhcvec(kt) = evhc
     end if

     if (column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' computing entrainment enhancem'//&
                          'ent factor '
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' ' 
       write (dpu,'(a,f14.7,a)')  ' slv(max(kt-2,1))/cp_air =  ',&
                                    slv(max(kt-2,1))/cp_air,' K'
       write (dpu,'(a,f14.7,a)')  ' slv(kt)/cp_air =  ', slv(kt)/&
                                    cp_air,' K'
       write (dpu,'(a,f14.7,a)')  ' jt2slv(kt)/cp_air =  ',jt2slv&
                                   /cp_air,' K'
       write (dpu,'(a,f14.7,a)')  ' ql(kt-1) = ',ql(kt-1)*   &
                                    1000.,' g/kg'
       write (dpu,'(a,f14.7,a)')  ' ql(kt) = ',ql(kt)*1000., &
                                  ' g/kg'
       write (dpu,'(a,f14.7,a)')  ' qltop(kt) = ',qltop(kt)* &
                                    1000.,' g/kg'
       write (dpu,'(a,f14.7,a)')  ' qa(kt-1) = ',qa(kt-1)
       write (dpu,'(a,f14.7,a)')  ' qa(kt) = ',  qa(kt)                    
       write (dpu,'(a,f14.7,a)')  ' evhc = ',evhc
     end if
    
     !-----------------------------------------------------------
     ! 
     ! Radiative forcing at the upper inversion if at a cloud top

     if (cloudtop) then
     
     !------------------------------------------------------
     !  estimate longwave opt depth in layer kt

       lwp = ql(kt) * (pi(kt+1) - pi(kt)) / grav
       opt_depth = 156*lwp

     !------------------------------------------------------
     ! Approximation to LW cooling frac at inversion
     ! The following formula is a polynomial approximation
     ! to exact solution which is
     ! radinvfrac = 1 - 2/opt_depth + 2/(exp(opt_depth)-1))
 
       radinvfrac  = opt_depth*(4.+opt_depth) / &
                     (6.*(4.+opt_depth) + opt_depth**2)

     !------------------------------------------------------
     ! units of radf = (m*m)/(s*s*s)
       radf = max(-radinvfrac*qrl(kt)*(zi(kt)-zi(kt+1)),0.)* &
              cp_air * chs(kt)

     else
       lwp        = 0.
       opt_depth  = 0.
       radinvfrac = 0.
       radf       = 0.
     end if

     !-----------------------------------------------------------
     ! 
     ! Solve cubic equation
     !   r^3 + trmp*r + trmq = 0,   r = sqrt<e>
     ! to estimate <e> for multilayer convection. Note, that if 
     ! the CL goes to the surface, vyb and vub are zero, and ebrk
     ! and lbrk have already incorporated the surface interfacial
     ! layer, so the formulas below still apply.  For a SRCL, 
     ! there are no interior interfaces so ebrk = lbrk = 0.
     ! The cases are:
     !    (1) no cloudtop cooling (radf=0) -- trmq = 0, 
     !        r = sqrt(-trmp)
     !    (2) radf > 0 but no interior CL interface -- trmp = 0, 
     !        trmq < 0
     !    (3) radf > 0 and interior CL interface(s) -- trmp < 0, 
     !        trmq < 0

     trma = 1. - (b1*a1l/lbulk)*(evhc*(-vys+vus)*(zi(kt)-zm(kt))& 
            +  (-vyb+vub)*(zm(kb-1)-zi(kb)) )
     trma = max(trma,0.5)  ! Prevents runaway entrainment instab.
     trmp = -ebrk(ncv) *(lbrk(ncv)/lbulk)/trma
     trmq = -b1*radf*(leng(kt)/lbulk)*(zi(kt)-zm(kt))/trma

     qq = (trmp/3.)**3+(trmq/2.)**2
     if (trmq .lt. 0.) then
       if (qq .gt. 0.) then 
         rootp = (-trmq/2.+sqrt(qq))**(1./3.)
         if (trmp .lt. 0.) then       
         !--------------------------------------------
         ! case 3 (in case 2, added term is zero)
           rootp = rootp  + (-trmq/2.-sqrt(qq))**(1./3.)
         end if
       else  
       !-------------------------------------------------
       ! also part of case 3
         angle = acos(-trmq/2./sqrt(-(trmp/3)**3))
         rootp = 2.*sqrt(-trmp/3.)*cos(angle/3.)
       end if
     else
     !------------------------------------------------------
     !  case 1: radf = 0, so trmq = 0
       rootp = sqrt(-trmp)
     endif
          
     if (column_match) then
       write (dpu,'(a,f14.7)')  ' trma = ',trma
       write (dpu,'(a,f14.7,a)')  ' trmp = ',trmp,' m2/s2'
       write (dpu,'(a,f14.7,a)')  ' trmq = ',trmq,' m3/s3'
       write (dpu,'(a,f14.7,a)')  ' qq   = ',  qq,' m6/s6'
       write (dpu,'(a,f14.7)')  ' b1   = ', b1
       write (dpu,'(a,f14.7)')  ' a1l  = ', a1l
       write (dpu,'(a,f14.7)')  ' vys  = ', vys
       write (dpu,'(a,f14.7)')  ' vus  = ', vus
       write (dpu,'(a,f14.7)')  ' vyb  = ', vyb
       write (dpu,'(a,f14.7)')  ' vub  = ', vub
       write (dpu,'(a,f14.7,a)') ' cwp  = ', 1000.*lwp,     &
                                 ' g/m2'
       write (dpu,'(a,f14.7)')  ' opt_depth  = ', opt_depth
       write (dpu,'(a,f14.7)')  ' radinvfrac = ',radinvfrac
       write (dpu,'(a,f14.7,a)')  ' radf = ', radf,' m2/s3' 
     end if
    
     !-----------------------------------------------------------
     ! limit CL-avg TKE used for entrainment
    
     ebrk(ncv) = rootp**2    
     ebrk(ncv) = max(min(ebrk(ncv),tkemax),tkemin) 
     wbrk(ncv) = ebrk(ncv)/b1

     if (column_match) then
       write (dpu,'(a,f14.7,a)')  ' rootp**2 = ',  rootp**2, ' m2/s2'
       write (dpu,'(a,f14.7,a)')  ' ebrk     = ', ebrk(ncv), ' m2/s2'
       write (dpu,'(a,f14.7,a)')  ' wbrk     = ', wbrk(ncv), ' m2/s2'
     end if
     
    !-----------------------------------------------------------
    ! ebrk should be greater than zero so if it is not a FATAL
    ! call is implemented

     if (ebrk(ncv) .eq. 0.) then
       call error_mesg ('edt_mod', &
                        'convective layer average tke = 0',  FATAL)
     end if
    
    !-----------------------------------------------------------
    ! Compute adjustment time for layer equal to lbulk / <e>
    ! Should this be divided by "c" = b1/mu, which for
    ! default value = 5.8/70 = 0.083
    
     do k = kb, kt, -1
       adj_time_inv(k) = sqrt(ebrk(ncv))/lbulk
     enddo
    
    !-----------------------------------------------------------
    ! We approximate TKE = <e> at entrainment interfaces con-
    ! sistent with entrainment closure.
            
     rcap(kt) = 1.   
     rcap(kb) = 1.   

    !-----------------------------------------------------------
    ! Calculate ratio rcap = e/<e> in convective layer interior. 
    ! Bound it by limits rmin = 0.1 to rmax = 2.0 to take care 
    ! of some pathological cases.

     if ((kb-kt).gt.1) then
       do k = kb-1, kt+1, -1
         rcap(k) = (mu*leng(k)/lbulk + wcap(k)/wbrk(ncv)) /    &
                   (mu*leng(k)/lbulk + 1.               )
         rcap(k) = min(max(rcap(k),rmin), rmax)
       end do
     end if
    
    !-----------------------------------------------------------
    ! Compute TKE throughout CL, and bound by tkemin & tkemax.
    !
    ! Question by Steve Klein:
    !
    ! Does tke(kb) properly account if kb is the top interface
    ! of another convective layer? 
    !
    
     do k = kb, kt, -1
       tke(k) = max(min(ebrk(ncv)*rcap(k),tkemax),tkemin) 
     end do 
                        
    !-----------------------------------------------------------
    ! Compute CL interior diffusivities, buoyancy and shear 
    ! production
    
     if ((kb-kt).gt.1) then

       do k = kb-1, kt+1, -1

         kvh(k)        = leng(k)*sqrt(tke(k))*shcl(ncv)
         kvm(k)        = leng(k)*sqrt(tke(k))*smcl(ncv)
         bprod(k)      = - kvh(k)*n2(k)
         sprod(k)      =   kvm(k)*s2(k)
         trans(k)      = mu*(ebrk(ncv)-tke(k))*           &
                         adj_time_inv(k)/b1 
         diss(k)       = sqrt(tke(k)*tke(k)*tke(k))/b1/   &
                         leng(k)
         shvec(k)      = shcl(ncv)
         smvec(k)      = smcl(ncv)
         ghvec(k)      = ghcl(ncv)

         turbtype(2,k) = 1
         isturb(k)     = 1.
         isturb(k-1)   = 1.
  
      !-------------------------------------------------
      ! compute sdevs
      ! 
      ! The approximation used is that qt and qs(Tl) are 
      ! correlated with strength kappa. The sign of the 
      ! correlation is positive if the vertical 
      ! gradients to qt and sli are of the same sign, 
      ! and negative otherwise.
 
         if (k .eq. kdim+1) call error_mesg ('edt_mod',   &
              'trying to compute sdevs at the surface',   &
                            FATAL)
         if (k .eq. 1)  call error_mesg ('edt_mod',       &
              'trying to compute sdevs at the model top', &
                            FATAL)
     
         qtsltmp    = (qt(k-1) - qt(k))/(zm(k-1) - zm(k))
         slsltmp    = (sl(k-1) - sl(k))/(zm(k-1) - zm(k))  
         qsltmp     = 0.5 * ( qsl    (k-1) + qsl    (k) )
         dqsldtltmp = 0.5 * ( dqsldtl(k-1) + dqsldtl(k) )   
         hlefftmp   = 0.5 * ( hleff  (k-1) + hleff  (k) )

         sdevs(k) = ( (mesovar*qsltmp)**2.0) +            &
                    ( ( (kvh(k)/sqrt(tke(k))) *           &
                    (qtsltmp-(kappa*slsltmp*dqsldtltmp/cp_air)))**2.0)
         sdevs(k) = sqrt(sdevs(k)) /                      &
                    (1.+hlefftmp*dqsldtltmp/cp_air)   

         if (column_match) then
           write (dpu,'(a)')  ' ' 
           write (dpu,'(a)')  ' ' 
           write (dpu,'(a,i4)')  ' sigmas for level ',k
           write (dpu,'(a)')  ' ' 
           write (dpu,'(a,f14.7,a)')  ' sigmas = ', 1000.*  &
                                        sdevs(k), ' g/kg'
           write (dpu,'(a,f14.7)')  ' acoef = ', 1. /       &
                                     ( 1. + hlefftmp*dqsldtltmp/cp_air )
           write (dpu,'(a,f14.7,a)')  ' sigmas/a = ',       &
                                        1000.*sdevs(k)*(1. +hlefftmp*dqsldtltmp/cp_air),&
                                      ' g/kg'
           write (dpu,'(a,f14.7)'  )  ' mesovar = ', mesovar
           write (dpu,'(a,f14.7,a)')  ' mesovar*qsl = ',    &
                                        mesovar*qsltmp*1000.,' g/kg'
           write (dpu,'(a,f14.7,a)')  ' turb.fluct = ',1000.&
                                      *(kvh(k)/sqrt(tke(k)))*(qtsltmp-            &
                                       (kappa*slsltmp*dqsldtltmp/cp_air)) ,' g/kg'
           write (dpu,'(a,f14.7,a)')  ' (kvh/sqrt(tke))* '//&
                                      ' qtsltmp = ',1000.*(kvh(k)/sqrt(tke(k)))*  &
                                        qtsltmp ,' g/kg'
           write (dpu,'(a,f14.7,a)')  ' (kvh/sqrt(tke))* '//&
                                      ' (kappa*slsltmp*dqsldtltmp/cp_air) = ',1000.*  &
                                        (kvh(k)/sqrt(tke(k)))*(kappa*slsltmp*       &
                                         dqsldtltmp/cp_air) ,' g/kg'
           write (dpu,'(a)')  ' '
           write (dpu,'(a)')  ' '
         end if
 
       end do
    
    
       !------------------------------------------------------
       !
       ! set sdevs at kt and kb equal to kt+1 and kb-1 values 
       ! respectively
       !
   
       sdevs(kt) = sdevs(kt+1)
       sdevs(kb) = sdevs(kb-1)

     end if

     !-----------------------------------------------------------
     ! Compute diffusivity we*dz and some diagnostics at the 
     ! upper inversion. Limit entrainment rate below the free 
     ! entrainment limit a1l * sqrt(e)

     kentr          = jtzm * a1l * sqrt(ebrk(ncv)) * &
                      min(evhc * ebrk(ncv)/(jtbu*leng(kt)),1.)
     kvh(kt)        = kentr
     kvm(kt)        = kentr
     bprod(kt)      = -kentr*n2h+radf
     sprod(kt)      =  kentr*s2(kt)
     trans(kt)      = mu*(ebrk(ncv)-tke(kt))*adj_time_inv(kt)/b1
     diss(kt)       = sqrt(tke(kt)*tke(kt)*tke(kt))/b1/leng(kt)
     radfvec(kt)    = radf
     turbtype(4,kt) = 1
     isturb(kt)     = 1.
    
     !-----------------------------------------------------------
     ! set isturb to 1 in the ambiguous layer
    
     isturb(kt-1) = 1.

     if (column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' at upper inversion: '
       write (dpu,'(a,f14.7,a)')  ' kentr = ', kentr, ' m2/s'
       write (dpu,'(a,f14.7)')  ' mu    = ', mu
       write (dpu,'(a,f14.7,a)')  ' 1/adj_time_inv = ',           &
                                    1./adj_time_inv(kt), ' sec'
       write (dpu,'(a)')  ' '
     end if
      
     !-----------------------------------------------------------
     ! compute sdevs at CL top
     !
     ! FOR SOME REASON THIS DERIVATION DOES NOT WORK
     !
     ! note that a water perturbation q* and buoyancy 
     ! perturbation b* are computed from the below equations:
     !
     !      q*   =    (w'qt') / sqrt(tke)
     !           =   - kentr * (dqt/dz) / sqrt(tke)
     !
     !      b*   =    (w'b')  / sqrt(tke)
     !           =    bprod   / sqrt(tke)
     !
     !
            
     qtsltmp     = (qt(kt-1) - qt(kt))/(zm(kt-1) - zm(kt)) 
     qstartmp    =  - kentr * qtsltmp / sqrt(tke(kt))
     bstartmp    =    bprod(kt)       / sqrt(tke(kt)) 
     temperature = pm(kt)/density(kt)/rdgas
    
     !sdevs(kt)   =  ((mesovar*qsl(kt))**2.0) +       &
     !     ((qstartmp-dqsldtl(kt)*bstartmp*temperature/grav)**2.0)
     !sdevs(kt) = sqrt(sdevs(kt))/(1.+hleff(kt)*dqsldtl(kt)/cp_air)
    
    
     if ((kb-kt).le.1) sdevs(kt) = mesovar*qsl(kt)/             &
                      (1.+hleff(kt)*dqsldtl(kt)/cp_air)
   

     if (column_match) then
       write (dpu,'(a)')  ' ' 
       write (dpu,'(a)')  ' ' 
       write (dpu,'(a,i4)')  ' sigmas at kt level # ',kt
       write (dpu,'(a)')  ' ' 
       write (dpu,'(a,f14.7,a)')  ' sigmas = ', 1000.*sdevs(kt),  &
                                  ' g/kg'
       write (dpu,'(a,f14.7)')    ' acoef = ', 1. /(1.+hleff(kt)* &
                                    dqsldtl(kt)/cp_air )                      
       write (dpu,'(a,f14.7,a)')  ' sigmas/a = ', 1000.*sdevs(kt)*&
                                  ( 1. + hleff(kt)*dqsldtl(kt)/cp_air ), ' g/kg'
       write (dpu,'(a,f14.7)'  )  ' mesovar = ', mesovar
       write (dpu,'(a,f14.7,a)')  ' mesovar*qsl = ',mesovar*      &
                                    qsl(kt)*1000.,' g/kg'
       write (dpu,'(a,f14.7,a)')  ' turb.fluct = ',1000.* &
                                   (qstartmp-dqsldtl(kt)*bstartmp*temperature/grav) ,    &
                                  ' g/kg'
       write (dpu,'(a,f14.7,a)')  ' temperature = ',temperature,  ' K'
       write (dpu,'(a,f14.7,a)')  ' qstartmp = ',1000.*qstartmp , ' g/kg'
       write (dpu,'(a,f14.7,a)')  ' bstartmp = ',bstartmp ,' m/s2'
       write (dpu,'(a,f14.7,a)')  ' jtbu = ',jtbu ,' m/s2'
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' '
     end if
                          
    !-----------------------------------------------------------
    ! Compute Kh, Km and some diagnostics at the lower inversion

     if (kb .lt. kdim+1) then 
          
       kentr        = jbzm * a1l * sqrt(ebrk(ncv)) * &
                                min(ebrk(ncv)/(jbbu*leng(kb)),1.)

    !------------------------------------------------------
    ! The addition of the previous value of kvh and khm 
    ! handles the case of 2 CLs entraining into each other
       kvh(kb)      = kentr + kvh(kb)   
       kvm(kb)      = kentr + kvm(kb)   
 
       bprod(kb)    = bprod(kb) - kvh(kb)*n2(kb)
       sprod(kb)    = sprod(kb) + kvm(kb)*s2(kb)
       trans(kb)    = trans(kb) + mu*(ebrk(ncv)-tke(kb))*    &
                      adj_time_inv(kb)/b1
       diss(kb)     = diss(kb)+sqrt(tke(kb)*tke(kb)*tke(kb)) &
                      /b1/leng(kb)
       turbtype(3,kb) = 1
       isturb(kb-1)   = 1.

    !------------------------------------------------------
    ! set isturb to 1 in the ambiguous layer
         
       isturb(kb) = 1.

       if (column_match) then
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' at lower inversion: '
         write (dpu,'(a,f14.7,a)')  ' kentr = ', kentr, ' m2/s'
         write (dpu,'(a,f14.7)')  ' mu    = ', mu
         write (dpu,'(a,f14.7,a)')  ' adj_time_inv = ',        &
                                      1./adj_time_inv(kb), ' sec'
         write (dpu,'(a)')  ' '
       end if
    
    !------------------------------------------------------
    ! compute sdevs at CL bottom
    !
    ! note that same formulas are used here as for CL top
   
       qtsltmp     = (qt(kb-1) - qt(kb))/(zm(kb-1) - zm(kb)) 
       qstartmp    =  - kentr * qtsltmp / sqrt(tke(kb))
       bstartmp    =    bprod(kb)       / sqrt(tke(kb)) 
       temperature = pm(kb-1)/density(kb-1)/rdgas
 
         !sdevs(kb)   = ((mesovar*qsl(kb-1))**2.0) + ((qstartmp-&
 !              dqsldtl(kb-1)*bstartmp*temperature/grav)&
 !       **2.0)
                 !sdevs(kb)   = sqrt(sdevs(kb))/(1.+hleff(kb-1)*        &
 !                               dqsldtl(kb-1)/cp_air)
                         
       if ((kb-kt).le.1) sdevs(kb) = mesovar*qsl(kb-1)/       &
                                     (1.+hleff(kb-1)*dqsldtl(kb-1)/cp_air)

       if (column_match) then
         write (dpu,'(a)')  ' ' 
         write (dpu,'(a)')  ' ' 
         write (dpu,'(a,i5)')  ' sigmas at lower inversion l'//&
                               'evel # ',kb
         write (dpu,'(a)')  ' ' 
         write (dpu,'(a,f14.7,a)')  ' sigmas = ', 1000.*       &
                                      sdevs(kb), ' g/kg'
         write (dpu,'(a,f14.7)  ')  ' acoef = ', 1. /(1.+      &
                                      hleff(kb-1)*dqsldtl(kb-1)/cp_air )                         
         write (dpu,'(a,f14.7,a)')  ' sigmas/a = ', 1000.*     &
                                      sdevs(kb)*( 1. + hleff(kb-1)*dqsldtl(kb-1)/cp_air ), &
                                    ' g/kg'
         write (dpu,'(a,f14.7)'  )  ' mesovar = ', mesovar
         write (dpu,'(a,f14.7,a)')  ' mesovar*qsl = ',mesovar* &
                                      qsl(kb-1)*1000.,' g/kg'
         write (dpu,'(a,f14.7,a)')  ' turb.fluct = ',1000.*    &
                                     (qstartmp-dqsldtl(kb-1)*bstartmp*temperature/    &
                                      grav) , ' g/kg'
         write (dpu,'(a,f14.7,a)')  ' temperature = ',         &
                                      temperature,' K'
         write (dpu,'(a,f14.7,a)')  ' qstartmp = ',1000.*      &
                                      qstartmp ,' g/kg'
         write (dpu,'(a,f14.7,a)')  ' bstartmp = ',bstartmp ,  ' m/s2'
         write (dpu,'(a,f14.7,a)')  ' jbbu = ',jbbu ,' m/s2'
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' '
       end if

     end if

    !-----------------------------------------------------------
    ! put minimum threshold on TKE to prevent possible division 
    ! by zero.

     if (kb .lt. kdim+1) then
       wcap(kb) = bprod(kb)*leng(kb)/sqrt(max(tke(kb),tkemin))
     else
       wcap(kb) = tkes / b1
     end if
     wcap(kt) = bprod(kt)*leng(kt) / sqrt(max(tke(kt),tkemin))

     if (column_match) then
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  'k     leng     rcap     wcap      tk'//&
            'e'
       write (dpu,'(a)')  '      (m)              (m2/s2)   (m2'//&
            '/s2)'
       write (dpu,'(a)')  '------------------------------------'//&
            '----'
       write (dpu,'(a)')  ' '
       do k = kb,kt,-1
         write(dpu,34) k,leng(k),rcap(k),wcap(k),tke(k)
       enddo  
34     format(1X,i2,1X,4(f8.4,1X))
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' k     kvh      kvm      bprod      '//&
            '  sprod         trans         diss'
       write (dpu,'(a)')  '      (m2/s)   (m2/s)   (m2/s3)     '//&
            ' (m2/s3)       (m2/s3)      (m2/s3)'
       write (dpu,'(a)')  '------------------------------------'//&
            '-------------------------------------'
       write (dpu,'(a)')  ' '
       do k = kb,kt,-1
         write(dpu,35) k,kvh(k),kvm(k),bprod(k),sprod(k),trans(k),  &
         diss(k)
       enddo  
35     format(1X,i2,1X,2(f8.4,1X),4(f12.9,1X))
       write (dpu,'(a)')  ' '
       write (dpu,'(a)')  ' '
     end if
                 
  !----------------------------------------------------------------
  ! End big loop over ncv
 
   end do      

!-----------------------------------------------------------------------
!
!      If the lowest CL reaches the surface, define the PBL 
!      depth as the CL top.
!
   if (ncvsurf .gt. 0) then
     pblh   = zi(ktop(ncvsurf))
     if(bflxs.ge.0.) then
       turbtype(2,kdim+1) = 1
     else
       turbtype(3,kdim+1) = 1
     end if
   else
     pblh = 0.
   end if
                        
       
!-----------------------------------------------------------------------
!
!      STABLE TURBULENT LAYERS
!

!----------------------------------------------------------------
! find turbulent lengthscales in all stable turbulent layers

   belongst(1) = .false.   ! k = 1 assumed nonturbulent
       
   any_stable  = .false.
       
   do k = 2, kdim
         
     belongst(k)  = (ri(k) .lt. ricrit) .and. (.not. belongcv(k))
     if (belongst(k)) any_stable = .true.
     if (belongst(k) .and. (.not.belongst(k-1)) ) then
       kt    = k     ! Top of stable turb layer
     else if (.not. belongst(k) .and. belongst(k-1)) then
       kb    = k-1   ! Base of stable turb layer
       lbulk = zm(kt-1) - zm(kb)
       do ks = kt, kb
         leng(ks)=lengthscale(zi(ks),lbulk)
         adj_time_inv(ks) = 1. / lbulk
       end do
     end if
  
   end do ! k

   !----------------------------------------------------------------
   ! Now look whether stable turb layer extends to ground. Note that 
   ! interface kdim+1 is assumed to always be stable-turbulent if it 
   ! is not convective. Note that if it is convective, kdim will 
   ! also be convective, so the above loop will have finished 
   ! finding all elevated turbulent layers.

   belongst(kdim+1) = .not. belongcv(kdim+1)
   
   if (belongst(kdim+1)) then  
   
     turbtype(1,kdim+1) = 1
     if (belongst(kdim)) then
     !------------------------------------------------------
     ! surface stable layer includes interior stable turb-
     ! ulent interface kt already defined above. Note that
     ! zm(kb) = 0.
       lbulk = zm(kt-1)     
     else                     
     !------------------------------------------------------
     ! surface stable BL with no interior turbulence            
       kt = kdim+1 
     end if
     lbulk  = zm(kt-1)
     pblh   = lbulk   ! PBL Height <-> lowest stable turb. layer            
     do ks = kt,kdim+1
       leng(ks) = lengthscale(zi(ks),lbulk)
       adj_time_inv(ks) = 1./ lbulk
     end do 
     adj_time_inv(kdim+1) = adj_time_inv(kdim+1)*sqrt(tkes)
   end if  ! for belongst if

   !----------------------------------------------------------------
   ! Calculate tke, kvh, kvm

   if (column_match .and. any_stable) then
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' STABLE LAYERS '
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' k     kvh      kvm      leng      t'//&
                        'ke      bprod        sprod        diss'
     write (dpu,'(a)')  '      (m2/s)    (m2/s)    (m)    (m2'//&
                        '/s2)   (m2/s3)      (m2/s3)      (m2/s3) '
     write (dpu,'(a)')  '------------------------------------'//&
                        '-----------------------------------------'
     write (dpu,'(a)')  ' '            
   end if

   do k = 2, kdim
     if (belongst(k)) then    
       turbtype(1,k) = 1
       isturb(k)     = 1.
       isturb(k-1)   = 1.
       call galperin(ri(k),gh,ckh,ckm)
 
 ! note that original cb code did not limit the maximum
 ! gh to 0.0233 so that we should check that if the  
 ! limiter was invoked
!!RSH drop this fail for now (maybe ok during early spinup ?? )
!                  if (gh .eq. 0.0233) call error_mesg ('edt_mod',&
!                    'galperin stability fn = 0.0233',&
!                                      FATAL)
       
                 ! ri(k) should be less than ricrit     
       if (ri(k) .gt. ricrit) call error_mesg ('edt_mod',&
                                  'ri(k) > ricrit, but belongst(k) = T',FATAL)
 
       tke(k) = b1*(leng(k)**2)*(-ckh*n2(k)+ckm*s2(k))
       tke(k) = max(min(tke(k),tkemax),tkemin)
       kvh(k) = leng(k) * sqrt(tke(k)) * ckh
       kvm(k) = leng(k) * sqrt(tke(k)) * ckm
       bprod(k)  = - kvh(k)*n2(k)
       sprod(k)  =   kvm(k)*s2(k)
       diss (k)  = sqrt(tke(k)*tke(k)*tke(k))/b1/leng(k)
       adj_time_inv(k) = adj_time_inv(k)*sqrt(tke(k))
       shvec(k)  = ckh
       smvec(k)  = ckm
       ghvec(k)  = gh
   
       if (column_match) then         
         write(dpu,37) k,kvh(k),kvm(k),leng(k),tke(k),bprod(k),&
                      sprod(k),diss(k)
37       format(1X,i2,1X,3(f8.4,1X),4(f12.9,1X))
       end if

                 !------------------------------------------------------
 ! compute sdevs 
 ! (note same formulas as used in CL calculations)
 ! 

       if (k .eq. kdim+1)  call error_mesg ('edt_mod',       &
                          'trying to compute sdevs at the surface', FATAL)
       if (k .eq. 1)       call error_mesg ('edt_mod',       &
                          'trying to compute sdevs at the model top', FATAL)
     
       qtsltmp    = (qt(k-1) - qt(k))/(zm(k-1) - zm(k))
       slsltmp    = (sl(k-1) - sl(k))/(zm(k-1) - zm(k))  
       qsltmp     = 0.5 * ( qsl    (k-1) + qsl    (k) )
       dqsldtltmp = 0.5 * ( dqsldtl(k-1) + dqsldtl(k) )   
       hlefftmp   = 0.5 * ( hleff  (k-1) + hleff  (k) )
     
       sdevs(k) = ( (mesovar*qsltmp)**2.0) +                 &
                  ( (  (kvh(k)/sqrt(tke(k))) *               &
                  (qtsltmp-(kappa*slsltmp*dqsldtltmp/cp_air)) )&
                   **2.0) 
       sdevs(k) = sqrt(sdevs(k)) / (1.+hlefftmp*dqsldtltmp/cp_air)   

       if (column_match) then
         write (dpu,'(a)')  ' ' 
         write (dpu,'(a)')  ' ' 
         write (dpu,'(a)')  ' sigmas  .... '
         write (dpu,'(a)')  ' ' 
         write (dpu,'(a,f14.7,a)')  ' sigmas = ',1000.*sdevs(k), ' g/kg'
         write (dpu,'(a,f14.7)  ')  ' acoef = ', 1. / ( 1. +   &
                                      hlefftmp*dqsldtltmp/cp_air )
         write (dpu,'(a,f14.7,a)')  ' sigmas/a = ', 1000.*     &
                                      sdevs(k)* ( 1. + hlefftmp*dqsldtltmp/cp_air ), ' g/kg'
         write (dpu,'(a,f14.7)'  )  ' mesovar = ', mesovar
         write (dpu,'(a,f14.7,a)')  ' mesovar*qsl = ',mesovar* &
                                      qsltmp*1000.,' g/kg'
         write (dpu,'(a,f14.7,a)')  ' turb.fluct = ',1000.*    &
                                     (kvh(k)/sqrt(tke(k)))*(qtsltmp-(kappa*slsltmp*   &
                                      dqsldtltmp/cp_air)) ,' g/kg'
         write (dpu,'(a,f14.7,a)')  ' (kvh/sqrt(tke))*qtsltm'//&
                                    'p = ',1000.*(kvh(k)/sqrt(tke(k)))*qtsltmp ,' g/kg'
         write (dpu,'(a,f14.7,a)')  ' (kvh/sqrt(tke))*(kappa'//&
                                    '*slsltmp*dqsldtltmp/cp_air) = ',1000.*(kvh(k)/sqrt  &
                                    (tke(k)))*(kappa*slsltmp*dqsldtltmp/cp_air),' g/kg'
         write (dpu,'(a)')  ' '
         write (dpu,'(a)')  ' '
       end if
  

     end if  ! belongs to stable layer
            
   end do  ! k loop

!-----------------------------------------------------------------------
! 
!      set tke at surface equal to diagnostic variable tkes

   tke(kdim+1) = tkes

   if (column_match) then
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7,a)')  ' pblh = ', pblh, ' m'
     write (dpu,'(a,f14.7,a)')  ' tkes = ', tkes, ' m2/s2'
     write (dpu,'(a)')  ' '
     do k = kdim-n_print_levels,kdim
       write (dpu,'(a,i4,a,f14.7)') 'k = ',k,'; isturb =',isturb(k)
     enddo
     write (dpu,'(a)')  ' '
   end if
       
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine caleddy

!
!======================================================================= 


!======================================================================= 
!
!      subroutine galperin
!        
      
subroutine galperin(ricl,gh,sh,sm)
        
!
!-----------------------------------------------------------------------
!
!      Given a Richardson number ricl, return the stability functions
!      sh and sm calculated according Galperin (1982).
!
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!
!
  
real, intent(in)  :: ricl
real, intent(out) :: gh, sh, sm
  
! internal variables

real              :: ri, trma, trmb, trmc, det

!-----------------------------------------------------------------------
!
!      code
!

   ri   = min(ricl,0.163)
   trma = alph3*alph4*ri+2.*b1*(alph2-alph4*alph5*ri)
   trmb = ri*(alph3+alph4)+2.*b1*(-alph5*ri+alph1)
   trmc = ri
   det = max(trmb*trmb-4.*trma*trmc,0.)
   gh = (-trmb + sqrt(det))/2./trma
   gh = max(gh,-0.28)
   gh = min(gh,0.0233)
   sh = alph5 / (1.+alph3*gh)
   sm = (alph1 + alph2*gh)/(1.+alph3*gh)/(1.+alph4*gh)
          
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine galperin

!
!======================================================================= 


!======================================================================= 
!
!      function lengthscale
!        
      
function lengthscale(height,depth)
        
!
!-----------------------------------------------------------------------
!
!  Calculate the turbulent length scale given the depth of the
!  turbulent layer.  Near the surface, the lengthscale asymptotes
!  to vonkarm*height.
!
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!
!
 
real              :: lengthscale 
real, intent(in)  :: height,depth

!-----------------------------------------------------------------------
!
!  code
!

   lengthscale = vonkarm*height/(1.+(vonkarm*height/(tunl*depth)))       
          
!-----------------------------------------------------------------------
! 
!      function end
!

end function lengthscale

!
!=======================================================================

!======================================================================= 
!
!  subroutine gaussian_cloud
!    
!
!  this subroutine computes a new cloud fraction and cloud conden-
!  sate using the Gaussian cloud modell of Mellor (1977) and
!  Sommeria and Deardorff (1977). 
!    
!  In this derivation a background mesoscale variability to total
!  water variability equivalent to a fraction, mesovar, of qsl is 
!  assumed.
!

subroutine gaussian_cloud (qxtop, qxmid, qxbot, acoef, sigmasf, qalyr, &
                           qclyr, sfuh,  sflh)
 
!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
!
!      qxtop      saturation excess at the top of the layer 
!                 (kg wat/kg air)
!      qxmid      saturation excess at the midpoint of the layer 
!                 (kg wat/kg air)
!      qxbot      saturation excess at the bottom of the layer 
!                 (kg wat/kg air)
!      sigmasf    standardard devations of water perturbation s
!                 (kg water / kg air)
!      acoef      thermo coefficient = 1./(1.+ L*dqsdT/cp)
!
!      ------
!      output
!      ------
!
!      qalyr      layer mean cloud fraction  (fraction)
!      qclyr      layer mean cloud condensate (kg condensate/kg air)
!      sfuh       saturated fraction of the upper half of the layer
!      sflh       saturated fraction of the lower half of the layer                  
!
!      --------
!      internal
!      --------
!
!
!-----------------------------------------------------------------------

real,    intent(in) :: qxtop, qxmid, qxbot, acoef, sigmasf
real,    intent(out):: qalyr, qclyr, sfuh, sflh

real :: q1top,q1mid,q1bot,cftop,cfmid,cfbot,qctop,qcmid,qcbot
real :: qcuh,qclh

!-----------------------------------------------------------------------
!
!      initialize variables

   qalyr    = 0.0
   qclyr    = 0.0
   sfuh     = 0.0
   sflh     = 0.0
   q1top    = 0.0
   q1mid    = 0.0
   q1bot    = 0.0
   cftop    = 0.0
   cfmid    = 0.0
   cfbot    = 0.0
   qctop    = 0.0
   qcmid    = 0.0
   qcbot    = 0.0
   qcuh     = 0.0
   qclh     = 0.0

!-----------------------------------------------------------------------
!
!      compute cloud fraction and cloud condensate at layer top, 
!      midpoint and bottom

   q1top = acoef * qxtop / sigmasf 
   cftop = max ( 0.5 * ( 2. -  erfcc( q1top / sqrt(2.) ) ) , 0.0 )
   qctop = cftop * q1top + ( exp(-0.5*q1top*q1top) / sqrt(2*fpi) )
   qctop = sigmasf * max ( qctop, 0. )  
   
   q1mid = acoef * qxmid / sigmasf  
   cfmid = max ( 0.5 * ( 2. -  erfcc( q1mid / sqrt(2.) ) ) , 0.0 )
   qcmid = cfmid * q1mid + ( exp(-0.5*q1mid*q1mid) / sqrt(2*fpi) )
   qcmid = sigmasf * max ( qcmid, 0. )   
   
   q1bot = acoef * qxbot / sigmasf  
   cfbot = max ( 0.5 * ( 2. -  erfcc( q1bot / sqrt(2.) ) ) , 0.0 )
   qcbot = cfbot * q1bot + ( exp(-0.5*q1bot*q1bot) / sqrt(2*fpi) )
   qcbot = sigmasf * max ( qcbot, 0. ) 
   
   sfuh  = 0.5 * (cftop + cfmid)
   sflh  = 0.5 * (cfmid + cfbot)
   qcuh  = 0.5 * (qctop + qcmid)
   qclh  = 0.5 * (qcmid + qcbot)
   qalyr = 0.5 * ( sfuh + sflh )
   qclyr = 0.5 * ( qcuh + qclh )

   if (qalyr .lt. qcminfrac) then
     qalyr = 0.
     qclyr = 0.
   end if

   if (column_match) then

     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' ========================================'//&
                        '======= '
     write (dpu,'(a)')  '             GAUSSIAN CLOUD MODEL         '
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7,a)')  ' sigmas = ', 1000.*sigmasf, ' g/kg'
     write (dpu,'(a,f14.7)'  )  ' acoef  = ', acoef
     write (dpu,'(a,f14.7,a)')  ' qxtop = ', qxtop*1000.,' g/kg'
     write (dpu,'(a,f14.7)')    ' q1top = ', q1top
     write (dpu,'(a,f14.7)')    ' cftop = ', cftop
     write (dpu,'(a,f14.7,a)')  ' qctop = ', 1000.*qctop,' g/kg'
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7,a)')  ' qxmid = ', qxmid*1000.,' g/kg'
     write (dpu,'(a,f14.7)')    ' q1mid = ', q1mid
     write (dpu,'(a,f14.7)')    ' cfmid = ', cfmid
     write (dpu,'(a,f14.7,a)')  ' qcmid = ', 1000.*qcmid,' g/kg'
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7,a)')  ' qxbot = ', qxbot*1000.,' g/kg'
     write (dpu,'(a,f14.7)')    ' q1bot = ', q1bot
     write (dpu,'(a,f14.7)')    ' cfbot = ', cfbot
     write (dpu,'(a,f14.7,a)')  ' qcbot = ', 1000.*qcbot,' g/kg'
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7)')    ' sfuh = ', sfuh
     write (dpu,'(a,f14.7)')    ' sflh = ', sflh
     write (dpu,'(a,f14.7,a)')  ' qcuh = ', 1000.*qcuh, ' g/kg'
     write (dpu,'(a,f14.7,a)')  ' qclh = ', 1000.*qclh, ' g/kg'
     write (dpu,'(a)')  ' '
     write (dpu,'(a,f14.7)')    ' qalyr = ', qalyr
     write (dpu,'(a,f14.7,a)')  ' qclyr = ', 1000.*qclyr, ' g/kg'
     write (dpu,'(a)')  ' '
     write (dpu,'(a)')  ' '
   end if    
      
       
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine gaussian_cloud

!
!======================================================================= 


!======================================================================= 
!
!  function erfcc
!        
      
function erfcc(x)
        
!
!-----------------------------------------------------------------------
!
!  This numerical recipes routine calculates the complementary
!  error function.
!
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!
!
 
real :: erfcc
real, intent(in) :: x 
real :: t,z


!-----------------------------------------------------------------------
!
!  code
!

   z=abs(x)      
   t=1./(1.+0.5*z)
  
   erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*       &
        (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*    &
        (1.48851587+t*(-.82215223+t*.17087277)))))))))
  
   if (x.lt.0.) erfcc=2.-erfcc
      
!-----------------------------------------------------------------------
! 
!  function end
!

end function erfcc

!
!=======================================================================

end module edt_mod
