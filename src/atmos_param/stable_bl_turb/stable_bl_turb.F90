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

  MODULE STABLE_BL_TURB_MOD

!=======================================================================
 use           mpp_mod, only: input_nml_file
 use           fms_Mod, ONLY: FILE_EXIST, OPEN_NAMELIST_FILE,          &
                              ERROR_MESG, FATAL, mpp_pe, mpp_root_pe,  &
                              CLOSE_FILE,                              &
                              check_nml_error, write_version_number,   &
                              stdlog
 use  Diag_Manager_Mod, ONLY: register_diag_field, send_data
 use  Time_Manager_Mod, ONLY: time_type
 use     Constants_Mod, ONLY: cp_air, hlv, hls, grav, vonkarm, tfreeze,&
                              rdgas, rvgas, omega
 use Monin_Obukhov_Mod, ONLY: stable_mix

 implicit none
 private
 public :: STABLE_BL_TURB, STABLE_BL_TURB_INIT, STABLE_BL_TURB_END


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  character(len=128) :: version = '$Id: stable_bl_turb.F90,v 19.0 2012/01/06 20:26:08 fms Exp $'
  character(len=128) :: tagname = '$Name: siena_201211 $'
  logical            :: module_is_initialized = .false.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------
! --- CONSTANTS
!---------------------------------------------------------------------

  real :: oalsm, oalsh
  real :: d608 = 0.
  
!---------------------------------------------------------------------
! --- NAMELIST
!---------------------------------------------------------------------

  real    :: akmax        = 1.e4 ! maximum diffusion coefficient value 
                                 !  (m2/s)
  real    :: alpha        = 0.5
  real    :: alsm         = 150.0
  real    :: alsh         = 150.0
  real    :: fmin         = 5.0e-5
  real    :: hpbl_cap     = 1000.
  real    :: ri_crit      = 0.2
  real    :: diff_min     = 0.001
  real    :: winddifmin   = 0.01
  real    :: small        = 1.e-5
  real    :: b_louis      = 9.4
  real    :: cmstar_louis = 7.4
  real    :: chstar_louis = 5.3
       
  NAMELIST / stable_bl_turb_nml / akmax, alpha, alsm, alsh, fmin, &
                                  hpbl_cap, diff_min, ri_crit

!---------------------------------------------------------------------
! DIAGNOSTICS FIELDS 
!---------------------------------------------------------------------

integer :: id_z_sbl, id_f_sbl

character(len=14) :: mod_name = 'stable_bl_turb'

real :: missing_value = -999.

!---------------------------------------------------------------------
 contains

!#######################################################################

! <SUBROUTINE NAME="STABLE_BL_TURB">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  STABLE_BL_TURB(is, js, Time, temp, qv,  ql,  qi,  um,  vm,
!                zhalf, zfull, u_star, b_star, lat, 
!                akm, akh, vspblcap, kbot )
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!       Starting integer for longitude window (used for diagnostics)
!  </IN>
!  <IN NAME="js" TYPE="integer">
!       Starting integer for latitude window (used for diagnostics)
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!       Time type variable (for diagnostics)
!  </IN>
!  <IN NAME="temp" TYPE="real">
!       Temperature (K)
!  </IN>
!  <IN NAME="qv" TYPE="real">
!       Water vapor specific humidity (kg/kg)
!  </IN>
!  <IN NAME="ql" TYPE="real">
!       Cloud liquid water specific humidity (kg/kg)
!  </IN>
!  <IN NAME="qi" TYPE="real">
!       Cloud ice water specific humidity (kg/kg)
!  </IN>
!  <IN NAME="um" TYPE="real">
!       Zonal wind velocity (m/s)
!  </IN>
!  <IN NAME="vm" TYPE="real">
!       Meridional wind velocity (m/s)
!  </IN>
!  <IN NAME="zhalf" TYPE="real">
!       Geopotential height of half levels (m)
!  </IN>
!  <IN NAME="zfull" TYPE="real">
!       Geopotential height of full levels (m)
!  </IN>
!  <IN NAME="u_star" TYPE="real">
!       Surface friction velocity (m/s)
!  </IN>
!  <IN NAME="b_star" TYPE="real">
!       Surface buoyancy scale (m/s2)
!  </IN>
!  <IN NAME="lat" TYPE="real">
!       Latitude (radians)
!  </IN>
!  <OUT NAME="akm" TYPE="real">
!       Momentum vertical diffusion coefficient (m2/s)
!  </OUT>
!  <OUT NAME="akh" TYPE="real">
!       Heat/Tracer vertical diffusion coefficient (m2/s)
!  </OUT>
!  <IN NAME="vspblcap" TYPE="real">
!       Cap to height of very stable enhanced PBL mixing, coming
!       from any other module (m)
!
!       In usual application this might be entrain_mod. This is also
!       an optional argument. 
!  </IN>
!  <IN NAME="kbot" TYPE="integer">
!       Integer indicating the lowest level above ground (integer)
!       
!       This optional argument is used only for eta coordinate model.
!  </IN>
! </SUBROUTINE>
!
 subroutine STABLE_BL_TURB( is, js, Time, temp, qv,  ql,  qi,  um,  vm,&
                            zhalf, zfull, u_star, b_star, lat,         &
                            akm, akh, vspblcap, kbot )

!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!
!       is, js   -  Starting indices for window
!       Time     -  Time used for diagnostics [time_type]
!       temp     -  temperature (K)
!       qv       -  specific humidity of water vapor (kg/kg)
!       ql       -  specific humidity of cloud liquid (kg/kg)
!       qi       -  specific humidity of cloud ice (kg/kg)
!       um, vm   -  Wind components (m/s)
!       u_star   -  surface friction velocity (m/s)
!       b_star   -  surface buoyancy (m/s2)
!       lat      -  latitude in radians
!       zhalf    -  Height at half levels (m)
!       zfull    -  Height at full levels (m)
!
!      --------------
!      optional input
!      --------------
!
!      vspblcap  cap to height of very stable PBL mixing, coming
!                from any other module (m) (in usual application
!                this might be entrain_mod or edt mod)
!
!      kbot      integer indicating the lowest true layer of atmosphere
!                this is used only for eta coordinate model
!
!---------------------------------------------------------------------
! Arguments (Intent out)
!       akm  -  mixing coefficient for momentum
!       akh  -  mixing coefficient for heat and moisture
!---------------------------------------------------------------------

  type(time_type), intent(in)                    :: Time
  integer,         intent(in)                    :: is, js
  real,            intent(in),  dimension(:,:)   :: u_star, b_star, lat
  real,            intent(in),  dimension(:,:,:) :: temp, qv, ql, qi
  real,            intent(in),  dimension(:,:,:) :: um, vm 
  real,            intent(in),  dimension(:,:,:) :: zhalf,  zfull
  real,            intent(out), dimension(:,:,:) :: akm,    akh

  real,     intent(in),   dimension(:,:), optional :: vspblcap
  integer,  intent(in),   dimension(:,:), optional :: kbot
  
!---------------------------------------------------------------------

  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3)-1) ::             &
        dsdzh, shear, buoync, Ri, Ritmp, fm, fh, lm, lh, xxm1, xxm2, &
        phi, zfunc, cmtmp, chtmp, fmtmp, fhtmp

  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3))  :: mask, hleff
  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3))  :: zfull_ag, slv
  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3)+1):: zhalf_ag

  real, dimension(SIZE(um,1),SIZE(um,2)) ::  zsurf,  &
        fcor, hpbl, z_sbl, f_sbl

  integer :: ix, jx, kx, i, j, k,  kxm
  integer :: shape1(1), shape3(3)
  logical :: used

!=======================================================================
! --- Initalize
!=======================================================================

! --- Check to see if STABLE_BL_TURB has been initialized
  if( .not. module_is_initialized ) CALL ERROR_MESG( ' STABLE_BL_TURB',                          &
       ' STABLE_BL_TURB_INIT has not been called', FATAL)

! --- Zero out output arrays
    akm(:,:,:) = 0.0
    akh(:,:,:) = 0.0
  z_sbl(:,:)   = 0.0
  f_sbl(:,:)   = 0.0
  
! --- Set dimensions etc
  ix  = SIZE( um, 1 )
  jx  = SIZE( um, 2 )
  kx  = SIZE( um, 3 )
  kxm = kx - 1

  shape1 =    ix * jx * kxm
  shape3 = (/ ix,  jx,  kxm /)

!====================================================================
! --- COMPUTE HEIGHT ABOVE SURFACE            
!====================================================================


       mask = 1.0
                   
       if (present(kbot)) then
            do j=1,jx
            do i=1,ix
                 zsurf(i,j) = zhalf(i,j,kbot(i,j)+1)
                 if (kbot(i,j).lt.kx) then
                    do k = kbot(i,j)+1,kx
                       mask(i,j,k) = 0.0
                    enddo
                 end if      
            enddo
            enddo
       else
            zsurf(:,:) = zhalf(:,:,kx+1)
       end if

       do k = 1, kx
            zfull_ag(:,:,k) = zfull(:,:,k) - zsurf(:,:)
            zhalf_ag(:,:,k) = zhalf(:,:,k) - zsurf(:,:)
       end do
       zhalf_ag(:,:,kx+1) = zhalf(:,:,kx+1) - zsurf(:,:)
       
!====================================================================
! --- DYNAMIC HEIGHT - also height relative to the surface     
!====================================================================

  fcor(:,:) = 2.0 * omega * SIN( lat(:,:) )
  fcor(:,:) = ABS( fcor(:,:) )
  fcor(:,:) = MAX( fcor(:,:), fmin )
  hpbl(:,:) = alpha * u_star(:,:) / fcor(:,:)

! --- bound
  hpbl(:,:) = MIN( hpbl(:,:), hpbl_cap )
  
! --- cap from entrainment turbulence
  if (present(vspblcap)) hpbl(:,:) = MIN( hpbl(:,:), vspblcap )
  
! --- height relative to the surface
! --- zfunc = zhalf_ag / hpbl, where stable conditions exist
! ---       = 1.0 for unstable conditions

  zfunc = 1.0
  do k = 1, kxm
      where( b_star(:,:) < 0.0) 
           zfunc(:,:,k) = min(1.,max(0.,zhalf_ag(:,:,k+1)/             &
                                     max(0.1,hpbl(:,:))))
      endwhere
  enddo
    
!====================================================================
! --- COMPUTE LIQUID WATER VIRTUAL STATIC ENERGY             
!====================================================================

   hleff   = (min(1.,max(0.,0.05*(temp   -tfreeze+20.)))*hlv + &
              min(1.,max(0.,0.05*(tfreeze -temp      )))*hls)
     
   slv     = cp_air*temp + grav*zfull_ag - hleff*(ql + qi)
   slv     = slv*(1+d608*(qv+ql+qi))
       
!====================================================================
! --- COMPUTE RICHARDSON NUMBER                 
!====================================================================

! --- D( )/DZ OPERATOR  
  
  dsdzh(:,:,1:kxm) = 1.0 / (zfull_ag(:,:,1:kxm) - zfull_ag(:,:,2:kx))

! --- WIND SHEAR SQUARED

  xxm1(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( um(:,:,1:kxm) - um(:,:,2:kx) )
  xxm2(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( vm(:,:,1:kxm) - vm(:,:,2:kx) )



  shear(:,:,:) = xxm1(:,:,:)*xxm1(:,:,:) + xxm2(:,:,:)*xxm2(:,:,:)

  where (shear .lt. (dsdzh*winddifmin*dsdzh*winddifmin)) 
         shear = dsdzh*winddifmin*dsdzh*winddifmin
  end where         

! --- BUOYANCY 
  xxm1(:,:,1:kxm) =       slv(:,:,1:kxm) - slv(:,:,2:kx) 
  xxm2(:,:,1:kxm) = 0.5*( slv(:,:,1:kxm) + slv(:,:,2:kx) )
 
 
  buoync(:,:,:) = grav * dsdzh(:,:,:) * xxm1(:,:,:) / xxm2(:,:,:)

! --- RICHARDSON NUMBER

  Ri(:,:,:) = buoync(:,:,:) / shear(:,:,:)   

!====================================================================
! --- MASK OUT UNDERGROUND VALUES FOR ETA COORDINATE
!====================================================================

  if( PRESENT( kbot ) ) then
     shear(:,:,1:kxm) =  shear(:,:,1:kxm) * mask(:,:,2:kx) 
    buoync(:,:,1:kxm) = buoync(:,:,1:kxm) * mask(:,:,2:kx) 
        Ri(:,:,1:kxm) =     Ri(:,:,1:kxm) * mask(:,:,2:kx) 
  endif

!====================================================================
! --- MIXING LENGTHS                 
!====================================================================

 do k = 1,kxm
   xxm1(:,:,k) = 1.0 / (vonkarm*zhalf_ag(:,:,k+1))
 end do 

  lm(:,:,:) = 1.0 / ( xxm1(:,:,1:kxm) + oalsm )
  lh(:,:,:) = 1.0 / ( xxm1(:,:,1:kxm) + oalsh )

!====================================================================
! --- STABILITY FUNCTIONS : STABLE SIDE       
!
! Note the very stable form of stability function acquired from 
! monin obukhov is weighted with the traditional stable form
! (phi = 1 + zeta/zeta_crit  or fm = (1 - Ri/Ricrit)**2)
! For Ricrit   = 0.2, phi retains the usual 1 + 5*zeta form.
!
! For heights greater than hpbl, the usual form is used.  For heights
! less than hpbl, the weight of the traditional form is given by
! zfunc which is linear in z/hpbl (see code above).
!====================================================================

  Ritmp = Ri
  where (Ritmp .lt. small) Ritmp = small

  CALL STABLE_MIX( Ritmp, phi)

  phi = (1-zfunc)*phi + zfunc* ((1-min(1.,(Ritmp/ri_crit)))**2.)
  
  fm(:,:,:) = phi(:,:,:)
  fh(:,:,:) =  fm(:,:,:) 
  
!====================================================================
! --- STABILITY FUNCTIONS : UNSTABLE SIDE (Louis 1979)
!
! f = 1.  - b * Ri / (1 + c*sqrt(-Ri))
!
! where b = 9.4 and
!
!              l * l * b * ( (  (1+(dz/z))**(1/3) - 1 )**(3/2))
! c = C_star * --------------------------------------------------
!              sqrt(z) * (dz**3/2)  
!
! where C_star(momentum) = 7.4, and C_star(heat) = 5.3
!     
!====================================================================
 
  Ritmp = Ri
  where (Ri .gt. 0.) Ritmp = 0.  
  
  zfunc(:,:,1:kxm) = 1.+(1./(dsdzh(:,:,1:kxm)*zhalf_ag(:,:,2:(kxm+1))))
  zfunc = zfunc **(1./3.) - 1.
  zfunc = zfunc **1.5
  zfunc = zfunc /  sqrt(zhalf_ag(:,:,2:(kxm+1))) 
  zfunc = zfunc * ( dsdzh(:,:,1:kxm) ** 1.5 )
  
  cmtmp =  cmstar_louis*lm(:,:,:)*lm(:,:,:)*b_louis*zfunc(:,:,:)
  chtmp =  chstar_louis*lh(:,:,:)*lh(:,:,:)*b_louis*zfunc(:,:,:)
  fmtmp(:,:,:) = 1. - (b_louis*Ritmp/(1.+cmtmp*sqrt(-1.*Ritmp)))
  fhtmp(:,:,:) = 1. - (b_louis*Ritmp/(1.+chtmp*sqrt(-1.*Ritmp)))
  
  where (Ri .lt. small)
      fm = fmtmp
      fh = fhtmp
  end where      
 
!====================================================================
! --- MIXING COEFFICENTS                 
!====================================================================

  shear(:,:,:) = SQRT( shear(:,:,:) )

! --- Momentum
  xxm1(:,:,:)    = lm(:,:,:) * lm(:,:,:) * fm(:,:,:)
   akm(:,:,2:kx) = xxm1(:,:,1:kxm) * shear(:,:,1:kxm) 
  where (akm .lt. diff_min) akm = 0.0
  where (akm .gt. akmax) akm = akmax
  
! --- Heat and Moisture
  xxm1(:,:,:)    = lh(:,:,:) * lh(:,:,:) * fh(:,:,:)
   akh(:,:,2:kx) = xxm1(:,:,1:kxm) * shear(:,:,1:kxm)
  where (akh .lt. diff_min) akh = 0.0
  where (akh .gt. akmax) akh = akmax

!====================================================================
! --- Extra diagnostics
!====================================================================

  where( b_star(:,:)  < 0.0 .and. hpbl (:,:) > 0.0 )
          z_sbl(:,:) = hpbl(:,:)
          f_sbl(:,:) = 1.0
  endwhere
  
  if ( id_z_sbl > 0 ) then
!     used = send_data ( id_z_sbl, z_sbl, Time, is, js )
     used = send_data ( id_z_sbl, z_sbl, Time)
  endif
  if ( id_f_sbl > 0 ) then
!     used = send_data ( id_f_sbl, f_sbl, Time, is, js )
     used = send_data ( id_f_sbl, f_sbl, Time)
  endif
  
!=======================================================================
  end subroutine STABLE_BL_TURB

!#######################################################################

! <SUBROUTINE NAME="STABLE_BL_TURB_INIT">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!     Initializes stable_bl_turb_mod: Reads and records namelist, 
!     sets up netcdf output if desired.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  STABLE_BL_TURB_INIT ( axes, Time )
!
!  </TEMPLATE>
!  <IN NAME=" axes" TYPE="integer">
!   Vector of axes integers
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   Time variable 
!  </IN>
! </SUBROUTINE>
!
 subroutine STABLE_BL_TURB_INIT ( axes, Time )
!=======================================================================
                   
 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time

 integer :: unit, io, ierr

!=======================================================================

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=stable_bl_turb_nml, iostat=io)
  ierr = check_nml_error(io,'stable_bl_turb_nml')
#else   
! -------------------------------------
  if( FILE_EXIST( 'input.nml' ) ) then
   unit = OPEN_NAMELIST_FILE ( file = 'input.nml')
   ierr = 1
   do while( ierr .ne. 0 )
   READ ( unit,  nml = stable_bl_turb_nml, iostat = io, end = 10 ) 
   ierr = check_nml_error (io, 'stable_bl_turb_nml')
   end do
10 continue
   CALL CLOSE_FILE( unit )
! -------------------------------------
  end if
#endif

!---------------------------------------------------------------------
! --- Output version
!---------------------------------------------------------------------

  if ( mpp_pe() == mpp_root_pe() ) then
       call write_version_number(version, tagname)
       unit = stdlog()
       WRITE( unit, nml = stable_bl_turb_nml ) 
  endif

!---------------------------------------------------------------------
! --- CONSTANTS
!---------------------------------------------------------------------
  d608 = (rvgas-rdgas)/rdgas
  oalsm = 1.0 / alsm
  oalsh = 1.0 / alsh

!---------------------------------------------------------------------
! --- initialize quantities for diagnostics output
!---------------------------------------------------------------------

   id_z_sbl = register_diag_field ( mod_name, &
     'z_sbl', axes(1:2), Time, &
     'Depth of stable boundary layer',              'm' )

   id_f_sbl = register_diag_field ( mod_name, &
     'f_sbl', axes(1:2), Time, &
     'Frequency of stable boundary layer',          ' ' )

!---------------------------------------------------------------------
 module_is_initialized = .true.
!=======================================================================
 end subroutine STABLE_BL_TURB_INIT

!#######################################################################

! <SUBROUTINE NAME="STABLE_BL_TURB_END">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    Closes down stable_bl_turb.  
!  </DESCRIPTION>
!  <TEMPLATE>
!   call STABLE_BL_TURB_END
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine STABLE_BL_TURB_END

!---------------------------------------------------------------------
 module_is_initialized = .false.
!=======================================================================
 end subroutine STABLE_BL_TURB_END

!#######################################################################
  end MODULE STABLE_BL_TURB_MOD
