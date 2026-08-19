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

module vert_turb_driver_mod

!-----------------------------------------------------------------------
!
!       driver for compuing vertical diffusion coefficients
!
!         choose either:
!              1) mellor-yamada 2.5 (with tke)
!              2) non-local K scheme
!              3) entrainment and diagnostic turbulence (edt) from
!                 Bretherton and Grenier
!
!-----------------------------------------------------------------------
!---------------- modules ---------------------


use      my25_turb_mod, only: my25_turb_init, my25_turb_end,  &
                              my25_turb, tke_surf, get_tke,   &
                              my25_turb_restart

use    diffusivity_mod, only: diffusivity, molecular_diff

use            edt_mod, only: edt_init, edt, edt_end

use    strat_cloud_mod, only: strat_cloud_on

use   shallow_conv_mod, only: shallow_conv_init, shallow_conv

use stable_bl_turb_mod, only: stable_bl_turb_init, stable_bl_turb

use        entrain_mod, only: entrain_init, entrain, entrain_end

use   diag_manager_mod, only: register_diag_field, send_data

use   time_manager_mod, only: time_type, get_time, operator(-)

use      constants_mod, only: rdgas, rvgas, kappa
 
use            mpp_mod, only: input_nml_file
use            fms_mod, only: mpp_pe, mpp_root_pe, stdlog, &
                              error_mesg, open_namelist_file, file_exist, &
                              check_nml_error, close_file, FATAL, &
                              write_version_number
 

use  field_manager_mod, only: MODEL_ATMOS

use tracer_manager_mod, only: get_tracer_index

implicit none
private

!---------------- interfaces ---------------------

public   vert_turb_driver_init, vert_turb_driver_end, vert_turb_driver
public   vert_turb_driver_restart


!-----------------------------------------------------------------------
!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: vert_turb_driver.F90,v 19.0 2012/01/06 20:27:33 fms Exp $'
character(len=128) :: tagname = '$Name:  $'
logical            :: module_is_initialized = .false.

!-----------------------------------------------------------------------
 real, parameter :: p00    = 1000.0E2
 real, parameter :: p00inv = 1./p00
 real :: d622   = 0.
 real :: d378   = 0.
 real :: d608   = 0.

!---------------- private data -------------------

 real :: gust_zi = 1000.   ! constant for computed gustiness (meters)

 integer :: nql, nqi, nqa    !  tracer indices for stratiform clouds

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

 logical :: do_shallow_conv  = .false.
 logical :: do_mellor_yamada = .true.
 logical :: do_diffusivity         = .false.
 logical :: do_molecular_diffusion = .false.
 logical :: do_edt                 = .false.
 logical :: do_stable_bl     = .false.
 logical :: use_tau          = .true.
 logical :: do_entrain    = .false.
 logical :: do_simple = .false. 

 character(len=24) :: gust_scheme  = 'constant' ! valid schemes are:
                                                !   => 'constant'
                                                !   => 'beljaars'
 real              :: constant_gust = 1.0
 real              :: gust_factor   = 1.0
 
 namelist /vert_turb_driver_nml/ do_shallow_conv, do_mellor_yamada, &
                                 gust_scheme, constant_gust, use_tau, &
                                 do_molecular_diffusion, do_stable_bl, &
                                 do_diffusivity, do_edt, do_entrain, &
                                 gust_factor, do_simple

!-------------------- diagnostics fields -------------------------------

integer :: id_tke,    id_lscale, id_lscale_0, id_z_pbl, id_gust,  &
           id_diff_t, id_diff_m, id_diff_sc, id_z_full, id_z_half,&
           id_uwnd,   id_vwnd,   id_diff_t_stab, id_diff_m_stab,  &
           id_diff_t_entr, id_diff_m_entr    

real :: missing_value = -999.

character(len=9) :: mod_name = 'vert_turb'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine vert_turb_driver (is, js, Time, Time_next, dt, tdtlw,     &
                             frac_land,   &
                             p_half, p_full, z_half, z_full, u_star,   &
                             b_star, q_star, rough, lat, convect,      &
                             u, v, t, q, r, um, vm, tm, qm, rm,        &
                             udt, vdt, tdt, qdt, rdt, ind_lcl, do_lcl_diffusivity_depth, diff_t, diff_m,  &
                             gust, z_pbl, mask, kbot           )

!-----------------------------------------------------------------------
integer,         intent(in)         :: is, js
type(time_type), intent(in)         :: Time, Time_next
   real,         intent(in)         :: dt
   real, intent(in), dimension(:,:) :: frac_land, u_star, b_star,  &
                                       q_star, rough, lat
logical, intent(in), dimension(:,:) :: convect       
   real, intent(in), dimension(:,:,:) :: tdtlw, p_half, p_full, &
                                         z_half, z_full, &
                                         u, v, t, q, um, vm, tm, qm, &
                                         udt, vdt, tdt, qdt
   real, intent(in) ,   dimension(:,:,:,:) :: r, rm, rdt
   integer, intent(in), dimension(:,:) :: ind_lcl
   logical, intent(in) :: do_lcl_diffusivity_depth
   real, intent(out),   dimension(:,:,:) :: diff_t, diff_m
   real, intent(out),   dimension(:,:)   :: gust, z_pbl 
   real, intent(in),optional, dimension(:,:,:) :: mask
integer, intent(in),optional, dimension(:,:) :: kbot
!-----------------------------------------------------------------------
real   , dimension(size(t,1),size(t,2),size(t,3))   :: ape, thv
logical, dimension(size(t,1),size(t,2),size(t,3)+1) :: lmask
real   , dimension(size(t,1),size(t,2),size(t,3)+1) :: el, diag3
real   , dimension(size(t,1),size(t,2),size(t,3)+1) :: tke
real   , dimension(size(t,1),size(t,2))             :: stbltop
real   , dimension(size(t,1),size(t,2))             :: el0, vspblcap
real   , dimension(size(diff_t,1),size(diff_t,2), &
                                  size(diff_t,3))   :: diff_sc,     &
                                                       diff_t_stab, &
                                                       diff_m_stab, &
       diff_t_entr, &
       diff_m_entr, &
       use_entr
real   , dimension(size(t,1),size(t,2),size(t,3))   :: tt, qq, uu, vv
real   , dimension(size(t,1),size(t,2),size(t,3))   :: qlin, qiin, qain
real    :: dt_tke
integer :: ie, je, nlev, sec, day, nt
logical :: used
!-----------------------------------------------------------------------
!----------------------- vertical turbulence ---------------------------
!-----------------------------------------------------------------------

      if (.not. module_is_initialized)  call error_mesg  &
                     ('vert_turb_driver in vert_turb_driver_mod',  &
                      'initialization has not been called', FATAL)

     nlev = size(p_full,3)
     ie = is + size(p_full,1) - 1
     je = js + size(p_full,2) - 1

!-----------------------------------------------------------------------
!---- set up state variable used by this module ----

      if (use_tau) then
      !-- variables at time tau
          uu = u
          vv = v
          tt = t
          qq = q
      else
      !-- variables at time tau+1
          uu = um + dt*udt
          vv = vm + dt*vdt
          tt = tm + dt*tdt
          qq = qm + dt*qdt
      endif

      !------ setup cloud variables: ql & qi & qa -----
      if (strat_cloud_on) then
           nt=size(r,4)
           if (nt == 0 .or. nt < max(nql,nqi,nqa))                    &
        call error_mesg ('vert_turb_driver',                  & 
                     'number of tracers less than nql or nqi or nqa', &
      FATAL) 
           if (use_tau) then
                qlin (:,:,:)=r(:,:,:,nql)
                qiin (:,:,:)=r(:,:,:,nqi)
                qain (:,:,:)=r(:,:,:,nqa)
           else
                qlin (:,:,:)=rm(:,:,:,nql)+rdt(:,:,:,nql)*dt
                qiin (:,:,:)=rm(:,:,:,nqi)+rdt(:,:,:,nqi)*dt
                qain (:,:,:)=rm(:,:,:,nqa)+rdt(:,:,:,nqa)*dt
           endif
      else
           qlin = 0.0
           qiin = 0.0
           qain = 0.0
      end if

!--------------------------------------------------------------------

!--------------------------------------------------------------------
! initialize output

   diff_t = 0.0
   diff_m = 0.0
   z_pbl = -999.0
   
!-------------------------------------------------------------------
! initiallize variables   
   vspblcap = 0.0   
   
!-----------------------------------------------------------------------
if (do_mellor_yamada) then

!    ----- virtual temp ----------
     ape(:,:,:)=(p_full(:,:,:)*p00inv)**(-kappa)
     if(do_simple) then 
       thv(:,:,:)=tt(:,:,:)*ape(:,:,:)
     else
       thv(:,:,:)=tt(:,:,:)*(qq(:,:,:)*d608+1.0)*ape(:,:,:)
     endif  
     if (present(mask)) where (mask < 0.5) thv = 200.

 endif

!---------------------------
 if (do_mellor_yamada) then
!---------------------------

!    ----- time step for prognostic tke calculation -----
     call get_time (Time_next-Time, sec, day)
     dt_tke = real(sec+day*86400)

!    --------------------- update tke-----------------------------------
!    ---- compute surface tke --------
!    ---- compute tke, master length scale (el0),  -------------
!    ---- length scale (el), and vert mix coeffs (diff_t,diff_m) ----

     call tke_surf  (is, js, u_star, kbot=kbot)



     if ( id_z_pbl > 0 ) then
     !------ compute pbl depth from k_profile if diagnostic needed -----
     call my25_turb (is, js, dt_tke, frac_land, p_half, p_full, thv, uu, vv, &
                     z_half, z_full, rough,   &
                     el0, el, diff_m, diff_t, &
                     mask=mask, kbot=kbot, &
                     ustar=u_star,bstar=b_star,h=z_pbl)
     else
     call my25_turb (is, js, dt_tke, frac_land, p_half, p_full, thv, uu, vv, &
                     z_half, z_full, rough,   &
                     el0, el, diff_m, diff_t, &
                     mask=mask, kbot=kbot)
     end if

!---------------------------
 else if (do_diffusivity) then
!--------------------------------------------------------------------
!----------- compute molecular diffusion, if desired  ---------------

    if (do_molecular_diffusion) then
      call molecular_diff (tt, p_half, diff_m, diff_t)
    else
      diff_m = 0.0
      diff_t = 0.0
    endif

!---------------------------
!------------------- non-local K scheme --------------

    if (do_lcl_diffusivity_depth) then 
      call diffusivity ( tt, qq, uu, vv, p_full, p_half, z_full, z_half,   &
                         u_star, b_star, z_pbl, diff_m, diff_t, &
                         ind_lcl = ind_lcl, kbot = kbot)
    else 
      call diffusivity ( tt, qq, uu, vv, p_full, p_half, z_full, z_half,   &
                         u_star, b_star, z_pbl, diff_m, diff_t, &
                         kbot = kbot)
    endif 
!---------------------------
else if (do_edt) then
!----------------------------

!    ----- time step for prognostic tke calculation -----
      call get_time (Time_next-Time, sec, day)
      dt_tke = real(sec+day*86400)
 

      tke = 0.0

    call edt(is,ie,js,je,dt_tke,Time_next,tdtlw, u_star,b_star,q_star, &
             tt,qq,  &
             qlin,qiin,qain,uu,vv,z_full,p_full,z_half,p_half,stbltop, &
             diff_m,diff_t,z_pbl,kbot=kbot,tke=tke)


 endif
 


 
!------------------------------------------------------------------
! --- boundary layer entrainment parameterization

   if( do_entrain ) then

       call entrain(is,ie,js,je,Time_next,tdtlw, convect,u_star,b_star,&
                    tt,qq, &
            qlin,qiin,qain,uu,vv,z_full,p_full,z_half,p_half,diff_m,   &
    diff_t,diff_m_entr,diff_t_entr,use_entr,z_pbl,vspblcap,    &
    kbot=kbot)
   
   endif

!-----------------------------------------------------------------------
! --- stable boundary layer parameterization

   if( do_stable_bl ) then

        if (do_entrain) then

CALL STABLE_BL_TURB( is, js, Time_next, tt, qq, qlin, qiin, uu,&
                     vv, z_half, z_full, u_star, b_star, lat,  &
     diff_m_stab, diff_t_stab,                 &
     vspblcap = vspblcap, kbot=kbot)
     
            diff_m = use_entr*diff_m_entr + (1-use_entr)*diff_m_stab
            diff_t = use_entr*diff_t_entr + (1-use_entr)*diff_t_stab
    
            !for diagnostic purposes only, save the stable_bl_turb
            !coefficient only where it was used
    
            diff_m_stab = (1-use_entr)*diff_m_stab
            diff_t_stab = (1-use_entr)*diff_t_stab    
         
else

CALL STABLE_BL_TURB( is, js, Time_next, tt, qq, qlin, qiin, uu,&
                     vv, z_half, z_full, u_star, b_star, lat,  &
     diff_m_stab, diff_t_stab,kbot=kbot)
     
            diff_m = diff_m +  MAX( diff_m_stab - diff_m, 0.0 )
            diff_t = diff_t +  MAX( diff_t_stab - diff_t, 0.0 )
    
end if
        
    endif
   
!-----------------------------------------------------------------------
!------------------ shallow convection ???? ----------------------------

   if (do_shallow_conv) then
        call shallow_conv (tt, qq, p_full, p_half, diff_sc, kbot)
        diff_t = diff_t + diff_sc
   endif

!-----------------------------------------------------------------------
!------------- define gustiness ------------

     if ( trim(gust_scheme) == 'constant' ) then
          gust = constant_gust
     else if ( trim(gust_scheme) == 'beljaars' ) then
!    --- from Beljaars (1994) and Beljaars and Viterbo (1999) ---
          where (b_star > 0.)
             gust = gust_factor * (u_star*b_star*gust_zi)**(1./3.)
          elsewhere
             gust = 0.
          endwhere
     endif

!-----------------------------------------------------------------------
!------------------------ diagnostics section --------------------------

if (do_mellor_yamada) then

!     --- set up local mask for fields with surface data ---
      if ( present(mask) ) then
         lmask(:,:,1)        = .true.
         lmask(:,:,2:nlev+1) = mask(:,:,1:nlev) > 0.5
      else
         lmask = .true.
      endif

!------- tke --------------------------------
      if ( id_tke > 0 ) then
         call get_tke(is,ie,js,je,tke)
         used = send_data ( id_tke, tke, Time_next, is, js, 1, &
                            mask=lmask )
      endif

!------- length scale (at half levels) ------
      if ( id_lscale > 0 ) then
         used = send_data ( id_lscale, el, Time_next, is, js, 1,  &
                            mask=lmask )
      endif

!------- master length scale -------
      if ( id_lscale_0 > 0 ) then
         used = send_data ( id_lscale_0, el0, Time_next, is, js )
      endif

end if

if (do_edt) then 
    
!     --- set up local mask for fields with surface data ---
    if ( present(mask) ) then
          lmask(:,:,1)        = .true.
          lmask(:,:,2:nlev+1) = mask(:,:,1:nlev) > 0.5
     else   
        lmask = .true.
       endif

!------- tke --------------------------------
      if ( id_tke > 0 ) then
        used = send_data ( id_tke, tke, Time_next, is, js, 1,     &
                          mask=lmask )
      endif
 
end if

!------- boundary layer depth -------
      if ( id_z_pbl > 0 ) then
         used = send_data ( id_z_pbl, z_pbl, Time_next, is, js )
      endif

!------- gustiness -------
      if ( id_gust > 0 ) then
         used = send_data ( id_gust, gust, Time_next, is, js )
      endif


!------- output diffusion coefficients ---------

  if ( id_diff_t > 0 .or. id_diff_m > 0 .or. id_diff_sc > 0 .or. &
       id_diff_t_stab > 0 .or. id_diff_m_stab > 0 .or.           &
       id_diff_t_entr > 0 .or. id_diff_m_entr > 0  ) then
!       --- set up local mask for fields without surface data ---
        if (present(mask)) then
            lmask(:,:,1:nlev) = mask(:,:,1:nlev) > 0.5
            lmask(:,:,nlev+1) = .false.
        else
            lmask(:,:,1:nlev) = .true.
            lmask(:,:,nlev+1) = .false.
        endif
!       -- dummy data at surface --
        diag3(:,:,nlev+1)=0.0
  endif

!------- diffusion coefficient for heat/moisture -------
   if ( id_diff_t > 0 ) then
      diag3(:,:,1:nlev) = diff_t(:,:,1:nlev)
      used = send_data ( id_diff_t, diag3, Time_next, is, js, 1, mask=lmask )
   endif

!------- diffusion coefficient for momentum -------
   if ( id_diff_m > 0 ) then
      diag3(:,:,1:nlev) = diff_m(:,:,1:nlev)
      used = send_data ( id_diff_m, diag3, Time_next, is, js, 1, mask=lmask )
   endif

!------- diffusion coefficient for shallow conv -------
 if (do_shallow_conv) then
   if ( id_diff_sc > 0 ) then
      diag3(:,:,1:nlev) = diff_sc(:,:,1:nlev)
      used = send_data ( id_diff_sc, diag3, Time_next, is, js, 1, mask=lmask)
   endif
 endif

!------- diffusion coefficients for stable boudary layer -------
   if (do_stable_bl) then
!------- for heat/moisture -------
    if ( id_diff_t_stab > 0 ) then
       diag3(:,:,1:nlev) = diff_t_stab(:,:,1:nlev)
      used = send_data ( id_diff_t_stab, diag3, Time_next, is, js, 1, mask=lmask )
  endif
!------- for momentum -------
    if ( id_diff_m_stab > 0 ) then
       diag3(:,:,1:nlev) = diff_m_stab(:,:,1:nlev)
     used = send_data ( id_diff_m_stab, diag3, Time_next, is, js, 1, mask=lmask )
    endif
 endif

!------- diffusion coefficients for entrainment module -------
 if (do_entrain) then
      if ( id_diff_t_entr > 0 ) then
       diag3(:,:,1:nlev) = diff_t_entr(:,:,1:nlev)
      used = send_data ( id_diff_t_entr, diag3, Time_next, is, js, 1, mask=lmask )
      endif
      if ( id_diff_m_entr > 0 ) then
       diag3(:,:,1:nlev) = diff_m_entr(:,:,1:nlev)
      used = send_data ( id_diff_m_entr, diag3, Time_next, is, js, 1, mask=lmask )
      endif
 endif

!--- geopotential height relative to the surface on full and half levels ----

   if ( id_z_half > 0 ) then
      !--- set up local mask for fields with surface data ---
      if ( present(mask) ) then
         lmask(:,:,1)        = .true.
         lmask(:,:,2:nlev+1) = mask(:,:,1:nlev) > 0.5
      else
         lmask = .true.
      endif
      used = send_data ( id_z_half, z_half, Time_next, is, js, 1, mask=lmask )
   endif
   
   if ( id_z_full > 0 ) then
      used = send_data ( id_z_full, z_full, Time_next, is, js, 1, rmask=mask)
   endif
   
!--- zonal and meridional wind on mass grid -------

   if ( id_uwnd > 0 ) then
      used = send_data ( id_uwnd, uu, Time_next, is, js, 1, rmask=mask)
   endif
  
   if ( id_vwnd > 0 ) then
      used = send_data ( id_vwnd, vv, Time_next, is, js, 1, rmask=mask)
   endif
  
 
   
!-----------------------------------------------------------------------

end subroutine vert_turb_driver

!#######################################################################

subroutine vert_turb_driver_init (lonb, latb, id, jd, kd, axes, Time, &
                                  doing_edt, doing_entrain)

!-----------------------------------------------------------------------
   real, dimension(:,:), intent(in) :: lonb, latb
   integer,         intent(in) :: id, jd, kd, axes(4)
   type(time_type), intent(in) :: Time
   logical,         intent(out) :: doing_edt, doing_entrain
!-----------------------------------------------------------------------
   integer, dimension(3) :: full = (/1,2,3/), half = (/1,2,4/)
   integer :: ierr, unit, io, logunit

      if (module_is_initialized)  &
          call error_mesg  &
                   ('vert_turb_driver_init in vert_turb_driver_mod',  &
                    'attempting to call initialization twice', FATAL)

!-----------------------------------------------------------------------
!--------------- read namelist ------------------

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=vert_turb_driver_nml, iostat=io)
   ierr = check_nml_error(io,'vert_turb_driver_nml')
#else   
      if (file_exist('input.nml')) then
         unit = open_namelist_file (file='input.nml')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=vert_turb_driver_nml, iostat=io, end=10)
            ierr = check_nml_error (io, 'vert_turb_driver_nml')
         enddo
  10     call close_file (unit)
      endif
#endif

!---------- output namelist --------------------------------------------

      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           write (logunit,nml=vert_turb_driver_nml)
      endif

!     --- check namelist option ---
      if ( trim(gust_scheme) /= 'constant' .and. &
           trim(gust_scheme) /= 'beljaars' ) call error_mesg &
         ('vert_turb_driver_mod', 'invalid value for namelist '//&
          'variable GUST_SCHEME', FATAL)

      if (do_molecular_diffusion .and. do_mellor_yamada)  &
         call error_mesg ( 'vert_turb_driver_mod', 'cannot activate '//&
              'molecular diffusion with mellor_yamada', FATAL)
 
       if (do_molecular_diffusion .and. do_edt)  &
         call error_mesg ( 'vert_turb_driver_mod', 'cannot activate '//&
           'molecular diffusion with EDT', FATAL)

!-----------------------------------------------------------------------
!s initialise these here as rdgas no longer a parameter
d622   = rdgas/rvgas
d378   = 1.-d622
d608   = d378/d622
!-----------------------------------------------------------------------
         
       if (strat_cloud_on) then
! get tracer indices for stratiform cloud variables
          nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
          nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
          nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
          if (mpp_pe() == mpp_root_pe()) &
                 write (logunit,'(a,3i4)') 'Stratiform cloud tracer indices: nql,nqi,nqa =',nql,nqi,nqa
          if (min(nql,nqi,nqa) <= 0) call error_mesg ('moist_processes', &
                         'stratiform cloud tracer(s) not found', FATAL)
          if (nql == nqi .or. nqa == nqi .or. nql == nqa) call error_mesg ('moist_processes',  &
       'tracers indices cannot be the same (i.e., nql=nqi=nqa).', FATAL)
      endif

!----------------------------------------------------------------------

      if (do_mellor_yamada) call my25_turb_init (id, jd, kd)

      if (do_shallow_conv)  call shallow_conv_init (kd)

      if (do_stable_bl)     call stable_bl_turb_init ( axes, Time )

      if (do_edt)           call edt_init (lonb, latb, axes,Time,id,jd,kd)

      if (do_entrain)       call entrain_init (lonb, latb, axes,Time,id,jd,kd)
      
!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----

   id_uwnd = register_diag_field ( mod_name, 'uwnd', axes(full), Time, &
        'zonal wind on mass grid', 'meters/second' ,                   &
         missing_value=missing_value    )

   id_vwnd = register_diag_field ( mod_name, 'vwnd', axes(full), Time, &
        'meridional wind on mass grid', 'meters/second' ,              &
        missing_value=missing_value    )

   id_z_full = &
   register_diag_field ( mod_name, 'z_full', axes(full), Time,    &
        'geopotential height relative to surface at full levels', &
         'meters' , missing_value=missing_value    )

   id_z_half = &
   register_diag_field ( mod_name, 'z_half', axes(half), Time,    &
        'geopotential height relative to surface at half levels', &
        'meters' , missing_value=missing_value    )

if (do_mellor_yamada) then

   id_tke = &
   register_diag_field ( mod_name, 'tke', axes(half), Time,      &
                        'turbulent kinetic energy',  'm2/s2'   , &
                        missing_value=missing_value               )

   id_lscale = &
   register_diag_field ( mod_name, 'lscale', axes(half), Time,    &
                        'turbulent length scale',  'm'   ,        &
                        missing_value=missing_value               )

   id_lscale_0 = &
   register_diag_field ( mod_name, 'lscale_0', axes(1:2), Time,   &
                        'master length scale',  'm'               )
endif

 if (do_edt) then
 
   id_tke = &
   register_diag_field ( mod_name, 'tke', axes(half), Time,      &
                         'turbulent kinetic energy',  'm2/s2'   , &
                         missing_value=missing_value               )
 
  end if

   id_z_pbl = &
   register_diag_field ( mod_name, 'z_pbl', axes(1:2), Time,       &
                        'depth of planetary boundary layer',  'm'  )

   id_gust = &
   register_diag_field ( mod_name, 'gust', axes(1:2), Time,        &
                        'wind gustiness in surface layer',  'm/s'  )

   id_diff_t = &
   register_diag_field ( mod_name, 'diff_t', axes(half), Time,    &
                        'vert diff coeff for temp',  'm2/s'   ,   &
                        missing_value=missing_value               )

   id_diff_m = &
   register_diag_field ( mod_name, 'diff_m', axes(half), Time,      &
                        'vert diff coeff for momentum',  'm2/s'   , &
                        missing_value=missing_value               )

if (do_shallow_conv) then

   id_diff_sc = &
   register_diag_field ( mod_name, 'diff_sc', axes(half), Time,      &
                        'vert diff coeff for shallow conv', 'm2/s' , &
                        missing_value=missing_value               )
endif

if (do_stable_bl) then
  id_diff_t_stab = &
    register_diag_field ( mod_name, 'diff_t_stab', axes(half), Time,       &
                       'vert diff coeff for temp',  'm2/s',                &
                        missing_value=missing_value               )

  id_diff_m_stab = &
    register_diag_field ( mod_name, 'diff_m_stab', axes(half), Time,       &
                       'vert diff coeff for momentum',  'm2/s',            &
                       missing_value=missing_value               )
 endif


if (do_entrain) then
  id_diff_m_entr = &
    register_diag_field ( mod_name, 'diff_m_entr', axes(half), Time,        &
            'momentum vert diff coeff from entrainment module',  'm2/s',    &
                        missing_value=missing_value               )

  id_diff_t_entr = &
    register_diag_field ( mod_name, 'diff_t_entr', axes(half), Time,        &
            'heat vert diff coeff from entrainment module',  'm2/s',        &
                        missing_value=missing_value               )

 endif


!-----------------------------------------------------------------------

   doing_edt = do_edt
   doing_entrain = do_entrain
   module_is_initialized =.true.

!-----------------------------------------------------------------------

end subroutine vert_turb_driver_init


!#######################################################################

subroutine vert_turb_driver_end

!-----------------------------------------------------------------------
      if (do_mellor_yamada) call my25_turb_end
      if (do_edt) call edt_end
      if (do_entrain) call entrain_end
      module_is_initialized =.false.

!-----------------------------------------------------------------------

end subroutine vert_turb_driver_end

!#######################################################################
! <SUBROUTINE NAME="vert_turb_driver_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine vert_turb_driver_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

   if (do_mellor_yamada) call my25_turb_restart(timestamp)
end subroutine vert_turb_driver_restart
! </SUBROUTINE> NAME="vert_turb_driver_restart"


end module vert_turb_driver_mod

