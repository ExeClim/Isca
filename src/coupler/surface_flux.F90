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

module surface_flux_mod
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">GFDL </CONTACT>
!
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
!
! <OVERVIEW>
!  Driver program for the calculation of fluxes on the exchange grids.
! </OVERVIEW>
!
! <DESCRIPTION>
!
! </DESCRIPTION>
!
! ============================================================================

use             fms_mod, only: FATAL, close_file, mpp_pe, mpp_root_pe, write_version_number
use             fms_mod, only: file_exist, check_nml_error, open_namelist_file, stdlog
use             fms_mod, only: error_mesg ! RUIZHI
use   monin_obukhov_mod, only: mo_drag, mo_profile
use  sat_vapor_pres_mod, only: escomp, descomp
use       constants_mod, only: cp_air, hlv, stefan, rdgas, rvgas, grav, vonkarm, dens_h2o
use             mpp_mod, only: input_nml_file

implicit none
private

! ==== public interface ======================================================
public  surface_flux, gp_surface_flux
! ==== end of public interface ===============================================

! <INTERFACE NAME="surface_flux">
!   <OVERVIEW>
!   For the calculation of fluxes on the exchange grids.
!   </OVERVIEW>
!   <DESCRIPTION>
!   For the calculation of fluxes on the exchange grids.
!   </DESCRIPTION>
!
!  <IN NAME="t_atm" TYPE="real, dimension(:)" UNITS="Kelvin">
!  Air temp lowest atmospheric level.
!  </IN>
!  <IN NAME="q_atm" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Mixing ratio at lowest atmospheric level (kg/kg).
!  </IN>
!  <IN NAME="u_atm" TYPE="real, dimension(:)" UNITS="m/s">
!  Zonal wind velocity at lowest atmospheric level.
!  </IN>
!  <IN NAME="v_atm" TYPE="real, dimension(:)" UNITS="m/s">
!  Meridional wind velocity at lowest atmospheric level.
!  </IN>
!  <IN NAME="p_atm" TYPE="real, dimension(:)" UNITS="Pascal">
!  Pressure lowest atmospheric level.
!  </IN>
!  <IN NAME="z_atm" TYPE="real, dimension(:)" UNITS="m" >
!  Height lowest atmospheric level.
!  </IN>
!  <IN NAME="p_surf" TYPE="real, dimension(:)" UNITS="Pascal">
!   Pressure at the earth's surface
!  </IN>
!  <IN NAME="t_surf" TYPE="real, dimension(:)" UNITS="Kelvin">
!   Temp at the earth's surface
!  </IN>
!  <IN NAME="t_ca" TYPE="real, dimension(:)" UNITS="Kelvin">
!   Air temp at the canopy
!  </IN>
!  <OUT NAME="q_surf" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Mixing ratio at earth surface (kg/kg).
!  </OUT>
!  <IN NAME="u_surf" TYPE="real, dimension(:)" UNITS="m/s">
!  Zonal wind velocity at earth surface.
!  </IN>
!  <IN NAME="v_surf" TYPE="real, dimension(:)" UNITS="m/s">
!  Meridional wind velocity at earth surface.
!  </IN>
!  <IN NAME="rough_mom" TYPE="real, dimension(:)" UNITS="m">
!  Momentum roughness length
!  </IN>
!  <IN NAME="rough_heat" TYPE="real, dimension(:)" UNITS="m">
!  Heat roughness length
!  </IN>
!  <IN NAME="rough_moist" TYPE="real, dimension(:)" UNITS="m">
! <Moisture roughness length
!  </IN>
!  <IN NAME="rough_scale" TYPE="real, dimension(:)" UNITS="dimensionless" >
!  Scale factor used to topographic roughness calculation
!  </IN>
!  <IN NAME="gust" TYPE="real, dimension(:)"  UNITS="m/s">
!   Gustiness factor
!  </IN>
!  <OUT NAME="flux_t" TYPE="real, dimension(:)" UNITS="W/m^2">
!  Sensible heat flux
!  </OUT>
!  <OUT NAME="flux_q" TYPE="real, dimension(:)" UNITS="kg/(m^2 s)">
!  Evaporative water flux
!  </OUT>
!  <OUT NAME="flux_r" TYPE="real, dimension(:)" UNITS="W/m^2">
!  Radiative energy flux
!  </OUT>
!  <OUT NAME="flux_u" TYPE="real, dimension(:)" UNITS="Pa">
!  Zonal momentum flux
!  </OUT>
!  <OUT NAME="flux_v" TYPE="real, dimension(:)" UNITS="Pa">
! Meridional momentum flux
!  </OUT>
!  <OUT NAME="cd_m" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Momentum exchange coefficient
!  </OUT>
!  <OUT NAME="cd_t" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Heat exchange coefficient
!  </OUT>
!  <OUT NAME="cd_q" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Moisture exchange coefficient
!  </OUT>
!  <OUT NAME="w_atm" TYPE="real, dimension(:)" UNITS="m/s">
!  Absolute wind at the lowest atmospheric level
!  </OUT>
!  <OUT NAME="u_star" TYPE="real, dimension(:)" UNITS="m/s">
!  Turbulent velocity scale
!  </OUT>
!  <OUT NAME="b_star" TYPE="real, dimension(:)" UNITS="m/s^2">
!  Turbulent buoyant scale
!  </OUT>
!  <OUT NAME="q_star" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Turbulent moisture scale
!  </OUT>
!  <OUT NAME="dhdt_surf" TYPE="real, dimension(:)" UNITS="(W/m^2)/K">
!  Sensible heat flux temperature sensitivity
!  </OUT>
!  <OUT NAME="dedt_surf" TYPE="real, dimension(:)" UNITS="1/K">
!   Moisture flux temperature sensitivity
!  </OUT>
!  <OUT NAME="dedq_surf" TYPE="real, dimension(:)" UNITS="(kg/m^2)/s">
!  Moisture flux humidity sensitivity
!  </OUT>
!  <OUT NAME="drdt_surf" TYPE="real, dimension(:)" UNITS="(W/m^2)/K">
!  Radiative energy flux temperature sensitivity
!  </OUT>
!  <OUT NAME="dhdt_atm" TYPE="real, dimension(:)" UNITS="(W/m^2)/K">
!  Derivative of sensible heat flux over temp at the lowest atmos level.
!  </OUT>
!  <OUT NAME="dedq_atm" TYPE="real, dimension(:)" UNITS="(kg/m^2/sec)/K">
!  Derivative of water vapor flux over temp at the lowest atmos level.
!  </OUT>
!  <OUT NAME="dtaudu_atm" TYPE="real, dimension(:)" UNITS="Pa/(m/s)">
!  Derivative of zonal wind stress w.r.t the lowest level zonal
!  wind speed of the atmos
!  </OUT>
!  <OUT NAME="dtaudv_atm" TYPE="real, dimension(:)" UNITS="Pa/(m/s)">
!  Derivative of meridional wind stress w.r.t the lowest level meridional
!  wind speed of the atmos
!  </OUT>
!  <OUT NAME="dt" TYPE="real">
!  Time step (it is not used presently)
!  </OUT>
!  <IN NAME="land" TYPE="logical, dimension(:)">
!  Indicates where land exists (true if exchange cell is on land).
!  </IN>
!  <IN NAME="seawater" TYPE="logical, dimension(:)">
!  Indicates where liquid ocean water exists
!  (true if exchange cell is on liquid ocean water).
!  </IN>
!  <IN NAME="avail" TYPE="logical, dimension(:)">
!  True where the exchange cell is active.
!  </IN>


interface surface_flux
!    module procedure surface_flux_0d
    module procedure surface_flux_1d
    module procedure surface_flux_2d
end interface
! </INTERFACE>

!-----------------------------------------------------------------------

character(len=*), parameter :: version = '$Id: surface_flux.F90,v 19.0.6.1 2013/01/24 20:23:30 pjp Exp $'
character(len=*), parameter :: tagname = '$Name:  $'

logical :: do_init = .true.

!jp As grav is no longer a `parameter`, initialisation of these variables
!   now happens in surface_flux_init
real :: d622, d378, hlars, gcp, kappa, d608


! ---- namelist with default values ------------------------------------------
! <NAMELIST NAME="surface_flux_nml">
!   <DATA NAME="no_neg_q"  TYPE="logical"  DEFAULT=".false.">
!    If q_atm_in (specific humidity) is negative (because of numerical truncation),
!    then override with 0.
!   </DATA>
!   <DATA NAME="use_virtual_temp"  TYPE="logical"  DEFAULT=".true.">
!    If true, use virtual potential temp to calculate the stability of the surface layer.
!    if false, use potential temp.
!   </DATA>
!   <DATA NAME="alt_gustiness"  TYPE="logical"  DEFAULT=".false.">
!   An alternative formulation for gustiness calculation.
!   A minimum bound on the wind speed used influx calculations, with the bound
!   equal to gust_const
!   </DATA>
!   <DATA NAME="old_dtaudv"  TYPE="logical"  DEFAULT=".false.">
!   The derivative of surface wind stress w.r.t. the zonal wind and
!   meridional wind are approximated by the same tendency.
!   </DATA>
!   <DATA NAME="use_mixing_ratio"  TYPE="logical"  DEFAULT=".false.">
!   An option to provide capability to run the Manabe Climate form of the
!   surface flux (coded for legacy purposes).
!   </DATA>
!   <DATA NAME="gust_const"  TYPE=""  DEFAULT="1.0">
!    Constant for alternative gustiness calculation.
!   </DATA>
!   <DATA NAME="gust_min"  TYPE=""  DEFAULT="0.0">
!    Minimum gustiness used when alt_gustiness = false.
!   </DATA>
!   <DATA NAME="ncar_ocean_flux"  TYPE="logical"  DEFAULT=".false.">
!    Use NCAR climate model turbulent flux calculation described by
!    Large and Yeager, NCAR Technical Document, 2004
!   </DATA>
!   <DATA NAME="ncar_ocean_flux_orig"  TYPE="logical"  DEFAULT=".false.">
!    Use NCAR climate model turbulent flux calculation described by
!    Large and Yeager, NCAR Technical Document, 2004, using the original
!    GFDL implementation, which contains a bug in the specification of
!    the exchange coefficient for the sensible heat.  This option is available
!    for legacy purposes, and is not recommended for new experiments.
!   </DATA>
!   <DATA NAME="raoult_sat_vap"  TYPE="logical"  DEFAULT=".false.">
!    Reduce saturation vapor pressures to account for seawater salinity.
!   </DATA>
! </NAMELIST>

logical :: no_neg_q              = .false.  ! for backwards compatibility
logical :: use_virtual_temp      = .true.
logical :: alt_gustiness         = .false.
logical :: old_dtaudv            = .false.
logical :: use_mixing_ratio      = .false.
real    :: gust_const            =  1.0
real    :: gust_min              =  0.0
logical :: ncar_ocean_flux       = .false.
logical :: ncar_ocean_flux_orig  = .false. ! for backwards compatibility
logical :: raoult_sat_vap        = .false.
logical :: do_simple             = .false.

real    :: land_humidity_prefactor  =  1.0    ! Default is that land makes no difference to evaporative fluxes
real    :: land_evap_prefactor  =  1.0    ! Default is that land makes no difference to evaporative fluxes

! RUIZHI: from Stephen
logical :: use_actual_surface_temperatures = .true. 
!Always true, apart from when running a dry model, where you can set this to false so that 
! escomp is called with temperatures of 200k always, preventing bad temperature errors.

logical :: use_frierson_mo_drag = .false. 
!Isca by default uses the full monin-obukhov formulae. 
!When true we switch to using the simplified formulae from 10.1175/JAS3753.1 eq 12-14.
! RUIZHI

real    :: flux_heat_gp  =  5.7    ! Default value for Jupiter of 5.7 Wm^-2
real    :: diabatic_acce =  1.0    ! Diabatic acceleration??


namelist /surface_flux_nml/ no_neg_q,             &
                            use_virtual_temp,     &
                            alt_gustiness,        &
                            gust_const,           &
                            gust_min,             &
                            old_dtaudv,           &
                            use_mixing_ratio,     &
                            ncar_ocean_flux,      &
                            ncar_ocean_flux_orig, &
                            raoult_sat_vap,       &
                            do_simple,            &
                            land_humidity_prefactor, & ! Added to make land 'dry', i.e. to decrease the evaporative heat flux in areas of land.
                            land_evap_prefactor, & ! Added to make land 'dry', i.e. to decrease the evaporative heat flux in areas of land.
                            flux_heat_gp,         &    ! prescribed lower boundary heat flux on a giant planet
                            diabatic_acce,       &
                            use_actual_surface_temperatures, & ! RUIZHI
                            use_frierson_mo_drag ! RUIZHI



contains


! ============================================================================
! <SUBROUTINE NAME="surface_flux_1d" INTERFACE="surface_flux">
!  <IN NAME="t_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="q_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="u_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="v_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="p_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="z_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="p_surf" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="t_surf" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="t_ca" TYPE="real, dimension(:)"> </IN>
!  <OUT NAME="q_surf" TYPE="real, dimension(:)"> </OUT>
!  <IN NAME="u_surf" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="v_surf" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="rough_mom" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="rough_heat" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="rough_moist" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="rough_scale" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="gust" TYPE="real, dimension(:)"> </IN>
!  <OUT NAME="flux_t" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="flux_q" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="flux_r" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="flux_u" TYPE="real, dimension(:)"></OUT>
!  <OUT NAME="flux_v" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="cd_m" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="cd_t" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="cd_q" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="w_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="u_star" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="b_star" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="q_star" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dhdt_surf" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dedt_surf" TYPE="real, dimension(:)"></OUT>
!  <OUT NAME="dedq_surf" TYPE="real, dimension(:)"></OUT>
!  <OUT NAME="drdt_surf" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dhdt_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dedq_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dtaudu_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dtaudv_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dt" TYPE="real"> </OUT>
!  <IN NAME="land" TYPE="logical, dimension(:)"> </IN>
!  <IN NAME="seawater" TYPE="logical, dimension(:)"> </IN>
!  <IN NAME="avail" TYPE="logical, dimension(:)"> </IN>


!<PUBLICROUTINE INTERFACE="surface_flux">
subroutine surface_flux_1d (                                           &
     t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     bucket, bucket_depth, max_bucket_depth_land,                      &
     depth_change_lh_1d, depth_change_conv_1d, depth_change_cond_1d,   &
     u_surf,    v_surf,                                                &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
     flux_t, flux_q, flux_r, flux_u, flux_v,                           &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     ex_del_m, ex_del_h, ex_del_q,                                     &
     temp_2m, u_10m, v_10m,                                            &
     q_2m, rh_2m,                                                      &
     dt,        land,      seawater,     avail  )
!</PUBLICROUTINE>
!  slm Mar 28 2002 -- remove agument drag_q since it is just cd_q*wind
! ============================================================================
  ! ---- arguments -----------------------------------------------------------
  logical, intent(in), dimension(:) :: land,  seawater, avail
  logical, intent(in) :: bucket         ! Add bucket model
  real, intent(in),  dimension(:) :: &
       t_atm,     q_atm_in,   u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,  &
       rough_mom, rough_heat, rough_moist,  rough_scale, gust
  real, intent(out), dimension(:) :: &
       flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q,                          & 
       ex_del_m, ex_del_h, ex_del_q,                         &
       temp_2m, u_10m, v_10m,                                &
       q_2m, rh_2m


  real, intent(inout), dimension(:) :: q_surf
  real, intent(inout), dimension(:) :: bucket_depth
  real, intent(inout), dimension(:) :: depth_change_lh_1d
  real, intent(in), dimension(:) :: depth_change_conv_1d, depth_change_cond_1d
  real, intent(in) :: max_bucket_depth_land
  real, intent(in) :: dt

  ! ---- local constants -----------------------------------------------------
  ! temperature increment and its reciprocal value for comp. of derivatives
  real, parameter:: del_temp=0.1, del_temp_inv=1.0/del_temp
  real:: zrefm, zrefh


  ! ---- local vars ----------------------------------------------------------
  real, dimension(size(t_atm(:))) ::                          &
       thv_atm,  th_atm,   tv_atm,    thv_surf,            &
       e_sat,    e_sat1,   q_sat,     q_sat1,    p_ratio,  &
       t_surf0,  t_surf1,  u_dif,     v_dif,               &
       rho_drag, drag_t,    drag_m,   drag_q,    rho,      &
       q_atm,    q_surf0,  dw_atmdu,  dw_atmdv,  w_gust,   &
       e_sat_2m, q_sat_2m

  integer :: i, nbad


  if (do_init) call surface_flux_init

  !---- use local value of surf temp ----
! RUIZHI: from Stephen
  t_surf0 = 200.   !  avoids out-of-bounds in es lookup
if (use_actual_surface_temperatures) then
  where (avail)
     where (land)
        t_surf0 = t_ca
     elsewhere
        t_surf0 = t_surf
     endwhere
  endwhere
! else
!   if ((maxval(q_atm_in) > 0.).and. ( mpp_pe() == mpp_root_pe() )) then
      ! print *, 'surface_flux_mod: Note that you are using dry models.'
!      call error_mesg('surface_flux_mod','Note that you are using dry models.', FATAL)
!   endif
endif

  t_surf1 = t_surf0 + del_temp

  call escomp ( t_surf0, e_sat  )  ! saturation vapor pressure
  call escomp ( t_surf1, e_sat1 )  ! perturbed  vapor pressure

if (.not. use_actual_surface_temperatures) then
  where (avail)
     where (land)
        t_surf0 = t_ca
     elsewhere
        t_surf0 = t_surf
     endwhere
  endwhere
endif
! end RUIZHI

  if(use_mixing_ratio) then
    ! surface mixing ratio at saturation
    q_sat   = d622*e_sat /(p_surf-e_sat )
    q_sat1  = d622*e_sat1/(p_surf-e_sat1)
  elseif(do_simple) then                  !rif:(09/02/09)
    q_sat   = d622*e_sat / p_surf
    q_sat1  = d622*e_sat1/ p_surf
  else
    ! surface specific humidity at saturation
    q_sat   = d622*e_sat /(p_surf-d378*e_sat )
    q_sat1  = d622*e_sat1/(p_surf-d378*e_sat1)
  endif

  ! initilaize surface air humidity according to surface type
  where (land)
!     q_surf0 = q_surf ! land calculates it
     q_surf0 = q_sat ! our simplified land evaporation model does not calculate q_surf, so we specify it as q_sat.
  elsewhere
     q_surf0 = q_sat  ! everything else assumes saturated sfc humidity
  endwhere

  if (raoult_sat_vap) where (seawater) q_surf0 = 0.98 * q_surf0

  ! check for negative atmospheric humidities
  where(avail) q_atm = q_atm_in
  if(no_neg_q) then
     where(avail .and. q_atm_in < 0.0) q_atm = 0.0
  endif

  ! initialize surface air humidity depending on whether surface is dry or wet (bucket empty or not)
  if (bucket) then
  where (bucket_depth <= 0.0)
      q_surf0 = q_atm
  elsewhere
      q_surf0 = q_sat    ! everything else assumes saturated sfc humidity
  end where
  endif

  ! generate information needed by monin_obukhov
  where (avail)
     p_ratio = (p_surf/p_atm)**kappa

     tv_atm  = t_atm  * (1.0 + d608*q_atm)     ! virtual temperature
     th_atm  = t_atm  * p_ratio                ! potential T, using p_surf as refernce
     thv_atm = tv_atm * p_ratio                ! virt. potential T, using p_surf as reference
     thv_surf= t_surf0 * (1.0 + d608*q_surf0 ) ! surface virtual (potential) T
!     thv_surf= t_surf0                        ! surface virtual (potential) T -- just for testing tun off the q_surf

     u_dif = u_surf - u_atm                    ! velocity components relative to surface
     v_dif = v_surf - v_atm
  endwhere

  if(alt_gustiness) then
     do i = 1, size(avail)
        if (.not.avail(i)) cycle
        w_atm(i) = max(sqrt(u_dif(i)**2 + v_dif(i)**2), gust_const)
        ! derivatives of surface wind w.r.t. atm. wind components
        if(w_atm(i) > gust_const) then
           dw_atmdu(i) = u_dif(i)/w_atm(i)
           dw_atmdv(i) = v_dif(i)/w_atm(i)
        else
           dw_atmdu(i) = 0.0
           dw_atmdv(i) = 0.0
        endif
     enddo
  else
     if (gust_min > 0.0) then
       where(avail)
         w_gust = max(gust,gust_min) ! minimum gustiness
       end where
     else
       where(avail)
         w_gust = gust
       end where
     endif

     where(avail)
        w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + w_gust*w_gust)
        ! derivatives of surface wind w.r.t. atm. wind components
        dw_atmdu = u_dif/w_atm
        dw_atmdv = v_dif/w_atm
     endwhere
  endif

  !  monin-obukhov similarity theory
  call mo_drag (thv_atm, thv_surf, z_atm,                  &
       rough_mom, rough_heat, rough_moist, w_atm,          &
       cd_m, cd_t, cd_q, u_star, b_star, avail             )

! added for 10m winds and 2m temperature add mo_profile()

  zrefm = 10. !want winds at 10m
  zrefh = 2.  !want temp and q at 2m

  call mo_profile( zrefm, zrefh, z_atm,   rough_mom, &
       rough_heat, rough_moist,          &
       u_star, b_star, q_star,        &
       ex_del_m, ex_del_h, ex_del_q, avail  )

! adapted from https://github.com/mom-ocean/MOM5/blob/3702ad86f9653f4e315b98613eb824a47d89cf00/src/coupler/flux_exchange.F90#L1932

     !    ------- reference temp -----------
        where (avail) &
           temp_2m = t_surf + (t_atm - t_surf) * ex_del_h !t_ca = canopy temperature, assuming that there is no canopy (no difference between land and ocean), t_ca = t_surf

     !    ------- reference u comp -----------
        where (avail) &
           u_10m = u_atm * ex_del_m ! setting u at surface to 0.

     !    ------- reference v comp -----------
       where (avail) &
           v_10m = v_atm * ex_del_m ! setting v at surface to 0.

! end of low level wind additions


      ! Add 2m q and RH

      ! ------- reference specific humidity -----------
      where (avail) &
           q_2m = q_surf + (q_atm - q_surf) * ex_del_q

      call escomp ( temp_2m, e_sat_2m )

      if(use_mixing_ratio) then
         ! surface mixing ratio at saturation
         q_sat_2m   = d622 * e_sat_2m / (p_surf - e_sat_2m)
      elseif(do_simple) then
         q_sat_2m   = d622 * e_sat_2m / p_surf
      else
         ! surface specific humidity at saturation
         q_sat_2m   = d622 * e_sat_2m / (p_surf - d378*e_sat)
      endif

      ! ------- reference relative humidity -----------
      where (avail) &
         rh_2m = q_2m / q_sat_2m


  ! override with ocean fluxes from NCAR calculation
  if (ncar_ocean_flux .or. ncar_ocean_flux_orig) then
    call  ncar_ocean_fluxes (w_atm, th_atm, t_surf0, q_atm, q_surf0, z_atm, &
                             seawater, cd_m, cd_t, cd_q, u_star, b_star     )
  end if

  where (avail)
     ! scale momentum drag coefficient on orographic roughness
     cd_m = cd_m*(log(z_atm/rough_mom+1)/log(z_atm/rough_scale+1))**2
     ! surface layer drag coefficients
     drag_t = cd_t * w_atm
     drag_q = cd_q * w_atm
     drag_m = cd_m * w_atm

     ! density
     rho = p_atm / (rdgas * tv_atm)

     ! sensible heat flux
     rho_drag = cp_air * drag_t * rho
     flux_t = rho_drag * (t_surf0 - th_atm)  ! flux of sensible heat (W/m**2)
     dhdt_surf =  rho_drag                   ! d(sensible heat flux)/d(surface temperature)
     dhdt_atm  = -rho_drag*p_ratio           ! d(sensible heat flux)/d(atmos temperature)

     ! evaporation
     rho_drag  =  drag_q * rho
  end where

! Add bucket - if bucket is on evaluate fluxes based on moisture availability.
! Note changes to avail statements to allow bucket to be switched on or off
  if (bucket) then
	  where (avail)
	      ! begin LJJ addition
  		where(land)
			where (bucket_depth >= max_bucket_depth_land*0.75)
				flux_q    =  rho_drag * (q_surf0 - q_atm)
			elsewhere
                flux_q    =  bucket_depth/(max_bucket_depth_land*0.75) * rho_drag * (q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))
			end where
		elsewhere
	        flux_q    =  rho_drag * (q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))
		end where

	    depth_change_lh_1d  = flux_q * dt/dens_h2o
	    where (flux_q > 0.0 .and. bucket_depth < depth_change_lh_1d) ! where more evaporation than what's in bucket, empty bucket
	        flux_q = bucket_depth * dens_h2o / dt
	        depth_change_lh_1d = flux_q * dt / dens_h2o
	    end where

	    where (bucket_depth <= 0.0)
	      dedt_surf = 0.
	      dedq_surf = 0.
	      dedq_atm = 0.
	    elsewhere
	      dedq_surf = 0.
	      dedq_atm = -rho_drag ! d(latent heat flux)/d(atmospheric mixing ratio)
		  where(land)
			  where (bucket_depth >= max_bucket_depth_land*0.75)
				  dedt_surf =  rho_drag * (q_sat1 - q_sat) *del_temp_inv
			  elsewhere
      	          dedt_surf =  bucket_depth/(max_bucket_depth_land*0.75) * rho_drag * (q_sat1 - q_sat) *del_temp_inv
			  end where
		  elsewhere
 	          dedt_surf =  rho_drag * (q_sat1 - q_sat) *del_temp_inv
		  end where

	    end where
	  end where
  else

! otherwise revert to simple land model
  where (avail)
     where (land)
!       Simplified land model uses simple prefactor in front of qsurf0. Land is therefore basically the same as sea, but with this prefactor, hence the changes to dedq_surf and dedt_surf also.
        flux_q    =  rho_drag * land_evap_prefactor * (land_humidity_prefactor*q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))
        dedq_surf = 0
        dedt_surf =  rho_drag * land_evap_prefactor * (land_humidity_prefactor*q_sat1 - q_sat) *del_temp_inv
!        dedq_surf = rho_drag 
!        dedt_surf = 0
     elsewhere
        flux_q    =  rho_drag * (q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))
        dedq_surf = 0
        dedt_surf =  rho_drag * (q_sat1 - q_sat) *del_temp_inv
     end where
     dedq_atm  = -rho_drag   ! d(latent heat flux)/d(atmospheric mixing ratio)
   end where
  endif

! end of Add bucket changes

  where (avail)

     q_star = flux_q / (u_star * rho)             ! moisture scale
     ! ask Chris and Steve K if we still want to keep this for diagnostics
     q_surf = q_atm + flux_q / (rho*cd_q*w_atm)   ! surface specific humidity

     ! upward long wave radiation
     flux_r    =   stefan*t_surf**4               ! (W/m**2)
     drdt_surf = 4*stefan*t_surf**3               ! d(upward longwave)/d(surface temperature)

     ! stresses
     rho_drag   = drag_m * rho
     flux_u     = rho_drag * u_dif   ! zonal      component of stress (Nt/m**2)
     flux_v     = rho_drag * v_dif   ! meridional component of stress

  elsewhere
     ! zero-out un-available data in output only fields
     flux_t     = 0.0
     flux_q     = 0.0
     flux_r     = 0.0
     flux_u     = 0.0
     flux_v     = 0.0
     dhdt_surf  = 0.0
     dedt_surf  = 0.0
     dedq_surf  = 0.0
     drdt_surf  = 0.0
     dhdt_atm   = 0.0
     dedq_atm   = 0.0
     u_star     = 0.0
     b_star     = 0.0
     q_star     = 0.0
     q_surf     = 0.0
     w_atm      = 0.0
  end where

  ! calculate d(stress component)/d(atmos wind component)
  dtaudu_atm = 0.0
  dtaudv_atm = 0.0
  if (old_dtaudv) then
     where(avail)
        dtaudv_atm = -rho_drag
        dtaudu_atm = -rho_drag
     endwhere
  else
     where(avail)
        dtaudu_atm = -cd_m*rho*(dw_atmdu*u_dif + w_atm)
        dtaudv_atm = -cd_m*rho*(dw_atmdv*v_dif + w_atm)
     endwhere
  endif

end subroutine surface_flux_1d
! </SUBROUTINE>

!#######################################################################

subroutine surface_flux_0d (                                                 &
     t_atm_0,     q_atm_0,      u_atm_0,     v_atm_0,   p_atm_0, z_atm_0,    &
     p_surf_0,    t_surf_0,     t_ca_0,      q_surf_0,                       &
     u_surf_0,    v_surf_0,                                                  &
     rough_mom_0, rough_heat_0, rough_moist_0, rough_scale_0, gust_0,        &
     flux_t_0,    flux_q_0,     flux_r_0,    flux_u_0,  flux_v_0,            &
     cd_m_0,      cd_t_0,       cd_q_0,                                      &
     w_atm_0,     u_star_0,     b_star_0,     q_star_0,                      &
     dhdt_surf_0, dedt_surf_0,  dedq_surf_0,  drdt_surf_0,                   &
     dhdt_atm_0,  dedq_atm_0,   dtaudu_atm_0, dtaudv_atm_0,                  &
     ex_del_m_0, ex_del_h_0, ex_del_q_0,                                     &
     temp_2m_0, u_10m_0, v_10m_0,                                            &
     q_2m_0, rh_2m_0,                                                        &
     dt,          land_0,       seawater_0,  avail_0  )

  ! ---- arguments -----------------------------------------------------------
  logical, intent(in) :: land_0,  seawater_0, avail_0
  real, intent(in) ::                                                  &
       t_atm_0,     q_atm_0,      u_atm_0,     v_atm_0,                &
       p_atm_0,     z_atm_0,      t_ca_0,                              &
       p_surf_0,    t_surf_0,     u_surf_0,    v_surf_0,               &
       rough_mom_0, rough_heat_0, rough_moist_0, rough_scale_0, gust_0
  real, intent(out) ::                                                 &
       flux_t_0,    flux_q_0,     flux_r_0,    flux_u_0,  flux_v_0,    &
       dhdt_surf_0, dedt_surf_0,  dedq_surf_0, drdt_surf_0,            &
       dhdt_atm_0,  dedq_atm_0,   dtaudu_atm_0,dtaudv_atm_0,           &
       w_atm_0,     u_star_0,     b_star_0,    q_star_0,               &
       cd_m_0,      cd_t_0,       cd_q_0,                              &
       ex_del_m_0, ex_del_h_0, ex_del_q_0,                             &
       temp_2m_0, u_10m_0, v_10m_0,                                    &
       q_2m_0, rh_2m_0
  real, intent(inout) :: q_surf_0
  real, intent(in)    :: dt

  ! ---- local vars ----------------------------------------------------------
  logical, dimension(1) :: land,  seawater, avail
  logical :: bucket
  real, dimension(1) :: &
       t_atm,     q_atm,      u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,             &
       rough_mom, rough_heat, rough_moist,  rough_scale, gust
  real, dimension(1) :: &
       flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q,                          &
       ex_del_m, ex_del_h, ex_del_q,                         &
       temp_2m, u_10m, v_10m,                                &
       q_2m, rh_2m

  real, dimension(1) :: q_surf
  real, dimension(1) :: bucket_depth
  real, dimension(1) :: depth_change_lh_1d
  real, dimension(1) :: depth_change_conv_1d, depth_change_cond_1d
  real :: max_bucket_depth_land

  avail = .true.

  t_atm(1)       = t_atm_0
  q_atm(1)       = q_atm_0
  u_atm(1)       = u_atm_0
  v_atm(1)       = v_atm_0
  p_atm(1)       = p_atm_0
  z_atm(1)       = z_atm_0
  t_ca(1)        = t_ca_0
  p_surf(1)      = p_surf_0
  t_surf(1)      = t_surf_0
  u_surf(1)      = u_surf_0
  v_surf(1)      = v_surf_0
  rough_mom(1)   = rough_mom_0
  rough_heat(1)  = rough_heat_0
  rough_moist(1) = rough_moist_0
  rough_scale(1) = rough_scale_0
  gust(1)        = gust_0
  q_surf(1)      = q_surf_0
  land(1)        = land_0
  seawater(1)    = seawater_0
  avail(1)       = avail_0

  call surface_flux_1d (                                                 &
       t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
       p_surf,    t_surf,     t_ca,      q_surf,                         &
       bucket, bucket_depth, max_bucket_depth_land,                      &
       depth_change_lh_1d, depth_change_conv_1d, depth_change_cond_1d,   &
       u_surf,    v_surf,                                                &
       rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
       flux_t, flux_q, flux_r, flux_u, flux_v,                           &
       cd_m,      cd_t,       cd_q,                                      &
       w_atm,     u_star,     b_star,     q_star,                        &
       dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
       dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
       ex_del_m, ex_del_h, ex_del_q,                                     &
       temp_2m, u_10m, v_10m,                                            &
       q_2m, rh_2m,                                                      &
       dt,        land,      seawater, avail  )

  flux_t_0     = flux_t(1)
  flux_q_0     = flux_q(1)
  flux_r_0     = flux_r(1)
  flux_u_0     = flux_u(1)
  flux_v_0     = flux_v(1)
  dhdt_surf_0  = dhdt_surf(1)
  dedt_surf_0  = dedt_surf(1)
  dedq_surf_0  = dedq_surf(1)
  drdt_surf_0  = drdt_surf(1)
  dhdt_atm_0   = dhdt_atm(1)
  dedq_atm_0   = dedq_atm(1)
  dtaudu_atm_0 = dtaudu_atm(1)
  dtaudv_atm_0 = dtaudv_atm(1)
  w_atm_0      = w_atm(1)
  u_star_0     = u_star(1)
  b_star_0     = b_star(1)
  q_star_0     = q_star(1)
  q_surf_0     = q_surf(1)
  cd_m_0       = cd_m(1)
  cd_t_0       = cd_t(1)
  cd_q_0       = cd_q(1)
  ex_del_m_0   = ex_del_m(1)
  ex_del_h_0   = ex_del_h(1)
  ex_del_q_0   = ex_del_q(1)
  temp_2m_0    = temp_2m(1)
  u_10m_0      = u_10m(1)
  v_10m_0      = v_10m(1)
  q_2m_0       = q_2m(1)
  rh_2m_0      = rh_2m(1)

end subroutine surface_flux_0d

subroutine surface_flux_2d (                                           &
     t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     bucket, bucket_depth, max_bucket_depth_land,                      &
     depth_change_lh,   depth_change_conv,   depth_change_cond,        &
     u_surf,    v_surf,                                                &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
     flux_t,    flux_q,     flux_r,    flux_u,    flux_v,              &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     ex_del_m, ex_del_h, ex_del_q,                                     &
     temp_2m, u_10m, v_10m,                                            &
     q_2m, rh_2m,                                                      &
     dt,        land,       seawater,  avail  )

  ! ---- arguments -----------------------------------------------------------
  logical, intent(in), dimension(:,:) :: land,  seawater, avail
  real, intent(in),  dimension(:,:) :: &
       t_atm,     q_atm_in,   u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,             &
       rough_mom, rough_heat, rough_moist, rough_scale, gust
  real, intent(out), dimension(:,:) :: &
       flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q,                          &
       ex_del_m, ex_del_h, ex_del_q,                         &
       temp_2m, u_10m, v_10m,                                &
       q_2m, rh_2m

  real, intent(inout), dimension(:,:) :: q_surf
  logical, intent(in) :: bucket
  real, intent(inout), dimension(:,:) :: bucket_depth
  real, intent(inout), dimension(:,:) :: depth_change_lh
  real, intent(in), dimension(:,:)    :: depth_change_conv, depth_change_cond
  real, intent(in) :: max_bucket_depth_land
  real, intent(in) :: dt

  ! ---- local vars -----------------------------------------------------------
  integer :: j

  do j = 1, size(t_atm,2)
     call surface_flux_1d (                                           &
          t_atm(:,j),     q_atm_in(:,j),   u_atm(:,j),     v_atm(:,j),     p_atm(:,j),     z_atm(:,j),    &
          p_surf(:,j),    t_surf(:,j),     t_ca(:,j),      q_surf(:,j),                                   &
          bucket, bucket_depth(:,j), max_bucket_depth_land,                                               &
          depth_change_lh(:,j), depth_change_conv(:,j), depth_change_cond(:,j),                           &
          u_surf(:,j),    v_surf(:,j),                                                                    &
          rough_mom(:,j), rough_heat(:,j), rough_moist(:,j), rough_scale(:,j), gust(:,j),                 &
          flux_t(:,j),    flux_q(:,j),     flux_r(:,j),    flux_u(:,j),    flux_v(:,j),                   &
          cd_m(:,j),      cd_t(:,j),       cd_q(:,j),                                                     &
          w_atm(:,j),     u_star(:,j),     b_star(:,j),     q_star(:,j),                                  &
          dhdt_surf(:,j), dedt_surf(:,j),  dedq_surf(:,j),  drdt_surf(:,j),                               &
          dhdt_atm(:,j),  dedq_atm(:,j),   dtaudu_atm(:,j), dtaudv_atm(:,j),                              &
          ex_del_m(:,j), ex_del_h(:,j), ex_del_q(:,j),                                                    &
          temp_2m(:,j), u_10m(:,j), v_10m(:,j),                                                           &
          q_2m(:,j), rh_2m(:,j),                                                                          &
          dt,             land(:,j),       seawater(:,j),  avail(:,j)  )
  end do
end subroutine surface_flux_2d


! ============================================================================
!  Initialization of the surface flux module--reads the nml.
!
subroutine surface_flux_init

! ---- local vars ----------------------------------------------------------
  integer :: unit, ierr, io

  ! read namelist
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, surface_flux_nml, iostat=io)
      ierr = check_nml_error(io,'surface_flux_nml')
#else
  if ( file_exist('input.nml')) then
     unit = open_namelist_file ()
     ierr=1;
     do while (ierr /= 0)
        read  (unit, nml=surface_flux_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'surface_flux_nml')
     enddo
10   call close_file (unit)
  endif
#endif

  ! write version number
  call write_version_number(version, tagname)

  unit = stdlog()
  if ( mpp_pe() == mpp_root_pe() )  write (unit, nml=surface_flux_nml)

  d622   = rdgas/rvgas
  d378   = 1.-d622
  hlars  = hlv/rvgas
  gcp    = grav/cp_air
  kappa  = rdgas/cp_air
  d608   = d378/d622

  ! d608 set to zero  if the use of
  ! virtual temperatures is turned off in namelist
  if(.not. use_virtual_temp) d608 = 0.0

  do_init = .false.

end subroutine surface_flux_init



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Original  code: GFDL.Climate.Model.Info@noaa.gov
! Update Jul2007: GFDL.Climate.Model.Info@noaa.gov (ch and ce exchange coeff bugfix)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes (u_del, t, ts, q, qs, z, avail, &
                              cd, ch, ce, ustar, bstar       )
real   , intent(in)   , dimension(:) :: u_del, t, ts, q, qs, z
logical, intent(in)   , dimension(:) :: avail
real   , intent(inout), dimension(:) :: cd, ch, ce, ustar, bstar

  real :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real :: cd_rt                                ! full drag coefficients @ z
  real :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real :: u, u10, tv, tstar, qstar, z0, xx, stab
  integer, parameter :: n_itts = 2
  integer               i, j

  if(ncar_ocean_flux_orig) then

      do i=1,size(u_del(:))
         if (avail(i)) then
             tv = t(i)*(1+0.608*q(i));
             u = max(u_del(i), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
             u10 = u;                                                ! first guess 10m wind

             cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
             cd_n10_rt = sqrt(cd_n10);
             ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
             stab = 0.5 + sign(0.5,t(i)-ts(i))
             ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c

             cd(i) = cd_n10;                                         ! first guess for exchange coeff's at z
             ch(i) = ch_n10;
             ce(i) = ce_n10;
             do j=1,n_itts                                           ! Monin-Obukhov iteration
                cd_rt = sqrt(cd(i));
                ustar(i) = cd_rt*u;                                   ! L-Y eqn. 7a
                tstar    = (ch(i)/cd_rt)*(t(i)-ts(i));                ! L-Y eqn. 7b
                qstar    = (ce(i)/cd_rt)*(q(i)-qs(i));                ! L-Y eqn. 7c
                bstar(i) = grav*(tstar/tv+qstar/(q(i)+1/0.608));
                zeta     = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
                zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
                x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
                x2 = max(x2, 1.0);                                    ! undocumented NCAR
                x = sqrt(x2);

                if (zeta > 0) then
                    psi_m = -5*zeta;                                    ! L-Y eqn. 8c
                    psi_h = -5*zeta;                                    ! L-Y eqn. 8c
                else
                    psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
                    psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
                end if

                u10 = u/(1+cd_n10_rt*(log(z(i)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9
                cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
                cd_n10_rt = sqrt(cd_n10);
                ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
                stab = 0.5 + sign(0.5,zeta)
                ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
                z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic

                xx = (log(z(i)/10)-psi_m)/vonkarm;
                cd(i) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
                xx = (log(z(i)/10)-psi_h)/vonkarm;
                ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2;                !     10b (this code is wrong)
                ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)**2;                !     10c (this code is wrong)
             end do
         end if
      end do

  else

      do i=1,size(u_del(:))
         if (avail(i)) then
             tv = t(i)*(1+0.608*q(i));
             u = max(u_del(i), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
             u10 = u;                                                ! first guess 10m wind

             cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
             cd_n10_rt = sqrt(cd_n10);
             ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
             stab = 0.5 + sign(0.5,t(i)-ts(i))
             ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c

             cd(i) = cd_n10;                                         ! first guess for exchange coeff's at z
             ch(i) = ch_n10;
             ce(i) = ce_n10;
             do j=1,n_itts                                           ! Monin-Obukhov iteration
                cd_rt = sqrt(cd(i));
                ustar(i) = cd_rt*u;                                   ! L-Y eqn. 7a
                tstar    = (ch(i)/cd_rt)*(t(i)-ts(i));                ! L-Y eqn. 7b
                qstar    = (ce(i)/cd_rt)*(q(i)-qs(i));                ! L-Y eqn. 7c
                bstar(i) = grav*(tstar/tv+qstar/(q(i)+1/0.608));
                zeta     = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
                zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
                x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
                x2 = max(x2, 1.0);                                    ! undocumented NCAR
                x = sqrt(x2);

                if (zeta > 0) then
                    psi_m = -5*zeta;                                    ! L-Y eqn. 8c
                    psi_h = -5*zeta;                                    ! L-Y eqn. 8c
                else
                    psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
                    psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
                end if

                u10 = u/(1+cd_n10_rt*(log(z(i)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9
                cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
                cd_n10_rt = sqrt(cd_n10);
                ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
                stab = 0.5 + sign(0.5,zeta)
                ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
                z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic

                xx = (log(z(i)/10)-psi_m)/vonkarm;
                cd(i) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
                xx = (log(z(i)/10)-psi_h)/vonkarm;
                ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10b (corrected code)
                ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10c (corrected code)
             end do
         end if
      end do

  endif

end subroutine ncar_ocean_fluxes

subroutine gp_surface_flux (dt_tg, p_half, num_levels)

real   , intent(inout), dimension(:,:,:) :: dt_tg
real   , intent(in), dimension(:,:,:) :: p_half
integer   , intent(in) :: num_levels


! add the internal heat flux
dt_tg(:,:,num_levels) = dt_tg(:,:,num_levels)                          &
  + diabatic_acce*grav*flux_heat_gp/(cp_air*(p_half(:,:,num_levels+1)    &
  - p_half(:,:,num_levels)))


end subroutine gp_surface_flux


end module surface_flux_mod

