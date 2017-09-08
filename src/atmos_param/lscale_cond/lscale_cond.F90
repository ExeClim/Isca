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

module lscale_cond_mod

!-----------------------------------------------------------------------
use            mpp_mod, only:  input_nml_file
use            fms_mod, only:  file_exist, error_mesg, open_namelist_file,  &
                               check_nml_error, mpp_pe, mpp_root_pe, FATAL,  &
                               close_file, write_version_number, stdlog
use sat_vapor_pres_mod, only:  compute_qs
use      constants_mod, only:  HLv,HLs,Cp_Air,Grav,rdgas,rvgas

implicit none
private
!-----------------------------------------------------------------------
!  ---- public interfaces ----

   public  lscale_cond, lscale_cond_init, lscale_cond_end

!-----------------------------------------------------------------------
!   ---- version number ----

 character(len=128) :: version = '$Id: lscale_cond.F90,v 19.0 2012/01/06 20:10:05 fms Exp $'
 character(len=128) :: tagname = '$Name:  $'
 logical            :: module_is_initialized=.false.

!-----------------------------------------------------------------------
!   ---- local/private data ----

    real :: d622 = 0.
    real :: d378 = 0.


!-----------------------------------------------------------------------
!   --- namelist ----

real    :: hc=1.00
logical :: do_evap=.false.
logical :: do_simple =.false.

namelist /lscale_cond_nml/  hc, do_evap, do_simple


!-----------------------------------------------------------------------
!           description of namelist variables
!
!  hc        =  relative humidity at which large scale condensation
!               occurs, where 0 <= hc <= 1 (default: hc=1.)
!
!  do_evap   =  flag for the re-evaporation of moisture in
!               sub-saturated layers below, if do_evap=.true. then
!               re-evaporation is performed (default: do_evap=.false.)
!
!-----------------------------------------------------------------------

contains

!#######################################################################

   subroutine lscale_cond (tin, qin, pfull, phalf, coldT, &
                           rain, snow, tdel, qdel, mask, conv)

!-----------------------------------------------------------------------
!
!                      large scale condensation
!
!-----------------------------------------------------------------------
!
!   input:  tin      temperature at full model levels
!           qin      specific humidity of water vapor at full
!                      model levels
!           pfull    pressure at full model levels
!           phalf    pressure at half (interface) model levels
!           coldT    should precipitation be snow at this point?
!   optional:
!           mask     optional mask (0 or 1.) 
!           conv     logical flag; if true then no large-scale
!                       adjustment is performed at that grid-point or
!                       model level
!
!  output:  rain     liquid precipitation (kg/m2)
!           snow     frozen precipitation (kg/m2)
!           tdel     temperature tendency at full model levels
!           qdel     specific humidity tendency (of water vapor) at
!                      full model levels
!
!-----------------------------------------------------------------------
!--------------------- interface arguments -----------------------------

   real   , intent(in) , dimension(:,:,:) :: tin, qin, pfull, phalf
   logical   , intent(in) , dimension(:,:):: coldT
   real   , intent(out), dimension(:,:)   :: rain,snow
   real   , intent(out), dimension(:,:,:) :: tdel, qdel
   real   , intent(in) , dimension(:,:,:), optional :: mask
   logical, intent(in) , dimension(:,:,:), optional :: conv
!-----------------------------------------------------------------------
!---------------------- local data -------------------------------------

logical,dimension(size(tin,1),size(tin,2),size(tin,3)) :: do_adjust
   real,dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &
                             qsat, dqsat, pmass
   real,dimension(size(tin,1),size(tin,2))             :: hlcp, precip
integer  k, kx
!-----------------------------------------------------------------------
!     computation of precipitation by condensation processes
!-----------------------------------------------------------------------

      if (.not. module_is_initialized) call error_mesg ('lscale_cond',  &
                         'lscale_cond_init has not been called.', FATAL)

      kx=size(tin,3)

!----- compute proper latent heat --------------------------------------
      if(do_simple) then
             hlcp = HLv/Cp_Air
      else
        WHERE (coldT)
             hlcp = HLs/Cp_Air
        ELSEWHERE
             hlcp = HLv/Cp_Air
        END WHERE
      endif

!--- saturation specific humidity (qsat) and deriv wrt temp (dqsat) ---

     call compute_qs (tin, pfull,qsat, hc = hc, dqsdT=dqsat) 

!--------- do adjustment where greater than saturated value ------------

   if (present(conv)) then
!!!!  do_adjust(:,:,:)=(.not.conv(:,:,:) .and. qin(:,:,:) > qsat(:,:,:))
      do_adjust(:,:,:)=(.not.conv(:,:,:) .and.   &
                         (qin(:,:,:) - qsat(:,:,:))*qsat(:,:,:) > 0.0)
   else
!!!!  do_adjust(:,:,:)=(qin(:,:,:) > qsat(:,:,:))
      do_adjust(:,:,:)=( (qin(:,:,:) - qsat(:,:,:))*qsat(:,:,:) > 0.0)
   endif

   if (present(mask)) then
      do_adjust(:,:,:)=do_adjust(:,:,:) .and. (mask(:,:,:) > 0.5)
   end if

!----------- compute adjustments to temp and spec humidity -------------
   do k = 1,kx
   where (do_adjust(:,:,k))
      qdel(:,:,k)=(qsat(:,:,k)-qin(:,:,k))/(1.0+hlcp(:,:)*dqsat(:,:,k))
      tdel(:,:,k)=-hlcp(:,:)*qdel(:,:,k)
   elsewhere
      qdel(:,:,k)=0.0
      tdel(:,:,k)=0.0
   endwhere
   end do
!------------ pressure mass of each layer ------------------------------

   do k=1,kx
      pmass(:,:,k)=(phalf(:,:,k+1)-phalf(:,:,k))/Grav
   enddo

!------------ re-evaporation of precipitation in dry layer below -------

   if (do_evap) then
      if (present(mask)) then
         call precip_evap (pmass,tin,qin,qsat,dqsat,hlcp,tdel,qdel,mask)
      else
         call precip_evap (pmass,tin,qin,qsat,dqsat,hlcp,tdel,qdel)
      endif
   endif

!------------ integrate precip -----------------------------------------

      precip(:,:)=0.0
   do k=1,kx
      precip(:,:)=precip(:,:)-pmass(:,:,k)*qdel(:,:,k)
   enddo
      precip(:,:)=max(precip(:,:),0.0)

   !assign precip to snow or rain
   if(do_simple) then !no snow!
        rain = precip
        snow = 0.
   else
     WHERE (coldT)
        snow = precip
        rain = 0.
     ELSEWHERE
        rain = precip
        snow = 0.
     END WHERE
   endif

!-----------------------------------------------------------------------

   end subroutine lscale_cond

!#######################################################################

subroutine precip_evap (pmass, tin, qin, qsat, dqsat, hlcp, &
                        tdel, qdel, mask)

!-----------------------------------------------------------------------
!        performs re-evaporation of falling precipitation
!-----------------------------------------------------------------------
   real, intent(in),    dimension(:,:,:) :: pmass, tin, qin, qsat, dqsat
   real, intent(in),    dimension(:,:)   :: hlcp
   real, intent(inout), dimension(:,:,:) :: tdel, qdel
   real, intent(in), dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
   real, dimension(size(tin,1),size(tin,2)) :: exq, def

   integer  k
!-----------------------------------------------------------------------
    exq(:,:)=0.0

    do k=1,size(tin,3)

        where (qdel(:,:,k) < 0.0)  exq(:,:) = exq(:,:) -  &
                                               qdel(:,:,k)*pmass(:,:,k)

        if (present(mask)) exq(:,:) = exq(:,:)*mask(:,:,k)

!  ---- evaporate precip where needed ------

        where ( (qdel(:,:,k) >= 0.0) .and. (exq(:,:) > 0.0) )
            exq(:,:) = exq(:,:) / pmass(:,:,k)
            def(:,:) = (qsat(:,:,k)-qin(:,:,k))/(1.+hlcp(:,:)*dqsat(:,:,k))
            def(:,:) = min(max(def(:,:),0.0),exq(:,:))
            qdel(:,:,k) = qdel(:,:,k) + def(:,:)
            tdel(:,:,k) = tdel(:,:,k) - def(:,:)*hlcp(:,:)
            exq(:,:) = (exq(:,:)-def(:,:))*pmass(:,:,k)
        endwhere

    enddo

!-----------------------------------------------------------------------

   end subroutine precip_evap

!#######################################################################

   subroutine lscale_cond_init ()

!-----------------------------------------------------------------------
!
!        initialization for large scale condensation
!
!-----------------------------------------------------------------------

  integer  unit,io,ierr, logunit

!----------- read namelist ---------------------------------------------

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=lscale_cond_nml, iostat=io)
      ierr = check_nml_error(io,"lscale_cond_nml")
#else
      if (file_exist('input.nml')) then
         unit = open_namelist_file ()
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=lscale_cond_nml, iostat=io, end=10)
            ierr = check_nml_error (io,'lscale_cond_nml')
         enddo
  10     call close_file (unit)
      endif
#endif

!---------- output namelist --------------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           logunit = stdlog()
           write (logunit,nml=lscale_cond_nml)
      endif

   !s initialise here as rdgas no longer a parameter
   d622 = rdgas/rvgas
   d378 = 1.-d622

      module_is_initialized=.true.

   end subroutine lscale_cond_init

!#######################################################################
   subroutine lscale_cond_end

      module_is_initialized=.false.

!---------------------------------------------------------------------

   end subroutine lscale_cond_end

!#######################################################################

end module lscale_cond_mod

