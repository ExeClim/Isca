
module moisture_convert_mod

!-----------------------------------------------------------------------

use  plev_constants_mod, only: RDGAS, RVGAS, TFREEZE

implicit none
private

public  sphum_to_mixrat, mixrat_to_sphum,  &
            e_to_mixrat,      e_to_sphum,  &
           mixrat_to_rh,     sphum_to_rh,  &
           rh_to_mixrat,     rh_to_sphum

!-----------------------------------------------------------------------

interface sphum_to_mixrat
    module procedure sphum_to_mixrat_0d, sphum_to_mixrat_3d
end interface

interface mixrat_to_sphum
    module procedure mixrat_to_sphum_0d, mixrat_to_sphum_3d
end interface

interface e_to_mixrat
    module procedure e_to_mixrat_0d, e_to_mixrat_3d
end interface

interface e_to_sphum
    module procedure e_to_sphum_0d, e_to_sphum_3d
end interface

interface mixrat_to_rh
    module procedure mixrat_to_rh_0d, mixrat_to_rh_3d
end interface

interface sphum_to_rh
    module procedure sphum_to_rh_0d, sphum_to_rh_3d
end interface

interface rh_to_mixrat
    module procedure rh_to_mixrat_0d, rh_to_mixrat_3d
end interface

interface rh_to_sphum
    module procedure rh_to_sphum_0d, rh_to_sphum_3d
end interface

interface lookup_es
    module procedure lookup_es_0d, lookup_es_3d
end interface

!-----------------------------------------------------------------------

real, parameter :: eps = RDGAS/RVGAS,  one_eps = 1.-eps

!-----------------------------------------------------------------------
!  parameters for es table size and resolution

 integer, parameter :: tcmin = -160  ! minimum temperature (degC) in lookup table
 integer, parameter :: tcmax =  100  ! maximum temperature (degC) in lookup table
 integer, parameter :: esres =  10   ! table resolution (increments per degree)
 integer, parameter :: nsize = (tcmax-tcmin)*esres+1    !  lookup table size
 integer, parameter :: nlim  = nsize-1

 real    :: tmin, tmax          !  lookup table limits in degK
 real    :: dtres, dtinv, teps

 real ::   TABLE(nsize)    !  sat vapor pres (es)
 real ::  DTABLE(nsize)    !  first derivative of es
 real :: D2TABLE(nsize)    ! second derivative of es

 logical :: table_is_initialized = .FALSE.

contains

!#######################################################################

  function sphum_to_mixrat_0d (sphum) result (mixrat)

    real, intent(in) :: sphum
    real             :: mixrat

    mixrat = sphum / (1.0 - sphum)

  end function sphum_to_mixrat_0d

!-----------------------------------------------------------------------

  function sphum_to_mixrat_3d (sphum) result (mixrat)

    real, intent(in) :: sphum(:,:,:)
    real, dimension(size(sphum,1),size(sphum,2),size(sphum,3)) :: mixrat

    mixrat = sphum / (1.0 - sphum)

  end function sphum_to_mixrat_3d

!#######################################################################

  function mixrat_to_sphum_0d (mixrat) result (sphum)

    real, intent(in) :: mixrat
    real             :: sphum

    sphum = mixrat / (1.0 + mixrat)

  end function mixrat_to_sphum_0d

!-----------------------------------------------------------------------

  function mixrat_to_sphum_3d (mixrat) result (sphum)

    real, intent(in) :: mixrat(:,:,:)
    real, dimension(size(mixrat,1),size(mixrat,2),size(mixrat,3)) ::  &
            sphum

    sphum = mixrat / (1.0 + mixrat)

  end function mixrat_to_sphum_3d

!#######################################################################

  function e_to_mixrat_0d (e, pres) result (mixrat)

    real, intent(in) :: e, pres
    real             :: mixrat

    mixrat = eps * e / (pres - e )

  end function e_to_mixrat_0d

!-----------------------------------------------------------------------

  function e_to_mixrat_3d (e, pres) result (mixrat)

    real, intent(in), dimension(:,:,:) :: e, pres
    real, dimension(size(e,1),size(e,2),size(e,3)) :: mixrat

    mixrat = eps * e / (pres - e )

  end function e_to_mixrat_3d

!#######################################################################

  function e_to_sphum_0d (e, pres) result (sphum)

    real, intent(in) :: e, pres
    real             :: sphum

    sphum = eps * e / (pres - one_eps * e )

  end function e_to_sphum_0d

!-----------------------------------------------------------------------

  function e_to_sphum_3d (e, pres) result (sphum)

    real, intent(in), dimension(:,:,:) :: e, pres
    real, dimension(size(e,1),size(e,2),size(e,3)) :: sphum

    sphum = eps * e / (pres - one_eps * e )

  end function e_to_sphum_3d

!#######################################################################

  function mixrat_to_rh_0d (mixrat, temp, pres) result (rh)

    real, intent(in) :: mixrat, temp, pres
    real             :: rh, esat, mixrat_sat

    call lookup_es (temp, esat)
    mixrat_sat = e_to_mixrat (esat, pres)

    rh = mixrat / mixrat_sat

  end function mixrat_to_rh_0d

!-----------------------------------------------------------------------

  function mixrat_to_rh_3d (mixrat, temp, pres) result (rh)

    real, intent(in), dimension(:,:,:) :: mixrat, temp, pres
    real, dimension(size(temp,1),size(temp,2),size(temp,3)) ::  &
            rh, esat, mixrat_sat

    call lookup_es (temp, esat)
    mixrat_sat = e_to_mixrat (esat, pres)

    rh = mixrat / mixrat_sat

  end function mixrat_to_rh_3d

!#######################################################################

  function sphum_to_rh_0d (sphum, temp, pres) result (rh)

    real, intent(in) :: sphum, temp, pres
    real  ::  rh, mixrat

    mixrat = sphum_to_mixrat (sphum)
    rh = mixrat_to_rh (mixrat, temp, pres)

  end function sphum_to_rh_0d

!-----------------------------------------------------------------------

  function sphum_to_rh_3d (sphum, temp, pres) result (rh)

    real, intent(in), dimension(:,:,:) :: sphum, temp, pres
    real, dimension(size(temp,1),size(temp,2),size(temp,3)) ::  &
            rh, mixrat

    mixrat = sphum_to_mixrat (sphum)
    rh = mixrat_to_rh (mixrat, temp, pres)

  end function sphum_to_rh_3d

!#######################################################################

  function rh_to_mixrat_0d (rh, temp, pres) result (mixrat)

    real, intent(in) :: rh, temp, pres
    real             ::  mixrat, esat, mixrat_sat

    call lookup_es (temp, esat)
    mixrat_sat = e_to_mixrat (esat, pres)

    mixrat = mixrat_sat * rh

  end function rh_to_mixrat_0d

!-----------------------------------------------------------------------

  function rh_to_mixrat_3d (rh, temp, pres) result (mixrat)

    real, intent(in), dimension(:,:,:) :: rh, temp, pres
    real, dimension(size(temp,1),size(temp,2),size(temp,3)) ::  &
            mixrat, esat, mixrat_sat

    call lookup_es (temp, esat)
    mixrat_sat = e_to_mixrat (esat, pres)

    mixrat = mixrat_sat * rh

  end function rh_to_mixrat_3d

!#######################################################################

  function rh_to_sphum_0d (rh, temp, pres) result (sphum)

    real, intent(in) :: rh, temp, pres
    real             ::  sphum, mixrat

    mixrat = rh_to_mixrat (rh, temp, pres)
    sphum  = mixrat_to_sphum (mixrat)

  end function rh_to_sphum_0d

!-----------------------------------------------------------------------

  function rh_to_sphum_3d (rh, temp, pres) result (sphum)

    real, intent(in), dimension(:,:,:) :: rh, temp, pres
    real, dimension(size(temp,1),size(temp,2),size(temp,3)) ::  &
            sphum, mixrat

    mixrat = rh_to_mixrat (rh, temp, pres)
    sphum  = mixrat_to_sphum (mixrat)

  end function rh_to_sphum_3d

!#######################################################################
!           saturation vapor pressure lookup table
!#######################################################################

  subroutine lookup_es_0d (temp, esat)
    real, intent(in)  :: temp
    real, intent(out) :: esat
    real    :: tmp, del
    integer :: ind

    call sat_vapor_pres_init

    tmp = temp-tmin
    ind = int(dtinv*(tmp+teps))
    del = tmp-dtres*real(ind)
    esat = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))

    if (ind < 0 .or. ind > nlim) then
       print *, 'ERROR: saturation vapor pressure error'
       stop 111
    endif

  end subroutine lookup_es_0d

!-----------------------------------------------------------------------

  subroutine lookup_es_3d (temp, esat)
    real, intent(in)  :: temp(:,:,:)
    real, intent(out) :: esat(:,:,:)
    real    :: tmp, del
    integer :: ind
    integer :: i, j, k

    call sat_vapor_pres_init

    do  k = 1, size(temp,3)
    do  j = 1, size(temp,2)
    do  i = 1, size(temp,1)
       tmp = temp(i,j,k)-tmin
       ind = int(dtinv*(tmp+teps))
       del = tmp-dtres*real(ind)
       esat(i,j,k) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
       if (ind < 0 .or. ind > nlim) then
          print *, 'ERROR: saturation vapor pressure error'
          stop 111
       endif
    enddo
    enddo
    enddo

  end subroutine lookup_es_3d

!-----------------------------------------------------------------------
!                Computation of the es values
!
!   Saturation vapor pressure (es) values are computed from
!   equations in the Smithsonian meteorological tables page 350.
!   For temperatures < 0C, sat vapor pres is computed over ice.
!   For temperatures > -20C, sat vapor pres is computed over water.
!   Between -20C and 0C the returned value is blended (over water
!   and over ice).  All sat vapor pres values are returned in pascals.
!
!   Reference:  Smithsonian meteorological tables, page 350.

function compute_es (tem) result (es)
 real, intent(in) :: tem(:)
 real :: es(size(tem,1))

 real, parameter :: TBASW = TFREEZE+100.
 real, parameter :: TBASI = TFREEZE
 real, parameter :: ESBASW = 101324.60
 real, parameter :: ESBASI =    610.71

 real    :: x, esice, esh2o
 integer :: i

   do i = 1, size(tem(:))

!  compute es over ice 

     if (tem(i) < TBASI) then
         x = -9.09718*(TBASI/tem(i)-1.0) - 3.56654*log10(TBASI/tem(i)) &
             +0.876793*(1.0-tem(i)/TBASI) + log10(ESBASI)
         esice =10.**(x)
     else
         esice = 0.
     endif

!  compute es over water greater than -20 c.
!  values over 100 c may not be valid
!  see smithsonian meteorological tables page 350.

     if (tem(i) > -20.+TBASI) then
         x = -7.90298*(TBASW/tem(i)-1) + 5.02808*log10(TBASW/tem(i)) &
             -1.3816e-07*(10**((1-tem(i)/TBASW)*11.344)-1)        &
             +8.1328e-03*(10**((TBASW/tem(i)-1)*(-3.49149))-1)    &
             +log10(ESBASW)
         esh2o = 10.**(x)
     else
         esh2o = 0.
     endif

!  derive blended es over ice and supercooled water between -20c and 0c

     if (tem(i) <= -20.+TBASI) then
         es(i) = esice
     else if (tem(i) >= TBASI) then
         es(i) = esh2o
     else
         es(i) = 0.05*((TBASI-tem(i))*esice + (tem(i)-TBASI+20.)*esh2o)
     endif

   enddo

end function compute_es

!-----------------------------------------------------------------------

subroutine sat_vapor_pres_init

  real    :: tem(3), es(3), hdtinv
  integer :: i

! increment used to generate derivative table
  real, parameter :: tinrc = .01           
  real, parameter :: tfact = 1./(2.*tinrc)

! return silently if this routine has already been called
  if (table_is_initialized) return

! global variables
      tmin = real(tcmin)+TFREEZE   ! minimum valid temp in table
      tmax = real(tcmax)+TFREEZE   ! maximum valid temp in table
      dtinv = real(esres)
      dtres = 1./dtinv
      teps = 1./real(2*esres)
! local variables
      hdtinv = dtinv*0.5! global variables
      tmin = real(tcmin)+TFREEZE   ! minimum valid temp in table
      tmax = real(tcmax)+TFREEZE   ! maximum valid temp in table
      dtinv = real(esres)
      dtres = 1./dtinv
      teps = 1./real(2*esres)
! local variables
      hdtinv = dtinv*0.5

! compute es tables from tcmin to tcmax
! estimate es derivative with small +/- difference

      do i = 1, nsize
         tem(1) = tmin + dtres*real(i-1)
         tem(2) = tem(1)-tinrc
         tem(3) = tem(1)+tinrc
         es = compute_es (tem)
          TABLE(i) = es(1)
         DTABLE(i) = (es(3)-es(2))*tfact
      enddo

! compute one-half second derivative using centered differences
! differencing des values in the table

      do i = 2, nsize-1
         D2TABLE(i) = 0.25*dtinv*(DTABLE(i+1)-DTABLE(i-1))
      enddo
    ! one-sided derivatives at boundaries
         D2TABLE(1)     = 0.50*dtinv*(DTABLE(2)    -DTABLE(1))
         D2TABLE(nsize) = 0.50*dtinv*(DTABLE(nsize)-DTABLE(nsize-1))

   table_is_initialized = .true.

end subroutine sat_vapor_pres_init

end module moisture_convert_mod

