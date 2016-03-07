
 module sat_vapor_pres_k_mod

! This module is what I (pjp) think a kernel should be.
! There have been many proposals as to what a kernel should look like.
! If fact, so many different ideas have been expressed that the lack
! of agreement has greatly hampered progress.
! The only way to move forward is to limit the requirments for a kernel
! to only what is widely agreeded upon.
! I believe that there are only two things widely agreeded upon.

! 1) A kernel should be independent of the rest of FMS so that it can
!    easily be ported into another programming system.
!    This requires that a kernel does not access anything by use association.
!    The one exception is this kernel, because it is not practical for physics
!    modules to avoid using a module that computes the saturation vapor
!    pressure of water vapor.

! 2) For the sake of thread safety, module globals should be written only at initialization.
!    In this case, the module globals are the tables and a handful of scalars.

! 3) A kernel should not read from an external file.

! One of the things that was not widely agreeded upon is that a kernel should
! not be a fortran module. This complicates things greatly for questionable
! benefit and could be done as a second step anyway, if necessary.

 implicit none
 private

 character(len=128), parameter :: version = '$Id: sat_vapor_pres_k.F90,v 18.0 2010/03/02 23:58:26 fms Exp $'
 character(len=128), parameter :: tagname = '$Name: testing $'

 public :: sat_vapor_pres_init_k
 public :: lookup_es_k
 public :: lookup_des_k
 public :: lookup_es_des_k
 public :: lookup_es2_k
 public :: lookup_des2_k
 public :: lookup_es2_des2_k
 public :: lookup_es3_k
 public :: lookup_des3_k
 public :: lookup_es3_des3_k
 public :: compute_qs_k
 public :: compute_mrs_k

 interface lookup_es_k
   module procedure lookup_es_k_0d
   module procedure lookup_es_k_1d
   module procedure lookup_es_k_2d
   module procedure lookup_es_k_3d
 end interface

 interface lookup_des_k
   module procedure lookup_des_k_0d
   module procedure lookup_des_k_1d
   module procedure lookup_des_k_2d
   module procedure lookup_des_k_3d
 end interface

 interface lookup_es_des_k
   module procedure lookup_es_des_k_0d
   module procedure lookup_es_des_k_1d
   module procedure lookup_es_des_k_2d
   module procedure lookup_es_des_k_3d
 end interface

 interface lookup_es2_k
   module procedure lookup_es2_k_0d
   module procedure lookup_es2_k_1d
   module procedure lookup_es2_k_2d
   module procedure lookup_es2_k_3d
 end interface

 interface lookup_des2_k
   module procedure lookup_des2_k_0d
   module procedure lookup_des2_k_1d
   module procedure lookup_des2_k_2d
   module procedure lookup_des2_k_3d
 end interface

 interface lookup_es2_des2_k
   module procedure lookup_es2_des2_k_0d
   module procedure lookup_es2_des2_k_1d
   module procedure lookup_es2_des2_k_2d
   module procedure lookup_es2_des2_k_3d
 end interface

 interface lookup_es3_k
   module procedure lookup_es3_k_0d
   module procedure lookup_es3_k_1d
   module procedure lookup_es3_k_2d
   module procedure lookup_es3_k_3d
 end interface

 interface lookup_des3_k
   module procedure lookup_des3_k_0d
   module procedure lookup_des3_k_1d
   module procedure lookup_des3_k_2d
   module procedure lookup_des3_k_3d
 end interface

 interface lookup_es3_des3_k
   module procedure lookup_es3_des3_k_0d
   module procedure lookup_es3_des3_k_1d
   module procedure lookup_es3_des3_k_2d
   module procedure lookup_es3_des3_k_3d
 end interface

 interface compute_qs_k
   module procedure compute_qs_k_0d
   module procedure compute_qs_k_1d
   module procedure compute_qs_k_2d
   module procedure compute_qs_k_3d
 end interface

 interface compute_mrs_k
   module procedure compute_mrs_k_0d
   module procedure compute_mrs_k_1d
   module procedure compute_mrs_k_2d
   module procedure compute_mrs_k_3d
 end interface

 real :: dtres, tepsl, tminl, dtinvl
 integer :: table_siz
 real, dimension(:), allocatable :: TABLE   !  sat vapor pres (es)
 real, dimension(:), allocatable :: DTABLE  !  first derivative of es
 real, dimension(:), allocatable :: D2TABLE ! second derivative of es
 real, dimension(:), allocatable :: TABLE2  !  sat vapor pres (es)
 real, dimension(:), allocatable :: DTABLE2 !  first derivative of es
 real, dimension(:), allocatable :: D2TABLE2 ! second derivative of es
 real, dimension(:), allocatable :: TABLE3  !  sat vapor pres (es)
 real, dimension(:), allocatable :: DTABLE3 !  first derivative of es
 real, dimension(:), allocatable :: D2TABLE3 ! second derivative of es

 logical  :: use_exact_qs
 logical  :: module_is_initialized = .false.

 contains

 subroutine sat_vapor_pres_init_k(table_size, tcmin, tcmax, TFREEZE, HLV, RVGAS, ES0, err_msg, &
                                  use_exact_qs_input, do_simple,  &
                                  construct_table_wrt_liq, &
                                  construct_table_wrt_liq_and_ice, &
                                  teps, tmin, dtinv)

! This routine has been generalized to return tables for any temperature range and resolution

 integer, intent(in) :: table_size
 real, intent(in) :: tcmin ! TABLE(1)          = sat vapor pressure at temperature tcmin (deg C)
 real, intent(in) :: tcmax ! TABLE(table_size) = sat vapor pressure at temperature tcmax (deg C)
 real, intent(in) :: TFREEZE, HLV, RVGAS, ES0
 logical, intent(in)  :: use_exact_qs_input, do_simple
 logical, intent(in)  :: construct_table_wrt_liq
 logical, intent(in)  :: construct_table_wrt_liq_and_ice
 character(len=*), intent(out) :: err_msg
 real, intent(out), optional :: teps, tmin, dtinv

! increment used to generate derivative table
  real, dimension(3) :: tem(3), es(3)
  real :: hdtinv, tinrc, tfact
  integer :: i

      err_msg = ''

      if (module_is_initialized) return

      if(allocated(TABLE) .or. allocated(DTABLE) .or. allocated(D2TABLE)) then
        err_msg = 'Attempt to allocate sat vapor pressure tables when already allocated'
        return
      else
        allocate(TABLE(table_size), DTABLE(table_size), D2TABLE(table_size))
      endif
      
   if (construct_table_wrt_liq) then
      if(allocated(TABLE2) .or. allocated(DTABLE2) .or. allocated(D2TABLE2)) then
        err_msg = 'Attempt to allocate sat vapor pressure table2s when already allocated'
        return
      else
        allocate(TABLE2(table_size), DTABLE2(table_size), D2TABLE2(table_size))
      endif
   endif

   if (construct_table_wrt_liq_and_ice) then
      if(allocated(TABLE3) .or. allocated(DTABLE3) .or. allocated(D2TABLE3)) then
        err_msg = 'Attempt to allocate sat vapor pressure table2s when already allocated'
        return
      else
        allocate(TABLE3(table_size), DTABLE3(table_size), D2TABLE3(table_size))
      endif
   endif

      table_siz = table_size
      dtres = (tcmax - tcmin)/(table_size-1)
      tminl = real(tcmin)+TFREEZE  ! minimum valid temp in table
      dtinvl = 1./dtres
      tepsl = .5*dtres
      tinrc = .1*dtres
      if(present(teps )) teps =tepsl
      if(present(tmin )) tmin =tminl
      if(present(dtinv)) dtinv=dtinvl

! To be able to compute tables for any temperature range and resolution,
! and at the same time exactly reproduce answers from memphis revision,
! it is necessary to compute ftact differently than it is in memphis.
      tfact = 5*dtinvl

      hdtinv = dtinvl*0.5

! compute es tables from tcmin to tcmax
! estimate es derivative with small +/- difference

      if (do_simple) then

        do i = 1, table_size
          tem(1) = tminl + dtres*real(i-1)
          TABLE(i) = ES0*610.78*exp(-hlv/rvgas*(1./tem(1) - 1./tfreeze)) 
          DTABLE(i) = hlv*TABLE(i)/rvgas/tem(1)**2.
        enddo

      else

        do i = 1, table_size
          tem(1) = tminl + dtres*real(i-1)
          tem(2) = tem(1)-tinrc
          tem(3) = tem(1)+tinrc
          es = compute_es_k (tem, TFREEZE)
          TABLE(i) = es(1)
          DTABLE(i) = (es(3)-es(2))*tfact
        enddo

      endif !if (do_simple)

! compute one-half second derivative using centered differences
! differencing des values in the table

      do i = 2, table_size-1
         D2TABLE(i) = 0.25*dtinvl*(DTABLE(i+1)-DTABLE(i-1))
      enddo
    ! one-sided derivatives at boundaries

         D2TABLE(1) = 0.50*dtinvl*(DTABLE(2)-DTABLE(1))

         D2TABLE(table_size) = 0.50*dtinvl*&
              (DTABLE(table_size)-DTABLE(table_size-1))
      
   if (construct_table_wrt_liq) then
! compute es tables from tcmin to tcmax
! estimate es derivative with small +/- difference
 
      do i = 1, table_size
        tem(1) = tminl + dtres*real(i-1)
        tem(2) = tem(1)-tinrc
        tem(3) = tem(1)+tinrc
!   pass in flag to force all values to be wrt liquid
        es = compute_es_liq_k (tem, TFREEZE)
        TABLE2(i) = es(1)
        DTABLE2(i) = (es(3)-es(2))*tfact
      enddo
 
! compute one-half second derivative using centered differences
! differencing des values in the table

     do i = 2, table_size-1
       D2TABLE2(i) = 0.25*dtinvl*(DTABLE2(i+1)-DTABLE2(i-1))
     enddo
! one-sided derivatives at boundaries

     D2TABLE2(1) = 0.50*dtinvl*(DTABLE2(2)-DTABLE2(1))

     D2TABLE2(table_size) = 0.50*dtinvl*&
          (DTABLE2(table_size)-DTABLE2(table_size-1))
   endif


   if (construct_table_wrt_liq_and_ice) then
! compute es tables from tcmin to tcmax
! estimate es derivative with small +/- difference
 
      do i = 1, table_size
        tem(1) = tminl + dtres*real(i-1)
        tem(2) = tem(1)-tinrc
        tem(3) = tem(1)+tinrc
!   pass in flag to force all values to be wrt liquid
        es = compute_es_liq_ice_k (tem, TFREEZE)
        TABLE3(i) = es(1)
        DTABLE3(i) = (es(3)-es(2))*tfact
      enddo
 
! compute one-half second derivative using centered differences
! differencing des values in the table

     do i = 2, table_size-1
       D2TABLE3(i) = 0.25*dtinvl*(DTABLE3(i+1)-DTABLE3(i-1))
     enddo
! one-sided derivatives at boundaries

     D2TABLE3(1) = 0.50*dtinvl*(DTABLE3(2)-DTABLE3(1))

     D2TABLE3(table_size) = 0.50*dtinvl*&
          (DTABLE3(table_size)-DTABLE3(table_size-1))
   endif

      use_exact_qs = use_exact_qs_input
      module_is_initialized = .true.

 end subroutine sat_vapor_pres_init_k

!#######################################################################

 function compute_es_k(tem, TFREEZE) result (es)
 real, intent(in) :: tem(:), TFREEZE
 real :: es(size(tem,1))
         
 real    :: x, esice, esh2o, TBASW, TBASI
 integer :: i
 real, parameter :: ESBASW = 101324.60
 real, parameter :: ESBASI =    610.71

   TBASW = TFREEZE+100.
   TBASI = TFREEZE

   do i = 1, size(tem)

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

 end function compute_es_k

!#######################################################################

 function compute_es_liq_k(tem, TFREEZE) result (es)
 real, intent(in) :: tem(:), TFREEZE
 real :: es(size(tem,1))
         
 real    :: x, esh2o, TBASW
 integer :: i
 real, parameter :: ESBASW = 101324.60

   TBASW = TFREEZE+100.

   do i = 1, size(tem)


!  compute es over water for all temps.
!  values over 100 c may not be valid
!  see smithsonian meteorological tables page 350.

         x = -7.90298*(TBASW/tem(i)-1) + 5.02808*log10(TBASW/tem(i)) &
             -1.3816e-07*(10**((1-tem(i)/TBASW)*11.344)-1)        &
             +8.1328e-03*(10**((TBASW/tem(i)-1)*(-3.49149))-1)    &
             +log10(ESBASW)
         esh2o = 10.**(x)


         es(i) = esh2o

   enddo

 end function compute_es_liq_k

!#######################################################################

 function compute_es_liq_ice_k(tem, TFREEZE) result (es)
 real, intent(in) :: tem(:), TFREEZE
 real :: es(size(tem,1))
         
 real    :: x, TBASW, TBASI
 integer :: i
 real, parameter :: ESBASW = 101324.60
 real, parameter :: ESBASI =    610.71

   TBASW = TFREEZE+100.
   TBASI = TFREEZE

   do i = 1, size(tem)

     if (tem(i) < TBASI) then

!  compute es over ice 

         x = -9.09718*(TBASI/tem(i)-1.0) - 3.56654*log10(TBASI/tem(i)) &
             +0.876793*(1.0-tem(i)/TBASI) + log10(ESBASI)
         es(i) =10.**(x)
     else

!  compute es over water 
!  values over 100 c may not be valid
!  see smithsonian meteorological tables page 350.

         x = -7.90298*(TBASW/tem(i)-1) + 5.02808*log10(TBASW/tem(i)) &
             -1.3816e-07*(10**((1-tem(i)/TBASW)*11.344)-1)        &
             +8.1328e-03*(10**((TBASW/tem(i)-1)*(-3.49149))-1)    &
             +log10(ESBASW)
         es(i) = 10.**(x)
     endif

   enddo

 end function compute_es_liq_ice_k

!#######################################################################

 subroutine compute_qs_k_3d (temp, press, eps, zvir, qs, nbad, q, hc, &
                          dqsdT, esat, es_over_liq, es_over_liq_and_ice)

 real, intent(in),  dimension(:,:,:)           :: temp, press   
 real, intent(in)                              :: eps, zvir
 real, intent(out), dimension(:,:,:)           :: qs   
 integer, intent(out)                          :: nbad
 real, intent(in),  dimension(:,:,:), optional :: q
 real, intent(in),                    optional :: hc
 real, intent(out), dimension(:,:,:), optional :: dqsdT, esat
 logical,intent(in),                  optional :: es_over_liq
 logical,intent(in),                  optional :: es_over_liq_and_ice

 real, dimension(size(temp,1), size(temp,2), size(temp,3)) ::   &
                                                  esloc, desat, denom
 integer :: i, j, k
 real    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0
   endif
 if (present(es_over_liq)) then
   if (present (dqsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dqsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dqsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif 
   if (nbad == 0) then
     if (present (q) .and. use_exact_qs) then
       qs = (1.0 + zvir*q)*eps*esloc/press
       if (present (dqsdT)) then
         dqsdT = (1.0 + zvir*q)*eps*desat/press
       endif
     else  ! (present(q))
       denom = press - (1.0 - eps)*esloc
       do k=1,size(qs,3)
         do j=1,size(qs,2)
           do i=1,size(qs,1)
             if (denom(i,j,k) > 0.0) then
               qs(i,j,k) = eps*esloc(i,j,k)/denom(i,j,k)         
             else
               qs(i,j,k) = eps
             endif
           end do
         end do
       end do
       if (present (dqsdT)) then
         dqsdT = eps*press*desat/denom**2
       endif
     endif ! (present(q))
   else ! (nbad = 0)
     qs = -999.
     if (present (dqsdT)) then
       dqsdT = -999.
     endif
     if (present (esat)) then
       esat = -999.
     endif 
   endif ! (nbad = 0)

     
 end subroutine compute_qs_k_3d

!#######################################################################

 subroutine compute_qs_k_2d (temp, press, eps, zvir, qs, nbad, q, hc, &
                          dqsdT, esat, es_over_liq, es_over_liq_and_ice)

 real, intent(in),  dimension(:,:)           :: temp, press   
 real, intent(in)                            :: eps, zvir
 real, intent(out), dimension(:,:)           :: qs   
 integer, intent(out)                        :: nbad
 real, intent(in),  dimension(:,:), optional :: q
 real, intent(in),                  optional :: hc
 real, intent(out), dimension(:,:), optional :: dqsdT, esat
 logical,intent(in),                optional :: es_over_liq
 logical,intent(in),                  optional :: es_over_liq_and_ice

 real, dimension(size(temp,1), size(temp,2)) :: esloc, desat, denom
 integer :: i, j
 real    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0
   endif

 if (present(es_over_liq)) then
   if (present (dqsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dqsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dqsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif 
   if (nbad == 0) then
     if (present (q) .and. use_exact_qs) then
       qs = (1.0 + zvir*q)*eps*esloc/press
       if (present (dqsdT)) then
         dqsdT = (1.0 + zvir*q)*eps*desat/press
       endif
     else  ! (present(q))
       denom = press - (1.0 - eps)*esloc
      do j=1,size(qs,2)
        do i=1,size(qs,1)
          if (denom(i,j) > 0.0) then
            qs(i,j) = eps*esloc(i,j)/denom(i,j)
          else
            qs(i,j) = eps
          endif
        end do
      end do
      if (present (dqsdT)) then
        dqsdT = eps*press*desat/denom**2
      endif
    endif ! (present(q))
   else ! (nbad = 0)
     qs = -999.
     if (present (dqsdT)) then
       dqsdT = -999.
     endif
     if (present (esat)) then
       esat = -999.
     endif 
   endif ! (nbad = 0)


 end subroutine compute_qs_k_2d

!#######################################################################

 subroutine compute_qs_k_1d (temp, press, eps, zvir, qs, nbad, q, hc, &
                          dqsdT, esat, es_over_liq, es_over_liq_and_ice)

 real, intent(in),  dimension(:)           :: temp, press   
 real, intent(in)                          :: eps, zvir
 real, intent(out), dimension(:)           :: qs   
 integer, intent(out)                      :: nbad
 real, intent(in),  dimension(:), optional :: q
 real, intent(in),                optional :: hc
 real, intent(out), dimension(:), optional :: dqsdT, esat
 logical,intent(in),              optional :: es_over_liq
 logical,intent(in),                  optional :: es_over_liq_and_ice

 real, dimension(size(temp,1)) :: esloc, desat, denom
 integer :: i
 real    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0
   endif

 if (present(es_over_liq)) then
   if (present (dqsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dqsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dqsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif 
   if (nbad == 0) then
     if (present (q) .and. use_exact_qs) then
       qs = (1.0 + zvir*q)*eps*esloc/press
       if (present (dqsdT)) then
         dqsdT = (1.0 + zvir*q)*eps*desat/press
       endif
     else  ! (present(q))
       denom = press - (1.0 - eps)*esloc
       do i=1,size(qs,1)
         if (denom(i) >  0.0) then
           qs(i) = eps*esloc(i)/denom(i)
         else
           qs(i) = eps
         endif
       end do
       if (present (dqsdT)) then
         dqsdT = eps*press*desat/denom**2
       endif
     endif ! (present(q))
   else ! (nbad = 0)
     qs = -999.
     if (present (dqsdT)) then
       dqsdT = -999.
     endif
     if (present (esat)) then
       esat = -999.
     endif 
   endif ! (nbad = 0)


 end subroutine compute_qs_k_1d

!#######################################################################

 subroutine compute_qs_k_0d (temp, press, eps, zvir, qs, nbad, q, hc, &
                          dqsdT, esat, es_over_liq, es_over_liq_and_ice)

 real, intent(in)                :: temp, press   
 real, intent(in)                :: eps, zvir
 real, intent(out)               :: qs   
 integer, intent(out)            :: nbad
 real, intent(in),      optional :: q
 real, intent(in),      optional :: hc
 real, intent(out),     optional :: dqsdT, esat
 logical,intent(in),    optional :: es_over_liq
 logical,intent(in),                  optional :: es_over_liq_and_ice

 real    :: esloc, desat, denom
 real    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0
   endif

 if (present(es_over_liq)) then
   if (present (dqsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dqsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dqsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif 
   if (nbad == 0) then
     if (present (q) .and. use_exact_qs) then
       qs = (1.0 + zvir*q)*eps*esloc/press
       if (present (dqsdT)) then
         dqsdT = (1.0 + zvir*q)*eps*desat/press
       endif
     else  ! (present(q))
       denom = press - (1.0 - eps)*esloc
       if (denom > 0.0) then
         qs = eps*esloc/denom
       else
         qs = eps
       endif
       if (present (dqsdT)) then
         dqsdT = eps*press*desat/denom**2
       endif
     endif ! (present(q))
   else ! (nbad = 0)
     qs = -999.
     if (present (dqsdT)) then
       dqsdT = -999.
     endif
     if (present (esat)) then
       esat = -999.
     endif 
   endif ! (nbad = 0)


 end subroutine compute_qs_k_0d

!#######################################################################

!#######################################################################

 subroutine compute_mrs_k_3d (temp, press, eps, zvir, mrs, nbad,   &
                 mr, hc, dmrsdT, esat,es_over_liq, es_over_liq_and_ice)

 real, intent(in),  dimension(:,:,:)           :: temp, press
 real, intent(in)                              :: eps, zvir
 real, intent(out), dimension(:,:,:)           :: mrs   
 integer, intent(out)                          :: nbad
 real, intent(in),  dimension(:,:,:), optional :: mr
 real, intent(in),                    optional :: hc
 real, intent(out), dimension(:,:,:), optional :: dmrsdT, esat
 logical,intent(in),                  optional :: es_over_liq
 logical,intent(in),                  optional :: es_over_liq_and_ice

 real, dimension(size(temp,1), size(temp,2), size(temp,3)) ::    &
                                                    esloc, desat, denom
 integer :: i, j, k
 real    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0
   endif

 if (present (es_over_liq)) then
   if (present (dmrsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dmrsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dmrsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif 
   if (nbad == 0) then
     if (present (mr) .and. use_exact_qs) then
       mrs = (eps + mr)*esloc/press
       if (present (dmrsdT)) then
         dmrsdT =  (eps + mr)*desat/press
       endif
     else ! (present (mr))
       denom = press - esloc
       do k=1,size(mrs,3)
         do j=1,size(mrs,2)
           do i=1,size(mrs,1)
             if (denom(i,j,k) > 0.0) then
               mrs(i,j,k) = eps*esloc(i,j,k)/denom(i,j,k) 
             else
               mrs(i,j,k) = eps
             endif
           end do
         end do
       end do
       if (present (dmrsdT)) then
         dmrsdT = eps*press*desat/denom**2
       endif
     endif !(present (mr))
   else
     mrs = -999.
     if (present (dmrsdT)) then
       dmrsdT = -999.
     endif
     if (present (esat)) then
       esat = -999.
     endif 
   endif

     
 end subroutine compute_mrs_k_3d

!#######################################################################

 subroutine compute_mrs_k_2d (temp, press, eps, zvir, mrs, nbad,  &
                 mr, hc, dmrsdT, esat,es_over_liq, es_over_liq_and_ice)

 real, intent(in),  dimension(:,:)           :: temp, press
 real, intent(in)                            :: eps, zvir
 real, intent(out), dimension(:,:)           :: mrs   
 integer, intent(out)                        :: nbad
 real, intent(in), dimension(:,:), optional  :: mr
 real, intent(in),                 optional :: hc
 real, intent(out), dimension(:,:), optional :: dmrsdT, esat
 logical,intent(in),               optional :: es_over_liq
 logical,intent(in),                  optional :: es_over_liq_and_ice

 real, dimension(size(temp,1), size(temp,2)) :: esloc, desat, denom
 integer :: i, j
 real    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0
   endif

 if (present (es_over_liq)) then
   if (present (dmrsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dmrsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dmrsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif 
   if (nbad == 0) then
     if (present (mr) .and. use_exact_qs) then
       mrs = (eps + mr)*esloc/press
       if (present (dmrsdT)) then
         dmrsdT = (eps + mr)*desat/press
       endif
     else ! (present (mr))
       denom = press - esloc
       do j=1,size(mrs,2)
         do i=1,size(mrs,1)
           if (denom(i,j) > 0.0) then
             mrs(i,j) = eps*esloc(i,j)/denom(i,j) 
           else
             mrs(i,j) = eps
           endif
         end do
       end do
       if (present (dmrsdT)) then
         dmrsdT = eps*press*desat/denom**2
       endif
     endif !(present (mr))
   else
     mrs = -999.
     if (present (dmrsdT)) then
       dmrsdT = -999.
     endif
     if (present (esat)) then
       esat = -999.
     endif 
   endif


 end subroutine compute_mrs_k_2d

!#######################################################################

 subroutine compute_mrs_k_1d (temp, press, eps, zvir, mrs, nbad,  &
                 mr, hc, dmrsdT, esat,es_over_liq, es_over_liq_and_ice)

 real, intent(in),  dimension(:)           :: temp, press
 real, intent(in)                          :: eps, zvir
 real, intent(out), dimension(:)           :: mrs   
 integer, intent(out)                      :: nbad
 real, intent(in),  dimension(:), optional :: mr
 real, intent(in),                optional :: hc
 real, intent(out), dimension(:), optional :: dmrsdT, esat
 logical,intent(in),              optional :: es_over_liq
 logical,intent(in),                  optional :: es_over_liq_and_ice

 real, dimension(size(temp,1)) :: esloc, desat, denom
 integer :: i
 real    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0
   endif

 if (present (es_over_liq)) then
   if (present (dmrsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dmrsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dmrsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif 
   if (nbad == 0) then
     if (present (mr) .and. use_exact_qs) then
       mrs = (eps + mr)*esloc/press
       if (present (dmrsdT)) then
         dmrsdT =  (eps + mr)*desat/press
       endif
     else ! (present (mr))
       denom = press - esloc
       do i=1,size(mrs,1)
         if (denom(i) > 0.0) then
           mrs(i) = eps*esloc(i)/denom(i) 
         else
           mrs(i) = eps
         endif
       end do
       if (present (dmrsdT)) then
         dmrsdT = eps*press*desat/denom**2
       endif
     endif !(present (mr))
   else
     mrs = -999.
     if (present (dmrsdT)) then
       dmrsdT = -999.
     endif
     if (present (esat)) then
       esat = -999.
     endif 
   endif


 end subroutine compute_mrs_k_1d

!#######################################################################

 subroutine compute_mrs_k_0d (temp, press, eps, zvir, mrs, nbad,   &
                 mr, hc, dmrsdT, esat,es_over_liq, es_over_liq_and_ice)

 real, intent(in)                              :: temp, press
 real, intent(in)                              :: eps, zvir
 real, intent(out)                             :: mrs   
 integer, intent(out)                          :: nbad
 real, intent(in),                    optional :: mr
 real, intent(in),                    optional :: hc
 real, intent(out),                   optional :: dmrsdT, esat
 logical,intent(in),                  optional :: es_over_liq
 logical,intent(in),                  optional :: es_over_liq_and_ice

 real    :: esloc, desat, denom
 real    :: hc_loc

   if (present(hc)) then
     hc_loc = hc
   else
     hc_loc = 1.0
   endif

 if (present (es_over_liq)) then
   if (present (dmrsdT)) then
     call lookup_es2_des2_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es2_k (temp, esloc, nbad)
   endif
 else if (present(es_over_liq_and_ice)) then
   if (present (dmrsdT)) then
     call lookup_es3_des3_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es3_k (temp, esloc, nbad)
   endif
 else
   if (present (dmrsdT)) then
     call lookup_es_des_k (temp, esloc, desat, nbad)
     desat = desat*hc_loc
   else
     call lookup_es_k (temp, esloc, nbad)
   endif
 endif
   esloc = esloc*hc_loc
   if (present (esat)) then
     esat = esloc
   endif 
   if (nbad == 0) then
     if (present (mr) .and. use_exact_qs) then
       mrs = (eps + mr)*esloc/press
       if (present (dmrsdT)) then
         dmrsdT = (eps + mr)*desat/press
       endif
     else ! (present (mr))
       denom = press - esloc
       if (denom > 0.0) then
         mrs = eps*esloc/denom 
       else
         mrs = eps       
       endif
       if (present (dmrsdT)) then
         dmrsdT = eps*press*desat/denom**2
       endif
     endif !(present (mr))
   else
     mrs = -999.
     if (present (dmrsdT)) then
       dmrsdT = -999.
     endif
     if (present (esat)) then
       esat = -999.
     endif 
   endif


 end subroutine compute_mrs_k_0d



!#######################################################################

 subroutine lookup_es_des_k_3d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: esat, desat
 integer, intent(out)                 :: nbad

 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz) then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j,k) = TABLE(ind+1) +  &
                     del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
       desat(i,j,k) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_es_des_k_3d

!#######################################################################

 subroutine lookup_es_des_k_2d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: esat, desat
 integer, intent(out)               :: nbad

 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j) = TABLE(ind+1) + &
                   del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
       desat(i,j) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     endif
   enddo
   enddo

 end subroutine lookup_es_des_k_2d

!#######################################################################

 subroutine lookup_es_des_k_1d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: esat, desat
 integer, intent(out)             :: nbad

 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i) = TABLE(ind+1) + &
                   del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
       desat(i) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     endif
   enddo

 end subroutine lookup_es_des_k_1d

!#######################################################################

 subroutine lookup_es_des_k_0d (temp, esat, desat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: esat, desat
 integer, intent(out) :: nbad

 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     esat = TABLE(ind+1) + &
            del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
     desat = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
   endif

 end subroutine lookup_es_des_k_0d

!#######################################################################

 subroutine lookup_es_k_3d(temp, esat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j,k) = TABLE(ind+1) + &
                     del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_es_k_3d

!#######################################################################

 subroutine lookup_des_k_3d(temp, desat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i,j,k) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_des_k_3d

!#######################################################################
 subroutine lookup_des_k_2d(temp, desat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i,j) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     endif
   enddo
   enddo

 end subroutine lookup_des_k_2d
!#######################################################################
 subroutine lookup_es_k_2d(temp, esat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j) = TABLE(ind+1) + del*(DTABLE(ind+1) +   &
                                                  del*D2TABLE(ind+1))
     endif
   enddo
   enddo

 end subroutine lookup_es_k_2d
!#######################################################################
 subroutine lookup_des_k_1d(temp, desat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     endif 
   enddo

 end subroutine lookup_des_k_1d
!#######################################################################
 subroutine lookup_es_k_1d(temp, esat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
     endif
   enddo

 end subroutine lookup_es_k_1d
!#######################################################################
 subroutine lookup_des_k_0d(temp, desat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     desat = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
   endif 

 end subroutine lookup_des_k_0d
!#######################################################################
 subroutine lookup_es_k_0d(temp, esat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     esat = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
   endif 

 end subroutine lookup_es_k_0d
!#######################################################################

 subroutine lookup_es2_des2_k_3d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: esat, desat
 integer, intent(out)                 :: nbad

 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz) then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j,k) = TABLE2(ind+1) +  &
                     del*(DTABLE2(ind+1) + del*D2TABLE2(ind+1))
       desat(i,j,k) = DTABLE2(ind+1) + 2.*del*D2TABLE2(ind+1)
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_es2_des2_k_3d

!#######################################################################

 subroutine lookup_es2_des2_k_2d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: esat, desat
 integer, intent(out)               :: nbad

 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j) = TABLE2(ind+1) + &
                   del*(DTABLE2(ind+1) + del*D2TABLE2(ind+1))
       desat(i,j) = DTABLE2(ind+1) + 2.*del*D2TABLE2(ind+1)
     endif
   enddo
   enddo

 end subroutine lookup_es2_des2_k_2d

!#######################################################################

 subroutine lookup_es2_des2_k_1d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: esat, desat
 integer, intent(out)             :: nbad

 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i) = TABLE2(ind+1) + &
                   del*(DTABLE2(ind+1) + del*D2TABLE2(ind+1))
       desat(i) = DTABLE2(ind+1) + 2.*del*D2TABLE2(ind+1)
     endif
   enddo

 end subroutine lookup_es2_des2_k_1d

!#######################################################################

 subroutine lookup_es2_des2_k_0d (temp, esat, desat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: esat, desat
 integer, intent(out) :: nbad

 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     esat = TABLE2(ind+1) + &
            del*(DTABLE2(ind+1) + del*D2TABLE2(ind+1))
     desat = DTABLE2(ind+1) + 2.*del*D2TABLE2(ind+1)
   endif

 end subroutine lookup_es2_des2_k_0d

!#######################################################################

 subroutine lookup_es2_k_3d(temp, esat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j,k) = TABLE2(ind+1) + &
                     del*(DTABLE2(ind+1) + del*D2TABLE2(ind+1))
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_es2_k_3d

!#######################################################################

 subroutine lookup_des2_k_3d(temp, desat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i,j,k) = DTABLE2(ind+1) + 2.*del*D2TABLE2(ind+1)
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_des2_k_3d

!#######################################################################
 subroutine lookup_des2_k_2d(temp, desat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i,j) = DTABLE2(ind+1) + 2.*del*D2TABLE2(ind+1)
     endif
   enddo
   enddo

 end subroutine lookup_des2_k_2d
!#######################################################################
 subroutine lookup_es2_k_2d(temp, esat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j) = TABLE2(ind+1) + del*(DTABLE2(ind+1) +   &
                                                  del*D2TABLE2(ind+1))
     endif
   enddo
   enddo

 end subroutine lookup_es2_k_2d
!#######################################################################
 subroutine lookup_des2_k_1d(temp, desat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i) = DTABLE2(ind+1) + 2.*del*D2TABLE2(ind+1)
     endif 
   enddo

 end subroutine lookup_des2_k_1d
!#######################################################################
 subroutine lookup_es2_k_1d(temp, esat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i) = TABLE2(ind+1) + del*(DTABLE2(ind+1) + del*D2TABLE2(ind+1))
     endif
   enddo

 end subroutine lookup_es2_k_1d
!#######################################################################
 subroutine lookup_des2_k_0d(temp, desat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     desat = DTABLE2(ind+1) + 2.*del*D2TABLE2(ind+1)
   endif 

 end subroutine lookup_des2_k_0d
!#######################################################################
 subroutine lookup_es2_k_0d(temp, esat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     esat = TABLE2(ind+1) + del*(DTABLE2(ind+1) + del*D2TABLE2(ind+1))
   endif 

 end subroutine lookup_es2_k_0d
!#######################################################################

!#######################################################################

 subroutine lookup_es3_des3_k_3d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: esat, desat
 integer, intent(out)                 :: nbad

 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz) then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j,k) = TABLE3(ind+1) +  &
                     del*(DTABLE3(ind+1) + del*D2TABLE3(ind+1))
       desat(i,j,k) = DTABLE3(ind+1) + 2.*del*D2TABLE3(ind+1)
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_es3_des3_k_3d

!#######################################################################

 subroutine lookup_es3_des3_k_2d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: esat, desat
 integer, intent(out)               :: nbad

 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j) = TABLE3(ind+1) + &
                   del*(DTABLE3(ind+1) + del*D2TABLE3(ind+1))
       desat(i,j) = DTABLE3(ind+1) + 2.*del*D2TABLE3(ind+1)
     endif
   enddo
   enddo

 end subroutine lookup_es3_des3_k_2d

!#######################################################################

 subroutine lookup_es3_des3_k_1d (temp, esat, desat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: esat, desat
 integer, intent(out)             :: nbad

 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i) = TABLE3(ind+1) + &
                   del*(DTABLE3(ind+1) + del*D2TABLE3(ind+1))
       desat(i) = DTABLE3(ind+1) + 2.*del*D2TABLE3(ind+1)
     endif
   enddo

 end subroutine lookup_es3_des3_k_1d

!#######################################################################

 subroutine lookup_es3_des3_k_0d (temp, esat, desat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: esat, desat
 integer, intent(out) :: nbad

 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     esat = TABLE3(ind+1) + &
            del*(DTABLE3(ind+1) + del*D2TABLE3(ind+1))
     desat = DTABLE3(ind+1) + 2.*del*D2TABLE3(ind+1)
   endif

 end subroutine lookup_es3_des3_k_0d

!#######################################################################

 subroutine lookup_es3_k_3d(temp, esat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j,k) = TABLE3(ind+1) + &
                     del*(DTABLE3(ind+1) + del*D2TABLE3(ind+1))
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_es3_k_3d

!#######################################################################

 subroutine lookup_des3_k_3d(temp, desat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i,j,k) = DTABLE3(ind+1) + 2.*del*D2TABLE3(ind+1)
     endif
   enddo
   enddo
   enddo

 end subroutine lookup_des3_k_3d

!#######################################################################
 subroutine lookup_des3_k_2d(temp, desat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i,j) = DTABLE3(ind+1) + 2.*del*D2TABLE3(ind+1)
     endif
   enddo
   enddo

 end subroutine lookup_des3_k_2d
!#######################################################################
 subroutine lookup_es3_k_2d(temp, esat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i,j) = TABLE3(ind+1) + del*(DTABLE3(ind+1) +   &
                                                  del*D2TABLE3(ind+1))
     endif
   enddo
   enddo

 end subroutine lookup_es3_k_2d
!#######################################################################
 subroutine lookup_des3_k_1d(temp, desat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       desat(i) = DTABLE3(ind+1) + 2.*del*D2TABLE3(ind+1)
     endif 
   enddo

 end subroutine lookup_des3_k_1d
!#######################################################################
 subroutine lookup_es3_k_1d(temp, esat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     if (ind < 0 .or. ind >= table_siz)  then
       nbad = nbad+1
     else
       del = tmp-dtres*real(ind)
       esat(i) = TABLE3(ind+1) + del*(DTABLE3(ind+1) + del*D2TABLE3(ind+1))
     endif
   enddo

 end subroutine lookup_es3_k_1d
!#######################################################################
 subroutine lookup_des3_k_0d(temp, desat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     desat = DTABLE3(ind+1) + 2.*del*D2TABLE3(ind+1)
   endif 

 end subroutine lookup_des3_k_0d
!#######################################################################
 subroutine lookup_es3_k_0d(temp, esat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   if (ind < 0 .or. ind >= table_siz)  then
     nbad = nbad+1
   else
     del = tmp-dtres*real(ind)
     esat = TABLE3(ind+1) + del*(DTABLE3(ind+1) + del*D2TABLE3(ind+1))
   endif 

 end subroutine lookup_es3_k_0d
!#######################################################################
 end module sat_vapor_pres_k_mod

