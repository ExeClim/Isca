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

module horiz_interp_bilinear_mod

  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

  ! <OVERVIEW>
  !   Performs spatial interpolation between grids using bilinear interpolation
  ! </OVERVIEW>

  ! <DESCRIPTION>
  !     This module can interpolate data from regular rectangular grid
  !     to rectangular/tripolar grid. The interpolation scheme is bilinear interpolation.
  !     There is an optional mask field for missing input data.
  !     An optional output mask field may be used in conjunction with
  !     the input mask to show where output data exists.
  ! </DESCRIPTION>

  use mpp_mod,               only: mpp_error, FATAL, stdout, mpp_pe, mpp_root_pe
  use fms_mod,               only: write_version_number
  use constants_mod,         only: PI
  use horiz_interp_type_mod, only: horiz_interp_type, stats

  implicit none
  private


  public :: horiz_interp_bilinear_new, horiz_interp_bilinear, horiz_interp_bilinear_del
  public :: horiz_interp_bilinear_init

  !--- public interface
  interface horiz_interp_bilinear_new
    module procedure horiz_interp_bilinear_new_1d
    module procedure horiz_interp_bilinear_new_2d
  end interface


  real, parameter :: epsln=1.e-10

  !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: horiz_interp_bilinear.F90,v 14.0 2007/03/15 22:39:57 fms Exp $'
  character(len=128) :: tagname = '$Name: siena_201211 $'
  logical            :: module_is_initialized = .FALSE.

contains

  !#######################################################################
  !  <SUBROUTINE NAME="horiz_interp_bilinear_init">
  !  <OVERVIEW>
  !     writes version number and tag name to logfile.out
  !  </OVERVIEW>
  !  <DESCRIPTION>       
  !     writes version number and tag name to logfile.out
  !  </DESCRIPTION>

  subroutine horiz_interp_bilinear_init

    if(module_is_initialized) return
    call write_version_number (version, tagname)
    module_is_initialized = .true.

  end subroutine horiz_interp_bilinear_init

  !  </SUBROUTINE>

  !########################################################################

  subroutine horiz_interp_bilinear_new_1d ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in),  dimension(:)        :: lon_in , lat_in
    real, intent(in),  dimension(:,:)      :: lon_out, lat_out
    integer, intent(in),          optional :: verbose
    logical, intent(in),          optional :: src_modulo

    logical :: src_is_modulo
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, n, m
    integer :: ie, is, je, js, ln_err, lt_err, warns, unit
    real    :: wtw, wte, wts, wtn, lon, lat, tpi, hpi
    real    :: glt_min, glt_max, gln_min, gln_max, min_lon, max_lon

    warns = 0
    if(present(verbose)) warns = verbose
    src_is_modulo = .true. 
    if (present(src_modulo)) src_is_modulo = src_modulo

    hpi = 0.5*pi
    tpi = 4.0*hpi
    glt_min = hpi
    glt_max = -hpi
    gln_min = tpi
    gln_max = -tpi
    min_lon = 0.0
    max_lon = tpi
    ln_err = 0
    lt_err = 0
    !-----------------------------------------------------------------------

    allocate ( Interp % wti (size(lon_out,1),size(lon_out,2),2),   &
               Interp % wtj (size(lon_out,1),size(lon_out,2),2),   &
               Interp % i_lon (size(lon_out,1),size(lon_out,2),2), &
               Interp % j_lat (size(lon_out,1),size(lon_out,2),2))
    !-----------------------------------------------------------------------

    nlon_in = size(lon_in(:))  ; nlat_in = size(lat_in(:))
    nlon_out = size(lon_out, 1); nlat_out = size(lon_out, 2)
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out

    if(src_is_modulo) then
       if(lon_in(nlon_in) - lon_in(1) .gt. tpi + epsln) &
            call mpp_error(FATAL,'horiz_interp_bilinear_mod: '// & 
            'The range of source grid longitude should be no larger than tpi')

       if(lon_in(1) .lt. 0.0 .OR. lon_in(nlon_in) > tpi ) then
          min_lon = lon_in(1)
          max_lon = lon_in(nlon_in)
       endif
    endif

    do n = 1, nlat_out
       do m = 1, nlon_out
          lon = lon_out(m,n)
          lat = lat_out(m,n)

          if(src_is_modulo) then
             if(lon .lt. min_lon) then
                lon = lon + tpi
             else if(lon .gt. max_lon) then
                lon = lon - tpi
             endif
          else  ! when the input grid is in not cyclic, the output grid should located inside
             ! the input grid
             if((lon .lt. lon_in(1)) .or. (lon .gt. lon_in(nlon_in))) &
                  call mpp_error(FATAL,'horiz_interp_bilinear_mod: ' //&
                  'when input grid is not modulo, output grid should locate inside input grid')
          endif

          glt_min = min(lat,glt_min);  glt_max = max(lat,glt_max)
          gln_min = min(lon,gln_min);  gln_max = max(lon,gln_max)

          is = indp(lon, lon_in ) 
          if( lon_in(is) .gt. lon ) is = max(is-1,1)
          if( lon_in(is) .eq. lon .and. is .eq. nlon_in) is = max(is - 1,1)
          ie = min(is+1,nlon_in)
          if(lon_in(is) .ne. lon_in(ie) .and. lon_in(is) .le. lon) then
             wtw = ( lon_in(ie) - lon) / (lon_in(ie) - lon_in(is) )
          else
             !     east or west of the last data value. this could be because a
             !     cyclic condition is needed or the dataset is too small. 
             ln_err = 1
             ie = 1
             is = nlon_in
             if (lon_in(ie) .ge. lon ) then
                wtw = (lon_in(ie) -lon)/(lon_in(ie)-lon_in(is)+tpi+epsln)
             else
                wtw = (lon_in(ie) -lon+tpi+epsln)/(lon_in(ie)-lon_in(is)+tpi+epsln)
             endif
          endif
          wte = 1. - wtw

          js = indp(lat, lat_in ) 

          if( lat_in(js) .gt. lat ) js = max(js - 1, 1)
          if( lat_in(js) .eq. lat .and. js .eq. nlat_in) js = max(js - 1, 1)
          je = min(js + 1, nlat_in)

          if ( lat_in(js) .ne. lat_in(je) .and. lat_in(js) .le. lat) then
             wts = ( lat_in(je) - lat )/(lat_in(je)-lat_in(js))
          else
             !     north or south of the last data value. this could be because a
             !     pole is not included in the data set or the dataset is too small.
             !     in either case extrapolate north or south
             lt_err = 1
             wts = 1.
          endif

          wtn = 1. - wts

          Interp % i_lon (m,n,1) = is; Interp % i_lon (m,n,2) = ie
          Interp % j_lat (m,n,1) = js; Interp % j_lat (m,n,2) = je
          Interp % wti   (m,n,1) = wtw
          Interp % wti   (m,n,2) = wte
          Interp % wtj   (m,n,1) = wts
          Interp % wtj   (m,n,2) = wtn

       enddo
    enddo

    unit = stdout()

    if (ln_err .eq. 1 .and. warns > 0) then
       write (unit,'(/,(1x,a))')                                      &
            '==> Warning: the geographic data set does not extend far   ', &
            '             enough east or west - a cyclic boundary       ', &
            '             condition was applied. check if appropriate   '
       write (unit,'(/,(1x,a,2f8.4))')                                &
            '    data required between longitudes:', gln_min, gln_max,     &
            '      data set is between longitudes:', lon_in(1), lon_in(nlon_in)
       warns = warns - 1
    endif

    if (lt_err .eq. 1 .and. warns > 0) then
       write (unit,'(/,(1x,a))')                                     &
            '==> Warning: the geographic data set does not extend far   ',&
            '             enough north or south - extrapolation from    ',&
            '             the nearest data was applied. this may create ',&
            '             artificial gradients near a geographic pole   ' 
       write (unit,'(/,(1x,a,2f8.4))')                             &
            '    data required between latitudes:', glt_min, glt_max,   &
            '      data set is between latitudes:', lat_in(1), lat_in(nlat_in)
    endif

    return

  end subroutine horiz_interp_bilinear_new_1d

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_bilinear_new">

  !   <OVERVIEW>
  !      Initialization routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !      Allocates space and initializes a derived-type variable
  !      that contains pre-computed interpolation indices and weights.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, src_modulo )

  !   </TEMPLATE>
  !   
  !   <IN NAME="lon_in" TYPE="real, dimension(:,:)" UNITS="radians">
  !      Longitude (in radians) for source data grid. 
  !   </IN>

  !   <IN NAME="lat_in" TYPE="real, dimension(:,:)" UNITS="radians">
  !      Latitude (in radians) for source data grid.
  !   </IN>

  !   <IN NAME="lon_out" TYPE="real, dimension(:,:)" UNITS="radians" >
  !      Longitude (in radians) for source data grid. 
  !   </IN>

  !   <IN NAME="lat_out" TYPE="real, dimension(:,:)" UNITS="radians" >
  !      Latitude (in radians) for source data grid. 
  !   </IN>

  !   <IN NAME="src_modulo" TYPE="logical, optional">
  !      logical variable to indicate if the boundary condition along zonal boundary
  !      is cyclic or not. When true, the zonal boundary condition is cyclic.
  !   </IN>

  !   <IN NAME="verbose" TYPE="integer, optional" >
  !      flag for the amount of print output.
  !   </IN>

  !   <INOUT NAME="Interp" TYPE="type(horiz_interp_type)" >
  !      A derived-type variable containing indices and weights used for subsequent 
  !      interpolations. To reinitialize this variable for a different grid-to-grid 
  !      interpolation you must first use the "horiz_interp_del" interface.
  !   </INOUT>

  subroutine horiz_interp_bilinear_new_2d ( Interp, lon_in, lat_in, lon_out, lat_out, &
       verbose, src_modulo )

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in),  dimension(:,:)      :: lon_in , lat_in
    real, intent(in),  dimension(:,:)      :: lon_out, lat_out
    integer, intent(in),          optional :: verbose
    logical, intent(in),          optional :: src_modulo
    integer                                :: warns 
    logical                                :: src_is_modulo
    integer                                :: nlon_in, nlat_in, nlon_out, nlat_out
    integer                                :: m, n, is, ie, js, je, num_solution
    real                                   :: lon, lat, quadra, x, y, y1, y2
    real                                   :: a1, b1, c1, d1, a2, b2, c2, d2, a, b, c
    real                                   :: lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4
    real                                   :: tpi, lon_min, lon_max

    tpi = 2.0*pi

    warns = 0
    if(present(verbose)) warns = verbose
    src_is_modulo = .true. 
    if (present(src_modulo)) src_is_modulo = src_modulo

    ! make sure lon and lat has the same dimension
    if(size(lon_out,1) /= size(lat_out,1) .or. size(lon_out,2) /= size(lat_out,2) ) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: when using bilinear ' // &
         'interplation, the output grids should be geographical grids')    

    if(size(lon_in,1) /= size(lat_in,1) .or. size(lon_in,2) /= size(lat_in,2) ) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: when using bilinear '// &
         'interplation, the input grids should be geographical grids')  

    !--- get the grid size 
    nlon_in  = size(lon_in,1) ; nlat_in  = size(lat_in,2)
    nlon_out = size(lon_out,1); nlat_out = size(lon_out,2)
    Interp%nlon_src = nlon_in;  Interp%nlat_src = nlat_in
    Interp%nlon_dst = nlon_out; Interp%nlat_dst = nlat_out

    allocate ( Interp % wti (size(lon_out,1),size(lon_out,2),2),   &
               Interp % wtj (size(lon_out,1),size(lon_out,2),2),   &
               Interp % i_lon (size(lon_out,1),size(lon_out,2),2), &
               Interp % j_lat (size(lon_out,1),size(lon_out,2),2))

    !--- first fine the neighbor points for the destination points.
    call find_neighbor(Interp, lon_in, lat_in, lon_out, lat_out, src_is_modulo)

    !***************************************************************************
    !         Algorithm explanation (from disscussion with Steve Garner )      *
    !                                                                          *
    !    lon(x,y) = a1*x + b1*y + c1*x*y + d1         (1)                      *
    !    lat(x,y) = a2*x + b2*y + c2*x*y + d2         (2)                      *
    !    f (x,y) = a3*x + b3*y + c3*x*y + d3          (3)                      *
    !    with x and y is between 0 and 1.                                      *
    !    lon1 = lon(0,0) = d1,          lat1 = lat(0,0) = d2                   *
    !    lon2 = lon(1,0) = a1+d1,       lat2 = lat(1,0) = a2+d2                *
    !    lon3 = lon(1,1) = a1+b1+c1+d1, lat3 = lat(1,1) = a2+b2+c2+d2          *
    !    lon4 = lon(0,1) = b1+d1,       lat4 = lat(0,1) = b2+d2                *
    !    where (lon1,lat1),(lon2,lat2),(lon3,lat3),(lon4,lat4) represents      *
    !    the four corners starting from the left lower corner of grid box      *
    !    that encloses a destination grid ( the rotation direction is          *
    !    counterclockwise ). With these conditions, we get                     *
    !    a1 = lon2-lon1,           a2 = lat2-lat1                              *
    !    b1 = lon4-lon1,           b2 = lat4-lat1                              *
    !    c1 = lon3-lon2-lon4+lon1, c2 = lat3-lat2-lat4+lat1                    *
    !    d1 = lon1                 d2 = lat1                                   *
    !    So given any point (lon,lat), from equation (1) and (2) we can        *
    !    solve (x,y).                                                          *
    !    From equation (3)                                                     *
    !    f1 = f(0,0) = d3,          f2 = f(1,0) = a3+d3                        *
    !    f3 = f(1,1) = a3+b3+c3+d3, f4 = f(0,1) = b3+d3                        *
    !    we obtain                                                             *
    !    a3 = f2-f1,       b3 = f4-f1                                          *
    !    c3 = f3-f2-f4+f1, d3 = f1                                             *
    !    at point (lon,lat) ---> (x,y)                                         *
    !    f(x,y) = (f2-f1)x + (f4-f1)y + (f3-f2-f4+f1)xy + f1                   *
    !           = f1*(1-x)*(1-y) + f2*x*(1-y) + f3*x*y + f4*y*(1-x)            *
    !    wtw=1-x; wte=x; wts=1-y; xtn=y                                        *
    !                                                                          *
    !***************************************************************************

    lon_min = minval(lon_in);
    lon_max = maxval(lon_in);
    !--- calculate the weight
    do n = 1, nlat_out
       do m = 1, nlon_out
          lon = lon_out(m,n)
          lat = lat_out(m,n)
          if(lon .lt. lon_min) then
             lon = lon + tpi
          else if(lon .gt. lon_max) then
             lon = lon - tpi
          endif
          is = Interp%i_lon(m,n,1); ie = Interp%i_lon(m,n,2)
          js = Interp%j_lat(m,n,1); je = Interp%j_lat(m,n,2)
          lon1 = lon_in(is,js); lat1 = lat_in(is,js);
          lon2 = lon_in(ie,js); lat2 = lat_in(ie,js);
          lon3 = lon_in(ie,je); lat3 = lat_in(ie,je);
          lon4 = lon_in(is,je); lat4 = lat_in(is,je); 
          if(lon .lt. lon_min) then
             lon1 = lon1 -tpi; lon4 = lon4 - tpi
          else if(lon .gt. lon_max) then
             lon2 = lon2 +tpi; lon3 = lon3 + tpi
          endif                      
          a1 = lon2-lon1
          b1 = lon4-lon1
          c1 = lon1+lon3-lon4-lon2
          d1 = lon1
          a2 = lat2-lat1
          b2 = lat4-lat1
          c2 = lat1+lat3-lat4-lat2
          d2 = lat1
          !--- the coefficient of the quadratic equation
          a  = b2*c1-b1*c2
          b  = a1*b2-a2*b1+c1*d2-c2*d1+c2*lon-c1*lat
          c  = a2*lon-a1*lat+a1*d2-a2*d1
          quadra = b*b-4*a*c
          if(abs(quadra) < epsln) quadra = 0.0
          if(quadra < 0.0) call mpp_error(FATAL, &
               "horiz_interp_bilinear_mod: No solution existed for this quadratic equation")
          if ( abs(a) .lt. epsln) then  ! a = 0 is a linear equation
             if( abs(b) .lt. epsln) call mpp_error(FATAL, &
                  "horiz_interp_bilinear_mod: no unique solution existed for this linear equation")
             y = -c/b
          else
             y1 = 0.5*(-b+sqrt(quadra))/a
             y2 = 0.5*(-b-sqrt(quadra))/a
             if(abs(y1) < epsln) y1 = 0.0
             if(abs(y2) < epsln) y2 = 0.0
             if(abs(1-y1) < epsln) y1 = 1.0
             if(abs(1-y2) < epsln) y2 = 1.0
             num_solution = 0
             if(y1 .le. 1 .and. y1 .ge. 0) then
                y = y1
                num_solution = num_solution +1
             endif
             if(y2 .le. 1 .and. y2 .ge. 0) then
                y = y2
                num_solution = num_solution + 1
             endif
             if(num_solution == 0) then
                call mpp_error(FATAL, "horiz_interp_bilinear_mod: No solution found")
             else if(num_solution == 2) then
                call mpp_error(FATAL, "horiz_interp_bilinear_mod: Two solutions found")
             endif
           endif
           if(abs(a1+c1*y) < epsln) call mpp_error(FATAL, &
               "horiz_interp_bilinear_mod: the denomenator is 0")
           if(abs(y) < epsln) y = 0.0
           if(abs(1-y) < epsln) y = 1.0
           x = (lon-b1*y-d1)/(a1+c1*y)
           if(abs(x) < epsln) x = 0.0
           if(abs(1-x) < epsln) x = 1.0
           ! x and y should be between 0 and 1.
           if( x>1 .or. x<0 .or. y>1 .or. y < 0) call mpp_error(FATAL, &
               "horiz_interp_bilinear_mod: weight should be between 0 and 1")
           Interp % wti(m,n,1)=1-x; Interp % wti(m,n,2)=x   
           Interp % wtj(m,n,1)=1-y; Interp % wtj(m,n,2)=y          
       enddo
    enddo

  end subroutine horiz_interp_bilinear_new_2d
  ! </SUBROUTINE>

  !#######################################################################
  ! this routine will search the source grid to fine the grid box that encloses 
  ! each destination grid.
  subroutine find_neighbor( Interp, lon_in, lat_in, lon_out, lat_out, src_modulo )
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in),       dimension(:,:) :: lon_in , lat_in
    real, intent(in),       dimension(:,:) :: lon_out, lat_out
    logical,                 intent(in)    :: src_modulo
    integer                                :: nlon_in, nlat_in, nlon_out, nlat_out
    integer                                :: max_step, n, m, l, i, j, ip1, jp1, step
    integer                                :: is, js, jstart, jend, istart, iend, npts
    integer, allocatable, dimension(:)     :: ilon, jlat
    real                                   :: lon_min, lon_max, lon, lat, tpi
    logical                                :: found
    real                                   :: lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4

    tpi = 2.0*pi
    nlon_in  = size(lon_in,1) ; nlat_in  = size(lat_in,2)
    nlon_out = size(lon_out,1); nlat_out = size(lon_out,2)

    lon_min = minval(lon_in);
    lon_max = maxval(lon_in);

    max_step = min(nlon_in,nlat_in)/2 ! can be adjusted if needed
    allocate(ilon(8*max_step), jlat(8*max_step) )

    do n = 1, nlat_out
       do m = 1, nlon_out
          found = .false.
          lon = lon_out(m,n)
          lat = lat_out(m,n)

          if(src_modulo) then
             if(lon .lt. lon_min) then
                lon = lon + tpi
             else if(lon .gt. lon_max) then
                lon = lon - tpi
             endif
          else
             if(lon .lt. lon_min .or. lon .gt. lon_max ) &
             call mpp_error(FATAL,'horiz_interp_bilinear_mod: ' //&
                  'when input grid is not modulo, output grid should locate inside input grid')
          endif
          !--- search for the surrounding four points locatioon.
          if(m==1 .and. n==1) then
             J_LOOP: do j = 1, nlat_in-1
                do i = 1, nlon_in
                   ip1 = i+1
                   jp1 = j+1
                   if(i==nlon_in) then
                      if(src_modulo)then
                         ip1 = 1
                      else
                         cycle
                      endif
                   endif
                   lon1 = lon_in(i,  j);   lat1 = lat_in(i,j)
                   lon2 = lon_in(ip1,j);   lat2 = lat_in(ip1,j)
                   lon3 = lon_in(ip1,jp1); lat3 = lat_in(ip1,jp1)
                   lon4 = lon_in(i,  jp1); lat4 = lat_in(i,  jp1)  

                   if(lon .lt. lon_min .or. lon .gt. lon_max) then
                      if(i .ne. nlon_in) then
                         cycle
                      else
                         if(lon .lt. lon_min) then
                             lon1 = lon1 -tpi; lon4 = lon4 - tpi
                         else if(lon .gt. lon_max) then
                             lon2 = lon2 +tpi; lon3 = lon3 + tpi
                         endif
                      endif
                   endif

                   if(lat .ge. intersect(lon1,lat1,lon2,lat2,lon))then ! south
                      if(lon .le. intersect(lat2,lon2,lat3,lon3,lat))then ! east
                         if(lat .le. intersect(lon3,lat3,lon4,lat4,lon))then ! north
                            if(lon .ge. intersect(lat4,lon4,lat1,lon1,lat))then  ! west
                               found = .true.
                               Interp % i_lon (m,n,1) = i; Interp % i_lon (m,n,2) = ip1
                               Interp % j_lat (m,n,1) = j; Interp % j_lat (m,n,2) = jp1
                               exit J_LOOP
                            endif
                         endif
                      endif
                   endif
                enddo
             enddo J_LOOP
          else
             step = 0
             do while ( .not. found .and. step .lt. max_step )
                !--- take the adajcent point as the starting point
                if(m == 1) then
                   is = Interp % i_lon (m,n-1,1)
                   js = Interp % j_lat (m,n-1,1)
                else
                   is = Interp % i_lon (m-1,n,1)
                   js = Interp % j_lat (m-1,n,1)
                endif
                if(step==0) then
                   npts = 1
                   ilon(1) = is
                   jlat(1) = js
                else
                   npts = 0
                   !--- bottom and top boundary
                   jstart = max(js-step,1)
                   jend   = min(js+step,nlat_in)

                   do l = -step, step
                      i = is+l
                      if(src_modulo)then
                         if( i < 1) then
                            i = i + nlon_in
                         else if (i > nlon_in) then
                            i = i - nlon_in
                         endif
                         if( i < 1 .or. i > nlon_in) call mpp_error(FATAL, &
                              'horiz_interp_bilinear_mod: max_step is too big, decrease max_step' )
                      else
                         if( i < 1 .or. i > nlon_in) cycle
                      endif

                      npts       = npts + 1
                      ilon(npts) = i
                      jlat(npts) = jstart
                      npts       = npts + 1
                      ilon(npts) = i
                      jlat(npts) = jend                         
                   enddo

                   !--- right and left boundary -----------------------------------------------
                   istart = is - step
                   iend   = is + step
                   if(src_modulo) then
                      if( istart < 1)       istart = istart + nlon_in
                      if( iend   > nlon_in) iend   = iend   - nlon_in
                   else 
                      istart = max(istart,1)
                      iend   = min(iend, nlon_in)
                   endif
                   do l = -step, step
                      j = js+l
                         if( j < 1 .or. j > nlat_in) cycle
                         npts = npts+1
                         ilon(npts) = istart
                         jlat(npts) = j
                         npts = npts+1
                         ilon(npts) = iend
                         jlat(npts) = j
                  end do
                end if

                !--- find the surrouding points             
                do l = 1, npts
                   i = ilon(l)
                   j = jlat(l)
                   ip1 = i+1
                   if(ip1>nlon_in) then
                      if(src_modulo) then
                         ip1 = 1
                      else
                         cycle
                      endif
                   endif
                   jp1 = j+1
                   if(jp1>nlat_in) cycle
                   lon1 = lon_in(i,  j);   lat1 = lat_in(i,j)
                   lon2 = lon_in(ip1,j);   lat2 = lat_in(ip1,j)
                   lon3 = lon_in(ip1,jp1); lat3 = lat_in(ip1,jp1)
                   lon4 = lon_in(i,  jp1); lat4 = lat_in(i,  jp1)  

                   if(lon .lt. lon_min .or. lon .gt. lon_max) then
                      if(i .ne. nlon_in) then
                         cycle
                      else
                         if(lon .lt. lon_min) then
                             lon1 = lon1 -tpi; lon4 = lon4 - tpi
                         else if(lon .gt. lon_max) then
                             lon2 = lon2 +tpi; lon3 = lon3 + tpi
                         endif
                      endif
                   endif

                   if(lat .ge. intersect(lon1,lat1,lon2,lat2,lon))then ! south
                      if(lon .le. intersect(lat2,lon2,lat3,lon3,lat))then ! east
                         if(lat .le. intersect(lon3,lat3,lon4,lat4,lon))then !north
                            if(lon .ge. intersect(lat4,lon4,lat1,lon1,lat))then ! west
                               found = .true.
                               is=i; js=j  
                               Interp % i_lon (m,n,1) = i; Interp % i_lon (m,n,2) = ip1
                               Interp % j_lat (m,n,1) = j; Interp % j_lat (m,n,2) = jp1
                               exit
                            endif
                         endif
                      endif
                   endif
                enddo
                step = step + 1
             enddo
          endif
          if(.not.found) then
             call mpp_error(FATAL, &
                  'horiz_interp_bilinear_mod: the destination point is not inside the source grid' )
          endif
       enddo
    enddo

  end subroutine find_neighbor

  !#######################################################################
  function intersect(x1, y1, x2, y2, x)
     real, intent(in) :: x1, y1, x2, y2, x
     real             :: intersect

     intersect = (y2-y1)*(x-x1)/(x2-x1) + y1

  return

  end function intersect

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_bilinear">

  !   <OVERVIEW>
  !      Subroutine for performing the horizontal interpolation between two grids.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Subroutine for performing the horizontal interpolation between two grids. 
  !     horiz_interp_bilinear_new must be called before calling this routine.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_bilinear ( Interp, data_in, data_out, verbose, mask_in,mask_out, missing_value, missing_permit)
  !   </TEMPLATE>
  !   
  !   <IN NAME="Interp" TYPE="type(horiz_interp_type)">
  !     Derived-type variable containing interpolation indices and weights.
  !     Returned by a previous call to horiz_interp_bilinear_new.
  !   </IN>
  !   <IN NAME="data_in" TYPE="real, dimension(:,:)">
  !      Input data on source grid.
  !   </IN>
  !   <IN NAME="verbose" TYPE="integer, optional">
  !      flag for the amount of print output.
  !               verbose = 0, no output; = 1, min,max,means; = 2, still more
  !   </IN>
  !   <IN NAME="mask_in" TYPE="real, dimension(:,:),optional">
  !      Input mask, must be the same size as the input data. The real value of
  !      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points 
  !      that should not be used or have missing data. 
  !   </IN>
  !   <IN NAME="missing_value" TYPE="real, optional">
  !      Use the missing_value to indicate missing data.
  !   </IN>

  !   <IN NAME="missing_permit" TUPE="integer, optional">
  !      numbers of points allowed to miss for the bilinear interpolation. The value
  !      should be between 0 and 3.
  !   </IN>

  !   <OUT NAME="data_out" TYPE="real, dimension(:,:)">
  !      Output data on destination grid.
  !   </OUT>
  !   <OUT NAME="mask_out" TYPE="real, dimension(:,:),optional">
  !      Output mask that specifies whether data was computed.
  !   </OUT>

  subroutine horiz_interp_bilinear ( Interp, data_in, data_out, verbose, mask_in,mask_out, &
       missing_value, missing_permit)
    !-----------------------------------------------------------------------
    type (horiz_interp_type), intent(in)        :: Interp
    real, intent(in),  dimension(:,:)           :: data_in
    real, intent(out), dimension(:,:)           :: data_out
    integer, intent(in),               optional :: verbose
    real, intent(in), dimension(:,:),  optional :: mask_in
    real, intent(out), dimension(:,:), optional :: mask_out
    real, intent(in),                  optional :: missing_value
    integer, intent(in),               optional :: missing_permit
    !-----------------------------------------------------------------------
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, n, m,         &
         is, ie, js, je, iverbose, max_missing, num_missing, &
         miss_in, miss_out, unit
    real    :: dwtsum, wtsum, min_in, max_in, avg_in, &
         min_out, max_out, avg_out, wtw, wte, wts, wtn
    real    :: mask(size(data_in,1), size(data_in,2) )

    num_missing = 0

    nlon_in  = Interp%nlon_src;  nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    if(present(mask_in)) then
       mask = mask_in
    else
       mask = 1.0
    endif

    if (present(verbose)) then
       iverbose = verbose
    else
       iverbose = 0
    endif

    if(present(missing_permit)) then
       max_missing = missing_permit
    else
       max_missing = 0
    endif

    if(max_missing .gt. 3 .or. max_missing .lt. 0) call mpp_error(FATAL, &
         'horiz_interp_bilinear_mod: missing_permit should be between 0 and 3')

    if (size(data_in,1) /= nlon_in .or. size(data_in,2) /= nlat_in) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: size of input array incorrect')

    if (size(data_out,1) /= nlon_out .or. size(data_out,2) /= nlat_out) &
         call mpp_error(FATAL,'horiz_interp_bilinear_mod: size of output array incorrect')

    do n = 1, nlat_out
       do m = 1, nlon_out
          is = Interp % i_lon (m,n,1); ie = Interp % i_lon (m,n,2)
          js = Interp % j_lat (m,n,1); je = Interp % j_lat (m,n,2)
          wtw = Interp % wti   (m,n,1)
          wte = Interp % wti   (m,n,2)
          wts = Interp % wtj   (m,n,1)
          wtn = Interp % wtj   (m,n,2)

          if(present(missing_value) ) then
             num_missing = 0
             if(data_in(is,js) == missing_value) then
                num_missing = num_missing+1
                mask(is,js) = 0.0
             endif
             if(data_in(ie,js) == missing_value) then
                num_missing = num_missing+1
                mask(ie,js) = 0.0
             endif
             if(data_in(ie,je) == missing_value) then
                num_missing = num_missing+1
                mask(ie,je) = 0.0
             endif
             if(data_in(is,je) == missing_value) then
                num_missing = num_missing+1
                mask(is,je) = 0.0
             endif
          endif

          dwtsum = data_in(is,js)*mask(is,js)*wtw*wts &
               + data_in(ie,js)*mask(ie,js)*wte*wts &
               + data_in(ie,je)*mask(ie,je)*wte*wtn &
               + data_in(is,je)*mask(is,je)*wtw*wtn 
          wtsum  = mask(is,js)*wtw*wts + mask(ie,js)*wte*wts  &
               + mask(ie,je)*wte*wtn + mask(is,je)*wtw*wtn

          if(.not. present(mask_in) .and. .not. present(missing_value)) wtsum = 1.0

          if(num_missing .gt. max_missing ) then
             data_out(m,n) = missing_value
             if(present(mask_out)) mask_out(m,n) = 0.0
          else if(wtsum .lt. epsln) then 
             if(present(missing_value)) then
                data_out(m,n) = missing_value
             else
                data_out(m,n) = 0.0
             endif
             if(present(mask_out)) mask_out(m,n) = 0.0      
          else
             data_out(m,n) = dwtsum/wtsum
             if(present(mask_out)) mask_out(m,n) = wtsum
          endif
       enddo
    enddo
    !***********************************************************************
    ! compute statistics: minimum, maximum, and mean
    !-----------------------------------------------------------------------
    if (iverbose > 0) then

       ! compute statistics of input data

       call stats (data_in, min_in, max_in, avg_in, miss_in, missing_value, mask_in)

       ! compute statistics of output data
       call stats (data_out, min_out, max_out, avg_out, miss_out, missing_value, mask_out)

       !---- output statistics ----
       unit = stdout()
       write (unit,900)
       write (unit,901)  min_in ,max_in, avg_in
       if (present(mask_in))  write (unit,903)  miss_in
       write (unit,902)  min_out,max_out,avg_out
       if (present(mask_out)) write (unit,903)  miss_out

900    format (/,1x,10('-'),' output from horiz_interp ',10('-'))
901    format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
902    format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
903    format ('          number of missing points = ',i6)

    endif

    return

  end subroutine horiz_interp_bilinear
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="horiz_interp_bilinear_del">

  !   <OVERVIEW>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_bilinear_new.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "horiz_interp_type" variables.
  !     Must be called before reinitializing with horiz_interp_bilinear_new.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call horiz_interp_bilinear_del ( Interp )
  !   </TEMPLATE>

  !   <INOUT NAME="Interp" TYPE="horiz_interp_type">
  !     A derived-type variable returned by previous call
  !     to horiz_interp_bilinear_new. The input variable must have
  !     allocated arrays. The returned variable will contain
  !     deallocated arrays.
  !   </INOUT>

  subroutine horiz_interp_bilinear_del( Interp )

    type (horiz_interp_type), intent(inout) :: Interp

    if(associated(Interp%wti))   deallocate(Interp%wti)
    if(associated(Interp%wtj))   deallocate(Interp%wtj)
    if(associated(Interp%i_lon)) deallocate(Interp%i_lon)
    if(associated(Interp%j_lat)) deallocate(Interp%j_lat)

  end subroutine horiz_interp_bilinear_del
  ! </SUBROUTINE>

  !#######################################################################

  function indp (value, array)
    integer                        :: indp
    real, dimension(:), intent(in) :: array
    real, intent(in)               :: value
    !
    !=======================================================================
    !
    !     indp = index of nearest data point within "array" corresponding to
    !            "value".

    !     inputs:
    !     value  = arbitrary data...same units as elements in "array"
    !     array  = array of data points  (must be monotonically increasing)

    !     output:
    !     indp =  index of nearest data point to "value"
    !             if "value" is outside the domain of "array" then indp = 1
    !             or "ia" depending on whether array(1) or array(ia) is
    !             closest to "value"
    !=======================================================================
    !
    integer i, ia, unit
    logical keep_going
    !
    ia = size(array(:))
    do i=2,ia
       if (array(i) .lt. array(i-1)) then
          unit = stdout()
          write (unit,*) &
               ' => Error: array must be monotonically increasing in "indp"' , &
               '           when searching for nearest element to value=',value
          write (unit,*) '           array(i) < array(i-1) for i=',i 
          write (unit,*) '           array(i) for i=1..ia follows:'
          call abort()
       endif
    enddo
    if (value .lt. array(1) .or. value .gt. array(ia)) then
       if (value .lt. array(1))  indp = 1
       if (value .gt. array(ia)) indp = ia
    else
       i=1
       keep_going = .true.
       do while (i .le. ia .and. keep_going)
          i = i+1
          if (value .le. array(i)) then
             indp = i
             if (array(i)-value .gt. value-array(i-1)) indp = i-1
             keep_going = .false.
          endif
       enddo
    endif
    return
  end function indp

  !######################################################################

end module horiz_interp_bilinear_mod
