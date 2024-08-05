! The original code is from http://romps.berkeley.edu/papers/pubdata/2016/lcl/lcl.f90
!
! Version 1.0 released by David Romps on September 12, 2017.
! 
! When using this code, please cite:
! 
! @article{Romps2017lcl,
!   Title   = {Exact expression for the lifting condensation level},
!   Author  = {David M. Romps},
!   Journal = {Journal of the Atmospheric Sciences},
!   Year    = {2017},
!   Volume  = {74},
! }
!
! This lcl function returns the height of the lifting condensation level
! (LCL) in meters.  The inputs are:
! - p in Pascals
! - T in Kelvins
! - Exactly one of rh, rhl, and rhs (dimensionless, from 0 to 1):
!    * The value of rh is interpreted to be the relative humidity with
!      respect to liquid water if T >= 273.15 K and with respect to ice if
!      T < 273.15 K. 
!    * The value of rhl is interpreted to be the relative humidity with
!      respect to liquid water
!    * The value of rhs is interpreted to be the relative humidity with
!      respect to ice
! - ldl is an optional logical flag.  If true, the lifting deposition
!   level (LDL) is returned instead of the LCL. 
! - min_lcl_ldl is an optional logical flag.  If true, the minimum of the
!   LCL and LDL is returned.

module lcl_mod

contains

function lcl(p,T,rh,rhl,rhs,return_ldl,return_min_lcl_ldl)

   implicit none

   ! In
   double precision, intent(in) :: p, T
   double precision, intent(in), optional :: rh, rhl, rhs
   logical, intent(in), optional :: return_ldl, return_min_lcl_ldl

   ! Out
   double precision :: lcl

   logical :: liquid2, solid2
   integer(kind=4) :: error

   double precision :: Ttrip, ptrip, E0v, E0s, ggr, rgasa, &
      rgasv, cva, cvv, cvl, cvs, cpa, cpv, al, bl, cl, as, bs, cs, &
      pv, qv, cpm, rgasm, rh2, rhl2, rhs2, ldl
   integer :: rh_counter
   logical :: return_ldl2, return_min_lcl_ldl2

   ! Set defaults for options
   return_ldl2         = .false.
   return_min_lcl_ldl2 = .false.

   if (present(return_ldl)) return_ldl2 = return_ldl
   if (present(return_min_lcl_ldl)) return_min_lcl_ldl2 = return_min_lcl_ldl

   ! Parameters
   Ttrip = 273.16     ! K
   ptrip = 611.65     ! Pa
   E0v   = 2.3740e6   ! J/kg
   E0s   = 0.3337e6   ! J/kg
   ggr   = 9.81       ! m/s^2
   rgasa = 287.04     ! J/kg/K 
   rgasv = 461.       ! J/kg/K 
   cva   = 719.       ! J/kg/K
   cvv   = 1418.      ! J/kg/K 
   cvl   = 4119.      ! J/kg/K 
   cvs   = 1861.      ! J/kg/K 
   cpa   = cva + rgasa
   cpv   = cvv + rgasv

   ! Ensure that, at most, one of liquid and solid are true
   if (return_ldl2 .and. return_min_lcl_ldl2) then
      print *,'Error in lcl: At most, only one of return_ldl and return_min_lcl_ldl can be true'
      stop
   end if

   ! Calculate pv from rh, rhl, or rhs
   rh_counter = 0
   if (present(rh )) rh_counter = rh_counter + 1
   if (present(rhl)) rh_counter = rh_counter + 1
   if (present(rhs)) rh_counter = rh_counter + 1
   if (rh_counter .ne. 1) then
      print *,'Error in lcl: Exactly one of rh, rhl, and rhs must be specified'
      stop
   end if

   if (present(rh)) then
      rh2 = rh
      ! The variable rh is assumed to be 
      ! with respect to liquid if T > Ttrip and 
      ! with respect to solid if T < Ttrip
      if (T > Ttrip) then
         pv = rh2 * pvstarl(T)
      else
         pv = rh2 * pvstars(T)
      end if
      rhl2 = pv / pvstarl(T)
      rhs2 = pv / pvstars(T)
   else if (present(rhl)) then
      rhl2 = rhl
      pv = rhl2 * pvstarl(T)
      rhs2 = pv / pvstars(T)
      if (T > Ttrip) then
         rh2 = rhl2
      else
         rh2 = rhs2
      end if
   else if (present(rhs)) then
      rhs2 = rhs
      pv = rhs2 * pvstars(T)
      rhl2 = pv / pvstarl(T)
      if (T > Ttrip) then
         rh2 = rhl2
      else
         rh2 = rhs2
      end if
   end if
   if (pv > p) then
      print *,'Warning in lcl: pv exceeds p, returning -9999'
      lcl = -9999.
      return
   end if

   ! Calculate lcl_liquid and lcl_solid
   qv = rgasa*pv / (rgasv*p + (rgasa-rgasv)*pv)
   rgasm = (1.-qv)*rgasa + qv*rgasv
   cpm = (1.-qv)*cpa + qv*cpv
   if (rh2==0) then
      lcl = cpm*T/ggr
      return
   end if
   al = -(cpv-cvl)/rgasv + cpm/rgasm
   bl = -(E0v-(cvv-cvl)*Ttrip)/(rgasv*T)
   !cl = pv/pvstarl(T)*exp(-(E0v-(cvv-cvl)*Ttrip)/(rgasv*T))
   cl = pv/pvstarl(T)*exp(bl)

   as = -(cpv-cvs)/rgasv + cpm/rgasm
   bs = -(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T)
   !cs = pv/pvstars(T)*exp(-(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T))
   cs = pv/pvstars(T)*exp(bs)

   lcl = cpm*T/ggr*( 1. - bl/(al*wapr(bl/al*cl**(1./al),-1,error,0)) )
   ldl = cpm*T/ggr*( 1. - bs/(as*wapr(bs/as*cs**(1./as),-1,error,0)) )

   ! Return either the LCL or the LDL
   if (return_ldl2 .and. return_min_lcl_ldl2) then
      print *,'Error in lcl: return_ldl and return_min_lcl_ldl cannot both be true'
      stop
   else if (return_ldl2) then
      lcl = ldl
   else if (return_min_lcl_ldl2) then
      lcl = min(lcl,ldl)
   end if

   contains

      ! The saturation vapor pressure over liquid water
      function pvstarl(T)

        implicit none

        double precision :: T        ! In
        double precision :: pvstarl  ! Out
  
        pvstarl = ptrip * (T/Ttrip)**((cpv-cvl)/rgasv) * &
          exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1./Ttrip - 1./T) )
      
      end function pvstarl
      
      ! The saturation vapor pressure over solid ice
      function pvstars(T)
         
        implicit none

        double precision :: T        ! In
        double precision :: pvstars  ! Out
  
        pvstars = ptrip * (T/Ttrip)**((cpv-cvs)/rgasv) * &
          exp( (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1./Ttrip - 1./T) )
    
      end function pvstars

end function lcl

function bisect ( xx, nb, ner, l )

!*****************************************************************************80
!
!! BISECT approximates the W function using bisection.
!
!  Discussion:
!
!    The parameter TOL, which determines the accuracy of the bisection
!    method, is calculated using NBITS (assuming the final bit is lost
!    due to rounding error).
!
!    N0 is the maximum number of iterations used in the bisection
!    method.
!
!    For XX close to 0 for Wp, the exponential approximation is used.
!    The approximation is exact to O(XX^8) so, depending on the value
!    of NBITS, the range of application of this formula varies. Outside
!    this range, the usual bisection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Parameters:
!
!    Input, double precision :: XX, the argument.
!
!    Input, integer ( kind = 4 ) NB, indicates the branch of the W function.
!    0, the upper branch;
!    nonzero, the lower branch.
!
!    Output, integer ( kind = 4 ) NER, the error flag.
!    0, success;
!    1, the routine did not converge.  Perhaps reduce NBITS and try again.
!
!    Input, integer ( kind = 4 ) L, the offset indicator.
!    1, XX represents the offset of the argument from -exp(-1).
!    not 1, XX is the actual argument.
!
!    Output, double precision :: BISECT, the value of W(X), as determined
!
  implicit none

  double precision :: bisect
!  double precision :: crude
  double precision :: d
  double precision :: f
  double precision :: fd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n0
  parameter ( n0 = 500 )
  integer ( kind = 4 ) nb
  integer ( kind = 4 ), save :: nbits = 0
  integer ( kind = 4 ) ner
  double precision :: r
  double precision :: test
  double precision :: tol
  double precision :: u
  double precision :: x
  double precision :: xx

  bisect = 0.0D+00
  ner = 0

  if ( nbits == 0 ) then
    call nbits_compute ( nbits )
  end if

  if ( l == 1 ) then
    x = xx - exp ( -1.0D+00 )
  else
    x = xx
  end if

  if ( nb == 0 ) then

    test = 1.0D+00 / ( 2.0D+00 ** nbits ) ** ( 1.0D+00 / 7.0D+00 )

    if ( abs ( x ) < test ) then

      bisect = x &
        * exp ( - x &
        * exp ( - x &
        * exp ( - x &
        * exp ( - x &
        * exp ( - x &
        * exp ( - x ))))))

      return

    else

      u = crude ( x, nb ) + 1.0D-03
      tol = abs ( u ) / 2.0D+00 ** nbits
      d = max ( u - 2.0D-03, -1.0D+00 )

      do i = 1, n0

        r = 0.5D+00 * ( u - d )
        bisect = d + r
!
!  Find root using w*exp(w)-x to avoid ln(0) error.
!
        if ( x < exp ( 1.0D+00 ) ) then

          f = bisect * exp ( bisect ) - x
          fd = d * exp ( d ) - x
!
!  Find root using ln(w/x)+w to avoid overflow error.
!
        else

          f = log ( bisect / x ) + bisect
          fd = log ( d / x ) + d

        end if

        if ( f == 0.0D+00 ) then
          return
        end if

        if ( abs ( r ) <= tol ) then
          return
        end if

        if ( 0.0D+00 < fd * f ) then
          d = bisect
        else
          u = bisect
        end if

      end do

    end if

  else

    d = crude ( x, nb ) - 1.0D-03
    u = min ( d + 2.0D-03, -1.0D+00 )
    tol = abs ( u ) / 2.0D+00 ** nbits

    do i = 1, n0

      r = 0.5D+00 * ( u - d )
      bisect = d + r
      f = bisect * exp ( bisect ) - x

      if ( f == 0.0D+00 ) then
        return
      end if

      if ( abs ( r ) <= tol ) then
        return
      end if

      fd = d * exp ( d ) - x

      if ( 0.0D+00 < fd * f ) then
        d = bisect
      else
        u = bisect
      end if

    end do

  end if
!
!  The iteration did not converge.
!
  ner = 1

  return
end function bisect

function crude ( xx, nb )

!*****************************************************************************80
!
!! CRUDE returns a crude approximation for the W function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Parameters:
!
!    Input, double precision :: XX, the argument.
!
!    Input, integer ( kind = 4 ) NB, indicates the desired branch.
!    * 0, the upper branch;
!    * nonzero, the lower branch.
!
!    Output, double precision :: CRUDE, the crude approximation to W at XX.
!
  implicit none

  double precision :: an2
  double precision :: c13
  double precision :: crude
  double precision :: em
  double precision :: em2
  double precision :: em9
  double precision :: eta
  integer ( kind = 4 ) init
  integer ( kind = 4 ) nb
  double precision :: reta
  double precision :: s2
  double precision :: s21
  double precision :: s22
  double precision :: s23
  double precision :: t
  double precision :: ts
  double precision :: xx
  double precision :: zl

  save c13
  save em
  save em2
  save em9
  save init
  save s2
  save s21
  save s22
  save s23

  data init / 0 /

  crude = 0.0D+00
!
!  Various mathematical constants.
!
  if ( init == 0 ) then
    init = 1
    em = - exp ( -1.0D+00 )
    em9 = - exp ( -9.0D+00 )
    c13 = 1.0D+00 / 3.0D+00
    em2 = 2.0D+00 / em
    s2 = sqrt ( 2.0D+00 )
    s21 = 2.0D+00 * s2 - 3.0D+00
    s22 = 4.0D+00 - 3.0D+00 * s2
    s23 = s2 - 2.0D+00
  end if
!
!  Crude Wp.
!
  if ( nb == 0 ) then

    if ( xx <= 20.0D+00 ) then
      reta = s2 * sqrt ( 1.0D+00 - xx / em )
      an2 = 4.612634277343749D+00 * sqrt ( sqrt ( reta + &
        1.09556884765625D+00 ) )
      crude = reta / ( 1.0D+00 + reta / ( 3.0D+00 &
        + ( s21 * an2 + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0D+00
    else
      zl = log ( xx )
      crude = log ( xx / log ( xx &
        / zl ** exp ( -1.124491989777808D+00 / &
        ( 0.4225028202459761D+00 + zl ))))
    end if

  else
!
!  Crude Wm.
!
    if ( xx <= em9 ) then
      zl = log ( -xx )
      t = -1.0D+00 - zl
      ts = sqrt ( t )
      crude = zl - ( 2.0D+00 * ts ) / ( s2 + ( c13 - t &
        / ( 2.7D+02 + ts * 127.0471381349219D+00 ) ) * ts )
    else
      zl = log ( -xx )
      eta = 2.0D+00 - em2 * xx
      crude = log ( xx / log ( - xx / ( ( 1.0D+00 &
        - 0.5043921323068457D+00 * ( zl + 1.0D+00 ) ) &
        * ( sqrt ( eta ) + eta / 3.0D+00 ) + 1.0D+00 ) ) )
     end if

  end if

  return
end function crude

subroutine nbits_compute ( nbits )

!*****************************************************************************80
!
!! NBITS_COMPUTE computes the mantissa length minus one.
!
!  Discussion:
!
!    NBITS is the number of bits (less 1) in the mantissa of the
!    floating point number number representation of your machine.
!    It is used to determine the level of accuracy to which the W
!    function should be calculated.
!
!    Most machines use a 24-bit matissa for single precision and
!    53-56 bits for double precision ::. The IEEE standard is 53
!    bits. The Fujitsu VP2200 uses 56 bits. Long word length
!    machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
!    single precision.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NBITS, the mantissa length, in bits, 
!    minus one.
!
  implicit none

  double precision :: b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbits
  double precision :: v

  nbits = 0

  b = 1.0D+00

  do

    b = b / 2.0D+00
    v = b + 1.0D+00

    if ( v == 1.0D+00 ) then
      return
    end if

    nbits = nbits + 1

  end do

  return
end subroutine nbits_compute

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp

function wapr ( x, nb, nerror, l )

!*****************************************************************************80
!
!! WAPR approximates the W function.
!
!  Discussion:
!
!    The call will fail if the input value X is out of range.
!    The range requirement for the upper branch is:
!      -exp(-1) <= X.
!    The range requirement for the lower branch is:
!      -exp(-1) < X < 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Parameters:
!
!    Input, double precision :: X, the argument.
!
!    Input, integer ( kind = 4 ) NB, indicates the desired branch.
!    * 0, the upper branch;
!    * nonzero, the lower branch.
!
!    Output, integer ( kind = 4 ) NERROR, the error flag.
!    * 0, successful call.
!    * 1, failure, the input X is out of range.
!
!    Input, integer ( kind = 4 ) L, indicates the interpretation of X.
!    * 1, X is actually the offset from -(exp-1), so compute W(X-exp(-1)).
!    * not 1, X is the argument; compute W(X);
!
!    Output, double precision :: WAPR, the approximate value of W(X).
!
  implicit none

  double precision :: an2
  double precision :: an3
  double precision :: an4
  double precision :: an5
  double precision :: an6
  double precision :: c13
  double precision :: c23
  double precision :: d12
  double precision :: delx
  double precision :: em
  double precision :: em2
  double precision :: em9
  double precision :: eta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) init
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nbits
  integer ( kind = 4 ) nerror
  integer ( kind = 4 ) niter
  double precision :: reta
  double precision :: s2
  double precision :: s21
  double precision :: s22
  double precision :: s23
  double precision :: t
  double precision :: tb
  double precision :: tb2
  double precision :: temp
  double precision :: temp2
  double precision :: ts
  double precision :: wapr
  double precision :: x
  double precision :: x0
  double precision :: x1
  double precision :: xx
  double precision :: zl
  double precision :: zn

  save an3
  save an4
  save an5
  save an6
  save c13
  save c23
  save d12
  save em
  save em2
  save em9
  save init
  save nbits
  save niter
  save s2
  save s21
  save s22
  save s23
  save tb
  save tb2
  save x0
  save x1

  data init / 0 /
  data niter / 1 /

  wapr = 0.0D+00
  nerror = 0

  if ( init == 0 ) then

    init = 1

    call nbits_compute ( nbits )

    if ( 56 <= nbits ) then
      niter = 2
    end if
!
!  Various mathematical constants.
!
    em = -exp ( -1.0D+00 )
    em9 = -exp ( -9.0D+00 )
    c13 = 1.0D+00 / 3.0D+00
    c23 = 2.0D+00 * c13
    em2 = 2.0D+00 / em
    d12 = -em2
    tb = 0.5D+00 ** nbits
    tb2 = sqrt ( tb )
    x0 = tb ** ( 1.0D+00 / 6.0D+00 ) * 0.5D+00
    x1 = ( 1.0D+00 - 17.0D+00 * tb ** ( 2.0D+00 / 7.0D+00 ) ) * em
    an3 = 8.0D+00 / 3.0D+00
    an4 = 135.0D+00 / 83.0D+00
    an5 = 166.0D+00 / 39.0D+00
    an6 = 3167.0D+00 / 3549.0D+00
    s2 = sqrt ( 2.0D+00 )
    s21 = 2.0D+00 * s2 - 3.0D+00
    s22 = 4.0D+00 - 3.0D+00 * s2
    s23 = s2 - 2.0D+00

  end if

  if ( l == 1 ) then

    delx = x

    if ( delx < 0.0D+00 ) then
      nerror = 1
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'WAPR - Fatal error!'
      write ( *, '(a)' ) '  The offset X is negative.'
      write ( *, '(a)' ) '  It must be nonnegative.'
      stop 1
    end if

    xx = x + em

  else

    if ( x < em ) then
      nerror = 1
      return
    else if ( x == em ) then
      wapr = -1.0D+00
      return
    end if

    xx = x
    delx = xx - em

  end if

  if ( nb == 0 ) then
!
!  Calculations for Wp.
!
    if ( abs ( xx ) <= x0 ) then
      wapr = xx / ( 1.0D+00 + xx / ( 1.0D+00 + xx &
        / ( 2.0D+00 + xx / ( 0.6D+00 + 0.34D+00 * xx ))))
      return
    else if ( xx <= x1 ) then
      reta = sqrt ( d12 * delx )
      wapr = reta / ( 1.0D+00 + reta / ( 3.0D+00 + reta / ( reta &
        / ( an4 + reta / ( reta * an6 + an5 ) ) + an3 ) ) ) &
        - 1.0D+00
      return
    else if ( xx <= 20.0D+00 ) then
      reta = s2 * sqrt ( 1.0D+00 - xx / em )
      an2 = 4.612634277343749D+00 * sqrt ( sqrt ( reta + &
        1.09556884765625D+00 ))
      wapr = reta / ( 1.0D+00 + reta / ( 3.0D+00 + ( s21 * an2 &
        + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0D+00
    else
      zl = log ( xx )
      wapr = log ( xx / log ( xx &
        / zl ** exp ( -1.124491989777808D+00 / &
        ( 0.4225028202459761D+00 + zl ))))
    end if
!
!  Calculations for Wm.
!
  else

    if ( 0.0D+00 <= xx ) then
      nerror = 1
      return
    else if ( xx <= x1 ) then
      reta = sqrt ( d12 * delx )
      wapr = reta / ( reta / ( 3.0D+00 + reta / ( reta / ( an4 &
        + reta / ( reta * an6 - an5 ) ) - an3 ) ) - 1.0D+00 ) - 1.0D+00
      return
    else if ( xx <= em9 ) then
      zl = log ( -xx )
      t = -1.0D+00 - zl
      ts = sqrt ( t )
      wapr = zl - ( 2.0D+00 * ts ) / ( s2 + ( c13 - t &
        / ( 270.0D+00 + ts * 127.0471381349219D+00 )) * ts )
    else
      zl = log ( -xx )
      eta = 2.0D+00 - em2 * xx
      wapr = log ( xx / log ( -xx / ( ( 1.0D+00 &
        - 0.5043921323068457D+00 * ( zl + 1.0D+00 ) ) &
        * ( sqrt ( eta ) + eta / 3.0D+00 ) + 1.0D+00 )))
    end if

  end if

  do i = 1, niter
    zn = log ( xx / wapr ) - wapr
    temp = 1.0D+00 + wapr
    temp2 = temp + c23 * zn
    temp2 = 2.0D+00 * temp * temp2
    wapr = wapr * ( 1.0D+00 + ( zn / temp ) * ( temp2 - zn ) &
      / ( temp2 - 2.0D+00 * zn ) )
  end do

  return
end function wapr

end module lcl_mod
