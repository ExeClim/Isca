module statistical_cloud_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif
  use            fms_mod, only: stdlog, FATAL, WARNING, NOTE, error_mesg, &
                                uppercase, check_nml_error
  use   time_manager_mod, only: time_type
  use sat_vapor_pres_mod, only: compute_qs, lookup_es
  use   diag_manager_mod, only: register_diag_field, send_data
  use      constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, HLV, KAPPA, TFREEZE, PI

  implicit none
  
  character(len=20), parameter :: mod_name = "statistical_cloud"
  integer, parameter  :: B_GAUSSIAN = 1
  integer, private    :: pdf_formula = B_GAUSSIAN
  character(len=32)   :: pdf_formula_name = 'gaussian'
  
  ! parameters to control the shape of sigma_s
  real :: s0 = 0.01
  real :: s1 = 0.09
  real :: s2 = 1.1    ! peak position
  real :: s3 = -0.31  ! left shift...

  namelist /statistical_cloud_nml/ pdf_formula_name, s0, s1, s2, s3

  contains

  subroutine statistical_cloud_init(axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    integer :: io, ierr, nml_unit, stdlog_unit
    character(len=32) :: pdf_str = ''
    
#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=statistical_cloud_nml, iostat=io)
    ierr = check_nml_error(io, 'statistical_cloud_nml')
#else
    if (file_exist('input.nml')) then
      nml_unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
        read(nml_unit, nml=statistical_cloud_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'statistical_cloud_nml')
      enddo
10    call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, statistical_cloud_nml)

    ! Select distribution formula
    pdf_str = uppercase(trim(pdf_formula_name))

    !if(pdf_str == 'GAUSSIAN' .or. pdf_str == 'GAUSS') then
    if (index(pdf_str, 'GAUSS') .ge. 0) then
      pdf_formula = B_GAUSSIAN
      call error_mesg(mod_name, 'Using default Gaussian probability density function.', NOTE)
    else
      call error_mesg(mod_name, '"' // trim(pdf_formula_name) // &
                '" is not a valid distribution choice.', FATAL)
    endif

  end subroutine statistical_cloud_init

  subroutine statistical_cloud(temp, sphum, pfull, ql, cf)
    real, intent(in),  dimension(:,:,:) :: temp, sphum, pfull
    real, intent(out), dimension(:,:,:) :: ql, cf
    real, dimension(size(temp,1),size(temp,2),size(temp,3)) :: qc
    
    select case(pdf_formula)
      case(B_GAUSSIAN)
        call gaussian_cloud(temp, sphum, pfull, ql, cf, qc)
        ! update new condensate
        !ql = qc
      case default
        call error_mesg(mod_name, 'invalid PDF formula.', FATAL)
    end select

  end subroutine statistical_cloud

  subroutine gaussian_cloud(temp, sphum, pfull, ql_in, cf, qc)
    ! Refer to the following paper:
    ! Mellor, The Gaussian cloud model relations, Jorunal of Atmospheric Sciences, 1977 (34)
    implicit none
    real, intent(in),  dimension(:,:,:) :: temp, sphum, pfull, ql_in
    real, intent(out), dimension(:,:,:) :: &
        cf,     & ! layer mean cloud fraction (fraction)
        qc        ! layer mean cloud condensate (kg condensate/kg air)
    real, dimension(size(temp,1), size(temp,2), size(temp,3)) :: &
        Tl,     & ! liquid water temperature
        qsl,    & ! saturated specific humidity
        dqsldT, & !
        deltaq, & ! saturation excess (kg water / kg air)
        sigmas, & ! standardard devations of water perturbation 's' (kg water / kg air)
        acoef,  & ! thermo coefficient = 1./(1.+ L*dqsldT/cp)
        Q1
    integer :: i, j, k

    Tl = temp - HLV / CP_AIR * ql_in
    !call compute_qs(Tl, pfull, qsl)
    !dqsldT = HLV * qsl / RVGAS * Tl**2
    call compute_qs(Tl, pfull, qsl, dqsdT=dqsldT) ! esat = esl

    acoef = 1.0 / (1.0 + HLV / dqsldT / CP_AIR)

    deltaq = sphum - qsl

    call calc_sigmas(pfull, sigmas)

    Q1 = acoef * deltaq / sigmas / 2.0

    ! Note that erfc(x) = 1 - erf(x) and erfcc is a concise approximation of erfc
    do i=1,size(temp,1)
      do j=1,size(temp,2)
        do k=1,size(temp,3)
          cf(i,j,k) = max(0.5 * (2.0 - erfcc( Q1(i,j,k) / sqrt(2.0) )), 0.0)
          ! qc(i,j,k) = cf(i,j,k) * Q1(i,j,k) + exp( -0.5 * Q1(i,j,k)**2 ) / sqrt(2.0*PI)
          ! qc(i,j,k) = 2.0 * sigmas(i,j,k) * max(qc(i,j,k), 0.0)  
          !if (i==1 .and. j==1) write(*,*) 'QL cf=', cf(i,j,k), ql_in(i,j,k), qc(i,j,k)
        enddo
      enddo
    enddo
    qc = cf * Q1 + exp( -0.5 * Q1**2 ) / sqrt(2.0*PI)
    qc = 2.0 * sigmas * max(qc, 0.0)

  end subroutine gaussian_cloud

  subroutine calc_sigmas(pfull, sigmas)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(out), dimension(:,:,:) :: sigmas
    real, dimension(size(pfull,1), size(pfull,2), size(pfull,3)) :: s
    real :: p0 = 1.0e5 ! Units: Pa

    s = pfull / p0
    sigmas = max(0.0, ( s0*(1.0-s) + s1*s + s*(1.0-s)*(s2*s+s3) ) * 1.0e-3)

  end subroutine calc_sigmas

  subroutine statistical_cloud_end()
    ! nothing yet
  end subroutine statistical_cloud_end

  function erfcc(x)
    ! This numerical recipes routine calculates the complementary error function.
    ! Original code is from src/atmos_param/edt/edt.F90, lines 4753-4795
    !
    ! Returns the complementary error function erfc(x) = 1-erf(x) with fractional
    ! error everywhere less than 1.2e-7
    implicit none

    real :: erfcc
    real, intent(in) :: x 
    real :: t, z

    z = abs(x)      
    t = 1.0 / (1.0 + 0.5 * z)

    erfcc = t * exp( -z * z   -1.26551223 &
                      + t * (  1.00002368 &
                      + t * (  0.37409196 &
                      + t * (  0.09678418 &
                      + t * ( -0.18628806 &
                      + t * (  0.27886807 &
                      + t * ( -1.13520398 &
                      + t * (  1.48851587 &
                      + t * ( -0.82215223 &
                      + t *    0.17087277 ) ) ) ) ) ) ) ) )

    if (x .lt. 0.0) erfcc = 2.0 - erfcc    

  end function erfcc

  ! subroutine log_normal_pdf(x, mu, sigma, pdf)
  !   ! Modified from https://people.sc.fsu.edu/~jburkardt/f_src/log_normal/log_normal.f90

  !   !*****************************************************************************80
  !   !
  !   !! LOG_NORMAL_PDF evaluates the Log Normal PDF.
  !   !
  !   !  Discussion:
  !   !
  !   !    PDF(A,B;X)
  !   !      = exp ( - 0.5 * ( ( log ( X ) - MU ) / SIGMA )^2 )
  !   !        / ( SIGMA * X * sqrt ( 2 * PI ) )
  !   !
  !   !    The Log Normal PDF is also known as the Cobb-Douglas PDF,
  !   !    and as the Antilog_normal PDF.
  !   !
  !   !    The Log Normal PDF describes a variable X whose logarithm
  !   !    is normally distributed.
  !   !
  !   !    The special case MU = 0, SIGMA = 1 is known as Gilbrat's PDF.
  !   !
  !   !  Licensing:
  !   !    This code is distributed under the GNU LGPL license.
  !   !
  !   !  Modified:
  !   !    10 February 1999
  !   !
  !   !  Author:
  !   !    John Burkardt
  !   !
  !   !  Parameters:
  !   !    Input, real X, the argument of the PDF.  (0.0 < X)
  !   !    Input, real MU, SIGMA, the parameters of the PDF. (0.0 < SIGMA)
  !   !    Output, real PDF, the value of the PDF.
    
  !     implicit none
    
  !     real :: mu, pdf, sigma, x
    
  !     if (x <= 0.0) then
  !       pdf = 0.0
  !     else
  !       pdf = exp(-0.5 * ( (log(x) - mu) / sigma) ** 2) &
  !         / ( sigma * x * sqrt (2.0 * PI) )
  !     end if

  ! end subroutine log_normal_pdf
   
end module statistical_cloud_mod