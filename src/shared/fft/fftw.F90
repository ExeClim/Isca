module fftw3
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
    public :: fftw3_init
    public :: fftw3_end
    public :: grid_to_fourier_fftw
    public :: fourier_to_grid_fftw

    interface grid_to_fourier_fftw
      module procedure grid_to_fourier_double_2d_fftw
      module procedure grid_to_fourier_float_2d_fftw
    end interface

    interface fourier_to_grid_fftw 
      module procedure fourier_to_grid_double_2d_fftw
      module procedure fourier_to_grid_float_2d_fftw
    end interface

    private

    ! real -> complex
    real(C_DOUBLE), pointer  :: real_input(:)
    real(C_DOUBLE), pointer  :: real_output(:)
    type(C_PTR) :: real_input_pointer, real_output_pointer

    ! complex -> real
    complex(C_DOUBLE_COMPLEX), pointer  :: complex_input(:)
    complex(C_DOUBLE_COMPLEX), pointer  :: complex_output(:)
    type(C_PTR) :: complex_input_pointer, complex_output_pointer

    logical     :: module_is_initialized = .false.

contains

!
! Real -> Complex
!
subroutine grid_to_fourier_double_2d_fftw(num, leng, lenc, grid, fourier)

integer(kind=4),            intent(in)    :: num    ! number of transformations
integer(kind=4),            intent(in)    :: leng   ! length of each transformation
real(C_DOUBLE),             intent(in)    :: grid(leng, num)    ! input grid
complex(C_DOUBLE_COMPLEX),  intent(out)   :: fourier(lenc, num) ! output fourier domain
real                                      :: fact ! factor by which to scale results
integer                                   :: i, j

fact = 1.0 / (leng - 1)

do j = 1, num
  do i = 1, leng - 1
    real_input(i) = grid(i,j)
  enddo

  call dfftw_execute_dft_r2c(real_input_pointer, real_input, complex_output)

  do i = 1, lenc
    fourier(i, j) = complex_output(i) * fact
  enddo
enddo
return
end subroutine grid_to_fourier_double_2d_fftw


! 
! Complex -> Real
!
subroutine fourier_to_grid_double_2d_fftw(num, leng, lenc, fourier, grid)

integer(kind=4),            intent(in)    :: num    ! number of transformations
integer(kind=4),            intent(in)    :: leng   ! length of each transformation
complex(C_DOUBLE_COMPLEX),  intent(in)    :: fourier(lenc, num) ! output fourier domain
real(C_DOUBLE),             intent(out)   :: grid(leng, num)    ! input grid

integer :: i, j

do j = 1, num
  do i = 1, lenc
    complex_input(i) = fourier(i, j)
  enddo

  call dfftw_execute_dft_c2r(complex_input_pointer, complex_input, real_output)

  do i = 1, leng - 1
    grid(i, j) = real_output(i)
  enddo

enddo
return

end subroutine fourier_to_grid_double_2d_fftw


subroutine fourier_to_grid_float_2d_fftw()
  !! TODO
end subroutine fourier_to_grid_float_2d_fftw

subroutine grid_to_fourier_float_2d_fftw()
  !! TODO
end subroutine grid_to_fourier_float_2d_fftw

subroutine fftw3_init(leng, lenc)
    integer, intent(in) :: leng, lenc

    ! setup real -> complex transform
    real_input_pointer = fftw_alloc_real(int(leng, C_SIZE_T))
    call c_f_pointer(real_input_pointer, real_input, [leng])

    complex_output_pointer = fftw_alloc_complex(int(lenc, C_SIZE_T))
    call c_f_pointer(complex_output_pointer, complex_output, [leng])

    call dfftw_plan_dft_r2c_1d(real_input_pointer, leng, real_input, complex_output, FFTW_MEASURE)

    ! setup complex -> real transform
    complex_input_pointer = fftw_alloc_complex(int((lenc), C_SIZE_T))
    call c_f_pointer(complex_input_pointer, complex_input, [lenc])

    real_output_pointer = fftw_alloc_real(int(leng, C_SIZE_T))
    call c_f_pointer(real_output_pointer, real_output, [leng])

    call dfftw_plan_dft_c2r_1d(complex_input_pointer, leng, complex_input, real_output, FFTW_MEASURE)

    module_is_initialized = .true.
end subroutine fftw3_init


subroutine fftw3_end()
  call fftw_cleanup()
  module_is_initialized = .false.
end subroutine fftw3_end


end module fftw3