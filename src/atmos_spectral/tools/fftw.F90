  module FFTW3
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

    logical :: module_is_initialized=.false.

contains

!
! Real -> Complex
!
subroutine grid_to_fourier_double_2d_fftw(num, leng, lenc, grid, fourier)

integer(kind=8),            intent(in)    :: num    ! number of transformations
integer(kind=4),            intent(in)    :: leng   ! length of each transformation
real(C_DOUBLE),             intent(in)    :: grid(leng, num)    ! input grid
complex(C_DOUBLE_COMPLEX),  intent(out)   :: fourier(lenc, num) ! output fourier domain

real    :: fact ! factor by which to scale results
integer :: i, j

complex(C_DOUBLE_COMPLEX) :: aout
real(C_DOUBLE)            :: ina, inb, ina_copy
type(C_PTR)               :: p, p2
parameter (N=10)
dimension ina(leng), inb(leng), ina_copy(leng), aout((leng/2+1))

call dfftw_plan_dft_r2c_1d(p, leng, ina, aout, FFTW_ESTIMATE)

do j = 1, num

  do i=1, leng
        ina(i) = grid(i, j)
        ina_copy(i) = grid(i, j)
  enddo

  call dfftw_execute_dft_r2c(p, ina, aout)

  ! do i=1, leng
  !       print *, 'ina_copy: ', ina_copy(i)
  ! enddo 

  do i = 1, lenc
    fourier(i, j) = aout(i)
  enddo

enddo

! call exit

return
end subroutine grid_to_fourier_double_2d_fftw

! 
! Complex -> Real
!
subroutine fourier_to_grid_double_2d_fftw(num, leng, lenc, fourier, grid)

integer(kind=8),            intent(in)    :: num    ! number of transformations
integer(kind=4),            intent(in)    :: leng   ! length of each transformation
complex(C_DOUBLE_COMPLEX),  intent(in)    :: fourier(lenc, num) ! output fourier domain
real(C_DOUBLE),             intent(out)   :: grid(leng, num)    ! input grid

real :: fact ! factor by which to scale results

complex(C_DOUBLE_COMPLEX) :: aout, aout_copy
real(C_DOUBLE)            :: ina, inb
type(C_PTR)               :: p, p2
parameter (N=10)
dimension ina(leng), inb(leng), ina_copy(leng), aout((leng/2+1)), aout_copy((leng/2+1))
integer :: i, j

fact = 1.0 / leng

call dfftw_plan_dft_c2r_1d(p2, leng, aout, inb, FFTW_ESTIMATE)


do j=1, num
  do i=1, lenc
    aout(i) = fourier(i, j)
    aout_copy(i) = fourier(i, j)
  enddo

  call dfftw_execute_dft_c2r(p2, aout, inb)

  do i = 1, leng
    inb(i) = inb(i)*fact
    grid(i, j) = inb(i)
  enddo

  ! do i=1,lenc
  !   print *, 'four: ', (inb(i))!, 'inb: ',inb(i)
  ! enddo
enddo

! call exit
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
    print *, 'leng, lenc', leng, lenc
    module_is_initialized = .true.
end subroutine fftw3_init


subroutine fftw3_end()
  call fftw_cleanup()
  module_is_initialized = .false.
end subroutine fftw3_end



end module FFTW3