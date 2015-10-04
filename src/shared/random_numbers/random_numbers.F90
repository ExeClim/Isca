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

module random_numbers_mod 
  ! Generic module to wrap random number generators. 
  !   The module defines a type that identifies the particular stream of random 
  !   numbers, and has procedures for initializing it and getting real numbers 
  !   in the range 0 to 1. 
  ! This version uses the Mersenne Twister to generate random numbers on [0, 1]. 
  !
  use MersenneTwister_mod, only: randomNumberSequence, & ! The random number engine.
                                 new_RandomNumberSequence, getRandomReal 
  use time_manager_mod, only: time_type, get_date
  implicit none
  private
  
  type randomNumberStream
    type(randomNumberSequence) :: theNumbers
  end type randomNumberStream
  
  interface getRandomNumbers
    module procedure getRandomNumber_Scalar, getRandomNumber_1D, getRandomNumber_2D
  end interface getRandomNumbers
  
  interface initializeRandomNumberStream
    module procedure initializeRandomNumberStream_S, initializeRandomNumberStream_V
  end interface initializeRandomNumberStream

  public :: randomNumberStream,                             &
            initializeRandomNumberStream, getRandomNumbers, &
            constructSeed
contains
  ! ---------------------------------------------------------
  ! Initialization
  ! ---------------------------------------------------------
  function initializeRandomNumberStream_S(seed) result(new) 
    integer, intent( in)     :: seed
    type(randomNumberStream) :: new
    
    new%theNumbers = new_RandomNumberSequence(seed)
    
  end function initializeRandomNumberStream_S
  ! ---------------------------------------------------------
  function initializeRandomNumberStream_V(seed) result(new) 
    integer, dimension(:), intent( in) :: seed
    type(randomNumberStream)           :: new
    
    new%theNumbers = new_RandomNumberSequence(seed)
    
  end function initializeRandomNumberStream_V
  ! ---------------------------------------------------------
  ! Procedures for drawing random numbers
  ! ---------------------------------------------------------
  subroutine getRandomNumber_Scalar(stream, number)
    type(randomNumberStream), intent(inout) :: stream
    real,                     intent(  out) :: number
    
    number = getRandomReal(stream%theNumbers)
  end subroutine getRandomNumber_Scalar
  ! ---------------------------------------------------------
  subroutine getRandomNumber_1D(stream, numbers)
    type(randomNumberStream), intent(inout) :: stream
    real, dimension(:),       intent(  out) :: numbers
    
    ! Local variables
    integer :: i
    
    do i = 1, size(numbers)
      numbers(i) = getRandomReal(stream%theNumbers)
    end do
  end subroutine getRandomNumber_1D
  ! ---------------------------------------------------------
  subroutine getRandomNumber_2D(stream, numbers)
    type(randomNumberStream), intent(inout) :: stream
    real, dimension(:, :),    intent(  out) :: numbers
    
    ! Local variables
    integer :: i
    
    do i = 1, size(numbers, 2)
      call getRandomNumber_1D(stream, numbers(:, i))
    end do
  end subroutine getRandomNumber_2D
  ! ---------------------------------------------------------
  ! Constructs a unique seed from grid cell index and model date/time
  !   The perm is supplied we generate a different seed by 
  !   circularly shifting the bits of the seed - this is useful 
  !   if we want to create more than one seed for a given 
  !   column and model date/time. 
  !   Note that abs(perm) must be <= the number of bits used 
  !   to represent the default integer (likely 32) 
  ! ---------------------------------------------------------
  function constructSeed(i, j, time, perm) result(seed)
    integer,           intent( in)  :: i, j
    type(time_type),   intent( in) :: time
    integer, optional, intent( in) :: perm
    integer, dimension(8) :: seed
    
    ! Local variables
    integer :: year, month, day, hour, minute, second
    
    
    call get_date(time, year, month, day, hour, minute, second)
    seed = (/ i, j, year, month, day, hour, minute, second /)
    if(present(perm)) seed = ishftc(seed, perm)
  end function constructSeed
end module random_numbers_mod

