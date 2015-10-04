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

! -*-f90-*-
! $Id: mpp_sum.inc,v 18.0 2010/03/02 23:57:50 fms Exp $

!#######################################################################

    subroutine MPP_SUM_SCALAR_( a, pelist )
!sums array a when only first element is passed: this routine just converts to a call to MPP_SUM_
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call MPP_SUM_( b, 1, pelist )
      a = b(1)
      return
    end subroutine MPP_SUM_SCALAR_

!#######################################################################
    subroutine MPP_SUM_2D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine MPP_SUM_2D_

!#######################################################################
    subroutine MPP_SUM_3D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine MPP_SUM_3D_

!#######################################################################
    subroutine MPP_SUM_4D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine MPP_SUM_4D_

!#######################################################################
    subroutine MPP_SUM_5D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine MPP_SUM_5D_
