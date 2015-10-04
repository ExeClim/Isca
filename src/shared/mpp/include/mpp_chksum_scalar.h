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

    function MPP_CHKSUM_( var, pelist )
!mold is a dummy array to be used by TRANSFER()
!must be same TYPE as result
!result is LONG_KIND, which will actually be int ifdef no_8byte_integers
      integer(LONG_KIND) :: MPP_CHKSUM_
      MPP_TYPE_, intent(in) :: var
      integer, intent(in), optional :: pelist(:)
      integer(LONG_KIND) :: mold(1)
      pointer( p, mold )

      p = LOC(var)
      MPP_CHKSUM_ = mpp_chksum( mold, pelist )
      return
    end function MPP_CHKSUM_
