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

  subroutine MPP_GLOBAL_SUM_TL_( domain, field, field_tl, gsum, gsum_tl, flags, position, tile_count )
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(inout) :: field(:,: MPP_EXTRA_INDICES_ )
    MPP_TYPE_, intent(inout) :: field_tl(:,: MPP_EXTRA_INDICES_ )
    MPP_TYPE_, intent(inout) :: gsum
    MPP_TYPE_, intent(inout) :: gsum_tl
    integer, intent(in), optional :: position
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: tile_count

    gsum = mpp_global_sum(domain, field, flags, position, tile_count )
    gsum_tl = mpp_global_sum(domain, field_tl, flags, position, tile_count )

    return
  end subroutine MPP_GLOBAL_SUM_TL_
