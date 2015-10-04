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
! $Id: mpp_transmit.inc,v 19.0 2012/01/06 22:04:03 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_TRANSMIT_SCALAR_( put_data, to_pe, get_data, from_pe, plen, glen, block, tag, recv_request, send_request)
      integer, intent(in) :: to_pe, from_pe
      MPP_TYPE_, intent(in)  :: put_data
      MPP_TYPE_, intent(out) :: get_data
      integer, optional,  intent(in) :: plen, glen
      logical, intent(in),  optional :: block
      integer, intent(in),  optional :: tag
      integer, intent(out), optional :: recv_request, send_request
      integer                       :: put_len, get_len
      MPP_TYPE_ :: put_data1D(1), get_data1D(1)
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      get_len=1; if(PRESENT(glen))get_len=glen
      call MPP_TRANSMIT_ ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                           recv_request=recv_request, send_request=send_request )

      return
    end subroutine MPP_TRANSMIT_SCALAR_

    subroutine MPP_TRANSMIT_2D_( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: recv_request, send_request
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine MPP_TRANSMIT_2D_

    subroutine MPP_TRANSMIT_3D_( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: recv_request, send_request
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine MPP_TRANSMIT_3D_

    subroutine MPP_TRANSMIT_4D_( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: recv_request, send_request
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine MPP_TRANSMIT_4D_

    subroutine MPP_TRANSMIT_5D_( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: recv_request, send_request
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine MPP_TRANSMIT_5D_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_RECV_( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request

      MPP_TYPE_ :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine MPP_RECV_

    subroutine MPP_SEND_( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(*)
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request
      MPP_TYPE_ :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag=tag, send_request=request )
    end subroutine MPP_SEND_

    subroutine MPP_RECV_SCALAR_( get_data, from_pe, glen, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      MPP_TYPE_, intent(out) :: get_data
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request

      integer, optional, intent(in) :: glen
      integer                       :: get_len
      MPP_TYPE_ :: get_data1D(1)
      MPP_TYPE_ :: dummy(1)

      pointer( ptr, get_data1D )
      ptr = LOC(get_data)
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe, block, tag, recv_request=request)

    end subroutine MPP_RECV_SCALAR_

    subroutine MPP_SEND_SCALAR_( put_data, to_pe, plen, tag, request)
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      MPP_TYPE_, intent(in) :: put_data
      integer, optional, intent(in) :: plen
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request
      integer                       :: put_len
      MPP_TYPE_ :: put_data1D(1)
      MPP_TYPE_ :: dummy(1)

      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )

    end subroutine MPP_SEND_SCALAR_

    subroutine MPP_RECV_2D_( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request

      MPP_TYPE_ :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine MPP_RECV_2D_

    subroutine MPP_SEND_2D_( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:)
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request
      MPP_TYPE_ :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine MPP_SEND_2D_

    subroutine MPP_RECV_3D_( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request

      MPP_TYPE_ :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request  )
    end subroutine MPP_RECV_3D_

    subroutine MPP_SEND_3D_( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:)
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request
      MPP_TYPE_ :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine MPP_SEND_3D_

    subroutine MPP_RECV_4D_( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request

      MPP_TYPE_ :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request  )
    end subroutine MPP_RECV_4D_

    subroutine MPP_SEND_4D_( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:,:)
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request
      MPP_TYPE_ :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine MPP_SEND_4D_

    subroutine MPP_RECV_5D_( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request

      MPP_TYPE_ :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request  )
    end subroutine MPP_RECV_5D_

    subroutine MPP_SEND_5D_( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:,:,:)
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: request
      MPP_TYPE_ :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine MPP_SEND_5D_
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_BROADCAST_SCALAR_( data, from_pe, pelist )
      MPP_TYPE_, intent(inout) :: data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(1)

      pointer( ptr, data1D )

      ptr = LOC(data)
      call MPP_BROADCAST_( data1D, 1, from_pe, pelist )

      return
    end subroutine MPP_BROADCAST_SCALAR_
    
    subroutine MPP_BROADCAST_2D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine MPP_BROADCAST_2D_
    
    subroutine MPP_BROADCAST_3D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
   end subroutine MPP_BROADCAST_3D_
    
    subroutine MPP_BROADCAST_4D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine MPP_BROADCAST_4D_
    
    subroutine MPP_BROADCAST_5D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine MPP_BROADCAST_5D_
