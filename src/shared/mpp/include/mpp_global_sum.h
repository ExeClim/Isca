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

  function MPP_GLOBAL_SUM_( domain, field, flags, position, tile_count )
    MPP_TYPE_ :: MPP_GLOBAL_SUM_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,: MPP_EXTRA_INDICES_ )
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: position
    integer, intent(in), optional :: tile_count

    MPP_TYPE_, dimension(:,:),       allocatable :: field2D
    MPP_TYPE_, dimension(:,:),       allocatable :: global2D
    MPP_TYPE_, dimension(MAX_TILES), save        :: gsum, nbrgsum, mygsum
    
    integer :: i,j, ioff,joff, isc, iec, jsc, jec, is, ie, js, je, ishift, jshift, ioffset, joffset
    integer :: gxsize, gysize
    integer :: global_flag, tile, ntile, nlist, n, list, m

    if( domain%max_ntile_pe > MAX_TILES ) call mpp_error(FATAL, "MPP_GLOBAL_SUM: number of tiles is exceed MAX_TILES")
    ntile     = size(domain%x(:))
    nlist     = size(domain%list(:))
    tile = 1
    if(present(tile_count)) tile = tile_count
    global_flag = NON_BITWISE_EXACT_SUM
    if(present(flags)) global_flag = flags

    call mpp_get_domain_shift(domain, ishift, jshift, position)

    if( size(field,1).EQ.domain%x(tile)%compute%size+ishift .AND. size(field,2).EQ.domain%y(tile)%compute%size+jshift )then
!field is on compute domain
        ioff = -domain%x(tile)%compute%begin + 1 
        joff = -domain%y(tile)%compute%begin + 1
    else if( size(field,1).EQ.domain%x(tile)%memory%size+ishift .AND. size(field,2).EQ.domain%y(tile)%memory%size+jshift )then
!field is on data domain
        ioff = -domain%x(tile)%data%begin + 1
        joff = -domain%y(tile)%data%begin + 1
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_SUM_: incoming field array must match either compute domain or data domain.' )
    end if

    if(domain%ntiles > MAX_TILES)  call mpp_error( FATAL,  &
         'MPP_GLOBAL_SUM_: number of tiles on this mosaic is greater than MAXTILES')

    call mpp_get_compute_domain( domain, is,  ie,  js,  je,  tile_count = tile_count )
    isc = is; iec = ie + ishift; jsc = js; jec = je + jshift

    call mpp_get_global_domain(domain, xsize = gxsize, ysize = gysize )
    MPP_GLOBAL_SUM_ = 0
    if( global_flag == BITWISE_EXACT_SUM )then
        !this is bitwise exact across different PE counts.

       allocate( field2D (isc:iec,jsc:jec) )
       do j = jsc, jec
          do i = isc, iec
             field2D(i,j) = sum( field(i+ioff:i+ioff,j+joff:j+joff MPP_EXTRA_INDICES_) )
          end do
       end do
       allocate( global2D( gxsize+ishift, gysize+jshift ) )
       global2D = 0.

       !call mpp_global_field( domain, field2D, global2D, position=position, tile_count=tile_count )
       
       if ( present( tile_count ) ) then
           call mpp_global_field( domain, field2D, global2D, position=position, tile_count=tile_count )
       else    
           call mpp_global_field( domain, field2D, global2D, position=position )
       endif
       
       ioffset = domain%x(tile)%goffset*ishift; joffset = domain%y(tile)%goffset*jshift
       mygsum(tile) = sum(global2D(1:gxsize+ioffset,1:gysize+joffset))
       deallocate(global2D, field2d)
       if( tile == ntile) then 
          if(domain%ntiles == 1 ) then
             MPP_GLOBAL_SUM_ = mygsum(tile)
          else if( nlist == 1) then
             MPP_GLOBAL_SUM_ = sum(mygsum(1:ntile))
          else ! need to sum by the order of tile_count
             ! first fill the global sum on current pe.
             do n = 1, ntile
                gsum(domain%tile_id(n)) = mygsum(n)
             end do
             !--- send the data to other pe if the current pe is the root pe of any tile
             if( mpp_domain_is_tile_root_pe(domain) ) then
                do list = 1, nlist - 1
                   m = mod( domain%pos+list, nlist )
                   call mpp_send( mygsum(1), plen=ntile, to_pe=domain%list(m)%pe, tag=COMM_TAG_1 )
                end do
             end if
             call mpp_sync_self()
             !--- receive data from root_pe of each tile
             do list = 1, nlist - 1
                m = mod( domain%pos+nlist-list, nlist )
                if( domain%list(m)%pe == domain%list(m)%tile_root_pe ) then
                    call mpp_recv( nbrgsum(1), glen=size(domain%list(m)%x(:)), from_pe=domain%list(m)%pe, tag=COMM_TAG_1)
                    do n = 1, size(domain%list(m)%x(:))
                       gsum(domain%list(m)%tile_id(n)) = nbrgsum(n)
                    end do
                end if
             end do

             MPP_GLOBAL_SUM_ = sum(gsum(1:domain%ntiles))
          end if
       end if
    else  !this is not bitwise-exact across different PE counts
       ioffset = domain%x(tile)%loffset*ishift; joffset = domain%y(tile)%loffset*jshift
       mygsum(tile) = sum( field(is+ioff:ie+ioff+ioffset, js+joff:je+joff+joffset MPP_EXTRA_INDICES_) )
       if(tile == ntile) then
          MPP_GLOBAL_SUM_ = sum(mygsum)
          call mpp_sum( MPP_GLOBAL_SUM_, domain%list(:)%pe )
       end if
    end if

    return
  end function MPP_GLOBAL_SUM_
