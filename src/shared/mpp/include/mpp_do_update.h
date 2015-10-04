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
    subroutine MPP_DO_UPDATE_3D_( f_addrs, domain, update, d_type, ke, b_addrs, b_size, flags)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND),         intent(in) :: f_addrs(:,:)
      type(domain2D),             intent(in) :: domain
      type(overlapSpec),          intent(in) :: update
      MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
      integer,                    intent(in) :: ke
      integer(LONG_KIND),         intent(in) :: b_addrs(:,:)
      integer,                    intent(in) :: b_size
      integer, optional,          intent(in) :: flags

      MPP_TYPE_ :: field(update%xbegin:update%xend, update%ybegin:update%yend,ke)
      MPP_TYPE_ :: fillbuffer(b_size)
      pointer(ptr_field, field)
      pointer(ptr_buffer, fillbuffer)
      integer                     :: update_flags
      type(overlap_type), pointer :: overPtr => NULL()      
      character(len=8)            :: text

!equate to mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer( ptr, buffer )
      integer :: buffer_pos
     
!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      integer,    allocatable :: msg1(:), msg2(:)
      logical :: send(8), recv(8), update_edge_only
      integer :: to_pe, from_pe, pos, msgsize
      integer :: n, l_size, l, m, i, j, k
      integer :: is, ie, js, je, tMe, dir
      integer :: start, start1, start2, index, is1, ie1, js1, je1, ni, nj, total
      integer :: buffer_recv_size, nlist, outunit

      outunit = stdout()
      ptr = LOC(mpp_domains_stack)
      l_size = size(f_addrs,1)

      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) )update_flags = flags

      update_edge_only = BTEST(update_flags, EDGEONLY)
      recv(1) = BTEST(update_flags,EAST)
      recv(3) = BTEST(update_flags,SOUTH)
      recv(5) = BTEST(update_flags,WEST)
      recv(7) = BTEST(update_flags,NORTH)
      if( update_edge_only ) then
         if( .NOT. (recv(1) .OR. recv(3) .OR. recv(5) .OR. recv(7)) ) then
            recv(1) = .true.
            recv(3) = .true.
            recv(5) = .true.
            recv(7) = .true.
         endif
      else
         recv(2) = recv(1) .AND. recv(3)
         recv(4) = recv(3) .AND. recv(5)
         recv(6) = recv(5) .AND. recv(7)
         recv(8) = recv(7) .AND. recv(1)
      endif
      send    = recv

      if(debug_message_passing) then
         nlist = size(domain%list(:))  
         allocate(msg1(0:nlist-1), msg2(0:nlist-1) )
         msg1 = 0
         msg2 = 0
         do m = 1, update%nrecv
            overPtr => update%recv(m)
            msgsize = 0
            do n = 1, overPtr%count
               dir = overPtr%dir(n)
               if(recv(dir)) then
                  is = overPtr%is(n); ie = overPtr%ie(n)
                  js = overPtr%js(n); je = overPtr%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
            end do
            from_pe = update%recv(m)%pe
            l = from_pe-mpp_root_pe()
            call mpp_recv( msg1(l), glen=1, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
            msg2(l) = msgsize
         enddo

         do m = 1, update%nsend
            overPtr => update%send(m)
            msgsize = 0
            do n = 1, overPtr%count
               dir = overPtr%dir(n)
               if(send(dir)) then
                  is = overPtr%is(n); ie = overPtr%ie(n)
                  js = overPtr%js(n); je = overPtr%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
            end do
            call mpp_send( msgsize, plen=1, to_pe=overPtr%pe, tag=COMM_TAG_1 )
         enddo
         call mpp_sync_self(check=EVENT_RECV)

         do m = 0, nlist-1
            if(msg1(m) .NE. msg2(m)) then
               print*, "My pe = ", mpp_pe(), ",domain name =", trim(domain%name), ",from pe=", &
                    domain%list(m)%pe, ":send size = ", msg1(m), ", recv size = ", msg2(m)
               call mpp_error(FATAL, "mpp_do_update: mismatch on send and recv size")
            endif
         enddo
         call mpp_sync_self()
         write(outunit,*)"NOTE from mpp_do_update: message sizes are matched between send and recv for domain " &
                          //trim(domain%name)
         deallocate(msg1, msg2)
      endif

      !recv
      buffer_pos = 0  
      do m = 1, update%nrecv
         overPtr => update%recv(m)
         if( overPtr%count == 0 )cycle
         call mpp_clock_begin(recv_clock)
         msgsize = 0
         do n = 1, overPtr%count
            dir = overPtr%dir(n)
            if(recv(dir)) then
               is = overPtr%is(n); ie = overPtr%ie(n)
               js = overPtr%js(n); je = overPtr%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end if
         end do

         msgsize = msgsize*ke*l_size
         if( msgsize.GT.0 )then
            from_pe = overPtr%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_2 )
            buffer_pos = buffer_pos + msgsize
         end if
         call mpp_clock_end(recv_clock)
      end do ! end do m = 1, update%nrecv
      buffer_recv_size = buffer_pos

      ! send
      do m = 1, update%nsend
         overPtr => update%send(m)
         if( overPtr%count == 0 )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         msgsize = 0
         do n = 1, overPtr%count
            dir = overPtr%dir(n)
            if( send(dir) )  msgsize = msgsize + overPtr%msgsize(n)
         enddo
         if( msgsize.GT.0 )then
            msgsize = msgsize*ke*l_size
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_START_UPDATE_DOMAINS: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
         end if

         do n = 1, overPtr%count
            dir = overPtr%dir(n)
            if( send(dir) ) then
               tMe = overPtr%tileMe(n)
               is = overPtr%is(n); ie = overPtr%ie(n)
               js = overPtr%js(n); je = overPtr%je(n)
               if( overptr%is_refined(n) ) then
                  do l=1,l_size  ! loop over number of fields
                     ptr_field = f_addrs(l, tMe)
                     do k = 1,ke  
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = field(i,j,k)
                           end do
                        end do
                     end do
                  end do
               else
                  select case( overPtr%rotation(n) )
                  case(ZERO)
                     do l=1,l_size  ! loop over number of fields
                        ptr_field = f_addrs(l, tMe)
                        do k = 1,ke  
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 buffer(pos) = field(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  case( MINUS_NINETY ) 
                     do l=1,l_size  ! loop over number of fields
                        ptr_field = f_addrs(l, tMe)
                        do k = 1,ke  
                           do i = is, ie
                              do j = je, js, -1
                                 pos = pos + 1
                                 buffer(pos) = field(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  case( NINETY ) 
                     do l=1,l_size  ! loop over number of fields
                        ptr_field = f_addrs(l, tMe)
                        do k = 1,ke  
                           do i = ie, is, -1
                              do j = js, je
                                 pos = pos + 1
                                 buffer(pos) = field(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  case( ONE_HUNDRED_EIGHTY ) 
                     do l=1,l_size  ! loop over number of fields
                        ptr_field = f_addrs(l, tMe)
                        do k = 1,ke  
                           do j = je, js, -1
                              do i = ie, is, -1
                                 pos = pos + 1
                                 buffer(pos) = field(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  end select
               end if
            endif
         end do ! do n = 1, overPtr%count

         call mpp_clock_end(pack_clock)
         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = overPtr%pe
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_2 )
            buffer_pos = pos
         end if
         call mpp_clock_end(send_clock)
      end do ! end do ist = 0,nlist-1

      !unpack recv
      !unpack halos in reverse order
!     ptr_rfield = f_addrs(1)
      call mpp_clock_begin(wait_clock)
      call mpp_sync_self(check=EVENT_RECV)
      call mpp_clock_end(wait_clock)
      buffer_pos = buffer_recv_size      

      do m = update%nrecv, 1, -1
         overPtr => update%recv(m)
         if( overPtr%count == 0 )cycle
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         do n = overPtr%count, 1, -1
            dir = overPtr%dir(n)
            if( recv(dir) ) then
               tMe = overPtr%tileMe(n)
               is = overPtr%is(n); ie = overPtr%ie(n)
               js = overPtr%js(n); je = overPtr%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               if(OverPtr%is_refined(n)) then
                  index = overPtr%index(n)
                  is1 = update%rSpec(tMe)%isNbr(index); ie1 = update%rSpec(tMe)%ieNbr(index)
                  js1 = update%rSpec(tMe)%jsNbr(index); je1 = update%rSpec(tMe)%jeNbr(index)
                  ni = ie1 - is1 + 1
                  nj = je1 - js1 + 1
                  total = ni*nj
                  start = (update%rSpec(tMe)%start(index)-1)*ke
                  if(start+total*ke>b_size ) call mpp_error(FATAL, &
                       "MPP_DO_UPDATE: b_size is less than the size of the data to be filled.")
                  msgsize = ie - is + 1
                  do l=1, l_size  ! loop over number of fields
                     ptr_buffer = b_addrs(l, tMe)
                     start1 = start + (js-js1)*ni + is - is1
                     do k = 1, ke
                        start2 = start1
                        do j = js, je
                           fillbuffer(start2+1:start2+msgsize) = buffer(pos+1:pos+msgsize)
                           start2 = start2 + ni
                           pos   = pos + msgsize
                        end do
                        start1 = start1 + total
                     end do
                  end do
               else
                  do l=1,l_size  ! loop over number of fields
                     ptr_field = f_addrs(l, tMe)
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              field(i,j,k) = buffer(pos)
                           end do
                        end do
                     end do
                  end do
               endif
            end if
         end do ! do n = 1, overPtr%count
         call mpp_clock_end(unpk_clock)
      end do

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)
      return
    end subroutine MPP_DO_UPDATE_3D_
