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
subroutine MPP_DO_UPDATE_NEST_FINE_3D_(f_addrs, nest_domain, update, d_type, ke, wb_addrs, eb_addrs, &
                                   sb_addrs, nb_addrs, flags)
  integer(LONG_KIND),         intent(in) :: f_addrs(:)
  type(nest_domain_type),     intent(in) :: nest_domain
  type(nestSpec),          intent(in) :: update
  MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
  integer,                    intent(in) :: ke
  integer(LONG_KIND),         intent(in) :: wb_addrs(:)
  integer(LONG_KIND),         intent(in) :: eb_addrs(:)
  integer(LONG_KIND),         intent(in) :: sb_addrs(:)
  integer(LONG_KIND),         intent(in) :: nb_addrs(:)
  integer,                    intent(in) :: flags

  character(len=8)            :: text
  type(overlap_type), pointer :: overPtr => NULL()
  logical   :: send(8), recv(8), update_edge_only
  integer   :: from_pe, to_pe, dir
  integer   :: m, n, l, i, j, k
  integer   :: is, ie, js, je, l_size
  integer   :: buffer_pos, msgsize
  integer   :: buffer_recv_size, pos
  MPP_TYPE_ :: field(update%xbegin:update%xend, update%ybegin:update%yend,ke)  
  MPP_TYPE_ :: wbuffer(update%west%is_you :update%west%ie_you,  update%west%js_you :update%west%je_you, ke)
  MPP_TYPE_ :: ebuffer(update%east%is_you :update%east%ie_you,  update%east%js_you :update%east%je_you, ke)
  MPP_TYPE_ :: sbuffer(update%south%is_you:update%south%ie_you, update%south%js_you:update%south%je_you,ke)
  MPP_TYPE_ :: nbuffer(update%north%is_you:update%north%ie_you, update%north%js_you:update%north%je_you,ke)
  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

  pointer(ptr_field, field)
  pointer(ptr_buffer, buffer )   
  pointer(ptr_wbuffer, wbuffer)
  pointer(ptr_ebuffer, ebuffer)
  pointer(ptr_sbuffer, sbuffer)
  pointer(ptr_nbuffer, nbuffer)

  update_edge_only = BTEST(flags, EDGEONLY)
  recv(1) = BTEST(flags,EAST)
  recv(3) = BTEST(flags,SOUTH)
  recv(5) = BTEST(flags,WEST)
  recv(7) = BTEST(flags,NORTH)
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

  ptr_buffer = LOC(mpp_domains_stack)
  l_size = size(f_addrs(:))

  !--- pre-post receiving
  buffer_pos = 0  
  do m = 1, update%nrecv
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_recv_clock)
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
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_FINE_3D_: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
     call mpp_clock_end(nest_recv_clock)
  end do ! end do m = 1, update%nrecv
  buffer_recv_size = buffer_pos

  !--- pack and send the data
  do m = 1, update%nsend
     overPtr => update%send(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_pack_clock)
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
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_FINE_3D_: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if

     do n = 1, overPtr%count
        dir = overPtr%dir(n)
        if( send(dir) ) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           do l=1,l_size  ! loop over number of fields
              ptr_field = f_addrs(l)
              do k = 1,ke  
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do
           end do
        endif
     end do ! do n = 1, overPtr%count

     call mpp_clock_end(nest_pack_clock)
     call mpp_clock_begin(nest_send_clock)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        to_pe = overPtr%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1 )
        buffer_pos = pos
     end if
     call mpp_clock_end(nest_send_clock)
  end do ! end do list = 0,nlist-1

  !unpack buffer
  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(nest_wait_clock)

  buffer_pos = buffer_recv_size      

  call mpp_clock_begin(nest_unpk_clock)
  do m = update%nrecv, 1, -1
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle

     pos = buffer_pos
     do n = overPtr%count, 1, -1
        dir = overPtr%dir(n)
        if( recv(dir) ) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           select case (dir)
           case ( 1 ) ! east
              do l=1,l_size  ! loop over number of fields
                 ptr_ebuffer = eb_addrs(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          ebuffer(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
              end do
           case ( 3 ) ! south
              do l=1,l_size  ! loop over number of fields
                 ptr_sbuffer = sb_addrs(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbuffer(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
              end do
           case ( 5 ) ! west
              do l=1,l_size  ! loop over number of fields
                 ptr_wbuffer = wb_addrs(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          wbuffer(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
              end do
           case ( 7 ) ! north
              do l=1,l_size  ! loop over number of fields
                 ptr_nbuffer = nb_addrs(l)
                 do k = 1,ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          nbuffer(i,j,k) = buffer(pos)
                       end do
                    end do
                 end do
              end do
           end select
        endif
     end do ! do n = 1, overPtr%count
  end do
  call mpp_clock_end(nest_unpk_clock)

  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self( )
  call mpp_clock_end(nest_wait_clock)
      return  

end subroutine MPP_DO_UPDATE_NEST_FINE_3D_



!###############################################################################
subroutine MPP_DO_UPDATE_NEST_COARSE_3D_(f_addrs, nest_domain, update, d_type, ke, b_addrs)
  integer(LONG_KIND),         intent(in) :: f_addrs(:)
  type(nest_domain_type),     intent(in) :: nest_domain
  type(nestSpec),             intent(in) :: update
  MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
  integer,                    intent(in) :: ke
  integer(LONG_KIND),         intent(in) :: b_addrs(:)

  character(len=8)            :: text
  type(overlap_type), pointer :: overPtr => NULL()
  integer   :: from_pe, to_pe
  integer   :: m, n, l, i, j, k
  integer   :: is, ie, js, je, l_size
  integer   :: buffer_pos, msgsize
  integer   :: buffer_recv_size, pos
  MPP_TYPE_ :: field(update%xbegin:update%xend, update%ybegin:update%yend,ke)  
  MPP_TYPE_ :: fillbuffer(update%center%is_you:update%center%ie_you,  update%center%js_you :update%center%je_you, ke)
  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

  pointer(ptr_field, field)
  pointer(ptr_buffer, buffer )   
  pointer(ptr_fillbuffer, fillbuffer)

  ptr_buffer = LOC(mpp_domains_stack)
  l_size = size(f_addrs(:))

  !--- pre-post receiving
  buffer_pos = 0  
  do m = 1, update%nrecv
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_recv_clock)
     msgsize = 0
     do n = 1, overPtr%count
        is = overPtr%is(n); ie = overPtr%ie(n)
        js = overPtr%js(n); je = overPtr%je(n)
        msgsize = msgsize + (ie-is+1)*(je-js+1)
     end do

     msgsize = msgsize*ke*l_size
     if( msgsize.GT.0 )then
        from_pe = overPtr%pe
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_COARSE_3D_: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_2 )
        buffer_pos = buffer_pos + msgsize
     end if
     call mpp_clock_end(nest_recv_clock)
  end do ! end do m = 1, update%nrecv
  buffer_recv_size = buffer_pos

  !--- pack and send the data
  do m = 1, update%nsend
     overPtr => update%send(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_pack_clock)
     pos = buffer_pos
     msgsize = 0
     do n = 1, overPtr%count
        msgsize = msgsize + overPtr%msgsize(n)
     enddo
     if( msgsize.GT.0 )then
        msgsize = msgsize*ke*l_size
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_COARSE_3D_: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if

     do n = 1, overPtr%count
        is = overPtr%is(n); ie = overPtr%ie(n)
        js = overPtr%js(n); je = overPtr%je(n)
        do l=1,l_size  ! loop over number of fields
           ptr_field = f_addrs(l)
           do k = 1,ke  
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    buffer(pos) = field(i,j,k)
                 end do
              end do
           end do
        end do
     end do ! do n = 1, overPtr%count

     call mpp_clock_end(nest_pack_clock)
     call mpp_clock_begin(nest_send_clock)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        to_pe = overPtr%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_2 )
        buffer_pos = pos
     end if
     call mpp_clock_end(nest_send_clock)
  end do ! end do list = 0,nlist-1

  !unpack buffer
  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(nest_wait_clock)

  buffer_pos = buffer_recv_size      

  call mpp_clock_begin(nest_unpk_clock)
  do m = update%nrecv, 1, -1
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle

     pos = buffer_pos
     do n = overPtr%count, 1, -1
        is = overPtr%is(n); ie = overPtr%ie(n)
        js = overPtr%js(n); je = overPtr%je(n)
        msgsize = (ie-is+1)*(je-js+1)*ke*l_size
        pos = buffer_pos - msgsize
        buffer_pos = pos
        do l=1,l_size  ! loop over number of fields
           ptr_fillbuffer = b_addrs(l)
           do k = 1,ke
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    fillbuffer(i,j,k) = buffer(pos)
                 end do
              end do
           end do
        end do
     end do ! do n = 1, overPtr%count
  end do
  call mpp_clock_end(nest_unpk_clock)

  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self( )
  call mpp_clock_end(nest_wait_clock)
      return  

end subroutine MPP_DO_UPDATE_NEST_COARSE_3D_
