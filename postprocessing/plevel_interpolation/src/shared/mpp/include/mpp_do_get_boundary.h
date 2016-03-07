! -*-f90-*- 
subroutine MPP_DO_GET_BOUNDARY_3D_( f_addrs, domain, bound, b_addrs, bsize, ke, d_type, flags)
  type(domain2D), intent(in)      :: domain
  type(overlapSpec),  intent(in)  :: bound
  integer(LONG_KIND), intent(in)  :: f_addrs(:,:)
  integer(LONG_KIND), intent(in)  :: b_addrs(:,:,:)
  integer,            intent(in)  :: bsize(:), ke
  MPP_TYPE_, intent(in)           :: d_type  ! creates unique interface
  integer, intent(in)             :: flags

  MPP_TYPE_ :: field(bound%xbegin:bound%xend, bound%ybegin:bound%yend,ke)
  MPP_TYPE_ :: ebuffer(bsize(1), ke), sbuffer(bsize(2), ke), wbuffer(bsize(3), ke), nbuffer(bsize(4), ke)
  pointer(ptr_field, field)
  pointer(ptr_ebuffer, ebuffer)
  pointer(ptr_sbuffer, sbuffer)  
  pointer(ptr_wbuffer, wbuffer)
  pointer(ptr_nbuffer, nbuffer)

  integer,    allocatable :: msg1(:), msg2(:)
  logical                 :: recv(4), send(4)
  integer                 :: nlist, buffer_pos, pos, tMe, from_pe
  integer                 :: i, j, k, l, m, n, index, buffer_recv_size
  integer                 :: is, ie, js, je, msgsize, l_size
  character(len=8)        :: text

  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))

  pointer( ptr, buffer )
  ptr = LOC(mpp_domains_stack)

  l_size = size(f_addrs,1)
  recv(1) = BTEST(flags,EAST)
  recv(2) = BTEST(flags,SOUTH)
  recv(3) = BTEST(flags,WEST)
  recv(4) = BTEST(flags,NORTH)

  send = recv

  nlist = size(domain%list(:))  

  if(debug_message_passing) then
      allocate(msg1(0:nlist-1), msg2(0:nlist-1) )
      msg1 = 0
      msg2 = 0

      do m = 1, bound%nrecv
         msgsize = 0
         do n = 1, bound%recv(m)%count
            if(recv(bound%recv(m)%dir(n))) then
               is = bound%recv(m)%is(n); ie = bound%recv(m)%ie(n)
               js = bound%recv(m)%js(n); je = bound%recv(m)%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end if
         end do
         from_pe = bound%recv(m)%pe
         l = from_pe-mpp_root_pe()
         call mpp_recv( msg1(l), glen=1, from_pe=from_pe, block=.FALSE.)
         msg2(l) = msgsize
      enddo

      do m = 1, bound%nsend
         msgsize = 0
         do n = 1, bound%send(m)%count
            if(recv(bound%send(m)%dir(n))) then
               is = bound%send(m)%is(n); ie = bound%send(m)%ie(n)
               js = bound%send(m)%js(n); je = bound%send(m)%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end if
         end do
         call mpp_send( msgsize, plen=1, to_pe=bound%send(m)%pe)
      enddo

      call mpp_sync_self(check=EVENT_RECV)

      do m = 0, nlist-1
         if(msg1(m) .NE. msg2(m)) then
            print*, "My pe = ", mpp_pe(), ",domain name =", trim(domain%name), ",from pe=", &
               domain%list(m)%pe, ":send size = ", msg1(m), ", recv size = ", msg2(m)
            call mpp_error(FATAL, "mpp_do_get_boundary: mismatch on send and recv size")
         endif
      enddo
      call mpp_sync_self()
      write(stdout(),*)"NOTE from mpp_do_get_boundary: message sizes are matched between send and recv for domain " &
                       //trim(domain%name)
      deallocate(msg1, msg2)
  endif
  !recv
  buffer_pos = 0     
  do m = 1, bound%nrecv
     msgsize = 0
     do n = 1, bound%recv(m)%count
        if(recv(bound%recv(m)%dir(n))) then
           is = bound%recv(m)%is(n); ie = bound%recv(m)%ie(n)
           js = bound%recv(m)%js(n); je = bound%recv(m)%je(n)
           msgsize = msgsize + (ie-is+1)*(je-js+1)
        end if
     end do
     msgsize = msgsize*ke*l_size
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_OLD: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=bound%recv(m)%pe, block=.false. )
        buffer_pos = buffer_pos + msgsize
     end if
  end do
  buffer_recv_size = buffer_pos

  ! send
  do m = 1, bound%nsend
     pos = buffer_pos
     do n = 1, bound%send(m)%count
        if(send(bound%send(m)%dir(n))) then
           is = bound%send(m)%is(n); ie = bound%send(m)%ie(n)
           js = bound%send(m)%js(n); je = bound%send(m)%je(n)
           tMe = bound%send(m)%tileMe(n)
           select case( bound%send(m)%rotation(n) )
           case(ZERO)
              do l=1,l_size
                 ptr_field = f_addrs(l, tMe)
                 do k = 1, ke
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( MINUS_NINETY )
              do l=1,l_size
                 ptr_field = f_addrs(l, tMe)
                 do k = 1, ke
                    do j = je, js, -1
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case( NINETY )
              do l=1,l_size
                 ptr_field = f_addrs(l, tMe)
                 do k = 1, ke
                    do j = js, je
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           case (ONE_HUNDRED_EIGHTY) 
              do l=1,l_size
                 ptr_field = f_addrs(l, tMe)
                 do k = 1, ke
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              end do
           end select
        end if ! if(send(bound%dir(n)))
     end do ! do n = 1, bound%count
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then  
        !--- maybe we do not need the following stack size check.
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_OLD: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=bound%send(m)%pe )
        buffer_pos = pos
     end if
  end do

  call mpp_clock_begin(wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(wait_clock)
  buffer_pos = buffer_recv_size  

  !unpack recv
  !unpack buffer in reverse order.
  do m = bound%nrecv, 1, -1
     do n = bound%recv(m)%count, 1, -1
        if(recv(bound%recv(m)%dir(n))) then
           is = bound%recv(m)%is(n); ie = bound%recv(m)%ie(n)
           js = bound%recv(m)%js(n); je = bound%recv(m)%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           tMe = bound%recv(m)%tileMe(n)
           select case( bound%recv(m)%dir(n) )
           case ( 1 ) ! EAST
              do l=1,l_size
                 ptr_ebuffer = b_addrs(1, l, tMe)              
                 do k = 1, ke
                    index = bound%recv(m)%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          ebuffer(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 2 ) ! SOUTH
              do l=1,l_size
                 ptr_sbuffer = b_addrs(2, l, tMe)   
                 do k = 1, ke
                    index = bound%recv(m)%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          sbuffer(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 3 ) ! WEST
              do l=1,l_size
                 ptr_wbuffer = b_addrs(3, l, tMe)   
                 do k = 1, ke
                    index = bound%recv(m)%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          wbuffer(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           case ( 4 ) ! norTH
              do l=1,l_size
                 ptr_nbuffer = b_addrs(4, l, tMe)   
                 do k = 1, ke
                    index = bound%recv(m)%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          nbuffer(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           end select
        end if
     end do
  end do

  call mpp_sync_self( )


end subroutine MPP_DO_GET_BOUNDARY_3D_


subroutine MPP_DO_GET_BOUNDARY_3D_V_(f_addrsx, f_addrsy, domain, boundx, boundy, b_addrsx, b_addrsy, &
                                        bsizex, bsizey, ke, d_type, flags)
  type(domain2D),     intent(in)  :: domain
  type(overlapSpec),  intent(in)  :: boundx, boundy
  integer(LONG_KIND), intent(in)  :: f_addrsx(:,:), f_addrsy(:,:)
  integer(LONG_KIND), intent(in)  :: b_addrsx(:,:,:), b_addrsy(:,:,:)
  integer,            intent(in)  :: bsizex(:), bsizey(:), ke
  MPP_TYPE_, intent(in)           :: d_type  ! creates unique interface
  integer, intent(in)             :: flags

  MPP_TYPE_ :: fieldx(boundx%xbegin:boundx%xend, boundx%ybegin:boundx%yend,ke)
  MPP_TYPE_ :: fieldy(boundy%xbegin:boundy%xend, boundy%ybegin:boundy%yend,ke)
  MPP_TYPE_ :: ebufferx(bsizex(1), ke), sbufferx(bsizex(2), ke), wbufferx(bsizex(3), ke), nbufferx(bsizex(4), ke)
  MPP_TYPE_ :: ebuffery(bsizey(1), ke), sbuffery(bsizey(2), ke), wbuffery(bsizey(3), ke), nbuffery(bsizey(4), ke)
  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)
  pointer(ptr_ebufferx, ebufferx)
  pointer(ptr_sbufferx, sbufferx)  
  pointer(ptr_wbufferx, wbufferx)
  pointer(ptr_nbufferx, nbufferx)
  pointer(ptr_ebuffery, ebuffery)
  pointer(ptr_sbuffery, sbuffery)  
  pointer(ptr_wbuffery, wbuffery)
  pointer(ptr_nbuffery, nbuffery)

  integer,    allocatable :: msg1(:), msg2(:)
  logical                 :: recv(4), send(4)
  integer                 :: nlist, buffer_pos, pos, tMe, m
  integer                 :: is, ie, js, je, msgsize, l_size, buffer_recv_size
  integer                 :: i, j, k, l, n, index, to_pe, from_pe
  integer                 :: rank_x, rank_y, cur_rank, ind_x, ind_y
  integer                 :: nsend_x, nsend_y, nrecv_x, nrecv_y
  character(len=8)        :: text

  MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
  pointer( ptr, buffer )
  ptr = LOC(mpp_domains_stack)

  l_size = size(f_addrsx,1)
  recv(1) = BTEST(flags,EAST)
  recv(2) = BTEST(flags,SOUTH)
  recv(3) = BTEST(flags,WEST)
  recv(4) = BTEST(flags,NORTH)
  send = recv

  nlist = size(domain%list(:))  

  nsend_x = boundx%nsend
  nsend_y = boundy%nsend
  nrecv_x = boundx%nrecv
  nrecv_y = boundy%nrecv

  if(debug_message_passing) then
     allocate(msg1(0:nlist-1), msg2(0:nlist-1) )
     msg1 = 0
     msg2 = 0

     cur_rank = get_rank_recv(domain, boundx, boundy, rank_x, rank_y, ind_x, ind_y) 

     do while ( ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y )
        msgsize = 0
        if(cur_rank == rank_x) then
           from_pe = boundx%recv(ind_x)%pe
           do n = 1, boundx%recv(ind_x)%count
              if(recv(boundx%recv(ind_x)%dir(n))) then
                 is = boundx%recv(ind_x)%is(n); ie = boundx%recv(ind_x)%ie(n)
                 js = boundx%recv(ind_x)%js(n); je = boundx%recv(ind_x)%je(n)
                 msgsize = msgsize + (ie-is+1)*(je-js+1)
              end if
           end do
           ind_x = ind_x+1
           if(ind_x .LE. nrecv_x) then
              rank_x = boundx%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = -1
           endif
        endif

        if(cur_rank == rank_y) then
           from_pe = boundy%recv(ind_y)%pe
           do n = 1, boundy%recv(ind_y)%count
              if(recv(boundy%recv(ind_y)%dir(n))) then
                 is = boundy%recv(ind_y)%is(n); ie = boundy%recv(ind_y)%ie(n)
                 js = boundy%recv(ind_y)%js(n); je = boundy%recv(ind_y)%je(n)
                 msgsize = msgsize + (ie-is+1)*(je-js+1)
              end if
           end do
           ind_y = ind_y+1
           if(ind_y .LE. nrecv_y) then
              rank_y = boundy%recv(ind_y)%pe - domain%pe 
              if(rank_y .LE.0) rank_y = rank_y + nlist
           else
              rank_y = -1
           endif
        endif
        cur_rank = max(rank_x, rank_y)
        m = from_pe-mpp_root_pe()
        call mpp_recv( msg1(m), glen=1, from_pe=from_pe, block=.FALSE.)
        msg2(m) = msgsize
     end do

     cur_rank = get_rank_send(domain, boundx, boundy, rank_x, rank_y, ind_x, ind_y) 
     do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
        msgsize = 0
        if(cur_rank == rank_x) then
           to_pe = boundx%send(ind_x)%pe
           do n = 1, boundx%send(ind_x)%count
              if(send(boundx%send(ind_x)%dir(n))) then
                 is = boundx%send(ind_x)%is(n); ie = boundx%send(ind_x)%ie(n)
                 js = boundx%send(ind_x)%js(n); je = boundx%send(ind_x)%je(n)
                 msgsize = msgsize + (ie-is+1)*(je-js+1)
              endif
           enddo
           ind_x = ind_x+1
           if(ind_x .LE. nsend_x) then
              rank_x = boundx%send(ind_x)%pe - domain%pe 
              if(rank_x .LT.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1
           endif
        endif

        if(cur_rank == rank_y) then
           to_pe = boundy%send(ind_y)%pe
           do n = 1, boundy%send(ind_y)%count
              if(send(boundy%send(ind_y)%dir(n))) then
                 is = boundy%send(ind_y)%is(n); ie = boundy%send(ind_y)%ie(n)
                 js = boundy%send(ind_y)%js(n); je = boundy%send(ind_y)%je(n)
                 msgsize = msgsize + (ie-is+1)*(je-js+1)
              end if
           end do
           ind_y = ind_y+1
           if(ind_y .LE. nsend_y) then
              rank_y = boundy%send(ind_y)%pe - domain%pe 
              if(rank_y .LT.0) rank_y = rank_y + nlist
           else
              rank_y = nlist+1
           endif
        endif
        cur_rank = min(rank_x, rank_y)
        call mpp_send( msgsize, plen=1, to_pe=to_pe)        
     enddo

      call mpp_sync_self(check=EVENT_RECV)
      do m = 0, nlist-1
         if(msg1(m) .NE. msg2(m)) then
            print*, "My pe = ", mpp_pe(), ",domain name =", trim(domain%name), ",from pe=", &
               domain%list(m)%pe, ":send size = ", msg1(m), ", recv size = ", msg2(m)
            call mpp_error(FATAL, "mpp_do_get_boundaryV: mismatch on send and recv size")
         endif
      enddo

      call mpp_sync_self()
      write(stdout(),*)"NOTE from mpp_do_get_boundary_V: message sizes are matched between send and recv for domain " &
                       //trim(domain%name)
      deallocate(msg1, msg2)
  endif

  !recv
  buffer_pos = 0     
  cur_rank = get_rank_recv(domain, boundx, boundy, rank_x, rank_y, ind_x, ind_y)   

  do while ( ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y )
     msgsize = 0
     if(cur_rank == rank_x) then
        from_pe = boundx%recv(ind_x)%pe
        do n = 1, boundx%recv(ind_x)%count
           if(recv(boundx%recv(ind_x)%dir(n))) then
              is = boundx%recv(ind_x)%is(n); ie = boundx%recv(ind_x)%ie(n)
              js = boundx%recv(ind_x)%js(n); je = boundx%recv(ind_x)%je(n)
              msgsize = msgsize + (ie-is+1)*(je-js+1)
           end if
        end do
        ind_x = ind_x+1
        if(ind_x .LE. nrecv_x) then
           rank_x = boundx%recv(ind_x)%pe - domain%pe 
           if(rank_x .LE.0) rank_x = rank_x + nlist
        else
           rank_x = -1
        endif
     endif

     if(cur_rank == rank_y) then
        from_pe = boundy%recv(ind_y)%pe
        do n = 1, boundy%recv(ind_y)%count
           if(recv(boundy%recv(ind_y)%dir(n))) then
              is = boundy%recv(ind_y)%is(n); ie = boundy%recv(ind_y)%ie(n)
              js = boundy%recv(ind_y)%js(n); je = boundy%recv(ind_y)%je(n)
              msgsize = msgsize + (ie-is+1)*(je-js+1)
           end if
        end do
        ind_y = ind_y+1
        if(ind_y .LE. nrecv_y) then
           rank_y = boundy%recv(ind_y)%pe - domain%pe 
           if(rank_y .LE.0) rank_y = rank_y + nlist
        else
           rank_y = -1
        endif
     endif
     cur_rank = max(rank_x, rank_y)
     msgsize = msgsize*ke*l_size
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_V_: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE. )
        buffer_pos = buffer_pos + msgsize
     end if
  end do
  buffer_recv_size = buffer_pos

  ! send
  cur_rank = get_rank_send(domain, boundx, boundy, rank_x, rank_y, ind_x, ind_y) 

  do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
     pos = buffer_pos
     if(cur_rank == rank_x) then
        to_pe = boundx%send(ind_x)%pe
        do n = 1, boundx%send(ind_x)%count
           if(send(boundx%send(ind_x)%dir(n))) then
              is = boundx%send(ind_x)%is(n); ie = boundx%send(ind_x)%ie(n)
              js = boundx%send(ind_x)%js(n); je = boundx%send(ind_x)%je(n)
              tMe = boundx%send(ind_x)%tileMe(n)
              select case( boundx%send(ind_x)%rotation(n) )
              case(ZERO)
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, tMe)
                    do k = 1, ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             buffer(pos) = fieldx(i,j,k)
                          end do
                       end do
                    end do
                 end do
              case( MINUS_NINETY )
                 if( BTEST(flags,SCALAR_BIT) ) then
                    do l=1,l_size
                       ptr_fieldy = f_addrsy(l, tMe)
                       do k = 1, ke
                          do j = je, js, -1
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 else
                    do l=1,l_size
                       ptr_fieldy = f_addrsy(l, tMe)
                       do k = 1, ke
                          do j = je, js, -1
                             do i = is, ie
                                pos = pos + 1
                                buffer(pos) = -fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 end if
              case( NINETY )
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, tMe)
                    do k = 1, ke
                       do j = js, je
                          do i = ie, is, -1
                             pos = pos + 1
                             buffer(pos) = fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end do
              case (ONE_HUNDRED_EIGHTY) 
                 if( BTEST(flags,SCALAR_BIT) ) then
                    do l=1,l_size
                       ptr_fieldx = f_addrsx(l, tMe)
                       do k = 1, ke
                          do j = je, js, -1
                             do i = ie, is, -1
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 else     
                    do l=1,l_size
                       ptr_fieldx = f_addrsx(l, tMe)
                       do k = 1, ke
                          do j = je, js, -1
                             do i = ie, is, -1
                                pos = pos + 1
                                buffer(pos) = -fieldx(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 end if
              end select
           end if ! if(send(boundx%dir(n)))
        end do  !do n = 1, boundx%count
        ind_x = ind_x+1
        if(ind_x .LE. nsend_x) then
           rank_x = boundx%send(ind_x)%pe - domain%pe 
           if(rank_x .LT.0) rank_x = rank_x + nlist
        else
           rank_x = nlist+1
        endif
     endif

     if(cur_rank == rank_y) then
        to_pe = boundy%send(ind_y)%pe
        do n = 1, boundy%send(ind_y)%count
           if(send(boundy%send(ind_y)%dir(n))) then
              is = boundy%send(ind_y)%is(n); ie = boundy%send(ind_y)%ie(n)
              js = boundy%send(ind_y)%js(n); je = boundy%send(ind_y)%je(n)
              tMe = boundy%send(ind_y)%tileMe(n)
              select case( boundy%send(ind_y)%rotation(n) )
              case(ZERO)
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, tMe)
                    do k = 1, ke
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             buffer(pos) = fieldy(i,j,k)
                          end do
                       end do
                    end do
                 end do
              case( MINUS_NINETY )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, tMe)
                    do k = 1, ke
                       do j = je, js, -1
                          do i = is, ie
                             pos = pos + 1
                             buffer(pos) = fieldx(i,j,k)
                          end do
                       end do
                    end do
                 end do
              case( NINETY )
                 if( BTEST(flags,SCALAR_BIT) ) then
                    do l=1,l_size
                       ptr_fieldx = f_addrsx(l, tMe)
                       do k = 1, ke
                          do j = js, je
                             do i = ie, is, -1
                                pos = pos + 1
                                buffer(pos) = fieldx(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 else
                    do l=1,l_size
                       ptr_fieldx = f_addrsx(l, tMe)
                       do k = 1, ke
                          do j = js, je
                             do i = ie, is, -1
                                pos = pos + 1
                                buffer(pos) = -fieldx(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 end if
              case (ONE_HUNDRED_EIGHTY) 
                 if( BTEST(flags,SCALAR_BIT) ) then
                    do l=1,l_size
                       ptr_fieldy = f_addrsy(l, tMe)
                       do k = 1, ke
                          do j = je, js, -1
                             do i = ie, is, -1
                                pos = pos + 1
                                buffer(pos) = fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 else     
                    do l=1,l_size
                       ptr_fieldy = f_addrsy(l, tMe)
                       do k = 1, ke
                          do j = je, js, -1
                             do i = ie, is, -1
                                pos = pos + 1
                                buffer(pos) = -fieldy(i,j,k)
                             end do
                          end do
                       end do
                    end do
                 end if
              end select
           end if ! if(send(boundy%dir(n)))
        end do    ! do n = 1, boundy%count
        ind_y = ind_y+1
        if(ind_y .LE. nsend_y) then
           rank_y = boundy%send(ind_y)%pe - domain%pe 
           if(rank_y .LT.0) rank_y = rank_y + nlist
        else
           rank_y = nlist+1
        endif
     endif
     cur_rank = min(rank_x, rank_y)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then  
        !--- maybe we do not need the following stack size check.
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_V_: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
        buffer_pos = pos
     end if

  end do       

  call mpp_sync_self(check=EVENT_RECV)

  !unpack recv
  !unpack buffer in reverse order.
  buffer_pos = buffer_recv_size  
  cur_rank = get_rank_unpack(domain, boundx, boundy, rank_x, rank_y, ind_x, ind_y) 

  do while(ind_x >0 .OR. ind_y >0)
     if(cur_rank == rank_y) then
        do n = boundy%recv(ind_y)%count, 1, -1
           if(recv(boundy%recv(ind_y)%dir(n))) then
              is = boundy%recv(ind_y)%is(n); ie = boundy%recv(ind_y)%ie(n)
              js = boundy%recv(ind_y)%js(n); je = boundy%recv(ind_y)%je(n)
              msgsize = (ie-is+1)*(je-js+1)*ke*l_size
              pos = buffer_pos - msgsize
              buffer_pos = pos
              tMe = boundy%recv(ind_y)%tileMe(n)
              select case( boundy%recv(ind_y)%dir(n) )
              case ( 1 ) ! EAST
                 do l=1,l_size
                    ptr_ebuffery = b_addrsy(1, l, tMe)     
                    do k = 1, ke
                       index = boundy%recv(ind_y)%index(n)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             ebuffery(index,k) = buffer(pos)
                             index = index + 1
                          end do
                       end do
                    end do
                 end do
              case ( 2 ) ! SOUTH
                 do l=1,l_size
                    ptr_sbuffery = b_addrsy(2, l, tMe)     
                    do k = 1, ke
                       index = boundy%recv(ind_y)%index(n)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             sbuffery(index,k) = buffer(pos)
                             index = index + 1
                          end do
                       end do
                    end do
                 end do
              case ( 3 ) ! WEST
                 do l=1,l_size
                    ptr_wbuffery = b_addrsy(3, l, tMe)     
                    do k = 1, ke
                       index = boundy%recv(ind_y)%index(n)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             wbuffery(index,k) = buffer(pos)
                             index = index + 1
                          end do
                       end do
                    end do
                 end do
              case ( 4 ) ! norTH
                 do l=1,l_size
                    ptr_nbuffery = b_addrsy(4, l, tMe)     
                    do k = 1, ke
                       index = boundy%recv(ind_y)%index(n)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             nbuffery(index,k) = buffer(pos)
                             index = index + 1
                          end do
                       end do
                    end do
                 end do
              end select
           end if
        end do
        ind_y = ind_y-1
        if(ind_y .GT. 0) then
           rank_y = boundy%recv(ind_y)%pe - domain%pe 
           if(rank_y .LE.0) rank_y = rank_y + nlist
        else
           rank_y = nlist+1
        endif
     endif

     if(cur_rank == rank_x) then
        do n = boundx%recv(ind_x)%count, 1, -1
           if(recv(boundx%recv(ind_x)%dir(n))) then
              is = boundx%recv(ind_x)%is(n); ie = boundx%recv(ind_x)%ie(n)
              js = boundx%recv(ind_x)%js(n); je = boundx%recv(ind_x)%je(n)
              msgsize = (ie-is+1)*(je-js+1)*ke*l_size
              pos = buffer_pos - msgsize
              buffer_pos = pos
              tMe = boundx%recv(ind_x)%tileMe(n)
              select case( boundx%recv(ind_x)%dir(n) )
              case ( 1 ) ! EAST
                 do l=1,l_size
                    ptr_ebufferx = b_addrsx(1, l, tMe)     
                    do k = 1, ke
                       index = boundx%recv(ind_x)%index(n)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             ebufferx(index,k) = buffer(pos)
                             index = index + 1
                          end do
                       end do
                    end do
                 end do
              case ( 2 ) ! SOUTH
                 do l=1,l_size
                    ptr_sbufferx = b_addrsx(2, l, tMe)     
                    do k = 1, ke
                       index = boundx%recv(ind_x)%index(n)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             sbufferx(index,k) = buffer(pos)
                             index = index + 1
                          end do
                       end do
                    end do
                 end do
              case ( 3 ) ! WEST
                 do l=1,l_size
                    ptr_wbufferx = b_addrsx(3, l, tMe)     
                    do k = 1, ke
                       index = boundx%recv(ind_x)%index(n)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             wbufferx(index,k) = buffer(pos)
                             index = index + 1
                          end do
                       end do
                    end do
                 end do
              case ( 4 ) ! norTH
                 do l=1,l_size
                    ptr_nbufferx = b_addrsx(4, l, tMe)     
                    do k = 1, ke
                       index = boundx%recv(ind_x)%index(n)
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             nbufferx(index,k) = buffer(pos)
                             index = index + 1
                          end do
                       end do
                    end do
                 end do
              end select
           end if
        end do
        ind_x = ind_x-1
        if(ind_x .GT. 0) then
           rank_x = boundx%recv(ind_x)%pe - domain%pe 
           if(rank_x .LE.0) rank_x = rank_x + nlist
        else
           rank_x = nlist+1
        endif
     endif
     cur_rank = min(rank_x, rank_y)
  end do

  call mpp_sync_self( )



end subroutine MPP_DO_GET_BOUNDARY_3D_V_
