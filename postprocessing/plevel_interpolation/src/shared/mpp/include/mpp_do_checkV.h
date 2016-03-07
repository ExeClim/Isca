! -*-f90-*- 
    subroutine MPP_DO_CHECK_3D_V_(f_addrsx,f_addrsy, domain, check_x, check_y, &
                                   d_type, ke, flags, name)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND),  intent(in)        :: f_addrsx(:,:), f_addrsy(:,:)
      type(domain2d),      intent(in)        :: domain
      type(overlapSpec),   intent(in)        :: check_x, check_y
      integer,             intent(in)        :: ke
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface
      integer, intent(in),          optional :: flags
      character(len=*), intent(in), optional :: name

      MPP_TYPE_ :: fieldx(check_x%xbegin:check_x%xend, check_x%ybegin:check_x%yend,ke)
      MPP_TYPE_ :: fieldy(check_y%xbegin:check_y%xend, check_y%ybegin:check_y%yend,ke)
      pointer(ptr_fieldx, fieldx)
      pointer(ptr_fieldy, fieldy)
      integer,    allocatable :: msg1(:), msg2(:)
      integer :: update_flags
      integer :: l_size, l, i, j, k, is, ie, js, je, n, m
      integer :: pos, nlist, msgsize
      integer :: to_pe, from_pe
      integer :: tMe
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer(ptr,buffer )
      integer :: buffer_pos
      character(len=8) :: text
      character(len=64) :: field_name      
      integer :: buffer_recv_size
      integer :: rank_x, rank_y, ind_x, ind_y, cur_rank
      integer :: nsend_x, nsend_y, nrecv_x, nrecv_y

      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) ) update_flags = flags

      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
      l_size = size(f_addrsx,1)
      nlist = size(domain%list(:))
      ptr = LOC(mpp_domains_stack)

      !--- if debug_update_level is not NO_DEBUG, check the consistency on the bounds 
      !--- (domain is symmetry or folded north edge). North bound will be checked when north edge is folded.
      !--- when domain is symmetry, For data on T-cell, no check is needed; for data on E-cell, 
      !--- data on East and West boundary will be checked ; For data on N-cell, data on North and South 
      !--- boundary will be checked; For data on C-cell, data on West, East, South, North will be checked.
      !--- The check will be done in the following way: Western boundary data sent to Eastern boundary to check
      !--- and Southern boundary to check

      if(present(name)) then
         field_name = name
      else
         field_name = "un-named"
      end if

      nsend_x = check_x%nsend
      nsend_y = check_y%nsend
      nrecv_x = check_x%nrecv
      nrecv_y = check_y%nrecv

      if(debug_message_passing) then
         allocate(msg1(0:nlist-1), msg2(0:nlist-1) )
         msg1 = 0
         msg2 = 0
         cur_rank = get_rank_recv(domain, check_x, check_y, rank_x, rank_y, ind_x, ind_y) 
         do while ( ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y )
            msgsize = 0
            if(cur_rank == rank_x) then
               from_pe = check_x%recv(ind_x)%pe
               do n = 1, check_x%recv(ind_x)%count
                  is = check_x%recv(ind_x)%is(n); ie = check_x%recv(ind_x)%ie(n)
                  js = check_x%recv(ind_x)%js(n); je = check_x%recv(ind_x)%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end do
               ind_x = ind_x+1
               if(ind_x .LE. nrecv_x) then
                  rank_x = check_x%recv(ind_x)%pe - domain%pe 
                  if(rank_x .LE.0) rank_x = rank_x + nlist
               else
                  rank_x = -1
               endif
            endif
            if(cur_rank == rank_y) then
               from_pe = check_y%recv(ind_y)%pe
               do n = 1, check_y%recv(ind_y)%count
                  is = check_y%recv(ind_y)%is(n); ie = check_y%recv(ind_y)%ie(n)
                  js = check_y%recv(ind_y)%js(n); je = check_y%recv(ind_y)%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end do
               ind_y = ind_y+1
               if(ind_y .LE. nrecv_y) then
                  rank_y = check_y%recv(ind_y)%pe - domain%pe 
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

         cur_rank = get_rank_send(domain, check_x, check_y, rank_x, rank_y, ind_x, ind_y) 
         do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
            msgsize = 0
            if(cur_rank == rank_x) then
               to_pe = check_x%send(ind_x)%pe
               do n = 1, check_x%send(ind_x)%count
                  is = check_x%send(ind_x)%is(n); ie = check_x%send(ind_x)%ie(n)
                  js = check_x%send(ind_x)%js(n); je = check_x%send(ind_x)%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               enddo
               ind_x = ind_x+1
               if(ind_x .LE. nsend_x) then
                  rank_x = check_x%send(ind_x)%pe - domain%pe 
                  if(rank_x .LT.0) rank_x = rank_x + nlist
               else
                  rank_x = nlist+1
               endif
            endif

            if(cur_rank == rank_y) then
               to_pe = check_y%send(ind_y)%pe
               do n = 1, check_y%send(ind_y)%count
                  is = check_y%send(ind_y)%is(n); ie = check_y%send(ind_y)%ie(n)
                  js = check_y%send(ind_y)%js(n); je = check_y%send(ind_y)%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end do
               ind_y = ind_y+1
               if(ind_y .LE. nsend_y) then
                  rank_y = check_y%send(ind_y)%pe - domain%pe 
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
               call mpp_error(FATAL, "mpp_do_checkV: mismatch on send and recv size")
            endif
         enddo

         call mpp_sync_self()
         write(stdout(),*)"NOTE from mpp_do_checkV: message sizes are matched between send and recv for domain " &
              //trim(domain%name)
         deallocate(msg1, msg2)
      endif

      !--- recv the data       
      cur_rank = get_rank_recv(domain, check_x, check_y, rank_x, rank_y, ind_x, ind_y) 

      do while ( ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y )
         msgsize = 0
         if(cur_rank == rank_x) then
            from_pe = check_x%recv(ind_x)%pe
            do n = 1, check_x%recv(ind_x)%count
               is = check_x%recv(ind_x)%is(n); ie = check_x%recv(ind_x)%ie(n)
               js = check_x%recv(ind_x)%js(n); je = check_x%recv(ind_x)%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            ind_x = ind_x+1
            if(ind_x .LE. nrecv_x) then
               rank_x = check_x%recv(ind_x)%pe - domain%pe 
               if(rank_x .LE.0) rank_x = rank_x + nlist
            else
               rank_x = -1
            endif
         endif
         if(cur_rank == rank_y) then
            from_pe = check_y%recv(ind_y)%pe
            do n = 1, check_y%recv(ind_y)%count
               is = check_y%recv(ind_y)%is(n); ie = check_y%recv(ind_y)%ie(n)
               js = check_y%recv(ind_y)%js(n); je = check_y%recv(ind_y)%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            ind_y = ind_y+1
            if(ind_y .LE. nrecv_y) then
               rank_y = check_y%recv(ind_y)%pe - domain%pe 
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
               call mpp_error( FATAL, 'MPP_DO_CHECK_V: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false. )
            buffer_pos = buffer_pos + msgsize
         end if
      enddo
      buffer_recv_size = buffer_pos

      !--- send the data
      cur_rank = get_rank_send(domain, check_x, check_y, rank_x, rank_y, ind_x, ind_y) 

      do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
         pos = buffer_pos
         if(cur_rank == rank_x) then
            to_pe = check_x%send(ind_x)%pe
            do n = 1, check_x%send(ind_x)%count
               is = check_x%send(ind_x)%is(n); ie = check_x%send(ind_x)%ie(n)
               js = check_x%send(ind_x)%js(n); je = check_x%send(ind_x)%je(n)
               tMe = check_x%send(ind_x)%tileMe(n)
               select case( check_x%send(ind_x)%rotation(n) )
               case(ZERO)
                  do l = 1, l_size ! loop over number of fields
                     ptr_fieldx = f_addrsx(l, tMe)
                     ptr_fieldy = f_addrsy(l, tMe)
                     do k = 1,ke  
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldx(i,j,k)
                           end do
                        end do
                     end do
                  end do
               case(MINUS_NINETY)
                  if( BTEST(update_flags,SCALAR_BIT) ) then
                     do l = 1, l_size ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
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
                     do l = 1, l_size ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
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
               case(NINETY)
                  do l = 1, l_size ! loop over number of fields
                     ptr_fieldx = f_addrsx(l, tMe)
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
               case(ONE_HUNDRED_EIGHTY) 
                  if( BTEST(update_flags,SCALAR_BIT) ) then
                     do l = 1, l_size ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
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
                     do l = 1, l_size ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
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
            end do
            ind_x = ind_x+1
            if(ind_x .LE. nsend_x) then
               rank_x = check_x%send(ind_x)%pe - domain%pe 
               if(rank_x .LT.0) rank_x = rank_x + nlist
            else
               rank_x = nlist+1
            endif
         endif

         if(cur_rank == rank_y) then
            to_pe = check_y%send(ind_y)%pe            
            do n = 1, check_y%send(ind_y)%count
               is = check_y%send(ind_y)%is(n); ie = check_y%send(ind_y)%ie(n)
               js = check_y%send(ind_y)%js(n); je = check_y%send(ind_y)%je(n)
               tMe = check_y%send(ind_y)%tileMe(n)
               select case( check_y%send(ind_y)%rotation(n) )
               case(ZERO)
                  do l = 1, l_size ! loop over number of fields
                     ptr_fieldx = f_addrsx(l, tMe)
                     ptr_fieldy = f_addrsy(l, tMe)
                     do k = 1,ke  
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldy(i,j,k)
                           end do
                        end do
                     end do
                  end do
               case(MINUS_NINETY)
                  do l = 1, l_size ! loop over number of fields
                     ptr_fieldx = f_addrsx(l, tMe)
                     ptr_fieldy = f_addrsy(l, tMe)
                     do k = 1, ke
                        do j = je, js, -1
                           do i = is, ie
                              pos = pos + 1
                              buffer(pos) = fieldx(i,j,k)
                           end do
                        end do
                     end do
                  end do
               case(NINETY)
                  if( BTEST(update_flags,SCALAR_BIT) ) then
                     do l = 1, l_size ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
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
                     do l = 1, l_size ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
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
               case(ONE_HUNDRED_EIGHTY) 
                  if( BTEST(update_flags,SCALAR_BIT) ) then
                     do l = 1, l_size ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
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
                     do l = 1, l_size ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
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
            end do
            ind_y = ind_y+1
            if(ind_y .LE. nsend_y) then
               rank_y = check_y%send(ind_y)%pe - domain%pe 
               if(rank_y .LT.0) rank_y = rank_y + nlist
            else
               rank_y = nlist+1
            endif
         endif
         cur_rank = min(rank_x, rank_y)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_CHECK_V: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
            buffer_pos = pos
         end if
      end do ! end do list = 0,nlist-1

      call mpp_sync_self(check=EVENT_RECV) ! To ensure recv is completed.
      buffer_pos = buffer_recv_size

      !--- compare the data in reverse order
      cur_rank = get_rank_unpack(domain, check_x, check_y, rank_x, rank_y, ind_x, ind_y) 

      CHECK_LOOP: do while(ind_x >0 .OR. ind_y >0)
         if(cur_rank == rank_y) then
            do n = check_y%recv(ind_y)%count, 1, -1
               is = check_y%recv(ind_y)%is(n); ie = check_y%recv(ind_y)%ie(n)
               js = check_y%recv(ind_y)%js(n); je = check_y%recv(ind_y)%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               tMe = check_y%recv(ind_y)%tileMe(n)
               do l=1,l_size  ! loop over number of fields
                  ptr_fieldx = f_addrsx(l, tMe)
                  ptr_fieldy = f_addrsy(l, tMe)
                  do k = 1,ke
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           if( fieldy(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_CHECK_V on pe = ", mpp_pe(), ": y component of vector ", &
                                   trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldy(i,j,k), &
                                   " does not equal to the value = ", buffer(pos), " on pe ", check_y%recv(ind_y)%pe
                              call mpp_error(debug_update_level, "MPP_DO_CHECK_V: mismatch on the boundary for symmetry point")
                              exit CHECK_LOOP
                           end if
                        end do
                     end do
                  end do
               end do
            end do
            ind_y = ind_y-1
            if(ind_y .GT. 0) then
               rank_y = check_y%recv(ind_y)%pe - domain%pe 
               if(rank_y .LE.0) rank_y = rank_y + nlist
            else
               rank_y = nlist+1
            endif
         endif

         if(cur_rank == rank_x) then
            do n = check_x%recv(ind_x)%count, 1, -1
               is = check_x%recv(ind_x)%is(n); ie = check_x%recv(ind_x)%ie(n)
               js = check_x%recv(ind_x)%js(n); je = check_x%recv(ind_x)%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               tMe = check_x%recv(ind_x)%tileMe(n)
               do l=1,l_size  ! loop over number of fields
                  ptr_fieldx = f_addrsx(l, tMe)
                  ptr_fieldy = f_addrsy(l, tMe)
                  do k = 1,ke
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           if( fieldx(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_CHECK_V on pe = ", mpp_pe(), ": x-component of vector ", &
                                   trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldx(i,j,k), &
                                   " does not equal to the value = ", buffer(pos), " on pe ", check_x%recv(ind_x)%pe
                              call mpp_error(debug_update_level, "MPP_DO_CHECK_V: mismatch on the boundary for symmetry point")
                              exit CHECK_LOOP
                           end if
                        end do
                     end do
                  end do
               end do
            end do
            ind_x = ind_x-1
            if(ind_x .GT. 0) then
               rank_x = check_x%recv(ind_x)%pe - domain%pe 
               if(rank_x .LE.0) rank_x = rank_x + nlist
            else
               rank_x = nlist+1
            endif
         endif
         cur_rank = min(rank_x, rank_y)
      end do CHECK_LOOP ! end do list = nlist-1,0,-1
      call mpp_sync_self()

      return

    end subroutine MPP_DO_CHECK_3D_V_
