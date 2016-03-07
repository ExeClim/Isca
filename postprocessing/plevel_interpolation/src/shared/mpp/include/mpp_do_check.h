! -*-f90-*- 
    subroutine MPP_DO_CHECK_3D_( f_addrs, domain, check, d_type, ke, flags, name)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND),         intent(in) :: f_addrs(:,:)
      type(domain2D),             intent(in) :: domain
      type(overlapSpec),          intent(in) :: check
      MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
      integer,                    intent(in) :: ke
      integer, optional,          intent(in) :: flags
      character(len=*), optional, intent(in) :: name

      MPP_TYPE_ :: field(check%xbegin:check%xend, check%ybegin:check%yend,ke)
      pointer(ptr_field, field)
      integer                     :: update_flags
      character(len=8)            :: text
      character(len=64)           :: field_name      

!equate to mpp_domains_stack
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer( ptr, buffer )
      integer :: buffer_pos
      integer,    allocatable :: msg1(:), msg2(:)     
!receive domains saved here for unpacking
!for non-blocking version, could be recomputed
      integer :: to_pe, from_pe, pos, msgsize
      integer :: n, l_size, l, m, i, j, k
      integer :: is, ie, js, je, tMe
      integer :: buffer_recv_size, nlist

      ptr = LOC(mpp_domains_stack)
      l_size = size(f_addrs,1)

      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) )update_flags = flags

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

      if(debug_message_passing) then
         nlist = size(domain%list(:))  
         allocate(msg1(0:nlist-1), msg2(0:nlist-1) )
         msg1 = 0
         msg2 = 0
         do m = 1, check%nrecv
            msgsize = 0
            do n = 1, check%recv(m)%count
               is = check%recv(m)%is(n); ie = check%recv(m)%ie(n)
               js = check%recv(m)%js(n); je = check%recv(m)%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            from_pe = check%recv(m)%pe
            l = from_pe-mpp_root_pe()
            call mpp_recv( msg1(l), glen=1, from_pe=from_pe, block=.FALSE.)
            msg2(l) = msgsize
         enddo

         do m = 1, check%nsend
            msgsize = 0
            do n = 1, check%send(m)%count
               is = check%send(m)%is(n); ie = check%send(m)%ie(n)
               js = check%send(m)%js(n); je = check%send(m)%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            call mpp_send(msgsize, plen=1, to_pe=check%send(m)%pe)
         enddo
         call mpp_sync_self(check=EVENT_RECV)

         do m = 0, nlist-1
            if(msg1(m) .NE. msg2(m)) then
               print*, "My pe = ", mpp_pe(), ",domain name =", trim(domain%name), ",from pe=", &
                    domain%list(m)%pe, ":send size = ", msg1(m), ", recv size = ", msg2(m)
               call mpp_error(FATAL, "mpp_do_check: mismatch on send and recv size")
            endif
         enddo
         call mpp_sync_self()
         write(stdout(),*)"NOTE from mpp_do_check: message sizes are matched between send and recv for domain " &
              //trim(domain%name)
         deallocate(msg1, msg2)
      endif

      buffer_pos = 0        
      !--- pre-post recv the data 
      do m = 1, check%nrecv
         msgsize = 0
         do n = 1, check%recv(m)%count
            is = check%recv(m)%is(n); ie = check%recv(m)%ie(n)
            js = check%recv(m)%js(n); je = check%recv(m)%je(n)
            msgsize = msgsize + (ie-is+1)*(je-js+1)
         end do
         msgsize = msgsize*ke*l_size

         if( msgsize.GT.0 )then
            from_pe = check%recv(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_CHECK: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE. )
            buffer_pos = buffer_pos + msgsize
         end if
      end do
      buffer_recv_size = buffer_pos

      !--- send the data
      do m = 1, check%nsend
         pos = buffer_pos
         do n = 1, check%recv(m)%count
            is = check%recv(m)%is(n); ie = check%recv(m)%ie(n)
            js = check%recv(m)%js(n); je = check%recv(m)%je(n)
            tMe = check%recv(m)%tileMe(n)
            select case( check%recv(m)%rotation(n) )
            case(ZERO)
               do l = 1, l_size ! loop over number of fields
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
            case(MINUS_NINETY)
               do l = 1, l_size ! loop over number of fields
                  ptr_field = f_addrs(l, tMe)
                  do k = 1,ke  
                     do j = je, js, -1
                        do i = is, ie

                           pos = pos + 1
                           buffer(pos) = field(i,j,k)
                        end do
                     end do
                  end do
               end do
            case(NINETY)
               do l = 1, l_size ! loop over number of fields
                  ptr_field = f_addrs(l, tMe)
                  do k = 1,ke  
                     do j = js, je
                        do i = ie, is, -1

                           pos = pos + 1
                           buffer(pos) = field(i,j,k)
                        end do
                     end do
                  end do
               end do
            case(ONE_HUNDRED_EIGHTY)
               do l = 1, l_size ! loop over number of fields
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
         end do
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = check%recv(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_CHECK: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
            buffer_pos = pos
         end if
      end do ! end do list = 0,nlist-1

      call mpp_sync_self(check=EVENT_RECV) ! To ensure recv is completed.
      buffer_pos = buffer_recv_size
      !--- compare the data in reverse order
      CHECK_LOOP: do m = check%nrecv, 1, -1
         do n = check%recv(m)%count, 1, -1
            is = check%recv(m)%is(n); ie = check%recv(m)%ie(n)
            js = check%recv(m)%js(n); je = check%recv(m)%je(n)
            msgsize = (ie-is+1)*(je-js+1)*ke*l_size
            pos = buffer_pos - msgsize
            buffer_pos = pos
            tMe = check%recv(m)%tileMe(n)
            do l=1, l_size  ! loop over number of fields
               ptr_field = f_addrs(l, tMe)
               do k = 1,ke
                  do j = js, je
                     do i = is, ie
                        pos = pos + 1
                        if( field(i,j,k) .NE. buffer(pos) ) then
                           print*,"Error from MPP_DO_CHECK on pe = ", mpp_pe(), ": field ", &
                                trim(field_name), " at point (", i, ",", j, ",", k, ") = ", field(i,j,k), &
                                " does not equal to the value = ", buffer(pos), " on pe ", check%recv(m)%pe
                           call mpp_error(debug_update_level, "MPP_DO_CHECK: mismatch on the boundary for symmetry point")
                           exit CHECK_LOOP
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do CHECK_LOOP ! end do list = nlist-1,0,-1
      call mpp_sync_self()

      return
    end subroutine MPP_DO_CHECK_3D_
