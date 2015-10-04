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
subroutine MPP_START_DO_UPDATE_3D_(id_update, f_addrs, domain, update, d_type, ke_max, ke_list, flags, reuse_id_update, name)
  integer,                    intent(in) :: id_update
  integer(LONG_KIND),         intent(in) :: f_addrs(:,:)
  type(domain2D),             intent(in) :: domain
  type(overlapSpec),          intent(in) :: update
  MPP_TYPE_,                  intent(in) :: d_type  ! creates unique interface
  integer,                    intent(in) :: ke_max
  integer,                    intent(in) :: ke_list(:,:)
  logical,                    intent(in) :: reuse_id_update
  character(len=*),           intent(in) :: name
  integer,                    intent(in) :: flags
  
  !--- local variables
  integer                     :: i, j, k, m, n, l, dir, count, tMe
  integer                     :: buffer_pos, msgsize, from_pe, to_pe, pos
  integer                     :: is, ie, js, je
  logical                     :: send(8), recv(8), update_edge_only
  integer                     :: l_size, ke_sum
  character(len=128)          :: text
  type(overlap_type), pointer :: overPtr => NULL()    
  MPP_TYPE_                   :: buffer(size(mpp_domains_stack_nonblock(:)))
  MPP_TYPE_                   :: field(update%xbegin:update%xend, update%ybegin:update%yend,ke_max)
  pointer( ptr, buffer )
  pointer(ptr_field, field)

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

  l_size = size(f_addrs,1)
  ke_sum = sum(ke_list)
  ptr = LOC(mpp_domains_stack_nonblock)

  buffer_pos = nonblock_data(id_update)%recv_pos

  ! pre-postrecv
  do m = 1, update%nrecv
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(recv_clock_nonblock)
     msgsize = 0
     !--- make sure the domain stack size is big enough.
     do n = 1, overPtr%count
        dir = overPtr%dir(n)
        if(recv(dir)) then
           msgsize = msgsize + overPtr%msgsize(n)
        end if
     end do

     msgsize = msgsize*ke_sum
     if( msgsize.GT.0 )then
        from_pe = overPtr%pe
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_START_DO_UPDATE: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        count = nonblock_data(id_update)%request_recv_count + 1
        if( count > MAX_REQUEST ) then
           write( text,'(a,i8,a,i8)' ) 'recv request count =', count, ' greater than MAX_REQEUST =', MAX_REQUEST
           call mpp_error(FATAL,'MPP_START_DO_UPDATE: '//trim(text))
        endif
        nonblock_data(id_update)%request_recv_count = count
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., &
             tag=id_update, request=nonblock_data(id_update)%request_recv(count))
        nonblock_data(id_update)%size_recv(count) = msgsize
#ifdef use_libMPI
        nonblock_data(id_update)%type_recv(count) = MPI_TYPE_
#endif
        buffer_pos = buffer_pos + msgsize
     end if
     call mpp_clock_end(recv_clock_nonblock)
  end do ! end do m = 1, update%nrecv  

  msgsize = buffer_pos - nonblock_data(id_update)%recv_pos
  if( reuse_id_update ) then
     if(msgsize .NE. nonblock_data(id_update)%recv_msgsize) then
        call mpp_error(FATAL,'MPP_START_DO_UPDATE: mismatch of recv msgsize for field '//trim(name) )
     endif
  else
     nonblock_data(id_update)%recv_msgsize = msgsize
     nonblock_data(id_update)%send_pos = buffer_pos
     nonblock_buffer_pos = nonblock_buffer_pos + msgsize
  endif

  ! send
  do m = 1, update%nsend
     overPtr => update%send(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(pack_clock_nonblock)
     pos = buffer_pos

     ! make sure the stacksize is big enough
     msgsize = 0
     do n = 1, overPtr%count
        dir = overPtr%dir(n)
        if( send(dir) )  msgsize = msgsize + overPtr%msgsize(n)
     enddo
     if( msgsize.GT.0 )then
        msgsize = msgsize*ke_sum
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_START_DO_UPDATE: mpp_domains_stack overflow, ' // &
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
                 do k = 1,ke_list(l,tMe)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          buffer(pos) = field(i,j,k)
                       end do
                    end do
                 end do
              enddo
           else
              select case( overPtr%rotation(n) )
              case(ZERO)
                 do l=1,l_size  ! loop over number of fields
                    ptr_field = f_addrs(l, tMe)
                    do k = 1,ke_list(l,tMe)  
                       do j = js, je
                          do i = is, ie
                             pos = pos + 1
                             buffer(pos) = field(i,j,k)
                          end do
                       end do
                    end do
                 enddo
              case( MINUS_NINETY ) 
                 do l=1,l_size  ! loop over number of fields
                    ptr_field = f_addrs(l, tMe)
                    do k = 1,ke_list(l,tMe)  
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

                    do k = 1,ke_list(l,tMe)  
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

                    do k = 1,ke_list(l,tMe)  
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

     call mpp_clock_end(pack_clock_nonblock)
     call mpp_clock_begin(send_clock_nonblock)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        to_pe = overPtr%pe
        count = nonblock_data(id_update)%request_send_count + 1
        if( count > MAX_REQUEST ) then
           write( text,'(a,i8,a,i8)' ) 'send request count =', count, ' greater than MAX_REQEUST =', MAX_REQUEST
           call mpp_error(FATAL,'MPP_START_DO_UPDATE: '//trim(text))
        endif
        nonblock_data(id_update)%request_send_count = count
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, &
                       tag=id_update, request=nonblock_data(id_update)%request_send(count))
        buffer_pos = pos
     end if
     call mpp_clock_end(send_clock_nonblock)
  end do ! end do ist = 0,nlist-1

  msgsize = buffer_pos - nonblock_data(id_update)%send_pos
  if( reuse_id_update ) then
     if(msgsize .NE. nonblock_data(id_update)%send_msgsize) then
        call mpp_error(FATAL,'MPP_START_DO_UPDATE: mismatch of send msgsize for field '//trim(name) )
     endif
  else
     nonblock_buffer_pos = nonblock_buffer_pos + msgsize
     nonblock_data(id_update)%send_msgsize = msgsize
  endif
  
  overPtr => NULL()

  return


end subroutine MPP_START_DO_UPDATE_3D_

!###############################################################################

subroutine MPP_COMPLETE_DO_UPDATE_3D_(id_update, f_addrs, domain, update, d_type, ke_max, ke_list, &
                                      b_addrs, b_size, flags) 
  integer,             intent(in) :: id_update
  integer(LONG_KIND),  intent(in) :: f_addrs(:,:)
  type(domain2d),      intent(in) :: domain
  type(overlapSpec),   intent(in) :: update
  integer,             intent(in) :: ke_max
  integer,             intent(in) :: ke_list(:,:)
  MPP_TYPE_,           intent(in) :: d_type  ! creates unique interface
  integer(LONG_KIND),  intent(in) :: b_addrs(:,:)
  integer,             intent(in) :: b_size
  integer,             intent(in) :: flags

  !--- local variables
  integer                     :: i, j, k, m, n, l, dir, count, tMe
  integer                     :: buffer_pos, msgsize, from_pe, pos
  integer                     :: is, ie, js, je
  integer                     :: start, start1, start2, index
  integer                     :: is1, ie1, js1, je1, ni, nj, total
  logical                     :: send(8), recv(8), update_edge_only
  integer                     :: l_size, ke_sum
  character(len=128)          :: text
  type(overlap_type), pointer :: overPtr => NULL()
  MPP_TYPE_                   :: recv_buffer(size(mpp_domains_stack_nonblock(:)))
  MPP_TYPE_                   :: field(update%xbegin:update%xend, update%ybegin:update%yend,ke_max)
  MPP_TYPE_                   :: buffer(b_size)
  pointer( ptr, recv_buffer )
  pointer(ptr_field, field)
  pointer(ptr_buffer, buffer) 


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

  ke_sum = sum(ke_list)
  l_size = size(f_addrs,1)
  ptr = LOC(mpp_domains_stack_nonblock)

  count = nonblock_data(id_update)%request_recv_count
  if(count > 0) then
     call mpp_clock_begin(wait_clock_nonblock)
     call mpp_sync_self(check=EVENT_RECV, request=nonblock_data(id_update)%request_recv(1:count), &
                        msg_size=nonblock_data(id_update)%size_recv(1:count),                     &
                        msg_type=nonblock_data(id_update)%type_recv(1:count) )
     call mpp_clock_end(wait_clock_nonblock)
     nonblock_data(id_update)%request_recv_count = 0
#ifdef use_libMPI
     nonblock_data(id_update)%request_recv(:)    = MPI_REQUEST_NULL
#else
     nonblock_data(id_update)%request_recv(:)    = 0
#endif
     nonblock_data(id_update)%size_recv(:) = 0
     nonblock_data(id_update)%type_recv(:) = 0
  endif 

  buffer_pos = nonblock_data(id_update)%recv_pos + nonblock_data(id_update)%recv_msgsize
  !--unpack the data
  call mpp_clock_begin(unpk_clock_nonblock)
  do m = update%nrecv, 1, -1
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle

     pos = buffer_pos
     do n = overPtr%count, 1, -1
        dir = overPtr%dir(n)
        if( recv(dir) ) then
           tMe = overPtr%tileMe(n)
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke_sum
           pos = buffer_pos - msgsize
           buffer_pos = pos
           if(OverPtr%is_refined(n)) then
              index = overPtr%index(n)
              is1 = update%rSpec(tMe)%isNbr(index); ie1 = update%rSpec(tMe)%ieNbr(index)
              js1 = update%rSpec(tMe)%jsNbr(index); je1 = update%rSpec(tMe)%jeNbr(index)
              ni = ie1 - is1 + 1
              nj = je1 - js1 + 1
              total = ni*nj
              start = (update%rSpec(tMe)%start(index)-1)*ke_max
              if(start+total*ke_max>size(buffer) ) call mpp_error(FATAL, &
                   "MPP_COMPETE_UPDATE_DOMAINS: size of buffer is less than the size of the data to be filled.")
              msgsize = ie - is + 1
              do l=1, l_size  ! loop over number of fields
                 ptr_buffer = b_addrs(l, tMe)
                 if(l==1) start = (update%rSpec(tMe)%start(index)-1)*ke_list(l,tMe)
                 start1 = start + (js-js1)*ni + is - is1
                 do k = 1, ke_list(l,tMe)
                    start2 = start1
                    do j = js, je
                       buffer(start2+1:start2+msgsize) = recv_buffer(pos+1:pos+msgsize)
                       start2 = start2 + ni
                       pos   = pos + msgsize
                    end do
                    start1 = start1 + total
                 end do
              enddo
           else
              do l=1, l_size  ! loop over number of fields
                 ptr_field = f_addrs(l, tMe)
                 do k = 1,ke_list(l,tMe)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          field(i,j,k) = recv_buffer(pos)
                       end do
                    end do
                 end do
              end do
           endif
        end if
     end do ! do n = 1, overPtr%count
  end do

  call mpp_clock_end(unpk_clock_nonblock)

  count = nonblock_data(id_update)%request_send_count
  if(count > 0) then
     call mpp_clock_begin(wait_clock_nonblock)
     call mpp_sync_self(check=EVENT_SEND, request=nonblock_data(id_update)%request_send(1:count))
     call mpp_clock_end(wait_clock_nonblock)
     nonblock_data(id_update)%request_send_count = 0
#ifdef use_libMPI
     nonblock_data(id_update)%request_send(:)    = MPI_REQUEST_NULL
#else
     nonblock_data(id_update)%request_send(:)    = 0
#endif
  endif 
  
!  call init_nonblock_type(nonblock_data(id_update))

  return

end subroutine MPP_COMPLETE_DO_UPDATE_3D_
