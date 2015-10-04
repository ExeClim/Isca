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
subroutine MPP_START_DO_UPDATE_3D_V_(id_update, f_addrsx, f_addrsy, domain, update_x, update_y,     &
                                     d_type, ke_max, ke_list, gridtype, flags, reuse_id_update, name)
  integer,             intent(in) :: id_update
  integer(LONG_KIND),  intent(in) :: f_addrsx(:,:), f_addrsy(:,:)
  type(domain2d),      intent(in) :: domain
  type(overlapSpec),   intent(in) :: update_x, update_y
  integer,             intent(in) :: ke_max
  integer,             intent(in) :: ke_list(:,:)
  MPP_TYPE_,           intent(in) :: d_type  ! creates unique interface
  integer,             intent(in) :: gridtype
  logical,             intent(in) :: reuse_id_update
  character(len=*),    intent(in) :: name
  integer,             intent(in) :: flags

  !---local variable ------------------------------------------
  integer            :: i, j, k, l, is, ie, js, je, n
  integer            :: pos, nlist, msgsize, tile, l_size
  integer            :: to_pe, from_pe, buffer_pos
  integer            :: tMe, dir, count, ke_sum
  logical            :: send(8), recv(8), update_edge_only
  character(len=128) :: text
  integer            :: rank_x, rank_y, ind_x, ind_y, cur_rank
  integer            :: nsend_x, nsend_y, nrecv_x, nrecv_y
  MPP_TYPE_          :: fieldx(update_x%xbegin:update_x%xend, update_x%ybegin:update_x%yend,ke_max)
  MPP_TYPE_          :: fieldy(update_y%xbegin:update_y%xend, update_y%ybegin:update_y%yend,ke_max)
  MPP_TYPE_          :: buffer(size(mpp_domains_stack_nonblock(:)))

  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)
  pointer( ptr, buffer )

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
  l_size = size(f_addrsx,1)
  nlist  = size(domain%list(:))
  ptr    = LOC(mpp_domains_stack_nonblock)

  !recv
  nsend_x = update_x%nsend
  nsend_y = update_y%nsend
  nrecv_x = update_x%nrecv
  nrecv_y = update_y%nrecv

  !--- recv
  cur_rank = get_rank_recv(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 

  buffer_pos = nonblock_data(id_update)%recv_pos
  call mpp_clock_begin(recv_clock_nonblock)
  do while (ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y)
     msgsize = 0
     select case(gridtype)
     case(BGRID_NE, BGRID_SW, AGRID)
        if(cur_rank == rank_x) then
           from_pe = update_x%recv(ind_x)%pe
           do n = 1, update_x%recv(ind_x)%count
              dir = update_x%recv(ind_x)%dir(n)
              if(recv(dir)) then
                 msgsize = msgsize + update_x%recv(ind_x)%msgsize(n)
              end if
           end do
           msgsize = msgsize*2
           ind_x = ind_x+1
           ind_y = ind_x
           if(ind_x .LE. nrecv_x) then
              rank_x = update_x%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = -1
           endif
           rank_y = rank_x
        endif
     case(CGRID_NE, CGRID_SW)
        if(cur_rank == rank_x) then
           from_pe = update_x%recv(ind_x)%pe
           do n = 1, update_x%recv(ind_x)%count
              dir = update_x%recv(ind_x)%dir(n)
              if(recv(dir)) then
                 msgsize = msgsize + update_x%recv(ind_x)%msgsize(n)
              end if
           end do
           ind_x = ind_x+1
           if(ind_x .LE. nrecv_x) then
              rank_x = update_x%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = -1
           endif
        endif
        if(cur_rank == rank_y) then
           from_pe = update_y%recv(ind_y)%pe
           do n = 1, update_y%recv(ind_y)%count
              dir = update_y%recv(ind_y)%dir(n)
              if(recv(dir)) then
                 msgsize = msgsize + update_y%recv(ind_y)%msgsize(n)
              end if
           end do
           ind_y = ind_y+1
           if(ind_y .LE. nrecv_y) then
              rank_y = update_y%recv(ind_y)%pe - domain%pe 
              if(rank_y .LE.0) rank_y = rank_y + nlist
           else
              rank_y = -1
           endif
        endif
     end select
     cur_rank = max(rank_x, rank_y)
     msgsize = msgsize*ke_sum

     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, buffer_pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_START_DO_UPDATE_V: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        count = nonblock_data(id_update)%request_recv_count + 1
        if( count > MAX_REQUEST ) then
           write( text,'(a,i8,a,i8)' ) 'request count =', count, ' greater than MAX_REQEUST =', MAX_REQUEST
           call mpp_error(FATAL,'MPP_START_DO_UPDATE_V: '//trim(text))
        endif
        nonblock_data(id_update)%request_recv_count = count

        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., &
                   tag=id_update, request=nonblock_data(id_update)%request_recv(count))
        nonblock_data(id_update)%size_recv(count) = msgsize
#ifdef use_libMPI
        nonblock_data(id_update)%type_recv(count) = MPI_TYPE_
#endif
        buffer_pos = buffer_pos + msgsize
     end if
  end do
  call mpp_clock_end(recv_clock_nonblock)
  msgsize = buffer_pos - nonblock_data(id_update)%recv_pos
  if( reuse_id_update ) then
     if(msgsize .NE. nonblock_data(id_update)%recv_msgsize) then
        call mpp_error(FATAL,'MPP_START_UPDATE_DOMAINS_V: mismatch of recv msgsize for field '//trim(name) )
     endif
  else
     nonblock_data(id_update)%recv_msgsize = msgsize
     nonblock_data(id_update)%send_pos = buffer_pos
     nonblock_buffer_pos = nonblock_buffer_pos + msgsize
  endif

  !--- send
  cur_rank = get_rank_send(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 

  do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
     call mpp_clock_begin(pack_clock_nonblock)
     pos = buffer_pos
     !--- make sure the domain stack size is big enough
     msgsize = 0
     if(cur_rank == rank_x) then
        do n = 1, update_x%send(ind_x)%count
           dir = update_x%send(ind_x)%dir(n)
           if( send(dir) ) msgsize = msgsize +  update_x%send(ind_x)%msgsize(n)
        enddo
     endif
     if(cur_rank == rank_y) then
        do n = 1, update_y%send(ind_y)%count
           dir = update_y%send(ind_y)%dir(n)
           if( send(dir) ) msgsize = msgsize +  update_y%send(ind_y)%msgsize(n)
        enddo
     endif

     if( msgsize.GT.0 )then
        msgsize = msgsize*ke_sum
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_START_DO_UPDATE_V: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if
     
     select case( gridtype )
     case(BGRID_NE, BGRID_SW, AGRID)
        if(cur_rank == rank_x) then
           to_pe = update_x%send(ind_x)%pe
           do n = 1, update_x%send(ind_x)%count
              dir = update_x%send(ind_x)%dir(n)
              if( send(dir) ) then 
                 tMe = update_x%send(ind_x)%tileMe(n)
                 is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                 js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                 if( update_x%send(ind_x)%is_refined(n) ) then
                    select case( update_x%send(ind_x)%rotation(n) )
                    case(ZERO)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 2
                                   buffer(pos-1) = fieldx(i,j,k)
                                   buffer(pos)   = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 2
                                      buffer(pos-1) =  fieldy(i,j,k)
                                      buffer(pos)   =  fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 2
                                      buffer(pos-1) = -fieldy(i,j,k)
                                      buffer(pos)   =  fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    case(NINETY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 2
                                      buffer(pos-1) =  fieldy(i,j,k)
                                      buffer(pos)   =  fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 2
                                      buffer(pos-1) =  fieldy(i,j,k)
                                      buffer(pos)   = -fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    case(ONE_HUNDRED_EIGHTY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 2
                                      buffer(pos-1) = fieldx(i,j,k)
                                      buffer(pos)   = fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 2
                                      buffer(pos-1) = -fieldx(i,j,k)
                                      buffer(pos)   = -fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    end select  ! select case( rotation(n) )
                 else   ! if( is_refined(n) )
                    select case( update_x%send(ind_x)%rotation(n) )
                    case(ZERO)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 2
                                   buffer(pos-1) = fieldx(i,j,k)
                                   buffer(pos)   = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case( MINUS_NINETY ) 
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do i = is, ie
                                   do j = je, js, -1
                                      pos = pos + 2
                                      buffer(pos-1) = fieldy(i,j,k)
                                      buffer(pos)   = fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do i = is, ie
                                   do j = je, js, -1
                                      pos = pos + 2
                                      buffer(pos-1) = -fieldy(i,j,k)
                                      buffer(pos)   =  fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    case( NINETY )
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do i = ie, is, -1
                                   do j = js, je
                                      pos = pos + 2
                                      buffer(pos-1) = fieldy(i,j,k)
                                      buffer(pos)   = fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do i = ie, is, -1
                                   do j = js, je
                                      pos = pos + 2
                                      buffer(pos-1) = fieldy(i,j,k)
                                      buffer(pos)   = -fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    case( ONE_HUNDRED_EIGHTY )
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = je, js, -1
                                   do i = ie, is, -1
                                      pos = pos + 2
                                      buffer(pos-1) =  fieldx(i,j,k)
                                      buffer(pos)   =  fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = je, js, -1
                                   do i = ie, is, -1
                                      pos = pos + 2
                                      buffer(pos-1) =  -fieldx(i,j,k)
                                      buffer(pos)   =  -fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    end select ! select case( rotation(n) )
                 end if ! if( is_refined(n) )
              end if ! if( send(dir) ) 
           end do ! do n = 1, update_x%send(ind_x)%count
           ind_x = ind_x+1
           ind_y = ind_x
           if(ind_x .LE. nsend_x) then
              rank_x = update_x%send(ind_x)%pe - domain%pe 
              if(rank_x .LT.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1
           endif
           rank_y = rank_x
        endif
     case(CGRID_NE, CGRID_SW)
        if(cur_rank == rank_x) then
           to_pe = update_x%send(ind_x)%pe
           do n = 1, update_x%send(ind_x)%count
              dir = update_x%send(ind_x)%dir(n)
              if( send(dir) ) then
                 tMe = update_x%send(ind_x)%tileMe(n)
                 is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                 js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                 if( update_x%send(ind_x)%is_refined(n) ) then
                    select case( update_x%send(ind_x)%rotation(n) )
                    case(ZERO)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 1
                                      buffer(pos) = fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 1
                                      buffer(pos) = -fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    case(NINETY)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1, ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case(ONE_HUNDRED_EIGHTY)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end select
                 else
                    select case( update_x%send(ind_x)%rotation(n) )
                    case(ZERO)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do i = is, ie
                                   do j = je, js, -1
                                      pos = pos + 1
                                      buffer(pos) = fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do i = is, ie
                                   do j = je, js, -1
                                      pos = pos + 1
                                      buffer(pos) = -fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    case(NINETY)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1, ke_list(l,tMe)
                             do i = ie, is, -1
                                do j = js, je
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case(ONE_HUNDRED_EIGHTY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = je, js, -1
                                   do i = ie, is, -1
                                      pos = pos + 1
                                      buffer(pos) = fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
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
                 end if
              end if
           end do
           ind_x = ind_x+1
           if(ind_x .LE. nsend_x) then
              rank_x = update_x%send(ind_x)%pe - domain%pe 
              if(rank_x .LT.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1 
           endif
        endif
        if(cur_rank == rank_y) then
           to_pe = update_y%send(ind_y)%pe
           do n = 1, update_y%send(ind_y)%count
              dir = update_y%send(ind_y)%dir(n)
              if( send(dir) ) then
                 tMe = update_y%send(ind_y)%tileMe(n)
                 is = update_y%send(ind_y)%is(n); ie = update_y%send(ind_y)%ie(n)
                 js = update_y%send(ind_y)%js(n); je = update_y%send(ind_y)%je(n)
                 if( update_y%send(ind_y)%is_refined(n) ) then
                    select case( update_y%send(ind_y)%rotation(n) )
                    case(ZERO)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       enddo
                    case(MINUS_NINETY)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case(NINETY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 1
                                      buffer(pos) = fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = js, je
                                   do i = is, ie
                                      pos = pos + 1
                                      buffer(pos) = -fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    case(ONE_HUNDRED_EIGHTY)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    end select
                 else
                    select case( update_y%send(ind_y)%rotation(n) )
                    case(ZERO)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do j = js, je
                                do i = is, ie
                                   pos = pos + 1
                                   buffer(pos) = fieldy(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case(MINUS_NINETY)
                       do l=1,l_size  ! loop over number of fields
                          ptr_fieldx = f_addrsx(l,tMe)
                          ptr_fieldy = f_addrsy(l,tMe)
                          do k = 1,ke_list(l,tMe)
                             do i = is, ie
                                do j = je, js, -1
                                   pos = pos + 1
                                   buffer(pos) = fieldx(i,j,k)
                                end do
                             end do
                          end do
                       end do
                    case(NINETY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do i = ie, is, -1
                                   do j = js, je
                                      pos = pos + 1
                                      buffer(pos) = fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do i = ie, is, -1
                                   do j = js, je
                                      pos = pos + 1
                                      buffer(pos) = -fieldx(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       end if
                    case(ONE_HUNDRED_EIGHTY)
                       if( BTEST(flags,SCALAR_BIT) ) then
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
                                do j = je, js, -1
                                   do i = ie, is, -1
                                      pos = pos + 1
                                      buffer(pos) = fieldy(i,j,k)
                                   end do
                                end do
                             end do
                          end do
                       else
                          do l=1,l_size  ! loop over number of fields
                             ptr_fieldx = f_addrsx(l,tMe)
                             ptr_fieldy = f_addrsy(l,tMe)
                             do k = 1,ke_list(l,tMe)
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
                 end if
              endif
           enddo
           ind_y = ind_y+1
           if(ind_y .LE. nsend_y) then
              rank_y = update_y%send(ind_y)%pe - domain%pe 
              if(rank_y .LT.0) rank_y = rank_y + nlist
           else
              rank_y = nlist+1
           endif
        endif
     end select
     call mpp_clock_end(pack_clock_nonblock)
     call mpp_clock_begin(send_clock_nonblock)
     cur_rank = min(rank_x, rank_y)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        count = nonblock_data(id_update)%request_send_count + 1
        if( count > MAX_REQUEST ) then
           write( text,'(a,i8,a,i8)' ) 'send request count =', count, ' greater than MAX_REQEUST =', MAX_REQUEST
           call mpp_error(FATAL,'MPP_START_DO_UPDATE_V: '//trim(text))
        endif
        nonblock_data(id_update)%request_send_count = count
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, &
                       tag=id_update, request=nonblock_data(id_update)%request_send(count) )
        buffer_pos = pos
     end if
     call mpp_clock_end(send_clock_nonblock)
  end do

  msgsize = buffer_pos - nonblock_data(id_update)%send_pos
  if( reuse_id_update ) then
     if(msgsize .NE. nonblock_data(id_update)%send_msgsize) then
        call mpp_error(FATAL,'MPP_START_DO_UPDATE_V: mismatch of send msgsize for field '//trim(name) )
     endif
  else
     nonblock_buffer_pos = nonblock_buffer_pos + msgsize
     nonblock_data(id_update)%send_msgsize = msgsize
  endif


end subroutine MPP_START_DO_UPDATE_3D_V_

!###############################################################################
subroutine MPP_COMPLETE_DO_UPDATE_3D_V_(id_update, f_addrsx, f_addrsy, domain, update_x, update_y,     &
                                        d_type, ke_max, ke_list, b_addrsx, b_addrsy, b_sizex, b_sizey, &
                                        gridtype, flags) 
  integer,             intent(in) :: id_update
  integer(LONG_KIND),  intent(in) :: f_addrsx(:,:), f_addrsy(:,:)
  type(domain2d),      intent(in) :: domain
  type(overlapSpec),   intent(in) :: update_x, update_y
  integer,             intent(in) :: ke_max
  integer,             intent(in) :: ke_list(:,:)
  MPP_TYPE_,           intent(in) :: d_type  ! creates unique interface
  integer(LONG_KIND),  intent(in) :: b_addrsx(:,:), b_addrsy(:,:)
  integer,             intent(in) :: b_sizex, b_sizey
  integer,             intent(in) :: gridtype
  integer,             intent(in) :: flags


  !--- local variables
  MPP_TYPE_ :: fieldx(update_x%xbegin:update_x%xend, update_x%ybegin:update_x%yend,ke_max)
  MPP_TYPE_ :: fieldy(update_y%xbegin:update_y%xend, update_y%ybegin:update_y%yend,ke_max)
  MPP_TYPE_ :: bufferx(b_sizex)
  MPP_TYPE_ :: buffery(b_sizey)
  pointer(ptr_fieldx, fieldx)
  pointer(ptr_fieldy, fieldy)
  pointer(ptr_bufferx, bufferx)
  pointer(ptr_buffery, buffery)


  MPP_TYPE_ :: recv_buffer(size(mpp_domains_stack_nonblock(:)))
  pointer( ptr, recv_buffer )

  integer :: i, j, k, l, is, ie, js, je, n, ke_sum, l_size
  integer :: pos, nlist, msgsize, tile, buffer_pos
  integer :: rank_x, rank_y, ind_x, ind_y, cur_rank
  integer :: index, is1, ie1, js1, je1, ni, nj, total, start1, start, start2
  logical :: recv(8), send(8), update_edge_only
  integer :: shift, midpoint
  integer :: tMe, dir, count

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
  l_size = size(f_addrsx,1)
  nlist  = size(domain%list(:))
  ptr = LOC(mpp_domains_stack_nonblock)

  buffer_pos = nonblock_data(id_update)%recv_pos + nonblock_data(id_update)%recv_msgsize
  cur_rank = get_rank_unpack(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 

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
     nonblock_data(id_update)%size_recv(:)    = 0
     nonblock_data(id_update)%type_recv(:)    = 0
  endif 

  do while (ind_x > 0 .OR. ind_y > 0)
     call mpp_clock_begin(unpk_clock_nonblock)
     pos = buffer_pos
     select case ( gridtype )
     case(BGRID_NE, BGRID_SW, AGRID)
        if(cur_rank == rank_x) then
           do n = update_x%recv(ind_x)%count, 1, -1    
              dir = update_x%recv(ind_x)%dir(n)
              if( recv(dir) ) then
                 tMe = update_x%recv(ind_x)%tileMe(n)
                 is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                 js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n) 
                 msgsize = (ie-is+1)*(je-js+1)*ke_sum*2
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 if(update_x%recv(ind_x)%is_refined(n)) then
                    index = update_x%recv(ind_x)%index(n)
                    is1 = update_x%rSpec(tMe)%isNbr(index); ie1 = update_x%rSpec(tMe)%ieNbr(index)
                    js1 = update_x%rSpec(tMe)%jsNbr(index); je1 = update_x%rSpec(tMe)%jeNbr(index)
                    ni = ie1 - is1 + 1
                    nj = je1 - js1 + 1
                    total = ni*nj
                    start = (update_x%rSpec(tMe)%start(index)-1)*ke_max

                    if(start+total*ke_max>size(bufferx) .or. start+total*ke_max>size(buffery) ) call mpp_error(FATAL, &
                         "MPP_COMPLETE_DO_UPDATE_3D_V: size of bufferx or buffery is less than the size of the data to be filled.")
                    msgsize = ie - is + 1
                    do l=1, l_size  ! loop over number of fields
                       ptr_bufferx = b_addrsx(l, tMe)
                       ptr_buffery = b_addrsy(l, tMe)
                       if(l==1) start = (update_x%rSpec(tMe)%start(index)-1)*ke_list(l,tMe)
                       start1 = start + (js-js1)*ni + is - is1 
                       do k = 1, ke_list(l,tMe)
                          start2 = start1
                          do j = js, je
                             do i = start2+1, start2+msgsize
                                pos = pos + 2
                                bufferx(i) = recv_buffer(pos-1)
                                buffery(i) = recv_buffer(pos)
                             end do
                             start2 = start2 + ni
                          end do
                          start1 = start1 + total
                       end do
                    enddo
                 else
                    do l=1, l_size  ! loop over number of fields
                       ptr_fieldx = f_addrsx(l, tMe)
                       ptr_fieldy = f_addrsy(l, tMe)
                       do k = 1,ke_list(l,tMe)
                          do j = js, je
                             do i = is, ie
                                pos = pos + 2
                                fieldx(i,j,k) = recv_buffer(pos-1)
                                fieldy(i,j,k) = recv_buffer(pos)
                             end do
                          end do
                       enddo
                    end do
                 end if
              end if ! end if( recv(dir) )
           end do  ! do dir=8,1,-1 
           ind_x = ind_x-1
           ind_y = ind_x
           if(ind_x .GT. 0) then
              rank_x = update_x%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1
           endif
           rank_y = rank_x
        endif
     case(CGRID_NE, CGRID_SW)
        if(cur_rank == rank_y) then
           do n = update_y%recv(ind_y)%count, 1, -1

              dir = update_y%recv(ind_y)%dir(n)
              if( recv(dir) ) then
                 tMe = update_y%recv(ind_y)%tileMe(n)
                 is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
                 js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
                 msgsize = (ie-is+1)*(je-js+1)*ke_sum
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 if(update_y%recv(ind_y)%is_refined(n)) then
                    index = update_y%recv(ind_y)%index(n)
                    is1 = update_y%rSpec(tMe)%isNbr(index); ie1 = update_y%rSpec(tMe)%ieNbr(index)
                    js1 = update_y%rSpec(tMe)%jsNbr(index); je1 = update_y%rSpec(tMe)%jeNbr(index)
                    ni = ie1 - is1 + 1
                    nj = je1 - js1 + 1
                    total = ni*nj
                    start = (update_y%rSpec(tMe)%start(index)-1)*ke_max
                    if(start+total*ke_max>size(buffery) ) call mpp_error(FATAL, &
                      "MPP_COMPLETE_DO_UPDATE_3D_V: size of buffery is less than the size of the data to be filled.")

                    msgsize = ie - is + 1
                    do l = 1, l_size
                       ptr_buffery = b_addrsy(l, tMe)
                       if(l==1) start = (update_y%rSpec(tMe)%start(index)-1)*ke_list(l,tMe)
                       start1 = start + (js-js1)*ni + is - is1 
                       do k = 1, ke_list(l,tMe)
                          start2 = start1
                          do j = js, je
                             do i = start2+1, start2+msgsize
                                pos = pos + 1
                                buffery(i) = recv_buffer(pos)
                             end do
                             start2 = start2 + ni
                          end do
                          start1 = start1 + total
                       end do
                    enddo
                 else
                    do l=1, l_size  ! loop over number of fields
                       ptr_fieldy = f_addrsy(l, tMe)
                       do k = 1,ke_list(l,tMe)
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                fieldy(i,j,k) = recv_buffer(pos)
                             end do
                          end do
                       end do 
                    end do
                 end if
              end if
           end do
           ind_y = ind_y-1
           if(ind_y .GT. 0) then
              rank_y = update_y%recv(ind_y)%pe - domain%pe
              if(rank_y .LE.0) rank_y = rank_y + nlist
           else
              rank_y = nlist+1
           endif
        endif
        if(cur_rank == rank_x) then
           do n = update_x%recv(ind_x)%count, 1, -1
              dir = update_x%recv(ind_x)%dir(n)
              if( recv(dir) ) then
                 tMe = update_x%recv(ind_x)%tileMe(n)
                 is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                 js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n) 
                 msgsize = (ie-is+1)*(je-js+1)*ke_sum
                 pos = buffer_pos - msgsize
                 buffer_pos = pos
                 if(update_x%recv(ind_x)%is_refined(n)) then
                    index = update_x%recv(ind_x)%index(n)
                    is1 = update_x%rSpec(tMe)%isNbr(index); ie1 = update_x%rSpec(tMe)%ieNbr(index)
                    js1 = update_x%rSpec(tMe)%jsNbr(index); je1 = update_x%rSpec(tMe)%jeNbr(index)
                    ni = ie1 - is1 + 1
                    nj = je1 - js1 + 1
                    total = ni*nj
                    start = (update_x%rSpec(tMe)%start(index)-1)*ke_max
                    if(start+total*ke_max>size(bufferx) ) call mpp_error(FATAL, &
                      "MPP_COMPLETE_DO_UPDATE_3D_V: size of bufferx is less than the size of the data to be filled.")

                    msgsize = ie - is + 1
                    do l=1, l_size  ! loop over number of fields
                       ptr_bufferx = b_addrsx(l, tMe)
                       if(l==1) start = (update_x%rSpec(tMe)%start(index)-1)*ke_list(l,tMe)
                       start1 = start + (js-js1)*ni + is - is1 
                       do k = 1, ke_list(l,tMe)
                          start2 = start1
                          do j = js, je
                             do i = start2+1, start2+msgsize
                                pos = pos + 1
                                bufferx(i) = recv_buffer(pos)
                             end do
                             start2 = start2 + ni
                          end do
                          start1 = start1 + total
                       end do
                    end do
                 else
                    do l=1, l_size  ! loop over number of fields
                       ptr_fieldx = f_addrsx(l, tMe)
                       do k = 1,ke_list(l,tMe)
                          do j = js, je
                             do i = is, ie
                                pos = pos + 1
                                fieldx(i,j,k) = recv_buffer(pos)
                             end do
                          end do
                       end do
                    end do
                 end if
              end if
           end do
           ind_x = ind_x-1
           if(ind_x .GT. 0) then
              rank_x = update_x%recv(ind_x)%pe - domain%pe 
              if(rank_x .LE.0) rank_x = rank_x + nlist
           else
              rank_x = nlist+1
           endif
        endif
     end select
     cur_rank = min(rank_x, rank_y)     
     call mpp_clock_end(unpk_clock_nonblock)
  end do

  ! ---northern boundary fold
  shift = 0
  if(domain%symmetry) shift = 1
  if( BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then
     j = domain%y(1)%global%end+shift
     if( domain%y(1)%data%begin.LE.j .AND. j.LE.domain%y(1)%data%end+shift )then !fold is within domain
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
           j  = domain%y(1)%global%end+shift
           is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
           if( .NOT. domain%symmetry ) is = is - 1
           do i = is ,ie, midpoint
              if( domain%x(1)%data%begin.LE.i .AND. i.LE. domain%x(1)%data%end+shift )then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 enddo
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off west edge, when update north or west direction
        j = domain%y(1)%global%end+shift 
        if ( recv(7) .OR. recv(5) ) then
           select case(gridtype)
           case(BGRID_NE)
              if(domain%symmetry) then
                 is = domain%x(1)%global%begin
              else
                 is = domain%x(1)%global%begin - 1
              end if
              if( is.GT.domain%x(1)%data%begin )then

                 if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-north BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)

                        do k = 1,ke_list(l,tMe)
                           do i = domain%x(1)%data%begin,is-1
                              fieldx(i,j,k) = fieldx(2*is-i,j,k)
                              fieldy(i,j,k) = fieldy(2*is-i,j,k)
                           end do
                        end do
                     end do
              end if
           case(CGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-north CGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do i = domain%x(1)%data%begin,is-1
                          fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if

        !off east edge
        is = domain%x(1)%global%end
        if(domain%x(1)%cyclic .AND. is.LT.domain%x(1)%data%end )then
           ie = domain%x(1)%data%end
           is = is + 1
           select case(gridtype)
           case(BGRID_NE)
              is = is + shift
              ie = ie + shift
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do i = is,ie
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,l_size
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do i = is, ie
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,SOUTH) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then      ! ---southern boundary fold
     ! NOTE: symmetry is assumed for fold-south boundary
     j = domain%y(1)%global%begin
     if( domain%y(1)%data%begin.LE.j .AND. j.LE.domain%y(1)%data%end+shift )then !fold is within domain
        midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           j  = domain%y(1)%global%begin
           is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
           do i = is ,ie, midpoint
              if( domain%x(1)%data%begin.LE.i .AND. i.LE. domain%x(1)%data%end+shift )then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off west edge, when update north or west direction
        j = domain%y(1)%global%begin
        if ( recv(3) .OR. recv(5) ) then
           select case(gridtype)
           case(BGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-south BGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do i = domain%x(1)%data%begin,is-1
                          fieldx(i,j,k) = fieldx(2*is-i,j,k)
                          fieldy(i,j,k) = fieldy(2*is-i,j,k)
                       end do
                    end do
                 end do
              end if
           case(CGRID_NE)
              is = domain%x(1)%global%begin
              if( is.GT.domain%x(1)%data%begin )then
                 if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO_UPDATE_V: folded-south CGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do i = domain%x(1)%data%begin,is-1
                          fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if

        !off east edge
        is = domain%x(1)%global%end
        if(domain%x(1)%cyclic .AND. is.LT.domain%x(1)%data%end )then
           ie = domain%x(1)%data%end
           is = is + 1
           select case(gridtype)
           case(BGRID_NE)
              is = is + shift
              ie = ie + shift
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do i = is,ie
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,l_size
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do i = is, ie
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,WEST) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then      ! ---eastern boundary fold
     ! NOTE: symmetry is assumed for fold-west boundary
     i = domain%x(1)%global%begin
     if( domain%x(1)%data%begin.LE.i .AND. i.LE.domain%x(1)%data%end+shift )then !fold is within domain
        midpoint = (domain%y(1)%global%begin+domain%y(1)%global%end-1+shift)/2
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           i  = domain%x(1)%global%begin
           js = domain%y(1)%global%begin; je = domain%y(1)%global%end+shift
           do j = js ,je, midpoint
              if( domain%y(1)%data%begin.LE.j .AND. j.LE. domain%y(1)%data%end+shift )then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off south edge, when update south or west direction
        i = domain%x(1)%global%begin
        if ( recv(3) .OR. recv(5) ) then
           select case(gridtype)
           case(BGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then

                 if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-west BGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do j = domain%y(1)%data%begin,js-1
                          fieldx(i,j,k) = fieldx(i,2*js-j,k)
                          fieldy(i,j,k) = fieldy(i,2*js-j,k)
                       end do
                    end do
                 end do
              end if
           case(CGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then
                 if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-west CGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do j = domain%y(1)%data%begin,js-1
                          fieldx(i,j,k) = fieldx(i, 2*js-j-1,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if

        !off north edge
        js = domain%y(1)%global%end
        if(domain%y(1)%cyclic .AND. js.LT.domain%y(1)%data%end )then
           je = domain%y(1)%data%end
           js = js + 1
           select case(gridtype)
           case(BGRID_NE)
              js = js + shift
              je = je + shift
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do j = js,je
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do j = js, je
                       fieldx(i,j,k) = -fieldx(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  else if( BTEST(domain%fold,EAST) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then      ! ---eastern boundary fold
     ! NOTE: symmetry is assumed for fold-west boundary
     i = domain%x(1)%global%end+shift
     if( domain%x(1)%data%begin.LE.i .AND. i.LE.domain%x(1)%data%end+shift )then !fold is within domain
        midpoint = (domain%y(1)%global%begin+domain%y(1)%global%end-1+shift)/2
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           i  = domain%x(1)%global%end+shift
           js = domain%y(1)%global%begin; je = domain%y(1)%global%end+shift
           do j = js ,je, midpoint
              if( domain%y(1)%data%begin.LE.j .AND. j.LE. domain%y(1)%data%end+shift )then
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       fieldx(i,j,k) = 0.
                       fieldy(i,j,k) = 0.
                    end do
                 end do
              end if
           end do
        endif

        ! the following code code block correct an error where the data in your halo coming from 
        ! other half may have the wrong sign
        !off south edge, when update south or west direction
        i = domain%x(1)%global%end+shift
        if ( recv(3) .OR. recv(1) ) then
           select case(gridtype)
           case(BGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then

                 if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-east BGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    ptr_fieldy = f_addrsy(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do j = domain%y(1)%data%begin,js-1
                          fieldx(i,j,k) = fieldx(i,2*js-j,k)
                          fieldy(i,j,k) = fieldy(i,2*js-j,k)
                       end do
                    end do
                 end do
              end if
           case(CGRID_NE)
              js = domain%y(1)%global%begin
              if( js.GT.domain%y(1)%data%begin )then
                 if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                      call mpp_error( FATAL, 'MPP_COMPLETE_DO__UPDATE_V: folded-east CGRID_NE west edge ubound error.' )
                 do l=1,l_size
                    ptr_fieldx = f_addrsx(l, 1)
                    do k = 1,ke_list(l,tMe)
                       do j = domain%y(1)%data%begin,js-1
                          fieldx(i,j,k) = fieldx(i, 2*js-j-1,k)
                       end do
                    end do
                 end do
              end if
           end select
        end if

        !off north edge
        js = domain%y(1)%global%end
        if(domain%y(1)%cyclic .AND. js.LT.domain%y(1)%data%end )then
           je = domain%y(1)%data%end
           js = js + 1
           select case(gridtype)
           case(BGRID_NE)
              js = js + shift
              je = je + shift
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 ptr_fieldy = f_addrsy(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do j = js,je
                       fieldx(i,j,k) = -fieldx(i,j,k)
                       fieldy(i,j,k) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           case(CGRID_NE)
              do l=1,l_size
                 ptr_fieldx = f_addrsx(l, 1)
                 do k = 1,ke_list(l,tMe)
                    do j = js, je
                       fieldx(i,j,k) = -fieldx(i,j,k)
                    end do
                 end do
              end do
           end select
        end if
     end if
  end if

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

  return

end subroutine MPP_COMPLETE_DO_UPDATE_3D_V_
