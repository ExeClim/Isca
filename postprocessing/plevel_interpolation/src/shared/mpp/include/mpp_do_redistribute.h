    subroutine MPP_DO_REDISTRIBUTE_3D_( f_in, f_out, d_comm, d_type )
      integer(LONG_KIND), intent(in)         :: f_in(:), f_out(:)
      type(DomainCommunicator2D), intent(in) :: d_comm
      MPP_TYPE_, intent(in)                  :: d_type
      MPP_TYPE_ :: field_in(d_comm%domain_in%x(1)%data%begin:d_comm%domain_in%x(1)%data%end, &
                            d_comm%domain_in%y(1)%data%begin:d_comm%domain_in%y(1)%data%end,d_comm%ke)
      pointer( ptr_field_in, field_in)
      MPP_TYPE_ :: field_out(d_comm%domain_out%x(1)%data%begin:d_comm%domain_out%x(1)%data%end, &
                             d_comm%domain_out%y(1)%data%begin:d_comm%domain_out%y(1)%data%end,d_comm%ke)
      pointer( ptr_field_out, field_out)
      type(domain2D), pointer :: domain_in, domain_out
      integer :: i, j, k, l, n, l_size
      integer :: is, ie, js, je
      integer :: ke
      integer :: list, pos, msgsize
      integer :: to_pe, from_pe
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer( ptr, buffer )
      integer :: buffer_pos, wordlen

!fix ke
      l_size = size(f_out(:))  ! equal to size(f_in(:))
      ke = d_comm%ke
      domain_in =>d_comm%domain_in; domain_out =>d_comm%domain_out

      buffer_pos = 0
      ptr = LOC(mpp_domains_stack)
      wordlen = size(TRANSFER(buffer(1),mpp_domains_stack))
!send
      n = d_comm%Slist_size
      do list = 0,n-1
         if( .NOT. d_comm%S_do_buf(list) )cycle
         to_pe = d_comm%cto_pe(list)
         is=d_comm%sendis(1,list); ie=d_comm%sendie(1,list)
         js=d_comm%sendjs(1,list); je=d_comm%sendje(1,list)
         pos = buffer_pos
         do l=1,l_size  ! loop over number of fields
           ptr_field_in = f_in(l)
           do k = 1,ke
              do j = js,je
                 do i = is,ie
                    pos = pos+1
                    buffer(pos) = field_in(i,j,k)
                 end do
              end do
           end do
         end do
         if( debug )write( stderr(),* )'PE', pe, ' to PE ', to_pe, 'is,ie,js,je=', is, ie, js, je
         msgsize = pos - buffer_pos
         call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
         buffer_pos = pos
      end do
!recv
      n = d_comm%Rlist_size
      do list = 0,n-1
         if( .NOT. d_comm%R_do_buf(list) )cycle
         from_pe = d_comm%cfrom_pe(list)
         is=d_comm%recvis(1,list); ie=d_comm%recvie(1,list)
         js=d_comm%recvjs(1,list); je=d_comm%recvje(1,list)
         msgsize = d_comm%R_msize(list)*l_size
         if( debug )write( stderr(),* )'PE', pe, ' from PE ', from_pe, 'is,ie,js,je=', is, ie, js, je
         call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
         pos = buffer_pos
         do l=1,l_size  ! loop over number of in/out fields
           ptr_field_out = f_out(l)
           do k = 1,ke
              do j = js,je
                 do i = is,ie
                    pos = pos+1
                    field_out(i,j,k) = buffer(pos)
                 end do
              end do
           end do
         end do
         buffer_pos = pos
      end do
      call mpp_sync_self()
    end subroutine MPP_DO_REDISTRIBUTE_3D_
