! -*-f90-*- 
    subroutine MPP_DO_UPDATE_AD_3D_V_(f_addrsx,f_addrsy, domainx, domainy, d_type, ke, gridtype, flags, name)
!updates data domain of 3D field whose computational domains have been computed
      integer(LONG_KIND), intent(in)         :: f_addrsx(:,:), f_addrsy(:,:)
      type(domain2d),     intent(in)         :: domainx, domainy
      integer,            intent(in)         :: ke
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface
      integer, intent(in),          optional :: gridtype
      integer, intent(in),          optional :: flags
      character(len=*), intent(in), optional :: name

      MPP_TYPE_ :: fieldx(domainx%x(1)%memory%begin:domainx%x(1)%memory%end, &
                          domainx%y(1)%memory%begin:domainx%y(1)%memory%end,ke)
      MPP_TYPE_ :: fieldy(domainy%x(1)%memory%begin:domainy%x(1)%memory%end, &
                          domainy%y(1)%memory%begin:domainy%y(1)%memory%end,ke)
      pointer(ptr_fieldx, fieldx)
      pointer(ptr_fieldy, fieldy)

      integer :: update_flags
      integer :: l_size, l, i, j, k, is, ie, js, je, m, n
      integer :: pos, nlist, list, msgsize
      integer :: to_pe, from_pe, midpoint
      integer :: tMe, dir
      logical :: send(8), recv(8)
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer(ptr,buffer )
      integer :: buffer_pos
      character(len=8) :: text
      character(len=64) :: field_name      
      type(overlapSpec), pointer :: overPtrx => NULL()
      type(overlapSpec), pointer :: overPtry => NULL()
      type(boundary),    pointer :: check_x  => NULL()
      type(boundary),    pointer :: check_y  => NULL()

      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) ) then 
          update_flags = flags
          ! The following test is so that SCALAR_PAIR can be used alone with the
          ! same default update pattern as without.
          if (BTEST(update_flags,SCALAR_BIT)) then
            if (.NOT.(BTEST(update_flags,WEST) .OR. BTEST(update_flags,EAST) &
                 .OR. BTEST(update_flags,NORTH) .OR. BTEST(update_flags,SOUTH))) &
              update_flags = update_flags + XUPDATE+YUPDATE   !default with SCALAR_PAIR
          end if 
      end if  

      if( BTEST(update_flags,NORTH) .AND. BTEST(domainx%fold,NORTH) .AND. BTEST(gridtype,SOUTH) ) &
           call mpp_error( FATAL, 'MPP_DO_UPDATE_V_AD: Incompatible grid offset and fold.' )

      recv(1) = BTEST(update_flags,EAST)
      recv(3) = BTEST(update_flags,SOUTH)
      recv(5) = BTEST(update_flags,WEST)
      recv(7) = BTEST(update_flags,NORTH)
      recv(2) = recv(1) .AND. recv(3)
      recv(4) = recv(3) .AND. recv(5)
      recv(6) = recv(5) .AND. recv(7)
      recv(8) = recv(7) .AND. recv(1)
      send    = recv

      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
      l_size = size(f_addrsx,1)
      nlist = size(domainx%list(:))
      ptr = LOC(mpp_domains_stack)

      !--- if debug_update_level is not NO_DEBUG, check the consistency on the bounds 
      !--- (domain is symmetry or folded north edge). North bound will be checked when north edge is folded.
      !--- when domain is symmetry, For data on T-cell, no check is needed; for data on E-cell, 
      !--- data on East and West boundary will be checked ; For data on N-cell, data on North and South 
      !--- boundary will be checked; For data on C-cell, data on West, East, South, North will be checked.
      !--- The check will be done in the following way: Western boundary data sent to Eastern boundary to check
      !--- and Southern boundary to check

      if(  debug_update_level .NE. NO_CHECK ) then      
         if(present(name)) then
            field_name = name
         else
            field_name = "un-named"
         end if

         !--- send the data
         do list = 0,nlist-1
            m = mod( domainx%pos+list, nlist )
            check_x => domainx%check%send(m)
            check_y => domainy%check%send(m)
            pos = buffer_pos
            do n = 1, check_x%count
               is = check_x%is(n); ie = check_x%ie(n)
               js = check_x%js(n); je = check_x%je(n)
               tMe = check_x%tileMe(n)
               select case( check_x%rotation(n) )
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
            do n = 1, check_y%count
               is = check_y%is(n); ie = check_y%ie(n)
               js = check_y%js(n); je = check_y%je(n)
               tMe = check_y%tileMe(n)
               select case( check_y%rotation(n) )
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

            msgsize = pos - buffer_pos
            if( msgsize.GT.0 )then
               to_pe = domainx%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_DO_UPDATE_V_AD: mpp_domains_stack overflow, ' // &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
               end if
               call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
               buffer_pos = pos
            end if
         end do ! end do list = 0,nlist-1

         !--- recv the data 
         do list = 0,nlist-1
            m = mod( domainx%pos+nlist-list, nlist )
            check_x=>domainx%check%recv(m)
            check_y=>domainy%check%recv(m)
            msgsize = 0
            do n = 1, check_x%count
               is = check_x%is(n); ie = check_x%ie(n)
               js = check_x%js(n); je = check_x%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            do n = 1, check_y%count
               is = check_y%is(n); ie = check_y%ie(n)
               js = check_y%js(n); je = check_y%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            msgsize = msgsize*ke*l_size

            if( msgsize.GT.0 )then
               from_pe = domainx%list(m)%pe
               mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
               if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                  write( text,'(i8)' )mpp_domains_stack_hwm
                  call mpp_error( FATAL, 'MPP_DO_UPDATE_V_AD: mpp_domains_stack overflow, '// &
                       'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
               end if
               call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
               buffer_pos = buffer_pos + msgsize
            end if
         end do

         !--- compare the data in reverse order
         CHECK_LOOP: do list = nlist-1,0,-1
            m = mod( domainx%pos+nlist-list, nlist )
            check_x=>domainx%check%recv(m)
            check_y=>domainy%check%recv(m)
            do n = check_y%count, 1, -1
               is = check_y%is(n); ie = check_y%ie(n)
               js = check_y%js(n); je = check_y%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               tMe = check_y%tileMe(n)
               do l=1,l_size  ! loop over number of fields
                  ptr_fieldx = f_addrsx(l, tMe)
                  ptr_fieldy = f_addrsy(l, tMe)
                  do k = 1,ke
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           if( fieldy(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_UPDATE_V_AD on pe = ", mpp_pe(), ": y component of vector ", &
                                   trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldy(i,j,k), &
                                   " does not equal to the value = ", buffer(pos), " on pe ", domainx%list(m)%pe
                              call mpp_error(debug_update_level, &
					"MPP_DO_UPDATE_V_AD: mismatch on the boundary for symmetry point")
                              exit CHECK_LOOP
                           end if
                        end do
                     end do
                  end do
               end do
            end do
            do n = check_x%count, 1, -1
               is = check_x%is(n); ie = check_x%ie(n)
               js = check_x%js(n); je = check_x%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               tMe = check_x%tileMe(n)
               do l=1,l_size  ! loop over number of fields
                  ptr_fieldx = f_addrsx(l, tMe)
                  ptr_fieldy = f_addrsy(l, tMe)
                  do k = 1,ke
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           if( fieldx(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_UPDATE_V_AD on pe = ", mpp_pe(), ": x-component of vector ", &
                                   trim(field_name), " at point (", i, ",", j, ",", k, ") = ", fieldx(i,j,k), &
                                   " does not equal to the value = ", buffer(pos), " on pe ", domainx%list(m)%pe
                              call mpp_error(debug_update_level, &
                                             "MPP_DO_UPDATE_V_AD: mismatch on the boundary for symmetry point")
                              exit CHECK_LOOP
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do CHECK_LOOP ! end do list = nlist-1,0,-1
      end if

      !--- send
      do list = 0,nlist-1
         m = mod( domainx%pos+list, nlist )
         overPtrx => domainx%update_recv(m)
         overPtry => domainy%update_recv(m)
         if( overPtrx%count == 0 .AND. overPtry%count == 0 )cycle
         call mpp_clock_begin(pack_clock)
         pos = buffer_pos
         select case( gridtype )
         case(BGRID_NE, BGRID_SW, AGRID)
            do n = 1, overPtrx%count
               dir = overPtrx%dir(n)
               if( recv(dir) ) then
                  tMe = overPtrx%tileMe(n)

                  is = overPtrx%is(n); ie = overPtrx%ie(n)
                  js = overPtrx%js(n); je = overPtrx%je(n)
                  select case( overPtrx%rotation(n) )
                  case(ZERO)
                     do l=1,l_size  ! loop over number of fields
                        ptr_fieldx = f_addrsx(l,tMe)
                        ptr_fieldy = f_addrsy(l,tMe)
                        do k = 1,ke
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
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l,tMe)
                           ptr_fieldy = f_addrsy(l,tMe)
                           do k = 1,ke
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
                           do k = 1,ke
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
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l,tMe)
                           ptr_fieldy = f_addrsy(l,tMe)
                           do k = 1,ke
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
                           do k = 1,ke
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
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l,tMe)
                           ptr_fieldy = f_addrsy(l,tMe)
                           do k = 1,ke
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
                           do k = 1,ke
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
                  end select ! select case( overPtrx%rotation(n) )
               end if ! if( send(dir) ) 
            end do ! do n = 1, overPtrx%count
         case(CGRID_NE, CGRID_SW)
            do n = 1, overPtrx%count
               dir = overPtrx%dir(n)
               if( recv(dir) ) then
                  tMe = overPtrx%tileMe(n)
                  is = overPtrx%is(n); ie = overPtrx%ie(n)
                  js = overPtrx%js(n); je = overPtrx%je(n)
                  select case( overPtrx%rotation(n) )
                  case(ZERO)
                     do l=1,l_size  ! loop over number of fields
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
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
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
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
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
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
                        do k = 1, ke
                           do i = ie, is, -1
                              do j = js, je
                                 pos = pos + 1
                                 buffer(pos) = fieldy(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  case(ONE_HUNDRED_EIGHTY)
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
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
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
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
            end do
            do n = 1, overPtry%count
               dir = overPtry%dir(n)
               if( recv(dir) ) then
                  tMe = overPtry%tileMe(n)
                  is = overPtry%is(n); ie = overPtry%ie(n)
                  js = overPtry%js(n); je = overPtry%je(n)
                  select case( overPtry%rotation(n) )
                  case(ZERO)
                     do l=1,l_size  ! loop over number of fields
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
                     do l=1,l_size  ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
                        do k = 1,ke
                           do i = is, ie
                              do j = je, js, -1
                                 pos = pos + 1
                                 buffer(pos) = fieldx(i,j,k)
                              end do
                           end do
                        end do
                     end do
                  case(NINETY)
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
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
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
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
                     if( BTEST(update_flags,SCALAR_BIT) ) then
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
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
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
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
            end do ! do n = 1, overPtry%count
         end select
         call mpp_clock_end(pack_clock)
         call mpp_clock_begin(send_clock)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = domainx%list(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE_V_AD: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe )
            buffer_pos = pos
         end if
         call mpp_clock_end(send_clock)
         !     write(stdout(),*) 'Update send checksum=',mpp_chksum(rbuffer(buffer_pos))
      end do ! do list = 0,nlist-1

!recv
      do list = 0,nlist-1
         m = mod( domainx%pos+nlist-list, nlist )
         overPtrx => domainx%update_send(m)
         overPtry => domainy%update_send(m)
         if( overPtrx%count == 0 .AND. overPtry%count == 0 )cycle
         call mpp_clock_begin(recv_clock)
         msgsize = 0
         select case(gridtype)
         case(BGRID_NE, BGRID_SW, AGRID)
            do n = 1, overPtrx%count
               dir = overPtrx%dir(n)
               if(send(dir)) then
                  tMe = overPtrx%tileMe(n)
                  is = overPtrx%is(n); ie = overPtrx%ie(n)
                  js = overPtrx%js(n); je = overPtrx%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
            end do
            msgsize = msgsize*2
         case(CGRID_NE, CGRID_SW)
            do n = 1, overPtrx%count
               dir = overPtrx%dir(n)
               if(send(dir)) then
                  is = overPtrx%is(n); ie = overPtrx%ie(n)
                  js = overPtrx%js(n); je = overPtrx%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
            end do
            do n = 1, overPtry%count
               dir = overPtry%dir(n)
               if(send(dir)) then
                  is = overPtry%is(n); ie = overPtry%ie(n)
                  js = overPtry%js(n); je = overPtry%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
            end do
         end select
         msgsize = msgsize*ke*l_size

         if( msgsize.GT.0 )then
            from_pe = domainx%list(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, buffer_pos+msgsize )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE_V_AD: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe )
            buffer_pos = buffer_pos + msgsize
         end if
         call mpp_clock_end(recv_clock)
      end do

      !unpack recv
      !unpack halos in reverse order

      do list = nlist-1,0,-1
         m = mod( domainx%pos+nlist-list, nlist )
         overPtrx => domainx%update_send(m)
         overPtry => domainy%update_send(m)
         if( overPtrx%count == 0 .AND. overPtry%count == 0 )cycle
         call mpp_clock_begin(unpk_clock)
         pos = buffer_pos
         select case ( gridtype )
         case(BGRID_NE, BGRID_SW, AGRID)
            do n = overPtrx%count, 1, -1
               dir = overPtrx%dir(n)
               if( send(dir) ) then
                  tMe = overPtrx%tileMe(n)
                  is = overPtrx%is(n); ie = overPtrx%ie(n)
                  js = overPtrx%js(n); je = overPtrx%je(n) 
                  msgsize = (ie-is+1)*(je-js+1)*ke*2*l_size
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  do l=1, l_size  ! loop over number of fields
                     ptr_fieldx = f_addrsx(l, tMe)
                     ptr_fieldy = f_addrsy(l, tMe)
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 2
                              fieldx(i,j,k) = buffer(pos-1)
                              fieldy(i,j,k) = buffer(pos)
                           end do
                        end do
                     end do
                  end do
               end if ! end if( send(dir) )
            end do  ! do dir=8,1,-1 
         case(CGRID_NE, CGRID_SW)
            do n = overPtry%count, 1, -1
               dir = overPtry%dir(n)
               if( send(dir) ) then
                  tMe = overPtry%tileMe(n)
                  is = overPtry%is(n); ie = overPtry%ie(n)
                  js = overPtry%js(n); je = overPtry%je(n)
                  msgsize = (ie-is+1)*(je-js+1)*ke*l_size
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  do l=1,l_size  ! loop over number of fields
                     ptr_fieldx = f_addrsx(l, tMe)
                     ptr_fieldy = f_addrsy(l, tMe)
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              fieldy(i,j,k) = buffer(pos)
                           end do
                        end do
                     end do
                  end do
               end if
            end do
            do n = overPtrx%count, 1, -1
               dir = overPtrx%dir(n)
               if( send(dir) ) then
                  tMe = overPtrx%tileMe(n)
                  is = overPtrx%is(n); ie = overPtrx%ie(n)
                  js = overPtrx%js(n); je = overPtrx%je(n) 
                  msgsize = (ie-is+1)*(je-js+1)*ke*l_size
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  do l=1,l_size  ! loop over number of fields
                     ptr_fieldx = f_addrsx(l, tMe)
                     ptr_fieldy = f_addrsy(l, tMe)
                     do k = 1,ke
                        do j = js, je
                           do i = is, ie
                              pos = pos + 1
                              fieldx(i,j,k) = buffer(pos)
                           end do
                        end do
                     end do
                  end do
               end if
            end do ! do n = 1, overPtrx%count 
         end select
         call mpp_clock_end(unpk_clock)
      end do

      ! ---northern boundary fold
      if( BTEST(domainx%fold,NORTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then

         j = domainy%y(1)%global%end
         if( domainy%y(1)%data%begin.LE.j .AND. j.LE.domainy%y(1)%data%end )then !fold is within domain
            midpoint = (domainy%x(1)%global%begin+domainy%x(1)%global%end-1)/2
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               j  = domainx%y(1)%global%end
               is = domainx%x(1)%global%begin; ie = domainx%x(1)%global%end
               if( .NOT. domainx%symmetry ) is = is - 1
               do i = is ,ie, midpoint
                  if( domainx%x(1)%data%begin.LE.i .AND. i.LE. domainx%x(1)%data%end )then
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)   
                        do k = 1,ke
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
            j = domainy%y(1)%global%end 
            if ( BTEST(update_flags,NORTH) .OR. BTEST(update_flags,WEST) ) then
               select case(gridtype)
               case(BGRID_NE)
                  is = domainx%x(1)%global%begin - 1
                  if(domainx%symmetry) then
                     is = domainx%x(1)%global%begin
                  else
                     is = domainx%x(1)%global%begin - 1
                  end if
                  if( is.GT.domainx%x(1)%data%begin )then

                     if( 2*is-domainx%x(1)%data%begin.GT.domainx%x(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V_AD: BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)   
                        do k = 1,ke
                           do i = domainx%x(1)%data%begin,is-1
                              fieldx(i,j,k) = fieldx(2*is-i,j,k)
                              fieldy(i,j,k) = fieldy(2*is-i,j,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  is = domainy%x(1)%global%begin
                  if( is.GT.domainy%x(1)%data%begin )then
                     if( 2*is-domainy%x(1)%data%begin-1.GT.domainy%x(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V_AD: CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldy = f_addrsy(l, 1)   
                        do k = 1,ke
                           do i = domainy%x(1)%data%begin,is-1
                              fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off east edge
            is = domainy%x(1)%global%end
            if(domainy%x(1)%cyclic .AND. is.LT.domainy%x(1)%data%end )then
               ie = domainy%x(1)%data%end
               is = is + 1
               select case(gridtype)
               case(BGRID_NE)
                  do l=1,l_size
                     ptr_fieldx = f_addrsx(l, 1)
                     ptr_fieldy = f_addrsy(l, 1)   
                     do k = 1,ke
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     ptr_fieldy = f_addrsy(l, 1)   
                     do k = 1,ke
                        do i = is, ie
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      end if

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)

      return

    end subroutine MPP_DO_UPDATE_AD_3D_V_
