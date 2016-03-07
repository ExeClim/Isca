! -*-f90-*- 
! this routine is used to retrieve scalar boundary data for symmetric domain. 

subroutine MPP_GET_BOUNDARY_2D_(field, domain, ebuffer, sbuffer, wbuffer, nbuffer, flags, &
                                position, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: field(:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffer(:), sbuffer(:), wbuffer(:), nbuffer(:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete
  
  MPP_TYPE_                             :: field3D(size(field,1),size(field,2),1)
  MPP_TYPE_, allocatable, dimension(:,:) :: ebuffer2D, sbuffer2D, wbuffer2D, nbuffer2D
  integer                               :: xcount, ycount

  pointer( ptr, field3D )
  ptr  = LOC(field)

  !--- We require wbuffer and ebuffer should coexist, sbuffer and nbuffer should coexist.
  xcount =0; ycount = 0
  if(present(ebuffer)) xcount = xcount + 1
  if(present(wbuffer)) xcount = xcount + 1
  if(present(sbuffer)) ycount = ycount + 1
  if(present(nbuffer)) ycount = ycount + 1
  if(xcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D: ebuffer and wbuffer should be paired together")
  if(ycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D: sbuffer and nbuffer should be paired together")
  if(xcount>0) then
     allocate(ebuffer2D(size(ebuffer(:)), 1), wbuffer2D(size(wbuffer(:)), 1) )
     ebuffer2D = RESHAPE( ebuffer, SHAPE(ebuffer2D) )
     wbuffer2D = RESHAPE( wbuffer, SHAPE(wbuffer2D) )
  end if

  if(ycount>0) then
     allocate(sbuffer2D(size(sbuffer(:)), 1), nbuffer2D(size(nbuffer(:)), 1) )
     sbuffer2D = RESHAPE( sbuffer, SHAPE(sbuffer2D) )
     nbuffer2D = RESHAPE( nbuffer, SHAPE(nbuffer2D) )
  end if

  if(xcount>0 .AND. ycount>0 ) then
     call mpp_get_boundary(field3D, domain, ebuffer=ebuffer2D, sbuffer=sbuffer2D, wbuffer=wbuffer2D, nbuffer=nbuffer2D, &
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  else if(xcount>0) then
     call mpp_get_boundary(field3D, domain, ebuffer=ebuffer2D, wbuffer=wbuffer2D, & 
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  else if(ycount>0) then
     call mpp_get_boundary(field3D, domain, sbuffer=sbuffer2D, nbuffer=nbuffer2D, &
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  end if

  if(xcount>0) then
     ebuffer = RESHAPE( ebuffer2D, SHAPE(ebuffer) )
     wbuffer = RESHAPE( wbuffer2D, SHAPE(wbuffer) )
     deallocate(ebuffer2D, wbuffer2D)
  end if
  if(ycount>0) then
     sbuffer = RESHAPE( sbuffer2D, SHAPE(sbuffer) )
     nbuffer = RESHAPE( nbuffer2D, SHAPE(nbuffer) )
     deallocate(sbuffer2D, nbuffer2D)
  end if

  return

end subroutine MPP_GET_BOUNDARY_2D_


!###############################################################################################
subroutine MPP_GET_BOUNDARY_3D_(field, domain, ebuffer, sbuffer, wbuffer, nbuffer, flags, &
                                position, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: field(:,:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffer(:,:), sbuffer(:,:), wbuffer(:,:), nbuffer(:,:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete

  integer                  :: update_flags, ntile
  logical                  :: need_ebuffer, need_sbuffer, need_wbuffer, need_nbuffer
  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrs=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrs=-9999
  integer, save    :: bsize(4)=0, isize=0, jsize=0, ksize=0, pos, list=0, l_size=0, upflags
  integer          :: buffer_size(4)
  integer          :: max_ntile, tile, update_position, ishift, jshift
  logical          :: do_update, is_complete, set_mismatch
  character(len=3) :: text
  MPP_TYPE_        :: d_type
  type(overlapSpec), pointer :: bound => NULL()

  ntile = size(domain%x(:))

  update_flags = XUPDATE+YUPDATE
  if(present(flags)) update_flags = flags
  update_position = CENTER
  if(present(position)) update_position = position

  !--- check if the suitable buffer are present
  need_ebuffer=.false.; need_sbuffer=.false.; need_wbuffer=.false.; need_nbuffer=.false.
  if( domain%symmetry .AND. PRESENT(position) ) then
     select case(position)
     case(CORNER)
       need_ebuffer=.true.; need_sbuffer=.true.; need_wbuffer=.true.; need_nbuffer=.true.
     case(NORTH)
       need_sbuffer=.true.; need_nbuffer=.true.
     case(EAST)
       need_ebuffer=.true.; need_wbuffer=.true.
     end select
   end if
   need_ebuffer = need_ebuffer .AND. BTEST(update_flags, EAST)
   need_sbuffer = need_sbuffer .AND. BTEST(update_flags, SOUTH)
   need_wbuffer = need_wbuffer .AND. BTEST(update_flags, WEST)
   need_nbuffer = need_nbuffer .AND. BTEST(update_flags, NORTH)
   if(need_ebuffer) then
      if(.NOT. PRESENT(ebuffer) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D: optional argument ebuffer should be presented')
   end if  
   if(need_sbuffer) then
      if(.NOT. PRESENT(sbuffer) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D: optional argument sbuffer should be presented')
   end if  
   if(need_wbuffer) then
      if(.NOT. PRESENT(wbuffer) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D: optional argument wbuffer should be presented')
   end if  
   if(need_nbuffer) then
      if(.NOT. PRESENT(nbuffer) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D: optional argument nbuffer should be presented')
   end if  

  tile = 1
  max_ntile = domain%max_ntile_pe
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete        
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  f_addrs(list, tile) = LOC(field)
  buffer_size = 0
  if(present(ebuffer)) then
     b_addrs(1, list, tile) = LOC(ebuffer)
     buffer_size(1) = size(ebuffer,1)
  end if

  if(present(sbuffer)) then
     b_addrs(2, list, tile) = LOC(sbuffer)
     buffer_size(2) = size(sbuffer,1)
  end if

  if(present(wbuffer)) then
     b_addrs(3, list, tile) = LOC(wbuffer)
     buffer_size(3) = size(wbuffer,1)
  end if

  if(present(nbuffer)) then
     b_addrs(4, list, tile) = LOC(nbuffer)
     buffer_size(4) = size(nbuffer,1)
  end if

  if(list == 1 .AND. tile == 1 )then
     isize=size(field,1); jsize=size(field,2); ksize = size(field,3); pos = update_position
     bsize = buffer_size; upflags = update_flags
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize .NE. size(field,1))
     set_mismatch = set_mismatch .OR. (jsize .NE. size(field,2))
     set_mismatch = set_mismatch .OR. (ksize .NE. size(field,3))
     set_mismatch = set_mismatch .OR. ANY( bsize .NE. buffer_size )
     set_mismatch = set_mismatch .OR. (update_position .NE. pos)
     set_mismatch = set_mismatch .OR. (upflags .NE. update_flags)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update )then 
     !--- only non-center data in symmetry domain will be retrieved.
     if(position == CENTER .OR. (.NOT. domain%symmetry) ) return 
     bound => search_bound_overlap(domain, update_position)
     call mpp_get_domain_shift(domain, ishift, jshift, update_position)
     if(size(field,1) .NE. domain%x(1)%memory%size+ishift .OR. size(field,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D: field is not on memory domain")
    if(ASSOCIATED(bound)) then
        call mpp_do_get_boundary(f_addrs(1:l_size,1:ntile), domain, bound, b_addrs(:,1:l_size,1:ntile), &
             bsize, ksize, d_type, update_flags)
     endif
     l_size=0; f_addrs=-9999; bsize=0; b_addrs=-9999; isize=0;  jsize=0;  ksize=0
  end if

end subroutine MPP_GET_BOUNDARY_3D_

!#########################################################################################
subroutine MPP_GET_BOUNDARY_4D_(field, domain, ebuffer, sbuffer, wbuffer, nbuffer, flags, &
                                position, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: field(:,:,:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffer(:,:,:), sbuffer(:,:,:), wbuffer(:,:,:), nbuffer(:,:,:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete
  
  MPP_TYPE_                              :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
  MPP_TYPE_, allocatable, dimension(:,:) :: ebuffer2D, sbuffer2D, wbuffer2D, nbuffer2D
  integer                                :: xcount, ycount
  pointer( ptr, field3D )
  ptr = LOC(field)

  !--- We require wbuffer and ebuffer should coexist, sbuffer and nbuffer should coexist.
  xcount = 0; ycount = 0
  if(present(ebuffer)) xcount = xcount + 1
  if(present(wbuffer)) xcount = xcount + 1
  if(present(sbuffer)) ycount = ycount + 1
  if(present(nbuffer)) ycount = ycount + 1
  if(xcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_4D: ebuffer and wbuffer should be paired together")
  if(ycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_4D: sbuffer and nbuffer should be paired together")
  if(xcount>0) then
     allocate(ebuffer2D(size(ebuffer,1), size(ebuffer,2)*size(ebuffer,3)) )
     allocate(ebuffer2D(size(wbuffer,1), size(wbuffer,2)*size(wbuffer,3)) )
     ebuffer2D = RESHAPE( ebuffer, SHAPE(ebuffer2D) )
     wbuffer2D = RESHAPE( wbuffer, SHAPE(wbuffer2D) )
  end if

  if(ycount>0) then
     allocate(sbuffer2D(size(sbuffer,1), size(sbuffer,2)*size(sbuffer,3)) )
     allocate(nbuffer2D(size(nbuffer,1), size(nbuffer,2)*size(nbuffer,3)) )
     sbuffer2D = RESHAPE( sbuffer, SHAPE(sbuffer2D) )
     nbuffer2D = RESHAPE( nbuffer, SHAPE(nbuffer2D) )
  end if

  if(xcount>0 .AND. ycount>0 ) then
     call mpp_get_boundary(field3D, domain, ebuffer=ebuffer2D, sbuffer=sbuffer2D, wbuffer=wbuffer2D, nbuffer=nbuffer2D, &
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  else if(xcount>0) then
     call mpp_get_boundary(field3D, domain, ebuffer=ebuffer2D, wbuffer=wbuffer2D, & 
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  else if(ycount>0) then
     call mpp_get_boundary(field3D, domain, sbuffer=sbuffer2D, nbuffer=nbuffer2D, &
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  end if

  if(xcount>0) then
     ebuffer = RESHAPE( ebuffer2D, SHAPE(ebuffer) )
     wbuffer = RESHAPE( wbuffer2D, SHAPE(wbuffer) )
     deallocate(ebuffer2D, wbuffer2D)
  end if
  if(ycount>0) then
     sbuffer = RESHAPE( sbuffer2D, SHAPE(sbuffer) )
     nbuffer = RESHAPE( nbuffer2D, SHAPE(nbuffer) )
     deallocate(sbuffer2D, nbuffer2D)
  end if

  return

end subroutine MPP_GET_BOUNDARY_4D_

!###############################################################################
subroutine MPP_GET_BOUNDARY_5D_(field, domain, ebuffer, sbuffer, wbuffer, nbuffer, flags, &
                                position, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: field(:,:,:,:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffer(:,:,:,:), sbuffer(:,:,:,:), wbuffer(:,:,:,:), nbuffer(:,:,:,:)
  integer,      intent(in), optional :: flags, position, tile_count
  logical,      intent(in), optional :: complete
  
  MPP_TYPE_                              :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
  MPP_TYPE_, allocatable, dimension(:,:) :: ebuffer2D, sbuffer2D, wbuffer2D, nbuffer2D
  integer                                :: xcount, ycount
  pointer( ptr, field3D )
  ptr = LOC(field)

  !--- We require wbuffer and ebuffer should coexist, sbuffer and nbuffer should coexist.
  xcount = 0; ycount = 0
  if(present(ebuffer)) xcount = xcount + 1
  if(present(wbuffer)) xcount = xcount + 1
  if(present(sbuffer)) ycount = ycount + 1
  if(present(nbuffer)) ycount = ycount + 1
  if(xcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_5D: ebuffer and wbuffer should be paired together")
  if(ycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_5D: sbuffer and nbuffer should be paired together")
  if(xcount>0) then
     allocate(ebuffer2D(size(ebuffer,1), size(ebuffer,2)*size(ebuffer,3)*size(ebuffer,4)) )
     allocate(ebuffer2D(size(wbuffer,1), size(wbuffer,2)*size(wbuffer,3)*size(wbuffer,4)) )
     ebuffer2D = RESHAPE( ebuffer, SHAPE(ebuffer2D) )
     wbuffer2D = RESHAPE( wbuffer, SHAPE(wbuffer2D) )
  end if

  if(ycount>0) then
     allocate(sbuffer2D(size(sbuffer,1), size(sbuffer,2)*size(sbuffer,3)*size(sbuffer,4)) )
     allocate(nbuffer2D(size(nbuffer,1), size(nbuffer,2)*size(nbuffer,3)*size(nbuffer,4)) )
     sbuffer2D = RESHAPE( sbuffer, SHAPE(sbuffer2D) )
     nbuffer2D = RESHAPE( nbuffer, SHAPE(nbuffer2D) )
  end if

  if(xcount>0 .AND. ycount>0 ) then
     call mpp_get_boundary(field3D, domain, ebuffer=ebuffer2D, sbuffer=sbuffer2D, wbuffer=wbuffer2D, nbuffer=nbuffer2D, &
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  else if(xcount>0) then
     call mpp_get_boundary(field3D, domain, ebuffer=ebuffer2D, wbuffer=wbuffer2D, & 
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  else if(ycount>0) then
     call mpp_get_boundary(field3D, domain, sbuffer=sbuffer2D, nbuffer=nbuffer2D, &
                           flags=flags, position=position, complete=complete, tile_count=tile_count)
  end if

  if(xcount>0) then
     ebuffer = RESHAPE( ebuffer2D, SHAPE(ebuffer) )
     wbuffer = RESHAPE( wbuffer2D, SHAPE(wbuffer) )
     deallocate(ebuffer2D, wbuffer2D)
  end if
  if(ycount>0) then
     sbuffer = RESHAPE( sbuffer2D, SHAPE(sbuffer) )
     nbuffer = RESHAPE( nbuffer2D, SHAPE(nbuffer) )
     deallocate(sbuffer2D, nbuffer2D)
  end if

  return

end subroutine MPP_GET_BOUNDARY_5D_


!####################################################################
! vector update
subroutine MPP_GET_BOUNDARY_2D_V_(fieldx, fieldy, domain, ebufferx, sbufferx, wbufferx, nbufferx, &
                                  ebuffery, sbuffery, wbuffery, nbuffery, flags, gridtype, &
                                  complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: fieldx(:,:), fieldy(:,:)
  MPP_TYPE_, intent(inout), optional :: ebufferx(:), sbufferx(:), wbufferx(:), nbufferx(:)
  MPP_TYPE_, intent(inout), optional :: ebuffery(:), sbuffery(:), wbuffery(:), nbuffery(:)
  integer,      intent(in), optional :: flags, gridtype, tile_count
  logical,      intent(in), optional :: complete
  
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)

  MPP_TYPE_, allocatable, dimension(:,:) :: ebufferx2D, sbufferx2D, wbufferx2D, nbufferx2D
  MPP_TYPE_, allocatable, dimension(:,:) :: ebuffery2D, sbuffery2D, wbuffery2D, nbuffery2D
  integer                                :: xxcount, xycount, yycount, yxcount

  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx  = LOC(fieldx)
  ptry  = LOC(fieldy)

  !--- We require wbuffex and ebufferx should coexist, sbufferx and nbufferx should coexist.
  !---            wbuffey and ebuffery should coexist, sbuffery and nbuffery should coexist.
  xxcount = 0; xycount = 0; yxcount = 0; yycount = 0
  if(present(ebufferx)) xxcount = xxcount + 1
  if(present(wbufferx)) xxcount = xxcount + 1
  if(present(ebuffery)) xycount = xycount + 1
  if(present(wbuffery)) xycount = xycount + 1
  if(present(sbufferx)) yxcount = yxcount + 1
  if(present(nbufferx)) yxcount = yxcount + 1
  if(present(sbuffery)) yycount = yycount + 1
  if(present(nbuffery)) yycount = yycount + 1

  if(xxcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: ebufferx and wbufferx should be paired together")
  if(xycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: ebuffery and wbuffery should be paired together")
  if(yxcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: sbufferx and nbufferx should be paired together")
  if(yycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: sbuffery and nbuffery should be paired together")

  if(xxcount>0) then
     allocate(ebufferx2D(size(ebufferx(:)), 1), wbufferx2D(size(wbufferx(:)), 1) )
     ebufferx2D = RESHAPE( ebufferx, SHAPE(ebufferx2D) )
     wbufferx2D = RESHAPE( wbufferx, SHAPE(wbufferx2D) )
  end if

  if(xycount>0) then
     allocate(ebuffery2D(size(ebuffery(:)), 1), wbuffery2D(size(wbuffery(:)), 1) )
     ebuffery2D = RESHAPE( ebuffery, SHAPE(ebuffery2D) )
     wbuffery2D = RESHAPE( wbuffery, SHAPE(wbuffery2D) )
  end if

  if(yxcount>0) then
     allocate(sbufferx2D(size(sbufferx(:)), 1), nbufferx2D(size(nbufferx(:)), 1) )
     sbufferx2D = RESHAPE( sbufferx, SHAPE(sbufferx2D) )
     nbufferx2D = RESHAPE( nbufferx, SHAPE(nbufferx2D) )
  end if

  if(yycount>0) then
     allocate(sbuffery2D(size(sbuffery(:)), 1), nbuffery2D(size(nbuffery(:)), 1) )
     sbuffery2D = RESHAPE( sbuffery, SHAPE(sbuffery2D) )
     nbuffery2D = RESHAPE( nbuffery, SHAPE(nbuffery2D) )
  end if

  !--- We are assuming flags will be always XUPDATE+YUPDATE, so there are three possible 
  if( xxcount>0 .AND. xycount>0 .AND. yxcount>0 .AND. yycount>0 ) then  ! BGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  ebufferx=ebufferx2D, sbufferx=sbufferx2D,              &
                           wbufferx=wbufferx2D, nbufferx=nbufferx2D, ebuffery=ebuffery2D, sbuffery=sbuffery2D, &
                           wbuffery=wbuffery2D, nbuffery=nbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  else if( xxcount>0 .AND. yycount>0 ) then  ! CGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  ebufferx=ebufferx2D, wbufferx=wbufferx2D,              &
                           sbuffery=sbuffery2D, nbuffery=nbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  else if( xycount>0 .AND. yxcount>0 ) then  ! DGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  sbufferx=sbufferx2D, nbufferx=nbufferx2D,              &
                           ebuffery=ebuffery2D, wbuffery=wbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  end if

  if(xxcount>0) then
     ebufferx = RESHAPE( ebufferx2D, SHAPE(ebufferx) )
     wbufferx = RESHAPE( wbufferx2D, SHAPE(wbufferx) )
     deallocate(ebufferx2D, wbufferx2D)
  end if
  if(xycount>0) then
     ebuffery = RESHAPE( ebuffery2D, SHAPE(ebuffery) )
     wbuffery = RESHAPE( wbuffery2D, SHAPE(wbuffery) )
     deallocate(ebuffery2D, wbuffery2D)
  end if

  if(yxcount>0) then
     sbufferx = RESHAPE( sbufferx2D, SHAPE(sbufferx) )
     nbufferx = RESHAPE( nbufferx2D, SHAPE(nbufferx) )
     deallocate(sbufferx2D, nbufferx2D)
  end if
  if(yycount>0) then
     sbuffery = RESHAPE( sbuffery2D, SHAPE(sbuffery) )
     nbuffery = RESHAPE( nbuffery2D, SHAPE(nbuffery) )
     deallocate(sbuffery2D, nbuffery2D)
  end if

  return

end subroutine MPP_GET_BOUNDARY_2D_V_


!###############################################################################################
subroutine MPP_GET_BOUNDARY_3D_V_(fieldx, fieldy, domain, ebufferx, sbufferx, wbufferx, nbufferx, &
                                  ebuffery, sbuffery, wbuffery, nbuffery, flags, gridtype, &
                                  complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: fieldx(domain%x(1)%memory%begin:,domain%y(1)%memory%begin:,:)
  MPP_TYPE_,            intent(in)   :: fieldy(domain%x(1)%memory%begin:,domain%y(1)%memory%begin:,:)
  MPP_TYPE_, intent(inout), optional :: ebufferx(:,:), sbufferx(:,:), wbufferx(:,:), nbufferx(:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffery(:,:), sbuffery(:,:), wbuffery(:,:), nbuffery(:,:)
  integer,      intent(in), optional :: flags, gridtype, tile_count
  logical,      intent(in), optional :: complete

  integer                 :: ntile, update_flags
  logical                 :: need_ebufferx, need_sbufferx, need_wbufferx, need_nbufferx
  logical                 :: need_ebuffery, need_sbuffery, need_wbuffery, need_nbuffery

  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrsx=-9999
  integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),  save :: f_addrsy=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrsx=-9999
  integer(LONG_KIND),dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_addrsy=-9999
  integer, save    :: bsizex(4)=0, bsizey(4)=0, isize(2)=0, jsize(2)=0, ksize=0, l_size=0, list=0
  integer, save    :: offset_type, upflags
  integer          :: bufferx_size(4), buffery_size(4)
  integer          :: max_ntile, tile, grid_offset_type
  logical          :: do_update, is_complete, set_mismatch
  character(len=3) :: text
  MPP_TYPE_        :: d_type
  type(overlapSpec), pointer :: boundx=>NULL()
  type(overlapSpec), pointer :: boundy=>NULL()
  integer                     :: position_x, position_y, ishift, jshift

  ntile = size(domain%x(:))
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

  !--- check if the suitable buffer are present
  need_ebufferx=.FALSE.; need_sbufferx=.FALSE.
  need_wbufferx=.FALSE.; need_nbufferx=.FALSE.
  need_ebuffery=.FALSE.; need_sbuffery=.FALSE.
  need_wbuffery=.FALSE.; need_nbuffery=.FALSE.
  if( domain%symmetry .AND. PRESENT(gridtype) ) then
     select case(gridtype)
     case(BGRID_NE, BGRID_SW)
       need_ebufferx=.true.; need_sbufferx=.true.; need_wbufferx=.true.; need_nbufferx=.true.
       need_ebuffery=.true.; need_sbuffery=.true.; need_wbuffery=.true.; need_nbuffery=.true.
     case(CGRID_NE, CGRID_SW)
       need_ebufferx=.true.; need_wbufferx=.true.; need_sbuffery=.true.; need_nbuffery=.true.
     case(DGRID_NE, DGRID_SW)
       need_ebuffery=.true.; need_wbuffery=.true.; need_sbufferx=.true.; need_nbufferx=.true.
     end select
   end if

   need_ebufferx = need_ebufferx .AND. BTEST(update_flags, EAST)
   need_sbufferx = need_sbufferx .AND. BTEST(update_flags, SOUTH)
   need_wbufferx = need_wbufferx .AND. BTEST(update_flags, WEST)
   need_nbufferx = need_nbufferx .AND. BTEST(update_flags, NORTH)
   need_ebuffery = need_ebuffery .AND. BTEST(update_flags, EAST)
   need_sbuffery = need_sbuffery .AND. BTEST(update_flags, SOUTH)
   need_wbuffery = need_wbuffery .AND. BTEST(update_flags, WEST)
   need_nbuffery = need_nbuffery .AND. BTEST(update_flags, NORTH)
   if(need_ebufferx) then
      if(.NOT. PRESENT(ebufferx) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D_V: optional argument ebufferx should be presented')
   end if  
   if(need_sbufferx) then
      if(.NOT. PRESENT(sbufferx) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D_V: optional argument sbufferx should be presented')
   end if  
   if(need_wbufferx) then
      if(.NOT. PRESENT(wbufferx) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D_V: optional argument wbufferx should be presented')
   end if  
   if(need_nbufferx) then
      if(.NOT. PRESENT(nbufferx) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D_V: optional argument nbufferx should be presented')
   end if  
   if(need_ebuffery) then
      if(.NOT. PRESENT(ebuffery) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D_V: optional argument ebuffery should be presented')
   end if  
   if(need_sbuffery) then
      if(.NOT. PRESENT(sbuffery) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D_V: optional argument sbuffery should be presented')
   end if  
   if(need_wbuffery) then
      if(.NOT. PRESENT(wbuffery) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D_V: optional argument wbuffery should be presented')
   end if  
   if(need_nbuffery) then
      if(.NOT. PRESENT(nbuffery) ) call mpp_error( FATAL,'MPP_GET_BOUNDARY_3D_V: optional argument nbuffery should be presented')
   end if  

  tile = 1
  max_ntile = domain%max_ntile_pe
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete        
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  f_addrsx(list, tile) = LOC(fieldx)
  f_addrsy(list, tile) = LOC(fieldy)

  bufferx_size = 0; buffery_size = 0     
  if(present(ebufferx)) then
     b_addrsx(1, list, tile) = LOC(ebufferx)
     bufferx_size(1) = size(ebufferx,1)
  end if

  if(present(sbufferx)) then
     b_addrsx(2, list, tile) = LOC(sbufferx)
     bufferx_size(2) = size(sbufferx,1)
  end if

  if(present(wbufferx)) then
     b_addrsx(3, list, tile) = LOC(wbufferx)
     bufferx_size(3) = size(wbufferx,1)
  end if

  if(present(nbufferx)) then
     b_addrsx(4, list, tile) = LOC(nbufferx)
     bufferx_size(4) = size(nbufferx,1)
  end if

  if(present(ebuffery)) then
     b_addrsy(1, list, tile) = LOC(ebuffery)
     buffery_size(1) = size(ebuffery,1)
  end if

  if(present(sbuffery)) then
     b_addrsy(2, list, tile) = LOC(sbuffery)
     buffery_size(2) = size(sbuffery,1)
  end if

  if(present(wbuffery)) then
     b_addrsy(3, list, tile) = LOC(wbuffery)
     buffery_size(3) = size(wbuffery,1)
  end if

  if(present(nbuffery)) then
     b_addrsy(4, list, tile) = LOC(nbuffery)
     buffery_size(4) = size(nbuffery,1)
  end if

  grid_offset_type = AGRID
  if(present(gridtype)) grid_offset_type = gridtype
  if(list == 1 .AND. tile == 1 )then
     isize(1)=size(fieldx,1); jsize(1)=size(fieldx,2); isize(2)=size(fieldy,1); jsize(2)=size(fieldy,2)
     ksize = size(fieldx,3); offset_type = grid_offset_type
     bsizex = bufferx_size; bsizey = buffery_size; upflags = update_flags
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize(1) .NE. size(fieldx,1))
     set_mismatch = set_mismatch .OR. (jsize(1) .NE. size(fieldx,2))
     set_mismatch = set_mismatch .OR. (ksize    .NE. size(fieldx,3))
     set_mismatch = set_mismatch .OR. (isize(2) .NE. size(fieldy,1))
     set_mismatch = set_mismatch .OR. (jsize(2) .NE. size(fieldy,2))
     set_mismatch = set_mismatch .OR. (ksize    .NE. size(fieldy,3))
     set_mismatch = set_mismatch .OR. ANY( bsizex .NE. bufferx_size )
     set_mismatch = set_mismatch .OR. ANY( bsizey .NE. buffery_size )
     set_mismatch = set_mismatch .OR. (offset_type .NE. grid_offset_type)
     set_mismatch = set_mismatch .OR. (upflags .NE. update_flags)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_GET_BOUNDARY_3D_V: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update )then
     select case(grid_offset_type)
     case (AGRID)
        position_x = CENTER
        position_y = CENTER
     case (BGRID_NE, BGRID_SW)
        position_x = CORNER
        position_y = CORNER
     case (CGRID_NE, CGRID_SW)
        position_x = EAST
        position_y = NORTH
     case (DGRID_NE, DGRID_SW)
        position_x = NORTH
        position_y = EAST
     case default
        call mpp_error(FATAL, "mpp_get_boundary.h: invalid value of grid_offset_type")
     end select

     boundx => search_bound_overlap(domain, position_x)
     boundy => search_bound_overlap(domain, position_y)  

     call mpp_get_domain_shift(domain, ishift, jshift, position_x)
     if(size(fieldx,1) .NE. domain%x(1)%memory%size+ishift .OR. size(fieldx,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D_V: fieldx is not on memory domain")
     call mpp_get_domain_shift(domain, ishift, jshift, position_y)
     if(size(fieldy,1) .NE. domain%x(1)%memory%size+ishift .OR. size(fieldy,2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "MPP_GET_BOUNDARY_3D_V: fieldy is not on memory domain")
     if(ASSOCIATED(boundx) ) then
        call mpp_do_get_boundary(f_addrsx(1:l_size,1:ntile), f_addrsy(1:l_size,1:ntile), domain, boundx, boundy, &
             b_addrsx(:,1:l_size,1:ntile), b_addrsy(:,1:l_size,1:ntile), bsizex, &
             bsizey, ksize, d_type, update_flags)
     endif
     l_size=0; f_addrsx=-9999; f_addrsy=-9999; bsizex=0; bsizey=0; 
     b_addrsx=-9999; b_addrsy=-9999; isize=0;  jsize=0;  ksize=0
  end if

end subroutine MPP_GET_BOUNDARY_3D_V_

!##########################################################################
subroutine MPP_GET_BOUNDARY_4D_V_(fieldx, fieldy, domain, ebufferx, sbufferx, wbufferx, nbufferx, &
                                  ebuffery, sbuffery, wbuffery, nbuffery, flags, gridtype, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  MPP_TYPE_, intent(inout), optional :: ebufferx(:,:,:), sbufferx(:,:,:), wbufferx(:,:,:), nbufferx(:,:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffery(:,:,:), sbuffery(:,:,:), wbuffery(:,:,:), nbuffery(:,:,:)
  integer,      intent(in), optional :: flags, gridtype, tile_count
  logical,      intent(in), optional :: complete
  
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))
  MPP_TYPE_, allocatable, dimension(:,:) :: ebufferx2D, sbufferx2D, wbufferx2D, nbufferx2D
  MPP_TYPE_, allocatable, dimension(:,:) :: ebuffery2D, sbuffery2D, wbuffery2D, nbuffery2D
  integer                                :: xxcount, xycount, yycount, yxcount

  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx  = LOC(fieldx)
  ptry  = LOC(fieldy)

  !--- We require wbuffex and ebufferx should coexist, sbufferx and nbufferx should coexist.
  !---            wbuffey and ebuffery should coexist, sbuffery and nbuffery should coexist.
  xxcount = 0; xycount = 0; yxcount = 0; yycount = 0
  if(present(ebufferx)) xxcount = xxcount + 1
  if(present(wbufferx)) xxcount = xxcount + 1
  if(present(ebuffery)) xycount = xycount + 1
  if(present(wbuffery)) xycount = xycount + 1
  if(present(sbufferx)) yxcount = yxcount + 1
  if(present(nbufferx)) yxcount = yxcount + 1
  if(present(sbuffery)) yycount = yycount + 1
  if(present(nbuffery)) yycount = yycount + 1

  if(xxcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: ebufferx and wbufferx should be paired together")
  if(xycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: ebuffery and wbuffery should be paired together")
  if(yxcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: sbufferx and nbufferx should be paired together")
  if(yycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: sbuffery and nbuffery should be paired together")

  if(xxcount>0) then
     allocate(ebufferx2D(size(ebufferx,1), size(ebufferx,2)*size(ebufferx,3)))
     allocate(wbufferx2D(size(wbufferx,1), size(wbufferx,2)*size(wbufferx,3)))
     ebufferx2D = RESHAPE( ebufferx, SHAPE(ebufferx2D) )
     wbufferx2D = RESHAPE( wbufferx, SHAPE(wbufferx2D) )
  end if

  if(xycount>0) then
     allocate(ebuffery2D(size(ebuffery,1), size(ebuffery,2)*size(ebuffery,3)))
     allocate(wbuffery2D(size(wbuffery,1), size(wbuffery,2)*size(wbuffery,3)))
     ebuffery2D = RESHAPE( ebuffery, SHAPE(ebuffery2D) )
     wbuffery2D = RESHAPE( wbuffery, SHAPE(wbuffery2D) )
  end if

  if(yxcount>0) then
     allocate(sbufferx2D(size(sbufferx,1), size(sbufferx,2)*size(sbufferx,3)))
     allocate(nbufferx2D(size(nbufferx,1), size(nbufferx,2)*size(nbufferx,3)))
     sbufferx2D = RESHAPE( sbufferx, SHAPE(sbufferx2D) )
     nbufferx2D = RESHAPE( nbufferx, SHAPE(nbufferx2D) )
  end if

  if(yycount>0) then
     allocate(sbuffery2D(size(sbuffery,1), size(sbuffery,2)*size(sbuffery,3)))
     allocate(nbuffery2D(size(nbuffery,1), size(nbuffery,2)*size(nbuffery,3)))
     sbuffery2D = RESHAPE( sbuffery, SHAPE(sbuffery2D) )
     nbuffery2D = RESHAPE( nbuffery, SHAPE(nbuffery2D) )
  end if

  !--- We are assuming flags will be always XUPDATE+YUPDATE, so there are three possible 
  if( xxcount>0 .AND. xycount>0 .AND. yxcount>0 .AND. yycount>0 ) then  ! BGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  ebufferx=ebufferx2D, sbufferx=sbufferx2D,              &
                           wbufferx=wbufferx2D, nbufferx=nbufferx2D, ebuffery=ebuffery2D, sbuffery=sbuffery2D, &
                           wbuffery=wbuffery2D, nbuffery=nbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  else if( xxcount>0 .AND. yycount>0 ) then  ! CGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  ebufferx=ebufferx2D, wbufferx=wbufferx2D,              &
                           sbuffery=sbuffery2D, nbuffery=nbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  else if( xycount>0 .AND. yxcount>0 ) then  ! DGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  sbufferx=sbufferx2D, nbufferx=nbufferx2D,              &
                           ebuffery=ebuffery2D, wbuffery=wbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  end if

  if(xxcount>0) then
     ebufferx = RESHAPE( ebufferx2D, SHAPE(ebufferx) )
     wbufferx = RESHAPE( wbufferx2D, SHAPE(wbufferx) )
     deallocate(ebufferx2D, wbufferx2D)
  end if
  if(xycount>0) then
     ebuffery = RESHAPE( ebuffery2D, SHAPE(ebuffery) )
     wbuffery = RESHAPE( wbuffery2D, SHAPE(wbuffery) )
     deallocate(ebuffery2D, wbuffery2D)
  end if

  if(yxcount>0) then
     sbufferx = RESHAPE( sbufferx2D, SHAPE(sbufferx) )
     nbufferx = RESHAPE( nbufferx2D, SHAPE(nbufferx) )
     deallocate(sbufferx2D, nbufferx2D)
  end if
  if(yycount>0) then
     sbuffery = RESHAPE( sbuffery2D, SHAPE(sbuffery) )
     nbuffery = RESHAPE( nbuffery2D, SHAPE(nbuffery) )
     deallocate(sbufferx2D, nbuffery2D)
  end if

  return

end subroutine MPP_GET_BOUNDARY_4D_V_

!###############################################################################
subroutine MPP_GET_BOUNDARY_5D_V_(fieldx, fieldy, domain, ebufferx, sbufferx, wbufferx, nbufferx, &
                                  ebuffery, sbuffery, wbuffery, nbuffery, flags, gridtype, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  MPP_TYPE_,            intent(in)   :: fieldx(:,:,:,:,:), fieldy(:,:,:,:,:)
  MPP_TYPE_, intent(inout), optional :: ebufferx(:,:,:,:), sbufferx(:,:,:,:), wbufferx(:,:,:,:), nbufferx(:,:,:,:)
  MPP_TYPE_, intent(inout), optional :: ebuffery(:,:,:,:), sbuffery(:,:,:,:), wbuffery(:,:,:,:), nbuffery(:,:,:,:)
  integer,      intent(in), optional :: flags, gridtype, tile_count
  logical,      intent(in), optional :: complete
  
  MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
  MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))
  MPP_TYPE_, allocatable, dimension(:,:) :: ebufferx2D, sbufferx2D, wbufferx2D, nbufferx2D
  MPP_TYPE_, allocatable, dimension(:,:) :: ebuffery2D, sbuffery2D, wbuffery2D, nbuffery2D
  integer                                :: xxcount, xycount, yycount, yxcount

  pointer( ptrx, field3Dx )
  pointer( ptry, field3Dy )
  ptrx  = LOC(fieldx)
  ptry  = LOC(fieldy)

  !--- We require wbuffex and ebufferx should coexist, sbufferx and nbufferx should coexist.
  !---            wbuffey and ebuffery should coexist, sbuffery and nbuffery should coexist.
  xxcount = 0; xycount = 0; yxcount = 0; yycount = 0
  if(present(ebufferx)) xxcount = xxcount + 1
  if(present(wbufferx)) xxcount = xxcount + 1
  if(present(ebuffery)) xycount = xycount + 1
  if(present(wbuffery)) xycount = xycount + 1
  if(present(sbufferx)) yxcount = yxcount + 1
  if(present(nbufferx)) yxcount = yxcount + 1
  if(present(sbuffery)) yycount = yycount + 1
  if(present(nbuffery)) yycount = yycount + 1

  if(xxcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: ebufferx and wbufferx should be paired together")
  if(xycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: ebuffery and wbuffery should be paired together")
  if(yxcount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: sbufferx and nbufferx should be paired together")
  if(yycount==1) call mpp_error(FATAL, "MPP_GET_BOUNDARY_2D_V: sbuffery and nbuffery should be paired together")

  if(xxcount>0) then
     allocate(ebufferx2D(size(ebufferx,1), size(ebufferx,2)*size(ebufferx,3)*size(ebufferx,4)))
     allocate(wbufferx2D(size(wbufferx,1), size(wbufferx,2)*size(wbufferx,3)*size(wbufferx,4)))
     ebufferx2D = RESHAPE( ebufferx, SHAPE(ebufferx2D) )
     wbufferx2D = RESHAPE( wbufferx, SHAPE(wbufferx2D) )
  end if

  if(xycount>0) then
     allocate(ebuffery2D(size(ebuffery,1), size(ebuffery,2)*size(ebuffery,3)*size(ebuffery,4)))
     allocate(wbuffery2D(size(wbuffery,1), size(wbuffery,2)*size(wbuffery,3)*size(wbuffery,4)))
     ebuffery2D = RESHAPE( ebuffery, SHAPE(ebuffery2D) )
     wbuffery2D = RESHAPE( wbuffery, SHAPE(wbuffery2D) )
  end if

  if(yxcount>0) then
     allocate(sbufferx2D(size(sbufferx,1), size(sbufferx,2)*size(sbufferx,3)*size(sbufferx,4)))
     allocate(nbufferx2D(size(nbufferx,1), size(nbufferx,2)*size(nbufferx,3)*size(nbufferx,4)))
     sbufferx2D = RESHAPE( sbufferx, SHAPE(sbufferx2D) )
     nbufferx2D = RESHAPE( nbufferx, SHAPE(nbufferx2D) )
  end if

  if(yycount>0) then
     allocate(sbuffery2D(size(sbuffery,1), size(sbuffery,2)*size(sbuffery,3)*size(sbuffery,4)))
     allocate(nbuffery2D(size(nbuffery,1), size(nbuffery,2)*size(nbuffery,3)*size(nbuffery,4)))
     sbuffery2D = RESHAPE( sbuffery, SHAPE(sbuffery2D) )
     nbuffery2D = RESHAPE( nbuffery, SHAPE(nbuffery2D) )
  end if

  !--- We are assuming flags will be always XUPDATE+YUPDATE, so there are three possible 
  if( xxcount>0 .AND. xycount>0 .AND. yxcount>0 .AND. yycount>0 ) then  ! BGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  ebufferx=ebufferx2D, sbufferx=sbufferx2D,              &
                           wbufferx=wbufferx2D, nbufferx=nbufferx2D, ebuffery=ebuffery2D, sbuffery=sbuffery2D, &
                           wbuffery=wbuffery2D, nbuffery=nbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  else if( xxcount>0 .AND. yycount>0 ) then  ! CGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  ebufferx=ebufferx2D, wbufferx=wbufferx2D,              &
                           sbuffery=sbuffery2D, nbuffery=nbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  else if( xycount>0 .AND. yxcount>0 ) then  ! DGRID
     call mpp_get_boundary(field3Dx, field3Dy, domain,  sbufferx=sbufferx2D, nbufferx=nbufferx2D,              &
                           ebuffery=ebuffery2D, wbuffery=wbuffery2D, flags=flags, gridtype=gridtype,           &
                           complete=complete, tile_count=tile_count)
  end if

  if(xxcount>0) then
     ebufferx = RESHAPE( ebufferx2D, SHAPE(ebufferx) )
     wbufferx = RESHAPE( wbufferx2D, SHAPE(wbufferx) )
     deallocate(ebufferx2D, wbufferx2D)
  end if
  if(xycount>0) then
     ebuffery = RESHAPE( ebuffery2D, SHAPE(ebuffery) )
     wbuffery = RESHAPE( wbuffery2D, SHAPE(wbuffery) )
     deallocate(ebuffery2D, wbuffery2D)
  end if

  if(yxcount>0) then
     sbufferx = RESHAPE( sbufferx2D, SHAPE(sbufferx) )
     nbufferx = RESHAPE( nbufferx2D, SHAPE(nbufferx) )
     deallocate(sbufferx2D, nbufferx2D)
  end if
  if(yycount>0) then
     sbuffery = RESHAPE( sbuffery2D, SHAPE(sbuffery) )
     nbuffery = RESHAPE( nbuffery2D, SHAPE(nbuffery) )
     deallocate(sbufferx2D, nbuffery2D)
  end if

  return

end subroutine MPP_GET_BOUNDARY_5D_V_
