! -*-f90-*-
    subroutine MPP_UPDATE_DOMAINS_AD_2D_( field, domain, flags, complete, position, &
                                          whalo, ehalo, shalo, nhalo, name, tile_count)
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_,        intent(inout)        :: field(:,:)
      type(domain2D),   intent(inout)        :: domain  ! Must be definable in mpp_update_init_comm
      integer,          intent(in), optional :: flags
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains_ad( field3D, domain, flags, complete, position, &
                                  whalo, ehalo, shalo, nhalo, name, tile_count)

      return
    end subroutine MPP_UPDATE_DOMAINS_AD_2D_

    subroutine MPP_UPDATE_DOMAINS_AD_3D_( field, domain, flags, complete, position, &
                                          whalo, ehalo, shalo, nhalo, name, tile_count)
!updates data domain of 3D field whose computational domains have been computed
      MPP_TYPE_,        intent(inout)        :: field(:,:,:)
      type(domain2D),   intent(inout)        :: domain  ! Must be definable in mpp_update_init_comm
      integer,          intent(in), optional :: flags
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      type(domain2d), pointer :: Dom => NULL()
      integer                 :: update_position, update_whalo, update_ehalo, update_shalo, update_nhalo, ntile

      integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: f_addrs=-9999
      integer          :: tile, max_ntile
      character(len=3) :: text
      logical          :: set_mismatch, is_complete
      logical          :: do_update
      integer,    save :: isize=0, jsize=0, ke=0, l_size=0, list=0
      integer,    save :: pos, whalosz, ehalosz, shalosz, nhalosz
      MPP_TYPE_        :: d_type

      if(present(whalo)) then
         update_whalo = whalo
         if(abs(update_whalo) > domain%whalo ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D: "// &
                "optional argument whalo should not be larger than the whalo when define domain.")
      else
         update_whalo = domain%whalo
      end if
      if(present(ehalo)) then
         update_ehalo = ehalo
         if(abs(update_ehalo) > domain%ehalo ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D: "// &
                "optional argument ehalo should not be larger than the ehalo when define domain.")
      else
         update_ehalo = domain%ehalo
      end if
      if(present(shalo)) then
         update_shalo = shalo
         if(abs(update_shalo) > domain%shalo ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D: "// &
                "optional argument shalo should not be larger than the shalo when define domain.")
      else
         update_shalo = domain%shalo
      end if
      if(present(nhalo)) then
         update_nhalo = nhalo
         if(abs(update_nhalo) > domain%nhalo ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D: "// &
                "optional argument nhalo should not be larger than the nhalo when define domain.")
      else
         update_nhalo = domain%nhalo
      end if

      !--- when there is NINETY or MINUS_NINETY rotation for some contact, the salar data can not be on E or N-cell,
      if(present(position)) then
         if(domain%rotated_ninety .AND. ( position == EAST .OR. position == NORTH ) )  &
              call mpp_error(FATAL, 'MPP_UPDATE_3D: hen there is NINETY or MINUS_NINETY rotation, ' // &
              'can not use scalar version update_domain for data on E or N-cell' )
      end if

      ntile = size(domain%x(:))

      max_ntile = domain%max_ntile_pe
      is_complete = .true.
      if(PRESENT(complete)) then
         is_complete = complete
      end if
      tile = 1

      if(max_ntile>1) then
         if(ntile>MAX_TILES) then
            write( text,'(i3)' ) MAX_TILES
            call mpp_error(FATAL,'MPP_UPDATE_AD_3D: MAX_TILES='//text//' exceeded number of tiles on this pe.' )
         endif
         if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D: "// &
             "optional argument tile_count should be present when number of tiles on this pe is more than 1")
         tile = tile_count
      end if
      do_update = (tile == ntile) .AND. is_complete
      list = list+1
      if(list > MAX_DOMAIN_FIELDS)then
         write( text,'(i2)' ) MAX_DOMAIN_FIELDS
         call mpp_error(FATAL,'MPP_UPDATE_AD_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
      endif
      f_addrs(list, tile) = LOC(field)
      update_position = CENTER
      if(present(position)) update_position = position
      if(list == 1 .AND. tile == 1 )then
         isize=size(field,1); jsize=size(field,2); ke = size(field,3); pos = update_position
         whalosz = update_whalo; ehalosz = update_ehalo; shalosz = update_shalo; nhalosz = update_nhalo
      else
         set_mismatch = .false.
         set_mismatch = set_mismatch .OR. (isize /= size(field,1))
         set_mismatch = set_mismatch .OR. (jsize /= size(field,2))
         set_mismatch = set_mismatch .OR. (ke /= size(field,3))
         set_mismatch = set_mismatch .OR. (update_position /= pos)
         set_mismatch = set_mismatch .OR. (update_whalo /= whalosz)
         set_mismatch = set_mismatch .OR. (update_ehalo /= ehalosz)
         set_mismatch = set_mismatch .OR. (update_shalo /= shalosz)
         set_mismatch = set_mismatch .OR. (update_nhalo /= nhalosz)
         if(set_mismatch)then
            write( text,'(i2)' ) list
            call mpp_error(FATAL,'MPP_UPDATE_AD_3D: Incompatible field at count '//text//' for group update.' )
         endif
      endif
      if(is_complete) then
         l_size = list
         list = 0
      end if
      if(do_update )then
         if( domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) )then
            Dom => search_domain(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, update_position)
            call mpp_do_update_ad( f_addrs(1:l_size,1:ntile), Dom, d_type, ke, flags, name )
         end if
         l_size=0; f_addrs=-9999; isize=0;  jsize=0;  ke=0
      endif

      return

    end subroutine MPP_UPDATE_DOMAINS_AD_3D_

    subroutine MPP_UPDATE_DOMAINS_AD_4D_( field, domain, flags, complete, position, &
                                          whalo, ehalo, shalo, nhalo, name, tile_count)
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_,        intent(inout)        :: field(:,:,:,:)
      type(domain2D),   intent(inout)        :: domain  ! Must be definable in mpp_update_init_comm
      integer,          intent(in), optional :: flags
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))

      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains_ad( field3D, domain, flags, complete, position, &
                                  whalo, ehalo, shalo, nhalo, name, tile_count)

      return
    end subroutine MPP_UPDATE_DOMAINS_AD_4D_

    subroutine MPP_UPDATE_DOMAINS_AD_5D_( field, domain, flags, complete, position, &
                                          whalo, ehalo, shalo, nhalo, name, tile_count )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_,        intent(inout)        :: field(:,:,:,:,:)
      type(domain2D),   intent(inout)        :: domain  ! Must be definable in mpp_update_init_comm
      integer,          intent(in), optional :: flags
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
      pointer( ptr, field3D )
      ptr = LOC(field)
      call mpp_update_domains_ad( field3D, domain, flags, complete, position, &
                                  whalo, ehalo, shalo, nhalo, name, tile_count )
      return
    end subroutine MPP_UPDATE_DOMAINS_AD_5D_

#ifdef VECTOR_FIELD_

!VECTOR_FIELD_ is set to false for MPP_TYPE_ integer or logical.
!vector fields
    subroutine MPP_UPDATE_DOMAINS_AD_2D_V_( fieldx, fieldy, domain, flags, gridtype, complete, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count)
!updates data domain of 2D field whose computational domains have been computed
      MPP_TYPE_,        intent(inout)        :: fieldx(:,:), fieldy(:,:)
      type(domain2D),   intent(inout)        :: domain
      integer,          intent(in), optional :: flags, gridtype
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),1)
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),1)

      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains_ad( field3Dx, field3Dy, domain, flags, gridtype, complete, &
                                  whalo, ehalo, shalo, nhalo, name, tile_count )
      return
    end subroutine MPP_UPDATE_DOMAINS_AD_2D_V_


    subroutine MPP_UPDATE_DOMAINS_AD_3D_V_( fieldx, fieldy, domain, flags, gridtype, complete, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count)
!updates data domain of 3D field whose computational domains have been computed
      MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
      type(domain2D),   intent(inout)        :: domain
      integer,          intent(in), optional :: flags, gridtype
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      type(domain2d),                pointer :: domainx => NULL()
      type(domain2d),                pointer :: domainy => NULL()
      integer                                :: update_whalo, update_ehalo, update_shalo, update_nhalo, ntile
      integer                                :: grid_offset_type
      logical                                :: exchange_uv
      
      integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: f_addrsx=-9999, f_addrsy=-9999
      logical          :: do_update, is_complete
      integer, save    :: isize(2)=0,jsize(2)=0,ke=0,l_size=0, offset_type=0, list=0
      integer, save    :: whalosz, ehalosz, shalosz, nhalosz
      integer          :: tile, max_ntile
      logical          :: set_mismatch
      character(len=3) :: text
      MPP_TYPE_        :: d_type

      if(present(whalo)) then
         update_whalo = whalo
         if(abs(update_whalo) > domain%whalo ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D_V: "// &
                "optional argument whalo should not be larger than the whalo when define domain.")
      else
         update_whalo = domain%whalo
      end if
      if(present(ehalo)) then
         update_ehalo = ehalo
         if(abs(update_ehalo) > domain%ehalo ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D_V: "// &
                "optional argument ehalo should not be larger than the ehalo when define domain.")
      else
         update_ehalo = domain%ehalo
      end if
      if(present(shalo)) then
         update_shalo = shalo
         if(abs(update_shalo) > domain%shalo ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D_V: "// &
                "optional argument shalo should not be larger than the shalo when define domain.")
      else
         update_shalo = domain%shalo
      end if
      if(present(nhalo)) then
         update_nhalo = nhalo
         if(abs(update_nhalo) > domain%nhalo ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D_V: "// &
                "optional argument nhalo should not be larger than the nhalo when define domain.")
      else
         update_nhalo = domain%nhalo
      end if

      grid_offset_type = AGRID
      if( PRESENT(gridtype) ) grid_offset_type = gridtype

      exchange_uv = .false.
      if(grid_offset_type == DGRID_NE) then
         exchange_uv = .true.
         grid_offset_type = CGRID_NE
      else if( grid_offset_type == DGRID_SW ) then
         exchange_uv = .true.
         grid_offset_type = CGRID_SW
      end if

      ntile = size(domain%x(:))

      max_ntile = domain%max_ntile_pe
      is_complete = .true.
      if(PRESENT(complete)) then
         is_complete = complete
      end if
      tile = 1

      if(max_ntile>1) then
         if(ntile>MAX_TILES) then
            write( text,'(i2)' ) MAX_TILES
            call mpp_error(FATAL,'MPP_UPDATE_AD_3D_V: MAX_TILES='//text//' exceeded number of tiles on this pe.' )
         endif
         if(.NOT. present(tile_count) ) call mpp_error(FATAL, "MPP_UPDATE_AD_3D_V: "// &
             "optional argument tile_count should be present when number of tiles on this pe is more than 1")
         tile = tile_count
      end if
      do_update = (tile == ntile) .AND. is_complete
      list = list+1
      if(list > MAX_DOMAIN_FIELDS)then
         write( text,'(i2)' ) MAX_DOMAIN_FIELDS
         call mpp_error(FATAL,'MPP_UPDATE_AD_3D_V: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
      endif

      f_addrsx(list, tile) = LOC(fieldx)
      f_addrsy(list, tile) = LOC(fieldy)
      if(list == 1 .AND. tile == 1)then
         isize(1)=size(fieldx,1); jsize(1)=size(fieldx,2); ke = size(fieldx,3)
         isize(2)=size(fieldy,1); jsize(2)=size(fieldy,2)
         offset_type = grid_offset_type
         whalosz = update_whalo; ehalosz = update_ehalo; shalosz = update_shalo; nhalosz = update_nhalo
      else
         set_mismatch = .false.
         set_mismatch = set_mismatch .OR. (isize(1) /= size(fieldx,1))
         set_mismatch = set_mismatch .OR. (jsize(1) /= size(fieldx,2))
         set_mismatch = set_mismatch .OR. (ke /= size(fieldx,3))
         set_mismatch = set_mismatch .OR. (isize(2) /= size(fieldy,1))
         set_mismatch = set_mismatch .OR. (jsize(2) /= size(fieldy,2))
         set_mismatch = set_mismatch .OR. (ke /= size(fieldy,3))
         set_mismatch = set_mismatch .OR. (grid_offset_type /= offset_type)
         set_mismatch = set_mismatch .OR. (update_whalo /= whalosz)
         set_mismatch = set_mismatch .OR. (update_ehalo /= ehalosz)
         set_mismatch = set_mismatch .OR. (update_shalo /= shalosz)
         set_mismatch = set_mismatch .OR. (update_nhalo /= nhalosz)
         if(set_mismatch)then
            write( text,'(i2)' ) list
            call mpp_error(FATAL,'MPP_UPDATE_AD_3D_V: Incompatible field at count '//text//' for group vector update.' )
         end if
      end if
      if(is_complete) then
         l_size = list
         list = 0
      end if
      if(do_update)then
         if( domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) )then
            domainx => search_domain(domain, update_whalo, update_ehalo, update_shalo, update_nhalo,    &
                 gridtype=grid_offset_type, direction='x')
            domainy => search_domain(domain, update_whalo, update_ehalo, update_shalo, update_nhalo,    &
                 gridtype=grid_offset_type, direction='y')
            if(exchange_uv) then
               call mpp_do_update_ad(f_addrsx(1:l_size,1:ntile),f_addrsy(1:l_size,1:ntile), domainy, domainx, &
                    d_type, ke, grid_offset_type, flags,name)
            else
               call mpp_do_update_ad(f_addrsx(1:l_size,1:ntile),f_addrsy(1:l_size,1:ntile), domainx, domainy, &
                    d_type, ke, grid_offset_type, flags,name)
            end if
         end if
         l_size=0; f_addrsx=-9999; f_addrsy=-9999; isize=0;  jsize=0;  ke=0
      end if

      return
    end subroutine MPP_UPDATE_DOMAINS_AD_3D_V_


    subroutine MPP_UPDATE_DOMAINS_AD_4D_V_( fieldx, fieldy, domain, flags, gridtype, complete, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count )
!updates data domain of 4D field whose computational domains have been computed
      MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:), fieldy(:,:,:,:)
      type(domain2D),   intent(inout)        :: domain
      integer,          intent(in), optional :: flags, gridtype
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4))
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4))

      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains_ad( field3Dx, field3Dy, domain, flags, gridtype, complete, &
                                  whalo, ehalo, shalo, nhalo, name, tile_count )

      return
    end subroutine MPP_UPDATE_DOMAINS_AD_4D_V_

    subroutine MPP_UPDATE_DOMAINS_AD_5D_V_( fieldx, fieldy, domain, flags, gridtype, complete, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count )
!updates data domain of 5D field whose computational domains have been computed
      MPP_TYPE_,        intent(inout)        :: fieldx(:,:,:,:,:), fieldy(:,:,:,:,:)
      type(domain2D),   intent(inout)        :: domain
      integer,          intent(in), optional :: flags, gridtype
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3Dx(size(fieldx,1),size(fieldx,2),size(fieldx,3)*size(fieldx,4)*size(fieldx,5))
      MPP_TYPE_ :: field3Dy(size(fieldy,1),size(fieldy,2),size(fieldy,3)*size(fieldy,4)*size(fieldy,5))

      pointer( ptrx, field3Dx )
      pointer( ptry, field3Dy )
      ptrx = LOC(fieldx)
      ptry = LOC(fieldy)
      call mpp_update_domains_ad( field3Dx, field3Dy, domain, flags, gridtype, complete, &
                                  whalo, ehalo, shalo, nhalo, name, tile_count )
      return
    end subroutine MPP_UPDATE_DOMAINS_AD_5D_V_
#endif /* VECTOR_FIELD_ */
