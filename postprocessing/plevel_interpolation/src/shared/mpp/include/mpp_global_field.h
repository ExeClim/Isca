    subroutine MPP_GLOBAL_FIELD_2D_( domain, local, global, flags, position,tile_count)
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:)
      MPP_TYPE_, intent(out) :: global(:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position
      integer, intent(in), optional :: tile_count

      MPP_TYPE_ :: local3D (size( local,1),size( local,2),1)
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),1)
      pointer( lptr,  local3D )
      pointer( gptr, global3D )
      lptr = LOC( local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, position,tile_count )

    end subroutine MPP_GLOBAL_FIELD_2D_

    subroutine MPP_GLOBAL_FIELD_3D_( domain, local, global, flags, position, tile_count)
!get a global field from a local field
!local field may be on compute OR data domain
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position
      integer, intent(in), optional :: tile_count

      integer :: ishift, jshift
      integer :: tile

      tile = 1; if(PRESENT(tile_count)) tile = tile_count

      call mpp_get_domain_shift(domain, ishift, jshift, position)
      call mpp_do_global_field( domain, local, global, tile, ishift, jshift, flags)

    end subroutine MPP_GLOBAL_FIELD_3D_

    subroutine MPP_GLOBAL_FIELD_4D_( domain, local, global, flags, position,tile_count )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position
      integer, intent(in), optional :: tile_count

      MPP_TYPE_ :: local3D (size( local,1),size( local,2),size( local,3)*size(local,4))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),size(global,3)*size(local,4))
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, position,tile_count )
    end subroutine MPP_GLOBAL_FIELD_4D_

    subroutine MPP_GLOBAL_FIELD_5D_( domain, local, global, flags, position,tile_count )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position
      integer, intent(in), optional :: tile_count

      MPP_TYPE_ :: local3D (size( local,1),size( local,2),size( local,3)*size( local,4)*size(local,5))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),size(global,3)*size(global,4)*size(local,5))
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, position,tile_count )
    end subroutine MPP_GLOBAL_FIELD_5D_
