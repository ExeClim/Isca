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
subroutine MPP_UPDATE_NEST_FINE_2D_(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
                                    flags, complete, position, extra_halo, name, tile_count) 
      MPP_TYPE_,             intent(in)      :: field(:,:)
      type(nest_domain_type), intent(inout)  :: nest_domain
      MPP_TYPE_,             intent(inout)   :: wbuffer(:,:)
      MPP_TYPE_,             intent(inout)   :: ebuffer(:,:)
      MPP_TYPE_,             intent(inout)   :: sbuffer(:,:)
      MPP_TYPE_,             intent(inout)   :: nbuffer(:,:)
      integer,          intent(in), optional :: flags
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      integer,          intent(in), optional :: extra_halo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
      MPP_TYPE_ :: wbuffer3D(size(wbuffer,1),size(wbuffer,2),1)
      MPP_TYPE_ :: ebuffer3D(size(ebuffer,1),size(ebuffer,2),1)
      MPP_TYPE_ :: sbuffer3D(size(sbuffer,1),size(sbuffer,2),1)
      MPP_TYPE_ :: nbuffer3D(size(nbuffer,1),size(nbuffer,2),1)
      pointer( ptr, field3D )
      pointer( ptr_w, wbuffer3D)
      pointer( ptr_e, ebuffer3D)
      pointer( ptr_s, sbuffer3D)
      pointer( ptr_n, nbuffer3D)
      ptr = LOC(field)
      ptr_w = LOC(wbuffer)
      ptr_e = LOC(ebuffer)
      ptr_s = LOC(sbuffer)
      ptr_n = LOC(nbuffer)
      call mpp_update_nest_fine( field3D, nest_domain, wbuffer3D, ebuffer3D, sbuffer3D, nbuffer3D, &
                                 flags, complete, position, extra_halo, name, tile_count) 

      return


end subroutine MPP_UPDATE_NEST_FINE_2D_

subroutine MPP_UPDATE_NEST_FINE_3D_(field, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, &
                                    flags, complete, position, extra_halo, name, tile_count) 
    MPP_TYPE_,             intent(in)      :: field(:,:,:)
    type(nest_domain_type), intent(inout)  :: nest_domain
    MPP_TYPE_,             intent(inout)   :: wbuffer(:,:,:)
    MPP_TYPE_,             intent(inout)   :: ebuffer(:,:,:)
    MPP_TYPE_,             intent(inout)   :: sbuffer(:,:,:)
    MPP_TYPE_,             intent(inout)   :: nbuffer(:,:,:)
    integer,          intent(in), optional :: flags
    logical,          intent(in), optional :: complete
    integer,          intent(in), optional :: position
    integer,          intent(in), optional :: extra_halo
    character(len=*), intent(in), optional :: name
    integer,          intent(in), optional :: tile_count

   MPP_TYPE_        :: d_type
   type(nestSpec), pointer :: update=>NULL()
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: f_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: wb_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: eb_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: sb_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: nb_addrs=-9999
   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile
   integer          :: add_halo, update_flags, update_position
   integer          :: wbuffersz, ebuffersz, sbuffersz, nbuffersz
   integer          :: isize, jsize, ksize, l_size
   integer, save    :: isize_save, jsize_save, ksize_save
   integer          :: wbuffersz_save, ebuffersz_save, sbuffersz_save, nbuffersz_save
   integer, save    :: add_halo_save, update_flags_save, update_position_save
   integer, save    :: list=0 

   add_halo = 0
   if(present(extra_halo)) add_halo = add_halo
   update_position = CENTER
   if(present(position)) update_position = position
   update_flags = XUPDATE+YUPDATE   !default
   if( PRESENT(flags) )update_flags = flags


   is_complete = .true.
   if(PRESENT(complete)) then
      is_complete = complete
   end if
   tile = 1
   if(present(tile_count)) tile = tile_count
   if( tile > 1 ) then
      call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D: currently do not support multiple tile per pe')
   endif

   list = list+1
   if(list > MAX_DOMAIN_FIELDS)then
      write( text,'(i2)' ) MAX_DOMAIN_FIELDS
      call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
   endif

   f_addrs(list) = LOC(field)
   wb_addrs(list) = LOC(wbuffer)
   eb_addrs(list) = LOC(ebuffer)
   sb_addrs(list) = LOC(sbuffer)
   nb_addrs(list) = LOC(nbuffer)

   wbuffersz = size(wbuffer); ebuffersz = size(ebuffer)
   sbuffersz = size(sbuffer); nbuffersz = size(nbuffer)
   isize=size(field,1); jsize=size(field,2); ksize = size(field,3)
   if(list == 1)then
      isize_save = isize; jsize_save = jsize; ksize_save = ksize
      update_position_save = update_position
      update_flags_save    = update_flags
      wbuffersz_save = wbuffersz; ebuffersz_save = ebuffersz
      sbuffersz_save = sbuffersz; nbuffersz_save = nbuffersz
      add_halo_save = add_halo
   else
      set_mismatch = .false.
      set_mismatch = set_mismatch .OR. (isize_save /= isize)
      set_mismatch = set_mismatch .OR. (jsize_save /= jsize)
      set_mismatch = set_mismatch .OR. (ksize_save /= ksize)
      set_mismatch = set_mismatch .OR. (update_position_save /= update_position)
      set_mismatch = set_mismatch .OR. (wbuffersz_save /= wbuffersz)
      set_mismatch = set_mismatch .OR. (ebuffersz_save /= ebuffersz)
      set_mismatch = set_mismatch .OR. (sbuffersz_save /= sbuffersz)
      set_mismatch = set_mismatch .OR. (nbuffersz_save /= nbuffersz)
      set_mismatch = set_mismatch .OR. (update_flags_save /= update_flags)
      set_mismatch = set_mismatch .OR. (add_halo_save /= add_halo)

      if(set_mismatch)then
         write( text,'(i2)' ) list
         call mpp_error(FATAL,'MPP_UPDATE_NEST_FINE_3D_: Incompatible field at count '//text//' for group update.' )
      endif
   endif

   if(is_complete) then
      l_size = list
      list = 0
   end if
      
   if(is_complete)then
      update => search_C2F_nest_overlap(nest_domain, add_halo, update_position)
      call mpp_do_update_nest_fine(f_addrs(1:l_size), nest_domain, update, d_type, ksize, &
            wb_addrs(1:l_size), eb_addrs(1:l_size), sb_addrs(1:l_size), nb_addrs(1:l_size), update_flags )

   endif


end subroutine MPP_UPDATE_NEST_FINE_3D_


!###############################################################################
subroutine MPP_UPDATE_NEST_FINE_4D_(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
                                    flags, complete, position, extra_halo, name, tile_count) 
      MPP_TYPE_,             intent(in)      :: field(:,:,:,:)
      type(nest_domain_type), intent(inout)  :: nest_domain
      MPP_TYPE_,             intent(inout)   :: wbuffer(:,:,:,:)
      MPP_TYPE_,             intent(inout)   :: ebuffer(:,:,:,:)
      MPP_TYPE_,             intent(inout)   :: sbuffer(:,:,:,:)
      MPP_TYPE_,             intent(inout)   :: nbuffer(:,:,:,:)
      integer,          intent(in), optional :: flags
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      integer,          intent(in), optional :: extra_halo
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
      MPP_TYPE_ :: wbuffer3D(size(wbuffer,1),size(wbuffer,2),size(wbuffer,3)*size(wbuffer,4))
      MPP_TYPE_ :: ebuffer3D(size(ebuffer,1),size(ebuffer,2),size(ebuffer,3)*size(ebuffer,4))
      MPP_TYPE_ :: sbuffer3D(size(sbuffer,1),size(sbuffer,2),size(sbuffer,3)*size(sbuffer,4))
      MPP_TYPE_ :: nbuffer3D(size(nbuffer,1),size(nbuffer,2),size(nbuffer,3)*size(nbuffer,4))

      pointer( ptr, field3D )
      pointer( ptr_w, wbuffer3D)
      pointer( ptr_e, ebuffer3D)
      pointer( ptr_s, sbuffer3D)
      pointer( ptr_n, nbuffer3D)
      ptr = LOC(field)
      ptr_w = LOC(wbuffer)
      ptr_e = LOC(ebuffer)
      ptr_s = LOC(sbuffer)
      ptr_n = LOC(nbuffer)
      call mpp_update_nest_fine( field3D, nest_domain, wbuffer3D, ebuffer3D, sbuffer3D, nbuffer3D, &
                                 flags, complete, position, extra_halo, name, tile_count) 

      return


end subroutine MPP_UPDATE_NEST_FINE_4D_



subroutine MPP_UPDATE_NEST_COARSE_2D_(field, nest_domain, buffer, complete, position, name, tile_count) 
      MPP_TYPE_,             intent(in)      :: field(:,:)
      type(nest_domain_type), intent(inout)  :: nest_domain
      MPP_TYPE_,             intent(inout)   :: buffer(:,:)
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
      MPP_TYPE_ :: buffer3D(size(buffer,1),size(buffer,2),1)
      pointer( ptr, field3D )
      pointer( ptr_b, buffer3D)
      ptr = LOC(field)
      ptr_b = LOC(buffer)
      call mpp_update_nest_coarse( field3D, nest_domain, buffer3D, complete, position, name, tile_count) 

      return


end subroutine MPP_UPDATE_NEST_COARSE_2D_


subroutine MPP_UPDATE_NEST_COARSE_3D_(field, nest_domain, buffer, complete, position, name, tile_count) 
   MPP_TYPE_,             intent(in)      :: field(:,:,:)
   type(nest_domain_type), intent(inout)  :: nest_domain
   MPP_TYPE_,             intent(inout)   :: buffer(:,:,:)
   logical,          intent(in), optional :: complete
   integer,          intent(in), optional :: position
   character(len=*), intent(in), optional :: name
   integer,          intent(in), optional :: tile_count

   MPP_TYPE_        :: d_type
   type(nestSpec), pointer :: update=>NULL()
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: f_addrs=-9999
   integer(LONG_KIND),dimension(MAX_DOMAIN_FIELDS),save :: b_addrs=-9999
   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile
   integer          :: update_position
   integer          :: buffersz, buffersz_save
   integer          :: isize, jsize, ksize, l_size
   integer, save    :: isize_save, jsize_save, ksize_save
   integer, save    :: update_position_save
   integer, save    :: list=0 

   update_position = CENTER
   if(present(position)) update_position = position

   is_complete = .true.
   if(PRESENT(complete)) then
      is_complete = complete
   end if
   tile = 1
   if(present(tile_count)) tile = tile_count
   if( tile > 1 ) then
      call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D: currently do not support multiple tile per pe')
   endif

   list = list+1
   if(list > MAX_DOMAIN_FIELDS)then
      write( text,'(i2)' ) MAX_DOMAIN_FIELDS
      call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
   endif

   f_addrs(list) = LOC(field)
   b_addrs(list) = LOC(buffer)

   buffersz = size(buffer)
   isize=size(field,1); jsize=size(field,2); ksize = size(field,3)
   if(list == 1)then
      isize_save = isize; jsize_save = jsize; ksize_save = ksize
      update_position_save = update_position
      buffersz_save = buffersz
   else
      set_mismatch = .false.
      set_mismatch = set_mismatch .OR. (isize_save /= isize)
      set_mismatch = set_mismatch .OR. (jsize_save /= jsize)
      set_mismatch = set_mismatch .OR. (ksize_save /= ksize)
      set_mismatch = set_mismatch .OR. (update_position_save /= update_position)
      set_mismatch = set_mismatch .OR. (buffersz_save /= buffersz)

      if(set_mismatch)then
         write( text,'(i2)' ) list
         call mpp_error(FATAL,'MPP_UPDATE_NEST_COARSE_3D_: Incompatible field at count '//text//' for group update.' )
      endif
   endif

   if(is_complete) then
      l_size = list
      list = 0
   end if
      
   if(is_complete)then
      update => search_F2C_nest_overlap(nest_domain, update_position)
      call mpp_do_update_nest_coarse(f_addrs(1:l_size), nest_domain, update, d_type, ksize, &
            b_addrs(1:l_size))
   endif

end subroutine MPP_UPDATE_NEST_COARSE_3D_

!###############################################################################
subroutine MPP_UPDATE_NEST_COARSE_4D_(field, nest_domain, buffer, complete, position, name, tile_count) 
      MPP_TYPE_,             intent(in)      :: field(:,:,:,:)
      type(nest_domain_type), intent(inout)  :: nest_domain
      MPP_TYPE_,             intent(inout)   :: buffer(:,:,:,:)
      logical,          intent(in), optional :: complete
      integer,          intent(in), optional :: position
      character(len=*), intent(in), optional :: name
      integer,          intent(in), optional :: tile_count

      MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
      MPP_TYPE_ :: buffer3D(size(buffer,1),size(buffer,2),size(buffer,3)*size(buffer,4))

      pointer( ptr, field3D )
      pointer( ptr_b, buffer3D)
      ptr = LOC(field)
      ptr_b = LOC(buffer)
      call mpp_update_nest_coarse( field3D, nest_domain, buffer3D, complete, position, name, tile_count) 

      return


end subroutine MPP_UPDATE_NEST_COARSE_4D_
