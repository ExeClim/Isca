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

module spectral_init_cond_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use               fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, field_size, stdlog, file_exist, &
                                 write_version_number, close_file, check_nml_error, read_data

use               mpp_mod, only: mpp_chksum

use       mpp_domains_mod, only: mpp_get_global_domain

use         constants_mod, only: grav, pi

use   vert_coordinate_mod, only: compute_vert_coord

use        transforms_mod, only: get_grid_boundaries, get_deg_lon, get_deg_lat, trans_grid_to_spherical, &
                                 trans_spherical_to_grid, grid_domain, spectral_domain, get_grid_domain, &
                                 get_spec_domain

use  press_and_geopot_mod, only: press_and_geopot_init, pressure_variables

use spectral_initialize_fields_mod, only: spectral_initialize_fields

use      ic_from_external_file_mod, only: ic_from_external_file

use       topog_regularization_mod, only: compute_lambda, regularize

use                 topography_mod, only: gaussian_topog_init, get_topog_mean, get_ocean_mask

use                tracer_type_mod, only: tracer_type

implicit none
private

character(len=128), parameter :: version = &
'$Id: spectral_init_cond.F90,v 19.0 2012/01/06 20:35:39 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: siena_201211 $'

public :: spectral_init_cond

real :: initial_temperature=264.
character(len=64) :: topography_option = 'flat'  ! realistic topography computed from high resolution raw data
character(len=64) :: topog_file_name  = 'topography.data.nc'
character(len=64) :: topog_field_name = 'zsurf'
character(len=256) :: land_field_name = 'land_mask'

namelist / spectral_init_cond_nml / initial_temperature, topography_option, topog_file_name, topog_field_name, land_field_name

contains

!=========================================================================================================================

subroutine spectral_init_cond(initial_state_option, tracer_attributes, reference_sea_level_press, triang_trunc, use_virtual_temperature, &
                              vert_coord_option, vert_difference_option, scale_heights, surf_res,    &
                              p_press, p_sigma, exponent, ocean_topog_smoothing, pk, bk, vors, divs, &
                              ts, ln_ps, ug, vg, tg, psg, vorg, divg, grid_tracers, surf_geopotential, ocean_mask)

character(len=*), intent(in) :: initial_state_option
type(tracer_type), intent(inout), dimension(:) :: tracer_attributes
real,    intent(in) :: reference_sea_level_press
logical, intent(in) :: triang_trunc, use_virtual_temperature
character(len=*), intent(in) :: vert_coord_option, vert_difference_option
real,    intent(in) :: scale_heights, surf_res, p_press, p_sigma, exponent, ocean_topog_smoothing
real,    intent(out), dimension(:)       :: pk, bk
complex, intent(out), dimension(:,:,:)   :: vors, divs, ts
complex, intent(out), dimension(:,:  )   :: ln_ps
real,    intent(out), dimension(:,:,:)   :: ug, vg, tg
real,    intent(out), dimension(:,:  )   :: psg
real,    intent(out), dimension(:,:,:)   :: vorg, divg
real,    intent(out), dimension(:,:,:,:) :: grid_tracers
real,    intent(out), dimension(:,:  )   :: surf_geopotential
logical, optional, intent(in), dimension(:,:) :: ocean_mask

integer :: unit, ierr, io

!------------------------------------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=spectral_init_cond_nml, iostat=io)
    ierr = check_nml_error(io, 'spectral_init_cond_nml')
#else
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=spectral_init_cond_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'spectral_init_cond_nml')
    enddo
20  call close_file (unit)
#endif
call write_version_number(version, tagname)
if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=spectral_init_cond_nml)

call compute_vert_coord(vert_coord_option, scale_heights, surf_res, exponent, p_press, p_sigma, reference_sea_level_press, pk,bk)

call get_topography(topography_option, ocean_topog_smoothing, surf_geopotential, ocean_mask)
call press_and_geopot_init(pk, bk, use_virtual_temperature, vert_difference_option)

if(initial_state_option == 'quiescent') then
  call spectral_initialize_fields(reference_sea_level_press, triang_trunc, initial_temperature, &
                  surf_geopotential, ln_ps, vors, divs, ts, psg, ug, vg, tg, vorg, divg)
else if(initial_state_option == 'input') then
  call ic_from_external_file(triang_trunc, tracer_attributes, vors, divs, ts, ln_ps, ug, &
                             vg, tg, psg, vorg, divg, grid_tracers)
endif

call check_vert_coord(size(ug,3), psg)

return
end subroutine spectral_init_cond

!================================================================================

subroutine check_vert_coord(num_levels, psg)
integer, intent(in) :: num_levels
real, intent(in), dimension(:,:) :: psg
real, dimension(size(psg,1), size(psg,2), num_levels  ) :: p_full, ln_p_full
real, dimension(size(psg,1), size(psg,2), num_levels+1) :: p_half, ln_p_half
integer :: i,j,k

call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, psg)
do k=1,size(p_full,3)
  do j=1,size(p_full,2)
    do i=1,size(p_full,1)
      if(p_half(i,j,k+1) < p_half(i,j,k)) then
        call error_mesg('check_vert_coord','Pressure levels intersect.',FATAL)
      endif
    enddo
  enddo
enddo

return
end subroutine check_vert_coord
!================================================================================

subroutine get_topography(topography_option, ocean_topog_smoothing, surf_geopotential, ocean_mask_in)

character(len=*), intent(in) :: topography_option
real,    intent(in) :: ocean_topog_smoothing
real,    intent(out), dimension(:,:) :: surf_geopotential
logical, intent(in), optional, dimension(:,:) :: ocean_mask_in
real,    dimension(size(surf_geopotential,1)  ) :: deg_lon
real,    dimension(size(surf_geopotential,2)  ) :: deg_lat
real,    dimension(size(surf_geopotential,1),size(surf_geopotential,2)) :: surf_height
real,    dimension(size(surf_geopotential,1),size(surf_geopotential,2)) :: land_ones
logical, dimension(size(surf_geopotential,1),size(surf_geopotential,2)) :: ocean_mask
complex, allocatable, dimension(:,:) :: spec_tmp
real :: fraction_smoothed, lambda
integer :: is, ie, js, je, ms, me, ns, ne, global_num_lon, global_num_lat
real, allocatable, dimension(:) :: blon, blat
logical :: topo_file_exists, water_file_exists
integer, dimension(4) :: siz
character(len=12) :: ctmp1='     by     ', ctmp2='     by     '

if(trim(topography_option) == 'flat') then
   surf_geopotential = 0.

else if(trim(topography_option) == 'input') then
   if(file_exist(trim('INPUT/'//topog_file_name))) then
     call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat) 
     call field_size(trim('INPUT/'//topog_file_name), trim(topog_field_name), siz)
     if ( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) then
       call read_data(trim('INPUT/'//topog_file_name), trim(topog_field_name), surf_height, grid_domain)
     else
       write(ctmp1(1: 4),'(i4)') siz(1)
       write(ctmp1(9:12),'(i4)') siz(2)
       write(ctmp2(1: 4),'(i4)') global_num_lon
       write(ctmp2(9:12),'(i4)') global_num_lat
       call error_mesg ('get_topography','Topography file contains data on a '// &
              ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
     endif
   else
     call error_mesg('get_topography','topography_option="'//trim(topography_option)//'"'// &
                     ' but '//trim('INPUT/'//topog_file_name)//' does not exist', FATAL)
   endif

!  RG add lines to read land mask as well as topography, then convert to an ocean_mask for smoothing
   if(file_exist(trim('INPUT/'//topog_file_name))) then
     call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat) 
     call field_size(trim('INPUT/'//topog_file_name), trim(land_field_name), siz)
     if ( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) then
       call read_data(trim('INPUT/'//topog_file_name), trim(land_field_name), land_ones, grid_domain)
     else
       write(ctmp1(1: 4),'(i4)') siz(1)
       write(ctmp1(9:12),'(i4)') siz(2)
       write(ctmp2(1: 4),'(i4)') global_num_lon
       write(ctmp2(9:12),'(i4)') global_num_lat
       call error_mesg ('get_topography','Land file contains data on a '// &
              ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
     endif
   else
     call error_mesg('get_topography','topography_option="'//trim(topography_option)//'"'// &
                     ' but '//trim('INPUT/'//topog_file_name)//' does not exist', FATAL)
   endif
    
   where(land_ones > 0.)
     ocean_mask = .false.
   elsewhere
     ocean_mask = .true.
   end where

   surf_geopotential = grav*surf_height

   if(ocean_topog_smoothing == 0.) then
!    Spectrally truncate the topography
     call get_spec_domain(ms, me, ns, ne)
     allocate(spec_tmp(ms:me, ns:ne))
     call trans_grid_to_spherical(surf_geopotential,spec_tmp)
     call trans_spherical_to_grid(spec_tmp,surf_geopotential)
     deallocate(spec_tmp)
   else
     call compute_lambda(ocean_topog_smoothing, ocean_mask, surf_geopotential, lambda, fraction_smoothed)

!  Note that the array surf_height is used here for the smoothed surf_geopotential,
!  then immediately loaded back into surf_geopotential
     call regularize(lambda, ocean_mask, surf_geopotential, surf_height, fraction_smoothed)
     surf_geopotential = surf_height

     if(mpp_pe() == mpp_root_pe()) then
       print '(/,"Message from subroutine get_topography:")'
       print '("lambda=",1pe16.8,"  fraction_smoothed=",1pe16.8,/)',lambda,fraction_smoothed
     endif
   endif




else if(trim(topography_option) == 'interpolated') then

!  Get realistic topography
   call get_grid_domain(is, ie, js, je)
   allocate(blon(is:ie+1), blat(js:je+1))
   call get_grid_boundaries(blon, blat)
   topo_file_exists = get_topog_mean(blon, blat, surf_height)
   if(.not.topo_file_exists) then
     call error_mesg('get_topography','topography_option="'//trim(topography_option)//'"'// &
                     ' but topography data file does not exist', FATAL)
   endif
   surf_geopotential = grav*surf_height

   if(ocean_topog_smoothing == 0.) then
!    Spectrally truncate the realistic topography 
     call get_spec_domain(ms, me, ns, ne)
     allocate(spec_tmp(ms:me, ns:ne))
     call trans_grid_to_spherical(surf_geopotential,spec_tmp)
     call trans_spherical_to_grid(spec_tmp,surf_geopotential)
     deallocate(spec_tmp)
   else
!    Do topography regularization
     if(present(ocean_mask_in)) then
       ocean_mask = ocean_mask_in
     else
       water_file_exists = get_ocean_mask(blon, blat, ocean_mask)
       if(.not.water_file_exists) then
         call error_mesg('get_topography','topography_option="'//trim(topography_option)//'"'// &
                         ' and ocean_mask is not present but water data file does not exist', FATAL)
       endif
     endif
     call compute_lambda(ocean_topog_smoothing, ocean_mask, surf_geopotential, lambda, fraction_smoothed)

!  Note that the array surf_height is used here for the smoothed surf_geopotential,
!  then immediately loaded back into surf_geopotential
     call regularize(lambda, ocean_mask, surf_geopotential, surf_height, fraction_smoothed)
     surf_geopotential = surf_height

     if(mpp_pe() == mpp_root_pe()) then
       print '(/,"Message from subroutine get_topography:")'
       print '("lambda=",1pe16.8,"  fraction_smoothed=",1pe16.8,/)',lambda,fraction_smoothed
     endif
   endif
   deallocate(blon, blat)

else if(trim(topography_option) == 'gaussian') then
   call get_deg_lon(deg_lon)
   call get_deg_lat(deg_lat)
   call gaussian_topog_init(deg_lon*pi/180, deg_lat*pi/180, surf_height)
   surf_geopotential = grav*surf_height
else 
   call error_mesg('get_topography','"'//trim(topography_option)//'" is an invalid value for topography_option.', FATAL)
endif

return
end subroutine get_topography
!=======================================================================================================
subroutine print_chksum(text, vors, divs, ts, ln_ps, ug, vg, tg, psg, vorg, divg, surf_geopotential)
character(len=*), intent(in) :: text
complex, intent(in), dimension(:,:,:) :: vors, divs, ts
complex, intent(in), dimension(:,:  ) :: ln_ps
real,    intent(in), dimension(:,:,:) :: ug, vg, tg, vorg, divg
real,    intent(in), dimension(:,:  ) :: psg, surf_geopotential

integer(kind=kind(ug)) :: chksum_vors, chksum_divs, chksum_ts, chksum_ln_ps, chksum_ug, chksum_vg
integer(kind=kind(ug)) :: chksum_tg, chksum_psg, chksum_vorg, chksum_divg, chksum_surf_geo

if (mpp_pe() == mpp_root_pe()) print '(/,a)',text

chksum_vors    = mpp_chksum(vors)
chksum_divs    = mpp_chksum(divs)
chksum_ts      = mpp_chksum(ts)
chksum_ln_ps   = mpp_chksum(ln_ps)
chksum_ug      = mpp_chksum(ug)
chksum_vg      = mpp_chksum(vg)
chksum_tg      = mpp_chksum(tg)
chksum_psg     = mpp_chksum(psg)
chksum_vorg    = mpp_chksum(vorg)
chksum_divg    = mpp_chksum(divg)
chksum_surf_geo = mpp_chksum(surf_geopotential)

if (mpp_pe() == mpp_root_pe()) then
  print '("mpp_chksum(vors   )=",z17)',chksum_vors
  print '("mpp_chksum(divs   )=",z17)',chksum_divs
  print '("mpp_chksum(ts     )=",z17)',chksum_ts
  print '("mpp_chksum(ln_ps  )=",z17)',chksum_ln_ps
  print '("mpp_chksum(ug     )=",z17)',chksum_ug
  print '("mpp_chksum(vg     )=",z17)',chksum_vg
  print '("mpp_chksum(tg     )=",z17)',chksum_tg
  print '("mpp_chksum(psg    )=",z17)',chksum_psg
  print '("mpp_chksum(vorg   )=",z17)',chksum_vorg
  print '("mpp_chksum(divg   )=",z17)',chksum_divg
  print '("mpp_chksum(surf_geopotential)=",z17)',chksum_surf_geo
endif

return
end subroutine print_chksum
!================================================================================
end module spectral_init_cond_mod
