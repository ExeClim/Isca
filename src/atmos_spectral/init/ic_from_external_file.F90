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

module ic_from_external_file_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use               fms_mod, only: error_mesg, FATAL, field_size, stdlog, file_exist, field_exist, &
                                 write_version_number, close_file, check_nml_error, read_data

use        transforms_mod, only: get_grid_boundaries, get_deg_lon, get_deg_lat, trans_grid_to_spherical, &
                                 trans_spherical_to_grid, grid_domain, spectral_domain, get_grid_domain, &
                                 get_spec_domain, get_lat_max, get_lon_max, vor_div_from_uv_grid, uv_grid_from_vor_div

use       tracer_type_mod, only: tracer_type

use     field_manager_mod, only: MODEL_ATMOS

use    tracer_manager_mod, only: get_number_tracers

implicit none
private

character(len=128), parameter :: version = &
'$Id: ic_from_external_file.F90,v 19.0 2012/01/06 20:35:31 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: siena_201211 $'

public :: ic_from_external_file

character(len=128) :: file_name = 'INPUT/init_cond.nc'
character(len=128) :: u_name    = 'u'
character(len=128) :: v_name    = 'v'
character(len=128) :: t_name    = 't'
character(len=128) :: ps_name   = 'ps'

namelist / ic_from_external_file_nml / file_name, u_name, v_name, t_name, ps_name

contains

!=========================================================================================================================

subroutine ic_from_external_file(triang_trunc, tracer_attributes, vors, divs, ts, ln_ps, ug, &
                                 vg, tg, psg, vorg, divg, grid_tracers)

logical, intent(in) :: triang_trunc
type(tracer_type),intent(inout), dimension(:) :: tracer_attributes
complex, intent(out), dimension(:,:,:)   :: vors, divs, ts
complex, intent(out), dimension(:,:  )   :: ln_ps
real,    intent(out), dimension(:,:,:)   :: ug, vg, tg
real,    intent(out), dimension(:,:  )   :: psg
real,    intent(out), dimension(:,:,:)   :: vorg, divg
real,    intent(out), dimension(:,:,:,:) :: grid_tracers

integer :: unit, ierr, io, lon_max, lat_max, siz(4), num_tracers, ntr
character(len=128) :: tr_name, mesg1, mesg2
real, dimension(size(psg,1), size(psg,2)) :: ln_psg

!-------------------------------------------------------------------------------
! Read the namelist
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ic_from_external_file_nml, iostat=io)
    ierr = check_nml_error(io, 'ic_from_external_file_nml')
#else
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=ic_from_external_file_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'ic_from_external_file_nml')
    enddo
20  call close_file (unit)
#endif
call write_version_number(version, tagname)
write (stdlog(), nml=ic_from_external_file_nml)

!-------------------------------------------------------------------------------
! Error checks
if(.not.file_exist(file_name)) then
  call error_mesg('ic_from_external_file',trim(file_name)//' does not exist', FATAL)
endif
if( .not.field_exist(file_name,u_name) .or. .not.field_exist(file_name,v_name)  .or. &
    .not.field_exist(file_name,t_name) .or. .not.field_exist(file_name,ps_name)) then
  call error_mesg('ic_from_external_file','One or more of '//trim(u_name)//', '// &
      trim(v_name)//', '//trim(t_name)//', '//trim(ps_name)// &
      ' does not exist in '//trim(file_name), FATAL)
endif 
call field_size(file_name, u_name, siz) ! Checking the size of just one field should suffice.
call get_lon_max(lon_max)
call get_lat_max(lat_max)

if(any(siz(1:3) /= (/lon_max,lat_max,size(ug,3)/))) then
  write(mesg1,'(3i5)') siz(1:3) 
  write(mesg2,'(3i5)') (/lon_max,lat_max,size(ug,3)/)
  call error_mesg('ic_from_external_file',trim(u_name)//' in file '//trim(file_name)//' has wrong dimensions'// &
  ' Its shape is '//trim(mesg1)//'  Should be '//trim(mesg2), FATAL)
endif

!-------------------------------------------------------------------------------
! Read the initial state

call read_data(file_name, u_name,   ug, domain=grid_domain)
call read_data(file_name, v_name,   vg, domain=grid_domain)
call read_data(file_name, t_name,   tg, domain=grid_domain)
call read_data(file_name, ps_name, psg, domain=grid_domain)

! Read the initial tracers
call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
do ntr=1,num_tracers
  tr_name = trim(tracer_attributes(ntr)%name)
  if(.not.field_exist(file_name,tr_name)) then
    call error_mesg('ic_from_external_file', &
                    trim(tr_name)//' is listed in field_table but does not exist in '//trim(file_name),FATAL)
  endif
  call read_data(file_name, tr_name, grid_tracers(:,:,:,ntr), domain=grid_domain)
enddo

!-------------------------------------------------------------------------------
! initial spectral fields (and spectrally-filtered) grid fields

call trans_grid_to_spherical(tg, ts)
call trans_spherical_to_grid(ts, tg)

ln_psg = alog(psg)
call trans_grid_to_spherical(ln_psg, ln_ps)
call trans_spherical_to_grid(ln_ps,  ln_psg)
psg = exp(ln_psg)

call vor_div_from_uv_grid(ug, vg, vors, divs, triang=triang_trunc)
call uv_grid_from_vor_div(vors, divs, ug, vg)
call trans_spherical_to_grid(vors, vorg)
call trans_spherical_to_grid(divs, divg)

return
end subroutine ic_from_external_file

!================================================================================

end module ic_from_external_file_mod
