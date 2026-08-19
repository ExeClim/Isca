module column_init_cond_mod 


#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use                constants_mod, only: grav, pi
use column_initialize_fields_mod, only: column_initialize_fields
use                      fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, field_size, stdlog, file_exist, &
                                        write_version_number, close_file, check_nml_error, read_data
use              mpp_domains_mod, only: mpp_get_global_domain     
use              mpp_domains_mod, only: mpp_get_global_domain
use                 spec_mpp_mod, only: grid_domain, get_grid_domain
use          vert_coordinate_mod, only: compute_vert_coord
use         press_and_geopot_mod, only: press_and_geopot_init, pressure_variables
use              tracer_type_mod, only: tracer_type

implicit none
private 

public :: column_init_cond

character(len=128), parameter :: version = '$Id: column_init_cond.F90,v 0.1 2018/14/11 HH:MM:SS isca Exp $'
character(len=128), parameter :: tagname = '$Name: isca_201811 $'

real :: initial_temperature = 264.
real :: surf_geopotential = 0.0
real :: surface_wind = 5. 

namelist / column_init_cond_nml / initial_temperature, surf_geopotential, surface_wind

contains

subroutine column_init_cond(initial_state_option, tracer_attributes, reference_sea_level_press, use_virtual_temperature, &
    vert_coord_option, vert_difference_option, scale_heights, surf_res,    &
    p_press, p_sigma, exponent, pk, bk, ug, vg, tg, psg, &
    grid_tracers, surf_geopotential_out)! NTL START HERE, lon_boundaries, lat_boundaries)

character(len=*), intent(in) :: initial_state_option
type(tracer_type), intent(inout), dimension(:) :: tracer_attributes
real,    intent(in) :: reference_sea_level_press
logical, intent(in) :: use_virtual_temperature
character(len=*), intent(in) :: vert_coord_option, vert_difference_option
real,    intent(in) :: scale_heights, surf_res, p_press, p_sigma, exponent
!real,    intent(in), dimension(:) :: lon_boundaries, lat_boundaries NTL START HERE 
real,    intent(out), dimension(:)       :: pk, bk
real,    intent(out), dimension(:,:,:)   :: ug, vg, tg
real,    intent(out), dimension(:,:  )   :: psg
real,    intent(out), dimension(:,:,:,:) :: grid_tracers
real,    intent(out), dimension(:,:  )   :: surf_geopotential_out

integer :: unit, ierr, io

!------------------------------------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=column_init_cond_nml, iostat=io)
    ierr = check_nml_error(io, 'column_init_cond_nml')
#else
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=column_init_cond_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'column_init_cond_nml')
    enddo
20  call close_file (unit)
#endif
call write_version_number(version, tagname)
if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=column_init_cond_nml)



call compute_vert_coord(vert_coord_option, scale_heights, surf_res, exponent, p_press, p_sigma, reference_sea_level_press, pk,bk)
surf_geopotential_out = surf_geopotential ! only option to set topography to uniform surface geopotential 
call press_and_geopot_init(pk, bk, use_virtual_temperature, vert_difference_option)

if(initial_state_option == 'default') then
  call column_initialize_fields(reference_sea_level_press, initial_temperature, surface_wind, &
  surf_geopotential_out,  psg, ug, vg, tg)
else 
!!!!!!!! NTL: intial condition from file not yet configured !!!!!!!!!!
!else if(initial_state_option == 'input') then 
!call ic_from_external_file(triang_trunc, tracer_attributes, vors, divs, ts, ln_ps, ug, &
!vg, tg, psg, vorg, divg, grid_tracers)
  call error_mesg('column_init_cond','invalid initial state, can only choose "default"', FATAL)
endif

call check_vert_coord(size(ug,3), psg)

return
end subroutine column_init_cond

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
          call error_mesg('column_init_cond: check_vert_coord','Pressure levels intersect.',FATAL)
        endif
      enddo
    enddo
  enddo
  
  return
  end subroutine check_vert_coord



end module column_init_cond_mod