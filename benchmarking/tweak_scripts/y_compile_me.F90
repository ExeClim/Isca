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

module main
#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif
!use mpp_domains_mod, only: domain2D
use               fms_mod, only: error_mesg, FATAL, field_size, stdlog, file_exist, field_exist, &
                                 write_version_number, close_file, check_nml_error, read_data
use y_intermediate_mod , only : grid_domain
!use        transforms_mod, only: get_grid_boundaries, get_deg_lon, get_deg_lat, trans_grid_to_spherical, &
!                                  trans_spherical_to_grid, get_grid_domain, &
!                                  get_spec_domain, get_lat_max, get_lon_max, vor_div_from_uv_grid, uv_grid_from_vor_div

use       tracer_type_mod, only: tracer_type
use     field_manager_mod, only: MODEL_ATMOS
use    tracer_manager_mod, only: get_number_tracers










contains 
subroutine  name
    print *, 'hello, world'
    
end subroutine name



end module