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

module spec_mpp_mod
!Balaji (GFDL.Climate.Model.Info@noaa.gov)
!This module holds the data for the domains used by the spectral transform module

!This is the version for the transpose method
  use fms_mod,         only: mpp_pe, mpp_root_pe, mpp_npes, write_version_number, mpp_error, FATAL

  use mpp_domains_mod, only: mpp_domains_init, domain1D, domain2D, GLOBAL_DATA_DOMAIN, &
                             mpp_define_domains, mpp_get_compute_domain, mpp_get_compute_domains, &
                             mpp_get_domain_components, mpp_get_pelist

  implicit none
  private
  
  character(len=128), private :: version = '$Id: spec_mpp.F90,v 14.0 2007/03/15 22:12:44 fms Exp $'
  character(len=128), private :: tagname = '$Name: siena_201211 $'
  type(domain2D), save, public :: grid_domain,  spectral_domain,  global_spectral_domain
  logical, private :: module_is_initialized=.FALSE.
  integer, private :: pe, npes

  public :: spec_mpp_init, get_grid_domain, get_spec_domain, spec_mpp_end, atmosphere_domain

  contains

!=======================================================================================================================

    subroutine spec_mpp_init( num_fourier, num_spherical, num_lon, lat_max, grid_layout, spectral_layout )
      integer, intent(in) ::  num_fourier, num_spherical, num_lon, lat_max
      integer, intent(in), optional :: grid_layout(2), spectral_layout(2)
      integer :: layout(2)
      character(len=4) :: chtmp1, chtmp2

      if( module_is_initialized ) return
      call mpp_domains_init()
      pe = mpp_pe()
      npes = mpp_npes()

      call write_version_number(version, tagname)

!grid domain: by default, 1D decomposition along Y
      layout = (/1,npes/)
      if( PRESENT(grid_layout) ) layout = grid_layout
      call mpp_define_domains( (/1,num_lon,1,lat_max/), layout, grid_domain )
      if(pe == mpp_root_pe()) call print_decomp (npes, layout, grid_domain )

!requirement of equal domains: can be generalized to retain mirror symmetry between N/S if unequal.
!the equal-domains requirement permits us to eliminate one buffer/unbuffer in the transpose_fourier routines.
      if( mod(lat_max,layout(2)).NE.0 ) then
!       call mpp_error( FATAL, 'SPEC_MPP_INIT: currently requires equal grid domains on all PEs.' )
        write(chtmp1,'(i4)') layout(2)
        write(chtmp2,'(i4)') lat_max
        call mpp_error( FATAL, 'SPEC_MPP_INIT:Requires num_lat_rows/num_pes=int;num_pes='&
       &//chtmp1//';num_lat_rows='//chtmp2 )
      endif

!spectral domain: by default, 1D decomposition along M
      layout=(/npes,1/)
      if( PRESENT(spectral_layout) ) layout = spectral_layout
      call mpp_define_domains( (/0,num_fourier,0,num_spherical/), layout, spectral_domain )

!global spectral domains (may be used for I/O) are the same as spectral domains, with global data boundaries
      call mpp_define_domains( (/0,num_fourier,0,num_spherical/), layout, global_spectral_domain, &
           xflags=GLOBAL_DATA_DOMAIN, yflags=GLOBAL_DATA_DOMAIN )

      module_is_initialized=.TRUE.
      return
    end subroutine spec_mpp_init
!=======================================================================================================================

subroutine print_decomp (npes, layout, Domain)
integer, intent(in) :: npes, layout(2)
type(domain2d), intent(in) :: Domain
integer, dimension(0:npes-1) :: xsize, ysize
integer :: i, j, xlist(layout(1)), ylist(layout(2))
type (domain1D) :: Xdom, Ydom

call mpp_get_compute_domains   ( Domain, xsize=xsize, ysize=ysize )
call mpp_get_domain_components ( Domain, Xdom, Ydom )
call mpp_get_pelist ( Xdom, xlist )
call mpp_get_pelist ( Ydom, ylist )

write (*,100)
write (*,110) (xsize(xlist(i)),i=1,layout(1))
write (*,120) (ysize(ylist(j)),j=1,layout(2))

100 format ('ATMOS MODEL DOMAIN DECOMPOSITION')
110 format ('  X-AXIS = ',24i4,/,(11x,24i4))
120 format ('  Y-AXIS = ',24i4,/,(11x,24i4))

end subroutine print_decomp
!=======================================================================================================================

subroutine get_grid_domain(is, ie, js, je)
integer, intent(out) :: is, ie, js, je

if(.not.module_is_initialized) call mpp_error( FATAL, 'subroutine get_grid_domain: spec_mpp is not initialized')

call mpp_get_compute_domain(grid_domain, is, ie, js, je)

return
end subroutine get_grid_domain
!=======================================================================================================================

subroutine get_spec_domain(ms, me, ns, ne)
integer, intent(out) :: ms, me, ns, ne

if(.not.module_is_initialized) call mpp_error( FATAL, 'subroutine get_spec_domain: spec_mpp is not initialized')

call mpp_get_compute_domain(spectral_domain, ms, me, ns, ne)

return
end subroutine get_spec_domain
!=======================================================================================================================

subroutine spec_mpp_end

module_is_initialized = .false.

return
end subroutine spec_mpp_end
!=======================================================================================================================

subroutine atmosphere_domain(Domain)
type(domain2d), intent(inout) :: Domain

Domain = grid_domain

end subroutine atmosphere_domain
!=======================================================================================================================

end module spec_mpp_mod
