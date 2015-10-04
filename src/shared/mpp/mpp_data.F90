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

module mpp_data_mod
#include <fms_platform.h>

#if defined(use_libMPI) && defined(sgi_mipspro)
  use mpi
#endif

  use mpp_parameter_mod, only : MAXPES

  implicit none
  private

  character(len=128), public :: version= &
       '$Id mpp_data.F90 $'
  character(len=128), public :: tagname= &
       '$Name: siena_201211 $'

#if defined(use_libSMA) || defined(use_MPI_SMA)
#include <mpp/shmem.fh>
#endif

#if defined(use_libMPI) && !defined(sgi_mipspro)
#include <mpif.h>  
!sgi_mipspro gets this from 'use mpi'
#endif

  !--- public data is used by mpp_mod
  public :: stat, mpp_stack, ptr_stack, status, ptr_status, sync, ptr_sync  
  public :: mpp_from_pe, ptr_from, remote_data_loc, ptr_remote

  !--- public data which is used by mpp_domains_mod. 
  !--- All othere modules should import these parameters from mpp_domains_mod. 
  public :: mpp_domains_stack, ptr_domains_stack
  public :: mpp_domains_stack_nonblock, ptr_domains_stack_nonblock

  !-------------------------------------------------------------------------------!
  ! The following data included in the .inc file are diffrent for sma or mpi case !
  !-------------------------------------------------------------------------------!

#ifdef use_libSMA
#include <mpp_data_sma.inc>
#else
#ifdef use_libMPI
#include <mpp_data_mpi.inc>
#else
#include <mpp_data_nocomm.inc>
#endif
#endif

end module mpp_data_mod
