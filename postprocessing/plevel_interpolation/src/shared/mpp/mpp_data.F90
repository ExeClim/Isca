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
       '$Name: testing $'

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
