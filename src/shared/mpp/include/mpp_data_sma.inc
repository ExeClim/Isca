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
! $Id: mpp_data_sma.inc,v 19.0 2012/01/06 22:01:19 fms Exp $

!----------------------------------------------------------------!
! The following data is used in mpp_mod and its components       !
!----------------------------------------------------------------!
real(DOUBLE_KIND)                       :: mpp_stack(1)
integer                                 :: mpp_from_pe
integer(LONG_KIND), dimension(0:MAXPES) :: status, remote_data_loc
integer :: sync(SHMEM_REDUCE_SYNC_SIZE+SHMEM_BCAST_SYNC_SIZE+SHMEM_BARRIER_SYNC_SIZE)


!we call shpalloc in mpp_init() to ensure all these are remotely accessible
!on PVP where shpalloc doesn't exist, module variables are automatically
!guaranteed to be remotely accessible
pointer( ptr_stack, mpp_stack )
pointer( ptr_sync, sync )
pointer( ptr_status, status )
pointer( ptr_from, mpp_from_pe )
pointer( ptr_remote, remote_data_loc )


!--- these are some dummy data
integer, parameter :: stat = -999

!-------------------------------------------------------------------!
! The following data is used in mpp_domains_mod and its components  !
!-------------------------------------------------------------------!
real(DOUBLE_KIND) :: mpp_domains_stack(1)
real(DOUBLE_KIND) :: mpp_domains_stack_nonblock(1)
pointer( ptr_domains_stack, mpp_domains_stack )  
pointer( ptr_domains_stack_nonblock, mpp_domains_stack_nonblock )  
