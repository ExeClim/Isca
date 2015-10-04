!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    subroutine MPP_TRANSMIT_( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, send_request )
!a message-passing routine intended to be reminiscent equally of both MPI and SHMEM

!put_data and get_data are contiguous MPP_TYPE_ arrays

!at each call, your put_data array is put to   to_pe's get_data
!              your get_data array is got from from_pe's put_data
!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get

!special PE designations:
!      NULL_PE: to disable a put or a get (e.g at boundaries)
!      ANY_PE:  if remote PE for the put or get is to be unspecific
!      ALL_PES: broadcast and collect operations (collect not yet implemented)

!ideally we would not pass length, but this f77-style call performs better (arrays passed by address, not descriptor)
!further, this permits <length> contiguous words from an array of any rank to be passed (avoiding f90 rank conformance check)

!caller is responsible for completion checks (mpp_sync_self) before and after

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(*)
      MPP_TYPE_, intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      integer, intent(out), optional :: recv_request, send_request

      integer :: i, outunit
      MPP_TYPE_, allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)
      integer(LONG_KIND),     save :: get_data_addr=-9999
      MPP_TYPE_                    :: get_data_local(get_len_nocomm)
      pointer(ptr_get_data, get_data_local)


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return
      
      outunit = stdout()
      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
          if( allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been deallocated by prior receive.' )
          if( get_len_nocomm > 0) then  ! pre-post recv
             ptr_get_data = get_data_addr
             do i = 1,get_len_nocomm
                get_data_local(i) = put_data(i)
             end do
             get_len_nocomm = 0
             get_data_addr = -9999
          else
             allocate( local_data(put_len) )
             do i = 1,put_len
                local_data(i) = put_data(i)
             end do
          endif
      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len )call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
              end if
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets

      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get: for libSMA, a get means do a wait to ensure put on remote PE is complete
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
          if( .NOT.allocated(local_data) ) then
             get_data_addr = LOC(get_data)
             get_len_nocomm  = get_len
          else
             do i = 1,get_len
                get_data(i) = local_data(i)
             end do
             deallocate(local_data)
          endif
      else if( from_pe.EQ.ANY_PE )then

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning, and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if
      return
    end subroutine MPP_TRANSMIT_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_BROADCAST_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      return
    end subroutine MPP_BROADCAST_

!####################################################################################
#include <mpp_transmit.inc>

