
! -*-f90-*-
! $Id: mpp_util_nocomm.inc,v 18.0 2010/03/02 23:58:09 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!         MISCELLANEOUS UTILITIES: mpp_error                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mpp_error_basic( errortype, errormsg )
  !a very basic error handler
  !uses ABORT and FLUSH calls, may need to use cpp to rename
  integer,                    intent(in) :: errortype
  character(len=*), intent(in), optional :: errormsg
  character(len=256)                     :: text
  logical                                :: opened
  integer                                :: istat, errunit, outunit

  if( .NOT.module_is_initialized )call ABORT()

  select case( errortype )
  case(NOTE)
     text = 'NOTE'         !just FYI
  case(WARNING)
     text = 'WARNING'      !probable error
  case(FATAL)
     text = 'FATAL'        !fatal error
  case default
     text = 'WARNING: non-existent errortype (must be NOTE|WARNING|FATAL)'
  end select

  if( npes.GT.1 )write( text,'(a,i5)' )trim(text)//' from PE', pe   !this is the mpp part
  if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)
  
  errunit = stderr()
  outunit = stdout()

  select case( errortype )
  case(NOTE)
     write( outunit,'(a)' )trim(text)
  case default
     write( errunit,'(/a/)' )trim(text)
     write( outunit,'(/a/)' )trim(text)
     if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
        call FLUSH(outunit)
#ifdef sgi_mipspro
        call TRACE_BACK_STACK_AND_PRINT()
#endif
        call ABORT() !automatically calls traceback on Cray systems
     end if
  end select

  error_state = errortype
  return
end subroutine mpp_error_basic

!#####################################################################
!--- makes a PE set out of a PE list. A PE list is an ordered list of PEs
!--- a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI
!--- if stride is non-uniform or not a power of 2, 
!--- will return error (not required for MPI but enforced for uniformity)
function get_peset(pelist)
  integer                       :: get_peset
  integer, intent(in), optional :: pelist(:)

  if( .NOT.PRESENT(pelist) )then !set it to current_peset_num
     get_peset = current_peset_num; return
  end if

  get_peset = 0

  return

end function get_peset

!#######################################################################
 !synchronize PEs in list
subroutine mpp_sync( pelist, check )
  integer, intent(in), optional :: pelist(:)
  integer, intent(in), optional :: check

  return
end subroutine mpp_sync

!#######################################################################
  !this is to check if current PE's outstanding puts are complete
  !but we can't use shmem_fence because we are actually waiting for
  !a remote PE to complete its get
subroutine mpp_sync_self( pelist, check )
  integer, intent(in), optional :: pelist(:)
  integer, intent(in), optional :: check
  return
end subroutine mpp_sync_self

