      module parkind

      implicit none
      save

!------------------------------------------------------------------
! rrtmg kinds
! Define integer and real kinds for various types.
!
! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!
!     integer kinds
!     -------------
!

#ifdef OVERLOAD_R4
          integer, parameter :: kind_ib = selected_int_kind(6)  ! 4 byte integer
#else
          integer, parameter :: kind_ib = selected_int_kind(13)  ! 8 byte integer
#endif
      integer, parameter :: kind_im = selected_int_kind(6)   ! 4 byte integer
      integer, parameter :: kind_in = kind(1)                ! native integer

!
!     real kinds
!     ----------
!
#ifdef OVERLOAD_R4
          integer, parameter :: kind_rb = selected_real_kind(6)  ! 4 byte real
#else
          integer, parameter :: kind_rb = selected_real_kind(13)  ! 8 byte real
#endif
      integer, parameter :: kind_rm = selected_real_kind(6)  ! 4 byte real
      integer, parameter :: kind_rn = kind(1.0)              ! native real

      end module parkind
