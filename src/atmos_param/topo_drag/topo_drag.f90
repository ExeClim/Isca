module topo_drag_mod

!==========================================================================
! TOPOGRAPHIC DRAG CLOSURE FOR GENERAL CIRCULATION MODELS -- Garner (2001)
!==========================================================================

!--------------------------------------------------------------------------
!  Calculates horizontal velocity tendency due to topographic drag
!--------------------------------------------------------------------------

  use       Fms_Mod, only: FILE_EXIST, OPEN_NAMELIST_FILE, ERROR_MESG, FATAL, &
                           READ_DATA, WRITE_DATA, CLOSE_FILE, mpp_pe, mpp_root_pe, &
                           write_version_number, stdlog
  use Constants_Mod, only: Grav,Cp_Air,Rdgas,Radius,Pi,Radian

  implicit none

  private

  character(len=128) :: version = '$Id: topo_drag.f90,v 11.0 2004/09/28 19:24:57 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'
  logical            :: module_is_initialized = .false.

  public topo_drag, topo_drag_init, topo_drag_end

contains

!#############################################################################      

  subroutine topo_drag(is,js,uwnd,vwnd,atmp,pfull,phalf,zfull,zhalf,    &
    z_pbl, taux,tauy,dtaux,dtauy,taus)

    integer,intent(in) :: is,js

!   INPUT
!   -----

!   UWND     Zonal wind (dimensioned IDIM x JDIM x KDIM)
!   VWND     Meridional wind (dimensioned IDIM x JDIM x KDIM)
!   ATMP     Temperature at full model levels (IDIM x JDIM x KDIM)
!   PFULL    Pressure at full model levels (IDIM x JDIM x KDIM)
!   PHALF    Pressure at half model levels (IDIM x JDIM x KDIM+1)
!   ZFULL    Height at full model levels (IDIM x JDIM x KDIM)
!   ZHALF    Height at half model levels (IDIM x JDIM x KDIM+1)

    real,intent(in),dimension(:,:,:) :: uwnd,vwnd,atmp
    real,intent(in),dimension(:,:,:) :: pfull,phalf,zfull,zhalf
    real,intent(in),dimension(:,:  ) :: z_pbl                    

!   OUTPUT
!   ------

!   TAUX,TAUY    Base momentum flux in kg/m/s^2 (IDIM x JDIM) for diagnostics
!   DTAUX,DTAUY  Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
!   TAUS         normalized, "clipped" saturation momentum flux at 1/2 levels

    real,intent(out),dimension(:,:)   :: taux,tauy
    real,intent(out),dimension(:,:,:) :: dtaux,dtauy,taus

!---------------------------------------------------------------------

      call error_mesg('topo_drag', &
      'This module is not supported as part of the public release', FATAL)

  end subroutine topo_drag

  !=====================================================================

  subroutine topo_drag_init(lonb,latb,ierr)

    real,    intent(in), dimension(:) :: lonb,latb
    integer, intent(out)              :: ierr


!------- write version number and namelist ---------

    if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
    endif

    module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('topo_drag_init', &
      'This module is not supported as part of the public release', FATAL)

  end subroutine topo_drag_init

  !=====================================================================

  subroutine topo_drag_end

      module_is_initialized = .false.
 
!---------------------------------------------------------------------

      call error_mesg('topo_drag_end', &
      'This module is not supported as part of the public release', FATAL)

  end subroutine topo_drag_end

!#############################################################################      

end module topo_drag_mod
