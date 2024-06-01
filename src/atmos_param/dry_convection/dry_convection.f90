module dry_convection_mod
  
  use           fms_mod, only: error_mesg, FATAL, stdlog, mpp_pe,             &
                               mpp_root_pe, open_namelist_file, close_file,   &
                               check_nml_error
  use     constants_mod, only: rdgas, cp_air, PSTD_MKS
  use  time_manager_mod, only: time_type
  use  diag_manager_mod, only: register_diag_field, send_data

!-----------------------------------------------------------------------------
  
  implicit none
  private

!                             ---  namelist ---
  real    :: small = 0.001
  real    :: tau = 120
  namelist /dry_convection_nml/ small,tau

 !                        ---  global variables --- 
  integer :: id_dt
  character(len=14), parameter :: mod_name='dry_convection'
  real :: missing_value = -1.e-10
  real :: Kappa
  public :: dry_convection_init, dry_convection

  contains

!-----------------------------------------------------------------------------

    subroutine dry_convection_init(axes, Time)

!                         ---  input arguments ---

      integer, intent(in) :: axes(4)
      type(time_type), intent(in) :: Time

!                         ---  local variables ---

      integer :: unit, ierr, io

!                         ---  executable code ---


      unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
         read(unit, nml=dry_convection_nml, iostat=io, end=20)
         ierr = check_nml_error (io, 'dry_convection_nml')
      enddo

20    call close_file (unit)

      if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=dry_convection_nml)

      id_dt = register_diag_field ( mod_name, 'dt_tg', axes(1:3), &
           Time, 'Temperature tendency', 'K/s', missing_value=missing_value)
      
      Kappa = rdgas/cp_air

    end subroutine dry_convection_init

!-----------------------------------------------------------------------------
  
  subroutine dry_convection(Time, temp0, p_full, p_half, dt_tg, cape, cin, lzb, lcl)

    ! DKOLL: do dry adjustment while conserving dry enthalpy
    ! For this to be exact, pressure level spacing needs to be fixed!
    ! --> does not work for rapidly changing surface pressures??
    ! Require:
    ! 1) cp/g [ TdP(k+1)+TdP(k) ](t=n+1) == cp/g [ TdP(k+1)+TdP(k) ](t=n)
    ! 2) theta(k+1,t=n+1) == theta(k,t=n+1)
    ! ==> (T*P**kappa) (k+1,t=n+1) == (T*P**kappa) (k,t=n+1)

!                         ---  input arguments ---
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(:,:,:) ::                                     &
         temp0,               &   ! temperature
         p_full,              &   ! pressure on full model levels
         p_half                   ! pressure on half model levels
    real, intent(out), dimension(size(temp0,1),size(temp0,2)) ::                    &
         cape,              &   !< convectively available potential energy
         cin                    !< convective inhibition
    integer, intent(out), dimension(:,:)   ::                                 &
         lcl,               &   !< lifting condensation level (index)
         lzb                    !< level of zero buoyancy     

!                         ---  output arguments ---
    real, intent(out), dimension(:,:,:) ::                                    &
         dt_tg                 ! temperature tendency

!                         ---  local variables ---
    real, dimension(size(temp0,1),size(temp0,2),size(temp0,3)) ::     &
         temp, pi, theta, dpres
    
    real,    dimension(size(temp0,1),size(temp0,2)) :: store
    logical, dimension(size(temp0,1),size(temp0,2)) :: do_daa
    
    integer :: kmax, iter, k
    logical :: do_any, did_adj

    logical :: used

!                         ---  executable code ---

! *********************************************************
    ! --- Set dimensions
    kmax  = size( temp0, 3 )

    ! --- Compute pressure thickness of layers
    do k=1, kmax
       dpres(:,:,k) = p_half(:,:,k+1) - p_half(:,:,k)
    end do

    ! --- Copy input temperature
    temp = temp0
    
    ! --- Compute exner function
    pi = ( p_full / Pstd_mks ) ** Kappa                                   

    ! --- Compute potential temperature
    theta = temp / pi                  

    did_adj = .false.
    
    do k = kmax - 1, 1, -1
       ! ----------------------------------------------
       
       ! --- Flag layers needing adjustment
       do_daa(:,:) = ( theta(:,:,k+1) - theta(:,:,k) ) > small
       do_any = ANY( do_daa(:,:) )
       
       ! --- Do adjustment
       if ( do_any ) then
          where ( do_daa )
             temp(:,:,k+1) = ( temp(:,:,k+1)*dpres(:,:,k+1) + temp(:,:,k)*dpres(:,:,k) )/ &
                  ( dpres(:,:,k+1) + (p_full(:,:,k)/p_full(:,:,k+1)) ** Kappa &
                  * dpres(:,:,k) )
             temp(:,:,k  ) = temp(:,:,k+1) * (p_full(:,:,k)/p_full(:,:,k+1)) ** Kappa
          end where
          did_adj = .true.
          ! print out message?
       end if
       
       ! ----------------------------------------------
    end do
           
    ! --- Compute change in temperature
    ! --- [Note: this is absolute change, need to divide by sec
    ! --- for dt_tg in atmosphere.f90]
    dt_tg = (temp - temp0)/tau

    if (id_dt   > 0) used = send_data ( id_dt, dt_tg, Time, 1, 1)

  end subroutine dry_convection

end module dry_convection_mod
