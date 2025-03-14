MODULE DRY_ADJ_MOD

    !=======================================================================
    !          DRY ADIABATIC ADJUSTMENT       
    !=======================================================================
    


     use       fms_mod, ONLY: file_exist, error_mesg, open_file, &
                              check_nml_error, mpp_pe, FATAL,    & 
                              NOTE, close_file, uppercase 

     use constants_mod, ONLY: grav, kappa
    !---------------------------------------------------------------------
     implicit none
     private
    
     public :: dry_adj, dry_adj_init, dry_adj_end, dry_adj_bdgt
    
    !---------------------------------------------------------------------
    
     character(len=128) :: version = '$Id: dry_adj.F90,v 17.0.4.1 2010/08/30 20:33:34 wfc Exp $'
     character(len=128) :: tag = '$Name: testing $'
     logical            :: module_is_initialized = .false.
    
    !---------------------------------------------------------------------
    
      real,    parameter :: p00     = 1000.0E2
      real               :: kappa_use = 0.0 ! modified kappa 
                                            ! (kappa_use = kappa * alpha) 
                                            ! initialised in dry_adj_init
    
    !---------------------------------------------------------------------
    ! --- NAMELIST
    !---------------------------------------------------------------------
    !     itermax - Max number of iterations
    !---------------------------------------------------------------------
    
      integer :: itermax = 5
      real    :: small = 0.001
      logical :: do_mcm_dry_adj = .true.
      logical :: do_relax = .false.
      real    :: tau_relax = 7200. 
      real    :: alpha = 1.0 ! values less than one achieve stable column. 
        
      namelist /dry_adj_nml/ itermax, small, do_mcm_dry_adj, do_relax, tau_relax!, alpha
    
    !---------------------------------------------------------------------
    
      contains
    
    !#######################################################################
    !#######################################################################
    
      SUBROUTINE DRY_ADJ ( temp0, pres, pres_int, dtemp, timestep, mask )
    
    !=======================================================================
    !  DRY ADIABATIC ADJUSTMENT
    !=======================================================================
    !---------------------------------------------------------------------
    ! Arguments (Intent in)
    !     temp0    - Temperature
    !     pres     - Pressure
    !     pres_int - Pressure at layer interface
    !     mask     -  OPTIONAL; floating point mask (0. or 1.) designating 
    !                 where data is present
    !---------------------------------------------------------------------
      real, intent(in), dimension(:,:,:) :: temp0, pres, pres_int
      real, intent(in) :: timestep 
      real, intent(in), OPTIONAL, dimension(:,:,:) :: mask
    
    !---------------------------------------------------------------------
    ! Arguments (Intent out)
    !     dtemp - Change in temperature
    !---------------------------------------------------------------------
      real, intent(out), dimension(:,:,:) :: dtemp
    
    !---------------------------------------------------------------------
    !  (Intent local)
    !---------------------------------------------------------------------
     
      real, dimension(size(temp0,1),size(temp0,2),size(temp0,3)) ::     &
            temp, pi, theta, pixdp, dpres, ppp
     
      real,    dimension(size(temp0,1),size(temp0,2)) :: store, xxx
      logical, dimension(size(temp0,1),size(temp0,2)) :: do_daa
     
      integer :: kmax, iter, k, i, j
      logical :: do_any, did_adj
    
    !====================================================================
    
    ! --- Check to see if dry_adj has been initialized
      if(.not. module_is_initialized ) CALL error_mesg( 'DRY_ADJ', &
                               'dry_adj_init has not been called', FATAL )
    
    ! --- Set dimensions
      kmax  = size( temp0, 3 )
    
    ! --- Compute pressure thickness of layers
      dpres(:,:,1:kmax) = pres_int(:,:,2:kmax+1) - pres_int(:,:,1:kmax)
    
    ! --- Copy input temperature
      temp = temp0
    
    ! --- Compute exner function
    !do i = 1,size(temp0,1)
    !   do j = 1,size(temp0,2)
    !      pi(i,j,:) = (pres(i,j,:) / pres(i,j,size(temp0,3))) ** kappa_use
    !   end do
    !end do
    ! --- Compute exner function
      pi = ( pres / p00 ) ** kappa_use                                   
    
    ! --- Compute product of pi and dpres
      pixdp = pi * dpres
    
    ! --- Compute potential temperature
      theta = temp / pi                  
    
      if(do_mcm_dry_adj) then
        do k = 2,kmax
          xxx = 0.5*kappa_use*(pres(:,:,k) - pres(:,:,k-1))/pres_int(:,:,k)
          ppp(:,:,k) = (1.0 + xxx)/(1.0 - xxx)
        enddo
      endif
        
    !-----------------------------------------------------------------
    ! iteration loop starts           
    !-----------------------------------------------------------------
      do iter = 1,itermax
    !-----------------------------------------------------------------           
    
      did_adj = .false.
    
      do k = 1,kmax - 1
    ! ----------------------------------------------
    
    ! --- Flag layers needing adjustment
      if(do_mcm_dry_adj) then
        do_daa(:,:) = temp(:,:,k+1) > ( temp(:,:,k)*ppp(:,:,k+1) + small )
      else
        do_daa(:,:) = ( theta(:,:,k+1) - theta(:,:,k) ) > small
      endif
      
      if( PRESENT( mask ) ) then
      do_daa(:,:) = do_daa(:,:) .and. ( mask(:,:,k+1) > 0.5 )
      endif
      do_any = ANY( do_daa(:,:) )
    
    ! --- Do adjustment
     if ( do_any ) then
       if(do_mcm_dry_adj) then
         where ( do_daa )
           temp(:,:,k) = (temp(:,:,k)  * dpres(:,:,k  )   + &
                          temp(:,:,k+1)* dpres(:,:,k+1) )   &
                       /(dpres(:,:,k) + ppp(:,:,k+1)*dpres(:,:,k+1))
           temp(:,:,k+1) = temp(:,:,k)*ppp(:,:,k+1)
         end where
         did_adj = .true.
       else
         where ( do_daa )
           store(:,:) = ( theta(:,:,k  ) * pixdp(:,:,k  )     &
                       +  theta(:,:,k+1) * pixdp(:,:,k+1) )   &
                      / ( pixdp(:,:,k  ) + pixdp(:,:,k+1) )
           theta(:,:,k  ) = store(:,:)
           theta(:,:,k+1) = store(:,:)
            temp(:,:,k  ) = pi(:,:,k  ) * theta(:,:,k  ) 
            temp(:,:,k+1) = pi(:,:,k+1) * theta(:,:,k+1)
         end where
         did_adj = .true.
       endif
     end if
    
    ! ----------------------------------------------
      end do
    
    ! --- If no adjusment made this pass, exit iteration loop.
      if ( .not. did_adj ) go to 900
    
    !-----------------------------------------------------------------
      end do
    !-----------------------------------------------------------------
    ! iteration loop ends           
    !-----------------------------------------------------------------
     if(.not.do_mcm_dry_adj) then
         call error_mesg ('DRY_ADJ', 'Non-convergence in dry_adj', NOTE)
     endif
     900 continue
        
    ! --- Compute change in temperature
     if (do_relax) then 
        dtemp = (temp - temp0) * (1 - exp(-timestep / tau_relax))
     else
        dtemp = temp - temp0
     endif 
        
    !=======================================================================
      end SUBROUTINE DRY_ADJ
    
    !#####################################################################
    !#####################################################################
    
      SUBROUTINE DRY_ADJ_INIT(alpha_in)
    
    !=======================================================================
    ! ***** INITIALIZE RAS
    !=======================================================================
    
      real, intent(in) :: alpha_in
    !---------------------------------------------------------------------
    !  (Intent local)
    !---------------------------------------------------------------------
    
      integer :: unit, io, ierr
    
    !=====================================================================
    
    !---------------------------------------------------------------------
    ! --- READ NAMELIST
    !---------------------------------------------------------------------
    
      if (file_exist('input.nml')) then
        unit = open_file (file='input.nml', action='read')
        ierr = 1
        do while (ierr /= 0)
           read  (unit, nml=dry_adj_nml, iostat=io, end=10)
           ierr = check_nml_error(io, 'dry_adj_nml')
        end do
 10     call close_file(unit)
     endif
    
    !------- write version number and namelist ---------
    
     unit = open_file (file='logfile.out', action='append')
     if ( mpp_pe() == 0 ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit,nml=dry_adj_nml)
     endif
     call close_file(unit)
    
    !-------------------------------------------------------------------

     alpha = alpha_in

     kappa_use = kappa * alpha 
    
      module_is_initialized = .true.
    
    !=====================================================================
      end SUBROUTINE DRY_ADJ_INIT
    
    
    !#######################################################################
    !#######################################################################
    
      SUBROUTINE DRY_ADJ_END()
    
    !-------------------------------------------------------------------
    
      module_is_initialized = .false.
    
    !=====================================================================
      end SUBROUTINE DRY_ADJ_END
    
    
    !#######################################################################
    !#######################################################################
    
      SUBROUTINE DRY_ADJ_BDGT ( dtemp, pres_int )
    
    !=======================================================================
    ! Budget check for dry adiabatic adjustment - a debugging tool
    !=======================================================================
    
    !---------------------------------------------------------------------
    ! Arguments (Intent in)
    !     dtemp    - Temperature change 
    !     pres_int - Pressure at layer interface
    !---------------------------------------------------------------------
      real, intent(in), dimension(:,:,:) :: dtemp, pres_int
    
    !---------------------------------------------------------------------
    !  (Intent local)
    !---------------------------------------------------------------------
    
     real, dimension(size(dtemp,1),size(dtemp,2),size(dtemp,3)) ::  dpres
     real    :: sum_dtemp
     integer :: imax, jmax, kmax, i, j, k
    
    !=======================================================================
    
      imax = size ( dtemp, 1 )
      jmax = size ( dtemp, 2 )
      kmax = size ( dtemp, 3 )
    
    ! --- Compute pressure thickness of layers
      dpres(:,:,1:kmax) = pres_int(:,:,2:kmax+1) - pres_int(:,:,1:kmax)
    
    ! --- Check budget
    
      do j = 1,jmax
      do i = 1,imax
    
        sum_dtemp = 0.                                                          
    
      do k = 1,kmax
        sum_dtemp = sum_dtemp + dtemp(i,j,k)*dpres(i,j,k) / Grav                                   
      end do
    
      if ( abs( sum_dtemp ) > 1.0e-4 ) then
        print *
        print *, ' DRY ADIABATIC ADJUSTMENT BUDGET CHECK AT i,j = ', i,j
        print *, ' sum_dtemp  = ',  sum_dtemp                                                                  
        print *, 'STOP'
    !    STOP
      endif
    
      end do
      end do
    
    !=======================================================================
      end SUBROUTINE DRY_ADJ_BDGT
    
    !#######################################################################
    !#######################################################################
      end MODULE DRY_ADJ_MOD