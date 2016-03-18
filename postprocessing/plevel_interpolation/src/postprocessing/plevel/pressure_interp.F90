
module pressure_interp_mod

!-----------------------------------------------------------------------

use moisture_convert_mod, only:  sphum_to_rh, rh_to_sphum
use   plev_constants_mod, only:  GRAV, RDGAS, RVGAS

implicit none
private

public pres_interp_type
public wind_interp, pres_interp, pres_interp_init
public pres_interp_free

type pres_interp_type
     integer          :: nx, ny, nz, kx
     integer, pointer :: index(:,:,:)
     real,    pointer :: factr(:,:,:)
     logical, pointer :: mask (:,:,:)
     logical          :: use_extrap
end type pres_interp_type

contains

!#######################################################################

   subroutine wind_interp ( pfull, u, v, pout, kbot, up, vp, mask )

     real, intent(in),  dimension(:,:,:) :: pfull, u, v
     real, intent(in),  dimension(:)     :: pout
  integer, intent(in),  dimension(:,:)   :: kbot
     real, intent(out), dimension(:,:,:) :: up, vp
  logical, intent(out), dimension(:,:,:), optional :: mask

  integer, dimension(size(pfull,1),size(pfull,2)) :: index
     real, dimension(size(pfull,1),size(pfull,2)) :: factr
     real, dimension(size(pfull,1),size(pfull,2),  &
                                   size(pfull,3)) :: log_pfull
     real :: log_pout
  integer :: i, j, k, n, idim, jdim

      idim = size(pfull,1);  jdim = size(pfull,2)

      log_pfull = log(pfull)

!     ----- initialize optional mask -----
!     (true for extrap below bottom level)

      if (present(mask)) mask = .true.

!  ---- loop through output levels ----

   do n = 1, size(pout)

      log_pout = log(pout(n))

!     ---- set up indexing (find full p-level just below pout) -----

      index = kbot

      do j = 1, jdim
      do i = 1, idim
         do k = 2, kbot(i,j)
            if (log_pout <= log_pfull(i,j,k)) then
                index(i,j) = k
                if (present(mask)) mask(i,j,n) = .false.
                exit
            endif
         enddo
      enddo
      enddo

    
!     ---- do interpolation ----
!     ---- (note: wind is held fixed for extrapolation) ----

      do j = 1, jdim
      do i = 1, idim
         k = index(i,j)
         factr(i,j) = (log_pout          -log_pfull(i,j,k)) /   &
                      (log_pfull(i,j,k-1)-log_pfull(i,j,k))
        !factr(i,j) = min(1.,factr(i,j)) ! no extrap above top level
        !factr(i,j) = max(0.,factr(i,j)) ! no extrap below bottom level
         factr(i,j) = min(+1.5,factr(i,j)) ! limit extrap above top level
         factr(i,j) = max(-0.5,factr(i,j)) ! limit extrap below bottom level
         up(i,j,n) = u(i,j,k) + factr(i,j) * (u(i,j,k-1) - u(i,j,k))
         vp(i,j,n) = v(i,j,k) + factr(i,j) * (v(i,j,k-1) - v(i,j,k))
      enddo
      enddo

   enddo

   end subroutine wind_interp

!#######################################################################

   subroutine pres_interp (Control, data, data_out)

      type(pres_interp_type), intent(in)  :: Control
                        real, intent(in)  :: data(:,:,:)
                        real, intent(out) :: data_out(:,:,:)

      integer  i,j,k,kd,ku,n

      do n = 1, Control%nz

         do j = 1, Control%ny
         do i = 1, Control%nx
             k  = Control % index(i,j,n)
             data_out(i,j,n) = data(i,j,k) + Control % factr(i,j,n) *  &
                              (data(i,j,k-1) - data(i,j,k))
         enddo
         enddo

      enddo

   end subroutine pres_interp

!#######################################################################

   function pres_interp_init ( pfull, pout, kbot,               &
                               tin, tout, zin, zout, qin, qout, &
                               use_extrap ) result ( Control )

      real, intent(in), dimension(:,:,:) :: pfull
      real, intent(in), dimension(:)     :: pout
   integer, intent(in), dimension(:,:)   :: kbot
      real, intent(in), dimension(:,:,:), optional :: tin, zin, qin
      real, intent(out),dimension(:,:,:), optional :: tout, zout, qout
   logical, intent(in),                   optional :: use_extrap
      type(pres_interp_type) :: Control

      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) ::  &
              log_pfull, tvin
      real, dimension(size(pfull,1),size(pfull,2)) :: tphalf
      real, dimension(size(pout,1)) :: log_pout
      real, dimension(size(pfull,1),size(pfull,2),size(pout)) :: tvout

      integer :: i, j, k, idim, jdim, n
      real    :: virtc, rdgrav2

!-----------------------------------------------------------------------

      idim = size(pfull,1);  jdim = size(pfull,2)

      log_pfull = log(pfull)
      log_pout  = log(pout)

!-----------------------------------------------------------------------
!------------------ setup control indexing -----------------------------

      Control % nx = idim
      Control % ny = jdim
      Control % nz = size(pout,1)
      Control % kx = size(pfull,3)

      Control % use_extrap = .true.
      if (present(use_extrap)) Control % use_extrap = use_extrap

      allocate (Control%index (Control%nx,Control%ny,Control%nz))
      allocate (Control%factr (Control%nx,Control%ny,Control%nz))
      allocate (Control%mask  (Control%nx,Control%ny,Control%nz))

      do n = 1, Control%nz

         Control % index(:,:,n) = kbot(:,:)
         Control % mask (:,:,n) = .true.

         do j = 1, jdim
         do i = 1, idim
             do k = 2, kbot(i,j)
                if (log_pout(n) <= log_pfull(i,j,k)) then
                    Control % index (i,j,n) = k
                    Control % mask  (i,j,n) = .false.
                    exit
                endif
             enddo
         enddo
         enddo

         do j = 1, jdim
         do i = 1, idim
             k = Control % index(i,j,n)
             Control % factr(i,j,n) =  &
                               (log_pout(n) - log_pfull(i,j,k)) / &
                        (log_pfull(i,j,k-1) - log_pfull(i,j,k))
            !Control % factr(i,j,n) = min(1.,Control%factr(i,j,n)) ! no extrap above top level
            !Control % factr(i,j,n) = max(0.,Control%factr(i,j,n)) ! no extrap below bottom level
             Control % factr(i,j,n) = min(+1.5,Control%factr(i,j,n)) ! limit extrap above top level
             Control % factr(i,j,n) = max(-0.5,Control%factr(i,j,n)) ! limit extrap below bottom level
         enddo
         enddo 

      enddo 

!-----------------------------------------------------------------------
    !--- temperature ---
      if (present(tin) .and. present(tout)) then
          call pres_interp (Control, tin, tout)
          if (Control%use_extrap) &
              call temp_extrap (Control, tin, log_pout, log_pfull, tout)
      endif
    !--- specific humidity ---
      if (present(qin) .and. present(qout)) then
          call pres_interp (Control, qin, qout)
          if (Control%use_extrap .and. present(tin) .and. present(tout)) &
              call q_extrap (Control, pfull, tin, qin, pout, tout, qout)
      endif

     !--- geopotential height ---
      if (present(zin) .and. present(zout) .and. &
          present(tin) .and. present(tout)) then
         !--- with moisture effects ---
          if (present(qin) .and. present(qout)) then
              virtc = (RVGAS-RDGAS)/RDGAS
              tvin  = tin  * (1.+virtc*qin)
              tvout = tout * (1.+virtc*qout)
          else
              tvin  = tin
              tvout = tout
          endif

          rdgrav2 = 0.5*RDGAS/GRAV
         ! interpolate between levels using hydrostatic equation
          do n = 1, Control%nz
          do j = 1, jdim
          do i = 1, idim
             k = Control % index(i,j,n)
             zout(i,j,n) = (log_pfull(i,j,k)-log_pout(n)) * &
                           (tvout(i,j,n)+tvin(i,j,k))*rdgrav2 + &
                           zin(i,j,k)
          enddo
          enddo
          enddo
      endif

!-----------------------------------------------------------------------

   end function pres_interp_init

!#######################################################################

   subroutine temp_extrap (Control, tin, log_pout, log_pfull, tout)

      type(pres_interp_type), intent(in) :: Control
      real, intent(in)    :: tin(:,:,:), log_pout(:), log_pfull(:,:,:)
      real, intent(inout) :: tout(:,:,:)

      integer :: i,j,k,n
      real :: rlapse, rrlaps, rglp21

      rlapse = 6.5e-3
      rrlaps = 1./rlapse

      do n = 1, Control%nz
      do j = 1, Control%ny
      do i = 1, Control%nx
         if ( Control % mask(i,j,n) ) then
            k = Control % index(i,j,n) - 1
            rglp21 = 0.5*(RDGAS/GRAV) * (log_pout(n)-log_pfull(i,j,k))
            tout(i,j,n) = tin(i,j,k) * (rrlaps+rglp21)/(rrlaps-rglp21)
         endif
      enddo
      enddo
      enddo

   end subroutine temp_extrap

!#######################################################################

   subroutine q_extrap (Control, pfull, tin, qin, pout, tout, qout)

      type(pres_interp_type), intent(in) :: Control
      real, intent(in)   , dimension(:,:,:) :: pfull, tin, qin, tout
      real, intent(in)   , dimension(:)     :: pout
      real, intent(inout), dimension(:,:,:) :: qout

      integer :: i, j, n, kb
      real    :: rh


      do n = 1, Control%nz
      do j = 1, Control%ny
      do i = 1, Control%nx
         if ( Control % mask(i,j,n) ) then
            kb = Control % index(i,j,n) - 1
            rh = sphum_to_rh ( qin(i,j,kb), tin(i,j,kb), pfull(i,j,kb) )
            qout(i,j,n) = rh_to_sphum ( rh, tout(i,j,n), pout(n) )
         endif
      enddo
      enddo
      enddo

   end subroutine q_extrap

!#######################################################################

   subroutine pres_interp_free (Control)
      type(pres_interp_type), intent(inout) :: Control

     deallocate (Control%index, Control%factr, Control%mask)
     Control%nx = 0; Control%ny = 0; Control%nz = 0; Control%kx = 0

   end subroutine pres_interp_free

!#######################################################################

end module pressure_interp_mod

