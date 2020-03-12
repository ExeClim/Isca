module convective_cloud_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  use                fms_mod, only: stdlog, FATAL, WARNING, NOTE, error_mesg, &
                                    uppercase, check_nml_error
  use       time_manager_mod, only: time_type
  use     sat_vapor_pres_mod, only: compute_qs, lookup_es
  use       diag_manager_mod, only: register_diag_field, send_data
  use          constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, HLV, KAPPA, RADIUS, TFREEZE
  use                lcl_mod, only: lcl
  use marine_strat_cloud_mod, only: calc_lcls

  implicit none

  character(len=14), parameter :: mod_name = "conv_cloud"

  real :: conv_rain_min = 0.14 ! mm/day, threshold to produce conv cld
  real :: pshallow = 7.5e4     ! copy from GFDL/AM4 diag_cloud.F90
  integer :: id_conv_cf

  namelist /convective_cloud_nml/ conv_rain_min, pshallow

  contains

  subroutine convective_cloud_init(axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    integer :: io, ierr, nml_unit, stdlog_unit

#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=convective_cloud_nml, iostat=io)
    ierr = check_nml_error(io, 'convective_cloud_nml')
#else
    if (file_exist('input.nml')) then
      nml_unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
          read(nml_unit, nml=convective_cloud_nml, iostat=io, end=10)
          ierr = check_nml_error(io, 'convective_cloud_nml')
      enddo
10    call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, convective_cloud_nml)

    id_conv_cf = register_diag_field (mod_name, 'conv_cf', axes(1:3), Time, &
        'Convective cloud fraction for the simple cloud scheme', 'unitless: values 0-1')

  end subroutine convective_cloud_init

  subroutine convective_cloud_diag(pfull, precip, klcls, klzbs, conv_cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(in),  dimension(:,:)   :: precip
    integer, intent(in), dimension(:,:) :: klcls, klzbs
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: conv_cf
    real :: used

    call calc_convective_cf2(pfull, precip, klcls, klzbs, conv_cf, Time)

  end subroutine convective_cloud_diag

  subroutine calc_convective_cf(pfull, precip, klcls, conv_cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(in),  dimension(:,:)   :: precip
    integer, intent(in), dimension(:,:) :: klcls
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: conv_cf
    real, dimension(size(pfull,1), size(pfull,2)) :: precip_mm_per_day, plcl
    real    :: a, b, tower_scale_coeff, used
    integer :: k

    call calc_lcls(klcls, pfull=pfull, plcls=plcl)

    precip_mm_per_day = precip * 24.0 * 3600.0  ! change units to mm/day

    a = -0.125 * log(0.14)  !0.246
    b = 0.125
    tower_scale_coeff = 0.25

    conv_cf = 0.0

    do k=1, size(pfull,3)
      where (precip_mm_per_day<0.14)
        conv_cf(:,:,k) = 0.0
      elsewhere (precip_mm_per_day>85.0)
        conv_cf(:,:,k) = 0.8
      elsewhere
        conv_cf(:,:,k) = a + b * log(precip_mm_per_day)
      end where
      ! Convective clouds only exist when above LCL
      where(pfull(:,:,k) > plcl)
        conv_cf(:,:,k) = 0.0
      end where
    end do

    where (pfull < pshallow)
      conv_cf = conv_cf * tower_scale_coeff
    end where
    !!!!!!conv_cf = conv_cf * tower_scale_coeff

    ! Output the diagnostics
    if (id_conv_cf > 0) then
      used = send_data(id_conv_cf, conv_cf, Time)
    endif

  end subroutine calc_convective_cf

  subroutine calc_convective_cf2(pfull, precip, klcls, klzbs, conv_cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(in),  dimension(:,:)   :: precip
    integer, intent(in), dimension(:,:) :: klcls, klzbs
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: conv_cf
    real, dimension(size(pfull,1), size(pfull,2)) :: precip_mm_per_day, plcl2d, conv_cf_tmp
    real    :: convcld_a, convcld_b, used, tower_scale_coeff
    integer :: i, j, k, klcl, ktop, nlayers, k_tower

    call calc_lcls(klcls, pfull=pfull, plcls=plcl2d)

    ! conv_rain_min = 4 mm/day
    precip_mm_per_day = precip * 24.0 * 3600.0  ! change units to mm/day

    !convcld_a = -0.125 * log(0.14)  !0.246, ! 0.001
    convcld_a = -0.125 * log(conv_rain_min)
    convcld_b =  0.125 !0.0418 !
    tower_scale_coeff = 0.25

    conv_cf = 0.0
    !conv_cf_tmp = convcld_a + convcld_b * log(1.0 + precip_mm_per_day)

    where (precip_mm_per_day<conv_rain_min)
        conv_cf_tmp = 0.0
    !elsewhere(precip_mm_per_day>85.0)
    !    conv_cf_tmp = 0.8
    elsewhere
        conv_cf_tmp = convcld_a + convcld_b * log(precip_mm_per_day)
        conv_cf_tmp = min(0.8, max(0.0, conv_cf_tmp))
    end where

    !ktop = 3
    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        klcl = klcls(i,j)
        ktop = klzbs(i,j)
        !if (ktop.ne.0 .and. klcl.ne.ktop) then
        if (ktop.ne.0 .and. pfull(i,j,ktop)<4e4 .and. klcl.ne.ktop) then
          !nlayers = klcl - ktop
          ! Random overlap assumption
          ! conv_cf(i,j,ktop:klcl) = 1.0 - (1.0 - conv_cf_tmp(i,j))**(1.0 / nlayers)
          ! Maxmimum overlap
          !conv_cf(i,j,ktop:klcl) = conv_cf_tmp(i,j)*0.25
          ! Third try:
          conv_cf(i,j,ktop:klcl) = conv_cf_tmp(i,j) !*0.5 !*0.5
          !write(*,*) 'QL', i,j, klcl, ktop, nlayers, conv_cf(i,j,klcl), pfull(i,j,ktop)
          ! k_tower = minloc(abs(pfull(i,j,:) - pshallow), 1)
          ! if (k_tower>ktop .and. k_tower<klcl) then
          !   conv_cf(i,j,ktop:k_tower) = conv_cf_tmp(i,j) * tower_scale_coeff
          ! end if
        end if
      end do
    end do

    !where (pfull < pshallow)
    !  conv_cf = conv_cf * tower_scale_coeff
    !end where

    !conv_cf = min(0.6, max(0.0, conv_cf))

    ! Output the diagnostics
    if (id_conv_cf > 0) then
      used = send_data(id_conv_cf, conv_cf, Time)
    endif

  end subroutine calc_convective_cf2

  subroutine add_anvil_clouds(pfull, precip, klcls, klzbs, conv_cf, Time)
    real, intent(in),  dimension(:,:,:) :: pfull
    real, intent(in),  dimension(:,:)   :: precip
    integer, intent(in), dimension(:,:) :: klcls, klzbs
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:,:) :: conv_cf
    real, dimension(size(pfull,1), size(pfull,2)) :: precip_mm_per_day, plcl2d, conv_cf_tmp
    real    :: convcld_a, convcld_b, used, tower_scale_coeff
    integer :: i, j, k, klcl, ktop, nlayers, k_tower

    call calc_lcls(klcls, pfull=pfull, plcls=plcl2d)

    precip_mm_per_day = precip * 24.0 * 3600.0  ! change units to mm/day

    convcld_a = -0.125 * log(0.14)  !0.246, ! 0.001
    convcld_b =  0.125 !0.0418 !
    tower_scale_coeff = 0.25

    conv_cf = 0.0
    conv_cf_tmp = convcld_a + convcld_b * log(1.0 + precip_mm_per_day)

    !ktop = 3
    do i=1, size(pfull,1)
      do j=1, size(pfull,2)
        klcl = klcls(i,j)
        ktop = klzbs(i,j)
        if (ktop.ne.0 .and. klcl.ne.ktop) then
          !conv_cf(i,j,ktop:klcl) = conv_cf_tmp(i,j) !*0.5 !*0.5
          k_tower = minloc(abs(pfull(i,j,:) - 4.e4), 1)
          if (k_tower>ktop .and. k_tower<klcl) then
            !conv_cf(i,j,ktop:k_tower) = min(0.6, max(0.0, (conv_cf_tmp(i,j)-0.3)*2))
            conv_cf(i,j,ktop:k_tower) = conv_cf_tmp(i,j) * 0.25 !*0.5
            !write(*,*) 'QL anvil', conv_cf_tmp(i,j)*0.25
          end if
        end if
      end do
    end do

    !conv_cf = min(0.6, max(0.0, conv_cf))

    ! Output the diagnostics
    if (id_conv_cf > 0) then
      used = send_data(id_conv_cf, conv_cf, Time)
    endif

  end subroutine add_anvil_clouds

  subroutine convective_cloud_end()
    ! nothing
  end subroutine convective_cloud_end

end module convective_cloud_mod
