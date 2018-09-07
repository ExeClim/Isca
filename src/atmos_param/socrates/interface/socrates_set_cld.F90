! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates cloud type
!
!------------------------------------------------------------------------------
module socrates_set_cld
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_CLD'
contains

subroutine set_cld(cld, control, dimen, spectrum, n_profile, n_layer, &
  cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim, &
  dp_corr_strat, dp_corr_conv)

use def_cld,      only: StrCld, allocate_cld, allocate_cld_prsc
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use soc_constants_mod,   only: i_def, r_def
use rad_pcf,      only: &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat, ip_cloud_csiw, &
  ip_clcmp_st_water, ip_clcmp_st_ice, ip_clcmp_cnv_water, ip_clcmp_cnv_ice, &
  ip_phase_water, ip_phase_ice, ip_cloud_type_homogen, &
  ip_cloud_type_water, ip_cloud_type_ice, &
  ip_cloud_type_strat, ip_cloud_type_conv, &
  ip_cloud_type_sw, ip_cloud_type_si, ip_cloud_type_cw, ip_cloud_type_ci, &
  ip_drop_unparametrized, ip_ice_unparametrized, i_normal, i_err_fatal

use fms_mod, only: error_mesg, FATAL

implicit none


! Cloud properties:
type(StrCld),      intent(out) :: cld

! Control options:
type(StrCtrl),     intent(in)  :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

integer(i_def), intent(in) :: n_profile
integer(i_def), intent(in) :: n_layer

real(r_def), intent(in), optional :: &
  cloud_frac(:,:), conv_frac(:,:), &
  liq_frac(:,:), ice_frac(:,:), liq_conv_frac(:,:), ice_conv_frac(:,:), &
  liq_mmr(:,:), ice_mmr(:,:), liq_conv_mmr(:,:), ice_conv_mmr(:,:), &
  liq_dim(:,:), ice_dim(:,:), liq_conv_dim(:,:), ice_conv_dim(:,:)
!   Liquid and ice cloud fractions, gridbox mean mixing ratios, and
!   effective dimensions

real(r_def), intent(in), optional :: dp_corr_strat, dp_corr_conv
!   Decorrelation pressure scales for cloud vertical overlap


! Local variables
integer :: i, j, k, l
!   Loop variables
integer :: i_phase, i_param_type, n_cloud_parameter
!   Working variables
integer :: i_cloud_type(dimen%nd_cloud_component)
!   Types of cloud to which each component contributes

real(r_def) :: condensed_min_dim
real(r_def) :: condensed_max_dim
!   Minimum and maximum dimensions of each condensed component

real(r_def) :: eps = EPSILON(1.0)
real(r_def) :: min_cloud_fraction = 0.0001

integer                      :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_CLD'
character (len=128) :: cmessage

! Functions called
integer, external :: set_n_cloud_parameter


! Allocate structure for the core radiation code interface
call allocate_cld(cld, dimen, spectrum)
call allocate_cld_prsc(cld, dimen, spectrum)

if (.not.control%l_cloud) then
  return
end if

!------------------------------------------------------------------------------
! Set properties of condensed components
!------------------------------------------------------------------------------

if (control%l_ice .and. control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_st_ice
  case (ip_cloud_conv_strat, ip_cloud_csiw)
    cld%n_condensed = 4
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_st_ice
    cld%type_condensed(3) = ip_clcmp_cnv_water
    cld%type_condensed(4) = ip_clcmp_cnv_ice
  end select
else if (control%l_ice .and. .not.control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water)
    cld%n_condensed = 1
    cld%type_condensed(1) = ip_clcmp_st_ice
  case (ip_cloud_conv_strat, ip_cloud_csiw)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_ice
    cld%type_condensed(2) = ip_clcmp_cnv_ice
  end select
else if (.not.control%l_ice .and. control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water)
    cld%n_condensed = 1
    cld%type_condensed(1) = ip_clcmp_st_water
  case (ip_cloud_conv_strat, ip_cloud_csiw)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_cnv_water
  end select
else
  cmessage = 'Cloud on, but no condensed components included.'
  ierr=i_err_fatal
!   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  call error_mesg(ModuleName,cmessage, FATAL)
end if

do i=1, cld%n_condensed
  select case (cld%type_condensed(i))
  case (ip_clcmp_st_water)
    i_phase = ip_phase_water
    i_param_type = control%i_st_water
  case (ip_clcmp_st_ice)
    i_phase = ip_phase_ice
    i_param_type = control%i_st_ice
  case (ip_clcmp_cnv_water)
    i_phase = ip_phase_water
    i_param_type = control%i_cnv_water
  case (ip_clcmp_cnv_ice)
    i_phase = ip_phase_ice
    i_param_type = control%i_cnv_ice
  end select

  select case (i_phase)
  case (ip_phase_water)
    if (i_param_type <= 0) then
      cld%i_condensed_param(i) = ip_drop_unparametrized
      cmessage = 'Prescribed liquid cloud not yet implemented.'
      ierr=i_err_fatal
    !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      call error_mesg(ModuleName,cmessage, FATAL)      
    else if (i_param_type > spectrum%dim%nd_drop_type) then
      cmessage = 'Liquid cloud type outside allowed range.'
      ierr=i_err_fatal
    !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      call error_mesg(ModuleName,cmessage, FATAL)      
    else if (spectrum%drop%l_drop_type(i_param_type)) then
      ! Take parametrisation from spectral file
      cld%i_condensed_param(i) = spectrum%drop%i_drop_parm(i_param_type)
      cld%condensed_n_phf(i) = spectrum%drop%n_phf(i_param_type)
      ! DEPENDS ON: set_n_cloud_parameter
      n_cloud_parameter = set_n_cloud_parameter( cld%i_condensed_param(i), &
        cld%type_condensed(i), cld%condensed_n_phf(i) )
      do j=1, spectrum%basic%n_band
        do k=1, n_cloud_parameter
          cld%condensed_param_list(k, i, j) &
            = spectrum%drop%parm_list(k, j, i_param_type)
        end do
      end do
      ! Assign droplet mass mixing ratio and effective radius
      condensed_min_dim = spectrum%drop%parm_min_dim(i_param_type)
      condensed_max_dim = spectrum%drop%parm_max_dim(i_param_type)
      select case (cld%type_condensed(i))
      case (ip_clcmp_st_water)
        if (present(liq_mmr).and.present(liq_dim)) then
          do k = dimen%id_cloud_top, n_layer
            do l = 1, n_profile
              cld%condensed_mix_ratio(l, k, i) = liq_mmr(l, k)
              cld%condensed_dim_char(l, k, i) = min( max( liq_dim(l, k), &
                condensed_min_dim ), condensed_max_dim )
            end do
          end do
        else
          cmessage = 'Liquid MMR and effective radius not provided.'
          ierr=i_err_fatal
        !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
          call error_mesg(ModuleName,cmessage, FATAL)
        end if
      case (ip_clcmp_cnv_water)
        if (present(liq_conv_mmr).and.present(liq_conv_dim)) then
          do k = dimen%id_cloud_top, n_layer
            do l = 1, n_profile
              cld%condensed_mix_ratio(l, k, i) = liq_conv_mmr(l, k)
              cld%condensed_dim_char(l, k, i) = min( max( liq_conv_dim(l, k), &
                condensed_min_dim ), condensed_max_dim )
            end do
          end do
        else
          cmessage = 'Convective liquid MMR and effective radius not provided.'
          ierr=i_err_fatal
        !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
          call error_mesg(ModuleName,cmessage, FATAL)          
        end if        
      end select
    else
      cmessage = 'Liquid cloud type not in spectral file.'
      ierr=i_err_fatal
    !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      call error_mesg(ModuleName,cmessage, FATAL)      
    end if
  case (ip_phase_ice)
    if (i_param_type <= 0) then
      cld%i_condensed_param(i) = ip_ice_unparametrized
      cmessage = 'Prescribed ice cloud not yet implemented.'
      ierr=i_err_fatal
    !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      call error_mesg(ModuleName,cmessage, FATAL)      
    else if (i_param_type > spectrum%dim%nd_ice_type) then
      cmessage = 'Ice cloud type outside allowed range.'
      ierr=i_err_fatal
    !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      call error_mesg(ModuleName,cmessage, FATAL)      
    else if (spectrum%ice%l_ice_type(i_param_type)) then
      ! Take parametrisation from spectral file
      cld%i_condensed_param(i) = spectrum%ice%i_ice_parm(i_param_type)
      cld%condensed_n_phf(i) = spectrum%ice%n_phf(i_param_type)
      n_cloud_parameter = set_n_cloud_parameter( cld%i_condensed_param(i), &
        cld%type_condensed(i), cld%condensed_n_phf(i) )
      do j=1, spectrum%basic%n_band
        do k=1, n_cloud_parameter
          cld%condensed_param_list(k, i, j) &
            = spectrum%ice%parm_list(k, j, i_param_type)
        end do
      end do
      ! Assign ice mass mixing ratio and effective dimension
      condensed_min_dim = spectrum%ice%parm_min_dim(i_param_type)
      condensed_max_dim = spectrum%ice%parm_max_dim(i_param_type)
      select case (cld%type_condensed(i))
      case (ip_clcmp_st_ice)
        if (present(ice_mmr).and.present(ice_dim)) then
          do k = dimen%id_cloud_top, n_layer
            do l = 1, n_profile
              cld%condensed_mix_ratio(l, k, i) = ice_mmr(l, k)
              cld%condensed_dim_char(l, k, i) = min( max( ice_dim(l, k), &
                condensed_min_dim ), condensed_max_dim )
            end do
          end do
        else
          cmessage = 'Ice MMR and effective radius not provided.'
          ierr=i_err_fatal
        !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
          call error_mesg(ModuleName,cmessage, FATAL)          
        end if        
      case (ip_clcmp_cnv_ice)
        if (present(ice_conv_mmr).and.present(ice_conv_dim)) then
          do k = dimen%id_cloud_top, n_layer
            do l = 1, n_profile
              cld%condensed_mix_ratio(l, k, i) = ice_conv_mmr(l, k)
              cld%condensed_dim_char(l, k, i) = min( max( ice_conv_dim(l, k), &
                condensed_min_dim ), condensed_max_dim )
            end do
          end do
        else
          cmessage = 'Convective ice MMR and effective radius not provided.'
          ierr=i_err_fatal
        !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
          call error_mesg(ModuleName,cmessage, FATAL)          
        end if        
      end select
    else
      cmessage = 'Ice cloud type not in spectral file.'
      ierr=i_err_fatal
    !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      call error_mesg(ModuleName,cmessage, FATAL)      
    end if
  end select
end do

! Set the decorrelation scalings for cloud vertical overlap
if (present(dp_corr_strat)) then
  cld%dp_corr_strat = dp_corr_strat
else
  cld%dp_corr_strat = 0.0_r_def
end if
if (present(dp_corr_conv)) then
  cld%dp_corr_conv  = dp_corr_conv
else
  cld%dp_corr_conv = 0.0_r_def
end if


!------------------------------------------------------------------------------
! Set cloud amounts and convert mixing ratios to in-cloud values
!------------------------------------------------------------------------------

! Set cloud fractions
select case (control%i_cloud_representation)
case (ip_cloud_homogen)
  cld%n_cloud_type = 1
  do i = 1, cld%n_condensed
    i_cloud_type(i) = ip_cloud_type_homogen
  end do
  if (present(cloud_frac)) then
    do k = dimen%id_cloud_top, n_layer
      do l = 1, n_profile
        cld%frac_cloud(l, k, ip_cloud_type_homogen) = cloud_frac(l, k)
      end do
    end do      
  else
    cmessage = 'Cloud fraction not provided.'
    ierr=i_err_fatal
    ! CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    call error_mesg(ModuleName,cmessage, FATAL)    
  end if
case (ip_cloud_ice_water)
  cld%n_cloud_type = 2
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      i_cloud_type(i) = ip_cloud_type_water
    case (ip_clcmp_st_ice)
      i_cloud_type(i) = ip_cloud_type_ice
    end select
  end do
  if (present(liq_frac).and.present(ice_frac).and.present(cloud_frac)) then
    do k = dimen%id_cloud_top, n_layer
      do l = 1, n_profile
        if (liq_frac(l, k) + ice_frac(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          cld%frac_cloud(l, k, ip_cloud_type_water) = &
            cloud_frac(l, k)*liq_frac(l, k) / (liq_frac(l, k)+ice_frac(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_ice) = &
            cloud_frac(l, k)*ice_frac(l, k) / (liq_frac(l, k)+ice_frac(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_water) = 0.0_r_def
          cld%frac_cloud(l, k, ip_cloud_type_ice) = 0.0_r_def
        end if
      end do
    end do
  else if (present(liq_frac).and.present(ice_frac)) then
    do k = dimen%id_cloud_top, n_layer
      do l = 1, n_profile
        cld%frac_cloud(l, k, ip_cloud_type_water) = liq_frac(l, k)
        cld%frac_cloud(l, k, ip_cloud_type_ice) = ice_frac(l, k)
      end do
    end do
  else
    cmessage = 'Liquid and ice cloud fractions not provided.'
    ierr=i_err_fatal
    ! CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    call error_mesg(ModuleName,cmessage, FATAL)    
  end if
case (ip_cloud_conv_strat)
  cld%n_cloud_type = 2
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      i_cloud_type(i) = ip_cloud_type_strat
    case (ip_clcmp_st_ice)
      i_cloud_type(i) = ip_cloud_type_strat
    case (ip_clcmp_cnv_water)
      i_cloud_type(i) = ip_cloud_type_conv
    case (ip_clcmp_cnv_ice)
      i_cloud_type(i) = ip_cloud_type_conv
    end select
  end do
  if (present(cloud_frac).and.present(conv_frac)) then
  do k = dimen%id_cloud_top, n_layer
    do l = 1, n_profile
      cld%frac_cloud(l, k, ip_cloud_type_strat) = cloud_frac(l, k)
      cld%frac_cloud(l, k, ip_cloud_type_conv) = conv_frac(l, k)
    end do
    end do
  else
    cmessage = 'Cloud fraction and convective cloud fraction not provided.'
    ierr=i_err_fatal
    ! CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    call error_mesg(ModuleName,cmessage, FATAL)
  end if
case (ip_cloud_csiw)
  cld%n_cloud_type = 4
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      i_cloud_type(i) = ip_cloud_type_sw
    case (ip_clcmp_st_ice)
      i_cloud_type(i) = ip_cloud_type_si
    case (ip_clcmp_cnv_water)
      i_cloud_type(i) = ip_cloud_type_cw
    case (ip_clcmp_cnv_ice)
      i_cloud_type(i) = ip_cloud_type_ci
    end select
  end do
  if (present(liq_frac).and.present(ice_frac).and.present(cloud_frac)) then
    do k = dimen%id_cloud_top, n_layer
      do l = 1, n_profile
        if (liq_frac(l, k) + ice_frac(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          cld%frac_cloud(l, k, ip_cloud_type_sw) = &
            cloud_frac(l, k)*liq_frac(l, k) / (liq_frac(l, k)+ice_frac(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_si) = &
            cloud_frac(l, k)*ice_frac(l, k) / (liq_frac(l, k)+ice_frac(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_sw) = 0.0_r_def
          cld%frac_cloud(l, k, ip_cloud_type_si) = 0.0_r_def
        end if
      end do
    end do
  else if (present(liq_frac).and.present(ice_frac)) then
    do k = dimen%id_cloud_top, n_layer
      do l = 1, n_profile
        cld%frac_cloud(l, k, ip_cloud_type_sw) = liq_frac(l, k)
        cld%frac_cloud(l, k, ip_cloud_type_si) = ice_frac(l, k)
      end do
    end do
  else
    cmessage = 'Liquid and ice cloud fractions not provided.'
    ierr=i_err_fatal
    ! CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    call error_mesg(ModuleName,cmessage, FATAL)
  end if
  if (present(liq_conv_frac).and.present(ice_conv_frac).and. &
      present(conv_frac)) then
    do k = dimen%id_cloud_top, n_layer
      do l = 1, n_profile
        if (liq_conv_frac(l, k) + ice_conv_frac(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          cld%frac_cloud(l, k, ip_cloud_type_cw) = conv_frac(l, k) &
            *liq_conv_frac(l, k) / (liq_conv_frac(l, k)+ice_conv_frac(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_ci) = conv_frac(l, k) &
            *ice_conv_frac(l, k) / (liq_conv_frac(l, k)+ice_conv_frac(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_cw) = 0.0_r_def
          cld%frac_cloud(l, k, ip_cloud_type_ci) = 0.0_r_def
        end if
      end do
    end do
  else if (present(liq_conv_frac).and.present(ice_conv_frac)) then
    do k = dimen%id_cloud_top, n_layer
      do l = 1, n_profile
        cld%frac_cloud(l, k, ip_cloud_type_cw) = liq_conv_frac(l, k)
        cld%frac_cloud(l, k, ip_cloud_type_ci) = ice_conv_frac(l, k)
      end do
    end do
  else
    cmessage = 'Liquid and ice convective cloud fractions not provided.'
    ierr=i_err_fatal
    ! CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    call error_mesg(ModuleName,cmessage, FATAL)
  end if
end select

! Convert mass mixing ratios to in-cloud values
do i = 1, cld%n_condensed
  do k = dimen%id_cloud_top, n_layer
    do l = 1, n_profile
      cld%condensed_mix_ratio(l, k, i) = cld%condensed_mix_ratio(l, k, i) &
        / max(cld%frac_cloud(l, k, i_cloud_type(i)), eps)
    end do
  end do
end do

! Normalise the cloud fractions
do k = dimen%id_cloud_top, n_layer
  do l = 1, n_profile
    cld%w_cloud(l, k) = sum(cld%frac_cloud(l, k, 1:cld%n_cloud_type))
    if (cld%w_cloud(l, k) > min_cloud_fraction) then
      do j=1, cld%n_cloud_type
        cld%frac_cloud(l, k, j) = cld%frac_cloud(l, k, j) / cld%w_cloud(l, k)
        ! write(6,*) cld%frac_cloud(l, k, j) / cld%w_cloud(l, k), cld%frac_cloud(l, k, j) ,  cld%w_cloud(l, k), 'div?', cld%n_cloud_type
      end do
    else
      cld%w_cloud(l, k) = 0.0_r_def
      cld%frac_cloud(l, k, 1:cld%n_cloud_type) = 0.0_r_def
    end if
    if (cld%w_cloud(l, k) > 1.0_r_def + min_cloud_fraction) then
      cmessage = 'Cloud fraction greater than 1.'
      ierr=i_err_fatal
    !   CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      call error_mesg(ModuleName,cmessage, FATAL)
    else if (cld%w_cloud(l, k) > 1.0_r_def) then
      cld%w_cloud(l, k) = 1.0_r_def
    end if
  end do
end do

end subroutine set_cld
end module socrates_set_cld
