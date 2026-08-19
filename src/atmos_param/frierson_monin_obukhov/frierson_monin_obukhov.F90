
module frierson_monin_obukhov_mod

    !======================================================================
    !
    !                         MONIN-OBUKHOV MODULE
    !
    !          Routines for computing surface drag coefficients from
    !               data at the lowest model level using
    !                          Monin-Obukhov scaling 
    !
    !======================================================================
    
    
    use constants_mod, only : grav, vonkarm
    
    use       fms_mod, only : file_exist, check_nml_error,              &
                              open_namelist_file, write_version_number, &
                              mpp_pe, mpp_root_pe, stdlog, close_file,  &
                              error_mesg, FATAL
    
    implicit none
    private
    
    ! public interfaces
    !=======================================================================
     public frierson_mo_drag
     public mo_profile
     public mo_diff
     public stable_mix
    !=======================================================================
    
    
    ! form of interfaces
    !=======================================================================
    ! call  mo_drag (pt, pt0, z, z0, zt, zq, speed, &
    !                drag_m, drag_t, drag_q,            &
    !                u_star, b_star, [mask])
    !
    !   (In the following the phrase "dimension(:) or (:,:)" means
    !      that this routine can be called either with all of the 
    !      variables so designated being either 1d or 2d arrays.
    !      All of these arrays should conform exactly.)
    !
    !   input: 
    !
    !      pt,  real, dimension(:) or (:,:)
    !         virtual potential temperature at lowest model level
    !         degrees Kelvin
    !
    !      pt0, real, dimension(:) or (:,:)
    !         virtual potential temperature at surface
    !         degrees Kelvin
    !
    !      z, real, dimension(:) or (:,:) 
    !         height above surface of lowest model layer 
    !         meters
    !
    !      z0, real, dimension(:) or (:,:)
    !          surface roughness for momentum 
    !          meters
    !
    !      zt, real, dimension(:) or (:,:)
    !          surface roughness for temperature
    !          meters
    !
    !      zq, real, dimension(:) or (:,:)
    !          surface roughness for moisture
    !          meters
    !
    !      speed, real, dimension(:) or (:,:) 
    !          wind speed at lowest model level with respect to surface 
    !             (any "gustiness" factor should be included in speed)
    !          meters/sec
    !
    !   inout:
    !
    !        drag_m, real, dimension(:) or (:,:)
    !              non-dimensional drag coefficient for momentum
    !           
    !        drag_t, real, dimension(:) or (:,:)
    !             non-dimensional drag coefficient for temperature
    !           
    !        drag_q, real, dimension(:) or (:,:)
    !             non-dimensional drag coefficient for moisture
    !
    !          (the input values are used only if the time-smoothing 
    !            option is turned on)
    !
    !   output:
    !
    !        u_star, real, dimension(:) or (:,:)
    !           friction velocity 
    !           meters/sec
    !
    !        b_star, real, dimension(:) or (:,:)
    !           buoyancy scale 
    !           (meters/sec)**2
    !
    !
    !            The magnitude of the wind stress is 
    !                 density*(ustar**2)
    !            The drag coefficient for momentum is 
    !                 (u_star/speed)**2
    !            The buoyancy flux is
    !                 density*ustar*bstar
    !            The drag coefficient for heat etc is
    !                 (u_star/speed)*(b_star/delta_b)
    !                 where delta_b is the buoyancy difference between
    !                  the surface and the lowest model level
    !
    !     
    !      
    !    optional:
    !       mask    : logical, dimension(:) 
    !                   computation performed only where mask = .true.
    !       NOTE(!) :  mask option is only available for 1d verison 
    !
    !==========================================================================
    !
    ! subroutine mo_profile(zref, z, z0, zt, zq, u_star, b_star, q_star, &
    !                          del_m, del_h, del_q, [mask])
    !
    !     (In the following the phrase "dimension(:) or (:,:)" means
    !      that this routine can be called either with all of the 
    !      variables so designated being either 1d or 2d arrays.
    !      All of these arrays should conform exactly.)
    !
    !  input:
    !
    !     zref, real
    !           height above surface to which interpolation is requested
    !           meters
    !
    !     z, real, dimension(:) or (:,:) 
    !        height of lowest model layer
    !        meters
    !
    !     z0, real, dimension(:) or (:,:)
    !         surface roughness for momentum 
    !         meters
    !
    !     zt, real, dimension(:) or (:,:)
    !         surface roughness for temperature
    !         meters
    !
    !     zq, real, dimension(:) or (:,:)
    !         surface roughness for moisture
    !         meters
    !
    !     u_star, real, dimension(:) or (:,:)
    !             friction velocity 
    !             meters/sec
    !
    !     b_star, real, dimension(:) or (:,:)
    !             buoyancy scale
    !             (meters/sec)**2
    !
    !     q_star, real, dimension(:) or (:,:)
    !             moisture scale
    !             kg/kg
    !
    !          (Note:  u_star and b_star are output from mo_drag,
    !                  q_star = flux_q/u_star/rho )
    !
    !    optional input:
    !
    !       mask, logical, dimension(:) 
    !                   computation performed only where mask = .true.
    !       NOTE(!):  mask option is only available for 1d verison 
    !
    !
    !    output:
    !
    !       del_m, real, dimension(:) or (:,:)
    !              dimensionless ratio, as defined below, for momentum
    !
    !       del_h, real, dimension(:) or (:,:)
    !              dimensionless ratio, as defined below, for temperature
    !
    !       del_q, real, dimension(:) or (:,:)
    !              dimensionless ratio, as defined below, for moisture
    !
    !          Ratios are  (f(zref) - f_surf)/(f(z) - f_surf)
    !
    !==========================================================================
    !               
    !  subroutine mo_diff(z, u_star, b_star, k_m, k_h, [mask])
    !
    !    input: 
    !        z, real, dimension(see below)
    !                height above surface of at which diffusivities 
    !                are desired
    !                meters
    !
    !        u_star, real, dimension(see below)
    !                surface friction velocity
    !                meters/sec
    !
    !        b_star, real, dimension(see below)
    !                buoyancy scale 
    !                (meters/sec)**2
    !
    !          (Note:  u_star and b_star are output from mo_drag)
    !
    !    optional input:
    !       mask    : logical, dimension(:) 
    !                   computation performed only where mask = .true.
    !       NOTE:  mask option is only available for 1d verisons -- see below
    !
    !    output:
    !
    !        k_m   : real, dimension(see below)
    !                kinematic diffusivity for momentum 
    !                (meters**2/sec)
    !
    !        k_h   : real, dimension(see below)
    !                kinemtatic diffusivity for temperature
    !                (meters**2/sec)
    !
    !    dimensions:  any of the following four options can be used
    !
    !      1) diffusivities desired on multiple levels, with 2d (x,y) input
    !              z(:,:,:), k_m(:,:,:), k_h(:,:,:)
    !              u_star(:,:), b_star(:,:) corresponding to the 1st and 2nd
    !                  indices of z -- vertical level is third index
    !              mask option NOT available
    !
    !      2) as in 1), but with only 1d input
    !              z(:,:), k_m(:,:), k_h(:,:)
    !              u_star(:), b_star(:) corresponding to the 1st index of z
    !              vertical level is second index
    !              mask option available
    !
    !      3) diffusivities desired on one level only, with 2d input
    !              z(:,:), k_m(:,:), k_h(:,:),u_star(:,:), b_star(:,:) 
    !              mask option NOT available
    !
    !      4) diffusivities desired on one level only, with 1d input
    !              z(:), k_m(:), k_h(:), u_star(:), b_star(:) 
    !              mask option available
    !
    !=======================================================================
    
    interface frierson_mo_drag
        module procedure  mo_drag_1d, mo_drag_2d
    end interface
    
    interface mo_profile
        module procedure  mo_profile_1d, mo_profile_2d
    end interface
    
    interface mo_diff
        module procedure  mo_diff_0d,         mo_diff_1d,         mo_diff_2d
        module procedure  mo_diff_one_lev_0d, mo_diff_one_lev_1d, mo_diff_one_lev_2d
    end interface
    
    interface stable_mix
        module procedure  stable_mix_0d, stable_mix_1d
        module procedure  stable_mix_2d, stable_mix_3d
    end interface
    
    
    !--------------------- version number ---------------------------------
    
    character(len=128) :: version = '$Id: frierson_monin_obukhov.f90,v 1.5.4.1 2005/05/13 18:16:40 pjp Exp $'
    character(len=128) :: tagname = '$Name:  $'
    
    !=======================================================================
    
    !  DEFAULT VALUES OF NAMELIST PARAMETERS:
    
    real    :: rich_crit  = 2.0
    real    :: drag_min   = 1.e-05
    real    :: relax_time = 0.            
    logical :: neutral    = .false.
    integer :: stable_option  = 1
    real    :: zeta_trans     = 0.5
    
    namelist /frierson_monin_obukhov_nml/ rich_crit, neutral, drag_min, relax_time, &
                                 stable_option, zeta_trans
    
    !=======================================================================
    
    !  MODULE VARIABLES
    
    real, parameter    :: small  = 1.e-04
    real               :: b_stab, r_crit, sqrt_drag_min, lambda, rich_trans
    logical            :: init = .false.
    
    
    contains
    
    !=======================================================================
    
    subroutine monin_obukhov_init
    
    integer :: unit, ierr, io
    
    !------------------- read namelist input -------------------------------
    
          if (file_exist('input.nml')) then
             unit = open_namelist_file ( )
             ierr=1; do while (ierr /= 0)
                read  (unit, nml=frierson_monin_obukhov_nml, iostat=io, end=10)
                ierr = check_nml_error(io,'frierson_monin_obukhov_nml')
             enddo
      10     call close_file (unit)
          endif
    
    !---------- output namelist to log-------------------------------------
    
          call write_version_number ( version, tagname )
          if ( mpp_pe() == mpp_root_pe() ) then
               write (stdlog(), nml=frierson_monin_obukhov_nml)
          endif
    
    if(rich_crit.le.0.25)  call error_mesg( &
            'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
            'rich_crit in monin_obukhov_mod must be > 0.25', FATAL)
    
    if(drag_min.le.0.0)  call error_mesg( &
            'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
            'drag_min in monin_obukhov_mod must be >= 0.0', FATAL)
    
    b_stab = 1.0/rich_crit
    r_crit = 0.95*rich_crit
    
    sqrt_drag_min = 0.0
    if(drag_min.ne.0.0) sqrt_drag_min = sqrt(drag_min)
    
    lambda     = 1.0 + (5.0 - b_stab)*zeta_trans   ! used only if stable_option = 2
    rich_trans = zeta_trans/(1.0 + 5.0*zeta_trans) ! used only if stable_option = 2
    
    init = .true.
    
    return
    end subroutine monin_obukhov_init
    
    !=======================================================================
    
    subroutine mo_drag_1d &
             (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star, mask)
    
    real, intent(in)   , dimension(:) :: pt, pt0, z, z0, zt, zq, speed
    real, intent(out),   dimension(:) :: drag_m, drag_t, drag_q
    real, intent(out),   dimension(:) :: u_star, b_star
    logical, intent(in), optional, dimension(:) :: mask
    
    real :: vk2, r_crit1
    
    real   , dimension(size(pt)) :: rich, fm, fh, fq
    integer, dimension(size(pt)) :: stabl, unstabl
    logical, dimension(size(pt)) :: avail
    
    real   , dimension(size(pt)) :: &
        rich_1d, z_1d, z0_1d, zt_1d, zq_1d, fm_1d, fh_1d, fq_1d, zeta_1d, delta_b, us, bs, qs, &
        u_star_eq, b_star_eq, drag_m_eq, drag_t_eq, drag_q_eq
    
    real :: xi, xi_1, xi_2
    
    integer :: i,npts,nptu
    
    if(.not.init) call monin_obukhov_init
    
    avail = .true.
    if (present(mask)) avail = mask
    
    vk2        = vonkarm**2
    
    ! delta_b is the difference in buoyancy b/w first level and surface
    ! rich is the richardson number: (dB/dz)/(du/dz)^2
    where(avail) 
       delta_b = grav*(pt0 - pt)/pt0
       rich    = - z*delta_b/(speed*speed + small)
    else where 
       rich = 0.0
    end where
    
    if(neutral) then
    
      where(avail)
        fm   = log(z/z0)
        fh   = log(z/zt)
        fq   = log(z/zq)
        us   = vonkarm/fm
        bs   = vonkarm/fh
        qs   = vonkarm/fq
        drag_m    = us*us
        drag_t    = us*bs
        drag_q    = us*qs
        u_star = us*speed
        b_star = bs*delta_b
      end where
    
    else
    
    ! large richardson number -- small velocity or large buoyancy
      where(avail .and. rich >= r_crit) 
        drag_m_eq   = drag_min
        drag_t_eq   = drag_min
        drag_q_eq   = drag_min
        us          = sqrt_drag_min
        bs          = sqrt_drag_min
        qs          = sqrt_drag_min
      end where
    
    ! solve for zeta on stable and unstable points separately
    
    ! this subroutine takes in the richardson numbers, and puts into the arrays
    ! stabl and unstabl (which have size npts and nptu) based on the sign of rich
      call separate_stabilities(npts, nptu, stabl, unstabl, rich, avail, r_crit)
    ! unstable points -- now treated as the neutral case.  
      if(nptu > 0) then
    
    !    do i =1, nptu
    !       rich_1d(i) = rich(unstabl(i))
    !       z_1d   (i) = z   (unstabl(i))
    !       z0_1d  (i) = z0  (unstabl(i))
    !       zt_1d  (i) = zt  (unstabl(i))
    !       zq_1d  (i) = zq  (unstabl(i))
    !    end do
    !
    !    call solve_zeta  &
    !     (nptu, rich_1d, z_1d, z0_1d, zt_1d, zq_1d, fm_1d, fh_1d, fq_1d, &
    !      unstable = .true.)
    !
        do i = 1, nptu
    !      fm(unstabl(i)) = fm_1d(i)
    !      fh(unstabl(i)) = fh_1d(i)
    !      fq(unstabl(i)) = fq_1d(i)
           fm(unstabl(i))   = log(z(unstabl(i))/z0(unstabl(i)))
           fh(unstabl(i))   = log(z(unstabl(i))/zt(unstabl(i)))
           fq(unstabl(i))   = log(z(unstabl(i))/zq(unstabl(i)))
        end do
      end if
    
    ! stable points -- now treated with phi = 1 + zeta (so drag coeff goes to 
    ! zero near Ri = 1).  
    ! no iteration is necessary -- we have an analytic solution
      if(npts > 0) then
    !    do i =1, npts
    !       rich_1d(i) = rich(stabl(i))
    !       z_1d   (i) = z   (stabl(i))
    !       z0_1d  (i) = z0  (stabl(i))
    !       zt_1d  (i) = zt  (stabl(i))
    !       zq_1d  (i) = zq  (stabl(i))
    !    end do
    !
    !    call solve_zeta  &
    !     (npts, rich_1d, z_1d, z0_1d, zt_1d, zq_1d, fm_1d, fh_1d, fq_1d, &
    !      unstable = .false.)
    !
        
        do i = 1, npts
          fm(stabl(i)) = log(z(stabl(i))/z0(stabl(i)))/(1. - & 
             rich(stabl(i))*(1. - z0(stabl(i))/z(stabl(i))))
          fh(stabl(i)) = log(z(stabl(i))/zt(stabl(i)))/(1. - & 
             rich(stabl(i))*(1. - zt(stabl(i))/z(stabl(i))))
          fq(stabl(i)) = log(z(stabl(i))/zq(stabl(i)))/(1. - & 
             rich(stabl(i))*(1. - zq(stabl(i))/z(stabl(i))))
    !      fm(stabl(i)) = fm_1d(i)
    !      fh(stabl(i)) = fh_1d(i)
    !      fq(stabl(i)) = fq_1d(i)
        end do
      end if
    
      where (avail .and. rich < r_crit)
        us   = max(vonkarm/fm, sqrt_drag_min)
        bs   = max(vonkarm/fh, sqrt_drag_min)
        qs   = max(vonkarm/fq, sqrt_drag_min)
        drag_m_eq   = us*us
        drag_t_eq   = us*bs
        drag_q_eq   = us*qs
      end where
    
      if(relax_time.ne.0.0) then
         call error_mesg('mo_drag', &
         'Non-zero relax_time not allowed in this version of monin_obukhov',FATAL)
    !    xi   = dt/relax_time
    !    xi_1 = 1.0/(1.0 + xi)
    !    xi_2 =  xi/(1.0 + xi)
    !    where (avail)
    !      drag_m = xi_1*drag_m + xi_2*drag_m_eq
    !      drag_t = xi_1*drag_t + xi_2*drag_t_eq
    !      drag_q = xi_1*drag_q + xi_2*drag_q_eq
    !      us = sqrt(drag_m)
    !      bs = drag_t/us
    !      u_star = us*speed
    !      b_star = bs*delta_b
    !   end where
      else
         where (avail)
           drag_m = drag_m_eq
           drag_t = drag_t_eq
           drag_q = drag_q_eq
           u_star = us*speed
           b_star = bs*delta_b
         end where
      endif
    
    end if
    
    return
    end subroutine mo_drag_1d
    
    !=======================================================================
    
    subroutine solve_zeta(n, rich, z, z0, zt, zq, f_m, f_h, f_q, unstable)
    
    integer, intent(in) :: n
    real   , intent(in), dimension(n)  :: rich, z, z0, zt, zq
    logical, intent (in) :: unstable
    real   , intent(out), dimension(n) :: f_m, f_h, f_q
    
    real :: vk, vk2, error, max_cor, zeta_min
    real, dimension(n) ::   &
              d_rich, rich_1, correction,  &
              z_z0, z_zt, z_zq, ln_z_z0, ln_z_zt, ln_z_zq, zeta, zeta_q, corr
    
    integer :: iter, max_iter,nn
    
    vk2 = vonkarm**2
    error = 1.e-4
    max_iter = 20
    zeta_min = 1.e-6
    
    z_z0 = z/z0
    z_zt = z/zt
    z_zq = z/zq
    ln_z_z0 = log(z_z0)
    ln_z_zt = log(z_zt)
    ln_z_zq = log(z_zq)
    
    zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
    if(.not.unstable)  zeta = zeta/(1.0 - rich/rich_crit)
    
    iter_loop: do iter = 1, max_iter
    
      where (abs(zeta).lt.zeta_min) zeta = zeta_min
      call solve_zeta_internal   &
           (unstable, zeta,z_z0,z_zt,ln_z_z0,ln_z_zt,f_m,f_h,rich_1,d_rich)
      correction = (rich - rich_1)/d_rich  
      corr = min(abs(correction),abs(correction/zeta))
      max_cor= maxval(corr)
    
      if(max_cor .gt. error) then
        where(corr.gt. error) zeta = zeta + correction
        cycle iter_loop
      else
        zeta_q = zeta/z_zq
        if(unstable) then
          call mo_integral_uh(f_q, zeta, zeta_q, ln_z_zq)
        else
          call mo_integral_s (f_q, zeta, zeta_q, ln_z_zq)
        end if
        return
      end if
    
    end do iter_loop
    
    call error_mesg ('solve_zeta in monin_obukhov_mod',  &
                     'no convergence in surface drag iteration', FATAL)
    
    end subroutine solve_zeta
    
    !=======================================================================
    
    subroutine solve_zeta_internal(unstable, zeta, z_z0, z_zt, &
            ln_z_z0, ln_z_zt, f_m, f_h, rich, d_rich)
    
    logical, intent(in) :: unstable
    real, intent(in), dimension(:) :: zeta, z_z0, z_zt, ln_z_z0, ln_z_zt
    real, intent(out), dimension(:) :: f_m, f_h, rich, d_rich
    
    real, dimension(size(zeta)) ::   &
       x, x_0, x_h, x_h_0, rzeta, zeta_0, zeta_h, df_m, df_h
    
    rzeta   = 1.0/zeta
    zeta_0 = zeta/z_z0
    zeta_h = zeta/z_zt
    
    if(unstable) then
      call mo_derivative_um(x  ,  zeta  )
      call mo_derivative_um(x_0,  zeta_0)
    else
      call mo_derivative_s(x  ,  zeta  )
      call mo_derivative_s(x_0,  zeta_0)
    endif
    
    df_m  = (x   - x_0  )*rzeta
    
    if(unstable) then
      call mo_derivative_uh(x  ,  zeta  )
      call mo_derivative_uh(x_0,  zeta_h)
    else
      call mo_derivative_s(x  ,  zeta  )
      call mo_derivative_s(x_0,  zeta_h)
    endif
    
    df_h  = (x   - x_0  )*rzeta
    
    if(unstable) then
      call mo_integral_um(f_m, zeta, zeta_0, ln_z_z0)
      call mo_integral_uh(f_h, zeta, zeta_h, ln_z_zt)
    else
      call mo_integral_s(f_m, zeta, zeta_0, ln_z_z0)
      call mo_integral_s(f_h, zeta, zeta_h, ln_z_zt)
    end if
    
    rich   = zeta*f_h/(f_m*f_m)
    d_rich = rich*( rzeta +  df_h/f_h - 2.0 *df_m/f_m)
    
    return
    end subroutine solve_zeta_internal
    
    !=======================================================================
    
    subroutine mo_derivative_um(phi, zeta)
    
    real   , intent(out),  dimension(:) :: phi
    real   , intent(in),   dimension(:) :: zeta
    
    
    !     phi = (1 - 16.0*zeta)**(-0.25)
    
        phi = (1 - 16.0*zeta)**(-0.5)
        phi = sqrt(phi)
    
    return
    end subroutine mo_derivative_um
    !=======================================================================
    
    subroutine mo_derivative_uh(phi, zeta)
    
    real   , intent(out),  dimension(:) :: phi
    real   , intent(in),   dimension(:) :: zeta
    
    
         phi = (1 - 16.0*zeta)**(-0.5)
    
    
    return
    end subroutine mo_derivative_uh
    !=======================================================================
    
    subroutine mo_derivative_s(phi, zeta)
    
    real   , intent(out),  dimension(:) :: phi
    real   , intent(in),   dimension(:) :: zeta
    
    
      phi = 1.0 + zeta*(5.0 + b_stab*zeta)/(1.0 + zeta)
    
    return
    end subroutine mo_derivative_s
    
    !=======================================================================
    
    subroutine mo_integral_um (psi, zeta, zeta_0, ln_z_z0_in)
    
    real   , intent(out), dimension(:) :: psi
    real   , intent(in), dimension(:)  :: zeta, zeta_0
    real   , optional,   dimension(:)  :: ln_z_z0_in
    
    real, dimension(size(zeta)) :: x, x_0, x1, x1_0, ln_z_z0, num, denom, y
    
    if(present(ln_z_z0_in)) then
       ln_z_z0 = ln_z_z0_in
    else
       ln_z_z0 = log(zeta/zeta_0)
    end if
    
    !     x      = (1 - 16.0*zeta)**(0.25)
    !     x_0    = (1 - 16.0*zeta_0)**(0.25)
    
         x      = sqrt(sqrt(1 - 16.0*zeta))
         x_0    = sqrt(sqrt(1 - 16.0*zeta_0))
    
         x1     = 1.0 + x
         x1_0   = 1.0 + x_0
         num    = x1*x1*(1.0 + x*x)
         denom  = x1_0*x1_0*(1.0 + x_0*x_0)
         y = atan(x) - atan(x_0)
         psi    = ln_z_z0 - log(num/denom) + 2*y
    
    return
    end subroutine mo_integral_um
    
    !=======================================================================
    subroutine mo_integral_uh (psi, zeta, zeta_0, ln_z_z0_in)
    
    real   , intent(out), dimension(:) :: psi
    real   , intent(in), dimension(:)  :: zeta, zeta_0
    real   , optional,   dimension(:)  :: ln_z_z0_in
    
    real, dimension(size(zeta)) :: x, x_0, ln_z_z0
    
    if(present(ln_z_z0_in)) then
       ln_z_z0 = ln_z_z0_in
    else
       ln_z_z0 = log(zeta/zeta_0)
    end if
    
         x     = sqrt(1 - 16.0*zeta)
         x_0   = sqrt(1 - 16.0*zeta_0)
         psi   = ln_z_z0 - 2.0*log( (1.0 + x)/(1.0 + x_0) )
    
    return
    end subroutine mo_integral_uh
    
    !=======================================================================
    
    subroutine mo_integral_s (psi, zeta, zeta_0, ln_z_z0_in)
    
    real   , intent(out), dimension(:) :: psi
    real   , intent(in), dimension(:)  :: zeta, zeta_0
    real   , optional,   dimension(:)  :: ln_z_z0_in
    
    real, dimension(size(zeta)) :: ln_z_z0
    
    if(present(ln_z_z0_in)) then
       ln_z_z0 = ln_z_z0_in
    else
       ln_z_z0 = log(zeta/zeta_0)
    end if
    
    
      psi = ln_z_z0 + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_0)) &
           + b_stab*(zeta - zeta_0) 
    
    return
    end subroutine mo_integral_s
    
    !=======================================================================
    
    subroutine mo_profile_1d(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, del_m, del_h, del_q, mask)
    ! zref_t only for compatability with lima revision of flux_exchange.f90. Nothing is done with it.
    
    real, intent(in)                :: zref, zref_t
    real, intent(in) , dimension(:) :: z, z0, zt, zq, u_star, b_star, q_star
    real, intent(out), dimension(:) :: del_m, del_h, del_q
    logical, optional, dimension(:) :: mask
    
    integer, dimension(size(z)) :: stabl, unstabl
    real   , allocatable, dimension(:) :: &
         zeta_1d, zeta_0_1d, zeta_h_1d, zeta_q_1d, zeta_ref_1d,   &
         ln_z_z0_1d, ln_z_zt_1d, ln_z_zq_1d, ln_z_zref_1d,         &
         f_m_1d, f_h_1d, f_q_1d, f_m_ref_1d, f_h_ref_1d, f_q_ref_1d
    
    integer :: i,j, npts,nptu
    
    real, dimension(size(z)) :: zeta, zeta_0, zeta_h, zeta_q, zeta_ref,  &
                                ln_z_z0, ln_z_zt, ln_z_zq, ln_z_zref,     &
                                f_m_ref, f_m_0, f_h_ref, f_h_0, f_q_ref, f_q_0, &
                                mo_length_inv
    
    logical, dimension(size(z)) :: avail, avail1
    
    
    if(.not.init) call monin_obukhov_init
    
    !-- zero output arrays --
        del_m = 0.0
        del_h = 0.0
        del_q = 0.0
    
    avail = .true.
    if (present(mask)) avail = mask
    
    where(avail) 
      ln_z_z0   = log(z/z0)
      ln_z_zt   = log(z/zt)
      ln_z_zq   = log(z/zq)
      ln_z_zref = log(z/zref)
    end where
    
    if(neutral) then
    
      where(avail)
        del_m = 1.0 - ln_z_zref/ln_z_z0
        del_h = 1.0 - ln_z_zref/ln_z_zt
        del_q = 1.0 - ln_z_zref/ln_z_zq
      end where
    
    else
    
      avail1 = avail .and. (u_star.ne.0.0)
    
      where(u_star .eq. 0.0)
        del_m = 0.0
        del_h = 0.0
        del_q = 0.0
      end where
    
      where(avail1) 
        mo_length_inv = - vonkarm * b_star/(u_star*u_star)
        zeta      = z   *mo_length_inv
        zeta_0    = z0  *mo_length_inv
        zeta_h    = zt  *mo_length_inv
        zeta_q    = zq  *mo_length_inv
        zeta_ref  = zref*mo_length_inv
      elsewhere
        zeta      = 0.0
      end where
    
      call separate_stabilities(npts, nptu, stabl, unstabl, zeta, avail1)
    
      if(nptu > 0) then
    !
    !    allocate(zeta_1d      (nptu))
    !    allocate(zeta_0_1d    (nptu))
    !    allocate(zeta_h_1d    (nptu))
    !    allocate(zeta_q_1d    (nptu))
    !    allocate(zeta_ref_1d  (nptu))
    !    allocate(ln_z_z0_1d   (nptu))
    !    allocate(ln_z_zt_1d   (nptu))
    !    allocate(ln_z_zq_1d   (nptu))
    !    allocate(ln_z_zref_1d (nptu))
    !    allocate(f_m_1d       (nptu))
    !    allocate(f_h_1d       (nptu))
    !    allocate(f_q_1d       (nptu))
    !    allocate(f_m_ref_1d   (nptu))
    !    allocate(f_h_ref_1d   (nptu))
    !    allocate(f_q_ref_1d   (nptu))
    !
    !    do i =1, nptu
    !       zeta_1d   (i)   = zeta     (unstabl(i))
    !       zeta_0_1d (i)   = zeta_0   (unstabl(i))
    !       zeta_h_1d (i)   = zeta_h   (unstabl(i))
    !       zeta_q_1d (i)   = zeta_q   (unstabl(i))
    !       zeta_ref_1d (i) = zeta_ref (unstabl(i))
    !       ln_z_z0_1d(i)   = ln_z_z0  (unstabl(i))
    !       ln_z_zt_1d(i)   = ln_z_zt  (unstabl(i))
    !       ln_z_zq_1d(i)   = ln_z_zq  (unstabl(i))
    !       ln_z_zref_1d(i) = ln_z_zref(unstabl(i))
    !    end do
    !
    !    call mo_integral_um(f_m_1d     ,zeta_1d,  zeta_0_1d,   ln_z_z0_1d  )
    !    call mo_integral_uh(f_h_1d     ,zeta_1d,  zeta_h_1d,   ln_z_zt_1d  )
    !    call mo_integral_uh(f_q_1d     ,zeta_1d,  zeta_q_1d,   ln_z_zq_1d  )
    !    call mo_integral_um(f_m_ref_1d ,zeta_1d,  zeta_ref_1d, ln_z_zref_1d)
    !    call mo_integral_uh(f_h_ref_1d ,zeta_1d,  zeta_ref_1d, ln_z_zref_1d)
    !    f_q_ref_1d = f_h_ref_1d
    !
        do i = 1, nptu
          del_m(unstabl(i)) = 1.0 - ln_z_zref(unstabl(i))/ln_z_z0(unstabl(i))
          del_h(unstabl(i)) = 1.0 - ln_z_zref(unstabl(i))/ln_z_zt(unstabl(i))
          del_q(unstabl(i)) = 1.0 - ln_z_zref(unstabl(i))/ln_z_zq(unstabl(i))
    !      del_m(unstabl(i)) = 1.0 - f_m_ref_1d(i)/f_m_1d(i)
    !      del_h(unstabl(i)) = 1.0 - f_h_ref_1d(i)/f_h_1d(i)
    !      del_q(unstabl(i)) = 1.0 - f_q_ref_1d(i)/f_q_1d(i)
        end do
    
    !    deallocate(zeta_1d      )
    !    deallocate(zeta_0_1d    )
    !    deallocate(zeta_h_1d    )
    !    deallocate(zeta_q_1d    )
    !    deallocate(zeta_ref_1d  )
    !    deallocate(ln_z_z0_1d   )
    !    deallocate(ln_z_zt_1d   )
    !    deallocate(ln_z_zq_1d   )
    !    deallocate(ln_z_zref_1d )
    !    deallocate(f_m_1d       )
    !    deallocate(f_h_1d       )
    !    deallocate(f_q_1d       )
    !    deallocate(f_m_ref_1d   )
    !    deallocate(f_h_ref_1d   )
    !    deallocate(f_q_ref_1d   )
    
      end if
    
      if(npts > 0) then
    
    !    allocate(zeta_1d      (npts))
    !    allocate(zeta_0_1d    (npts))
    !    allocate(zeta_h_1d    (npts))
    !    allocate(zeta_q_1d    (npts))
    !    allocate(zeta_ref_1d  (npts))
    !    allocate(ln_z_z0_1d   (npts))
    !    allocate(ln_z_zt_1d   (npts))
    !    allocate(ln_z_zq_1d   (npts))
    !    allocate(ln_z_zref_1d (npts))
    !    allocate(f_m_1d       (npts))
    !    allocate(f_h_1d       (npts))
    !    allocate(f_q_1d       (npts))
    !    allocate(f_m_ref_1d   (npts))
    !    allocate(f_h_ref_1d   (npts))
    !    allocate(f_q_ref_1d   (npts))
    
    !    do i =1, npts
    !       zeta_1d   (i)   = zeta     (stabl(i))
    !       zeta_0_1d (i)   = zeta_0   (stabl(i))
    !       zeta_h_1d (i)   = zeta_h   (stabl(i))
    !       zeta_q_1d (i)   = zeta_q   (stabl(i))
    !       zeta_ref_1d (i) = zeta_ref (stabl(i))
    !       ln_z_z0_1d(i)   = ln_z_z0  (stabl(i))
    !       ln_z_zt_1d(i)   = ln_z_zt  (stabl(i))
    !       ln_z_zq_1d(i)   = ln_z_zq  (stabl(i))
    !       ln_z_zref_1d(i) = ln_z_zref(stabl(i))
    !    end do
    
    !    call mo_integral_s(f_m_1d     , zeta_1d, zeta_0_1d,   ln_z_z0_1d  )
    !    call mo_integral_s(f_h_1d     , zeta_1d, zeta_h_1d,   ln_z_zt_1d  )
    !    call mo_integral_s(f_q_1d     , zeta_1d, zeta_q_1d,   ln_z_zq_1d  )
    !    call mo_integral_s(f_m_ref_1d , zeta_1d, zeta_ref_1d, ln_z_zref_1d)
    !    f_h_ref_1d = f_m_ref_1d
    !    f_q_ref_1d = f_m_ref_1d
    
        do i = 1, npts
          del_m(stabl(i)) = ( ln_z_z0(stabl(i)) - ln_z_zref(stabl(i)) + &
                    zeta_ref(stabl(i)) - zeta_0(stabl(i)) )/ &
                    ( ln_z_z0(stabl(i)) + zeta(stabl(i)) - zeta_0(stabl(i)) )
          del_h(stabl(i)) = ( ln_z_zt(stabl(i)) - ln_z_zref(stabl(i)) + &
                    zeta_ref(stabl(i)) - zeta_h(stabl(i)) )/ &
                    ( ln_z_zt(stabl(i)) + zeta(stabl(i)) - zeta_h(stabl(i)) )
          del_q(stabl(i)) = ( ln_z_zq(stabl(i)) - ln_z_zref(stabl(i)) + &
                    zeta_ref(stabl(i)) - zeta_q(stabl(i)) )/ &
                    ( ln_z_zq(stabl(i)) + zeta(stabl(i)) - zeta_q(stabl(i)) )
    !      del_m(stabl(i)) = 1.0 - f_m_ref_1d(i)/f_m_1d(i)
    !      del_h(stabl(i)) = 1.0 - f_h_ref_1d(i)/f_h_1d(i)
    !      del_q(stabl(i)) = 1.0 - f_q_ref_1d(i)/f_q_1d(i)
        end do
    
    !    deallocate(zeta_1d      )
    !    deallocate(zeta_0_1d    )
    !    deallocate(zeta_h_1d    )
    !    deallocate(zeta_q_1d    )
    !    deallocate(zeta_ref_1d  )
    !    deallocate(ln_z_z0_1d   )
    !    deallocate(ln_z_zt_1d   )
    !    deallocate(ln_z_zq_1d   )
    !    deallocate(ln_z_zref_1d )
    !    deallocate(f_m_1d       )
    !    deallocate(f_h_1d       )
    !    deallocate(f_q_1d       )
    !    deallocate(f_m_ref_1d   )
    !    deallocate(f_h_ref_1d   )
    !    deallocate(f_q_ref_1d   )
    
      end if
    
    end if
    
    return
    end subroutine mo_profile_1d
    
    !=======================================================================
    
    subroutine mo_diff_1d(z, u_star, b_star, k_m, k_h, mask)
    
    real, intent(in),  dimension(:,:) :: z
    real, intent(in),  dimension(:)   :: u_star, b_star
    real, intent(out), dimension(:,:) :: k_m, k_h
    logical, optional, dimension(:)   :: mask
    
    
    real   , dimension(size(z,1),size(z,2)) :: phi_m , phi_h , zeta
    real   , allocatable, dimension(:)      :: phi_m1, phi_h1, zeta1
    real   , dimension(size(z,1))           :: mo_length_inv
    integer, dimension(size(z,1))           :: stab, unstab
    logical, dimension(size(z,1))           :: avail, avail1
    integer :: n, i, k, ns, nu
    
    
    if(.not.init) call monin_obukhov_init
    
    avail = .true.
    if (present(mask)) avail = mask
    
    if(neutral) then
    
      do k =1, size(z,2)
        where(avail) 
            k_m(:,k) = vonkarm * u_star * z(:,k)
            k_h(:,k) = k_m(:,k)
        end where
      end do
    
    else
    
      avail1 = avail .and. (u_star.ne.0.0)
    
      do k =1, size(z,2)
        where(u_star .eq. 0.0)
          k_m(:,k) = 0.0
          k_h(:,k) = 0.0
        end where
      end do
    
      do k =1, size(z,2)
        where(avail1) 
           mo_length_inv = - vonkarm * b_star/(u_star*u_star)
           zeta(:,k) = z(:,k)*mo_length_inv
        elsewhere
           mo_length_inv = 0.0
        endwhere
      end do 
    
      call separate_stabilities(ns, nu, stab, unstab, mo_length_inv, avail1)
    
      if(nu > 0) then
    
    
    !    allocate(zeta1 (nu))
    !    allocate(phi_m1(nu))
    !    allocate(phi_h1(nu))
    
    
        do k = 1, size(z,2)
    !      do i = 1, nu
    !        zeta1(i) = zeta(unstab(i),k)
    !      end do
    
    !      call mo_derivative_um(phi_m1, zeta1)
    !      call mo_derivative_uh(phi_h1, zeta1)
    
          do i = 1, nu
    !        phi_m(unstab(i),k) = phi_m1(i)
    !        phi_h(unstab(i),k) = phi_h1(i)
            phi_m(unstab(i),k) = 1.0
            phi_h(unstab(i),k) = 1.0
          end do
        end do
    
    !    deallocate(zeta1)
    !    deallocate(phi_m1)
    !    deallocate(phi_h1)
    
      end if
    
      if(ns > 0) then
    
    !    allocate(zeta1 (ns))
    !    allocate(phi_m1(ns))
    !    allocate(phi_h1(ns))
    
        do k = 1, size(z,2)
    !      do i =1, ns
    !        zeta1(i) = zeta(stab(i),k)
    !      end do
    
    !      call mo_derivative_s(phi_m1, zeta1)
    !      phi_h1 = phi_m1
    
          do i = 1, ns
            phi_m(stab(i),k) = 1.0 + zeta(stab(i),k)
            phi_h(stab(i),k) = 1.0 + zeta(stab(i),k)
    !        phi_m(stab(i),k) = phi_m1(i)
    !        phi_h(stab(i),k) = phi_h1(i)
          end do
        end do
    
    !    deallocate(zeta1)
    !    deallocate(phi_m1)
    !    deallocate(phi_h1)
    
      end if
    
      do k =1, size(z,2)
        where(avail1) 
          k_m(:,k) = vonkarm * u_star * z(:,k)/phi_m(:,k)
          k_h(:,k) = vonkarm * u_star * z(:,k)/phi_h(:,k)
        end where
      end do
    
    endif
    
    return
    end subroutine mo_diff_1d
    
    !=======================================================================
    
    subroutine separate_stabilities(ns, nu, stab, unstab, x, mask, crit)
    
    real   , intent(in), dimension(:) :: x
    logical, intent(in), dimension(:) :: mask
    real,    intent(in),  optional    :: crit
    integer, intent(out) :: ns ,nu
    integer, intent(out), dimension(:) :: stab, unstab
    
    real    :: xcrit
    integer :: i,j
    
    ns=0
    nu=0
    
    xcrit = 1.e10
    if(present(crit)) xcrit = crit
    
    do i=1,size(x)
      if (mask(i) .and. x(i) >= 0.0 .and. x(i) < xcrit) then
         ns=ns+1
         stab(ns)=i
      else if (mask(i) .and. x(i) < 0.0) then
         nu=nu+1
         unstab(nu)=i
      endif
    enddo
    
    return
    end subroutine separate_stabilities
    
    !=======================================================================
    ! The following routines allow some of the above subroutines to be called
    ! with different dimensions of the input and output
    !
    !=======================================================================
    
    subroutine mo_drag_2d &
             (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)
    
    real, intent(in),  dimension(:,:) :: z, speed, pt, pt0, z0, zt, zq
    real, intent(out), dimension(:,:) :: drag_m, drag_t, drag_q
    real, intent(out), dimension(:,:) :: u_star, b_star
    
    real, dimension(size(pt,1)*size(pt,2)) :: z_1d, speed_1d, pt_1d, pt0_1d, &
       z0_1d, zt_1d, zq_1d, drag_m_1d, drag_t_1d, drag_q_1d, u_star_1d, b_star_1d
    
    integer :: shape1(1), shape2(2)
    
    shape1 =    size(pt,1) * size(pt,2)
    shape2 = (/ size(pt,1),  size(pt,2) /)
    
              z_1d = reshape (z      , shape1)
          speed_1d = reshape (speed  , shape1)
             pt_1d = reshape (pt     , shape1)
            pt0_1d = reshape (pt0    , shape1)
             z0_1d = reshape (z0     , shape1)
             zt_1d = reshape (zt     , shape1)
             zq_1d = reshape (zq     , shape1)
         u_star_1d = reshape (u_star , shape1)
         b_star_1d = reshape (b_star , shape1)
    
    call mo_drag_1d &
             (pt_1d, pt0_1d, z_1d, z0_1d, zt_1d, zq_1d, speed_1d,&
              drag_m_1d, drag_t_1d, drag_q_1d, u_star_1d, b_star_1d)
    
          drag_m = reshape (drag_m_1d   , shape2)
          drag_t = reshape (drag_t_1d   , shape2)
          drag_q = reshape (drag_q_1d   , shape2)
          u_star = reshape (u_star_1d   , shape2)
          b_star = reshape (b_star_1d   , shape2)
    
    
    return
    end subroutine mo_drag_2d
    !=======================================================================
    
    subroutine mo_profile_2d(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, del_m, del_h, del_q)
    ! zref_t only for compatability with lima revision of flux_exchange.f90. Nothing is done with it.
    
    real, intent(in)                  :: zref, zref_t
    real, intent(in) , dimension(:,:) :: z, z0, zt, zq, u_star, b_star, q_star
    real, intent(out), dimension(:,:) :: del_m, del_h, del_q
    
    real, dimension(size(z,1)*size(z,2)) :: z_1d, z0_1d,  zt_1d, zq_1d, &
                                              u_star_1d, b_star_1d, q_star_1d, &
                                              del_m_1d, del_h_1d, del_q_1d
    
    integer :: shape1(1), shape2(2)
    
    shape1 =    size(z,1) * size(z,2)
    shape2 = (/ size(z,1),  size(z,2) /)
    
              z_1d = reshape (z        , shape1)
             z0_1d = reshape (z0       , shape1)
             zt_1d = reshape (zt       , shape1)
             zq_1d = reshape (zq       , shape1)
         u_star_1d = reshape (u_star   , shape1)
         b_star_1d = reshape (b_star   , shape1)
         q_star_1d = reshape (q_star   , shape1)
    
    call mo_profile_1d&
        (zref, zref_t, z_1d, z0_1d, zt_1d, zq_1d, u_star_1d, b_star_1d, q_star_1d, del_m_1d, del_h_1d, del_q_1d)
    
          del_m = reshape (del_m_1d , shape2)
          del_h = reshape (del_h_1d , shape2)
          del_q = reshape (del_q_1d , shape2)
    
    return
    end subroutine mo_profile_2d
    
    !=======================================================================
    
    subroutine mo_diff_one_lev_1d(z, u_star, b_star, k_m, k_h, mask)
    
    real, intent(in),  dimension(:) :: z
    real, intent(in),  dimension(:)   :: u_star, b_star
    real, intent(out), dimension(:) :: k_m, k_h
    logical, optional, dimension(:)   :: mask
    
    real, dimension(size(z),1) :: z_nlev, k_m_nlev, k_h_nlev
    logical, dimension(size(z))  :: avail
    
    avail = .true.
    if (present(mask)) avail = mask
    
    z_nlev(:,1) = z
    
    call mo_diff_1d(z_nlev, u_star, b_star, k_m_nlev, k_h_nlev, avail)
    
    where(avail)
      k_m = k_m_nlev(:,1)
      k_h = k_h_nlev(:,1)
    end where
    
     
    return
    end subroutine mo_diff_one_lev_1d
    
    !=======================================================================
    
    subroutine mo_diff_one_lev_0d(z, u_star, b_star, k_m, k_h)
    
    real, intent(in)  :: z
    real, intent(in)  :: u_star, b_star
    real, intent(out) :: k_m, k_h
    
    real, dimension(1,1) :: z1, k_m1, k_h1
    real, dimension(1) :: u_star1, b_star1
    
    z1      = z
    u_star1 = u_star
    b_star1 = b_star
    call mo_diff_1d(z1, u_star1, b_star1, k_m1, k_h1)
    k_m = k_m1(1,1)
    k_h = k_h1(1,1)
    
    return
    end subroutine mo_diff_one_lev_0d
    
    !=======================================================================
    
    subroutine mo_diff_0d(z, u_star, b_star, k_m, k_h)
    
    real, intent(in),  dimension(:) :: z
    real, intent(in)                :: u_star, b_star
    real, intent(out), dimension(:) :: k_m, k_h
    
    real, dimension(size(z,1),1) :: z1, k_m1, k_h1
    real, dimension(1) :: u_star1, b_star1
    
    z1(:,1) = z
    u_star1 = u_star
    b_star1 = b_star
    call mo_diff_1d(z1, u_star1, b_star1, k_m1, k_h1)
    k_m = k_m1(:,1)
    k_h = k_h1(:,1)
    
    return
    end subroutine mo_diff_0d
    !=======================================================================
    
    subroutine mo_diff_2d(z, u_star, b_star, k_m, k_h)
    
    real, intent(in),  dimension(:,:,:) :: z
    real, intent(in),  dimension(:,:)   :: u_star, b_star
    real, intent(out), dimension(:,:,:) :: k_m, k_h
    
    real, dimension(size(z,1)*size(z,2),size(z,3)) :: z1, k_m1, k_h1
    
    real, dimension(size(z,1)*size(z,2)) :: u_star1, b_star1
    
    integer :: i, j, k, ii
    
    do j = 1, size(z,2)
      do i = 1, size(z,1)
        ii = i + (j -1)*size(z,1)
        u_star1(ii) = u_star(i,j)
        b_star1(ii) = b_star(i,j)
        do k = 1, size(z,3)
          z1(ii,k) = z(i,j,k)
        end do
      end do
    end do
    
    call mo_diff_1d(z1,u_star1,b_star1,k_m1,k_h1)
    
    do j = 1, size(z,2)
      do i = 1, size(z,1)
        ii = i + (j -1)*size(z,1)
        do k = 1, size(z,3)
          k_m(i,j,k) = k_m1(ii,k)
          k_h(i,j,k) = k_h1(ii,k)
        end do
      end do
    end do
    
    return
    end subroutine mo_diff_2d
    
    !=======================================================================
    
    subroutine mo_diff_one_lev_2d(z, u_star, b_star, k_m, k_h)
    
    real, intent(in),  dimension(:,:) :: z
    real, intent(in),  dimension(:,:) :: u_star, b_star
    real, intent(out), dimension(:,:) :: k_m, k_h
    
    real, dimension(size(z,1),size(z,2),1) :: z_nlev, k_m_nlev, k_h_nlev
    
    z_nlev(:,:,1) = z
    
    call mo_diff_2d(z_nlev, u_star,b_star, k_m_nlev, k_h_nlev)
    
    k_m = k_m_nlev(:,:,1)
    k_h = k_h_nlev(:,:,1)
      
    return
    end subroutine mo_diff_one_lev_2d
    !=======================================================================
    
    subroutine stable_mix_3d(rich, mix)
    
    real, intent(in) , dimension(:,:,:)  :: rich
    real, intent(out), dimension(:,:,:)  :: mix
    
    real, dimension(size(rich,1),size(rich,2),size(rich,3)) :: &
          r, a, b, c, zeta, phi
        
    mix = 0.0 
      
    if(stable_option == 1) then
    
      where(rich > 0.0 .and. rich < rich_crit)
        r = 1.0/rich
        a = r - b_stab
        b = r - (1.0 + 5.0)
        c = - 1.0
        
        zeta = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a)
        phi = 1.0 + b_stab*zeta + (5.0 - b_stab)*zeta/(1.0 + zeta)
        mix = 1./(phi*phi)
      end where
    
    else if(stable_option == 2) then
    
      where(rich > 0.0 .and. rich <= rich_trans)
        mix = (1.0 - 5.0*rich)**2
      end where
      where(rich > rich_trans .and. rich < rich_crit)
        mix = ((1.0 - b_stab*rich)/lambda)**2
      end where
    
    end if
    
    return
    end subroutine stable_mix_3d
    !=======================================================================
    
    subroutine stable_mix_2d(rich, mix)
    
    real, intent(in) , dimension(:,:)  :: rich
    real, intent(out), dimension(:,:)  :: mix
    
    real, dimension(size(rich,1),size(rich,2),1) :: rich_3d, mix_3d
    
    rich_3d(:,:,1) = rich
    
    call stable_mix_3d(rich_3d, mix_3d)
    
    mix = mix_3d(:,:,1)
    
    return
    end subroutine stable_mix_2d
    
    !=======================================================================
    
    subroutine stable_mix_1d(rich, mix)
    
    real, intent(in) , dimension(:)  :: rich
    real, intent(out), dimension(:)  :: mix
    
    real, dimension(size(rich),1,1) :: rich_3d, mix_3d
    
    rich_3d(:,1,1) = rich
    
    call stable_mix_3d(rich_3d, mix_3d)
    
    mix = mix_3d(:,1,1)
    
    return
    end subroutine stable_mix_1d
    
    !=======================================================================
    
    subroutine stable_mix_0d(rich, mix)
    
    real, intent(in) :: rich
    real, intent(out) :: mix
    
    real, dimension(1,1,1) :: rich_3d, mix_3d
    
    rich_3d(1,1,1) = rich
    
    call stable_mix_3d(rich_3d, mix_3d)
    
    mix = mix_3d(1,1,1)
    
    return
    end subroutine stable_mix_0d
    !=======================================================================
    
    end module frierson_monin_obukhov_mod
    