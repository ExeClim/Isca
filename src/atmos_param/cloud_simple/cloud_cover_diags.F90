module cloud_cover_diags_mod

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  use           fms_mod, only: stdlog, FATAL, WARNING, NOTE, error_mesg, &
                               uppercase, check_nml_error
  use  time_manager_mod, only: time_type
  use  diag_manager_mod, only: register_diag_field, send_data

  implicit none

  character(len=14), parameter :: mod_name = "cloud_cover"

  ! Overlap assumptions include 'maximum-random', 'maximum', and 'random',
  ! and the default is 'maximum-random'.
  character(len=32) :: overlap_assumption = 'maximum-random'

  ! If true, output cloud covers based on all of the three assumptions
  logical :: do_test_overlap = .false.

  real :: cf_min = 0.0 ! The cloud free threshold

  real :: mid_btm  = 7.0e4 ! Bottom (Top) pressure of middle (low) clouds
  real :: high_btm = 4.0e4 ! Bottom (Top) pressure of high (middle) clouds

  ! ----- outputs for cloud amount diagnostics ----- !
  integer :: id_tot_cld_amt, id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt

  ! if do_test_overlap is .true., then output cloud covers based on all of the three assumptions
  integer :: id_tot_cld_amt_mxr, id_high_cld_amt_mxr, id_mid_cld_amt_mxr, id_low_cld_amt_mxr
  integer :: id_tot_cld_amt_max, id_high_cld_amt_max, id_mid_cld_amt_max, id_low_cld_amt_max
  integer :: id_tot_cld_amt_rnd, id_high_cld_amt_rnd, id_mid_cld_amt_rnd, id_low_cld_amt_rnd

  namelist /cloud_cover_diag_nml/ &
              overlap_assumption, cf_min, mid_btm, high_btm, do_test_overlap

  contains

  ! ===================================================
  !             cloud cover diags init
  ! ===================================================
  subroutine cloud_cover_diags_init(axes, Time)
    type(time_type), intent(in)       :: Time
    integer, intent(in), dimension(4) :: axes
    integer :: io, ierr, nml_unit, stdlog_unit

#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=cloud_cover_diag_nml, iostat=io)
    ierr = check_nml_error(io, 'cloud_cover_diag_nml')
#else
    if (file_exist('input.nml')) then
      nml_unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
        read(nml_unit, nml=cloud_cover_diag_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'cloud_cover_diag_nml')
      enddo
10    call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, cloud_cover_diag_nml)

    id_tot_cld_amt = register_diag_field (mod_name, 'tot_cld_amt', axes(1:2), Time, &
                          'total cloud amount (%)', 'percent')
    id_high_cld_amt = register_diag_field (mod_name, 'high_cld_amt', axes(1:2), Time, &
                          'high cloud amount (%)', 'percent')
    id_mid_cld_amt = register_diag_field (mod_name, 'mid_cld_amt', axes(1:2), Time, &
                          'mid cloud amount (%)', 'percent')
    id_low_cld_amt = register_diag_field (mod_name, 'low_cld_amt', axes(1:2), Time, &
                          'low cloud amount (%)', 'percent')

    if (do_test_overlap) then
      ! Maximum-random overlap
      id_tot_cld_amt_mxr = register_diag_field (mod_name, 'tot_cld_amt_mxr', axes(1:2), Time, &
              'total cloud amount (maximum-random assumption)', 'percent')
      id_high_cld_amt_mxr = register_diag_field (mod_name, 'high_cld_amt_mxr', axes(1:2), Time, &
              'high cloud amount (maximum-random assumption)', 'percent')
      id_mid_cld_amt_mxr = register_diag_field (mod_name, 'mid_cld_amt_mxr', axes(1:2), Time, &
              'mid cloud amount (maximum-random assumption)', 'percent')
      id_low_cld_amt_mxr = register_diag_field (mod_name, 'low_cld_amt_mxr', axes(1:2), Time, &
              'low cloud amount (maximum-random assumption)', 'percent')

      ! Maximum overlap
      id_tot_cld_amt_max = register_diag_field (mod_name, 'tot_cld_amt_max', axes(1:2), Time, &
              'total cloud amount (maximum assumption)', 'percent')
      id_high_cld_amt_max = register_diag_field (mod_name, 'high_cld_amt_max', axes(1:2), Time, &
              'high cloud amount (maximum assumption)', 'percent')
      id_mid_cld_amt_max = register_diag_field (mod_name, 'mid_cld_amt_max', axes(1:2), Time, &
              'mid cloud amount (maximum assumption)', 'percent')
      id_low_cld_amt_max = register_diag_field (mod_name, 'low_cld_amt_max', axes(1:2), Time, &
              'low cloud amount (maximum assumption)', 'percent')

      ! Random overlap
      id_tot_cld_amt_rnd = register_diag_field (mod_name, 'tot_cld_amt_rnd', axes(1:2), Time, &
              'total cloud amount (random assumption)', 'percent')
      id_high_cld_amt_rnd = register_diag_field (mod_name, 'high_cld_amt_rnd', axes(1:2), Time, &
              'high cloud amount (random assumption)', 'percent')
      id_mid_cld_amt_rnd = register_diag_field (mod_name, 'mid_cld_amt_rnd', axes(1:2), Time, &
              'mid cloud amount (random assumption)', 'percent')
      id_low_cld_amt_rnd = register_diag_field (mod_name, 'low_cld_amt_rnd', axes(1:2), Time, &
              'low cloud amount (random assumption)', 'percent')
    end if

  end subroutine cloud_cover_diags_init

  ! ===================================================
  !             cloud cover diags init
  ! ===================================================
  subroutine cloud_cover_diags(cf, p_full, p_half, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full, p_half
    type(time_type),   intent(in)       :: Time
    character(len=32) :: overlap_str = ''
    real, dimension(size(cf,1), size(cf,2)) :: tca, low_ca, mid_ca, high_ca

    overlap_str = uppercase(trim(overlap_assumption))

    if (overlap_str == 'MAXIMUM-RANDOM') then
      call diag_cldamt_maxrnd_overlap(cf, p_full, p_half, Time, tca, high_ca, mid_ca, low_ca)

    else if (overlap_str == 'MAXIMUM') then
      call diag_cldamt_max_overlap(cf, p_full, Time, tca, high_ca, mid_ca, low_ca)

    else if (overlap_str == 'RANDOM') then
      call diag_cldamt_random_overlap(cf, p_full, Time, tca, high_ca, mid_ca, low_ca)

    else
      call error_mesg(mod_name, '"'//trim(overlap_assumption)//'"'// &
        ' is not a valid cloud overlap assumption.', FATAL)
    end if

    ! Diagnostics output
    call output_cldamt(tca, high_ca, mid_ca, low_ca, Time)

    if (do_test_overlap) then
      call diag_cldamt_maxrnd_overlap(cf, p_full, p_half, Time, tca, high_ca, mid_ca, low_ca)
      call diag_cldamt_max_overlap(cf, p_full, Time, tca, high_ca, mid_ca, low_ca)
      call diag_cldamt_random_overlap(cf, p_full, Time, tca, high_ca, mid_ca, low_ca)
    end if

  end subroutine cloud_cover_diags

  subroutine cldovrlap(pint, cld, nmxrgn, pmxrgn)
    ! The original codes are from CAM.
    ! Please refer to:
    !  https://github.com/E3SM-Project/E3SM/blob/master/components/cam/src/physics/cam/pkg_cldoptics.F90#L136

    !-----------------------------------------------------------------------
    ! Purpose:
    ! Partitions each column into regions with clouds in neighboring layers.
    ! This information is used to implement maximum overlap in these regions
    ! with random overlap between them.
    ! On output,
    !    nmxrgn contains the number of regions in each column
    !    pmxrgn contains the interface pressures for the lower boundaries of
    !           each region!
    ! Author: W. Collins
    !-----------------------------------------------------------------------

    ! Input arguments
    real, intent(in), dimension(:,:,:) :: pint  ! Interface pressure
    real, intent(in), dimension(:,:,:) :: cld   ! Fractional cloud cover

    ! Output arguments
    ! Number of maximally overlapped regions
    integer, intent(out), dimension(size(cld,1), size(cld,2)) :: nmxrgn
    real, intent(out), dimension(size(pint,1), size(pint,2), size(pint,3)) :: pmxrgn

    !---------------------------Local variables-----------------------------
    integer :: i, j  ! Lat/Longitude index
    integer :: k     ! Level index
    integer :: n     ! Max-overlap region counter

    real, dimension(size(pint,1), size(pint,2), size(pint,3)) :: pnm  ! Interface pressure
    logical :: cld_found                          ! Flag for detection of cloud
    logical, dimension(size(cld,3)) :: cld_layer  ! Flag for cloud in layer
    integer :: pver, pverp

    pver = size(cld,3)
    pverp = pver + 1

    do i = 1,size(cld,1)
      do j = 1,size(cld,2)
        cld_found = .false.
        ! True if cloud fraction greater than cf_min
        cld_layer(:) = cld(i,j,:) > cf_min    ! 0.0
        pmxrgn(i,j,:) = 0.0
        pnm(i,j,:) = pint(i,j,:)
        n = 1
        do k = 1, pver
          if (cld_layer(k) .and. .not. cld_found) then
              cld_found = .true.
          else if (.not. cld_layer(k) .and. cld_found) then
            cld_found = .false.
            if (count(cld_layer(k:pver)) == 0) then
              exit
            endif
            pmxrgn(i,j,n) = pnm(i,j,k)
            n = n + 1
          endif
        end do
        pmxrgn(i,j,n) = pnm(i,j,pverp)
        nmxrgn(i,j) = n
      end do
    end do
  end subroutine cldovrlap

  subroutine diag_cldamt_maxrnd_overlap(cld, pmid, pint, Time, tca, high_ca, mid_ca, low_ca)
    ! Original codes are from CESM, refer to:
    ! https://github.com/E3SM-Project/E3SM/blob/master/components/cam/src/physics/cam/cloud_cover_diags.F90

    !!! pmid is pfull, and pint is phalf
    real, intent(in),  dimension(:,:,:) :: cld, pmid, pint
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:)   :: tca, high_ca, mid_ca, low_ca

    !---------------------------Local workspace-----------------------------
    ! Total, low, middle and high random overlap cloud cover
    real :: used

    integer :: i,j,k       ! lat/lon,level indices
    integer, dimension(size(cld,1), size(cld,2)) :: irgn    ! Max-overlap region index
    integer :: max_nmxrgn  ! maximum value of nmxrgn over columns
    integer :: ityp        ! Type counter
    real, dimension(size(cld,1), size(cld,2)) :: clrsky      ! Max-random clear sky fraction
    real, dimension(size(cld,1), size(cld,2)) :: clrskymax   ! Maximum overlap clear sky fraction
    integer, dimension(size(cld,1), size(cld,2)) :: nmxrgn
    real, dimension(size(pint,1), size(pint,2), size(pint,3)) :: pmxrgn

    !------------------------------Cloud Range Paramters-------------------------------
    real :: plowmax=1.2e5,  plowmin         ! Max/min prs for low cloud cover range
    real :: pmedmax,        pmedmin         ! Max/min prs for mid cloud cover range
    real :: phghmax,        phghmin=5.0e3   ! Max/min prs for hgh cloud cover range

    real, dimension(4) :: ptypmin
    real, dimension(4) :: ptypmax

    plowmin = mid_btm
    pmedmax = mid_btm
    pmedmin = high_btm
    phghmax = high_btm

    ptypmin = (/phghmin, plowmin, pmedmin, phghmin/)
    ptypmax = (/plowmax, plowmax, pmedmax, phghmax/)

    ! call the overlap subroutine to obtain the nmxrgn and pmxrgn
    call cldovrlap(pint, cld, nmxrgn, pmxrgn)

    ! Initialize region number
    max_nmxrgn = -1
    do i=1,size(cld,1)
      do j=1,size(cld,2)
        max_nmxrgn = max(max_nmxrgn, nmxrgn(i,j))
      end do
    end do

    do ityp= 1,4
      irgn = 1
      do k =1,max_nmxrgn-1
          do i=1,size(cld,1)
            do j=1,size(cld,2)
              if (pmxrgn(i,j,irgn(i,j)) < ptypmin(ityp) .and. irgn(i,j) < nmxrgn(i,j)) then
                  irgn(i,j) = irgn(i,j) + 1
              end if
            end do
          end do
      end do
      !
      ! Compute cloud amount by estimating clear-sky amounts
      !
      clrsky = 1.0
      clrskymax = 1.0

      do i=1,size(cld,1)
        do j=1,size(cld,2)
          do k=1,size(cld,3)
            if (pmid(i,j,k) >= ptypmin(ityp) .and. pmid(i,j,k) <= ptypmax(ityp)) then
                if (pmxrgn(i,j,irgn(i,j)) < pmid(i,j,k) .and. irgn(i,j) < nmxrgn(i,j)) then
                  irgn(i,j) = irgn(i,j) + 1
                  clrsky(i,j) = clrsky(i,j) * clrskymax(i,j)
                  clrskymax(i,j) = 1.0
                endif
                clrskymax(i,j) = min(clrskymax(i,j), 1.0-cld(i,j,k))
            endif
          end do
        end do
      end do

      if (ityp == 1) tca = 1.0 - (clrsky * clrskymax)
      if (ityp == 2) low_ca = 1.0 - (clrsky * clrskymax)
      if (ityp == 3) mid_ca = 1.0 - (clrsky * clrskymax)
      if (ityp == 4) high_ca = 1.0 - (clrsky * clrskymax)
    end do

    if (do_test_overlap) then
      if (id_tot_cld_amt_mxr > 0) then
        used = send_data (id_tot_cld_amt_mxr, tca*1.0e2, Time)
      endif
      if (id_high_cld_amt_mxr > 0) then
        used = send_data (id_high_cld_amt_mxr, high_ca*1.0e2, Time)
      endif
      if (id_mid_cld_amt_mxr > 0) then
        used = send_data (id_mid_cld_amt_mxr, mid_ca*1.0e2, Time)
      endif
      if (id_low_cld_amt_mxr > 0) then
        used = send_data (id_low_cld_amt_mxr, low_ca*1.0e2, Time)
      endif
    endif

  end subroutine diag_cldamt_maxrnd_overlap

  subroutine diag_cldamt_max_overlap(cf, p_full, Time, tca, high_ca, mid_ca, low_ca)
    real, intent(in),  dimension(:,:,:) :: cf, p_full
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:)   :: tca, high_ca, mid_ca, low_ca
    integer :: i, j, ks, ke
    logical, dimension(size(cf,3)) :: ind_mid
    real :: used

    tca = 1.0
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    ! total cld amount
    tca = maxval(cf, 3)

    do i=1,size(cf,1)
      do j=1,size(cf,2)
        ! low cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        low_ca(i,j) = maxval(cf(i,j,ks:ke), 1)

        ! middle cloud amount
        ind_mid = high_btm<=p_full(i,j,:) .and. p_full(i,j,:)<=mid_btm
        ks = minloc(p_full(i,j,:), 1, mask=ind_mid)
        ke = maxloc(p_full(i,j,:), 1, mask=ind_mid)
        mid_ca(i,j) = maxval(cf(i,j,ks:ke), 1)

        ! high cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        high_ca(i,j) = maxval(cf(i,j,ks:ke), 1)
      enddo
    enddo

    if (do_test_overlap) then
      if (id_tot_cld_amt_max > 0) then
        used = send_data (id_tot_cld_amt_max, tca*1.0e2, Time)
      endif
      if (id_high_cld_amt_max > 0) then
        used = send_data (id_high_cld_amt_max, high_ca*1.0e2, Time)
      endif
      if (id_mid_cld_amt_max > 0) then
        used = send_data (id_mid_cld_amt_max, mid_ca*1.0e2, Time)
      endif
      if (id_low_cld_amt_max > 0) then
        used = send_data (id_low_cld_amt_max, low_ca*1.0e2, Time)
      endif
    endif

  end subroutine diag_cldamt_max_overlap

  subroutine diag_cldamt_random_overlap(cf, p_full, Time, tca, high_ca, mid_ca, low_ca)
    real, intent(in),  dimension(:,:,:) :: cf, p_full
    type(time_type),   intent(in)       :: Time
    real, intent(out), dimension(:,:)   :: tca, high_ca, mid_ca, low_ca
    integer :: i, j, ks, ke
    logical, dimension(size(cf,3)) :: ind_mid
    real :: used

    tca = 1.0
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    do i=1,size(cf,1)
      do j=1,size(cf,2)
        ! low cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        call random_overlap_single_column(cf(i,j,ks:ke), low_ca(i,j))

        ! middle cloud amount
        ind_mid = high_btm<=p_full(i,j,:) .and. p_full(i,j,:)<=mid_btm
        ks = minloc(p_full(i,j,:), 1, mask=ind_mid)
        ke = maxloc(p_full(i,j,:), 1, mask=ind_mid)
        call random_overlap_single_column(cf(i,j,ks:ke), mid_ca(i,j))

        ! high cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        call random_overlap_single_column(cf(i,j,ks:ke), high_ca(i,j))
      enddo
    enddo
    tca = 1.0 - (1.0-high_ca)*(1.0-mid_ca)*(1.0-low_ca)

    if (do_test_overlap) then
      ! Diagnostics output
      if (id_tot_cld_amt_rnd > 0) then
        used = send_data (id_tot_cld_amt_rnd, tca*1.0e2, Time)
      endif
      if (id_high_cld_amt_rnd > 0) then
        used = send_data (id_high_cld_amt_rnd, high_ca*1.0e2, Time)
      endif
      if (id_mid_cld_amt_rnd > 0) then
        used = send_data (id_mid_cld_amt_rnd, mid_ca*1.0e2, Time)
      endif
      if (id_low_cld_amt_rnd > 0) then
        used = send_data (id_low_cld_amt_rnd, low_ca*1.0e2, Time)
      endif
    endif

  end subroutine diag_cldamt_random_overlap

  subroutine random_overlap_single_column(cf, cldamt)
    implicit none
    real, dimension(:), intent(in) :: cf
    real, intent(out) :: cldamt
    integer :: k

    ! Random overlap
    cldamt = 1.0
    do k=1,size(cf,1)
      cldamt = cldamt * (1-cf(k))
    end do
    cldamt = 1.0 - cldamt
  end subroutine random_overlap_single_column

  subroutine output_cldamt(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:) :: tca, high_ca, mid_ca, low_ca
    type(time_type),  intent(in)     :: Time
    real :: used

    if ( id_tot_cld_amt > 0 ) then
      used = send_data ( id_tot_cld_amt, tca, Time)
    endif
    if ( id_high_cld_amt > 0 ) then
      used = send_data ( id_high_cld_amt, high_ca, Time)
    endif
    if ( id_mid_cld_amt > 0 ) then
      used = send_data ( id_mid_cld_amt, mid_ca, Time)
    endif
    if ( id_low_cld_amt > 0 ) then
      used = send_data ( id_low_cld_amt, low_ca, Time)
    endif
  end subroutine output_cldamt

  ! ===================================================
  !             cloud cover diags end
  ! ===================================================
  subroutine cloud_cover_diags_end()

  end subroutine cloud_cover_diags_end

end module cloud_cover_diags_mod
