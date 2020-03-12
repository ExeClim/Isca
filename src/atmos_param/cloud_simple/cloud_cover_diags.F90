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

  logical :: adjust_top = .false.
  logical :: do_test_overlap = .false.
  logical :: do_cam_cld_cover_diag = .true.

  real :: cf_min = 1e-4

  ! ----- outputs for cloud amount diagnostics ----- !
  integer :: id_tot_cld_amt, id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt
  integer :: id_tot_cld_amt_mxr, id_high_cld_amt_mxr, id_mid_cld_amt_mxr, id_low_cld_amt_mxr
  integer :: id_tot_cld_amt_max, id_high_cld_amt_max, id_mid_cld_amt_max, id_low_cld_amt_max
  integer :: id_tot_cld_amt_rnd, id_high_cld_amt_rnd, id_mid_cld_amt_rnd, id_low_cld_amt_rnd
  integer :: id_tot_cld_amt_cam, id_high_cld_amt_cam, id_mid_cld_amt_cam, id_low_cld_amt_cam

  namelist /cloud_cover_diag_nml/ do_test_overlap, do_cam_cld_cover_diag, adjust_top, cf_min

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
                          'total cloud amount', 'percent')
    id_high_cld_amt = register_diag_field (mod_name, 'high_cld_amt', axes(1:2), Time, &
                          'high cloud amount', 'percent')
    id_mid_cld_amt = register_diag_field (mod_name, 'mid_cld_amt', axes(1:2), Time, &
                          'mid cloud amount', 'percent')
    id_low_cld_amt = register_diag_field (mod_name, 'low_cld_amt', axes(1:2), Time, &
                          'low cloud amount', 'percent')
    if (do_test_overlap) then
      id_tot_cld_amt_mxr = register_diag_field (mod_name, 'tot_cld_amt_mxr', axes(1:2), Time, &
              'total cloud amount', 'percent')
      id_high_cld_amt_mxr = register_diag_field (mod_name, 'high_cld_amt_mxr', axes(1:2), Time, &
              'high cloud amount', 'percent')
      id_mid_cld_amt_mxr = register_diag_field (mod_name, 'mid_cld_amt_mxr', axes(1:2), Time, &
              'mid cloud amount', 'percent')
      id_low_cld_amt_mxr = register_diag_field (mod_name, 'low_cld_amt_mxr', axes(1:2), Time, &
              'low cloud amount', 'percent')
      ! Max overlap
      id_tot_cld_amt_max = register_diag_field (mod_name, 'tot_cld_amt_max', axes(1:2), Time, &
              'total cloud amount', 'percent')
      id_high_cld_amt_max = register_diag_field (mod_name, 'high_cld_amt_max', axes(1:2), Time, &
              'high cloud amount', 'percent')
      id_mid_cld_amt_max = register_diag_field (mod_name, 'mid_cld_amt_max', axes(1:2), Time, &
              'mid cloud amount', 'percent')
      id_low_cld_amt_max = register_diag_field (mod_name, 'low_cld_amt_max', axes(1:2), Time, &
              'low cloud amount', 'percent')
      ! Random overlap
      id_tot_cld_amt_rnd = register_diag_field (mod_name, 'tot_cld_amt_rnd', axes(1:2), Time, &
              'total cloud amount', 'percent')
      id_high_cld_amt_rnd = register_diag_field (mod_name, 'high_cld_amt_rnd', axes(1:2), Time, &
              'high cloud amount', 'percent')
      id_mid_cld_amt_rnd = register_diag_field (mod_name, 'mid_cld_amt_rnd', axes(1:2), Time, &
              'mid cloud amount', 'percent')
      id_low_cld_amt_rnd = register_diag_field (mod_name, 'low_cld_amt_rnd', axes(1:2), Time, &
              'low cloud amount', 'percent')
    end if
    if (do_cam_cld_cover_diag) then
      id_tot_cld_amt_cam = register_diag_field (mod_name, 'tot_cld_amt_cam', axes(1:2), Time, &
              'total cloud amount', 'percent')
      id_high_cld_amt_cam = register_diag_field (mod_name, 'high_cld_amt_cam', axes(1:2), Time, &
              'high cloud amount', 'percent')
      id_mid_cld_amt_cam = register_diag_field (mod_name, 'mid_cld_amt_cam', axes(1:2), Time, &
              'mid cloud amount', 'percent')
      id_low_cld_amt_cam = register_diag_field (mod_name, 'low_cld_amt_cam', axes(1:2), Time, &
              'low cloud amount', 'percent')
    end if
  end subroutine cloud_cover_diags_init

  ! ===================================================
  !             cloud cover diags init
  ! ===================================================
  subroutine cloud_cover_diags(cf, p_full, p_half, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full, p_half
    type(time_type),   intent(in)       :: Time

    call diag_cloud_amount(cf, p_full, p_half, Time)
    if (do_test_overlap) then
      call diag_cldamt_maxrnd_overlap(cf, p_full, Time)
      call diag_cldamt_max_overlap(cf, p_full, Time)
      call diag_cldamt_random_overlap(cf, p_full, Time)
    end if
    if (do_cam_cld_cover_diag) then
      call cloud_cover_diags_cam(cf, p_full, p_half, Time)
    end if
    
  end subroutine cloud_cover_diags

  subroutine diag_cloud_amount(cf, p_full, p_half, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full, p_half
    type(time_type),   intent(in)       :: Time
    real,    dimension(size(cf,1),size(cf,2)) :: tca, high_ca, mid_ca, low_ca
    integer, dimension(size(cf,1),size(cf,2),size(cf,3)) :: ktop, kbot
    real,    dimension(size(cf,1),size(cf,2),size(cf,3)) :: cldamt, cloud
    integer, dimension(size(cf,1),size(cf,2)) :: nclds

    call max_rnd_overlap(cf, p_full, p_half, nclds, ktop, kbot, cldamt)
    call compute_tca_random(nclds, cldamt, tca)
    !call expand_cloud(nclds, ktop, kbot, cldamt, cloud)
    call compute_isccp_clds2(p_full, nclds, ktop, cldamt, high_ca, mid_ca, low_ca)

    ! Diagnostics output
    call output_cldamt(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cloud_amount

  subroutine max_rnd_overlap(cf, pfull, phalf, nclds, ktop, kbot, cldamt)
    !max_rnd_overlap returns various cloud specification properties
    !    obtained with the maximum-random overlap assumption.

    real,    dimension(:,:,:), intent(in)             :: cf, pfull, phalf
    integer, dimension(:,:),   intent(out)            :: nclds
    integer, dimension(:,:,:), intent(out)            :: ktop, kbot
    real,    dimension(:,:,:), intent(out)            :: cldamt

    ! local variables:
    real, dimension (size(cf,1), size(cf,2), size(cf,3))  :: cldamt_cs
    integer    :: kdim
    integer    :: top_t, bot_t
    integer    :: tmp_top, tmp_bot, nlev
    logical    :: already_in_cloud, cloud_bottom_reached
    real       :: maxcldfrac
    real       :: totcld_bot, max_bot
    real       :: totcld_top, max_top, tmp_val
    integer    :: i, j, k, kc, t

    kdim     = size(cf,3)
    nclds    = 0
    ktop     = 1
    kbot     = 1
    cldamt   = 0.0

    do j=1,size(cf,2)
      do i=1,size(cf,1)
        ! set a flag indicating that we are searching for the next cloud top.
        already_in_cloud  = .false.
        cloud_bottom_reached = .false.
        ! march down the column.
        do k=1,kdim
          ! find a layer containing cloud in the column.
          if (cf(i,j,k) .gt. cf_min) then
            if (.not. already_in_cloud)  then
              nclds(i,j) = nclds(i,j) + 1
              already_in_cloud = .true.
              cloud_bottom_reached = .false.
              ktop(i,j,nclds(i,j)) = k
              maxcldfrac = 0.0
            endif
            maxcldfrac = MAX(maxcldfrac, cf(i,j,k))
          endif

          if (cf(i,j,k) <= cf_min .and. already_in_cloud) then
            cloud_bottom_reached = .true.
            kbot(i,j,nclds(i,j)) = k - 1
          else if (already_in_cloud .and. k == kdim) then
            cloud_bottom_reached = .true.
            kbot(i,j,nclds(i,j)) = kdim
          endif
          !--------------------------------------------------------------------
          ! define the cloud fraction as the largest value of any layer in the cloud.
          !--------------------------------------------------------------------
          if (cloud_bottom_reached) then
            cldamt_cs(i,j,nclds(i,j)) = maxcldfrac
            !----------------------------------------------------------------------
            !    if adjust_top is true, the top and bottom indices of multi-layer
            !    clouds are adjusted to be those that are the most exposed to top
            !    and bottom view.
            !----------------------------------------------------------------------
            if (adjust_top) then
              ! define the cloud thickness.
              nlev = kbot(i,j,nclds(i,j)) - ktop(i,j,nclds(i,j)) + 1
              if (nlev > 1) then
                ! use the current top and bottom as the first guess for the new values.
                tmp_top = ktop(i,j,nclds(i,j))
                tmp_bot = kbot(i,j,nclds(i,j))
                ! initialize local search variables.
                totcld_bot = 0.
                totcld_top = 0.
                max_bot    = 0.
                max_top    = 0.

                do t=1,nlev
                  ! find adjusted cloud top.
                  top_t   = ktop(i,j,nclds(i,j)) + t - 1
                  tmp_val = MAX(0., cf(i,j,top_t) - totcld_top)
                  if (tmp_val > max_top) then
                    max_top = tmp_val
                    tmp_top = top_t
                  end if
                  totcld_top = totcld_top + tmp_val
                  ! find adjusted cloud base.
                  bot_t   = kbot(i,j,nclds(i,j)) - t + 1
                  tmp_val = MAX(0., cf(i,j,bot_t) - totcld_bot)
                  if (tmp_val > max_bot) then
                    max_bot = tmp_val
                    tmp_bot = bot_t
                  end if
                  totcld_bot = totcld_bot + tmp_val
                end do
                ! assign tmp_top and tmp_bot as the new ktop and kbot.
                ktop(i,j,nclds(i,j)) = tmp_top
                kbot(i,j,nclds(i,j)) = tmp_bot
              endif  !(nlev > 1)
            endif ! (adjust_top)
            !---------------------------------------------------------------------
            !    reset already_in_cloud and cloud_bottom_reached to indicate that
            !    the current cloud has been exited.
            !---------------------------------------------------------------------
            already_in_cloud     = .false.
            cloud_bottom_reached = .false.
          endif   ! (cloud_bottom_reached)
        end do
      end do
    end do
    !---------------------------------------------------------------------
    !    place cloud properties into physical-space arrays for return to
    !    calling routine. NOTE THAT ALL LEVELS IN A GIVEN CLOUD ARE
    !    ASSIGNED THE SAME PROPERTIES.
    !---------------------------------------------------------------------
    do j=1,size(cf,2)
      do i=1,size(cf,1)
        do kc=1, nclds(i,j)
          !cldamt(i,j,ktop(i,j,kc):kbot(i,j,kc)) = cldamt_cs(i,j,kc)
          cldamt(i,j,kc) = cldamt_cs(i,j,kc)
        end do
      end do
    end do
  end subroutine max_rnd_overlap

  subroutine compute_tca_random(nclds, cldamt, tca)
    ! This subroutine was adapted from AM4 src/atmos_param/clouds/clouds.F90
    integer, intent(in)  :: nclds (:,:)
    real,    intent(in)  :: cldamt(:,:,:)
    real,    intent(out) :: tca   (:,:)
    integer :: i, j, k

    !---- compute total cloud amount assuming that -----
    !       independent clouds overlap randomly
    tca = 1.0
    do i=1,size(cldamt,1)
      do j=1,size(cldamt,2)
        do k = 1,nclds(i,j)
          tca(i,j) = tca(i,j) * (1.0 - cldamt(i,j,k))
        enddo
      enddo
    enddo
    tca = (1.0 - tca) * 1.0e2 ! unit percent
  end subroutine compute_tca_random

  subroutine compute_isccp_clds2(pfull, nclds, ktop, cldamt, high_ca, mid_ca, low_ca)
    real,     dimension(:,:,:), intent(in)  :: pfull, cldamt
    integer,  dimension(:,:),   intent(in)  :: nclds
    integer,  dimension(:,:,:), intent(in)  :: ktop
    real,     dimension(:,:),   intent(out) :: high_ca, mid_ca, low_ca
    real,     parameter :: mid_btm = 6.8e4, high_btm = 4.4e4
    ! local array
    integer :: i, j, k, k_top

    high_ca = 1.0
    mid_ca  = 1.0
    low_ca  = 1.0

    do i=1,size(cldamt,1)
      do j=1,size(cldamt,2)
        do k = 1,nclds(i,j)
          k_top = ktop(i,j,k)
          if (pfull(i,j,k_top)>mid_btm) then
            low_ca(i,j) = low_ca(i,j) * (1.0 - cldamt(i,j,k))
          else if (pfull(i,j,k_top)<high_btm) then
            high_ca(i,j) = high_ca(i,j) * (1.0 - cldamt(i,j,k))
          else
            mid_ca(i,j) = mid_ca(i,j) * (1.0 - cldamt(i,j,k))
          end if
        enddo
      enddo
    enddo

    low_ca  = (1.0 - low_ca)  * 1.0e2 ! unit percent
    mid_ca  = (1.0 - mid_ca)  * 1.0e2
    high_ca = (1.0 - high_ca) * 1.0e2

  end subroutine compute_isccp_clds2

  subroutine cldovrlap(pint, cld, nmxrgn, pmxrgn)
    !subroutine cldovrlap(lchnk   ,ncol    ,pint    ,cld     ,nmxrgn  ,pmxrgn  )
    !  This code is borrowed from CESM.
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
    !
    ! Input arguments
    !
    !integer, intent(in) :: ncol                ! number of atmospheric columns
    real, intent(in), dimension(:,:,:) :: pint  ! Interface pressure
    real, intent(in), dimension(:,:,:) :: cld   ! Fractional cloud cover
    !
    ! Output arguments
    !
    ! Number of maximally overlapped regions
    integer,  intent(out), dimension(size(cld,1), size(cld,2)) :: nmxrgn
    real, intent(out), dimension(size(pint,1), size(pint,2), size(pint,3)) :: pmxrgn
    ! Maximum values of pressure for each maximally overlapped region.
    !    0->pmxrgn(i,1) is range of pressure for
    !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
    !    2nd region, etc
    !
    !---------------------------Local variables-----------------------------
    !
    integer :: i, j  ! Lat/Longitude index
    integer :: k     ! Level index
    integer :: n     ! Max-overlap region counter

    real, dimension(size(pint,1), size(pint,2), size(pint,3)) :: pnm  ! Interface pressure

    logical :: cld_found                          ! Flag for detection of cloud
    logical, dimension(size(cld,3)) :: cld_layer  ! Flag for cloud in layer
    !
    !------------------------------------------------------------------------
    !
    integer :: pver, pverp

    pver = size(cld,3)
    pverp = pver + 1

    do i = 1,size(cld,1)
      do j = 1,size(cld,2)
        cld_found = .false.
        cld_layer(:) = cld(i,j,:) > 0.0
        pmxrgn(i,j,:) = 0.0
        pnm(i,j,:) = pint(i,j,:) !* 10.0  ! why multiplied by 10?
        n = 1
        do k = 1, pver
          if (cld_layer(k) .and.  .not. cld_found) then
              cld_found = .true.
          else if ( .not. cld_layer(k) .and. cld_found) then
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

  subroutine cloud_cover_diags_cam(cld, pmid, pint, Time)
    !!! pmid is pfull, and pint is phalf
    real,    intent(in), dimension(:,:,:) :: cld, pmid, pint
    ! Total, low, middle and high random overlap cloud cover
    real, dimension(size(cld,1), size(cld,2)) :: cldtot, cldlow, cldmed, cldhgh

    type(time_type),  intent(in)  :: Time
    real :: used

    !---------------------------Local workspace-----------------------------
    integer :: i,j,k       ! lat/lon,level indices
    integer, dimension(size(cld,1), size(cld,2)) :: irgn    ! Max-overlap region index
    integer :: max_nmxrgn  ! maximum value of nmxrgn over columns
    integer :: ityp        ! Type counter
    real, dimension(size(cld,1), size(cld,2)) :: clrsky      ! Max-random clear sky fraction
    real, dimension(size(cld,1), size(cld,2)) :: clrskymax   ! Maximum overlap clear sky fraction

    !------------------------------Parameters-------------------------------
    real, parameter :: plowmax=1.2e5, plowmin=7.0e4  ! Max/min prs for low cloud cover range
    real, parameter :: pmedmax=7.0e4, pmedmin=4.0e4  ! Max/min prs for mid cloud cover range
    real, parameter :: phghmax=4.0e4, phghmin=5.0e3  ! Max/min prs for hgh cloud cover range

    real, dimension(4) :: ptypmin
    real, dimension(4) :: ptypmax

    data ptypmin /phghmin, plowmin, pmedmin, phghmin/
    data ptypmax /plowmax, plowmax, pmedmax, phghmax/

    integer, dimension(size(cld,1), size(cld,2)) :: nmxrgn
    real, dimension(size(pint,1), size(pint,2), size(pint,3)) :: pmxrgn

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
      do k=1,size(cld,3)
        do i=1,size(cld,1)
          do j=1,size(cld,2)
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

      if (ityp == 1) cldtot = 1.0 - (clrsky * clrskymax)
      if (ityp == 2) cldlow = 1.0 - (clrsky * clrskymax)
      if (ityp == 3) cldmed = 1.0 - (clrsky * clrskymax)
      if (ityp == 4) cldhgh = 1.0 - (clrsky * clrskymax)
    end do

    ! Write the output diagnostics
    if ( id_tot_cld_amt_cam > 0 ) then
      used = send_data ( id_tot_cld_amt_cam, cldtot*1e2, Time)
    endif
    if ( id_high_cld_amt_cam > 0 ) then
      used = send_data ( id_high_cld_amt_cam, cldlow*1e2, Time)
    endif
    if ( id_mid_cld_amt_cam > 0 ) then
      used = send_data ( id_mid_cld_amt_cam, cldmed*1e2, Time)
    endif
    if ( id_low_cld_amt_cam > 0 ) then
      used = send_data ( id_low_cld_amt_cam, cldhgh*1e2, Time)
    endif

  end subroutine cloud_cover_diags_cam

  subroutine diag_cldamt_maxrnd_overlap(cf, p_full, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full
    type(time_type),   intent(in)       :: Time
    real, dimension(size(cf,1), size(cf,2)) :: tca, high_ca, mid_ca, low_ca
    integer :: i, j, ks, ke !, ks_mid, ke_mid
    logical, dimension(size(cf,3)) :: ind_mid
    real :: mid_btm = 7e4, high_btm = 4e4

    tca = 1.0
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    do i=1,size(cf,1)
      do j=1,size(cf,2)
        ! total cloud amount
        !ks = 1
        !ke = size(cf,3)
        !call max_rnd_overlap_single_lev(cf(i,j,ks:ke), p_full(i,j,ks:ke), tca(i,j))

        ! high cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        call max_rnd_overlap_single_lev(cf(i,j,ks:ke), p_full(i,j,ks:ke), high_ca(i,j))
        !ks_mid = ke + 1

        ! low cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        ke = size(cf,3) !maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        call max_rnd_overlap_single_lev(cf(i,j,ks:ke), p_full(i,j,ks:ke), low_ca(i,j))
        !ke_mid = ks - 1

        ! middle cloud amount
        !ks = ks_mid
        !ke = ke_mid
        ind_mid = high_btm<=p_full(i,j,:) .and. p_full(i,j,:)<=mid_btm
        ks = minloc(p_full(i,j,:), 1, mask=ind_mid)
        ke = maxloc(p_full(i,j,:), 1, mask=ind_mid)
        call max_rnd_overlap_single_lev(cf(i,j,ks:ke), p_full(i,j,ks:ke), mid_ca(i,j))
      enddo
    enddo
    tca = 1.0 - (1.0-high_ca)*(1.0-mid_ca)*(1.0-low_ca)

    ! Diagnostics output
    call output_cldamt_max_random(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cldamt_maxrnd_overlap

  subroutine diag_cldamt_max_overlap(cf, p_full, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full
    type(time_type),   intent(in)       :: Time
    real, dimension(size(cf,1), size(cf,2)) :: tca, high_ca, mid_ca, low_ca
    integer :: i, j, ks, ke
    logical, dimension(size(cf,3)) :: ind_mid
    real :: mid_btm = 7e4, high_btm = 4e4
    tca = 1.0
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    ! total cld amount
    tca = maxval(cf, 3)

    do i=1,size(cf,1)
      do j=1,size(cf,2)
        ! high cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        high_ca(i,j) = maxval(cf(i,j,ks:ke), 1)

        ! low cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        low_ca(i,j) = maxval(cf(i,j,ks:ke), 1)

        ! middle cloud amount
        ind_mid = high_btm<=p_full(i,j,:) .and. p_full(i,j,:)<=mid_btm
        ks = minloc(p_full(i,j,:), 1, mask=ind_mid)
        ke = maxloc(p_full(i,j,:), 1, mask=ind_mid)
        mid_ca(i,j) = maxval(cf(i,j,ks:ke), 1)
      enddo
    enddo

    ! Diagnostics output
    call output_cldamt_max(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cldamt_max_overlap

  subroutine diag_cldamt_random_overlap(cf, p_full, Time)
    real, intent(in),  dimension(:,:,:) :: cf, p_full
    type(time_type),   intent(in)       :: Time
    real, dimension(size(cf,1), size(cf,2)) :: tca, high_ca, mid_ca, low_ca
    integer :: i, j, ks, ke
    logical, dimension(size(cf,3)) :: ind_mid
    real :: mid_btm = 7e4, high_btm = 4e4

    tca = 1.0
    high_ca = 1.0
    mid_ca = 1.0
    low_ca = 1.0

    do i=1,size(cf,1)
      do j=1,size(cf,2)
        ! total cloud amount
        !ks = 1
        !ke = size(cf,3)
        !call random_overlap_single_lev(cf(i,j,ks:ke), tca(i,j))

        ! high cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)<high_btm)
        call random_overlap_single_lev(cf(i,j,ks:ke), high_ca(i,j))

        ! low cloud amount
        ks = minloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        ke = maxloc(p_full(i,j,:), 1, mask=p_full(i,j,:)>mid_btm)
        call random_overlap_single_lev(cf(i,j,ks:ke), low_ca(i,j))

        ! middle cloud amount
        ind_mid = high_btm<=p_full(i,j,:) .and. p_full(i,j,:)<=mid_btm
        ks = minloc(p_full(i,j,:), 1, mask=ind_mid)
        ke = maxloc(p_full(i,j,:), 1, mask=ind_mid)
        call random_overlap_single_lev(cf(i,j,ks:ke), mid_ca(i,j))
      enddo
    enddo
    tca = 1.0 - (1.0-high_ca)*(1.0-mid_ca)*(1.0-low_ca)

    ! Diagnostics output
    call output_cldamt_random(tca, high_ca, mid_ca, low_ca, Time)
  end subroutine diag_cldamt_random_overlap

  subroutine max_rnd_overlap_single_lev(cf, pfull, cldamt)
    implicit none
    real, dimension(:), intent(in) :: cf, pfull
    real, intent(out) :: cldamt
    ! local variables:
    integer, dimension(size(cf,1)) :: ktop, kbot
    real,    dimension(size(cf,1)) :: cldamt_cs
    integer :: nclds, kdim, k
    logical :: already_in_cloud, cloud_bottom_reached
    real    :: maxcldfrac

    kdim  = size(cf,1)
    nclds = 0
    ktop  = 1
    kbot  = 1
    maxcldfrac = 0.0

    ! set a flag indicating that we are searching for the next cloud top.
    already_in_cloud = .false.
    cloud_bottom_reached = .false.
    ! march down the column.
    do k=1,kdim
      ! find a layer containing cloud in the column.
      if (cf(k) .gt. cf_min) then
        if (.not. already_in_cloud) then
          nclds = nclds + 1
          already_in_cloud = .true.
          cloud_bottom_reached = .false.
          ktop(nclds) = k
          maxcldfrac = 0.0
        endif
        maxcldfrac = MAX(maxcldfrac, cf(k))
      endif

      if (cf(k) <= cf_min .and. already_in_cloud) then
        cloud_bottom_reached = .true.
        kbot(nclds) = k - 1
      else if (already_in_cloud .and. k == kdim) then
        cloud_bottom_reached = .true.
        kbot(nclds) = kdim
      endif

      if (cloud_bottom_reached) then
        cldamt_cs(nclds) = maxcldfrac
        already_in_cloud = .false.
        cloud_bottom_reached = .false.
      endif ! (cloud_bottom_reached)
    end do

    ! Random overlap
    cldamt = 1.0
    do k=1, nclds
      cldamt = cldamt * (1-cldamt_cs(k))
    end do
    cldamt = 1.0 - cldamt
  end subroutine max_rnd_overlap_single_lev

  subroutine random_overlap_single_lev(cf, cldamt)
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
  end subroutine random_overlap_single_lev

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

  subroutine output_cldamt_max_random(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:) :: tca, high_ca, mid_ca, low_ca
    type(time_type),  intent(in)     :: Time
    real :: used

    if ( id_tot_cld_amt > 0 ) then
      used = send_data ( id_tot_cld_amt_mxr, tca*1e2, Time)
    endif
    if ( id_high_cld_amt > 0 ) then
      used = send_data ( id_high_cld_amt_mxr, high_ca*1e2, Time)
    endif
    if ( id_mid_cld_amt > 0 ) then
      used = send_data ( id_mid_cld_amt_mxr, mid_ca*1e2, Time)
    endif
    if ( id_low_cld_amt > 0 ) then
      used = send_data ( id_low_cld_amt_mxr, low_ca*1e2, Time)
    endif
  end subroutine output_cldamt_max_random

  subroutine output_cldamt_max(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:) :: tca, high_ca, mid_ca, low_ca
    type(time_type),  intent(in)     :: Time
    real :: used

    if ( id_tot_cld_amt_max > 0 ) then
      used = send_data ( id_tot_cld_amt_max, tca*1e2, Time)
    endif
    if ( id_high_cld_amt_max > 0 ) then
      used = send_data ( id_high_cld_amt_max, high_ca*1e2, Time)
    endif
    if ( id_mid_cld_amt_max > 0 ) then
      used = send_data ( id_mid_cld_amt_max, mid_ca*1e2, Time)
    endif
    if ( id_low_cld_amt_max > 0 ) then
      used = send_data ( id_low_cld_amt_max, low_ca*1e2, Time)
    endif
  end subroutine output_cldamt_max

  subroutine output_cldamt_random(tca, high_ca, mid_ca, low_ca, Time)
    real, intent(in), dimension(:,:) :: tca, high_ca, mid_ca, low_ca
    type(time_type),  intent(in)     :: Time
    real :: used

    if ( id_tot_cld_amt_rnd > 0 ) then
      used = send_data ( id_tot_cld_amt_rnd, tca*1e2, Time)
    endif
    if ( id_high_cld_amt_rnd > 0 ) then
      used = send_data ( id_high_cld_amt_rnd, high_ca*1e2, Time)
    endif
    if ( id_mid_cld_amt_rnd > 0 ) then
      used = send_data ( id_mid_cld_amt_rnd, mid_ca*1e2, Time)
    endif
    if ( id_low_cld_amt_rnd > 0 ) then
      used = send_data ( id_low_cld_amt_rnd, low_ca*1e2, Time)
    endif
  end subroutine output_cldamt_random

  ! ===================================================
  !             cloud cover diags end
  ! ===================================================
  subroutine cloud_cover_diags_end()
    ! nothing
  end subroutine cloud_cover_diags_end

end module cloud_cover_diags_mod
