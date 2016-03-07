
! -*-f90-*-
! $Id: mpp_domains_define.inc,v 18.0 2010/03/02 23:57:15 fms Exp $


  ! <SUBROUTINE NAME="mpp_define_layout2D" INTERFACE="mpp_define_layout">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"></IN>
  !  <IN NAME="ndivs" TYPE="integer"></IN>
  !  <OUT NAME="layout" TYPE="integer" DIM="(2)"></OUT>
  ! </SUBROUTINE>
  subroutine mpp_define_layout2D( global_indices, ndivs, layout )
    integer, intent(in) :: global_indices(:) !(/ isg, ieg, jsg, jeg /)
    integer, intent(in) :: ndivs !number of divisions to divide global domain
    integer, intent(out) :: layout(:)

    integer :: isg, ieg, jsg, jeg, isz, jsz, idiv, jdiv

    if(size(global_indices(:)) .NE. 4) call mpp_error(FATAL,"mpp_define_layout2D: size of global_indices should be 4")
    if(size(layout(:)) .NE. 2) call mpp_error(FATAL,"mpp_define_layout2D: size of layout should be 2")

    isg = global_indices(1)
    ieg = global_indices(2)
    jsg = global_indices(3)
    jeg = global_indices(4)

    isz = ieg - isg + 1
    jsz = jeg - jsg + 1
    !first try to divide ndivs in the domain aspect ratio: if imperfect aspect, reduce idiv till it divides ndivs
    idiv = nint( sqrt(float(ndivs*isz)/jsz) )
    idiv = max(idiv,1) !for isz=1 line above can give 0
    do while( mod(ndivs,idiv).NE.0 )
       idiv = idiv - 1
    end do                 !will terminate at idiv=1 if not before
    jdiv = ndivs/idiv

    layout = (/ idiv, jdiv /)
    return
  end subroutine mpp_define_layout2D

  !############################################################################
  ! <SUBROUTINE NAME="mpp_define_mosaic_pelist">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"></IN>
  !  <IN NAME="pelist" TYPE="integer" DIM="(0:)">  </IN>
  ! </SUBROUTINE>
  ! NOTE: The following routine may need to revised to improve the capability.
  !  It is very hard to make it balance for all the situation.
  !  Hopefully some smart idea will come up someday.
  subroutine mpp_define_mosaic_pelist( sizes, pe_start, pe_end, pelist, costpertile)
    integer, dimension(:), intent(in)           :: sizes
    integer, dimension(:), intent(inout)        :: pe_start, pe_end
    integer, dimension(:), intent(in), optional :: pelist, costpertile
    integer, dimension(size(sizes(:)))          :: costs
    integer, dimension(:), allocatable          :: pes
    integer                                     :: ntiles, npes, totcosts, avgcost
    integer                                     :: ntiles_left, npes_left, pos, n, tile
    integer                                     :: cost_on_tile, cost_on_pe, npes_used, errunit

    ntiles = size(sizes(:))
    if(size(pe_start(:)) .NE. ntiles .OR. size(pe_end(:)) .NE. ntiles ) then
       call mpp_error(FATAL, "mpp_define_mosaic_pelist: size mismatch between pe_start/pe_end and sizes")
    end if

    if(present(costpertile)) then
       if(size(costpertile(:)) .NE. ntiles ) then
          call mpp_error(FATAL, "mpp_define_mosaic_pelist: size mismatch between costpertile and sizes")
       end if
       costs = sizes*costpertile
    else
       costs = sizes
    end if

    if( PRESENT(pelist) )then
       if( .NOT.any(pelist.EQ.mpp_pe()) )then
          errunit = stderr()
          write( errunit,* )'pe=', mpp_pe(), ' pelist=', pelist
          call mpp_error( FATAL, 'mpp_define_mosaic_pelist: pe must be in pelist.' )
       end if
       npes = size(pelist(:))
       allocate( pes(0:npes-1) )
       pes(:) = pelist(:)
    else
       npes = mpp_npes()
       allocate( pes(0:npes-1) )
       call mpp_get_current_pelist(pes)
    end if

    ntiles_left = ntiles
    npes_left = npes
    pos = pes(0)

    do while( ntiles_left > 0 )
       if( npes_left == 1 ) then ! all left tiles will on the last processor, imbalance possibly.
          do n = 1, ntiles
             if(costs(n) > 0) then
                pe_start(n) = pos
                pe_end(n) = pos
                costs(n)  = 0
             end if
          end do
          ntiles_left = 0
          npes_left = 0
       else
          totcosts = sum(costs)   
          avgcost  = CEILING(real(totcosts)/npes_left )
          tile = minval(maxloc(costs))
          cost_on_tile = costs(tile)  
          pe_start(tile) = pos
          ntiles_left = ntiles_left - 1
          costs(tile) = 0
          totcosts = totcosts - cost_on_tile
          if(cost_on_tile .GE. avgcost ) then
             npes_used = min(ceiling(real(cost_on_tile)/avgcost), npes_left)
             if( ntiles_left > 0 .AND. npes_used == npes_left ) npes_used = npes_used - 1
             pe_end(tile) = pos + npes_used - 1
             npes_left = npes_left - npes_used
             pos = pos + npes_used
          else
             !--- find other tiles to share the pe
             pe_end(tile) = pos
             cost_on_pe = cost_on_tile
             do while(ntiles_left>npes_left)  ! make sure all the pes are used.
                tile = minval(minloc(costs, costs> 0 ))
                cost_on_tile = costs(tile)
                cost_on_pe = cost_on_pe + cost_on_tile
                if(cost_on_pe > avgcost ) exit
                pe_start(tile) = pos
                pe_end(tile) = pos
                ntiles_left = ntiles_left - 1
                costs(tile) = 0
                totcosts = totcosts - cost_on_tile
             end do
             npes_left = npes_left - 1
             pos = pos + 1
          end if
       end if
    end do

    if(npes_left .NE. 0 ) call mpp_error(FATAL, "mpp_define_mosaic_pelist: the left npes should be zero")
    deallocate(pes)

  end subroutine mpp_define_mosaic_pelist

  !#####################################################################
  subroutine mpp_compute_extent(isg,ieg,ndivs,ibegin,iend, extent )
    integer,                 intent(in)          :: isg, ieg, ndivs
    integer, dimension(0:), intent(out)          :: ibegin, iend
    integer, dimension(0:), intent(in), optional :: extent

    integer :: ndiv, imax, ndmax, ndmirror
    integer :: is, ie, n
    logical :: symmetrize, use_extent
    !statement functions
    logical :: even, odd
    even(n) = (mod(n,2).EQ.0)
    odd (n) = (mod(n,2).EQ.1)

    use_extent = .false.
    if(PRESENT(extent)) then
       if( size(extent(:)).NE.ndivs ) &
            call mpp_error( FATAL, 'mpp_compute_extent: extent array size must equal number of domain divisions.' )       
       use_extent = .true.
       if(ALL(extent ==0)) use_extent = .false.
    endif

    is = isg
    if(use_extent) then
       ibegin(0) = isg
       do ndiv = 0, ndivs-2
          if(extent(ndiv) .LE. 0) call mpp_error( FATAL, 'mpp_compute_extent: domain extents must be positive definite.' )
          iend(ndiv) = ibegin(ndiv) + extent(ndiv) - 1
          ibegin(ndiv+1) = iend(ndiv) + 1
       enddo
       iend(ndivs-1) = ibegin(ndivs-1) + extent(ndivs-1) - 1
       if(iend(ndivs-1) .NE. ieg) call mpp_error(FATAL, 'mpp_compute_extent: extent array limits do not match global domain.' )
    else
       do ndiv=0,ndivs-1
          !modified for mirror-symmetry
          !original line
          !                 ie = is + CEILING( float(ieg-is+1)/(ndivs-ndiv) ) - 1

          !problem of dividing nx points into n domains maintaining symmetry
          !i.e nx=18 n=4 4554 and 5445 are solutions but 4455 is not.
          !this will always work for nx even n even or odd
          !this will always work for nx odd, n odd
          !this will never  work for nx odd, n even: for this case we supersede the mirror calculation
          !                 symmetrize = .NOT. ( mod(ndivs,2).EQ.0 .AND. mod(ieg-isg+1,2).EQ.1 )
          !nx even n odd fails if n>nx/2
          symmetrize = ( even(ndivs) .AND. even(ieg-isg+1) ) .OR. &
               (  odd(ndivs) .AND.  odd(ieg-isg+1) ) .OR. &
               (  odd(ndivs) .AND. even(ieg-isg+1) .AND. ndivs.LT.(ieg-isg+1)/2 )

          !mirror domains are stored in the list and retrieved if required.
          if( ndiv.EQ.0 )then
             !initialize max points and max domains
             imax = ieg
             ndmax = ndivs
          end if
          !do bottom half of decomposition, going over the midpoint for odd ndivs
          if( ndiv.LT.(ndivs-1)/2+1 )then
             !domain is sized by dividing remaining points by remaining domains
             ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
             ndmirror = (ndivs-1) - ndiv !mirror domain
             if( ndmirror.GT.ndiv .AND. symmetrize )then !only for domains over the midpoint
                !mirror extents, the max(,) is to eliminate overlaps
                ibegin(ndmirror) = max( isg+ieg-ie, ie+1 )
                iend(ndmirror)   = max( isg+ieg-is, ie+1 )
                imax = ibegin(ndmirror) - 1
                ndmax = ndmax - 1
             end if
          else
             if( symmetrize )then
                !do top half of decomposition by retrieving saved values
                is = ibegin(ndiv)
                ie = iend(ndiv)
             else
                ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
             end if
          end if
          ibegin(ndiv) = is
          iend(ndiv) = ie
          if( ie.LT.is )call mpp_error( FATAL,  &
               'MPP_DEFINE_DOMAINS(mpp_compute_extent): domain extents must be positive definite.' )
          if( ndiv.EQ.ndivs-1 .AND. iend(ndiv).NE.ieg ) &
               call mpp_error( FATAL, 'mpp_compute_extent: domain extents do not span space completely.' )
          is = ie + 1
       end do
    endif


  end subroutine mpp_compute_extent

  !#####################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !              MPP_DEFINE_DOMAINS: define layout and decomposition            !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! <SUBROUTINE NAME="mpp_define_domains1D" INTERFACE="mpp_define_domains">>
  !   <IN NAME="global_indices" TYPE="integer" DIM="(2)"> </IN>
  !   <IN NAME="ndivs" TYPE="integer">  </IN>
  !   <INOUT NAME="domain" TYPE="type(domain1D)"> </INOUT>
  !   <IN NAME="pelist" TYPE="integer" DIM="(0:)">  </IN>
  !   <IN NAME="flags" TYPE="integer">  </IN>
  !   <IN NAME="halo" TYPE="integer">  </IN>
  !   <IN NAME="extent" TYPE="integer" DIM="(0:)">  </IN>
  !   <IN NAME="maskmap" TYPE="logical" DIM="(0:)"> </IN>
  ! </SUBROUTINE>
    !routine to divide global array indices among domains, and assign domains to PEs
    !domain is of type domain1D
    !ARGUMENTS:
    !      global_indices(2)=(isg,ieg) gives the extent of global domain
    !      ndivs is number of divisions of domain: even divisions unless extent is present.
    !      domain is the returned domain1D
    !      pelist (optional) list of PEs to which domains are to be assigned (default 0...npes-1)
    !                 size of pelist must correspond to number of mask=.TRUE. divisions
    !      flags define whether compute and data domains are global (undecomposed) and whether global domain has periodic boundaries
    !      halo (optional) defines halo width (currently the same on both sides)
    !      extent (optional) array defines width of each division (used for non-uniform domain decomp, for e.g load-balancing)
    !      maskmap (optional) a division whose maskmap=.FALSE. is not assigned to any domain
    !  By default we assume decomposition of compute and data domains, non-periodic boundaries, no halo, as close to uniform extents
    !  as the input parameters permit
  subroutine mpp_define_domains1D( global_indices, ndivs, domain, pelist, flags, halo, extent, maskmap, &
                                   memory_size, begin_halo, end_halo )
    integer,           intent(in) :: global_indices(:) !(/ isg, ieg /)
    integer,           intent(in) :: ndivs
    type(domain1D), intent(inout) :: domain !declared inout so that existing links, if any, can be nullified
    integer, intent(in), optional :: pelist(0:)
    integer, intent(in), optional :: flags, halo
    integer, intent(in), optional :: extent(0:)
    logical, intent(in), optional :: maskmap(0:)
    integer, intent(in), optional :: memory_size
    integer, intent(in), optional :: begin_halo, end_halo

    logical              :: compute_domain_is_global, data_domain_is_global
    integer              :: ndiv, n, isg, ieg, i
    integer, allocatable :: pes(:)
    integer              :: ibegin(0:ndivs-1), iend(0:ndivs-1)
    logical              :: mask(0:ndivs-1)
    integer              :: halosz, halobegin, haloend
    integer              :: errunit

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )
    if(size(global_indices(:)) .NE. 2) call mpp_error(FATAL,"mpp_define_domains1D: size of global_indices should be 2")
    !get global indices
    isg = global_indices(1)
    ieg = global_indices(2)
    if( ndivs.GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more divisions requested than rows available.' )
    !get the list of PEs on which to assign domains; if pelist is absent use 0..npes-1
    if( PRESENT(pelist) )then
       if( .NOT.any(pelist.EQ.mpp_pe()) )then
          errunit = stderr()
          write( errunit,* )'pe=', mpp_pe(), ' pelist=', pelist
          call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: pe must be in pelist.' )
       end if
       allocate( pes(0:size(pelist(:))-1) )
       pes(:) = pelist(:)
    else
       allocate( pes(0:mpp_npes()-1) )
       pes(:) = (/ (i,i=0,mpp_npes()-1) /)
    end if

    !get number of real domains: 1 mask domain per PE in pes
    mask = .TRUE.                 !default mask
    if( PRESENT(maskmap) )then
       if( size(maskmap(:)).NE.ndivs ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: maskmap array size must equal number of domain divisions.' )
       mask(:) = maskmap(:)
    end if
    if( count(mask).NE.size(pes(:)) ) &
         call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: number of TRUEs in maskmap array must match PE count.' )

    !get halosize
    halosz = 0
    if( PRESENT(halo) ) then
       halosz = halo
       !--- if halo is present, begin_halo and end_halo should not present
       if(present(begin_halo) .OR. present(end_halo) ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: when halo is present, begin_halo and end_halo should not present")
    end if
    halobegin = halosz; haloend = halosz
    if(present(begin_halo)) halobegin = begin_halo
    if(present(end_halo))   haloend   = end_halo
    halosz = max(halobegin, haloend)
    !get flags
    compute_domain_is_global = .FALSE.
    data_domain_is_global    = .FALSE.
    domain%cyclic = .FALSE.
    domain%goffset = 1
    domain%loffset = 1
    if( PRESENT(flags) )then
       !NEW: obsolete flag global_compute_domain, since ndivs is non-optional and you cannot have global compute and ndivs.NE.1 
       compute_domain_is_global = ndivs.EQ.1
       !if compute domain is global, data domain must also be
       data_domain_is_global    = BTEST(flags,GLOBAL) .OR. compute_domain_is_global
       domain%cyclic  = BTEST(flags,CYCLIC) .AND. halosz.NE.0
       if(BTEST(flags,CYCLIC)) domain%goffset = 0
    end if

    !set up links list
    allocate( domain%list(0:ndivs-1) )

    !set global domain
    domain%list(:)%global%begin     = isg
    domain%list(:)%global%end       = ieg
    domain%list(:)%global%size      = ieg-isg+1
    domain%list(:)%global%max_size  = ieg-isg+1
    domain%list(:)%global%is_global = .TRUE. !always

    !get compute domain
    if( compute_domain_is_global )then
       domain%list(:)%compute%begin = isg
       domain%list(:)%compute%end   = ieg
       domain%list(:)%compute%is_global = .TRUE.
       domain%list(:)%pe = pes(:)
       domain%pos = 0
    else
       domain%list(:)%compute%is_global = .FALSE.
       n = 0
       call mpp_compute_extent(isg, ieg, ndivs, ibegin, iend, extent)
       do ndiv=0,ndivs-1
          domain%list(ndiv)%compute%begin = ibegin(ndiv)
          domain%list(ndiv)%compute%end   = iend(ndiv)
          if( mask(ndiv) )then
             domain%list(ndiv)%pe = pes(n)
             if( mpp_pe().EQ.pes(n) )domain%pos = ndiv
             n = n + 1
          else
             domain%list(ndiv)%pe = NULL_PE
          end if
       end do
    end if

    domain%list(:)%compute%size  = domain%list(:)%compute%end - domain%list(:)%compute%begin + 1

    !get data domain
    !data domain is at least equal to compute domain
    domain%list(:)%data%begin = domain%list(:)%compute%begin
    domain%list(:)%data%end   = domain%list(:)%compute%end
    domain%list(:)%data%is_global = .FALSE.
    !apply global flags
    if( data_domain_is_global )then
       domain%list(:)%data%begin  = isg
       domain%list(:)%data%end    = ieg
       domain%list(:)%data%is_global = .TRUE.
    end if
    !apply margins
    domain%list(:)%data%begin = domain%list(:)%data%begin - halobegin
    domain%list(:)%data%end   = domain%list(:)%data%end   + haloend
    domain%list(:)%data%size  = domain%list(:)%data%end - domain%list(:)%data%begin + 1

    !--- define memory domain, if memory_size is not present or memory size is 0, memory domain size
    !--- will be the same as data domain size. if momory_size is present, memory_size should greater than
    !--- or equal to data size. The begin of memory domain will be always the same as data domain.
    domain%list(:)%memory%begin = domain%list(:)%data%begin
    domain%list(:)%memory%end   = domain%list(:)%data%end
    if( present(memory_size) ) then
       if(memory_size > 0) then
          if( domain%list(domain%pos)%data%size > memory_size ) call mpp_error(FATAL, &
               "mpp_domains_define.inc: data domain size is larger than memory domain size on this pe")
          domain%list(:)%memory%end   = domain%list(:)%memory%begin + memory_size - 1
       end if
    end if
    domain%list(:)%memory%size = domain%list(:)%memory%end - domain%list(:)%memory%begin + 1
    domain%list(:)%memory%is_global = domain%list(:)%data%is_global

    domain%compute = domain%list(domain%pos)%compute
    domain%data    = domain%list(domain%pos)%data
    domain%global  = domain%list(domain%pos)%global
    domain%memory  = domain%list(domain%pos)%memory
    domain%compute%max_size = MAXVAL( domain%list(:)%compute%size )
    domain%data%max_size    = MAXVAL( domain%list(:)%data%size )
    domain%global%max_size  = domain%global%size
    domain%memory%max_size  = domain%memory%size

    !PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
    deallocate( pes )
    return

  end subroutine mpp_define_domains1D

  !################################################################################
  !--- define the IO domain. 
  subroutine mpp_define_io_domain(domain, io_layout)
    type(domain2D), intent(inout) :: domain
    integer,        intent(in   ) :: io_layout(2)
    integer                       :: layout(2)
    integer                       :: npes_in_group
    type(domain2D), pointer       :: io_domain=>NULL() 
    integer                       :: i, j, n, m
    integer                       :: ipos, jpos, igroup, jgroup
    integer                       :: ipos_beg, ipos_end, jpos_beg, jpos_end
    integer                       :: whalo, ehalo, shalo, nhalo
    integer                       :: npes_x, npes_y

    if(io_layout(1) * io_layout(2) .LE. 0) then
       call mpp_error(NOTE,  &
         "mpp_domains_define.inc(mpp_define_io_domain): io domain will not be defined for "//trim(domain%name)// &
         " when one or both entry of io_layout is not positive")
       return
    endif

    layout(1) = size(domain%x(1)%list(:))
    layout(2) = size(domain%y(1)%list(:)) 

    if(ASSOCIATED(domain%io_domain)) call mpp_error(FATAL, &
       "mpp_domains_define.inc(mpp_define_io_domain): io_domain is already defined")

    if(mod(layout(1), io_layout(1)) .NE. 0) call mpp_error(FATAL, &
       "mpp_domains_define.inc(mpp_define_io_domain): "//trim(domain%name)//" domain layout(1) must be divided by io_layout(1)")    
    if(mod(layout(2), io_layout(2)) .NE. 0) call mpp_error(FATAL, &
       "mpp_domains_define.inc(mpp_define_io_domain): "//trim(domain%name)//" domain layout(2) must be divided by io_layout(2)")
    if(size(domain%x(:)) > 1) call mpp_error(FATAL, &
       "mpp_domains_define.inc(mpp_define_io_domain): "//trim(domain%name)// &
       ": multiple tile per pe is not supported yet for this routine")
 
    allocate(domain%io_domain)
    domain%io_layout = io_layout
    io_domain => domain%io_domain
    npes_x = layout(1)/io_layout(1)
    npes_y = layout(2)/io_layout(2)
    npes_in_group = npes_x*npes_y

    io_domain%whalo    = domain%whalo
    io_domain%ehalo    = domain%ehalo
    io_domain%shalo    = domain%shalo
    io_domain%nhalo    = domain%nhalo
    io_domain%ntiles   = 1
    io_domain%pe       = domain%pe
    io_domain%symmetry = domain%symmetry
    allocate(io_domain%list(0:npes_in_group-1))
    do i = 0, npes_in_group-1
       allocate( io_domain%list(i)%x(1), io_domain%list(i)%y(1), io_domain%list(i)%tile_id(1) )
    enddo
    ipos = mod(domain%x(1)%pos, npes_x)
    jpos = mod(domain%y(1)%pos, npes_y)
    igroup  = domain%x(1)%pos/npes_x
    jgroup  = domain%y(1)%pos/npes_y
    ipos_beg = igroup*npes_x; ipos_end = ipos_beg + npes_x - 1
    jpos_beg = jgroup*npes_y; jpos_end = jpos_beg + npes_y - 1
    n = 0
    do j = jpos_beg, jpos_end
       do i = ipos_beg, ipos_end
          io_domain%list(n)%pe = domain%pearray(i,j)
          m = io_domain%list(n)%pe - mpp_root_pe()
          io_domain%list(n)%x(1)%compute = domain%list(m)%x(1)%compute
          io_domain%list(n)%y(1)%compute = domain%list(m)%y(1)%compute
          igroup                         = domain%list(m)%x(1)%pos/npes_x
          jgroup                         = domain%list(m)%y(1)%pos/npes_y
          io_domain%list(n)%tile_id(1)   = jgroup*io_layout(1) + igroup
          n = n + 1
       enddo
    enddo

    allocate(io_domain%x(1), io_domain%y(1), io_domain%tile_id(1) )
    allocate(io_domain%x(1)%list(0:npes_x-1), io_domain%y(1)%list(0:npes_y-1) )
    n    = jpos*npes_x + ipos
    io_domain%pos          = n
    io_domain%x(1)%compute = domain%x(1)%compute
    io_domain%x(1)%data    = domain%x(1)%data
    io_domain%x(1)%memory  = domain%x(1)%memory
    io_domain%y(1)%compute = domain%y(1)%compute
    io_domain%y(1)%data    = domain%y(1)%data
    io_domain%y(1)%memory  = domain%y(1)%memory
    io_domain%x(1)%global%begin     = io_domain%list(0)%x(1)%compute%begin
    io_domain%x(1)%global%end       = io_domain%list(npes_in_group-1)%x(1)%compute%end
    io_domain%x(1)%global%size      = io_domain%x(1)%global%end - io_domain%x(1)%global%begin + 1
    io_domain%x(1)%global%max_size  = io_domain%x(1)%global%size
    io_domain%y(1)%global%begin     = io_domain%list(0)%y(1)%compute%begin
    io_domain%y(1)%global%end       = io_domain%list(npes_in_group-1)%y(1)%compute%end
    io_domain%y(1)%global%size      = io_domain%y(1)%global%end - io_domain%y(1)%global%begin + 1
    io_domain%y(1)%global%max_size  = io_domain%y(1)%global%size
    io_domain%x(1)%pos     = ipos
    io_domain%y(1)%pos     = jpos
    io_domain%tile_id(1)   = io_domain%list(n)%tile_id(1)
    io_domain%tile_root_pe = io_domain%list(0)%pe

    !z1l
!!$    do j = 0, npes_y - 1
!!$       n = j*npes_x + ipos
!!$       io_domain%y(1)%list(j) = io_domain%list(n)%y(1)
!!$    enddo
!!$    do i = 0, npes_x - 1
!!$       n = jpos*npes_x + i
!!$       io_domain%x(1)%list(i) = io_domain%list(n)%x(1)
!!$    enddo

    whalo = domain%whalo
    ehalo = domain%ehalo
    shalo = domain%shalo
    nhalo = domain%nhalo
    
    io_domain=>NULL()


  end subroutine mpp_define_io_domain

  ! <SUBROUTINE NAME="mpp_define_domains2D" INTERFACE="mpp_define_domains">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"> </IN>
  !  <IN NAME="layout" TYPE="integer" DIM="(2)"></IN>
  !  <INOUT NAME="domain" TYPE="type(domain2D)"></INOUT>
  !  <IN NAME="pelist" TYPE="integer" DIM="(0:)"></IN>
  !  <IN NAME="xflags, yflags" TYPE="integer"></IN>
  !  <IN NAME="xhalo, yhalo" TYPE="integer"></IN>
  !  <IN NAME="xextent, yextent" TYPE="integer" DIM="(0:)"></IN>
  !  <IN NAME="maskmap" TYPE="logical" DIM="(:,:)"></IN>
  !  <IN NAME="name" TYPE="character(len=*)"></IN>
  ! </SUBROUTINE>
  subroutine mpp_define_domains2D( global_indices, layout, domain, pelist, xflags, yflags,    &
         xhalo, yhalo, xextent, yextent, maskmap, name, symmetry,  memory_size,               &
         whalo, ehalo, shalo, nhalo, is_mosaic, tile_count, tile_id, complete, x_cyclic_offset, y_cyclic_offset )
    !define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
    integer,          intent(in)           :: global_indices(:) !(/ isg, ieg, jsg, jeg /)
    integer,          intent(in)           :: layout(:)
    type(domain2D),   intent(inout)        :: domain
    integer,          intent(in), optional :: pelist(0:)
    integer,          intent(in), optional :: xflags, yflags, xhalo, yhalo
    integer,          intent(in), optional :: xextent(0:), yextent(0:)
    logical,          intent(in), optional :: maskmap(0:,0:)
    character(len=*), intent(in), optional :: name
    logical,          intent(in), optional :: symmetry
    logical,          intent(in), optional :: is_mosaic                  ! indicate if calling mpp_define_domains from mpp_define_mosaic.
    integer,          intent(in), optional :: memory_size(:)
    integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! halo size for West, East, South and North direction.
                                                                         ! if whalo and ehalo is not present, 
                                                                         ! will take the value of xhalo
                                                                         ! if shalo and nhalo is not present, 
                                                                         ! will take the value of yhalo
    integer, intent(in),          optional :: tile_count                 ! tile number on current pe, default value is 1
                                                                         ! this is for the situation that multiple tiles on one processor.
    integer, intent(in),          optional :: tile_id                    ! tile id
    logical, intent(in),          optional :: complete                   ! true indicate mpp_define_domain is completed for mosaic definition.    
    integer, intent(in),          optional :: x_cyclic_offset            ! offset for x-cyclic boundary condition, 
                                                                         ! (0,j) = (ni, mod(j+x_cyclic_offset,nj))
                                                                         ! (ni+1,j) = ( 1, mod(j+nj-x_cyclic_offset,nj) )
    integer, intent(in),          optional :: y_cyclic_offset            ! offset for y-cyclic boundary condition
                                                                         ! (i,0)    = (mod(i+y_cyclic_offset,ni), nj))
                                                                         ! (i,nj+1) = (mod(mod(i+ni-y_cyclic_offset,ni), 1) )

    integer              :: i, j, m, n, xhalosz, yhalosz, memory_xsize, memory_ysize
    integer              :: whalosz, ehalosz, shalosz, nhalosz
    integer              :: ipos, jpos, pos, tile, nlist, cur_tile_id
    integer              :: ndivx, ndivy, isg, ieg, jsg, jeg, ishift, jshift, errunit, logunit
    integer              :: x_offset, y_offset, start_pos, nfold
    logical              :: from_mosaic, is_complete
    logical              :: mask(0:layout(1)-1,0:layout(2)-1)
    integer, allocatable :: pes(:), pesall(:)
    integer              :: pearray(0:layout(1)-1,0:layout(2)-1)
    integer              :: ibegin(0:layout(1)-1), iend(0:layout(1)-1)
    integer              :: jbegin(0:layout(2)-1), jend(0:layout(2)-1)
    character(len=8)     :: text

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: You must first call mpp_domains_init.' )
    if(PRESENT(name)) then
       if(len_trim(name) > NAME_LENGTH) call mpp_error(FATAL,  &
            "mpp_domains_define.inc(mpp_define_domains2D): the len_trim of optional argument name ="//trim(name)// &
            " is greater than NAME_LENGTH, change the argument name or increase NAME_LENGTH")  
       domain%name = name
    endif
    if(size(global_indices(:)) .NE. 4) call mpp_error(FATAL,   &
       "mpp_define_domains2D: size of global_indices should be 4 for "//trim(domain%name) )
    if(size(layout(:)) .NE. 2) call mpp_error(FATAL,"mpp_define_domains2D: size of layout should be 2 for "//trim(domain%name) )

    ndivx = layout(1); ndivy = layout(2)
    isg = global_indices(1); ieg = global_indices(2); jsg = global_indices(3); jeg = global_indices(4)

    from_mosaic = .false.
    if(present(is_mosaic)) from_mosaic = is_mosaic
    is_complete = .true.
    if(present(complete)) is_complete = complete
    tile = 1
    if(present(tile_count)) tile = tile_count
    cur_tile_id = 1
    if(present(tile_id)) cur_tile_id = tile_id
  
    if( PRESENT(pelist) )then
       allocate( pes(0:size(pelist(:))-1) )
       pes = pelist
       if(from_mosaic) then
          allocate( pesall(0:mpp_npes()-1) )
          call mpp_get_current_pelist(pesall)   
       else
          allocate( pesall(0:size(pes(:))-1) )
          pesall = pes
       end if
    else
       allocate( pes(0:mpp_npes()-1) )
       allocate( pesall(0:mpp_npes()-1) )
       call mpp_get_current_pelist(pes)
       pesall = pes
    end if

    !--- at least of one of x_cyclic_offset and y_cyclic_offset must be zero
    !--- folded boundary condition is not supported when either x_cyclic_offset or  y_cyclic_offset is nonzero.
    !--- Since we only implemented Folded-north boundary condition currently, we only consider y-flags. 
    x_offset = 0; y_offset = 0
    if(PRESENT(x_cyclic_offset)) x_offset = x_cyclic_offset
    if(PRESENT(y_cyclic_offset)) y_offset = y_cyclic_offset
    if(x_offset*y_offset .NE. 0) call mpp_error(FATAL, &
       'MPP_DEFINE_DOMAINS2D: At least one of x_cyclic_offset and y_cyclic_offset must be zero for '//trim(domain%name))

    !--- x_cyclic_offset and y_cyclic_offset should no larger than the global grid size.
    if(abs(x_offset) > jeg-jsg+1) call mpp_error(FATAL, &
         'MPP_DEFINE_DOMAINS2D: absolute value of x_cyclic_offset is greater than jeg-jsg+1 for '//trim(domain%name))
    if(abs(y_offset) > ieg-isg+1) call mpp_error(FATAL, &
         'MPP_DEFINE_DOMAINS2D: absolute value of y_cyclic_offset is greater than ieg-isg+1 for '//trim(domain%name))

    !--- when there is more than one tile on one processor, all the tile will limited on this processor
    if( tile > 1 .AND. size(pes(:)) > 1) call mpp_error(FATAL, &
         'MPP_DEFINE_DOMAINS2D: there are more than one tile on this pe, '// &
         'all the tile should be limited on this pe for '//trim(domain%name))

    !--- the position of current pe is changed due to mosaic, because  pes 
    !--- is only part of the pelist in mosaic (pesall). We assume the pe 
    !--- distribution are contious in mosaic.
    pos = mpp_pe() - mpp_root_pe()

    domain%symmetry = .FALSE.
    if(present(symmetry)) domain%symmetry = symmetry
    if(domain%symmetry) then
       ishift = 1; jshift = 1
    else
       ishift = 0; jshift = 0
    end if

    !--- first compute domain decomposition.
    call mpp_compute_extent(isg, ieg, ndivx, ibegin, iend, xextent)    
    call mpp_compute_extent(jsg, jeg, ndivy, jbegin, jend, yextent)    

    xhalosz = 0; yhalosz = 0
    if(present(xhalo)) xhalosz = xhalo
    if(present(yhalo)) yhalosz = yhalo
    whalosz = xhalosz; ehalosz = xhalosz
    shalosz = yhalosz; nhalosz = yhalosz
    if(present(whalo)) whalosz = whalo
    if(present(ehalo)) ehalosz = ehalo
    if(present(shalo)) shalosz = shalo
    if(present(nhalo)) nhalosz = nhalo

    !--- configure maskmap
    mask = .TRUE.
    if( PRESENT(maskmap) )then
       if( size(maskmap,1).NE.ndivx .OR. size(maskmap,2).NE.ndivy ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: maskmap array does not match layout for '//trim(domain%name) )
       mask(:,:) = maskmap(:,:)
    end if
    !number of unmask domains in layout must equal number of PEs assigned
    n = count(mask)
    if( n.NE.size(pes(:)) )then
       write( text,'(i8)' )n
       call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: incorrect number of PEs assigned for ' // &
            'this layout and maskmap. Use '//text//' PEs for this domain decomposition for '//trim(domain%name) )
    end if

    memory_xsize = 0; memory_ysize = 0
    if(present(memory_size)) then
       if(size(memory_size(:)) .NE. 2) call mpp_error(FATAL,  &
             "mpp_define_domains2D: size of memory_size should be 2 for "//trim(domain%name))
       memory_xsize = memory_size(1)
       memory_ysize = memory_size(2)
    end if

    !--- set up domain%list.
    !--- set up 2-D domain decomposition for T, E, C, N and computing overlapping
    !--- when current tile is the last tile in the mosaic. 
    nlist = size(pesall(:))
    if( .NOT. Associated(domain%x) ) then
       allocate(domain%tileList(1))
       domain%tileList(1)%xbegin = global_indices(1)
       domain%tileList(1)%xend   = global_indices(2)
       domain%tileList(1)%ybegin = global_indices(3)
       domain%tileList(1)%yend   = global_indices(4)
       allocate(domain%x(1), domain%y(1) )
       allocate(domain%tile_id(1))
       domain%tile_id        = cur_tile_id
       domain%ntiles         = 1
       domain%max_ntile_pe   = 1
       domain%ncontacts      = 0
       domain%rotated_ninety = .FALSE.
       allocate( domain%list(0:nlist-1) )        
       do i = 0, nlist-1
          allocate( domain%list(i)%x(1), domain%list(i)%y(1), domain%list(i)%tile_id(1) )
       end do
    end if

    start_pos = 0
    do n = 0, nlist-1
       if(pesall(n) == pes(0)) then
          start_pos = n
          exit
       endif
    enddo
    
    !place on PE array; need flag to assign them to j first and then i
    pearray(:,:) = NULL_PE
    ipos = NULL_PE; jpos = NULL_PE
    n = 0
    m = start_pos
    do j = 0,ndivy-1
       do i = 0,ndivx-1
          if( mask(i,j) )then
             pearray(i,j) = pes(n)
             domain%list(m)%x(tile)%compute%begin = ibegin(i)
             domain%list(m)%x(tile)%compute%end   = iend(i)
             domain%list(m)%y(tile)%compute%begin = jbegin(j)
             domain%list(m)%y(tile)%compute%end   = jend(j)
             domain%list(m)%x(tile)%compute%size  = domain%list(m)%x(tile)%compute%end - domain%list(m)%x(tile)%compute%begin + 1
             domain%list(m)%y(tile)%compute%size  = domain%list(m)%y(tile)%compute%end - domain%list(m)%y(tile)%compute%begin + 1
             domain%list(m)%tile_id(tile)         = cur_tile_id
             domain%list(m)%x(tile)%pos           = i
             domain%list(m)%y(tile)%pos           = j
             domain%list(m)%tile_root_pe          = pes(0)  
             domain%list(m)%pe                    = pesall(m)  

             if( pes(n).EQ.mpp_pe() )then
                ipos = i
                jpos = j
             end if
             n = n + 1
             m = m + 1
          end if
       end do
    end do

    !Considering mosaic, the following will only be done on the pe in the pelist
    !when there is only one tile, all the current pe will be in the pelist.
    if( ANY(pes == mpp_pe()) ) then
       domain%io_layout = layout
       domain%tile_root_pe  = pes(0)
       if( ipos.EQ.NULL_PE .OR. jpos.EQ.NULL_PE ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: pelist must include this PE for '//trim(domain%name) )
       if( debug ) then
         errunit = stderr()
         write( errunit, * )'pe, tile, ipos, jpos=', mpp_pe(), tile, ipos, jpos, ' pearray(:,jpos)=', &
                  pearray(:,jpos), ' pearray(ipos,:)=', pearray(ipos,:)
       endif

       !--- when tile is not equal to 1, the layout for that tile always ( 1, 1), so no need for pearray in domain
       if( tile == 1 ) then
          allocate( domain%pearray(0:ndivx-1,0:ndivy-1) )
          domain%pearray = pearray
       end if

       domain%pe  = mpp_pe()
       domain%pos  = pos
       domain_cnt = domain_cnt + INT(1,KIND=LONG_KIND)
       domain%id = domain_cnt*DOMAIN_ID_BASE  ! Must be LONG_KIND arithmetic

       !do domain decomposition using 1D versions in X and Y, 
       call mpp_define_domains( global_indices(1:2), ndivx, domain%x(tile), &
            pack(pearray(:,jpos),mask(:,jpos)), xflags, xhalo, xextent, mask(:,jpos), memory_xsize, whalo, ehalo )
       call mpp_define_domains( global_indices(3:4), ndivy, domain%y(tile), &
            pack(pearray(ipos,:),mask(ipos,:)), yflags, yhalo, yextent, mask(ipos,:), memory_ysize, shalo, nhalo  )
       if( domain%x(tile)%list(ipos)%pe.NE.domain%y(tile)%list(jpos)%pe ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain%x%list(ipos)%pe.NE.domain%y%list(jpos)%pe.' ) 

       !set up fold, when the boundary is folded, there is only one tile.
       domain%fold = 0
       nfold = 0
       if( PRESENT(xflags) )then
          if( BTEST(xflags,WEST) ) then
             !--- make sure no cross-domain in y-direction
             if(domain%x(tile)%data%begin .LE. domain%x(tile)%global%begin .AND. &
                domain%x(tile)%compute%begin > domain%x(tile)%global%begin ) then
                call mpp_error(FATAL, &
                 'MPP_DEFINE_DOMAINS: the domain could not be crossed when west is folded')
             endif
             if( domain%x(tile)%cyclic )call mpp_error( FATAL, &
                   'MPP_DEFINE_DOMAINS: an axis cannot be both folded west and cyclic for '//trim(domain%name) )
             domain%fold = domain%fold + FOLD_WEST_EDGE
             nfold = nfold+1
          endif
          if( BTEST(xflags,EAST) ) then
             !--- make sure no cross-domain in y-direction
             if(domain%x(tile)%data%end .GE. domain%x(tile)%global%end .AND. &
                domain%x(tile)%compute%end < domain%x(tile)%global%end ) then
                call mpp_error(FATAL, &
                 'MPP_DEFINE_DOMAINS: the domain could not be crossed when north is folded')
             endif
             if( domain%x(tile)%cyclic )call mpp_error( FATAL, &
                  'MPP_DEFINE_DOMAINS: an axis cannot be both folded east and cyclic for '//trim(domain%name) )
             domain%fold = domain%fold + FOLD_EAST_EDGE
             nfold = nfold+1
          endif
       endif
       if( PRESENT(yflags) )then
          if( BTEST(yflags,SOUTH) ) then
             !--- make sure no cross-domain in y-direction
             if(domain%y(tile)%data%begin .LE. domain%y(tile)%global%begin .AND. &
                domain%y(tile)%compute%begin > domain%y(tile)%global%begin ) then
                call mpp_error(FATAL, &
                 'MPP_DEFINE_DOMAINS: the domain could not be crossed when south is folded')
             endif
             if( domain%y(tile)%cyclic )call mpp_error( FATAL, &
                 'MPP_DEFINE_DOMAINS: an axis cannot be both folded north and cyclic for '//trim(domain%name))
             domain%fold = domain%fold + FOLD_SOUTH_EDGE
             nfold = nfold+1
          endif
          if( BTEST(yflags,NORTH) ) then
             !--- make sure no cross-domain in y-direction
             if(domain%y(tile)%data%end .GE. domain%y(tile)%global%end .AND. &
                domain%y(tile)%compute%end < domain%y(tile)%global%end ) then
                write(*,*)'WARNING from mpp_domains_define.inc(mpp_define_domains2D): '// &
                'the halo data may be not correct because the domain is crossed and ' // &
                   'the halo lies on the north folded region'
             endif
             if( domain%y(tile)%cyclic )call mpp_error( FATAL, &
                'MPP_DEFINE_DOMAINS: an axis cannot be both folded south and cyclic for '//trim(domain%name) )
             domain%fold = domain%fold + FOLD_NORTH_EDGE
             nfold = nfold+1
          endif
       endif
       if(nfold > 1) call mpp_error(FATAL, &
           'MPP_DEFINE_DOMAINS2D: number of folded edge is greater than 1 for '//trim(domain%name) )
 
       if(nfold == 1) then
           if( x_offset .NE. 0 .OR. y_offset .NE. 0) call mpp_error(FATAL, &
                  'MPP_DEFINE_DOMAINS2D: For the foled_north/folded_south/fold_east/folded_west boundary condition,  '// &
                  'x_cyclic_offset and y_cyclic_offset must be zero for '//trim(domain%name))
       endif
       if( BTEST(domain%fold,SOUTH) .OR. BTEST(domain%fold,NORTH) )then
          if( domain%y(tile)%cyclic )call mpp_error( FATAL, &
              'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic for '//trim(domain%name) )
          if( modulo(domain%x(tile)%global%size,2).NE.0 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in X must be even ' // &
                  'when there is a fold in Y for '//trim(domain%name) )
          !check if folded domain boundaries line up in X: compute domains lining up is a sufficient condition for symmetry
          n = ndivx - 1
          do i = 0,n/2
             if( domain%x(tile)%list(i)%compute%size.NE.domain%x(tile)%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries ' // &
                                         'must line up (mirror-symmetric extents) for '//trim(domain%name) )
          end do
       end if
       if( BTEST(domain%fold,WEST) .OR. BTEST(domain%fold,EAST) )then
          if( domain%x(tile)%cyclic )call mpp_error( FATAL, &
             'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic for '//trim(domain%name) )
          if( modulo(domain%y(tile)%global%size,2).NE.0 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in Y must be even '//&
                   'when there is a fold in X for '//trim(domain%name) )
          !check if folded domain boundaries line up in Y: compute domains lining up is a sufficient condition for symmetry
          n = ndivy - 1
          do i = 0,n/2
             if( domain%y(tile)%list(i)%compute%size.NE.domain%y(tile)%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must '//&
                        'line up (mirror-symmetric extents) for '//trim(domain%name) )
          end do
       end if

       !set up domain%list
       if( mpp_pe().EQ.pes(0) .AND. PRESENT(name) )then
          logunit = stdlog()
          write( logunit, '(/a,i3,a,i3)' )trim(name)//' domain decomposition: ', ndivx, ' X', ndivy
          write( logunit, '(3x,a)' )'pe,   is,  ie,  js,  je,    isd, ied, jsd, jed'
       end if
    end if   !     if( ANY(pes == mpp_pe()) ) 

    if(is_complete) then
       domain%whalo = whalosz; domain%ehalo = ehalosz
       domain%shalo = shalosz; domain%nhalo = nhalosz
       allocate(domain%update_T, domain%update_E, domain%update_C, domain%update_N)
       allocate(domain%check_E,  domain%check_C,  domain%check_N )
       domain%update_T%nsend = 0
       domain%update_T%nrecv = 0
       domain%update_C%nsend = 0
       domain%update_C%nrecv = 0
       domain%update_E%nsend = 0
       domain%update_E%nrecv = 0
       domain%update_N%nsend = 0
       domain%update_N%nrecv = 0

       if( BTEST(domain%fold,SOUTH) ) then
          call compute_overlaps_fold_south(domain, CENTER,      0,      0)
          call compute_overlaps_fold_south(domain, CORNER, ishift, jshift)
          call compute_overlaps_fold_south(domain, EAST,   ishift,      0)
          call compute_overlaps_fold_south(domain, NORTH,  0,      jshift)
       else if( BTEST(domain%fold,WEST) ) then
          call compute_overlaps_fold_west(domain, CENTER,      0,      0)
          call compute_overlaps_fold_west(domain, CORNER, ishift, jshift)
          call compute_overlaps_fold_west(domain, EAST,   ishift,      0)
          call compute_overlaps_fold_west(domain, NORTH,       0, jshift)
       else if( BTEST(domain%fold,EAST) ) then
          call compute_overlaps_fold_east(domain, CENTER,      0,      0)
          call compute_overlaps_fold_east(domain, CORNER, ishift, jshift)
          call compute_overlaps_fold_east(domain, EAST,   ishift,      0)
          call compute_overlaps_fold_east(domain, NORTH,       0, jshift)
       else
          call compute_overlaps(domain, CENTER,      0,      0, x_offset, y_offset)
          call compute_overlaps(domain, CORNER, ishift, jshift, x_offset, y_offset)
          call compute_overlaps(domain, EAST,   ishift,      0, x_offset, y_offset)
          call compute_overlaps(domain, NORTH,       0, jshift, x_offset, y_offset)
       endif
       call check_overlap_pe_order(domain, domain%update_T, trim(domain%name)//" update_T in mpp_define_domains")
       call check_overlap_pe_order(domain, domain%update_C, trim(domain%name)//" update_C in mpp_define_domains")
       call check_overlap_pe_order(domain, domain%update_E, trim(domain%name)//" update_E in mpp_define_domains")
       call check_overlap_pe_order(domain, domain%update_N, trim(domain%name)//" update_N in mpp_define_domains")


       !--- when ncontacts is nonzero, set_check_overlap will be called in mpp_define
       if(domain%symmetry .AND. domain%ncontacts == 0 ) then
          call set_check_overlap( domain, CORNER )
          call set_check_overlap( domain, EAST   )
          call set_check_overlap( domain, NORTH  )
          allocate(domain%bound_E,  domain%bound_C,  domain%bound_N )
          call set_bound_overlap( domain, CORNER )
          call set_bound_overlap( domain, EAST   )
          call set_bound_overlap( domain, NORTH  )
       end if
    end if

    !print out decomposition, this didn't consider maskmap.
    if( mpp_pe() .EQ. pes(0) .AND. PRESENT(name) )then
       write(*,*) trim(name)//' domain decomposition'
       write(*,'(a,i4,a,i4,a,i4,a,i4)')'whalo = ', whalosz, ", ehalo = ", ehalosz, ", shalo = ", shalosz, ", nhalo = ", nhalosz
       write (*,110) (domain%x(1)%list(i)%compute%size, i= 0, layout(1)-1)
       write (*,120) (domain%y(1)%list(i)%compute%size, i= 0, layout(2)-1)
110    format ('  X-AXIS = ',24i4,/,(11x,24i4))
120    format ('  Y-AXIS = ',24i4,/,(11x,24i4))
    endif

    deallocate( pes, pesall)

    return
  end subroutine mpp_define_domains2D

  !#####################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_define_mosaic: define mosaic domain                        !
!  NOTE: xflags and yflags is not in mpp_define_mosaic, because such relation !
!        are already defined in the mosaic relation.                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!??? do we need optional argument xextent and yextent
!??? how to specify pelist, we may use two dimensional variable pelist to represent.
!z1l: We assume the tilelist are in always limited to 1, 2, ... num_tile. If we want 
!     to remove this limitation, we need to add one more argument tilelist.
  subroutine mpp_define_mosaic( global_indices, layout, domain, num_tile, num_contact, tile1, tile2, &
                                istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2, pe_start,      &
                                pe_end, pelist, whalo, ehalo, shalo, nhalo, xextent, yextent,                  &
                                maskmap, name, memory_size, symmetry, xflags, yflags )
    integer,          intent(in)           :: global_indices(:,:)  ! The size of first indice is 4, (/ isg, ieg, jsg, jeg /)
                                                                   ! The size of second indice is number of tiles in mosaic.
    integer,          intent(in)           :: layout(:,:)
    type(domain2D),   intent(inout)        :: domain
    integer,          intent(in)           :: num_tile             ! number of tiles in the mosaic
    integer,          intent(in)           :: num_contact          ! number of contact region between tiles.
    integer,          intent(in)           :: tile1(:), tile2(:)   ! tile number
    integer,          intent(in)           :: istart1(:), iend1(:) ! i-index in tile_1 of contact region
    integer,          intent(in)           :: jstart1(:), jend1(:) ! j-index in tile_1 of contact region
    integer,          intent(in)           :: istart2(:), iend2(:) ! i-index in tile_2 of contact region
    integer,          intent(in)           :: jstart2(:), jend2(:) ! j-index in tile_2 of contact region
    integer,          intent(in)           :: pe_start(:)          ! start pe of the pelist used in each tile
    integer,          intent(in)           :: pe_end(:)            ! end pe of the pelist used in each tile
    integer,          intent(in), optional :: pelist(:)            ! list of processors used in mosaic
    integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
    integer,          intent(in), optional :: xextent(:,:), yextent(:,:) 
    logical,          intent(in), optional :: maskmap(:,:,:)         
    character(len=*), intent(in), optional :: name
    integer,          intent(in), optional :: memory_size(2)
    logical,          intent(in), optional :: symmetry
    integer,          intent(in), optional :: xflags, yflags

    integer              :: n, m, ndivx, ndivy, nc, nlist, nt, pos, n1, n2
    integer              :: whalosz, ehalosz, shalosz, nhalosz, xhalosz, yhalosz, t1, t2, tile
    integer              :: flags_x, flags_y
    logical, allocatable :: mask(:,:)
    integer, allocatable :: pes(:), xext(:), yext(:), pelist_tile(:), ntile_per_pe(:), tile_count(:)
    logical              :: is_symmetry
    integer, allocatable :: align1(:), align2(:), is1(:), ie1(:), js1(:), je1(:), is2(:), ie2(:), js2(:), je2(:)
    integer, allocatable :: isgList(:), iegList(:), jsgList(:), jegList(:)
    real,    allocatable :: refine1(:), refine2(:)

    mosaic_defined = .true.
    !--- the size of first indice of global_indices must be 4.
    if(size(global_indices, 1) .NE. 4) call mpp_error(FATAL, &
         'mpp_domains_define.inc: The size of first dimension of global_indices is not 4')
    !--- the size of second indice of global_indices must be num_tile
    if(size(global_indices, 2) .NE. num_tile) call mpp_error(FATAL, &
         'mpp_domains_define.inc: The size of second dimension of global_indices is not equal num_tile')    
    !--- the size of first indice of layout must be 2. The second dimension size of layout must equal num_tile.
    if(size(layout, 1) .NE. 2) call mpp_error(FATAL, &
         'mpp_domains_define.inc: The size of first dimension of layout is not 2')
    if(size(layout,2) .NE. num_tile)  call mpp_error(FATAL, &
         'mpp_domains_define.inc: The size of second dimension of layout is not equal num_tile')

    !--- setup pelist for the mosaic ---------------------
    nlist = mpp_npes()
    allocate(pes(0:nlist-1))
    if(present(pelist)) then
       if( nlist .NE. size(pelist(:))) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size of pelist is not equal mpp_npes')
       pes = pelist
    else
       call mpp_get_current_pelist(pes)
    end if
    !--- pelist should be monotonic increasing by 1.
    do n = 1, nlist-1
       if(pes(n) - pes(n-1) .NE. 1) call mpp_error(FATAL, &
           'mpp_domains_define.inc: pelist is not monotonic increasing by 1')
    end do

    is_symmetry        = .FALSE.
    if(present(symmetry)) is_symmetry = symmetry

    if(size(pe_start(:)) .NE. num_tile .OR. size(pe_end(:)) .NE. num_tile ) call mpp_error(FATAL, &
         'mpp_domains_define.inc: size of pe_start and/or pe_end is not equal num_tile')
    !--- make sure pe_start and pe_end is in the pelist.   
    if( ANY( pe_start < pes(0) ) ) call mpp_error(FATAL, 'mpp_domains_define.inc: not all the pe_start are in the pelist')
    if( ANY( pe_end > pes(nlist-1)) ) call mpp_error(FATAL, 'mpp_domains_define.inc: not all the pe_end are in the pelist')

    !--- calculate number of tiles on each pe.
    allocate( ntile_per_pe(0:nlist-1) )
    ntile_per_pe = 0
    do n = 1, num_tile
       do m = pe_start(n) - mpp_root_pe(), pe_end(n) - mpp_root_pe()
          ntile_per_pe(m) = ntile_per_pe(m) + 1
       end do
    end do
    if(ANY(ntile_per_pe == 0)) call mpp_error(FATAL, &
         'mpp_domains_define.inc: At least one pe in pelist is not used by any tile in the mosaic')

    !--- check the size comformable of xextent and yextent
    if( PRESENT(xextent) ) then
       if(size(xextent,1) .GT. maxval(layout(1,:)) ) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size mismatch between xextent and layout') 
       if(size(xextent,2) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size of xextent is not eqaul num_tile')    
    end if
    if( PRESENT(yextent) ) then
       if(size(yextent,1) .GT. maxval(layout(2,:)) ) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size mismatch between yextent and layout') 
       if(size(yextent,2) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size of yextent is not eqaul num_tile')    
    end if

    !--- check the size comformable of maskmap
    !--- since the layout is different between tiles, so the actual size of maskmap for each tile is 
    !--- not diffrent. When define maskmap for multiple tiles, user can choose the maximum value
    !--- of layout of all tiles to the first and second dimension of maskmap. 
    if(present(maskmap)) then
       if(size(maskmap,1) .GT. maxval(layout(1,:)) .or. size(maskmap,2) .GT. maxval(layout(2,:))) &
            call mpp_error(FATAL, 'mpp_domains_define.inc: size mismatch between maskmap and layout')  
       if(size(maskmap,3) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define.inc: the third dimension of maskmap is not equal num_tile')  
    end if

    allocate(domain%tileList(num_tile))
    do n = 1, num_tile
       domain%tileList(n)%xbegin = global_indices(1,n)
       domain%tileList(n)%xend   = global_indices(2,n)
       domain%tileList(n)%ybegin = global_indices(3,n)
       domain%tileList(n)%yend   = global_indices(4,n)
    enddo
    !--- define some mosaic information in domain type
    nt = ntile_per_pe(mpp_pe()-mpp_root_pe())
    allocate(domain%tile_id(nt), domain%x(nt), domain%y(nt) )
    allocate(domain%list(0:nlist-1))

    do n = 0, nlist-1
       nt = ntile_per_pe(mpp_root_pe()+n)
       allocate(domain%list(n)%x(nt), domain%list(n)%y(nt), domain%list(n)%tile_id(nt) )
    end do

    pe = mpp_pe()
    pos = 0
    do n = 1, num_tile
       if( pe .GE. pe_start(n) .AND. pe .LE. pe_end(n)) then
          pos = pos + 1
          domain%tile_id(pos) = n
       end if
    end do

    domain%rotated_ninety = .FALSE.
    domain%ntiles         = num_tile      
    domain%max_ntile_pe   = maxval(ntile_per_pe)
    domain%ncontacts      = num_contact
    
    deallocate(ntile_per_pe)
    !---call mpp_define_domain to define domain decomposition for each tile.
    allocate(tile_count(pes(0):pes(0)+nlist-1))
    tile_count = 0  ! tile number on current pe    

    do n = 1, num_tile
       allocate(mask(layout(1,n), layout(2,n)))
       allocate(pelist_tile(pe_start(n):pe_end(n)) )
       tile_count(pe_start(n)) = tile_count(pe_start(n)) + 1
       do m = pe_start(n), pe_end(n)
          pelist_tile(m) = m
       end do
       mask = .TRUE.
       if(present(maskmap))  mask = maskmap(1:layout(1,n), 1:layout(2,n), n)
       ndivx = layout(1,n); ndivy = layout(2,n)
       allocate(xext(ndivx), yext(ndivy))
       xext = 0; yext = 0
       if(present(xextent)) xext = xextent(1:ndivx,n)
       if(present(yextent)) yext = yextent(1:ndivy,n)
       ! when num_tile is one, we assume only folded_north and cyclic_x, cyclic_y boundary condition is the possible
       ! z1l: when we decide to support multiple-tile tripolar grid, we will redesign the following part.
       if(num_tile == 1) then
          flags_x = 0
          flags_y = 0
          if(PRESENT(xflags)) flags_x = xflags
          if(PRESENT(yflags)) flags_y = yflags
          do m = 1, num_contact
             if(istart1(m) == iend1(m) ) then  ! x-direction contact, possible cyclic, folded-west or folded-east
                if(istart2(m) .NE. iend2(m) ) call mpp_error(FATAL,  &
                   "mpp_domains_define: for one tile mosaic, when istart1=iend1, istart2 must equal iend2")
                if(istart1(m) == istart2(m) ) then ! folded west or folded east
                   if(istart1(m) == global_indices(1,n) ) then
                      if(.NOT. BTEST(flags_x,WEST) )  flags_x = flags_x + FOLD_WEST_EDGE
                   else if(istart1(m) == global_indices(2,n) ) then
                      if(.NOT. BTEST(flags_x,EAST) )  flags_x = flags_x + FOLD_EAST_EDGE
                   else
                       call mpp_error(FATAL, "mpp_domains_define: when istart1=iend1,jstart1=jend1, "//&
                         "istart1 should equal global_indices(1) or global_indices(2)")
                   endif
                else
                if(.NOT. BTEST(flags_x,CYCLIC))  flags_x = flags_x + CYCLIC_GLOBAL_DOMAIN
                endif
             else if( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic, folded-south or folded-north
                if(jstart2(m) .NE. jend2(m) ) call mpp_error(FATAL,  &
                   "mpp_domains_define: for one tile mosaic, when jstart1=jend1, jstart2 must equal jend2")
                if(jstart1(m) == jstart2(m) ) then ! folded south or folded north
                   if(jstart1(m) == global_indices(3,n) ) then
                      if(.NOT. BTEST(flags_y,SOUTH) )  flags_y = flags_y + FOLD_SOUTH_EDGE
                   else if(jstart1(m) == global_indices(4,n) ) then
                   if(.NOT. BTEST(flags_y,NORTH) )  flags_y = flags_y + FOLD_NORTH_EDGE
                else 
                       call mpp_error(FATAL, "mpp_domains_define: when istart1=iend1,jstart1=jend1, "//&
                         "istart1 should equal global_indices(1) or global_indices(2)")
                   endif
                else
                   if(.NOT. BTEST(flags_y,CYCLIC))  flags_y = flags_y + CYCLIC_GLOBAL_DOMAIN
                end if
             else 
               call mpp_error(FATAL,  &
                   "mpp_domains_define: for one tile mosaic, invalid boundary contact")
             end if
          end do
          call mpp_define_domains(global_indices(:,n), layout(:,n), domain, pelist=pelist_tile, xflags = flags_x, &
                                  yflags = flags_y, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo,           &
                                  xextent=xext, yextent=yext, maskmap=mask, name=name, symmetry=is_symmetry,      &
                                  memory_size = memory_size, is_mosaic = .true.)
       else
          call mpp_define_domains(global_indices(:,n), layout(:,n), domain, pelist=pelist_tile,                   &
                                  whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, xextent=xext, yextent=yext, &
                                  maskmap=mask, name=name, symmetry=is_symmetry, memory_size = memory_size,       &
                                  is_mosaic = .true., tile_count = tile_count(pe_start(n)), tile_id=n, complete = n==num_tile)
       end if
       deallocate(mask, xext, yext, pelist_tile)
    end do

    deallocate(pes, tile_count)

    if(num_contact == 0 .OR. num_tile == 1) return

    !--- loop through each contact region and find the contact for each tile ( including alignment )
    !--- we assume the tiles list is continuous and starting from 1.
    allocate(is1(num_contact), ie1(num_contact), js1(num_contact), je1(num_contact) )
    allocate(is2(num_contact), ie2(num_contact), js2(num_contact), je2(num_contact) )
    allocate(isgList(num_tile), iegList(num_tile), jsgList(num_tile), jegList(num_tile) )
    allocate(align1(num_contact), align2(num_contact), refine1(num_contact), refine2(num_contact))
    !--- get the global domain for each tile
    do n = 1, num_tile
       isgList(n) = domain%tileList(n)%xbegin; iegList(n) = domain%tileList(n)%xend
       jsgList(n) = domain%tileList(n)%ybegin; jegList(n) = domain%tileList(n)%yend
    end do

    !--- transfer the contact index to domain index.
    nc = 0
    do n = 1, num_contact
       t1 = tile1(n)
       t2 = tile2(n)
       is1(n) = istart1(n) + isgList(t1) - 1; ie1(n) = iend1(n) + isgList(t1) - 1
       js1(n) = jstart1(n) + jsgList(t1) - 1; je1(n) = jend1(n) + jsgList(t1) - 1
       is2(n) = istart2(n) + isgList(t2) - 1; ie2(n) = iend2(n) + isgList(t2) - 1
       js2(n) = jstart2(n) + jsgList(t2) - 1; je2(n) = jend2(n) + jsgList(t2) - 1
       call check_alignment( is1(n), ie1(n), js1(n), je1(n), isgList(t1), iegList(t1), jsgList(t1), jegList(t1), align1(n))
       call check_alignment( is2(n), ie2(n), js2(n), je2(n), isgList(t2), iegList(t2), jsgList(t2), jegList(t2), align2(n))
       if( (align1(n) == WEST .or. align1(n) == EAST ) .NEQV. (align2(n) == WEST .or. align2(n) == EAST ) )&
             domain%rotated_ninety=.true.   
    end do

    !--- calculate the refinement ratio between tiles
    do n = 1, num_contact
       n1 = max(abs(iend1(n) - istart1(n)), abs(jend1(n) - jstart1(n)) ) + 1
       n2 = max(abs(iend2(n) - istart2(n)), abs(jend2(n) - jstart2(n)) ) + 1
       refine1(n) = real(n2)/n1
       refine2(n) = real(n1)/n2
    end do    

    whalosz = 0; ehalosz = 0; shalosz = 0; nhalosz = 0
    if(present(whalo)) whalosz = whalo
    if(present(ehalo)) ehalosz = ehalo
    if(present(shalo)) shalosz = shalo
    if(present(nhalo)) nhalosz = nhalo
    xhalosz = max(whalosz, ehalosz)
    yhalosz = max(shalosz, nhalosz)

    !--- computing the overlap for the contact region with halo size xhalosz and yhalosz
    call define_contact_point( domain, CENTER, num_contact, tile1, tile2, align1, align2, refine1, refine2, &
                               is1, ie1, js1, je1, is2, ie2, js2, je2, isgList, iegList, jsgList, jegList )

    call set_contact_point( domain, CORNER )
    call set_contact_point( domain, EAST   )
    call set_contact_point( domain, NORTH  )

    !--- goffset setting is needed for exact global sum
    do m = 1, size(domain%tile_id(:))
       tile = domain%tile_id(m)
       do n = 1, num_contact
          if( tile1(n) == tile ) then
             if(align1(n) == EAST ) domain%x(m)%goffset = 0
             if(align1(n) == NORTH) domain%y(m)%goffset = 0
          end if
          if( tile2(n) == tile ) then
             if(align2(n) == EAST ) domain%x(m)%goffset = 0
             if(align2(n) == NORTH) domain%y(m)%goffset = 0
          end if
       end do
    end do
    call check_overlap_pe_order(domain, domain%update_T, trim(domain%name)//" update_T in mpp_define_mosaic")
    call check_overlap_pe_order(domain, domain%update_C, trim(domain%name)//" update_C in mpp_define_mosaic")
    call check_overlap_pe_order(domain, domain%update_E, trim(domain%name)//" update_E in mpp_define_mosaic")
    call check_overlap_pe_order(domain, domain%update_N, trim(domain%name)//" update_N in mpp_define_mosaic")

    !--- set the overlapping for boundary check if domain is symmetry
    if(debug_update_level .NE. NO_CHECK) then
       call set_check_overlap( domain, CORNER )
       call set_check_overlap( domain, EAST   )
       call set_check_overlap( domain, NORTH  )
    endif
    if(domain%symmetry) then
       allocate(domain%bound_E,  domain%bound_C,  domain%bound_N )
       call set_bound_overlap( domain, CORNER )
       call set_bound_overlap( domain, EAST   )
       call set_bound_overlap( domain, NORTH  )
       call check_overlap_pe_order(domain, domain%bound_C, trim(domain%name)//" bound_C")
       call check_overlap_pe_order(domain, domain%bound_E, trim(domain%name)//" bound_E")
       call check_overlap_pe_order(domain, domain%bound_N, trim(domain%name)//" bound_N")
    end if

    !--- release memory
    deallocate(align1, align2, is1, ie1, js1, je1, is2, ie2, js2, je2  )
    deallocate(isgList, iegList, jsgList, jegList, refine1, refine2 )
  end subroutine mpp_define_mosaic

!#####################################################################
  logical function mpp_mosaic_defined()
  ! Accessor function for value of mosaic_defined
    mpp_mosaic_defined = mosaic_defined
  end function mpp_mosaic_defined
!#####################################################################

  subroutine compute_overlaps( domain, position, ishift, jshift, x_cyclic_offset, y_cyclic_offset )
    !computes remote domain overlaps
    !assumes only one in each direction
    !will calculate the overlapping for T,E,C,N-cell seperately.
    type(domain2D), intent(inout) :: domain
    integer, intent(in)           :: position, ishift, jshift
    integer, intent(in)           :: x_cyclic_offset, y_cyclic_offset

    integer                          :: i, m, n, nlist, tMe, tNbr, dir
    integer                          :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd, jed
    integer                          :: isg, ieg, jsg, jeg, ioff, joff
    integer                          :: list, middle, ni, nj, isgd, iegd, jsgd, jegd
    integer                          :: ism, iem, jsm, jem, whalo, ehalo, shalo, nhalo
    logical                          :: folded, need_adjust_1, need_adjust_2, need_adjust_3, folded_north
    type(overlap_type)               :: overlap
    type(overlapSpec),   pointer     :: update=>NULL()
    type(overlap_type),  pointer     :: overlapList(:)=>NULL()
    type(overlap_type),  pointer     :: checkList(:)=>NULL()
    type(overlapSpec),   pointer     :: check =>NULL()
    integer                          :: nsend, nrecv
    integer                          :: nsend_check, nrecv_check

    !--- since we restrict that if multiple tiles on one pe, all the tiles are limited to this pe. 
    !--- In this case, if ntiles on this pe is greater than 1, no overlapping between processor within each tile
    !--- In this case the overlapping exist only for tMe=1 and tNbr=1
    if(size(domain%x(:)) > 1) return

    !--- if there is no halo, no need to compute overlaps.
    if(domain%whalo==0 .AND. domain%ehalo==0 .AND. domain%shalo==0 .AND. domain%nhalo==0) return

    !--- when there is only one tile, n will equal to np
    nlist = size(domain%list(:))

    select case(position)
    case (CENTER)
       update => domain%update_T
       check  => NULL()
    case (CORNER)
       update => domain%update_C
       check  => domain%check_C
    case (EAST)
       update => domain%update_E
       check  => domain%check_E
    case (NORTH)
       update => domain%update_N
       check  => domain%check_N
    case default
       call mpp_error(FATAL, &
        "mpp_domains_define.inc(compute_overlaps): the value of position should be CENTER, EAST, CORNER or NORTH")
    end select

    allocate(overlapList(MAXLIST) )
    allocate(checkList(MAXLIST)   )

    !--- overlap is used to store the overlapping temporarily.
    call allocate_update_overlap( overlap, MAXOVERLAP)
    !send
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec, position=position )
    call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ni, ysize=nj, position=position ) !cyclic offsets
    call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )

    update%xbegin = ism; update%xend = iem
    update%ybegin = jsm; update%yend = jem
    if(ASSOCIATED(check)) then
       check%xbegin  = ism; check%xend  = iem
       check%ybegin  = jsm; check%yend  = jem
    endif
    update%whalo  = domain%whalo; update%ehalo = domain%ehalo
    update%shalo  = domain%shalo; update%nhalo = domain%nhalo
    whalo         = domain%whalo; ehalo        = domain%ehalo
    shalo         = domain%shalo; nhalo        = domain%nhalo

    ioff = ni - ishift
    joff = nj - jshift
    middle = (isg+ieg)/2+1
    tMe = 1; tNbr = 1
    folded_north = BTEST(domain%fold,NORTH)
    if( BTEST(domain%fold,SOUTH) .OR. BTEST(domain%fold,EAST) .OR. BTEST(domain%fold,WEST) ) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps): folded south, east or west boundary condition " // &
            "is not supported, please use other version of compute_overlaps for "//trim(domain%name))
    endif

    nsend = 0
    nsend_check = 0

    do list = 0,nlist-1
       m = mod( domain%pos+list, nlist )
       if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
          !to_pe's eastern halo
          dir = 1
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
          !--- to make sure the consistence between pes
          if( domain%symmetry .AND. (position == NORTH .OR. position == CORNER ) &
               .AND. ( jsc == je .or. jec == js ) ) then
             !--- do nothing, this point will come from other pe
          else 
             if( ie.GT.ieg )then
                if( domain%x(tMe)%cyclic .AND. iec.LT.is )then !try cyclic offset
                   is = is-ioff; ie = ie-ioff
                   call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jeg, nj)
                end if
             end if
             !--- when the south face is folded, the east halo point at right side domain will be folded.
             !--- the position should be on CORNER or NORTH
             if( je == jeg .AND. folded_north .AND. (position == CORNER .OR. position == NORTH) & 
                  .AND. is .GE. middle .AND. domain%list(m)%x(tNbr)%compute%end+ehalo+ishift .LE. ieg ) then
                call insert_update_overlap( overlap, domain%list(m)%pe,              &
                                            is, ie, js, je-1, isc, iec, jsc, jec, dir)
                is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
                js = je
                select case (position)
                case(NORTH)
                   i=is; is = isg+ieg-ie; ie = isg+ieg-i  
                case(CORNER)
                   i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift      
                end select
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, .true.)
             else 
                if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
                   ! (js, jeg) and (jsg, je).
                   call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                               is, ie, jsg, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
                   call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                               is, ie, js, jeg, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
                else
                   call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                               is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
                end if
             end if
          end if

          !to_pe's SE halo
          dir = 2
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%begin-shalo;    je = domain%list(m)%y(tNbr)%compute%begin-1
          need_adjust_1 = .true.; need_adjust_2 = .true.; need_adjust_3 = .true.
          if( ie.GT.ieg )then
             if( domain%x(tMe)%cyclic .AND. iec.LT.is )then !try cyclic offset
                is = is-ioff; ie = ie-ioff
                need_adjust_1 = .false.
                if(jsg .GT. js) then
                   if( domain%y(tMe)%cyclic .AND. je.LT.jsc )then !try cyclic offset
                      js = js+joff; je = je+joff
                      need_adjust_2 = .false.
                      if(x_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jeg, nj)
                      else if(y_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, ieg, ni)
                      end if
                   end if
                else
                   call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jeg, nj)
                   need_adjust_3 = .false.
                end if
             end if
          end if
          if( need_adjust_3 .AND. jsg.GT.js )then
             if( need_adjust_2 .AND. domain%y(tMe)%cyclic .AND. je.LT.jsc )then !try cyclic offset
                js = js+joff; je = je+joff
                if(need_adjust_1 .AND. ie.LE.ieg) then
                   call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, ieg, ni)  
                end if
             end if
          end if
          if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
             ! (js, jeg) and (jsg, je).   
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, jeg, isc, iec, jsc, jec, dir)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, jsg, je, isc, iec, jsc, jec, dir)        
          else if(is > ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, ieg) and (isg, ie).            
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ieg, js, je, isc, iec, jsc, jec, dir)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         isg, ie, js, je, isc, iec, jsc, jec, dir)
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir)
          end if
          !to_pe's southern halo
          dir = 3
          is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
          js = domain%list(m)%y(tNbr)%compute%begin-shalo; je = domain%list(m)%y(tNbr)%compute%begin-1
          if( jsg.GT.js )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsc )then !try cyclic offset
                js = js+joff; je = je+joff
                call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, ieg, ni)
             end if
          end if
          if(is>ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, ieg) and (isg, ie). 
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         isg, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ieg, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
          end if

          !to_pe's SW halo
          dir = 4
          is = domain%list(m)%x(tNbr)%compute%begin-whalo; ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%begin-shalo; je = domain%list(m)%y(tNbr)%compute%begin-1
          need_adjust_1 = .true.; need_adjust_2 = .true.; need_adjust_3 = .true.
          if( isg.GT.is )then
             if( domain%x(tMe)%cyclic .AND. ie.LT.isc )then !try cyclic offset
                is = is+ioff; ie = ie+ioff
                need_adjust_1 = .false.
                if(jsg .GT. js) then
                   if( domain%y(tMe)%cyclic .AND. je.LT.jsc )then !try cyclic offset
                      js = js+joff; je = je+joff
                      need_adjust_2 = .false.
                      if(x_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jeg, nj)
                      else if(y_cyclic_offset .NE. 0) then
                         call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, ieg, ni)
                      end if
                   end if
                else                  
                   call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jeg, nj)
                   need_adjust_3 = .false.
                end if
             end if
          end if
          if( need_adjust_3 .AND. jsg.GT.js )then
             if( need_adjust_2 .AND. domain%y(tMe)%cyclic .AND. je.LT.jsc )then !try cyclic offset
                js = js+joff; je = je+joff
                if(need_adjust_1 .AND. isg.LE.is  )then
                   call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, ieg, ni)  
                end if
             end if
          end if
          if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
             ! (js, jeg) and (jsg, je).   
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, jeg, isc, iec, jsc, jec, dir)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, jsg, je, isc, iec, jsc, jec, dir)        
          else if(is > ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, ieg) and (isg, ie).            
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ieg, js, je, isc, iec, jsc, jec, dir)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         isg, ie, js, je, isc, iec, jsc, jec, dir)
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir)
          end if

          !to_pe's western halo
          dir = 5
          is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
          if( isg.GT.is )then
             if( domain%x(tMe)%cyclic .AND. ie.LT.isc )then !try cyclic offset
                is = is+ioff; ie = ie+ioff
                call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jeg, nj)
             end if
          end if
          !--- when the north face is folded, some point at j=nj will be folded.
          !--- the position should be on CORNER or NORTH
          if( je == jeg .AND. folded_north .AND. (position == CORNER .OR. position == NORTH) &
               .AND. ( (domain%x(tMe)%cyclic .AND. domain%list(m)%x(tNbr)%compute%begin == isg)            &
               .OR. ( domain%list(m)%x(tNbr)%compute%begin-1 .GE. middle ) ) ) then
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je-1, isc, iec, jsc, jec, dir)
             !--- consider at j = jeg for west edge.
             !--- when the data is at corner and not symmetry, i = isg -1 will get from cyclic condition
             if(position == CORNER .AND. .NOT. domain%symmetry .AND. domain%list(m)%x(tNbr)%compute%begin == isg) then
                call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                           ie, ie, je, je, isc, iec, jsc, jec, dir, .true.)
             end if
             is = domain%list(m)%x(tNbr)%compute%begin-whalo;  ie = domain%list(m)%x(tNbr)%compute%begin-1
             je = domain%list(m)%y(tNbr)%compute%end+jshift;   js = je  
             if ( domain%list(m)%x(tNbr)%compute%begin == isg ) then
                select case (position)
                case(NORTH)
                   i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                case(CORNER)
                   i=is; is = 2*isg-ie-2+2*ishift; ie = 2*isg-i-2+2*ishift  
                end select
                if(ie .GT. domain%x(tMe)%compute%end+ishift) call mpp_error( FATAL, &
                     'mpp_domains_define.inc(compute_overlaps): west edge ubound error send.' )     
             else 
                select case (position)
                case(NORTH)
                   i=is; is = isg+ieg-ie; ie = isg+ieg-i
                case(CORNER)
                   i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                end select
             end if
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir, .true.)
          else
             if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
                ! (js, jeg) and (jsg, je).
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, jsg, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, jeg, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             end if
          end if

          !to_pe's NW halo
          dir = 6
          is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          need_adjust_1 = .true.; need_adjust_2 = .true.; need_adjust_3 = .true.
          if( isg.GT.is )then
             if( domain%x(tMe)%cyclic .AND. ie.LT.isc )then !try cyclic offset
                is = is+ioff; ie = ie+ioff
                need_adjust_1 = .false.
                if(je .GT. jeg) then
                   if( domain%y(tMe)%cyclic .AND. jec.LT.js )then !try cyclic offset
                      js = js-joff; je = je-joff
                      need_adjust_2 = .false.
                      if(x_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jeg, nj)
                      else if(y_cyclic_offset .NE. 0) then
                         call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, ieg, ni)
                      end if
                   end if
                else                  
                   call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jeg, nj)
                   need_adjust_3 = .false.
                end if
             end if
          end if
          folded = .FALSE.
          if( need_adjust_3 .AND. je.GT.jeg )then
             if( need_adjust_2 .AND. domain%y(tMe)%cyclic .AND. jec.LT.js )then !try cyclic offset
                js = js-joff; je = je-joff
                if( need_adjust_1  .AND. isg.LE.is)then
                   call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, ieg, ni)   
                end if
             else if( folded_north )then
                folded = .TRUE.
                call get_fold_index_north(isg, ieg, jeg, ishift, position, is, ie, js, je)
             end if
          end if
          if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
             ! (js, jeg) and (jsg, je).   
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, jeg, isc, iec, jsc, jec, dir, folded)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, jsg, je, isc, iec, jsc, jec, dir, folded)        
          else if(is > ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, ieg) and (isg, ie).            
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ieg, js, je, isc, iec, jsc, jec, dir, folded)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         isg, ie, js, je, isc, iec, jsc, jec, dir, folded)
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir, folded)
          end if
          !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
          if(is .LT. isg .AND. domain%x(tMe)%cyclic) then
             is = is + ioff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, is, js, je, isc, iec, jsc, jec, dir, folded)
          endif

          !to_pe's northern halo
          dir = 7
          folded = .FALSE. 
          is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          if( je.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. jec.LT.js )then !try cyclic offset
                js = js-joff; je = je-joff
                call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, ieg, ni)
             else if( folded_north )then
                folded = .TRUE.
                call get_fold_index_north(isg, ieg, jeg, ishift, position, is, ie, js, je)
             end if
          end if
          !--- when domain symmetry and position is EAST or CORNER, the point when isc == ie, 
          !--- no need to send, because the data on that point will come from other pe. 
          !--- come from two pe ( there will be only one point on one pe. ).
          if( domain%symmetry .AND. (position == EAST .OR. position == CORNER ) &
               .AND. ( isc == ie .or. iec == is ) ) then
             !--- do nothing, this point will come from other pe
          else
             if(is>ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
                ! (is, ieg) and (isg, ie). 
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            isg, ie, js, je, isc, iec, jsc, jec, dir, folded, symmetry=domain%symmetry)
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ieg, js, je, isc, iec, jsc, jec, dir, folded, symmetry=domain%symmetry)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, folded, symmetry=domain%symmetry)
             end if
          end if
          !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
          if(is .LT. isg .AND. domain%x(tMe)%cyclic) then
             is = is + ioff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, is, js, je, isc, iec, jsc, jec, dir, folded)
          endif

          !to_pe's NE halo
          dir = 8
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          need_adjust_1 = .true.; need_adjust_2 = .true.; need_adjust_3 = .true.
          if( ie.GT.ieg )then
             if( domain%x(tMe)%cyclic .AND. iec.LT.is )then !try cyclic offset
                is = is-ioff; ie = ie-ioff
                need_adjust_1 = .false.
                if(je .GT. jeg) then
                   if( domain%y(tMe)%cyclic .AND. jec.LT.js )then !try cyclic offset
                      js = js-joff; je = je-joff
                      need_adjust_2 = .false.
                      if(x_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jeg, nj)
                      else if(y_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, ieg, ni)
                      end if
                   end if
                else                  
                   call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jeg, nj)
                   need_adjust_3 = .false.
                end if
             end if
          end if
          folded = .false.
          if( need_adjust_3 .AND. je.GT.jeg )then
             if( need_adjust_2 .AND. domain%y(tMe)%cyclic .AND. jec.LT.js )then !try cyclic offset
                js = js-joff; je = je-joff
                if( need_adjust_1 .AND. ie.LE.ieg)then
                   call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, ieg, ni)  
                end if
             else if( folded_north )then
                folded = .TRUE.
                call get_fold_index_north(isg, ieg, jeg, ishift, position, is, ie, js, je)
             end if
          end if
          if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
             ! (js, jeg) and (jsg, je).   
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, jeg, isc, iec, jsc, jec, dir, folded)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, jsg, je, isc, iec, jsc, jec, dir, folded)        
          else if(is > ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, ieg) and (isg, ie).            
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ieg, js, je, isc, iec, jsc, jec, dir, folded)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         isg, ie, js, je, isc, iec, jsc, jec, dir, folded)
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir, folded)
          end if

          !--- Now calculate the overlapping for fold-edge. Currently we only consider about folded-north
          !--- for folded-north-edge, only need to consider to_pe's north(7) direction 
          !--- only position at NORTH and CORNER need to be considered
          if( folded_north .AND. ( position == NORTH .OR. position == CORNER) ) then
             if( domain%y(tMe)%data%begin .LE. jeg .AND. jeg .LE. domain%y(tMe)%data%end+jshift )then !fold is within domain
                dir = 7
                !--- calculate the overlapping for sending
                if( domain%x(tMe)%pos .LT. (size(domain%x(tMe)%list(:))+1)/2 )then 
                   js = domain%list(m)%y(tNbr)%compute%end+jshift;   je = js
                   if( je == jeg )then   ! fold is within domain.
                      is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
                      select case (position)
                      case(NORTH)
                         is = max(is, middle)
                         i=is; is = isg+ieg-ie; ie = isg+ieg-i
                      case(CORNER)
                         is = max(is, middle)
                         i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                      end select
                      call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                                 is, ie, js, je, isc, iec, jsc, jec, dir, .true.)           
                      is = max(is, isc); ie = min(ie, iec)
                      js = max(js, jsc); je = min(je, jec)
                      if(debug_update_level .NE. NO_CHECK .AND. ie.GE.is .AND. je.GE.js )then
                         nsend_check = nsend_check+1
                         if(nsend_check > size(checkList(:)) ) then
                            call expand_check_overlap_list(checkList, nlist)
                         endif
                         call allocate_check_overlap(checkList(nsend_check), 1)
                         call insert_check_overlap(checkList(nsend_check), domain%list(m)%pe, &
                                                   tMe, 4, ONE_HUNDRED_EIGHTY, is, ie, js, je)
                      end if
                   end if
                end if
             end if
          end if
       end if
       !--- copy the overlapping information 
       if( overlap%count > 0) then
         nsend = nsend + 1
         if(nsend > size(overlapList(:)) ) then
            call mpp_error(NOTE, 'mpp_domains_define.inc(compute_overlaps): overlapList for send is expanded')
            call expand_update_overlap_list(overlapList, nlist)
         endif
         call add_update_overlap( overlapList(nsend), overlap)
         call init_overlap_type(overlap)
       endif
    end do  ! end of send set up.

    ! copy the overlapping information into domain data structure
    if(nsend>0) then
       allocate(update%send(nsend))
       update%nsend = nsend
       do m = 1, nsend
          call add_update_overlap( update%send(m), overlapList(m) )
       enddo
    endif

    if(nsend_check>0) then
       check%nsend = nsend_check
       allocate(check%send(nsend_check))
       do m = 1, nsend_check
          call add_check_overlap( check%send(m), checkList(m) )
       enddo
    endif

    do m = 1,size(overlapList(:))
       call deallocate_overlap_type(overlapList(m))
    enddo

    if(debug_update_level .NE. NO_CHECK) then
       do m = 1,size(checkList(:))
          call deallocate_overlap_type(checkList(m))
       enddo
    endif

    isgd = isg - domain%whalo
    iegd = ieg + domain%ehalo
    jsgd = jsg - domain%shalo
    jegd = jeg + domain%nhalo

    ! begin setting up recv 
    nrecv = 0      
    nrecv_check = 0
    do list = 0,nlist-1
       m = mod( domain%pos+nlist-list, nlist )
       if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
          isc = domain%list(m)%x(1)%compute%begin; iec = domain%list(m)%x(1)%compute%end+ishift
          jsc = domain%list(m)%y(1)%compute%begin; jec = domain%list(m)%y(1)%compute%end+jshift
          !recv_e  
          dir = 1
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( domain%symmetry .AND. (position == NORTH .OR. position == CORNER ) &
               .AND. ( jsd == je .or. jed == js ) ) then
             ! --- do nothing, this point will come from other pe
          else
             if( ied.GT.ieg )then
                if( domain%x(tMe)%cyclic .AND. ie.LT.isd )then !try cyclic offset
                   is = is+ioff; ie = ie+ioff
                   call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jeg, nj)
                end if
             end if

             !--- when the north face is folded, the east halo point at right side domain will be folded.
             !--- the position should be on CORNER or NORTH
             if( jed == jeg .AND. folded_north .AND. (position == CORNER .OR. position == NORTH) &        
                  .AND. isd .GE. middle .AND. ied .LE. ieg ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                                     is, ie, js, je, isd, ied, jsd, jed-1, dir)
                is=isc; ie=iec; js=jsc; je=jec
                jsd = jed
                select case (position)
                case(NORTH)
                   i=is; is = isg+ieg-ie; ie = isg+ieg-i  
                case(CORNER)
                   i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift      
                end select
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd, jed, dir, .TRUE.)
             else 
                if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
                   ! (js, jeg) and (jsg, je).
                   call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                               is, ie, js, jeg, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
                   call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                               is, ie, jsg, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
                else
                   call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                               is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
                end if
             end if
          end if
          !recv_se     
          dir = 2 
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          need_adjust_1 = .true.; need_adjust_2 = .true.; need_adjust_3 = .true.
          if( jsd.LT.jsg )then
             if( domain%y(tMe)%cyclic .AND. js.GT.jed )then !try cyclic offset
                js = js-joff; je = je-joff
                need_adjust_1 = .false.
                if( ied.GT.ieg )then
                   if( domain%x(tMe)%cyclic .AND. ie.LT.isd )then !try cyclic offset
                      is = is+ioff; ie = ie+ioff
                      need_adjust_2 = .false.
                      if(x_cyclic_offset .NE. 0) then
                         call apply_cyclic_offset(js, je, x_cyclic_offset, jsgd, jeg, nj)
                      else if(y_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, iegd, ni)
                      end if
                   end if
                else
                   call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, ieg, ni)
                   need_adjust_3 = .false.
                end if
             end if
          end if
          if( need_adjust_3 .AND. ied.GT.ieg )then
             if( need_adjust_2 .AND. domain%x(tMe)%cyclic .AND. ie.LT.isd )then !try cyclic offset
                is = is+ioff; ie = ie+ioff
                if( need_adjust_1 .AND. jsd.GE.jsg )then
                   call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jeg, nj)
                end if
             end if
          end if
          if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
             ! (js, jeg) and (jsgd, je).  
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, jeg, isd, ied, jsd, jed, dir)
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, jsgd, je, isd, ied, jsd, jed, dir)
          else if(is > ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, iegd) and (isg, ie).            
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, iegd, js, je, isd, ied, jsd, jed, dir)
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        isg, ie, js, je, isd, ied, jsd, jed, dir)
          else
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, je, isd, ied, jsd, jed, dir)
          end if

          !recv_s      
          dir = 3
          isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
          jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          if( jsd.LT.jsg )then
             if( domain%y(tMe)%cyclic .AND. js.GT.jed )then !try cyclic offset
                js = js-joff; je = je-joff
                call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, ieg, ni)
             end if
          end if
          if(is>ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, ieg) and (isg, ie). 
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ieg, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        isg, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
          else
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
          end if

          !recv_sw 
          dir = 4
          isd = domain%x(tMe)%data%begin; ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%data%begin; jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          need_adjust_1 = .true.; need_adjust_2 = .true.; need_adjust_3 = .true.
          if( jsd.LT.jsg )then
             if( domain%y(tMe)%cyclic .AND. js.GT.jed )then !try cyclic offset
                js = js-joff; je = je-joff
                need_adjust_1 = .false.
                if( isd.LT.isg )then
                   if( domain%x(tMe)%cyclic .AND. is.GT.ied )then !try cyclic offset
                      is = is-ioff; ie = ie-ioff
                      need_adjust_2 = .false.
                      if(x_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(js, je, -x_cyclic_offset, jsgd, jeg, nj)
                      else if(y_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(is, ie, -y_cyclic_offset, isgd, ieg, ni)
                      end if
                   end if
                else            
                   call apply_cyclic_offset(is, ie, -y_cyclic_offset, isg, ieg, ni) 
                   need_adjust_3 = .false.  
                end if
             end if
          end if
          if( need_adjust_3 .AND. isd.LT.isg )then
             if( need_adjust_2 .AND. domain%x(tMe)%cyclic .AND. is.GT.ied )then !try cyclic offset
                is = is-ioff; ie = ie-ioff
                if(need_adjust_1 .AND. jsd.GE.jsg) then
                   call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jeg, nj)
                end if
             end if
          end if
          if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
             ! (js, jeg) and (jsgd, je). 
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, jeg, isd, ied, jsd, jed, dir)
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, jsgd, je, isd, ied, jsd, jed, dir)
          else if(is > ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, ieg) and (isgd, ie).  
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ieg, js, je, isd, ied, jsd, jed, dir)
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        isgd, ie, js, je, isd, ied, jsd, jed, dir)
          else
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, je, isd, ied, jsd, jed, dir)
          end if
          !recv_w      
          dir = 5
          isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( isd.LT.isg )then
             if( domain%x(tMe)%cyclic .AND. is.GT.ied )then !try cyclic offset
                is = is-ioff; ie = ie-ioff
                call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jeg, nj)
             end if
          end if
          !--- when the north face is folded, some point at j=nj will be folded.
          !--- the position should be on CORNER or NORTH
          if( jed == jeg .AND. folded_north .AND. (position == CORNER .OR. position == NORTH) & 
               .AND. ( (domain%x(tMe)%cyclic .AND. isd < isg ) .OR. ( ied .GE. middle ) ) )  then
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, je, isd, ied, jsd, jed-1, dir)
             !--- consider at j = jeg for west edge.
             !--- when the data is at corner and not symmetry, i = isg -1 will get from cyclic condition
             if( position == CORNER .AND. .NOT. domain%symmetry .AND. isd < isg ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, ied, ied, jed, jed, dir)
             end if
             is=isc; ie=iec; js=jsc; je=jec
             if(isd < isg) then
                select case (position)
                case(NORTH)
                   i=is; is = 2*isg-ie-1; ie = 2*isg-i-1     
                case(CORNER)
                   ied = ied -1 + ishift
                   i=is; is = 2*isg-ie-2+2*ishift; ie = 2*isg-i-2+2*ishift 
                end select
                if(ie .GT. domain%x(tMe)%compute%end+ishift) call mpp_error( FATAL, &
                     'mpp_domains_define.inc(compute_overlaps): west edge ubound error recv.' )     
             else 
                select case (position)
                case(NORTH)
                   i=is; is = isg+ieg-ie; ie = isg+ieg-i
                case(CORNER)
                   i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                end select
             end if
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, je, isd, ied, jed, jed, dir, .TRUE.)
          else
             if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
                ! (js, jeg) and (jsg, je).
                call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                           is, ie, js, jeg, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
                call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                           is, ie, jsg, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             else
                call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                           is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             end if
          end if

          !recv_nw     
          dir = 6
          folded = .false.
          isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          need_adjust_1 = .true.; need_adjust_2 = .true.; need_adjust_3 = .true.
          if( jed.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsd )then !try cyclic offset
                js = js+joff; je = je+joff
                need_adjust_1 = .false.
                if( isd.LT.isg )then
                   if( domain%x(tMe)%cyclic .AND. is.GE.ied )then !try cyclic offset
                      is = is-ioff; ie = ie-ioff
                      need_adjust_2 = .false.
                      if(x_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jegd, nj)
                      else if(y_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(is, ie, y_cyclic_offset, isgd, ieg, ni)
                      end if
                   end if
                else
                   call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, ieg, ni) 
                   need_adjust_3 = .false. 
                end if
             else if( folded_north )then
                folded = .TRUE.
                call get_fold_index_north(isg, ieg, jeg, ishift, position, is, ie, js, je)
             end if
          end if
          if( need_adjust_3 .AND. isd.LT.isg )then
             if( need_adjust_2 .AND. domain%x(tMe)%cyclic .AND. is.GE.ied )then !try cyclic offset
                is = is-ioff; ie = ie-ioff
                if( need_adjust_1 .AND. jed.LE.jeg )then
                   call apply_cyclic_offset(js, je, -x_cyclic_offset, jsg, jeg, nj)
                end if
             end if
          end if

          if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
             ! (js, jegd) and (jsg, je).  
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, jegd, isd, ied, jsd, jed, dir, folded)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, jsg, je, isd, ied, jsd, jed, dir, folded)
          else if(is > ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, ieg) and (isgd, ie). 
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ieg, js, je, isd, ied, jsd, jed, dir, folded)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         isgd, ie, js, je, isd, ied, jsd, jed, dir, folded)
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isd, ied, jsd, jed, dir, folded)
          end if

          !--- when north edge is folded, is will be less than isg when position is EAST and CORNER
          if(is .LT. isg .AND. domain%x(tMe)%cyclic) then
             is = is + ioff
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, is, js, je, isd, ied, jsd, jed, dir, folded )
          endif

          !recv_n      
          dir = 7
          folded = .false.
          isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( jed.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsd )then !try cyclic offset
                js = js+joff; je = je+joff
                call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, ieg, ni)
             else if( folded_north )then
                folded = .TRUE.
                call get_fold_index_north(isg, ieg, jeg, ishift, position, is, ie, js, je)
             end if
          end if
          !--- when domain symmetry and position is EAST or CORNER, the point at i=isd will 
          !--- come from two pe ( there will be only one point on one pe. ).
          if( domain%symmetry .AND. (position == EAST .OR. position == CORNER ) &
               .AND. (isd == ie .or. ied == is ) ) then
             !--- do nothing, this point will come from other pe
          else
             if(is>ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
                ! (is, ieg) and (isg, ie).
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ieg, js, je, isd, ied, jsd, jed, dir, folded, symmetry=domain%symmetry)
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            isg, ie, js, je, isd, ied, jsd, jed, dir, folded, symmetry=domain%symmetry)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd, jed, dir, folded, symmetry=domain%symmetry)
             end if
          end if
          !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
          if(is .LT. isg .AND. domain%x(tMe)%cyclic) then
             is = is + ioff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, is, js, je, isd, ied, jsd, jed, dir, folded)
          endif
          !recv_ne     
          dir = 8
          folded = .false.
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          need_adjust_1 = .true.; need_adjust_2 = .true.; need_adjust_3 = .true.
          if( jed.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsd )then !try cyclic offset
                js = js+joff; je = je+joff
                need_adjust_1 = .false.
                if( ied.GT.ieg )then
                   if( domain%x(tMe)%cyclic .AND. ie.LT.isd )then !try cyclic offset
                      is = is+ioff; ie = ie+ioff
                      need_adjust_2 = .false.
                      if(x_cyclic_offset .NE. 0) then  
                         call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jegd, nj)
                      else if(y_cyclic_offset .NE. 0) then 
                         call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, iegd, ni)
                      end if
                   end if
                else           
                   call apply_cyclic_offset(is, ie, y_cyclic_offset, isg, ieg, ni) 
                   need_adjust_3 = .false.         
                end if
             else if( folded_north )then
                folded = .TRUE.
                call get_fold_index_north(isg, ieg, jeg, ishift, position, is, ie, js, je)
             end if
          end if
          if(  need_adjust_3 .AND. ied.GT.ieg )then
             if(  need_adjust_2 .AND. domain%x(tMe)%cyclic .AND. ie.LT.isd )then !try cyclic offset
                is = is+ioff; ie = ie+ioff
                if( need_adjust_1 .AND. jed.LE.jeg)then
                   call apply_cyclic_offset(js, je, x_cyclic_offset, jsg, jeg, nj)
                end if
             end if
          end if
          if(js > je) then ! seperate into two regions due to x_cyclic_offset is nonzero, the two region are
             ! (js, jegd) and (jsg, je). 
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, jsg, je, isd, ied, jsd, jed, dir, folded)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, jegd, isd, ied, jsd, jed, dir, folded)
          else if(is > ie) then ! seperate into two regions due to y_cyclic_offset is nonzero, the two region are
             ! (is, iegd) and (isg, ie).    
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, iegd, js, je, isd, ied, jsd, jed, dir, folded)
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         isg, ie, js, je, isd, ied, jsd, jed, dir, folded)
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isd, ied, jsd, jed, dir, folded)
          end if
       endif
       !--- Now calculate the overlapping for fold-edge. Currently we only consider about folded-north
       !--- for folded-north-edge, only need to consider to_pe's north(7) direction 
       !--- only position at NORTH and CORNER need to be considered
       if( folded_north .AND. ( position == NORTH .OR. position == CORNER) ) then
          if( domain%y(tMe)%data%begin .LE. jeg .AND. jeg .LE. domain%y(tMe)%data%end )then !fold is within domain
             dir = 7
             !--- calculating overlapping for receving on north
             if( domain%x(tMe)%pos .GE. size(domain%x(tMe)%list(:))/2 )then 
                jsd = domain%y(tMe)%compute%end+jshift;   jed = jsd
                if( jed == jeg )then   ! fold is within domain.
                   isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
                   is=isc; ie=iec; js = jsc; je = jec
                   select case (position)
                   case(NORTH)
                      isd = max(isd, middle)
                      i=is; is = isg+ieg-ie; ie = isg+ieg-i
                   case(CORNER)
                      isd = max(isd, middle)
                      i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                   end select
                   call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                              is, ie, js, je, isd, ied, jsd, jed, dir, .TRUE.) 
                   is = max(is, isd); ie = min(ie, ied)
                   js = max(js, jsd); je = min(je, jed)
                   if(debug_update_level .NE. NO_CHECK .AND. ie.GE.is .AND. je.GE.js )then
                      nrecv_check = nrecv_check+1          
                      if(nrecv_check > size(checkList(:)) ) then
                         call expand_check_overlap_list(checkList, nlist)
                      endif
                      call allocate_check_overlap(checkList(nrecv_check), 1)
                      call insert_check_overlap(checkList(nrecv_check), domain%list(m)%pe, &
                                                tMe, 4, ONE_HUNDRED_EIGHTY, is, ie, js, je)
                   end if
                end if
             end if
          end if
       end if
       !--- copy the overlapping information
       if( overlap%count > 0) then
          nrecv = nrecv + 1
          if(nrecv > size(overlapList(:)) )then
            call mpp_error(NOTE, 'mpp_domains_define.inc(compute_overlaps): overlapList for recv is expanded')
            call expand_update_overlap_list(overlapList, nlist)
          endif
          call add_update_overlap( overlapList(nrecv), overlap)
          call init_overlap_type(overlap)
       endif
    enddo ! end of recv do loop

    ! copy the overlapping information into domain
    if(nrecv>0) then
       allocate(update%recv(nrecv))
       update%nrecv = nrecv
       do m = 1, nrecv
          call add_update_overlap( update%recv(m), overlapList(m) )
          do n = 1, update%recv(m)%count
             if(update%recv(m)%tileNbr(n) == domain%tile_id(tMe)) then
                if(update%recv(m)%dir(n) == 1) domain%x(tMe)%loffset = 0
                if(update%recv(m)%dir(n) == 7) domain%y(tMe)%loffset = 0
             endif
          enddo
       enddo
    endif

    if(nrecv_check>0) then
       check%nrecv = nrecv_check
       allocate(check%recv(nrecv_check))
       do m = 1, nrecv_check
          call add_check_overlap( check%recv(m), checkList(m) )
       enddo
    endif

    call deallocate_overlap_type(overlap)
    do m = 1,size(overlapList(:))
       call deallocate_overlap_type(overlapList(m))
    enddo

    if(debug_update_level .NE. NO_CHECK) then
       do m = 1,size(checkList(:))
          call deallocate_overlap_type(checkList(m))
       enddo
    endif

    deallocate(overlapList)
    deallocate(checkList)
    update=>NULL()
    check=>NULL()
    domain%initialized = .true.

  end subroutine compute_overlaps

  !####################################################################################
  subroutine compute_overlaps_fold_south( domain, position, ishift, jshift)
    !computes remote domain overlaps
    !assumes only one in each direction
    !will calculate the overlapping for T,E,C,N-cell seperately.
    type(domain2D), intent(inout) :: domain
    integer, intent(in)           :: position, ishift, jshift

    integer                          :: i, m, n, nlist, tMe, tNbr, dir
    integer                          :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd, jed
    integer                          :: isg, ieg, jsg, jeg, ioff, joff
    integer                          :: list, middle, ni, nj, isgd, iegd, jsgd, jegd
    integer                          :: ism, iem, jsm, jem, whalo, ehalo, shalo, nhalo
    logical                          :: folded
    type(overlap_type)               :: overlap
    type(overlapSpec),   pointer     :: update=>NULL()
    type(overlap_type)               :: overlapList(MAXLIST)
    type(overlap_type)               :: checkList(MAXLIST)
    type(overlapSpec),   pointer     :: check =>NULL()
    integer                          :: nsend, nrecv
    integer                          :: nsend_check, nrecv_check

    !--- since we restrict that if multiple tiles on one pe, all the tiles are limited to this pe. 
    !--- In this case, if ntiles on this pe is greater than 1, no overlapping between processor within each tile
    !--- In this case the overlapping exist only for tMe=1 and tNbr=1
    if(size(domain%x(:)) > 1) return

    !--- if there is no halo, no need to compute overlaps.
    if(domain%whalo==0 .AND. domain%ehalo==0 .AND. domain%shalo==0 .AND. domain%nhalo==0) return

    !--- when there is only one tile, n will equal to np
    nlist = size(domain%list(:))

    select case(position)
    case (CENTER)
       update => domain%update_T
       check  => NULL()
    case (CORNER)
       update => domain%update_C
       check  => domain%check_C
    case (EAST)
       update => domain%update_E
       check  => domain%check_E
    case (NORTH)
       update => domain%update_N
       check  => domain%check_N
    case default
       call mpp_error(FATAL, &
        "mpp_domains_define.inc(compute_overlaps_fold_south): the value of position should be CENTER, EAST, CORNER or NORTH")
    end select

    !--- overlap is used to store the overlapping temporarily.
    call allocate_update_overlap( overlap, MAXOVERLAP)

    !send
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec, position=position  )
    call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ni, ysize=nj, position=position  ) !cyclic offsets
    call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )
    update%xbegin = ism; update%xend = iem
    update%ybegin = jsm; update%yend = jem
    if(ASSOCIATED(check)) then
       check%xbegin  = ism; check%xend  = iem
       check%ybegin  = jsm; check%yend  = jem
    endif
    update%whalo  = domain%whalo; update%ehalo = domain%ehalo
    update%shalo  = domain%shalo; update%nhalo = domain%nhalo
    whalo         = domain%whalo; ehalo        = domain%ehalo
    shalo         = domain%shalo; nhalo        = domain%nhalo


    ioff = ni - ishift
    joff = nj - jshift
    middle = (isg+ieg)/2+1
    tMe = 1; tNbr = 1

    if(.NOT. BTEST(domain%fold,SOUTH)) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_south): "//&
            "boundary condition in y-direction should be folded-south for "//trim(domain%name))
    endif
    if(.NOT. domain%x(tMe)%cyclic) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_south): "//&
            "boundary condition in x-direction should be cyclic for "//trim(domain%name))
    endif

    if(.not. domain%symmetry) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_south): "//&
            "when south boundary is folded, the domain must be symmetry for "//trim(domain%name))
    endif

    nsend = 0
    nsend_check = 0
    do list = 0,nlist-1
       m = mod( domain%pos+list, nlist )
       if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
          !to_pe's eastern halo
          dir = 1
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
          !--- to make sure the consistence between pes
          if( (position == NORTH .OR. position == CORNER ) .AND. ( jsc == je .or. jec == js ) ) then
             !--- do nothing, this point will come from other pe
          else 
             if( ie.GT.ieg .AND. iec.LT.is )then ! cyclic is assumed
                is = is-ioff; ie = ie-ioff
             end if
             !--- when the south face is folded, the east halo point at right side domain will be folded.
             !--- the position should be on CORNER or NORTH
             if( js == jsg .AND. (position == CORNER .OR. position == NORTH) & 
                  .AND. is .GE. middle .AND. domain%list(m)%x(tNbr)%compute%end+ehalo+jshift .LE. ieg ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js+1, je, isc, iec, jsc, jec, dir)
                is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
                je = js
                select case (position)
                case(NORTH)
                   i=is; is = isg+ieg-ie; ie = isg+ieg-i  
                case(CORNER)
                   i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift      
                end select
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, .true.)
             else 
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             end if
          end if

          !to_pe's SE halo
          dir = 2
          folded = .false.
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%begin-shalo;    je = domain%list(m)%y(tNbr)%compute%begin-1
          if( ie.GT.ieg .AND. iec.LT.is )then ! cyclic is assumed
             is = is-ioff; ie = ie-ioff
          end if
          if( js.LT.jsg )then
             folded = .TRUE.
             call get_fold_index_south(isg, ieg, jsg, ishift, position, is, ie, js, je)
          end if

          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, folded)

          !to_pe's southern halo
          dir = 3
          folded = .FALSE.
          is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
          js = domain%list(m)%y(tNbr)%compute%begin-shalo; je = domain%list(m)%y(tNbr)%compute%begin-1
          folded = .FALSE. 
          if( js.LT.jsg )then
             folded = .TRUE.
             call get_fold_index_south(isg, ieg, jsg, ishift, position, is, ie, js, je)
          end if
          !--- when domain symmetry and position is EAST or CORNER, the point when isc == ie, 
          !--- no need to send, because the data on that point will come from other pe. 
          !--- come from two pe ( there will be only one point on one pe. ).
          if( (position == EAST .OR. position == CORNER ) .AND. ( isc == ie .or. iec == is ) ) then
             !--- do nothing, this point will come from other pe
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir, folded, symmetry=domain%symmetry)
          endif
          !--- when south edge is folded, ie will be less than isg when position is EAST and CORNER
          if(is .LT. isg) then
             is = is + ioff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, is, js, je, isc, iec, jsc, jec, dir, folded)
          endif

          !to_pe's SW halo
          dir = 4
          folded = .false.
          is = domain%list(m)%x(tNbr)%compute%begin-whalo; ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%begin-shalo; je = domain%list(m)%y(tNbr)%compute%begin-1
          if( isg.GT.is .AND. ie.LT.isc )then !cyclic offset
             is = is+ioff; ie = ie+ioff
          end if
          if( js.LT.jsg )then
             folded = .TRUE.
             call get_fold_index_south(isg, ieg, jsg, ishift, position, is, ie, js, je)
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, folded)
          !--- when south edge is folded, is will be less than isg when position is EAST and CORNER
          if(is .LT. isg) then
             is = is + ioff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, is, js, je, isc, iec, jsc, jec, dir, folded)
          endif

          !to_pe's western halo
          dir = 5
          is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift

          !--- to make sure the consistence between pes
          if( (position == NORTH .OR. position == CORNER ) .AND. ( jsc == je .or. jec == js ) ) then
             !--- do nothing, this point will come from other pe
          else 
             if( isg.GT.is .AND. ie.LT.isc )then ! cyclic offset
                is = is+ioff; ie = ie+ioff
             end if
             !--- when the south face is folded, some point at j=nj will be folded.
             !--- the position should be on CORNER or NORTH
             if( js == jsg .AND. (position == CORNER .OR. position == NORTH) &
                  .AND. ( domain%list(m)%x(tNbr)%compute%begin == isg .OR. domain%list(m)%x(tNbr)%compute%begin-1 .GE. middle  ) ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js+1, je, isc, iec, jsc, jec, dir)
                is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
                js = domain%list(m)%y(tNbr)%compute%begin;   je = js  
                if ( domain%list(m)%x(tNbr)%compute%begin == isg ) then
                   select case (position)
                   case(NORTH)
                      i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                   case(CORNER)
                      i=is; is = 2*isg-ie-2+2*ishift; ie = 2*isg-i-2+2*ishift  
                   end select
                   if(ie .GT. domain%x(tMe)%compute%end+ishift) call mpp_error( FATAL, &
                        'mpp_domains_define.inc(compute_overlaps_fold_south): west edge ubound error send.' )     
                else 
                   select case (position)
                   case(NORTH)
                      i=is; is = isg+ieg-ie; ie = isg+ieg-i
                   case(CORNER)
                      i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                   end select
                end if
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             end if
          endif

          !to_pe's NW halo
          dir = 6
          is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          if( isg.GT.is .AND. ie.LT.isc )then ! cyclic offset
             is = is+ioff; ie = ie+ioff
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir)

          !to_pe's northern halo
          dir = 7
          is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)

          !to_pe's NE halo
          dir = 8
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          if( ie.GT.ieg .AND. iec.LT.is )then !cyclic offset
             is = is-ioff; ie = ie-ioff
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir)

          !--- Now calculate the overlapping for fold-edge. 
          !--- only position at NORTH and CORNER need to be considered
          if( ( position == NORTH .OR. position == CORNER) ) then
             if( domain%y(tMe)%data%begin .LE. jsg .AND. jsg .LE. domain%y(tMe)%data%end+jshift )then !fold is within domain
                dir = 3
                !--- calculate the overlapping for sending
                if( domain%x(tMe)%pos .LT. (size(domain%x(tMe)%list(:))+1)/2 )then 
                   js = domain%list(m)%y(tNbr)%compute%begin;   je = js
                   if( js == jsg )then   ! fold is within domain.
                      is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
                      select case (position)
                      case(NORTH)
                         is = max(is, middle)
                         i=is; is = isg+ieg-ie; ie = isg+ieg-i
                      case(CORNER)
                         is = max(is, middle)
                         i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                      end select
                      call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                                 is, ie, js, je, isc, iec, jsc, jec, dir, .true.)           
                      is = max(is, isc); ie = min(ie, iec)
                      js = max(js, jsc); je = min(je, jec)
                      if(debug_update_level .NE. NO_CHECK .AND. ie.GE.is .AND. je.GE.js )then
                         nsend_check = nsend_check+1
                         call allocate_check_overlap(checkList(nsend_check), 1)
                         call insert_check_overlap(checkList(nsend_check), domain%list(m)%pe, &
                                                   tMe, 2, ONE_HUNDRED_EIGHTY, is, ie, js, je)
                      end if
                   end if
                end if
             end if
          end if
       end if
       !--- copy the overlapping information 
       if( overlap%count > 0) then
         nsend = nsend + 1
         if(nsend > MAXLIST) call mpp_error(FATAL,  &
             "mpp_domains_define.inc(compute_overlaps_south): nsend is greater than MAXLIST, increase MAXLIST")
         call add_update_overlap(overlapList(nsend), overlap)
         call init_overlap_type(overlap)
       endif
    end do  ! end of send set up.

    ! copy the overlapping information into domain data structure
    if(nsend>0) then
       allocate(update%send(nsend))
       update%nsend = nsend   
       do m = 1, nsend
          call add_update_overlap( update%send(m), overlapList(m) )
       enddo
    endif

    if(nsend_check>0) then
       allocate(check%send(nsend_check))
       check%nsend = nsend_check
       do m = 1, nsend_check
          call add_check_overlap( check%send(m), checkList(m) )
       enddo
    endif

    do m = 1, MAXLIST
       call deallocate_overlap_type(overlapList(m))
       if(debug_update_level .NE. NO_CHECK) call deallocate_overlap_type(checkList(m))
    enddo

    isgd = isg - domain%whalo
    iegd = ieg + domain%ehalo
    jsgd = jsg - domain%shalo
    jegd = jeg + domain%nhalo

    ! begin setting up recv 
    nrecv = 0    
    nrecv_check = 0      
    do list = 0,nlist-1
       m = mod( domain%pos+nlist-list, nlist )
       if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
          isc = domain%list(m)%x(1)%compute%begin; iec = domain%list(m)%x(1)%compute%end+ishift
          jsc = domain%list(m)%y(1)%compute%begin; jec = domain%list(m)%y(1)%compute%end+jshift
          !recv_e  
          dir = 1
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( (position == NORTH .OR. position == CORNER ) .AND. ( jsd == je .or. jed == js ) ) then
             ! --- do nothing, this point will come from other pe
          else
             if( ied.GT.ieg .AND. ie.LT.isd )then !cyclic offset
                is = is+ioff; ie = ie+ioff
             end if

             !--- when the south face is folded, the east halo point at right side domain will be folded.
             !--- the position should be on CORNER or NORTH
             if( jsd == jsg .AND. (position == CORNER .OR. position == NORTH) &        
                  .AND. isd .GE. middle .AND. ied .LE. ieg ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd+1, jed, dir)
                is=isc; ie=iec; js=jsc; je=jec
                jed = jsd
                select case (position)
                case(NORTH)
                   i=is; is = isg+ieg-ie; ie = isg+ieg-i  
                case(CORNER)
                   i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift      
                end select
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd, jed, dir, .TRUE.)
             else 
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             end if
          end if

          !recv_se     
          dir = 2 
          folded = .false.
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          if( jsd.LT.jsg )then
             folded = .true.
             call get_fold_index_south(isg, ieg, jsg, ishift, position, is, ie, js, je)
          end if
          if( ied.GT.ieg .AND. ie.LT.isd )then !cyclic offset
             is = is+ioff; ie = ie+ioff
          endif
          call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                     is, ie, js, je, isd, ied, jsd, jed, dir, folded)

          !recv_s      
          dir = 3
          folded = .false.
          isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
          jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          if( jsd.LT.jsg )then
             folded = .true.
             call get_fold_index_south(isg, ieg, jsg, ishift, position, is, ie, js, je)              
          end if
          if( (position == EAST .OR. position == CORNER ) .AND. (isd == ie .or. ied == is ) ) then
             !--- do nothing, this point will come from other pe
          else
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, je, isd, ied, jsd, jed, dir, folded, symmetry=domain%symmetry)
          end if
          !--- when south edge is folded, is will be less than isg when position is EAST and CORNER
          if(is .LT. isg ) then
             is = is + ioff
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, is, js, je, isd, ied, jsd, jed, dir, folded)
          endif

          !recv_sw 
          dir = 4
          folded = .false.
          isd = domain%x(tMe)%data%begin; ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%data%begin; jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          if( jsd.LT.jsg )then
             folded = .true.
             call get_fold_index_south(isg, ieg, jsg, ishift, position, is, ie, js, je)        
          end if
          if( isd.LT.isg .AND. is.GT.ied ) then ! cyclic offset
             is = is-ioff; ie = ie-ioff
          end if
          call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                     is, ie, js, je, isd, ied, jsd, jed, dir, folded)
          !--- when southth edge is folded, is will be less than isg when position is EAST and CORNER
          if(is .LT. isg ) then
             is = is + ioff
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, is, js, je, isd, ied, jsd, jed, dir, folded )
          endif

          !recv_w      
          dir = 5
          isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( (position == NORTH .OR. position == CORNER ) .AND. ( jsd == je .or. jed == js ) ) then
             ! --- do nothing, this point will come from other pe
          else
             if( isd.LT.isg .AND. is.GT.ied )then ! cyclic offset
                is = is-ioff; ie = ie-ioff
             end if
             !--- when the south face is folded, some point at j=nj will be folded.
             !--- the position should be on CORNER or NORTH
             if( jsd == jsg .AND. (position == CORNER .OR. position == NORTH) & 
                  .AND. ( isd < isg  .OR. ied .GE. middle ) )  then
                call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                           is, ie, js, je, isd, ied, jsd+1, jed, dir)
                is=isc; ie=iec; js=jsc; je=jec
                if(isd < isg) then
                   select case (position)
                   case(NORTH)
                      i=is; is = 2*isg-ie-1; ie = 2*isg-i-1     
                   case(CORNER)
                      ied = ied -1 + ishift
                      i=is; is = 2*isg-ie-2+2*ishift; ie = 2*isg-i-2+2*ishift 
                   end select
                   if(ie .GT. domain%x(tMe)%compute%end+ishift) call mpp_error( FATAL, &
                        'mpp_domains_define.inc(compute_overlaps): west edge ubound error recv.' )     
                else 
                   select case (position)
                   case(NORTH)
                      i=is; is = isg+ieg-ie; ie = isg+ieg-i
                   case(CORNER)
                      i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                   end select
                end if
                call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                           is, ie, js, je, isd, ied, jsd, jsd, dir, .TRUE.)
             else
                call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                           is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             end if
          endif

          !recv_nw     
          dir = 6
          isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( isd.LT.isg .AND. is.GE.ied )then !cyclic offset
             is = is-ioff; ie = ie-ioff
          endif

          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isd, ied, jsd, jed, dir)

          !recv_n      
          dir = 7
          isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)

          !recv_ne     
          dir = 8
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if(  ied.GT.ieg .AND. ie.LT.isd )then ! cyclic offset
             is = is+ioff; ie = ie+ioff
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isd, ied, jsd, jed, dir)

          !--- Now calculate the overlapping for fold-edge. 
          !--- for folded-south-edge, only need to consider to_pe's south(3) direction 
          !--- only position at NORTH and CORNER need to be considered
          if( ( position == NORTH .OR. position == CORNER) ) then
             if( domain%y(tMe)%data%begin .LE. jsg .AND. jsg .LE. domain%y(tMe)%data%end+jshift )then !fold is within domain
                dir = 3
                !--- calculating overlapping for receving on north
                if( domain%x(tMe)%pos .GE. size(domain%x(tMe)%list(:))/2 )then 
                   jsd = domain%y(tMe)%compute%begin;   jed = jsd
                   if( jsd == jsg )then   ! fold is within domain.
                      isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
                      is=isc; ie=iec; js = jsc; je = jec
                      select case (position)
                      case(NORTH)
                         isd = max(isd, middle)
                         i=is; is = isg+ieg-ie; ie = isg+ieg-i
                      case(CORNER)
                         isd = max(isd, middle)
                         i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                      end select
                      call insert_update_overlap(overlap, domain%list(m)%pe, &
                                                 is, ie, js, je, isd, ied, jsd, jed, dir, .TRUE.) 
                      is = max(is, isd); ie = min(ie, ied)
                      js = max(js, jsd); je = min(je, jed)
                      if(debug_update_level .NE. NO_CHECK .AND. ie.GE.is .AND. je.GE.js )then
                         nrecv_check = nrecv_check+1                      
                         call allocate_check_overlap(checkList(nrecv_check), 1)
                         call insert_check_overlap(checkList(nrecv_check), domain%list(m)%pe, &
                                                   tMe, 2, ONE_HUNDRED_EIGHTY, is, ie, js, je)
                      endif
                   endif
                endif
             endif
          endif
       endif
       !--- copy the overlapping information
       if( overlap%count > 0) then
          nrecv = nrecv + 1
          if(nrecv > MAXLIST) call mpp_error(FATAL,  &
             "mpp_domains_define.inc(compute_overlaps_south): nrecv is greater than MAXLIST, increase MAXLIST")
          call add_update_overlap( overlapList(nrecv), overlap)
          call init_overlap_type(overlap)
       endif
    enddo ! end of recv do loop

    ! copy the overlapping information into domain
    if(nrecv>0) then
       update%nrecv = nrecv
       allocate(update%recv(nrecv))
       do m = 1, nrecv
          call add_update_overlap( update%recv(m), overlapList(m) )
          do n = 1, update%recv(m)%count
             if(update%recv(m)%tileNbr(n) == domain%tile_id(tMe)) then
                if(update%recv(m)%dir(n) == 1) domain%x(tMe)%loffset = 0
                if(update%recv(m)%dir(n) == 7) domain%y(tMe)%loffset = 0
             endif
          enddo
       enddo
    endif

    if(nrecv_check>0) then
       check%nrecv = nrecv_check
       allocate(check%recv(nrecv_check))
       do m = 1, nrecv_check
          call add_check_overlap( check%recv(m), checkList(m) )
       enddo
    endif

    call deallocate_overlap_type(overlap)
    do m = 1, MAXLIST
       call deallocate_overlap_type(overlapList(m))
       if(debug_update_level .NE. NO_CHECK) call deallocate_overlap_type(checkList(m))
    enddo

    update => NULL()
    check=>NULL()
    domain%initialized = .true.

  end subroutine compute_overlaps_fold_south

  !####################################################################################
  subroutine compute_overlaps_fold_west( domain, position, ishift, jshift)
    !computes remote domain overlaps
    !assumes only one in each direction
    !will calculate the overlapping for T,E,C,N-cell seperately.
    type(domain2D), intent(inout) :: domain
    integer, intent(in)           :: position, ishift, jshift

    integer                          :: j, m, n, nlist, tMe, tNbr, dir
    integer                          :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd, jed
    integer                          :: isg, ieg, jsg, jeg, ioff, joff
    integer                          :: list, middle, ni, nj, isgd, iegd, jsgd, jegd
    integer                          :: ism, iem, jsm, jem, whalo, ehalo, shalo, nhalo
    logical                          :: folded
    type(overlap_type)               :: overlap
    type(overlapSpec),   pointer     :: update=>NULL()
    type(overlap_type)               :: overlapList(MAXLIST)
    type(overlap_type)               :: checkList(MAXLIST)
    type(overlapSpec),   pointer     :: check =>NULL()
    integer                          :: nsend, nrecv
    integer                          :: nsend_check, nrecv_check

    !--- since we restrict that if multiple tiles on one pe, all the tiles are limited to this pe. 
    !--- In this case, if ntiles on this pe is greater than 1, no overlapping between processor within each tile
    !--- In this case the overlapping exist only for tMe=1 and tNbr=1
    if(size(domain%x(:)) > 1) return

    !--- if there is no halo, no need to compute overlaps.
    if(domain%whalo==0 .AND. domain%ehalo==0 .AND. domain%shalo==0 .AND. domain%nhalo==0) return

    !--- when there is only one tile, n will equal to np
    nlist = size(domain%list(:))

    select case(position)
    case (CENTER)
       update => domain%update_T
       check  => NULL()
    case (CORNER)
       update => domain%update_C
       check  => domain%check_C
    case (EAST)
       update => domain%update_E
       check  => domain%check_E
    case (NORTH)
       update => domain%update_N
       check  => domain%check_N
    case default
       call mpp_error(FATAL, &
        "mpp_domains_define.inc(compute_overlaps_fold_west): the value of position should be CENTER, EAST, CORNER or NORTH")
    end select

    !--- overlap is used to store the overlapping temporarily.
    call allocate_update_overlap( overlap, MAXOVERLAP)

    !send
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec, position=position )
    call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ni, ysize=nj, position=position ) !cyclic offsets
    call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )
    update%xbegin = ism; update%xend = iem
    update%ybegin = jsm; update%yend = jem
    if(ASSOCIATED(check)) then
       check%xbegin  = ism; check%xend  = iem
       check%ybegin  = jsm; check%yend  = jem
    endif
    update%whalo  = domain%whalo; update%ehalo = domain%ehalo
    update%shalo  = domain%shalo; update%nhalo = domain%nhalo
    whalo         = domain%whalo; ehalo        = domain%ehalo
    shalo         = domain%shalo; nhalo        = domain%nhalo

    ioff = ni - ishift
    joff = nj - jshift
    middle = (jsg+jeg)/2+1
    tMe = 1; tNbr = 1

    if(.NOT. BTEST(domain%fold,WEST)) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_west): "//&
            "boundary condition in y-direction should be folded-west for "//trim(domain%name))
    endif
    if(.NOT. domain%y(tMe)%cyclic) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_west): "//&
            "boundary condition in y-direction should be cyclic for "//trim(domain%name))
    endif

    if(.not. domain%symmetry) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_west): "//&
            "when west boundary is folded, the domain must be symmetry for "//trim(domain%name))
    endif

    nsend = 0
    nsend_check = 0
    do list = 0,nlist-1
       m = mod( domain%pos+list, nlist )
       if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
          !to_pe's eastern halo
          dir = 1
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)

          !to_pe's SE halo
          dir = 2
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%begin-shalo;    je = domain%list(m)%y(tNbr)%compute%begin-1
          if( js.LT.jsg .AND. jsc.GT.je )then ! cyclic is assumed
             js = js+joff; je = je+joff
          end if

          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                       is, ie, js, je, isc, iec, jsc, jec, dir)

          !to_pe's southern halo
          dir = 3
          is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
          js = domain%list(m)%y(tNbr)%compute%begin-shalo; je = domain%list(m)%y(tNbr)%compute%begin-1
          !--- to make sure the consistence between pes
          if( (position == EAST .OR. position == CORNER ) .AND. ( isc == ie .or. iec == is ) ) then
             !--- do nothing, this point will come from other pe
          else 
             if( js.LT.jsg .AND. jsc.GT.je) then ! cyclic offset
                js = js+joff; je = je+joff
             endif

             !--- when the west face is folded, the south halo points at 
             !--- the position should be on CORNER or EAST
             if( is == isg .AND. (position == CORNER .OR. position == EAST) &
                  .AND. ( domain%list(m)%y(tNbr)%compute%begin == jsg .OR. domain%list(m)%y(tNbr)%compute%begin-1 .GE. middle  ) ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, &
                                            is+1, ie, js, je, isc, iec, jsc, jec, dir)
                is = domain%list(m)%x(tNbr)%compute%begin; ie = is
                js = domain%list(m)%y(tNbr)%compute%begin-shalo;    je = domain%list(m)%y(tNbr)%compute%begin-1
                if ( domain%list(m)%y(tNbr)%compute%begin == jsg ) then
                   select case (position)
                   case(EAST)
                      j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                   case(CORNER)
                      j=js; js = 2*jsg-je-2+2*jshift; je = 2*jsg-j-2+2*jshift  
                   end select
                   if(je .GT. domain%y(tMe)%compute%end+jshift) call mpp_error( FATAL, &
                        'mpp_domains_define.inc(compute_overlaps_fold_west: south edge ubound error send.' )     
                else 
                   select case (position)
                   case(EAST)
                      j=js; js = jsg+jeg-je; je = jsg+jeg-j
                   case(CORNER)
                      j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                   end select
                end if
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             end if
          endif

          !to_pe's SW halo
          dir = 4
          folded = .false.
          is = domain%list(m)%x(tNbr)%compute%begin-whalo; ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%begin-shalo; je = domain%list(m)%y(tNbr)%compute%begin-1
          if( jsg.GT.js .AND. je.LT.jsc )then !cyclic offset
             js = js+joff; je = je+joff
          end if
          if( is.LT.isg )then
             folded = .TRUE.
             call get_fold_index_west(jsg, jeg, isg, jshift, position, is, ie, js, je)
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, folded)
          !--- when south edge is folded, js will be less than jsg when position is EAST and CORNER
          if(js .LT. jsg) then
             js = js + joff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, js, isc, iec, jsc, jec, dir, folded)
          endif

          !to_pe's western halo
          dir = 5
          folded = .FALSE.
          is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
          if( isg.GT.is )then
             folded = .true.
             call get_fold_index_west(jsg, jeg, isg, jshift, position, is, ie, js, je)
          end if
          !--- when domain symmetry and position is EAST or CORNER, the point when isc == ie, 
          !--- no need to send, because the data on that point will come from other pe. 
          !--- come from two pe ( there will be only one point on one pe. ).
          if( (position == EAST .OR. position == CORNER ) .AND. ( jsc == je .or. jec == js ) ) then
             !--- do nothing, this point will come from other pe
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir, folded, symmetry=domain%symmetry)
          endif
          !--- when south edge is folded, ie will be less than isg when position is EAST and CORNER
          if(js .LT. jsg) then
             js = js + ioff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, js, isc, iec, jsc, jec, dir, folded)
          endif

          !to_pe's NW halo
          dir = 6
          folded = .false.
          is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          if( je.GT.jeg .AND. jec.LT.js )then ! cyclic offset
             js = js-joff; je = je-joff
          end if
          if( is.LT.isg )then
             folded = .TRUE.
             call get_fold_index_west(jsg, jeg, isg, jshift, position, is, ie, js, je)
          end if

          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, folded)

          !to_pe's northern halo
          dir = 7
          is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          !--- to make sure the consistence between pes
          if( (position == EAST .OR. position == CORNER ) .AND. ( isc == ie .or. iec == is ) ) then
             !--- do nothing, this point will come from other pe
          else 
             if( je.GT.jeg .AND. jec.LT.js) then ! cyclic offset
                js = js-joff; je = je-joff
             endif
             !--- when the west face is folded, the south halo points at 
             !--- the position should be on CORNER or EAST
             if( is == isg .AND. (position == CORNER .OR. position == EAST) &
                  .AND. ( js .GE. middle .AND. domain%list(m)%y(tNbr)%compute%end+nhalo+jshift .LE. jeg ) ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is+1, ie, js, je, isc, iec, jsc, jec, dir)
                is = domain%list(m)%x(tNbr)%compute%begin; ie = is
                js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
                select case (position)
                case(EAST)
                   j=js; js = jsg+jeg-je; je = jsg+jeg-j
                case(CORNER)
                   j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                end select
                call insert_update_overlap( overlap, domain%list(m)%pe, &
                                            is, ie, js, je, isc, iec, jsc, jec, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             end if
          endif

          !to_pe's NE halo
          dir = 8
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          if( je.GT.jeg .AND. jec.LT.js )then !cyclic offset
             js = js-joff; je = je-joff
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir)

          !--- Now calculate the overlapping for fold-edge. 
          !--- only position at EAST and CORNER need to be considered
          if( ( position == EAST .OR. position == CORNER) ) then
             if( domain%x(tMe)%compute%begin-whalo .LE. isg .AND. isg .LE. domain%x(tMe)%data%end+ishift )then !fold is within domain
                dir = 5
                !--- calculate the overlapping for sending
                if( domain%y(tMe)%pos .LT. (size(domain%y(tMe)%list(:))+1)/2 )then 
                   is = domain%list(m)%x(tNbr)%compute%begin;   ie = is
                   if( is == isg )then   ! fold is within domain.
                      js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
                      select case (position)
                      case(EAST)
                         js = max(js, middle)
                         j=js; js = jsg+jeg-je; je = jsg+jeg-j
                      case(CORNER)
                         js = max(js, middle)
                         j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                      end select
                      call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                                 is, ie, js, je, isc, iec, jsc, jec, dir, .true.)           
                      is = max(is, isc); ie = min(ie, iec)
                      js = max(js, jsc); je = min(je, jec)
                      if(debug_update_level .NE. NO_CHECK .AND. ie.GE.is .AND. je.GE.js )then
                         nsend_check = nsend_check+1
                         call allocate_check_overlap(checkList(nsend_check), 1)
                         call insert_check_overlap(checkList(nsend_check), domain%list(m)%pe, &
                                                   tMe, 3, ONE_HUNDRED_EIGHTY, is, ie, js, je)
                      end if
                   end if
                end if
             end if
          end if
       end if
       !--- copy the overlapping information 
       if( overlap%count > 0) then
         nsend = nsend + 1
         if(nsend > MAXLIST) call mpp_error(FATAL,  &
             "mpp_domains_define.inc(compute_overlaps_west): nsend is greater than MAXLIST, increase MAXLIST")
         call add_update_overlap(overlapList(nsend), overlap)
         call init_overlap_type(overlap)
       endif
    end do  ! end of send set up.

   ! copy the overlapping information into domain data structure
    if(nsend>0) then
       update%nsend = nsend  
       allocate(update%send(nsend))
       do m = 1, nsend
          call add_update_overlap( update%send(m), overlapList(m) )
       enddo
    endif

    if(nsend_check>0) then
       check%nsend = nsend_check
       allocate(check%send(nsend_check))
       do m = 1, nsend_check
          call add_check_overlap( check%send(m), checkList(m) )
       enddo
    endif
 
    do m = 1, MAXLIST
       call deallocate_overlap_type(overlapList(m))
       if(debug_update_level .NE. NO_CHECK) call deallocate_overlap_type(checkList(m))
    enddo

    isgd = isg - domain%whalo
    iegd = ieg + domain%ehalo
    jsgd = jsg - domain%shalo
    jegd = jeg + domain%nhalo

    ! begin setting up recv     
    nrecv = 0  
    nrecv_check = 0
    do list = 0,nlist-1
       m = mod( domain%pos+nlist-list, nlist )
       if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
          isc = domain%list(m)%x(1)%compute%begin; iec = domain%list(m)%x(1)%compute%end+ishift
          jsc = domain%list(m)%y(1)%compute%begin; jec = domain%list(m)%y(1)%compute%end+jshift
          !recv_e  
          dir = 1
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)

          !recv_se     
          dir = 2 
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          if( jsd.LT.jsg .AND. js.GE.jed )then ! cyclic is assumed
             js = js-joff; je = je-joff
          end if
          call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                     is, ie, js, je, isd, ied, jsd, jed, dir)

          !recv_s      
          dir = 3
          folded = .false.
          isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
          jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec

          if( (position == EAST .OR. position == CORNER ) .AND. ( isd == ie .or. ied == is ) ) then
             !--- do nothing, this point will come from other pe
          else
             if( jsd.LT.jsg .AND. js .GT. jed)then
                js = js-joff; je = je-joff
             end if
             !--- when the west face is folded, the south halo points at 
             !--- the position should be on CORNER or EAST
             if( isd == isg .AND. (position == CORNER .OR. position == EAST) &
                  .AND. ( jsd < jsg .OR. jed .GE. middle ) ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd+1, ied, jsd, jed, dir)
                is=isc; ie=iec; js=jsc; je=jec
                if(jsd<jsg) then
                   select case (position)
                   case(EAST)
                      j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                   case(CORNER)
                      j=js; js = 2*jsg-je-2+2*jshift; je = 2*jsg-j-2+2*jshift  
                   end select
                   if(je .GT. domain%y(tMe)%compute%end+jshift) call mpp_error( FATAL, &
                        'mpp_domains_define.inc(compute_overlaps_fold_west: south edge ubound error recv.' )
                else    
                   select case (position)
                   case(EAST)
                      j=js; js = jsg+jeg-je; je = jsg+jeg-j
                   case(CORNER)
                      j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                   end select
                end if
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, isd, jsd, jed, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             end if
          endif

          !recv_sw 
          dir = 4
          folded = .false.
          isd = domain%x(tMe)%data%begin; ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%data%begin; jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          if( isd.LT.isg )then
             folded = .true.
             call get_fold_index_west(jsg, jeg, isg, jshift, position, is, ie, js, je)        
          end if
          if( jsd.LT.jsg .AND. js.GT.jed ) then ! cyclic offset
             js = js-joff; je = je-joff
          end if
          call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                     is, ie, js, je, isd, ied, jsd, jed, dir, folded)
          !--- when west edge is folded, js will be less than jsg when position is EAST and CORNER
          if(js .LT. jsg ) then
             js = js + joff
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, js, isd, ied, jsd, jed, dir, folded )
          endif

          !recv_w      
          dir = 5
          folded = .false.
          isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( isd.LT.isg )then
             folded = .true.
             call get_fold_index_west(jsg, jeg, isg, jshift, position, is, ie, js, je)        
          end if
          if( (position == EAST .OR. position == CORNER ) .AND. (jsd == je .or. jed == js ) ) then
             !--- do nothing, this point will come from other pe
          else
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, je, isd, ied, jsd, jed, dir, folded, symmetry=domain%symmetry)
          end if
          !--- when west edge is folded, js will be less than jsg when position is EAST and CORNER
          if(js .LT. jsg ) then
             js = js + joff
             call insert_update_overlap(overlap, domain%list(m)%pe, &
                                        is, ie, js, js, isd, ied, jsd, jed, dir, folded)
          endif

          !recv_nw     
          dir = 6
          folded = .false.
          isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( isd.LT.isg) then
             folded = .true.
             call get_fold_index_west(jsg, jeg, isg, jshift, position, is, ie, js, je)    
          end if
          if( jed.GT.jeg .AND. je.LT.jsd )then !cyclic offset
             js = js+joff; je = je+joff
          endif

          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isd, ied, jsd, jed, dir)

          !recv_n      
          dir = 7
          folded = .false.
          isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( (position == EAST .OR. position == CORNER ) .AND. ( isd == ie .or. ied == is ) ) then
             !--- do nothing, this point will come from other pe
          else
             if( jed.GT.jeg .AND. je.LT.jsd)then
                js = js+joff; je = je+joff
             end if
             !--- when the west face is folded, the south halo points at 
             !--- the position should be on CORNER or EAST
             if( isd == isg .AND. (position == CORNER .OR. position == EAST) &
                  .AND. jsd .GE. middle .AND. jed .LE. jeg ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd+1, ied, jsd, jed, dir)
                is=isc; ie=iec; js=jsc; je=jec
                select case (position)
                case(EAST)
                   j=js; js = jsg+jeg-je; je = jsg+jeg-j
                case(CORNER)
                   j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                end select
                call insert_update_overlap( overlap, domain%list(m)%pe, &
                                            is, ie, js, je, isd, isd, jsd, jed, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             end if
          endif

          !recv_ne     
          dir = 8
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if(  jed.GT.jeg .AND. je.LT.jsd )then ! cyclic offset
             js = js+joff; je = je+joff
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isd, ied, jsd, jed, dir)

          !--- Now calculate the overlapping for fold-edge. 
          !--- for folded-south-edge, only need to consider to_pe's south(3) direction 
          !--- only position at EAST and CORNER need to be considered
          if( ( position == EAST .OR. position == CORNER) ) then
             if( domain%x(tMe)%data%begin .LE. isg .AND. isg .LE. domain%x(tMe)%data%end+ishift )then !fold is within domain
                dir = 5
                !--- calculating overlapping for receving on north
                if( domain%y(tMe)%pos .GE. size(domain%y(tMe)%list(:))/2 )then 
                   isd = domain%x(tMe)%compute%begin;   ied = isd
                   if( isd == isg )then   ! fold is within domain.
                      jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
                      is=isc; ie=iec; js = jsc; je = jec
                      select case (position)
                      case(EAST)
                         jsd = max(jsd, middle)
                         j=js; js = jsg+jeg-je; je = jsg+jeg-j
                      case(CORNER)
                         jsd = max(jsd, middle)
                         j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                      end select
                      call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                                 is, ie, js, je, isd, ied, jsd, jed, dir, .TRUE.) 
                      is = max(is, isd); ie = min(ie, ied)
                      js = max(js, jsd); je = min(je, jed)
                      if(debug_update_level .NE. NO_CHECK .AND. ie.GE.is .AND. je.GE.js )then
                         nrecv_check = nrecv_check+1   
                         call allocate_check_overlap(checkList(nrecv_check), 1)
                         call insert_check_overlap(checkList(nrecv_check), domain%list(m)%pe, &
                                                   tMe, 3, ONE_HUNDRED_EIGHTY, is, ie, js, je)                   
                      endif
                   endif
                endif
             endif
          endif
       endif
       !--- copy the overlapping information
       if( overlap%count > 0) then
          nrecv = nrecv + 1
         if(nrecv > MAXLIST) call mpp_error(FATAL,  &
             "mpp_domains_define.inc(compute_overlaps_west): nrecv is greater than MAXLIST, increase MAXLIST")
          call add_update_overlap( overlapList(nrecv), overlap)
          call init_overlap_type(overlap)
       endif
    enddo ! end of recv do loop

    ! copy the overlapping information into domain
    if(nrecv>0) then
       update%nrecv = nrecv
       allocate(update%recv(nrecv))
       do m = 1, nrecv
          call add_update_overlap( update%recv(m), overlapList(m) )
          do n = 1, update%recv(m)%count
             if(update%recv(m)%tileNbr(n) == domain%tile_id(tMe)) then
                if(update%recv(m)%dir(n) == 1) domain%x(tMe)%loffset = 0
                if(update%recv(m)%dir(n) == 7) domain%y(tMe)%loffset = 0
             endif
          enddo
       enddo
    endif

    if(nrecv_check>0) then
       check%nrecv = nrecv_check
       allocate(check%recv(nrecv_check))
       do m = 1, nrecv_check
          call add_check_overlap( check%recv(m), checkList(m) )
       enddo
    endif

    call deallocate_overlap_type(overlap)
    do m = 1, MAXLIST
       call deallocate_overlap_type(overlapList(m))
       if(debug_update_level .NE. NO_CHECK) call deallocate_overlap_type(checkList(m))
    enddo

    update=>NULL()
    check=>NULL()
    domain%initialized = .true.

  end subroutine compute_overlaps_fold_west

  !###############################################################################
  subroutine compute_overlaps_fold_east( domain, position, ishift, jshift )
    !computes remote domain overlaps
    !assumes only one in each direction
    !will calculate the overlapping for T,E,C,N-cell seperately.
    !here assume fold-east and y-cyclic boundary condition
    type(domain2D), intent(inout) :: domain
    integer, intent(in)           :: position, ishift, jshift

    integer                          :: j, m, n, nlist, tMe, tNbr, dir
    integer                          :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd
    integer                          :: jed, isg, ieg, jsg, jeg, ioff, joff
    integer                          :: list, middle, ni, nj, isgd, iegd, jsgd, jegd
    integer                          :: ism, iem, jsm, jem, whalo, ehalo, shalo, nhalo
    logical                          :: folded
    type(overlap_type)               :: overlap
    type(overlapSpec),   pointer     :: update=>NULL()
    type(overlap_type)               :: overlapList(MAXLIST)
    type(overlap_type)               :: checkList(MAXLIST)
    type(overlapSpec),   pointer     :: check =>NULL()
    integer                          :: nsend, nrecv
    integer                          :: nsend_check, nrecv_check

    !--- since we restrict that if multiple tiles on one pe, all the tiles are limited to this pe. 
    !--- In this case, if ntiles on this pe is greater than 1, no overlapping between processor within each tile
    !--- In this case the overlapping exist only for tMe=1 and tNbr=1
    if(size(domain%x(:)) > 1) return

    !--- if there is no halo, no need to compute overlaps.
    if(domain%whalo==0 .AND. domain%ehalo==0 .AND. domain%shalo==0 .AND. domain%nhalo==0) return

    !--- when there is only one tile, n will equal to np
    nlist = size(domain%list(:))

    select case(position)
    case (CENTER)
       update => domain%update_T
    case (CORNER)
       update => domain%update_C
       check  => domain%check_C
    case (EAST)
       update => domain%update_E
       check  => domain%check_E
    case (NORTH)
       update => domain%update_N
       check  => domain%check_N
    case default
       call mpp_error(FATAL, &
        "mpp_domains_define.inc(compute_overlaps_fold_east): the value of position should be CENTER, EAST, CORNER or NORTH")
    end select

    !--- overlap is used to store the overlapping temporarily.
    call allocate_update_overlap( overlap, MAXOVERLAP)

    !send
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec, position=position )
    call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ni, ysize=nj, position=position ) !cyclic offsets
    call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )
    update%xbegin = ism; update%xend = iem
    update%ybegin = jsm; update%yend = jem
    if(ASSOCIATED(check)) then
       check%xbegin  = ism; check%xend  = iem
       check%ybegin  = jsm; check%yend  = jem
    endif
    update%whalo  = domain%whalo; update%ehalo = domain%ehalo
    update%shalo  = domain%shalo; update%nhalo = domain%nhalo
    whalo         = domain%whalo; ehalo        = domain%ehalo
    shalo         = domain%shalo; nhalo        = domain%nhalo

    ioff = ni - ishift
    joff = nj - jshift
    middle = (jsg+jeg)/2+1
    tMe = 1; tNbr = 1

    if(.NOT. BTEST(domain%fold,EAST)) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_east): "//&
            "boundary condition in y-direction should be folded-east for "//trim(domain%name))
    endif
    if(.NOT. domain%y(tMe)%cyclic) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_east): "//&
            "boundary condition in y-direction should be cyclic for "//trim(domain%name))
    endif
    if(.not. domain%symmetry) then
       call mpp_error(FATAL, "mpp_domains_define.inc(compute_overlaps_fold_east): "//&
            "when east boundary is folded, the domain must be symmetry for "//trim(domain%name))
    endif

    nsend = 0
    nsend_check = 0
    do list = 0,nlist-1
       m = mod( domain%pos+list, nlist )
       if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
          !to_pe's eastern halo
          dir = 1
          folded = .false.
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
          if( ie.GT.ieg )then
             folded = .true.
             call get_fold_index_east(jsg, jeg, ieg, jshift, position, is, ie, js, je)
          end if
          !--- when domain symmetry and position is EAST or CORNER, the point when jsc == je, 
          !--- no need to send, because the data on that point will come from other pe. 
          !--- come from two pe ( there will be only one point on one pe. ).
          if( (position == EAST .OR. position == CORNER ) .AND. ( jsc == je .or. jec == js ) ) then
             !--- do nothing, this point will come from other pe
          else
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, je, isc, iec, jsc, jec, dir, folded, symmetry=domain%symmetry)
          endif
          !--- when east edge is folded, js .LT. jsg
          if(js .LT. jsg) then
             js = js + ioff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, js, isc, iec, jsc, jec, dir, folded)
          endif

          !to_pe's SE halo
          dir = 2
          folded = .FALSE.
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%begin-shalo;    je = domain%list(m)%y(tNbr)%compute%begin-1
          if( jsg.GT.js .AND. je.LT.jsc )then !try cyclic offset
             js = js+joff; je = je+joff
          end if

          if( ie.GT.ieg )then
             folded = .TRUE.
             call get_fold_index_east(jsg, jeg, ieg, jshift, position, is, ie, js, je)
          end if

          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, folded)
          !--- when east edge is folded, 
          if(js .LT. jsg) then
             js = js + joff
             call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                         is, ie, js, js, isc, iec, jsc, jec, dir, folded)
          endif

          !to_pe's southern halo
          dir = 3
          is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
          js = domain%list(m)%y(tNbr)%compute%begin-shalo; je = domain%list(m)%y(tNbr)%compute%begin-1
          !--- to make sure the consistence between pes
          if( (position == EAST .OR. position == CORNER ) .AND. ( isc == ie .or. iec == is ) ) then
             !--- do nothing, this point will come from other pe
          else 
             if( js.LT.jsg .AND. jsc.GT.je) then ! cyclic offset
                js = js+joff; je = je+joff
             endif
             !--- when the east face is folded, the south halo points at 
             !--- the position should be on CORNER or EAST
             if( ie == ieg .AND. (position == CORNER .OR. position == EAST) &
                  .AND. ( domain%list(m)%y(tNbr)%compute%begin == jsg .OR. &
                        domain%list(m)%y(tNbr)%compute%begin-1 .GE. middle  ) ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie-1, js, je, isc, iec, jsc, jec, dir)
                !--- consider at i = ieg for east edge.
                !--- when the data is at corner and not symmetry, j = jsg -1 will get from cyclic condition
                if(position == CORNER .AND. .NOT. domain%symmetry .AND. domain%list(m)%y(tNbr)%compute%begin == jsg) then
                   call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                              ie, ie, je, je, isc, iec, jsc, jec, dir, .true.)
                end if

                ie = domain%list(m)%x(tNbr)%compute%end+ishift; is = ie
                js = domain%list(m)%y(tNbr)%compute%begin-shalo;  je = domain%list(m)%y(tNbr)%compute%begin-1
                if ( domain%list(m)%y(tNbr)%compute%begin == jsg ) then
                   select case (position)
                   case(EAST)
                      j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                   case(CORNER)
                      j=js; js = 2*jsg-je-2+2*jshift; je = 2*jsg-j-2+2*jshift  
                   end select
                   if(je .GT. domain%y(tMe)%compute%end+jshift) call mpp_error( FATAL, &
                        'mpp_domains_define.inc(compute_overlaps_fold_east: south edge ubound error send.' )     
                else 
                   select case (position)
                   case(EAST)
                      j=js; js = jsg+jeg-je; je = jsg+jeg-j
                   case(CORNER)
                      j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                   end select
                end if
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             end if
          endif

          !to_pe's SW halo
          dir = 4
          is = domain%list(m)%x(tNbr)%compute%begin-whalo; ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%begin-shalo; je = domain%list(m)%y(tNbr)%compute%begin-1
          if( js.LT.jsg .AND. jsc.GT.je )then ! cyclic is assumed
             js = js+joff; je = je+joff
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir)

          !to_pe's western halo
          dir = 5
          is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)

          !to_pe's NW halo
          dir = 6
          is = domain%list(m)%x(tNbr)%compute%begin-whalo;    ie = domain%list(m)%x(tNbr)%compute%begin-1
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          if( je.GT.jeg .AND. jec.LT.js )then !cyclic offset
             js = js-joff; je = je-joff
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir)

          !to_pe's northern halo
          dir = 7
          folded = .FALSE. 
          is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end+ishift
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          !--- to make sure the consistence between pes
          if( (position == EAST .OR. position == CORNER ) .AND. ( isc == ie .or. iec == is ) ) then
             !--- do nothing, this point will come from other pe
          else 
             if( je.GT.jeg .AND. jec.LT.js) then ! cyclic offset
                js = js-joff; je = je-joff
             endif
             !--- when the east face is folded, the north halo points at 
             !--- the position should be on CORNER or EAST
             if( ie == ieg .AND. (position == CORNER .OR. position == EAST) &
                  .AND. ( js .GE. middle .AND. domain%list(m)%y(tNbr)%compute%end+nhalo+jshift .LE. jeg ) ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie-1, js, je, isc, iec, jsc, jec, dir)
                ie = domain%list(m)%x(tNbr)%compute%end+ishift; is = ie
                js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
                select case (position)
                case(EAST)
                   j=js; js = jsg+jeg-je; je = jsg+jeg-j
                case(CORNER)
                   j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                end select
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isc, iec, jsc, jec, dir, symmetry=domain%symmetry)
             end if
          endif

          !to_pe's NE halo
          dir = 8
          folded = .false.
          is = domain%list(m)%x(tNbr)%compute%end+1+ishift; ie = domain%list(m)%x(tNbr)%compute%end+ehalo+ishift
          js = domain%list(m)%y(tNbr)%compute%end+1+jshift; je = domain%list(m)%y(tNbr)%compute%end+nhalo+jshift
          if( je.GT.jeg .AND. jec.LT.js )then ! cyclic offset
             js = js-joff; je = je-joff
          end if
          if( ie.GT.ieg )then
             folded = .TRUE.
             call get_fold_index_east(jsg, jeg, ieg, jshift, position, is, ie, js, je)
          end if

          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isc, iec, jsc, jec, dir, folded)

          !--- Now calculate the overlapping for fold-edge. 
          !--- only position at EAST and CORNER need to be considered
          if( ( position == EAST .OR. position == CORNER) ) then
             if( domain%x(tMe)%data%begin .LE. ieg .AND. ieg .LE. domain%x(tMe)%data%end+ishift )then !fold is within domain
                dir = 1
                !--- calculate the overlapping for sending
                if( domain%y(tMe)%pos .LT. (size(domain%y(tMe)%list(:))+1)/2 )then 
                   ie = domain%list(m)%x(tNbr)%compute%end+ishift;   is = ie
                   if( ie == ieg )then   ! fold is within domain.
                      js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end+jshift
                      select case (position)
                      case(EAST)
                         js = max(js, middle)
                         j=js; js = jsg+jeg-je; je = jsg+jeg-j
                      case(CORNER)
                         js = max(js, middle)
                         j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                      end select
                      call insert_update_overlap(overlap, domain%list(m)%pe, &
                                                 is, ie, js, je, isc, iec, jsc, jec, dir, .true.)           
                      is = max(is, isc); ie = min(ie, iec)
                      js = max(js, jsc); je = min(je, jec)
                      if(debug_update_level .NE. NO_CHECK .AND. ie.GE.is .AND. je.GE.js )then
                         nsend_check = nsend_check+1
                         call allocate_check_overlap(checkList(nsend_check), 1)
                         call insert_check_overlap(checkList(nsend_check), domain%list(m)%pe, &
                                                   tMe, 1, ONE_HUNDRED_EIGHTY, is, ie, js, je)
                      end if
                   end if
                end if
             end if
          end if
       end if
       !--- copy the overlapping information 
       if( overlap%count > 0) then
         nsend = nsend + 1
         if(nsend > MAXLIST) call mpp_error(FATAL,  &
             "mpp_domains_define.inc(compute_overlaps_east): nsend is greater than MAXLIST, increase MAXLIST")
         call add_update_overlap(overlapList(nsend), overlap)
         call init_overlap_type(overlap)
       endif
    end do  ! end of send set up.

    ! copy the overlapping information into domain data structure
    if(nsend>0) then
       update%nsend = nsend   
       allocate(update%send(nsend))
       do m = 1, nsend
          call add_update_overlap( update%send(m), overlapList(m) )
       enddo
    endif

   if(nsend_check>0) then
       check%nsend = nsend_check
       allocate(check%send(nsend_check))
       do m = 1, nsend_check
          call add_check_overlap( check%send(m), checkList(m) )
       enddo
    endif
 
    do m = 1, MAXLIST
       call deallocate_overlap_type(overlapList(m))
       if(debug_update_level .NE. NO_CHECK) call deallocate_overlap_type(checkList(m))
    enddo

    isgd = isg - domain%whalo
    iegd = ieg + domain%ehalo
    jsgd = jsg - domain%shalo
    jegd = jeg + domain%nhalo

    ! begin setting up recv 
    nrecv = 0 
    nrecv_check = 0             
    do list = 0,nlist-1
       m = mod( domain%pos+nlist-list, nlist )
       if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
          isc = domain%list(m)%x(1)%compute%begin; iec = domain%list(m)%x(1)%compute%end+ishift
          jsc = domain%list(m)%y(1)%compute%begin; jec = domain%list(m)%y(1)%compute%end+jshift
          !recv_e  
          dir = 1
          folded = .false.
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( ied.GT.ieg )then
             folded = .true.
             call get_fold_index_east(jsg, jeg, ieg, jshift, position, is, ie, js, je)        
          end if
          if( (position == EAST .OR. position == CORNER ) .AND. (jsd == je .or. jed == js ) ) then
             !--- do nothing, this point will come from other pe
          else
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, je, isd, ied, jsd, jed, dir, folded, symmetry=domain%symmetry)
          end if
          !--- when west edge is folded, js will be less than jsg when position is EAST and CORNER
          if(js .LT. jsg ) then
             js = js + joff
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, js, isd, ied, jsd, jed, dir, folded)
          endif

          !recv_se     
          dir = 2 
          folded = .false.
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          if( ied.GT.ieg )then
             folded = .true.
             call get_fold_index_east(jsg, jeg, ieg, jshift, position, is, ie, js, je)        
          end if
          if( jsd.LT.jsg .AND. js.GT.jed ) then ! cyclic offset
             js = js-joff; je = je-joff
          end if
          call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                     is, ie, js, je, isd, ied, jsd, jed, dir, folded)
          !--- when west edge is folded, js will be less than jsg when position is EAST and CORNER
          if(js .LT. jsg ) then
             js = js + joff
             call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                        is, ie, js, js, isd, ied, jsd, jed, dir, folded )
          endif

          !recv_s      
          dir = 3
          folded = .false.
          isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
          jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec

          if( (position == EAST .OR. position == CORNER ) .AND. ( isd == ie .or. ied == is ) ) then
             !--- do nothing, this point will come from other pe
          else
             if( jsd.LT.jsg .AND. js .GT. jed)then
                js = js-joff; je = je-joff
             end if
             !--- when the east face is folded, the south halo points at 
             !--- the position should be on CORNER or EAST
             if( ied == ieg .AND. (position == CORNER .OR. position == EAST) &
                  .AND. ( jsd < jsg .OR. jed .GE. middle ) ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied-1, jsd, jed, dir)
                is=isc; ie=iec; js=jsc; je=jec
                if(jsd<jsg) then
                   select case (position)
                   case(EAST)
                      j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                   case(CORNER)
                      j=js; js = 2*jsg-je-2+2*jshift; je = 2*jsg-j-2+2*jshift  
                   end select
                   if(je .GT. domain%y(tMe)%compute%end+jshift) call mpp_error( FATAL, &
                        'mpp_domains_define.inc(compute_overlaps_fold_west: south edge ubound error recv.' )
                else    
                   select case (position)
                   case(EAST)
                      j=js; js = jsg+jeg-je; je = jsg+jeg-j
                   case(CORNER)
                      j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                   end select
                end if
                call insert_update_overlap( overlap, domain%list(m)%pe, &
                                            is, ie, js, je, ied, ied, jsd, jed, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             end if
          endif

          !recv_sw 
          dir = 4
          isd = domain%x(tMe)%data%begin; ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%data%begin; jed = domain%y(tMe)%compute%begin-1
          is=isc; ie=iec; js=jsc; je=jec
          if( jsd.LT.jsg .AND. js.GE.jed )then ! cyclic is assumed
             js = js-joff; je = je-joff
          end if
          call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                     is, ie, js, je, isd, ied, jsd, jed, dir)

          !recv_w      
          dir = 5
          isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          call insert_update_overlap( overlap, domain%list(m)%pe, &
                                      is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)

          !recv_nw     
          dir = 6
          folded = .false.
          isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if(  jed.GT.jeg .AND. je.LT.jsd )then ! cyclic offset
             js = js+joff; je = je+joff
          end if
          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isd, ied, jsd, jed, dir)

          !recv_n      
          dir = 7
          folded = .false.
          isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end+ishift
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( (position == EAST .OR. position == CORNER ) .AND. ( isd == ie .or. ied == is ) ) then
             !--- do nothing, this point will come from other pe
          else
             if( jed.GT.jeg .AND. je.LT.jsd)then
                js = js+joff; je = je+joff
             end if
             !--- when the east face is folded, the south halo points at 
             !--- the position should be on CORNER or EAST
             if( ied == ieg .AND. (position == CORNER .OR. position == EAST) &
                  .AND. jsd .GE. middle .AND. jed .LE. jeg ) then
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied-1, jsd, jed, dir)
                is=isc; ie=iec; js=jsc; je=jec
                select case (position)
                case(EAST)
                   j=js; js = jsg+jeg-je; je = jsg+jeg-j
                case(CORNER)
                   j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                end select
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, ied, ied, jsd, jed, dir, .true.)
             else
                call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                            is, ie, js, je, isd, ied, jsd, jed, dir, symmetry=domain%symmetry)
             end if
          endif

          !recv_ne     
          dir = 8
          folded = .false.
          isd = domain%x(tMe)%compute%end+1+ishift; ied = domain%x(tMe)%data%end+ishift
          jsd = domain%y(tMe)%compute%end+1+jshift; jed = domain%y(tMe)%data%end+jshift
          is=isc; ie=iec; js=jsc; je=jec
          if( ied.GT.ieg) then
             folded = .true.
             call get_fold_index_east(jsg, jeg, ieg, jshift, position, is, ie, js, je)    
          end if
          if( jed.GT.jeg .AND. je.LT.jsd )then !cyclic offset
             js = js+joff; je = je+joff
          endif

          call insert_update_overlap( overlap, domain%list(m)%pe, & 
                                      is, ie, js, je, isd, ied, jsd, jed, dir)
          !--- Now calculate the overlapping for fold-edge. 
          !--- for folded-south-edge, only need to consider to_pe's south(3) direction 
          !--- only position at EAST and CORNER need to be considered
          if( ( position == EAST .OR. position == CORNER) ) then
             if( domain%x(tMe)%data%begin .LE. ieg .AND. ieg .LE. domain%x(tMe)%data%end+ishift )then !fold is within domain
                dir = 1
                !--- calculating overlapping for receving on north
                if( domain%y(tMe)%pos .GE. size(domain%y(tMe)%list(:))/2 )then 
                   ied = domain%x(tMe)%compute%end+ishift;   isd = ied
                   if( ied == ieg )then   ! fold is within domain.
                      jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end+jshift
                      is=isc; ie=iec; js = jsc; je = jec
                      select case (position)
                      case(EAST)
                         jsd = max(jsd, middle)
                         j=js; js = jsg+jeg-je; je = jsg+jeg-j
                      case(CORNER)
                         jsd = max(jsd, middle)
                         j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
                      end select
                      call insert_update_overlap(overlap, domain%list(m)%pe, & 
                                                 is, ie, js, je, isd, ied, jsd, jed, dir, .TRUE.) 
                      is = max(is, isd); ie = min(ie, ied)
                      js = max(js, jsd); je = min(je, jed)
                      if(debug_update_level .NE. NO_CHECK .AND. ie.GE.is .AND. je.GE.js )then
                         nrecv_check = nrecv_check+1                      
                         call allocate_check_overlap(checkList(nrecv_check), 1)
                         call insert_check_overlap(checkList(nrecv_check), domain%list(m)%pe, &
                                                   tMe, 3, ONE_HUNDRED_EIGHTY, is, ie, js, je)  
                      endif
                   endif
                endif
             endif
          endif
       endif
       !--- copy the overlapping information
       if( overlap%count > 0) then
          nrecv = nrecv + 1
         if(nrecv > MAXLIST) call mpp_error(FATAL,  &
             "mpp_domains_define.inc(compute_overlaps_east): nrecv is greater than MAXLIST, increase MAXLIST")
          call add_update_overlap( overlapList(nrecv), overlap)
          call init_overlap_type(overlap)
       endif
    enddo ! end of recv do loop

    ! copy the overlapping information into domain
    if(nrecv>0) then
       update%nrecv = nrecv
       allocate(update%recv(nrecv))
       do m = 1, nrecv
          call add_update_overlap( update%recv(m), overlapList(m) )
          do n = 1, update%recv(m)%count
             if(update%recv(m)%tileNbr(n) == domain%tile_id(tMe)) then
                if(update%recv(m)%dir(n) == 1) domain%x(tMe)%loffset = 0
                if(update%recv(m)%dir(n) == 7) domain%y(tMe)%loffset = 0
             endif
          enddo
       enddo
    endif

    if(nrecv_check>0) then
       check%nrecv = nrecv_check
       allocate(check%recv(nrecv_check))
       do m = 1, nrecv_check
          call add_check_overlap( check%recv(m), checkList(m) )
       enddo
    endif

    call deallocate_overlap_type(overlap)
    do m = 1, MAXLIST
       call deallocate_overlap_type(overlapList(m))
       if(debug_update_level .NE. NO_CHECK) call deallocate_overlap_type(checkList(m))
    enddo

    update=>NULL()
    check=>NULL()

    domain%initialized = .true.

  end subroutine compute_overlaps_fold_east

  !#####################################################################################
  subroutine get_fold_index_west(jsg, jeg, isg, jshift, position, is, ie, js, je)
    integer,    intent(in) :: jsg, jeg, isg, jshift, position
    integer, intent(inout) :: is, ie, js, je
    integer                :: i, j

    select case(position)
    case(CENTER)
       j=js; js = jsg+jeg-je; je = jsg+jeg-j
       i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
    case(EAST)
       j=js; js = jsg+jeg-je; je = jsg+jeg-j
       i=is; is = 2*isg-ie; ie = 2*isg-i
    case(NORTH)
       j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
       i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
    case(CORNER)
       j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
       i=is; is = 2*isg-ie; ie = 2*isg-i
    end select

  end subroutine get_fold_index_west

  !#####################################################################################
  subroutine get_fold_index_east(jsg, jeg, ieg, jshift, position, is, ie, js, je)
    integer,    intent(in) :: jsg, jeg, ieg, jshift, position
    integer, intent(inout) :: is, ie, js, je
    integer                :: i, j

    select case(position)
    case(CENTER)
       j=js; js = jsg+jeg-je; je = jsg+jeg-j
       i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
    case(EAST)
       j=js; js = jsg+jeg-je; je = jsg+jeg-j
       i=is; is = 2*ieg-ie; ie = 2*ieg-i
    case(NORTH)
       j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
       i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
    case(CORNER)
       j=js; js = jsg+jeg-je-1+jshift; je = jsg+jeg-j-1+jshift
       i=is; is = 2*ieg-ie; ie = 2*ieg-i
    end select

  end subroutine get_fold_index_east

  !#####################################################################################
  subroutine get_fold_index_south(isg, ieg, jsg, ishift, position, is, ie, js, je)
    integer,    intent(in) :: isg, ieg, jsg, ishift, position
    integer, intent(inout) :: is, ie, js, je
    integer                :: i, j

    select case(position)
    case(CENTER)
       i=is; is = isg+ieg-ie; ie = isg+ieg-i
       j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
    case(EAST)
       i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
       j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
    case(NORTH)
       i=is; is = isg+ieg-ie; ie = isg+ieg-i
       j=js; js = 2*jsg-je; je = 2*jsg-j
    case(CORNER)
       i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
       j=js; js = 2*jsg-je; je = 2*jsg-j
    end select

  end subroutine get_fold_index_south
  !#####################################################################################
  subroutine get_fold_index_north(isg, ieg, jeg, ishift, position, is, ie, js, je)
    integer,    intent(in) :: isg, ieg, jeg, ishift, position
    integer, intent(inout) :: is, ie, js, je
    integer                :: i, j

    select case(position)
    case(CENTER)
       i=is; is = isg+ieg-ie; ie = isg+ieg-i
       j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
    case(EAST)
       i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
       j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
    case(NORTH)
       i=is; is = isg+ieg-ie; ie = isg+ieg-i
       j=js; js = 2*jeg-je; je = 2*jeg-j
    case(CORNER)
       i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
       j=js; js = 2*jeg-je; je = 2*jeg-j
    end select

  end subroutine get_fold_index_north


  !#####################################################################################
  ! add offset to the index
  subroutine apply_cyclic_offset(lstart, lend, offset, gstart, gend, gsize)
    integer, intent(inout) :: lstart, lend
    integer, intent(in   ) :: offset, gstart, gend, gsize

    lstart = lstart + offset
    if(lstart > gend)   lstart = lstart - gsize
    if(lstart < gstart) lstart = lstart + gsize
    lend = lend + offset
    if(lend > gend)     lend   = lend   - gsize
    if(lend < gstart)   lend   = lend   + gsize

    return

  end subroutine apply_cyclic_offset

  !###################################################################################
  ! this routine setup the overlapping for mpp_update_domains for arbitrary halo update.
  ! should be the halo size defined in mpp_define_domains.
  ! xhalo_out, yhalo_out should not be exactly the same as xhalo_in, yhalo_in
  ! currently we didn't consider about tripolar grid situation, because in the folded north
  ! region, the overlapping is specified through list of points, not through rectangular.
  ! But will return back to solve this problem in the future.
  subroutine set_overlaps(domain, overlap_in, overlap_out, whalo_out, ehalo_out, shalo_out, nhalo_out)
    type(domain2d),    intent(in)    :: domain
    type(overlapSpec), intent(in)    :: overlap_in
    type(overlapSpec), intent(inout) :: overlap_out
    integer,           intent(in)    :: whalo_out, ehalo_out, shalo_out, nhalo_out
    integer                          :: nlist, m, n, isoff, ieoff, jsoff, jeoff, rotation
    integer                          :: whalo_in, ehalo_in, shalo_in, nhalo_in
    integer                          :: dir
    type(overlap_type)               :: overlap 
    type(overlap_type), allocatable  :: send(:), recv(:)
    type(overlap_type), pointer      :: ptrIn  => NULL()
    integer                          :: nsend, nrecv, nsend_in, nrecv_in

    if( domain%fold .NE. 0) call mpp_error(FATAL, &
         "mpp_domains_define.inc(set_overlaps): folded domain is not implemented for arbitrary halo update, contact developer")

    whalo_in = domain%whalo
    ehalo_in = domain%ehalo
    shalo_in = domain%shalo
    nhalo_in = domain%nhalo

    if( .NOT. domain%initialized) call mpp_error(FATAL, &
         "mpp_domains_define.inc: domain is not defined yet")

    nlist = size(domain%list(:))
    isoff = whalo_in - abs(whalo_out)
    ieoff = ehalo_in - abs(ehalo_out)
    jsoff = shalo_in - abs(shalo_out)
    jeoff = nhalo_in - abs(nhalo_out)

    nsend = 0
    nsend_in = overlap_in%nsend
    nrecv_in = overlap_in%nrecv
    if(nsend_in>0) allocate(send(nsend_in))
    if(nrecv_in>0) allocate(recv(nrecv_in))
    call allocate_update_overlap(overlap, MAXOVERLAP)

    overlap_out%whalo  = whalo_out
    overlap_out%ehalo  = ehalo_out
    overlap_out%shalo  = shalo_out
    overlap_out%nhalo  = nhalo_out
    overlap_out%xbegin = overlap_in%xbegin
    overlap_out%xend   = overlap_in%xend
    overlap_out%ybegin = overlap_in%ybegin
    overlap_out%yend   = overlap_in%yend

    !--- setting up overlap.
    do m = 1, nsend_in
       ptrIn  => overlap_in%send(m)
       if(ptrIn%count .LE. 0) call mpp_error(FATAL, &
          "mpp_domains_define.inc(set_overlaps): number of overlap for send should be a positive number for"//trim(domain%name) )
       do n = 1, ptrIn%count
          dir = ptrIn%dir(n)
          rotation = ptrIn%rotation(n)
          select case(dir)
          case(1) ! to_pe's eastern halo
             if(ehalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, 0, 0, n, dir, rotation)
             else if(ehalo_out<0) then
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, 0, 0, n, dir, rotation)
             end if
          case(2) ! to_pe's southeast halo
             if(ehalo_out>0 .AND. shalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, jsoff, 0, n, dir, rotation)
             else if(ehalo_out<0 .AND. shalo_out < 0) then ! three parts: southeast, south and east.
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, 0, shalo_out, n, dir, rotation)
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, jsoff, 0, n, dir-1, rotation)
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, 0, shalo_out, n, dir+1, rotation)
             end if
          case(3) ! to_pe's southern halo
             if(shalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, 0, jsoff, 0, n, dir, rotation)
             else if(shalo_out<0) then
                call set_single_overlap(ptrIn, overlap, 0, 0, 0, shalo_out, n, dir, rotation)
             end if
          case(4) ! to_pe's southwest halo
             if(whalo_out>0 .AND. shalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, isoff, 0, jsoff, 0, n, dir, rotation)
             else if(whalo_out<0 .AND. shalo_out < 0) then
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, 0, shalo_out, n, dir, rotation)
                call set_single_overlap(ptrIn, overlap, isoff, 0, 0, shalo_out, n, dir-1, rotation)
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, jsoff, 0, n, dir+1, rotation)
             end if
          case(5) ! to_pe's western halo
             if(whalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, isoff, 0, 0, 0, n, dir, rotation)
             else if(whalo_out<0) then
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, 0, 0, n, dir, rotation)
             end if
          case(6) ! to_pe's northwest halo
             if(whalo_out>0 .AND. nhalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, isoff, 0, 0, -jeoff, n, dir, rotation)
             else if(whalo_out<0 .AND. nhalo_out < 0) then
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, -nhalo_out, 0, n, dir, rotation)
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, 0, -jeoff, n, dir-1, rotation)
                call set_single_overlap(ptrIn, overlap, isoff, 0, -nhalo_out, 0, n, dir+1, rotation)
             end if
          case(7) ! to_pe's northern halo
             if(nhalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, 0, 0, -jeoff, n, dir, rotation)
             else if(nhalo_out<0) then
                call set_single_overlap(ptrIn, overlap, 0, 0, -nhalo_out, 0, n, dir, rotation)
             end if
          case(8) ! to_pe's northeast halo
             if(ehalo_out>0 .AND. nhalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, 0, -jeoff, n, dir, rotation)
             else if(ehalo_out<0 .AND. nhalo_out < 0) then
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, -nhalo_out, 0, n, dir, rotation)
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, -nhalo_out, 0, n, dir-1, rotation)
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, 0, -jeoff, n, 1, rotation)
             end if
          end select
       end do ! do n = 1, ptrIn%count
       if(overlap%count>0) then       
          nsend = nsend+1
          call add_update_overlap(send(nsend), overlap)
          call init_overlap_type(overlap)
       endif
    end do ! end do list = 0, nlist-1

    if(nsend>0) then
       overlap_out%nsend = nsend
       allocate(overlap_out%send(nsend));
       do n = 1, nsend
          call add_update_overlap(overlap_out%send(n), send(n) )
       enddo
    endif

    !--------------------------------------------------
    !                     recving
    !---------------------------------------------------
    overlap%count = 0
    nrecv = 0
    do m = 1, nrecv_in
       ptrIn  => overlap_in%recv(m)
       if(ptrIn%count .LE. 0) call mpp_error(FATAL, &
          "mpp_domains_define.inc(set_overlaps): number of overlap for recv should be a positive number")
       overlap%count = 0
       do n = 1, ptrIn%count
          dir = ptrIn%dir(n)
          rotation = ptrIn%rotation(n)  
          select case(dir)
          case(1) ! eastern halo
             if(ehalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, 0, 0, n, dir)
             else if(ehalo_out<0) then
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, 0, 0, n, dir)
             end if
          case(2) ! southeast halo
             if(ehalo_out>0 .AND. shalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, jsoff, 0, n, dir)
             else if(ehalo_out<0 .AND. shalo_out < 0) then
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, 0, shalo_out, n, dir)
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, jsoff, 0, n, dir-1)
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, 0, shalo_out, n, dir+1)
             end if
          case(3) ! southern halo
             if(shalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, 0, jsoff, 0, n, dir)
             else if(shalo_out<0) then
                call set_single_overlap(ptrIn, overlap, 0, 0, 0, shalo_out, n, dir)
             end if
          case(4) ! southwest halo
             if(whalo_out>0 .AND. shalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, isoff, 0, jsoff, 0, n, dir)
             else if(whalo_out<0 .AND. shalo_out < 0) then
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, 0, shalo_out, n, dir)
                call set_single_overlap(ptrIn, overlap, isoff, 0, 0, shalo_out, n, dir-1)
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, jsoff, 0, n, dir+1)
             end if
          case(5) ! western halo
             if(whalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, isoff, 0, 0, 0, n, dir)
             else if(whalo_out<0) then
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, 0, 0, n, dir)
             end if
          case(6) ! northwest halo
             if(whalo_out>0 .AND. nhalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, isoff, 0, 0, -jeoff, n, dir)
             else if(whalo_out<0 .AND. nhalo_out < 0) then
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, -nhalo_out, 0, n, dir)
                call set_single_overlap(ptrIn, overlap, 0, whalo_out, 0, -jeoff, n, dir-1)
                call set_single_overlap(ptrIn, overlap, isoff, 0, -nhalo_out, 0, n, dir+1)
             end if
          case(7) ! northern halo
             if(nhalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, 0, 0, -jeoff, n, dir)
             else if(nhalo_out<0) then
                call set_single_overlap(ptrIn, overlap, 0, 0, -nhalo_out, 0, n, dir)
             end if
          case(8) ! northeast halo
             if(ehalo_out>0 .AND. nhalo_out > 0) then
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, 0, -jeoff, n, dir)
             else if(ehalo_out<0 .AND. nhalo_out < 0) then
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, -nhalo_out, 0, n, dir)
                call set_single_overlap(ptrIn, overlap, 0, -ieoff, -nhalo_out, 0, n, dir-1)
                call set_single_overlap(ptrIn, overlap, -ehalo_out, 0, 0, -jeoff, n, 1)
             end if
          end select
       end do ! do n = 1, ptrIn%count
       if(overlap%count>0) then       
          nrecv = nrecv+1
          call add_update_overlap(recv(nrecv), overlap)
          call init_overlap_type(overlap)
       endif
    end do ! end do list = 0, nlist-1

    if(nrecv>0) then
       overlap_out%nrecv = nrecv
       allocate(overlap_out%recv(nrecv));
       do n = 1, nrecv
          call add_update_overlap(overlap_out%recv(n), recv(n) )
       enddo
    endif

    call deallocate_overlap_type(overlap)
    do n = 1, nsend_in
       call deallocate_overlap_type(send(n))
    enddo
    do n = 1, nrecv_in
       call deallocate_overlap_type(recv(n))
    enddo
    if(allocated(send)) deallocate(send)
    if(allocated(recv)) deallocate(recv)
    ptrIn => NULL()


  end subroutine set_overlaps

 !##############################################################################
  subroutine set_single_overlap(overlap_in, overlap_out, isoff, ieoff, jsoff, jeoff, index, dir, rotation)
    type(overlap_type),    intent(in) :: overlap_in
    type(overlap_type), intent(inout) :: overlap_out
    integer,               intent(in) :: isoff, jsoff, ieoff, jeoff       
    integer,               intent(in) :: index      
    integer,               intent(in) :: dir
    integer, optional,     intent(in) :: rotation
    integer                           :: rotate
    integer                           :: count

    if( overlap_out%pe == NULL_PE ) then
       overlap_out%pe  = overlap_in%pe
    else
       if(overlap_out%pe .NE. overlap_in%pe) call mpp_error(FATAL,  &
            "mpp_domains_define.inc(set_single_overlap): mismatch of pe between overlap_in and overlap_out")
    endif

    overlap_out%count = overlap_out%count + 1
    count = overlap_out%count 
    if(count > MAXOVERLAP) call mpp_error(FATAL, &
         "set_single_overlap: number of overlap is greater than MAXOVERLAP, increase MAXOVERLAP")
    rotate = ZERO
    if(present(rotation)) rotate = rotation
    overlap_out%rotation  (count) = overlap_in%rotation(index)
    overlap_out%is_refined(count) = overlap_in%is_refined(index)
    overlap_out%dir       (count) = dir
    overlap_out%tileMe    (count) = overlap_in%tileMe(index)   
    overlap_out%tileNbr   (count) = overlap_in%tileNbr(index)   

    select case(rotate)
    case(ZERO)
       overlap_out%is(count)      = overlap_in%is(index) + isoff
       overlap_out%ie(count)      = overlap_in%ie(index) + ieoff
       overlap_out%js(count)      = overlap_in%js(index) + jsoff
       overlap_out%je(count)      = overlap_in%je(index) + jeoff
    case(NINETY)
       overlap_out%is(count)      = overlap_in%is(index) - jeoff
       overlap_out%ie(count)      = overlap_in%ie(index) - jsoff
       overlap_out%js(count)      = overlap_in%js(index) + isoff
       overlap_out%je(count)      = overlap_in%je(index) + ieoff
    case(MINUS_NINETY)
       overlap_out%is(count)      = overlap_in%is(index) + jsoff
       overlap_out%ie(count)      = overlap_in%ie(index) + jeoff
       overlap_out%js(count)      = overlap_in%js(index) - ieoff
       overlap_out%je(count)      = overlap_in%je(index) - isoff
    case default
       call mpp_error(FATAL, "mpp_domains_define.inc: the value of rotation should be ZERO, NINETY or MINUS_NINETY")
    end select

  end subroutine set_single_overlap

  !###################################################################################
  !--- compute the overlapping between tiles for the T-cell. 
  subroutine define_contact_point( domain, position, num_contact, tile1, tile2, align1, align2, &
       refine1, refine2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,        &
       isgList, iegList, jsgList, jegList )
    type(domain2D),     intent(inout) :: domain
    integer,               intent(in) :: position          
    integer,               intent(in) :: num_contact      ! number of contact regions
    integer, dimension(:), intent(in) :: tile1, tile2     ! tile number
    integer, dimension(:), intent(in) :: align1, align2   ! align direction of contact region
    real,    dimension(:), intent(in) :: refine1, refine2 ! refinement between tiles
    integer, dimension(:), intent(in) :: istart1, iend1   ! i-index in tile_1 of contact region
    integer, dimension(:), intent(in) :: jstart1, jend1   ! j-index in tile_1 of contact region
    integer, dimension(:), intent(in) :: istart2, iend2   ! i-index in tile_2 of contact region
    integer, dimension(:), intent(in) :: jstart2, jend2   ! j-index in tile_2 of contact region
    integer, dimension(:), intent(in) :: isgList, iegList ! i-global domain of each tile
    integer, dimension(:), intent(in) :: jsgList, jegList ! j-global domain of each tile

    integer              :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer              :: isc1, iec1, jsc1, jec1, isc2, iec2, jsc2, jec2
    integer              :: isd1, ied1, jsd1, jed1, isd2, ied2, jsd2, jed2
    integer              :: is, ie, js, je, ioff, joff, isoff, ieoff, jsoff, jeoff
    integer              :: ntiles, max_contact
    integer              :: nlist, list, m, n, l, count, numS, numR
    integer              :: whalo, ehalo, shalo, nhalo
    integer              :: t1, t2, tt, pos
    integer              :: ntileMe, ntileNbr, tMe, tNbr, tileMe, dir
    integer              :: nxd, nyd, nxc, nyc, ism, iem, jsm, jem
    integer              :: dirlist(8)
    !--- is2Send and is1Send will figure out the overlapping for sending from current pe.
    !--- is1Recv and iscREcv will figure out the overlapping for recving onto current pe.
    integer, dimension(4*num_contact) :: is1Send, ie1Send, js1Send, je1Send
    integer, dimension(4*num_contact) :: is2Send, ie2Send, js2Send, je2Send
    integer, dimension(4*num_contact) :: is2Recv, ie2Recv, js2Recv, je2Recv
    integer, dimension(4*num_contact) :: is1Recv, ie1Recv, js1Recv, je1Recv
    integer, dimension(4*num_contact) :: align1Recv, align2Recv, align1Send, align2Send
    real,    dimension(4*num_contact) :: refineRecv, refineSend
    integer, dimension(4*num_contact) :: rotateSend, rotateRecv, tileSend, tileRecv
    integer                           :: nsend, nrecv, nsend2, nrecv2
    type(contact_type), dimension(domain%ntiles)            :: eCont, wCont, sCont, nCont
    type(overlap_type), dimension(0:size(domain%list(:))-1) :: overlapSend, overlapRecv

    if( position .NE. CENTER )  call mpp_error(FATAL,  "mpp_domains_define.inc: " //&
         "routine define_contact_point can only be used to calculate overlapping for cell center.")

    ntiles  = domain%ntiles

    eCont(:)%ncontact = 0; 

    do n = 1, ntiles
       eCont(n)%ncontact = 0; sCont(n)%ncontact = 0; wCont(n)%ncontact = 0; nCont(n)%ncontact = 0;
       allocate(eCont(n)%tile(num_contact), wCont(n)%tile(num_contact) )
       allocate(nCont(n)%tile(num_contact), sCont(n)%tile(num_contact) )
       allocate(eCont(n)%align1(num_contact), eCont(n)%align2(num_contact) )
       allocate(wCont(n)%align1(num_contact), wCont(n)%align2(num_contact) )
       allocate(sCont(n)%align1(num_contact), sCont(n)%align2(num_contact) )
       allocate(nCont(n)%align1(num_contact), nCont(n)%align2(num_contact) )
       allocate(eCont(n)%refine1(num_contact), eCont(n)%refine2(num_contact) )
       allocate(wCont(n)%refine1(num_contact), wCont(n)%refine2(num_contact) )
       allocate(sCont(n)%refine1(num_contact), sCont(n)%refine2(num_contact) )
       allocate(nCont(n)%refine1(num_contact), nCont(n)%refine2(num_contact) )
       allocate(eCont(n)%is1(num_contact), eCont(n)%ie1(num_contact), eCont(n)%js1(num_contact), eCont(n)%je1(num_contact))
       allocate(eCont(n)%is2(num_contact), eCont(n)%ie2(num_contact), eCont(n)%js2(num_contact), eCont(n)%je2(num_contact))
       allocate(wCont(n)%is1(num_contact), wCont(n)%ie1(num_contact), wCont(n)%js1(num_contact), wCont(n)%je1(num_contact))
       allocate(wCont(n)%is2(num_contact), wCont(n)%ie2(num_contact), wCont(n)%js2(num_contact), wCont(n)%je2(num_contact))
       allocate(sCont(n)%is1(num_contact), sCont(n)%ie1(num_contact), sCont(n)%js1(num_contact), sCont(n)%je1(num_contact))
       allocate(sCont(n)%is2(num_contact), sCont(n)%ie2(num_contact), sCont(n)%js2(num_contact), sCont(n)%je2(num_contact))
       allocate(nCont(n)%is1(num_contact), nCont(n)%ie1(num_contact), nCont(n)%js1(num_contact), nCont(n)%je1(num_contact))
       allocate(nCont(n)%is2(num_contact), nCont(n)%ie2(num_contact), nCont(n)%js2(num_contact), nCont(n)%je2(num_contact))
    end do

    !--- set up the east, south, west and north contact for each tile.
    do n = 1, num_contact
       t1 = tile1(n)
       t2 = tile2(n)
       select case(align1(n))
       case (EAST)
          call fill_contact( eCont(t1), t2, istart1(n), iend1(n), jstart1(n), jend1(n), istart2(n), iend2(n), &
               jstart2(n), jend2(n), align1(n), align2(n), refine1(n), refine2(n))
       case (WEST)
          call fill_contact( wCont(t1), t2, istart1(n), iend1(n), jstart1(n), jend1(n), istart2(n), iend2(n), &
               jstart2(n), jend2(n), align1(n), align2(n), refine1(n), refine2(n))
       case (SOUTH)
          call fill_contact( sCont(t1), t2, istart1(n), iend1(n), jstart1(n), jend1(n), istart2(n), iend2(n), &
               jstart2(n), jend2(n), align1(n), align2(n), refine1(n), refine2(n))
       case (NORTH)
          call fill_contact( nCont(t1), t2, istart1(n), iend1(n), jstart1(n), jend1(n), istart2(n), iend2(n), &
               jstart2(n), jend2(n), align1(n), align2(n), refine1(n), refine2(n))
       end select
       select case(align2(n))
       case (EAST)
          call fill_contact( eCont(t2), t1, istart2(n), iend2(n), jstart2(n), jend2(n), istart1(n), iend1(n), &
               jstart1(n), jend1(n), align2(n), align1(n), refine2(n), refine1(n))
       case (WEST)
          call fill_contact( wCont(t2), t1, istart2(n), iend2(n), jstart2(n), jend2(n), istart1(n), iend1(n), &
               jstart1(n), jend1(n), align2(n), align1(n), refine2(n), refine1(n))
       case (SOUTH)
          call fill_contact( sCont(t2), t1, istart2(n), iend2(n), jstart2(n), jend2(n), istart1(n), iend1(n), &
               jstart1(n), jend1(n), align2(n), align1(n), refine2(n), refine1(n))
       case (NORTH)
          call fill_contact( nCont(t2), t1, istart2(n), iend2(n), jstart2(n), jend2(n), istart1(n), iend1(n), &
               jstart1(n), jend1(n), align2(n), align1(n), refine2(n), refine1(n))
       end select
    end do

    !--- the tile number of current pe, halo size 
    whalo = domain%whalo
    ehalo = domain%ehalo
    shalo = domain%shalo
    nhalo = domain%nhalo

    !--- find if there is an extra point in x and y direction depending on position
    nlist = size(domain%list(:))

    max_contact = 4*num_contact ! should be enough

    ntileMe = size(domain%x(:))
    refineSend = 1; refineRecv = 1

    !--------------------------------------------------------------------------------------------------
    !    loop over each tile on current domain to set up the overlapping for each tile
    !--------------------------------------------------------------------------------------------------
    !--- first check the overlap within the tiles.
    do n = 1, domain%update_T%nsend
       pos = domain%update_T%send(n)%pe - mpp_root_pe() 
       call add_update_overlap(overlapSend(pos), domain%update_T%send(n) )
    enddo
    do n = 1, domain%update_T%nrecv
       pos = domain%update_T%recv(n)%pe - mpp_root_pe() 
       call add_update_overlap(overlapRecv(pos), domain%update_T%recv(n) )
    enddo

    call mpp_get_memory_domain(domain, ism, iem, jsm, jem)
    domain%update_T%xbegin = ism;   domain%update_T%xend  = iem
    domain%update_T%ybegin = jsm;   domain%update_T%yend  = jem
    domain%update_T%whalo  = whalo; domain%update_T%ehalo = ehalo
    domain%update_T%shalo  = shalo; domain%update_T%nhalo = nhalo

    do tMe = 1, ntileMe
       tileMe = domain%tile_id(tMe)
       rotateSend = ZERO; rotateRecv = ZERO

       !--- loop over all the contact region to figure out the index for overlapping region.
       count = 0
       do n = 1, eCont(tileMe)%ncontact  ! east contact
          count = count+1
          tileRecv(count)   = eCont(tileMe)%tile(n);    tileSend(count)   = eCont(tileMe)%tile(n)
          align1Recv(count) = eCont(tileMe)%align1(n);  align2Recv(count) = eCont(tileMe)%align2(n)
          align1Send(count) = eCont(tileMe)%align1(n);  align2Send(count) = eCont(tileMe)%align2(n)
          refineSend(count) = eCont(tileMe)%refine2(n); refineRecv(count) = eCont(tileMe)%refine1(n)
          is1Recv(count)    = eCont(tileMe)%is1(n) + 1; ie1Recv(count)    = is1Recv(count) + ehalo - 1  
          js1Recv(count)    = eCont(tileMe)%js1(n);     je1Recv(count)    = eCont(tileMe)%je1(n)
          select case(eCont(tileMe)%align2(n))
          case ( WEST )  ! w <-> e 
             is2Recv(count) = eCont(tileMe)%is2(n);     ie2Recv(count) = is2Recv(count) + ehalo - 1
             js2Recv(count) = eCont(tileMe)%js2(n);     je2Recv(count) = eCont(tileMe)%je2(n)
             ie1Send(count) = eCont(tileMe)%is1(n);     is1Send(count) = ie1Send(count) - whalo + 1  
             js1Send(count) = eCont(tileMe)%js1(n);     je1Send(count) = eCont(tileMe)%je1(n)
             ie2Send(count) = eCont(tileMe)%is2(n) - 1; is2Send(count) = ie2Send(count) - whalo + 1  
             js2Send(count) = eCont(tileMe)%js2(n);     je2Send(count) = eCont(tileMe)%je2(n)
          case ( SOUTH ) ! s <-> e
             rotateRecv(count) = NINETY;                rotateSend(count) = MINUS_NINETY
             js2Recv(count) = eCont(tileMe)%js2(n);     je2Recv(count) = js2Recv(count) + ehalo -1
             is2Recv(count) = eCont(tileMe)%is2(n);     ie2Recv(count) = eCont(tileMe)%ie2(n)
             ie1Send(count) = eCont(tileMe)%is1(n);     is1Send(count) = ie1Send(count) - shalo + 1  
             js1Send(count) = eCont(tileMe)%js1(n);     je1Send(count) = eCont(tileMe)%je1(n)
             is2Send(count) = eCont(tileMe)%is2(n);     ie2Send(count) = eCont(tileMe)%ie2(n)
             je2Send(count) = eCont(tileMe)%js2(n) - 1; js2Send(count) = je2Send(count) - shalo + 1  
          end select
       end do

       do n = 1, sCont(tileMe)%ncontact  ! south contact 
          count = count+1
          tileRecv(count)   = sCont(tileMe)%tile(n);    tileSend(count)   = sCont(tileMe)%tile(n)
          align1Recv(count) = sCont(tileMe)%align1(n);  align2Recv(count) = sCont(tileMe)%align2(n);
          align1Send(count) = sCont(tileMe)%align1(n);  align2Send(count) = sCont(tileMe)%align2(n);
          refineSend(count) = sCont(tileMe)%refine2(n); refineRecv(count) = sCont(tileMe)%refine1(n)
          is1Recv(count)    = sCont(tileMe)%is1(n);     ie1Recv(count)    = sCont(tileMe)%ie1(n)
          je1Recv(count)    = sCont(tileMe)%js1(n) - 1; js1Recv(count)    = je1Recv(count) - shalo + 1
          select case(sCont(tileMe)%align2(n))
          case ( NORTH ) ! n <-> s
             is2Recv(count) = sCont(tileMe)%is2(n);   ie2Recv(count) = sCont(tileMe)%ie2(n)
             je2Recv(count) = sCont(tileMe)%je2(n);   js2Recv(count) = je2Recv(count) - shalo + 1
             is1Send(count) = sCont(tileMe)%is1(n);   ie1Send(count) = sCont(tileMe)%ie1(n)
             js1Send(count) = sCont(tileMe)%js1(n);   je1Send(count) = js1Send(count) + nhalo -1
             is2Send(count) = sCont(tileMe)%is2(n);   ie2Send(count) = sCont(tileMe)%ie2(n)
             js2Send(count) = sCont(tileMe)%je2(n)+1; je2Send(count) = js2Send(count) + nhalo - 1
          case ( EAST )  ! e <-> s
             rotateRecv(count) = MINUS_NINETY;        rotateSend(count) = NINETY
             ie2Recv(count) = sCont(tileMe)%ie2(n);   is2Recv(count) = ie2Recv(count) - shalo + 1
             js2Recv(count) = sCont(tileMe)%js2(n);   je2Recv(count) = sCont(tileMe)%je2(n)
             is1Send(count) = sCont(tileMe)%is1(n);   ie1Send(count) = sCont(tileMe)%ie1(n)
             js1Send(count) = sCont(tileMe)%js1(n);   je1Send(count) = js1Send(count) + ehalo - 1
             is2Send(count) = sCont(tileMe)%ie2(n)+1; ie2Send(count) = is2Send(count) + ehalo - 1
             js2Send(count) = sCont(tileMe)%js2(n);   je2Send(count) = sCont(tileMe)%je2(n)
          end select
       end do

       do n = 1, wCont(tileMe)%ncontact  ! west contact
          count = count+1
          tileRecv(count)   = wCont(tileMe)%tile(n);    tileSend(count)   = wCont(tileMe)%tile(n)
          align1Recv(count) = wCont(tileMe)%align1(n);  align2Recv(count) = wCont(tileMe)%align2(n);
          align1Send(count) = wCont(tileMe)%align1(n);  align2Send(count) = wCont(tileMe)%align2(n);
          refineSend(count) = wCont(tileMe)%refine2(n); refineRecv(count) = wCont(tileMe)%refine1(n)
          ie1Recv(count)    = wCont(tileMe)%is1(n) - 1; is1Recv(count)    = ie1Recv(count) - whalo + 1
          js1Recv(count)    = wCont(tileMe)%js1(n);     je1Recv(count)    = wCont(tileMe)%je1(n)
          select case(wCont(tileMe)%align2(n))
          case ( EAST )  ! e <-> w
             ie2Recv(count) = wCont(tileMe)%ie2(n);   is2Recv(count) = ie2Recv(count) - whalo + 1
             js2Recv(count) = wCont(tileMe)%js2(n);   je2Recv(count) = wCont(tileMe)%je2(n) 
             is1Send(count) = wCont(tileMe)%is1(n);   ie1Send(count) = is1Send(count) + ehalo - 1 
             js1Send(count) = wCont(tileMe)%js1(n);   je1Send(count) = wCont(tileMe)%je1(n) 
             is2Send(count) = wCont(tileMe)%ie2(n)+1; ie2Send(count) = is2Send(count) + ehalo - 1 
             js2Send(count) = wCont(tileMe)%js2(n);   je2Send(count) = wCont(tileMe)%je2(n) 
          case ( NORTH ) ! n <-> w
             rotateRecv(count) = NINETY;              rotateSend(count) = MINUS_NINETY
             je2Recv(count) = wCont(tileMe)%je2(n);   js2Recv(count) = je2Recv(count) - whalo + 1
             is2Recv(count) = wCont(tileMe)%is2(n);   ie2Recv(count) = wCont(tileMe)%ie2(n)
             is1Send(count) = wCont(tileMe)%is1(n);   ie1Send(count) = is1Send(count) + nhalo - 1
             js1Send(count) = wCont(tileMe)%js1(n);   je1Send(count) = wCont(tileMe)%je1(n)
             js2Send(count) = wCont(tileMe)%je2(n)+1; je2Send(count) = js2Send(count) + nhalo - 1
             is2Send(count) = wCont(tileMe)%is2(n);   ie2Send(count) = wCont(tileMe)%ie2(n)
          end select
       end do

       do n = 1, nCont(tileMe)%ncontact  ! north contact 
          count = count+1
          tileRecv(count)   = nCont(tileMe)%tile(n);   tileSend(count)   = nCont(tileMe)%tile(n)
          align1Recv(count) = nCont(tileMe)%align1(n); align2Recv(count) = nCont(tileMe)%align2(n);
          align1Send(count) = nCont(tileMe)%align1(n); align2Send(count) = nCont(tileMe)%align2(n);
          refineSend(count) = nCont(tileMe)%refine2(n); refineRecv(count) = nCont(tileMe)%refine1(n)
          is1Recv(count)    = nCont(tileMe)%is1(n);    ie1Recv(count)    = nCont(tileMe)%ie1(n)
          js1Recv(count)    = nCont(tileMe)%je1(n)+1;  je1Recv(count)    = js1Recv(count) + nhalo - 1
          select case(nCont(tileMe)%align2(n))
          case ( SOUTH ) ! s <-> n
             is2Recv(count) = nCont(tileMe)%is2(n);   ie2Recv(count) = nCont(tileMe)%ie2(n)
             js2Recv(count) = nCont(tileMe)%js2(n);   je2Recv(count) = js2Recv(count) + nhalo - 1
             is1Send(count) = nCont(tileMe)%is1(n);   ie1Send(count) = nCont(tileMe)%ie1(n)
             je1Send(count) = nCont(tileMe)%je1(n);   js1Send(count) = je1Send(count) - shalo + 1
             is2Send(count) = nCont(tileMe)%is2(n);   ie2Send(count) = nCont(tileMe)%ie2(n)
             je2Send(count) = nCont(tileMe)%js2(n)-1; js2Send(count) = je2Send(count) - shalo + 1
          case ( WEST )  ! w <-> n
             rotateRecv(count) = MINUS_NINETY;        rotateSend(count) = NINETY
             is2Recv(count) = nCont(tileMe)%ie2(n);   ie2Recv(count) = is2Recv(count) + nhalo - 1
             js2Recv(count) = nCont(tileMe)%js2(n);   je2Recv(count) = nCont(tileMe)%je2(n)
             is1Send(count) = nCont(tileMe)%is1(n);   ie1Send(count) = nCont(tileMe)%ie1(n)
             je1Send(count) = nCont(tileMe)%je1(n);   js1Send(count) = je1Send(count) - whalo + 1
             ie2Send(count) = nCont(tileMe)%is2(n)-1; is2Send(count) = ie2Send(count) - whalo + 1
             js2Send(count) = nCont(tileMe)%js2(n);   je2Send(count) = nCont(tileMe)%je2(n)
          end select
       end do

       numS = count
       numR = count
       !--- figure out the index for corner overlapping,
       !--- fill_corner_contact will be updated to deal with the situation that there are multiple tiles on 
       !--- each side of six sides of cubic grid.
       if(.NOT. domain%rotated_ninety)  then
          call fill_corner_contact(eCont, sCont, wCont, nCont, isgList, iegList, jsgList, jegList, numR, numS, &
               tileRecv, tileSend, is1Recv, ie1Recv, js1Recv, je1Recv, is2Recv, ie2Recv,   &
               js2Recv, je2Recv, is1Send, ie1Send, js1Send, je1Send, is2Send, ie2Send,     &
               js2Send, je2Send, align1Recv, align2Recv, align1Send, align2Send,           &
               whalo, ehalo, shalo, nhalo, tileMe )
       end if

       isc = domain%x(tMe)%compute%begin; iec = domain%x(tMe)%compute%end
       jsc = domain%y(tMe)%compute%begin; jec = domain%y(tMe)%compute%end

       !--- compute the overlapping for send.
       do n = 1, numS
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             ntileNbr = size(domain%list(m)%x(:))
             do tNbr = 1, ntileNbr
                if( domain%list(m)%tile_id(tNbr) .NE. tileSend(n) ) cycle
                isc1 = max(isc, is1Send(n)); iec1 = min(iec, ie1Send(n))
                jsc1 = max(jsc, js1Send(n)); jec1 = min(jec, je1Send(n))  
                if( isc1 > iec1 .OR. jsc1 > jec1 ) cycle
                !--- loop over 8 direction to get the overlapping starting from east with clockwise.
                do dir = 1, 8
                   !--- get the to_pe's data domain.
                   select case ( dir )
                   case ( 1 )  ! eastern halo
                      if( align2Send(n) .NE. EAST ) cycle
                      isd = domain%list(m)%x(tNbr)%compute%end+1; ied = domain%list(m)%x(tNbr)%compute%end+ehalo
                      jsd = domain%list(m)%y(tNbr)%compute%begin; jed = domain%list(m)%y(tNbr)%compute%end
                   case ( 2 )  ! southeast halo
                      isd = domain%list(m)%x(tNbr)%compute%end+1; ied = domain%list(m)%x(tNbr)%compute%end+ehalo
                      jsd = domain%list(m)%y(tNbr)%compute%begin-shalo;    jed = domain%list(m)%y(tNbr)%compute%begin-1
                   case ( 3 )  ! southern halo
                      if( align2Send(n) .NE. SOUTH ) cycle
                      isd = domain%list(m)%x(tNbr)%compute%begin; ied = domain%list(m)%x(tNbr)%compute%end
                      jsd = domain%list(m)%y(tNbr)%compute%begin-shalo;    jed = domain%list(m)%y(tNbr)%compute%begin-1
                   case ( 4 )  ! southwest halo
                      isd = domain%list(m)%x(tNbr)%compute%begin-whalo;    ied = domain%list(m)%x(tNbr)%compute%begin-1
                      jsd = domain%list(m)%y(tNbr)%compute%begin-shalo;    jed = domain%list(m)%y(tNbr)%compute%begin-1
                   case ( 5 )  ! western halo
                      if( align2Send(n) .NE. WEST ) cycle
                      isd = domain%list(m)%x(tNbr)%compute%begin-whalo;    ied = domain%list(m)%x(tNbr)%compute%begin-1
                      jsd = domain%list(m)%y(tNbr)%compute%begin; jed = domain%list(m)%y(tNbr)%compute%end
                   case ( 6 )  ! northwest halo
                      isd = domain%list(m)%x(tNbr)%compute%begin-whalo;    ied = domain%list(m)%x(tNbr)%compute%begin-1
                      jsd = domain%list(m)%y(tNbr)%compute%end+1; jed = domain%list(m)%y(tNbr)%compute%end+nhalo
                   case ( 7 )  ! northern halo
                      if( align2Send(n) .NE. NORTH ) cycle
                      isd = domain%list(m)%x(tNbr)%compute%begin; ied = domain%list(m)%x(tNbr)%compute%end
                      jsd = domain%list(m)%y(tNbr)%compute%end+1; jed = domain%list(m)%y(tNbr)%compute%end+nhalo
                   case ( 8 )  ! northeast halo
                      isd = domain%list(m)%x(tNbr)%compute%end+1; ied = domain%list(m)%x(tNbr)%compute%end+ehalo
                      jsd = domain%list(m)%y(tNbr)%compute%end+1; jed = domain%list(m)%y(tNbr)%compute%end+nhalo
                   end select
                   isd = max(isd, is2Send(n)); ied = min(ied, ie2Send(n))
                   jsd = max(jsd, js2Send(n)); jed = min(jed, je2Send(n))  
                   if( isd > ied .OR. jsd > jed ) cycle
                   ioff = 0; joff = 0
                   nxd = ied - isd + 1
                   nyd = jed - jsd + 1
                   select case ( align2Send(n) )
                   case ( WEST, EAST )
                      joff = jsd - js2Send(n)
                      joff = get_refine_count( joff, refineSend(n), .true., dir )
                      nyd  = get_refine_count( nyd, refineSend(n), .false., dir  )
                   case ( SOUTH, NORTH )
                      ioff = isd - is2Send(n) 
                      ioff = get_refine_count( ioff, refineSend(n), .true., dir  )
                      nxd  = get_refine_count( nxd, refineSend(n), .false., dir  )
                   end select

                   !--- get the index in current pe.
                   select case ( rotateSend(n) )
                   case ( ZERO )
                      isc2 = is1Send(n) + ioff; iec2 = isc2 + nxd - 1
                      jsc2 = js1Send(n) + joff; jec2 = jsc2 + nyd - 1
                   case ( NINETY )                     ! N -> W or S -> E
                      iec2 = ie1Send(n) - joff; isc2 = iec2 - nyd + 1 
                      jsc2 = js1Send(n) + ioff; jec2 = jsc2 + nxd - 1
                   case ( MINUS_NINETY )               ! W -> N or E -> S
                      isc2 = is1Send(n) + joff; iec2 = isc2 + nyd - 1
                      jec2 = je1Send(n) - ioff; jsc2 = jec2 - nxd + 1 
                   end select
                   is = max(isc1,isc2); ie = min(iec1,iec2)
                   js = max(jsc1,jsc2); je = min(jec1,jec2)
                   if(ie.GE.is .AND. je.GE.js )then
                      if(.not. associated(overlapSend(m)%tileMe)) call allocate_update_overlap(overlapSend(m), MAXOVERLAP)
                      call insert_mosaic_update_overlap(overlapSend(m), domain%list(m)%pe, tMe, tNbr, &
                                                        is, ie, js, je, dir, rotateSend(n), .true., refineSend(n) .NE. 1 )
                   endif
                end do ! end do dir = 1, 8
             end do ! end do tNbr = 1, ntileNbr
          end do ! end do list = 0, nlist-1
       end do ! end do n = 1, numS

       !--- compute the overlapping for recv.
       do n = 1, numR
          do list = 0, nlist-1
             m = mod( domain%pos+nlist-list, nlist )
             ntileNbr = size(domain%list(m)%x(:))
             do tNbr = 1, ntileNbr
                if( domain%list(m)%tile_id(tNbr) .NE. tileRecv(n) ) cycle
                isc = domain%list(m)%x(tNbr)%compute%begin; iec = domain%list(m)%x(tNbr)%compute%end
                jsc = domain%list(m)%y(tNbr)%compute%begin; jec = domain%list(m)%y(tNbr)%compute%end
                isc = max(isc, is2Recv(n)); iec = min(iec, ie2Recv(n))
                jsc = max(jsc, js2Recv(n)); jec = min(jec, je2Recv(n))      
                if( isc > iec .OR. jsc > jec ) cycle
                !--- find the offset for this overlapping.
                ioff = 0; joff = 0
                nxc = iec - isc + 1; nyc = jec - jsc + 1
                select case ( align2Recv(n) )
                case ( WEST, EAST )
                   joff = jsc - js2Recv(n)
                   joff = get_refine_count( joff, 1/refineRecv(n), .true., dir  )
                   nyc  = get_refine_count( nyc, 1/refineRecv(n), .false., dir  )
                case ( NORTH, SOUTH )
                   ioff = isc - is2Recv(n) 
                   ioff = get_refine_count( ioff, 1/refineRecv(n), .true., dir  )
                   nxc  = get_refine_count( nxc, 1/refineRecv(n), .false., dir  )
                end select

                !--- get the index in current pe.
                select case ( rotateRecv(n) )
                case ( ZERO )
                   isd1 = is1Recv(n) + ioff; ied1 = isd1 + nxc - 1
                   jsd1 = js1Recv(n) + joff; jed1 = jsd1 + nyc - 1
                case ( NINETY )                      ! N -> W or S -> E
                   isd1 = is1Recv(n);        ied1 = isd1 + nyc - 1
                   jed1 = je1Recv(n) - ioff; jsd1 = jed1 - nxc + 1
                case ( MINUS_NINETY )                ! W -> N or E -> S
                   ied1 = ie1Recv(n) - joff; isd1 = ied1 - nyc + 1
                   jsd1 = js1Recv(n);        jed1 = jsd1 + nxc - 1
                end select

                !--- loop over 8 direction to get the overlapping starting from east with clockwise.
                do dir = 1, 8
                   select case ( dir )
                   case ( 1 )  ! eastern halo
                      if( align1Recv(n) .NE. EAST ) cycle
                      isd2 = domain%x(tMe)%compute%end+1; ied2 = domain%x(tMe)%data%end
                      jsd2 = domain%y(tMe)%compute%begin; jed2 = domain%y(tMe)%compute%end
                   case ( 2 )  ! southeast halo
                      isd2 = domain%x(tMe)%compute%end+1; ied2 = domain%x(tMe)%data%end
                      jsd2 = domain%y(tMe)%data%begin;    jed2 = domain%y(tMe)%compute%begin-1
                   case ( 3 )  ! southern halo
                      if( align1Recv(n) .NE. SOUTH ) cycle
                      isd2 = domain%x(tMe)%compute%begin; ied2 = domain%x(tMe)%compute%end
                      jsd2 = domain%y(tMe)%data%begin;    jed2 = domain%y(tMe)%compute%begin-1
                   case ( 4 )  ! southwest halo
                      isd2 = domain%x(tMe)%data%begin;    ied2 = domain%x(tMe)%compute%begin-1
                      jsd2 = domain%y(tMe)%data%begin;    jed2 = domain%y(tMe)%compute%begin-1
                   case ( 5 )  ! western halo
                      if( align1Recv(n) .NE. WEST ) cycle
                      isd2 = domain%x(tMe)%data%begin;    ied2 = domain%x(tMe)%compute%begin-1
                      jsd2 = domain%y(tMe)%compute%begin; jed2 = domain%y(tMe)%compute%end
                   case ( 6 )  ! northwest halo
                      isd2 = domain%x(tMe)%data%begin;    ied2 = domain%x(tMe)%compute%begin-1
                      jsd2 = domain%y(tMe)%compute%end+1; jed2 = domain%y(tMe)%data%end
                   case ( 7 )  ! northern halo
                      if( align1Recv(n) .NE. NORTH ) cycle
                      isd2 = domain%x(tMe)%compute%begin; ied2 = domain%x(tMe)%compute%end
                      jsd2 = domain%y(tMe)%compute%end+1; jed2 = domain%y(tMe)%data%end
                   case ( 8 )  ! northeast halo
                      isd2 = domain%x(tMe)%compute%end+1; ied2 = domain%x(tMe)%data%end
                      jsd2 = domain%y(tMe)%compute%end+1; jed2 = domain%y(tMe)%data%end
                   end select
                   is = max(isd1,isd2); ie = min(ied1,ied2)
                   js = max(jsd1,jsd2); je = min(jed1,jed2)
                   if(ie.GE.is .AND. je.GE.js )then
                      if(.not. associated(overlapRecv(m)%tileMe)) call allocate_update_overlap(overlapRecv(m), MAXOVERLAP)
                      call insert_mosaic_update_overlap(overlapRecv(m), domain%list(m)%pe, tMe, tNbr, &
                                                        is, ie, js, je, dir, rotateRecv(n), .true., refineSend(n) .NE. 1 )
                      count                              = overlapRecv(m)%count
                      overlapRecv(m)%isMe        (count) = is
                      overlapRecv(m)%ieMe        (count) = ie
                      overlapRecv(m)%jsMe        (count) = js
                      overlapRecv(m)%jeMe        (count) = je
                      !--- get the index on the pe sending data when there is refinement
                      if( overlapRecv(m)%is_refined(count) ) then
                         select case ( align1Recv(n) )
                         case ( WEST, EAST )
                            isoff = 0;               ieoff = ie - is + 1
                            jsoff = js - js1Recv(n); jeoff = je - js1Recv(n) + 1
                            jsoff = get_refine_count( jsoff, refineRecv(n), .true., dir  )
                            jeoff = get_refine_count( jeoff, refineRecv(n), .false., dir  ) - 1
                         case ( NORTH, SOUTH )
                            isoff = is - is1Recv(n); ieoff = ie - is1Recv(n) + 1
                            jsoff = 0;               jeoff = je - js
                            isoff = get_refine_count( isoff, refineRecv(n), .true., dir  )
                            ieoff = get_refine_count( ieoff, refineRecv(n), .false., dir  ) - 1
                         end select

                         !--- get the index in neighbor pe.
                         select case ( rotateRecv(n) )
                         case ( ZERO )
                            isd2 = is2Recv(n) + isoff; ied2 = is2Recv(n) + ieoff
                            jsd2 = js2Recv(n) + jsoff; jed2 = js2Recv(n) + jeoff
                         case ( NINETY )                      ! N -> W or S -> E
                            isd2 = ie2Recv(n) - jeoff; ied2 = ie2Recv(n) - jsoff
                            jsd2 = js2Recv(n);         jed2 = js2Recv(n) + ieoff - isoff
                         case ( MINUS_NINETY )                ! W -> N or E -> S
                            isd2 = is2Recv(n);         ied2 = is2Recv(n) + jeoff - jsoff
                            jsd2 = je2Recv(n) - ieoff; jed2 = je2Recv(n) - isoff
                         end select
                         isd2= max(isd2,isc); ied2 = min(ied2, iec); jsd2= max(jsd2,jsc); jed2 = min(jed2, jec)
                         if( isd2 > ied2 .OR. jsd2 > jed2 ) call mpp_error(FATAL, &
                              "mpp_domains_define.inc:  the neighbor index is out of range")
                         overlapRecv(m)%is(count) = isd2
                         overlapRecv(m)%ie(count) = ied2
                         overlapRecv(m)%js(count) = jsd2
                         overlapRecv(m)%je(count) = jed2
                      end if
                   endif
                end do ! end do dir = 1, 8
             end do ! end do tNbr = 1, ntileNbr
          end do ! end do list = 0, nlist-1
       end do ! end do n = 1, numR
    end do   ! end do tMe = 1, ntileMe

    !--- copy the overlapping information into domain data
    nsend = 0; nsend2 = 0
    do list = 0, nlist-1
       m = mod( domain%pos+list, nlist )
       if(overlapSend(m)%count>0) nsend = nsend + 1
    enddo

    dirlist(1) = 1; dirlist(2) = 3; dirlist(3) = 5; dirlist(4) = 7
    dirlist(5) = 2; dirlist(6) = 4; dirlist(7) = 6; dirlist(8) = 8
    
    ! copy the overlap information into domain.
    if(nsend >0) then
       if(associated(domain%update_T%send)) then
          do m = 1, domain%update_T%nsend
             call deallocate_overlap_type(domain%update_T%send(m))
          enddo
          deallocate(domain%update_T%send)
       endif
       domain%update_T%nsend = nsend
       allocate(domain%update_T%send(nsend))
       do list = 0, nlist-1
          m = mod( domain%pos+list, nlist )
          ntileNbr = size(domain%list(m)%x(:))  
          !--- for the send, the list should be in tileNbr order and dir order to be consistent with Recv
          if(overlapSend(m)%count > 0) then
             nsend2 = nsend2+1
             if(nsend2>nsend) call mpp_error(FATAL, &
                  "mpp_domains_define.inc(define_contact_point): nsend2 is greater than nsend")
             call allocate_update_overlap(domain%update_T%send(nsend2), overlapSend(m)%count)

             do tNbr = 1, ntileNbr
                do tt = 1, ntileMe
                   if(domain%list(m)%pe == domain%pe) then ! own processor
                      tMe = tNbr+tt-1
                      if(tMe > ntileMe) tMe = tMe - ntileMe
                   else
                      tMe = tt
                   end if
                   do n = 1, 8  ! loop over 8 direction
                      do l = 1, overlapSend(m)%count
                         if(overlapSend(m)%tileMe(l) .NE. tMe) cycle
                         if(overlapSend(m)%tileNbr(l) .NE. tNbr) cycle
                         if(overlapSend(m)%dir(l) .NE. dirlist(n) ) cycle
                         call insert_mosaic_update_overlap(domain%update_T%send(nsend2), overlapSend(m)%pe,                     &
                              overlapSend(m)%tileMe(l), overlapSend(m)%tileNbr(l), overlapSend(m)%is(l), overlapSend(m)%ie(l),  &
                              overlapSend(m)%js(l), overlapSend(m)%je(l), overlapSend(m)%dir(l), overlapSend(m)%rotation(l),    &
                              overlapSend(m)%from_contact(l),  overlapSend(m)%is_refined(l) )                   
                      end do
                   end do
                end do
             end do
          end if
       enddo
    endif

    if(nsend2 .NE. nsend) call mpp_error(FATAL, &
         "mpp_domains_define.inc(define_contact_point): nsend2 does not equal to nsend")

    nrecv = 0; nrecv2 = 0
    do list = 0, nlist-1
       m = mod( domain%pos+list, nlist )
       if(overlapRecv(m)%count>0) nrecv = nrecv + 1
    enddo

    if(nrecv >0) then
       if(associated(domain%update_T%recv)) then
          do m = 1, domain%update_T%nrecv
             call deallocate_overlap_type(domain%update_T%recv(m))
          enddo
          deallocate(domain%update_T%recv)
       endif
       domain%update_T%nrecv = nrecv
       allocate(domain%update_T%recv(nrecv))

       do list = 0, nlist-1
          m = mod( domain%pos+nlist-list, nlist )
          ntileNbr = size(domain%list(m)%x(:))   
          if(overlapRecv(m)%count > 0) then
             nrecv2 = nrecv2 + 1
             if(nrecv2>nrecv) call mpp_error(FATAL, &
                   "mpp_domains_define.inc(define_contact_point): nrecv2 is greater than nrecv")
             call allocate_update_overlap(domain%update_T%recv(nrecv2), overlapRecv(m)%count)
             do tMe = 1, ntileMe
                do tt = 1, ntileNbr
                   !--- make sure the same order tile for different pe count
                   if(domain%list(m)%pe == domain%pe) then ! own processor
                      tNbr = tMe+tt-1
                      if(tNbr>ntileNbr) tNbr = tNbr - ntileNbr
                   else
                      tNbr = tt
                   end if
                   do n = 1, 8  ! loop over 8 direction
                      do l = 1, overlapRecv(m)%count
                         if(overlapRecv(m)%tileMe(l) .NE. tMe) cycle
                         if(overlapRecv(m)%tileNbr(l) .NE. tNbr) cycle
                         if(overlapRecv(m)%dir(l) .NE. dirlist(n) ) cycle
                         call insert_mosaic_update_overlap(domain%update_T%recv(nrecv2), overlapRecv(m)%pe, &
                              overlapRecv(m)%tileMe(l), overlapRecv(m)%tileNbr(l), overlapRecv(m)%is(l), overlapRecv(m)%ie(l),  &
                              overlapRecv(m)%js(l), overlapRecv(m)%je(l), overlapRecv(m)%dir(l), overlapRecv(m)%rotation(l),    &
                              overlapRecv(m)%from_contact(l),  overlaprecv(m)%is_refined(l) )                   
                         count = domain%update_T%recv(nrecv2)%count
                         domain%update_T%recv(nrecv2)%isMe        (count) = overlapRecv(m)%isMe        (l)
                         domain%update_T%recv(nrecv2)%ieMe        (count) = overlapRecv(m)%ieMe        (l)
                         domain%update_T%recv(nrecv2)%jsMe        (count) = overlapRecv(m)%jsMe        (l)
                         domain%update_T%recv(nrecv2)%jeMe        (count) = overlapRecv(m)%jeMe        (l)
                         domain%update_T%recv(nrecv2)%index       (count) = overlapRecv(m)%index       (l)
                      end do
                   end do
                end do
             end do
          end if
       end do
    endif

    if(nrecv2 .NE. nrecv) call mpp_error(FATAL, &
    "mpp_domains_define.inc(define_contact_point): nrecv2 does not equal to nrecv")

    do m = 0,nlist-1
       call deallocate_overlap_type(overlapSend(m))
       call deallocate_overlap_type(overlapRecv(m))
    enddo
    !--- release memory
    do n = 1, ntiles
       deallocate(eCont(n)%tile, wCont(n)%tile, sCont(n)%tile, nCont(n)%tile )
       deallocate(eCont(n)%align1, wCont(n)%align1, sCont(n)%align1, nCont(n)%align1)
       deallocate(eCont(n)%align2, wCont(n)%align2, sCont(n)%align2, nCont(n)%align2)
       deallocate(eCont(n)%refine1, wCont(n)%refine1, sCont(n)%refine1, nCont(n)%refine1)
       deallocate(eCont(n)%refine2, wCont(n)%refine2, sCont(n)%refine2, nCont(n)%refine2)
       deallocate(eCont(n)%is1, eCont(n)%ie1, eCont(n)%js1, eCont(n)%je1 )
       deallocate(eCont(n)%is2, eCont(n)%ie2, eCont(n)%js2, eCont(n)%je2 )
       deallocate(wCont(n)%is1, wCont(n)%ie1, wCont(n)%js1, wCont(n)%je1 )
       deallocate(wCont(n)%is2, wCont(n)%ie2, wCont(n)%js2, wCont(n)%je2 )
       deallocate(sCont(n)%is1, sCont(n)%ie1, sCont(n)%js1, sCont(n)%je1 )
       deallocate(sCont(n)%is2, sCont(n)%ie2, sCont(n)%js2, sCont(n)%je2 )
       deallocate(nCont(n)%is1, nCont(n)%ie1, nCont(n)%js1, nCont(n)%je1 )
       deallocate(nCont(n)%is2, nCont(n)%ie2, nCont(n)%js2, nCont(n)%je2 )
    end do

    !--- define the refined overlapping and group the overlapping 
    !--- between same tiles( in different direction or different pe ). 
    !--- set up the index position in the return buffer.
    call define_refine_overlap( domain, CENTER )  

    domain%initialized = .true.

  end subroutine define_contact_point

  !--- currently assume no grid is divided over two pes.
  !--- we may remove this restriction in the future. 
  function get_refine_count(num , refine, use_floor, dir )
    integer, intent(in) :: num
    real,    intent(in) :: refine
    logical, intent(in) :: use_floor
    integer, intent(in) :: dir
    integer             :: get_refine_count
    real, parameter     :: EPSLN = 1.e-10
    real                :: rval

    rval = num*refine

    if(mod(dir,2) == 1) then ! WEST, EAST, SOUTH, NORTH
       get_refine_count = nint(rval)
       if(abs(rval-get_refine_count) > EPSLN) call mpp_error(FATAL, &
            "get_refine_count: The two adjacent tiles with refinement does not align well in E/W/S/N direction")
    else
       get_refine_count = floor(rval)
       if(use_floor) then
          if( rval - get_refine_count > 1 - EPSLN ) get_refine_count = get_refine_count + 1 ! Deal with truncation error
       else
          if(rval - get_refine_count  > EPSLN) get_refine_count = get_refine_count + 1 ! Deal with truncation error
       end if
    end if


    return

  end function get_refine_count

  !####################################################################################
  !--- The following will define the refined overlapping and group the overlapping 
  !--- between same tiles( in different direction or different pe ). 
  !--- set up the index position in the return buffer.
  subroutine define_refine_overlap(domain, position)
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: position

    logical                        :: new_overlap
    integer,           parameter   :: MAXOVERLAP = 100
    integer, dimension(MAXOVERLAP) :: isMe, ieMe, jsMe, jeMe, isNbr, ieNbr, jsNbr, jeNbr
    integer, dimension(MAXOVERLAP) :: rotation, tileList, dirList
    integer                        :: count, tMe, ntileMe, nlist, pos
    integer                        :: m, n, l, dir, istart, iend, dirdiff
    type(overlap_type), pointer    :: OverPtr => NULL()
    type(overlapSpec), pointer     :: update => NULL()

    select case(position)
    case (CENTER)
       update => domain%update_T
    case (CORNER)
       update => domain%update_C
    case (EAST)
       update => domain%update_E
    case (NORTH)
       update => domain%update_N
    case default
       call mpp_error(FATAL, &
         "mpp_domains_define.inc(define_refine_overlap): position should be CENTER, CORNER, EAST or NORTH")
    end select

    ntileMe = size(domain%x(:))
    nlist = size(domain%list(:))    

    allocate(update%rSpec(ntileMe) )

    do tMe = 1, ntileMe
       count = 0
       do m = 1, update%nrecv
          overPtr => update%recv(m)
          if( overptr%count == 0)cycle
          pos = overPtr%pe - mpp_root_pe()
          do l = 1, overPtr%count
             if(OverPtr%is_refined(l) .AND. overPtr%tileMe(l) == tMe ) then
                dir = overPtr%dir(l)
                new_overlap = .true.
                do n = 1, count
                   if( tileList(n) == domain%list(pos)%tile_id(overPtr%tileNbr(l)) ) then ! same tile
                      dirdiff = dir - dirList(n)   ! dir always no less than dirList(n).
                      if( (dir == 8 .AND. dirList(n) == 1) .OR. (dir == 1 .AND. dirList(n) == 8) ) dirdiff = 1
                      dirdiff = abs(dirdiff)                      
                      if( dirdiff .LE. 1 ) then    ! same direction
                         if( mod(dirList(n),2) == 0 ) then
                            if(mod(dir,2) == 0 ) call mpp_error(FATAL, &
                                 "mpp_get_mosaic_refine_overlap: corner overlap should be unique.")
                            dirList(n) = dir
                         end if
                         if( dirList(n) == 1 .OR. dirList(n) == 5 ) then ! EAST, or WEST
                            if( isMe(n) .NE. overPtr%isMe(l) .OR. ieMe(n) .NE. overPtr%ieMe(l) ) call mpp_error(FATAL, &
                                 "mpp_get_mosaic_refine_overlap: mismatch on my index in the overlapping 1" )
                            jeMe(n) = max(jeMe(n),  overPtr%jeMe(l) )
                            jsMe(n) = min(jsMe(n),  overPtr%jsMe(l) )
                         else if(dirList(n) == 3 .OR. dirList(n) == 7 ) then ! SOUTH or NORTH
                            if( jsMe(n) .NE. overPtr%jsMe(l) .OR. jeMe(n) .NE. overPtr%jeMe(l) ) call mpp_error(FATAL, &
                                 "mpp_get_mosaic_refine_overlap: mismatch on my index in the overlapping 2" )
                            ieMe(n) = max(ieMe(n),  overPtr%ieMe(l) )
                            isMe(n) = min(isMe(n),  overPtr%isMe(l) )
                         end if
                         isNbr(n) = min(isNbr(n), overPtr%is(l))
                         ieNbr(n) = max(ieNbr(n), overPtr%ie(l))
                         jsNbr(n) = min(jsNbr(n), overPtr%js(l))
                         jeNbr(n) = max(jeNbr(n), overPtr%je(l))
                         OverPtr%index(l) = n
                         new_overlap = .false.
                         exit
                      end if
                   end if
                end do

                if(new_overlap) then
                   count = count + 1
                   if(count > MAXOVERLAP) call mpp_error(FATAL,  &
                        "mpp_get_mosaic_refine_overlap: number of overlapping is larger than MAXOVERLAP")
                   tileList(count) = domain%list(pos)%tile_id(overPtr%tileNbr(l))
                   dirList(count)  = dir
                   rotation(count) = overPtr%rotation(l)
                   OverPtr%index(l)= count
                   isMe(count) = overPtr%isMe(l);    ieMe(count) = overPtr%ieMe(l)
                   jsMe(count) = overPtr%jsMe(l);    jeMe(count) = overPtr%jeMe(l)
                   isNbr(count) = overPtr%is(l);     ieNbr(count) = overPtr%ie(l)
                   jsNbr(count) = overPtr%js(l);     jeNbr(count) = overPtr%je(l)
                end if
             end if
          end do
       end do
       if(count > 0) then
          allocate(update%rSpec(tMe)%isMe(count),  update%rSpec(tMe)%ieMe(count) )
          allocate(update%rSpec(tMe)%jsMe(count),  update%rSpec(tMe)%jeMe(count) )
          allocate(update%rSpec(tMe)%isNbr(count), update%rSpec(tMe)%ieNbr(count) )
          allocate(update%rSpec(tMe)%jsNbr(count), update%rSpec(tMe)%jeNbr(count) )
          allocate(update%rSpec(tMe)%start(count), update%rSpec(tMe)%end  (count) )
          allocate(update%rSpec(tMe)%dir  (count), update%rSpec(tMe)%rotation(count) )
          update%rSpec(tMe)%isMe = isMe(1:count); update%rSpec(tMe)%isNbr = isNbr(1:count)
          update%rSpec(tMe)%ieMe = ieMe(1:count); update%rSpec(tMe)%ieNbr = ieNbr(1:count)
          update%rSpec(tMe)%jsMe = jsMe(1:count); update%rSpec(tMe)%jsNbr = jsNbr(1:count)
          update%rSpec(tMe)%jeMe = jeMe(1:count); update%rSpec(tMe)%jeNbr = jeNbr(1:count)
          update%rSpec(tMe)%count = count
          update%rSpec(tMe)%dir   = dirList(1:count)
          update%rSpec(tMe)%rotation = rotation(1:count)
          update%rSpec(tMe)%total = sum( (update%rSpec(tMe)%ieNbr-update%rSpec(tMe)%isNbr+1)* &
               (update%rSpec(tMe)%jeNbr-update%rSpec(tMe)%jsNbr+1)  )
          iend = 0
          do n = 1, count
             istart = iend + 1
             iend = istart + (ieNbr(n) - isNbr(n) + 1)*(jeNbr(n) - jsNbr(n) + 1) - 1
             update%rSpec(tMe)%start(n) = istart
             update%rSpec(tMe)%end(n)   = iend
          enddo
       end if
  end do

end subroutine define_refine_overlap

!##############################################################################
!--- always fill the contact according to index order.
subroutine fill_contact(Contact, tile, is1, ie1, js1, je1, is2, ie2, js2, je2, align1, align2, refine1, refine2 )
  type(contact_type), intent(inout) :: Contact
  integer,            intent(in)    :: tile
  integer,            intent(in)    :: is1, ie1, js1, je1
  integer,            intent(in)    :: is2, ie2, js2, je2
  integer,            intent(in)    :: align1, align2
  real,               intent(in)    :: refine1, refine2 
  integer                           :: pos, n

  do pos = 1, Contact%ncontact
     select case(align1)
     case(WEST, EAST)  
        if( js1 < Contact%js1(pos) ) exit
     case(SOUTH, NORTH)  
        if( is1 < Contact%is1(pos) ) exit
     end select
  end do

  Contact%ncontact = Contact%ncontact + 1
  do n = Contact%ncontact, pos+1, -1  ! shift the data if needed.
     Contact%tile(n)   = Contact%tile(n-1)
     Contact%align1(n) = Contact%align1(n-1)
     Contact%align2(n) = Contact%align2(n-1)
     Contact%is1(n) = Contact%is1(n-1); Contact%ie1(n) = Contact%ie1(n-1)
     Contact%js1(n) = Contact%js1(n-1); Contact%je1(n) = Contact%je1(n-1)
     Contact%is2(n) = Contact%is2(n-1); Contact%ie2(n) = Contact%ie2(n-1)
     Contact%js2(n) = Contact%js2(n-1); Contact%je2(n) = Contact%je2(n-1)
  end do

  Contact%tile(pos)    = tile
  Contact%align1(pos)  = align1
  Contact%align2(pos)  = align2
  Contact%refine1(pos) = refine1
  Contact%refine2(pos) = refine2    
  Contact%is1(pos) = is1; Contact%ie1(pos) = ie1
  Contact%js1(pos) = js1; Contact%je1(pos) = je1
  Contact%is2(pos) = is2; Contact%ie2(pos) = ie2
  Contact%js2(pos) = js2; Contact%je2(pos) = je2

end subroutine fill_contact

!############################################################################
! this routine sets the overlapping between tiles for E,C,N-cell based on T-cell overlapping
subroutine set_contact_point(domain, position)
  type(domain2d), intent(inout) :: domain
  integer,           intent(in) :: position

  integer                         :: ishift, jshift, nlist, list, m, n
  integer                         :: ntileMe, tMe, dir, count, pos, nsend, nrecv
  integer                         :: isoff1, ieoff1, isoff2, ieoff2, jsoff1, jeoff1, jsoff2, jeoff2
  type(overlap_type),  pointer    :: ptrIn  => NULL()
  type(overlapSpec),   pointer    :: update_in  => NULL()
  type(overlapSpec),   pointer    :: update_out => NULL()
  type(overlap_type)              :: overlapList(0:size(domain%list(:))-1)
  type(overlap_type)              :: overlap

  call mpp_get_domain_shift(domain, ishift, jshift, position)
  update_in => domain%update_T
  select case(position)
  case (CORNER)
     update_out => domain%update_C
  case (EAST)
     update_out => domain%update_E
  case (NORTH)
     update_out => domain%update_N
  case default
     call mpp_error(FATAL, "mpp_domains_define.inc(set_contact_point): the position should be CORNER, EAST or NORTH")
  end select

  update_out%xbegin = update_in%xbegin; update_out%xend = update_in%xend + ishift
  update_out%ybegin = update_in%ybegin; update_out%yend = update_in%yend + jshift
  update_out%whalo  = update_in%whalo;  update_out%ehalo  = update_in%ehalo
  update_out%shalo  = update_in%shalo;  update_out%nhalo  = update_in%nhalo

  nlist = size(domain%list(:))
  ntileMe = size(domain%x(:))
  call allocate_update_overlap(overlap, MAXOVERLAP)
  do m = 0, nlist-1
     call init_overlap_type(overlapList(m))
  enddo

  !--- first copy the send information in update_out to send
  nsend = update_out%nsend
  do m = 1, nsend
     pos = update_out%send(m)%pe - mpp_root_pe()
     call add_update_overlap(overlapList(pos), update_out%send(m))  
     call deallocate_overlap_type(update_out%send(m))    
  enddo
  if(ASSOCIATED(update_out%send) )deallocate(update_out%send)

  !--- loop over the list of overlapping.
  nsend = update_in%nsend
  do m = 1, nsend
     ptrIn  => update_in%send(m)
     pos  = PtrIn%pe - mpp_root_pe()
     do n = 1, ptrIn%count
        dir = ptrIn%dir(n)
        ! only set overlapping between tiles for send ( ptrOut%overlap(1) is false )
        if(ptrIn%from_contact(n)) then  
           select case ( dir )               
           case ( 1 ) ! to_pe's eastern halo
              select case(ptrIn%rotation(n))
              case (ZERO)  !  W -> E
                 isoff1 = ishift; ieoff1 = ishift; jsoff1 = 0;      jeoff1 = jshift
              case (NINETY) ! S -> E
                 isoff1 = 0;      ieoff1 = jshift; jsoff1 = ishift; jeoff1 = ishift
              end select
           case ( 2 ) ! to_pe's south-eastearn halo
              select case(ptrIn%rotation(n))
              case (ZERO)  
                 isoff1 = ishift; ieoff1 = ishift; jsoff1 = 0;      jeoff1 = 0
              case (NINETY) 
                 isoff1 = jshift; ieoff1 = jshift; jsoff1 = ishift; jeoff1 = ishift
              case (MINUS_NINETY) 
                 isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = 0
              end select
           case ( 3 ) ! to_pe's southern halo
              select case(ptrIn%rotation(n))
              case (ZERO)  !  N -> S
                 isoff1 = 0;      ieoff1 = ishift; jsoff1 = 0;      jeoff1 = 0
              case (MiNUS_NINETY) ! E -> S
                 isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = ishift
              end select
           case ( 4 ) ! to_pe's south-westearn halo
              select case(ptrIn%rotation(n))
              case (ZERO)  
                 isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = 0
              case (NINETY) 
                 isoff1 = jshift; ieoff1 = jshift; jsoff1 = 0;      jeoff1 = 0
              case (MINUS_NINETY) 
                 isoff1 = 0;      ieoff1 = 0;      jsoff1 = ishift; jeoff1 = ishift
              end select
           case ( 5 ) ! to_pe's western halo
              select case(ptrIn%rotation(n))
              case (ZERO)  !  E -> W
                 isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = jshift
              case (NINETY) ! N -> W
                 isoff1 = 0;      ieoff1 = jshift; jsoff1 = 0;      jeoff1 = 0
              end select
           case ( 6 ) ! to_pe's north-westearn halo
              select case(ptrIn%rotation(n))
              case (ZERO)  
                 isoff1 = 0;      ieoff1 = 0;      jsoff1 = jshift; jeoff1 = jshift
              case (NINETY) 
                 isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = 0
              case (MINUS_NINETY) 
                 isoff1 = jshift; ieoff1 = jshift; jsoff1 = ishift; jeoff1 = ishift
              end select
           case ( 7 ) ! to_pe's northern halo
              select case(ptrIn%rotation(n))
              case (ZERO)  !  S -> N
                 isoff1 = 0;      ieoff1 = ishift; jsoff1 = jshift; jeoff1 = jshift
              case (MINUS_NINETY) ! W -> N
                 isoff1 = jshift; ieoff1 = jshift; jsoff1 = 0;      jeoff1 = ishift
              end select
           case ( 8 ) ! to_pe's north-eastearn halo
              select case(ptrIn%rotation(n))
              case (ZERO)  
                 isoff1 = ishift; ieoff1 = ishift; jsoff1 = jshift; jeoff1 = jshift
              case (NINETY) 
                 isoff1 = 0;      ieoff1 = 0;      jsoff1 = ishift; jeoff1 = ishift
              case (MINUS_NINETY) 
                 isoff1 = jshift; ieoff1 = jshift; jsoff1 = 0;      jeoff1 = 0
              end select
           end select
           call insert_mosaic_update_overlap(overlap, PtrIn%pe, PtrIn%tileMe(n), PtrIn%tileNbr(n), &
                                             Ptrin%is(n) + isoff1, Ptrin%ie(n) + ieoff1, Ptrin%js(n) + jsoff1, &
                                             Ptrin%je(n) + jeoff1, PtrIn%dir(n), PtrIn%rotation(n), PtrIn%from_contact(n), &
                                             PtrIn%is_refined(n) )
        end if
     end do ! do n = 1, prtIn%count
     if(overlap%count > 0) then
        call add_update_overlap(overlapList(pos), overlap)
        call init_overlap_type(overlap)
     endif
  end do  ! do list = 0, nlist-1

  nsend = 0
  do list = 0, nlist-1
     m = mod( domain%pos+list, nlist )
     if(overlapList(m)%count>0) nsend = nsend+1
  enddo

  update_out%nsend = nsend
  if(nsend>0) then
     allocate(update_out%send(nsend))
     pos = 0
     do list = 0, nlist-1
        m = mod( domain%pos+list, nlist )
        if(overlapList(m)%count>0) then
           pos = pos+1
           if(pos>nsend) call mpp_error(FATAL, &
                "mpp_domains_define.inc(set_contact_point): pos should be no larger than nsend")
           call add_update_overlap(update_out%send(pos), overlapList(m))
           call deallocate_overlap_type(overlapList(m))
        endif
     enddo
     if(pos .NE. nsend) call mpp_error(FATAL, &
          "mpp_domains_define.inc(set_contact_point): pos should equal to nsend")
  endif



  !--- first copy the recv information in update_out to recv
  nrecv = update_out%nrecv
  do m = 1, nrecv
     pos = update_out%recv(m)%pe - mpp_root_pe()
     call add_update_overlap(overlapList(pos), update_out%recv(m))  
     call deallocate_overlap_type(update_out%recv(m))     
  enddo
  if(ASSOCIATED(update_out%recv) )deallocate(update_out%recv)

  !--- loop over the list of overlapping.
  nrecv = update_in%nrecv
  do m=1,nrecv
     ptrIn  => update_in%recv(m)
     pos  = PtrIn%pe - mpp_root_pe()
     do n = 1, ptrIn%count
        dir = ptrIn%dir(n)
        ! only set overlapping between tiles for recv ( ptrOut%overlap(1) is false )
        if(ptrIn%from_contact(n)) then  
           select case ( dir )               
           case ( 1 ) ! E
              isoff1 = ishift; ieoff1 = ishift; jsoff1 = 0;      jeoff1 = jshift
           case ( 2 ) ! SE
              isoff1 = ishift; ieoff1 = ishift; jsoff1 = 0;      jeoff1 = 0
           case ( 3 ) ! S
              isoff1 = 0;      ieoff1 = ishift; jsoff1 = 0;      jeoff1 = 0
           case ( 4 ) ! SW
              isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = 0
           case ( 5 ) ! W 
              isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = jshift
           case ( 6 ) ! NW
              isoff1 = 0;       ieoff1 = 0;      jsoff1 = jshift; jeoff1 = jshift
           case ( 7 ) ! N
              isoff1 = 0;      ieoff1 = ishift; jsoff1 = jshift; jeoff1 = jshift
           case ( 8 ) ! NE
              isoff1 = ishift; ieoff1 = ishift; jsoff1 = jshift; jeoff1 = jshift
           end select
           call insert_mosaic_update_overlap(overlap, PtrIn%pe, PtrIn%tileMe(n), PtrIn%tileNbr(n), &
                                             Ptrin%is(n) + isoff1, Ptrin%ie(n) + ieoff1, Ptrin%js(n) + jsoff1, &
                                             Ptrin%je(n) + jeoff1, PtrIn%dir(n), PtrIn%rotation(n), PtrIn%from_contact(n), &
                                             PtrIn%is_refined(n) )
           count                     = overlap%count
           overlap%isMe      (count) = Ptrin%isMe(n) + isoff1
           overlap%ieMe      (count) = Ptrin%ieMe(n) + ieoff1
           overlap%jsMe      (count) = Ptrin%jsMe(n) + jsoff1
           overlap%jeMe      (count) = Ptrin%jeMe(n) + jeoff1
           !--- figure out neighbor offset when there is refinemnet
           if( PtrIn%is_refined(n) ) then
              isoff2 = isoff1; ieoff2 = ieoff1; jsoff2 = jsoff1; jeoff2 = jeoff1
              if( ptrIn%rotation(n) == NINETY) then
                 select case ( dir )                                       
                 case ( 1, 5 ) !S -> E, N -> W
                    isoff2 = jsoff1; ieoff2 = jeoff1; jsoff2 = isoff1; jeoff2 = ieoff1
                 case ( 2 ) 
                    isoff2 = jshift; ieoff2 = jshift; jsoff2 = ishift; jeoff2 = ishift
                 case ( 4 ) 
                    isoff2 = jshift; ieoff2 = jshift; jsoff2 = 0;      jeoff2 = 0
                 case ( 6 ) 
                    isoff2 = 0;      ieoff2 = 0;      jsoff2 = 0;      jeoff2 = 0
                 case ( 8 ) 
                    isoff2 = 0;      ieoff2 = 0;      jsoff2 = ishift; jeoff2 = ishift
                 end select
              else if( ptrIn%rotation(n) == MINUS_NINETY) then
                 select case ( dir )                                       
                 case ( 3, 7 ) !E -> S, W -> N
                    isoff2 = jsoff1; ieoff2 = jeoff1; jsoff2 = isoff1; jeoff2 = ieoff1
                 case ( 2 )
                    isoff2 = 0;      ieoff2 = 0;      jsoff2 = 0;      jeoff2 = 0 
                 case ( 4 ) 
                    isoff2 = 0;      ieoff2 = 0;      jsoff2 = ishift; jeoff2 = ishift
                 case ( 6 ) 
                    isoff2 = jshift; ieoff2 = jshift; jsoff2 = ishift; jeoff2 = ishift
                 case ( 8 ) 
                    isoff2 = jshift; ieoff2 = jshift; jsoff2 = 0;      jeoff2 = 0
                 end select
              end if
              !--- calculate overlapping
              overlap%is(count) = Ptrin%is(n) + isoff2
              overlap%ie(count) = Ptrin%ie(n) + ieoff2
              overlap%js(count) = Ptrin%js(n) + jsoff2
              overlap%je(count) = Ptrin%je(n) + jeoff2
           end if
        end if
     end do ! do n = 1, ptrIn%count
     if(overlap%count > 0) then
        call add_update_overlap(overlapList(pos), overlap)
        call init_overlap_type(overlap)
     endif
     do tMe = 1, size(domain%x(:))
        do n = 1, overlap%count
           if(overlap%tileMe(n) == tMe) then
              if(overlap%dir(n) == 1 ) domain%x(tMe)%loffset = 0  
              if(overlap%dir(n) == 7 ) domain%y(tMe)%loffset = 0  
           end if
        end do
     end do
  end do ! do list = 0, nlist-1

  nrecv = 0
  do list = 0, nlist-1
     m = mod( domain%pos+nlist-list, nlist )
     if(overlapList(m)%count>0) nrecv = nrecv+1
  enddo

  update_out%nrecv = nrecv
  if(nrecv>0) then
     allocate(update_out%recv(nrecv))
     pos = 0
     do list = 0, nlist-1
        m = mod( domain%pos+nlist-list, nlist )
        if(overlapList(m)%count>0) then
           pos = pos+1
           if(pos>nrecv) call mpp_error(FATAL, &
                "mpp_domains_define.inc(set_contact_point): pos should be no larger than nrecv")
           call  add_update_overlap(update_out%recv(pos), overlapList(m))
           call deallocate_overlap_type(overlapList(m))
        endif
     enddo
     if(pos .NE. nrecv) call mpp_error(FATAL, &
          "mpp_domains_define.inc(set_contact_point): pos should equal to nrecv")
  endif



  !--- define the refined overlapping and group the overlapping 
  !--- between same tiles( in different direction or different pe ). 
  !--- set up the index position in the return buffer.
  call define_refine_overlap( domain, position )  

  call deallocate_overlap_type(overlap)

end subroutine set_contact_point

!--- set up the overlapping for boundary check if the domain is symmetry. The check will be 
!--- done on current pe for east boundary for E-cell, north boundary for N-cell, 
!--- East and North boundary for C-cell
subroutine set_check_overlap( domain, position )
type(domain2d),    intent(in) :: domain
integer,           intent(in) :: position
integer                       :: nlist, m, n
integer, parameter            :: MAXCOUNT = 100
integer                       :: is, ie, js, je
integer                       :: nsend, nrecv, pos, maxsize, rotation
type(overlap_type)            :: overlap
type(overlapSpec),   pointer  :: update  => NULL()
type(overlapSpec),   pointer  :: check   => NULL()

select case(position)
case (CORNER)
   update => domain%update_C
   check  => domain%check_C
case (EAST)
   update => domain%update_E
   check  => domain%check_E
case (NORTH)
   update => domain%update_N
   check  => domain%check_N
case default
   call mpp_error(FATAL, "mpp_domains_define.inc(set_check_overlap): position should be CORNER, EAST or NORTH")
end select

check%xbegin = update%xbegin; check%xend = update%xend
check%ybegin = update%ybegin; check%yend = update%yend
check%nsend  = 0
check%nrecv  = 0
if( .NOT. domain%symmetry ) return

nsend = 0
maxsize = 0
do m = 1, update%nsend
    do n = 1, update%send(m)%count
       if( update%send(m)%is_refined(n) ) cycle 
       if( update%send(m)%rotation(n) == ONE_HUNDRED_EIGHTY ) cycle
       if( ( (position == EAST .OR. position == CORNER) .AND. update%send(m)%dir(n) == 1 ) .OR. &
            ( (position == NORTH .OR. position == CORNER) .AND. update%send(m)%dir(n) == 7 ) ) then
            maxsize = max(maxsize, update%send(m)%count)
            nsend = nsend + 1
            exit
       endif      
    enddo
enddo

if(nsend>0) then
   allocate(check%send(nsend))
   call allocate_check_overlap(overlap, maxsize)
endif


nlist = size(domain%list(:))
!--- loop over the list of domains to find the boundary overlap for send
pos = 0
do m = 1, update%nsend
 do n = 1, update%send(m)%count
    if( update%send(m)%is_refined(n) ) cycle
    if( update%send(m)%rotation(n) == ONE_HUNDRED_EIGHTY ) cycle
    ! comparing east direction on currently pe
    if( (position == EAST .OR. position == CORNER) .AND. update%send(m)%dir(n) == 1 ) then  
       rotation = update%send(m)%rotation(n)
       select case( rotation ) 
       case( ZERO ) ! W -> E
          is = update%send(m)%is(n) - 1
          ie = is
          js = update%send(m)%js(n)
          je = update%send(m)%je(n)
       case( NINETY ) ! S -> E
          is = update%send(m)%is(n)
          ie = update%send(m)%ie(n)
          js = update%send(m)%js(n) - 1
          je = js
       end select
       call insert_check_overlap(overlap, update%send(m)%pe, &
                                 update%send(m)%tileMe(n), 1, rotation, is, ie, js, je)
    end if

    ! comparing north direction on currently pe
    if( (position == NORTH .OR. position == CORNER) .AND. update%send(m)%dir(n) == 7 ) then  
       rotation = update%send(m)%rotation(n)
       select case( rotation ) 
       case( ZERO ) ! S->N
          is = update%send(m)%is(n)
          ie = update%send(m)%ie(n)
          js = update%send(m)%js(n) - 1
          je = js
       case( MINUS_NINETY ) ! W->N
          is = update%send(m)%is(n) - 1
          ie = is
          js = update%send(m)%js(n)
          je = update%send(m)%je(n)
       end select
       call insert_check_overlap(overlap, update%send(m)%pe, &
            update%send(m)%tileMe(n), 4, rotation, is, ie, js, je)
    end if
 end do ! do n =1, update%send(m)%count
 if(overlap%count>0) then
   pos = pos+1
   if(pos>nsend)call mpp_error(FATAL, "mpp_domains_define.inc(set_check_overlap): pos is greater than nsend")
   call add_check_overlap(check%send(pos), overlap)
   call init_overlap_type(overlap)
 endif
end do ! end  do list = 0, nlist

if(pos .NE. nsend)call mpp_error(FATAL, "mpp_domains_define.inc(set_check_overlap): pos is greater than nsend")

nrecv = 0
maxsize = 0
do m = 1, update%nrecv
    do n = 1, update%recv(m)%count
       if( update%recv(m)%is_refined(n) ) cycle 
       if( update%recv(m)%rotation(n) == ONE_HUNDRED_EIGHTY ) cycle
       if( ( (position == EAST .OR. position == CORNER) .AND. update%recv(m)%dir(n) == 1 ) .OR. &
            ( (position == NORTH .OR. position == CORNER) .AND. update%recv(m)%dir(n) == 7 ) ) then
            maxsize = max(maxsize, update%recv(m)%count)
            nrecv = nrecv + 1
            exit
       endif      
    enddo
enddo

if(nsend>0) call deallocate_overlap_type(overlap)

if(nrecv>0) then
   allocate(check%recv(nrecv))
   call allocate_check_overlap(overlap, maxsize)
endif

pos = 0
do m = 1, update%nrecv
 do n = 1, update%recv(m)%count
    if( update%recv(m)%is_refined(n) ) cycle
    if( update%recv(m)%rotation(n) == ONE_HUNDRED_EIGHTY ) cycle
    if( (position == EAST .OR. position == CORNER) .AND. update%recv(m)%dir(n) == 1 ) then  
       is = update%recv(m)%is(n) - 1
       ie = is
       js = update%recv(m)%js(n)
       je = update%recv(m)%je(n)
       call insert_check_overlap(overlap, update%recv(m)%pe, &
                                 update%recv(m)%tileMe(n), 1, update%recv(m)%rotation(n), is, ie, js, je)
    end if
    if( (position == NORTH .OR. position == CORNER) .AND. update%recv(m)%dir(n) == 7 ) then  
       is = update%recv(m)%is(n)
       ie = update%recv(m)%ie(n)
       js = update%recv(m)%js(n) - 1
       je = js
       call insert_check_overlap(overlap, update%recv(m)%pe, &
            update%recv(m)%tileMe(n), 3, update%recv(m)%rotation(n), is, ie, js, je)
    end if
 end do ! n = 1, overlap%count
 if(overlap%count>0) then
   pos = pos+1
   if(pos>nrecv)call mpp_error(FATAL, "mpp_domains_define.inc(set_check_overlap): pos is greater than nrecv")    
   call add_check_overlap(check%recv(pos), overlap)
   call init_overlap_type(overlap)
 endif
end do ! end  do list = 0, nlist

if(pos .NE. nrecv)call mpp_error(FATAL, "mpp_domains_define.inc(set_check_overlap): pos is greater than nrecv")
if(nrecv>0) call deallocate_overlap_type(overlap)

end subroutine set_check_overlap

!#############################################################################
!--- set up the overlapping for boundary if the domain is symmetry.
subroutine set_bound_overlap( domain, position )
  type(domain2d), intent(inout) :: domain
  integer,        intent(in)    :: position
  integer                       :: m, n, l, count, dr, tMe
  integer, parameter            :: MAXCOUNT = 100
  integer, dimension(MAXCOUNT)  :: dir, rotation, is, ie, js, je, isMe, ieMe, jsMe, jeMe, tileMe
  integer, dimension(size(domain%x(:)), 4)           :: nrecvl
  integer, dimension(size(domain%x(:)), 4, MAXCOUNT) :: isl, iel, jsl, jel, islMe, ielMe, jslMe, jelMe
  type(overlap_type),  pointer   :: overlap => NULL()
  type(overlapSpec),   pointer   :: update => NULL()
  type(overlapSpec),   pointer   :: bound => NULL()
  integer                        :: nlist_send, nlist_recv, ishift, jshift
  integer                        :: ism, iem, jsm, jem, nsend, nrecv

  if( position == CENTER .OR. .NOT. domain%symmetry ) return
  call mpp_get_domain_shift(domain, ishift, jshift, position)
  call mpp_get_memory_domain ( domain, ism, iem, jsm, jem )

  select case(position)
  case (CORNER)
     update => domain%update_C
     bound  => domain%bound_C
  case (EAST)
     update => domain%update_E
     bound  => domain%bound_E
  case (NORTH)
     update => domain%update_N
     bound  => domain%bound_N
  case default
     call mpp_error( FATAL, "mpp_domains_mod(set_bound_overlap): invalid option of position")
  end select

  bound%xbegin = ism; bound%xend = iem + ishift
  bound%ybegin = jsm; bound%yend = jem + jshift

  nlist_send = update%nsend
  nlist_recv = update%nrecv
  bound%nsend = nlist_send
  bound%nrecv = nlist_recv
  if(nlist_send >0) then
     allocate(bound%send(nlist_send))
     bound%send(:)%count = 0
  endif
  if(nlist_recv >0) then
     allocate(bound%recv(nlist_recv))
     bound%recv(:)%count = 0
  endif
  !--- loop over the list of domains to find the boundary overlap for send

  nsend = 0
  do m = 1, nlist_send
     overlap => update%send(m)
     if( overlap%count == 0 ) cycle
     count = 0
     do n = 1, overlap%count
        if( (position == EAST .OR. position == CORNER) .AND. overlap%dir(n) == 1) then ! east
           count=count+1
           dir(count) = 1
           rotation(count) = overlap%rotation(n)
           tileMe(count)   = overlap%tileMe(n)
           select case( rotation(count) ) 
           case( ZERO ) ! W -> E
              is(count) = overlap%is(n) - 1
              ie(count) = is(count)
              js(count) = overlap%js(n)
              je(count) = overlap%je(n)
           case( NINETY ) ! S -> E
              is(count) = overlap%is(n)
              ie(count) = overlap%ie(n)
              js(count) = overlap%js(n) - 1
              je(count) = js(count)
           end select
        end if
        if( (position == NORTH .OR. position == CORNER) .AND. overlap%dir(n) == 3 ) then ! south
           count=count+1
           dir(count) = 2
           rotation(count) = overlap%rotation(n)
           tileMe(count)   = overlap%tileMe(n)
           select case( rotation(count) ) 
           case( ZERO ) ! N->S
              is(count) = overlap%is(n)
              ie(count) = overlap%ie(n)
              js(count) = overlap%je(n) + 1
              je(count) = js(count)
           case( MINUS_NINETY ) ! E->S
              is(count) = overlap%ie(n) + 1
              ie(count) = is(count)
              js(count) = overlap%js(n)
              je(count) = overlap%je(n)
           end select
        end if
        if( (position == EAST .OR. position == CORNER) .AND. overlap%dir(n) == 5 ) then ! west
           count=count+1
           dir(count) = 3
           rotation(count) = overlap%rotation(n)
           tileMe(count) = overlap%tileMe(n)
           select case( rotation(count) ) 
           case( ZERO ) ! E->W
              is(count) = overlap%ie(n) + 1
              ie(count) = is(count)
              js(count) = overlap%js(n)
              je(count) = overlap%je(n)
           case( NINETY ) ! N->W
              is(count) = overlap%is(n)
              ie(count) = overlap%ie(n)
              js(count) = overlap%je(n) + 1
              je(count) = js(count)
           end select
        end if
        if( (position == NORTH .OR. position == CORNER) .AND. overlap%dir(n) == 7 ) then ! north
           count=count+1
           dir(count) = 4
           rotation(count) = overlap%rotation(n)
           tileMe(count) = overlap%tileMe(n)
           select case( rotation(count) ) 
           case( ZERO ) ! S->N
              is(count) = overlap%is(n)
              ie(count) = overlap%ie(n)
              js(count) = overlap%js(n) - 1
              je(count) = js(count)
           case( MINUS_NINETY ) ! W->N
              is(count) = overlap%is(n) - 1
              ie(count) = is(count)
              js(count) = overlap%js(n)
              je(count) = overlap%je(n)
           end select
        end if
     end do ! do n =1, overlap%count
     if(count>0) then
        nsend = nsend + 1
        bound%send(nsend)%count = count
        bound%send(nsend)%pe    = overlap%pe
        allocate(bound%send(nsend)%is(count),  bound%send(nsend)%ie(count) )
        allocate(bound%send(nsend)%js(count),  bound%send(nsend)%je(count) )
        allocate(bound%send(nsend)%dir(count), bound%send(nsend)%rotation(count) )
        allocate(bound%send(nsend)%tileMe(count))
        bound%send(nsend)%is(:)       = is(1:count)
        bound%send(nsend)%ie(:)       = ie(1:count)
        bound%send(nsend)%js(:)       = js(1:count)
        bound%send(nsend)%je(:)       = je(1:count)
        bound%send(nsend)%dir(:)      = dir(1:count)
        bound%send(nsend)%tileMe(:)   = tileMe(1:count)
        bound%send(nsend)%rotation(:) = rotation(1:count)
     end if
  end do ! end  do list = 0, nlist

  !--- loop over the list of domains to find the boundary overlap for recv
  bound%nsend = nsend
  nrecvl(:,:) = 0    
  nrecv       = 0
  
  do m = 1, nlist_recv
     overlap => update%recv(m)
     if( overlap%count == 0 ) cycle
     count = 0
     do n = 1, overlap%count
        if( (position == EAST .OR. position == CORNER) .AND. overlap%dir(n) == 1) then ! east
           count=count+1
           dir(count) = 1
           rotation(count) = overlap%rotation(n)
           isMe(count)     = overlap%isMe(n) - 1
           ieMe(count) = isMe(count)
           jsMe(count)     = overlap%jsMe(n)
           jeMe(count)     = overlap%jeMe(n)
           tileMe(count)   = overlap%tileMe(n)
           if( overlap%is_refined(n)) then
              select case( rotation(count) ) 
              case( ZERO ) ! W -> E
                 is(count) = overlap%is(n) - 1
                 ie(count) = is(count)
                 js(count) = overlap%js(n)
                 je(count) = overlap%je(n)
              case( NINETY ) ! S -> E
                 is(count) = overlap%is(n)
                 ie(count) = overlap%ie(n)
                 js(count) = overlap%js(n) - 1
                 je(count) = js(count)
              end select
           else
              is(count) = overlap%is(n) - 1
              ie(count) = is(count)
              js(count) = overlap%js(n)
              je(count) = overlap%je(n)
           end if
           tMe                        = tileMe(count)
           nrecvl(tMe, 1) = nrecvl(tMe,1) + 1
           islMe(tMe,1,nrecvl(tMe, 1)) = isMe(count)
           ielMe(tMe,1,nrecvl(tMe, 1)) = ieMe(count)
           jslMe(tMe,1,nrecvl(tMe, 1)) = jsMe(count)
           jelMe(tMe,1,nrecvl(tMe, 1)) = jeMe(count)
           isl  (tMe,1,nrecvl(tMe, 1)) = is  (count)
           iel  (tMe,1,nrecvl(tMe, 1)) = ie  (count)
           jsl  (tMe,1,nrecvl(tMe, 1)) = js  (count)
           jel  (tMe,1,nrecvl(tMe, 1)) = je  (count)
        end if

        if( (position == NORTH .OR. position == CORNER) .AND. overlap%dir(n) == 3) then ! south
           count=count+1
           dir(count) = 2
           rotation(count) = overlap%rotation(n)
           isMe(count)     = overlap%isMe(n)
           ieMe(count)     = overlap%ieMe(n)
           jsMe(count)     = overlap%jeMe(n) + 1
           jeMe(count) = jsMe(count)
           tileMe(count)   = overlap%tileMe(n)
           if( overlap%is_refined(n)) then
              select case( rotation(count) ) 
              case( ZERO ) ! N->S
                 is(count) = overlap%is(n)
                 ie(count) = overlap%ie(n)
                 js(count) = overlap%je(n) + 1
                 je(count) = js(count)
              case( MINUS_NINETY ) ! E->S
                 is(count) = overlap%ie(n) + 1
                 ie(count) = is(count)
                 js(count) = overlap%js(n)
                 je(count) = overlap%je(n)
              end select
           else
              is(count) = overlap%is(n)
              ie(count) = overlap%ie(n)
              js(count) = overlap%je(n) + 1
              je(count) = js(count)
           end if
           tMe                        = tileMe(count)
           nrecvl(tMe, 2) = nrecvl(tMe,2) + 1
           islMe(tMe,2,nrecvl(tMe, 2)) = isMe(count)
           ielMe(tMe,2,nrecvl(tMe, 2)) = ieMe(count)
           jslMe(tMe,2,nrecvl(tMe, 2)) = jsMe(count)
           jelMe(tMe,2,nrecvl(tMe, 2)) = jeMe(count)
           isl  (tMe,2,nrecvl(tMe, 2)) = is  (count)
           iel  (tMe,2,nrecvl(tMe, 2)) = ie  (count)
           jsl  (tMe,2,nrecvl(tMe, 2)) = js  (count)
           jel  (tMe,2,nrecvl(tMe, 2)) = je  (count)
        end if

        if( (position == EAST .OR. position == CORNER) .AND. overlap%dir(n) == 5) then ! west
           count=count+1
           dir(count) = 3
           rotation(count) = overlap%rotation(n)
           isMe(count)     = overlap%ieMe(n) + 1 
           ieMe(count) = isMe(count)
           jsMe(count)     = overlap%jsMe(n)
           jeMe(count)     = overlap%jeMe(n)
           tileMe(count)   = overlap%tileMe(n)
           if( overlap%is_refined(n)) then
              select case( rotation(count) ) 
              case( ZERO ) ! E->W
                 is(count) = overlap%ie(n) + 1
                 ie(count) = is(count)
                 js(count) = overlap%js(n)
                 je(count) = overlap%je(n)
              case( NINETY ) ! S -> E
                 is(count) = overlap%is(n)
                 ie(count) = overlap%ie(n)
                 js(count) = overlap%je(n) + 1
                 je(count) = js(count)
              end select
           else
              is(count) = overlap%ie(n) + 1
              ie(count) = is(count)
              js(count) = overlap%js(n)
              je(count) = overlap%je(n)
           end if
           tMe                        = tileMe(count)
           nrecvl(tMe, 3) = nrecvl(tMe,3) + 1
           islMe(tMe,3,nrecvl(tMe, 3)) = isMe(count)
           ielMe(tMe,3,nrecvl(tMe, 3)) = ieMe(count)
           jslMe(tMe,3,nrecvl(tMe, 3)) = jsMe(count)
           jelMe(tMe,3,nrecvl(tMe, 3)) = jeMe(count)
           isl  (tMe,3,nrecvl(tMe, 3)) = is  (count)
           iel  (tMe,3,nrecvl(tMe, 3)) = ie  (count)
           jsl  (tMe,3,nrecvl(tMe, 3)) = js  (count)
           jel  (tMe,3,nrecvl(tMe, 3)) = je  (count)
        end if

        if( (position == NORTH .OR. position == CORNER) .AND. overlap%dir(n) == 7) then ! north
           count=count+1
           dir(count) = 4
           rotation(count) = overlap%rotation(n)
           isMe(count)     = overlap%isMe(n)
           ieMe(count)     = overlap%ieMe(n)
           jsMe(count)     = overlap%jsMe(n) - 1
           jeMe(count) = jsMe(count)
           tileMe(count) = overlap%tileMe(n)
           if( overlap%is_refined(n)) then
              select case( rotation(count) ) 
              case( ZERO ) ! S->N
                 is(count) = overlap%is(n)
                 ie(count) = overlap%ie(n)
                 js(count) = overlap%js(n) - 1
                 je(count) = js(count)
              case( MINUS_NINETY ) ! W->N
                 is(count) = overlap%is(n) - 1
                 ie(count) = is(count)
                 js(count) = overlap%js(n)
                 je(count) = overlap%je(n)
              end select
           else
              is(count) = overlap%is(n)
              ie(count) = overlap%ie(n)
              js(count) = overlap%js(n) - 1
              je(count) = js(count)
           end if
           tMe                        = tileMe(count)
           nrecvl(tMe, 4) = nrecvl(tMe,4) + 1
           islMe(tMe,4,nrecvl(tMe, 4)) = isMe(count)
           ielMe(tMe,4,nrecvl(tMe, 4)) = ieMe(count)
           jslMe(tMe,4,nrecvl(tMe, 4)) = jsMe(count)
           jelMe(tMe,4,nrecvl(tMe, 4)) = jeMe(count)
           isl  (tMe,4,nrecvl(tMe, 4)) = is  (count)
           iel  (tMe,4,nrecvl(tMe, 4)) = ie  (count)
           jsl  (tMe,4,nrecvl(tMe, 4)) = js  (count)
           jel  (tMe,4,nrecvl(tMe, 4)) = je  (count)
        end if
     end do ! do n = 1, overlap%count
     if(count>0) then
        nrecv = nrecv + 1
        bound%recv(nrecv)%count = count
        bound%recv(nrecv)%pe    = overlap%pe
        allocate(bound%recv(nrecv)%isMe(count),   bound%recv(nrecv)%ieMe(count) )
        allocate(bound%recv(nrecv)%jsMe(count),   bound%recv(nrecv)%jeMe(count) )
        allocate(bound%recv(nrecv)%is(count),     bound%recv(nrecv)%ie(count) )
        allocate(bound%recv(nrecv)%js(count),     bound%recv(nrecv)%je(count) )
        allocate(bound%recv(nrecv)%dir(count),    bound%recv(nrecv)%index(count)  )
        allocate(bound%recv(nrecv)%tileMe(count), bound%recv(nrecv)%rotation(count) )
        bound%recv(nrecv)%isMe(:)     = isMe(1:count)
        bound%recv(nrecv)%ieMe(:)     = ieMe(1:count)
        bound%recv(nrecv)%jsMe(:)     = jsMe(1:count)
        bound%recv(nrecv)%jeMe(:)     = jeMe(1:count)
        bound%recv(nrecv)%is(:)       = is(1:count)
        bound%recv(nrecv)%ie(:)       = ie(1:count)
        bound%recv(nrecv)%js(:)       = js(1:count)
        bound%recv(nrecv)%je(:)       = je(1:count)
        bound%recv(nrecv)%dir(:)      = dir(1:count)
        bound%recv(nrecv)%tileMe(:)   = tileMe(1:count)
        bound%recv(nrecv)%rotation(:) = rotation(1:count)
     end if
  end do ! end  do list = 0, nlist

  bound%nrecv = nrecv

  !--- find the boundary index for each contact within the east boundary
  do m = 1, nrecv
     do n = 1, bound%recv(m)%count
        tMe = bound%recv(m)%tileMe(n)
        dr = bound%recv(m)%dir(n)
        bound%recv(m)%index(n) = 1
        do l = 1, nrecvl(tMe,dr)
           if(dr == 1 .OR. dr == 3) then  ! EAST, WEST
              if( bound%recv(m)%jsMe(n) > jslMe(tMe, dr, l) ) then
                 bound%recv(m)%index(n) = bound%recv(m)%index(n)                  + &
                      max(abs(jel(tMe, dr, l)-jsl(tMe, dr, l)), &
                      abs(iel(tMe, dr, l)-isl(tMe, dr, l))) + 1 - jshift
              end if
           else                             ! South, North
              if( bound%recv(m)%isMe(n) > islMe(tMe, dr, l) ) then
                 bound%recv(m)%index(n) = bound%recv(m)%index(n)                  + &
                      max(abs(jel(tMe, dr, l)-jsl(tMe, dr, l)), &
                      abs(iel(tMe, dr, l)-isl(tMe, dr, l))) + 1 - ishift
              end if
           end if
        end do
     end do
  end do


end subroutine set_bound_overlap


!#############################################################################

subroutine fill_corner_contact(eCont, sCont, wCont, nCont, isg, ieg, jsg, jeg, numR, numS, tileRecv, tileSend, &
 is1Recv, ie1Recv, js1Recv, je1Recv, is2Recv, ie2Recv, js2Recv, je2Recv,         &
 is1Send, ie1Send, js1Send, je1Send, is2Send, ie2Send, js2Send, je2Send,         &
 align1Recv, align2Recv, align1Send, align2Send,                                 &
 whalo, ehalo, shalo, nhalo, tileMe)
type(contact_type), dimension(:), intent(in) :: eCont, sCont, wCont, nCont
integer, dimension(:),            intent(in) :: isg, ieg, jsg, jeg
integer,                       intent(inout) :: numR, numS
integer, dimension(:),         intent(inout) :: tileRecv, tileSend
integer, dimension(:),         intent(inout) :: is1Recv, ie1Recv, js1Recv, je1Recv
integer, dimension(:),         intent(inout) :: is2Recv, ie2Recv, js2Recv, je2Recv
integer, dimension(:),         intent(inout) :: is1Send, ie1Send, js1Send, je1Send
integer, dimension(:),         intent(inout) :: is2Send, ie2Send, js2Send, je2Send
integer, dimension(:),         intent(inout) :: align1Recv, align2Recv, align1Send, align2Send
integer,                          intent(in) :: tileMe, whalo, ehalo, shalo, nhalo
integer                                      :: is1, ie1, js1, je1, is2, ie2, js2, je2
integer                                      :: tn, tc, n, m
logical                                      :: found_corner

found_corner = .false.
!--- southeast for recving
if(eCont(tileMe)%ncontact > 0) then     
 if(eCont(tileMe)%js1(1) == jsg(tileMe) ) then
    tn = eCont(tileMe)%tile(1)
    if(econt(tileMe)%js2(1) > jsg(tn) ) then  ! the corner tile is tn.
       if( econt(tileMe)%js2(1) - jsg(tn) < shalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: southeast tile for recv 1 is not tiled properly")
       found_corner = .true.; tc = tn
       is1 = eCont(tileMe)%ie1(1) + 1; je1 = eCont(tileMe)%js1(1) - 1
       is2 = eCont(tileMe)%is2(1);     je2 = eCont(tileMe)%js2(1) - 1
    else if(sCont(tn)%ncontact >0) then ! the corner tile may be south tile of tn.
       if(sCont(tn)%is1(1) == isg(tn)) then ! corner is nc.
          found_corner = .true.; tc = sCont(tn)%tile(1)
          is1 = eCont(tileMe)%ie1(1) + 1; je1 = eCont(tileMe)%js1(1) - 1
          is2 = sCont(tn)%is2(1);         je2 = sCont(tn)%je2(1)
       end if
    end if
 end if
end if
if( .not. found_corner ) then  ! not found,
 n = sCont(tileMe)%ncontact
 if( n > 0) then
    if( sCont(tileMe)%ie1(n) == ieg(tileMe)) then
       tn = sCont(tileMe)%tile(n)
       if(scont(tileMe)%ie2(n) < ieg(tn) ) then  ! the corner tile is tn.
          if(ieg(tn) - scont(tileMe)%ie2(n) < ehalo ) call mpp_error(FATAL, &
               "mpp_domains_define.inc: southeast tile for recv 2 is not tiled properly")
          found_corner = .true.; tc = tn
          is1 = sCont(tileMe)%ie1(n) + 1; je1 = sCont(tileMe)%js1(n) - 1
          is2 = sCont(tileMe)%ie2(n) + 1; je2 = sCont(tileMe)%je2(n)
       else if(eCont(tn)%ncontact >0) then ! the corner tile may be east tile of tn.
          m = eCont(tn)%ncontact
          if(eCont(tn)%je1(m) == jeg(tn)) then ! corner is nc.
             found_corner = .true.; tc = eCont(tn)%tile(m)
             is1 = sCont(tileMe)%ie1(n) + 1; je1 = sCont(tileMe)%js1(n) - 1
             is2 = eCont(tn)%is2(m);         je2 = eCont(tn)%je2(m)
          end if
       end if
    end if
 end if
end if
if(found_corner) then
 numR = numR + 1
 tileRecv(numR) = tc; align1Recv(numR) = SOUTH_EAST;  align2Recv(numR) = NORTH_WEST
 is1Recv(numR) = is1;             ie1Recv(numR) = is1 + ehalo - 1
 js1Recv(numR) = je1 - shalo + 1; je1Recv(numR) = je1
 is2Recv(numR) = is2;             ie2Recv(numR) = is2 + ehalo - 1
 js2Recv(numR) = je2 - shalo + 1; je2Recv(numR) = je2  
end if

!--- southwest for recving
found_corner = .false.
if(wCont(tileMe)%ncontact > 0) then
 if(wCont(tileMe)%js1(1) == jsg(tileMe) ) then
    tn = wCont(tileMe)%tile(1)
    if(wcont(tileMe)%js2(1) > jsg(tn) ) then  ! the corner tile is tn.
       if( wcont(tileMe)%js2(1) - jsg(tn) < shalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: southwest tile for recv 1 is not tiled properly")
       found_corner = .true.; tc = tn
       ie1 = wCont(tileMe)%is1(1) - 1; je1 = wCont(tileMe)%js1(1) - 1
       ie2 = wCont(tileMe)%is2(1);     je2 = wCont(tileMe)%js2(1) - 1
    else if(sCont(tn)%ncontact >0) then ! the corner tile may be south tile of tn.
       n = sCont(tn)%ncontact
       if(sCont(tn)%ie1(n) == ieg(tn)) then ! corner is nc.
          found_corner = .true.;  tc = sCont(tn)%tile(n)
          ie1 = wCont(tileMe)%is1(1) - 1; je1 = wCont(tileMe)%js1(1) - 1
          ie2 = sCont(tn)%ie2(1);         je2 = sCont(tn)%je2(1)
       end if
    end if
 end if
end if
if( .not. found_corner ) then  ! not found,
 n = sCont(tileMe)%ncontact
 if( n > 0) then
    if( sCont(tileMe)%is1(1) == isg(tileMe)) then
       tn = sCont(tileMe)%tile(1)
       if(sCont(tileMe)%is2(1) > isg(tn) ) then  ! the corner tile is tn.
          if( scont(tileMe)%is2(1)-isg(tn) < whalo ) call mpp_error(FATAL, &
               "mpp_domains_define.inc: southwest tile for recv 1 is not tiled properly")
          found_corner = .true.; tc = tn
          ie1 = sCont(tileMe)%is1(1) - 1; je1 = sCont(tileMe)%js1(1) - 1
          ie2 = sCont(tileMe)%is2(1) - 1; je2 = sCont(tileMe)%js2(1)
       else if(wCont(tn)%ncontact >0) then ! the corner tile may be west tile of tn.
          m = wCont(tn)%ncontact
          if(wCont(tn)%je1(m) == jeg(tn)) then ! corner is nc.
             found_corner = .true.; tc = wCont(tn)%tile(m)
             ie1 = sCont(tileMe)%is1(1) - 1; je1 = sCont(tileMe)%js1(1) - 1
             ie2 = wCont(tn)%ie2(m);         je2 = wCont(tn)%je2(m)
          end if
       end if
    end if
 end if
end if
if(found_corner) then
 numR = numR + 1
 tileRecv(numR) = tc; align1Recv(numR) = SOUTH_WEST; align2Recv(numR) = NORTH_EAST
 is1Recv(numR) = ie1 - whalo + 1; ie1Recv(numR) = ie1
 js1Recv(numR) = je1 - shalo + 1; je1Recv(numR) = je1
 is2Recv(numR) = ie2 - whalo + 1; ie2Recv(numR) = ie2
 js2Recv(numR) = je2 - shalo + 1; je2Recv(numR) = je2 
end if

!--- northwest for recving
found_corner = .false.
n = wCont(tileMe)%ncontact
if( n > 0) then
 if(wCont(tileMe)%je1(n) == jeg(tileMe) ) then
    tn = wCont(tileMe)%tile(n)
    if(wcont(tileMe)%je2(n) < jeg(tn) ) then  ! the corner tile is tn.
       if( jeg(tn) - wcont(tileMe)%je2(n) < nhalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: northwest tile for recv 1 is not tiled properly")
       found_corner = .true.; tc = tn
       ie1 = wCont(tileMe)%is1(n) - 1; js1 = wCont(tileMe)%je1(n) + 1
       ie2 = wCont(tileMe)%is2(n);     js2 = wCont(tileMe)%je2(n) + 1
    else if(nCont(tn)%ncontact >0) then ! the corner tile may be south tile of tn.
       m = nCont(tn)%ncontact
       if(nCont(tn)%ie1(m) == ieg(tn)) then ! corner is nc.
          found_corner = .true.; tc = nCont(tn)%tile(m)
          ie1 = wCont(tileMe)%is1(n) - 1; js1 = wCont(tileMe)%je1(n) + 1
          ie2 = nCont(tn)%ie2(m);         js2 = nCont(tn)%js2(m)
       end if
    endif
 endif
end if
if( .not. found_corner ) then  ! not found,
 if( nCont(tileMe)%ncontact > 0) then
    if( nCont(tileMe)%is1(1) == isg(tileMe)) then
       tn = nCont(tileMe)%tile(1)
       if(nCont(tileMe)%is2(1) > isg(tn) ) then  ! the corner tile is tn.
          if( ncont(tileMe)%is2(1)-isg(tn) < whalo ) call mpp_error(FATAL, &
               "mpp_domains_define.inc: northwest tile for recv 2 is not tiled properly")
          found_corner = .true.; tc = tn
          ie1 = nCont(tileMe)%is1(1) - 1; js1 = nCont(tileMe)%je1(1) + 1
          ie2 = nCont(tileMe)%is2(1) - 1; js2 = nCont(tileMe)%js2(1)
       else if(wCont(tn)%ncontact >0) then ! the corner tile may be west tile of tn.
          if(wCont(tn)%js1(1) == jsg(tn)) then ! corner is nc.
             found_corner = .true.; tc = wCont(tn)%tile(1)
             ie1 = nCont(tileMe)%is1(1) - 1; js1 = nCont(tileMe)%je1(1) + 1
             ie2 = wCont(tn)%ie2(1);         js2 = wCont(tn)%js2(1)
          end if
       end if
    end if
 end if
end if
if(found_corner) then
 numR = numR + 1
 tileRecv(numR) = tc; align1Recv(numR) =NORTH_WEST;  align2Recv(numR) = SOUTH_EAST
 is1Recv(numR) = ie1 - whalo + 1; ie1Recv(numR) = ie1
 js1Recv(numR) = js1;             je1Recv(numR) = js1 + nhalo - 1
 is2Recv(numR) = ie2 - whalo + 1; ie2Recv(numR) = ie2
 js2Recv(numR) = js2;             je2Recv(numR) = js2 + nhalo - 1  
end if

!--- northeast for recving
found_corner = .false.
n = eCont(tileMe)%ncontact
if( n > 0) then
 if(eCont(tileMe)%je1(n) == jeg(tileMe) ) then
    tn = eCont(tileMe)%tile(n)
    if(econt(tileMe)%je2(n) < jeg(tn) ) then  ! the corner tile is tn.
       if( jeg(tn) - econt(tileMe)%je2(n) < nhalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: northeast tile for recv 1 is not tiled properly")
       found_corner = .true.; tc = tn
       is1 = eCont(tileMe)%ie1(n) + 1; js1 = eCont(tileMe)%je1(n) + 1
       is2 = eCont(tileMe)%is2(1);     js2 = eCont(tileMe)%je2(1) + 1
    else if(nCont(tn)%ncontact >0) then ! the corner tile may be south tile of tn.
       if(nCont(tn)%is1(1) == isg(tn)) then ! corner is nc.
          found_corner = .true.; tc = nCont(tn)%tile(1)
          is1 = eCont(tileMe)%ie1(n) + 1; js1 = eCont(tileMe)%je1(n) + 1
          is2 = nCont(tn)%is2(1);         js2 = nCont(tn)%js2(1)
       end if
    end if
 end if
end if
if( .not. found_corner ) then  ! not found,
 n = nCont(tileMe)%ncontact
 if( n > 0) then
    if( nCont(tileMe)%ie1(n) == ieg(tileMe)) then
       tn = nCont(tileMe)%tile(n)
       if(nCont(tileMe)%ie2(n) < ieg(tn) ) then  ! the corner tile is tn.
          if(ieg(tn) - sCont(tileMe)%ie2(n) < ehalo ) call mpp_error(FATAL, &
               "mpp_domains_define.inc: northeast tile for recv 2 is not tiled properly")
          found_corner = .true.; tc = tn
          is1 = sCont(tileMe)%ie1(n) + 1; js1 = sCont(tileMe)%je1(n) + 1
          is2 = sCont(tileMe)%ie2(n) + 1; js2 = sCont(tileMe)%js2(n)
       else if(eCont(tn)%ncontact >0) then ! the corner tile may be east tile of tn.
          if(eCont(tn)%js1(1) == jsg(tn)) then ! corner is nc.
             found_corner = .true.; tc = eCont(tn)%tile(1)
             is1 = sCont(tileMe)%ie1(n) + 1; js1 = sCont(tileMe)%je1(n) + 1
             is2 = eCont(tn)%is2(m);         js2 = eCont(tn)%js2(m)
          end if
       end if
    end if
 end if
end if
if(found_corner) then
 numR = numR + 1
 tileRecv(numR) = tc; align1Recv(numR) =NORTH_EAST;  align2Recv(numR) = SOUTH_WEST
 is1Recv(numR) = is1; ie1Recv(numR) = is1 + ehalo - 1
 js1Recv(numR) = js1; je1Recv(numR) = js1 + nhalo - 1
 is2Recv(numR) = is2; ie2Recv(numR) = is2 + ehalo - 1
 js2Recv(numR) = js2; je2Recv(numR) = js2 + nhalo - 1  
end if

!--- to_pe's southeast for sending
do n = 1, wCont(tileMe)%ncontact
 tn = wCont(tileMe)%tile(n)
 if(wCont(tileMe)%js2(n) == jsg(tn) ) then
    if(wcont(tileMe)%js1(n) > jsg(tileMe) ) then  ! send to  tile tn.
       if( wcont(tileMe)%js1(n) - jsg(tileMe) < shalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: southeast tile for send 1 is not tiled properly")
       numS = numS+1; tileSend(numS) = tn
       align1Send(numS) = NORTH_WEST;  align2Send(numS) = SOUTH_EAST
       is1Send(numS) = wCont(tileMe)%is1(n);     ie1Send(numS) = is1Send(numS) + ehalo - 1  
       je1Send(numS) = wCont(tileMe)%js1(n) - 1; js1Send(numS) = je1Send(numS) - shalo + 1
       is2Send(numS) = wCont(tileMe)%ie2(n) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
       je2Send(numS) = wCont(tileMe)%js2(n) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
    end if
 end if
end do
do n = 1, nCont(tileMe)%ncontact
 tn = nCont(tileMe)%tile(n)
 if(nCont(tileMe)%ie2(n) == ieg(tn) ) then
    if(nCont(tileMe)%ie1(n) < ieg(tileMe) ) then  ! send to  tile tn.
       if( ieg(tileMe) - nCont(tileMe)%ie1(n) < ehalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: southeast tile for send 2 is not tiled properly")
       numS = numS+1; tileSend(numS) = tn
       align1Send(numS) = NORTH_WEST;  align2Send(numS) = SOUTH_EAST
       is1Send(numS) = nCont(tileMe)%ie1(n) + 1; ie1Send(numS) = is1Send(numS) + ehalo - 1  
       je1Send(numS) = nCont(tileMe)%je1(n)    ; js1Send(numS) = je1Send(numS) - shalo + 1
       is2Send(numS) = nCont(tileMe)%ie2(n) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
       je2Send(numS) = nCont(tileMe)%je2(n) - 1; js2Send(numS) = je2Send(numS) - shalo + 1
    end if
 end if
end do

!--- found the corner overlap that is not specified through contact line.
n = wCont(tileMe)%ncontact
found_corner = .false.
if( n > 0) then
 tn = wCont(tileMe)%tile(n)
 if( wCont(tileMe)%je1(n) == jeg(tileMe) .AND. wCont(tileMe)%je2(n) == jeg(tn) ) then
    m = nCont(tn)%ncontact
    if(m >0) then
       tc = nCont(tn)%tile(m)
       if( nCont(tn)%ie1(m) == ieg(tn) .AND. nCont(tn)%ie2(m) == ieg(tc) ) found_corner = .true.
    end if
 end if
end if
if( .not. found_corner ) then  ! not found, then starting from north contact
 if( nCont(tileMe)%ncontact > 0) then
    tn = nCont(tileMe)%tile(1)
    if( nCont(tileMe)%is1(1) == isg(tileMe) .AND. nCont(tileMe)%is2(1) == isg(tn) ) then
       if(wCont(tn)%ncontact >0) then 
          tc = wCont(tn)%tile(1)
          if( wCont(tn)%js1(1) == jsg(tn) .AND. wCont(tn)%js2(1) == jsg(tc) ) found_corner = .true.
       end if
    end if
 end if
end if

if(found_corner) then
 numS = numS+1; tileSend(numS) = tc
 align1Send(numS) = NORTH_WEST;  align2Send(numS) = SOUTH_EAST
 is1Send(numS) = isg(tileMe); ie1Send(numS) = is1Send(numS) + ehalo - 1  
 je1Send(numS) = jeg(tileMe); js1Send(numS) = je1Send(numS) - shalo + 1
 is2Send(numS) = ieg(tc) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
 je2Send(numS) = jsg(tc) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
end if

!--- to_pe's southwest for sending
do n = 1, eCont(tileMe)%ncontact
 tn = eCont(tileMe)%tile(n)
 if(eCont(tileMe)%js2(n) == jsg(tn) ) then
    if(econt(tileMe)%js1(n) > jsg(tileMe) ) then  ! send to  tile tn.
       if( econt(tileMe)%js1(n) - jsg(tileMe) < shalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: southwest tile for send 1 is not tiled properly")
       numS = numS+1; tileSend(numS) = tn
       align1Send(numS) = NORTH_EAST;  align2Send(numS) = SOUTH_WEST
       ie1Send(numS) = eCont(tileMe)%ie1(n);     is1Send(numS) = ie1Send(numS) - whalo + 1  
       je1Send(numS) = eCont(tileMe)%js1(n) - 1; js1Send(numS) = je1Send(numS) - shalo + 1
       ie2Send(numS) = eCont(tileMe)%is2(n) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
       je2Send(numS) = eCont(tileMe)%js2(n) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
    end if
 end if
end do
do n = 1, nCont(tileMe)%ncontact
 tn = nCont(tileMe)%tile(n)
 if(nCont(tileMe)%is2(n) == isg(tn) ) then
    if(ncont(tileMe)%is1(n) > isg(tileMe) ) then  ! send to  tile tn.
       if( ncont(tileMe)%is1(n) - isg(tileMe) < whalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: southwest tile for send 2 is not tiled properly")
       numS = numS+1; tileSend(numS) = tn
       align1Send(numS) = NORTH_EAST;  align2Send(numS) = SOUTH_WEST
       ie1Send(numS) = nCont(tileMe)%is1(n) - 1; is1Send(numS) = ie1Send(numS) - whalo + 1  
       ie1Send(numS) = nCont(tileMe)%je1(n)    ; js1Send(numS) = je1Send(numS) - shalo + 1
       ie2Send(numS) = nCont(tileMe)%is2(n) - 1; is2Send(numS) = je2Send(numS) - whalo + 1
       je2Send(numS) = nCont(tileMe)%js2(n) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
    end if
 end if
end do

!--- found the corner overlap that is not specified through contact line.
n = eCont(tileMe)%ncontact
found_corner = .false.
if( n > 0) then
 tn = eCont(tileMe)%tile(n)
 if( eCont(tileMe)%je1(n) == jeg(tileMe) .AND. eCont(tileMe)%je2(n) == jeg(tn) ) then
    if(nCont(tn)%ncontact >0) then
       tc = nCont(tn)%tile(1)
       if( nCont(tn)%is1(1) == isg(tn) .AND. nCont(tn)%is2(n) == isg(tc) ) found_corner = .true.
    end if
 end if
end if
if( .not. found_corner ) then  ! not found, then starting from north contact
 n = nCont(tileMe)%ncontact
 if( n > 0) then
    tn = nCont(tileMe)%tile(n)
    if( nCont(tileMe)%ie1(n) == ieg(tileMe) .AND. nCont(tileMe)%ie2(n) == ieg(tn) ) then
       if(eCont(tn)%ncontact >0) then
          tc = eCont(tn)%tile(1)
          if( eCont(tn)%js1(1) == jsg(tn) .AND. eCont(tn)%js2(n) == jsg(tc) ) found_corner = .true.
       end if
    end if
 end if
end if

if(found_corner) then
 numS = numS+1; tileSend(numS) = tc
 align1Send(numS) = NORTH_EAST;  align2Send(numS) = SOUTH_WEST
 ie1Send(numS) = ieg(tileMe); is1Send(numS) = ie1Send(numS) - whalo + 1  
 je1Send(numS) = jeg(tileMe); js1Send(numS) = je1Send(numS) - shalo + 1
 ie2Send(numS) = isg(tc) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
 je2Send(numS) = jsg(tc) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
end if

!--- to_pe's northwest for sending
do n = 1, eCont(tileMe)%ncontact
 tn = eCont(tileMe)%tile(n)
 if(eCont(tileMe)%je2(n) == jeg(tn) ) then
    if(econt(tileMe)%je1(n) < jeg(tileMe) ) then  ! send to  tile tn.
       if( jeg(tileMe) - econt(tileMe)%je1(n) < nhalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: northwest tile for send 1 is not tiled properly")
       numS = numS+1; tileSend(numS) = tn
       align1Send(numS) = SOUTH_EAST;  align2Send(numS) = NORTH_WEST
       ie1Send(numS) = eCont(tileMe)%ie1(n)    ; is1Send(numS) = ie1Send(numS) - whalo + 1  
       js1Send(numS) = eCont(tileMe)%je1(n) + 1; je1Send(numS) = js1Send(numS) + nhalo - 1
       ie2Send(numS) = eCont(tileMe)%is2(n) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
       js2Send(numS) = eCont(tileMe)%je2(n) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
    end if
 end if
end do

do n = 1, sCont(tileMe)%ncontact
 tn = sCont(tileMe)%tile(n)
 if(sCont(tileMe)%is2(n) == isg(tn) ) then
    if(scont(tileMe)%is1(n) > isg(tileMe) ) then  ! send to  tile tn.
       if( scont(tileMe)%is1(n) - isg(tileMe) < whalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: southwest tile for send 2 is not tiled properly")
       numS = numS+1; tileSend(numS) = tn
       align1Send(numS) = SOUTH_EAST;  align2Send(numS) = NORTH_WEST
       ie1Send(numS) = nCont(tileMe)%is1(n) - 1; is1Send(numS) = ie1Send(numS) - whalo + 1  
       js1Send(numS) = nCont(tileMe)%je1(n)    ; je1Send(numS) = js1Send(numS) + nhalo - 1
       ie2Send(numS) = nCont(tileMe)%is2(n) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
       js2Send(numS) = nCont(tileMe)%je2(n) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
    end if
 end if
end do

!--- found the corner overlap that is not specified through contact line.
n = eCont(tileMe)%ncontact
found_corner = .false.
if( n > 0) then
 tn = eCont(tileMe)%tile(1)
 if( eCont(tileMe)%js1(1) == jsg(tileMe) .AND. eCont(tileMe)%js2(1) == jsg(tn) ) then
    if(sCont(tn)%ncontact >0) then
       tc = sCont(tn)%tile(1)
       if( sCont(tn)%is1(1) == isg(tn) .AND. sCont(tn)%is2(1) == isg(tc) ) found_corner = .true.
    end if
 end if
end if
if( .not. found_corner ) then  ! not found, then starting from north contact
 n = sCont(tileMe)%ncontact
 found_corner = .false.
 if( n > 0) then
    tn = sCont(tileMe)%tile(n)
    if( sCont(tileMe)%ie1(n) == ieg(tileMe) .AND. sCont(tileMe)%ie2(n) == ieg(tn) ) then
       if(eCont(tn)%ncontact >0) then
          tc = eCont(tn)%tile(n)
          if( eCont(tn)%je1(n) == jeg(tn) .AND. eCont(tn)%je2(n) == jeg(tc) ) found_corner = .true.
       end if
    end if
 end if
end if

if(found_corner) then
 numS = numS+1; tileSend(numS) = tc
 align1Send(numS) = SOUTH_EAST;  align2Send(numS) = NORTH_WEST
 ie1Send(numS) = ieg(tileMe); is1Send(numS) = ie1Send(numS) - whalo + 1  
 js1Send(numS) = jsg(tileMe); je1Send(numS) = js1Send(numS) + nhalo - 1
 ie2Send(numS) = isg(tc) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
 js2Send(numS) = jeg(tc) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
end if

!--- to_pe's northeast for sending
do n = 1, wCont(tileMe)%ncontact
 tn = wCont(tileMe)%tile(n)
 if(wCont(tileMe)%je2(n) == jeg(tn) ) then
    if(wcont(tileMe)%je1(n) < jeg(tileMe) ) then  ! send to  tile tn.
       if( jeg(tileMe) - wcont(tileMe)%je1(n) < nhalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: northeast tile for send 1 is not tiled properly")
       numS = numS+1; tileSend(numS) = tn
       align1Send(numS) = SOUTH_WEST;  align2Send(numS) = NORTH_EAST
       is1Send(numS) = wCont(tileMe)%is1(n)    ; ie1Send(numS) = is1Send(numS) + ehalo - 1  
       js1Send(numS) = wCont(tileMe)%je1(n) + 1; je1Send(numS) = js1Send(numS) + nhalo - 1
       is2Send(numS) = wCont(tileMe)%ie2(n) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
       js2Send(numS) = wCont(tileMe)%je2(n) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
    end if
 end if
end do

do n = 1, sCont(tileMe)%ncontact
 tn = sCont(tileMe)%tile(n)
 if(sCont(tileMe)%ie2(n) == ieg(tn) ) then
    if(sCont(tileMe)%ie1(n) < ieg(tileMe) ) then  ! send to  tile tn.
       if( ieg(tileMe) - sCont(tileMe)%ie1(n) < ehalo ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: southeast tile for send 2 is not tiled properly")
       numS = numS+1; tileSend(numS) = tn
       align1Send(numS) = SOUTH_WEST;  align2Send(numS) = NORTH_EAST
       is1Send(numS) = sCont(tileMe)%ie1(n) + 1; ie1Send(numS) = is1Send(numS) + ehalo - 1  
       js1Send(numS) = sCont(tileMe)%js1(n)    ; je1Send(numS) = js1Send(numS) + nhalo - 1
       is2Send(numS) = sCont(tileMe)%ie2(n) + 1; ie2Send(numS) = is1Send(numS) + ehalo - 1
       js2Send(numS) = sCont(tileMe)%je2(n) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
    end if
 end if
end do

!--- found the corner overlap that is not specified through contact line.
n = wCont(tileMe)%ncontact
found_corner = .false.
if( n > 0) then
 tn = wCont(tileMe)%tile(1)
 if( wCont(tileMe)%js1(n) == jsg(tileMe) .AND. wCont(tileMe)%js2(n) == jsg(tn) ) then
    m = sCont(tn)%ncontact
    if(m >0) then
       tc = sCont(tn)%tile(m)
       if( sCont(tn)%ie1(m) == ieg(tn) .AND. sCont(tn)%ie2(m) == ieg(tc) ) found_corner = .true.
    end if
 end if
end if
if( .not. found_corner ) then  ! not found, then starting from north contact
 n = sCont(tileMe)%ncontact
 found_corner = .false.
 if( n > 0) then
    tn = sCont(tileMe)%tile(1)
    if( sCont(tileMe)%is1(1) == isg(tileMe) .AND. sCont(tileMe)%is2(1) == isg(tn) ) then
       m = wCont(tn)%ncontact
       if( m > 0 ) then
          tc = wCont(tn)%tile(m)
          if( wCont(tn)%je1(m) == jeg(tn) .AND. wCont(tn)%je2(m) == jeg(tc) ) found_corner = .true.
       end if
    end if
 end if
end if
if(found_corner) then
 numS = numS+1; tileSend(numS) = tc
 align1Send(numS) = SOUTH_WEST;  align2Send(numS) = NORTH_EAST
 is1Send(numS) = isg(tileMe); ie1Send(numS) = is1Send(numS) + ehalo - 1  
 js1Send(numS) = jsg(tileMe); je1Send(numS) = js1Send(numS) + nhalo - 1
 is2Send(numS) = ieg(tc) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
 js2Send(numS) = jeg(tc) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
end if

end subroutine fill_corner_contact

!--- find the alignment direction, check if index is reversed, if reversed, exchange index.
subroutine check_alignment( is, ie, js, je, isg, ieg, jsg, jeg, alignment )
integer, intent(inout) :: is, ie, js, je, isg, ieg, jsg, jeg
integer, intent(out)   :: alignment

integer :: i, j

if ( is == ie ) then      ! x-alignment
 if ( is == isg ) then
    alignment = WEST
 else if ( is == ieg ) then
    alignment = EAST
 else
    call mpp_error(FATAL, 'mpp_domains_define.inc: The contact region is not on the x-boundary of the tile')
 end if
 if ( js > je ) then
    j = js; js = je; je = j
 end if
else if ( js == je ) then ! y-alignment
 if ( js == jsg ) then
    alignment = SOUTH
 else if ( js == jeg ) then
    alignment = NORTH
 else
    call mpp_error(FATAL, 'mpp_domains_define.inc: The contact region is not on the y-boundary of the tile')
 end if
 if ( is > ie ) then
    i = is; is = ie; ie = i
 end if
else
 call mpp_error(FATAL, 'mpp_domains_define.inc: The contact region should be line contact' )
end if

end subroutine check_alignment
!#####################################################################  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_MODIFY_DOMAIN: modify extent of domain                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <SUBROUTINE NAME="mpp_modify_domain1D" INTERFACE="mpp_modify_domain">
!   <IN NAME="domain_in" TYPE="type(domain1D)" > </IN>
!   <IN NAME="hbegin,hend" TYPE="integer,optional" > </IN>
!   <IN NAME="cbegin,cend" TYPE="integer,optional" > </IN>
!   <IN NAME="gbegin,gend" TYPE="integer,optional" > </IN>
!   <INOUT NAME="domain_out" TYPE="type(domain1D)" > </INOUT>

! <PUBLICROUTINE>
subroutine mpp_modify_domain1D(domain_in,domain_out,cbegin,cend,gbegin,gend, hbegin, hend)
  ! </PUBLICROUTINE>
type(domain1D), intent(in)    :: domain_in
type(domain1D), intent(inout) :: domain_out
integer, intent(in), optional :: hbegin, hend             ! halo size 
integer, intent(in), optional :: cbegin, cend             ! extent of compute_domain
integer, intent(in), optional :: gbegin, gend             ! extent of global domain
integer :: ndivs, global_indices(2) !(/ isg, ieg /)
integer                       :: flag
! get the global indices of the input domain
global_indices(1) = domain_in%global%begin;  global_indices(2) = domain_in%global%end

! get the layout
ndivs = size(domain_in%list(:))

! get the flag
flag = 0
if(domain_in%cyclic) flag = flag + CYCLIC_GLOBAL_DOMAIN
if(domain_in%data%is_global) flag = flag + GLOBAL_DATA_DOMAIN

call mpp_define_domains( global_indices, ndivs, domain_out, pelist = domain_in%list(:)%pe, &
   flags = flag, begin_halo = hbegin, end_halo = hend, extent = domain_in%list(:)%compute%size )

if(present(cbegin)) domain_out%compute%begin = cbegin
if(present(cend))   domain_out%compute%end = cend
domain_out%compute%size = domain_out%compute%end - domain_out%compute%begin + 1
if(present(gbegin)) domain_out%global%begin = gbegin
if(present(gend))   domain_out%global%end = gend
domain_out%global%size = domain_out%global%end - domain_out%global%begin + 1

end subroutine mpp_modify_domain1D
! </SUBROUTINE>

!#######################################################################
!----------------------------------------------------------------------------------
! <SUBROUTINE NAME="mpp_modify_domain2D" INTERFACE="mpp_modify_domain">
!   <IN NAME="domain_in" TYPE="type(domain2D)" > </IN>
!   <IN NAME="isc,iec" TYPE="integer,optional" > </IN>
!   <IN NAME="jsc,jec" TYPE="integer,optional" > </IN>
!   <IN NAME="isg,ieg" TYPE="integer,optional" > </IN>
!   <IN NAME="jsg,jeg" TYPE="integer,optional" > </IN>
!   <IN NAME="whalo,ehalo" TYPE="integer,optional" > </IN>
!   <IN NAME="shalo,nhalo" TYPE="integer,optional" > </IN>
!   <INOUT NAME="domain_out" TYPE="type(domain2D)" > </INOUT>

! <PUBLICROUTINE>
subroutine mpp_modify_domain2D(domain_in, domain_out, isc, iec, jsc, jec, isg, ieg, jsg, jeg, whalo, ehalo, shalo, nhalo)
  ! </PUBLICROUTINE>
type(domain2D), intent(in)    :: domain_in
type(domain2D), intent(inout) :: domain_out
integer, intent(in), optional :: isc, iec, jsc, jec
integer, intent(in), optional :: isg, ieg, jsg, jeg
integer, intent(in), optional :: whalo, ehalo, shalo, nhalo
integer                       :: global_indices(4), layout(2)
integer                       :: xflag, yflag, nlist, i

if(present(whalo) .or. present(ehalo) .or. present(shalo) .or. present(nhalo) ) then
   ! get the global indices of the input domain
 global_indices(1) = domain_in%x(1)%global%begin;  global_indices(2) = domain_in%x(1)%global%end
 global_indices(3) = domain_in%y(1)%global%begin;  global_indices(4) = domain_in%y(1)%global%end

 ! get the layout
 layout(1) = size(domain_in%x(1)%list(:)); layout(2) = size(domain_in%y(1)%list(:))

 ! get the flag
 xflag = 0; yflag = 0
 if(domain_in%x(1)%cyclic) xflag = xflag + CYCLIC_GLOBAL_DOMAIN
 if(domain_in%x(1)%data%is_global) xflag = xflag + GLOBAL_DATA_DOMAIN
 if(domain_in%y(1)%cyclic) yflag = yflag + CYCLIC_GLOBAL_DOMAIN
 if(domain_in%y(1)%data%is_global) yflag = yflag + GLOBAL_DATA_DOMAIN

 call mpp_define_domains( global_indices, layout, domain_out, pelist = domain_in%list(:)%pe, &
      xflags = xflag, yflags = yflag,  whalo = whalo, ehalo = ehalo,     &
      shalo = shalo, nhalo = nhalo,                                      &
      xextent = domain_in%x(1)%list(:)%compute%size,                     &
      yextent = domain_in%y(1)%list(:)%compute%size,                     &
      symmetry=domain_in%symmetry,                                       &
      maskmap = domain_in%pearray .NE. NULL_PE )

else    
 call mpp_define_null_domain(domain_out)
 nlist = size(domain_in%list(:))
 allocate(domain_out%list(0:nlist-1) )
 do i = 0, nlist-1
    allocate(domain_out%list(i)%tile_id(1))
    domain_out%list(i)%tile_id(1) = 1
 enddo
 call mpp_modify_domain(domain_in%x(1), domain_out%x(1), isc, iec, isg, ieg)
 call mpp_modify_domain(domain_in%y(1), domain_out%y(1), jsc, jec, jsg, jeg)
endif

end subroutine mpp_modify_domain2D
! </SUBROUTINE>

!#####################################################################


subroutine mpp_define_null_domain1D(domain)
type(domain1D), intent(inout) :: domain

domain%global%begin  = -1; domain%global%end  = -1; domain%global%size = 0
domain%data%begin    = -1; domain%data%end    = -1; domain%data%size = 0
domain%compute%begin = -1; domain%compute%end = -1; domain%compute%size = 0
domain%pe = NULL_PE

end subroutine mpp_define_null_domain1D

!#####################################################################


subroutine mpp_define_null_domain2D(domain)
type(domain2D), intent(inout) :: domain

allocate(domain%x(1), domain%y(1), domain%tile_id(1))
call mpp_define_null_domain(domain%x(1))
call mpp_define_null_domain(domain%y(1))
domain%pe = NULL_PE
domain%tile_id(1)   = 1
domain%ntiles       = 1
domain%max_ntile_pe = 1
domain%ncontacts    = 0

end subroutine mpp_define_null_domain2D

!####################################################################

subroutine mpp_deallocate_domain1D(domain)
  type(domain1D), intent(inout) :: domain

  if(ASSOCIATED(domain%list)) deallocate(domain%list)

end subroutine mpp_deallocate_domain1D

!####################################################################

subroutine mpp_deallocate_domain2D(domain)
  type(domain2D), intent(inout) :: domain

  call deallocate_domain2D_local(domain)
  if(ASSOCIATED(domain%io_domain) ) then
     call deallocate_domain2D_local(domain%io_domain) 
     deallocate(domain%io_domain)
  endif

end subroutine mpp_deallocate_domain2D

!##################################################################

subroutine deallocate_domain2D_local(domain)
type(domain2D), intent(inout) :: domain
integer :: i, nlist, ntileMe

ntileMe = size(domain%x(:))

if(ASSOCIATED(domain%pearray))deallocate(domain%pearray)
do i = 1, ntileMe
   call mpp_deallocate_domain1D(domain%x(i))
   call mpp_deallocate_domain1D(domain%y(i))
enddo
deallocate(domain%x, domain%y, domain%tile_id)

if(ASSOCIATED(domain%list)) then
 do i = 0, size(domain%list(:))-1
    deallocate(domain%list(i)%x, domain%list(i)%y, domain%list(i)%tile_id)
 enddo
 deallocate(domain%list)
endif

if(ASSOCIATED(domain%check_C))  call deallocate_overlapSpec(domain%check_C)
if(ASSOCIATED(domain%check_E))  call deallocate_overlapSpec(domain%check_E)
if(ASSOCIATED(domain%check_N))  call deallocate_overlapSpec(domain%check_N)
if(ASSOCIATED(domain%bound_C))  call deallocate_overlapSpec(domain%bound_C)
if(ASSOCIATED(domain%bound_E))  call deallocate_overlapSpec(domain%bound_E)
if(ASSOCIATED(domain%bound_N))  call deallocate_overlapSpec(domain%bound_N)
if(ASSOCIATED(domain%update_T)) call deallocate_overlapSpec(domain%update_T)
if(ASSOCIATED(domain%update_E)) call deallocate_overlapSpec(domain%update_E)
if(ASSOCIATED(domain%update_C)) call deallocate_overlapSpec(domain%update_C)
if(ASSOCIATED(domain%update_N)) call deallocate_overlapSpec(domain%update_N)

end subroutine deallocate_domain2D_local

!####################################################################

subroutine allocate_check_overlap(overlap, count)
  type(overlap_type), intent(inout) :: overlap
  integer,            intent(in   ) :: count

  overlap%count = 0
  overlap%pe    = NULL_PE
  if(associated(overlap%tileMe)) call mpp_error(FATAL, &
       "allocate_check_overlap(mpp_domains_define): overlap is already been allocated")
  if(count < 1)  call mpp_error(FATAL, &
       "allocate_check_overlap(mpp_domains_define): count should be a positive integer")
  allocate(overlap%tileMe  (count), overlap%dir(count) )
  allocate(overlap%is      (count), overlap%ie (count) )
  allocate(overlap%js      (count), overlap%je (count) )
  allocate(overlap%rotation(count)                     )
  overlap%rotation     =  ZERO

end subroutine allocate_check_overlap

!#######################################################################
subroutine insert_check_overlap(overlap, pe, tileMe, dir, rotation, is, ie, js, je)
  type(overlap_type), intent(inout) :: overlap
  integer,            intent(in   ) :: pe
  integer,            intent(in   ) :: tileMe, dir, rotation
  integer,            intent(in   ) :: is, ie, js, je
  integer                           :: count

  overlap%count = overlap%count + 1
  count = overlap%count
  if(.NOT. associated(overlap%tileMe)) call mpp_error(FATAL, &
     "mpp_domains_define.inc(insert_check_overlap): overlap is not assigned any memory")
  if(count > size(overlap%tileMe(:)) ) call mpp_error(FATAL, &
     "mpp_domains_define.inc(insert_check_overlap): overlap%count is greater than size(overlap%tileMe)")
  if( overlap%pe == NULL_PE ) then
     overlap%pe = pe
  else
     if(overlap%pe .NE. pe) call mpp_error(FATAL,  &
           "mpp_domains_define.inc(insert_check_overlap): mismatch on pe")
  endif
  overlap%tileMe  (count) = tileMe
  overlap%dir     (count) = dir
  overlap%rotation(count) = rotation
  overlap%is      (count) = is
  overlap%ie      (count) = ie
  overlap%js      (count) = js
  overlap%je      (count) = je

end subroutine insert_check_overlap

!#######################################################################
!--- this routine add the overlap_in into overlap_out
subroutine add_check_overlap( overlap_out, overlap_in)
  type(overlap_type), intent(inout) :: overlap_out
  type(overlap_type), intent(in   ) :: overlap_in
  type(overlap_type)                :: overlap
  integer                          :: count, count_in, count_out

  ! if overlap_out%count == 0, then just copy overlap_in to overlap_out
  count_in  = overlap_in %count
  count_out = overlap_out%count
  count     = count_in+count_out
  if(count_in == 0) call mpp_error(FATAL, &
       "add_check_overlap(mpp_domains_define): overlap_in%count is zero")

  if(count_out == 0) then
     if(associated(overlap_out%tileMe)) call mpp_error(FATAL, &
          "add_check_overlap(mpp_domains_define): overlap is already been allocated but count=0")
     call allocate_check_overlap(overlap_out, count_in)
     overlap_out%pe  = overlap_in%pe
  else  ! need to expand the dimension size of overlap
     call allocate_check_overlap(overlap, count_out)
     if(overlap_out%pe .NE. overlap_in%pe) call mpp_error(FATAL, &
         "mpp_domains_define.inc(add_check_overlap): mismatch of pe between overlap_in and overlap_out")
     overlap%tileMe      (1:count_out) = overlap_out%tileMe      (1:count_out)
     overlap%is          (1:count_out) = overlap_out%is          (1:count_out)
     overlap%ie          (1:count_out) = overlap_out%ie          (1:count_out)
     overlap%js          (1:count_out) = overlap_out%js          (1:count_out)
     overlap%je          (1:count_out) = overlap_out%je          (1:count_out)
     overlap%dir         (1:count_out) = overlap_out%dir         (1:count_out)
     overlap%rotation    (1:count_out) = overlap_out%rotation    (1:count_out)
     call deallocate_overlap_type(overlap_out)
     call allocate_check_overlap(overlap_out, count)
     overlap_out%tileMe      (1:count_out) = overlap%tileMe      (1:count_out)
     overlap_out%is          (1:count_out) = overlap%is          (1:count_out)
     overlap_out%ie          (1:count_out) = overlap%ie          (1:count_out)
     overlap_out%js          (1:count_out) = overlap%js          (1:count_out)
     overlap_out%je          (1:count_out) = overlap%je          (1:count_out)
     overlap_out%dir         (1:count_out) = overlap%dir         (1:count_out)
     overlap_out%rotation    (1:count_out) = overlap%rotation    (1:count_out)
     call deallocate_overlap_type(overlap)
  end if
  overlap_out%count                           = count
  overlap_out%tileMe      (count_out+1:count) = overlap_in%tileMe      (1:count_in)
  overlap_out%is          (count_out+1:count) = overlap_in%is          (1:count_in)
  overlap_out%ie          (count_out+1:count) = overlap_in%ie          (1:count_in)
  overlap_out%js          (count_out+1:count) = overlap_in%js          (1:count_in)
  overlap_out%je          (count_out+1:count) = overlap_in%je          (1:count_in)
  overlap_out%dir         (count_out+1:count) = overlap_in%dir         (1:count_in)
  overlap_out%rotation    (count_out+1:count) = overlap_in%rotation    (1:count_in)

end subroutine add_check_overlap

!####################################################################
subroutine init_overlap_type(overlap)
  type(overlap_type), intent(inout) :: overlap

  overlap%count = 0
  overlap%pe    = NULL_PE  

end subroutine init_overlap_type

!####################################################################

subroutine allocate_update_overlap( overlap, count)
  type(overlap_type), intent(inout) :: overlap
  integer,           intent(in   ) :: count

  overlap%count = 0
  overlap%pe    = NULL_PE
  if(associated(overlap%tileMe)) call mpp_error(FATAL, &
       "allocate_update_overlap(mpp_domains_define): overlap is already been allocated")
  if(count < 1)  call mpp_error(FATAL, &
       "allocate_update_overlap(mpp_domains_define): count should be a positive integer")
  allocate(overlap%tileMe      (count), overlap%tileNbr (count) )
  allocate(overlap%is          (count), overlap%ie      (count) )
  allocate(overlap%js          (count), overlap%je      (count) )
  allocate(overlap%isMe        (count), overlap%ieMe    (count) )
  allocate(overlap%jsMe        (count), overlap%jeMe    (count) )
  allocate(overlap%dir         (count), overlap%rotation(count) )
  allocate(overlap%is_refined  (count), overlap%index   (count) )
  allocate(overlap%from_contact(count)                          )
  overlap%is_refined   = .FALSE.
  overlap%rotation     =  ZERO
  overlap%from_contact = .FALSE.
  overlap%index        = -1   

end subroutine allocate_update_overlap

  !#####################################################################################
  subroutine insert_update_overlap(overlap, pe, is1, ie1, js1, je1, is2, ie2, js2, je2, dir, reverse, symmetry)
    type(overlap_type), intent(inout) :: overlap
    integer,            intent(in   ) :: pe
    integer,            intent(in   ) :: is1, ie1, js1, je1, is2, ie2, js2, je2
    integer,            intent(in   ) :: dir
    logical, optional,  intent(in   ) :: reverse, symmetry

    logical :: is_reverse, is_symmetry, is_overlapped
    integer :: is, ie, js, je, count

    is_reverse = .FALSE.
    if(PRESENT(reverse)) is_reverse = reverse
    is_symmetry = .FALSE.
    if(PRESENT(symmetry)) is_symmetry = symmetry

    is = max(is1,is2); ie = min(ie1,ie2)
    js = max(js1,js2); je = min(je1,je2)
    is_overlapped = .false.
    !--- to avoid unnecessary ( duplicate overlap ) for symmetry domain
    if(is_symmetry .AND. (dir == 1 .OR. dir == 5)) then  ! x-direction
       if( ie .GE. is .AND. je .GT. js ) is_overlapped = .true.
    else if(is_symmetry .AND. (dir == 3 .OR. dir == 7)) then ! y-direction
       if( ie .GT. is .AND. je .GE. js ) is_overlapped = .true.
    else if(ie.GE.is .AND. je.GE.js )then
       is_overlapped = .true.
    endif

    if(is_overlapped) then
       if( overlap%pe == NULL_PE ) then
          overlap%pe = pe
       else
          if(overlap%pe .NE. pe) call mpp_error(FATAL,  &
               "mpp_domains_define.inc(insert_update_overlap): mismatch on pe")
       endif
       overlap%count = overlap%count+1
       count = overlap%count
       if(count > MAXOVERLAP) call mpp_error(FATAL, &
            "mpp_domains_define.inc(insert_update_overlap): number of overlap is greater than MAXOVERLAP, increase MAXOVERLAP")
       overlap%is(count) = is
       overlap%ie(count) = ie
       overlap%js(count) = js
       overlap%je(count) = je
       overlap%tileMe (count) = 1
       overlap%tileNbr(count) = 1
       overlap%dir(count) = dir
       if(is_reverse) then
          overlap%rotation(count) = ONE_HUNDRED_EIGHTY
       else
          overlap%rotation(count) = ZERO
       end if
    end if

  end subroutine insert_update_overlap

  !#####################################################################################
  subroutine insert_mosaic_update_overlap(overlap, pe, tileMe, tileNbr, is, ie, js, je, dir, &
                                          rotation, from_contact, refined)
    type(overlap_type), intent(inout) :: overlap
    integer,            intent(in   ) :: tileMe, tileNbr, pe
    integer,            intent(in   ) :: is, ie, js, je
    integer,            intent(in   ) :: dir, rotation
    logical,            intent(in   ) :: from_contact, refined
    integer                           :: count

    if( overlap%pe == NULL_PE ) then
       overlap%pe = pe
    else
       if(overlap%pe .NE. pe) call mpp_error(FATAL,  &
            "mpp_domains_define.inc(insert_mosaic_update_overlap): mismatch on pe")
    endif
    overlap%count = overlap%count+1
    count = overlap%count
    if(count > MAXOVERLAP) call mpp_error(FATAL, &
         "mpp_domains_define.inc(insert_mosaic_update_overlap): number of overlap is greater than MAXOVERLAP, increase MAXOVERLAP")
    overlap%tileMe      (count) = tileMe
    overlap%tileNbr     (count) = tileNbr
    overlap%is          (count) = is
    overlap%ie          (count) = ie
    overlap%js          (count) = js
    overlap%je          (count) = je
    overlap%dir         (count) = dir
    overlap%rotation    (count) = rotation
    overlap%from_contact(count) = from_contact
    overlap%is_refined  (count) = refined

  end subroutine insert_mosaic_update_overlap

!#######################################################################
subroutine deallocate_overlap_type( overlap)
  type(overlap_type), intent(inout) :: overlap

  if(overlap%count == 0) then
     if( .NOT. associated(overlap%tileMe)) return
  else
     if( .NOT. associated(overlap%tileMe)) call mpp_error(FATAL, &
          "deallocate_overlap_type(mpp_domains_define): overlap is not been allocated")
  endif
  if(ASSOCIATED(overlap%tileMe)) deallocate(overlap%tileMe)
  if(ASSOCIATED(overlap%tileNbr)) deallocate(overlap%tileNbr)
  if(ASSOCIATED(overlap%is)) deallocate(overlap%is)
  if(ASSOCIATED(overlap%ie)) deallocate(overlap%ie)
  if(ASSOCIATED(overlap%js)) deallocate(overlap%js)
  if(ASSOCIATED(overlap%je)) deallocate(overlap%je)
  if(ASSOCIATED(overlap%isMe)) deallocate(overlap%isMe)
  if(ASSOCIATED(overlap%ieMe)) deallocate(overlap%ieMe)
  if(ASSOCIATED(overlap%jsMe)) deallocate(overlap%jsMe)
  if(ASSOCIATED(overlap%jeMe)) deallocate(overlap%jeMe)
  if(ASSOCIATED(overlap%dir)) deallocate(overlap%dir)
  if(ASSOCIATED(overlap%rotation)) deallocate(overlap%rotation)
  if(ASSOCIATED(overlap%is_refined)) deallocate(overlap%is_refined)
  if(ASSOCIATED(overlap%index)) deallocate(overlap%index)
  if(ASSOCIATED(overlap%from_contact)) deallocate(overlap%from_contact)
  overlap%count = 0

end subroutine deallocate_overlap_type

!#######################################################################
subroutine deallocate_overlapSpec(overlap)
type(overlapSpec), intent(inout) :: overlap
integer                          :: n

   if(ASSOCIATED(overlap%send)) then
      do n = 1, size(overlap%send(:))
         call deallocate_overlap_type(overlap%send(n))
      enddo
      deallocate(overlap%send)
   endif
   if(ASSOCIATED(overlap%recv)) then
      do n = 1, size(overlap%recv(:))
         call deallocate_overlap_type(overlap%recv(n))
      enddo
      deallocate(overlap%recv)
   endif


end subroutine deallocate_overlapSpec

!#######################################################################
!--- this routine add the overlap_in into overlap_out
subroutine add_update_overlap( overlap_out, overlap_in)
  type(overlap_type), intent(inout) :: overlap_out
  type(overlap_type), intent(in   ) :: overlap_in
  type(overlap_type)                :: overlap
  integer                          :: count, count_in, count_out

  ! if overlap_out%count == 0, then just copy overlap_in to overlap_out
  count_in  = overlap_in %count
  count_out = overlap_out%count
  count     = count_in+count_out
  if(count_in == 0) call mpp_error(FATAL, &
       "mpp_domains_define.inc(add_update_overlap): overlap_in%count is zero")

  if(count_out == 0) then
     if(associated(overlap_out%tileMe)) call mpp_error(FATAL, &
          "mpp_domains_define.inc(add_update_overlap): overlap is already been allocated but count=0")
     call allocate_update_overlap(overlap_out, count_in)
     overlap_out%pe  = overlap_in%pe
  else  ! need to expand the dimension size of overlap
     if(overlap_in%pe .NE. overlap_out%pe) call mpp_error(FATAL, &
          "mpp_domains_define.inc(add_update_overlap): mismatch of pe between overlap_in and overlap_out")

     call allocate_update_overlap(overlap, count_out)
     overlap%tileMe      (1:count_out) = overlap_out%tileMe      (1:count_out)
     overlap%tileNbr     (1:count_out) = overlap_out%tileNbr     (1:count_out)
     overlap%is          (1:count_out) = overlap_out%is          (1:count_out)
     overlap%ie          (1:count_out) = overlap_out%ie          (1:count_out)
     overlap%js          (1:count_out) = overlap_out%js          (1:count_out)
     overlap%je          (1:count_out) = overlap_out%je          (1:count_out)
     overlap%isMe        (1:count_out) = overlap_out%isMe        (1:count_out)
     overlap%ieMe        (1:count_out) = overlap_out%ieMe        (1:count_out)
     overlap%jsMe        (1:count_out) = overlap_out%jsMe        (1:count_out)
     overlap%jeMe        (1:count_out) = overlap_out%jeMe        (1:count_out)
     overlap%dir         (1:count_out) = overlap_out%dir         (1:count_out)
     overlap%rotation    (1:count_out) = overlap_out%rotation    (1:count_out)
     overlap%is_refined  (1:count_out) = overlap_out%is_refined  (1:count_out)
     overlap%index       (1:count_out) = overlap_out%index       (1:count_out)
     overlap%from_contact(1:count_out) = overlap_out%from_contact(1:count_out)
     call deallocate_overlap_type(overlap_out)
     call allocate_update_overlap(overlap_out, count)
     overlap_out%tileMe      (1:count_out) = overlap%tileMe      (1:count_out)
     overlap_out%tileNbr     (1:count_out) = overlap%tileNbr     (1:count_out)
     overlap_out%is          (1:count_out) = overlap%is          (1:count_out)
     overlap_out%ie          (1:count_out) = overlap%ie          (1:count_out)
     overlap_out%js          (1:count_out) = overlap%js          (1:count_out)
     overlap_out%je          (1:count_out) = overlap%je          (1:count_out)
     overlap_out%isMe        (1:count_out) = overlap%isMe        (1:count_out)
     overlap_out%ieMe        (1:count_out) = overlap%ieMe        (1:count_out)
     overlap_out%jsMe        (1:count_out) = overlap%jsMe        (1:count_out)
     overlap_out%jeMe        (1:count_out) = overlap%jeMe        (1:count_out)
     overlap_out%dir         (1:count_out) = overlap%dir         (1:count_out)
     overlap_out%rotation    (1:count_out) = overlap%rotation    (1:count_out)
     overlap_out%is_refined  (1:count_out) = overlap%is_refined  (1:count_out)
     overlap_out%index       (1:count_out) = overlap%index       (1:count_out)
     overlap_out%from_contact(1:count_out) = overlap%from_contact(1:count_out)        
     call deallocate_overlap_type(overlap)
  end if
  overlap_out%count                           = count
  overlap_out%tileMe      (count_out+1:count) = overlap_in%tileMe      (1:count_in)
  overlap_out%tileNbr     (count_out+1:count) = overlap_in%tileNbr     (1:count_in)
  overlap_out%is          (count_out+1:count) = overlap_in%is          (1:count_in)
  overlap_out%ie          (count_out+1:count) = overlap_in%ie          (1:count_in)
  overlap_out%js          (count_out+1:count) = overlap_in%js          (1:count_in)
  overlap_out%je          (count_out+1:count) = overlap_in%je          (1:count_in)
  overlap_out%isMe        (count_out+1:count) = overlap_in%isMe        (1:count_in)
  overlap_out%ieMe        (count_out+1:count) = overlap_in%ieMe        (1:count_in)
  overlap_out%jsMe        (count_out+1:count) = overlap_in%jsMe        (1:count_in)
  overlap_out%jeMe        (count_out+1:count) = overlap_in%jeMe        (1:count_in)
  overlap_out%dir         (count_out+1:count) = overlap_in%dir         (1:count_in)
  overlap_out%rotation    (count_out+1:count) = overlap_in%rotation    (1:count_in)
  overlap_out%is_refined  (count_out+1:count) = overlap_in%is_refined  (1:count_in)
  overlap_out%index       (count_out+1:count) = overlap_in%index       (1:count_in)
  overlap_out%from_contact(count_out+1:count) = overlap_in%from_contact(1:count_in)

end subroutine add_update_overlap

!##############################################################################
subroutine expand_update_overlap_list(overlapList, npes)
  type(overlap_type), pointer :: overlapList(:)
  integer,      intent(in   ) :: npes
  type(overlap_type), pointer,save           :: newlist(:) => NULL()
  integer                                    :: nlist_old, nlist, m

  nlist_old = size(overlaplist(:))
  if(nlist_old .GE. npes) call mpp_error(FATAL, &
     'mpp_domains_define.inc(expand_update_overlap_list): size of overlaplist should be smaller than npes')
  nlist = min(npes, 2*nlist_old)
  allocate(newlist(nlist))
  do m = 1, nlist_old
     call add_update_overlap(newlist(m), overlaplist(m))
     call deallocate_overlap_type(overlapList(m))
  enddo

  deallocate(overlapList)
  overlaplist => newlist
  newlist => NULL()

  return

end subroutine expand_update_overlap_list

!##################################################################################
subroutine expand_check_overlap_list(overlaplist, npes)
  type(overlap_type), pointer :: overlaplist(:)
  integer,         intent(in) :: npes
  type(overlap_type), pointer,save           :: newlist(:) => NULL()
  integer                                    :: nlist_old, nlist, m

  nlist_old = size(overlaplist(:))
  if(nlist_old .GE. npes) call mpp_error(FATAL, &
     'mpp_domains_define.inc(expand_check_overlap_list): size of overlaplist should be smaller than npes')
  nlist = min(npes, 2*nlist_old)
  allocate(newlist(nlist))
  do m = 1,size(overlaplist(:))
     call add_check_overlap(newlist(m), overlaplist(m))
     call deallocate_overlap_type(overlapList(m))
  enddo
  deallocate(overlapList)
  overlaplist => newlist


  return

end subroutine expand_check_overlap_list


!###############################################################################
subroutine check_overlap_pe_order(domain, overlap, name)
  type(domain2d),    intent(in) :: domain
  type(overlapSpec), intent(in) :: overlap
  character(len=*),  intent(in) :: name
  integer                       :: nlist, m
  integer                       :: pos1, pos2

  nlist = size(domain%list(:))

  do m = 2, overlap%nsend
     pos1 = overlap%send(m-1)%pe - domain%pe
     if(pos1 .LT.0) pos1 = pos1 + nlist
     pos2 = overlap%send(m)%pe - domain%pe
     if(pos2 .LT.0) pos2 = pos2 + nlist
     if(pos2 .LE. pos1) then
        print*, trim(name)//" at pe = ", domain%pe, &
                ": send pe is ", overlap%send(m-1)%pe, overlap%send(m)%pe
        call mpp_error(FATAL, &
             "mpp_domains_define.inc(check_overlap_pe_order): pe is not in right order for send")
     endif
  enddo


  do m = 2, overlap%nrecv
     pos1 = overlap%recv(m-1)%pe - domain%pe
     if(pos1 .LE.0) pos1 = pos1 + nlist
     pos2 = overlap%recv(m)%pe - domain%pe
     if(pos2 .LE.0) pos2 = pos2 + nlist
     if(pos2 .GE. pos1) then
        print*, trim(name)//" at pe = ", domain%pe, ": recv pe is ", overlap%recv(m-1)%pe, overlap%recv(m)%pe
        call mpp_error(FATAL, &
             "mpp_domains_define.inc(check_overlap_pe_order): pe is not in right order for recv")
     endif
  enddo


end subroutine check_overlap_pe_order
