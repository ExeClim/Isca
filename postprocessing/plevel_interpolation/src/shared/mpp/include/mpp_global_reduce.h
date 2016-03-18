  function MPP_GLOBAL_REDUCE_2D_( domain, field, locus, position )
    MPP_TYPE_ :: MPP_GLOBAL_REDUCE_2D_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,:)
    integer, intent(out), optional :: locus(2)
    integer, intent(in),  optional :: position
    MPP_TYPE_ :: field3D(size(field,1),size(field,2),1)
    integer :: locus3D(3)
    pointer( ptr, field3D )
    ptr = LOC(field)
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_2D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D, position )
        locus = locus3D(1:2)
    else
        MPP_GLOBAL_REDUCE_2D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, position = position )
    end if
    return
  end function MPP_GLOBAL_REDUCE_2D_

  function MPP_GLOBAL_REDUCE_3D_( domain, field, locus, position )
    MPP_TYPE_ :: MPP_GLOBAL_REDUCE_3D_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(0:,0:,:)
    integer, intent(out), optional :: locus(3)
    integer, intent(in),  optional :: position

    MPP_TYPE_ :: local
    integer, save :: l_locus(3)
    MPP_TYPE_, save :: g_val  ! need storage class w/ global address; not sure whether fn result has required class
    integer, save   :: here   ! need storage class w/ global address
    integer :: ioff, joff, isc, iec, jsc, jec, ishift, jshift

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_REDUCE: You must first call mpp_domains_init.' )

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec )
    call mpp_get_domain_shift(domain, ishift, jshift, position)
    iec = iec + ishift
    jec = jec + jshift

    if( size(field,1).EQ. iec-isc+1  .AND. size(field,2).EQ. jec-jsc+1 )then
!field is on compute domain
        ioff = isc
        joff = jsc
    else if( size(field,1).EQ.domain%x(1)%memory%size+ishift .AND. size(field,2).EQ.domain%y(1)%memory%size+jshift )then
!field is on data domain
        ioff = domain%x(1)%data%begin
        joff = domain%y(1)%data%begin
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_REDUCE_: incoming field array must match either compute domain or data domain.' )
    end if

!get your local max/min
    local = REDUCE_VAL_(field(isc-ioff: iec-ioff, jsc-joff: jec-joff,:))
!find the global
    g_val = local
    call MPP_REDUCE_( g_val, domain%list(:)%pe )
!find locus of the global max/min
    if( PRESENT(locus) )then
!which PE is it on? min of all the PEs that have it
        here = mpp_npes()+1
        if( g_val == local )here = pe
        call mpp_min( here, domain%list(:)%pe )
!find the locus here
        if( pe.EQ.here )l_locus = REDUCE_LOC_(field(isc-ioff: iec-ioff, jsc-joff: jec-joff,:))
        l_locus(1) = l_locus(1) + ioff
        l_locus(2) = l_locus(2) + joff
        call mpp_broadcast( l_locus, 3, here, domain%list(:)%pe )
        locus = l_locus
    end if
    MPP_GLOBAL_REDUCE_3D_ = g_val
    return
  end function MPP_GLOBAL_REDUCE_3D_

  function MPP_GLOBAL_REDUCE_4D_( domain, field, locus, position )
    MPP_TYPE_ :: MPP_GLOBAL_REDUCE_4D_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,:,:,:)
    integer, intent(out), optional :: locus(4)
    integer, intent(in),  optional :: position

    MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
    integer :: locus3D(3)
    pointer( ptr, field3D )
    ptr = LOC(field)
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_4D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D, position )
        locus(1:2) = locus3D(1:2)
        locus(3) = modulo(locus3D(3),size(field,3))
        locus(4) = (locus3D(3)-locus(3))/size(field,3) + 1
        if( locus(3).EQ.0 )then
            locus(3) = size(field,3)
            locus(4) = locus(4) - 1
        end if
    else
        MPP_GLOBAL_REDUCE_4D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, position = position  )
    end if
    return
  end function MPP_GLOBAL_REDUCE_4D_

  function MPP_GLOBAL_REDUCE_5D_( domain, field, locus, position )
    MPP_TYPE_ :: MPP_GLOBAL_REDUCE_5D_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,:,:,:,:)
    integer, intent(out), optional :: locus(5)
    integer, intent(in),  optional :: position

    MPP_TYPE_ :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4)*size(field,5))
    integer :: locus3D(3)
    pointer( ptr, field3D )
    ptr = LOC(field)
    if( PRESENT(locus) )then
        MPP_GLOBAL_REDUCE_5D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, locus3D, position )
        locus(1:2) = locus3D(1:2)
        locus(3) = modulo(locus3D(3),size(field,3))
        locus(4) = modulo(locus3D(3),size(field,3)*size(field,4))
        locus(5) = (locus3D(3)-locus(4))/size(field,3)/size(field,4) + 1
        if( locus(3).EQ.0 )then
            locus(3) = size(field,3)
            locus(4) = locus(4) - 1
        end if
    else
        MPP_GLOBAL_REDUCE_5D_ = MPP_GLOBAL_REDUCE_3D_( domain, field3D, position = position  )
    end if
    return
  end function MPP_GLOBAL_REDUCE_5D_
