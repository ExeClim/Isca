! -*-f90-*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                 MPP_DOMAINS: initialization and termination                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <SUBROUTINE NAME="mpp_domains_init">
!  <OVERVIEW>
!    Initialize domain decomp package.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Called to initialize the <TT>mpp_domains_mod</TT> package.
!    
!    <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to have
!    <TT>mpp_domains_mod</TT> keep you informed of what it's up
!    to. <TT>MPP_DEBUG</TT> returns even more information for debugging.
!    
!    <TT>mpp_domains_init</TT> will call <TT>mpp_init</TT>, to make sure
!    <LINK SRC="mpp.html"><TT>mpp_mod</TT></LINK> is initialized. (Repeated
!    calls to <TT>mpp_init</TT> do no harm, so don't worry if you already
!    called it).
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_init(flags)
!  </TEMPLATE>
!  <IN NAME="flags" TYPE="integer"></IN>
! </SUBROUTINE>
    subroutine mpp_domains_init(flags)
      integer, intent(in), optional :: flags
#ifdef use_libSMA
      integer                       :: l=0
#endif
      integer                       :: unit_begin, unit_end, unit_nml, io_status, unit
      logical                       :: opened
#ifdef use_libSMA
      integer(POINTER_KIND) :: ptr_info_var(16)
      pointer( ptr_info, ptr_info_var )
#endif 

      if( module_is_initialized )return
      call mpp_init(flags)           !this is a no-op if already initialized
      call mpp_pset_init             !this is a no-op if already initialized
      module_is_initialized = .TRUE.
      pe = mpp_root_pe()
      unit = stdlog()
      if( mpp_pe() .EQ.mpp_root_pe() ) write( unit,'(/a)' )'MPP_DOMAINS module '//trim(version)//trim(tagname)

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
          domain_clocks_on = flags.EQ.MPP_DOMAIN_TIME
      end if

      !--- namelist 
      unit_begin = 103
      unit_end   = 512
      do unit_nml = unit_begin, unit_end
         inquire( unit_nml,OPENED=opened )
         if( .NOT.opened )exit
      end do

      open(unit_nml,file='input.nml', iostat=io_status)
      read(unit_nml,mpp_domains_nml,iostat=io_status)
      close(unit_nml)

      select case(lowercase(trim(debug_update_domain)))
      case("none")
         debug_update_level = NO_CHECK
      case("fatal")
         debug_update_level = FATAL
      case("warning")
         debug_update_level = WARNING
      case("note")
         debug_update_level = NOTe
      case default
         call mpp_error(FATAL, "mpp_domains_init: debug_update_level should be 'none', 'fatal', 'warning', or 'note'") 
      end select

      call mpp_domains_set_stack_size(32768) !default, pretty arbitrary
#ifdef use_libSMA
      call mpp_malloc( ptr_info, 16, l )
#endif 

!NULL_DOMAIN is a domaintype that can be used to initialize to undef
      call mpp_define_null_domain(NULL_DOMAIN1d);
      call mpp_define_null_domain(NULL_DOMAIN2d);

      if( domain_clocks_on )then
          pack_clock      = mpp_clock_id( 'Halo pack' )
          pack_loop_clock = mpp_clock_id( 'Halo pack loop' )
          send_clock      = mpp_clock_id( 'Halo send' )
          recv_clock      = mpp_clock_id( 'Halo recv' )
          unpk_clock      = mpp_clock_id( 'Halo unpk' )
          wait_clock      = mpp_clock_id( 'Halo wait' )
      end if
      return
    end subroutine mpp_domains_init

!#####################################################################
! <SUBROUTINE NAME="mpp_domains_exit">
!  <OVERVIEW>
!    Exit <TT>mpp_domains_mod</TT>.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Serves no particular purpose, but is provided should you require to
!    re-initialize <TT>mpp_domains_mod</TT>, for some odd reason.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_exit()
!  </TEMPLATE>
! </SUBROUTINE>
    subroutine mpp_domains_exit()
      integer :: unit 
      if( .NOT.module_is_initialized )return
      call mpp_max(mpp_domains_stack_hwm)
      unit = stdout()
      if( mpp_pe().EQ.mpp_root_pe() )write( unit,* )'MPP_DOMAINS_STACK high water mark=', mpp_domains_stack_hwm
      module_is_initialized = .FALSE.
      return
    end subroutine mpp_domains_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_CHECK_FIELD: Check parallel                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <SUBROUTINE NAME="mpp_check_field_3D" INTERFACE="mpp_check_field">
!   <IN NAME="field_in" TYPE="real, dimension(:,:,:)" > </IN>
!   <IN NAME="pelist1, pelist2" TYPE="integer, dimension(:)" > </IN>
!   <IN NAME="domain" TYPE="type(domain2d)" > </IN>
!   <IN NAME="mesg" TYPE="character(len=*)" > </IN>
!   <IN NAME="w_halo, s_halo, e_halo, n_halo" TYPE="integer, optional" > </IN>
!   <IN NAME="force_abort" TYPE="logical,optional" > </IN>
! </SUBROUTINE>

  subroutine mpp_check_field_3D(field_in, pelist1, pelist2, domain, mesg, &
                                 w_halo, s_halo, e_halo, n_halo, force_abort, position  )
!  This routine is used to do parallel checking for 3d data between n and m pe. The comparison is
!  is done on pelist2. When size of pelist2 is 1, we can check the halo; otherwise,
!  halo can not be checked.

  real, dimension(:,:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d),        intent(in) :: domain   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found      
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when true, call mpp_error if any difference
                                                    ! found. default value is false.  
  integer, intent(in), optional     :: position    ! when domain is symmetry, only value = CENTER is 
                                                   ! implemented. 
            
    integer :: k
    character(len=256) :: temp_mesg


    do k = 1, size(field_in,3)
       write(temp_mesg, '(a, i3)') trim(mesg)//" at level " , k
       call mpp_check_field_2d(field_in(:,:,k), pelist1, pelist2, domain, temp_mesg, &
                                 w_halo, s_halo, e_halo, n_halo, force_abort, position )
    enddo
    
  end subroutine mpp_check_field_3D
  
  
!#####################################################################################
! <SUBROUTINE NAME="mpp_check_field_2D" INTERFACE="mpp_check_field">
!   <IN NAME="field_in" TYPE="real, dimension(:,:)" > </IN>
! </SUBROUTINE>
  subroutine mpp_check_field_2d(field_in, pelist1, pelist2, domain, mesg, &
                                 w_halo, s_halo, e_halo, n_halo,force_abort, position  )
!  This routine is used to do parallel checking for 2d data between n and m pe. The comparison is
!  is done on pelist2. When size of pelist2 is 1, we can check the halo; otherwise,
!  halo can not be checked.

  real, dimension(:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d),        intent(in) :: domain   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.
  integer, intent(in), optional     :: position    ! when domain is symmetry, only value = CENTER is 
                                                   ! implemented. 

     if(present(position)) then
        if(position .NE. CENTER .AND. domain%symmetry) call mpp_error(FATAL,  &
           'mpp_check_field: when domain is symmetry, only value CENTER is implemented, contact author')
     endif

     if(size(pelist2(:)) == 1) then                 
        call mpp_check_field_2d_type1(field_in, pelist1, pelist2, domain, mesg, &
                                      w_halo, s_halo, e_halo, n_halo, force_abort )
     else if(size(pelist1(:)) == 1) then
        call mpp_check_field_2d_type1(field_in, pelist2, pelist1, domain, mesg, &
                                      w_halo, s_halo, e_halo, n_halo, force_abort )
     else if(size(pelist1(:)) .gt. 1 .and. size(pelist2(:)) .gt. 1) then
        call mpp_check_field_2d_type2(field_in, pelist1, pelist2, domain, mesg, force_abort )
     else
        call mpp_error(FATAL, 'mpp_check_field: size of both pelists should be greater than 0')
     endif

  end subroutine mpp_check_field_2D


!####################################################################################

  subroutine mpp_check_field_2d_type1(field_in, pelist1, pelist2, domain, mesg, &
                                 w_halo, s_halo, e_halo, n_halo,force_abort  )
!  This routine is used to check field between running on 1 pe (pelist2) and
!  n pe(pelist1). The need_to_be_checked data is sent to the pelist2 and All the
!  comparison is done on pelist2.

  real, dimension(:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d),        intent(in) :: domain   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found      
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.         
! some local data

  integer                :: pe,npes, p
  integer                :: hwest, hsouth, heast, hnorth, isg, ieg, jsg, jeg, xhalo, yhalo
  integer                :: i,j,im,jm,l,is,ie,js,je,isc,iec,jsc,jec,isd,ied,jsd,jed
  real,dimension(:,:), allocatable :: field1,field2
  real,dimension(:),   allocatable :: send_buffer
  integer, dimension(4)  ::  ibounds
  logical                :: check_success,  error_exit

  check_success = .TRUE.
  error_exit    = .FALSE.
  if(present(force_abort)) error_exit = force_abort
  hwest  = 0; if(present(w_halo)) hwest  = w_halo
  heast  = 0; if(present(e_halo)) heast  = e_halo
  hsouth = 0; if(present(s_halo)) hsouth = s_halo
  hnorth = 0; if(present(n_halo)) hnorth = n_halo

  pe = mpp_pe ()
  npes = mpp_npes()

  call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
  xhalo = isc - isd
  yhalo = jsc - jsd
  !--- need to checked halo size should not be bigger than x_halo or y_halo
  if(hwest .gt. xhalo .or. heast .gt. xhalo .or. hsouth .gt. yhalo .or. hnorth .gt. yhalo) &
     call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': The halo size is not correct')

  is = isc - hwest; ie = iec + heast; js = jsc - hsouth; je = jec + hnorth
  allocate(field2(is:ie,js:je))

  ! check if the field_in is on compute domain or data domain
  if((size(field_in,1) .eq. iec-isc+1) .and. (size(field_in,2) .eq. jec-jsc+1)) then
     !if field_in on compute domain, you can not check halo points
     if( hwest .ne. 0 .or. heast .ne. 0 .or. hsouth .ne. 0 .or. hnorth .ne. 0 ) &
         call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': field is on compute domain, can not check halo')
     field2(:,:) = field_in(:,:)
  else if((size(field_in,1) .eq. ied-isd+1) .and. (size(field_in,2) .eq. jed-jsd+1)) then
     field2(is:ie,js:je) = field_in(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1)
  else if((size(field_in,1) .eq. ieg-isg+1) .and. (size(field_in,2) .eq. jeg-jsg+1)) then
     if( hwest .ne. 0 .or. heast .ne. 0 .or. hsouth .ne. 0 .or. hnorth .ne. 0 ) &
         call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': field is on compute domain, can not check halo')
     field2(is:ie,js:je) = field_in(1:ie-is+1,1:je-js+1)
  else if((size(field_in,1) .eq. ieg-isg+1+2*xhalo) .and. (size(field_in,2) .eq. jeg-jsg+1+2*yhalo)) then
     field2(is:ie,js:je) = field_in(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1)
  else
     print*, 'on pe ', pe, 'domain: ', isc, iec, jsc, jec, isd, ied, jsd, jed, 'size of field: ', size(field_in,1), size(field_in,2)
     call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//':field is not on compute, data or global domain')
  endif
     
  call mpp_sync_self()
  
  if(any(pelist1 == pe)) then  ! send data to root pe
  
     im = ie-is+1; jm=je-js+1
     allocate(send_buffer(im*jm))
     
     ibounds(1) = is; ibounds(2) = ie; ibounds(3) = js; ibounds(4) = je
     l = 0
     do i = is,ie
     do j = js,je
        l = l+1
        send_buffer(l) = field2(i,j)
     enddo
     enddo
!  send the check bounds and data to the root pe
     ! Force use of "scalar", integer pointer mpp interface
     call mpp_send(ibounds(1), plen=4, to_pe=pelist2(1))
     call mpp_send(send_buffer(1),plen=im*jm, to_pe=pelist2(1))
     deallocate(send_buffer)
     
   else if(pelist2(1) == pe) then        ! receive data and compare
     do p = pelist1(1), pelist1(size(pelist1(:)))
     ! Force use of "scalar", integer pointer mpp interface
        call mpp_recv(ibounds(1), glen=4,from_pe=p)
        is = ibounds(1); ie = ibounds(2); js=ibounds(3); je=ibounds(4)
        im = ie-is+1; jm=je-js+1
        if(allocated(field1)) deallocate(field1)
        if(allocated(send_buffer)) deallocate(send_buffer)
        allocate(field1(is:ie,js:je),send_buffer(im*jm))
     ! Force use of "scalar", integer pointer mpp interface
        call mpp_recv(send_buffer(1),glen=im*jm,from_pe=p)
        l = 0
        
!  compare here, the comparison criteria can be changed according to need
        do i = is,ie
        do j = js,je
           l = l+1
           field1(i,j) = send_buffer(l)
           if(field1(i,j) .ne. field2(i,j)) then
         !   write to standard output
             print*,trim(mesg)//": ", i, j, field1(i,j), field2(i,j), field1(i,j) - field2(i,j)
!             write(stdout(),'(a,2i,2f)') trim(mesg), i, j, pass_field(i,j), field_check(i,j)
             check_success = .FALSE.
             if(error_exit) call mpp_error(FATAL,"mpp_check_field: can not reproduce at this point")
           endif
        enddo
        enddo
      enddo
        
      if(check_success) then
         print*, trim(mesg)//": ", 'comparison between 1 pe and ', npes-1, ' pes is ok'
      endif
  ! release memery
      deallocate(field1, send_buffer)
    endif
      
    deallocate(field2)
    
    call mpp_sync()
    
  end subroutine mpp_check_field_2d_type1
  
!####################################################################

  subroutine mpp_check_field_2d_type2(field_in, pelist1, pelist2, domain, mesg,force_abort)
!  This routine is used to check field between running on m pe (root pe) and
!  n pe. This routine can not check halo.
 
  real, dimension(:,:),  intent(in) :: field_in
  type(domain2d),        intent(in) :: domain
  integer, dimension(:), intent(in) :: pelist1
  integer, dimension(:), intent(in) :: pelist2
  character(len=*),      intent(in) :: mesg
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.
! some local variables                              
  logical                :: check_success, error_exit
  real, dimension(:,:), allocatable :: field1, field2
  integer :: i, j, pe, npes, isd,ied,jsd,jed, is, ie, js, je
  type(domain2d) :: domain1, domain2
  
    check_success = .TRUE.
    error_exit    = .FALSE.
    if(present(force_abort)) error_exit = force_abort
    pe = mpp_pe()
    npes = mpp_npes()
    call mpp_sync_self()
    if(any(pelist1 == pe)) domain1 = domain
    if(any(pelist2 == pe)) domain2 = domain
    
!  Comparison is made on pelist2.
    if(any(pelist2 == pe)) then
       call mpp_get_data_domain(domain2, isd, ied, jsd, jed)
       call mpp_get_compute_domain(domain2, is, ie, js, je)
       allocate(field1(isd:ied, jsd:jed),field2(isd:ied, jsd:jed))
      if((size(field_in,1) .ne. ied-isd+1) .or. (size(field_in,2) .ne. jed-jsd+1)) &
         call mpp_error(FATAL,'mpp_check_field: input field is not on the data domain')
      field2(isd:ied, jsd:jed) = field_in(:,:)
    endif
      
!  broadcast domain
    call mpp_broadcast_domain(domain1)
    call mpp_broadcast_domain(domain2)
    
    call mpp_redistribute(domain1,field_in,domain2,field1)
    
    if(any(pelist2 == pe)) then
        do i =is,ie
        do j =js,je
          if(field1(i,j) .ne. field2(i,j)) then
             print*, trim(mesg)//": ", i, j, field1(i,j), field2(i,j), field1(i,j) - field2(i,j)
!             write(stdout(),'(a,2i,2f)') trim(mesg), i, j, field_check(i,j), field_out(i,j)
             check_success = .FALSE.
             if(error_exit) call mpp_error(FATAL,"mpp_check_field: can not reproduce at this point")
          endif
        enddo
        enddo
        if(check_success) &
             print*, trim(mesg)//": ", 'comparison between ', size(pelist1(:)), ' pes and ', &
                   size(pelist2(:)), ' pe on', pe, ' pes is ok'
    endif          
                   
    if(any(pelist2 == pe))    deallocate(field1, field2)
    
    call mpp_sync()
    
    return
    
  end subroutine mpp_check_field_2d_type2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_BROADCAST_DOMAIN                                           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_domain( domain )
!broadcast domain (useful only outside the context of its own pelist)
      type(domain2D), intent(inout) :: domain
      integer, allocatable :: pes(:)
      logical :: native         !true if I'm on the pelist of this domain
      integer :: listsize, listpos
      integer :: n
      integer, dimension(5) :: msg, info         !pe and compute domain of each item in list
      
      if( .NOT.module_is_initialized ) &
                 call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN: You must first call mpp_domains_init.' ) 
      
!get the current pelist
      allocate( pes(0:mpp_npes()-1) )
      call mpp_get_current_pelist(pes)
      
!am I part of this domain?
      native = ASSOCIATED(domain%list)
      
!set local list size
      if( native )then
          listsize = size(domain%list(:))
      else
          listsize = 0
      end if
      call mpp_max(listsize)
      
      if( .NOT.native )then
!initialize domain%list and set null values in message
          allocate( domain%list(0:listsize-1) )
          domain%pe = NULL_PE
          domain%pos = -1
          allocate(domain%x(1), domain%y(1))
          do n = 0, listsize-1
             allocate(domain%list(n)%x(1), domain%list(n)%y(1) )
          end do
          domain%x%compute%begin =  1
          domain%x%compute%end   = -1
          domain%y%compute%begin =  1
          domain%y%compute%end   = -1
      end if
!initialize values in info
      info(1) = domain%pe
      call mpp_get_compute_domain( domain, info(2), info(3), info(4), info(5) )
      
!broadcast your info across current pelist and unpack if needed
      listpos = 0
      do n = 0,mpp_npes()-1
         msg = info
         if( mpp_pe().EQ.pes(n) .AND. debug )write( stderr(),* )'PE ', mpp_pe(), 'broadcasting msg ', msg
         call mpp_broadcast( msg, 5, pes(n) )
!no need to unpack message if native
!no need to unpack message from non-native PE
         if( .NOT.native .AND. msg(1).NE.NULL_PE )then
             domain%list(listpos)%pe = msg(1)
             domain%list(listpos)%x%compute%begin = msg(2)
             domain%list(listpos)%x%compute%end   = msg(3)
             domain%list(listpos)%y%compute%begin = msg(4)
             domain%list(listpos)%y%compute%end   = msg(5)
             listpos = listpos + 1
             if( debug )write( stderr(),* )'PE ', mpp_pe(), 'received domain from PE ', msg(1), 'is,ie,js,je=', msg(2:5)
         end if
      end do 
    end subroutine mpp_broadcast_domain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_UPDATE_DOMAINS: fill halos for 2D decomposition            !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef VECTOR_FIELD_
#define VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_UPDATE_DOMAINS_2D_
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r8_2D
#undef MPP_UPDATE_DOMAINS_3D_
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r8_3D
#undef MPP_UPDATE_DOMAINS_4D_
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r8_4D
#undef MPP_UPDATE_DOMAINS_5D_
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r8_5D
#ifdef  VECTOR_FIELD_
#undef MPP_UPDATE_DOMAINS_2D_V_
#define MPP_UPDATE_DOMAINS_2D_V_ mpp_update_domain2D_r8_2Dv
#undef MPP_UPDATE_DOMAINS_3D_V_
#define MPP_UPDATE_DOMAINS_3D_V_ mpp_update_domain2D_r8_3Dv
#undef MPP_UPDATE_DOMAINS_4D_V_
#define MPP_UPDATE_DOMAINS_4D_V_ mpp_update_domain2D_r8_4Dv
#undef MPP_UPDATE_DOMAINS_5D_V_
#define MPP_UPDATE_DOMAINS_5D_V_ mpp_update_domain2D_r8_5Dv
#endif
#undef MPP_REDISTRIBUTE_2D_
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_r8_2D
#undef MPP_REDISTRIBUTE_3D_
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_r8_3D
#undef MPP_REDISTRIBUTE_4D_
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_r8_4D
#undef MPP_REDISTRIBUTE_5D_
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_r8_5D
#include <mpp_update_domains2D.h>
#undef  VECTOR_FIELD_

#ifdef OVERLOAD_C8
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_UPDATE_DOMAINS_2D_
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c8_2D
#undef MPP_UPDATE_DOMAINS_3D_
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c8_3D
#undef MPP_UPDATE_DOMAINS_4D_
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c8_4D
#undef MPP_UPDATE_DOMAINS_5D_
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c8_5D
#undef MPP_REDISTRIBUTE_2D_
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_c8_2D
#undef MPP_REDISTRIBUTE_3D_
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_c8_3D
#undef MPP_REDISTRIBUTE_4D_
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_c8_4D
#undef MPP_REDISTRIBUTE_5D_
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_c8_5D
#include <mpp_update_domains2D.h>
#endif

#ifndef no_8byte_integers
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_UPDATE_DOMAINS_2D_
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i8_2D
#undef MPP_UPDATE_DOMAINS_3D_
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i8_3D
#undef MPP_UPDATE_DOMAINS_4D_
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i8_4D
#undef MPP_UPDATE_DOMAINS_5D_
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i8_5D
#undef MPP_REDISTRIBUTE_2D_
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_i8_2D
#undef MPP_REDISTRIBUTE_3D_
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_i8_3D
#undef MPP_REDISTRIBUTE_4D_
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_i8_4D
#undef MPP_REDISTRIBUTE_5D_
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_i8_5D
#include <mpp_update_domains2D.h>
#endif

#ifdef OVERLOAD_R4
#undef VECTOR_FIELD_
#define VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_UPDATE_DOMAINS_2D_
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r4_2D
#undef MPP_UPDATE_DOMAINS_3D_
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r4_3D
#undef MPP_UPDATE_DOMAINS_4D_
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r4_4D
#undef MPP_UPDATE_DOMAINS_5D_
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r4_5D
#ifdef  VECTOR_FIELD_
#undef MPP_UPDATE_DOMAINS_2D_V_
#define MPP_UPDATE_DOMAINS_2D_V_ mpp_update_domain2D_r4_2Dv
#undef MPP_UPDATE_DOMAINS_3D_V_
#define MPP_UPDATE_DOMAINS_3D_V_ mpp_update_domain2D_r4_3Dv
#undef MPP_UPDATE_DOMAINS_4D_V_
#define MPP_UPDATE_DOMAINS_4D_V_ mpp_update_domain2D_r4_4Dv
#undef MPP_UPDATE_DOMAINS_5D_V_
#define MPP_UPDATE_DOMAINS_5D_V_ mpp_update_domain2D_r4_5Dv
#endif
#undef MPP_REDISTRIBUTE_2D_
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_r4_2D
#undef MPP_REDISTRIBUTE_3D_
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_r4_3D
#undef MPP_REDISTRIBUTE_4D_
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_r4_4D
#undef MPP_REDISTRIBUTE_5D_
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_r4_5D
#include <mpp_update_domains2D.h>
#undef  VECTOR_FIELD_
#endif

#ifdef OVERLOAD_C4
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_UPDATE_DOMAINS_2D_
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c4_2D
#undef MPP_UPDATE_DOMAINS_3D_
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c4_3D
#undef MPP_UPDATE_DOMAINS_4D_
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c4_4D
#undef MPP_UPDATE_DOMAINS_5D_
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c4_5D
#undef MPP_REDISTRIBUTE_2D_
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_c4_2D
#undef MPP_REDISTRIBUTE_3D_
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_c4_3D
#undef MPP_REDISTRIBUTE_4D_
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_c4_4D
#undef MPP_REDISTRIBUTE_5D_
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_c4_5D
#include <mpp_update_domains2D.h>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_UPDATE_DOMAINS_2D_
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i4_2D
#undef MPP_UPDATE_DOMAINS_3D_
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i4_3D
#undef MPP_UPDATE_DOMAINS_4D_
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i4_4D
#undef MPP_UPDATE_DOMAINS_5D_
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i4_5D
#undef MPP_REDISTRIBUTE_2D_
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_i4_2D
#undef MPP_REDISTRIBUTE_3D_
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_i4_3D
#undef MPP_REDISTRIBUTE_4D_
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_i4_4D
#undef MPP_REDISTRIBUTE_5D_
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_i4_5D
#include <mpp_update_domains2D.h>

!*******************************************************
#undef VECTOR_FIELD_
#define VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_DO_UPDATE_3D_
#define MPP_DO_UPDATE_3D_ mpp_do_update_r8_3d
#ifdef  VECTOR_FIELD_
#undef MPP_DO_UPDATE_3D_V_
#define MPP_DO_UPDATE_3D_V_ mpp_do_update_r8_3dv
#endif
#include <mpp_do_update.h>
#include <mpp_do_updateV.h>

#ifdef OVERLOAD_C8
#undef  VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_DO_UPDATE_3D_
#define MPP_DO_UPDATE_3D_ mpp_do_update_c8_3d
#include <mpp_do_update.h>
#define VECTOR_FIELD_
#endif

#ifndef no_8byte_integers
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_DO_UPDATE_3D_
#define MPP_DO_UPDATE_3D_ mpp_do_update_i8_3d
#include <mpp_do_update.h>
#endif

#ifdef OVERLOAD_R4
#undef VECTOR_FIELD_
#define VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_DO_UPDATE_3D_
#define MPP_DO_UPDATE_3D_ mpp_do_update_r4_3d
#ifdef  VECTOR_FIELD_
#undef MPP_DO_UPDATE_3D_V_
#define MPP_DO_UPDATE_3D_V_ mpp_do_update_r4_3dv
#endif
#include <mpp_do_update.h>
#include <mpp_do_updateV.h>
#endif

#ifdef OVERLOAD_C4
#undef  VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_DO_UPDATE_3D_
#define MPP_DO_UPDATE_3D_ mpp_do_update_c4_3d
#include <mpp_do_update.h>
#define VECTOR_FIELD_
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_DO_UPDATE_3D_
#define MPP_DO_UPDATE_3D_ mpp_do_update_i4_3d
#include <mpp_do_update.h>


#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_DO_CHECK_3D_
#define MPP_DO_CHECK_3D_ mpp_do_check_r8_3d
#ifdef  VECTOR_FIELD_
#undef MPP_DO_CHECK_3D_V_
#define MPP_DO_CHECK_3D_V_ mpp_do_check_r8_3dv
#endif
#include <mpp_do_check.h>
#include <mpp_do_checkV.h>

#ifdef OVERLOAD_C8
#undef  VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_DO_CHECK_3D_
#define MPP_DO_CHECK_3D_ mpp_do_check_c8_3d
#include <mpp_do_check.h>
#define VECTOR_FIELD_
#endif

#ifndef no_8byte_integers
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_DO_CHECK_3D_
#define MPP_DO_CHECK_3D_ mpp_do_check_i8_3d
#include <mpp_do_check.h>
#endif

#ifdef OVERLOAD_R4
#undef VECTOR_FIELD_
#define VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_DO_CHECK_3D_
#define MPP_DO_CHECK_3D_ mpp_do_check_r4_3d
#ifdef  VECTOR_FIELD_
#undef MPP_DO_CHECK_3D_V_
#define MPP_DO_CHECK_3D_V_ mpp_do_check_r4_3dv
#endif
#include <mpp_do_check.h>
#include <mpp_do_checkV.h>
#endif

#ifdef OVERLOAD_C4
#undef  VECTOR_FIELD_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_DO_CHECK_3D_
#define MPP_DO_CHECK_3D_ mpp_do_check_c4_3d
#include <mpp_do_check.h>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_DO_CHECK_3D_
#define MPP_DO_CHECK_3D_ mpp_do_check_i4_3d
#include <mpp_do_check.h>


!bnc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_UPDATE_DOMAINS_AD: adjoint fill halos for 2D decomposition !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$#undef VECTOR_FIELD_
!!$#define VECTOR_FIELD_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(DOUBLE_KIND)
!!$#undef MPP_UPDATE_DOMAINS_AD_2D_
!!$#define MPP_UPDATE_DOMAINS_AD_2D_ mpp_update_domain2D_ad_r8_2D
!!$#undef MPP_UPDATE_DOMAINS_AD_3D_
!!$#define MPP_UPDATE_DOMAINS_AD_3D_ mpp_update_domain2D_ad_r8_3D
!!$#undef MPP_UPDATE_DOMAINS_AD_4D_
!!$#define MPP_UPDATE_DOMAINS_AD_4D_ mpp_update_domain2D_ad_r8_4D
!!$#undef MPP_UPDATE_DOMAINS_AD_5D_
!!$#define MPP_UPDATE_DOMAINS_AD_5D_ mpp_update_domain2D_ad_r8_5D
!!$#ifdef  VECTOR_FIELD_
!!$#undef MPP_UPDATE_DOMAINS_AD_2D_V_
!!$#define MPP_UPDATE_DOMAINS_AD_2D_V_ mpp_update_domain2D_ad_r8_2Dv
!!$#undef MPP_UPDATE_DOMAINS_AD_3D_V_
!!$#define MPP_UPDATE_DOMAINS_AD_3D_V_ mpp_update_domain2D_ad_r8_3Dv
!!$#undef MPP_UPDATE_DOMAINS_AD_4D_V_
!!$#define MPP_UPDATE_DOMAINS_AD_4D_V_ mpp_update_domain2D_ad_r8_4Dv
!!$#undef MPP_UPDATE_DOMAINS_AD_5D_V_
!!$#define MPP_UPDATE_DOMAINS_AD_5D_V_ mpp_update_domain2D_ad_r8_5Dv
!!$#endif
!!$#include <mpp_update_domains2D_ad.h>
!!$#undef  VECTOR_FIELD_
!!$
!!$#ifdef OVERLOAD_C8
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(DOUBLE_KIND)
!!$#undef MPP_UPDATE_DOMAINS_AD_2D_
!!$#define MPP_UPDATE_DOMAINS_AD_2D_ mpp_update_domain2D_ad_c8_2D
!!$#undef MPP_UPDATE_DOMAINS_AD_3D_
!!$#define MPP_UPDATE_DOMAINS_AD_3D_ mpp_update_domain2D_ad_c8_3D
!!$#undef MPP_UPDATE_DOMAINS_AD_4D_
!!$#define MPP_UPDATE_DOMAINS_AD_4D_ mpp_update_domain2D_ad_c8_4D
!!$#undef MPP_UPDATE_DOMAINS_AD_5D_
!!$#define MPP_UPDATE_DOMAINS_AD_5D_ mpp_update_domain2D_ad_c8_5D
!!$#include <mpp_update_domains2D_ad.h>
!!$#endif
!!$
!!$#ifndef no_8byte_integers
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(LONG_KIND)
!!$#undef MPP_UPDATE_DOMAINS_AD_2D_
!!$#define MPP_UPDATE_DOMAINS_AD_2D_ mpp_update_domain2D_ad_i8_2D
!!$#undef MPP_UPDATE_DOMAINS_AD_3D_
!!$#define MPP_UPDATE_DOMAINS_AD_3D_ mpp_update_domain2D_ad_i8_3D
!!$#undef MPP_UPDATE_DOMAINS_AD_4D_
!!$#define MPP_UPDATE_DOMAINS_AD_4D_ mpp_update_domain2D_ad_i8_4D
!!$#undef MPP_UPDATE_DOMAINS_AD_5D_
!!$#define MPP_UPDATE_DOMAINS_AD_5D_ mpp_update_domain2D_ad_i8_5D
!!$#include <mpp_update_domains2D_ad.h>
!!$#endif
!!$
!!$#ifdef OVERLOAD_R4
!!$#undef VECTOR_FIELD_
!!$#define VECTOR_FIELD_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(FLOAT_KIND)
!!$#undef MPP_UPDATE_DOMAINS_AD_2D_
!!$#define MPP_UPDATE_DOMAINS_AD_2D_ mpp_update_domain2D_ad_r4_2D
!!$#undef MPP_UPDATE_DOMAINS_AD_3D_
!!$#define MPP_UPDATE_DOMAINS_AD_3D_ mpp_update_domain2D_ad_r4_3D
!!$#undef MPP_UPDATE_DOMAINS_AD_4D_
!!$#define MPP_UPDATE_DOMAINS_AD_4D_ mpp_update_domain2D_ad_r4_4D
!!$#undef MPP_UPDATE_DOMAINS_AD_5D_
!!$#define MPP_UPDATE_DOMAINS_AD_5D_ mpp_update_domain2D_ad_r4_5D
!!$#ifdef  VECTOR_FIELD_
!!$#undef MPP_UPDATE_DOMAINS_AD_2D_V_
!!$#define MPP_UPDATE_DOMAINS_AD_2D_V_ mpp_update_domain2D_ad_r4_2Dv
!!$#undef MPP_UPDATE_DOMAINS_AD_3D_V_
!!$#define MPP_UPDATE_DOMAINS_AD_3D_V_ mpp_update_domain2D_ad_r4_3Dv
!!$#undef MPP_UPDATE_DOMAINS_AD_4D_V_
!!$#define MPP_UPDATE_DOMAINS_AD_4D_V_ mpp_update_domain2D_ad_r4_4Dv
!!$#undef MPP_UPDATE_DOMAINS_AD_5D_V_
!!$#define MPP_UPDATE_DOMAINS_AD_5D_V_ mpp_update_domain2D_ad_r4_5Dv
!!$#endif
!!$#include <mpp_update_domains2D_ad.h>
!!$#endif
!!$
!!$#ifdef OVERLOAD_C4
!!$#undef  VECTOR_FIELD_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(FLOAT_KIND)
!!$#undef MPP_UPDATE_DOMAINS_AD_2D_
!!$#define MPP_UPDATE_DOMAINS_AD_2D_ mpp_update_domain2D_ad_c4_2D
!!$#undef MPP_UPDATE_DOMAINS_AD_3D_
!!$#define MPP_UPDATE_DOMAINS_AD_3D_ mpp_update_domain2D_ad_c4_3D
!!$#undef MPP_UPDATE_DOMAINS_AD_4D_
!!$#define MPP_UPDATE_DOMAINS_AD_4D_ mpp_update_domain2D_ad_c4_4D
!!$#undef MPP_UPDATE_DOMAINS_AD_5D_
!!$#define MPP_UPDATE_DOMAINS_AD_5D_ mpp_update_domain2D_ad_c4_5D
!!$#include <mpp_update_domains2D_ad.h>
!!$#endif
!!$
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(INT_KIND)
!!$#undef MPP_UPDATE_DOMAINS_AD_2D_
!!$#define MPP_UPDATE_DOMAINS_AD_2D_ mpp_update_domain2D_ad_i4_2D
!!$#undef MPP_UPDATE_DOMAINS_AD_3D_
!!$#define MPP_UPDATE_DOMAINS_AD_3D_ mpp_update_domain2D_ad_i4_3D
!!$#undef MPP_UPDATE_DOMAINS_AD_4D_
!!$#define MPP_UPDATE_DOMAINS_AD_4D_ mpp_update_domain2D_ad_i4_4D
!!$#undef MPP_UPDATE_DOMAINS_AD_5D_
!!$#define MPP_UPDATE_DOMAINS_AD_5D_ mpp_update_domain2D_ad_i4_5D
!!$#include <mpp_update_domains2D_ad.h>
!!$
!!$!*******************************************************
!!$#undef VECTOR_FIELD_
!!$#define VECTOR_FIELD_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(DOUBLE_KIND)
!!$#undef MPP_DO_UPDATE_AD_3D_
!!$#define MPP_DO_UPDATE_AD_3D_ mpp_do_update_ad_r8_3d
!!$#ifdef  VECTOR_FIELD_
!!$#undef MPP_DO_UPDATE_AD_3D_V_
!!$#define MPP_DO_UPDATE_AD_3D_V_ mpp_do_update_ad_r8_3dv
!!$#endif
!!$#include <mpp_do_update_ad.h>
!!$#include <mpp_do_updateV_ad.h>
!!$#undef  VECTOR_FIELD_
!!$
!!$#ifdef OVERLOAD_C8
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(DOUBLE_KIND)
!!$#undef MPP_DO_UPDATE_AD_3D_
!!$#define MPP_DO_UPDATE_AD_3D_ mpp_do_update_ad_c8_3d
!!$#include <mpp_do_update_ad.h>
!!$#endif
!!$                                                                                  
!!$#ifndef no_8byte_integers
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(LONG_KIND)
!!$#undef MPP_DO_UPDATE_AD_3D_
!!$#define MPP_DO_UPDATE_AD_3D_ mpp_do_update_ad_i8_3d
!!$#include <mpp_do_update_ad.h>
!!$#endif
!!$
!!$#ifdef OVERLOAD_R4
!!$#undef  VECTOR_FIELD_
!!$#define VECTOR_FIELD_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(FLOAT_KIND)
!!$#undef MPP_DO_UPDATE_AD_3D_
!!$#define MPP_DO_UPDATE_AD_3D_ mpp_do_update_ad_r4_3d
!!$#ifdef  VECTOR_FIELD_
!!$#undef MPP_DO_UPDATE_AD_3D_V_
!!$#define MPP_DO_UPDATE_AD_3D_V_ mpp_do_update_ad_r4_3dv
!!$#endif
!!$#include <mpp_do_update_ad.h>
!!$#include <mpp_do_updateV_ad.h>
!!$#endif
!!$
!!$#ifdef OVERLOAD_C4
!!$#undef  VECTOR_FIELD_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(FLOAT_KIND)
!!$#undef MPP_DO_UPDATE_AD_3D_
!!$#define MPP_DO_UPDATE_AD_3D_ mpp_do_update_ad_c4_3d
!!$#include <mpp_do_update_ad.h>
!!$#endif
!!$
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(INT_KIND)
!!$#undef MPP_DO_UPDATE_AD_3D_
!!$#define MPP_DO_UPDATE_AD_3D_ mpp_do_update_ad_i4_3d
!!$#include <mpp_do_update_ad.h>

!bnc


!********************************************************
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_DO_REDISTRIBUTE_3D_
#define MPP_DO_REDISTRIBUTE_3D_ mpp_do_redistribute_r8_3D
#include <mpp_do_redistribute.h>
#undef  VECTOR_FIELD_

#ifdef OVERLOAD_C8
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_DO_REDISTRIBUTE_3D_
#define MPP_DO_REDISTRIBUTE_3D_ mpp_do_redistribute_c8_3D
#include <mpp_do_redistribute.h>
#endif

#ifndef no_8byte_integers
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_DO_REDISTRIBUTE_3D_
#define MPP_DO_REDISTRIBUTE_3D_ mpp_do_redistribute_i8_3D
#include <mpp_do_redistribute.h>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(LONG_KIND)
#undef MPP_DO_REDISTRIBUTE_3D_
#define MPP_DO_REDISTRIBUTE_3D_ mpp_do_redistribute_l8_3D
#include <mpp_do_redistribute.h>
#endif

#ifdef OVERLOAD_R4
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_DO_REDISTRIBUTE_3D_
#define MPP_DO_REDISTRIBUTE_3D_ mpp_do_redistribute_r4_3D
#include <mpp_do_redistribute.h>
#undef  VECTOR_FIELD_
#endif

#ifdef OVERLOAD_C4
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_DO_REDISTRIBUTE_3D_
#define MPP_DO_REDISTRIBUTE_3D_ mpp_do_redistribute_c4_3D
#include <mpp_do_redistribute.h>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_DO_REDISTRIBUTE_3D_
#define MPP_DO_REDISTRIBUTE_3D_ mpp_do_redistribute_i4_3D
#include <mpp_do_redistribute.h>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(INT_KIND)
#undef MPP_DO_REDISTRIBUTE_3D_
#define MPP_DO_REDISTRIBUTE_3D_ mpp_do_redistribute_l4_3D
#include <mpp_do_redistribute.h>

#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_GET_BOUNDARY_2D_
#define MPP_GET_BOUNDARY_2D_ mpp_get_boundary_r8_2d
#undef MPP_GET_BOUNDARY_3D_
#define MPP_GET_BOUNDARY_3D_ mpp_get_boundary_r8_3d
#undef MPP_GET_BOUNDARY_4D_
#define MPP_GET_BOUNDARY_4D_ mpp_get_boundary_r8_4d
#undef MPP_GET_BOUNDARY_5D_
#define MPP_GET_BOUNDARY_5D_ mpp_get_boundary_r8_5d
#undef MPP_GET_BOUNDARY_2D_V_
#define MPP_GET_BOUNDARY_2D_V_ mpp_get_boundary_r8_2dv
#undef MPP_GET_BOUNDARY_3D_V_
#define MPP_GET_BOUNDARY_3D_V_ mpp_get_boundary_r8_3dv
#undef MPP_GET_BOUNDARY_4D_V_
#define MPP_GET_BOUNDARY_4D_V_ mpp_get_boundary_r8_4dv
#undef MPP_GET_BOUNDARY_5D_V_
#define MPP_GET_BOUNDARY_5D_V_ mpp_get_boundary_r8_5dv
#include <mpp_get_boundary.h>

#ifdef OVERLOAD_R4
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_GET_BOUNDARY_2D_
#define MPP_GET_BOUNDARY_2D_ mpp_get_boundary_r4_2d
#undef MPP_GET_BOUNDARY_3D_
#define MPP_GET_BOUNDARY_3D_ mpp_get_boundary_r4_3d
#undef MPP_GET_BOUNDARY_4D_
#define MPP_GET_BOUNDARY_4D_ mpp_get_boundary_r4_4d
#undef MPP_GET_BOUNDARY_5D_
#define MPP_GET_BOUNDARY_5D_ mpp_get_boundary_r4_5d
#undef MPP_GET_BOUNDARY_2D_V_
#define MPP_GET_BOUNDARY_2D_V_ mpp_get_boundary_r4_2dv
#undef MPP_GET_BOUNDARY_3D_V_
#define MPP_GET_BOUNDARY_3D_V_ mpp_get_boundary_r4_3dv
#undef MPP_GET_BOUNDARY_4D_V_
#define MPP_GET_BOUNDARY_4D_V_ mpp_get_boundary_r4_4dv
#undef MPP_GET_BOUNDARY_5D_V_
#define MPP_GET_BOUNDARY_5D_V_ mpp_get_boundary_r4_5dv
#include <mpp_get_boundary.h>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_DO_GET_BOUNDARY_3D_
#define MPP_DO_GET_BOUNDARY_3D_ mpp_do_get_boundary_r8_3d
#undef MPP_DO_GET_BOUNDARY_3DV_
#define MPP_DO_GET_BOUNDARY_3D_V_ mpp_do_get_boundary_r8_3dv
#include <mpp_do_get_boundary.h>

#ifdef OVERLOAD_R4
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_DO_GET_BOUNDARY_3D_
#define MPP_DO_GET_BOUNDARY_3D_ mpp_do_get_boundary_r4_3d
#undef MPP_DO_GET_BOUNDARY_3D_V_
#define MPP_DO_GET_BOUNDARY_3D_V_ mpp_do_get_boundary_r4_3dv
#include <mpp_do_get_boundary.h>
#endif
