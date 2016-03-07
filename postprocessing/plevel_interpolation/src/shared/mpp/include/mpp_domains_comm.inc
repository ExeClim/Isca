
! -*-f90-*-
! $Id: mpp_domains_comm.inc,v 18.0 2010/03/02 23:57:13 fms Exp $


    function mpp_redistribute_init_comm(domain_in,l_addrs_in, domain_out,l_addrs_out, &
                                        isize_in,jsize_in,ksize_in,isize_out,jsize_out,ksize_out) RESULT(d_comm)
      type(DomainCommunicator2D), pointer       :: d_comm
      type(domain2D),target,      intent(in)    :: domain_in
      integer(LONG_KIND),         intent(in)    :: l_addrs_in(:)
      type(domain2D),target,      intent(in)    :: domain_out
      integer(LONG_KIND),         intent(in)    :: l_addrs_out(:)
      integer,                    intent(in)    :: isize_in
      integer,                    intent(in)    :: jsize_in
      integer,                    intent(in)    :: ksize_in
      integer,                    intent(in)    :: isize_out 
      integer,                    intent(in)    :: jsize_out
      integer,                    intent(in)    :: ksize_out 

      integer(LONG_KIND) :: domain_id
      integer :: m, list
      integer :: is, ie, js, je, ke, ioff, joff, list_size
      integer :: isc, iec, jsc, jec
      integer :: lsize,rsize,msgsize,to_pe,from_pe
      integer,           allocatable,dimension(:) :: isL, jsL
      integer(LONG_KIND),allocatable,dimension(:,:) :: slist_addr
      character(len=8) :: text


    ! This test determines whether input fields are from allocated memory (LOC gets global
    ! address) or "static" memory (need shmem_ptr). This probably needs to be generalized
    ! to determine appropriate mechanism for each incoming address.

    ! "Concurrent" run mode may leave field_in or field_out unassociated if pe does not
    ! contain in/out data. Use of STATIC option for ocean complicates this as ocean component
    ! always defined. Field_out is always a boundary structure and so is always allocated or
    ! not depending on whether it's used. If field out is defined (>0), then it is used otherwise
    ! field in must be defined.

!fix ke
      ke = 0
      if( domain_in%pe /= NULL_PE )ke = ksize_in
      if( domain_out%pe /= NULL_PE )then
          if( ke /= 0 .AND. ke /= ksize_out ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: mismatch between field_in and field_out.' )
          ke = ksize_out
      end if
      if( ke == 0 )call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: either domain_in or domain_out must be native.' )
!check sizes
      if( domain_in%pe /= NULL_PE )then
          if( isize_in /= domain_in%x(1)%data%size .OR. jsize_in /= domain_in%y(1)%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: field_in must be on data domain of domain_in.' )
      end if
      if( domain_out%pe /= NULL_PE )then
          if( isize_out /= domain_out%x(1)%data%size .OR. jsize_out /= domain_out%y(1)%data%size ) &
               call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: field_out must be on data domain of domain_out.' )
      end if


    ! Create unique domain identifier
      list_size = size(l_addrs_in(:))
      if(l_addrs_out(1) > 0)then
        domain_id = set_domain_id(domain_out%id,ke+list_size)
      else
        domain_id = set_domain_id(domain_in%id,ke+list_size)
      endif

      d_comm =>get_comm(domain_id,l_addrs_in(1),l_addrs_out(1))

      if(d_comm%initialized)return   ! Found existing field/domain communicator

      d_comm%l_addr = l_addrs_in(1)
      d_comm%domain_in =>domain_in
      d_comm%Slist_size = size(domain_out%list(:))
      d_comm%isize_in = isize_in
      d_comm%jsize_in = jsize_in
      d_comm%ke = ke

!send 
      lsize = d_comm%Slist_size-1
      allocate(d_comm%sendis(1,0:lsize), d_comm%sendie(1,0:lsize), &
               d_comm%sendjs(1,0:lsize), d_comm%sendje(1,0:lsize), &
               d_comm%S_msize(0:lsize),isL(0:lsize),jsL(0:lsize))
      allocate(slist_addr(list_size,0:lsize))
      allocate(d_comm%cto_pe(0:lsize), d_comm%S_do_buf(0:lsize))

      isL=0;jsL=0
      slist_addr = -9999
      d_comm%cto_pe=-1
      d_comm%sendis=0; d_comm%sendie=0
      d_comm%sendjs=0; d_comm%sendje=0;
      d_comm%S_msize=0
      d_comm%S_do_buf=.false.

      ioff = domain_in%x(1)%data%begin
      joff = domain_in%y(1)%data%begin

      call mpp_get_compute_domain( domain_in, isc, iec, jsc, jec )
      do list = 0,lsize
         m = mod( domain_out%pos+list+lsize+1, lsize+1 )
         d_comm%cto_pe(list) = domain_out%list(m)%pe
         to_pe =  d_comm%cto_pe(list)
         is = domain_out%list(m)%x(1)%compute%begin
         ie = domain_out%list(m)%x(1)%compute%end
         js = domain_out%list(m)%y(1)%compute%begin
         je = domain_out%list(m)%y(1)%compute%end
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie >= is .AND. je >= js )then
           d_comm%S_do_buf(list) = .true.
           d_comm%sendis(1,list)=is; d_comm%sendie(1,list)=ie
           d_comm%sendjs(1,list)=js; d_comm%sendje(1,list)=je
           d_comm%S_msize(list) = (ie-is+1)*(je-js+1)*ke
           isL(list) = is-ioff+1; jsL(list) = js-joff+1
         end if
      end do 

      call mpp_sync_self()
!recv 
      d_comm%domain_out =>domain_out
      d_comm%Rlist_size = size(domain_in%list(:))
      d_comm%isize_out = isize_out
      d_comm%jsize_out = jsize_out

      rsize = d_comm%Rlist_size-1
      allocate(d_comm%recvis(1,0:rsize), d_comm%recvie(1,0:rsize), &
               d_comm%recvjs(1,0:rsize), d_comm%recvje(1,0:rsize), &
               d_comm%R_msize(0:rsize))
      allocate(d_comm%cfrom_pe(0:rsize), d_comm%R_do_buf(0:rsize))
      allocate(d_comm%isizeR(0:rsize), d_comm%jsizeR(0:rsize))
      allocate(d_comm%sendisR(1,0:rsize), d_comm%sendjsR(1,0:rsize))
      allocate(d_comm%rem_addrl(list_size,0:rsize))
      d_comm%rem_addrl=-9999
      d_comm%cfrom_pe=-1
      d_comm%recvis=0; d_comm%recvie=0
      d_comm%recvjs=0; d_comm%recvje=0;
      d_comm%R_msize=0
      d_comm%R_do_buf=.false.
      d_comm%isizeR=0; d_comm%jsizeR=0
      d_comm%sendisR=0; d_comm%sendjsR=0

      call mpp_get_compute_domain( domain_out, isc, iec, jsc, jec )
      do list = 0,rsize
         m = mod( domain_in%pos+rsize+1-list, rsize+1 )
         d_comm%cfrom_pe(list) = domain_in%list(m)%pe
         from_pe = d_comm%cfrom_pe(list)
         is = domain_in%list(m)%x(1)%compute%begin
         ie = domain_in%list(m)%x(1)%compute%end
         js = domain_in%list(m)%y(1)%compute%begin
         je = domain_in%list(m)%y(1)%compute%end
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie >= is .AND. je >= js )then
           d_comm%R_do_buf(list) = .true.
           d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
           d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
           d_comm%R_msize(list) = (ie-is+1)*(je-js+1)*ke
         end if
      end do 

      d_comm%isize_max = isize_in; call mpp_max(d_comm%isize_max)
      d_comm%jsize_max = jsize_in; call mpp_max(d_comm%jsize_max)

    ! Handles case where S_msize and/or R_msize are 0 size array
      msgsize = ( MAXVAL( (/0,sum(d_comm%S_msize(:))/) ) + MAXVAL( (/0,sum(d_comm%R_msize(:))/) ) ) * list_size
      if(msgsize>0)then
         mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, msgsize )
         if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
             write( text,'(i8)' )mpp_domains_stack_hwm
             call mpp_error( FATAL, 'MPP_REDISTRIBUTE_INIT_COMM: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                  //trim(text)//') from all PEs.' )
         end if
      end if

      DEALLOCATE(slist_addr,isL,jsL)

      d_comm%initialized = .true.

    end function mpp_redistribute_init_comm


    function mpp_global_field_init_comm(domain,l_addr,isize_g,jsize_g,isize_l, &
                                        jsize_l, ksize,l_addr2,flags, position) RESULT(d_comm)
      type(DomainCommunicator2D), pointer       :: d_comm
      type(domain2D),target,      intent(in)    :: domain
      integer(LONG_KIND),         intent(in)    :: l_addr
      integer,                    intent(in)    :: isize_g
      integer,                    intent(in)    :: jsize_g
      integer,                    intent(in)    :: isize_l
      integer,                    intent(in)    :: jsize_l 
      integer,                    intent(in)    :: ksize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2
      integer, optional,          intent(in)    :: flags
      integer, optional,          intent(in)    :: position

      integer(LONG_KIND) :: domain_id
      integer :: n, lpos, rpos, list, nlist, tile_id
      integer :: update_flags
      logical :: xonly, yonly
      integer :: is, ie, js, je, ioff, joff, ishift, jshift
      integer :: lsize,msgsize,from_pe
      integer,           allocatable,dimension(:) :: isL, jsL
      integer(LONG_KIND),allocatable,dimension(:,:) :: slist_addr
      integer(LONG_KIND),save       ,dimension(2)   :: rem_addr
      character(len=8) :: text

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: must first call mpp_domains_init.' )
      update_flags=XUPDATE+YUPDATE; xonly = .FALSE.;  yonly = .FALSE.
      if( PRESENT(flags) )then
          update_flags = flags
          xonly = BTEST(flags,EAST)
          yonly = BTEST(flags,SOUTH)
          if( .NOT.xonly .AND. .NOT.yonly )call mpp_error( WARNING,  &
                    'MPP_GLOBAL_FIELD: you must have flags=XUPDATE, YUPDATE or XUPDATE+YUPDATE' )
         if(xonly .AND. yonly) then
            xonly = .false.; yonly = .false.
         endif
      end if

      call mpp_get_domain_shift(domain, ishift, jshift, position=position)

      if( isize_g /= (domain%x(1)%global%size+ishift) .OR. jsize_g /= (domain%y(1)%global%size+jshift) ) &
           call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_INIT_COMM: incoming arrays do not match domain.' )

      if( isize_l == (domain%x(1)%compute%size+ishift) .AND. jsize_l == (domain%y(1)%compute%size+jshift) )then
!local is on compute domain                                                                      
        ioff = -domain%x(1)%compute%begin + 1                                                       
        joff = -domain%y(1)%compute%begin + 1
      elseif( isize_l == (domain%x(1)%memory%size+ishift) .AND. jsize_l == (domain%y(1)%memory%size+jshift) )then
!local is on data domain                                                                        
        ioff = -domain%x(1)%data%begin + 1                                                         
        joff = -domain%y(1)%data%begin + 1
      else
        call mpp_error(FATAL,'MPP_GLOBAL_FIELD_INIT_COMM: incoming field array must match either compute domain or data domain.')
      endif

    ! Create unique domain identifier
      domain_id=set_domain_id(domain%id,ksize,update_flags, position=position)
      d_comm =>get_comm(domain_id,l_addr,l_addr2)
      
      if(d_comm%initialized)return   ! Found existing field/domain communicator
      d_comm%domain => domain
      d_comm%isize_in = isize_l; d_comm%isize_out = isize_g
      d_comm%jsize_in = jsize_l; d_comm%jsize_out = jsize_g
      d_comm%ke = ksize
      d_comm%gf_ioff=ioff; d_comm%gf_joff=joff

!fill off-domains (note loops begin at an offset of 1)
      if( xonly )then
        lsize = size(domain%x(1)%list(:))
!send
        allocate(d_comm%cto_pe(0:lsize-1))
        d_comm%cto_pe=-1
        do list = 0,lsize-1
          lpos = mod(domain%x(1)%pos+lsize-list,lsize)
          d_comm%cto_pe(list) = domain%x(1)%list(lpos)%pe
        end do
!recv
        allocate(d_comm%cfrom_pe(0:lsize-1))
        allocate(d_comm%recvis(1,0:lsize-1), d_comm%recvie(1,0:lsize-1), &
                 d_comm%recvjs(1,0:lsize-1), d_comm%recvje(1,0:lsize-1), &
                 d_comm%R_msize(0:lsize-1))                          
        d_comm%cfrom_pe=-1
        d_comm%recvis=0; d_comm%recvie=0
        d_comm%recvjs=0; d_comm%recvje=0;
        d_comm%R_msize=0
        do list = 0,lsize-1
          rpos = mod(domain%x(1)%pos+list,lsize)
          from_pe = domain%x(1)%list(rpos)%pe
          d_comm%cfrom_pe(list) = from_pe
          is = domain%list(from_pe)%x(1)%compute%begin; ie = domain%list(from_pe)%x(1)%compute%end+ishift
          js = domain%y(1)%compute%begin; je = domain%y(1)%compute%end+jshift
          d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
          d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
          d_comm%R_msize(list) = (ie-is+1) * (je-js+1) * ksize
        end do

      elseif( yonly )then
        lsize =  size(domain%y(1)%list(:))
!send  
        allocate(d_comm%cto_pe(0:lsize-1))
        d_comm%cto_pe=-1
        do list = 0,lsize
          lpos = mod(domain%y(1)%pos+lsize-list,lsize)
          d_comm%cto_pe(list) = domain%y(1)%list(lpos)%pe
        end do
!recv    
        allocate(d_comm%cfrom_pe(0:lsize-1))
        allocate(d_comm%recvis(1,0:lsize-1), d_comm%recvie(1,0:lsize-1), &
                 d_comm%recvjs(1,0:lsize-1), d_comm%recvje(1,0:lsize-1), &
                 d_comm%R_msize(0:lsize-1))                          
        d_comm%cfrom_pe=-1
        d_comm%recvis=0; d_comm%recvie=0
        d_comm%recvjs=0; d_comm%recvje=0;
        d_comm%R_msize=0
        do list = 0,lsize-1
          rpos = mod(domain%y(1)%pos+list,lsize)
          from_pe = domain%y(1)%list(rpos)%pe
          d_comm%cfrom_pe(list) = from_pe
          is = domain%x(1)%compute%begin; ie = domain%x(1)%compute%end+ishift
          js = domain%list(from_pe)%y(1)%compute%begin; je = domain%list(from_pe)%y(1)%compute%end+jshift
          d_comm%recvis(1,list)=is; d_comm%recvie(1,list)=ie
          d_comm%recvjs(1,list)=js; d_comm%recvje(1,list)=je
          d_comm%R_msize(list) = (ie-is+1) * (je-js+1) * ksize
        end do

      else
         nlist = size(domain%list(:))  
         tile_id = domain%tile_id(1)

         lsize = 0
         do list = 0,nlist-1
            if( domain%list(list)%tile_id(1) .NE. tile_id ) cycle
            lsize = lsize+1
         end do

         !send
         allocate(d_comm%cto_pe(0:lsize-1))
         d_comm%cto_pe=-1
         n = 0
         do list = 0,nlist-1
            lpos = mod(domain%pos+nlist-list,nlist)
            if( domain%list(lpos)%tile_id(1) .NE. tile_id ) cycle
            d_comm%cto_pe(n) = domain%list(lpos)%pe
            n = n + 1
         end do
         !recv
         allocate(d_comm%cfrom_pe(0:lsize-1))
         allocate(d_comm%recvis(1,0:lsize-1), d_comm%recvie(1,0:lsize-1), &
              d_comm%recvjs(1,0:lsize-1), d_comm%recvje(1,0:lsize-1), &
              d_comm%R_msize(0:lsize-1))
         d_comm%cfrom_pe=-1
         d_comm%recvis=0; d_comm%recvie=0
         d_comm%recvjs=0; d_comm%recvje=0;
         d_comm%R_msize=0
         n = 0
         do list = 0,nlist-1
            rpos = mod(domain%pos+list,nlist)
            if( domain%list(rpos)%tile_id(1) .NE. tile_id ) cycle
            d_comm%cfrom_pe(n) = domain%list(rpos)%pe
            is = domain%list(rpos)%x(1)%compute%begin; ie = domain%list(rpos)%x(1)%compute%end+ishift
            js = domain%list(rpos)%y(1)%compute%begin; je = domain%list(rpos)%y(1)%compute%end+jshift
            d_comm%recvis(1,n)=is; d_comm%recvie(1,n)=ie
            d_comm%recvjs(1,n)=js; d_comm%recvje(1,n)=je
            d_comm%R_msize(n) = (je-js+1) * (ie-is+1) * ksize
            n = n+1
         end do

      endif

      d_comm%Slist_size = lsize
      d_comm%Rlist_size = lsize

!send 
      allocate(d_comm%sendis(1,0:lsize-1), d_comm%sendie(1,0:lsize-1), &
               d_comm%sendjs(1,0:lsize-1), d_comm%sendje(1,0:lsize-1), &
               d_comm%S_msize(0:lsize-1),isL(0:lsize-1),jsL(0:lsize-1))
      allocate(slist_addr(2,0:lsize-1))
      isL=0; jsL=0
      slist_addr = -9999
      d_comm%sendis=0; d_comm%sendie=0
      d_comm%sendjs=0; d_comm%sendje=0;
      d_comm%S_msize=0
      do list = 0,lsize-1
        is=domain%x(1)%compute%begin; ie=domain%x(1)%compute%end+ishift
        js=domain%y(1)%compute%begin; je=domain%y(1)%compute%end+jshift
        d_comm%sendis(1,list)=is; d_comm%sendie(1,list)=ie
        d_comm%sendjs(1,list)=js; d_comm%sendje(1,list)=je
        d_comm%S_msize(list) = (je-js+1) * (ie-is+1) * ksize
        isL(list) = ioff+domain%x(1)%compute%begin; jsL(list) = joff+domain%y(1)%compute%begin
      end do

!recv
      allocate(d_comm%isizeR(0:lsize-1), d_comm%jsizeR(0:lsize-1))
      allocate(d_comm%sendisR(1,0:lsize-1), d_comm%sendjsR(1,0:lsize-1))
      if(.not.PRESENT(l_addr2))then
        allocate(d_comm%rem_addr(0:lsize-1))
        d_comm%rem_addr=-9999
      else
        allocate(d_comm%rem_addrx(0:lsize-1),d_comm%rem_addry(0:lsize-1))
        d_comm%rem_addrx=-9999; d_comm%rem_addry=-9999
      endif
      d_comm%isizeR=0; d_comm%jsizeR=0
      d_comm%sendisR=0; d_comm%sendjsR=0
      rem_addr = -9999

    ! Handles case where S_msize and/or R_msize are 0 size array
      msgsize = MAXVAL( (/0,sum(d_comm%S_msize(:))/) ) + MAXVAL( (/0,sum(d_comm%R_msize(:))/) )
      if(msgsize>0)then
         mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, msgsize )
         if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
             write( text,'(i8)' )mpp_domains_stack_hwm
             call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_INIT_COMM: mpp_domains_stack overflow, call mpp_domains_set_stack_size(' &
                  //trim(text)//') from all PEs.' )
         end if
      end if
            
      DEALLOCATE(slist_addr,isL,jsL)

      d_comm%initialized = .true.

    end function mpp_global_field_init_comm


    subroutine mpp_redistribute_free_comm(domain_in,l_addr,domain_out,l_addr2,ksize,lsize)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      type(domain2D),             intent(in)    :: domain_in
      integer(LONG_KIND),         intent(in)    :: l_addr
      type(domain2D),             intent(in)    :: domain_out  
      integer(LONG_KIND),         intent(in)    :: l_addr2
      integer,                    intent(in)    :: ksize,lsize

      integer(LONG_KIND) :: domain_id

      if(l_addr2 > 0)then
        domain_id = set_domain_id(domain_out%id,ksize+lsize)
      else
        domain_id = set_domain_id(domain_in%id,ksize+lsize)
      endif
      call free_comm(domain_id,l_addr,l_addr2)
    end subroutine mpp_redistribute_free_comm


    subroutine mpp_global_field_free_comm(domain,l_addr,ksize,l_addr2,flags)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      type(domain2D),             intent(in)    :: domain
      integer(LONG_KIND),         intent(in)    :: l_addr
      integer,                    intent(in)    :: ksize
      integer(LONG_KIND),optional,intent(in)    :: l_addr2
      integer,           optional,intent(in)    :: flags

      integer :: update_flags
      integer(LONG_KIND) :: domain_id

      update_flags=0; if(PRESENT(flags))update_flags=flags
      domain_id=set_domain_id(domain%id,ksize,update_flags)
      call free_comm(domain_id,l_addr,l_addr2)
    end subroutine mpp_global_field_free_comm


    subroutine free_comm(domain_id,l_addr,l_addr2)
    ! Since initialization of the d_comm type is expensive, freeing should be a rare
    ! event. Thus no attempt is made to salvage freed d_comm's.
      integer(LONG_KIND),         intent(in) :: domain_id
      integer(LONG_KIND),         intent(in) :: l_addr
      integer(LONG_KIND),optional,intent(in) :: l_addr2

      integer(LONG_KIND) :: dc_key,a_key
      integer :: dc_idx,a_idx,i_idx,insert,insert_a,insert_i
      integer :: a2_idx,insert_a2


      i_idx = find_key(domain_id,ids_sorted(1:n_ids),insert_i)
      a_idx = find_key(l_addr,addrs_sorted(1:a_sort_len),insert_a)
      a_key = int(addrs_idx(a_idx),KIND(LONG_KIND))
      if(PRESENT(l_addr2))then
        a2_idx = find_key(l_addr2,addrs2_sorted(1:a2_sort_len),insert_a2)
        a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(LONG_KIND))
      endif
      dc_key = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
      dc_idx = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)

      if(dc_idx < 0)then
        call mpp_error(FATAL,'FREE_COMM: attempt to remove nonexistent domains communicator key')
      endif
      call deallocate_comm(d_comm(dc_idx))
      call pop_key(dcKey_sorted,d_comm_idx,dc_sort_len,dc_idx)
      call pop_key(addrs_sorted,addrs_idx,a_sort_len,a_idx)
      if(PRESENT(l_addr2))call pop_key(addrs2_sorted,addrs2_idx,a2_sort_len,a2_idx)
    end subroutine free_comm


    function get_comm(domain_id,l_addr,l_addr2)
      integer(LONG_KIND),intent(in)       :: domain_id
      integer(LONG_KIND),intent(in)       :: l_addr
      integer(LONG_KIND),intent(in),optional :: l_addr2
      type(DomainCommunicator2D), pointer :: get_comm

      integer(LONG_KIND) :: dc_key,a_key
      integer :: i,dc_idx,a_idx,i_idx,insert,insert_a,insert_i
      integer :: a2_idx,insert_a2

      if(.not.ALLOCATED(d_comm))ALLOCATE(d_comm(MAX_FIELDS))
      i_idx   = find_key(domain_id,ids_sorted(1:n_ids),insert_i)
      a_idx = find_key(l_addr,addrs_sorted(1:a_sort_len),insert_a)
      a_key = int(addrs_idx(a_idx),KIND(LONG_KIND))
      if(PRESENT(l_addr2))then
        a2_idx = find_key(l_addr2,addrs2_sorted(1:a2_sort_len),insert_a2)
        a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(LONG_KIND))
      endif
      dc_key   = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
      dc_idx   = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)
      if(dc_idx > 0)then
        get_comm =>d_comm(d_comm_idx(dc_idx))
      else
        if(i_idx<0)then
          if(n_ids == MAX_DOM_IDS)then
            call mpp_error(FATAL,'GET_COMM: Maximum number of domains exceeded')
          endif
          n_ids = n_ids+1
          i_idx = push_key(ids_sorted,ids_idx,i_sort_len,insert_i,domain_id,n_ids)
        endif
        if(a_idx<0)then
          if(n_addrs == MAX_ADDRS)then
             call mpp_error(FATAL,'GET_COMM: Maximum number of memory addresses exceeded')
          endif
          n_addrs = n_addrs + 1
          a_idx = push_key(addrs_sorted,addrs_idx,a_sort_len,insert_a,l_addr,n_addrs)
        endif
        if(PRESENT(l_addr2))then
          if(a2_idx<0)then
            if(n_addrs2 == MAX_ADDRS2)then
               call mpp_error(FATAL,'GET_COMM: Maximum number of 2nd memory addresses exceeded')
            endif
            n_addrs2 = n_addrs2 + 1
            a2_idx = push_key(addrs2_sorted,addrs2_idx,a2_sort_len,insert_a2,l_addr2,n_addrs2)
          endif
        endif
        if(n_comm == MAX_FIELDS)then
          call mpp_error(FATAL,'GET_COMM: Maximum number of fields exceeded')
        endif
        a_key = int(addrs_idx(a_idx),KIND(8))
        if(PRESENT(l_addr2))a_key = a_key + ADDR2_BASE*int(addrs2_idx(a2_idx),KIND(8))
        dc_key   = DOMAIN_ID_BASE*int(ids_idx(i_idx),KIND(LONG_KIND)) + a_key
        dc_idx   = find_key(dc_key,dcKey_sorted(1:dc_sort_len),insert)
        if(dc_idx /= -1)call mpp_error(FATAL,'GET_COMM: attempt to insert existing key')
        n_comm = n_comm + 1
        i = push_key(dcKey_sorted,d_comm_idx,dc_sort_len,insert,dc_key,n_comm)
        d_comm_idx(insert) = n_comm
        if(PRESENT(l_addr2))then
          d_comm(n_comm)%l_addrx = l_addr
          d_comm(n_comm)%l_addry = l_addr2
        else
          d_comm(n_comm)%l_addr = l_addr
        endif
        get_comm =>d_comm(n_comm)
      endif
    end function get_comm


    function push_key(sorted,idx,n_idx,insert,key,ival)
      integer(LONG_KIND),intent(inout),dimension(:)   :: sorted
      integer,           intent(inout),dimension(-1:) :: idx  ! Start -1 to simplify first call logic in get_comm
      integer,           intent(inout)                :: n_idx
      integer,           intent(in)                   :: insert
      integer(LONG_KIND),intent(in)                   :: key
      integer,           intent(in)                   :: ival

      integer                                         :: push_key,i

      do i=n_idx,insert,-1
        sorted(i+1) = sorted(i)
        idx(i+1)    = idx(i)
      end do
      sorted(insert) = key
      n_idx = n_idx + 1
      idx(insert) = ival
      push_key = insert
    end function push_key


    subroutine pop_key(sorted,idx,n_idx,key_idx)
      integer(LONG_KIND),intent(inout),dimension(:)   :: sorted
      integer,           intent(inout),dimension(-1:) :: idx  ! Start -1 to simplify first call logic in get_comm
      integer,           intent(inout)                :: n_idx
      integer,           intent(in)                   :: key_idx

      integer                                         :: i

      do i=key_idx,n_idx-1
        sorted(i) = sorted(i+1)
        idx(i)    = idx(i+1)  
      end do                  
      sorted(n_idx) = -9999 
      idx(n_idx) = -9999
      n_idx = n_idx - 1  
    end subroutine pop_key


    function find_key(key,sorted,insert) RESULT(n)
    ! The algorithm used here requires monotonic keys w/out repetition.
      integer(LONG_KIND),intent(in)              :: key        ! new address to be found in list
      integer(LONG_KIND),dimension(:),intent(in) :: sorted  ! list of sorted local addrs
      integer,                        intent(out) :: insert
      integer :: n, n_max, n_min, n_key
      logical :: not_found

      n_key = size(sorted(:))
      insert = 1
      n = -1  ! value not in list
      if(n_key == 0)return  ! first call
      
      if(key < sorted(1))then
        insert = 1; return
      elseif(key > sorted(n_key))then
        insert = n_key+1; return
      endif
        
      if(key == sorted(1))then
        n = 1; return
      elseif(key == sorted(n_key))then
        n = n_key; return
      endif
        
      not_found = .true.
      n = n_key/2 + 1
      n_min=1; n_max=n_key
      do while(not_found)
        if(key == sorted(n))then
          not_found = .false.
        elseif(key > sorted(n))then
          if(key < sorted(n+1))then
            insert = n+1; exit
          endif
          n_min = n
          n = (n+1+n_max)/2
        else
          if(key > sorted(n-1))then
            insert = n; exit
          endif
          n_max = n
          n = (n+n_min)/2
        endif
        if(n==1 .or. n==n_key)exit
      end do
      if(not_found)n = -1  ! value not in list
    end function find_key


    subroutine deallocate_comm(d_comm)
      type(DomainCommunicator2D), intent(inout) :: d_comm

      d_comm%domain  =>NULL()
      d_comm%domain_in  =>NULL()
      d_comm%domain_out =>NULL()

      d_comm%initialized=.false.
      d_comm%id=-9999
      d_comm%l_addr  =-9999
      d_comm%l_addrx =-9999
      d_comm%l_addry =-9999

      if( _ALLOCATED(d_comm%sendis) )   DEALLOCATE(d_comm%sendis);    !!d_comm%sendis =>NULL()
      if( _ALLOCATED(d_comm%sendie) )   DEALLOCATE(d_comm%sendie);    !!d_comm%sendie =>NULL()
      if( _ALLOCATED(d_comm%sendjs) )   DEALLOCATE(d_comm%sendjs);    !!d_comm%sendjs =>NULL()
      if( _ALLOCATED(d_comm%sendje) )   DEALLOCATE(d_comm%sendje);    !!d_comm%sendje =>NULL()
      if( _ALLOCATED(d_comm%S_msize) )  DEALLOCATE(d_comm%S_msize);   !!d_comm%S_msize =>NULL()
      if( _ALLOCATED(d_comm%S_do_buf) ) DEALLOCATE(d_comm%S_do_buf);  !!d_comm%S_do_buf =>NULL()
      if( _ALLOCATED(d_comm%cto_pe) )   DEALLOCATE(d_comm%cto_pe);    !!d_comm%cto_pe  =>NULL()
      if( _ALLOCATED(d_comm%recvis) )   DEALLOCATE(d_comm%recvis);    !!d_comm%recvis =>NULL()
      if( _ALLOCATED(d_comm%recvie) )   DEALLOCATE(d_comm%recvie);    !!d_comm%recvie =>NULL()
      if( _ALLOCATED(d_comm%recvjs) )   DEALLOCATE(d_comm%recvjs);    !!d_comm%recvjs =>NULL()
      if( _ALLOCATED(d_comm%recvje) )   DEALLOCATE(d_comm%recvje);    !!d_comm%recvje =>NULL()
      if( _ALLOCATED(d_comm%R_msize) )  DEALLOCATE(d_comm%R_msize);   !!d_comm%R_msize =>NULL()
      if( _ALLOCATED(d_comm%R_do_buf) ) DEALLOCATE(d_comm%R_do_buf);  !!d_comm%R_do_buf =>NULL()
      if( _ALLOCATED(d_comm%cfrom_pe) ) DEALLOCATE(d_comm%cfrom_pe);  !!d_comm%cfrom_pe  =>NULL()
      d_comm%Slist_size=0; d_comm%Rlist_size=0
      d_comm%isize=0; d_comm%jsize=0; d_comm%ke=0
      d_comm%isize_in=0; d_comm%jsize_in=0
      d_comm%isize_out=0; d_comm%jsize_out=0
      d_comm%isize_max=0; d_comm%jsize_max=0
      d_comm%gf_ioff=0; d_comm%gf_joff=0
    ! Remote data
      if( _ALLOCATED(d_comm%isizeR) )   DEALLOCATE(d_comm%isizeR);    !!dd_comm%isizeR =>NULL()
      if( _ALLOCATED(d_comm%jsizeR) )   DEALLOCATE(d_comm%jsizeR);    !!dd_comm%jsizeR =>NULL()
      if( _ALLOCATED(d_comm%sendisR) )  DEALLOCATE(d_comm%sendisR);   !!dd_comm%sendisR =>NULL()
      if( _ALLOCATED(d_comm%sendjsR) )  DEALLOCATE(d_comm%sendjsR);   !!dd_comm%sendjsR =>NULL()
      if( _ALLOCATED(d_comm%rem_addr) ) DEALLOCATE(d_comm%rem_addr);  !!dd_comm%rem_addr  =>NULL()
      if( _ALLOCATED(d_comm%rem_addrx) )DEALLOCATE(d_comm%rem_addrx); !!dd_comm%rem_addrx =>NULL()
      if( _ALLOCATED(d_comm%rem_addry) )DEALLOCATE(d_comm%rem_addry); !!dd_comm%rem_addry =>NULL()
      if( _ALLOCATED(d_comm%rem_addrl) )DEALLOCATE(d_comm%rem_addrl); !!dd_comm%rem_addrl =>NULL()
    end subroutine deallocate_comm


    function set_domain_id(d_id,ksize,flags,gtype, position, whalo, ehalo, shalo, nhalo)
      integer(LONG_KIND), intent(in) :: d_id
      integer           , intent(in) :: ksize
      integer           , optional, intent(in) :: flags
      integer           , optional, intent(in) :: gtype
      integer           , optional, intent(in) :: position
      integer           , optional, intent(in) :: whalo, ehalo, shalo, nhalo

      integer(LONG_KIND)             :: set_domain_id
      
      set_domain_id=d_id + KE_BASE*int(ksize,KIND(d_id))
      if(PRESENT(flags))set_domain_id=set_domain_id+int(flags,KIND(d_id))
      if(PRESENT(gtype))set_domain_id=set_domain_id+GT_BASE*int(gtype,KIND(d_id))  ! Must be LONG_KIND arithmetic
      !--- gtype is never been used to set id. we need to add position to calculate id to seperate
      !--- BGRID and CGRID or scalar variable.
      if(present(position)) set_domain_id=set_domain_id+GT_BASE*int(2**position, KIND(d_id))
      !z1l ???? the following calculation may need to be revised
      if(present(whalo)) then
         if(whalo>=0) then
            set_domain_id=set_domain_id+GT_BASE*int(2**4*2**whalo, KIND(d_id))
         else
            set_domain_id=set_domain_id-GT_BASE*int(2**4*2**(-whalo), KIND(d_id))
         endif
      end if
      if(present(ehalo)) then
         if(ehalo>=0) then
            set_domain_id=set_domain_id+GT_BASE*int(2**4*2**ehalo, KIND(d_id))
         else
            set_domain_id=set_domain_id-GT_BASE*int(2**4*2**(-ehalo), KIND(d_id))
         endif
       end if
      if(present(shalo)) then
         if(shalo>=0) then
            set_domain_id=set_domain_id+GT_BASE*int(2**4*2**shalo, KIND(d_id))
         else
            set_domain_id=set_domain_id-GT_BASE*int(2**4*2**(-shalo), KIND(d_id))
         endif
      end if
      if(present(nhalo)) then
         if(nhalo>=0) then
            set_domain_id=set_domain_id+GT_BASE*int(2**4*2**nhalo, KIND(d_id))
         else
            set_domain_id=set_domain_id-GT_BASE*int(2**4*2**(-nhalo), KIND(d_id))
         endif
       end if
    end function set_domain_id


!#######################################################################

