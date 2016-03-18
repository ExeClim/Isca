    subroutine MPP_READ_2DDECOMP_2D_( unit, field, domain, data, tindex, tile_count )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:)
      integer, intent(in), optional :: tindex, tile_count
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),1)
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_read( unit, field, domain, data3D, tindex, tile_count)
      return
    end subroutine MPP_READ_2DDECOMP_2D_

    subroutine MPP_READ_2DDECOMP_3D_( unit, field, domain, data, tindex, tile_count )
!mpp_read reads <data> which has the domain decomposition <domain>
      integer,           intent(in) :: unit
      type(fieldtype),   intent(in) :: field
      type(domain2D),    intent(in) :: domain
      MPP_TYPE_,      intent(inout) :: data(:,:,:)
      integer, intent(in), optional :: tindex, tile_count

      MPP_TYPE_, allocatable :: cdata(:,:,:)
      MPP_TYPE_, allocatable :: gdata(:)
      integer :: len, lenx,leny,lenz,i,j,k,n
!NEW: data may be on compute OR data domain
      logical :: data_has_halos, halos_are_global, x_is_global, y_is_global
      integer :: is, ie, js, je, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ism, iem, jsm, jem
      integer :: ioff, joff, position

      call mpp_clock_begin(mpp_read_clock)
      
      if (.NOT. present(tindex) .AND. mpp_file(unit)%time_level .ne. -1) &
      call mpp_error(FATAL, 'MPP_READ: need to specify a time level for data with time axis')

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_READ: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_READ: invalid unit number.' )

      call mpp_get_compute_domain( domain, is,  ie,  js,  je, tile_count=tile_count )
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, x_is_global=x_is_global, &
                                   y_is_global=y_is_global, tile_count=tile_count )
      call mpp_get_memory_domain ( domain, ism, iem, jsm, jem )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, tile_count=tile_count )

      ! when domain is symmetry, extra point is needed for some data on x/y direction
      position = CENTER
      if(mpp_domain_is_symmetry(domain)) then
         if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+1 ) then  ! CENTER
            data_has_halos = .FALSE.
         else if( size(data,1).EQ.ie-is+2 .AND. size(data,2).EQ.je-js+1 ) then ! EAST
            data_has_halos = .FALSE.
            position = EAST
            ie = ie + 1
         else if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+2 ) then ! NORTH
            position = NORTH
            data_has_halos = .FALSE.
            je = je + 1
         else if( size(data,1).EQ.ie-is+2 .AND. size(data,2).EQ.je-js+2 ) then ! CORNER
            position = CORNER
            data_has_halos = .FALSE.
            ie = ie + 1;   je = je + 1
         else if( size(data,1).EQ.iem-ism+1 .AND. size(data,2).EQ.jem-jsm+1 )then ! CENTER
            data_has_halos = .TRUE.
         else if( size(data,1).EQ.iem-ism+2 .AND. size(data,2).EQ.jem-jsm+1 )then ! EAST
            position = EAST
            data_has_halos = .TRUE.
            ie = ie + 1
         else if( size(data,1).EQ.iem-ism+1 .AND. size(data,2).EQ.jem-jsm+2 )then ! NORTH
            position = NORTH
            data_has_halos = .TRUE.
            je = je + 1
         else if( size(data,1).EQ.iem-ism+2 .AND. size(data,2).EQ.jem-jsm+2 )then ! CORNER
            position = CORNER
            data_has_halos = .TRUE.
            ie = ie + 1;  je = je + 1
         else
            call mpp_error( FATAL, 'MPP_READ: when domain is symmetry, data must be either on ' &
                      //'compute domain or data domain with the consideration of shifting.' )
         end if
      else
         if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+1 )then
            data_has_halos = .FALSE.
         else if( size(data,1).EQ.iem-ism+1 .AND. size(data,2).EQ.jem-jsm+1 )then
            data_has_halos = .TRUE.
         else
            call mpp_error( FATAL, 'MPP_READ: data must be either on compute domain or data domain.' )
         end if
      endif
      halos_are_global = x_is_global .AND. y_is_global
      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( halos_are_global )then !you can read directly into data array
              if( pe.EQ.0 )call read_record( unit, field, size(data(:,:,:)), data, tindex )
          else
              lenx=size(data,1)
              leny=size(data,2)
              lenz=size(data,3)
              len=lenx*leny*lenz
              allocate(gdata(len))          
! read field on pe 0 and pass to all pes
              if( pe.EQ.0 ) call read_record( unit, field, len, gdata, tindex )
! broadcasting global array, this can be expensive!          
              call mpp_transmit( put_data=gdata(1), plen=len, to_pe=ALL_PES, &
                                 get_data=gdata(1), glen=len, from_pe=0 )
              ioff = is; joff = js
              if( data_has_halos )then
                  ioff = isd; joff = jsd
              end if
              do k=1,size(data,3)
                 do j=js,je
                    do i=is,ie
                       n=(i-isg+1) + (j-jsg)*lenx + (k-1)*lenx*leny
                       data(i-ioff+1,j-joff+1,k)=gdata(n)
                    enddo
                 enddo
              enddo
              deallocate(gdata)
          end if
      else if( data_has_halos )then
! for uniprocessor or multithreaded read
! read compute domain as contiguous data

          allocate( cdata(is:ie,js:je,size(data,3)) )
          call read_record(unit,field,size(cdata(:,:,:)),cdata,tindex,domain,position,tile_count)

          data(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1,:) = cdata(:,:,:)
          deallocate(cdata)
      else
          call read_record(unit,field,size(data(:,:,:)),data,tindex,domain,position,tile_count)
      end if

      call mpp_clock_end(mpp_read_clock)

      return
    end subroutine MPP_READ_2DDECOMP_3D_
