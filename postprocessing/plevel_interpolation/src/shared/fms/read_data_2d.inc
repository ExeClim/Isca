! -*-f90-*-
! $Id: read_data_2d.inc,v 13.0 2006/03/28 21:39:18 fms Exp $


   if((size(data,1)== ied-isd+1).and.(size(data,2)==jed-jsd+1)) then
     no_halo = .false.
   elseif ((size(data,1)== ie-is+1).and.(size(data,2)==je-js+1))  then
     no_halo = .true.
   else
     call mpp_error(FATAL,'fms_io read_data_2d, data must be either data or compute domain')
   endif
  

   if (.not.module_is_initialized) call fms_io_init ( )

   if (.not.associated(Current_domain)) call mpp_error &
        (FATAL,'set_domain not called')

   if (present(end)) end = .false. 

   if (do_read()) read (unit,end=10) gdata
   if (.not.read_all_pe) then
      ! did not read on all PEs need to broadcast data
      len = size(gdata,1)*size(gdata,2)
!           call mpp_transmit ( gdata, len, ALL_PES, &
!                               gdata, len, 0        ) 
      call mpp_broadcast ( gdata, len, mpp_root_pe() )      
   endif
   if(.not. no_halo .or. read_data_bug) then
! data dimensioned on data domain
      data(is:ie,js:je) = gdata(is:ie,js:je)
   else
! data defined on compute domain
      data(:,:) = gdata(is:ie,js:je)
   endif
   return  

10 call read_eof (end)

