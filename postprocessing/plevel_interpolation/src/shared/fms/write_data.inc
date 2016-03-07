
   if (.not.module_is_initialized) call fms_io_init ( )

   if (.not.associated(Current_domain)) call mpp_error &
        (FATAL, 'set_domain not called')

 ! get the global data and write only on root pe
   call mpp_global_field ( Current_domain, data, gdata ) 
   if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata

