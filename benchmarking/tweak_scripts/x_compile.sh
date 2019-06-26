echo compiling netcdf test
mpiifort -Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8 -DRRTM_NO_COMPILE -I/usr/local/include -fpp -stack_temps -safe_cray_ptr -ftz -assume byterecl -shared-intel -i4 -r8 -g -O2 -diag-disable 6843 -L/home/qv18258/files/local/lib/ -I/home/qv18258/files/local/include -lnetcdff test.f90 -o nettest.x

echo compiling originator.F90
mpiifort -Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8 -DRRTM_NO_COMPILE -I/usr/local/include -I/home/qv18258/files/local/include -fpp -stack_temps -safe_cray_ptr -ftz -assume byterecl  -shared-intel -i4 -r8 -g -O2 -diag-disable 6843 -L/home/qv18258/files/local/lib/  -c y_originator.f90

echo compiling y_intermediate.F90
mpiifort -Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8 -DRRTM_NO_COMPILE -I/usr/local/include -I/home/qv18258/files/local/include -fpp -stack_temps -safe_cray_ptr -ftz -assume byterecl  -shared-intel -i4 -r8 -g -O2 -diag-disable 6843 -L/home/qv18258/files/local/lib/  -c y_intermediate.f90

echo compiling y_compile_me.F90
mpiifort -Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8 -DRRTM_NO_COMPILE -I/usr/local/include -I/home/qv18258/files/local/include -fpp -stack_temps -safe_cray_ptr -ftz -assume byterecl  -shared-intel -i4 -r8 -g -O2 -diag-disable 6843 -L/home/qv18258/files/local/lib/  -c y_compile_me.F90
