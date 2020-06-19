# template for the Intel fortran compiler
# typical use with mkmf
# mkmf -t template.ifc -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include
CPPFLAGS = -I/usr/local/include
NETCDF_LIBS = `nc-config --libs`

# FFLAGS:
#  -fpp: Use the fortran preprocessor
#  -stack_temps:  Put temporary runtime arrays on the stack, not heap.
#  -safe_cray_ptr: Cray pointers don't alias other variables.
#  -ftz: Denormal numbers are flushed to zero.
#  -assume byterecl: Specifies the units for the OPEN statement as bytes.
#  -shared-intel:  Load intel libraries dynamically
#  -i4: 4 byte integers
#  -r8: 8 byte reals
#  -g: Generate symbolic debugging info in code
#  -O2: Level 2 speed optimisations
#  -diag-disable 6843:
#       This suppresses the warning: `warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.` of which
#       there are a lot of instances in the GFDL codebase.
FFLAGS = -Duse_netCDF3 $(CPPFLAGS) -fpp -stack_temps -safe_cray_ptr -ftz -assume byterecl -shared-intel -i4 -r8 -g -O2 -diag-disable 6843 -mcmodel large
#FFLAGS = $(CPPFLAGS) -fltconsistency -stack_temps -safe_cray_ptr -ftz -shared-intel -assume byterecl -g -O0 -i4 -r8 -check -warn -warn noerrors -debug variable_locations -inline_debug_info -traceback
FC = $(F90)
LD = $(F90) $(NETCDF_LIBS)
#CC = mpicc

LDFLAGS = -lnetcdff -lnetcdf -lmpi -shared-intel
CFLAGS = -D__IFC
