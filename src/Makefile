# Makefile created by mkmf $Id: mkmf,v 16.1 2010/05/19 18:49:19 fms Exp $ 

SRCROOT = /scratch/sit204/FMS2013/GFDLmoistModel/src/

CPPDEFS = -Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8


include /scratch/sit204/FMS2013/GFDLmoistModel/bin/mkmf.template.ia64


.DEFAULT:
	-echo $@ does not exist.
all: fms_moist.x
MersenneTwister.o: $(SRCROOT)shared/random_numbers/MersenneTwister.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/random_numbers/MersenneTwister.F90
affinity.o: $(SRCROOT)shared/mpp/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/affinity.c
atmos_model.o: $(SRCROOT)atmos_solo/atmos_model.F90 mpp.o fms.o atmosphere.o time_manager.o fms_io.o mpp_domains.o mpp_io.o diag_manager.o field_manager.o tracer_manager.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_solo/atmos_model.F90
atmosphere.o: $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90 mpp.o fms.o constants.o transforms.o time_manager.o press_and_geopot.o spectral_dynamics.o tracer_type.o hs_forcing.o field_manager.o tracer_manager.o idealized_moist_phys.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90
axis_utils.o: $(SRCROOT)shared/axis_utils/axis_utils.F90 mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/axis_utils/axis_utils.F90
constants.o: $(SRCROOT)shared/constants/constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/constants/constants.F90
create_xgrid.o: $(SRCROOT)shared/mosaic/create_xgrid.c /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/create_xgrid.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic	$(SRCROOT)shared/mosaic/create_xgrid.c
diag_axis.o: $(SRCROOT)shared/diag_manager/diag_axis.F90 mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)shared/diag_manager/diag_data.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)shared/diag_manager/diag_grid.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/diag_manager/diag_grid.F90
diag_manager.o: $(SRCROOT)shared/diag_manager/diag_manager.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h time_manager.o mpp_io.o mpp.o fms.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o constants.o mpp_domains.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/diag_manager/diag_manager.F90
diag_output.o: $(SRCROOT)shared/diag_manager/diag_output.F90 mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)shared/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)shared/diag_manager/diag_util.F90 diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_util.F90
diffusivity.o: $(SRCROOT)atmos_param/diffusivity/diffusivity.F90 constants.o mpp.o fms.o monin_obukhov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/diffusivity/diffusivity.F90
edt.o: $(SRCROOT)atmos_param/edt/edt.F90 constants.o mpp.o fms.o fms_io.o diag_manager.o time_manager.o monin_obukhov.o sat_vapor_pres.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/edt/edt.F90
entrain.o: $(SRCROOT)atmos_param/entrain/entrain.F90 constants.o mpp.o fms.o diag_manager.o time_manager.o sat_vapor_pres.o monin_obukhov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/entrain/entrain.F90
every_step_diagnostics.o: $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90 fms.o time_manager.o field_manager.o tracer_manager.o diag_manager.o press_and_geopot.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90
fft.o: $(SRCROOT)shared/fft/fft.F90 platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fft/fft.F90
fft99.o: $(SRCROOT)shared/fft/fft99.F90 constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fft/fft99.F90
field_manager.o: $(SRCROOT)shared/field_manager/field_manager.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/field_manager	$(SRCROOT)shared/field_manager/field_manager.F90
fm_util.o: $(SRCROOT)shared/field_manager/fm_util.F90 field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/field_manager/fm_util.F90
fms.o: $(SRCROOT)shared/fms/fms.F90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fms/fms.F90
fms_io.o: $(SRCROOT)shared/fms/fms_io.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_2d.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_3d.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_4d.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms	$(SRCROOT)shared/fms/fms_io.F90
fv_advection.o: $(SRCROOT)atmos_spectral/model/fv_advection.F90 fms.o constants.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/fv_advection.F90
gauss_and_legendre.o: $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90
gaussian_topog.o: $(SRCROOT)shared/topography/gaussian_topog.F90 fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/topography/gaussian_topog.F90
get_cal_time.o: $(SRCROOT)shared/time_manager/get_cal_time.F90 fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/get_cal_time.F90
global_integral.o: $(SRCROOT)atmos_spectral/model/global_integral.F90 fms.o press_and_geopot.o transforms.o constants.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/global_integral.F90
gradient.o: $(SRCROOT)shared/mosaic/gradient.F90 mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)shared/mosaic/gradient_c2l.c /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/constant.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic	$(SRCROOT)shared/mosaic/gradient_c2l.c
grid.o: $(SRCROOT)shared/mosaic/grid.F90 mpp.o constants.o fms.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/grid.F90
grid_fourier.o: $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 fms.o constants.o fft.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/grid_fourier.F90
horiz_interp.o: $(SRCROOT)shared/horiz_interp/horiz_interp.F90 fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
hs_forcing.o: $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90 mpp.o fms.o constants.o time_manager.o diag_manager.o field_manager.o tracer_manager.o interpolator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90
ic_from_external_file.o: $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90 mpp.o fms.o transforms.o tracer_type.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/ic_from_external_file.F90
idealized_moist_phys.o: $(SRCROOT)atmos_spectral/driver/solo/idealized_moist_phys.F90 mpp.o fms.o constants.o time_manager.o vert_turb_driver.o vert_diff.o two_stream_gray_rad.o mixed_layer.o lscale_cond.o qe_moist_convection.o diag_manager.o transforms.o spectral_dynamics.o surface_flux.o sat_vapor_pres.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/driver/solo/idealized_moist_phys.F90
implicit.o: $(SRCROOT)atmos_spectral/model/implicit.F90 fms.o constants.o press_and_geopot.o matrix_invert.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/implicit.F90
interp.o: $(SRCROOT)shared/mosaic/interp.c /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/interp.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic	$(SRCROOT)shared/mosaic/interp.c
interpolator.o: $(SRCROOT)atmos_shared/interpolator/interpolator.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h mpp.o mpp_io.o mpp_domains.o diag_manager.o fms.o horiz_interp.o time_manager.o time_interp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)atmos_shared/interpolator/interpolator.F90
jablonowski_2006.o: $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90 mpp.o fms.o constants.o vert_coordinate.o transforms.o press_and_geopot.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/jablonowski_2006.F90
leapfrog.o: $(SRCROOT)atmos_spectral/model/leapfrog.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/leapfrog.F90
lscale_cond.o: $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90 mpp.o fms.o sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90
matrix_invert.o: $(SRCROOT)atmos_spectral/model/matrix_invert.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/matrix_invert.F90
memuse.o: $(SRCROOT)shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/memutils/memuse.c
memutils.o: $(SRCROOT)shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/memutils/memutils.F90
mixed_layer.o: $(SRCROOT)atmos_spectral/driver/solo/mixed_layer.F90 fms.o field_manager.o tracer_manager.o constants.o diag_manager.o time_manager.o transforms.o vert_diff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/driver/solo/mixed_layer.F90
monin_obukhov.o: /scratch/sit204/old_FMS_model_work/FMS2006/MoistModel/src/atmos_param/monin_obukhov/monin_obukhov.f90 constants.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	/scratch/sit204/old_FMS_model_work/FMS2006/MoistModel/src/atmos_param/monin_obukhov/monin_obukhov.f90
mosaic.o: $(SRCROOT)shared/mosaic/mosaic.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/mosaic.F90
mosaic_util.o: $(SRCROOT)shared/mosaic/mosaic_util.c /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic	$(SRCROOT)shared/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)shared/mpp/mpp.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/system_clock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_sma.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_mpi.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_nocomm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_a.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_s.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_a.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_s.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_sma.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_sma.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_sma.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_sma.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_mpi.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_mpi.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_mpi.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_mpi.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_nocomm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_nocomm.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_nocomm.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_nocomm.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_int.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_scalar.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_gather.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp.F90
mpp_data.o: $(SRCROOT)shared/mpp/mpp_data.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_sma.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_mpi.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)shared/mpp/mpp_domains.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_define_nest_domains.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_util.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_comm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_define.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_misc.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D_nonblock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nonblock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV_nonblock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_check.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_checkV.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_nest_domains.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nest.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_redistribute.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_get_boundary.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_get_boundary.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_reduce.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_reduce.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum_tl.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_field.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_global_field.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_domains.F90
mpp_io.o: $(SRCROOT)shared/mpp/mpp_io.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_util.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_misc.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_connect.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_read.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_read_2Ddecomp.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_write.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write_2Ddecomp.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_io.F90
mpp_memutils.o: $(SRCROOT)shared/mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)shared/mpp/mpp_parameter.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)shared/mpp/mpp_pset.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)shared/mpp/mpp_utilities.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/mpp_utilities.F90
my25_turb.o: $(SRCROOT)atmos_param/my25_turb/my25_turb.F90 mpp.o fms.o fms_io.o tridiagonal.o constants.o monin_obukhov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/my25_turb/my25_turb.F90
nsclock.o: $(SRCROOT)shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/nsclock.c
platform.o: $(SRCROOT)shared/platform/platform.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/platform/platform.F90
polvani_2004.o: $(SRCROOT)atmos_spectral/init/polvani_2004.F90 mpp.o fms.o constants.o vert_coordinate.o transforms.o press_and_geopot.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/polvani_2004.F90
polvani_2007.o: $(SRCROOT)atmos_spectral/init/polvani_2007.F90 mpp.o fms.o constants.o vert_coordinate.o transforms.o press_and_geopot.o diag_manager.o tracer_type.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/polvani_2007.F90
press_and_geopot.o: $(SRCROOT)atmos_spectral/model/press_and_geopot.F90 fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/press_and_geopot.F90
qe_moist_convection.o: $(SRCROOT)atmos_param/qe_moist_convection/qe_moist_convection.F90 fms.o sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/qe_moist_convection/qe_moist_convection.F90
rad_utilities.o: $(SRCROOT)atmos_param/sea_esf_rad/null/rad_utilities.F90 fms.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/null/rad_utilities.F90
random_numbers.o: $(SRCROOT)shared/random_numbers/random_numbers.F90 MersenneTwister.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/random_numbers/random_numbers.F90
read_mosaic.o: $(SRCROOT)shared/mosaic/read_mosaic.c /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/read_mosaic.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/constant.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic	$(SRCROOT)shared/mosaic/read_mosaic.c
sat_vapor_pres.o: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 constants.o fms.o mpp_io.o mpp.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
shallow_conv.o: $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90 sat_vapor_pres.o mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90
spec_mpp.o: $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 fms.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spec_mpp.F90
spectral_damping.o: $(SRCROOT)atmos_spectral/model/spectral_damping.F90 fms.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/spectral_damping.F90
spectral_dynamics.o: $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90 mpp.o fms.o constants.o time_manager.o field_manager.o tracer_manager.o diag_manager.o transforms.o vert_advection.o implicit.o press_and_geopot.o spectral_damping.o leapfrog.o fv_advection.o water_borrowing.o global_integral.o spectral_init_cond.o tracer_type.o every_step_diagnostics.o mpp_domains.o polvani_2004.o polvani_2007.o jablonowski_2006.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/spectral_dynamics.F90
spectral_init_cond.o: $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90 mpp.o fms.o mpp_domains.o constants.o vert_coordinate.o transforms.o press_and_geopot.o spectral_initialize_fields.o ic_from_external_file.o topog_regularization.o topography.o tracer_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/spectral_init_cond.F90
spectral_initialize_fields.o: $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90 fms.o constants.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90
spherical.o: $(SRCROOT)atmos_spectral/tools/spherical.F90 fms.o spec_mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spherical.F90
spherical_fourier.o: $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 fms.o mpp.o mpp_domains.o constants.o spec_mpp.o spherical.o gauss_and_legendre.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spherical_fourier.F90
stable_bl_turb.o: $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90 mpp.o fms.o diag_manager.o time_manager.o constants.o monin_obukhov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90
strat_cloud.o: $(SRCROOT)atmos_param/strat_cloud/null/strat_cloud.F90 time_manager.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/null/strat_cloud.F90
surface_flux.o: $(SRCROOT)coupler/surface_flux.F90 fms.o monin_obukhov.o sat_vapor_pres.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)coupler/surface_flux.F90
test_fms_io.o: $(SRCROOT)shared/fms/test_fms_io.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/fms/test_fms_io.F90
test_mpp.o: $(SRCROOT)shared/mpp/test_mpp.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/mpp/test_mpp.F90
test_mpp_domains.o: $(SRCROOT)shared/mpp/test_mpp_domains.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/mpp/test_mpp_domains.F90
test_mpp_io.o: $(SRCROOT)shared/mpp/test_mpp_io.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/mpp/test_mpp_io.F90
test_mpp_pset.o: $(SRCROOT)shared/mpp/test_mpp_pset.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/mpp/test_mpp_pset.F90
threadloc.o: $(SRCROOT)shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/threadloc.c
time_interp.o: $(SRCROOT)shared/time_interp/time_interp.F90 time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_interp/time_interp.F90
time_interp_external.o: $(SRCROOT)shared/time_interp/time_interp_external.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o fms.o platform.o horiz_interp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include	$(SRCROOT)shared/time_interp/time_interp_external.F90
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/time_manager.F90
topog_regularization.o: $(SRCROOT)atmos_spectral/init/topog_regularization.F90 fms.o mpp.o constants.o transforms.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/topog_regularization.F90
topography.o: $(SRCROOT)shared/topography/topography.F90 gaussian_topog.o horiz_interp.o fms.o fms_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/topography/topography.F90
tracer_manager.o: $(SRCROOT)shared/tracer_manager/tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/tracer_manager/tracer_manager.F90
tracer_type.o: $(SRCROOT)atmos_spectral/model/tracer_type.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/tracer_type.F90
transforms.o: $(SRCROOT)atmos_spectral/tools/transforms.F90 mpp.o fms.o mpp_domains.o spec_mpp.o constants.o spherical_fourier.o grid_fourier.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/transforms.F90
tridiagonal.o: $(SRCROOT)shared/tridiagonal/tridiagonal.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/tridiagonal/tridiagonal.F90
two_stream_gray_rad.o: $(SRCROOT)atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90 fms.o constants.o diag_manager.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90
vert_advection.o: $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_shared/vert_advection/vert_advection.F90
vert_coordinate.o: $(SRCROOT)atmos_spectral/init/vert_coordinate.F90 mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/vert_coordinate.F90
vert_diff.o: $(SRCROOT)atmos_param/vert_diff/vert_diff.F90 constants.o fms.o field_manager.o tracer_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/vert_diff/vert_diff.F90
vert_turb_driver.o: $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90 my25_turb.o diffusivity.o edt.o strat_cloud.o shallow_conv.o stable_bl_turb.o entrain.o diag_manager.o time_manager.o constants.o mpp.o fms.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90
water_borrowing.o: $(SRCROOT)atmos_spectral/model/water_borrowing.F90 fms.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/water_borrowing.F90
./mpp_global_reduce.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_reduce.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_reduce.h .
./spectral_initialize_fields.F90: $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90
	cp $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90 .
./interpolator.F90: $(SRCROOT)atmos_shared/interpolator/interpolator.F90
	cp $(SRCROOT)atmos_shared/interpolator/interpolator.F90 .
./entrain.F90: $(SRCROOT)atmos_param/entrain/entrain.F90
	cp $(SRCROOT)atmos_param/entrain/entrain.F90 .
./mpp_update_domains2D.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D.h .
./mpp_pset.F90: $(SRCROOT)shared/mpp/mpp_pset.F90
	cp $(SRCROOT)shared/mpp/mpp_pset.F90 .
./test_mpp_pset.F90: $(SRCROOT)shared/mpp/test_mpp_pset.F90
	cp $(SRCROOT)shared/mpp/test_mpp_pset.F90 .
./polvani_2004.F90: $(SRCROOT)atmos_spectral/init/polvani_2004.F90
	cp $(SRCROOT)atmos_spectral/init/polvani_2004.F90 .
./fms_io.F90: $(SRCROOT)shared/fms/fms_io.F90
	cp $(SRCROOT)shared/fms/fms_io.F90 .
./two_stream_gray_rad.F90: $(SRCROOT)atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90
	cp $(SRCROOT)atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90 .
./write_data.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/write_data.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/write_data.inc .
./mpp_do_updateV.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV.h .
./read_data_2d.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_2d.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_2d.inc .
./mpp_transmit_mpi.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_mpi.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_mpi.h .
./vert_coordinate.F90: $(SRCROOT)atmos_spectral/init/vert_coordinate.F90
	cp $(SRCROOT)atmos_spectral/init/vert_coordinate.F90 .
./mpp_data_sma.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_sma.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_sma.inc .
./read_data_4d.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_4d.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_4d.inc .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./mosaic.F90: $(SRCROOT)shared/mosaic/mosaic.F90
	cp $(SRCROOT)shared/mosaic/mosaic.F90 .
./vert_diff.F90: $(SRCROOT)atmos_param/vert_diff/vert_diff.F90
	cp $(SRCROOT)atmos_param/vert_diff/vert_diff.F90 .
./gaussian_topog.F90: $(SRCROOT)shared/topography/gaussian_topog.F90
	cp $(SRCROOT)shared/topography/gaussian_topog.F90 .
./mpp_transmit.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit.inc .
./fft.F90: $(SRCROOT)shared/fft/fft.F90
	cp $(SRCROOT)shared/fft/fft.F90 .
./global_integral.F90: $(SRCROOT)atmos_spectral/model/global_integral.F90
	cp $(SRCROOT)atmos_spectral/model/global_integral.F90 .
./vert_turb_driver.F90: $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90
	cp $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90 .
./atmosphere.F90: $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90
	cp $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90 .
./spectral_dynamics.F90: $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90
	cp $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90 .
./polvani_2007.F90: $(SRCROOT)atmos_spectral/init/polvani_2007.F90
	cp $(SRCROOT)atmos_spectral/init/polvani_2007.F90 .
./mpp_io_connect.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_connect.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_connect.inc .
./matrix_invert.F90: $(SRCROOT)atmos_spectral/model/matrix_invert.F90
	cp $(SRCROOT)atmos_spectral/model/matrix_invert.F90 .
./random_numbers.F90: $(SRCROOT)shared/random_numbers/random_numbers.F90
	cp $(SRCROOT)shared/random_numbers/random_numbers.F90 .
./mpp_data.F90: $(SRCROOT)shared/mpp/mpp_data.F90
	cp $(SRCROOT)shared/mpp/mpp_data.F90 .
./lscale_cond.F90: $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90
	cp $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90 .
./gradient.F90: $(SRCROOT)shared/mosaic/gradient.F90
	cp $(SRCROOT)shared/mosaic/gradient.F90 .
./horiz_interp_type.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 .
./stable_bl_turb.F90: $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90
	cp $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90 .
./mpp_util_nocomm.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_nocomm.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_nocomm.inc .
./mpp_transmit_sma.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_sma.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_sma.h .
./mpp_global_field.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_field.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_field.h .
./mpp_global_sum_tl.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum_tl.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum_tl.h .
./read_mosaic.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/read_mosaic.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/read_mosaic.h .
./mpp_data_nocomm.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_nocomm.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_nocomm.inc .
./mpp_util_mpi.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_mpi.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_mpi.inc .
./transforms.F90: $(SRCROOT)atmos_spectral/tools/transforms.F90
	cp $(SRCROOT)atmos_spectral/tools/transforms.F90 .
./mpp_update_domains2D_nonblock.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D_nonblock.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D_nonblock.h .
./mpp_reduce_mpi.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_mpi.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_mpi.h .
./mpp_utilities.F90: $(SRCROOT)shared/mpp/mpp_utilities.F90
	cp $(SRCROOT)shared/mpp/mpp_utilities.F90 .
./mpp_do_check.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_check.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_check.h .
./axis_utils.F90: $(SRCROOT)shared/axis_utils/axis_utils.F90
	cp $(SRCROOT)shared/axis_utils/axis_utils.F90 .
./mixed_layer.F90: $(SRCROOT)atmos_spectral/driver/solo/mixed_layer.F90
	cp $(SRCROOT)atmos_spectral/driver/solo/mixed_layer.F90 .
./mpp_chksum.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum.h .
./mpp_chksum_int.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_int.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_int.h .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./mpp_write.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write.h .
./mpp_get_boundary.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_get_boundary.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_get_boundary.h .
./parse.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/field_manager/parse.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/field_manager/parse.inc .
./edt.F90: $(SRCROOT)atmos_param/edt/edt.F90
	cp $(SRCROOT)atmos_param/edt/edt.F90 .
./test_mpp_domains.F90: $(SRCROOT)shared/mpp/test_mpp_domains.F90
	cp $(SRCROOT)shared/mpp/test_mpp_domains.F90 .
./spherical_fourier.F90: $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90
	cp $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 .
./read_data_3d.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_3d.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_3d.inc .
./diag_axis.F90: $(SRCROOT)shared/diag_manager/diag_axis.F90
	cp $(SRCROOT)shared/diag_manager/diag_axis.F90 .
./fms_platform.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h .
./horiz_interp_conserve.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 .
./mpp_gather.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_gather.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_gather.h .
./fft99.F90: $(SRCROOT)shared/fft/fft99.F90
	cp $(SRCROOT)shared/fft/fft99.F90 .
./gradient_c2l.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/gradient_c2l.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/gradient_c2l.h .
./qe_moist_convection.F90: $(SRCROOT)atmos_param/qe_moist_convection/qe_moist_convection.F90
	cp $(SRCROOT)atmos_param/qe_moist_convection/qe_moist_convection.F90 .
./diag_util.F90: $(SRCROOT)shared/diag_manager/diag_util.F90
	cp $(SRCROOT)shared/diag_manager/diag_util.F90 .
./mpp_update_nest_domains.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_nest_domains.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_nest_domains.h .
./diag_manager.F90: $(SRCROOT)shared/diag_manager/diag_manager.F90
	cp $(SRCROOT)shared/diag_manager/diag_manager.F90 .
./mpp_error_s_a.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_a.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_a.h .
./mpp_error_a_s.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_s.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_s.h .
./mpp_define_nest_domains.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_define_nest_domains.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_define_nest_domains.inc .
./horiz_interp_bicubic.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 .
./water_borrowing.F90: $(SRCROOT)atmos_spectral/model/water_borrowing.F90
	cp $(SRCROOT)atmos_spectral/model/water_borrowing.F90 .
./fms.F90: $(SRCROOT)shared/fms/fms.F90
	cp $(SRCROOT)shared/fms/fms.F90 .
./mosaic_util.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h .
./leapfrog.F90: $(SRCROOT)atmos_spectral/model/leapfrog.F90
	cp $(SRCROOT)atmos_spectral/model/leapfrog.F90 .
./field_manager.F90: $(SRCROOT)shared/field_manager/field_manager.F90
	cp $(SRCROOT)shared/field_manager/field_manager.F90 .
./hs_forcing.F90: $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90
	cp $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90 .
./fm_util.F90: $(SRCROOT)shared/field_manager/fm_util.F90
	cp $(SRCROOT)shared/field_manager/fm_util.F90 .
./system_clock.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/system_clock.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/system_clock.h .
./mpp_comm_sma.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_sma.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_sma.inc .
./mpp_io_write.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_write.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_write.inc .
./sat_vapor_pres.F90: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90
	cp $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 .
./mpp_do_update.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update.h .
./mpp_error_s_s.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_s.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_s.h .
./mpp_io.F90: $(SRCROOT)shared/mpp/mpp_io.F90
	cp $(SRCROOT)shared/mpp/mpp_io.F90 .
./mpp_util.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util.inc .
./mpp_do_update_nonblock.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nonblock.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nonblock.h .
./mpp_comm_nocomm.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_nocomm.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_nocomm.inc .
./mpp_domains_reduce.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_reduce.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_reduce.inc .
./mpp.F90: $(SRCROOT)shared/mpp/mpp.F90
	cp $(SRCROOT)shared/mpp/mpp.F90 .
./mpp_io_misc.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_misc.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_misc.inc .
./tridiagonal.F90: $(SRCROOT)shared/tridiagonal/tridiagonal.F90
	cp $(SRCROOT)shared/tridiagonal/tridiagonal.F90 .
./memuse.c: $(SRCROOT)shared/memutils/memuse.c
	cp $(SRCROOT)shared/memutils/memuse.c .
./mpp_domains_comm.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_comm.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_comm.inc .
./mpp_data_mpi.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_mpi.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_mpi.inc .
./mpp_domains_misc.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_misc.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_misc.inc .
./my25_turb.F90: $(SRCROOT)atmos_param/my25_turb/my25_turb.F90
	cp $(SRCROOT)atmos_param/my25_turb/my25_turb.F90 .
./implicit.F90: $(SRCROOT)atmos_spectral/model/implicit.F90
	cp $(SRCROOT)atmos_spectral/model/implicit.F90 .
./diag_table.F90: $(SRCROOT)shared/diag_manager/diag_table.F90
	cp $(SRCROOT)shared/diag_manager/diag_table.F90 .
./mpp_sum_mpi.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_mpi.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_mpi.h .
./spherical.F90: $(SRCROOT)atmos_spectral/tools/spherical.F90
	cp $(SRCROOT)atmos_spectral/tools/spherical.F90 .
./mpp_chksum_scalar.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_scalar.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_scalar.h .
./mpp_do_update_nest.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nest.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nest.h .
./tracer_type.F90: $(SRCROOT)atmos_spectral/model/tracer_type.F90
	cp $(SRCROOT)atmos_spectral/model/tracer_type.F90 .
./sat_vapor_pres_k.F90: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
	cp $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 .
./gradient_c2l.c: $(SRCROOT)shared/mosaic/gradient_c2l.c
	cp $(SRCROOT)shared/mosaic/gradient_c2l.c .
./test_mpp_io.F90: $(SRCROOT)shared/mpp/test_mpp_io.F90
	cp $(SRCROOT)shared/mpp/test_mpp_io.F90 .
./diag_grid.F90: $(SRCROOT)shared/diag_manager/diag_grid.F90
	cp $(SRCROOT)shared/diag_manager/diag_grid.F90 .
./mpp_domains.F90: $(SRCROOT)shared/mpp/mpp_domains.F90
	cp $(SRCROOT)shared/mpp/mpp_domains.F90 .
./mpp_transmit_nocomm.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_nocomm.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_nocomm.h .
./rad_utilities.F90: $(SRCROOT)atmos_param/sea_esf_rad/null/rad_utilities.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/null/rad_utilities.F90 .
./spec_mpp.F90: $(SRCROOT)atmos_spectral/tools/spec_mpp.F90
	cp $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 .
./diag_data.F90: $(SRCROOT)shared/diag_manager/diag_data.F90
	cp $(SRCROOT)shared/diag_manager/diag_data.F90 .
./mpp_memutils.F90: $(SRCROOT)shared/mpp/mpp_memutils.F90
	cp $(SRCROOT)shared/mpp/mpp_memutils.F90 .
./mpp_reduce_nocomm.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_nocomm.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_nocomm.h .
./mpp_parameter.F90: $(SRCROOT)shared/mpp/mpp_parameter.F90
	cp $(SRCROOT)shared/mpp/mpp_parameter.F90 .
./mpp_io_read.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_read.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_read.inc .
./ic_from_external_file.F90: $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90
	cp $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90 .
./strat_cloud.F90: $(SRCROOT)atmos_param/strat_cloud/null/strat_cloud.F90
	cp $(SRCROOT)atmos_param/strat_cloud/null/strat_cloud.F90 .
./test_fms_io.F90: $(SRCROOT)shared/fms/test_fms_io.F90
	cp $(SRCROOT)shared/fms/test_fms_io.F90 .
./horiz_interp_spherical.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 .
./memutils.F90: $(SRCROOT)shared/memutils/memutils.F90
	cp $(SRCROOT)shared/memutils/memutils.F90 .
./mpp_util_sma.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_sma.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_sma.inc .
./MersenneTwister.F90: $(SRCROOT)shared/random_numbers/MersenneTwister.F90
	cp $(SRCROOT)shared/random_numbers/MersenneTwister.F90 .
./gauss_and_legendre.F90: $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90
	cp $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 .
./constants.F90: $(SRCROOT)shared/constants/constants.F90
	cp $(SRCROOT)shared/constants/constants.F90 .
./mpp_io_util.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_util.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_util.inc .
./tracer_manager.F90: $(SRCROOT)shared/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)shared/tracer_manager/tracer_manager.F90 .
./every_step_diagnostics.F90: $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90
	cp $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90 .
./mpp_do_get_boundary.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_get_boundary.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_get_boundary.h .
./mpp_write_2Ddecomp.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write_2Ddecomp.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write_2Ddecomp.h .
./vert_advection.F90: $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90
	cp $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 .
./topography.F90: $(SRCROOT)shared/topography/topography.F90
	cp $(SRCROOT)shared/topography/topography.F90 .
./test_mpp.F90: $(SRCROOT)shared/mpp/test_mpp.F90
	cp $(SRCROOT)shared/mpp/test_mpp.F90 .
./mpp_do_checkV.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_checkV.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_checkV.h .
./horiz_interp_bilinear.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 .
./mpp_sum.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum.inc .
./affinity.c: $(SRCROOT)shared/mpp/affinity.c
	cp $(SRCROOT)shared/mpp/affinity.c .
./mpp_do_updateV_nonblock.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV_nonblock.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV_nonblock.h .
./mosaic_util.c: $(SRCROOT)shared/mosaic/mosaic_util.c
	cp $(SRCROOT)shared/mosaic/mosaic_util.c .
./read_mosaic.c: $(SRCROOT)shared/mosaic/read_mosaic.c
	cp $(SRCROOT)shared/mosaic/read_mosaic.c .
./diffusivity.F90: $(SRCROOT)atmos_param/diffusivity/diffusivity.F90
	cp $(SRCROOT)atmos_param/diffusivity/diffusivity.F90 .
./mpp_comm_mpi.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_mpi.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_mpi.inc .
./interp.c: $(SRCROOT)shared/mosaic/interp.c
	cp $(SRCROOT)shared/mosaic/interp.c .
./mpp_comm.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm.inc .
./spectral_init_cond.F90: $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90
	cp $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90 .
./mpp_do_global_field.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_global_field.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_global_field.h .
./grid_fourier.F90: $(SRCROOT)atmos_spectral/tools/grid_fourier.F90
	cp $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 .
./mpp_sum_sma.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_sma.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_sma.h .
./topog_regularization.F90: $(SRCROOT)atmos_spectral/init/topog_regularization.F90
	cp $(SRCROOT)atmos_spectral/init/topog_regularization.F90 .
./spectral_damping.F90: $(SRCROOT)atmos_spectral/model/spectral_damping.F90
	cp $(SRCROOT)atmos_spectral/model/spectral_damping.F90 .
./idealized_moist_phys.F90: $(SRCROOT)atmos_spectral/driver/solo/idealized_moist_phys.F90
	cp $(SRCROOT)atmos_spectral/driver/solo/idealized_moist_phys.F90 .
./constant.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/constant.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/constant.h .
./mpp_sum_nocomm.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_nocomm.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_nocomm.h .
./threadloc.c: $(SRCROOT)shared/mpp/threadloc.c
	cp $(SRCROOT)shared/mpp/threadloc.c .
./surface_flux.F90: $(SRCROOT)coupler/surface_flux.F90
	cp $(SRCROOT)coupler/surface_flux.F90 .
./shallow_conv.F90: $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90
	cp $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90 .
./mpp_global_sum.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum.h .
./monin_obukhov.f90: /scratch/sit204/old_FMS_model_work/FMS2006/MoistModel/src/atmos_param/monin_obukhov/monin_obukhov.f90
	cp /scratch/sit204/old_FMS_model_work/FMS2006/MoistModel/src/atmos_param/monin_obukhov/monin_obukhov.f90 .
./get_cal_time.F90: $(SRCROOT)shared/time_manager/get_cal_time.F90
	cp $(SRCROOT)shared/time_manager/get_cal_time.F90 .
./interp.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/interp.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/interp.h .
./mpp_reduce_sma.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_sma.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_sma.h .
./mpp_domains_util.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_util.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_util.inc .
./time_interp.F90: $(SRCROOT)shared/time_interp/time_interp.F90
	cp $(SRCROOT)shared/time_interp/time_interp.F90 .
./press_and_geopot.F90: $(SRCROOT)atmos_spectral/model/press_and_geopot.F90
	cp $(SRCROOT)atmos_spectral/model/press_and_geopot.F90 .
./horiz_interp.F90: $(SRCROOT)shared/horiz_interp/horiz_interp.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp.F90 .
./create_xgrid.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/create_xgrid.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/create_xgrid.h .
./mpp_do_redistribute.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_redistribute.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_redistribute.h .
./diag_output.F90: $(SRCROOT)shared/diag_manager/diag_output.F90
	cp $(SRCROOT)shared/diag_manager/diag_output.F90 .
./mpp_domains_define.inc: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_define.inc
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_define.inc .
./nsclock.c: $(SRCROOT)shared/mpp/nsclock.c
	cp $(SRCROOT)shared/mpp/nsclock.c .
./jablonowski_2006.F90: $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90
	cp $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90 .
./grid.F90: $(SRCROOT)shared/mosaic/grid.F90
	cp $(SRCROOT)shared/mosaic/grid.F90 .
./mpp_error_a_a.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_a.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_a.h .
./fv_advection.F90: $(SRCROOT)atmos_spectral/model/fv_advection.F90
	cp $(SRCROOT)atmos_spectral/model/fv_advection.F90 .
./atmos_model.F90: $(SRCROOT)atmos_solo/atmos_model.F90
	cp $(SRCROOT)atmos_solo/atmos_model.F90 .
./create_xgrid.c: $(SRCROOT)shared/mosaic/create_xgrid.c
	cp $(SRCROOT)shared/mosaic/create_xgrid.c .
./time_interp_external.F90: $(SRCROOT)shared/time_interp/time_interp_external.F90
	cp $(SRCROOT)shared/time_interp/time_interp_external.F90 .
./mpp_read_2Ddecomp.h: /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_read_2Ddecomp.h
	cp /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_read_2Ddecomp.h .
SRC = $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90 $(SRCROOT)atmos_param/edt/edt.F90 $(SRCROOT)shared/time_interp/time_interp_external.F90 $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)atmos_param/strat_cloud/null/strat_cloud.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)shared/random_numbers/random_numbers.F90 $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)atmos_spectral/init/polvani_2004.F90 $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)shared/topography/gaussian_topog.F90 $(SRCROOT)atmos_param/vert_diff/vert_diff.F90 $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)atmos_spectral/model/implicit.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)atmos_param/qe_moist_convection/qe_moist_convection.F90 $(SRCROOT)atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90 $(SRCROOT)atmos_spectral/tools/transforms.F90 $(SRCROOT)atmos_spectral/init/vert_coordinate.F90 $(SRCROOT)shared/mpp/mpp_utilities.F90 $(SRCROOT)atmos_spectral/tools/spherical.F90 $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)atmos_spectral/driver/solo/idealized_moist_phys.F90 $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)shared/axis_utils/axis_utils.F90 $(SRCROOT)coupler/surface_flux.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)atmos_spectral/model/fv_advection.F90 $(SRCROOT)atmos_param/sea_esf_rad/null/rad_utilities.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)shared/diag_manager/diag_table.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)atmos_param/diffusivity/diffusivity.F90 $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)shared/tridiagonal/tridiagonal.F90 $(SRCROOT)atmos_solo/atmos_model.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90 $(SRCROOT)atmos_spectral/init/polvani_2007.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)atmos_spectral/model/spectral_damping.F90 /scratch/sit204/old_FMS_model_work/FMS2006/MoistModel/src/atmos_param/monin_obukhov/monin_obukhov.f90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90 $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)atmos_param/entrain/entrain.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)atmos_spectral/model/matrix_invert.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90 $(SRCROOT)atmos_param/my25_turb/my25_turb.F90 $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)atmos_spectral/model/tracer_type.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90 $(SRCROOT)atmos_shared/interpolator/interpolator.F90 $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90 $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)atmos_spectral/model/leapfrog.F90 $(SRCROOT)atmos_spectral/model/global_integral.F90 $(SRCROOT)atmos_spectral/init/topog_regularization.F90 $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)atmos_spectral/model/press_and_geopot.F90 $(SRCROOT)shared/fft/fft99.F90 $(SRCROOT)shared/random_numbers/MersenneTwister.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)atmos_spectral/driver/solo/mixed_layer.F90 $(SRCROOT)atmos_spectral/model/water_borrowing.F90 $(SRCROOT)shared/mpp/test_mpp.F90 $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_reduce.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_a.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_define_nest_domains.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_s.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_sma.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/system_clock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_write.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_s.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nonblock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_reduce.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_nocomm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_misc.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_comm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_misc.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_mpi.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/write_data.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_2d.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_mpi.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_sma.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_mpi.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_scalar.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_4d.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nest.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_nocomm.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_nocomm.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_read.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_sma.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_connect.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_util.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_get_boundary.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write_2Ddecomp.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_checkV.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_nocomm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_sma.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum_tl.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_field.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV_nonblock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/read_mosaic.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_nocomm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_mpi.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_mpi.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D_nonblock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_mpi.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_global_field.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_check.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_sma.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/constant.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_nocomm.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_int.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/field_manager/parse.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_get_boundary.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_3d.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/interp.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_sma.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_util.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/create_xgrid.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_redistribute.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_define.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/gradient_c2l.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_gather.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_a.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_nest_domains.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_read_2Ddecomp.h
OBJ = spectral_dynamics.o mpp_data.o affinity.o gradient_c2l.o atmosphere.o fms_io.o ic_from_external_file.o edt.o time_interp_external.o every_step_diagnostics.o tracer_manager.o mpp_memutils.o strat_cloud.o diag_util.o mosaic_util.o time_interp.o random_numbers.o mosaic.o memuse.o polvani_2004.o memutils.o jablonowski_2006.o threadloc.o horiz_interp_bicubic.o gaussian_topog.o vert_diff.o diag_axis.o implicit.o diag_manager.o qe_moist_convection.o two_stream_gray_rad.o transforms.o vert_coordinate.o mpp_utilities.o spherical.o mpp_pset.o fms.o test_mpp_io.o idealized_moist_phys.o create_xgrid.o interp.o diag_output.o horiz_interp_conserve.o mpp_parameter.o axis_utils.o surface_flux.o platform.o fv_advection.o rad_utilities.o mpp_domains.o diag_table.o diag_grid.o diffusivity.o time_manager.o spherical_fourier.o vert_advection.o tridiagonal.o atmos_model.o mpp_io.o sat_vapor_pres_k.o stable_bl_turb.o polvani_2007.o horiz_interp_spherical.o topography.o spectral_damping.o monin_obukhov.o horiz_interp.o gradient.o spectral_initialize_fields.o fft.o fm_util.o entrain.o horiz_interp_bilinear.o diag_data.o grid.o matrix_invert.o get_cal_time.o read_mosaic.o hs_forcing.o my25_turb.o constants.o tracer_type.o sat_vapor_pres.o test_mpp_pset.o spec_mpp.o test_fms_io.o shallow_conv.o interpolator.o grid_fourier.o field_manager.o vert_turb_driver.o lscale_cond.o gauss_and_legendre.o nsclock.o leapfrog.o global_integral.o topog_regularization.o mpp.o press_and_geopot.o fft99.o MersenneTwister.o horiz_interp_type.o test_mpp_domains.o mixed_layer.o water_borrowing.o test_mpp.o spectral_init_cond.o
OFF = /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_reduce.h $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90 $(SRCROOT)atmos_shared/interpolator/interpolator.F90 $(SRCROOT)atmos_param/entrain/entrain.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D.h $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)atmos_spectral/init/polvani_2004.F90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/write_data.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_2d.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_mpi.h $(SRCROOT)atmos_spectral/init/vert_coordinate.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_sma.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_4d.inc $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)atmos_param/vert_diff/vert_diff.F90 $(SRCROOT)shared/topography/gaussian_topog.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit.inc $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)atmos_spectral/model/global_integral.F90 $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90 $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90 $(SRCROOT)atmos_spectral/init/polvani_2007.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_connect.inc $(SRCROOT)atmos_spectral/model/matrix_invert.F90 $(SRCROOT)shared/random_numbers/random_numbers.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_nocomm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_sma.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_field.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum_tl.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/read_mosaic.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_nocomm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_mpi.inc $(SRCROOT)atmos_spectral/tools/transforms.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_domains2D_nonblock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_mpi.h $(SRCROOT)shared/mpp/mpp_utilities.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_check.h $(SRCROOT)shared/axis_utils/axis_utils.F90 $(SRCROOT)atmos_spectral/driver/solo/mixed_layer.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_int.h $(SRCROOT)shared/time_manager/time_manager.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_get_boundary.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/field_manager/parse.inc $(SRCROOT)atmos_param/edt/edt.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/fms/read_data_3d.inc $(SRCROOT)shared/diag_manager/diag_axis.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/include/fms_platform.h $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_gather.h $(SRCROOT)shared/fft/fft99.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/gradient_c2l.h $(SRCROOT)atmos_param/qe_moist_convection/qe_moist_convection.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_update_nest_domains.h $(SRCROOT)shared/diag_manager/diag_manager.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_a.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_s.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)atmos_spectral/model/water_borrowing.F90 $(SRCROOT)shared/fms/fms.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/mosaic_util.h $(SRCROOT)atmos_spectral/model/leapfrog.F90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90 $(SRCROOT)shared/field_manager/fm_util.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/system_clock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_sma.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_write.inc $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_s_s.h $(SRCROOT)shared/mpp/mpp_io.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nonblock.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_nocomm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_reduce.inc $(SRCROOT)shared/mpp/mpp.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_misc.inc $(SRCROOT)shared/tridiagonal/tridiagonal.F90 $(SRCROOT)shared/memutils/memuse.c /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_comm.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_data_mpi.inc /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_misc.inc $(SRCROOT)atmos_param/my25_turb/my25_turb.F90 $(SRCROOT)atmos_spectral/model/implicit.F90 $(SRCROOT)shared/diag_manager/diag_table.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)atmos_spectral/tools/spherical.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_chksum_scalar.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_update_nest.h $(SRCROOT)atmos_spectral/model/tracer_type.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)atmos_param/sea_esf_rad/null/rad_utilities.F90 $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)shared/mpp/mpp_parameter.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_read.inc $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90 $(SRCROOT)atmos_param/strat_cloud/null/strat_cloud.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/memutils/memutils.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_util_sma.inc $(SRCROOT)shared/random_numbers/MersenneTwister.F90 $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 $(SRCROOT)shared/constants/constants.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_io_util.inc $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_get_boundary.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)shared/mpp/test_mpp.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_checkV.h $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum.inc $(SRCROOT)shared/mpp/affinity.c /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)atmos_param/diffusivity/diffusivity.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm_mpi.inc $(SRCROOT)shared/mosaic/interp.c /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_comm.inc $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_global_field.h $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_sma.h $(SRCROOT)atmos_spectral/init/topog_regularization.F90 $(SRCROOT)atmos_spectral/model/spectral_damping.F90 $(SRCROOT)atmos_spectral/driver/solo/idealized_moist_phys.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/constant.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_sum_nocomm.h $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)coupler/surface_flux.F90 $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_global_sum.h /scratch/sit204/old_FMS_model_work/FMS2006/MoistModel/src/atmos_param/monin_obukhov/monin_obukhov.f90 $(SRCROOT)shared/time_manager/get_cal_time.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/interp.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_reduce_sma.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_util.inc $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)atmos_spectral/model/press_and_geopot.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mosaic/create_xgrid.h /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_do_redistribute.h $(SRCROOT)shared/diag_manager/diag_output.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_domains_define.inc $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90 $(SRCROOT)shared/mosaic/grid.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_error_a_a.h $(SRCROOT)atmos_spectral/model/fv_advection.F90 $(SRCROOT)atmos_solo/atmos_model.F90 $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)shared/time_interp/time_interp_external.F90 /scratch/sit204/FMS2013/GFDLmoistModel/src/shared/mpp/include/mpp_read_2Ddecomp.h
clean: neat
	-rm -f .fms_moist.x.cppdefs $(OBJ) fms_moist.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
fms_moist.x: $(OBJ) 
	$(LD) $(OBJ) -o fms_moist.x  $(LDFLAGS)
