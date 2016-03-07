
! -*-f90-*-
! $Id: mpp_comm.inc,v 16.0 2008/07/30 22:46:18 fms Exp $


#ifdef use_libSMA
#include <mpp_comm_sma.inc>
#elif defined(use_libMPI)
#include <mpp_comm_mpi.inc>
#else
#include <mpp_comm_nocomm.inc>
#endif

#ifndef no_8byte_integers
#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_1d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_2d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_4d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:,:,:)
#include <mpp_chksum_int.h>
#endif

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_1d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_2d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_4d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_INT_
#define MPP_CHKSUM_INT_ mpp_chksum_i4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_RANK_
#define MPP_RANK_  (:,:,:,:,:)
#include <mpp_chksum_int.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_0d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ !
#include <mpp_chksum_scalar.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_1d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_2d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_4d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#ifdef OVERLOAD_C8
#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_0d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ !
#include <mpp_chksum_scalar.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_1d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_2d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_4d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>
#endif

#ifdef OVERLOAD_R4
#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_0d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ !
#include <mpp_chksum_scalar.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_1d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_2d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_4d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_r4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>
#endif

#ifdef OVERLOAD_C4
#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_0d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ !
#include <mpp_chksum_scalar.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_1d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_2d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_4d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#undef MPP_CHKSUM_
#define MPP_CHKSUM_ mpp_chksum_c4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_RANK_
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>
#endif


