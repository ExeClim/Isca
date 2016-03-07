
! -*-f90-*-
! $Id: mpp_domains_reduce.inc,v 18.0 2010/03/02 23:57:19 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_REDUCE: get global max/min of field                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r8_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r8_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r8_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ maxval
#undef REDUCE_LOC_
#define REDUCE_LOC_ maxloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r8_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r8_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r8_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ minval
#undef REDUCE_LOC_
#define REDUCE_LOC_ minloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

#ifdef OVERLOAD_R4
#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r4_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r4_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r4_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ maxval
#undef REDUCE_LOC_
#define REDUCE_LOC_ maxloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r4_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r4_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r4_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ minval
#undef REDUCE_LOC_
#define REDUCE_LOC_ minloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#ifndef no_8byte_integers
#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i8_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i8_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i8_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ maxval
#undef REDUCE_LOC_
#define REDUCE_LOC_ maxloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i8_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i8_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i8_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ minval
#undef REDUCE_LOC_
#define REDUCE_LOC_ minloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i4_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i4_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i4_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ maxval
#undef REDUCE_LOC_
#define REDUCE_LOC_ maxloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i4_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i4_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i4_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ minval
#undef REDUCE_LOC_
#define REDUCE_LOC_ minloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                   MPP_GLOBAL_SUM: global sum of field                       !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#ifdef OVERLOAD_C8
#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>
#endif

#ifdef OVERLOAD_R4
#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>
#endif

#ifdef OVERLOAD_C4
#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>
#endif

#ifndef no_8byte_integers
#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>
#endif

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>
    
!gag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   MPP_GLOBAL_SUM_TL: global sum of forward and tangent-linear fields        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_r8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_r8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_r8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_r8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum_tl.h>

#ifdef OVERLOAD_C8
#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_c8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_c8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_c8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_c8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum_tl.h>
#endif

#ifdef OVERLOAD_R4
#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_r4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_r4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_r4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_r4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum_tl.h>
#endif

#ifdef OVERLOAD_C4
#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_c4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_c4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_c4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_c4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum_tl.h>
#endif

#ifndef no_8byte_integers
#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_i8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_i8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_i8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_i8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum_tl.h>
#endif

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_i4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_i4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_i4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum_tl.h>

#undef MPP_GLOBAL_SUM_TL_
#define MPP_GLOBAL_SUM_TL_ mpp_global_sum_tl_i4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum_tl.h>
!gag

!bnc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!         MPP_GLOBAL_SUM_AD: global adjoint sum of field                      !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_r8_2d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(DOUBLE_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_r8_3d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(DOUBLE_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_r8_4d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(DOUBLE_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_r8_5d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(DOUBLE_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#ifdef OVERLOAD_C8
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_c8_2d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(DOUBLE_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_c8_3d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(DOUBLE_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_c8_4d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(DOUBLE_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_c8_5d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(DOUBLE_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$#endif
!!$
!!$#ifdef OVERLOAD_R4
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_r4_2d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(FLOAT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_r4_3d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(FLOAT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_r4_4d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(FLOAT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_r4_5d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ real(FLOAT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$#endif
!!$
!!$#ifdef OVERLOAD_C4
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_c4_2d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(FLOAT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_c4_3d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(FLOAT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_c4_4d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(FLOAT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_c4_5d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ complex(FLOAT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$#endif
!!$
!!$#ifndef no_8byte_integers
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_i8_2d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(LONG_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_i8_3d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(LONG_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_i8_4d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(LONG_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_i8_5d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(LONG_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$#endif
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_i4_2d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(INT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_i4_3d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(INT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_i4_4d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(INT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$
!!$#undef MPP_GLOBAL_SUM_AD_
!!$#define MPP_GLOBAL_SUM_AD_ mpp_global_sum_ad_i4_5d
!!$#undef MPP_EXTRA_INDICES_
!!$#define MPP_EXTRA_INDICES_ ,:,:,:
!!$#undef MPP_TYPE_
!!$#define MPP_TYPE_ integer(INT_KIND)
!!$#include <mpp_global_sum_ad.h>
!!$!bnc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_FIELD: get global field from domain field           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r8_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r8_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r8_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_field.h>

#ifdef OVERLOAD_C8
#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c8_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c8_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c8_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_field.h>
#endif

#ifndef no_8byte_integers
#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i8_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i8_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i8_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_field.h>

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l8_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l8_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l8_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_global_field.h>
#endif

#ifdef OVERLOAD_R4
#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r4_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r4_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r4_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_field.h>
#endif

#ifdef OVERLOAD_C4
#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c4_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c4_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c4_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_field.h>
#endif

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i4_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i4_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i4_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_field.h>

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l4_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l4_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l4_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_global_field.h>

!****************************************************
#undef MPP_DO_GLOBAL_FIELD_3D_
#define MPP_DO_GLOBAL_FIELD_3D_ mpp_do_global_field2D_r8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_do_global_field.h>                                    

#ifdef OVERLOAD_C8
#undef MPP_DO_GLOBAL_FIELD_3D_
#define MPP_DO_GLOBAL_FIELD_3D_ mpp_do_global_field2D_c8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_do_global_field.h>                                    
#endif

#ifndef no_8byte_integers
#undef MPP_DO_GLOBAL_FIELD_3D_
#define MPP_DO_GLOBAL_FIELD_3D_ mpp_do_global_field2D_i8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_do_global_field.h>                                    

#undef MPP_DO_GLOBAL_FIELD_3D_
#define MPP_DO_GLOBAL_FIELD_3D_ mpp_do_global_field2D_l8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_do_global_field.h>                                    
#endif

#ifdef OVERLOAD_R4
#undef MPP_DO_GLOBAL_FIELD_3D_
#define MPP_DO_GLOBAL_FIELD_3D_ mpp_do_global_field2D_r4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_do_global_field.h>                                    
#endif

#ifdef OVERLOAD_C4
#undef MPP_DO_GLOBAL_FIELD_3D_
#define MPP_DO_GLOBAL_FIELD_3D_ mpp_do_global_field2D_c4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_do_global_field.h>                                    
#endif

#undef MPP_DO_GLOBAL_FIELD_3D_
#define MPP_DO_GLOBAL_FIELD_3D_ mpp_do_global_field2D_i4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_do_global_field.h>                                    

#undef MPP_DO_GLOBAL_FIELD_3D_
#define MPP_DO_GLOBAL_FIELD_3D_ mpp_do_global_field2D_l4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_do_global_field.h>

