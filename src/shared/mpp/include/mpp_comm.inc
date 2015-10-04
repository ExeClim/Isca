!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -*-f90-*-
! $Id: mpp_comm.inc,v 19.0 2012/01/06 22:00:37 fms Exp $


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

!#################################################
#undef MPP_GATHER_1D_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_GATHER_1D_ mpp_gather_int4_1d
#define MPP_GATHER_1DV_ mpp_gather_int4_1dv
#include <mpp_gather.h>

#undef MPP_GATHER_1D_
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_GATHER_1D_ mpp_gather_real4_1d
#define MPP_GATHER_1DV_ mpp_gather_real4_1dv
#include <mpp_gather.h>

#undef MPP_GATHER_1D_
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_GATHER_1D_ mpp_gather_real8_1d
#define MPP_GATHER_1DV_ mpp_gather_real8_1dv
#include <mpp_gather.h>

