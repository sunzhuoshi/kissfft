#pragma once
/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

#ifdef USE_SIMD
# include <xmmintrin.h>
# define kiss_fft_scalar __m128
#define KISS_FFT_MALLOC(nbytes) _mm_malloc(nbytes,16)
#define KISS_FFT_FREE _mm_free
#else	
#define KISS_FFT_MALLOC malloc
#define KISS_FFT_FREE free
#endif	


#ifdef FIXED_POINT
#ifdef ISPC
#include "ispc/stdint.h"
#else
#include <stdint.h>
#endif
# if (FIXED_POINT == 32)
#  define kiss_fft_scalar int32_t
# else	
#  define kiss_fft_scalar int16_t
# endif
#else
# ifndef kiss_fft_scalar
/*  default is float */
#   define kiss_fft_scalar float
# endif
#endif

typedef struct {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
} kiss_fft_cpx;

typedef struct kiss_fft_state* kiss_fft_cfg;

