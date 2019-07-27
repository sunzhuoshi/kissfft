/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */
#pragma once
#include "../kiss_fft_def.h"
#include <stddef.h>

kiss_fft_cfg kiss_fft_alloc(int nfft, int inverse_fft, void * mem, size_t * lenmem);
void kf_bfly2(kiss_fft_cpx * Fout, const size_t fstride, const kiss_fft_cfg st, const size_t m);
void kf_bfly4(kiss_fft_cpx * Fout, const size_t fstride, const kiss_fft_cfg st, const size_t m);
