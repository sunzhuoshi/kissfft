/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */
#include "orig.h"
#include "../_kiss_fft_guts.h"
#include <math.h>
#include <stdlib.h>

/*
 * original copies of ISPC methods
 */

 /*  facbuf is populated by p1,m1,p2,m2, ...
	 where
	 p[i] * m[i] = m[i-1]
	 m0 = n                  */
static
void kf_factor(int n, int * facbuf)
{
	int p = 4;
	double floor_sqrt;
	floor_sqrt = floor(sqrt((double)n));

	/*factor out powers of 4, powers of 2, then any remaining primes */
	do {
		while (n % p) {
			switch (p) {
			case 4: p = 2; break;
			case 2: p = 3; break;
			default: p += 2; break;
			}
			if (p > floor_sqrt)
				p = n;          /* no more factors, skip to end */
		}
		n /= p;
		*facbuf++ = p;
		*facbuf++ = n;
	} while (n > 1);
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
kiss_fft_cfg kiss_fft_alloc(int nfft, int inverse_fft, void * mem, size_t * lenmem)
{
	kiss_fft_cfg st = NULL;
	size_t memneeded = sizeof(struct kiss_fft_state)
		+ sizeof(kiss_fft_cpx)*(nfft - 1); /* twiddle factors*/

	if (lenmem == NULL) {
		st = (kiss_fft_cfg)KISS_FFT_MALLOC(memneeded);
	}
	else {
		if (mem != NULL && *lenmem >= memneeded)
			st = (kiss_fft_cfg)mem;
		*lenmem = memneeded;
	}
	if (st) {
		int i;
		st->nfft = nfft;
		st->inverse = inverse_fft;

		for (i = 0; i < nfft; ++i) {
			const double pi = 3.141592653589793238462643383279502884197169399375105820974944;
			double phase = -2 * pi*i / nfft;
			if (st->inverse)
				phase *= -1;
			kf_cexp(st->twiddles + i, phase);
		}

		kf_factor(nfft, st->factors);
	}
	return st;
}

void kf_bfly2(
	kiss_fft_cpx * Fout,
	const size_t fstride,
	const kiss_fft_cfg st,
	const size_t m
)
{
	kiss_fft_cpx * Fout2;
	kiss_fft_cpx * tw1 = st->twiddles;
	kiss_fft_cpx t;
	Fout2 = Fout + m;
	size_t k = m;
	do {
		C_FIXDIV(*Fout, 2); C_FIXDIV(*Fout2, 2);

		C_MUL(t, *Fout2, *tw1);
		tw1 += fstride;
		C_SUB(*Fout2, *Fout, t);
		C_ADDTO(*Fout, t);
		++Fout2;
		++Fout;
	} while (--k);
}

void kf_bfly4(
	kiss_fft_cpx * Fout,
	const size_t fstride,
	const kiss_fft_cfg st,
	const size_t m
)
{
	kiss_fft_cpx *tw1, *tw2, *tw3;
	kiss_fft_cpx scratch[6];
	size_t k = m;
	const size_t m2 = 2 * m;
	const size_t m3 = 3 * m;


	tw3 = tw2 = tw1 = st->twiddles;

	do {
		C_FIXDIV(*Fout, 4); C_FIXDIV(Fout[m], 4); C_FIXDIV(Fout[m2], 4); C_FIXDIV(Fout[m3], 4);

		C_MUL(scratch[0], Fout[m], *tw1);
		C_MUL(scratch[1], Fout[m2], *tw2);
		C_MUL(scratch[2], Fout[m3], *tw3);

		C_SUB(scratch[5], *Fout, scratch[1]);
		C_ADDTO(*Fout, scratch[1]);
		C_ADD(scratch[3], scratch[0], scratch[2]);
		C_SUB(scratch[4], scratch[0], scratch[2]);
		C_SUB(Fout[m2], *Fout, scratch[3]);
		tw1 += fstride;
		tw2 += fstride * 2;
		tw3 += fstride * 3;
		C_ADDTO(*Fout, scratch[3]);

		if (st->inverse) {
			Fout[m].r = scratch[5].r - scratch[4].i;
			Fout[m].i = scratch[5].i + scratch[4].r;
			Fout[m3].r = scratch[5].r + scratch[4].i;
			Fout[m3].i = scratch[5].i - scratch[4].r;
		}
		else {
			Fout[m].r = scratch[5].r + scratch[4].i;
			Fout[m].i = scratch[5].i - scratch[4].r;
			Fout[m3].r = scratch[5].r - scratch[4].i;
			Fout[m3].i = scratch[5].i + scratch[4].r;
		}
		++Fout;
	} while (--k);
}

void kf_bfly5(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        const size_t m
        )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    size_t u;
    kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);

        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r);

        scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
        scratch[6].i = -S_MUL(scratch[10].r,ya.i) - S_MUL(scratch[9].r,yb.i);

        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);

        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r);
        scratch[12].r = - S_MUL(scratch[10].i,yb.i) + S_MUL(scratch[9].i,ya.i);
        scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);

        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}