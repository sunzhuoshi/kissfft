/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */
#include "../_kiss_fft_guts.h"
#include "exlib.ispc"

struct ispc_cpx {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
};

#if 0
#define QUALIFIER uniform
#else 
#define QUALIFIER 
#endif

inline void mul(QUALIFIER ispc_cpx &m, const QUALIFIER ispc_cpx &a, const QUALIFIER ispc_cpx &b) {
    m.r = a.r * b.r - a.i * b.i;
    m.i = a.r * b.i + a.i * b.r;
}

inline void sub(QUALIFIER ispc_cpx &res, const QUALIFIER ispc_cpx &a, const QUALIFIER ispc_cpx &b) {
    res.r = a.r - b.r;
    res.i = a.i - b.i;
}

inline void add(QUALIFIER ispc_cpx &res, const QUALIFIER ispc_cpx &a, const QUALIFIER ispc_cpx &b) {
    res.r = a.r + b.r;
    res.i = a.i + b.i;
}

inline void addto(QUALIFIER ispc_cpx &res, const QUALIFIER ispc_cpx &a) {
    res.r += a.r;
    res.i += a.i;
}

inline void sum_sa(QUALIFIER ispc_cpx &res, const QUALIFIER ispc_cpx &a, const QUALIFIER ispc_cpx &b) {
    res.r = a.r - b.i;
    res.i = a.i + b.r;
}

inline void sum_as(QUALIFIER ispc_cpx &res, const QUALIFIER ispc_cpx &a, const QUALIFIER ispc_cpx &b) {
    res.r = a.r + b.i;
    res.i = a.i - b.r;
}

struct ispc_state {
    int nfft;
    int inverse;
    int factors[2*MAXFACTORS];
    ispc_cpx twiddles[1];
};

export void ispc_bfly2(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m
        )
{
#if 1
    ispc_cpx * uniform Fout2 = Fout + m;
    ispc_cpx * uniform tw = st->twiddles;

    uniform int chunk = 0;
    foreach (i = 0 ... m ) {
        ispc_cpx t_i;
        ispc_cpx t_out2;
        ispc_cpx t_out;
        ispc_cpx t;
        
        float *uniform twi = (float *uniform)&tw[chunk * fstride];
        float *uniform F = (float *uniform)&Fout[chunk];
        float *uniform Fm = (float *uniform)&Fout2[chunk];

        aos_to_soa2(twi, &t_i.r, &t_i.i);
        aos_to_soa2(F, &t_out.r, &t_out.i);
        aos_to_soa2(Fm, &t_out2.r, &t_out2.i);
        // TODO: support FIXED_POINT
        // C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);
        mul(t, t_out2, t_i);
        sub(t_out2, t_out, t);
        addto(t_out,  t);

        soa_to_aos2(t_out.r, t_out.i, F);
        soa_to_aos2(t_out2.r, t_out.i, Fm);
        chunk += programCount;
    }
#else
    ispc_cpx * uniform Fout2;
    ispc_cpx * uniform tw1 = st->twiddles;
    uniform ispc_cpx t;
    Fout2 = Fout + m;
    do{
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t, *Fout2, *tw1);
        tw1 += fstride;
        C_SUB(*Fout2,  *Fout, t);
        C_ADDTO(*Fout,  t);
        ++Fout2;
        ++Fout;
    }while (--m);
#endif 
}

export void ispc_bfly4(
        uniform ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        uniform ispc_state * uniform st,
        uniform int m
        )
{
#if 0
    uniform ispc_cpx * uniform tw = st->twiddles;
    ispc_cpx scratch[6];
    const uniform size_t m2=2*m;
    const uniform size_t m3=3*m;

    uniform int Chunk = 0;

    foreach (i = 0 ... m-1) 
    {
        uniform float *uniform F = (uniform float *uniform)&Fout[Chunk];
        uniform float *uniform Fm = (uniform float *uniform)&Fout[Chunk+m];
        uniform float *uniform Fm2 = (uniform float *uniform)&Fout[Chunk+m2];
        uniform float *uniform Fm3 = (uniform float *uniform)&Fout[Chunk+m3];

        uniform float *uniform t1 = (uniform float *uniform)&tw[fstride * Chunk];
        uniform float *uniform t2 = (uniform float *uniform)&tw[fstride * 2 * Chunk];
        uniform float *uniform t3 = (uniform float *uniform)&tw[fstride * 3 * Chunk];

        // TODO: support FIXED_POINT
        //C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);
        ispc_cpx t_out, t_out1, t_out2, t_out3, t_1, t_2, t_3;
        aos_to_soa2(F, &t_out.r, &t_out.i);
        aos_to_soa2(Fm, &t_out1.r, &t_out1.i);
        aos_to_soa2(Fm2, &t_out2.r, &t_out2.i);
        aos_to_soa2(Fm3, &t_out3.r, &t_out3.i);

        aos_to_soa2(t1, &t_1.r, &t_1.i);
        aos_to_soa2(t2, &t_2.r, &t_2.i);
        aos_to_soa2(t3, &t_3.r, &t_3.i);

        mul(scratch[0], t_out1, t_1);
        mul(scratch[1], t_out2, t_2);
        mul(scratch[2], t_out3, t_3);
        sub(scratch[5], t_out, scratch[1]);
        addto(t_out, scratch[1]);
        add(scratch[3], scratch[0], scratch[2]);
        sub(scratch[4], scratch[0], scratch[2]);
        sub(t_out2, t_out, scratch[3]);
        addto(t_out, scratch[3]);

        if(st->inverse) {
            sum_sa(t_out1, scratch[5], scratch[4]);
            sum_as(t_out3, scratch[5], scratch[4]);
        }else{
            sum_as(t_out1, scratch[5], scratch[4]);
            sum_sa(t_out3, scratch[5], scratch[4]);
        }

        soa_to_aos2(t_out.r, t_out.i, F);
        soa_to_aos2(t_out1.r, t_out1.i, Fm);
        soa_to_aos2(t_out2.r, t_out2.i, Fm2);
        soa_to_aos2(t_out3.r, t_out3.i, Fm3);

        Chunk += programCount;
    }
#else 
    /*
    ispc_cpx * uniform tw;
    ispc_cpx uniform scratch[6];
    const uniform size_t m2=2*m;
    const uniform size_t m3=3*m;

    tw = st->twiddles;

    for (uniform int i=0; i<m; ++i) {
        // TODO: support FIXED_POINT
        //C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);
        uniform ispc_cpx t_out = Fout[i];
        uniform ispc_cpx t_out1 = Fout[m + i];
        uniform ispc_cpx t_out2 = Fout[m2 + i];
        uniform ispc_cpx t_out3 = Fout[m3 + i];
        uniform ispc_cpx t_1 = tw[fstride * i];
        uniform ispc_cpx t_2 = tw[fstride * 2 * i];
        uniform ispc_cpx t_3 = tw[fstride * 3 * i];

        mul(scratch[0], t_out1, t_1);
        mul(scratch[1], t_out2, t_2);
        mul(scratch[2], t_out3, t_3);
        sub(scratch[5], t_out, scratch[1]);
        addto(t_out, scratch[1]);
        add(scratch[3], scratch[0], scratch[2]);
        sub(scratch[4], scratch[0], scratch[2]);
        sub(t_out2, t_out, scratch[3]);
        addto(t_out, scratch[3]);

        if(st->inverse) {
            sum_sa(t_out1, scratch[5], scratch[4]);
            sum_as(t_out3, scratch[5], scratch[4]);
        }else{
            sum_as(t_out1, scratch[5], scratch[4]);
            sum_sa(t_out3, scratch[5], scratch[4]);
        }
        Fout[i] = t_out;
        Fout[m + i] = t_out1;
        Fout[m2 + i] = t_out2;
        Fout[m3 + i] = t_out3;
    }
    */
    ispc_cpx * uniform tw1,* uniform tw2,* uniform tw3;
    ispc_cpx uniform scratch[6];
    uniform size_t k=m;
    const uniform size_t m2=2*m;
    const uniform size_t m3=3*m;


    tw3 = tw2 = tw1 = st->twiddles;

    do {
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);

        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );

        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
        ++Fout;
    }while(--k);
#endif
}

export void ispc_bfly5(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m
        )
{
    ispc_cpx * uniform Fout0,* uniform Fout1,* uniform Fout2,* uniform Fout3,* uniform Fout4;
    uniform int u;
    ispc_cpx uniform scratch[13];
    ispc_cpx * uniform twiddles = st->twiddles;
    ispc_cpx * uniform tw;
    uniform ispc_cpx ya,yb;
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