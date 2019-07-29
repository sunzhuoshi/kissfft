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

inline void mul(ispc_cpx &m, const ispc_cpx &a, const ispc_cpx &b) {
    m.r = a.r * b.r - a.i * b.i;
    m.i = a.r * b.i + a.i * b.r;
}

inline void mul(ispc_cpx &m, kiss_fft_scalar s) {
    m.r *= s;
    m.i *= s;
}

inline void sub(ispc_cpx &res, const ispc_cpx &a, const ispc_cpx &b) {
    res.r = a.r - b.r;
    res.i = a.i - b.i;
}

inline void add(ispc_cpx &res, const ispc_cpx &a, const ispc_cpx &b) {
    res.r = a.r + b.r;
    res.i = a.i + b.i;
}

inline void addto(ispc_cpx &res, const ispc_cpx &a) {
    res.r += a.r;
    res.i += a.i;
}

inline void sum_sa(ispc_cpx &res, const ispc_cpx &a, const ispc_cpx &b) {
    res.r = a.r - b.i;
    res.i = a.i + b.r;
}

inline void sum_as(ispc_cpx &res, const ispc_cpx &a, const ispc_cpx &b) {
    res.r = a.r + b.i;
    res.i = a.i - b.r;
}

struct ispc_state {
    uniform int nfft;
    uniform int inverse;
    uniform int factors[2*MAXFACTORS];
    uniform ispc_cpx twiddles[1];
};

inline void ispc_bfly2_scalar(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m,
        uniform int offset
        )
{
    ispc_cpx * uniform Fout2 = Fout + m;
    ispc_cpx * uniform tw = st->twiddles + offset;
    uniform ispc_cpx t;
    uniform int l = m - offset;

    while (l--) {
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t, *Fout2, *tw);
        tw += fstride;
        C_SUB(*Fout2,  *Fout, t);
        C_ADDTO(*Fout,  t);
        ++Fout2;
        ++Fout;
    }    
}

export void ispc_bfly2(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m
        )
{
#if USE_ISPC == 1
    uniform int chunk = 0;
    ispc_cpx * uniform Fout2 = Fout + m;
    ispc_cpx * uniform tw = st->twiddles;

    uniform int vl = m - m % programCount;

    foreach (i = 0 ... vl ) {
    // no better solution, just ignore perf warning
#pragma ignore warning(perf)	        
        ispc_cpx t_i = tw[i * fstride];
        ispc_cpx t_out2;
        ispc_cpx t_out;
        ispc_cpx t;
        
        kiss_fft_scalar *uniform F = (kiss_fft_scalar *uniform)&Fout[chunk];
        kiss_fft_scalar *uniform Fm = (kiss_fft_scalar *uniform)&Fout2[chunk];
        // TODO: support FIXED_POINT
        //C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        aos_to_soa2(F, &t_out.r, &t_out.i);
        aos_to_soa2(Fm, &t_out2.r, &t_out2.i);
        mul(t, t_out2, t_i);
        sub(t_out2, t_out, t);
        addto(t_out,  t);

        soa_to_aos2(t_out.r, t_out.i, F);
        soa_to_aos2(t_out2.r, t_out2.i, Fm);
        chunk += programCount;
    }
    // rear part that can't vectorized
    ispc_bfly2_scalar(Fout + chunk, fstride, st, m, vl);
#elif USE_ISPC == 2
    ispc_bfly2_scalar(Fout, fstride, st, m, 0);
#endif 
}

inline void ispc_bfly4_scalar(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m,
        uniform int offset
        )
{
    ispc_cpx * uniform tw1 = st->twiddles + offset
        ,* uniform tw2 = st->twiddles + offset * 2
        ,* uniform tw3 = st->twiddles + offset * 3;
    uniform ispc_cpx uniform scratch[6];
    const uniform size_t m2=2*m;
    const uniform size_t m3=3*m;
    uniform int l = m - offset;

    while (l--) {
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
    }
}

export void ispc_bfly4(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m
        )
{
#if USE_ISPC == 1
    ispc_cpx * uniform tw = st->twiddles;
    ispc_cpx * uniform Fout_m = Fout + m;
    ispc_cpx * uniform Fout_m2 = Fout + m * 2;
    ispc_cpx * uniform Fout_m3 = Fout + m * 3;
    ispc_cpx scratch[6];

    uniform int chunk = 0;
    uniform int vl = m - m % programCount;

    foreach (i = 0 ... vl) {
        kiss_fft_scalar *uniform F = (kiss_fft_scalar *uniform)&Fout[chunk];
        kiss_fft_scalar *uniform Fm = (kiss_fft_scalar *uniform)&Fout_m[chunk];
        kiss_fft_scalar *uniform Fm2 = (kiss_fft_scalar *uniform)&Fout_m2[chunk];
        kiss_fft_scalar *uniform Fm3 = (kiss_fft_scalar *uniform)&Fout_m3[chunk];

        // TODO: support FIXED_POINT
        //C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m*2],4); C_FIXDIV(Fout[m*3],4);
        ispc_cpx t_out, t_out1, t_out2, t_out3;
        // no better solution, just ignore perf warning
#pragma ignore warning(perf)	        
        ispc_cpx t_1 = tw[i * fstride], t_2 = tw[i * fstride * 2], t_3 = tw[i * fstride * 3];

        aos_to_soa2(F, &t_out.r, &t_out.i);
        aos_to_soa2(Fm, &t_out1.r, &t_out1.i);
        aos_to_soa2(Fm2, &t_out2.r, &t_out2.i);
        aos_to_soa2(Fm3, &t_out3.r, &t_out3.i);

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

        chunk += programCount;
    }
    ispc_bfly4_scalar(Fout + chunk, fstride, st, m, vl);
#elif USE_ISPC == 2
    ispc_bfly4_scalar(Fout, fstride, st, m, 0);
#endif
}

export void ispc_bfly3_scalar(
    ispc_cpx * uniform Fout,
    uniform const size_t fstride,
    ispc_state * uniform st,
    uniform int m,
    uniform int offset
) 
{
     const uniform size_t m2 = 2*m;
     ispc_cpx * uniform tw1 = st->twiddles + offset,
        * uniform tw2 = st->twiddles + offset * 2;
     uniform ispc_cpx scratch[5];
     uniform ispc_cpx epi3;
     uniform int l = m - offset;
     epi3 = st->twiddles[fstride*m];

     while (l--) {
         C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);

         C_MUL(scratch[1],Fout[m] , *tw1);
         C_MUL(scratch[2],Fout[m2] , *tw2);

         C_ADD(scratch[3],scratch[1],scratch[2]);
         C_SUB(scratch[0],scratch[1],scratch[2]);
         tw1 += fstride;
         tw2 += fstride*2;

         Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
         Fout[m].i = Fout->i - HALF_OF(scratch[3].i);

         C_MULBYSCALAR( scratch[0] , epi3.i );

         C_ADDTO(*Fout,scratch[3]);

         Fout[m2].r = Fout[m].r + scratch[0].i;
         Fout[m2].i = Fout[m].i - scratch[0].r;

         Fout[m].r -= scratch[0].i;
         Fout[m].i += scratch[0].r;

         ++Fout;
     }
}

export void ispc_bfly3(
    ispc_cpx * uniform Fout,
    uniform const size_t fstride,
    ispc_state * uniform st,
    uniform int m
) 
{
#if 0
    // TODO: add ISPC implementation
#else
    ispc_bfly3_scalar(Fout, fstride, st, m, 0);
#endif
}

inline void ispc_bfly5_scalar(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m,
        uniform int offset
        )
{
    ispc_cpx * uniform Fout0 = Fout,
        * uniform Fout1 = Fout + m,
        * uniform Fout2 = Fout + m * 2,
        * uniform Fout3 = Fout + m * 3,
        * uniform Fout4 = Fout + m * 4;
    uniform ispc_cpx uniform scratch[13];
    ispc_cpx * uniform tw = st->twiddles;
    uniform ispc_cpx ya = tw[fstride * m],
        yb = tw[fstride * m * 2];
    uniform int l = m - offset;

    for (uniform int u=0; u<l; ++u ) {
        uniform ispc_cpx t_1 = tw[(offset + u) * fstride],
            t_2 = tw[(offset + u) * fstride * 2],
            t_3 = tw[(offset + u) * fstride * 3],
            t_4 = tw[(offset + u) * fstride * 4];
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, t_1);
        C_MUL(scratch[2] ,*Fout2, t_2);
        C_MUL(scratch[3] ,*Fout3, t_3);
        C_MUL(scratch[4] ,*Fout4, t_4);

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

export void ispc_bfly5(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m
        )
{
#if USE_ISPC == 1
    ispc_cpx * uniform Fout0 = Fout,
        * uniform Fout1 = Fout + m,
        * uniform Fout2 = Fout + m * 2,
        * uniform Fout3 = Fout + m * 3,
        * uniform Fout4 = Fout + m * 4;
    ispc_cpx scratch[13];
    ispc_cpx * uniform tw = st->twiddles;
    ispc_cpx ya = tw[fstride * m], 
        yb = tw[fstride * m * 2];

    uniform int chunk = 0;
    uniform int ul = m % programCount;;
    uniform int vl = m - ul;

    foreach (i = 0 ... vl ) {
#pragma ignore warning(perf)	        
        ispc_cpx t_1 = tw[i * fstride], t_2 = tw[i * fstride * 2], t_3 = tw[i * fstride * 3], t_4 = tw[i * fstride * 4];
        ispc_cpx t_out0, t_out1, t_out2, t_out3, t_out4;

        kiss_fft_scalar *uniform F0 = (kiss_fft_scalar *uniform)&Fout0[chunk];
        kiss_fft_scalar *uniform F1 = (kiss_fft_scalar *uniform)&Fout1[chunk];
        kiss_fft_scalar *uniform F2 = (kiss_fft_scalar *uniform)&Fout2[chunk];
        kiss_fft_scalar *uniform F3 = (kiss_fft_scalar *uniform)&Fout3[chunk];
        kiss_fft_scalar *uniform F4 = (kiss_fft_scalar *uniform)&Fout4[chunk];

        aos_to_soa2(F0, &t_out0.r, &t_out0.i);
        aos_to_soa2(F1, &t_out1.r, &t_out1.i);
        aos_to_soa2(F2, &t_out2.r, &t_out2.i);
        aos_to_soa2(F3, &t_out3.r, &t_out3.i);
        aos_to_soa2(F4, &t_out4.r, &t_out4.i);

        scratch[0] = t_out0;

        mul(scratch[1], t_out1, t_1);
        mul(scratch[2], t_out2, t_2);
        mul(scratch[3], t_out3, t_3);
        mul(scratch[4], t_out4, t_4);

        add(scratch[7], scratch[1], scratch[4]);
        sub(scratch[10], scratch[1], scratch[4]);
        add(scratch[8], scratch[2], scratch[3]);
        sub(scratch[9], scratch[2], scratch[3]);

        t_out0.r += scratch[7].r + scratch[8].r;
        t_out0.i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + scratch[7].r * ya.r + scratch[8].r * yb.r;
        scratch[5].i = scratch[0].i + scratch[7].i * ya.r + scratch[8].i * yb.r;

        scratch[6].r = scratch[10].i * ya.i + scratch[9].i * yb.i;
        scratch[6].i = - scratch[10].r * ya.i - scratch[9].r * yb.i;

        sub(t_out1, scratch[5], scratch[6]);
        add(t_out4, scratch[5], scratch[6]);

        scratch[11].r = scratch[0].r + scratch[7].r * yb.r + scratch[8].r * ya.r;
        scratch[11].i = scratch[0].i + scratch[7].i * yb.r + scratch[8].i * ya.r;
        scratch[12].r = - scratch[10].i * yb.i + scratch[9].i * ya.i;
        scratch[12].i = scratch[10].r * yb.i - scratch[9].r * ya.i;

        add(t_out2, scratch[11], scratch[12]);
        sub(t_out3, scratch[11], scratch[12]);

        soa_to_aos2(t_out0.r, t_out0.i, F0);
        soa_to_aos2(t_out1.r, t_out1.i, F1);
        soa_to_aos2(t_out2.r, t_out2.i, F2);
        soa_to_aos2(t_out3.r, t_out3.i, F3);
        soa_to_aos2(t_out4.r, t_out4.i, F4);

        chunk += programCount;
    }
    ispc_bfly5_scalar(Fout + chunk, fstride, st, m, vl);
#elif USE_ISPC == 2    
    ispc_bfly5_scalar(Fout, fstride, st, m, 0);
#endif
}