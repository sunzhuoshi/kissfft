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

inline void ispc_bfly2_scalar(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m
        )
{
    ispc_cpx * uniform Fout2 = Fout + m;
    ispc_cpx * uniform tw = st->twiddles;
    uniform ispc_cpx t;

    while (m--) {
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
// it seems that no obvious performance gain vs. scalar ISPC code
// need to be double confirmed with more benchmark
#if 0
    uniform int chunk = 0;
    ispc_cpx * uniform Fout2 = Fout + m;
    ispc_cpx * uniform tw = st->twiddles;

    uniform int ul = m % programCount;;
    uniform int vl = m - ul;

    foreach (i = 0 ... vl ) {
// no better solution, just ignore perf warning
#pragma ignore warning(perf)	        
        ispc_cpx t_i = tw[i * fstride];
        ispc_cpx t_out2;
        ispc_cpx t_out;
        ispc_cpx t;
        
        float *uniform F = (float *uniform)&Fout[chunk];
        float *uniform Fm = (float *uniform)&Fout2[chunk];

        aos_to_soa2(F, &t_out.r, &t_out.i);
        aos_to_soa2(Fm, &t_out2.r, &t_out2.i);
        // TODO: support FIXED_POINT
        // C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);
        mul(t, t_out2, t_i);
        sub(t_out2, t_out, t);
        addto(t_out,  t);

        soa_to_aos2(t_out.r, t_out.i, F);
        soa_to_aos2(t_out2.r, t_out2.i, Fm);
        chunk += programCount;
    }

    // last part that can't vectorized
    ispc_bfly2_scalar(Fout + chunk, fstride, st, m);
#else
    ispc_bfly2_scalar(Fout, fstride, st, m);
#endif 
}

export void ispc_bfly4(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m
        )
{
#if 0
    ispc_cpx * uniform tw = st->twiddles;
    ispc_cpx * uniform Fout_m = Fout + m;
    ispc_cpx * uniform Fout_m2 = Fout + m * 2;
    ispc_cpx * uniform Fout_m3 = Fout + m * 3;
    ispc_cpx scratch[6];

    uniform int chunk = 0;

    foreach (i = 0 ... m) 
    {
        float *uniform F = (float *uniform)&Fout[chunk];
        float *uniform Fm = (float *uniform)&Fout_m[chunk];
        float *uniform Fm2 = (float *uniform)&Fout_m2[chunk];
        float *uniform Fm3 = (float *uniform)&Fout_m3[chunk];

        // TODO: support FIXED_POINT
        //C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);
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
#else 
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

export void ispc_bfly3(
    ispc_cpx * uniform Fout,
    uniform const size_t fstride,
    ispc_state * uniform st,
    uniform int m
) 
{
     uniform size_t k=m;
     const uniform size_t m2 = 2*m;
     ispc_cpx * uniform tw1,* uniform tw2;
     uniform ispc_cpx scratch[5];
     uniform ispc_cpx epi3;
     epi3 = st->twiddles[fstride*m];

     tw1=tw2=st->twiddles;

     do{
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
     }while(--k);    
}

export void ispc_bfly5(
        ispc_cpx * uniform Fout,
        uniform const size_t fstride,
        ispc_state * uniform st,
        uniform int m
        )
{
#if 0
    ispc_cpx * uniform Fout0 = Fout,
        * uniform Fout1 = Fout + m,
        * uniform Fout2 = Fout + m * 2,
        * uniform Fout3 = Fout + m * 3,
        * uniform Fout4 = Fout + m * 4;
    QUALIFIER ispc_cpx scratch[13];
    ispc_cpx * uniform tw = st->twiddles;
    QUALIFIER ispc_cpx ya = tw[fstride * m], 
        yb = tw[fstride * m * 2];

    uniform int chunk = 0;
    foreach (i = 0...m) {
    //for (uniform int i=0; i<m; ++i) {
        QUALIFIER ispc_cpx t_1 = tw[i * fstride], 
            t_2 = tw[i * fstride * 2], 
            t_3 = tw[i * fstride * 3], 
            t_4 = tw[i * fstride * 4];
        QUALIFIER ispc_cpx t_out0, t_out1, t_out2, t_out3, t_out4;

        float *uniform F0 = (float *uniform)&Fout0[chunk];
        float *uniform F1 = (float *uniform)&Fout1[chunk];
        float *uniform F2 = (float *uniform)&Fout2[chunk];
        float *uniform F3 = (float *uniform)&Fout3[chunk];
        float *uniform F4 = (float *uniform)&Fout4[chunk];

        aos_to_soa2(F0, &t_out0.r, &t_out0.i);
        aos_to_soa2(F1, &t_out1.r, &t_out1.i);
        aos_to_soa2(F2, &t_out2.r, &t_out2.i);
        aos_to_soa2(F3, &t_out3.r, &t_out3.i);
        aos_to_soa2(F4, &t_out4.r, &t_out4.i);

/*
        t_out0 = Fout0[i];
        t_out1 = Fout1[i];
        t_out2 = Fout2[i];
        t_out3 = Fout3[i];
        t_out4 = Fout4[i];
*/
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
/*
        Fout0[i] = t_out0;
        Fout1[i] = t_out1;
        Fout2[i] = t_out2;
        Fout3[i] = t_out3;
        Fout4[i] = t_out4;
*/
    }
#else    
    ispc_cpx * uniform Fout0 = Fout,
        * uniform Fout1 = Fout + m,
        * uniform Fout2 = Fout + m * 2,
        * uniform Fout3 = Fout + m * 3,
        * uniform Fout4 = Fout + m * 4;
    ispc_cpx uniform scratch[13];
    ispc_cpx * uniform tw = st->twiddles;
    uniform ispc_cpx ya = tw[fstride * m],
        yb = tw[fstride * m * 2];

    for (uniform int u=0; u<m; ++u ) {
        uniform ispc_cpx t_1 = tw[u * fstride],
            t_2 = tw[u * fstride * 2],
            t_3 = tw[u * fstride * 3],
            t_4 = tw[u * fstride * 4];
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
#endif
}