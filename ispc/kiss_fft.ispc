#include "../_kiss_fft_guts.h"

struct ispc_cpx {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
};

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
    ispc_cpx * uniform Fout2;
    ispc_cpx * uniform tw1 = st->twiddles;
    uniform ispc_cpx t;
    Fout2 = Fout + m;
    do{
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t,  *Fout2 , *tw1);
        tw1 += fstride;
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
        ++Fout2;
        ++Fout;
    }while (--m);
}