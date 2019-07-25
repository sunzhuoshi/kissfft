#include "../_kiss_fft_guts.h"
#include "exlib.ispc"

struct ispc_cpx {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
};

export void ispc_test_soa2aos2(
    ispc_cpx * uniform buf,
    uniform int m
) {
    uniform int chunk = 0;
    foreach (i=0 ... m-1) {
        ispc_cpx t;
        float *uniform p = (float *uniform)&buf[chunk];
        aos_to_soa2(p, &t.r, &t.i);
        soa_to_aos2(t.r, t.i, p);
    }
}