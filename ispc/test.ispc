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
    int nfft;
    int inverse;
    int factors[2*MAXFACTORS];
    ispc_cpx twiddles[1];
};

export void ispc_test_soa2aos2(
    ispc_cpx * uniform buf,
    uniform int m
) {
    uniform int chunk = 0;
    ispc_cpx t;

    foreach (i=0 ... m) {
        float *uniform p = (float *uniform)&buf[chunk];
        aos_to_soa2(p, &t.r, &t.i);
        soa_to_aos2(t.r, t.i, p);
		chunk += programCount;
    }
}
