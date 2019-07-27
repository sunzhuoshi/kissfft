#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <xmmintrin.h>
#include "../_kiss_fft_guts.h"

#if USE_ISPC
#include "../ispc/test_ispc.h"
#include "../ispc/kiss_fft_ispc.h"
#include "orig.h"

#define NFFT 240 //8 * 3 * 5

void print_debug(kiss_fft_cpx *buf, unsigned int m, unsigned int count) {
    for (unsigned int i=0; i<m && i<count; ++i) {
        printf("[%f, %f]\n", buf[i].r, buf[i].i);
    }
}

kiss_fft_scalar rand_scalar() {
#if 0
    static kiss_fft_scalar tmp = .0f;
    tmp += 1.0f;
    return tmp;
#else
    return rand();
#endif 
}

double snr_compare(kiss_fft_cpx * vec1, kiss_fft_cpx * vec2, int n)
{
	int k;
	double sigpow = 1e-10, noisepow = 1e-10, err, snr, scale = 0;

	for (k = 0; k < n; ++k) {
		sigpow += (double)vec1[k].r * (double)vec1[k].r +
			(double)vec1[k].i * (double)vec1[k].i;
		err = (double)vec1[k].r - (double)vec2[k].r;
		noisepow += err * err;
		err = (double)vec1[k].i - (double)vec2[k].i;
		noisepow += err * err;

		if (vec1[k].r)
			scale += (double)vec2[k].r / (double)vec1[k].r;
	}

	snr = 10 * log10(sigpow / noisepow);
	scale /= n;
	if (snr < 10) {
		printf("\npoor snr, try a scaling factor %f\n", scale);
		exit(1);
	}
	return snr;
}

void test_soa2aos2() 
{
    unsigned int m = 16;
    kiss_fft_cpx *buf = _mm_malloc(sizeof(kiss_fft_cpx) * m, 16);
    for (unsigned int i=0; i<m; ++i) {
        buf[i].r = rand_scalar();
        buf[i].i = rand_scalar();
    }
    printf("### original ###\n");
    print_debug(buf, m, m);
    ispc_test_soa2aos2((struct ispc_cpx *)buf, m);
    printf("### test ###\n");
    print_debug(buf, m, m);
    free(buf);
}

void test_bfly2()
{
	kiss_fft_cpx buf_c[NFFT];
	kiss_fft_cpx buf_ispc[NFFT];
	int fstride = 4;
	unsigned int m = 8;

	kiss_fft_cfg state_c = kiss_fft_alloc(NFFT, 0, 0, 0);
	kiss_fft_cfg state_ispc = kiss_fft_alloc(NFFT, 0, 0, 0);
	for (int i = 0; i < NFFT; ++i) {
		buf_c[i].r = buf_ispc[i].r = rand_scalar();
		buf_c[i].i = buf_ispc[i].i = rand_scalar();
	}
	kf_bfly2(buf_c, fstride, state_c, m);
	printf("### C ###\n");
	print_debug(buf_c, NFFT, 4);
	ispc_bfly2((struct ispc_cpx *)buf_ispc, fstride, (struct ispc_state *)state_ispc, m);
	printf("### ISPC ###\n");
	print_debug(buf_ispc, NFFT, 4);
	free(state_c);
	free(state_ispc);
}

void test_bfly4()
{
	kiss_fft_cpx buf_c[NFFT];
	kiss_fft_cpx buf_ispc[NFFT];
	int fstride = 1;
	unsigned int m = 30;

	kiss_fft_cfg state_c = kiss_fft_alloc(NFFT, 0, 0, 0);
	kiss_fft_cfg state_ispc = kiss_fft_alloc(NFFT, 0, 0, 0);
	for (int i = 0; i < NFFT; ++i) {
		buf_c[i].r = buf_ispc[i].r = rand_scalar();
		buf_c[i].i = buf_ispc[i].i = rand_scalar();
	}
	kf_bfly4(buf_c, fstride, state_c, m);
	printf("### C ###\n");
	print_debug(buf_c, NFFT, 16);
	ispc_bfly4((struct ispc_cpx *)buf_ispc, fstride, (struct ispc_state *)state_ispc, m);
	printf("### ISPC ###\n");
	print_debug(buf_ispc, NFFT, 16);
	free(state_c);
	free(state_ispc);
}

#endif

int main() {
    srand(0);
#if USE_ISPC
	//test_bfly2();
	test_bfly4();
	//test_soa2aos2();
#else 
    fprintf(stderr, "To test, please build with -DUSE_ISPC=1");
#endif
    return 0;
}