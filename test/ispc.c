#include <stdlib.h>
#include <stdio.h>
#include "../_kiss_fft_guts.h"

#if USE_ISPC
#include "../ispc/test_ispc.h"

void print_debug(kiss_fft_cpx *buf, int m);
kiss_fft_scalar rand_scalar();
void test_aos2();

void print_debug(kiss_fft_cpx *buf, int m) {
    for (int i=0; i<m; ++i) {
        printf("[%f, %f]\n", buf[i].r, buf[i].i);
    }
}

kiss_fft_scalar rand_scalar() {
#if 1
    static kiss_fft_scalar tmp = .0f;
    tmp += 1.0f;
    return tmp;
#else
    return rand();
#endif 
}

void test_aos2() 
{
    int m = 16;
    kiss_fft_cpx *buf = malloc(sizeof(kiss_fft_cpx) * m);
    for (int i=0; i<m; ++i) {
        buf[i].r = rand_scalar();
        buf[i].i = rand_scalar();
    }
    printf("### original ###\n");
    print_debug(buf, m);
    ispc_test_soa2aos2((struct ispc_cpx *)buf, m);
    printf("### test ###\n");
    print_debug(buf, m);
    free(buf);
}

#endif

int main() {
    srand(0);
#if USE_ISPC
    test_aos2();
#else 
    fprintf(stderr, "To test, please build with -DUSE_ISPC=1");
#endif
    return 0;
}