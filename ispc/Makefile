TARGET=--target=avx2-i32x16
#TARGET=--target=avx2-i64x4
#TARGET=--target=sse4-i32x8
#TARGET=--target=sse4-i32x4

all:
	ispc $(TARGET) $(TYPEFLAGS) kiss_fft.ispc -h kiss_fft_ispc.h -o kiss_fft_ispc.o 

asm:
	ispc $(TARGET) $(TYPEFLAGS) kiss_fft.ispc -o kiss_fft_ispc.asm --emit-asm

cpp:
	ispc $(TARGET) $(TYPEFLAGS) kiss_fft.ispc -o kiss_fft_ispc.cpp --emit-c++ --target=generic-4

clean:
	rm -rf kiss_fft_ispc.h 
	rm -rf *.o
