all:
	ispc $(TYPEFLAGS) kiss_fft.ispc -h kiss_fft_ispc.h -o kiss_fft_ispc.o 

asm:
	ispc $(TYPEFLAGS) kiss_fft.ispc -o kiss_fft_ispc.asm --emit-asm

cpp:
	ispc $(TYPEFLAGS) kiss_fft.ispc -o kiss_fft_ispc.cpp --emit-c++ --target=generic-4

clean:
	rm -rf kiss_fft_ispc.h 
	rm -rf *.o
