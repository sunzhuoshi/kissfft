# use "make testall USE_ISPC=1" to test ISPC version(experimental)
KFVER=131
TYPEFLAGS=-Dkiss_fft_scalar=float

ifeq "$(USE_ISPC)" ""
 USE_ISPC=0
endif

ISPC_MAKE=@echo ""
ISPC_OBJS=
MAKE_SIMD=make -C test DATATYPE=simd CFLAGADD="$(CFLAGADD)" test
ifeq ($(USE_ISPC), 1)
 ISPC_MAKE=make -C ispc TYPEFLAGS=$(TYPEFLAGS)
 ISPC_OBJS=ispc/kiss_fft_ispc.o
 MAKE_SIMD=@echo "NOTE: no simd datatype support for ISPC"
endif
	
ifeq ($(shell uname -s),Darwin)
	SHARED := -Wl,-install_name,libkissfft.dylib -o libkissfft.dylib
else
	SHARED := -Wl,-soname,libkissfft.so -o libkissfft.so
endif

all:
	gcc -Wall -fPIC -c *.c $(TYPEFLAGS) -o kiss_fft.o 
	$(ISPC_MAKE)
	ar crus libkissfft.a kiss_fft.o $(ISPC_OBJS)
	gcc -shared $(SHARED) kiss_fft.o $(ISPC_OBJS)

install: all
	cp libkissfft.so /usr/local/lib/

doc:
	@echo "Start by reading the README file.  If you want to build and test lots of stuff, do a 'make testall'"
	@echo "but be aware that 'make testall' has dependencies that the basic kissfft software does not."
	@echo "It is generally unneeded to run these tests yourself, unless you plan on changing the inner workings"
	@echo "of kissfft and would like to make use of its regression tests."

testall:
	# The simd and int32_t types may or may not work on your machine
	# USE_ISPC: $(USE_ISPC)
	make -C test testcpp && test/testcpp
	$(MAKE_SIMD)
	make -C test DATATYPE=int32_t CFLAGADD="$(CFLAGADD)" test
	make -C test DATATYPE=int16_t CFLAGADD="$(CFLAGADD)" test
	make -C test DATATYPE=float CFLAGADD="$(CFLAGADD)" test
	make -C test DATATYPE=double CFLAGADD="$(CFLAGADD)" test
	echo "all tests passed"

tarball: clean
	git archive --prefix=kissfft/ -o kissfft$(KFVER).tar.gz v$(KFVER)
	git archive --prefix=kissfft/ -o kissfft$(KFVER).zip v$(KFVER)

clean:
	cd test && make clean
	cd tools && make clean
	make -C ispc clean
	rm -f kiss_fft*.tar.gz *~ *.pyc kiss_fft*.zip 

asm: kiss_fft.s

kiss_fft.s: kiss_fft.c kiss_fft.h _kiss_fft_guts.h
	[ -e kiss_fft.s ] && mv kiss_fft.s kiss_fft.s~ || true
	gcc -S kiss_fft.c -O3 -mtune=native -ffast-math -fomit-frame-pointer -unroll-loops -dA -fverbose-asm 
	gcc -o kiss_fft_short.s -S kiss_fft.c -O3 -mtune=native -ffast-math -fomit-frame-pointer -dA -fverbose-asm -DFIXED_POINT
	[ -e kiss_fft.s~ ] && diff kiss_fft.s~ kiss_fft.s || true
