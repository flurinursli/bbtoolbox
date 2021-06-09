include ../../Makefile.inc

LIB_NAME = libfft.a

F90_OBJ = m_fft_real.o m_fft_complex.o

# default rule
fft : dirs $(LIB_NAME)

$(F90_OBJ) : %.o : %.f90
	$(FC) $(OPT) $(OMP) -cpp -DDOUBLE_PREC -I../common/include -I$(INCL_FFTW) -c -o $@ $<

$(LIB_NAME): $(F90_OBJ)
	ar qc $@ $(F90_OBJ)
	mv $(LIB_NAME) lib/
	cp *.mod include/

dirs:
	mkdir -p include lib

clean:
	rm -rf include
	rm -rf lib
	rm -rf *.o
	rm -rf *.mod
