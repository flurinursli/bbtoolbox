include ../../Makefile.inc

F90_OBJ = hyp_2F1.o quadpack.o m_rtt.o

LIB_NAME = librtt.a

rtt: dirs $(LIB_NAME)

$(F90_OBJ) : %.o : %.f90
	$(FC) $(OPT) -cpp -I../common/include -I../llsq/include64 -I../interpolation/include64 -I../fft/include -c -o $@ $<

$(LIB_NAME) : $(F90_OBJ)
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

.PHONY: clean
