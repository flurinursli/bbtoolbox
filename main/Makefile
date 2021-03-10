include ../Makefile.inc

PP_FLAGS =

ifeq ($(ERROR_TRAP),y)
	PP_FLAGS += -DERROR_TRAP
endif

ifeq ($(PERF),y)
	PP_FLAGS += -DPERF
endif

ifeq ($(MKL),y)
	PP_FLAGS += -DMKL
endif

ifeq ($(MPI),y)
	PP_FLAGS += -DMPI
endif

OBJ = m_toolbox.o m_roughness.o m_source.o main.o

#ODG = debug.o

main : $(OBJ)
	$(FC) $(OPT) $(OMP) -o ../toolbox.exe $(OBJ) -L../lib/common/lib -L../lib/llsq/lib -L../lib/fft/lib -L../lib/interpolation/lib -L../lib/filter/lib -L../lib/rtt/lib -L$(FFTW_PATH)/lib -L$(LAPACK_PATH) -L$(GSL_PATH) -lcommon -lllsqf -lfilterf -lfft -linterplf -lllsq -lrtt $(LINK_FLAGS) -lgsl -lgslcblas -linterpl

#dbg : $(ODG)
#	$(FC) $(OPT) $(OMP) -o ../debug.exe $(ODG) -L../lib/common/lib -L../lib/llsq/lib -L../lib/fft/lib -L../lib/interpolation/lib -L../lib/rtt/lib -L$(FFTW_PATH)/lib -L$(LAPACK_PATH) -L$(GSL_PATH) -lrtt -lfft -lllsq -lcommon $(LINK_FLAGS) -lgsl -lgslcblas -linterpl

# how to get main.o
$(OBJ) : %.o : %.f90
	$(FC) $(OPT) $(OMP) -cpp $(PP_FLAGS) -I../lib/common/include -I../lib/llsq/include32 -I../lib/interpolation/include32 -I../lib/fft/include -I../lib/filter/include -I../lib/llsq/include32 -I../lib/rtt/include -c -o $@ $<

# how to get debug.o
#$(ODG) : %.o : %.f90
#	$(FC) $(OPT) $(OMP) -cpp -I../lib/common/include -I../lib/rtt/include -c -o $@ $<

clean:
	rm -rf *.mod
	rm -rf *.o

.PHONY : clean
