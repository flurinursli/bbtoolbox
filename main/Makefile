include ../Makefile.inc

PP_FLAGS =

ifeq ($(MODE),debug)
	PP_FLAGS += -DPERF -DERROR_TRAP -DDEBUG
else
ifeq ($(ERROR_TRAP),y)
	PP_FLAGS += -DERROR_TRAP
endif
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

m_toolbox.o : $(CPP_OBJ)

# list C++ objects
CPP_OBJ = geo2utm.o

# compile C++ objects
$(CPP_OBJ) : %.o : %.cpp
	$(C++) -O3 -c -o $@ $<

OBJ = m_toolbox.o m_userdef.o m_timeseries.o m_roughness.o m_source.o m_noise.o m_rik.o main.o

main : $(OBJ) $(CPP_OBJ)
	$(FC) $(OPT) $(OMP) -o ../toolbox.exe $(CPP_OBJ) $(OBJ) -L../lib/common/lib -L../lib/llsq/lib -L../lib/fft/lib -L../lib/interpolation/lib -L../lib/filter/lib -L../lib/rtt/lib -L../lib/rng/lib -L$(PROJ4_PATH) -L$(TRNG4_PATH)/lib -L$(FFTW_PATH)/lib -L$(LAPACK_PATH) -L$(GSL_PATH) -lcommon -lllsqf -lfilterf -lfft -linterplf -lllsq -lrtt -lrngf $(LINK_FLAGS) -lgsl -lgslcblas -linterpl -ltrng4 -lstdc++ -lproj

# dbg	: $(OBJ_DBG) $(CPP_OBJ)
# 	$(FC) $(OPT) $(OMP) -o ../toolbox.exe $(CPP_OBJ) $(OBJ) -L../lib/common/lib -L../lib/llsq/lib -L../lib/fft/lib -L../lib/interpolation/lib -L../lib/filter/lib -L../lib/rtt/lib -L../lib/rng/lib -L$(PROJ4_PATH) -L$(TRNG4_PATH)/lib -L$(FFTW_PATH)/lib -L$(LAPACK_PATH) -L$(GSL_PATH) -lcommon -lllsqf -lfilterf -lfft -linterplf -lllsq -lrtt -lrngf $(LINK_FLAGS) -lgsl -lgslcblas -linterpl -ltrng4 -lstdc++ -lproj

# how to get main.o
$(OBJ) : %.o : %.f90
	$(FC) $(OPT) $(OMP) -cpp $(PP_FLAGS) -I../lib/common/include -I../lib/llsq/include32 -I../lib/interpolation/include32 -I../lib/fft/include -I../lib/filter/include -I../lib/llsq/include32 -I../lib/rtt/include -I../lib/rng/include -c -o $@ $<

clean:
	rm -rf *.mod
	rm -rf *.o

.PHONY : clean
