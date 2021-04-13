include ../../Makefile.inc

PRECISION = single

ifeq ($(PRECISION),double)
	LIB_NAME = libeikonal.a
else
	LIB_NAME = libeikonalf.a
endif

PP_FLAGS =

ifeq ($(PRECISION),double)
	PP_FLAGS += -DDOUBLE_PREC
endif

F90_OBJS = m_eikonal.o

eikonal: dirs $(LIB_NAME)

$(F90_OBJS) : %.o : %.f90
	$(FC) $(OPT) -cpp $(PP_FLAGS) -I../common/include -c -o $@ $<

# list C++ objects
CPP_OBJS = fmm.o

# compile C++ objects
$(CPP_OBJS) : %.o : %.cpp
	$(CXX) -O2 -cpp -DNDEBUG $(PP_FLAGS) -c -o $@ $<

# define list of Fortran + C objects
OBJS = $(F90_OBJS)
OBJS += $(CPP_OBJS)

$(LIB_NAME): $(OBJS)
	ar qc $@ $(OBJS)
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
