include ../../Makefile.inc

PRECISION = single

ifeq ($(PRECISION),double)
	LIB_NAME = libcompgeo.a
else
	LIB_NAME = libcompgeof.a
endif

PP_FLAGS =

ifeq ($(PRECISION),double)
	PP_FLAGS += -DDOUBLE_PREC
endif

F90_OBJS = m_compgeo.o

geometry: dirs $(LIB_NAME)

$(F90_OBJS) : %.o : %.f90
	$(FC) $(OPT) -cpp $(PP_FLAGS) -I../common/include -c -o $@ $<

# list C++ objects
CPP_OBJS = opttritri.o

# compile C++ objects
$(CPP_OBJS) : %.o : %.cpp
	$(CXX) -O2 -cpp $(PP_FLAGS) $(INCL) -c -o $@ $<

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
