include ../../Makefile.inc

INCL = -I$(TRNG4_PATH)/include

PRECISION = single

ifeq ($(PRECISION),double)
	LIB_NAME = librng.a
else
	LIB_NAME = librngf.a
endif

PP_FLAGS =

ifeq ($(PRECISION),double)
	PP_FLAGS += -DDOUBLE_PREC
endif

# list Fortran objects
F90_OBJS = m_random.o

# default rule
rng : dirs $(LIB_NAME)

# compile Fortran objects
$(F90_OBJS) : %.o : %.f90
	$(FC) $(OPT) $(OMP) -cpp $(PP_FLAGS) -I../common/include -c -o $@ $<

# list C++ objects
CPP_OBJS = trng.o

# compile C++ objects
$(CPP_OBJ) : %.o : %.cpp
	$(C++) -O3 -c -o $@ $<

# compile C++ objects
$(CPP_OBJS) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) -cpp $(PP_FLAGS) $(INCL) -c -o $@ $<

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
