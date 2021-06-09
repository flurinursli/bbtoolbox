include ../../Makefile.inc

PP_FLAGS =

ifeq ($(MKL),y)
	PP_FLAGS += -DMKL
endif

R32_OBJ = lib32/m_llsq.o
R64_OBJ = lib64/m_llsq.o

LIB64 = libllsq.a
LIB32 = libllsqf.a

llsq: dirs $(LIB64) $(LIB32)

$(R64_OBJ) : lib64/%.o : %.f90
	$(FC) $(OPT) $(OMP) -cpp $(PP_FLAGS) -DDOUBLE_PREC -I../common/include -c -o $@ $<

$(R32_OBJ) : lib32/%.o : %.f90
	$(FC) $(OPT) $(OMP) -cpp $(PP_FLAGS) -I../common/include -c -o $@ $<

$(LIB64) : $(R64_OBJ)
	ar qc $@ $(R64_OBJ)
	mv $(LIB64) lib/
	cp *.mod include64/

$(LIB32) : $(R32_OBJ)
	ar qc $@ $(R32_OBJ)
	mv $(LIB32) lib/
	cp *.mod include32/

dirs :
	mkdir -p lib lib32 lib64 include32 include64

clean :
	rm -rf lib
	rm -rf lib32
	rm -rf lib64
	rm -rf include32
	rm -rf include64
	rm -rf *.mod

.PHONY: clean dirs llsq
