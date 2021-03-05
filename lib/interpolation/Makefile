include ../../Makefile.inc

R32_OBJ = lib32/m_interpolation.o
R64_OBJ = lib64/m_interpolation.o

LIB64 = libinterpl.a
LIB32 = libinterplf.a

interpolation: dirs $(LIB64) $(LIB32)

$(R64_OBJ) : lib64/%.o : %.f90
	$(FC) $(OPT) -cpp -DDOUBLE_PREC -I../common/include -c -o $@ $<

$(R32_OBJ) : lib32/%.o : %.f90
	$(FC) $(OPT) -cpp -I../common/include -c -o $@ $<

$(LIB64) : $(R64_OBJ)
	ar qc $@ $(R64_OBJ)
	mv $(LIB64) lib/
	cp *r64.mod include64/

$(LIB32) : $(R32_OBJ)
	ar qc $@ $(R32_OBJ)
	mv $(LIB32) lib/
	cp *r32.mod include32/

dirs :
	mkdir -p lib lib32 lib64 include32 include64

clean :
	rm -rf lib
	rm -rf lib32
	rm -rf lib64
	rm -rf include32
	rm -rf include64
	rm -rf *.mod

.PHONY: clean
