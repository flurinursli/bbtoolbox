include ../../Makefile.inc

LIB_NAME = libcommon.a

F90_OBJ = m_precisions.o m_stat.o m_colors.o m_strings.o m_parser.o m_logfile.o

# default rule
common: dirs $(LIB_NAME)

$(F90_OBJ) : %.o : %.f90
	$(FC) $(OPT) -fopenmp -c -o $@ $<

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
