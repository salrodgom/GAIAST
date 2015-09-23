LINKFLAGS_FOR=-O3 -march=native
COMP_FOR = gfortran
all:
	${COMP_FOR} ${LINKFLAGS_FOR} -c iast.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -o iast iast.o
clean:;         @rm -f *.o *.mod iast
