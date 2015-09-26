LINKFLAGS_FOR=-O3 -march=native -Wall
COMP_FOR = gfortran
all:
	make install
	make execute
	make clean
install:
	${COMP_FOR} ${LINKFLAGS_FOR} -c iast.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -o iast iast.o
execute:
	./iast < input && gnuplot png.gp && display iast.png
clean:;         @rm -f *.o *.mod iast
