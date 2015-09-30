LINKFLAGS_FOR=-O3 -march=native
COMP_FOR = gfortran
all:
	make install
	make execute
	make clean
install:
	${COMP_FOR} ${LINKFLAGS_FOR} iast.f90 -o iast
execute:
	./iast < input && gnuplot png.gp && display iast.png
clean:;         @rm -f *.o *.mod iast isotermaN.dat
