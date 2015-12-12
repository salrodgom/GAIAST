LINKFLAGS_FOR=-O2 -march=native -fbounds-check
COMP_FOR = gfortran
install:
	${COMP_FOR} ${LINKFLAGS_FOR} gaiast.f90 -o gaiast
all:
	make install
	make execute
	make clean
execute:
	time ./gaiast < input
clean:;         @rm -f *.o *.mod gaiast isotermaN.dat
