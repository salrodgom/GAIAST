LINKFLAGS_FOR = -ffree-form -march=native -Wall -O3
COMP_FOR = gfortran
GNUPLOT_SCRIPT = png_ternary.gp
install:
	${COMP_FOR} ${LINKFLAGS_FOR} gaiast.f90 -o gaiast
all:
	make install
	make execute
	make clean
	make plot
execute:
	time ./gaiast < input
plot:
	gnuplot ${GNUPLOT_SCRIPT} ; display iast.png
clean:;         @rm -f *.o *.mod gaiast isotermaN.dat
