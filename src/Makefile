LINKFLAGS_FOR = -ffree-form -march=native
COMP_FOR = gfortran
GNUPLOT_SCRIPT = png.gp
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
