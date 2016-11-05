LINKFLAGS_FOR = -pedantic -O3 -march=native
COMP_FOR = gfortran
GNUPLOT_SCRIPT = png.gp
check:
	GFORTRANVERSION = $(shell gfortran -dumpversion)
	GFORTRANVERSIONGTEQ4 := $(shell expr `gfortran -dumpversion` \== 4.9.2)
	ifeq "$(GFORTRANVERSIONGTEQ4)" "1"
		$(error Version de gfortran no valida: use 4.9.2 o superior )
	endif
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
