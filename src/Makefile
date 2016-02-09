LINKFLAGS_FOR = -pedantic -O2 -march=native
COMP_FOR = gfortran
GCCVERSION = $(shell gfortran -dumpversion)
GCCVERSIONGTEQ4 := $(shell expr `${COMP_FOR} -dumpversion | cut -f1 -d.` \>= 4)
GNUPLOT_SCRIPT = png.gp
ifeq "$(GCCVERSIONGTEQ4)" "1"
    LINKFLAGS_FOR += -Wtype-limits
endif
install:
	${COMP_FOR} ${LINKFLAGS_FOR} gaiast.f90 -o gaiast
all:
	make install
	make execute
	make clean
	make plot
execute:
	./gaiast < input
plot:
	gnuplot ${GNUPLOT_SCRIPT} ; display iast.png
clean:;         @rm -f *.o *.mod gaiast isotermaN.dat
