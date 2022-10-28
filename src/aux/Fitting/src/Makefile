LINKFLAGS_FOR = -std=gnu -pedantic -ffree-form -march=native -O3 -Wall -frecursive
COMP_FOR = gfortran
all:
	make install
	make execute
	make clean
install:
	${COMP_FOR} ${LINKFLAGS_FOR} random.f90 -c
	${COMP_FOR} ${LINKFLAGS_FOR} sorting.f90 -c
	${COMP_FOR} ${LINKFLAGS_FOR} globals.f90 -c
	${COMP_FOR} ${LINKFLAGS_FOR} ga.f90 -c
	${COMP_FOR} ${LINKFLAGS_FOR} simplex.f90 -c
	${COMP_FOR} ${LINKFLAGS_FOR} -o fitting.exe main.f90 *.o
execute:
	./fitting < input
plot:
	gnuplot ${GNUPLOT_SCRIPT} ; display iast.png
clean:;         @rm -f *.o *.mod fitting
