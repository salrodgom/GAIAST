#!/usr/bin/gnuplot 
a=1.0
b=0.1
f(x)=a*x**b
fit f(x) 'isotermaN.dat' via a
fit f(x) 'isotermaN.dat' via b
fit f(x) 'isotermaN.dat' via a,b
set t postscript eps color enhanced blacktext 'Helvetica,16'
set o 'ajuste.eps'
set logscale x
set xtics 1e-6,100,1e14
plot 'isotermaN.dat' w p pt 6 t 'input',f(x) w l lt 1 lc rgb 'black' t 'output'
