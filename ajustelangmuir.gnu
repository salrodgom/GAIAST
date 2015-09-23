#!/usr/bin/gnuplot
alfa=1.0
f(x)=Nmax*alfa*x/(1+alfa*x)
fit f(x) 'isotermaN.dat' via Nmax
fit f(x) 'isotermaN.dat' via alfa
fit f(x) 'isotermaN.dat' via Nmax,alfa
set t postscript eps color enhanced blacktext 'Helvetica,16'
set o 'ajuste.eps'
set logscale x
set xtics 1e-6,100,1e14
plot 'isotermaN.dat' w p pt 6 t 'input',f(x) w l lt 1 lc rgb 'black' t 'output'
