#!/usr/bin/gnuplot
c=1.0
f(x)=K1*x/(1+(K1*x/(alfa*(1+k2*x))**c))**(1/c)
fit f(x) 'isotermaN.dat' via K1,alfa,k2
fit f(x) 'isotermaN.dat' via c
fit f(x) 'isotermaN.dat' via alfa,c
fit f(x) 'isotermaN.dat' via K1,k2
fit f(x) 'isotermaN.dat' via K1,k2,alfa,c
set t postscript eps color enhanced blacktext 'Helvetica,16'
set o 'ajuste.eps'
set logscale x
set xtics 1e-6,100,1e14
plot 'isotermaN.dat' w p pt 6 t 'input',f(x) w l lt 1 lc rgb 'black' t 'output'
