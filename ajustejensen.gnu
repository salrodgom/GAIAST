#!/usr/bin/gnuplot
set t postscript eps color enhanced
set o 'cosa.eps'
c=1.0
f(x)=K1*x/(1+(K1*x/(alfa*(1+k2*x))**c))**(1/c)
fit f(x) 'isotermaN.dat' via K1,alfa,k2
fit f(x) 'isotermaN.dat' via c
fit f(x) 'isotermaN.dat' via alfa,c
fit f(x) 'isotermaN.dat' via K1,k2
set logscale x
plot 'isotermaN.dat',f(x)
