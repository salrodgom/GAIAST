#!/usr/bin/gnuplot 
c=1.0
f(x)=Nmax*alfa*x/(1+(alfa*x)**c)**(1/c)
fit f(x) 'isotermaN.dat' via  Nmax,alfa
fit f(x) 'isotermaN.dat' via  c
fit f(x) 'isotermaN.dat' via  Nmax,alfa
fit f(x) 'isotermaN.dat' via  c,alfa
fit f(x) 'isotermaN.dat' via  Nmax,alfa,c
