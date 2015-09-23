#!/usr/bin/gnuplot 
f(x)=Nmax*alfa*x/(1+alfa*x)
fit f(x) 'isotermaN.dat' via Nmax,alfa
