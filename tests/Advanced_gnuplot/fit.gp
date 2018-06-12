set logscale x
iso1(x)=(a0*a1*x**a2)/(1+a1*x**a2) + (a3*a4*x**a5)/(1+a4*x**a5) 
iso2(x)=(b0*b1*x**b2)/(1+b1*x**b2) + (b3*b4*x**b5)/(1+b4*x**b5) 
a2=1
a5=1
b2=1
b5=1
fit iso1(x) 'isoterma1.dat' via a0,a1,a3,a4
fit iso2(x) 'isoterma2.dat' via b0,b1,b3,b4
fit iso1(x) 'isoterma1.dat' via a0,a1,a3,a4,a2,a5
fit iso2(x) 'isoterma2.dat' via b0,b1,b3,b4,b2,b5

