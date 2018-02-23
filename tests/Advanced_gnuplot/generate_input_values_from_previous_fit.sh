#!/bin/bash
function make_gp {
# langmuir
echo "set logscale x
iso1(x)=(a0*a1*x**a2)/(1+a1*x**a2) + (a3*a4*x**a5)/(1+a4*x**a5) 
iso2(x)=(b0*b1*x**b2)/(1+b1*x**b2) + (b3*b4*x**b5)/(1+b4*x**b5) 
a2=1
a5=1
b2=1
b5=1
fit iso1(x) 'isoterma1.dat' via a0,a1,a3,a4
fit iso2(x) 'isoterma2.dat' via b0,b1,b3,b4
#fit iso1(x) 'isoterma1.dat' via a0,a1,a3,a4,a2,a5
#fit iso2(x) 'isoterma2.dat' via b0,b1,b3,b4,b2,b5
" > fit.gp
}
# langmuir_freundlich_dualsite
cc real2binary_string_IEEE.c -o real2binary_string_IEEE -lm
make_gp
gnuplot  < fit.gp
if [ -f parameters.txt ] ; then rm parameters.txt ; touch parameters.txt ; fi
for parm in a0 a1 1.0 a3 a4 1.0 b0 b1 1.0 b3 b4 1.0 ; do
 parameter=$(grep ${parm} fit.log | grep "="  | grep "+/-" | awk '{print $3}' | tail -n1 )
 if [ -z "${parameter}" ]; then
  binary=$(./real2binary_string_IEEE ${parm})
  echo "$binary # $parm" >> parameters.txt
 else
  binary=$(./real2binary_string_IEEE ${parameter})
  echo "$binary # $parameter" >> parameters.txt
 fi
done
rm fit.log real2binary_string_IEEE
cat input.top > input
cat parameters.txt | awk '{print $1}' >> input
make all
exit 0
