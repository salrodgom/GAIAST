#!/bin/bash
declare -A parameter
function make_gp {
# langmuir
echo "set logscale x
iso1(x)=(a0*a1*x**a2)/(1+a1*x**a2) + (a3*a4*x**a5)/(1+a4*x**a5) # langmuir
iso2(x)=(b0*b1*x**b2)/(1+b1*x**b2) + (b3*b4*x**b5)/(1+b4*x**b5) # langmuir
a2=1
a5=1
b2=1
b5=1
fit iso1(x) 'isoterma1.dat' via a0,a1,a3,a4
fit iso2(x) 'isoterma2.dat' via b0,b1,b3,b4
fit iso1(x) 'isoterma1.dat' via a0,a1,a3,a4,a2,a5
fit iso2(x) 'isoterma2.dat' via b0,b1,b3,b4,b2,b5
" > fit.gp
}
make_gp
gnuplot  < fit.gp
if [ -f parameters.txt ] ; then rm parameters.txt ; touch parameters.txt ; fi
for parm in a0 a1 a2 a3 a4 a5 b0 b1 b2 b3 b4 b5 ; do
 parameter[${parm}]=$(grep ${parm} fit.log | grep "="  | grep "+/-" | awk '{print $3}' | tail -n1 )
done
rm fit.log 
cat input.top.norefit > input
echo "${parameter[a0]} ${parameter[a1]} ${parameter[a2]} ${parameter[a3]} ${parameter[a4]} ${parameter[a5]}" >> input
echo "${parameter[b0]} ${parameter[b1]} ${parameter[b2]} ${parameter[b3]} ${parameter[b4]} ${parameter[b5]}" >> input
make all
exit 0
