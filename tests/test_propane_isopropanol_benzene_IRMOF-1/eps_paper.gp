#!/usr/bin/gnuplot -persist
set terminal postscript eps color enhanced blacktext 'Helvetica,12'
set size ratio 1
set o 'isotermas.eps'
set k bottom
#set format x "%2.0t{/Symbol \264}10^{%L}"
set  format x "10^{%L}"
set xlabel 'Fugacity / Pa'
set ylabel 'Loading / mol/uc'
set logscale x
set yrange [0:]
set size 0.5,0.5
plot 'isoterma1.dat'    w p pt 6 lc rgb 'red'   t 'benzene',\
     'isoterma2.dat'    w p pt 6 lc rgb 'green' t 'propane',\
     'isoterma3.dat'    w p pt 6 lc rgb 'blue'  t 'isopropanol',\
     'curves.txt' u 1:2 w l lt 1 lc rgb 'red'      notitle,\
     'curves.txt' u 4:5 w l lt 1 lc rgb 'green'    notitle,\
     'curves.txt' u 7:8 w l lt 1 lc rgb 'blue'     notitle
