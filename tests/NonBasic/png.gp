#!/usr/bin/gnuplot -persist
set terminal pngcairo size 820,308 enhanced font 'Helvetica,10'
set size ratio 1
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 11 lc rgb '#808080' lt 1
set o 'iast.png'
set k top left
set format x "10^{%L}"
set xlabel "Fugacity [Pa]"
set ylabel "Loading [mol/kg]"
set multiplot layout 1,2
set logscale x
! sort -gk1 adsorcion.dat > c ; mv c adsorcion.dat
set xrange [1e-9:1e4]
plot 'isoterma1.dat' w p pt 7 lc rgb "red"   t "p-xylene",\
     'isoterma2.dat' w p pt 7 lc rgb "dark-green" t "o-xylene",\
     'isoterma3.dat' w p pt 7 lc rgb "blue" t "m-xylene",\
     'isoterma4.dat' w p pt 7 lc rgb "orange" t "ethylbenzene",\
     'curves.txt' u 1:2 w l lt 1 lc rgb "red"   notitle,\
     'curves.txt' u 4:5 w l lt 1 lc rgb "dark-green" notitle,\
     'curves.txt' u 7:8 w l lt 1 lc rgb "blue" notitle,\
     'curves.txt' u 10:11 w l lt 1 lc rgb "orange" notitle
set ytics nomirror
#set y2tics nomirror
#set y2label 'Selectivity / -' tc rgb "blue"
plot 'adsorcion.dat' u 1:2 w lp ps 0.5 pt 7 lc rgb "red" t "p-xylene",\
     'adsorcion.dat' u 1:3 w lp ps 0.5 pt 7 lc rgb "dark-green" t "o-xylene",\
     'adsorcion.dat' u 1:4 w lp ps 0.5 pt 7 lc rgb "blue" t "m-xylene",\
     'adsorcion.dat' u 1:5 w lp ps 0.5 pt 7 lc rgb "orange" t "ethylbenzene"
unset multiplot
