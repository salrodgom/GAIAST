#!/usr/bin/gnuplot -persist
set terminal pngcairo size 820,308 enhanced font 'Verdana,10'
set size ratio 1
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 11 lc rgb '#808080' lt 1
set o 'iast.png'
set k top left
set title "A, B, and C pure adsorption isotherms in X nanoporous material"
set format x "10^{%L}"
set xlabel "Fugacity / Pa"
set ylabel "Loading / mol/kg"
set multiplot layout 1,2
set logscale x
set xrange [1e-1:1e4]
! sort -gk1 adsorcion.dat > c ; mv c adsorcion.dat
plot 'isoterma1.dat' w p pt 7 lc rgb "red"   t "A",\
     'isoterma2.dat' w p pt 7 lc rgb "dark-green" t "B",\
     'isoterma2.dat' w p pt 7 lc rgb "orange" t "C",\
     'curves.txt' u 1:2 w l lt 1 lc rgb "red"   notitle,\
     'curves.txt' u 4:5 w l lt 1 lc rgb "dark-green" notitle,\
     'curves.txt' u 7:8 w l lt 1 lc rgb "orange" notitle
unset ytics
set autoscale
set key outside
set ytics nomirror
set y2tics nomirror
set y2label "Selectivity / -" tc rgb "blue"
set title "A/B/C ternay mixture in X nanoporous material"
plot 'adsorcion.dat' u 1:2 w l lt 1 lc rgb "red" t "A",\
     'adsorcion.dat' u 1:3 w l lt 1 lc rgb "dark-green" t "B",\
     'adsorcion.dat' u 1:4 w l lt 1 lc rgb "orange" t "C",\
     'adsorcion.dat' u 1:5 w l lt 1 lc rgb "black" t "A+B+C",\
     'adsorcion.dat' u 1:($3/$2) axes x1y2 w l lt 1 lc rgb "blue" t "B/A",1.0 axes x1y2 w l lt 2 lc rgb "gray" notitle,\
     'adsorcion.dat' u 1:($3/$4) axes x1y2 w l lt 1 lc rgb "dark-blue" t "B/C"
unset multiplot
