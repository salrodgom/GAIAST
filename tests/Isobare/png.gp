#!/usr/bin/gnuplot -persist
set terminal pngcairo size 820,308 enhanced font 'Verdana,10'
set size ratio 1
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 11 lc rgb '#808080' lt 1
set o 'iast.png'
set k top left
#set format x "%2.0t{/Symbol \264}10^{%L}"
set format x "10^{%L}"
set xlabel 'Pressure / a. u.'
set ylabel 'Loading / a. u.'
set multiplot layout 1,2
set logscale x
set yrange [0:]
plot 'isoterma1.dat' w p pt 6 lc rgb 'red'   t 'C_1',\
     'isoterma2.dat' w p pt 6 lc rgb 'green' t 'C_2',\
     'curves.txt' u 1:2 w l lt 1 lc rgb 'red'   notitle,\
     'curves.txt' u 4:5 w l lt 1 lc rgb 'green' notitle
plot 'adsorcion.dat' u 1:2 w p pt 7 ps 0.5 lt 1 lc rgb 'red' t 'C_1',\
     'adsorcion.dat' u 1:3 w p pt 7 ps 0.5 lt 1 lc rgb 'green' t 'C_2'
unset multiplot
