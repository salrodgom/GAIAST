#!/usr/bin/gnuplot -persist
set terminal pngcairo size 308,308 enhanced font 'Helvetica,10'
set size ratio 1
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 11 lc rgb '#808080' lt 1
set o 'iast.png'
set k top right
set xlabel 'Temperature / K'
set ylabel 'Loading / a. u.'
set yrange [0:]
plot 'isoterma1.dat' w p pt 6 lc rgb 'red'   t 'C_1',\
     'curves.txt' u 1:2 w l lt 1 lc rgb 'red'   notitle
