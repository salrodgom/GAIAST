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
set xlabel 'Fugacity / Pa'
set ylabel 'Loading / mol kg^-^1'
set multiplot layout 1,2
set logscale x
set yrange [0:18]
set xrange [1e2:1e9]
set title 'Pure Components'
plot 'isoterma1.dat' w p pt 6 lc rgb 'red'   t 'C_2H_6 GCMC Jurn',\
     'isoterma2.dat' w p pt 6 lc rgb 'blue'  t 'C_2H_4 GCMC Jurn',\
     'curves.txt' u 1:2 w l lt 1 lc rgb 'red'  notitle,\
     'curves.txt' u 4:5 w l lt 1 lc rgb 'blue' notitle,\
     'ethene-323K-Jurn-EXP.txt' w p pt 7 lc rgb 'black' t 'C_2H_4 EXP 323K'
set title 'Binary Mixture'
set yrange [0:18]
plot 'adsorcion.dat' u 1:2 w p pt 7 ps 0.5 lt 1 lc rgb 'red' t 'C_2H_6 Jurn IAST',\
     'adsorcion.dat' u 1:3 w p pt 7 ps 0.5 lt 1 lc rgb 'blue' t 'C_2H_4 Jurn IAST',\
     'jurn-ethane-mix.dat' u 1:2 w p pt 6 ps 1 lt 1 lc rgb 'red' t 'C_2H_6 GCMC Jurn',\
     'jurn-ethene-mix.dat' u 1:2 w p pt 6 ps 1 lt 1 lc rgb 'blue' t 'C_2H_4 GCMC Jurn'
unset multiplot
