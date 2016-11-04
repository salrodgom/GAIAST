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
set xlabel 'Fugacity / kPa'
set ylabel 'Loading / mol kg^-^1'
set multiplot layout 1,2
set logscale x
set yrange [0:]
set xrange [1e-1:]
plot 'isoterma1.dat' w p pt 6 lc rgb 'red'   t 'C_2H_6 Own GCMC',\
     'isoterma2.dat' w p pt 6 lc rgb 'blue' t 'C_2H_4 Own GCMC',\
     'curves.txt' u 1:2 w l lt 1 lc rgb 'red'   notitle,\
     'curves.txt' u 4:5 w l lt 1 lc rgb 'blue' notitle,\
     'ethene-303-jurn-dubbeldam-EXP.txt' u ($1/1000):2 w p pt 7 ps 1 lc rgb 'blue' t 'C_2H_4 Jurn EXP 303 K',\
     'exp-ethene-Cu-BTC-296K-2.dat'                    w p pt 9 ps 1 lc rgb 'blue' t 'C_2H_4 Jurn EXP 296 K'
plot 'adsorcion.dat' u 1:2 w p pt 7 ps 0.5 lt 1 lc rgb 'red' t 'ethane IAST',\
     'adsorcion.dat' u 1:3 w p pt 7 ps 0.5 lt 1 lc rgb 'blue' t 'ethene IAST',\
     'ethane-mix-simu-298K.dat' w p pt 6 ps 1 lt 1 lc rgb 'red'  t 'ethane Own GCMC',\
     'ethene-mix-simu-298K.dat' w p pt 6 ps 1 lt 1 lc rgb 'blue' t 'ethane Own GCMC'
#
#     'ghff-ethane-mix.dat' u ($1/1000):2 w p pt 6 ps 1 lt 1 lc rgb 'red' t 'ethane GCMC GHFF',\
#     'ghff-ethane-mix.dat' u ($1/1000):2 w p pt 6 ps 1 lt 1 lc rgb 'blue' t 'ethene GCMC GHFF'
unset multiplot
