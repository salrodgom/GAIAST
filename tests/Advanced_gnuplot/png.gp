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
set xlabel 'Pressure / kPa'
set ylabel 'Loading / mol kg^-^1'
set multiplot layout 1,2
set logscale x
set yrange [0:]
set k top left
#set title 'Pure Compounds'
plot 'isoterma1.dat' w p pt 7 lc rgb 'red'   t '3-MP {/Symbol m}VTMC',\
     'isoterma2.dat' w p pt 7 lc rgb 'blue' t '23-DMB {/Symbol m}VTMC',\
     'curves.txt' u 1:2 w l lt 1 lc rgb 'red' title 'G.A. Fit 3-MP',\
     'curves.txt' u 4:5 w l lt 1 lc rgb 'blue' title 'G.A. Fit 23-DMB'
#set title 'Binary Mixture'
plot 'adsorcion.dat' u 1:2 w l  lt 1 lc rgb 'red' t '3-MP IAST',\
     'adsorcion.dat' u 1:3 w l  lt 1 lc rgb 'blue' t '23-DMB IAST',\
     '3-methylpentane.2load'   u 1:2 w p pt 7 ps 1 lt 1 lc rgb 'red' title '3-MP {/Symbol m}VTMC',\
     '23-dimethylbutane.2load' u 1:2 w p pt 7 ps 1 lt 1 lc rgb 'blue' title '23-DMB {/Symbol m}VTMC'
unset multiplot
