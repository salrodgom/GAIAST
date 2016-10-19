#!/usr/bin/gnuplot -persist
#set terminal pngcairo size 300,300 enhanced font 'Helvetica,10'
set size ratio 1
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 11 lc rgb '#808080' lt 1
#set o 'isotermas_saturation.png'
set k bottom
#set format x "%2.0t{/Symbol \264}10^{%L}"
set format x "10^{%L}"
set xlabel 'Fugacity / Pa'
set ylabel '% saturation / -'
set logscale x
set yrange [0:101]
plot 'isoterma1.dat' u 1:(100*$2/14) w p pt 6 lc rgb 'red'   t 'benzene',\
     'isoterma2.dat' u 1:(100*$2/17) w p pt 6 lc rgb 'green' t 'propane',\
     'isoterma3.dat' u 1:(100*$2/16) w p pt 6 lc rgb 'blue' t 'isopropanol',\
     'curves.txt' u 1:(100*$2/14) w l lt 1 lc rgb 'red'   notitle,\
     'curves.txt' u 4:(100*$5/17) w l lt 1 lc rgb 'green' notitle,\
     'curves.txt' u 7:(100*$8/16) w l lt 1 lc rgb 'blue' notitle,20 w l lt -1 lc rgb 'gray' notitle,40 w l lt -1 lc rgb 'gray' notitle,60 w l lt -1 lc rgb 'gray' notitle,80 w l lt -1 lc rgb 'gray' notitle
