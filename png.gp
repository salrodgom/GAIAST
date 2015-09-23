#!/usr/bin/gnuplot -persist
#set t postscript eps color enhanced blacktext 'Helvetica,14'
set terminal pngcairo size 820,308 enhanced font 'Verdana,10'
set size ratio 1
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 11 lc rgb '#808080' lt 1
set o 'iast.png'
set k top left
#set xtics 1e-6,100,1e14
#set format x "%2.0t{/Symbol \264}10^{%L}"
set xlabel 'Pressure / kPa'
set ylabel 'Loading / mol kg^-^1'
set format x "10^{%L}"
#set size 1,0.5
set multiplot layout 1,2
#set size 0.5,0.5
set logscale x
plot './isoterma1.dat' pt 6 t 'C_1','./isoterma2.dat' w p pt 6 t 'C_2','iso1.dat' w l lt 1 lc rgb 'red' notitle,'iso2.dat' w l lt 1 lc rgb 'green' notitle
plot './adsorcion.dat' u 1:2 w l lt 1 lc rgb 'red' t 'C_1',\
     './adsorcion.dat' u 1:3 w l lt 1 lc rgb 'green' t 'C_2',\
     './adsorcion.dat' u 1:4 w l lt 1 lc rgb 'blue' t 'C_1+C_2'
unset multiplot
