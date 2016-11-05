#!/usr/bin/gnuplot -persist
set term postscript eps color enhanced blacktext 'Helvetica,12'
set out 'iast.eps'
set size ratio 1
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 11 lc rgb '#808080' lt 1
set k top left
#set format x "%2.0t{/Symbol \264}10^{%L}"
set format x "10^{%L}"
set xlabel 'Bulk fluid phase fugacity / kPa'
set ylabel 'Loading / mol kg^-^1'
set multiplot layout 1,2
set logscale x
set yrange [0:]
set xrange [1e-1:]
plot 'isoterma1.dat' w p pt 6 lc rgb 'red'   t 'C_2H_6 Own GCMC 298K',\
     'isoterma2.dat' w p pt 6 lc rgb 'blue' t 'C_2H_4 Own GCMC 298K',\
     'test_jurn/isoterma1.dat' u ($1/1000):2 w p pt 6 lc rgb '#FF1493'   t 'C_2H_6 Jurn GCMC 298K',\
     'test_jurn/isoterma2.dat' u ($1/1000):2 w p pt 6 lc rgb '#1E90FF'   t 'C_2H_4 Jurn GCMC 298K',\
     'curves.txt' u 1:2 w l lt 1 lc rgb 'red'   notitle,\
     'curves.txt' u 4:5 w l lt 1 lc rgb 'blue' notitle,\
     'ethene-303-jurn-dubbeldam-EXP.txt' u ($1/1000):2 w p pt 7 ps 1 lc rgb 'blue' t 'C_2H_4 Jurn EXP 303K',\
     'exp-ethane-Cu-BTC-296K-2.dat'                    w p pt 7 ps 1 lc rgb 'red'  t 'C_2H_6 Khrisna EXP 303K',\
     'exp-ethene-Cu-BTC-296K-2.dat'                    w p pt 9 ps 1 lc rgb 'blue' t 'C_2H_4 Khrisna EXP 296K'
plot 'adsorcion.dat' u 1:2 w p pt 7 ps 0.5 lt 1 lc rgb 'red' t 'C_2H_6 IAST',\
     'adsorcion.dat' u 1:3 w p pt 7 ps 0.5 lt 1 lc rgb 'blue' t 'C_2H_4 IAST',\
     'ethane-mix-simu-298K.dat' w p pt 6 ps 1 lt 1 lc rgb 'red'  t 'C_2H_6 Own GCMC',\
     'ethene-mix-simu-298K.dat' w p pt 6 ps 1 lt 1 lc rgb 'blue' t 'C_2H_4 Own GCMC'
unset multiplot
