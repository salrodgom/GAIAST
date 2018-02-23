#!/bin/bash
cc real2binary_string_IEEE.c -o real2binary_string_IEEE -lm
if [ -f parameters.txt ] ; then rm parameters.txt ; touch parameters.txt ; fi
for parm in 2.18176103       9.63111973       1.99902534       2.00305057       6.50844193       1.06016827 2.00119758       4.26060200      0.609405756       2.24413228       7.98052168       1.99992418  ; do
  binary=$(./real2binary_string_IEEE ${parm})
  echo "$binary # $parm" >> parameters.txt
done
rm real2binary_string_IEEE
