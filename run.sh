#!/bin/bash

rm sod.dat
gcc -o sod_exact sod_exact.c -lm
./sod_exact

(
    echo 'set xrange [0:1]'
    echo 'set yrange [0:1]'
    echo 'set xlabel "position (x)"'
    echo 'set ylabel "density"'
    echo 'plot "sod.dat" w l'
    echo 'pause -1'
) | gnuplot -persist