#!/bin/bash

gnuplot -persist <<-EOFMarker

set key box top left
set title noenhanced '$1'

set xlabel 'independent variable'
set ylabel 'dependent variable'

set xrange [0:]

set fit errorvariables
f(x) = a*x + b
fit f(x) '$1' u 1:2:3 yerrors via a, b
title_f(a,b) = sprintf('y = (%.2f±%.2f) x + (%.2f±%.2f)', a, a_err, b, b_err)

plot '$1' u 1:2:3 lc 1 w yerrorbars notitle, f(x) lc 1 title title_f(a,b)

EOFMarker
