#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'independent variable'
set ylabel 'entanglement energy'
plot '$1' u 2:3

EOFMarker
