#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'Φ / 2π'
set ylabel 'entanglement energy'
plot '$1' u 2:3

EOFMarker
