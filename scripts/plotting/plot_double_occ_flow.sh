#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'independent variable'
set ylabel 'double occupancy'
plot '$1'

EOFMarker
