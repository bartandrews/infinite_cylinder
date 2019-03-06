#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'Φ / 2π'
set ylabel '<Q_L>'
plot '$1'

EOFMarker
