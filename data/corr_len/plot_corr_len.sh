#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'V_1 / t_1'
set ylabel 'ξ / a'
plot '$1'

EOFMarker
