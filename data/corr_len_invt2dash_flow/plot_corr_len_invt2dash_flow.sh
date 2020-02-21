#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel '1 / t_2''
set ylabel 'ξ / a'
plot '$1'

EOFMarker
