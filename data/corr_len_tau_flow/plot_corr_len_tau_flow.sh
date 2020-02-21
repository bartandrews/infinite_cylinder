#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'τ'
set ylabel 'ξ / a'
plot '$1'

EOFMarker
