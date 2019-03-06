#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'bond'
set ylabel 'entanglement energy'
plot '$1'

EOFMarker
