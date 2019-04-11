#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'U/t'
set ylabel 'n_d'
plot '$1' u 1:2 w lp 

EOFMarker
