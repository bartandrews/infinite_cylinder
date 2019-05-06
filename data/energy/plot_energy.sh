#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'Js/J'
set ylabel 'E'
plot '$1' u 1:2 w lp 

EOFMarker
