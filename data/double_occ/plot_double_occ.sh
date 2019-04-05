#!/bin/bash

gnuplot -persist <<-EOFMarker

unset key
set title noenhanced '$1'

set xlabel 'U/t'
set ylabel 'n_d U'
plot '$1'

EOFMarker
