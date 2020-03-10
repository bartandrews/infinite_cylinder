#!/bin/bash

gnuplot -persist <<-EOFMarker

set title noenhanced '$1'
set xlabel 'sites along x'
set ylabel 'sites along y'
set cblabel 'correlation function'

plot '$1' matrix with image

EOFMarker
