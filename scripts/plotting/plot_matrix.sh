#!/bin/bash

gnuplot -persist <<-EOFMarker

set title noenhanced '$1'
set xlabel 'sites along x'
set ylabel 'sites along y'
set cblabel 'magnitude'
set size ratio -1

plot '$1' matrix u 2:1:3 with image notitle

EOFMarker
