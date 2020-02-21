#!/bin/bash

# sort by first column
sort -n $1 > $1.temp

# insert two blank lines between each block
awk -v i=1 'NR>1 && $i!=p { print "\n" }{ p=$i } 1' $1.temp > $1.temp2

# define the minimum charge number
min_charge=$(awk 'NR == 1 {print $1}' $1.temp2)

gnuplot -persist <<-EOFMarker

stats '$1.temp2' u 2:3 nooutput
blocks = STATS_blocks

set key out right box

set title noenhanced '$1'
set xlabel 'τ'
set ylabel 'entanglement energy'

#set yrange [:6]

plot for[i=0:blocks-1] '$1.temp2' index i u 2:3 pt i+1 title sprintf("%i", i+$min_charge)

EOFMarker

# remove temp files
rm $1.temp $1.temp2
