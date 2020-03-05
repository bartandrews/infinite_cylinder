#!/bin/bash

# $1 = the block of data to plot
# $2 = the data file

if (( $1 == 1 ))
then
    awk '/LylB=9.330721809337284/,/^$/' $2 > $2.temp
elif (( $1 == 2 ))
then
    awk '/LylB=13.996082714005926/,/^$/' $2 > $2.temp
elif (( $1 == 3 ))
then
    awk '/LylB=8.080642122531591/,/^$/' $2 > $2.temp 
elif (( $1 == 4 )) 
then
    awk '/LylB=12.120963183797386/,/^$/' $2 > $2.temp
elif (( $1 == 5 ))
then
    awk '/LylB=7.227546035131528/,/^$/' $2 > $2.temp
elif (( $1 == 6 ))
then
    awk '/LylB=10.841319052697292/,/^$/' $2 > $2.temp
fi

# remove first and last line
sed -i '1d;$d' $2.temp

# sort by first column
sort -n $2.temp > $2.temp2

# insert two blank lines between each block
awk -v i=1 'NR>1 && $i!=p { print "\n" }{ p=$i } 1' $2.temp2 > $2.temp3

# define the minimum charge number
min_charge=$(awk 'NR == 1 {print $1}' $2.temp3)

gnuplot -persist <<-EOFMarker

stats '$2.temp3' u 2:3 nooutput
blocks = STATS_blocks

set key out right box

set title noenhanced '$2'
set xlabel 'momentum / π'
set ylabel 'entanglement energy'

plot for[i=0:blocks-1] '$2.temp3' index i u 2:3 pt i+1 title sprintf("%i", i+$min_charge)

EOFMarker

# remove temp files
rm $2.temp $2.temp2 $2.temp3
