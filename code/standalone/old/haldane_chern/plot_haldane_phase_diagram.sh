#!/bin/bash

gnuplot -persist <<-EOFMarker

	set title '$1 (for the 0 band)' noenhanced
	set xlabel 'phi'
	set ylabel 'M/t_2'

	set palette maxcolors 100
	set palette defined (0 "blue", 50 "white", 99 "red")

	set xtics ('-π' -pi, 0, 'π' pi)
	set ytics ('-3√3' -3*sqrt(3), 0, '3√3' 3*sqrt(3))

	set pm3d map
	unset key

	set cblabel 'C'

	plot '$1' u 1:2:3 with image

EOFMarker
