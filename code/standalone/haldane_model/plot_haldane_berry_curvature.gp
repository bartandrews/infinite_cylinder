set pm3d map
set size square
unset key

set palette maxcolors 100
set palette defined (0 "blue", 50 "white", 99 "red")

set xtics ('-π' -pi, 0, 'π' pi)
set ytics ('-π' -pi, 0, 'π' pi)

set cblabel 'B'

set title 'Haldane Model Berry Curvature (for the 0 band)'
set xlabel 'kx'
set ylabel 'ky'

splot 'haldane_berry_curvature.txt' u 1:2:3
