set pm3d map
set size square
unset key

set xtics ('-π' -pi, 0, 'π' pi)
set ytics ('-π' -pi, 0, 'π' pi)

set title 'Berry Curvature (for a+ band)'
set xlabel 'kx'
set ylabel 'ky'

splot 'berry_curvature.txt' u 1:2:3
