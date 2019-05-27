set title 'Haldane Model Phase Diagram (for the 0 band)'
set xlabel 'phi'
set ylabel 'M/t2'

set palette maxcolors 100
set palette defined (0 "blue", 50 "white", 99 "red")

set xtics ('-π' -pi, 0, 'π' pi)
set ytics ('-3√3' -3*sqrt(3), 0, '3√3' 3*sqrt(3))

set pm3d map
unset key

set cblabel 'C'

splot 'haldane_phase_diagram.txt' u 1:2:3
