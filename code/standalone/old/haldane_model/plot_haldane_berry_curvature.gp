set title 'Haldane Model Berry Curvature (for the 0 band)'
set xlabel 'k_x'
set ylabel 'k_y'

set pm3d map
set size square
unset key

set palette maxcolors 100
set palette defined (0 "blue", 50 "white", 99 "red")

set xrange [-pi:pi]
set yrange [-pi:pi]

set xtics ('-2π/3' -2*pi/3, 0, '2π/3' 2*pi/3)
set ytics ('-2π/(3√3)' -2*pi/(3*sqrt(3)), 0, '2π/(3√3)' 2*pi/(3*sqrt(3)))

set cblabel 'B'

K1x = 2*pi/3
K1y = 2*pi/(3*sqrt(3))

set arrow 1 from K1x,K1y to 0,2*K1y nohead front
set arrow 2 from 0,2*K1y to -K1x,K1y nohead front
set arrow 3 from -K1x,K1y to -K1x,-K1y,0 nohead front
set arrow 4 from -K1x,-K1y to 0,-2*K1y nohead front
set arrow 5 from 0,-2*K1y to K1x,-K1y nohead front
set arrow 6 from K1x,0 to K1x,K1y nohead front

set arrow 7 from K1x,K1y to 0,0 nohead front dt 2 lw 2
set arrow 8 from 0,0 to K1x,0 nohead front dt 2 lw 2
set arrow 9 from K1x,0 to K1x,-K1y nohead front dt 2 lw 2

set label "Γ" at 0,0 front
set label "M" at K1x,0 front
set label "K" at K1x,K1y front
set label "K'" at K1x,-K1y front

splot 'haldane_berry_curvature.txt' u 1:2:3
