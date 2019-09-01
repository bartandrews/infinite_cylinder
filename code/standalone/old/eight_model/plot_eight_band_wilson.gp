### Start multiplot (1x2 layout)
set multiplot layout 1,2 rowsfirst
# --- GRAPH a
set title 'A/B character'
set ylabel 'Energy / meV'
set cblabel '|C_A|^2'

set size ratio 1
set ytics 2
#set yrange [-5:5]

set key out below

num_points = 3*30

set xrange [0: num_points]
set xtics ('K' 0, 'Γ' num_points/3, 'M' 2*num_points/3, "K\'" num_points)

set arrow from num_points/3,graph(0,0) to num_points/3,graph(1,1) nohead dt 2
set arrow from 2*num_points/3,graph(0,0) to 2*num_points/3,graph(1,1) nohead dt 2

#set label "C=-1" at 70,10.5
#set label "C=1" at 70,7.5
#set label "C=-1" at 70,4.5
#set label "C=1" at 70,1.5
#set label "C=-1" at 70,-1.5
#set label "C=1" at 70,-4.5
#set label "C=-1" at 70,-7.5
#set label "C=1" at 70,-11

plot '2D_eight_band_structure.dat' u 1:2:10 w p pt 4 lc palette title 'Ax↑',\
     '2D_eight_band_structure.dat' u 1:3:11 w p pt 6 lc palette title 'Bx↑',\
     '2D_eight_band_structure.dat' u 1:4:12 w p pt 8 lc palette title 'Ay↑',\
     '2D_eight_band_structure.dat' u 1:5:13 w p pt 12 lc palette title 'By↑',\
     '2D_eight_band_structure.dat' u 1:6:14 w p pt 5 lc palette title 'Ax↓' ,\
     '2D_eight_band_structure.dat' u 1:7:15 w p pt 7 lc palette title 'Bx↓',\
     '2D_eight_band_structure.dat' u 1:8:16 w p pt 9 lc palette title 'Ay↓',\
     '2D_eight_band_structure.dat' u 1:9:17 w p pt 13 lc palette title 'By↓'
# --- GRAPH b
unset title

unset yrange
set size ratio 1
set ytics autofreq 

set key out below

unset xrange
set xtics autofreq 

unset arrow

set ylabel 'Σ(hybrid Wannier charge centres) / 2π'
set xlabel 'k_x / b_1'

plot 'eight_wilson_loop.dat' i 0 u 1:2 pt 4 title 'Ax↑',\
     'eight_wilson_loop.dat' i 1 u 1:2 pt 6 title 'Bx↑',\
     'eight_wilson_loop.dat' i 2 u 1:2 pt 8 title 'Ay↑',\
     'eight_wilson_loop.dat' i 3 u 1:2 pt 12 title 'By↑',\
     'eight_wilson_loop.dat' i 4 u 1:2 pt 5 title 'Ax↓',\
     'eight_wilson_loop.dat' i 5 u 1:2 pt 7 title 'Bx↓',\
     'eight_wilson_loop.dat' i 6 u 1:2 pt 9 title 'Ay↓',\
     'eight_wilson_loop.dat' i 7 u 1:2 pt 13 title 'By↓'
unset multiplot
### End multiplot
