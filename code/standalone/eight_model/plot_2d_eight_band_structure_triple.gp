set ylabel 'Energy / meV'

set size ratio 0.2
set ytics 2
set cbtics 0.5
set cbrange [0:1]

set key out right

num_points = 3*30

set xrange [0: num_points]
set xtics ('K' 0, 'Γ' num_points/3, 'M' 2*num_points/3, "K\'" num_points)

set arrow from num_points/3,graph(0,0) to num_points/3,graph(1,1) nohead dt 2
set arrow from 2*num_points/3,graph(0,0) to 2*num_points/3,graph(1,1) nohead dt 2

### Start multiplot (2x2 layout)
set multiplot layout 3,1 rowsfirst
# --- GRAPH a
set title 'A/B character'
set cblabel '|C_A|^2'
plot '2D_eight_band_structure.dat' u 1:2:10 w p pt 4 lc palette title 'Ax↑',\
     '2D_eight_band_structure.dat' u 1:3:11 w p pt 6 lc palette title 'Bx↑',\
     '2D_eight_band_structure.dat' u 1:4:12 w p pt 8 lc palette title 'Ay↑',\
     '2D_eight_band_structure.dat' u 1:5:13 w p pt 12 lc palette title 'By↑',\
     '2D_eight_band_structure.dat' u 1:6:14 w p pt 5 lc palette title 'Ax↓' ,\
     '2D_eight_band_structure.dat' u 1:7:15 w p pt 7 lc palette title 'Bx↓',\
     '2D_eight_band_structure.dat' u 1:8:16 w p pt 9 lc palette title 'Ay↓',\
     '2D_eight_band_structure.dat' u 1:9:17 w p pt 13 lc palette title 'By↓'
# --- GRAPH b
set title 'x/y character'
set cblabel '|C_x|^2'
plot '2D_eight_band_structure.dat' u 1:2:18 w p pt 4 lc palette title 'Ax↑',\
     '2D_eight_band_structure.dat' u 1:3:19 w p pt 6 lc palette title 'Bx↑',\
     '2D_eight_band_structure.dat' u 1:4:20 w p pt 8 lc palette title 'Ay↑',\
     '2D_eight_band_structure.dat' u 1:5:21 w p pt 12 lc palette title 'By↑',\
     '2D_eight_band_structure.dat' u 1:6:22 w p pt 5 lc palette title 'Ax↓' ,\
     '2D_eight_band_structure.dat' u 1:7:23 w p pt 7 lc palette title 'Bx↓',\
     '2D_eight_band_structure.dat' u 1:8:24 w p pt 9 lc palette title 'Ay↓',\
     '2D_eight_band_structure.dat' u 1:9:25 w p pt 13 lc palette title 'By↓'
# --- GRAPH c
set title '↑/↓ character'
set cblabel '|C_↑|^2'
plot '2D_eight_band_structure.dat' u 1:2:26 w p pt 4 lc palette title 'Ax↑',\
     '2D_eight_band_structure.dat' u 1:3:27 w p pt 6 lc palette title 'Bx↑',\
     '2D_eight_band_structure.dat' u 1:4:28 w p pt 8 lc palette title 'Ay↑',\
     '2D_eight_band_structure.dat' u 1:5:29 w p pt 12 lc palette title 'By↑',\
     '2D_eight_band_structure.dat' u 1:6:30 w p pt 5 lc palette title 'Ax↓' ,\
     '2D_eight_band_structure.dat' u 1:7:31 w p pt 7 lc palette title 'Bx↓',\
     '2D_eight_band_structure.dat' u 1:8:32 w p pt 9 lc palette title 'Ay↓',\
     '2D_eight_band_structure.dat' u 1:9:33 w p pt 13 lc palette title 'By↓'
unset multiplot
### End multiplot
