set title 'A/B character'
set ylabel 'Energy / meV'
set cblabel '|C_A|^2'

set size ratio 0.4
set ytics 2

set key out below

num_points = 3*30

set xrange [0: num_points]
set xtics ('K' 0, 'Γ' num_points/3, 'M' 2*num_points/3, "K\'" num_points)

set arrow from num_points/3,graph(0,0) to num_points/3,graph(1,1) nohead dt 2
set arrow from 2*num_points/3,graph(0,0) to 2*num_points/3,graph(1,1) nohead dt 2

plot '2D_two_band_structure.dat' u 1:2:4 w p pt 4 lc palette title 'A',\
     '2D_two_band_structure.dat' u 1:3:5 w p pt 6 lc palette title 'B'
