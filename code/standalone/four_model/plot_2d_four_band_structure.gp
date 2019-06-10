set title 'Four Model Band Structure (2D)'
set ylabel 'Energy / meV'
set cblabel '|C_A|^2'

set size ratio 0.4
set ytics 2

set key out right

set xrange [0: 300]
set xtics ('K' 0, 'Γ' 100, 'M' 200, "K\'" 300)

set arrow from 100,graph(0,0) to 100,graph(1,1) nohead dt 2
set arrow from 200,graph(0,0) to 200,graph(1,1) nohead dt 2

plot '2D_four_band_structure.txt' u 1:2:6 w p lc palette title 'band 0',\
     '2D_four_band_structure.txt' u 1:3:7 w p lc palette title 'band 1',\
     '2D_four_band_structure.txt' u 1:4:8 w p lc palette title 'band 2',\
     '2D_four_band_structure.txt' u 1:5:9 w p lc palette title 'band 3'
