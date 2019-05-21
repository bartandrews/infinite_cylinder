set title '1.20 degrees'

#set xlabel 'path'
set ylabel 'Energy / meV'

unset key

set size ratio 0.25
set ytics 2

set cblabel '|C_A|^2'

set xrange [0: 606]
set xtics ('K' 0, 'Γ' 202, 'M' 404, "K\'" 606)

set arrow from 202,graph(0,0) to 202,graph(1,1) nohead dt 2
set arrow from 404,graph(0,0) to 404,graph(1,1) nohead dt 2

plot 'band_structure.txt' u 1:2:4 w p lc palette,\
     'band_structure.txt' u 1:3:5 w p lc palette
