set ylabel 'Berry flux along Wilson loop / 2π'
set xlabel 'k_x / b_1'

set key out right

plot 'two_wilson_loop.dat' i 0 u 1:2 pt 4 title 'A',\
     'two_wilson_loop.dat' i 1 u 1:2 pt 6 title 'B'
