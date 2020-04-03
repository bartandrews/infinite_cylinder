set ylabel 'Berry flux along Wilson loop / 2π'
set xlabel 'k_x / b_1'

set key out right

plot 'eight_wilson_loop.dat' i 0 u 1:2 pt 4 title 'Ax↑',\
     'eight_wilson_loop.dat' i 1 u 1:2 pt 6 title 'Bx↑',\
     'eight_wilson_loop.dat' i 2 u 1:2 pt 8 title 'Ay↑',\
     'eight_wilson_loop.dat' i 3 u 1:2 pt 12 title 'By↑',\
     'eight_wilson_loop.dat' i 4 u 1:2 pt 5 title 'Ax↓',\
     'eight_wilson_loop.dat' i 5 u 1:2 pt 7 title 'Bx↓',\
     'eight_wilson_loop.dat' i 6 u 1:2 pt 9 title 'Ay↓',\
     'eight_wilson_loop.dat' i 7 u 1:2 pt 13 title 'By↓'
