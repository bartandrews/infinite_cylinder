### Start multiplot (1x2 layout)
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
#set title 'A/B character'
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

plot '2D_two_band_structure.dat' u 1:2:4 w p pt 4 lc palette title 'A',\
     '2D_two_band_structure.dat' u 1:3:5 w p pt 6 lc palette title 'B'

# --- GRAPH b
unset title

unset yrange
set size ratio 1
set ytics autofreq 

set key out below

unset xrange
set xtics autofreq 

unset arrow

set ylabel 'Σ(HWCC) / 2π'
set xlabel 'k_x / b_1'

plot 'two_wilson_loop.dat' i 0 u 1:2 pt 4 title 'A',\
     'two_wilson_loop.dat' i 1 u 1:2 pt 6 title 'B'

# --- GRAPH c
set size ratio 0.25
set size 1,1
set origin 0, -0.25
unset key
set ylabel '# States'
set xlabel 'Energy / meV'

# Each bar is half the (visual) width of its x-range.
set boxwidth 0.5 absolute
set style fill solid 1.0 noborder

bin_width = 0.5;

bin_number(x) = floor(x/bin_width)

rounded(x) = bin_width * ( bin_number(x) + 0.5 )

plot 'two_energies.dat' using (rounded($1)):(1) smooth frequency with boxes

unset multiplot
### End multiplot
