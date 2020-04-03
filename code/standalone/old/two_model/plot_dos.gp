unset key
set ylabel '# States'
set xlabel 'Energy / meV'

set size ratio 0.1

# Each bar is half the (visual) width of its x-range.
set boxwidth 0.5 absolute
set style fill solid 1.0 noborder

bin_width = 0.5;

bin_number(x) = floor(x/bin_width)

rounded(x) = bin_width * ( bin_number(x) + 0.5 )

plot 'two_energies.dat' using (rounded($1)):(1) smooth frequency with boxes
