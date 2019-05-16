set xlabel 'path'
set ylabel 'E'

unset key

set size ratio 0.125
set ytics 3

set xrange [-3*pi/4: ( 3*pi/4 + 3*pi/(4*sqrt(3)) )]
set xtics ('-K' -3*pi/4, 'Γ' 0, 'M' 3*pi/4, 'K' ( 3*pi/4 + 3*pi/(4*sqrt(3)) ))

set arrow from 0,graph(0,0) to 0,graph(1,1) nohead dt 2
set arrow from 3*pi/4,graph(0,0) to 3*pi/4,graph(1,1) nohead dt 2

plot 'band_structure.txt' u ($1 > -3*pi/4 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.035 ? $3 : 1/0) notitle ls 1,\
     'band_structure.txt' u ($1 > 0 && $1 < 3*pi/4 ? $1 : 1/0):(abs($2) < 0.05 ? $3 : 1/0) notitle ls 1,\
     'band_structure.txt' u ($2 > 0 && $2 < 3*pi/(4*sqrt(3)) ? $2+3*pi/4 : 1/0):(abs($1-3*pi/4) < 0.031734 ? $3 : 1/0) notitle ls 1,\
     'band_structure.txt' u ($1 > -3*pi/4 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.035 ? $4 : 1/0) notitle ls 1,\
     'band_structure.txt' u ($1 > 0 && $1 < 3*pi/4 ? $1 : 1/0):(abs($2) < 0.05 ? $4 : 1/0) notitle ls 1,\
     'band_structure.txt' u ($2 > 0 && $2 < 3*pi/(4*sqrt(3)) ? $2+3*pi/4 : 1/0):(abs($1-3*pi/4) < 0.031734 ? $4 : 1/0) notitle ls 1,\
     'band_structure.txt' u ($1 > -3*pi/4 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.035 ? $5 : 1/0) notitle ls 2,\
     'band_structure.txt' u ($1 > 0 && $1 < 3*pi/4 ? $1 : 1/0):(abs($2) < 0.05 ? $5 : 1/0) notitle ls 2,\
     'band_structure.txt' u ($2 > 0 && $2 < 3*pi/(4*sqrt(3)) ? $2+3*pi/4 : 1/0):(abs($1-3*pi/4) < 0.031734 ? $5 : 1/0) notitle ls 2,\
     'band_structure.txt' u ($1 > -3*pi/4 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.035 ? $6 : 1/0) notitle ls 2,\
     'band_structure.txt' u ($1 > 0 && $1 < 3*pi/4 ? $1 : 1/0):(abs($2) < 0.05 ? $6 : 1/0) notitle ls 2,\
     'band_structure.txt' u ($2 > 0 && $2 < 3*pi/(4*sqrt(3)) ? $2+3*pi/4 : 1/0):(abs($1-3*pi/4) < 0.031734 ? $6 : 1/0) notitle ls 2,\
     'band_structure.txt' u ($1 > -3*pi/4 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.035 ? $7 : 1/0) notitle ls 3,\
     'band_structure.txt' u ($1 > 0 && $1 < 3*pi/4 ? $1 : 1/0):(abs($2) < 0.05 ? $7 : 1/0) notitle ls 3,\
     'band_structure.txt' u ($2 > 0 && $2 < 3*pi/(4*sqrt(3)) ? $2+3*pi/4 : 1/0):(abs($1-3*pi/4) < 0.031734 ? $7 : 1/0) notitle ls 3,\
     'band_structure.txt' u ($1 > -3*pi/4 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.035 ? $8 : 1/0) notitle ls 3,\
     'band_structure.txt' u ($1 > 0 && $1 < 3*pi/4 ? $1 : 1/0):(abs($2) < 0.05 ? $8 : 1/0) notitle ls 3,\
     'band_structure.txt' u ($2 > 0 && $2 < 3*pi/(4*sqrt(3)) ? $2+3*pi/4 : 1/0):(abs($1-3*pi/4) < 0.031734 ? $8 : 1/0) notitle ls 3,\
     'band_structure.txt' u ($1 > -3*pi/4 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.035 ? $9 : 1/0) notitle ls 4,\
     'band_structure.txt' u ($1 > 0 && $1 < 3*pi/4 ? $1 : 1/0):(abs($2) < 0.05 ? $9 : 1/0) notitle ls 4,\
     'band_structure.txt' u ($2 > 0 && $2 < 3*pi/(4*sqrt(3)) ? $2+3*pi/4 : 1/0):(abs($1-3*pi/4) < 0.031734 ? $9 : 1/0) notitle ls 4,\
     'band_structure.txt' u ($1 > -3*pi/4 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.035 ? $10 : 1/0) notitle ls 4,\
     'band_structure.txt' u ($1 > 0 && $1 < 3*pi/4 ? $1 : 1/0):(abs($2) < 0.05 ? $10 : 1/0) notitle ls 4,\
     'band_structure.txt' u ($2 > 0 && $2 < 3*pi/(4*sqrt(3)) ? $2+3*pi/4 : 1/0):(abs($1-3*pi/4) < 0.031734 ? $10 : 1/0) notitle ls 4
