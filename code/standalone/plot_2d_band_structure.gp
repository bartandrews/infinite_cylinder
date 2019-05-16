set title '2D Band Structure'

set xlabel 'path'
set ylabel 'E'

set xrange [-2*pi/3: 2*pi*((3+sqrt(3))/9)]
set xtics ('-K' -2*pi/3, 'Γ' 0, 'M' 2*pi/3, 'K' 2*pi*((3+sqrt(3))/9))

plot 'band_structure.txt' u ($1 > -2*pi/3 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.02 ? $3 : 1/0) title 'up x' ls 1,\
     'band_structure.txt' u ($1 > 0 && $1 < 2*pi/3 ? $1 : 1/0):(abs($2) < 0.05 ? $3 : 1/0) notitle ls 1,\
     'band_structure.txt' u ($2 > 0 && $2 < 2*pi/(3*sqrt(3)) ? $2+2*pi/3 : 1/0):(abs($1-2*pi/3) < 0.031734 ? $3 : 1/0) notitle ls 1,\
     'band_structure.txt' u ($1 > -2*pi/3 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.02 ? $5 : 1/0) title 'up y' ls 2,\
     'band_structure.txt' u ($1 > 0 && $1 < 2*pi/3 ? $1 : 1/0):(abs($2) < 0.05 ? $5 : 1/0) notitle ls 2,\
     'band_structure.txt' u ($2 > 0 && $2 < 2*pi/(3*sqrt(3)) ? $2+2*pi/3 : 1/0):(abs($1-2*pi/3) < 0.031734 ? $5 : 1/0) notitle ls 2,\
     'band_structure.txt' u ($1 > -2*pi/3 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.02 ? $7 : 1/0) title 'down x' ls 3,\
     'band_structure.txt' u ($1 > 0 && $1 < 2*pi/3 ? $1 : 1/0):(abs($2) < 0.05 ? $7 : 1/0) notitle ls 3,\
     'band_structure.txt' u ($2 > 0 && $2 < 2*pi/(3*sqrt(3)) ? $2+2*pi/3 : 1/0):(abs($1-2*pi/3) < 0.031734 ? $7 : 1/0) notitle ls 3,\
     'band_structure.txt' u ($1 > -2*pi/3 && $1 < 0 ? $1 : 1/0):(abs($2 - $1/sqrt(3)) < 0.02 ? $9 : 1/0) title 'down y' ls 4,\
     'band_structure.txt' u ($1 > 0 && $1 < 2*pi/3 ? $1 : 1/0):(abs($2) < 0.05 ? $9 : 1/0) notitle ls 4,\
     'band_structure.txt' u ($2 > 0 && $2 < 2*pi/(3*sqrt(3)) ? $2+2*pi/3 : 1/0):(abs($1-2*pi/3) < 0.031734 ? $9 : 1/0) notitle ls 4
