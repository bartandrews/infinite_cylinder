set title 'Band Structure'

set size square

set xtics ('-π' -pi, 0, 'π' pi)
set ytics ('-π' -pi, 0, 'π' pi)

set xlabel 'kx'
set ylabel 'ky'
set zlabel 'E'

splot 'band_structure.txt' u 1:2:3 title 'a+',\
      'band_structure.txt' u 1:2:4 title 'a-',\
      'band_structure.txt' u 1:2:5 title 'b+',\
      'band_structure.txt' u 1:2:6 title 'b-',\
      'band_structure.txt' u 1:2:7 title 'c+',\
      'band_structure.txt' u 1:2:8 title 'c-',\
      'band_structure.txt' u 1:2:9 title 'd+',\
      'band_structure.txt' u 1:2:10 title 'd-' 
