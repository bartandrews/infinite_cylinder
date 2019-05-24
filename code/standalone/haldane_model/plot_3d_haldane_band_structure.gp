set title 'Haldane Model Band Structure (3D)'
set xlabel 'kx'
set ylabel 'ky'
set zlabel 'Energy / meV'
set cblabel '|C_A|^2'

unset key

set size square

set xtics ('-2π/3' -2*pi/3*0.97, 0, '2π/3' 2*pi/3*0.97)
set ytics ('-2π/(3*sqrt(3))' -2*pi/(3*sqrt(3)), 0, '2π/(3*sqrt(3))' 2*pi/(3*sqrt(3)))

splot '3D_haldane_band_structure.txt' u 1:2:4:6 w p lc palette,\
      '3D_haldane_band_structure.txt' u 1:2:3:5 w p lc palette 
