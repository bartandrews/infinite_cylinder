set title 'Haldane Model Band Structure (3D)'
set xlabel 'kx'
set ylabel 'ky'
set zlabel 'Energy / meV'
set cblabel '|C_A|^2'

K1x = 2*pi/3
K1y = 2*pi/(3*sqrt(3))

set arrow 1 from K1x,K1y,0 to 0,2*K1y,0 nohead
set arrow 2 from 0,2*K1y,0 to -K1x,K1y,0 nohead
set arrow 3 from -K1x,K1y,0 to -K1x,-K1y,0 nohead
set arrow 4 from -K1x,-K1y,0 to 0,-2*K1y,0 nohead
set arrow 5 from 0,-2*K1y,0 to K1x,-K1y,0 nohead
set arrow 6 from K1x,-K1y,0 to K1x,K1y,0 nohead

set label "Γ" at 0,0,0
set label "M" at K1x,0,0
set label "K" at K1x,K1y,0
set label "K'" at K1x,-K1y,0
 
unset key

set size square

set xtics ('-2π/3' -2*pi/3*0.97, 0, '2π/3' 2*pi/3*0.97)
set ytics ('-2π/(3√3)' -2*pi/(3*sqrt(3)), 0, '2π/(3√3)' 2*pi/(3*sqrt(3)))

splot '3D_haldane_band_structure.txt' u 1:2:4:6 w p lc palette,\
      '3D_haldane_band_structure.txt' u 1:2:3:5 w p lc palette 
