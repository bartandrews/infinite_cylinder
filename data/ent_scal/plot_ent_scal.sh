#!/bin/bash

gnuplot -persist <<-EOFMarker

set key box out right
set title noenhanced '$1'

set xlabel 'cylinder circumference, L_y'
set ylabel 'entanglement entropy, S'

set xrange [0:]

f(x) = a*x + b
fit f(x) '$1' u 1:2 via a, b
title_f(a,b) = sprintf('S_{vN} = %.2f L_y + %.2f', a, b)

g(x) = c*x + d
fit g(x) '$1' u 1:3 via c, d
title_g(c,d) = sprintf('S_{inf} = %.2f L_y + %.2f', c, d)

#theoretical gamma in the Abelian case
sum_square_sites=4
gamma=-0.5*log(sum_square_sites)
title_gamma = sprintf('γ = -0.5 log(%i) = %.2f', sum_square_sites, gamma)

plot '$1' u 1:2 w p lc 1 title 'S_{vN}', f(x) lc 1 title title_f(a,b),\
'$2' u 1:3 w p pt 6 lc 2 title 'S_{inf}', g(x) lc 2 title title_g(c,d),\
gamma w l lc black dashtype 2 title title_gamma

EOFMarker
