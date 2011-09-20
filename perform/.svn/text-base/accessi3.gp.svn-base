# !/usr/bin/gnuplot

set terminal aqua enhanced

set multiplot
set size 1,.5
set xlabel "Taglia"

set ylabel "t [ns]"
set origin 0, .5
set title "Velocita' di accesso ai dati"
p './result/accessi_m3.dat' u 1:($4)*1e9 w l t "metodo", './result/accessi_f3.dat' u 1:($4)*1e9 w l t "friend"

set ylabel "t_{accesso}/t_{tot}"
set origin  0, 0
set title "Peso relativo dell'accesso"

p './result/accessi_m3.dat' u 1:((($3) - ($2))/($3)) w l t "metodo", './result/accessi_f3.dat' u 1:((($3)-($2))/($3)) w l t "friend"

pause -1
