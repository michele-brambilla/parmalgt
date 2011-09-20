#!/usr/bin/gnuplot
# Plot dei dati generati con il programma verifica

set terminal postscript color enhanced

set output "./result/operazioni.eps"
set ylabel "Gflops"
set xlabel "log_{10} numero operazioni"

set multiplot
set size 1, .5

set origin 0.0,0.0
set title "Prodotto matrici"
p [1:9] './result/res.dat' u 1:(($4)*10**(-9)) index 2 w linespoint t "scalari",\
  	'./result/res.dat' u 1:(($4)*10**(-9)) index 3 w linespoint t "perturbative"

set origin 0.0,0.5
set title "Somma matrici"
p [1:9] './result/res.dat' u 1:(($4)) index 0 w linespoint t "scalari",\
  	'./result/res.dat' u 1:(($4)) index 1 w linespoint t "perturbative"
unset multiplot




set output "./result/tempi.eps"

set log y
set ylabel "tempo [s]"
set key bottom
set multiplot
set size .5, .5

set origin 0.0,0.0
set title "Prodotto matrici scalari"
p [1:9] './result/res.dat' u 1:2 index 2 w linespoint t "clock()",\
  	'./result/res.dat' u 1:3 index 2 w linespoint t "gettimeofday()"

set origin 0.5,0.0
set title "Prodotto matrici pt"
p [1:9] './result/res.dat' u 1:2 index 3 w linespoint t "clock()",\
  	'./result/res.dat' u 1:3 index 3 w linespoint t "gettimeofday()"

set origin 0.0,0.5
set title "Somma matrici scalari"
p [1:9] './result/res.dat' u 1:2 index 0 w linespoint t "clock()",\
  	'./result/res.dat' u 1:3 index 0 w linespoint t "gettimeofday()"

set origin 0.5,0.5
set title "Somma matrici pt"
p [1:9] './result/res.dat' u 1:2 index 1 w linespoint t "clock()",\
  	'./result/res.dat' u 1:3 index 1 w linespoint t "gettimeofday()"

unset multiplot



set output "./result/confronto.eps"

unset log y
set ylabel "Gflops"
set key right bottom
set multiplot
set size .5, .5

set origin 0.0,0.0
set title "Prodotto matrici scalari"
p [1:9] './result/res.dat' u 1:(($4)*10**(-9)) index 2 w linespoint t "clock()",\
  	'./result/res.dat' u 1:(($5)*10**(-9)) index 2 w linespoint t "gettimeofday()"

set origin 0.5,0.0
set title "Prodotto matrici perturbative"
p [1:9] './result/res.dat' u 1:(($4)*10**(-9)) index 3 w linespoint t "clock()",\
  	'./result/res.dat' u 1:(($5)*10**(-9)) index 3 w linespoint t "gettimeofday()"

set origin 0.0,0.5
set title "Somma matrici scalari"
p [1:9] './result/res.dat' u 1:(($4)) index 0 w linespoint t "clock()",\
  	'./result/res.dat' u 1:(($5)) index 0 w linespoint t "gettimeofday()"

set origin 0.5,0.5
set title "Somma matrici perturbative"
p [1:9] './result/res.dat' u 1:(($4)) index 1 w linespoint t "clock()",\
  	'./result/res.dat' u 1:(($5)) index 1 w linespoint t "gettimeofday()"
unset multiplot


!ps2pdf ./result/operazioni.eps ./result/operazioni.pdf
!ps2pdf ./result/tempi.eps ./result/tempi.pdf
!ps2pdf ./result/confronto.eps ./result/confronto.pdf
!rm ./result/operazioni.eps
!rm ./result/tempi.eps
!rm ./result/confronto.eps

