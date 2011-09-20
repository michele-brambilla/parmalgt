#!/usr/bin/gnuplot
# Plot dei dati generati con il programma verifica

set terminal postscript color enhanced
set output "verifica-gc.eps"
set ylabel "Gflops"
set xlabel "log_{10} numero operazioni"

set multiplot
set size 1, .5

set origin 0.0,0.0
set title "Prodotto matrici"
p [4.5:8.5] 'res280109.dat' u 1:(($4)*10**(-9)) index 2 w linespoint t "scalari",\
  	'res280109.dat' u 1:(($4)*10**(-9)) index 3 w linespoint t "perturbative"

set origin 0.0,0.5
set title "Somma matrici"
p [4.5:8.5] 'res280109.dat' u 1:(($4)) index 0 w linespoint t "scalari",\
  	'res280109.dat' u 1:(($4)) index 1 w linespoint t "perturbative"
unset multiplot

set terminal x11

set terminal postscript color enhanced
set output "verifica-t64.eps"
set ylabel "Gflops"
set xlabel "log_{10} numero operazioni"

set multiplot
set size 1, .5

set origin 0.0,0.0
set title "Prodotto matrici"
p [4.5:8.5] 'res280109-t64.dat' u 1:(($4)*10**(-9)) index 2 w linespoint t "scalari",\
  	'res280109-t64.dat' u 1:(($4)*10**(-9)) index 3 w linespoint t "perturbative"

set origin 0.0,0.5
set title "Somma matrici"
p [4.5:8.5] 'res280109-t64.dat' u 1:(($4)) index 0 w linespoint t "scalari",\
  	'res280109-t64.dat' u 1:(($4)) index 1 w linespoint t "perturbative"
unset multiplot


!ps2pdf verifica-gc.eps verifica-gc.pdf
!rm verifica-gc.eps
!ps2pdf verifica-t64.eps verifica-t64.pdf
!rm verifica-t64.eps

