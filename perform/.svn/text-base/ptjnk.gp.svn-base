#! /usr/gnuplot
# Plot del calcolo dei flops nel programma

set terminal postscript color enhanced
set output "prlgt_flops.eps"

set title "Prestazione PTjnk.exe"
set xlabel "taglia reticolo"
set ylabel "Gflops"

p "ptjnk.dat" w l t "ordine 8", "ptjnk6.dat" w l t "ordine 6"
# "ptjnk.dat" w l t "ordine ###"
