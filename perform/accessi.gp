#! /bin/gnuplot

set terminal x11 enhanced

set xlabel "log_{10} numero ripetizioni"
set ylabel "tempo [s]"

set multiplot
set size .5, 1

set origin 0, 0
set title "Velocita' di accesso ai dati - Umu.get()"
p 'accessi_m.dat' index 0 w l t "taglia 2^4", \
  'accessi_m.dat' index 1 w l t "taglia 4^4", \
  'accessi_m.dat' index 2 w l t "taglia 6^4", \
  'accessi_m.dat' index 3 w l t "taglia 8^4"

set origin .5, 0
set title "Velocita' di accesso ai dati - U[get()]"
p 'accessi_f.dat' index 0 w l t "taglia 2^4", \
  'accessi_f.dat' index 1 w l t "taglia 4^4", \
  'accessi_f.dat' index 2 w l t "taglia 6^4", \
  'accessi_f.dat' index 3 w l t "taglia 8^4"

unset multiplot

pause -1

set terminal postscript color enhanced
set output "accessi.eps"

set multiplot
set size .5, 1

set origin 0, 0
set title "Velocita' di accesso ai dati - Umu.get()"
p 'accessi_m.dat' index 0 w l t "taglia 2^4", \
  'accessi_m.dat' index 1 w l t "taglia 4^4", \
  'accessi_m.dat' index 2 w l t "taglia 6^4", \
  'accessi_m.dat' index 3 w l t "taglia 8^4"

set origin .5, 0
set title "Velocita' di accesso ai dati - U[get()]"
p 'accessi_f.dat' index 0 w l t "taglia 2^4", \
  'accessi_f.dat' index 1 w l t "taglia 4^4", \
  'accessi_f.dat' index 2 w l t "taglia 6^4", \
  'accessi_f.dat' index 3 w l t "taglia 8^4"

unset multiplot

!ps2pdf accessi.eps
!rm accessi.eps
