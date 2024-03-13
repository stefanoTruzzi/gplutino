# **********************************************************************
# Produce a simple plot for the advection code.
# **********************************************************************
reset

set xlabel "x" font ", 18"
set ylabel "y" font ", 18"
set tics       font ", 18"
set yrange[-0.05:1.2]
set xrange[-0.1:0.1]

do for [n=0:30]{
  
  datafile = sprintf ("data.%04d.out",n)   # File name
  plot datafile u 1:5 w lp    title "Density"
  
  pause 1
}
