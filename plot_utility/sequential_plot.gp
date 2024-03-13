reset
dtype = "tab"

# -----------------------------------------------
# 1. Set font size
# -----------------------------------------------
set xlabel "x"  font ", 12"
set ylabel "y"  font ", 12"
set  tics       font ", 12"
set xtics       font ", 12"
set ytics       font ", 12"
set title       font ", 12"

# -----------------------------------------------
# 2. Set plot margins
# -----------------------------------------------
#set lmargin at screen 0.15
#set rmargin at screen 0.8
#set bmargin at screen 0.1
#set tmargin at screen 0.9

# -----------------------------------------------
# 3. Set Window size and aspect ratio
# -----------------------------------------------
#ywinsize = 960                      # Assume 640 pixel in y direction
#nx1 = 50
#nx2 = 50
#xwinsize = ywinsize*nx1*1.3/nx2     # extra factor 1.3 account for cbar
#set term qt size xwinsize, ywinsize

# -----------------------------------------------
# 4. Select the variable you want to plot, 
#    together with title
# -----------------------------------------------
# nvar      = prs         # = rho, vx1, vx2, vx3, prs, Bx1, Bx2, Bx3

rho = 3
vx1 = 4
vx2 = 5
vx3 = 6
prs = 7
bx1 = 8
bx2 = 9
bx3 = 10 

var_name  = "Pressure"

# -----------------------------------------------
# 5. Plot 
# -----------------------------------------------

set xrange[0:2*pi]
set yrange[0:2*pi]

do for [n=0:200]{
  
  datafile = sprintf ("data.%04d.out",n)   # File name
  #set table "contour_table.dat"
  set contour base
  set view 0,0
  unset table
  unset key
  plot datafile u 1:2:rho w image, #"contour_table.dat"# u 1:2:3 w l \
      lc rgb "gray"
  pause 0.5

  
}

