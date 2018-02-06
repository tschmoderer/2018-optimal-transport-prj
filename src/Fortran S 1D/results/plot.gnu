 set contour
 set cntrparam levels 30
 unset key
 set pm3d
 unset colorbox
 set hidden3d
 set title "Transport Optimal"
 set xlabel "x"
 set ylabel "t"
 set palette gray
 set view 0,0
 set dgrid3d           51 ,          51
 splot "results/transport.dat" with lines
