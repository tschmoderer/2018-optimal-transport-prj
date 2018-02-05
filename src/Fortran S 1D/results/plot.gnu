 set contour
 set cntrparam levels 30
 unset key
 set pm3d
 unset colorbox
 set hidden3d
 set title "Transport Optimal"
 set xlabel "x"
 set ylabel "t"
 set dgrid3d           30 ,          32
 splot "results/transport.dat" with lines
