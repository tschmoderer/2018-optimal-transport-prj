#!/bin/bash

files="$(ls -1v Transport/ | grep .dat)"

for file in $files; do 
postFile="Transport/${file/.dat/.png}"
file="Transport/${file}"

gnuplot <<- EOF 
    set term png
	set title "Transport Optimal"
	set xr [0:1]
	set yr [0:1]
	set xlabel "x"
	set ylabel "y"
 	set zlabel "f"
 	set contour
	set cntrparam levels 15
	unset key
	set pm3d
#	set palette gray
	unset colorbox
	set hidden3d
	load "plot.gnu"
    set output "${postFile}"
    splot "${file}" with lines
#	plot "${file}" matrix with image
EOF
done


ffmpeg -framerate 10 -pattern_type glob -i 'Transport/*.png' transport.mp4 -y

#rm Transport/*

gnuplot <<- EOF 
    set term png
	set title "Energie"
	set logscale
	set xlabel "Iteration"
	set ylabel "J(m,f)"
    set output "energie.png"
    unset key
    plot "data.dat" u 1:2 lt 4 with lines
EOF
