#!/bin/bash

files="$(ls -1v Transport/ | grep .dat)"

for file in $files; do 
postFile="Transport/${file/.dat/.png}"
file="Transport/${file}"

gnuplot <<- EOF 
    set term png
		set title "Transport Optimal"
		set xr [0:1]
		set yr [-0.1:0.5]
 		set xlabel "x"
 		set ylabel "f"
    set output "${postFile}"
    plot "${file}" u 1:2 with lines
EOF
done


ffmpeg -framerate 10 -pattern_type glob -i 'Transport/*.png' transport.mp4 -y

gnuplot <<- EOF 
#    set term png
    set multiplot layout 3,1 rowsfirst   
 	set xlabel "iter"
 	plot "data.dat" u 1:2 with lines title 'Energie' 	 
 	plot "data.dat" u 1:3 with lines title 'minF' 	
 	plot "data.dat" u 1:4 with lines title 'divV'
#    set output "data.png"
EOF




#rm Transport/*
