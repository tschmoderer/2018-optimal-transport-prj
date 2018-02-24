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
	set term png
	set output "energie.png"
	set xlabel "Iteration"
	set ylabel "J(m,f)"
	set logscale
 	plot "data.dat" u 1:2 lt 4 with lines title 'Energie' 	 
EOF

gnuplot <<- EOF
	set term png
	set output "f0.png"
	set xlabel "x"
	set ylabel "Densité"
 	plot "f0.dat" u 1:2 lt -1 with lines title 'Densité initiale' 	 
EOF

gnuplot <<- EOF
	set term png
	set output "f1.png"
	set xlabel "x"
	set ylabel "Densité"
 	plot "f1.dat" u 1:2 lt -1 with lines title 'Densité finale' 	 
EOF

gnuplot <<- EOF
	set term png
	set output "transport.png"
	load "plot.gnu"	 
	unset tics
	splot "transport.dat" with lines
EOF





























#rm Transport/*
