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
		unset key
 		set xlabel "x"
 		set ylabel "f"
    set output "${postFile}"
    plot "${file}" u 1:2 with lines
EOF
done

ffmpeg -framerate 1 -pattern_type glob -i 'Transport/*.png' transport.mp4 -y

gnuplot <<- EOF
	set term png
	set output "energie.png"
	set xlabel "Iteration"
	set ylabel "J(m,f)"
	set logscale
 	plot "data.dat" u 1:2 lt 3 with lines title 'Energie' 	 
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
	unset colorbox
	unset tics
	set xlabel "x"
	set ylabel "t"
	set title "Transport Optimal Généralisé"
	set pm3d map
	set palette rgb -21,-22,-23
	plot "transport.dat" matrix with image
EOF

gnuplot <<- EOF
	set term png
	set output "obstacle.png"
	load "plot.gnu"	 
	unset colorbox
	unset tics
	set xlabel "x"
	set ylabel "t"
	set title "Transport Optimal Généralisé"
	set pm3d map
	set palette model RGB defined(0 "white", 1 "royalblue")
	plot "obstacle.dat" matrix with image
EOF

convert obstacle.png -transparent white obstacle.png 
composite -gravity center obstacle.png transport.png resultat.png



#rm Transport/*
