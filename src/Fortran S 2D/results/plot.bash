#!/bin/bash

files="$(ls -1v Transport/ | grep .dat)"

for file in $files; do 
postFile="Transport/${file/.dat/.png}"
file="Transport/${file}"

gnuplot <<- EOF 
    set term png
		set title "Transport Optimal"
		set xr [0:1]
		set yr [-0.1:0.2]
 		set xlabel "x"
 		set ylabel "f"
 		load "plot.gnu"
    set output "${postFile}"
    splot "${file}"
EOF
done


ffmpeg -framerate 10 -pattern_type glob -i 'Transport/*.png' transport.mp4 -y

#rm Transport/*
