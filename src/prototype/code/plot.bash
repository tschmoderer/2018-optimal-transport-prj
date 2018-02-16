#!/bin/bash

files="$(ls -1v results/ | grep .dat)"

for file in $files; do 
postFile="results/${file/.dat/.png}"
file="results/${file}"

gnuplot <<- EOF 
    set term png
	set pm3d
	set palette gray
	unset colorbox
    set output "${postFile}"
    plot "${file}" matrix with image
EOF
done


ffmpeg -framerate 10 -pattern_type glob -i "results/*.png" transport.mp4 -y

#rm Transport/*
