#!/bin/bash

files="$(ls -1v results/ | grep .dat)"

for file in $files; do 
	postFile="results/${file/.dat/.png}"
	file="results/${file}"

	gnuplot <<- EOF 
		set term png
		
		unset colorbox
		unset tics
		
		
		set output "${postFile}"
		plot "${file}" matrix with image
	EOF
done

ffmpeg  -framerate 10 -i "results/%05d.png"  transport.mp4 -y

#rm Transport/*
