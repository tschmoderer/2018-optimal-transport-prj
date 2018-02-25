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
	set palette rgb -21,-22,-23
#	set palette gray
	unset colorbox
	set hidden3d
	load "plot1.gnu"
  set output "${postFile}"
  splot "${file}" with lines
#	plot "${file}" matrix with image
EOF
done

ffmpeg -framerate 10 -pattern_type glob -i 'Transport/*.png' transport.mp4 -y


#deuxième film
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
	unset ztics
 	set contour
	set cntrparam levels 15
	unset key
#	set pm3d
#	set palette rgb -21,-22,-23
#	set palette gray
#	unset colorbox
	set dgrid3d 100,100
#	set hidden3d
#	unset surface
	set view 0,0
	unset surface
  set output "${postFile}"
  splot "${file}" with lines
EOF
done

ffmpeg -framerate 10 -pattern_type glob -i 'Transport/*.png' transport_contour.mp4 -y


# Energie
gnuplot <<- EOF 
    set term png
		set title "Energie"
		set logscale
		set xlabel "Iteration"
		set ylabel "J(m,f)"
    set output "energie.png"
    unset key
    plot "data.dat" u 1:2 lt 3 with lines
EOF

#F0
gnuplot <<- EOF
	set term png
	set output "f0.png"
	set xlabel "x"
	set ylabel "y"
	set pm3d map
	set palette rgb -21,-22,-23
	set contour
	set cntrparam levels 15
	unset colorbox
	unset key
	set dgrid3d 100,100
	set hidden3d
	unset surface
#	set view 0,0
	set title 'Densité initiale' 	 
 	plot "f0.dat" matrix with image
EOF

#F1
gnuplot <<- EOF
	set term png
	set output "f1.png"
	set xlabel "x"
	set ylabel "y"
	set pm3d map
	set palette rgb -21,-22,-23
	set contour
	set cntrparam levels 15
	unset colorbox
	unset key
	set dgrid3d 32,32
	set hidden3d
	unset surface
	set  title 'Densité finale'
 	splot "f1.dat" with lines 
EOF
































