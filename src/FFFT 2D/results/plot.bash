#!/bin/bash

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
	load "plot.gnu"
	set xlabel "x"
	set ylabel "y"
	set xtics 0,0.1,1
	unset tics 
	set pm3d
	set palette rgb -21,-22,-23
	set palette gray
	set contour
	set cntrparam levels 15
	unset colorbox
	unset key
	set hidden3d
	unset surface
	set title 'Densité initiale' 	 
 	plot "f0.dat" matrix with image
EOF

#F1
gnuplot <<- EOF
	set term png
	set output "f1.png"
	load "plot.gnu"
	set xlabel "x"
	set ylabel "y"
	set xtics 0,0.1,1
	unset tics 
	set pm3d
	set palette rgb -21,-22,-23
	set palette gray
	set contour
	set cntrparam levels 15
	unset colorbox
	unset key
	set hidden3d
	unset surface
	set title 'Densité finale' 	 
 	plot "f1.dat" matrix with image
EOF

#traitement obstacle

files="$(ls -1v Obstacle/O_* | grep .dat)"

for file in $files; do 
postFile="${file/.dat/.png}"
file="${file}"
gnuplot <<- EOF 
  set term png
  set output "${postFile}"
  
  load "plot.gnu"
  load "plotC.gnu"
  
	set title "Transport Optimal Généralisé"
	set xlabel "x"
	set ylabel "y"
	unset tics
	set cbrange [0:1]
	set pm3d map
	set palette model RGB defined(0 "white", 1 "royalblue")
	unset colorbox
	unset key
	set hidden3d
	
  plot "${file}" matrix with image
EOF
convert ${postFile} -transparent white ${postFile}
done


#Morphing
files="$(ls -1v Transport/C_* | grep .dat)"

for file in $files; do 
postFile="${file/.dat/.png}"
file="${file}"

gnuplot <<- EOF 
  set term png
  set output "${postFile}"
  
  load "plot.gnu"
  load "plotC.gnu"
  
	set title "Transport Optimal Généralisé"
	set xlabel "x"
	set ylabel "y"
	unset tics

	set pm3d map
	set palette gray
	unset colorbox
	unset key

	set hidden3d
	
  plot "${file}" matrix with image
EOF
done

ffmpeg -framerate 10 -pattern_type glob -i 'Transport/C_*.png' transport_morphing.mp4 -y

#deuxième film
files="$(ls -1v Obstacle/T_* | grep .dat)"

for file in $files; do 
postFile="${file/.dat/.png}"
file="${file}"
#recuperer l'image d'obstacle qui a la meme valeur que nous 
number="${file//[^0-9]}"
obstacle="$(ls -1v Obstacle/O_* | grep ${number}.png)"

gnuplot <<- EOF 
  set term png
  set output "${postFile}"
  
  load "plot.gnu"
  load "plotC.gnu"
  
	set title "Transport Optimal Généralisé"
	set xlabel "x"
	set ylabel "y"
	unset tics
	
	set pm3d map
	set palette rgb -21,-22,-23
	unset colorbox
	unset key

	set hidden3d
	
  plot "${file}" matrix with image
EOF
composite -gravity center ${obstacle} ${postFile} ${postFile}
done

ffmpeg -framerate 10 -pattern_type glob -i 'Obstacle/T_*.png' transport_obstacle.mp4 -y

#Make films
files="$(ls -1v Transport/3D_* | grep .dat)"

for file in $files; do 
postFile="${file/.dat/.png}"
file="${file}"

gnuplot <<- EOF 
  set term png
  set output "${postFile}"
  load "plot.gnu"
  load "plot3D.gnu"
  load "plotC.gnu"
	set title "Transport Optimal Généralisé"
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
	unset colorbox
	set hidden3d  
  splot "${file}" with lines
EOF
done

ffmpeg -framerate 10 -pattern_type glob -i 'Transport/3D_*.png' transport.mp4 -y


#deuxième film
files="$(ls -1v Transport/C_* | grep .dat)"

for file in $files; do 
postFile="${file/.dat/.png}"
file="${file}"

gnuplot <<- EOF 
  set term png
  set output "${postFile}"
  
  load "plot.gnu"
  load "plotC.gnu"
  
	set title "Transport Optimal Généralisé"
	set xlabel "x"
	set ylabel "y"
	unset tics
	
	set pm3d map
	set palette rgb -21,-22,-23
	unset colorbox
	unset key

	set hidden3d
	
  plot "${file}" matrix with image
EOF
done

ffmpeg -framerate 10 -pattern_type glob -i 'Transport/C_*.png' transport_contour.mp4 -y























