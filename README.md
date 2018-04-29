# 2018-optimal-transport-prj

This repository contains my master thesis on Optimal Transport problem. 
Quickly, the optimal transport theory is the theory dealing with finding the shortest path between two probability density. 

For more detailed information see : http://www.gpeyre.com/ 

## Repository description 
- Rapport :
  - biblio.bib : The sources for this project
  - rapport.pdf : The result of this work (french)
  - soutenance.pdf : The presentation of this thesis
- src : 
  - CG 1D/2D : The algorithm implemented with conjugate gradient step for convergence, for 1D or 2D density (matlab)
  - FFFT 1D/2D : The algorithm implemented with the fast fourier transform for the projection on constraint step (Fortran)
  - Fortran 1D/2D : The algorithm implemented with conjugate gradient step for convergence, for 1D or 2D density (Fortran)
  - Fortran S 1D/2D : The algorithm implemented with conjugate gradient step for convergence and staggered grid discretization (Fortran)
  
## Example
For some example on 2D density, you can watch : https://www.youtube.com/playlist?list=PLJ92u2ph2rW-cRQQNv39aexr4nDht8Q_I
In this playlist, I apply optimal transport to finding shortest path in a labyrinth, morphing between image ..

## Credits
This work is under GNU GENERAL PUBLIC LICENSE, you can use this work with suitable citation. 

## Contact 
If any trouble detected, any help needed or any contribution please contact : timothee.schmoderer -at- netcourrier.com
