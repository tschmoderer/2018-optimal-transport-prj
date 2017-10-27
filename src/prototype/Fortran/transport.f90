program transport

implicit none
!! Les basiques !!
integer :: i,j,k;

!!! Initialisation !!!
!! L'espace temps !!
integer, parameter :: d = 1; ! Dimension du sytème
integer, parameter :: N = 50; ! Nb de points dans la direction x
integer, parameter :: Q = 50; ! Nb de points dans la direction t

!! Paramètres de la méthode !!
real :: alpha = 1.0; ! Doit etre dans ]0,2[
real :: beta = 1.0; ! Doit etre dans [0,1]
real :: gamma = 0.5; ! Doit etre positif

!! Grille centrée !!
real, dimension(0:N) :: GcX;
real, dimension(0:Q) :: GcT;

!! Grille décentrée !!
real, dimension(-1:N) :: GsX;
real, dimension(-1:Q) :: GsT;

!! Variables centrées !!
real, dimension(Q+1,N+1) :: m = 0,f = 0;

!! Variables décentrées !!
real, dimension(Q+1,N+2) :: mbar = 0;
real, dimension(Q+2,N+1) :: fbar = 0;

!! Les valeurs aux frontières !!
real, dimension(N+1) :: f0;
real, dimension(N+1) :: f1;



!!! Initialisation des variables !!!

!! Grille centrée !!
do i=0,N
    GcX(i) = i/(1.0*N);
end do
do j=0,Q
    GcT = j/(1.0*Q);
end do

!! Grille décentrée !!
do i=-1,N
    GsX(i) = (i+0.5)/(1.0*N);
end do
do j=-1,Q
    GsT(i) = (j+0.5)/(1.0*Q);
end do

!! Les fontières !!
call finitial(GcX,N,f0);
call ffinal(GcX,N,f1);

!! Ecrire les données dans les fichiers !!

open (unit=0,file="files/f0");
open (unit=1,file="files/f1");
open (unit=2,file="files/parameters");
write (0,*), f0;
write (1,*), f1;
write (2,*), N, Q;

! call system('cd .. && FreeFem++ poisson_2d.pde && cd Fortran');


!!! Zone de tests !!!

call RANDOM_NUMBER(m);
call RANDOM_NUMBER(f);
call RANDOM_NUMBER(mbar);
call RANDOM_NUMBER(fbar);



!!! Fin Zone de tests !!!
end program