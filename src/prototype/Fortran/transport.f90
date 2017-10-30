program transport
implicit none

!! Les basiques -  indices de boucles!!
integer :: i,j,k;

!!! Initialisation !!!
!! L'espace temps !!
integer, parameter :: N = 70; ! Nb de points dans la direction x
integer, parameter :: Q = 50; ! Nb de points dans la direction t

!   1
!   |
!   |
!   |
!   t  Q
!   |
!   |
!   |        N
!   0 ------ x ------ 1

!! Paramètres de la méthode !!
double precision :: alpha = 1; ! Doit etre dans ]0,2[
double precision :: beta = 1.0; ! Doit etre dans [0,1]
double precision :: gamma = 0.5; ! Doit etre positif

!! Grille centrée !!
double precision, dimension(0:N) :: GcX;

!! Variables centrées !!
double precision, dimension(N+1,Q+1) :: m,f; ! N+1 colonnes et Q+1 lignes

!! Variables décentrées !!
double precision, dimension(N+2,Q+1) :: mbar; ! N+2 colonnes et Q+1 lignes
double precision, dimension(N+1,Q+2) :: fbar; ! N+1 colonnes et Q+2 lignes

!! Les valeurs aux frontières !!
double precision, dimension(N+1) :: f0; ! La densité initiale --> en bas
double precision, dimension(N+1) :: f1; ! la densité finale --> en haut 

!!! Initialisation des variables !!!

!! Grille centrée !!
do i=0,N
    GcX(i) = i/(1.0*N);
end do


!! Les fontières !!
call finitial(GcX,N,f0);
call ffinal(GcX,N,f1);

!! Ecrire les données dans les fichiers !!

open (unit=0,file="files/f0"); write (0,*), f0; close(0);
open (unit=1,file="files/f1"); write (1,*), f1; close(1);
open (unit=2,file="files/parameters"); write (2,*), N, Q; close(2);

!! Générer la cnstante dans la méthode de poisson !! 

call system('FreeFem++ -v 0 poisson_2d_constante.pde');

!! call check(alpha,beta,gamma,N,Q);
print *, "Lancement de l'algorithme";
call DR(alpha,beta,gamma,N,Q);
print *, "Fin de l'algorithme";

end program