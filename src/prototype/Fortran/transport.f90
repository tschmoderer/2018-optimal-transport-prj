program transport
implicit none

!!! Initialisation !!!
integer :: i;
!! L'espace temps !!
integer, parameter :: N = 30; ! Nb de points dans la direction x i.e Nb colonnes
integer, parameter :: Q = 30; ! Nb de points dans la direction t i.e Nb lignes

!   1
!   |
!   | 
!   |  Q+1
!   t  
!   |  pts
!   |
!   |     N+1 pts
!   0 ------ x ------ 1

!! Paramètres de la méthode !!
double precision :: alpha = 1.983; ! Doit etre dans ]0,2[
double precision :: beta = 1.0; ! Doit etre dans [0,1]
double precision :: gamma = 1.0; ! Doit etre positif

!! Grille centrée selon la variable  x !!
double precision, dimension(0:N) :: GcX = (/(i/(1.0*N),i=0,N)/);

!! Les valeurs aux frontières !!
double precision, dimension(N+1) :: f0; ! La densité initiale --> en bas
double precision, dimension(N+1) :: f1; ! la densité finale --> en haut 

!! Variables centrées !!
double precision, dimension(Q+1,N+1) :: m,f; ! Q+1 lignes, N+1 colonnes

!! Variables décentrées !! 
double precision, dimension(Q+1,N+2) :: mbar; ! Q+1 lignes, N+2 colonnes
double precision, dimension(Q+2,N+1) :: fbar; ! Q+2 lignes, N+1 colonnes

!!! Initialisation des variables !!!
!! Les fontières !!
call finitial(GcX,N,f0);
call ffinal(GcX,N,f1);

!! Ecrire les données dans les fichiers !!

open (unit=0,file="files/f0"); write (0,*), f0; close(0);
open (unit=1,file="files/f1"); write (1,*), f1; close(1);
open (unit=2,file="files/parameters"); write (2,*), N, Q; close(2);

!! Générer la cnstante dans la méthode de poisson !! 

call system('FreeFem++ -v 0 poisson_2d_constante.pde');

!! Lancement de l'algorithme de résolution !!
call check(alpha,beta,gamma,N,Q);
stop
call check(alpha,beta,gamma,N,Q);
stop;
print *, "Lancement de l'algorithme";
call DR(alpha,beta,gamma,N,Q);
print *, "Fin de l'algorithme";

end program