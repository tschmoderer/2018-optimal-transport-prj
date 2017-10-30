program transport
implicit none

!! Les basiques -  indices de boucles!!
integer :: i,j,k;

!!! Initialisation !!!

!! L'espace temps !!
integer, parameter :: d = 1; ! Dimension du sytème
integer, parameter :: N = 10; ! Nb de points dans la direction x
integer, parameter :: Q = 10; ! Nb de points dans la direction t

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
double precision :: alpha = 1.0; ! Doit etre dans ]0,2[
double precision :: beta = 1.0; ! Doit etre dans [0,1]
double precision :: gamma = 0.5; ! Doit etre positif

!! Grille centrée !!
double precision, dimension(0:N) :: GcX;
double precision, dimension(0:Q) :: GcT;

!! Grille décentrée !!
double precision, dimension(-1:N) :: GsX;
double precision, dimension(-1:Q) :: GsT;

!! Variables centrées !!
double precision, dimension(N+1,Q+1) :: m,f; ! N+1 colonnes et Q+1 lignes

!! Variables décentrées !!
double precision, dimension(N+2,Q+1) :: mbar; ! N+2 colonnes et Q+1 lignes
double precision, dimension(N+1,Q+2) :: fbar; ! N+1 colonnes et Q+2 lignes

!! Les valeurs aux frontières !!
double precision, dimension(N+1) :: f0; ! La densité initiale --> en bas
double precision, dimension(N+1) :: f1; ! la densité finale --> en haut 

!! Variables dans la projection sur C !! 
double precision, dimension(N+2,Q+1) :: projCmbar;
double precision, dimension(N+1,Q+2) :: projCfbar;

!! variable prximité G2 !!
double precision, dimension(N+1,Q+1) :: mt,ft;
double precision, dimension(N+2,Q+1) :: mbart;
double precision, dimension(N+1,Q+2) :: fbart;

!! Le cout de la solution !!
double precision R;

!! Varaiable de la projection sur J !! 
double precision, dimension(N+1,Q+1) :: Pm,Pf;

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

call system('FreeFem++ poisson_2d_constante.pde');

close(0);
close(1);
close(2);

!!! Zone de tests !!!

call RANDOM_NUMBER(m);
call RANDOM_NUMBER(f);
call RANDOM_NUMBER(mbar);
call RANDOM_NUMBER(fbar);
!! call check_adjoint(N,Q); !! test les opérateurs et les adjoints
!! call cost(R,m,f,N,Q); print *, R; !! test la fonction cout
!! call proxJ(Pm,Pf,m,f,gamma,N,Q); print *, Pm; !! test l'opérateur de proximité de la fonction cout
!! call projC(projCmbar,projCfbar,mbar,fbar,N,Q); !! test projection sur C , inch allah ça marche
!! call proxG2(mbart,fbart,mt,ft,mbar,fbar,m,f,N,Q); !! projection sur la variable G2
!! call proxG1(mbar,fbar,m,f,gamma,N,Q); !! Proximité de G1
!!! Fin Zone de tests !!!


call DR(alpha,gamma,N,Q);

end program