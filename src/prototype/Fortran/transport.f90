program transport
implicit none

!! Les basiques -  indices de boucles!!
integer :: i,j,k;

!!! Initialisation !!!

!! L'espace temps !!
integer, parameter :: d = 1; ! Dimension du sytème
integer, parameter :: N = 80; ! Nb de points dans la direction x
integer, parameter :: Q = 30; ! Nb de points dans la direction t

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
real, dimension(N+1,Q+1) :: m,f; ! N+1 colonnes et Q+1 lignes

!! Variables décentrées !!
real, dimension(N+2,Q+1) :: mbar; ! N+2 colonnes et Q+1 lignes
real, dimension(N+1,Q+2) :: fbar; ! N+1 colonnes et Q+2 lignes

!! Les valeurs aux frontières !!
real, dimension(N+1) :: f0; ! La densité initiale --> en bas
real, dimension(N+1) :: f1; ! la densité finale --> en haut 

!! Variables dans la projection sur C !! 
real, dimension(N+2,Q+1) :: projCmbar;
real, dimension(N+1,Q+2) :: projCfbar;

!! variable prximité G2 !!
real, dimension(N+1,Q+1) :: mt,ft;
real, dimension(N+2,Q+1) :: mbart;
real, dimension(N+1,Q+2) :: fbart;

!! Le cout de la solution !!
real R;

!! Varaiable de la projection sur J !! 
real, dimension(N+1,Q+1) :: Pm,Pf;

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

close(unit=0);
close(unit=1);
close(unit=2);

!!! Zone de tests !!!

call RANDOM_NUMBER(m);
call RANDOM_NUMBER(f);
call RANDOM_NUMBER(mbar);
call RANDOM_NUMBER(fbar);
!! call check_adjoint(N,Q); !! test les opérateurs et les adjoints
!! call cost(R,m,f,N,Q); print *, R; !! test la fonction cout
!! call proxJ(Pm,Pf,m,f,gamma,N,Q); print *, Pm; !! test l'opérateur de proximité de la fonction cout
projCmbar = 0; projCfbar = 0;
call projC(projCmbar,projCfbar,mbar,fbar,N,Q); !! test projection sur C , inch allah ça marche
!! call proxG2(mbart,fbart,mt,ft,mbar,fbar,m,f,N,Q)
print *, "test"
!!! Fin Zone de tests !!!

end program