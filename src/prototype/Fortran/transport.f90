program transport
implicit none

include 'global.inc'

! 

!!! Initialisation !!!
integer :: i,j;
double precision :: mu = 0.1 , sigma = 0.005;
double precision, dimension(N+1) :: gauss;
double precision, dimension(N+1) :: f0, f1;
!!! Construction des matrices opérateurs !!!

!normalise(f) = f/sum(f(:));
!print*, gauss(mu,sigma,N)


!! Matrice opérateur B !! 
do i = 1,2*Q+2*N+4
    do j = 1,2*(Q+1)*(N+1)+N+Q+2
    end do
end do 
end program