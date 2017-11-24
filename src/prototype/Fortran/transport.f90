program transport
implicit none

include 'global.inc'

integer :: i,j,k;
double precision :: mu = 0.1 , sigma = 0.005 , minimum = 0.000001;

double precision, dimension(N+1) :: f0, f1;

! projection sur C !
double precision, dimension(2*Q+2*N+4+(N+1)*(Q+1),2*(Q+1)*(N+1)+Q+N+2) :: A, Atmp;
double precision, dimension(2*Q+2*N+4+(N+1)*(Q+1),2*Q+2*N+4+(N+1)*(Q+1)) :: delta;
double precision, dimension((N+1)*(Q+1)+2*(Q+1)+2*(N+1)) :: y;
double precision, dimension(2*(N+1)*(Q+1)+N+Q+2,2*(N+1)*(Q+1)+N+Q+2) :: IdPC;
double precision, dimension(2*(N+1)) :: tmpPC;

! Proximal of G2 !
double precision, dimension((N+1)*(Q+2)+(N+2)*(Q+1),(N+1)*(Q+2)+(N+2)*(Q+1)) :: tmpPG2,tmptmpPG2;

! Variables DR ! 
double precision, dimension((N+2)*(Q+1)+(N+1)*(Q+2)) :: wU0,wU1,zU;
double precision, dimension(2*(N+1)*(Q+1)) :: wV0,wV1,zV;

! variables export files 
character(10) :: charI;
character(200) :: cmd;
double precision, dimension(Q+1,N+1) :: tmpZV;

! variables construction des inverses !
integer INFO;
integer, dimension((N+1)*(Q+1)+2*(Q+1)+2*(N+1)) :: IPIVy;
integer, dimension(2*(Q+1)*(N+1)+Q+N+2) :: IPIVA;
integer, dimension((N+1)*(Q+2)+(N+2)*(Q+1)) :: IPIV;
double precision, dimension((N+1)*(Q+2)+(N+2)*(Q+1)) :: work;

external SGESV
external DGETRF
external DGETRI

!!! Construction des matrices opérateurs !!!

print *, 'Begin Initialisation :'

call gauss(mu,sigma,N,f0);
call gauss(mu+0.8,sigma,N,f1);

f0 = f0 + minimum;
f1 = f1 + minimum;

f0 = f0/sum(f0(:));
f1 = f1/sum(f1(:));

print *, 'Construction of matrix operators'

call boundary(B);
call interpolation(Interp);
call divergence(D);

print *, 'Construction of the projection on C'
! projection sur C ! 
A(1:(N+1)*(Q+1),1:2*(Q+1)*(N+1)+Q+N+2) = D;
A((N+1)*(Q+1)+1:2*Q+2*N+4+(N+1)*(Q+1),1:2*(Q+1)*(N+1)+Q+N+2) = B;
delta = matmul(A,transpose(A));

y(1:(N+1)*(Q+1)) = 0; ! condition de divergence
y((N+1)*(Q+1)+1:(N+1)*(Q+1)+2*(Q+1)) = 0; ! condition sur m
do i = 1,2*(N+1)
	if (modulo(i,2) .EQ. 1) then
		tmpPC(i) = f1((i+1)/2); ! on met f1
	else
		tmpPC(i) = f0(i/2); ! on met f0
	end if
end do
y((N+1)*(Q+1)+2*(Q+1)+1:(N+1)*(Q+1)+2*(Q+1)+2*(N+1)) = tmpPC;

do i = 1,2*Q+2*N+4+(N+1)*(Q+1)
	print *, delta(i,:), 'ENDL';
end do

!call SGESV(2*Q+2*N+4+(N+1)*(Q+1),(N+1)*(Q+1)+2*(Q+1)+2*(N+1),delta,2*Q+2*N+4+(N+1)*(Q+1),IPIVy,y,2*Q+2*N+4+(N+1)*(Q+1),INFO);
call SGESV(2*Q+2*N+4+(N+1)*(Q+1),1,delta,2*Q+2*N+4+(N+1)*(Q+1),IPIVy,y,(N+1)*(Q+1)+2*(Q+1)+2*(N+1),INFO);

Cst = matmul(transpose(A),y);

Atmp = A;
call SGESV(2*Q+2*N+4+(N+1)*(Q+1),2*(Q+1)*(N+1)+Q+N+2,delta,2*Q+2*N+4+(N+1)*(Q+1),IPIVA,A,2*Q+2*N+4+(N+1)*(Q+1),INFO);

IdPC = 0;
do i = 1,(N+1)*(Q+2)+(N+2)*(Q+1)
	IdPC(i,i) = 1;
end do
P = IdPC - matmul(transpose(Atmp),A);

print *, 'Construction of the projection on G2'

! projection sur G2 !
tmpPG2 = 0;
do i = 1,(N+1)*(Q+2)+(N+2)*(Q+1)
	tmpPG2(i,i) = 1;
end do 
tmptmpPG2 = tmpPG2 + matmul(transpose(Interp),Interp);
call DGETRF((N+1)*(Q+2)+(N+2)*(Q+1),(N+1)*(Q+2)+(N+2)*(Q+1),tmptmpPG2,(N+1)*(Q+2)+(N+2)*(Q+1),IPIV,INFO);
call DGETRI((N+1)*(Q+2)+(N+2)*(Q+1),tmptmpPG2,(N+1)*(Q+2)+(N+2)*(Q+1),IPIV,work,(N+1)*(Q+2)+(N+2)*(Q+1),INFO);
pG2 = tmptmpPG2;

!! CHECK !!

print *, 'tmpC :';
do i = 1,2*(N+1)*(Q+1)+N+Q+2
!print *, Cst(i) , 'ENDL';
end do
!call check(INFO);

stop 0
!! END CHECK !!

! Début ! 
print *, 'Begin DR algorithm'
wU0 = 0; wV0 = 0;
wU1 = 0; wV1 = 0;
zU = 0; zV = 0;

do i = 1,niter
	call proxG1(2*zU-wU0,2*zV-wV0,wU1,wV1);
	wU1 = wU0 + alpha*(wU1 - zU);
	wV1 = wV0 + alpha*(wV1 - zV);
	
	call proxG2(wU1,wV1,zU,zV);
	
	wU0 = wU1;
	wV0 = wV1;
	
	write(charI,'(I5.5)') i
	open(7,file='results/f_'//trim(charI)//'.dat'); 
	write(7,*) "# ", "X ", "T ", "Z";
	tmpZV = reshape(zV,(/Q+1,N+1/));
	do j = 1,N+1 
		do k = 1,Q+1 
			write(7,*) (j-1)/(1.0*N),(k-1)/(1.0*Q), tmpZV(k,j); 
		end do 
	end do
	close(7);
	write(cmd,*) 'set contour; set dgrid3d ', N, ',', Q, ';', 'splot "results/f_'//trim(charI), &
		'.dat" with lines; set title "Iteration Nb '//trim(charI),' "';
	open(8,file="plot.gnu"); write(8,*) cmd; close(8);
	call system("gnuplot -p plot.gnu");

end do


















end program
