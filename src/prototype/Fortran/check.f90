subroutine check(INFO)
	implicit none
	include 'global.inc'
	
	integer :: INFO;
	integer :: i,j;
	print *, 'Dans CHECK : '; 
	!print *, 'f0 : ', f0;
	!print *, 'f1 : ', f1;

	print *, 'B : ', B;
	!print *, 'Interp : ', Interp;
	!print *, 'D : ', D;

	!print *, 'y : ', y;
	!print *, 'A : ', A;
	!print *, 'delta : ', delta;
	!print *, 'Cst : ', Cst;
	!print *, 'P : ', P;

	!print *, 'PG2 : ', pG2;
	!print *, 'PG2 : ', pG2;
	
	print *, 'P : '; 
	do i =1,2*(N+1)*(Q+1)+N+Q+2
		print *, P(i,:);
		print *, 'ENDL'
	end do
end subroutine
