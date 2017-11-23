!! Etant donnée des valeurs décentrées bar et fbar
!! Calcul leur divergence d

subroutine divergence(D)
    implicit none
    include 'global.inc'
 
	double precision, dimension((N+1)*(Q+1),(N+2)*(Q+1)) :: Dm; 
	double precision, dimension((N+1)*(Q+1),(N+1)*(Q+2)) :: Df; 
	double precision, dimension(Q+1,Q+2) :: tmp;
	integer :: i,j;
	
	do i = 1,(N+1)*(Q+1)
		do j = 1,(N+2)*(Q+1)
			if (i .EQ. j) then
				Dm(i,j) = -1;
			else if (j .EQ. i+Q+1) then 
				Dm(i,j) = 1;
			end if
		end do
	end do
	
	do i = 1,Q+1
		do j = 1,Q+2
			if (i .EQ. j) then
				tmp(i,j) = -1;
			else if (j .EQ. i+1) then 
				tmp(i,j) = 1;
			end if
		end do
	end do	
	
	do i = 1,N+1
		Df((i-1)*(Q+1)+1:i*(Q+1),(i-1)*(Q+2)+1:i*(Q+2)) = tmp;
	end do
		
	D(1:(N+1)*(Q+1),1:(N+2)*(Q+1)) = N*Dm;
	D(1:(N+1)*(Q+1),(N+2)*(Q+1)+1:2*(N+1)*(Q+1)+N+Q+2) = Q*Df;	
end subroutine
