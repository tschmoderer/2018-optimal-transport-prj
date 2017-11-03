!! Etant donné m,f 
!! renvoi l'opérateur de proximité sur J

subroutine proxJ(Pm,Pf,m,f,g,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, intent(in) :: g;
    double precision, dimension(Q+1,N+1), intent(in) :: m,f; 
    double precision, dimension(Q+1,N+1), intent(out) :: Pm,Pf; 
    double precision, dimension(Q+1,N+1) :: x0,x1,P,dP; 
  
    integer :: i,j,k;
    
    x0 = 1000;
    x1 = 0;
    k = 0;
    do while ((sum(abs(x0-x1)) .GT. 1e-5 ) .AND. (k .LT. 50))
        x0 = x1;
        call polynome(P,x0,m,f,g,N,Q);
        call polynome_derivative(dP,x0,f,g,N,Q)
        x1 = x0 - P/(1.0*dP);
        k = k+1;
    end do
	
	do i = 1,Q+1
		do j = 1,N+1
			if (x1(i,j) .GT. 0) then 
				Pf(i,j) = x1(i,j);
   				Pm(i,j) = (x1(i,j)*m(i,j))/(1.0*(x1(i,j)+2.0*g));
			else 
				Pf(i,j) = 0;
				Pm(i,j) = 0;
			end if 
		end do
	end do
end subroutine
