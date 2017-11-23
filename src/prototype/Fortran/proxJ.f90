!! Etant donné m,f 
!! renvoi l'opérateur de proximité sur J

subroutine proxJ(V,pV)
    implicit none
    include 'global.inc'
    
    double precision, dimension((N+1)*(Q+1)+(N+1)*(Q+1)), intent(in) :: V;
	double precision, dimension((N+1)*(Q+1)+(N+1)*(Q+1)), intent(out) :: pV;
	
	double precision, dimension((Q+1)*(N+1)) :: mt,ft,Pf,Pm;
	double precision, dimension((Q+1)*(N+1)) :: x0,x1,poly,dpoly,ddpoly;
	integer :: i,j,k;
    
    mt = V(1:(N+1)*(Q+1));
    ft = V((N+1)*(Q+1)+1:2*(N+1)*(Q+1));
    
    x0 = 1000;
    x1 = 0;
    k = 0;
    do while ((sum(abs(x0-x1)) .GT. 1e-5 ) .AND. (k .LT. 1500))
        x0 = x1;
        poly = (x0-ft)*((x0+g)**2) - 0.5*g*mt**2;
        dpoly = 2*(x0+g)*(x0-ft)+(x0+g)**2;
        ddpoly = 2*(3*x0+2*g-ft);
        x1 = x0-2*poly*dpoly/(2*dpoly**2-poly*ddpoly);
      !  call polynome(P,x0,m,f,g,N,Q);
       ! call polynome_derivative(dP,x0,f,g,N,Q)
      !  x1 = x0 - P/(1.0*dP);
        k = k+1;
    end do
	
	do i = 1,(Q+1)*(N+1)
		if (x1(i) .GT. 0) then 
			Pf(i) = x1(i);
			Pm(i) = (x1(i)*mt(i))/(1.0*(x1(i)+2.0*g));
		else 
			Pf(i) = 0;
			Pm(i) = 0;
		end if 
	end do
	
	pV(1:(N+1)*(Q+1)) = Pm;
    pV((N+1)*(Q+1)+1:2*(N+1)*(Q+1)) = Pf;
end subroutine
