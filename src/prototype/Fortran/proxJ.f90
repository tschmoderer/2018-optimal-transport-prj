subroutine proxJ(Pm,Pf,m,f,g,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, intent(in) :: g;
    double precision, dimension(N+1,Q+1), intent(in) :: m,f; 
    double precision, dimension(N+1,Q+1), intent(out) :: Pm,Pf; 
    double precision, dimension(N+1,Q+1) :: x0,x1; 
    double precision, dimension(N+1,Q+1) :: P,dP;
    integer :: k;
    
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

    Pf = x1;
    Pm = (x1*m)/(1.0*(x1+2*g));
end subroutine

subroutine polynome(P,x,m,f,g,N,Q)
    integer, intent(in) :: N,Q;
    double precision, dimension(N+1,Q+1), intent(in) :: x,m,f;
    double precision, dimension(N+1,Q+1), intent(out) :: P;
    double precision :: g; 
    P = -g*m**2 + (x-f)*((x+2*g)**2);
end subroutine

subroutine polynome_derivative(dP,x,f,g,N,Q)
    integer, intent(in) :: N,Q;
    double precision, dimension(N+1,Q+1), intent(in) :: x,f;
    double precision, dimension(N+1,Q+1), intent(out) :: dP;
    double precision :: g;
    dP = 2*x*(x-f) + (x+2*g)**2;
end subroutine