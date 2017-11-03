subroutine polynome_derivative(dP,x,f,g,N,Q)
    integer, intent(in) :: N,Q;
    double precision, dimension(Q+1,N+1), intent(in) :: x,f;
    double precision, dimension(Q+1,N+1), intent(out) :: dP;
    double precision, intent(in) :: g;
    
    dP = (x + 2.0*g)*(3.0*x + 2.0*g - 2.0*f) + 1e-16; ! evite division par 0
end subroutine
