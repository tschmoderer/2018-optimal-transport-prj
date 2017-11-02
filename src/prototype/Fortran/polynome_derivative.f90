subroutine polynome_derivative(dP,x,f,g,N,Q)
    integer, intent(in) :: N,Q;
    double precision, dimension(Q+1,N+1), intent(in) :: x,f;
    double precision, dimension(Q+1,N+1), intent(out) :: dP;
    double precision :: g;
    dP = 2.0*(x+2.0*g)*(x-f) + (x+2.0*g)**2;
end subroutine