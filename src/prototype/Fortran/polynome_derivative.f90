subroutine polynome_derivative(dP,x,f,g,N,Q)
    integer, intent(in) :: N,Q;
    double precision, dimension(N+1,Q+1), intent(in) :: x,f;
    double precision, dimension(N+1,Q+1), intent(out) :: dP;
    double precision :: g;
    dP = 2*x*(x-f) + (x+2*g)**2;
end subroutine