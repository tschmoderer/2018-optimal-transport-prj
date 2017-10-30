subroutine polynome(P,x,m,f,g,N,Q)
    integer, intent(in) :: N,Q;
    double precision, dimension(N+1,Q+1), intent(in) :: x,m,f;
    double precision, dimension(N+1,Q+1), intent(out) :: P;
    double precision :: g; 
    P = -g*m**2 + (x-f)*((x+2*g)**2);
end subroutine
