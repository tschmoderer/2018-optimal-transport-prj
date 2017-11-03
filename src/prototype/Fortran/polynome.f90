subroutine polynome(P,x,m,f,g,N,Q)
    integer, intent(in) :: N,Q;
    double precision, dimension(Q+1,N+1), intent(in) :: x,m,f;
    double precision, dimension(Q+1,N+1), intent(out) :: P;
    double precision, intent(in) :: g; 
    
    P = -g*m*m + (x-f)*((x+2.0*g)*(x+2.0*g));
end subroutine
